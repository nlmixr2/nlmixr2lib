#!/usr/bin/env python3
"""Handle an upstream-popPK dependency for the current extraction task.

When the dispatched agent (running ``/extract-literature-model``)
discovers that the current paper's PK is fixed to a separate
publication NOT on disk, this script is the default action per
Phase 1 Step 6:

  1. Try to acquire the upstream PDF via the OA-PDF ladder
     (``acquire-paper.R``).
  2. If acquired -> place at ``papers/PMID_<id>/PMID_<id>.pdf`` and
     drop a ``trim_queue`` marker so the trim daemon processes it.
  3. If NOT acquired -> write
     ``papers/PMID_<id>/PMID_<id>_needs_acquisition.flag`` with clear
     "operator drops PDF here" instructions.
  4. Either way: create a new task YAML for the upstream paper in
     ``todo/`` (next free ``NNN-<author>_<year>_<drug>.yaml`` slot)
     that re-uses this same ``/extract-literature-model`` skill.
  5. Edit the current task's YAML to add
     ``depends_on: [<upstream_task_id>]``.
  6. Write a deferral report and exit cleanly.

The current task does NOT inline-extract the upstream model. It
pauses (via the dispatchable depends_on chain in
runner/orchestrator.py) until the upstream task commits, then
re-dispatches with the upstream PK parameters available on disk.

USAGE
-----
    python3 handle-upstream-dependency.py \\
        --queue-dir /home/.../nlmixr2lib_ingestion \\
        --current-task-id 1006-clozapine-pd \\
        --upstream-pmid 19349931 \\
        --upstream-citation 'Ng W et al. 2009, Ther Drug Monit 31:360-6' \\
        --upstream-drug clozapine

EXIT CODES
----------
    0  success (upstream queued; current task depends_on edited)
    1  current task YAML not found
    2  bad arguments
    3  filesystem error during file placement
    4  task-YAML edit failed (depends_on injection)
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

REPO_ROOT_ENV = "NLMIXR2LIB_REPO"  # optional override
DEFAULT_ACQUIRE_RSCRIPT = ".claude/skills/extract-literature-model/scripts/acquire-paper.R"


def _slugify_drug(drug: str) -> str:
    """Lowercase, strip non-alphanumeric; preserves camelCase / mixed-case where
    sensible by simply lowercasing. Used for the upstream task filename."""
    return re.sub(r"[^A-Za-z0-9]+", "_", drug).strip("_").lower()


def _slugify_author(author: str) -> str:
    """Author tokens: keep hyphens, strip diacritics best-effort, lowercase."""
    # Common diacritic stripping (covers ñ, é, etc.). We rely on simple ASCII
    # for ID slugs; the full citation lives in the YAML notes.
    import unicodedata

    nfkd = unicodedata.normalize("NFKD", author)
    ascii_only = "".join(c for c in nfkd if not unicodedata.combining(c))
    return re.sub(r"[^A-Za-z0-9\-]+", "", ascii_only)


def _parse_citation(citation: str) -> tuple[str, int | None]:
    """Extract first-author surname and 4-digit year from a citation string.

    Heuristic — robust enough for the common citation shapes in this corpus:
       "Ng W et al. 2009, Ther Drug Monit 31:360-6"  -> ("Ng", 2009)
       "Hamberg AK, Dahl ML, ... 2007;81(4):529-38"   -> ("Hamberg", 2007)
       "<bad string>"                                  -> ("upstream", None)
    """
    # First "word"-looking token = author surname
    m = re.match(r"\s*([A-Za-z][A-Za-z\-']{1,40})", citation)
    author = m.group(1) if m else "upstream"
    # First 4-digit year (1900-2099)
    y = re.search(r"\b(19|20)\d{2}\b", citation)
    year = int(y.group(0)) if y else None
    return author, year


def _next_task_id(todo_dir: Path, author: str, year: int | None, drug: str) -> str:
    """Compute the next free NNN-<author>_<year>_<drug> task ID.

    Scans every YAML in todo_dir + done/prefilter*/ for existing 3-digit
    prefixes; picks max+1. If `year` is None we fall back to "0000".
    """
    used: set[int] = set()
    for path in todo_dir.parent.glob("todo/*.yaml"):
        m = re.match(r"^(\d{3,4})-", path.stem)
        if m:
            used.add(int(m.group(1)))
    # Also consult done/ subdirectories
    for done_sub in (todo_dir.parent / "done").glob("*/"):
        for path in done_sub.glob("*.yaml"):
            m = re.match(r"^(\d{3,4})-", path.stem)
            if m:
                used.add(int(m.group(1)))
    next_n = max(used, default=0) + 1
    year_str = str(year) if year is not None else "0000"
    return f"{next_n:03d}-{author.lower()}_{year_str}_{_slugify_drug(drug)}"


def _try_acquire_pdf(
    *,
    queue_dir: Path,
    pmid: str | None,
    doi: str | None,
    out_pdf: Path,
    acquire_script: Path | None,
) -> bool:
    """Attempt to acquire the upstream PDF via acquire-paper.R.

    Returns True if a PDF landed at ``out_pdf``; False otherwise (paywall,
    script missing, network error, etc.).
    """
    if acquire_script is None or not acquire_script.exists():
        return False
    if doi is None and pmid is None:
        return False

    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    args = ["Rscript", str(acquire_script), "--out", str(out_pdf), "--retries", "1"]
    if doi is not None:
        args += ["--doi", doi]
    elif pmid is not None:
        args += ["--pmid", pmid]
    try:
        res = subprocess.run(args, capture_output=True, text=True, timeout=180)
    except (subprocess.SubprocessError, FileNotFoundError) as exc:
        sys.stderr.write(f"acquire-paper failed: {exc}\n")
        return False
    if res.returncode == 0 and out_pdf.exists() and out_pdf.stat().st_size > 1024:
        return True
    return False


def _build_task_yaml(
    *,
    task_id: str,
    upstream_drug: str,
    upstream_pmid: str | None,
    upstream_doi: str | None,
    upstream_citation: str,
    pdf_path: Path,
    queue_dir: Path,
    downstream_task_id: str,
) -> str:
    """Render the upstream task YAML body using the corpus's standard template."""
    expected_pdf = pdf_path
    expected_xml = pdf_path.with_suffix(".pmc.xml") if pdf_path.suffix == ".pdf" else None
    drug_label = upstream_drug.strip() or "(unknown)"
    pmid_line = f"PMID:          {upstream_pmid}" if upstream_pmid else "PMID:          (unknown)"
    doi_line = f"DOI:           {upstream_doi}" if upstream_doi else "DOI:           (unknown)"
    expected_filename = f"{task_id.split('-', 1)[-1].split('_', 1)[0].capitalize()}_{task_id.split('_')[1] if '_' in task_id else 'YYYY'}_{_slugify_drug(upstream_drug)}.R"

    worktree = f"/home/bill/github/nlmixr2/nlmixr2lib/.claude/worktrees/{task_id}"
    pmid_str = upstream_pmid or "(unknown)"
    doi_str = upstream_doi or "(unknown)"

    prompt_body = (
        "Use the /extract-literature-model skill to add a population PK\n"
        "model to the nlmixr2lib R package. This task was queued as a\n"
        "follow-up to the downstream task `" + downstream_task_id + "`;\n"
        "that task's PK is fixed to this paper. The downstream task has\n"
        "depends_on: [" + task_id + "] and will resume after this task\n"
        "commits.\n"
        "\n"
        "Paper metadata (from upstream citation):\n"
        "  Citation:  " + upstream_citation + "\n"
        "  " + pmid_line + "\n"
        "  " + doi_line + "\n"
        "  Drug:      " + drug_label + "\n"
        "  Target:    inst/modeldb/specificDrugs/\n"
        "\n"
        "Source files on disk (authoritative for parameter values + ODE form):\n"
        "  PDF: " + str(expected_pdf) + "\n"
        "\n"
        "If the PDF is missing at the path above, the pre-dispatch hook\n"
        "will DEFER this task — the operator drops the PDF and deletes\n"
        "the matching `_needs_acquisition.flag` companion when ready.\n"
        "\n"
        "Hard constraints:\n"
        "  - Follow the six phases of /extract-literature-model exactly,\n"
        "    including Phase 1 errata search and Phase 4 verification.\n"
        "  - Parameter VALUES and EQUATIONS must come from the source on\n"
        "    disk. Never substitute training-data values.\n"
        "  - Before committing: run `Rscript -e 'devtools::load_all(\\\".\\\");\n"
        "    devtools::check(error_on = \\\"error\\\", args = \\\"--no-build-vignettes\\\")'`\n"
        "    and expect 0 errors.\n"
        "  - Regenerate the model database with `nlmixr2lib:::buildModelDb()`.\n"
        "  - Push the branch. DO NOT run `gh pr create`. Print suggested PR\n"
        "    title + body at the end.\n"
        "\n"
        "Deliverables: write a one-paragraph summary to\n"
        "  " + str(queue_dir / "reports") + "/" + task_id + ".md\n"
        "covering PK structure, source locations, stop-and-ask decisions,\n"
        "deviations, branch name, and commit SHA.\n"
    )

    yaml = (
        "schema_version: 2\n"
        f"id: {task_id}\n"
        f"title: Extract upstream popPK model — {drug_label} ({upstream_citation[:60]})\n"
        "model: claude-opus-4-7\n"
        "effort: high\n"
        "priority: normal\n"
        "depends_on: []\n"
        "allowed_tools:\n"
        "- Read\n"
        "- Edit\n"
        "- Write\n"
        "- Bash\n"
        "- Grep\n"
        "- Glob\n"
        "- WebFetch\n"
        "- AskUserQuestion\n"
        f"working_dir: {worktree}\n"
        f"prompt: {json.dumps(prompt_body)}\n"
    )
    return yaml


def _write_needs_acquisition_flag(
    *,
    expected_pdf: Path,
    upstream_pmid: str | None,
    upstream_doi: str | None,
    upstream_citation: str,
    upstream_drug: str,
    downstream_task_id: str,
) -> Path:
    """Write the `_needs_acquisition.flag` companion that the trim daemon
    and pre-dispatch hook honour as "wait for operator-supplied PDF"."""
    flag_path = expected_pdf.with_name(expected_pdf.stem + "_needs_acquisition.flag")
    flag_path.parent.mkdir(parents=True, exist_ok=True)
    now = dt.datetime.now(dt.UTC).isoformat(timespec="seconds")
    body = (
        "schema_version: 1\n"
        f"created_at: {now}\n"
        f"flagged_by: extract-literature-model upstream-dependency helper\n"
        f"downstream_task: {downstream_task_id}\n"
        "reason: |\n"
        "  The OA-PDF ladder could not find an open-access copy.\n"
        "  Operator must acquire this paper via institutional access or\n"
        "  interlibrary loan and drop the PDF at the path below.\n"
        "paper_metadata:\n"
        f"  pmid:     {upstream_pmid or '(unknown)'}\n"
        f"  doi:      {upstream_doi or '(unknown)'}\n"
        f"  drug:     {upstream_drug}\n"
        f"  citation: {json.dumps(upstream_citation)}\n"
        "acquisition_instructions: |\n"
        f"  1. Acquire the upstream paper PDF.\n"
        f"  2. Save it at:\n"
        f"       {expected_pdf}\n"
        f"  3. Delete this `_needs_acquisition.flag` so the trim daemon\n"
        f"     picks it up.\n"
        f"  4. The trim_queue is monitored by `_scripts/trim_daemon.py`.\n"
        f"  5. Once trimmed, the upstream extraction task will dispatch\n"
        f"     automatically. The downstream task `{downstream_task_id}`\n"
        f"     will re-dispatch after the upstream commits, since it has\n"
        f"     depends_on pointing at the upstream task.\n"
    )
    flag_path.write_text(body, encoding="utf-8")
    return flag_path


def _add_depends_on(yaml_path: Path, upstream_task_id: str, downstream_task_id: str) -> None:
    """Edit the current task's YAML to add ``depends_on: [<upstream_task_id>]``.

    Replaces an existing ``depends_on: []`` line in-place; if the field already
    has entries, appends to the existing list. If the YAML can't be parsed,
    raises a clean error rather than corrupting the file.
    """
    if not yaml_path.exists():
        raise FileNotFoundError(f"task YAML not found: {yaml_path}")
    text = yaml_path.read_text(encoding="utf-8")

    # Use a minimal text-based edit to avoid breaking the prompt's quoting.
    new_block = (
        "depends_on:\n"
        f"- {upstream_task_id}\n"
        f"# {dt.datetime.now(dt.UTC).date()}: upstream-popPK dependency auto-handled —\n"
        "# the current task's PK is fixed to the paper above; on re-dispatch the\n"
        "# downstream task imports its PK parameters with provenance comments.\n"
        "# See .claude/skills/extract-literature-model/scripts/handle-upstream-dependency.py\n"
    )

    if "depends_on: []" in text:
        new_text = text.replace("depends_on: []", new_block.rstrip(), 1)
    else:
        # Add the upstream to an existing list. Find the depends_on: block
        # and inject the new entry at its head.
        m = re.search(r"^(depends_on:\s*\n((?:-\s.+\n)+))", text, re.MULTILINE)
        if m is None:
            raise RuntimeError(
                f"cannot find depends_on field in {yaml_path}; "
                "manual edit required to inject upstream task id"
            )
        old_block = m.group(0)
        injected = (
            "depends_on:\n"
            f"- {upstream_task_id}\n"
            + "".join(m.group(2).splitlines(keepends=True))
            + f"# {dt.datetime.now(dt.UTC).date()}: upstream dependency auto-added\n"
        )
        new_text = text.replace(old_block, injected, 1)

    yaml_path.write_text(new_text, encoding="utf-8")


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--queue-dir", required=True, type=Path,
                    help="Runner queue directory (contains todo/, papers/, .claude_task_runner/)")
    ap.add_argument("--current-task-id", required=True,
                    help="Task ID of the DOWNSTREAM task whose PK is fixed to the upstream paper.")
    ap.add_argument("--upstream-pmid", default=None,
                    help="Upstream PubMed ID (preferred).")
    ap.add_argument("--upstream-doi", default=None,
                    help="Upstream DOI (used if --upstream-pmid not provided).")
    ap.add_argument("--upstream-citation", required=True,
                    help="Verbatim citation string (used for metadata and slug parsing).")
    ap.add_argument("--upstream-drug", required=True,
                    help="Generic drug name (e.g. clozapine, warfarin).")
    ap.add_argument("--upstream-task-id-override", default=None,
                    help="Override the auto-computed upstream task id (advanced).")
    ap.add_argument("--no-acquire", action="store_true",
                    help="Skip the OA-PDF ladder; jump directly to needs_acquisition.flag.")
    ap.add_argument("--acquire-script", default=None, type=Path,
                    help="Path to acquire-paper.R (default: discover relative to this script).")
    args = ap.parse_args(argv)

    if not args.upstream_pmid and not args.upstream_doi:
        sys.stderr.write("error: --upstream-pmid or --upstream-doi required\n")
        return 2

    queue_dir: Path = args.queue_dir.resolve()
    if not queue_dir.is_dir():
        sys.stderr.write(f"error: queue dir not found: {queue_dir}\n")
        return 2

    # Verify the current task YAML exists.
    current_yaml = queue_dir / "todo" / f"{args.current_task_id}.yaml"
    if not current_yaml.exists():
        sys.stderr.write(f"error: current task YAML not found: {current_yaml}\n")
        return 1

    # Resolve acquire-paper.R relative to this script's location unless overridden.
    if args.acquire_script is not None:
        acquire_script: Path | None = args.acquire_script
    else:
        acquire_script = (
            Path(__file__).resolve().parent / "acquire-paper.R"
            if (Path(__file__).resolve().parent / "acquire-paper.R").exists()
            else None
        )

    # Compute the upstream task id + expected on-disk paths.
    author, year = _parse_citation(args.upstream_citation)
    author_slug = _slugify_author(author)
    if args.upstream_task_id_override:
        upstream_task_id = args.upstream_task_id_override
    else:
        upstream_task_id = _next_task_id(
            todo_dir=queue_dir / "todo",
            author=author_slug,
            year=year,
            drug=args.upstream_drug,
        )

    pmid_key = args.upstream_pmid or "noPMID"
    papers_subdir = queue_dir / "papers" / f"PMID_{pmid_key}"
    expected_pdf = papers_subdir / f"PMID_{pmid_key}.pdf"

    # Step 1: try OA-PDF acquisition unless suppressed.
    acquired = False
    if not args.no_acquire:
        acquired = _try_acquire_pdf(
            queue_dir=queue_dir,
            pmid=args.upstream_pmid,
            doi=args.upstream_doi,
            out_pdf=expected_pdf,
            acquire_script=acquire_script,
        )

    # Step 2: place trim_queue marker (if acquired) or needs_acquisition.flag.
    trim_queue_dir = queue_dir / "_scripts" / "trim_queue"
    if acquired:
        trim_queue_dir.mkdir(parents=True, exist_ok=True)
        marker = trim_queue_dir / (expected_pdf.name + ".req")
        marker.touch()
        sys.stderr.write(f"acquired upstream PDF -> {expected_pdf}\n")
    else:
        flag_path = _write_needs_acquisition_flag(
            expected_pdf=expected_pdf,
            upstream_pmid=args.upstream_pmid,
            upstream_doi=args.upstream_doi,
            upstream_citation=args.upstream_citation,
            upstream_drug=args.upstream_drug,
            downstream_task_id=args.current_task_id,
        )
        sys.stderr.write(f"upstream not acquired; placed flag -> {flag_path}\n")

        # Append to corpus-wide needs_acquisition.jsonl for audit visibility.
        jsonl_path = queue_dir / "needs_acquisition.jsonl"
        try:
            with jsonl_path.open("a", encoding="utf-8") as f:
                entry = {
                    "pmid": args.upstream_pmid,
                    "doi": args.upstream_doi,
                    "drug": args.upstream_drug,
                    "downstream_task": args.current_task_id,
                    "upstream_task": upstream_task_id,
                    "expected_path": str(expected_pdf),
                    "flagged_at": dt.datetime.now(dt.UTC).isoformat(timespec="seconds"),
                    "source": "extract-literature-model upstream-dependency helper",
                }
                f.write(json.dumps(entry) + "\n")
        except OSError as exc:
            sys.stderr.write(f"warn: could not append to needs_acquisition.jsonl: {exc}\n")

    # Step 3: write the upstream task YAML.
    upstream_yaml = queue_dir / "todo" / f"{upstream_task_id}.yaml"
    if upstream_yaml.exists():
        sys.stderr.write(
            f"warn: upstream task YAML already exists at {upstream_yaml}; "
            f"skipping creation (idempotent re-run)\n"
        )
    else:
        upstream_yaml.write_text(
            _build_task_yaml(
                task_id=upstream_task_id,
                upstream_drug=args.upstream_drug,
                upstream_pmid=args.upstream_pmid,
                upstream_doi=args.upstream_doi,
                upstream_citation=args.upstream_citation,
                pdf_path=expected_pdf,
                queue_dir=queue_dir,
                downstream_task_id=args.current_task_id,
            ),
            encoding="utf-8",
        )
        sys.stderr.write(f"wrote upstream task YAML -> {upstream_yaml}\n")

    # Step 4: edit the current task YAML to add depends_on.
    try:
        _add_depends_on(current_yaml, upstream_task_id, args.current_task_id)
    except Exception as exc:
        sys.stderr.write(f"error: failed to add depends_on to {current_yaml}: {exc}\n")
        return 4
    sys.stderr.write(f"added depends_on:[{upstream_task_id}] to {current_yaml}\n")

    # Step 5: emit a JSON summary to stdout for callers to parse.
    summary = {
        "downstream_task_id": args.current_task_id,
        "upstream_task_id": upstream_task_id,
        "upstream_pmid": args.upstream_pmid,
        "upstream_doi": args.upstream_doi,
        "upstream_drug": args.upstream_drug,
        "acquired": acquired,
        "expected_pdf_path": str(expected_pdf),
        "upstream_task_yaml": str(upstream_yaml),
        "downstream_task_yaml": str(current_yaml),
    }
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
