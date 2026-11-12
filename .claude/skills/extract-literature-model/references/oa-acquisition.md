# Open-access PDF acquisition ladder

Detailed source-by-source acquisition strategy used by `scripts/acquire-paper.R`. The script is the entry point — read this reference when the script reports a failure and you need to understand which step failed, or when adding a new source to the ladder.

## Lead-PDF source ladder

For each missing or invalid file (file does not exist, is < 10 KB, or whose first 4 bytes are not `%PDF`), attempt OA acquisition before giving up. Try sources in this order, stopping at the first that yields a valid `%PDF`-headed download of ≥ 10 KB:

1. **CrossRef link array.** `curl -sS "https://api.crossref.org/works/<DOI>" -A "literature-acquisition (mailto:wdenney@humanpredictions.com)"` → inspect the `message.link[]` array; download any URL whose `content-type` contains `pdf` or whose URL ends in `.pdf`.
2. **Unpaywall.** `curl -sS "https://api.unpaywall.org/v2/<DOI>?email=wdenney@humanpredictions.com"` → if `best_oa_location.url_for_pdf` is set, download it.
3. **Europe PMC.** `curl -sS "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:<DOI>&format=json"` → if any result has a `pmcid` field, download `https://europepmc.org/articles/<PMCID>?pdf=render`. Most effective for Wiley, Elsevier, and Springer papers in PMC after embargo.
4. **NCBI PMC ID converter.** `curl -sS "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids=<DOI>&format=json&tool=literature-acquisition&email=wdenney@humanpredictions.com"` → download `https://www.ncbi.nlm.nih.gov/pmc/articles/<PMCID>/pdf/`.
5. **Publisher-specific patterns** (only for known-OA outlets):
   - BMC: `https://<journal>.biomedcentral.com/track/pdf/<DOI>`
   - Frontiers: `https://www.frontiersin.org/articles/<DOI>/pdf`
   - Nature OA / Springer Nature OA: `https://www.nature.com/articles/<id>.pdf` (where `<id>` is the trailing segment of the DOI)
   - Hindawi: `https://onlinelibrary.wiley.com/doi/pdf/<DOI>` (post-Wiley acquisition)
   - Dovepress: `https://www.dovepress.com/getfile.php?fileID=<file_id>`
   - Static Springer Nature media (for Nature Comm SI files): `https://static-content.springer.com/esm/art%3A<doi-percent-encoded>/MediaObjects/<id>_MOESM1_ESM.pdf` (on 403, use the `media.springernature.com` form under Supplements below)

## Validation rules

- Size ≥ 10 KB AND first 4 bytes equal `%PDF`. Wiley/Elsevier non-OA endpoints typically return ~5 KB HTML challenge pages — reject (delete and try the next source). Never extract from a file whose head bytes are not `%PDF`.
- The size floor and `%PDF` check are for PDFs. Check any other supplement or model file with `file -b --mime-type <file>`: reject `text/html` or `application/xhtml+xml` (an error or challenge page; `file` sees through a BOM, leading whitespace and letter case). Otherwise the type must match the extension: `text/*` for code and text (`.txt`, `.ctl`, `.mod`, `.m`, `.R`, `.csv`); `application/zip` or an `application/vnd.openxmlformats-officedocument.*` type, i.e. the ZIP signature `PK`, for `.zip`, `.docx`, `.xlsx`. No size floor: model code can be well under 10 KB.
- **Title-content sanity check** for the lead PDF: after a successful download, `pdftotext -l 1 <file> -` and confirm the title / first-author line approximately matches what the task block expects. CrossRef DOIs occasionally point to a different paper than the task's metadata claims (publishers' DOI sequence is not always continuous, e.g. `10.1038/aps.2014.120` was a Zhang ROR review when the task expected Lu 2015 tacrolimus). If the downloaded PDF's title disagrees with the task expectation, don't extract from it — sidecar-ask the operator with the actual-vs-expected metadata.

## Reasonable-attempts cap

Don't loop forever. The 5-source ladder above with one retry per source is the cap; after that, record and defer (below). If a source returns a 5xx error, retry once after 5 seconds; otherwise treat as failed and move on.

## When to give up: record and defer

If all 5 source paths have been tried and none produced a valid PDF, this is not a stop-and-ask: add the paper to `needs_acquisition.jsonl` and defer the task (SKILL.md Phase 1 Step 0). List what was tried in the report: CrossRef link array (n URLs), Unpaywall (`oa_status`, PDF URL or none), Europe PMC and NCBI PMC (PMCID or not found), publisher URLs, and what each returned. An HTTP failure (403, 404, HTML challenge page) goes in the row's `last_transport_error`, never in `reason`: a 403 is a scraper block, not a paywall, and a `reason` marks the paper unobtainable. A PDF whose title disagrees with the task is the identity question under Validation rules, not an acquisition failure.

## Supplements and errata

The same ladder applies, with publisher-specific tweaks:

- **Wiley supplements**: try `https://onlinelibrary.wiley.com/action/downloadSupplement?doi=<DOI>&file=<filename>`. Filenames often follow the pattern `<journal-prefix><articleid>-sup-<NNNN>-<descriptor>.pdf`; if the exact filename isn't known, fetch the article landing page and grep for `downloadSupplement` URLs.
- **Springer Nature SI**: `https://static-content.springer.com/esm/art%3A<doi-percent-encoded>/MediaObjects/<id>_MOESM1_ESM.pdf`. Increment `MOESM2_ESM`, `MOESM3_ESM`, etc. for additional supplements. When `static-content.springer.com` returns 403, use `https://media.springernature.com/original/springer-static/esm/art%3A<DOI with / encoded as %2F>/MediaObjects/<filename>`. ESM files are often not PDFs (model code as `.txt`, `.docx`, `.xlsx`, `.zip`): take `<filename>` (`<journal-id>_<year>_<article-no>_MOESM<n>_ESM.<ext>`) from the article page's ESM links rather than assuming `.pdf`.
- **BMC**: PMC mirror often has supplements at `https://www.ncbi.nlm.nih.gov/pmc/articles/<PMCID>/bin/<filename>`.
- **Frontiers**: typically inline at the article landing page; fetch the page and grep for `Image_*.pdf` / `Table_*.docx` / `Presentation_*.pdf` in the HTML.
- **Errata**: search PubMed (`https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=<title>+AND+erratum&retmode=json`) for the original PMID + erratum link; download the erratum PDF via the same OA ladder.

If a supplement is unobtainable but the lead PDF is on disk, handle it as SKILL.md Phase 1 Step 0 says for a supplement that cannot be obtained.
