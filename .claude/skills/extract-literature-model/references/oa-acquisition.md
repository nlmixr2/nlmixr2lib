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
   - Static Springer Nature media (for Nature Comm SI files): `https://static-content.springer.com/esm/art%3A<doi-percent-encoded>/MediaObjects/<id>_MOESM1_ESM.pdf`

## Validation rules

- Size ≥ 10 KB AND first 4 bytes equal `%PDF`. Wiley/Elsevier non-OA endpoints typically return ~5 KB HTML challenge pages — reject (delete and try the next source). Never extract from a file whose head bytes are not `%PDF`.
- **Title-content sanity check** for the lead PDF: after a successful download, `pdftotext -l 1 <file> -` and confirm the title / first-author line approximately matches what the task block expects. CrossRef DOIs occasionally point to a different paper than the task's metadata claims (publishers' DOI sequence is not always continuous, e.g. `10.1038/aps.2014.120` was a Zhang ROR review when the task expected Lu 2015 tacrolimus). If the downloaded PDF's title disagrees with the task expectation, don't extract from it — sidecar-ask the operator with the actual-vs-expected metadata.

## Reasonable-attempts cap

Don't loop forever. The 5-source ladder above with one retry per source is the cap; after that, sidecar. If a source returns a 5xx error, retry once after 5 seconds; otherwise treat as failed and move on.

## When to give up and sidecar

If all 5 source paths have been tried and none produced a valid PDF whose title matches the task expectation, write a sidecar request describing what was attempted:

> Lead PDF for <paper> not on disk and OA acquisition failed. Tried: CrossRef link array (n URLs), Unpaywall (oa_status=<status>, pdf_url=<url-or-none>), Europe PMC (PMCID=<id-or-not-found>), NCBI PMC (PMCID=<id-or-not-found>), publisher landing (<list-of-tried-URLs>). Each returned <reason: HTML challenge page / 404 / not-PDF / title-mismatch>. Options: (A) operator drops the PDF on disk and re-dispatches, (B) operator emails corresponding author, (C) skip this task. Which applies?

## Supplements and errata

The same ladder applies, with publisher-specific tweaks:

- **Wiley supplements**: try `https://onlinelibrary.wiley.com/action/downloadSupplement?doi=<DOI>&file=<filename>`. Filenames often follow the pattern `<journal-prefix><articleid>-sup-<NNNN>-<descriptor>.pdf`; if the exact filename isn't known, fetch the article landing page and grep for `downloadSupplement` URLs.
- **Springer Nature SI**: `https://static-content.springer.com/esm/art%3A<doi-percent-encoded>/MediaObjects/<id>_MOESM1_ESM.pdf`. Increment `MOESM2_ESM`, `MOESM3_ESM`, etc. for additional supplements.
- **BMC**: PMC mirror often has supplements at `https://www.ncbi.nlm.nih.gov/pmc/articles/<PMCID>/bin/<filename>`.
- **Frontiers**: typically inline at the article landing page; fetch the page and grep for `Image_*.pdf` / `Table_*.docx` / `Presentation_*.pdf` in the HTML.
- **Errata**: search PubMed (`https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=<title>+AND+erratum&retmode=json`) for the original PMID + erratum link; download the erratum PDF via the same OA ladder.

If a supplement is unobtainable but the lead PDF is on disk, decide based on whether the missing supplement contains parameter values you need: if it does, sidecar-ask before extracting (see Phase 4 missing-parameter pathway); if it does not, document the gap in the vignette Errata and proceed.
