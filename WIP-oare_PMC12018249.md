# WIP — Xu 2025 checkpoint-BsAb QSP (PMC12018249) — paused on sidecar request-001

## State
- `inst/modeldb/specificDrugs/Xu_2025_checkpoint_bsab_qsp.R` — COMPLETE draft, 765 lines,
  59 ODE states, 152 ini() parameters, ASCII-clean, parses under rxode2.
- Vignette NOT yet written. NEWS.md not yet updated. Nothing committed.

## Verified so far
- PK units defect RESOLVED: Table S1 prints CL1 = 12 and Q1 = 1.725E-2 both as "L/h";
  literal reading gives a 0.15 s half-life. Reading CL1 as uL/h and Q1 as mL/h gives a
  simulated terminal t1/2 of 140.9 h, matching the ~142 h read off Supplementary Fig S3A
  (MSLN/CD3 3 mg/kg, blue "S" general-parameter curve). Encoded + commented in ini().
- lsoda / liblsoda / dop853 agree to all printed digits on the untreated trajectory
  (100 -> 3233.7 mm3 over 25 d). Use lsoda (2.6 s); dop853 takes 43 s.
- Drug binding cascade works: 10 mg/kg anti-TIGIT/PD-L1 ip biw drops
  PD1_tumor_PDL1_per_cell 1.26e-14 -> 7.10e-15 and TIGIT_tumor_CD155_per_cell
  1.52e-10 -> 2.42e-11; Erk_a 0.300 -> 0.337; NFkB 106 -> 129.
- km1 (1e-14), km2 (1e-10), km4 (4.1e-6), km21 (1000), km24/25/27 (1e-8) each land on the
  simulated operating point of the species they gate -> the cell/volume scaling is right.

## BLOCKER (sidecar request-001)
`n3` couples the immune module to tumour killing. Published range 4e-6 to 1.06e-5
(midpoint 6.512e-6) gives an immune multiplier of 1.66e-10 against a baseline of 1,
so TGI = 0.0% at every dose. Sweep (anti-TIGIT/PD-L1 ip biw, TGI at day 25):
    n3=6.51e-6 -> 0.0% at 1/3/10/18 mg/kg   (control V25 = 3234 mm3)
    n3=1e4     -> 3.8 / 7.0 / 11.1 / 13.0%  (control 2972)
    n3=3e4     -> 17.0 / 28.0 / 38.6 / 42.5% (control 2240)
    n3=1e5     -> 31.1 / 43.5 / 53.4 / 56.8% (control 622)
    n3=3e5     -> 31.0 / 45.9 / 58.9 / 64.0% (control 166)
Structure reproduces the paper's monotone-saturating dose-response only at n3 ~ 3e4-3e5.
Also internally inconsistent: the IFN-gamma Hill constant in killing_tumor_vivo is 1e-3
but max achievable IFNg is ~1e-6 (kf_IFNg_p=1.22e-7, kel_IFNg=0.022, promoter frac <= 2/11).

## Next steps once answered
1. Set n3 per the operator's choice; re-run the gate.
2. Write `vignettes/articles/Xu_2025_checkpoint_bsab_qsp.Rmd` (PK vs Fig S3A; binding
   cascade; untreated growth vs the Fig 5A day-25 scenario; TGI dose-response vs Fig 5B/6D;
   Errata: PK units, n3, IFN-g Hill constant, tumor_MHC_per_cell alias, drug_11 trimer typo,
   km ranges, SHP099 PK, dimer_CTLA4_10 = 0, mAb arms out of scope).
   No PKNCA: this is a mechanistic model (see references/endogenous-validation.md).
3. buildModelDb(), lint, ASCII gate, render gate, NEWS.md, push, PR text.
4. Report to /home/bill/gitlab/nlmixr2lib_ingestion/refharvest/reports/oare_PMC12018249.md
