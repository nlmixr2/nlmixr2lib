Ramakrishnan_2019_bortezomib_ncih929 <- function() {
  description <- "QSP. In vitro (NCI-H929 human multiple myeloma cell line). Pharmacodynamic systems model of bortezomib signaling from Ramakrishnan 2019: bortezomib in the culture medium (20 nM, first-order loss with a 144 h half-life) linearly stimulates the synthesis of pJNK (direct), pNFkB (two transit compartments), pAKT (one transit compartment) and p21, and a network of indirect-response turnover equations for pmTOR, pBAD (one transit compartment), Bcl-2, caspase 8 and caspase 9 (relative expression, 1 = untreated control at 0 h; no pp53 node) drives myeloma cell growth (pAKT x pmTOR / p21) and death (caspase 8 x caspase 9). Equations assembled from the MM.1S equations plus the published NCI-H929 differences (Data S1, Figure 1d). Naive-pooled maximum-likelihood fit (ADAPT 5); no between-well variability. Sibling cell-line models: Ramakrishnan_2019_bortezomib_u266, _rpmi8226, _mm1s."
  reference <- paste(
    "Ramakrishnan V, Mager DE. (2019). Pharmacodynamic Models of Differential",
    "Bortezomib Signaling Across Several Cell Lines of Multiple Myeloma.",
    "CPT Pharmacometrics Syst Pharmacol 8(3):146-157. doi:10.1002/psp4.12358.",
    "Equations: Supplementary Data S1 MM.1S Eqs. 20-34 modified as described in",
    "the Data S1 NCI-H929 section and Figure 1d; parameter values from Table 1;",
    "ADAPT 5 model code (Data S2, RPMI8226 version) for the variance model and",
    "the separate untreated-control growth equation."
  )
  vignette <- "Ramakrishnan_2019_bortezomib_myeloma"

  # Paper-mechanistic states with no canonical analogue in
  # inst/references/compartment-names.md: bortezomib in the culture medium,
  # the nine measured signaling proteins (relative expression), the transit
  # compartments that delay pNFkB / pAKT / pBAD (pattern below), and the
  # treated and untreated relative myeloma cell numbers.
  paper_specific_compartments <- c(
    "media",
    "pjnk",
    "pnfkb",
    "pakt",
    "p21",
    "pmtor",
    "cas8",
    "pbad",
    "bcl2",
    "cas9",
    "mm_cells",
    "mm_cells_control"
  )
  paper_specific_compartment_pattern <- "^(pnfkb|pakt|pbad)_transit[0-9]+$"

  dosing <- "media"

  units <- list(
    time = "h",
    dosing = "nM",
    concentration = "nM for media; protein states and cell numbers are relative to the untreated control at 0 h (unitless)"
  )

  compartmentData <- list(
    media = list(
      analyte = "bortezomib in the culture medium",
      units = "nM",
      specimen = "not applicable",
      verified = TRUE
    ),
    pjnk = list(
      analyte = "phosphorylated JNK (pJNK)",
      units = "unitless (relative expression)",
      specimen = "not applicable",
      verified = TRUE
    ),
    pnfkb_transit1 = list(
      analyte = "pNFkB signal transit compartment 1",
      units = "unitless (relative)",
      specimen = "not applicable",
      verified = TRUE
    ),
    pnfkb_transit2 = list(
      analyte = "pNFkB signal transit compartment 2",
      units = "unitless (relative)",
      specimen = "not applicable",
      verified = TRUE
    ),
    pnfkb = list(
      analyte = "phosphorylated NFkB (pNFkB)",
      units = "unitless (relative expression)",
      specimen = "not applicable",
      verified = TRUE
    ),
    pakt_transit1 = list(
      analyte = "pAKT signal transit compartment 1",
      units = "unitless (relative)",
      specimen = "not applicable",
      verified = TRUE
    ),
    pakt = list(
      analyte = "phosphorylated AKT (pAKT)",
      units = "unitless (relative expression)",
      specimen = "not applicable",
      verified = TRUE
    ),
    p21 = list(analyte = "p21", units = "unitless (relative expression)", specimen = "not applicable", verified = TRUE),
    pmtor = list(
      analyte = "phosphorylated mTOR (pmTOR)",
      units = "unitless (relative expression)",
      specimen = "not applicable",
      verified = TRUE
    ),
    cas8 = list(
      analyte = "caspase 8",
      units = "unitless (relative expression)",
      specimen = "not applicable",
      verified = TRUE
    ),
    pbad_transit1 = list(
      analyte = "Ser112-phosphorylated BAD, first (pJNK-inhibited) compartment",
      units = "unitless (relative)",
      specimen = "not applicable",
      verified = TRUE
    ),
    pbad = list(
      analyte = "Ser112-phosphorylated BAD (pBAD)",
      units = "unitless (relative expression)",
      specimen = "not applicable",
      verified = TRUE
    ),
    bcl2 = list(
      analyte = "Bcl-2",
      units = "unitless (relative expression)",
      specimen = "not applicable",
      verified = TRUE
    ),
    cas9 = list(
      analyte = "caspase 9",
      units = "unitless (relative expression)",
      specimen = "not applicable",
      verified = TRUE
    ),
    mm_cells = list(
      analyte = "viable NCI-H929 cells under bortezomib",
      units = "unitless (relative to untreated control at 0 h)",
      specimen = "not applicable",
      verified = TRUE
    ),
    mm_cells_control = list(
      analyte = "viable untreated NCI-H929 cells",
      units = "unitless (relative to untreated control at 0 h)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  # The paper's protein and cell outputs each share ONE estimated variance
  # coefficient; rxode2 cannot tie one residual parameter to several
  # endpoints, so each endpoint carries its own copy of the shared value.
  paper_specific_residual_sds <- c(
    "propSd_pjnk",
    "propSd_pnfkb",
    "propSd_pakt",
    "propSd_p21",
    "propSd_pmtor",
    "propSd_cas8",
    "propSd_pbad",
    "propSd_bcl2",
    "propSd_cas9",
    "propSd_mm_cells",
    "propSd_mm_cells_control"
  )

  covariateData <- list()

  population <- list(
    species = "in vitro (NCI-H929 human multiple myeloma cell line)",
    n_subjects = "Cell-culture time courses (naive pooled): 9 modelled signaling proteins by Western blot (pp53 measured but unchanged and not modelled) and viable-cell counts, untreated control and 20 nM bortezomib",
    n_studies = 1L,
    disease_state = "Multiple myeloma ('more sensitive' cell line, wild-type TP53)",
    dose_range = "Single bortezomib concentration of 20 nM in the culture medium",
    regions = "University at Buffalo, NY (USA)",
    notes = "Data from Ramakrishnan and Mager 2018 (J Pharmacol Exp Ther 365:734-751). Protein data truncated at 24 h for model development; cell proliferation followed to 96 h and fitted without holding any protein at its 24 h value."
  )

  ini({
    # Bortezomib loss from the medium -- main text Methods: 'mono-exponential
    # decay with a half-life of 144 hours'; Data S2 code KDEG = 0.0048135208.
    lkdeg_media <- fixed(log(log(2) / 144))
    label("First-order bortezomib degradation rate in culture medium (1/h)") # Methods: t1/2 = 144 h (Chudasama 2015)

    # Direct drug targets -- Table 1 'NCI-H929'.
    lkdeg_pjnk <- log(0.1394)
    label("pJNK turnover rate constant (1/h)") # Table 1 NCI-H929 'k tpJNK' = 0.1394 (no pJNK transit compartment in this cell line)
    slope_pjnk <- 0.2635
    label("Linear stimulation of pJNK synthesis by bortezomib SmpJNK (1/nM)") # Table 1 NCI-H929 'S mpJNK' = 0.2635
    lktr_pnfkb <- log(0.1123)
    label("pNFkB transit rate constant ktpNFkB (1/h)") # Table 1 NCI-H929 'k tpNFkB' = 0.1123
    slope_pnfkb <- 3.808
    label("Linear stimulation of pNFkB synthesis by bortezomib SmpNFkB (1/nM)") # Table 1 NCI-H929 'S mpNFkB' = 3.808
    lktr_pakt <- log(0.1633)
    label("pAKT transit rate constant ktpAKT (1/h)") # Table 1 NCI-H929 'k tpAKT' = 0.1633
    slope_pakt <- 0.4531
    label("Linear stimulation of pAKT synthesis by bortezomib SmpAKT (1/nM)") # Table 1 NCI-H929 'S mpAKT' = 0.4531
    lkdeg_p21 <- log(0.07804)
    label("p21 degradation rate constant kdp21 (1/h)") # Table 1 NCI-H929 'k dp21' = 0.07804
    slope_p21 <- 3.627
    label("Linear stimulation of p21 synthesis by bortezomib Smp21 (1/nM)") # Table 1 NCI-H929 'S mp21' = 3.627

    # Downstream proteins -- Table 1 'NCI-H929'.
    lkdeg_pmtor <- log(0.005520)
    label("pmTOR degradation rate constant kdpmTOR (1/h)") # Table 1 NCI-H929 'k dpmTOR' = 0.005520
    lkdeg_cas8 <- log(0.2219)
    label("Caspase 8 degradation rate constant kdCas8 (1/h)") # Table 1 NCI-H929 'k dCas8' = 0.2219
    gamma_cas8 <- 3.163
    label("Power coefficient on pJNK stimulation of caspase 8 synthesis (unitless)") # Table 1 NCI-H929 'gamma Cas8' = 3.163
    lktr_pbad <- log(0.3298)
    label("pBAD transit rate constant (1/h)") # Table 1 NCI-H929 'k dpBAD' (described as 'pBAD transit rate constant') = 0.3298
    gamma_pbad <- 1.103
    label("Power coefficient on pJNK-stimulated loss of pBAD (unitless)") # Table 1 NCI-H929 'gamma pBAD' = 1.103
    lkdeg_bcl2 <- log(0.1105)
    label("Bcl-2 degradation rate constant kdBcl2 (1/h)") # Table 1 NCI-H929 'k dBcl-2' = 0.1105
    gamma_bcl2 <- 6.167
    label("Power coefficient on pBAD stimulation of Bcl-2 synthesis (unitless)") # Table 1 NCI-H929 'gamma Bcl-2' = 6.167
    lkdeg_cas9 <- log(1.544)
    label("Caspase 9 degradation rate constant kdCas9 (1/h)") # Table 1 NCI-H929 'k dCas9' = 1.544

    # Myeloma cell growth and death -- main text Eq. 1.
    lkgrow <- log(0.02278)
    label("Myeloma cell growth rate constant kg (1/h)") # Table 1 NCI-H929 'k g' = 0.02278
    lkdeath <- log(0.004006)
    label("Myeloma cell death rate constant kd (1/h)") # Table 1 NCI-H929 'k d' = 0.004006

    # Residual error -- Data S1 Eq. 2, Var = (sigma1 + sigma2 * Y)^2. Table 1
    # reports only the slope sigma2 for each output type (proteins, cells).
    propSd_pjnk <- 0.7509
    label("Proportional residual error, pJNK (fraction)") # Table 1 NCI-H929 'sigma2 Protein' = 0.7509 (shared by all proteins)
    propSd_pnfkb <- 0.7509
    label("Proportional residual error, pNFkB (fraction)") # Table 1 NCI-H929 'sigma2 Protein' = 0.7509
    propSd_pakt <- 0.7509
    label("Proportional residual error, pAKT (fraction)") # Table 1 NCI-H929 'sigma2 Protein' = 0.7509
    propSd_p21 <- 0.7509
    label("Proportional residual error, p21 (fraction)") # Table 1 NCI-H929 'sigma2 Protein' = 0.7509
    propSd_pmtor <- 0.7509
    label("Proportional residual error, pmTOR (fraction)") # Table 1 NCI-H929 'sigma2 Protein' = 0.7509
    propSd_cas8 <- 0.7509
    label("Proportional residual error, caspase 8 (fraction)") # Table 1 NCI-H929 'sigma2 Protein' = 0.7509
    propSd_pbad <- 0.7509
    label("Proportional residual error, pBAD (fraction)") # Table 1 NCI-H929 'sigma2 Protein' = 0.7509
    propSd_bcl2 <- 0.7509
    label("Proportional residual error, Bcl-2 (fraction)") # Table 1 NCI-H929 'sigma2 Protein' = 0.7509
    propSd_cas9 <- 0.7509
    label("Proportional residual error, caspase 9 (fraction)") # Table 1 NCI-H929 'sigma2 Protein' = 0.7509
    propSd_mm_cells <- 0.2987
    label("Proportional residual error, treated myeloma cells (fraction)") # Table 1 NCI-H929 'sigma2 MMcell' = 0.2987
    propSd_mm_cells_control <- 0.2987
    label("Proportional residual error, untreated myeloma cells (fraction)") # Table 1 NCI-H929 'sigma2 MMcell' = 0.2987 (shared)
  })

  model({
    kdeg_media <- exp(lkdeg_media)
    kdeg_pjnk <- exp(lkdeg_pjnk)
    ktr_pnfkb <- exp(lktr_pnfkb)
    ktr_pakt <- exp(lktr_pakt)
    kdeg_p21 <- exp(lkdeg_p21)
    kdeg_pmtor <- exp(lkdeg_pmtor)
    kdeg_cas8 <- exp(lkdeg_cas8)
    ktr_pbad <- exp(lktr_pbad)
    kdeg_bcl2 <- exp(lkdeg_bcl2)
    kdeg_cas9 <- exp(lkdeg_cas9)
    kgrow <- exp(lkgrow)
    kdeath <- exp(lkdeath)

    # Data S1 Eq. 1: every protein starts at 1 (untreated control, 0 h), so
    # each zero-order synthesis rate equals its first-order loss rate.

    # Bortezomib in the medium (Cb, nM); a 20 nM dose into media starts treatment.
    d/dt(media) <- -kdeg_media * media

    # Data S1 NCI-H929: 'no transit compartment' for pJNK (Figure 1d).
    d/dt(pjnk) <- kdeg_pjnk * (1 + slope_pjnk * media) - kdeg_pjnk * pjnk

    # pNFkB and pAKT as MM.1S Eqs. 23-27 (Figure 1d: pNFkB 2TC, pAKT 1TC).
    d/dt(pnfkb_transit1) <- ktr_pnfkb * (1 + slope_pnfkb * media) - ktr_pnfkb * pnfkb_transit1
    d/dt(pnfkb_transit2) <- ktr_pnfkb * (pnfkb_transit1 - pnfkb_transit2)
    d/dt(pnfkb) <- ktr_pnfkb * (pnfkb_transit2 - pnfkb)
    d/dt(pakt_transit1) <- ktr_pakt * (1 + slope_pakt * media) - ktr_pakt * pakt_transit1
    d/dt(pakt) <- ktr_pakt * (pakt_transit1 - pakt)

    # MM.1S Eq. 28 with no power on bortezomib (no gamma p21 in Table 1) and
    # the pAKT inhibition of p21 'retained for NCI-H929' (Data S1), written
    # as the division by pAKT used for that interaction in U266 Eq. 13.
    d/dt(p21) <- kdeg_p21 / pakt * (1 + slope_p21 * media) - kdeg_p21 * p21

    # MM.1S Eq. 29: pmTOR stimulated by pAKT.
    d/dt(pmtor) <- kdeg_pmtor * pakt - kdeg_pmtor * pmtor

    # MM.1S Eq. 30: caspase 8 stimulated by pJNK^gamma, inhibited by pNFkB.
    d/dt(cas8) <- kdeg_cas8 / pnfkb * pjnk^gamma_cas8 - kdeg_cas8 * cas8

    # pBAD with one transit compartment (Figure 1d 'pBAD 1TC') and a power
    # on the pJNK-stimulated loss (Data S1), as in RPMI8226 code XP(10)-XP(13).
    d/dt(pbad_transit1) <- ktr_pbad - ktr_pbad * pjnk^gamma_pbad * pbad_transit1
    d/dt(pbad) <- ktr_pbad * (pbad_transit1 - pbad)

    # MM.1S Eq. 32 plus the pNFkB stimulation of Bcl-2 'retained for
    # NCI-H929' (Data S1), multiplicative as in U266 Eq. 17.
    d/dt(bcl2) <- kdeg_bcl2 * pbad^gamma_bcl2 * pnfkb - kdeg_bcl2 * bcl2

    # MM.1S Eq. 33: caspase 9 stimulated by caspase 8, inhibited by Bcl-2 and pNFkB.
    d/dt(cas9) <- kdeg_cas9 / (bcl2 * pnfkb) * cas8 - kdeg_cas9 * cas9

    # Main text Eq. 1: treated cells. Data S1: NCI-H929 needed no 24 h hold
    # of the protein drivers.
    d/dt(mm_cells) <- kgrow / p21 * pakt * pmtor * mm_cells -
      kdeath * cas8 * cas9 * mm_cells

    # Untreated control cells grow first-order at kg (Data S2 code XP(17)).
    d/dt(mm_cells_control) <- kgrow * mm_cells_control

    media(0) <- 0
    pjnk(0) <- 1
    pnfkb_transit1(0) <- 1
    pnfkb_transit2(0) <- 1
    pnfkb(0) <- 1
    pakt_transit1(0) <- 1
    pakt(0) <- 1
    p21(0) <- 1
    pmtor(0) <- 1
    cas8(0) <- 1
    pbad_transit1(0) <- 1
    pbad(0) <- 1
    bcl2(0) <- 1
    cas9(0) <- 1
    mm_cells(0) <- 1
    mm_cells_control(0) <- 1

    pjnk ~ prop(propSd_pjnk)
    pnfkb ~ prop(propSd_pnfkb)
    pakt ~ prop(propSd_pakt)
    p21 ~ prop(propSd_p21)
    pmtor ~ prop(propSd_pmtor)
    cas8 ~ prop(propSd_cas8)
    pbad ~ prop(propSd_pbad)
    bcl2 ~ prop(propSd_bcl2)
    cas9 ~ prop(propSd_cas9)
    mm_cells ~ prop(propSd_mm_cells)
    mm_cells_control ~ prop(propSd_mm_cells_control)
  })
}
