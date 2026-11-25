Ramakrishnan_2019_bortezomib_u266 <- function() {
  description <- "QSP. In vitro (U266 human multiple myeloma cell line). Pharmacodynamic systems model of bortezomib signaling from Ramakrishnan 2019: bortezomib in the culture medium (20 nM, first-order loss with a 144 h half-life) linearly stimulates the synthesis of pJNK, pNFkB and pAKT through transit compartments, and a network of indirect-response turnover equations for pmTOR, pp53, p21, pBAD, Bcl-2, caspase 8 and caspase 9 (relative expression, 1 = untreated control at 0 h) drives myeloma cell growth (pAKT x pmTOR / p21) and death (caspase 8 x caspase 9). Proteins driving the cell equation are held at their 24 h values thereafter, as in the final U266 fit. Naive-pooled maximum-likelihood fit (ADAPT 5); no between-well variability. Sibling cell-line models: Ramakrishnan_2019_bortezomib_rpmi8226, _mm1s, _ncih929."
  reference <- paste(
    "Ramakrishnan V, Mager DE. (2019). Pharmacodynamic Models of Differential",
    "Bortezomib Signaling Across Several Cell Lines of Multiple Myeloma.",
    "CPT Pharmacometrics Syst Pharmacol 8(3):146-157. doi:10.1002/psp4.12358.",
    "Equations from Supplementary Data S1 Eqs. 1-19 (U266); parameter values",
    "from Table 1; ADAPT 5 model code (Data S2, RPMI8226 version) for the",
    "variance model and the separate untreated-control growth equation."
  )
  vignette <- "Ramakrishnan_2019_bortezomib_myeloma"

  # Paper-mechanistic states with no canonical analogue in
  # inst/references/compartment-names.md: bortezomib in the culture medium,
  # the ten measured signaling proteins (relative expression), the transit
  # compartments that delay pJNK / pNFkB / pAKT / pBAD (pattern below), the
  # treated and untreated relative myeloma cell numbers, and the five
  # '_held' copies of the cell-equation drivers that stop changing at 24 h.
  paper_specific_compartments <- c(
    "media",
    "pjnk",
    "pnfkb",
    "pakt",
    "pmtor",
    "cas8",
    "pp53",
    "p21",
    "pbad",
    "bcl2",
    "cas9",
    "mm_cells",
    "mm_cells_control",
    "pakt_held",
    "pmtor_held",
    "p21_held",
    "cas8_held",
    "cas9_held"
  )
  paper_specific_compartment_pattern <- "^(pjnk|pnfkb|pakt|pbad)_transit[0-9]+$"

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
    pjnk_transit1 = list(
      analyte = "pJNK signal transit compartment 1",
      units = "unitless (relative)",
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
    pp53 = list(
      analyte = "phosphorylated p53 (pp53)",
      units = "unitless (relative expression)",
      specimen = "not applicable",
      verified = TRUE
    ),
    p21 = list(analyte = "p21", units = "unitless (relative expression)", specimen = "not applicable", verified = TRUE),
    pbad_transit1 = list(
      analyte = "Ser112-phosphorylated BAD, first (pJNK-inhibited) compartment",
      units = "unitless (relative)",
      specimen = "not applicable",
      verified = TRUE
    ),
    pbad_transit2 = list(
      analyte = "pBAD signal transit compartment 2",
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
      analyte = "viable U266 cells under bortezomib",
      units = "unitless (relative to untreated control at 0 h)",
      specimen = "not applicable",
      verified = TRUE
    ),
    mm_cells_control = list(
      analyte = "viable untreated U266 cells",
      units = "unitless (relative to untreated control at 0 h)",
      specimen = "not applicable",
      verified = TRUE
    ),
    pakt_held = list(
      analyte = "pAKT as seen by the cell equation (frozen after 24 h)",
      units = "unitless (relative expression)",
      specimen = "not applicable",
      verified = TRUE
    ),
    pmtor_held = list(
      analyte = "pmTOR as seen by the cell equation (frozen after 24 h)",
      units = "unitless (relative expression)",
      specimen = "not applicable",
      verified = TRUE
    ),
    p21_held = list(
      analyte = "p21 as seen by the cell equation (frozen after 24 h)",
      units = "unitless (relative expression)",
      specimen = "not applicable",
      verified = TRUE
    ),
    cas8_held = list(
      analyte = "caspase 8 as seen by the cell equation (frozen after 24 h)",
      units = "unitless (relative expression)",
      specimen = "not applicable",
      verified = TRUE
    ),
    cas9_held = list(
      analyte = "caspase 9 as seen by the cell equation (frozen after 24 h)",
      units = "unitless (relative expression)",
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
    "propSd_pmtor",
    "propSd_cas8",
    "propSd_pp53",
    "propSd_p21",
    "propSd_pbad",
    "propSd_bcl2",
    "propSd_cas9",
    "propSd_mm_cells",
    "propSd_mm_cells_control"
  )

  covariateData <- list()

  population <- list(
    species = "in vitro (U266 human multiple myeloma cell line)",
    n_subjects = "Cell-culture time courses (naive pooled): 10 signaling proteins by Western blot and viable-cell counts, untreated control and 20 nM bortezomib",
    n_studies = 1L,
    disease_state = "Multiple myeloma ('less sensitive' cell line, mutated TP53)",
    dose_range = "Single bortezomib concentration of 20 nM in the culture medium",
    regions = "University at Buffalo, NY (USA)",
    notes = "Data from Ramakrishnan and Mager 2018 (J Pharmacol Exp Ther 365:734-751). Protein data truncated at 24 h for model development (complete cell death by 24-30 h); cell proliferation followed to 96 h."
  )

  ini({
    # Bortezomib loss from the medium -- main text Methods: 'mono-exponential
    # decay with a half-life of 144 hours'; Data S2 code KDEG = 0.0048135208.
    lkdeg_media <- fixed(log(log(2) / 144))
    label("First-order bortezomib degradation rate in culture medium (1/h)") # Methods: t1/2 = 144 h (Chudasama 2015)

    # Upstream drug targets -- Data S1 Eqs. 3-9, Table 1 'U266'.
    lktr_pjnk <- log(0.07654)
    label("pJNK transit rate constant ktpJNK (1/h)") # Table 1 U266 'k tpJNK' = 0.07654
    slope_pjnk <- 2.757
    label("Linear stimulation of pJNK synthesis by bortezomib SmpJNK (1/nM)") # Table 1 U266 'S mpJNK' = 2.757
    lktr_pnfkb <- log(0.03528)
    label("pNFkB transit rate constant ktpNFkB (1/h)") # Table 1 U266 'k tpNFkB' = 0.03528
    slope_pnfkb <- 3.984
    label("Linear stimulation of pNFkB synthesis by bortezomib SmpNFkB (1/nM)") # Table 1 U266 'S mpNFkB' = 3.984
    lktr_pakt <- log(0.08517)
    label("pAKT transit rate constant ktpAKT (1/h)") # Table 1 U266 'k tpAKT' = 0.08517
    slope_pakt <- 0.5079
    label("Linear stimulation of pAKT synthesis by bortezomib SmpAKT (1/nM)") # Table 1 U266 'S mpAKT' = 0.5079

    # Downstream proteins -- Data S1 Eqs. 10-18, Table 1 'U266'.
    lkdeg_pmtor <- log(0.01434)
    label("pmTOR degradation rate constant kdpmTOR (1/h)") # Table 1 U266 'k dpmTOR' = 0.01434
    lkdeg_cas8 <- log(0.3287)
    label("Caspase 8 degradation rate constant kdCas8 (1/h)") # Table 1 U266 'k dCas8' = 0.3287
    lkdeg_pp53 <- log(0.3222)
    label("pp53 degradation rate constant kdpp53 (1/h)") # Table 1 U266 'k dpp53' = 0.3222
    lkdeg_p21 <- log(4.709)
    label("p21 degradation rate constant kdp21 (1/h)") # Table 1 U266 'k dp21' = 4.709 (289 %CV)
    gamma_p21 <- fixed(2.5)
    label("Power coefficient on pp53 stimulation of p21 synthesis (unitless)") # Table 1 U266 'gamma p21' = 2.500*
    lktr_pbad <- log(0.08116)
    label("pBAD transit rate constant ktpBAD (1/h)") # Table 1 U266 'k tpBAD' = 0.08116
    lkdeg_bcl2 <- log(0.01094)
    label("Bcl-2 degradation rate constant kdBcl2 (1/h)") # Table 1 U266 'k dBcl-2' = 0.01094
    gamma_bcl2 <- fixed(5)
    label("Power coefficient on pBAD stimulation of Bcl-2 synthesis (unitless)") # Table 1 U266 'gamma Bcl-2' = 5.000*
    lkdeg_cas9 <- log(0.9148)
    label("Caspase 9 degradation rate constant kdCas9 (1/h)") # Table 1 U266 'k dCas9' = 0.9148
    gamma_cas9 <- 1.278
    label("Power coefficient on caspase 8 stimulation of caspase 9 synthesis (unitless)") # Table 1 U266 'gamma Cas9' = 1.278

    # Myeloma cell growth and death -- Data S1 Eq. 19 (main text Eq. 1).
    lkgrow <- log(0.01882)
    label("Myeloma cell growth rate constant kg (1/h)") # Table 1 U266 'k g' = 0.01882
    lkdeath <- log(0.003776)
    label("Myeloma cell death rate constant kd (1/h)") # Table 1 U266 'k d' = 0.003776

    # Residual error -- Data S1 Eq. 2, Var = (sigma1 + sigma2 * Y)^2. Table 1
    # reports only the slope sigma2 for each output type (proteins, cells).
    propSd_pjnk <- 0.4724
    label("Proportional residual error, pJNK (fraction)") # Table 1 U266 'sigma2 Protein' = 0.4724 (shared by all proteins)
    propSd_pnfkb <- 0.4724
    label("Proportional residual error, pNFkB (fraction)") # Table 1 U266 'sigma2 Protein' = 0.4724
    propSd_pakt <- 0.4724
    label("Proportional residual error, pAKT (fraction)") # Table 1 U266 'sigma2 Protein' = 0.4724
    propSd_pmtor <- 0.4724
    label("Proportional residual error, pmTOR (fraction)") # Table 1 U266 'sigma2 Protein' = 0.4724
    propSd_cas8 <- 0.4724
    label("Proportional residual error, caspase 8 (fraction)") # Table 1 U266 'sigma2 Protein' = 0.4724
    propSd_pp53 <- 0.4724
    label("Proportional residual error, pp53 (fraction)") # Table 1 U266 'sigma2 Protein' = 0.4724
    propSd_p21 <- 0.4724
    label("Proportional residual error, p21 (fraction)") # Table 1 U266 'sigma2 Protein' = 0.4724
    propSd_pbad <- 0.4724
    label("Proportional residual error, pBAD (fraction)") # Table 1 U266 'sigma2 Protein' = 0.4724
    propSd_bcl2 <- 0.4724
    label("Proportional residual error, Bcl-2 (fraction)") # Table 1 U266 'sigma2 Protein' = 0.4724
    propSd_cas9 <- 0.4724
    label("Proportional residual error, caspase 9 (fraction)") # Table 1 U266 'sigma2 Protein' = 0.4724
    propSd_mm_cells <- 0.3645
    label("Proportional residual error, treated myeloma cells (fraction)") # Table 1 U266 'sigma2 MMcell' = 0.3645
    propSd_mm_cells_control <- 0.3645
    label("Proportional residual error, untreated myeloma cells (fraction)") # Table 1 U266 'sigma2 MMcell' = 0.3645 (shared)
  })

  model({
    kdeg_media <- exp(lkdeg_media)
    ktr_pjnk <- exp(lktr_pjnk)
    ktr_pnfkb <- exp(lktr_pnfkb)
    ktr_pakt <- exp(lktr_pakt)
    kdeg_pmtor <- exp(lkdeg_pmtor)
    kdeg_cas8 <- exp(lkdeg_cas8)
    kdeg_pp53 <- exp(lkdeg_pp53)
    kdeg_p21 <- exp(lkdeg_p21)
    ktr_pbad <- exp(lktr_pbad)
    kdeg_bcl2 <- exp(lkdeg_bcl2)
    kdeg_cas9 <- exp(lkdeg_cas9)
    kgrow <- exp(lkgrow)
    kdeath <- exp(lkdeath)

    # Data S1 Eq. 1: every protein starts at 1 (untreated control, 0 h), so
    # each zero-order synthesis rate equals its first-order loss rate.

    # Bortezomib in the medium (Cb, nM); a 20 nM dose into media starts treatment.
    d/dt(media) <- -kdeg_media * media

    # Eqs. 3-9: bortezomib stimulates the first transit compartment.
    d/dt(pjnk_transit1) <- ktr_pjnk * (1 + slope_pjnk * media) - ktr_pjnk * pjnk_transit1
    d/dt(pjnk) <- ktr_pjnk * (pjnk_transit1 - pjnk)
    d/dt(pnfkb_transit1) <- ktr_pnfkb * (1 + slope_pnfkb * media) - ktr_pnfkb * pnfkb_transit1
    d/dt(pnfkb_transit2) <- ktr_pnfkb * (pnfkb_transit1 - pnfkb_transit2)
    d/dt(pnfkb) <- ktr_pnfkb * (pnfkb_transit2 - pnfkb)
    d/dt(pakt_transit1) <- ktr_pakt * (1 + slope_pakt * media) - ktr_pakt * pakt_transit1
    d/dt(pakt) <- ktr_pakt * (pakt_transit1 - pakt)

    # Eq. 10: pmTOR stimulated by pAKT.
    r_pmtor <- kdeg_pmtor * pakt - kdeg_pmtor * pmtor
    d/dt(pmtor) <- r_pmtor

    # Eq. 11: caspase 8 stimulated by pJNK, inhibited by pNFkB.
    r_cas8 <- kdeg_cas8 / pnfkb * pjnk - kdeg_cas8 * cas8
    d/dt(cas8) <- r_cas8

    # Eqs. 12-13: pp53 stimulated by pJNK; p21 by pp53^gamma; both inhibited by pAKT.
    d/dt(pp53) <- kdeg_pp53 / pakt * pjnk - kdeg_pp53 * pp53
    r_p21 <- kdeg_p21 / pakt * pp53^gamma_p21 - kdeg_p21 * p21
    d/dt(p21) <- r_p21

    # Eqs. 14-16: pJNK accelerates loss of the first pBAD compartment.
    d/dt(pbad_transit1) <- ktr_pbad - ktr_pbad * pjnk * pbad_transit1
    d/dt(pbad_transit2) <- ktr_pbad * (pbad_transit1 - pbad_transit2)
    d/dt(pbad) <- ktr_pbad * (pbad_transit2 - pbad)

    # Eq. 17: Bcl-2 stimulated by pBAD^gamma and pNFkB, inhibited by pp53.
    d/dt(bcl2) <- kdeg_bcl2 / pp53 * pbad^gamma_bcl2 * pnfkb - kdeg_bcl2 * bcl2

    # Eq. 18: caspase 9 stimulated by caspase 8^gamma, inhibited by Bcl-2 and pNFkB.
    r_cas9 <- kdeg_cas9 / (bcl2 * pnfkb) * cas8^gamma_cas9 - kdeg_cas9 * cas9
    d/dt(cas9) <- r_cas9

    # Data S1 text after Eq. 19: for U266 the cell-dynamics parameters were
    # re-estimated with pAKT, pmTOR, p21, caspase 8 and caspase 9 fixed at
    # their 24 h values. The '_held' states follow those proteins until 24 h
    # and stay constant afterwards; time 0 is the start of treatment.
    pakt_rate <- ktr_pakt * (pakt_transit1 - pakt)
    d/dt(pakt_held) <- (t < 24) * pakt_rate
    d/dt(pmtor_held) <- (t < 24) * r_pmtor
    d/dt(p21_held) <- (t < 24) * r_p21
    d/dt(cas8_held) <- (t < 24) * r_cas8
    d/dt(cas9_held) <- (t < 24) * r_cas9

    # Eq. 19 (main text Eq. 1): treated cells.
    d/dt(mm_cells) <- kgrow / p21_held * pakt_held * pmtor_held * mm_cells -
      kdeath * cas8_held * cas9_held * mm_cells

    # Untreated control cells grow first-order at kg (Data S2 code XP(17)).
    d/dt(mm_cells_control) <- kgrow * mm_cells_control

    media(0) <- 0
    pjnk_transit1(0) <- 1
    pjnk(0) <- 1
    pnfkb_transit1(0) <- 1
    pnfkb_transit2(0) <- 1
    pnfkb(0) <- 1
    pakt_transit1(0) <- 1
    pakt(0) <- 1
    pmtor(0) <- 1
    cas8(0) <- 1
    pp53(0) <- 1
    p21(0) <- 1
    pbad_transit1(0) <- 1
    pbad_transit2(0) <- 1
    pbad(0) <- 1
    bcl2(0) <- 1
    cas9(0) <- 1
    pakt_held(0) <- 1
    pmtor_held(0) <- 1
    p21_held(0) <- 1
    cas8_held(0) <- 1
    cas9_held(0) <- 1
    mm_cells(0) <- 1
    mm_cells_control(0) <- 1

    pjnk ~ prop(propSd_pjnk)
    pnfkb ~ prop(propSd_pnfkb)
    pakt ~ prop(propSd_pakt)
    pmtor ~ prop(propSd_pmtor)
    cas8 ~ prop(propSd_cas8)
    pp53 ~ prop(propSd_pp53)
    p21 ~ prop(propSd_p21)
    pbad ~ prop(propSd_pbad)
    bcl2 ~ prop(propSd_bcl2)
    cas9 ~ prop(propSd_cas9)
    mm_cells ~ prop(propSd_mm_cells)
    mm_cells_control ~ prop(propSd_mm_cells_control)
  })
}
