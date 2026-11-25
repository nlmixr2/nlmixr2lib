Ramakrishnan_2019_bortezomib_rpmi8226 <- function() {
  description <- "QSP. In vitro (RPMI8226 human multiple myeloma cell line). Pharmacodynamic systems model of bortezomib signaling from Ramakrishnan 2019: bortezomib in the culture medium (20 nM, first-order loss with a 144 h half-life) linearly stimulates the synthesis of pJNK (direct), pNFkB and pAKT (one transit compartment each), and a network of indirect-response turnover equations for pmTOR, pp53, p21, pBAD (three transit compartments), Bcl-2, caspase 8 and caspase 9 (relative expression, 1 = untreated control at 0 h) drives myeloma cell growth (pAKT x pmTOR / p21) and death (caspase 8 x caspase 9). Structure from the authors' ADAPT 5 code. Naive-pooled maximum-likelihood fit (ADAPT 5); no between-well variability. Sibling cell-line models: Ramakrishnan_2019_bortezomib_u266, _mm1s, _ncih929."
  reference <- paste(
    "Ramakrishnan V, Mager DE. (2019). Pharmacodynamic Models of Differential",
    "Bortezomib Signaling Across Several Cell Lines of Multiple Myeloma.",
    "CPT Pharmacometrics Syst Pharmacol 8(3):146-157. doi:10.1002/psp4.12358.",
    "Equations from the ADAPT 5 model code for RPMI8226 (Supplementary Data S2)",
    "and the Data S1 description of the RPMI8226 differences from U266;",
    "parameter values from Table 1."
  )
  vignette <- "Ramakrishnan_2019_bortezomib_myeloma"

  # Paper-mechanistic states with no canonical analogue in
  # inst/references/compartment-names.md: bortezomib in the culture medium,
  # the ten measured signaling proteins (relative expression), the transit
  # compartments that delay pNFkB / pAKT / pBAD (pattern below), and the
  # treated and untreated relative myeloma cell numbers.
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
    cas8 = list(
      analyte = "caspase 8",
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
    pp53 = list(
      analyte = "phosphorylated p53 (pp53)",
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
    pbad_transit2 = list(
      analyte = "pBAD signal transit compartment 2",
      units = "unitless (relative)",
      specimen = "not applicable",
      verified = TRUE
    ),
    pbad_transit3 = list(
      analyte = "pBAD signal transit compartment 3",
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
    p21 = list(analyte = "p21", units = "unitless (relative expression)", specimen = "not applicable", verified = TRUE),
    mm_cells_control = list(
      analyte = "viable untreated RPMI8226 cells",
      units = "unitless (relative to untreated control at 0 h)",
      specimen = "not applicable",
      verified = TRUE
    ),
    mm_cells = list(
      analyte = "viable RPMI8226 cells under bortezomib",
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
    species = "in vitro (RPMI8226 human multiple myeloma cell line)",
    n_subjects = "Cell-culture time courses (naive pooled): 10 signaling proteins by Western blot and viable-cell counts, untreated control and 20 nM bortezomib",
    n_studies = 1L,
    disease_state = "Multiple myeloma ('less sensitive' cell line, mutated TP53)",
    dose_range = "Single bortezomib concentration of 20 nM in the culture medium",
    regions = "University at Buffalo, NY (USA)",
    notes = "Data from Ramakrishnan and Mager 2018 (J Pharmacol Exp Ther 365:734-751). Protein data truncated at 24 h for model development; cell proliferation followed to 96 h and fitted without holding any protein at its 24 h value."
  )

  ini({
    # Bortezomib loss from the medium -- main text Methods: 'mono-exponential
    # decay with a half-life of 144 hours'; Data S2 code KDEG = 0.0048135208.
    lkdeg_media <- fixed(log(log(2) / 144))
    label("First-order bortezomib degradation rate in culture medium (1/h)") # Methods: t1/2 = 144 h (Chudasama 2015)

    # Upstream drug targets -- Data S2 XP(2), XP(4)-XP(7); Table 1 'RPMI8226'.
    lkdeg_pjnk <- log(0.07269)
    label("pJNK degradation rate constant kdpJNK (1/h)") # Table 1 RPMI8226 'k dpJNK' = 0.07269
    slope_pjnk <- 0.2965
    label("Linear stimulation of pJNK synthesis by bortezomib SmpJNK (1/nM)") # Table 1 RPMI8226 'S mpJNK' = 0.2965
    lktr_pnfkb <- log(0.07685)
    label("pNFkB transit rate constant ktpNFkB (1/h)") # Table 1 RPMI8226 'k tpNFkB' = 0.07685
    slope_pnfkb <- 0.04046
    label("Linear stimulation of pNFkB synthesis by bortezomib SmpNFkB (1/nM)") # Table 1 RPMI8226 'S mpNFkB' = 0.04046
    lktr_pakt <- log(0.04528)
    label("pAKT transit rate constant ktpAKT (1/h)") # Table 1 RPMI8226 'k tpAKT' = 0.04528
    slope_pakt <- 1.011
    label("Linear stimulation of pAKT synthesis by bortezomib SmpAKT (1/nM)") # Table 1 RPMI8226 'S mpAKT' = 1.011

    # Downstream proteins -- Data S2 XP(3), XP(8)-XP(16); Table 1 'RPMI8226'.
    lkdeg_pmtor <- log(0.0007853)
    label("pmTOR degradation rate constant kdpmTOR (1/h)") # Table 1 RPMI8226 'k dpmTOR' = 0.0007853
    lkdeg_cas8 <- log(0.1656)
    label("Caspase 8 degradation rate constant kdCas8 (1/h)") # Table 1 RPMI8226 'k dCas8' = 0.1656
    gamma_cas8 <- 5.640
    label("Power coefficient on pNFkB inhibition of caspase 8 synthesis (unitless)") # Table 1 RPMI8226 'gamma Cas8' = 5.640
    lkdeg_pp53 <- log(0.04395)
    label("pp53 degradation rate constant kdpp53 (1/h)") # Table 1 RPMI8226 'k dpp53' = 0.04395
    gamma_pp53 <- 1.746
    label("Power coefficient on pJNK stimulation of pp53 synthesis (unitless)") # Table 1 RPMI8226 'gamma pp53' = 1.746
    lktr_pbad <- log(0.2257)
    label("pBAD transit rate constant ktpBAD (1/h)") # Table 1 RPMI8226 'k tpBAD' = 0.2257
    gamma_pbad <- 0.3705
    label("Power coefficient on pJNK-stimulated loss of pBAD (unitless)") # Table 1 RPMI8226 'gamma pBAD' = 0.3705
    lkdeg_bcl2 <- log(0.008890)
    label("Bcl-2 degradation rate constant kdBcl2 (1/h)") # Table 1 RPMI8226 'k dBcl-2' = 0.008890
    lkdeg_cas9 <- log(0.008027)
    label("Caspase 9 degradation rate constant kdCas9 (1/h)") # Table 1 RPMI8226 'k dCas9' = 0.008027
    gamma_cas9 <- fixed(5)
    label("Power coefficient on caspase 8 stimulation of caspase 9 synthesis (unitless)") # Table 1 RPMI8226 'gamma Cas9' = 5.000*
    lkdeg_p21 <- log(0.01834)
    label("p21 degradation rate constant kdp21 (1/h)") # Table 1 RPMI8226 'k dp21' = 0.01834
    gamma_p21 <- fixed(5)
    label("Power coefficient on pp53 stimulation of p21 synthesis (unitless)") # Table 1 RPMI8226 'gamma p21' = 5.000*

    # Myeloma cell growth and death -- Data S2 XP(17)-XP(18).
    lkgrow <- log(0.01680)
    label("Myeloma cell growth rate constant kg (1/h)") # Table 1 RPMI8226 'k g' = 0.01680
    lkdeath <- log(0.05166)
    label("Myeloma cell death rate constant kd (1/h)") # Table 1 RPMI8226 'k d' = 0.05166

    # Residual error -- Data S1 Eq. 2 / Data S2 VARMOD, Var = (sigma1 + sigma2 * Y)^2.
    # Table 1 reports only the slope sigma2 for each output type.
    propSd_pjnk <- 0.4324
    label("Proportional residual error, pJNK (fraction)") # Table 1 RPMI8226 'sigma2 Protein' = 0.4324 (shared by all proteins)
    propSd_pnfkb <- 0.4324
    label("Proportional residual error, pNFkB (fraction)") # Table 1 RPMI8226 'sigma2 Protein' = 0.4324
    propSd_pakt <- 0.4324
    label("Proportional residual error, pAKT (fraction)") # Table 1 RPMI8226 'sigma2 Protein' = 0.4324
    propSd_pmtor <- 0.4324
    label("Proportional residual error, pmTOR (fraction)") # Table 1 RPMI8226 'sigma2 Protein' = 0.4324
    propSd_cas8 <- 0.4324
    label("Proportional residual error, caspase 8 (fraction)") # Table 1 RPMI8226 'sigma2 Protein' = 0.4324
    propSd_pp53 <- 0.4324
    label("Proportional residual error, pp53 (fraction)") # Table 1 RPMI8226 'sigma2 Protein' = 0.4324
    propSd_p21 <- 0.4324
    label("Proportional residual error, p21 (fraction)") # Table 1 RPMI8226 'sigma2 Protein' = 0.4324
    propSd_pbad <- 0.4324
    label("Proportional residual error, pBAD (fraction)") # Table 1 RPMI8226 'sigma2 Protein' = 0.4324
    propSd_bcl2 <- 0.4324
    label("Proportional residual error, Bcl-2 (fraction)") # Table 1 RPMI8226 'sigma2 Protein' = 0.4324
    propSd_cas9 <- 0.4324
    label("Proportional residual error, caspase 9 (fraction)") # Table 1 RPMI8226 'sigma2 Protein' = 0.4324
    propSd_mm_cells <- 0.4118
    label("Proportional residual error, treated myeloma cells (fraction)") # Table 1 RPMI8226 'sigma2 MMcell' = 0.4118
    propSd_mm_cells_control <- 0.4118
    label("Proportional residual error, untreated myeloma cells (fraction)") # Table 1 RPMI8226 'sigma2 MMcell' = 0.4118 (shared)
  })

  model({
    kdeg_media <- exp(lkdeg_media)
    kdeg_pjnk <- exp(lkdeg_pjnk)
    ktr_pnfkb <- exp(lktr_pnfkb)
    ktr_pakt <- exp(lktr_pakt)
    kdeg_pmtor <- exp(lkdeg_pmtor)
    kdeg_cas8 <- exp(lkdeg_cas8)
    kdeg_pp53 <- exp(lkdeg_pp53)
    ktr_pbad <- exp(lktr_pbad)
    kdeg_bcl2 <- exp(lkdeg_bcl2)
    kdeg_cas9 <- exp(lkdeg_cas9)
    kdeg_p21 <- exp(lkdeg_p21)
    kgrow <- exp(lkgrow)
    kdeath <- exp(lkdeath)

    # Every protein starts at 1 (untreated control, 0 h), so each zero-order
    # synthesis rate equals its first-order loss rate (Data S1 Eq. 1; the
    # code sets ks = kd * 1).

    # XP(1): bortezomib in the medium (Cb, nM); a 20 nM dose into media starts treatment.
    d/dt(media) <- -kdeg_media * media

    # XP(2): pJNK stimulated directly by bortezomib (no transit compartment).
    d/dt(pjnk) <- kdeg_pjnk * (1 + slope_pjnk * media) - kdeg_pjnk * pjnk

    # XP(3): caspase 8 stimulated by pJNK, inhibited by pNFkB^gamma.
    d/dt(cas8) <- kdeg_cas8 / pnfkb^gamma_cas8 * pjnk - kdeg_cas8 * cas8

    # XP(4)-XP(5): pNFkB with one transit compartment. The code writes the
    # drug term as smpnfkb * X(1)**G; G is not among the Table 1 estimates
    # and the paper describes every bortezomib stimulation as linear, so the
    # term is linear here (G = 1).
    d/dt(pnfkb_transit1) <- ktr_pnfkb * (1 + slope_pnfkb * media) - ktr_pnfkb * pnfkb_transit1
    d/dt(pnfkb) <- ktr_pnfkb * (pnfkb_transit1 - pnfkb)

    # XP(6)-XP(7): pAKT with one transit compartment.
    d/dt(pakt_transit1) <- ktr_pakt * (1 + slope_pakt * media) - ktr_pakt * pakt_transit1
    d/dt(pakt) <- ktr_pakt * (pakt_transit1 - pakt)

    # XP(8): pmTOR stimulated by pAKT.
    d/dt(pmtor) <- kdeg_pmtor * pakt - kdeg_pmtor * pmtor

    # XP(9): pp53 stimulated by pJNK^gamma, inhibited by pAKT.
    d/dt(pp53) <- kdeg_pp53 / pakt * pjnk^gamma_pp53 - kdeg_pp53 * pp53

    # XP(10)-XP(13): pJNK^gamma accelerates loss of the first pBAD
    # compartment, followed by three transit compartments.
    d/dt(pbad_transit1) <- ktr_pbad - ktr_pbad * pjnk^gamma_pbad * pbad_transit1
    d/dt(pbad_transit2) <- ktr_pbad * (pbad_transit1 - pbad_transit2)
    d/dt(pbad_transit3) <- ktr_pbad * (pbad_transit2 - pbad_transit3)
    d/dt(pbad) <- ktr_pbad * (pbad_transit3 - pbad)

    # XP(14): Bcl-2 stimulated by pBAD and pNFkB, inhibited by pp53.
    d/dt(bcl2) <- kdeg_bcl2 / pp53 * pbad * pnfkb - kdeg_bcl2 * bcl2

    # XP(15): caspase 9 stimulated by caspase 8^gamma, inhibited by Bcl-2 and pNFkB.
    d/dt(cas9) <- kdeg_cas9 / (bcl2 * pnfkb) * cas8^gamma_cas9 - kdeg_cas9 * cas9

    # XP(16): p21 stimulated by pp53^gamma, inhibited by pAKT.
    d/dt(p21) <- kdeg_p21 / pakt * pp53^gamma_p21 - kdeg_p21 * p21

    # XP(17): untreated control cells grow first-order at kg.
    d/dt(mm_cells_control) <- kgrow * mm_cells_control

    # XP(18): treated cells (main text Eq. 1).
    d/dt(mm_cells) <- kgrow / p21 * pmtor * pakt * mm_cells -
      kdeath * cas8 * cas9 * mm_cells

    media(0) <- 0
    pjnk(0) <- 1
    cas8(0) <- 1
    pnfkb_transit1(0) <- 1
    pnfkb(0) <- 1
    pakt_transit1(0) <- 1
    pakt(0) <- 1
    pmtor(0) <- 1
    pp53(0) <- 1
    pbad_transit1(0) <- 1
    pbad_transit2(0) <- 1
    pbad_transit3(0) <- 1
    pbad(0) <- 1
    bcl2(0) <- 1
    cas9(0) <- 1
    p21(0) <- 1
    mm_cells_control(0) <- 1
    mm_cells(0) <- 1

    pjnk ~ prop(propSd_pjnk)
    cas8 ~ prop(propSd_cas8)
    pnfkb ~ prop(propSd_pnfkb)
    pakt ~ prop(propSd_pakt)
    pmtor ~ prop(propSd_pmtor)
    pp53 ~ prop(propSd_pp53)
    pbad ~ prop(propSd_pbad)
    bcl2 ~ prop(propSd_bcl2)
    cas9 ~ prop(propSd_cas9)
    p21 ~ prop(propSd_p21)
    mm_cells_control ~ prop(propSd_mm_cells_control)
    mm_cells ~ prop(propSd_mm_cells)
  })
}
