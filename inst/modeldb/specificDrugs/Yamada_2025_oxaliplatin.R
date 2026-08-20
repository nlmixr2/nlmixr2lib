Yamada_2025_oxaliplatin <- function() {
  description <- "Three-compartment population PK model describing free and total plasma platinum simultaneously after oxaliplatin administration, with zero-order IV input, first-order elimination and a mono-exponentially decaying time-dependent free fraction, in patients with locally advanced unresectable or metastatic HER2-negative, CLDN18.2-positive gastric or gastroesophageal junction (G/GEJ) adenocarcinoma receiving mFOLFOX6 with or without concomitant zolbetuximab (Yamada 2025, ILUSTRO Cohort 2). The disposition ODEs describe FREE platinum; total platinum is derived as the free concentration divided by the time-dependent free fraction. Concomitant zolbetuximab reduces the distribution volume of all three compartments by 12.3%."
  reference <- paste(
    "Yamada A, Yang J, Bonate PL, Heo N, Poondru S. Population",
    "pharmacokinetic analysis of fluorouracil and oxaliplatin in the",
    "absence or presence of zolbetuximab in locally advanced unresectable",
    "or metastatic gastric or gastroesophageal junction adenocarcinoma.",
    "Cancer Chemother Pharmacol. 2025;95:89.",
    "doi:10.1007/s00280-025-04808-2 (PMCID PMC12449379).",
    "Parameter estimates are from Supplementary Table 2 ('Parameter",
    "estimates of best base model for oxaliplatin and models with",
    "zolbetuximab effect on oxaliplatin PK'), 'Zolbetuximab as a covariate",
    "on V1-V3 parameter estimate' column -- the final model per Results,",
    "'Evaluation of zolbetuximab effect on oxaliplatin PK' -- in the online",
    "Supplementary Material (280_2025_4808_MOESM1_ESM.pdf).",
    "Companion fluorouracil model from the same paper:",
    "modellib('Yamada_2025_fluorouracil').",
    sep = " "
  )
  vignette <- "Yamada_2025_fluorouracil_oxaliplatin_zolbetuximab"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. The ODE states carry FREE (ultrafilterable) platinum;
  # total platinum is an algebraic observable, not a state. The paper reports
  # CL in L/h and V in L, so with dose in mg the natural concentration unit is
  # mg/L = ug/mL; the paper's own tables are printed in ng platinum/mL (a
  # factor of 1000 larger). See the vignette source-trace table.
  #
  # DOSE BASIS -- IMPORTANT. `amt` for this model is PLATINUM mass, not
  # oxaliplatin mass, because the assay (ICP-MS) measures platinum and the
  # disposition parameters were estimated against platinum concentrations.
  # An oxaliplatin dose must therefore be converted before use:
  #
  #     amt_platinum_mg = dose_oxaliplatin_mg * 195.08 / 397.29
  #                     = dose_oxaliplatin_mg * 0.4910
  #
  # (195.08 g/mol = atomic mass of platinum; 397.29 g/mol = molar mass of
  # oxaliplatin, C8H14N2O4Pt.) The paper does not state the dose basis it
  # used in NONMEM, but it is recoverable and unambiguous: simulating the
  # 85 mg/m^2 dose as oxaliplatin mass overpredicts every one of the four
  # dose-normalized exposure metrics in Table 2 by a factor of 2.01-2.07,
  # against an oxaliplatin:platinum molar-mass ratio of 2.0365, whereas the
  # platinum-mass basis reproduces all four to within 1.5%. Dose
  # normalization in the paper's Table 2 is nevertheless by the
  # ADMINISTERED OXALIPLATIN dose in mg (that is what "ng platinum/mL/mg"
  # means there). The vignette walks this check explicitly.
  compartmentData <- list(
    central     = list(analyte = "platinum (free)", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "platinum (free)", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "platinum (free)", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CONMED_ZOLBETUXIMAB = list(
      description        = "Concomitant zolbetuximab coadministration (1 = oxaliplatin given in the presence of zolbetuximab, 0 = oxaliplatin alone)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (oxaliplatin alone; ILUSTRO Cohort 2 Cycle 1, in which zolbetuximab was deliberately delayed to Day 3 so that oxaliplatin PK could be sampled without it)",
      notes              = "Enters as the fractional covariate model P = theta_p * (1 + theta_zolbetuximab * X_zolbetuximab) (Methods, 'Evaluation of zolbetuximab effect as a covariate'; Article Equation a), applied with a single shared coefficient to V1, V2 and V3. The estimate is -0.123 (RSE 32.0%), i.e. a 12.3% reduction in each distribution volume and hence in steady-state volume; the paper judges this statistically significant (dOFV -13.834, likelihood-ratio P = 0.0002) but not clinically relevant. Effects on CL and on the free fraction fp were also screened and rejected (dOFV -3.649 and -10.654 respectively, against the 10.83 backward-elimination criterion at alpha = 0.001). Time-varying in principle -- a patient contributes X = 0 in Cycle 1 and X = 1 in Cycle 2 -- so this column belongs on the observation/dosing records rather than being a subject-level baseline.",
      source_name        = "X_zolbetuximab"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 21L,
    n_studies      = 1L,
    disease_state  = "Adults with previously untreated, locally advanced unresectable or metastatic HER2-negative, CLDN18.2-positive gastric or gastroesophageal junction (G/GEJ) adenocarcinoma (CLDN18.2 positivity defined as >= 75% of tumor cells with moderate to strong membranous CLDN18 staining intensity). Cohort 2 of the phase 2 ILUSTRO study (NCT03505320).",
    dose_range     = "Oxaliplatin 85 mg/m^2 IV infused over 2 h concurrently with folinic acid 400 mg/m^2, on Days 1, 15 and 29 of 42-day cycles for up to 4 cycles, as part of mFOLFOX6 (followed by 5-FU 400 mg/m^2 IV bolus over 5-15 min then 5-FU 2400 mg/m^2 over 46-48 h). Zolbetuximab 800 mg/m^2 loading dose then 600 mg/m^2, each a minimum 2-h IV infusion, on Days 1 and 22 of each cycle except Cycle 1, where it was given on Day 3 so that oxaliplatin PK could be sampled in its absence.",
    notes          = "Baseline clinical and demographic characteristics are not tabulated in this paper; it cites the previously published ILUSTRO Cohort 2 report (reference 4) for them, which was not available on disk at extraction time. Analysis dataset: 232 non-BQL and 5 BQL free-platinum concentrations and 237 non-BQL total-platinum concentrations originally; the 5 BQL free-platinum values and 34 non-BQL values (17 free, 17 total) were excluded during base model development for duplicate sampling or measurement error. PK sampling: predose, end of infusion, and 0.5, 1, 2, 5 and 24 h after the start of dosing in Cycles 1 (chemotherapy alone) and 2 (chemotherapy plus zolbetuximab). Assay: validated ICP-MS, LLOQ 1 ng platinum/mL (total) and 3 ng platinum/mL (free). Estimation in NONMEM 7.5.0; the final model's covariance step used MATRIX=S and its condition number was 582.53."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural disposition parameters -- Supplementary Table 2,
    # 'Zolbetuximab as a covariate on V1-V3' column. That column is the
    # final model: Results, 'Evaluation of zolbetuximab effect on
    # oxaliplatin PK' states "The final model was a 3-compartment model
    # ... including time-dependent fp; IIV on CL, V1-V3, fp0_A, and
    # fp0_B; and zolbetuximab effect on V1-V3". The Discussion quotes
    # exactly these values (CL 8.09 L/h; V1 51.1 L, V2 171 L, V3 1400 L).
    #
    # These describe FREE (ultrafilterable) platinum disposition.
    #
    # NOTE on the half-lives quoted in the Discussion (0.19, 13.0 and
    # 446 h for the alpha, beta and gamma phases): those are computed
    # from the BEST BASE model column (CL 8.15, V1 44.7, Q2 118, V2 161,
    # Q3 3.35, V3 1510), which reproduces them to 0.194 / 12.9 / 450 h.
    # The final-model column above gives 0.23 / 14.3 / 427 h. See the
    # vignette, which checks both.
    # ------------------------------------------------------------------
    lcl  <- log(8.09); label("Systemic clearance of free platinum CL (L/h)")                                   # Supp Table 2, V1-V3 column: CL 8.09 L/h (RSE 21.9%)
    lvc  <- log(51.1); label("Central volume of distribution V1 (L)")                                          # Supp Table 2, V1-V3 column: V1 51.1 L (RSE 21.3%)
    lq   <- log(113);  label("Intercompartmental clearance to peripheral1, Q2 (L/h)")                          # Supp Table 2, V1-V3 column: Q2 113 L/h (RSE 15.4%)
    lvp  <- log(171);  label("First peripheral volume of distribution V2 (L)")                                 # Supp Table 2, V1-V3 column: V2 171 L (RSE 16.0%)
    lq2  <- log(3.21); label("Intercompartmental clearance to peripheral2, Q3 (L/h)")                          # Supp Table 2, V1-V3 column: Q3 3.21 L/h (RSE 25.1%)
    lvp2 <- log(1400); label("Second peripheral volume of distribution V3 (L)")                                # Supp Table 2, V1-V3 column: V3 1400 L (RSE 192.1%)

    # ------------------------------------------------------------------
    # Time-dependent free fraction of plasma platinum -- Article
    # Equation b (Results, 'Base model for free and total platinum
    # concentration-time profiles'):
    #
    #     fp = fp0_A * exp(-alpha * TALD) + fp0_B
    #
    # where TALD is time after the start of infusion of the last dose.
    # Protein binding of platinum is irreversible and time dependent
    # (Discussion), so the free fraction decays from fp0_A + fp0_B at
    # the start of a dose towards the time-independent floor fp0_B.
    # With the final-model values that is 0.226 immediately after a dose
    # decaying towards 0.0318.
    #
    # Canonical names: fuA / kfu / fuB (log-transformed as lfuA / lkfu /
    # lfuB); see inst/references/parameter-names.md, "Time-varying
    # protein binding".
    # ------------------------------------------------------------------
    lfuA <- log(0.194);  label("Amplitude of the time-decaying component of the free fraction, fp0_A (fraction)")            # Supp Table 2, V1-V3 column: fp0_A 0.194 (RSE 13.8%)
    lkfu <- log(0.0393); label("First-order rate constant of time-dependent platinum protein binding, alpha (1/h)")          # Supp Table 2, V1-V3 column: alpha 0.0393 /h (RSE 17.2%)
    lfuB <- log(0.0318); label("Time-independent (asymptotic) free fraction, fp0_B (fraction)")                              # Supp Table 2, V1-V3 column: fp0_B 0.0318 (RSE 37.7%)

    # ------------------------------------------------------------------
    # Covariate effect -- Supplementary Table 2 'Zolbetuximab effect
    # (%RSE)' block, 'V1-V3' row. A single shared coefficient scales all
    # three volumes through P = theta_p * (1 + theta_zolb * X_zolb).
    # 1 + (-0.123) = 0.877, i.e. the 12.3% decrease reported in the
    # Abstract and Discussion.
    # ------------------------------------------------------------------
    e_zolbetuximab_vc_vp_vp2 <- -0.123; label("Fractional effect of concomitant zolbetuximab on V1, V2 and V3 (unitless)")  # Supp Table 2, V1-V3 column: -0.123 (RSE 32.0%)

    # ------------------------------------------------------------------
    # Inter-individual variability -- Supplementary Table 2 'IIV
    # [shrinkage]' rows, reported as CV%. Converted to log-normal
    # variances with the library convention omega^2 = log(CV^2 + 1).
    #
    # etalvc is a SINGLE eta shared by V1, V2 and V3: the paper reports
    # one "omega V1-V3" and Methods, 'Base model for free platinum
    # concentration-time profile' describes "a common IIV on distribution
    # volume of all compartments (V1, V2, and V3)". It is therefore
    # applied to all three volumes in model() rather than being split
    # into separate etalvc / etalvp / etalvp2 terms.
    #
    # No IIV was retained on alpha in the final model (Supplementary
    # Table 2 has no omega_alpha row), so lkfu carries no eta.
    #
    # Reported shrinkage: 14.3% (CL), 24.3% (V1-V3), 13.8% (fp0_A),
    # 19.0% (fp0_B).
    # ------------------------------------------------------------------
    etalvc  ~ 0.020524                                                                # Supp Table 2: omega_V1-V3 = 14.4 % CV -> log(1 + 0.144^2) = 0.020524
    etalcl  ~ 0.041956                                                                # Supp Table 2: omega_CL    = 20.7 % CV -> log(1 + 0.207^2) = 0.041956
    etalfuA ~ 0.025906                                                                # Supp Table 2: omega_fp0_A = 16.2 % CV -> log(1 + 0.162^2) = 0.025906
    etalfuB ~ 0.144987                                                                # Supp Table 2: omega_fp0_B = 39.5 % CV -> log(1 + 0.395^2) = 0.144987

    # ------------------------------------------------------------------
    # Residual variability -- Supplementary Table 2 'Residual variability
    # (%RSE)' rows. Independent proportional errors for the free and the
    # total platinum observations. Both are printed with an explicit
    # percent sign in the source table, so they are read directly as
    # CV% and converted to fractions.
    # ------------------------------------------------------------------
    propSd         <- 0.385; label("Proportional residual error on free platinum concentration Cc (fraction)")     # Supp Table 2, V1-V3 column: 38.5 % CV (RSE 11.2%)
    propSd_totalPt <- 0.146; label("Proportional residual error on total platinum concentration totalPt (fraction)")  # Supp Table 2, V1-V3 column: 14.6 % CV (RSE 10.8%)
  })

  model({
    # 1. Zolbetuximab effect on the distribution volumes. A single
    #    fractional coefficient is shared by V1, V2 and V3 (Supp Table 2
    #    'V1-V3' row), so it is computed once here.
    zolb_v <- 1 + e_zolbetuximab_vc_vp_vp2 * CONMED_ZOLBETUXIMAB

    # 2. Individual parameters. etalvc is the paper's single common
    #    IIV term on V1-V3 and so appears on all three volumes.
    cl  <- exp(lcl  + etalcl)
    vc  <- exp(lvc  + etalvc) * zolb_v
    q   <- exp(lq)
    vp  <- exp(lvp  + etalvc) * zolb_v
    q2  <- exp(lq2)
    vp2 <- exp(lvp2 + etalvc) * zolb_v

    # 3. Individual free-fraction parameters
    fuA <- exp(lfuA + etalfuA)
    kfu <- exp(lkfu)
    fuB <- exp(lfuB + etalfuB)

    # 4. Micro-constants
    kel <- cl  / vc
    k12 <- q   / vc
    k21 <- q   / vp
    k13 <- q2  / vc
    k31 <- q2  / vp2

    # 5. ODE system for FREE platinum. Zero-order input: oxaliplatin is
    #    given as a 2-h IV infusion, encoded on the dosing records
    #    (rate / dur), so no absorption compartment is required.
    d/dt(central)     <- -(kel + k12 + k13) * central + k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    # 6. Time-dependent free fraction (Article Equation b). TALD is time
    #    after the start of the last infusion, which is rxode2's
    #    tad(central) for a model dosed into central.
    fp <- fuA * exp(-kfu * tad(central)) + fuB

    # 7. Observations. Cc is the free (ultrafilterable) platinum
    #    concentration -- the concentration of the central ODE state.
    #    Total platinum is algebraic: free divided by the free fraction.
    Cc      <- central / vc
    totalPt <- Cc / fp

    Cc      ~ prop(propSd)
    totalPt ~ prop(propSd_totalPt)
  })
}
