Li_2017_brentuximab <- function() {
  description <- "Semimechanistic coupled population PK model for brentuximab vedotin antibody-drug conjugate (ADC) and its released small-molecule payload monomethyl auristatin E (MMAE) in adults with CD30-expressing hematologic malignancies (Li 2017). ADC is described by a linear 3-compartment model with first-order elimination; MMAE by a linear 2-compartment model with first-order elimination. MMAE formation is driven by (1) proteolytic degradation of the ADC (scaled by a time-decaying drug-antibody ratio DAR(t) and a cycle-dependent fraction Fmc = Cycle^Fm) and (2) a first-order deconjugation flux proportional to the per-ADC MMAE payload above the minimum-detectable DAR. Modeled in molar units (amount nmol, volume L, concentration pmol/mL = nmol/L = nM) following the paper's convention."
  reference <- "Li H, Han TH, Hunder NN, Jang G, Zhao B. Population Pharmacokinetics of Brentuximab Vedotin in Patients With CD30-Expressing Hematologic Malignancies. J Clin Pharmacol. 2017;57(9):1148-1158. doi:10.1002/jcph.920. PMID 28513851."
  vignette <- "Li_2017_brentuximab"
  units <- list(
    time          = "day",
    dosing        = "nmol",
    concentration = "pmol/mL"
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. Power effects on ADC (CL, V1, V2, V3, Q2, Q3) and MMAE (CLM, V4, V5, Q5) parameters; reference 75 kg per Li 2017 Table 3/4 (pooled-population median 74.8 kg rounded to 75 kg in the reference-subject definition).",
      source_name        = "BW"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Time-fixed. Enters only on ADC V1 with a multiplicative ratio e_sexf_vc = 0.873 for females (males are the reference). Paper reports 136/314 = 43% female.",
      source_name        = "SEX"
    ),
    CYCLE = list(
      description        = "Treatment cycle number, 1 = first dosing cycle, 2 = second, ... (integer count, time-varying across a multi-cycle treatment course)",
      units              = "(count)",
      type               = "count",
      reference_category = NULL,
      notes              = "Required for the MMAE submodel only: Fmc = CYCLE^Fm with Fm = -0.261 captures the cycle-over-cycle decline in ADC-to-MMAE proteolytic conversion thought to reflect tumor-burden reduction. Increment CYCLE at the start of each new dosing cycle; leave CYCLE >= 1 throughout (CYCLE^Fm is undefined at 0). Does not affect ADC PK or MMAE deconjugation flux.",
      source_name        = "CYCLE"
    )
  )

  population <- list(
    n_subjects     = 314L,
    n_studies      = 5L,
    age_range      = "12-87 years (median 35)",
    age_median     = "35 years",
    weight_range   = "41.4-168 kg (median 74.8)",
    weight_median  = "75 kg (reference-subject value; pooled-population median 74.8 kg)",
    sex_female_pct = 43.3,
    race_ethnicity = c(White = 85, Other = 15),
    disease_state  = "Relapsed/refractory CD30-expressing hematologic malignancies: Hodgkin lymphoma (HL, 77%), systemic anaplastic large cell lymphoma (sALCL, 21%), other (2%).",
    dose_range     = "0.1-3.6 mg/kg IV every 3 weeks (Studies 1, 3, 4, 5) or 0.4-1.4 mg/kg IV weekly for 3 weeks of each 4-week cycle (Study 2); 30-minute or 2-hour infusion. Dose capped at 180 mg (equivalent to 1.8 mg/kg for a 100-kg subject). Licensed regimen is 1.8 mg/kg IV every 3 weeks.",
    regions        = "Global (North America, Europe, Canada).",
    study_phase    = "2 phase 1 dose-ranging studies, 2 pivotal phase 2 studies (HL, sALCL), 1 phase 1 clinical pharmacology study (5 studies total).",
    n_observations = "7081 ADC concentrations + 7452 MMAE concentrations across 314 subjects.",
    reference_subject = "75 kg male, Cycle 1 (WT / 75)^exponent and CYCLE^Fm both equal 1.",
    notes          = "Baseline characteristics reported in Li 2017 Table 2; data collected November 2006 - June 2011. Renal function (creatinine clearance 22-150 mL/min, median 122), hepatic function (ALT, AST, bilirubin, albumin), tumor size, age, race, and manufacturing process were evaluated as candidate covariates and were not statistically significant (only body weight on ADC CL/V1/V2/V3/Q2/Q3 and sex on ADC V1 were retained)."
  )

  ini({
    # ADC structural parameters -- reference: 75 kg male (Li 2017 Table 3).
    # Modeled in molar units: amount in nmol, volume in L, clearance in L/day
    # so that concentration Cc = A / V is in nmol/L = pmol/mL = nM (matching
    # the paper's reported "concentration in units of pmol/mL").
    lcl  <- log(1.56);  label("ADC clearance (CL, L/day)")                         # Li 2017 Table 3: CL 1.56 L/d
    lvc  <- log(4.29);  label("ADC central volume (V1, L)")                        # Li 2017 Table 3: V1 4.29 L
    lq   <- log(2.83);  label("ADC intercompartmental clearance to peripheral 1 (Q2, L/day)") # Li 2017 Table 3: Q2 2.83 L/d
    lvp  <- log(3.83);  label("ADC peripheral volume 1 (V2, L)")                   # Li 2017 Table 3: V2 3.83 L
    lq2  <- log(0.708); label("ADC intercompartmental clearance to peripheral 2 (Q3, L/day)") # Li 2017 Table 3: Q3 0.708 L/d
    lvp2 <- log(9.52);  label("ADC peripheral volume 2 (V3, L)")                   # Li 2017 Table 3: V3 9.52 L

    # MMAE structural parameters -- reference: 75 kg (Li 2017 Table 4). CLM and
    # V4 are apparent parameters (ADC -> MMAE stoichiometry is already absorbed
    # into the DAR-based formation term in the ODE).
    lclm  <- log(55.7);   label("MMAE apparent clearance (CLM, L/day)")            # Li 2017 Table 4: CLM 55.7 L/d
    lvcm  <- log(79.8);   label("MMAE apparent central volume (V4, L)")            # Li 2017 Table 4: V4 79.8 L
    lqm   <- log(65.0);   label("MMAE apparent intercompartmental clearance (Q5, L/day)") # Li 2017 Table 4: Q5 65.0 L/d
    lvpm  <- log(28.1);   label("MMAE apparent peripheral volume (V5, L)")         # Li 2017 Table 4: V5 28.1 L
    lbeta <- log(0.0785); label("Macro rate constant of DAR decline (beta, 1/day)") # Li 2017 Table 4: beta 0.0785 /d

    # Fm (exponent on Cycle for the MMAE formation fraction Fmc = Cycle^Fm).
    # Paper reports Fm = -0.261 with BSV 130 %CV (Li 2017 Table 4). Stored as
    # fm = -exp(lfm + etalfm) so the magnitude is log-normally distributed and
    # the sign is held negative.
    lfm <- log(0.261);    label("Log magnitude of Fm (Fm = -exp(lfm)); Fm is the exponent in Fmc = CYCLE^Fm") # Li 2017 Table 4: Fm = -0.261

    # Covariate effects on ADC -- Li 2017 Table 3.
    e_wt_adc_cl <- 0.698; label("Power exponent of WT on ADC CL, Q2, Q3 (unitless)") # Li 2017 Table 3: Theta_BW,CL,Q2,Q3 = 0.698
    e_wt_adc_vc <- 0.503; label("Power exponent of WT on ADC V1, V2, V3 (unitless)") # Li 2017 Table 3: Theta_BW,V1,V2,V3 = 0.503
    e_sexf_vc   <- 0.873; label("Multiplicative effect of female sex on ADC V1 (ratio female:male)") # Li 2017 Table 3: Theta_SEX,V1 = 0.873

    # Covariate effects on MMAE -- allometric exponents fixed by the paper
    # (Li 2017 Table 4 footnotes b and c).
    e_wt_mmae_cl <- fix(0.75); label("Power exponent of WT on MMAE CLM, Q5 (unitless, fixed)") # Li 2017 Table 4 fixed
    e_wt_mmae_vc <- fix(1.0);  label("Power exponent of WT on MMAE V4, V5 (unitless, fixed)")  # Li 2017 Table 4 fixed

    # IIV (log-normal; omega^2 = log(CV^2 + 1)). CL-V1 block with correlation
    # 0.229 on the log scale per Li 2017 Table 3 (Corr(CL, V1) row).
    # omega_cl^2 = log(0.469^2 + 1) = 0.1988
    # omega_v1^2 = log(0.135^2 + 1) = 0.01806
    # cov       = 0.229 * sqrt(0.1988 * 0.01806) = 0.01372
    etalcl + etalvc ~ c(0.1988,
                        0.01372, 0.01806)     # Li 2017 Table 3: CL 46.9%, V1 13.5%, Corr 0.229
    etalq   ~ fix(0.02225)                    # Li 2017 Table 3: Q2 15%  CV (fixed in paper); omega^2 = log(1.0225) = 0.02225
    etalvp  ~ fix(0.06062)                    # Li 2017 Table 3: V2 25%  CV (fixed in paper); omega^2 = log(1.0625) = 0.06062
    etalq2  ~ 0.1859                          # Li 2017 Table 3: Q3 45.2% CV; omega^2 = log(1.2043) = 0.1859
    etalvp2 ~ 0.7530                          # Li 2017 Table 3: V3 106%  CV; omega^2 = log(2.1236) = 0.7530

    # CLM-V4 block with correlation 0.634 on the log scale per Li 2017 Table 4.
    # omega_clm^2 = log(0.607^2 + 1) = 0.3137
    # omega_v4^2  = log(0.782^2 + 1) = 0.4773
    # cov         = 0.634 * sqrt(0.3137 * 0.4773) = 0.2454
    etalclm + etalvcm ~ c(0.3137,
                          0.2454, 0.4773)     # Li 2017 Table 4: CLM 60.7%, V4 78.2%, Corr 0.634
    etalbeta ~ 0.6741                         # Li 2017 Table 4: beta 98.1% CV; omega^2 = log(1.9624) = 0.6741
    etalfm   ~ 0.9895                         # Li 2017 Table 4: Fm 130% CV; omega^2 = log(2.69) = 0.9895
    # Q5 and V5 have BSV fixed to 0 in Li 2017 Table 4 -- no IIV terms.

    # Residual error -- converted from the paper's mixed units (ADC sigma1 in
    # ug/mL, MMAE sigma1 in ng/mL) to the native molar unit (pmol/mL) so the
    # SDs are directly applied to Cc and Cmmae as computed by the model.
    #   ADC:  0.0125 ug/mL * (1000 / MW_ADC_Da) = 0.0125 * 1000 / 153000
    #                                           = 8.17e-5 ng/pmol... re-derive:
    #         1 ug/mL ADC = 10^-6 g/mL / 153000 g/mol * 10^12 pmol/mol
    #                     * 1 mL = 10^-6 / 153000 * 10^12 = 6.536 pmol/mL.
    #         0.0125 ug/mL = 0.0817 pmol/mL. (Equivalent to LLOQ 12.5 ng/mL.)
    #   MMAE: 0.0119 ng/mL * (1 / MW_MMAE_kDa) = 0.0119 / 0.718 = 0.01658 pmol/mL.
    CcaddSd      <- fix(0.0817); label("Additive residual error on ADC Cc (pmol/mL; equivalent to 0.0125 ug/mL = LLOQ 12.5 ng/mL, fixed)") # Li 2017 Table 3: sigma1 0.0125 ug/mL fixed
    CcpropSd     <- 0.329;       label("Proportional residual error on ADC Cc (fraction)") # Li 2017 Table 3: sigma2 32.9% CV
    CmmaeaddSd   <- 0.01658;     label("Additive residual error on MMAE Cmmae (pmol/mL; equivalent to 0.0119 ng/mL)") # Li 2017 Table 4: sigma1 0.0119 ng/mL
    CmmaepropSd  <- 0.368;       label("Proportional residual error on MMAE Cmmae (fraction)") # Li 2017 Table 4: sigma2 36.8% CV
  })

  model({
    # Individual PK parameters with allometric/sex adjustment. WT reference
    # is 75 kg; SEXF is 0 for males (reference) and 1 for females, so
    # e_sexf_vc^SEXF equals 1 for males and 0.873 for females (Li 2017
    # Table 3 Eq. 2).
    cl_adc <- exp(lcl  + etalcl)  * (WT / 75)^e_wt_adc_cl
    v1     <- exp(lvc  + etalvc)  * (WT / 75)^e_wt_adc_vc * e_sexf_vc^SEXF
    q2     <- exp(lq   + etalq)   * (WT / 75)^e_wt_adc_cl
    v2     <- exp(lvp  + etalvp)  * (WT / 75)^e_wt_adc_vc
    q3     <- exp(lq2  + etalq2)  * (WT / 75)^e_wt_adc_cl
    v3     <- exp(lvp2 + etalvp2) * (WT / 75)^e_wt_adc_vc

    cl_mmae <- exp(lclm + etalclm) * (WT / 75)^e_wt_mmae_cl
    v4      <- exp(lvcm + etalvcm) * (WT / 75)^e_wt_mmae_vc
    q5      <- exp(lqm)            * (WT / 75)^e_wt_mmae_cl
    v5      <- exp(lvpm)           * (WT / 75)^e_wt_mmae_vc

    beta <- exp(lbeta + etalbeta)
    fm   <- -exp(lfm + etalfm)   # Enforce sign negative (paper point estimate Fm = -0.261)

    # Micro-constants (Li 2017 supplement M1)
    k10 <- cl_adc  / v1
    k12 <- q2      / v1
    k21 <- q2      / v2
    k13 <- q3      / v1
    k31 <- q3      / v3
    k40 <- cl_mmae / v4
    k45 <- q5      / v4
    k54 <- q5      / v5

    # Time-varying drug-antibody ratio within each dosing cycle (Li 2017
    # Eq. 5). DAR resets to DAR0 = 4 at each dose and decays exponentially
    # toward DAR0 * alpha = 1 (alpha = 0.25 fixed because the ELISA
    # quantifies ADC with at least one drug). tad() returns time since the
    # most recent dose (days).
    dar <- 4 * (0.25 + 0.75 * exp(-beta * tad()))   # Li 2017 Eq. 5 and supplement M1

    # Cycle-dependent proteolytic fraction (Li 2017 Eq. 6). CYCLE = 1 for the
    # first treatment cycle, incremented at each subsequent cycle.
    fmc <- CYCLE^fm

    # ODE system: 5 PK compartments + 2 bookkeeping compartments for the two
    # MMAE-formation pathways (Li 2017 supplement M1). Molar balance: A(1)-A(5)
    # are in nmol; V1/V4 are in L; concentrations are nmol/L = pmol/mL.
    d/dt(adc_central)     <- -(k10 + k12 + k13) * adc_central + k21 * adc_peripheral1 + k31 * adc_peripheral2
    d/dt(adc_peripheral1) <-  k12 * adc_central - k21 * adc_peripheral1
    d/dt(adc_peripheral2) <-  k13 * adc_central - k31 * adc_peripheral2

    # MMAE central: proteolytic flux (cycle-scaled) + deconjugation flux
    # - clearance - distribution to peripheral + return from peripheral.
    # The deconjugation term uses (DAR - DAR0 * alpha) so that at DAR = 1
    # (long after a dose) the deconjugation flux reaches zero.
    d/dt(mmae_central)    <-  k10 * adc_central * dar * fmc +
                               adc_central * (dar - 4 * 0.25) * beta -
                               (k40 + k45) * mmae_central + k54 * mmae_peripheral
    d/dt(mmae_peripheral) <-  k45 * mmae_central - k54 * mmae_peripheral

    # Cumulative pathway bookkeeping (Li 2017 supplement M1). These integrate
    # the total MMAE amount released along each pathway and are not used in
    # the observation equations; they enable reproducing the paper's "~13%
    # deconjugation at steady state" derivation by computing
    # deconjugation_pct = pathway_deconjugation / (pathway_proteolytic + pathway_deconjugation) * 100.
    d/dt(pathway_proteolytic)   <- k10 * adc_central * dar
    d/dt(pathway_deconjugation) <- adc_central * (dar - 4 * 0.25) * beta

    # Observations -- concentrations in pmol/mL (= nmol/L = nM).
    Cc    <- adc_central  / v1
    Cmmae <- mmae_central / v4

    Cc    ~ add(CcaddSd)    + prop(CcpropSd)
    Cmmae ~ add(CmmaeaddSd) + prop(CmmaepropSd)
  })
}
