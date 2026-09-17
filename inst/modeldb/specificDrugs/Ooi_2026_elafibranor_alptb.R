Ooi_2026_elafibranor_alptb <- function() {
  description <- paste0(
    "Joint indirect-response exposure-response model for alkaline ",
    "phosphatase (ALP) and total bilirubin in patients with primary biliary ",
    "cholangitis receiving elafibranor (Ooi 2026 final joint ALP-TB model, ",
    "Table 2). Each biomarker is a turnover compartment held at its ",
    "estimated baseline by a zero-order production rate and a first-order ",
    "degradation rate parameterised as a turnover half-life. A placebo ",
    "effect multiplies the production rate from the first dose onwards. The ",
    "drug effect also multiplies production and is driven by a STATIC ",
    "per-subject exposure metric, the sum of the steady-state dosing-interval ",
    "AUCs of elafibranor and its equipotent active metabolite GFT1007, ",
    "supplied as the covariates AUC_ELA and AUC_GFT1007: an inhibitory Emax ",
    "function with Hill coefficient fixed to 1 for ALP, and a linear slope ",
    "for total bilirubin. This is a SEQUENTIAL PK/PD model, so it contains ",
    "no PK compartments and takes no dose records; compute the two AUCs from ",
    "modellib('Ooi_2026_elafibranor') and ",
    "modellib('Ooi_2026_elafibranor_gft1007') as F * dose / CL, or supply ",
    "observed individual values. Baseline ALP rises 27.7% with any NCI ",
    "hepatic impairment and baseline total bilirubin follows a power ",
    "function of liver stiffness. The joint model estimates the correlation ",
    "between the two baseline random effects and a strong negative ",
    "correlation between the logit-scale Emax and the ALP placebo effect. ",
    "Residual error is combined additive plus proportional for each ",
    "biomarker, with inter-individual variability on the proportional ",
    "component. Fitted to 206 patients (1892 ALP and 1693 total-bilirubin ",
    "observations) from the phase II GFT505B-216-1 and phase III ELATIVE ",
    "(GFT505B-319-1) trials."
  )

  reference <- paste(
    "Ooi QX, Brendel K, van Beek S, Aguiar Zdovc J, Bardol M, Dehez M.",
    "Population Pharmacokinetics and Pharmacokinetics-Pharmacodynamics",
    "Analyses of Elafibranor to Support Dose Selection in Primary Biliary",
    "Cholangitis. CPT Pharmacometrics Syst Pharmacol. 2026;15(0):e70247.",
    "doi:10.1002/psp4.70247.",
    sep = " "
  )

  vignette <- "Ooi_2026_elafibranor"

  # Time is days since first dose (the NONMEM data item TSFDD is renamed to
  # TIME in Supplementary Datafile S6 $INPUT), matching the turnover
  # half-lives reported in days in Table 2. There are no dose records.
  units <- list(
    time = "day",
    dosing = "none (exposure enters as a static covariate)",
    concentration = "U/L (ALP); umol/L (total bilirubin)"
  )

  covariateData <- list(
    AUC_ELA = list(
      description = "Steady-state dosing-interval AUC of elafibranor",
      units = "umol*h/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Per-subject, time-fixed. Supplementary Datafile S6 $PK computes it ",
        "as AUCSSP = F * dose_mg / CL * 1e6 / 384.49 with CL in mL/h, i.e. ",
        "the individual relative bioavailability times the milligram dose ",
        "divided by the individual apparent clearance, converted to umol ",
        "using the elafibranor molecular weight 384.49 g/mol. Set to 0 for ",
        "placebo, which makes the drug term vanish exactly. Median of the ",
        "AUC sum was 32.3 umol*h/L on 80 mg/day and 39.3 umol*h/L on ",
        "120 mg/day (Figure 5 legend)."
      ),
      source_name = "AUCSSP"
    ),
    AUC_GFT1007 = list(
      description = "Steady-state dosing-interval AUC of GFT1007, the active metabolite of elafibranor",
      units = "umol*h/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Per-subject, time-fixed. Supplementary Datafile S6 $PK computes it ",
        "as AUCSSM = F * dose_mg / CL * 1e6 / 386.51, i.e. the same form as ",
        "AUC_ELA but converted with the GFT1007 molecular weight ",
        "386.51 g/mol. Enters the model only through the sum ",
        "AUC_ELA + AUC_GFT1007, because elafibranor and GFT1007 were shown ",
        "to be equipotent (Methods 2.3.2). Set to 0 for placebo. GFT1007 ",
        "contributes about five sixths of the sum at 80 mg/day."
      ),
      source_name = "AUCSSM"
    ),
    HEPIMP = list(
      description = "Any NCI hepatic impairment at baseline",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (NCI hepatic impairment score 0, normal hepatic function)",
      notes = paste0(
        "Supplementary Datafile S6 $PK codes the effect as ",
        "IF(NCIHISN.EQ.0) 1 else (1 + THETA), i.e. any score above 0 - mild ",
        "or worse - shares one coefficient, which is exactly the canonical ",
        "HEPIMP dichotomy. In the PKPD analysis set 37.4% had score 0, ",
        "62.1% score 1 and 0.5% score 2 (Table S2), so HEPIMP = 1 for 62.6%."
      ),
      source_name = "NCIHISN"
    ),
    LSM = list(
      description = "Baseline liver stiffness measured by transient elastography",
      units = "kPa",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Power model normalised to 8.1 kPa on the baseline total-bilirubin ",
        "typical value (Supplementary Datafile S6 $PK, ",
        "EBTBILLIVSTBL = (LIVSTBL/8.1)**THETA(12); 8.1 kPa is also the ",
        "reference-patient value in the Figure S4 forest plots). Mean ",
        "10.2 kPa (SD 8.2) in the PKPD analysis set, missing in 25.7% ",
        "(Table S2); the control stream sets the factor to 1 when the ",
        "measurement is missing, which is reproduced by supplying ",
        "LSM = 8.1 for such subjects."
      ),
      source_name = "LIVSTBL"
    )
  )

  compartmentData <- list(
    alp = list(
      analyte = "alkaline phosphatase",
      units = "U/L",
      specimen = "serum",
      verified = TRUE
    ),
    tbili = list(
      analyte = "total bilirubin",
      units = "umol/L",
      specimen = "serum",
      verified = TRUE
    )
  )

  # `alp` and `tbili` follow the register's bare clinical-biomarker
  # PD-output compartment pattern (`ast`, `cpk`, `ldl`, `hdl`, `urate`,
  # `Hba1c`), but the operator ruling is that a compartment canonical needs a
  # SECOND independent paper, so they are declared paper-specific here. A
  # future ALP or bilirubin turnover model is the trigger to promote them.
  paper_specific_compartments <- c("alp", "tbili")

  # The logit-scale Emax random effect and the random effects on the
  # proportional residual-error magnitudes.
  paper_specific_etas <- c(
    "etalogitemax_alp",
    "etapropSd_alp",
    "etapropSd_tbili"
  )

  population <- list(
    species = "human",
    n_subjects = 206L,
    n_studies = 2L,
    n_observations = 3585L,
    age_range = "not reported; mean 57.0 years (SD 8.6)",
    age_median = "mean 57.0 years (SD 8.6) (Table S2)",
    sex_female_pct = 95.6,
    race_ethnicity = c(
      White = 92.7,
      Asian = 1.9,
      `Black or African American` = 1.5,
      `American Indian or Alaska Native` = 0.5,
      `Multiple or other` = 2.4,
      `Unknown or not reported` = 1.0
    ),
    disease_state = "primary biliary cholangitis; 95.6% on prior ursodeoxycholic acid and 6.8% with prior obeticholic acid",
    dose_range = "placebo, elafibranor 80 mg/day or 120 mg/day for up to 12 weeks (phase II) or 52 weeks (phase III)",
    regions = "not reported",
    hepatic_function = "NCI hepatic impairment score 0 in 37.4%, 1 in 62.1%, 2 in 0.5%; mean baseline liver stiffness 10.2 kPa (SD 8.2), missing in 25.7%",
    notes = paste0(
      "Baseline ALP 308 U/L (SD 122) and baseline total bilirubin ",
      "9.79 umol/L (SD 5.1) (Table S2). The ALP data set held 206 patients ",
      "and 1892 observations, the total-bilirubin data set 205 patients and ",
      "1693 observations (Results 3.2.1); the joint model excluded one ",
      "patient (ID 3183) from both (Supplementary Datafile S6 $DATA and ",
      "Figure 4 legend). Studies GFT505B-216-1 (phase II, NCT03124108) and ",
      "GFT505B-319-1 (phase III ELATIVE, NCT04526665)."
    )
  )

  ini({
    # ==================================================================
    # BASELINES AND TURNOVER -- Table 2 (final joint ALP-TB model)
    # ==================================================================
    lrbase_alp <- log(251)
    label("Typical baseline alkaline phosphatase (U/L)")           # Table 2 'BASE ALP, typical (U/L)' 251 [RSE 5.81%]
    lrbase_tbili <- log(8.52)
    label("Typical baseline total bilirubin (umol/L)")             # Table 2 'BASE TB, typical (umol/L)' 8.52 [RSE 3.20%]
    lthalfrec_alp <- log(10.7)
    label("Alkaline phosphatase turnover half-life (day)")         # Table 2 'HL ALP (day)' 10.7 [RSE 4.43%]; footnote HL_ALP = ln(2)/K_out,ALP
    lthalfrec_tbili <- log(1240)
    label("Total bilirubin turnover half-life (day)")              # Table 2 'HL TB (day)' 1240 [RSE 29.0%]; footnote HL_TB = ln(2)/K_out,TB

    # ==================================================================
    # DRUG EFFECT -- Table 2
    # Both effects multiply the zero-order production rate; a negative Emax
    # and a negative slope therefore inhibit production.
    # ==================================================================
    emax_alp <- -0.731
    label("Maximum fractional inhibition of alkaline phosphatase production (fraction)")  # Table 2 'E max' -0.731 [RSE 8.68%]
    lec50_alp <- log(24.1)
    label("Sum of elafibranor and GFT1007 steady-state AUC at half maximum ALP effect (umol*h/L)")  # Table 2 'AUC tau,ss sum,50 (umol.h/L)' 24.1 [RSE 18.1%]
    hill_alp <- fixed(1.00)
    label("Hill coefficient of the ALP Emax model (unitless)")     # Table 2 'Hill' 1.00, fixed (Results 3.2.2 'The Hill coefficient of the Emax model was fixed to 1')
    slope_tbili <- -0.0100
    label("Linear slope of the drug effect on total bilirubin production (L/(umol*h))")   # Table 2 'Slope (L/[umol.h])' -0.0100 [RSE 51.3%]

    # ==================================================================
    # PLACEBO EFFECT -- Table 2
    # Fractional multiplier on the production rate, applied from the first
    # dose onwards (Supplementary Datafile S6 $DES, IF(TIME.GT.0)).
    # ==================================================================
    lpbo_alp <- log(0.992)
    label("Placebo multiplier on alkaline phosphatase production (fraction)")  # Table 2 'Placebo effect ALP' 0.992 [RSE 1.64%]
    lpbo_tbili <- log(0.958)
    label("Placebo multiplier on total bilirubin production (fraction)")       # Table 2 'Placebo effect TB' 0.958 [RSE 20.8%]

    # ==================================================================
    # COVARIATE EFFECTS -- Table 2
    # ==================================================================
    e_hepimp_rbase_alp <- 0.277
    label("Proportional increase in baseline ALP with any NCI hepatic impairment (fraction)")  # Table 2 'NCI hepatic impairment score on BASE ALP (proportional increase)' 0.277 [RSE 28.8%]
    e_lsm_rbase_tbili <- 0.333
    label("Power exponent of baseline liver stiffness on baseline total bilirubin (unitless)") # Table 2 'Liver stiffness on BASE TB (power)' 0.333 [RSE 20.6%]

    # ==================================================================
    # INTER-INDIVIDUAL VARIABILITY -- Table 2
    # Table 2 reports IIV on the approximate SD (CV) scale; variances are its
    # square and off-diagonals are correlation x SD x SD.
    # ==================================================================
    # Correlated baselines. Variances 0.320^2 and 0.384^2; covariance
    # 0.148 * 0.320 * 0.384 = 0.01818624.
    etalrbase_alp + etalrbase_tbili ~ c(
      0.10240000,
      0.01818624, 0.14745600
    )   # Table 2 'IIV BASE ALP (CV)' 0.320, 'IIV BASE TB (CV)' 0.384, 'Correlation between the IIV BASE ALP and BASE TB' 0.148

    # Correlated logit-Emax and ALP placebo effect. Table 2 prints the Emax
    # IIV as 0.231 on the BACK-TRANSFORMED Emax scale; the model's random
    # effect acts on the logit scale, where Supplementary Datafile S6
    # $OMEGA BLOCK(2) gives the variance 0.71883 (SD 0.8478). The two agree
    # by the delta method: 0.8478 * 0.731 * (1 - 0.731) / 0.731 = 0.228.
    # That stream block is the final one - its off-diagonal -0.087288 with
    # these two diagonals returns the printed correlation -0.979 exactly.
    # Covariance below is -0.979 * sqrt(0.71883) * 0.105 = -0.08715.
    etalogitemax_alp + etalpbo_alp ~ c(
      0.71883000,
      -0.08715450, 0.01102500
    )   # Datafile S6 $OMEGA BLOCK(2) 0.71883; Table 2 'IIV placebo effect ALP (CV)' 0.105 and 'Correlation between IIV on E max and placebo effect ALP' -0.979

    etalpbo_tbili ~ 0.952576   # Table 2 'IIV placebo effect TB (CV)' 0.976; 0.976^2
    etapropSd_alp ~ 0.864900   # Table 2 'IIV proportional RUV ALP (CV)' 0.930; 0.930^2
    etapropSd_tbili ~ 0.029584 # Table 2 'IIV proportional RUV TB (CV)' 0.172; 0.172^2

    # Table 2 carries no IIV on the turnover half-lives, on EC50 or on the
    # total-bilirubin slope; Supplementary Datafile S6 fixes those four
    # omegas to 0, so no eta is declared for them.

    # ==================================================================
    # RESIDUAL ERROR -- Table 2
    # Combined additive plus proportional per biomarker, with IIV on the
    # proportional component (Datafile S6 $ERROR,
    # Y = IPRED * (1 + EPS*EXP(ETA)) + EPS_add).
    # ==================================================================
    propSd_alp <- 0.0659
    label("Proportional residual error for alkaline phosphatase (fraction)")  # Table 2 'Proportional RUV ALP (CV)' 0.0659 [RSE 11.1%]
    addSd_alp <- 14.8
    label("Additive residual error for alkaline phosphatase (U/L)")           # Table 2 'Additive RUV ALP (SD, U/L)' 14.8 [RSE 4.65%]
    propSd_tbili <- 0.167
    label("Proportional residual error for total bilirubin (fraction)")       # Table 2 'Proportional RUV TB (CV)' 0.167 [RSE 4.59%]
    addSd_tbili <- 0.484
    label("Additive residual error for total bilirubin (umol/L)")             # Table 2 'Additive RUV TB (SD, umol/L)' 0.484 [RSE 26.7%]
  })

  model({
    # ------------------------------------------------------------------
    # 1. Drug exposure metric. Elafibranor and GFT1007 are equipotent, so
    #    only their summed steady-state AUC enters (Methods 2.3.2).
    #    Both are 0 on placebo, which makes both drug terms vanish exactly.
    # ------------------------------------------------------------------
    aucSum <- AUC_ELA + AUC_GFT1007

    # ------------------------------------------------------------------
    # 2. Baselines and turnover rates
    # ------------------------------------------------------------------
    rbase_alp <- exp(lrbase_alp + etalrbase_alp) *
      (1 + e_hepimp_rbase_alp * HEPIMP)
    rbase_tbili <- exp(lrbase_tbili + etalrbase_tbili) *
      (LSM / 8.1)^e_lsm_rbase_tbili

    kout_alp <- log(2) / exp(lthalfrec_alp)
    kout_tbili <- log(2) / exp(lthalfrec_tbili)
    kin_alp <- kout_alp * rbase_alp
    kin_tbili <- kout_tbili * rbase_tbili

    # ------------------------------------------------------------------
    # 3. Drug effects on production
    #    Emax is constrained to (-1, 0) through a logit transform of its
    #    absolute value, which is where its random effect acts
    #    (Datafile S6 $PK, 'Code for logit transformed EMAX').
    # ------------------------------------------------------------------
    posemax_alp <- -emax_alp
    emax_alp_i <- -expit(logit(posemax_alp) + etalogitemax_alp)
    ec50_alp <- exp(lec50_alp)

    drug_alp <- 1 + emax_alp_i * aucSum^hill_alp /
      (ec50_alp^hill_alp + aucSum^hill_alp)
    drug_tbili <- 1 + slope_tbili * aucSum

    # ------------------------------------------------------------------
    # 4. Placebo effect, switched on from the first dose onwards
    # ------------------------------------------------------------------
    pbo_alp <- exp(lpbo_alp + etalpbo_alp)
    pbo_tbili <- exp(lpbo_tbili + etalpbo_tbili)
    pbo_alp_t <- 1 + (pbo_alp - 1) * (t > 0)
    pbo_tbili_t <- 1 + (pbo_tbili - 1) * (t > 0)

    # ------------------------------------------------------------------
    # 5. Indirect-response system (Datafile S6 $DES)
    # ------------------------------------------------------------------
    d/dt(alp) <- kin_alp * pbo_alp_t * drug_alp - kout_alp * alp
    d/dt(tbili) <- kin_tbili * pbo_tbili_t * drug_tbili - kout_tbili * tbili

    alp(0) <- rbase_alp
    tbili(0) <- rbase_tbili

    # ------------------------------------------------------------------
    # 6. Observations and residual error
    # ------------------------------------------------------------------
    propSd_alp_i <- propSd_alp * exp(etapropSd_alp)
    propSd_tbili_i <- propSd_tbili * exp(etapropSd_tbili)

    alp ~ add(addSd_alp) + prop(propSd_alp_i)
    tbili ~ add(addSd_tbili) + prop(propSd_tbili_i)
  })
}
