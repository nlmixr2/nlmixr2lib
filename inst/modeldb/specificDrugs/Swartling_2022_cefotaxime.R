Swartling_2022_cefotaxime <- function() {
  description <- "Two-compartment intravenous population PK model for cefotaxime in critically ill adult ICU patients not on renal replacement therapy (Swartling 2022; ACCIS study, 51 patients at seven Swedish ICUs including one burn unit). Clearance and intercompartmental clearance scale allometrically with body weight at ICU admission (fixed exponent 0.75, reference 92 kg) and both volumes scale linearly with it. Clearance increases linearly with Cockcroft-Gault estimated creatinine clearance (centred at 94 mL/min) up to 120 mL/min and is flat above it. A single random effect is shared by clearance and central volume, scaled by an estimated factor on the central volume (IIV 49% on CL, 64% on Vc); residual error is proportional."
  reference <- "Swartling M, Smekal AK, Furebring M, Lipcsey M, Jonsson S, Nielsen EI. Population pharmacokinetics of cefotaxime in intensive care patients. Eur J Clin Pharmacol. 2022;78(2):251-258. doi:10.1007/s00228-021-03218-6"
  vignette <- "Swartling_2022_cefotaxime"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Cefotaxime is given as a 5-min IV infusion, so the dose
  # enters `central` directly and there is no depot state. `central` is
  # verified: Swartling 2022 'Cefotaxime analysis' states that TOTAL SERUM
  # cefotaxime was measured by LC-MS/MS, and Online Resource 2 (the final
  # NONMEM control stream) uses ADVAN3 TRANS4 with `S1 = V1`, so the observed
  # compartment is the central one. `peripheral1` is left unverified because
  # the paper does not state what matrix the peripheral distribution
  # compartment represents.
  compartmentData <- list(
    central = list(analyte = "cefotaxime", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "cefotaxime", units = "mg", specimen = "serum", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight at ICU admission",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Swartling 2022 enters body weight AT ICU ADMISSION (control-stream column BWAI) as a time-independent covariate, added before the other covariates on allometric principles ('Pharmacokinetic modelling' section): CL and Q scale as (WT/92)^0.75 and Vc and Vp as (WT/92). The reference 92 kg is the Table 1 cohort median (IQR 81-102, range 55-124). The paper chose ICU-admission weight deliberately because weight gained from resuscitation fluids may reflect volume of distribution. Missing values (7 patients) were imputed by linear regression on the weight before ICU admission (Table 1 median 87 kg, range 55-123), which is a separate column (BWBI) not used by the final model.",
      source_name = "BWAI"
    ),
    CRCL = list(
      description = "Estimated creatinine clearance by the Cockcroft-Gault formula computed with body weight at ICU admission; raw mL/min, NOT BSA-normalized",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Swartling 2022 tests eCLcr (Cockcroft-Gault, reference 18) as a TIME-DEPENDENT covariate on CL, carried forward to the next observation without interpolation (first value carried backward). Online Resource 2 names the column COGBWA (Cockcroft-Gault on body weight at ICU admission); the Table 1 footnote confirms the formula used body weight at ICU admission. Centring value 94 mL/min is the Table 1 median over treatment days 1-3 (IQR 65-138, range 5-258; day 1 median 85, day 2 97, day 3 105). The final model is piecewise linear with a slope fixed to zero above 120 mL/min, which the paper states is equivalent to truncating eCLcr at 120 mL/min. Creatinine was drawn routinely at 6 a.m., so the value paired with a cefotaxime sample may lag it by up to 24 h. Supply the raw (un-normalized) Cockcroft-Gault value; a BSA-normalized eGFR would rescale the effect.",
      source_name = "COGBWA"
    )
  )

  # Screened in the Swartling 2022 covariate analysis (forward inclusion
  # p < 0.05, backward elimination p < 0.001) but NOT retained, so these are
  # documentation only and are not referenced in model().
  covariatesDataExcluded <- list(
    DIS_BURN_RECENT = list(
      description = "Treated in the burn-unit ICU (1) versus a general ICU (0)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (general ICU)",
      notes = "Tested as a time-independent categorical covariate on the volume parameters because V in burn patients may differ from other ICU patients; not retained after eCLcr entered on CL ('Inclusion of the other covariates (e.g. burn patients) after inclusion of eCLcr as a covariate did not significantly improve the model'). 18 of 51 patients (35%) were treated in the burn unit. Control-stream column BURN.",
      source_name = "BURN"
    ),
    TRTDAY = list(
      description = "Day of antibiotic treatment (1, 2 or 3), categorical",
      units = "(day)",
      type = "categorical",
      reference_category = NULL,
      notes = "Tested as a time-dependent categorical covariate on the volume parameters and CL. The Discussion reports it was statistically significant (p < 0.05) before eCLcr entered on CL, but it was not retained in the final model. Only 34 of 51 patients had day-3 samples.",
      source_name = "day of treatment"
    ),
    SAPS3 = list(
      description = "Simplified Acute Physiology Score 3 at ICU admission",
      units = "(score)",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested as a continuous covariate on CL to capture non-renal factors; not retained. Table 1 median 55 (IQR 47-63, range 25-89). Two missing values were imputed with the population mean.",
      source_name = "SAPS3"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 51L,
    n_studies = 1L,
    n_sites = 7L,
    n_concentrations = 263L,
    age_range = "23-90 years",
    age_median = "64 years (IQR 50-73)",
    weight_range = "55-124 kg (body weight at ICU admission)",
    weight_median = "92 kg (IQR 81-102) at ICU admission",
    sex_female_pct = 35,
    disease_state = "Critically ill adults (>18 years) in the intensive care unit treated with cefotaxime for a proven or suspected infection, 18 in a burn-unit ICU and 33 in general ICUs. Most common infections: lower respiratory tract (59%), skin and soft tissue (25%), urinary tract (8%). SAPS3 at admission median 55 (range 25-89). Pregnant patients, patients with treatment restrictions and patients on renal replacement therapy were excluded.",
    dose_range = "Cefotaxime 1000-3000 mg as a 5-min IV infusion, 2-6 times daily, at physician discretion; most common 1000 mg t.i.d. followed by 2000 mg t.i.d. An extra dose in the middle of the first dosing interval was given in four cases.",
    regions = "Sweden (seven ICUs in five hospitals)",
    renal_function = "Cockcroft-Gault eCLcr median 94 mL/min over treatment days 1-3 (IQR 65-138, range 5-258).",
    notes = "Sub-study of the prospective observational multi-centre ACCIS study (Antibiotic Concentrations in Critical Ill ICU Patients in Sweden; ACTRN12616000167460), December 2015 to July 2017. Two samples per day (mid-interval and just before the next dose) for up to three consecutive days from the first day of treatment; median 6 samples per patient (range 2-6). Total serum cefotaxime by LC-MS/MS, quantification range 0.50-50 mg/L; 15 samples (6%) below the LLOQ were set to LLOQ/2. Four erroneously high supposed troughs were excluded. NONMEM 7.4, FOCE-INTERACTION. Only cefotaxime, not desacetylcefotaxime, was measured."
  )

  ini({
    # Structural parameters, final model (Swartling 2022 Table 2 'Estimate'
    # column). The typical individual has body weight at ICU admission 92 kg
    # and eCLcr 94 mL/min (Abstract; Table 2 footnotes a-d).
    lcl <- log(11.1); label("Clearance at WT = 92 kg and CRCL = 94 mL/min (L/h)")  # Swartling 2022 Table 2: CL 11.1 L/h (RSE 8.2%)
    lvc <- log(5.13); label("Central volume of distribution at WT = 92 kg (L)")    # Swartling 2022 Table 2: Vc 5.13 L (RSE 28%)
    lvp <- log(18.2); label("Peripheral volume of distribution at WT = 92 kg (L)") # Swartling 2022 Table 2: Vp 18.2 L (RSE 12%)
    lq <- log(14.5); label("Intercompartmental clearance at WT = 92 kg (L/h)")     # Swartling 2022 Table 2: Q 14.5 L/h (RSE 19%)

    # Allometric body-weight exponents. Hard-coded (not estimated) in Online
    # Resource 2: TVCL = THETA(1)*(BWAI/92)**0.75, TVV1 = THETA(2)*(BWAI/92),
    # TVV2 = THETA(3)*(BWAI/92), TVQ = THETA(4)*(BWAI/92)**0.75.
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent of WT/92 on CL and Q (unitless)") # Swartling 2022 Table 2 footnotes a and d; Online Resource 2 $PK
    e_wt_vc_vp <- fixed(1); label("Allometric exponent of WT/92 on Vc and Vp (unitless)")  # Swartling 2022 Table 2 footnotes b and c; Online Resource 2 $PK

    # eCLcr effect on CL. Table 2 footnote a and Online Resource 2:
    #   CLCOV = THETA(7)*((COGBWA-94)/1000) for COGBWA <= 120, and
    #   CLCOV = THETA(7)*((120-94)/1000) above it;
    #   TVCL = THETA(1)*(BWAI/92)**0.75*(1+CLCOV)
    # so the coefficient is a fractional change in CL per 1000 mL/min.
    e_crcl_cl <- 6.67; label("Fractional change in CL per 1000 mL/min of CRCL below the cap (unitless)") # Swartling 2022 Table 2: theta_cov (eCLcr <= 120) 6.67 (RSE 17%)
    crcl_cap <- fixed(120); label("Upper cap applied to CRCL before the CL effect (mL/min)")              # Swartling 2022 Results: slope fixed to zero above 120 mL/min, 'equals truncating eCLcr at high values'; Table 2 footnote a

    # Shared random effect. Online Resource 2: CL = TVCL*EXP(ETA(1)) and
    # V1 = TVV1*EXP(ETA(1)*THETA(6)), a single $OMEGA; the Results state IIV
    # on CL and Vc was highly positively correlated, so one eta with a scaling
    # factor on Vc was used.
    scale_etalvc <- 1.31; label("Scaling factor applied to the shared CL eta for Vc (unitless)") # Swartling 2022 Table 2: f_CL,Vc 1.31 (RSE 39%)

    # IIV on CL. Table 2 prints 'IIV CL (% CV)' 49 and 'IIV Vc (% CV)' 64,
    # with footnote b stating the Vc CV was 'derived as CV for CL times the
    # estimated f_CL,Vc'. That multiplication (49 x 1.31 = 64.2) is exact only
    # when the %CV is the standard deviation of eta, i.e. sqrt(omega): under
    # the exact log-normal form the Vc CV would be
    # sqrt(exp(1.31^2 * log(1 + 0.49^2)) - 1) = 66%. The proportional-error
    # row in the same table is likewise the SD directly (THETA(5) with $SIGMA
    # 1 FIX). So omega^2 = 0.49^2.
    etalcl ~ 0.2401 # 0.49^2; Swartling 2022 Table 2: IIV CL 49% CV (RSE 15%, shrinkage 0.8%)

    # Residual error. Online Resource 2: W = IPRED*THETA(5); Y = IPRED +
    # W*EPS(1) with $SIGMA 1 FIX, so THETA(5) is the proportional SD.
    propSd <- 0.333; label("Proportional residual error (fraction)") # Swartling 2022 Table 2: proportional residual error 33.3% CV (RSE 5.9%, shrinkage 7.8%)
  })
  model({
    # Covariate terms (Swartling 2022 Table 2 footnotes a-d; Online Resource 2).
    crcl_capped <- min(CRCL, crcl_cap)
    crcl_eff <- 1 + e_crcl_cl * (crcl_capped - 94) / 1000
    wt_cl_q <- (WT / 92)^e_wt_cl_q
    wt_vc_vp <- (WT / 92)^e_wt_vc_vp

    # Individual parameters. One eta is shared by CL and Vc, scaled on Vc.
    cl <- exp(lcl + etalcl) * wt_cl_q * crcl_eff
    vc <- exp(lvc + scale_etalvc * etalcl) * wt_vc_vp
    vp <- exp(lvp) * wt_vc_vp
    q <- exp(lq) * wt_cl_q

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg, volumes in L, so central/vc is mg/L: the TOTAL serum
    # cefotaxime concentration measured by the LC-MS/MS assay.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
