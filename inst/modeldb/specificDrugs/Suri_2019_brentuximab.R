Suri_2019_brentuximab <- function() {
  description <- "Coupled population PK model for the brentuximab vedotin antibody-drug conjugate (ADC) and its released payload monomethyl auristatin E (MMAE) in 661 adults with previously untreated stage III or IV classical Hodgkin lymphoma receiving brentuximab vedotin 1.2 mg/kg every 2 weeks with doxorubicin, vinblastine and dacarbazine (A+AVD) in the phase III ECHELON-1 study (Suri 2019). ADC is a linear 3-compartment model with zero-order input and first-order elimination; MMAE is a 2-compartment model with first-order elimination, formed through an intermediate lag compartment from (a) a one-time target-binding flux KD*Target*ADC (initial Target = 1, irreversibly depleted) and (b) a direct proteolytic flux FM*(1 - exp(-ALFM*tad))*K10*ADC with FM fixed to 1. The direct-flux time function is not printed in the source; the form used here is the one that reproduces the paper's own cycle-3 MMAE AUC and 49% dose-1-to-dose-5 decline (see the vignette). Amounts are in umol and concentrations in umol/L."
  reference <- "Suri A, Mould DR, Song G, Collins GP, Endres CJ, Gomez-Navarro J, Venkatakrishnan K. Population Pharmacokinetic Modeling and Exposure-Response Assessment for the Antibody-Drug Conjugate Brentuximab Vedotin in Hodgkin's Lymphoma in the Phase III ECHELON-1 Study. Clin Pharmacol Ther. 2019;106(6):1268-1279. doi:10.1002/cpt.1530. PMCID PMC6896233."
  vignette <- "Suri_2019_brentuximab"
  paper_specific_compartments <- c("lag")

  units <- list(
    time = "h",
    dosing = "umol",
    concentration = "umol/L"
  )

  compartmentData <- list(
    central = list(analyte = "brentuximab vedotin (ADC)", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "brentuximab vedotin (ADC)", units = "umol", specimen = "tissue", verified = TRUE),
    peripheral2 = list(analyte = "brentuximab vedotin (ADC)", units = "umol", specimen = "tissue", verified = TRUE),
    target = list(
      analyte = "hypothetical ADC binding target (unitless fraction remaining)",
      units = "(unitless)",
      specimen = "not applicable",
      verified = TRUE
    ),
    lag = list(
      analyte = "monomethyl auristatin E (in transit to plasma)",
      units = "umol",
      specimen = "not applicable",
      verified = TRUE
    ),
    central_mmae = list(analyte = "monomethyl auristatin E", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1_mmae = list(analyte = "monomethyl auristatin E", units = "umol", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effects normalised to the population mean 1.8 m^2 (Suri 2019 Table 1 mean, printed to one decimal place; the supplement states continuous covariates were 'normalized for the population mean'). Exponents: ADC CL 1.1, ADC Vc 0.893, ADC V3 1.47, MMAE CL 1.04. Cohort range 1.3-2.8 m^2. The one patient missing weight / BSA was assigned the population median.",
      source_name = "BSA"
    ),
    ALB = list(
      description = "Serum albumin concentration",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effects normalised to the population mean 39.1 g/L (Suri 2019 Table 1). Exponents: ADC CL -0.477, MMAE CL 0.0275. Cohort range 17-53 g/L.",
      source_name = "Albumin"
    ),
    CRCL = list(
      description = "Creatinine clearance by the Cockcroft-Gault equation (raw mL/min, NOT BSA-normalised)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on MMAE CL only (exponent 0.125), normalised to the population mean 134.1 mL/min (Suri 2019 Table 1, footnote b: Cockcroft-Gault). Cohort range 29.2-476.7 mL/min.",
      source_name = "CrCl"
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Categorical effect on ADC central volume only: vc *= 0.934^SEXF. The supplement's categorical form fixes theta = 1 for the reference subgroup coded 0 and gives males as the example reference; the Results statement that Vc is ~20% lower in females is consistent with 0.934 combined with the lower BSA of females under the 0.893 BSA exponent. 285 of 661 patients (43.1%) were female.",
      source_name = "Sex"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened; 'no discernible effect' on ADC or MMAE PK (Suri 2019 Results, data not shown)."
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator (1 = Asian, 0 = non-Asian)",
      units = "(binary)",
      type = "binary",
      notes = "Screened; no discernible effect on ADC or MMAE PK (Suri 2019 Results and Discussion)."
    ),
    ADA_POS = list(
      description = "Anti-drug antibody positive (including neutralising ADA)",
      units = "(binary)",
      type = "binary",
      notes = "Screened; immunogenicity status was not identified as a covariate on ADC clearance and had no effect on MMAE PK (Suri 2019 Results). Contrast Suri_2018_brentuximab, where ADA status was retained on ADC CL."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Reported in Suri 2019 Table 1 (mean 66.2 umol/L) but the final MMAE model carries Cockcroft-Gault creatinine clearance (CRCL) instead."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Reported in Suri 2019 Table 1 (mean 7.1 umol/L) but not retained in either final model (contrast Suri_2018_brentuximab, where bilirubin was on MMAE CL)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 661L,
    n_studies = 1L,
    age_range = "18-82 years (mean 38.7, SD 15.8)",
    weight_range = "40.8-165.5 kg (mean 73.5, SD 18.0)",
    bsa_range = "1.3-2.8 m^2 (mean 1.8, SD 0.3)",
    sex_female_pct = 43.1,
    race_ethnicity = c(White = 84.3, Black = 3.0, Asian = 8.5, Other = 2.7, NotReported = 1.5),
    disease_state = "Previously untreated, histologically confirmed stage III or IV classical Hodgkin lymphoma (ECHELON-1, NCT01712490, A+AVD arm).",
    dose_range = "Brentuximab vedotin 1.2 mg/kg as a 30-minute IV infusion on days 1 and 15 of 28-day cycles for up to 6 cycles, capped at 120 mg for patients above 100 kg, given within about 1 hour after doxorubicin, vinblastine and dacarbazine.",
    regions = "International (Americas, Europe, Asia).",
    study_phase = "Phase III",
    n_observations = "33,164 concentration records (16,536 ADC, 16,628 MMAE) and 7,209 dosing records; 347 post-dose BLQ records excluded (42 ADC, 305 MMAE).",
    albumin = "39.1 g/L mean (SD 5.3), range 17-53",
    renal_function = "Cockcroft-Gault CrCl 134.1 mL/min mean (SD 45.4), range 29.2-476.7",
    reference_subject = "BSA 1.8 m^2, ALB 39.1 g/L, CRCL 134.1 mL/min, male (SEXF = 0).",
    notes = "Demographics from Suri 2019 Table 1 (N = 661). Sparse sampling (pre-dose and end of infusion on days 1 and 15 of every cycle, plus 24 and 48 h in cycles 1 and 3) with an intensive subset of 59 patients (Suri 2019 supplement Table S4)."
  )

  ini({
    # ADC structural parameters (Suri 2019 Table 2, ADC columns).
    lcl  <- log(0.0615); label("ADC clearance (L/h)")                                     # Suri 2019 Table 2: 0.0615 (1.0% SE)
    lvc  <- log(3.58);   label("ADC central volume (L)")                                   # Suri 2019 Table 2: 3.58 (0.9% SE)
    lq   <- log(0.113);  label("ADC intercompartmental clearance to peripheral 1 (L/h)")  # Suri 2019 Table 2: Q2 0.113 (3.0% SE)
    lvp  <- log(3.26);   label("ADC peripheral volume 1 (L)")                              # Suri 2019 Table 2: V2 3.26 (1.9% SE)
    lq2  <- log(0.0239); label("ADC intercompartmental clearance to peripheral 2 (L/h)")  # Suri 2019 Table 2: Q3 0.0239 (2.3% SE)
    lvp2 <- log(15.7);   label("ADC peripheral volume 2 (L)")                              # Suri 2019 Table 2: V3 15.7 (4.0% SE)

    # ADC covariate effects (Suri 2019 Table 2; power form on the covariate
    # divided by its population mean, categorical form theta^SEXF, per the
    # supplement's covariate-model section).
    e_alb_cl  <- -0.477; label("Power exponent of (ALB / 39.1) on ADC CL (unitless)")   # Suri 2019 Table 2: 'Albumin on clearance' -0.477 (2.2% SE)
    e_bsa_cl  <- 1.1;    label("Power exponent of (BSA / 1.8) on ADC CL (unitless)")    # Suri 2019 Table 2: 'BSA on clearance' 1.1 (4.9% SE)
    e_bsa_vc  <- 0.893;  label("Power exponent of (BSA / 1.8) on ADC Vc (unitless)")    # Suri 2019 Table 2: 'BSA on central volume' 0.893 (6.3% SE)
    e_sexf_vc <- 0.934;  label("Multiplier on ADC Vc for female sex, applied as e_sexf_vc^SEXF (unitless)") # Suri 2019 Table 2: 'Sex on central volume' 0.934 (1.4% SE)
    e_bsa_vp2 <- 1.47;   label("Power exponent of (BSA / 1.8) on ADC peripheral volume 2 (unitless)")     # Suri 2019 Table 2: 'BSA on peripheral volume 2' 1.47 (8.5% SE)

    # MMAE structural parameters (Suri 2019 Table 2, MMAE columns). Apparent
    # values because the fraction metabolised FM is fixed to 1.
    lcl_mmae   <- log(1.45);   label("Apparent MMAE clearance (L/h)")                     # Suri 2019 Table 2: 1.45 (0.2% SE)
    lvc_mmae   <- log(35.5);   label("Apparent MMAE central volume (L)")                  # Suri 2019 Table 2: 35.5 (0.3% SE)
    lq_mmae    <- log(13.2);   label("Apparent MMAE intercompartmental clearance (L/h)")  # Suri 2019 Table 2: Q2 13.2 (0.3% SE)
    lvp_mmae   <- log(17.7);   label("Apparent MMAE peripheral volume (L)")               # Suri 2019 Table 2: V2 17.7 (0.3% SE)
    lkd_mmae   <- log(0.0376); label("ADC target-binding rate constant forming MMAE (1/h per umol of ADC)") # Suri 2019 Table 2: KD 0.0376 1/hour (0.3% SE)
    lalfm_mmae <- log(2.35);   label("Rate constant of the time-after-dose function on direct ADC to MMAE conversion (1/h)") # Suri 2019 Table 2: ALFM 2.35 1/hour (0.2% SE)
    lklag_mmae <- log(4.51);   label("Lag-compartment transfer rate constant (1/h)")      # Suri 2019 Table 2: Klag 4.51 1/hour (0.3% SE)
    # Suri 2019 Table 2: 'Fraction metabolized 1 FIX' (sparse data could not
    # inform it); multiplies the direct flux in model().
    fm_mmae    <- fixed(1);    label("Fraction of eliminated ADC converted directly to MMAE (fraction)")

    # MMAE covariate effects (Suri 2019 Table 2).
    e_bsa_cl_mmae  <- 1.04;   label("Power exponent of (BSA / 1.8) on MMAE CL (unitless)")      # Suri 2019 Table 2: 'BSA on clearance' 1.04 (1% SE)
    e_crcl_cl_mmae <- 0.125;  label("Power exponent of (CRCL / 134.1) on MMAE CL (unitless)")   # Suri 2019 Table 2: 'Creatinine clearance on clearance' 0.125 (3.3% SE)
    e_alb_cl_mmae  <- 0.0275; label("Power exponent of (ALB / 39.1) on MMAE CL (unitless)")     # Suri 2019 Table 2: 'Albumin concentration on clearance' 0.0275 (2.8% SE)

    # IIV. The supplement reports each IIV as sqrt(omega^2) and calls it an
    # approximate CV, so omega^2 = (CV/100)^2 (not log(CV^2 + 1)).
    etalcl  ~ 0.039204  # Suri 2019 Table 2: ADC CL 19.8% CV -> 0.198^2
    etalvc  ~ 0.0196    # Suri 2019 Table 2: ADC Vc 14.0% CV -> 0.140^2
    etalvp  ~ 0.065025  # Suri 2019 Table 2: ADC V2 25.5% CV -> 0.255^2
    etalq2  ~ 0.171396  # Suri 2019 Table 2: ADC Q3 41.4% CV -> 0.414^2
    etalvp2 ~ 0.603729  # Suri 2019 Table 2: ADC V3 77.7% CV -> 0.777^2

    etalcl_mmae ~ 0.148996  # Suri 2019 Table 2: MMAE CL 38.6% CV -> 0.386^2
    etalvc_mmae ~ 0.470596  # Suri 2019 Table 2: MMAE Vc 68.6% CV -> 0.686^2
    etalkd_mmae ~ 1.0201    # Suri 2019 Table 2: MMAE KD 101.0% CV -> 1.010^2

    # Residual error: log-transform-both-sides, homoscedastic additive on the
    # log scale (Suri 2019 supplement), encoded as lnorm().
    expSd      <- 0.181; label("Log-scale residual SD for ADC Cc (log units)")        # Suri 2019 Table 2: ADC residual variability 18.1% CV (0.4% SE)
    expSd_mmae <- 0.395; label("Log-scale residual SD for MMAE Cc_mmae (log units)")  # Suri 2019 Table 2: MMAE residual variability 39.5% CV (0.3% SE)
  })

  model({
    # 1. Covariates normalised to the Suri 2019 Table 1 population means.
    nbsa  <- BSA / 1.8
    nalb  <- ALB / 39.1
    ncrcl <- CRCL / 134.1

    # 2. Individual ADC parameters
    cl_adc <- exp(lcl + etalcl) * nbsa^e_bsa_cl * nalb^e_alb_cl
    v1_adc <- exp(lvc + etalvc) * nbsa^e_bsa_vc * e_sexf_vc^SEXF
    q2_adc <- exp(lq)
    v2_adc <- exp(lvp + etalvp)
    q3_adc <- exp(lq2 + etalq2)
    v3_adc <- exp(lvp2 + etalvp2) * nbsa^e_bsa_vp2

    # 3. Individual MMAE parameters
    cl_mmae <- exp(lcl_mmae + etalcl_mmae) * nbsa^e_bsa_cl_mmae * ncrcl^e_crcl_cl_mmae * nalb^e_alb_cl_mmae
    vc_mmae <- exp(lvc_mmae + etalvc_mmae)
    qm      <- exp(lq_mmae)
    vp_mmae <- exp(lvp_mmae)
    kd      <- exp(lkd_mmae + etalkd_mmae)
    alfm    <- exp(lalfm_mmae)
    klag    <- exp(lklag_mmae)

    # 4. Micro-constants
    k10 <- cl_adc / v1_adc
    k12 <- q2_adc / v1_adc
    k21 <- q2_adc / v2_adc
    k13 <- q3_adc / v1_adc
    k31 <- q3_adc / v3_adc
    k40 <- cl_mmae / vc_mmae
    k45 <- qm / vc_mmae
    k54 <- qm / vp_mmae

    # 5. Time-after-dose function on direct ADC -> MMAE conversion. Suri 2019
    # Figure 1b names FM and ALFM but prints no equation. The literal reading
    # of the prose, FM * exp(-ALFM * tad), with ALFM = 2.35 1/h shuts direct
    # conversion off within about an hour of each dose and predicts a cycle-3
    # MMAE AUC ~40-fold below the paper's own simulated 10-11 ng*day/mL. The
    # complementary form below reproduces both the paper's cycle-3 AUC and its
    # 49% dose-1-to-dose-5 decline without tuning; see the vignette.
    # tad() is assigned first: a fixed() theta multiplying an expression that
    # calls tad() inline fails rxode2's mu-reference parse.
    tafd <- tad()
    fmt  <- fm_mmae * (1 - exp(-alfm * tafd))

    # 6. ODE system (Suri 2019 Figure 1b). The target pool starts at 1 and is
    # irreversibly depleted by ADC; the binding flux does not remove ADC from
    # central, matching the Mould-lab ADC model also used by Suri_2018 and
    # Zhou_2025 brentuximab.
    target(0) <- 1
    d/dt(central)          <- -(k10 + k12 + k13) * central + k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1)      <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2)      <-  k13 * central - k31 * peripheral2
    d/dt(target)           <- -kd * target * central
    d/dt(lag)              <-  kd * target * central + fmt * k10 * central - klag * lag
    d/dt(central_mmae)     <-  klag * lag - (k40 + k45) * central_mmae + k54 * peripheral1_mmae
    d/dt(peripheral1_mmae) <-  k45 * central_mmae - k54 * peripheral1_mmae

    # 7. Observations. Suri 2019 Figure 3a and 3c plot ADC and MMAE in uM, so
    # amounts are umol (dose_umol = dose_mg / 153.4 for the ~153.4 kDa ADC)
    # and concentrations umol/L. Cc * 153.4 gives ug/mL; Cc_mmae * 718 gives
    # ng/mL (MMAE 718 g/mol).
    Cc      <- central / v1_adc
    Cc_mmae <- central_mmae / vc_mmae

    Cc      ~ lnorm(expSd)
    Cc_mmae ~ lnorm(expSd_mmae)
  })
}
