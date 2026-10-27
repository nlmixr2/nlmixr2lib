Baaz_2023_panitumumab_folfox <- function() {
  description <- paste(
    "Joint tumor growth inhibition (TGI) and non-target-progression (NTP)",
    "time-to-event model for FOLFOX with or without panitumumab in RAS exon 2",
    "wild-type metastatic colorectal cancer (Baaz 2023, PRIME study,",
    "n = 248). The sum of the largest target-lesion diameters (tumor_size, mm)",
    "grows exponentially at kge and is reduced by a dose-driven kill term,",
    "(kdrug_folfox * DOSE_FOLFOX_MGM2 + kdrug_panitumumab *",
    "DOSE_PANITUMUMAB_MGKG) * exp(-lambda * t), where the doses are the total",
    "dose per 14-day cycle and lambda is a resistance rate shared by both",
    "drugs. There is no PK model. The NTP hazard is linked to the tumor's",
    "rate of change, h = max(e_dsld_haz * d(tumor_size)/dt + base_haz, 0), and",
    "its cumulative hazard gives the NTP survivor function sur. Time is in",
    "months from the start of treatment. The separately calibrated Weibull",
    "dropout model is modellib('Baaz_2023_panitumumab_folfox_dropout').",
    sep = " "
  )
  reference <- paste(
    "Baaz M, Cardilin T, Jirstrand M. Model-based prediction of",
    "progression-free survival for combination therapies in oncology.",
    "CPT Pharmacometrics Syst Pharmacol. 2023;12(9):1227-1237.",
    "doi:10.1002/psp4.13003.",
    "Parameter estimates from Table 1; model code from the Monolix project in",
    "Data S1.",
    sep = " "
  )
  vignette <- "Baaz_2023_panitumumab_folfox"

  units <- list(
    time = "month",
    dosing = "n/a (no dose events; the regimen enters through the per-cycle dose covariates DOSE_FOLFOX_MGM2 and DOSE_PANITUMUMAB_MGKG)",
    concentration = "mm (the observed variable is the sum of the largest target-lesion diameters, not a drug concentration)"
  )

  covariateData <- list(
    DOSE_FOLFOX_MGM2 = list(
      description = "Total FOLFOX dose per 14-day treatment cycle, the summed mg/m^2 amounts of oxaliplatin, leucovorin and fluorouracil; 0 when FOLFOX is not given",
      units = "mg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Baaz 2023 drives the PD model with 'the total dose per cycle (14",
        "days), Di', combining the three FOLFOX drugs' doses, but does not",
        "print the resulting number. The value that reproduces the paper's",
        "own SLD VPC (Figure S1) is 1285 mg/m^2, the sum of the four",
        "per-administration amounts the Methods list for FOLFOX4:",
        "oxaliplatin 85 + leucovorin 200 + 5-FU bolus 400 + 5-FU infusion",
        "600 mg/m^2. The full two-day cycle total (2485 mg/m^2) doubles the",
        "FOLFOX kill rate and predicts a median 2-month SLD ratio of about",
        "0.5 against roughly 0.75-0.8 in Figure S1; see the vignette. Use",
        "1285 to reproduce the paper. Constant within a subject; the model",
        "assumes the protocol schedule for everyone."
      ),
      source_name = "D_FOLFOX (Monolix regressor Dose1)"
    ),
    DOSE_PANITUMUMAB_MGKG = list(
      description = "Panitumumab dose per 14-day treatment cycle; 0 when panitumumab is not given",
      units = "mg/kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "PRIME gave panitumumab 6 mg/kg once every two weeks, so",
        "DOSE_PANITUMUMAB_MGKG = 6 reproduces the paper (Baaz 2023",
        "Methods). The Monolix code in Data S1 uses a 0/1 regressor",
        "(Dose2) and an a2 term per month; Table 1 reports the same",
        "potency per mg/kg, so kdrug_panitumumab * 6 equals the Monolix a2.",
        "Constant within a subject."
      ),
      source_name = "D_Pani (Monolix regressor Dose2)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Exploratory post hoc covariate analysis only (Baaz 2023 Table 3): weight 1e-2 (RSE 29%) on kg and 4e-3 (RSE 37%) on aFOLFOX, in theta_i = theta_m * exp(eta_i + Z_i * w). Not encoded: the analysis does not report the covariate units or centering values, and the uncentered form as printed moves the typical baseline SLD far from the observed baseline in Figure S1. The paper's PFS predictions all use the covariate-free model."
    ),
    ALP = list(
      description = "Alkaline phosphatase",
      units = "U/L (assumed; not stated)",
      type = "continuous",
      notes = "Baaz 2023 Table 3: weight 5e-5 (RSE 27%) on SLD0 (Table S2 automated builder: same). Not encoded; units and centering are not reported (see AGE)."
    ),
    HGB = list(
      description = "Hemoglobin",
      units = "not stated",
      type = "continuous",
      notes = "Baaz 2023 Table 3: weight -9e-3 (RSE 5%) on SLD0 and -5e-3 (RSE 29%) on gamma. Not encoded; units and centering are not reported (see AGE)."
    ),
    LDH = list(
      description = "Lactate dehydrogenase",
      units = "U/L (assumed; not stated)",
      type = "continuous",
      notes = "Baaz 2023 Table 3: weight -2e-4 (RSE 23%) on SLD0. Not encoded; units and centering are not reported (see AGE)."
    ),
    NTARGET = list(
      description = "Number of target lesions at baseline",
      units = "(count)",
      type = "continuous",
      notes = "Baaz 2023 Table 3 'xLB' (lesions at baseline): weight 0.26 (RSE 5%) on SLD0 and 0.17 (RSE 25%) on gamma; the most important covariate for SLD0. Not encoded; centering is not reported (see AGE)."
    ),
    PLT = list(
      description = "Platelet count",
      units = "not stated",
      type = "continuous",
      notes = "Baaz 2023 Table 3: weight -2e-3 (RSE 37%) on kg. Not encoded; units and centering are not reported (see AGE)."
    ),
    WBC = list(
      description = "White blood cell count",
      units = "not stated",
      type = "continuous",
      notes = "Baaz 2023 Table 3: weight -3e-2 (RSE 40%) on aFOLFOX. Not encoded; units and centering are not reported (see AGE)."
    ),
    TUM_NRAS_EXON2_FAIL = list(
      description = "NRAS exon 2 mutation test failed (1 = failure)",
      units = "(binary)",
      type = "binary",
      notes = "Baaz 2023 Table 3 'xNRAS2,F': weight -2 (RSE 36%) on gamma. Not encoded; part of the exploratory covariate analysis (see AGE)."
    )
  )

  compartmentData <- list(
    tumor_size = list(
      analyte = "tumor size (sum of the largest diameters of the target lesions, RECIST 1.1)",
      units = "mm",
      specimen = "not applicable",
      verified = TRUE
    ),
    cumhaz = list(
      analyte = "cumulative hazard of non-target progression (new lesion, unequivocal non-target progression or death)",
      units = "(unitless)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 248L,
    n_studies = 1L,
    age_range = "27-85 years (whole PRIME population)",
    age_median = "61-62 years by arm (whole PRIME population)",
    sex_female_pct = 35,
    race_ethnicity = c(White = 92),
    disease_state = "metastatic colorectal cancer, RAS exon 2 wild-type, first-line treatment",
    dose_range = "FOLFOX4 every 2 weeks (oxaliplatin 85 mg/m^2 day 1; leucovorin 200 mg/m^2 days 1-2; 5-FU 400 mg/m^2 bolus then 600 mg/m^2 infusion days 1-2), alone or with panitumumab 6 mg/kg every 2 weeks",
    regions = "multinational (56-60% Western Europe, Canada, Australia)",
    notes = paste(
      "PRIME phase III study data obtained via Project Data Sphere.",
      "Patients with a RAS exon 2 mutation or only one SLD measurement were",
      "excluded, leaving 127 FOLFOX and 121 FOLFOX + panitumumab patients",
      "(Baaz 2023 Methods). SLD was assessed about every 2 months; the first",
      "radiologist's reading was used. Demographics (Table S1) describe the",
      "whole randomised PRIME arms (325 and 331 patients), not the analysed",
      "subset: 62-67% male, 91-93% White, ECOG 0-1 in 94%, colon primary in",
      "65-66%, liver metastases in 86-88%."
    )
  )

  ini({
    # Baaz 2023 Table 1, 'Estimated parameters for the TGI and TTE models'.
    # Monolix SAEM; parameters with BSV are log-normal, theta_i =
    # theta_m * exp(eta_i), with a diagonal Omega (Methods).

    # TGI parameters (time in months).
    lrbase <- log(123)
    label("Typical baseline sum of largest target diameters SLD0 (mm)") # Table 1: SLD0 = 123 mm (RSE 5%)

    lkge <- log(0.074)
    label("Tumor growth rate kg (1/month)") # Table 1: kg = 0.074 1/month (RSE 7%)

    lkdrug_folfox <- log(2.0e-4)
    label("FOLFOX potency aFOLFOX, kill rate per mg/m^2 of FOLFOX per cycle (m^2/(mg*month))") # Table 1: aFOLFOX = 2.0e-4 m2/(mg month) (RSE 4%)

    lkdrug_panitumumab <- log(7.3e-3)
    label("Panitumumab potency aPani, growth-rate reduction per mg/kg per cycle (kg/(mg*month))") # Table 1: aPani = 7.3e-3 kg/(mg month) (RSE 15%)

    llambda <- log(0.18)
    label("Resistance rate gamma shared by FOLFOX and panitumumab (1/month)") # Table 1: gamma = 0.18 1/month (RSE 10%)

    # NTP time-to-event parameters.
    e_dsld_haz <- 0.0041
    label("Slope linking the NTP hazard to d(tumor_size)/dt, alpha (1/mm)") # Table 1: alpha = 0.0041 1/mm (RSE 22%)

    lbase_haz <- log(0.035)
    label("Baseline NTP hazard beta (1/month)") # Table 1: beta = 0.035 1/month (RSE 19%)

    # Table 1 'BSV' column. Monolix reports omega as the standard deviation of
    # eta, so each value is squared to the variance nlmixr2 expects. Methods:
    # 'no correlation between random effects was assumed'; the random effect
    # on aPani was dropped (high correlation with aFOLFOX).
    etalrbase ~ 0.5776 # Table 1: BSV SLD0 = 0.76 (RSE 5%) -> 0.76^2
    etalkge ~ 0.2809 # Table 1: BSV kg = 0.53 (RSE 11%) -> 0.53^2
    etalkdrug_folfox ~ 0.1681 # Table 1: BSV aFOLFOX = 0.41 (RSE 7.66%) -> 0.41^2
    etallambda ~ 1.1236 # Table 1: BSV gamma = 1.06 (RSE 7%) -> 1.06^2

    # Residual error on SLD: 'an additive and a proportional term' (Methods).
    addSd <- 2.27
    label("Additive residual error on SLD (mm)") # Table 1: s1 = 2.27 mm (RSE 12%), footnote 'Standard deviation of additive error'
    propSd <- 0.15
    label("Proportional residual error on SLD (fraction)") # Table 1: s2 = 0.15 (RSE 5%), footnote 'Standard deviation of proportional error'
  })
  model({
    # Individual parameters.
    rbase <- exp(lrbase + etalrbase)
    kge <- exp(lkge + etalkge)
    kdrug_folfox <- exp(lkdrug_folfox + etalkdrug_folfox)
    kdrug_panitumumab <- exp(lkdrug_panitumumab)
    lambda <- exp(llambda + etallambda)
    base_haz <- exp(lbase_haz)

    # Dose-driven kill rates, k_s,i = a_i * D_i * exp(-gamma * t) (Methods).
    # Data S1 Monolix code: atot = (a1 + a2) * exp(-gamma * t) with a1, a2
    # the per-month rates; Table 1 reports them per unit of per-cycle dose.
    kkill <- (kdrug_folfox * DOSE_FOLFOX_MGM2 +
      kdrug_panitumumab * DOSE_PANITUMUMAB_MGKG) * exp(-lambda * t)

    # TGI: dSLD/dt = kg*SLD - ks,Pani*SLD - ks,FOLFOX*SLD. Data S1 caps the
    # growth at SLD > 1e6 mm to avoid overflow.
    dsld <- kge * tumor_size - kkill * tumor_size
    if (tumor_size > 1e6) {
      dsld <- 0
    }
    d/dt(tumor_size) <- dsld
    tumor_size(0) <- rbase

    # NTP hazard linked to the SLD rate of change: h1 = max(alpha*dSLD/dt +
    # beta, 0) (Methods; Data S1). The survivor function of NTP is sur.
    haz_ntp <- max(e_dsld_haz * dsld + base_haz, 0)
    d/dt(cumhaz) <- haz_ntp
    sur <- exp(-cumhaz)

    # Observation: Monolix combined error with additive s1 and proportional
    # s2; combined1 (SD = s1 + s2*f) is the Monolix default.
    tumor_size ~ add(addSd) + prop(propSd) + combined1()
  })
}
