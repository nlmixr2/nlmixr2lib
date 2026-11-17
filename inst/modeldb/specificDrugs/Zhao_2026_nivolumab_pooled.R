Zhao_2026_nivolumab_pooled <- function() {
  description <- "Two-compartment population PK model for subcutaneous and intravenous nivolumab (anti-PD-1 IgG4) in advanced solid tumors, pooled from historical intravenous studies and the subcutaneous CheckMate 8KX study, with first-order subcutaneous absorption, logit-scale bioavailability and time-varying clearance (sigmoid Emax); the Stage 1 pre-specified model that supplied the $PRIOR information for the CheckMate 67T analysis (Zhao 2026)"
  reference <- "Zhao Y, Hu Z, Gibson E, Zhu L, Vezina H, Kondic AG. Model-based clinical pharmacology profiling of subcutaneous nivolumab in renal cell carcinoma: bridging to alternative dosing regimens and indications. Clin Pharmacol Ther. 2026. doi:10.1002/cpt.70481"
  vignette <- "Zhao_2026_nivolumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  # Same structure as the Stage 2 model: two compartments, zero-order
  # intravenous infusion into the sampled serum pool and first-order
  # absorption from the subcutaneous injection site (Data S1).
  compartmentData <- list(
    depot = list(analyte = "nivolumab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "nivolumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "nivolumab", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power scaling on CL and Vc with reference weight 80 kg (Table S4 Note).",
      source_name = "WTB"
    ),
    CRCL = list(
      description = "Baseline estimated glomerular filtration rate",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Power scaling on CL with reference 90 mL/min/1.73 m^2 (Table S4 Note). Source column GFRB / EGFR; stored under the canonical CRCL (BSA-normalized renal function).",
      source_name = "GFRB"
    ),
    SEXF = list(
      description = "Biological sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "The reference subject is male (Table S4 Note). Female sex carries exponential effects on CL and Vc and a multiplicative factor on the typical subcutaneous bioavailability, the same covariate structure as the Stage 2 model (Data S1).",
      source_name = "SEX"
    ),
    ECOG_GE1 = list(
      description = "Baseline Eastern Cooperative Oncology Group (ECOG) performance-status indicator (1 if ECOG >= 1, else 0)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ECOG performance status = 0)",
      notes = "Exponential effect on CL. Table S4 also lists an F1_PS factor held at 1.00 (not estimated), carried here as a fixed() parameter so the effect is inert. Renamed from the source column PS to the canonical ECOG_GE1.",
      source_name = "PS"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 3540L,
    n_studies = 20L,
    disease_state = "Advanced solid tumors (NSCLC, melanoma, RCC, urothelial carcinoma, SCCHN, gastric cancer, colorectal cancer, hepatocellular carcinoma and others) and classical Hodgkin lymphoma",
    dose_range = "Intravenous nivolumab 0.1 - 20 mg/kg (and 240 mg flat) Q2W or Q3W; subcutaneous nivolumab with rHuPH20 720 or 960 mg single dose, 1200 mg Q4W or 600 mg Q2W (CheckMate 8KX)",
    regions = "Global (phase I / II / III studies)",
    n_subjects_sc = 82L,
    n_subjects_iv = 3458L,
    notes = "Stage 1 dataset: CheckMate 8KX (CA2098KX, 82 patients, 683 samples; subcutaneous nivolumab without rHuPH20 excluded) plus historical intravenous monotherapy studies. The main text says 17 historical intravenous studies; Table S2 lists 19 intravenous studies (CA209001, 003, 005, 009, 010, 017, 025, 026, 032, 037, 039, 051, 057, 063, 066, 067, 141, 205, 275), and n_subjects is the Table S2 sum over those 19 studies plus CheckMate 8KX. Baseline demographics of the Stage 1 dataset are not reported; the reference values (80 kg, eGFR 90 mL/min/1.73 m^2, 65-year-old male, PS 0) are described as the approximate dataset medians (Table S4 Note)."
  )

  ini({
    # Structural parameters, Table S4 ('Parameter estimates of the final
    # population pharmacokinetic model from stage 1'). Reference subject:
    # 80 kg, eGFR 90 mL/min/1.73 m^2, 65-year-old male, PS 0 (Table S4
    # Note). The source time unit is hours; rates are converted to days.
    lcl <- log(0.0108 * 24); label("Baseline clearance CL0_REF at the reference covariates (L/day)")     # Table S4: CL = 0.0108 L/hr; x 24 = 0.2592 L/day
    lvc <- log(4.29);        label("Central volume of distribution VC_REF (L)")                          # Table S4: VC = 4.29 L
    lq  <- log(0.0312 * 24); label("Intercompartmental clearance Q_REF (L/day)")                         # Table S4: Q = 0.0312 L/hr; x 24 = 0.7488 L/day
    lvp <- log(2.73);        label("Peripheral volume of distribution VP_REF (L)")                       # Table S4: VP = 2.73 L

    # Subcutaneous absorption.
    lka <- log(0.0129 * 24);                  label("First-order subcutaneous absorption rate constant KA_REF (1/day)")     # Table S4: KA = 0.0129 /hr; x 24 = 0.3096 /day
    logitfdepot <- log(0.815 / (1 - 0.815));  label("Logit of typical subcutaneous bioavailability F_REF (fraction)")       # Table S4: F1 = 0.815; logit(0.815) = 1.4828

    # Covariate effects on CL0: power on WT and eGFR; exponential on female
    # sex and PS >= 1.
    e_wt_cl       <-  0.607; label("Power exponent of WT on CL (unitless)")                                  # Table S4: CL_WTB = 0.607
    e_crcl_cl     <-  0.166; label("Power exponent of CRCL (eGFR) on CL (unitless)")                         # Table S4: CL_GFR = 0.166
    e_sex_cl      <- -0.169; label("Exponential coefficient of female sex on CL (unitless)")                 # Table S4: CL_SEX = -0.169
    e_ecog_ge1_cl <-  0.169; label("Exponential coefficient of ECOG_GE1 on CL (unitless)")                   # Table S4: CL_PS = 0.169

    # Covariate effects on Vc: power on WT; exponential on female sex.
    e_wt_vc  <-  0.621; label("Power exponent of WT on Vc (unitless)")                                       # Table S4: VC_WTB = 0.621
    e_sex_vc <- -0.136; label("Exponential coefficient of female sex on Vc (unitless)")                      # Table S4: VC_SEX = -0.136

    # Covariate effects on the typical subcutaneous bioavailability,
    # multiplicative on the fraction scale before the logit transform.
    e_sex_fdepot      <- 0.895;    label("Multiplicative factor of female sex on subcutaneous bioavailability (unitless)")  # Table S4: F1_SEX = 0.895
    e_ecog_ge1_fdepot <- fixed(1); label("Multiplicative factor of ECOG_GE1 on subcutaneous bioavailability (unitless)")    # Table S4: F1_PS = 1.00, held constant by the authors (superscript f)

    # Time-varying clearance, sigmoid Emax in time since first dose. CL
    # approaches CL0 x exp(-0.263) = 0.769 x baseline.
    cl_time_max   <- -0.263;       label("Maximal change in log CL, EMAX_REF (unitless)")                    # Table S4: CL_EMAX = -0.263
    lcl_t50       <- log(1.53e3 / 24); label("Log of the time to 50%% of the maximal change in CL, T50 (log days)")         # Table S4: CL_T50 = 1.53E+03 hr; / 24 = 63.75 days
    lcl_time_hill <- log(3.84);       label("Log of the Hill / sigmoidicity exponent of time on CL (log unitless)")  # Table S4: CL_HILL = 3.84

    # Inter-individual variability (Table S4 'Random effects'): variance (SD)
    # on the diagonal, covariance (correlation) off it; each parenthetical
    # reproduces (0.0402 / sqrt(0.122 x 0.0853) = 0.394; 0.180 /
    # sqrt(0.132 x 0.807) = 0.551).
    etalcl + etalvc ~ c(0.122,
                        0.0402, 0.0853)                                                                      # Table S4: ZCL 0.122 (0.349), ZCL:ZVC 0.0402 (0.395), ZVC 0.0853 (0.292)
    etalvp ~ 0.212                                                                                           # Table S4: ZVP 0.212 (0.461)
    etacl_time_max ~ 0.0537                                                                                  # Table S4: ZEMAX 0.0537 (0.232); additive on EMAX
    etalka + etalogitfdepot ~ c(0.132,
                                0.180, 0.807)                                                                # Table S4: ZKA 0.132 (0.364), ZKA:F1 0.180 (0.551), F1 0.807 (0.898)

    # Residual error: proportional standard deviation (Table S4 theta10).
    propSd <- 0.188; label("Proportional residual error (fraction)")                                         # Table S4: Proportional = 0.188 (theta10)
  })
  model({
    cl_base <- exp(lcl + etalcl) *
      (WT / 80)^e_wt_cl *
      (CRCL / 90)^e_crcl_cl *
      exp(e_sex_cl * SEXF) *
      exp(e_ecog_ge1_cl * ECOG_GE1)

    vc <- exp(lvc + etalvc) *
      (WT / 80)^e_wt_vc *
      exp(e_sex_vc * SEXF)

    # Time-varying CL with additive IIV on EMAX; t is time since first dose.
    cl_t50 <- exp(lcl_t50)
    cl_time_hill <- exp(lcl_time_hill)
    cl_time_max_i <- cl_time_max + etacl_time_max
    cl <- cl_base * exp(cl_time_max_i * t^cl_time_hill / (cl_t50^cl_time_hill + t^cl_time_hill))

    vp <- exp(lvp + etalvp)
    q <- exp(lq)
    ka <- exp(lka + etalka)

    fdepot_tv <- expit(logitfdepot) *
      e_sex_fdepot^SEXF *
      e_ecog_ge1_fdepot^ECOG_GE1
    fdepot <- expit(logit(fdepot_tv) + etalogitfdepot)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Subcutaneous doses go to depot; intravenous doses are infused into
    # central with F = 1.
    f(depot) <- fdepot

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
