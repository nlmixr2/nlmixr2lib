Takada_2025_vancomycin <- function() {
  description <- "Two-compartment IV population PK model for vancomycin in Japanese patients of advanced age (aged 75 years and older, body mass index below 25 kg/m^2) receiving therapeutic drug monitoring (Takada 2025). Clearance scales as a power function of Cockcroft-Gault creatinine clearance (exponent 0.63, reference 3.09 L/h = 51.5 mL/min) and of serum albumin (exponent 0.22, reference 2.3 g/dL); the albumin term is the novelty of this analysis, added because creatinine-based renal-function estimates underestimate clearance in low-muscle-mass patients of advanced age. Intercompartmental clearance and both volumes are covariate-free. Between-subject variability is on clearance only; the residual-error magnitudes were not reported by the source and are encoded as zero."
  reference <- "Takada K, Samura M, Igarashi Y, Suzuki A, Ishigo T, Fujii S, Ibe Y, Yoshida H, Tanaka H, Ebihara F, Maruyama T, Hamada Y, Komatsu T, Tomizawa A, Takuma A, Chiba H, Yagi Y, Nishi Y, Enoki Y, Taguchi K, Tanikawa K, Kunishima H, Matsumoto K. Development and validation of a population pharmacokinetic model of vancomycin for patients of advanced age. J Pharm Health Care Sci. 2025;11:22. doi:10.1186/s40780-025-00423-8"
  vignette <- "Takada_2025_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Vancomycin is given as an IV infusion, so the dose
  # enters `central` directly and there is no depot state. Takada 2025
  # Methods says only "blood concentrations of VCM", assayed by latex
  # immunoturbidimetry (Nanopia TDM Vancomycin, Sekisui Medical) in the
  # modeling cohort; the paper never states serum versus plasma, so
  # `specimen` is left unverified for both states.
  compartmentData <- list(
    central     = list(analyte = "vancomycin", units = "mg", specimen = "serum",  verified = FALSE),
    peripheral1 = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalized), from serum creatinine measured by the enzymatic method",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column CLcr. Takada 2025 Additional File 2: Table 2 (file MOESM4), Equations 1-2, give the estimating equation: CLcr (male, mL/min) = ([140 - age] * body weight) / (serum creatinine * 72); CLcr (female) = CLcr (male) * 0.85. Raw Cockcroft-Gault in mL/min, NOT BSA-normalized, following the Delattre_2010_amikacin.R precedent for storing raw Cockcroft-Gault under the canonical CRCL column. The source model equation is written in L/h with reference 3.09 L/h; model() converts CRCL from mL/min to L/h via * 0.06 (3.09 L/h == 51.5 mL/min, the Table 1 modeling-cohort median). Modeling-cohort median 51.5 mL/min (range 9.7-121.2) per Table 1; the Results narrative quotes 3.06 L/h [51.0 (9.7-121.0) mL/min] for the same quantity, a rounding inconsistency with Table 1. The paper's Limitations note that few patients sat near either end of that range.",
      source_name        = "CLcr"
    ),
    ALB = list(
      description        = "Serum albumin, measured by the modified bromocresol purple method",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column Alb, reported by Takada 2025 in US-convention g/dL (Table 1 modeling-cohort median 2.3 g/dL, range 1.2-4.2; 90.6% of the cohort below 3.0 g/dL). Stored under the canonical CRCL-style SI convention g/L, so model() applies the inline conversion ALB * 0.1 before entering the (Alb/2.3)^0.22 power term, which was calibrated in g/dL. Reference 2.3 g/dL == 23 g/L. All eight validation-cohort hospitals also used the modified bromocresol purple method. The paper's rationale (Discussion) is that albumin tracks muscle mass in patients of advanced age, so it corrects the creatinine-based clearance estimate when muscle mass is low.",
      source_name        = "Alb"
    )
  )

  # Screened in the Takada 2025 covariate search (Methods, "Development of a
  # population pharmacokinetic model") but NOT retained in the final model,
  # so they are documentation only and are not referenced in model(). The
  # stepwise search in Table 2 retained only CLcr and Alb, both on CL; no
  # covariate was retained on Vc.
  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as a candidate covariate on Vc (Methods: 'Age, sex, actual body weight, BMI, and Alb levels were used to explore the covariates for Vc'). Not retained; Table 2 shows no Vc covariate entering the stepwise search. Age also enters the model indirectly as an input to the Cockcroft-Gault CLcr equation.",
      source_name        = "Age"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "male",
      notes              = "Screened as a candidate covariate on Vc (Methods). Not retained. Sex also enters the model indirectly through the 0.85 female multiplier in the Cockcroft-Gault CLcr equation (Additional File 2: Table 2, file MOESM4, Equation 2). Modeling cohort 42.1% female (Table 1).",
      source_name        = "Sex"
    ),
    WT = list(
      description        = "Actual body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as a candidate covariate on Vc (Methods). Not retained: neither Vc nor Vp carries any covariate in the Table 3 final model, and the paper reports Vc and Vp in absolute litres rather than per-kilogram. Body weight enters indirectly through the Cockcroft-Gault CLcr equation. Modeling-cohort median 47 kg (range 26-70) per Table 1.",
      source_name        = "Body weight"
    ),
    BMI = list(
      description        = "Body mass index",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as a candidate covariate on both CL and Vc (Methods). Not retained; Table 2's stepwise search shows only renal-function and albumin terms entering CL. BMI below 25 kg/m^2 was an inclusion criterion, so the covariate range is truncated by design (modeling-cohort median 18.6 kg/m^2, range 11.0-24.8; 48.4% below 18.5).",
      source_name        = "BMI"
    ),
    CREAT = list(
      description        = "Serum creatinine, measured by the enzymatic method",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Not screened directly as a covariate, but it is the input to both renal-function estimates and it defines the validation subgroup the paper's conclusion rests on (patients with serum creatinine below 0.6 mg/dL, where the developed model beat every previous model). Modeling-cohort median 0.64 mg/dL (range 0.22-3.00); 42.1% of the modeling cohort had serum creatinine below 0.60 mg/dL (Table 1).",
      source_name        = "SCr"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 159L,
    n_studies        = 1L,
    n_sites          = 1L,
    n_concentrations = 417L,
    age_range        = "75-99 years (eligibility: 75 years and older)",
    age_median       = "84 years (range 75-99)",
    weight_range     = "26-70 kg",
    weight_median    = "47 kg (range 26-70)",
    height_median    = "158 cm (range 135-180)",
    bmi_median       = "18.6 kg/m^2 (range 11.0-24.8); eligibility required BMI below 25 kg/m^2",
    sex_female_pct   = 42.1,
    disease_state    = "Inpatients of advanced age receiving intravenous vancomycin under routine therapeutic drug monitoring for any indication. Inclusion: age 75 years and older and BMI below 25 kg/m^2. Exclusion: patients undergoing dialysis, patients who received vancomycin for fewer than 3 days, and patients with no vancomycin blood samples. Age 85 years and older accounted for 49.7% of the modeling cohort and BMI below 18.5 kg/m^2 for 48.4%.",
    dose_range       = "Daily maintenance dose median 1500 mg/day (range 250-3000); 67.9% received 1000-2000 mg/day, 27.7% below 1000 mg/day, 4.4% above 2000 mg/day. Doses chosen by the treating physician, not protocolized.",
    regions          = "Japan (Yokohama General Hospital, Kanagawa)",
    renal_function   = "Cockcroft-Gault creatinine clearance median 51.5 mL/min (range 9.7-121.2), with 47.2% at or below 50 mL/min; creatinine-based eGFR median 65.5 mL/min (range 12.4-159.1), with 28.3% at or below 50 mL/min; serum creatinine median 0.64 mg/dL (range 0.22-3.00), with 42.1% below 0.60 mg/dL",
    nutritional_status = "Serum albumin median 2.3 g/dL (range 1.2-4.2): below 2.0 g/dL in 19.5%, 2.0-2.9 in 71.1%, 3.0-3.5 in 5.7%, 3.6 and above in 3.8%",
    notes            = "Retrospective single-centre analysis of routine TDM data collected at Yokohama General Hospital between August 2016 and September 2024 (Takada 2025 Table 1). 417 concentrations from 159 patients: 65 peaks and 352 troughs, with no samples below the limit of quantitation (2.5 ug/mL for the modeling-cohort assay). Fit in Phoenix NLME with FOCE-ELS (equivalent to NONMEM FOCE with interaction); covariates selected by forward inclusion and backward elimination with a 3.85-point drop in the objective function taken as p < 0.05; evaluated by visual predictive check and a 1000-run bootstrap (Table 3). A separate 133-patient multicentre cohort drawn from eight hospitals between September 2020 and December 2023 was used only for external validation of predictive performance and did not contribute to the parameter estimates reproduced here."
  )

  ini({
    # Structural parameters, Takada 2025 Table 3 (final-model column). The
    # reference patient has a Cockcroft-Gault creatinine clearance of
    # 3.09 L/h (51.5 mL/min) and a serum albumin of 2.3 g/dL, both the
    # Table 1 modeling-cohort medians.
    lcl <- log(1.96);  label("Clearance at CLcr=3.09 L/h and Alb=2.3 g/dL (CL, L/h)")            # Takada 2025 Table 3 theta1: 1.96 L/h (SE 0.06, CV 3.11%, 95% CI 1.84-2.08)
    lq  <- log(4.86);  label("Intercompartmental clearance (Q, L/h)")                             # Takada 2025 Table 3 theta4: 4.86 L/h (SE 0.91, CV 18.63%, 95% CI 3.08-6.64); see vignette Errata for the conflicting 3.24 in the Results narrative
    lvc <- log(31.78); label("Central volume of distribution (Vc, L)")                            # Takada 2025 Table 3 theta5: 31.78 L (SE 3.86, CV 12.16%, 95% CI 24.19-39.38)
    lvp <- log(53.64); label("Peripheral volume of distribution (Vp, L)")                         # Takada 2025 Table 3 theta6: 53.64 L (SE 4.22, CV 7.86%, 95% CI 45.36-61.93)

    # Covariate effects on clearance, both estimated (SEs and 95% CIs
    # reported in Table 3), so neither is wrapped in fixed().
    e_crcl_cl <- 0.63; label("Power exponent on (CLcr/3.09 L/h) for CL (unitless)")               # Takada 2025 Table 3 theta2: 0.63 (SE 0.06, CV 9.72%, 95% CI 0.51-0.76)
    e_alb_cl  <- 0.22; label("Power exponent on (Alb/2.3 g/dL) for CL (unitless)")                # Takada 2025 Table 3 theta3: 0.22 (SE 0.09, CV 41.75%, 95% CI 0.03-0.40)

    # Between-subject variability on clearance only. The Table 3 footnote
    # defines eta as "normally distributed with mean 0 and variance omega^2,
    # etaCL = 0.11", i.e. 0.11 is the VARIANCE of etaCL (34.1% CV on the
    # exponential scale). No random effect is reported for Q, Vc, or Vp; the
    # Results narrative states that Vc "was a fixed value". The 7-parameter
    # base model of Table 2 is consistent with exactly one omega: four
    # structural thetas + one omega + two residual-error parameters.
    # NOTE: the paper's own dosing nomograms substitute the constant 0.11
    # into exp(etaCL) as if it were a fixed offset -- see vignette Errata.
    etalcl ~ 0.11  # Takada 2025 Table 3 footnote: etaCL = 0.11, stated as the variance omega^2

    # Residual variability. Takada 2025 Additional File 1: Table 1 (file MOESM3) lists the
    # three candidate Phoenix residual models that were evaluated, of which
    # the two-parameter "additive plus multiplicative" form
    #   Cobs = C + eps * sqrt[1 + C^2 * (Cmultstdev/sigma)^2]
    # is algebraically identical to nlmixr2's add(sigma) + prop(Cmultstdev),
    # and is the only one of the three consistent with the 7-parameter base
    # model of Table 2. The paper reports neither which model was selected
    # nor any residual-error VALUE anywhere in the article or its ten
    # supplementary files, so both magnitudes are encoded as zero rather than
    # invented. Simulations from this model are therefore residual-error-free.
    propSd <- fixed(0); label("Proportional residual SD (fraction; 0 -- not reported in the source)")  # Takada 2025: no residual-error estimate published
    addSd  <- fixed(0); label("Additive residual SD (ug/mL; 0 -- not reported in the source)")         # Takada 2025: no residual-error estimate published
  })
  model({
    # Unit alignment. The source model equation is written with creatinine
    # clearance in L/h and albumin in US-convention g/dL, whereas the
    # canonical covariate columns carry CRCL in mL/min and ALB in SI g/L.
    crcl_Lh <- CRCL * 0.06  # mL/min -> L/h (51.5 mL/min == 3.09 L/h, the reference)
    alb_gdL <- ALB * 0.1    # SI g/L -> US-convention g/dL (23 g/L == 2.3 g/dL, the reference)

    # Individual PK parameters. Takada 2025 Table 3 / Results:
    #   CL (L/h) = 1.96 * (CLcr/3.09)^0.63 * (Alb/2.3)^0.22 * exp(etaCL)
    #   Q  (L/h) = 4.86
    #   Vc (L)   = 31.78
    #   Vp (L)   = 53.64
    cl <- exp(lcl + etalcl) * (crcl_Lh / 3.09)^e_crcl_cl * (alb_gdL / 2.3)^e_alb_cl
    q  <- exp(lq)
    vc <- exp(lvc)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # Dose in mg, volumes in L, so central/vc is mg/L == ug/mL, the units the
    # source reports vancomycin concentrations and AUCss in.
    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
