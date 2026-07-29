Thomson_1989_lisinopril <- function() {
  description <- "One-compartment population PK model for oral lisinopril (an ACE inhibitor) at steady state in elderly and renal-disease hypertensive adults (Thomson 1989). First-order absorption with apparent clearance CL/F driven by body weight, serum creatinine, age, and a binary compensated-cardiac-failure indicator; apparent volume V/F and absorption rate ka are population means without retained covariate effects."
  reference   <- "Thomson AH, Kelly JG, Whiting B. Lisinopril population pharmacokinetics in elderly and renal disease patients with hypertension. Br J Clin Pharmacol 1989;27(1):57-65. doi:10.1111/j.1365-2125.1989.tb05335.x."
  vignette    <- "Thomson_1989_lisinopril"
  units       <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight at the steady-state PK profile (kg)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear multiplicative effect on CL/F: CL/F = 0.251 * WT * (other factors). Reference cohort mean 72 kg (Thomson 1989 'Data set'). The CL/F coefficient theta1 = 0.251 L/(h*kg) -- not a separate fixed allometric exponent. Thomson 1989 Discussion warns that extrapolation above 90 kg is unwise because only four obese patients were included.",
      source_name        = "wt"
    ),
    AGE = list(
      description        = "Age at the steady-state PK profile (years)",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL/F with reference 65 years: (AGE/65)^(-0.451) per Thomson 1989 Eq. and Table 3. Cohort mean 65 years; the renal-disease trial enrolled two patients under 40 years, so Thomson 1989 Discussion warns against extrapolation below 40 years.",
      source_name        = "Age"
    ),
    CREAT = list(
      description        = "Serum creatinine concentration at the steady-state PK profile (umol/L)",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL/F with reference 70 umol/L: (CREAT/70)^(-0.887) per Thomson 1989 Eq. and Table 3. Captures the dominant renal-function effect on lisinopril clearance; Thomson 1989 reports this non-linear creatinine-only encoding outperformed both the linear-on-Cr and the Jelliffe-nomogram creatinine-clearance covariate (Table 2 models 2-7).",
      source_name        = "Cr"
    ),
    DIS_CHF = list(
      description        = "Compensated congestive heart failure indicator at study entry",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no diagnosed cardiac failure)",
      notes              = "Multiplicative covariate on CL/F: CL/F is multiplied by 0.645 (~35% lower) when DIS_CHF = 1 per Thomson 1989 Table 2 model 11 and Eq. for cardiac failure. Encoded as theta6^DIS_CHF in the model() body so that DIS_CHF = 0 gives factor 1 and DIS_CHF = 1 gives factor 0.645. Time-fixed per subject; 13 of 60 patients in the analysis cohort had compensated cardiac failure. Thomson 1989 Discussion notes that all CHF patients were on background cardiac glycosides, and a likelihood comparison favoured the disease indicator over the concomitant-medication indicator (Table 2 model 11 LLD = 24 vs model 12 LLD = 15), so the effect is attributed to the disease, not the drug.",
      source_name        = "chf"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 60L,
    n_studies      = 2L,
    n_observations = 381L,
    age_range      = "Adults (cohort mean 65 years; two renal-disease patients under 40 years; elderly trial inclusion was 65+ years)",
    age_median     = "mean 65 years",
    weight_range   = "Adults (cohort mean 72 kg; only four obese patients > 90 kg)",
    weight_median  = "mean 72 kg",
    sex_female_pct = NA_real_,
    disease_state  = "Hypertensive adults from two multicentre trials: Trial I (n=40, elderly hypertension; age >= 65 years; mild-to-moderate systolic/diastolic or isolated systolic hypertension) and Trial II (n=20, renal-disease hypertension stratified by creatinine clearance 30-60, <30, or on haemodialysis). 13 patients had compensated cardiac failure on background cardiac glycosides; one patient was on haemodialysis.",
    dose_range     = "Oral lisinopril 2.5-40 mg daily (median 10 mg daily) at the steady-state profile.",
    regions        = "United Kingdom and Ireland (multicentre)",
    notes          = "Steady-state concentration-time profiles (samples at 0, 1, 2, 4, 6, 8, 12 h after the morning dose) collected after at least 2 weeks at a constant dose. Lisinopril measured by radioimmunoassay (Hichens et al. 1981). 79 of the initial 140 enrolled patients were excluded for missing dosing/sampling, missing biochemistry, or unreliable compliance; one outlier with tenfold-elevated concentrations was also excluded after preliminary investigation. Comorbidities and concomitant medications detailed in Thomson 1989 Table 1; 22 of 60 patients (37%) received no other drugs during the lisinopril study period. Demographic distributions (age, weight, serum creatinine) shown in Thomson 1989 Figure 1; dose and peak-concentration distributions in Figure 2."
  )

  ini({
    # Structural parameters -- Thomson 1989 Table 3 final-model estimates (column theta_i).
    # Reference subject for the covariate model: WT = 70 kg, CREAT = 70 umol/L, AGE = 65 years,
    # DIS_CHF = 0. The structural CL/F coefficient theta1 = 0.251 L/(h*kg) gives the per-kg
    # apparent clearance; CL/F at the reference subject is then 0.251 * 70 = 17.6 L/h.
    lcl <- log(0.251); label("Apparent clearance coefficient theta1 (L/(h*kg))")   # Thomson 1989 Table 3: theta1 = 0.251 (SE 0.029)
    lvc <- log(36.7);  label("Apparent volume of distribution V/F (L)")            # Thomson 1989 Table 3: theta2 = 36.7 (SE 3.9)
    lka <- log(0.104); label("First-order absorption rate constant ka (1/h)")      # Thomson 1989 Table 3: theta3 = 0.104 (SE 0.006)

    # Covariate-effect parameters on CL/F. Power-form exponents on (CREAT/70) and (AGE/65)
    # are estimated; the cardiac-failure factor is a multiplicative scalar raised to the
    # binary DIS_CHF indicator (encoded as theta6^DIS_CHF so DIS_CHF = 0 leaves CL/F
    # unchanged and DIS_CHF = 1 multiplies CL/F by 0.645).
    e_creat_cl <- -0.887; label("Power exponent of (CREAT/70) on CL/F (unitless)")   # Thomson 1989 Table 3: theta4 = -0.887 (SE 0.108)
    e_age_cl   <- -0.451; label("Power exponent of (AGE/65) on CL/F (unitless)")     # Thomson 1989 Table 3: theta5 = -0.451 (SE 0.172)
    e_chf_cl   <-  0.645; label("Multiplicative factor on CL/F when DIS_CHF = 1 (unitless)")   # Thomson 1989 Table 3: theta6 = 0.645 (SE 0.112)

    # Inter-individual variability. Thomson 1989 reports "log-additive" IIV, i.e. P_i =
    # P_typ * exp(eta_i), eta_i ~ N(0, omega^2). The values below are the omega^2 variances
    # taken directly from Table 3 (ii). The reported CV%s use the small-variance
    # approximation CV ~ sqrt(omega^2): sqrt(0.266) = 51.6% (rounded to 52% in the paper),
    # sqrt(1.61) = 127%, sqrt(0.538) = 73%.
    etalcl ~ 0.266   # Thomson 1989 Table 3 (ii): omega^2_CL = 0.266 (SE 0.059); CV ~ 52%
    etalvc ~ 1.61    # Thomson 1989 Table 3 (ii): omega^2_V  = 1.61 (SE 0.52);  CV ~ 127%
    etalka ~ 0.538   # Thomson 1989 Table 3 (ii): omega^2_ka = 0.538 (SE 0.198); CV ~ 73%

    # Residual error. Thomson 1989 used a NONMEM log-additive error model:
    # log(Cobs) = log(Cpred) + epsilon, epsilon ~ N(0, sigma^2 = 0.0772). For small sigma
    # this maps directly to a linear-space proportional residual with propSd = sqrt(sigma^2)
    # (CV ~ 28%, matching the paper's reported residual CV of 28%). The lognormal form
    # `Cc ~ lnorm(expSd)` with expSd = sqrt(0.0772) is an alternative exact mapping; the
    # proportional form is preferred here for consistency with the rest of the library and
    # because the variance is small enough that the two forms produce indistinguishable
    # simulations.
    propSd <- sqrt(0.0772); label("Proportional residual SD (fraction)")   # Thomson 1989 Table 3 (ii): sigma^2 = 0.0772 (SE 0.0184); CV ~ 28%
  })

  model({
    # Individual PK parameters with multiplicative log-normal IIV. The apparent-PK form
    # uses the paper's exact final equation:
    #   CL/F = theta1 * WT * (CREAT/70)^theta4 * (AGE/65)^theta5 * theta6^DIS_CHF
    # exp(lcl) = theta1 = 0.251 L/(h*kg); multiplying by WT gives the per-subject apparent
    # clearance in L/h. V/F is the population mean without retained covariates (Table 2
    # models 13-18 found no significant covariate effect on volume).
    cl <- exp(lcl + etalcl) * WT *
            (CREAT / 70)^e_creat_cl *
            (AGE   / 65)^e_age_cl *
            e_chf_cl^DIS_CHF
    vc <- exp(lvc + etalvc)
    ka <- exp(lka + etalka)

    kel <- cl / vc

    # One-compartment open model with first-order absorption (Thomson 1989 Methods: model
    # 19 retained over the zero-order alternative model 20 by both log-likelihood and the
    # observation that standard errors were unobtainable in the zero-order fit). Dose enters
    # the depot compartment in mg; central holds the systemic amount in mg.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Observed plasma lisinopril concentration. central / vc is mg/L = ug/mL; multiplying by
    # 1000 converts to ng/mL (the paper's reporting units; cohort peak range 6.4-343 ng/mL,
    # Thomson 1989 Figure 2b).
    Cc <- (central / vc) * 1000

    Cc ~ prop(propSd)
  })
}
