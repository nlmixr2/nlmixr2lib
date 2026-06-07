Sarashina_2005_epinastine <- function() {
  description <- "Two-compartment population PK model with first-order absorption for oral epinastine in healthy adults and paediatric atopic dermatitis patients (Sarashina 2005), with linear-in-WT CL/F and V1/F plus food-status and formulation covariate effects"
  reference <- "Sarashina A, Tatami S, Yamamura N, Tsuda Y, Igarashi T. Population pharmacokinetics of epinastine, a histamine H1 receptor antagonist, in adults and children. Br J Clin Pharmacol. 2005;59(1):43-53. doi:10.1111/j.1365-2125.2005.02250.x"
  vignette <- "Sarashina_2005_epinastine"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed. Sarashina 2005 enters WT linearly into CL/F and V1/F: CL/F = (theta1 + WT * theta10) * food * form and V1/F = (theta2 + WT * theta11) * food. The y-intercept (theta1, theta2) and per-kg slope (theta10, theta11) are both estimated.",
      source_name        = "WT"
    ),
    FED = list(
      description        = "Fed-vs-fasted dose-record indicator (1 = fed at dosing, 0 = fasted)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted)",
      notes              = "Per-record covariate. Adults: 607 fasted vs 724 non-fasted observations (Table 2). Paediatric: 10 fasting vs 169 non-fasting observations. Fed reduces Cmax and AUC (food/fasted ratios 0.67 and 0.62 respectively per the Discussion), and introduces an absorption lag.",
      source_name        = "FOOD"
    ),
    FORM_SYRUP = list(
      description        = "Dry-syrup vs tablet formulation indicator (1 = dry syrup, 0 = tablet reference)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (tablet)",
      notes              = "Per-subject formulation. Tablet is the reference (theta9 multiplier = 1); dry syrup increases CL/F by a factor of 1.06 (theta9). The covariate-columns register names tablet as one common comparator; document the reference here. Paediatric patients all received dry syrup; healthy adults received either tablet or dry syrup depending on the trial.",
      source_name        = "FORM"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 124,
    n_studies      = 6,
    age_range      = "2-26 years (paediatric 2-15 years; adult 20-26 years)",
    age_median     = "adult 22.3 years (mean); paediatric 10.2 years (mean)",
    weight_range   = "14.1-82 kg (paediatric 14.1-68 kg; adult 50-82 kg)",
    weight_median  = "adult 63.0 kg (mean); paediatric 36.9 kg (mean)",
    sex_female_pct = 19.4,
    race_ethnicity = c(Japanese = 100),
    disease_state  = "62 healthy adult volunteers (all male) plus 62 paediatric atopic dermatitis patients (38 male, 24 female).",
    dose_range     = "10, 20, or 40 mg oral epinastine once daily (paediatric: 10 mg if 14 kg to <24 kg, 20 mg if 24 kg or more).",
    regions        = "Japan",
    notes          = "Six clinical trials pooled (Table 1). 1510 plasma observations: 1331 from adults, 179 from paediatric atopic dermatitis patients. Adult sampling was rich after first dose; paediatric sampling was sparse (3 trough samples per patient at 2-6, 6-10, and 10-14 weeks of daily dosing). See Table 2 for full demographics."
  )

  ini({
    # Structural parameters - final estimates from Sarashina 2005 Table 4.
    # The paper reports CL/F and V1/F as LINEAR functions of body weight:
    #   CL/F = (theta1 + WT * theta10) * FOODCL_F * FORMCL_F
    #   V1/F = (theta2 + WT * theta11) * FOODV1_F
    # Q, V2/F, and Ka have no retained covariate effects (Table 3 / Table 4).
    # Reparameterised as exp(lcl) * (1 + e_wt_cl * WT) using the algebraic
    # identity theta1 + WT * theta10 = theta1 * (1 + (theta10/theta1) * WT)
    # so the canonical 'lcl' is log(theta1) (CL/F intercept at WT = 0) and
    # 'e_wt_cl' = theta10 / theta1 is the fractional WT slope (per kg).
    # 'lcl' and 'lvc' therefore back-transform to the paper's intercept
    # values 19.1 L/h and 174 L, not to a typical-subject CL/F.

    lcl <- log(19.1) ; label("CL/F intercept theta1 (L/h, log)")                              # Sarashina 2005 Table 4: theta1 = 19.1 L/h
    lvc <- log(174)  ; label("V1/F intercept theta2 (L, log)")                                # Sarashina 2005 Table 4: theta2 = 174 L
    lq  <- log(34.4) ; label("Inter-compartmental clearance Q/F theta3 (L/h, log)")           # Sarashina 2005 Table 4: theta3 = 34.4 L/h
    lvp <- log(452)  ; label("Peripheral volume V2/F theta4 (L, log)")                        # Sarashina 2005 Table 4: theta4 = 452 L
    lka <- log(1.18) ; label("First-order absorption rate constant Ka theta5 (1/h, log)")     # Sarashina 2005 Table 4: theta5 = 1.18 1/h

    # Covariate effects (Sarashina 2005 Table 4 and Discussion)
    e_wt_cl         <- 0.0421 ; label("Fractional per-kg WT slope on CL/F (= theta10/theta1, per kg)")   # Sarashina 2005 Table 4: theta10 = 0.805 L/h/kg; 0.805 / 19.1 = 0.04215
    e_wt_vc         <- 0.0227 ; label("Fractional per-kg WT slope on V1/F (= theta11/theta2, per kg)")   # Sarashina 2005 Table 4: theta11 = 3.95 L/kg; 3.95 / 174 = 0.02270
    e_fed_cl        <- 1.41   ; label("Fed/fasted ratio on CL/F (theta6, unitless)")                     # Sarashina 2005 Table 4: theta6 = 1.41
    e_fed_vc        <- 1.75   ; label("Fed/fasted ratio on V1/F (theta7, unitless)")                     # Sarashina 2005 Table 4: theta7 = 1.75
    e_fed_tlag      <- 0.234  ; label("Fed-state absorption lag time (additive, h)")                     # Sarashina 2005 Table 4: theta8 = 0.234 h
    e_form_syrup_cl <- 1.06   ; label("Dry-syrup/tablet ratio on CL/F (theta9, unitless)")               # Sarashina 2005 Table 4: theta9 = 1.06

    # Inter-individual variability. Sarashina 2005 reports the NONMEM
    # omega-squared values in Table 5 with the convention CV% ~ sqrt(omega^2);
    # sqrt(0.101) = 0.318 = 31.8% matches Table 4's reported CL/F CV%.
    etalcl ~ 0.101    # Sarashina 2005 Table 5: omega^2 CL/F = 0.101 (~31.8% CV)
    etalvc ~ 0.107    # Sarashina 2005 Table 5: omega^2 V1/F = 0.107 (~32.7% CV)
    etalq  ~ 0.226    # Sarashina 2005 Table 5: omega^2 Q    = 0.226 (~47.5% CV)
    etalvp ~ 1.43     # Sarashina 2005 Table 5: omega^2 V2/F = 1.43  (~119.6% CV)
    etalka ~ 0.323    # Sarashina 2005 Table 5: omega^2 Ka   = 0.323 (~56.8% CV)

    # Combined additive + proportional residual error (Sarashina 2005 Table 4)
    propSd <- 0.279 ; label("Proportional residual SD (fraction)")    # Sarashina 2005 Table 4: 27.9% (Table 5 sigma^2 = 0.0776; sqrt = 0.279)
    addSd  <- 0.425 ; label("Additive residual SD (ng/mL)")           # Sarashina 2005 Table 4: 0.425 ng/mL (Table 5 sigma^2 = 0.181; sqrt = 0.425)
  })
  model({
    # Covariate-effect multipliers / additive lag time (Sarashina 2005 final model equations)
    food_cl <- e_fed_cl ^ FED                   # 1 if fasted, 1.41 if fed
    food_vc <- e_fed_vc ^ FED                   # 1 if fasted, 1.75 if fed
    form_cl <- e_form_syrup_cl ^ FORM_SYRUP     # 1 if tablet, 1.06 if dry syrup
    tlag    <- e_fed_tlag * FED                 # 0 if fasted, 0.234 h if fed (additive)

    # Individual PK parameters. CL/F and V1/F retain the paper's LINEAR-in-WT form
    # via the (1 + e_wt_<p> * WT) factor: at WT=0 this collapses to exp(l<p>) (the
    # paper's intercept), and at WT=63 kg (adult mean) it yields the published
    # typical-value clearance.  Log-normal IIV applied multiplicatively.
    cl <- exp(lcl + etalcl) * (1 + e_wt_cl * WT) * food_cl * form_cl
    vc <- exp(lvc + etalvc) * (1 + e_wt_vc * WT) * food_vc
    q  <- exp(lq  + etalq)
    vp <- exp(lvp + etalvp)
    ka <- exp(lka + etalka)

    # Micro-rate constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment first-order absorption (Sarashina 2005 Methods Step 1)
    d/dt(depot)       = -ka * depot
    d/dt(central)     =  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) =                                k12 * central - k21 * peripheral1

    alag(depot) <- tlag

    # Observation. Dose in mg, volumes in L -> central/vc has units mg/L = ug/mL;
    # multiply by 1000 to match the paper's reported ng/mL plasma concentrations.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
