Xu_2011_sirukumab <- function() {
  description <- "Two-compartment population PK model for sirukumab (anti-IL-6 human IgG1 kappa monoclonal antibody, CNTO 136) in healthy adults following a single intravenous infusion, with first-order elimination from the central compartment and allometric body-weight scaling (Xu 2011)."
  reference <- "Xu Z, Bouman-Thio E, Comisar C, Frederick B, Van Hartingsveldt B, Marini JC, Davis HM, Zhou H. Pharmacokinetics, pharmacodynamics and safety of a human anti-IL-6 monoclonal antibody (sirukumab) in healthy subjects in a first-in-human study. Br J Clin Pharmacol. 2011;72(2):270-281. doi:10.1111/j.1365-2125.2011.03964.x"
  vignette <- "Xu_2011_sirukumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used for allometric scaling with reference weight 70 kg per Xu 2011 Table 4 and the population PK final-model equations. Exponents are 0.75 for CL and Q and 1 for V1 and V2 (fixed in the source).",
      source_name        = "WT"
    )
  )

  population <- list(
    n_subjects     = 34L,
    n_studies      = 1L,
    age_range      = "18-54 years",
    age_median     = "30 years",
    weight_range   = "49-99 kg",
    weight_median  = "71.3 kg",
    sex_female_pct = 16,
    race_ethnicity = c(White = 71, Black = 16, Asian = 9, Other = 4),
    disease_state  = "Healthy adult volunteers",
    dose_range     = "0.3, 1, 3, 6 (male and female cohorts), or 10 mg/kg as a single 10-15 min IV infusion",
    regions        = "Single-center study conducted in the United States",
    ada_positive_pct = 0,
    notes          = "Baseline demographics from Xu 2011 Table 1 and Section 'Subject demographic and baseline characteristics'. Double-blind, placebo-controlled, ascending single-dose first-in-human study (C0524T01). Forty-five subjects enrolled; 34 received sirukumab (6 per cohort for cohorts 1-5 in 0.3, 1, 3, 6 mg/kg male, and 6 mg/kg female groups; 4 in the 10 mg/kg cohort) and 11 received placebo. Population PK dataset was the 34 sirukumab-treated subjects. No subject developed antibodies to sirukumab."
  )

  ini({
    # Structural parameters - Xu 2011 Table 4 final-model estimates for a 70 kg reference subject.
    # A two-compartment model with zero-order IV input and first-order elimination from the
    # central compartment best fit the data versus one-, three-compartment, and Michaelis-Menten
    # alternatives (Xu 2011, Results / Population PK analysis).
    lcl <- log(0.364); label("Clearance CL for a 70 kg subject (L/day)")                         # Xu 2011 Table 4, theta1 (CL)
    lvc <- log(3.28);  label("Central volume of distribution V1 for a 70 kg subject (L)")        # Xu 2011 Table 4, theta2 (V1)
    lq  <- log(0.588); label("Inter-compartmental clearance Q for a 70 kg subject (L/day)")      # Xu 2011 Table 4, theta3 (Q)
    lvp <- log(4.97);  label("Peripheral volume of distribution V2 for a 70 kg subject (L)")     # Xu 2011 Table 4, theta4 (V2)

    # Allometric exponents fixed in the source: 0.75 for CL and Q, 1 for V1 and V2; reference
    # weight 70 kg. See Xu 2011 Results: "with the exponents for the allometric functions being
    # fixed to 0.75 for CL and Q and to 1 for V1 and V2".
    allo_cl <- fixed(0.75); label("Allometric exponent on CL (unitless)")                         # Xu 2011 Results, allometric scaling
    allo_q  <- fixed(0.75); label("Allometric exponent on Q (unitless)")                          # Xu 2011 Results, allometric scaling
    allo_v1 <- fixed(1);    label("Allometric exponent on V1 (unitless)")                         # Xu 2011 Results, allometric scaling
    allo_v2 <- fixed(1);    label("Allometric exponent on V2 (unitless)")                         # Xu 2011 Results, allometric scaling

    # IIV - Xu 2011 Table 4 reports "IIV (%)" under an exponential (log-normal) model for the
    # PK parameters. Converting CV% to the log-scale variance via omega^2 = log(CV^2 + 1):
    #   CL  IIV 24.3% -> omega^2 = log(0.243^2 + 1) = 0.057371
    #   V1  IIV 19.3% -> omega^2 = log(0.193^2 + 1) = 0.036572
    #   Q   IIV 53.4% -> omega^2 = log(0.534^2 + 1) = 0.250880
    #   V2  IIV 28.3% -> omega^2 = log(0.283^2 + 1) = 0.077043
    # The final-model table does not report IIV correlations, so IIV is treated as diagonal.
    etalcl ~ 0.057371  # Xu 2011 Table 4, IIV(%) on CL = 24.3
    etalvc ~ 0.036572  # Xu 2011 Table 4, IIV(%) on V1 = 19.3
    etalq  ~ 0.250880  # Xu 2011 Table 4, IIV(%) on Q  = 53.4
    etalvp ~ 0.077043  # Xu 2011 Table 4, IIV(%) on V2 = 28.3

    # Residual error - Xu 2011 Table 4 reports a combined additive + proportional model.
    # Proportional error variability (%): 21.7 -> propSd = 0.217 (fraction).
    # Additive error (ug/L): 0.0228. With observation units ug/mL (= mg/L), this converts to
    # 0.0228 ug/L x (1 mL / 1000 uL) ... use: 0.0228 ug/L x (1 L / 1000 mL) = 2.28e-5 ug/mL.
    propSd <- 0.217;   label("Proportional residual error (fraction)")                            # Xu 2011 Table 4, Proportional error variability (%) = 21.7
    addSd  <- 2.28e-5; label("Additive residual error (ug/mL)")                                   # Xu 2011 Table 4, Additive error (ug/L) = 0.0228
  })

  model({
    # Individual PK parameters with Xu 2011 allometric weight scaling (reference 70 kg).
    cl <- exp(lcl + etalcl) * (WT / 70)^allo_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^allo_v1
    q  <- exp(lq  + etalq)  * (WT / 70)^allo_q
    vp <- exp(lvp + etalvp) * (WT / 70)^allo_v2

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment model with IV dosing into the central compartment and first-order
    # elimination. Dose units mg, volume units L -> concentration mg/L = ug/mL.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
