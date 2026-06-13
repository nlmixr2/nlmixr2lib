Barau_2012_mycophenolic_acid <- function() {
  description <- "One-compartment population PK model for mycophenolic acid (MPA, active moiety of mycophenolate mofetil MMF) after oral MMF dosing in paediatric liver transplant recipients (Barau 2012). First-order absorption and first-order elimination, with diagonal (uncorrelated) inter-individual variability on ka, CL/F, and V/F and proportional residual error. Two covariates are retained in the final model: a linear-with-age effect on ka of the form ka_TV = 3.9 - 2.2 * (AGE / 8.65 years), so ka declines from 3.9 1/h at AGE = 0 to 1.7 1/h at the cohort median age of 8.65 years; and a power-on-binary effect on V/F of the form V/F = 64.7 L * 2.3^POSTTX_EARLY, where POSTTX_EARLY = 1 within the first 6 months post-transplant (POD <= 180 days) and 0 thereafter, so V/F is 64.7 L in the stable post-transplant period and 148.8 L in the immediate post-transplant period (paper attributes the volume increase to the higher unbound MPA fraction associated with low serum albumin in the immediate post-transplant period). Apparent clearance CL/F = 12.7 L/h carries no retained covariate effect in the final model. Enterohepatic recirculation, the MPAG metabolite compartment, and protein binding are not modelled here -- the paper attributes the absence of secondary peaks to surgical removal of the gallbladder in the liver-transplant recipients."
  reference <- paste(
    "Barau C, Furlan V, Debray D, Taburet AM, Barrail-Tran A.",
    "Population pharmacokinetics of mycophenolic acid and dose optimization",
    "with limited sampling strategy in liver transplant children.",
    "Br J Clin Pharmacol. 2012;74(3):515-524.",
    "doi:10.1111/j.1365-2125.2012.04213.x",
    sep = " "
  )
  vignette <- "Barau_2012_mycophenolic_acid"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    AGE = list(
      description        = "Subject age in years",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject (study-specific approximation; in a long follow-up the age is updated occasion-to-occasion, but in this paper's modelling window subjects are anchored at study entry). Centred at the cohort median 8.65 years (Barau 2012 Results: 'when the age was 8.65 years'); used in the additive ka covariate equation ka_TV = 3.9 - 2.2 * (AGE / 8.65). The cohort age range is 1.1-15.2 years (model-building set) and 1.1-18.0 years overall.",
      source_name        = "AGE"
    ),
    POD = list(
      description        = "Days post-transplantation",
      units              = "days",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying within subject; integer- or fractional-day valued. Paper's covariate equation collapses POD to a binary 6-month indicator: V/F is multiplied by 2.3 when the observation falls in the immediate post-transplant period (POD <= 180 days, equivalent to <= 6 months) and by 1 otherwise. The dichotomization is derived inside model() (posttx_early <- POD <= 180), so users supply POD directly in their dataset. The 180-day cutoff is the standard pharmacology approximation to the paper's '6 months' threshold; the paper does not provide a NONMEM control stream so the exact day-count translation is not author-stated. Median time since transplantation in the cohort was 17.2 months (range 0.2-188.5 months, Barau 2012 Methods 'Patients and study design').",
      source_name        = "POD"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 16L,
    n_studies      = 1L,
    age_range      = "1.1-15.2 years (model-building set; cohort overall 1.1-18.0 years across 28 patients)",
    age_median     = "8.7 years (model-building set; cohort overall median 8.65 years)",
    weight_range   = "9.3-49.2 kg (model-building set)",
    weight_median  = "23.8 kg (model-building set)",
    sex_female_pct = 50,
    race_ethnicity = "Not reported in source paper.",
    disease_state  = "Paediatric liver transplant recipients (n = 28 in the full cohort; n = 16 used for model building and n = 12 for model validation). Indications for transplantation: biliary atresia (14), fulminant hepatitis (8), progressive familial intrahepatic cholestasis (3), and Alagille syndrome / cystic fibrosis / Wilson disease (1 each).",
    dose_range     = "Oral mycophenolate mofetil (MMF) 186-594 mg/m^2 (= 6.1-22.2 mg/kg) twice daily; median starting dose 380 mg/m^2 (= 13.0 mg/kg) twice daily. MMF dose adjustments target MPA AUC(0,12 h) between 30 and 60 mg/L*h.",
    regions        = "France (Hopitaux Universitaires Paris-Sud, Le Kremlin Bicetre / Chatenay-Malabry).",
    co_medication  = "Tacrolimus (n = 23) or ciclosporin (n = 5); steroids (n = 14). Source paper found no MPA-exposure difference by calcineurin-inhibitor cotreatment (Discussion attributes the absence of a co-medication effect to negligible enterohepatic recirculation in the cohort after gall-bladder removal).",
    pod_range      = "Time since transplantation 0.2-188.5 months (median 17.2 months); 7 of the 16 model-building patients were sampled at <= 6 months post-transplant and 9 at > 6 months (Barau 2012 Table 3).",
    biochemistry   = "Cohort serum albumin median 31.7 g/L (range 17.2-35.0 g/L); ALAT median 93 IU/L; ASAT 64 IU/L; total bilirubin 20 umol/L; serum creatinine 37 umol/L; creatinine clearance 155 mL/min (Schwartz formula). Albumin is lower in the immediate post-transplant period (median 26.8 vs 33.0 g/L), associated with a higher unbound MPA fraction (median 2.8% vs 1.0%) and the 2.3-fold higher V/F. No biochemistry covariate is retained in the final model -- the post-transplant-period indicator captures the album/unbound-fraction effect parsimoniously.",
    notes          = "Patients underwent intensive PK sampling (predose plus 0.5, 1, 2, 4, 6, 8 h post-dose) during a 12-h dosing interval. The 16 model-building patients contributed one occasion each; the 12 model-validation patients contributed 26 intensive sets used for Bayesian individual-AUC estimation in Adapt II. The paper notes high inter-individual variability on ka (308% CV) likely reflects irregular individual concentration vs. time profiles and the absence of regular enterohepatic secondary peaks in the liver-transplant cohort."
  )

  ini({
    # Structural PK -- Barau 2012 Table 2 final-model estimates. Time in hours;
    # apparent clearance CL/F in L/h; apparent volume V/F in L; absorption rate
    # ka in 1/h. The reference subject for Table 2 is a paediatric liver
    # transplant recipient at the cohort median age 8.65 years in the stable
    # post-transplant period (> 6 months, POSTTX_EARLY = 0); the in-file
    # covariate equations below recreate individual-specific values inside
    # model().

    # ka structural parameter and age covariate. Paper's covariate equation
    # (Barau 2012 Results): ka_i = (kaTV - beta_age_ka * AGE / 8.65) * exp(eta_ka).
    # At AGE = 8.65 years (cohort median): ka_TV = 3.9 - 2.2 * 1 = 1.7 1/h (matches
    # the paper's prose 'ka, estimated at 1.7 h^-1 at age 8.7 years'). At AGE = 0
    # the typical-value intercept is 3.9 1/h. NB: the published linear-in-age
    # form goes structurally negative for AGE > 15.3 years (3.9 - 2.2 * 15.3 /
    # 8.65 = 0); within the model-building cohort (max age 15.2 years) the
    # typical-value ka remains positive, but the validation cohort extends to
    # 18.0 years where the typical-value ka is mathematically negative. See
    # vignette Assumptions and deviations for documentation of this limitation.
    lka       <- log(3.9)  ; label("Absorption rate intercept kaTV at AGE = 0 (1/h)")                                      # Barau 2012 Table 2 final ka = 3.9 1/h (at AGE = 0; combined with -beta_age_ka * AGE/8.65 yields 1.7 1/h at the cohort median age 8.65 y)
    e_age_ka  <- 2.2       ; label("Slope of ka decrease with AGE/8.65 (1/h per (AGE / 8.65 years))")                       # Barau 2012 Table 2 beta_age_ka = 2.2 (additive age effect; ka_TV = 3.9 - 2.2 * AGE/8.65)

    # CL/F (no retained covariates in the final model).
    lcl       <- log(12.7) ; label("Apparent oral clearance CL/F (L/h)")                                                    # Barau 2012 Table 2 final CL/F = 12.7 L/h

    # V/F structural parameter and post-transplant-period covariate. Paper's
    # covariate equation: V/F_i = V/F_TV * beta_time^POSTTX_EARLY * exp(eta_VF),
    # where POSTTX_EARLY = 1 if POD <= 180 days (= 6 months post-transplant) and
    # 0 otherwise. V/F_TV is the value in the stable post-transplant period
    # (> 6 months): 64.7 L. In the immediate post-transplant period the typical
    # V/F is 64.7 * 2.3 = 148.81 L.
    lvc                <- log(64.7) ; label("Apparent volume of distribution V/F in the stable post-transplant period (> 6 months) (L)")  # Barau 2012 Table 2 final V/F = 64.7 L (reference: > 6 months post-transplant)
    e_posttx_early_vc  <- 2.3       ; label("Fold-change in V/F during the immediate post-transplant period (<= 6 months) vs > 6 months reference (unitless)") # Barau 2012 Table 2 beta_time_V/F = 2.3 (V/F is 2.3-fold higher in the immediate period; 148.8 L vs 64.7 L)

    # Inter-individual variability -- Barau 2012 Table 2 reports IIV as %CV for
    # an exponential random-effect model q_i = q * exp(eta_i) with a diagonal
    # variance matrix (Methods: 'a diagonal variance matrix Omega was chosen').
    # Convert CV% to log-scale variance via omega^2 = log(CV^2 + 1):
    #   ka    CV 308.4% -> log(3.084^2 + 1) = 2.3524
    #   CL/F  CV  28.4% -> log(0.284^2 + 1) = 0.0776
    #   V/F   CV  41.8% -> log(0.418^2 + 1) = 0.1610
    etalka  ~ 2.3524 # Barau 2012 Table 2 IIV ka   = 308.4% CV -> log(1 + 3.084^2)
    etalcl  ~ 0.0776 # Barau 2012 Table 2 IIV CL/F =  28.4% CV -> log(1 + 0.284^2)
    etalvc  ~ 0.1610 # Barau 2012 Table 2 IIV V/F  =  41.8% CV -> log(1 + 0.418^2)

    # Residual unexplained variability -- Barau 2012 Table 2 reports a
    # proportional residual error model (Methods: 'Residual error was modelled
    # on a proportional error model'). sigma = 59.6% per Table 2 (s (%) =
    # 59.6); encoded directly as the SD of the relative error.
    propSd <- 0.596 ; label("Proportional residual error (fraction)")                                                       # Barau 2012 Table 2 final s = 59.6% (proportional error model)
  })

  model({
    # Binary post-transplant-period indicator at 6 months (180 days). Paper's
    # covariate equation uses time_post_transplantation = 1 in the immediate
    # post-transplant period (<= 6 months) and 0 thereafter; per the canonical
    # POD entry in covariate-columns.md, POD is supplied in days and the
    # dichotomisation is derived here.
    posttx_early <- (POD <= 180)

    # Individual PK parameters. ka follows the paper's additive age structure
    # (Barau 2012 Results): ka_i = (kaTV - beta_age_ka * AGE / 8.65) *
    # exp(eta_ka). V/F follows the paper's power-on-binary post-transplant
    # structure: V/F_i = V/F_TV * beta_time^POSTTX_EARLY * exp(eta_V/F). CL/F
    # carries no covariates in the final model.
    ka <- (exp(lka) - e_age_ka * (AGE / 8.65)) * exp(etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc) * e_posttx_early_vc ^ posttx_early

    # Micro-constants for the one-compartment first-order disposition.
    kel <- cl / vc

    # One-compartment oral PK with first-order absorption. Dose lands in
    # depot; bioavailability is absorbed into V/F and CL/F (apparent
    # parameters) as is standard for an oral-only popPK fit. MMF is hydrolysed
    # to MPA in the gut wall and plasma; the model treats the depot as MPA-
    # equivalent under the assumption of 1:1 molar conversion at the doses
    # used (the apparent V/F and CL/F already absorb the MMF-to-MPA conversion
    # factor through the bioavailability F).
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Plasma MPA concentration (mg/L). Dose is in mg of MMF -- the apparent
    # parameters V/F and CL/F implicitly fold in the MMF -> MPA molar
    # conversion factor (MPA MW 320.3, MMF MW 433.5) via the bioavailability
    # F, so no explicit unit conversion is applied here. See vignette
    # Assumptions and deviations.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
