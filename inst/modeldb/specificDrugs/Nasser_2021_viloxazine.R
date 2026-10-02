Nasser_2021_viloxazine <- function() {
  description <- "Joint parent (viloxazine) + metabolite (5-hydroxyviloxazine glucuronide, 5-HVLX-gluc) population PK model for once-daily oral viloxazine extended-release capsules in children (6-11 years) and adolescents (12-17 years) with ADHD (Nasser 2021, pooled phase 3 studies P301-P304). Viloxazine is one-compartment with first-order absorption (ka = 0.068 1/h, so the terminal phase is absorption-limited) and two parallel first-order elimination routes from the central compartment: formation of 5-HVLX-gluc (CLV) and all remaining viloxazine elimination (CLL). 5-HVLX-gluc is one-compartment with first-order elimination (CLM) and shares the viloxazine apparent volume (assumed by the authors for identifiability). Body weight enters as power functions centred on the 36.35 kg median on the shared volume, CLV and CLM; F is fixed to 1."
  reference <- "Nasser A, Gomeni R, Wang Z, Kosheleff AR, Xie L, Adeojo LW, Schwabe S. Population Pharmacokinetics of Viloxazine Extended-Release Capsules in Pediatric Subjects With Attention Deficit/Hyperactivity Disorder. J Clin Pharmacol. 2021;61(12):1626-1637. doi:10.1002/jcph.1940."
  vignette <- "Nasser_2021_viloxazine"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power function P = Pref * (WT / 36.35)^g on the shared apparent volume V2/F = V3/F (g = 0.78), on the viloxazine metabolic clearance CLV (g = 0.59) and on the 5-HVLX-gluc clearance CLM (g = 0.70); the viloxazine clearance CLL carries no weight effect (Methods 'Covariate analysis' equation; Table 2; Table S1 model 10). 36.35 kg is the population median used as typCov. Cohort mean 44.5 (SD 16.5) kg, range 20-92.5 kg (Table 1; Results 'Subject Characteristics').",
      source_name = "WT"
    )
  )

  # Covariates Nasser 2021 screened (Methods 'Covariate analysis') but did not
  # retain. Age was formally tested on V2 (Table S1 model 1b: not retained
  # alongside weight); the rest did not pass the graphical/EBE screen. BMI was
  # set aside because of its collinearity with weight (r = 0.87). Documentation
  # only -- none is referenced in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Formally tested jointly with weight on V2 (Table S1 model 1b, dOFV -0.29) and not retained; the effect of age was not significant after controlling for body weight (Discussion 'Body Weight'). Median 11 years used as typCov during testing. Cohort mean 11.2 (SD 3.1) years (Table 1)."
    ),
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      notes = "Screened, not retained. Cohort mean 149.0 (SD 17.9) cm (Table 1)."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Graphically associated with model parameters but set aside in favour of body weight because of collinearity (Pearson r = 0.87). Cohort mean 19.3 (SD 3.4) kg/m^2 (Table 1)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened, not retained. 158 of 495 subjects (31.9 pct) female (Table 1)."
    ),
    RACE_BLACK = list(
      description = "Black or African American race indicator",
      units = "(binary)",
      type = "binary",
      notes = "Race screened, not retained. 214 of 495 (43.2 pct) Black or African American, 258 (52.1 pct) White, 23 (4.6 pct) Other (Table 1)."
    ),
    RACE_HISPANIC = list(
      description = "Hispanic / Latino ethnicity indicator",
      units = "(binary)",
      type = "binary",
      notes = "Ethnicity screened, not retained. 112 of 495 (22.6 pct) Hispanic or Latino (Table 1)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened, not retained. Cohort mean 15 (SD 7) U/L (Table 1)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened, not retained. Cohort mean 24 (SD 6) U/L (Table 1)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened, not retained. Cohort mean 0.59 (SD 0.17) mg/dL (Table 1)."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "viloxazine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "viloxazine", units = "mg", specimen = "plasma", verified = TRUE),
    central_gluc = list(
      analyte = "5-hydroxyviloxazine glucuronide (5-HVLX-gluc)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 495L,
    n_studies = 4L,
    age_range = "6-17 years (mean 11.2, SD 3.1); children 6-11 years n = 263, adolescents 12-17 years n = 232",
    weight_range = "20-92.5 kg (mean 44.5, SD 16.5; median 36.35); children median 31.5 kg, adolescents median 57.25 kg",
    sex_female_pct = 31.9,
    race_ethnicity = "White 52.1%, Black or African American 43.2%, Other 4.6%; Hispanic or Latino 22.6%",
    disease_state = "Attention deficit/hyperactivity disorder (ADHD), viloxazine ER monotherapy",
    dose_range = "Viloxazine ER 100, 200, 400 or 600 mg orally once daily (children 100-400 mg/day in P301/P303; adolescents 200-600 mg/day in P302/P304), after 1-3 weeks of titration. 86 subjects at 100 mg, 197 at 200 mg, 164 at 400 mg, 48 at 600 mg.",
    regions = "United States (multicentre)",
    n_observations = "Sparse, flexible sampling at steady state: up to 5 samples per subject (predose and 1, 2, 4 and 6 h post-dose) of both viloxazine and 5-HVLX-gluc. LLOQ 0.01 ug/mL (viloxazine) and 0.005 ug/mL (5-HVLX-gluc).",
    notes = "Pooled phase 3 randomized, double-blind, placebo-controlled trials P301 (NCT03247530), P302 (NCT03247517), P303 (NCT03247543) and P304 (NCT03247556) (Methods 'Data Source'; Table 1). NONMEM 7.4, FOCE-I."
  )

  ini({
    # Final model (Table S1 model 10) estimates from Nasser 2021 Table 2,
    # column 'Model Estimate' (mean +/- SE (RSE)).
    lka <- log(0.068); label("Viloxazine first-order absorption rate constant ka (1/h)") # Table 2: ka = 0.068 +/- 0.0028 1/h (RSE 4.10)
    lvc <- log(14.6); label("Apparent volume of distribution of viloxazine, shared by 5-HVLX-gluc (V2/F = V3/F), at WT = 36.35 kg (L)") # Table 2: V2/F = 14.60 +/- 0.67 L (RSE 4.60); Methods: V2/F and V3/F assumed identical
    lcl_nonmet <- log(0.87); label("Apparent viloxazine clearance CLL through pathways other than 5-HVLX-gluc formation (L/h)") # Table 2: CLL = 0.87 +/- 0.18 L/h (RSE 21.30)
    lcl_met <- log(4.72); label("Apparent viloxazine metabolic clearance CLV forming 5-HVLX-gluc, at WT = 36.35 kg (L/h)") # Table 2: CLV = 4.72 +/- 0.21 L/h (RSE 4.40)
    lcl_gluc <- log(6.75); label("Apparent 5-HVLX-gluc elimination clearance CLM at WT = 36.35 kg (L/h)") # Table 2: CLM = 6.75 +/- 0.31 L/h (RSE 4.60)

    e_wt_vc <- 0.78; label("Power exponent of body weight (WT / 36.35 kg) on the shared volume V2/F = V3/F (unitless)") # Table 2: WT,V = 0.78 +/- 0.08 (RSE 10.40)
    e_wt_cl_met <- 0.59; label("Power exponent of body weight (WT / 36.35 kg) on CLV (unitless)") # Table 2: WT,CLV = 0.59 +/- 0.07 (RSE 12.00)
    e_wt_cl_gluc <- 0.70; label("Power exponent of body weight (WT / 36.35 kg) on CLM (unitless)") # Table 2: WT,CLM = 0.70 +/- 0.07 (RSE 9.40)

    # IIV: Table 2 'Random effect' rows, read as log-scale variances (NONMEM
    # OMEGA). The variance reading reproduces the paper's own Monte Carlo
    # medians in Tables S2 and S3; the SD reading overshoots them (see the
    # vignette). Lognormal IIV on every parameter (Methods).
    etalvc ~ 0.10 # Table 2 random effect V2/F = 0.10 +/- 0.03 (RSE 28.00)
    etalcl_nonmet ~ 3.03 # Table 2 random effect CLL = 3.03 +/- 0.50 (RSE 16.40)
    etalka ~ 0.17 # Table 2 random effect ka = 0.1700 +/- 0.0266 (RSE 15.60)
    etalcl_met ~ 0.11 # Table 2 random effect CLV = 0.11 +/- 0.01 (RSE 11.40)
    etalcl_gluc ~ 0.08 # Table 2 random effect CLM = 0.08 +/- 0.01 (RSE 14.90)

    # Residual error: combined additive + proportional (Methods). Table 2
    # prints a single Additive / Proportional pair for the joint fit; it is
    # applied to both analytes here. The rows are read as standard deviations
    # (the Figure 2 90 pct prediction interval is reproduced on that reading,
    # not on the variance reading; see the vignette).
    addSd <- 0.12; label("Additive residual SD for viloxazine (ug/mL)") # Table 2: residual Additive = 0.12 +/- 0.01 (RSE 5.50)
    propSd <- 0.29; label("Proportional residual SD for viloxazine (fraction)") # Table 2: residual Proportional = 0.29 +/- 0.01 (RSE 3.50)
    addSd_gluc <- 0.12; label("Additive residual SD for 5-HVLX-gluc (ug/mL)") # Table 2: residual Additive = 0.12 (single pair reported for both analytes)
    propSd_gluc <- 0.29; label("Proportional residual SD for 5-HVLX-gluc (fraction)") # Table 2: residual Proportional = 0.29 (single pair reported for both analytes)
  })

  model({
    # Covariate model (Methods 'Covariate analysis'):
    # P = Pref * (Cov / typCov)^g, typCov = median body weight 36.35 kg.
    wt_ratio <- WT / 36.35

    ka <- exp(lka + etalka)
    vc <- exp(lvc + etalvc) * wt_ratio^e_wt_vc
    cl_nonmet <- exp(lcl_nonmet + etalcl_nonmet)
    cl_met <- exp(lcl_met + etalcl_met) * wt_ratio^e_wt_cl_met
    cl_gluc <- exp(lcl_gluc + etalcl_gluc) * wt_ratio^e_wt_cl_gluc
    # V3/F (5-HVLX-gluc) identical to V2/F (viloxazine) (Methods)
    vc_gluc <- vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - (cl_met + cl_nonmet) / vc * central
    # 5-HVLX-gluc formed 1:1 on a mass basis; the paper states no
    # molecular-weight correction, so the metabolite state is in
    # viloxazine-mass equivalents and the fitted clearances absorb the ratio
    d/dt(central_gluc) <- cl_met / vc * central - cl_gluc / vc_gluc * central_gluc

    # F = 1 (Methods); dose mg / volume L = mg/L = ug/mL
    Cc <- central / vc
    Cc_gluc <- central_gluc / vc_gluc

    Cc ~ add(addSd) + prop(propSd)
    Cc_gluc ~ add(addSd_gluc) + prop(propSd_gluc)
  })
}
