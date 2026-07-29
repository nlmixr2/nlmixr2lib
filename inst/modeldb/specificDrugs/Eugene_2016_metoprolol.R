Eugene_2016_metoprolol <- function() {
  description <- "One-compartment population PK model for oral metoprolol tartrate with first-order absorption and lag time in elderly inpatients with multiple comorbidities; sex as the only covariate on apparent clearance (Eugene 2016)."
  reference <- "Eugene AR. Gender based Dosing of Metoprolol in the Elderly using Population Pharmacokinetic Modeling and Simulations. Int J Clin Pharmacol Toxicol. 2016;5(3):209-215. doi:10.19070/2167-910X-1600035"
  vignette <- "Eugene_2016_metoprolol"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    SEXF = list(
      description        = "Biological sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male) in the canonical column. The Eugene 2016 MONOLIX model uses female as its structural reference (see notes).",
      notes              = "Eugene 2016 reports the typical apparent clearance with female as the structural reference (CL_Female = 59.1 L/h, Table 1) and applies the male effect as a log-additive term betaCL_Male = 0.572 on the female baseline, so CL_Male = CL_Female * exp(0.572) = 105 L/h. To store under the canonical SEXF (1 = female, 0 = male) while preserving Eugene 2016's published female-reference CL, the effect is applied in model() as exp(e_sexf_cl * (1 - SEXF)); SEXF = 1 (female) yields factor 1 and SEXF = 0 (male) yields exp(0.572) ~= 1.77 (the male:female CL ratio).",
      source_name        = "SEX"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 7L,
    n_studies      = 1L,
    age_range      = "61-88 years (mean 75 +/- 9 years; males 65 +/- 6 years, females 83 +/- 13 years)",
    age_median     = "not reported (means reported instead)",
    weight_range   = "45-75 kg (mean 62 +/- 10 kg; males 65 +/- 6 kg, females 58 +/- 13 kg)",
    weight_median  = "not reported (means reported instead)",
    sex_female_pct = 42.9,
    race_ethnicity = c(White = 100),
    disease_state  = "Chronically ill elderly inpatients with multiple comorbidities",
    dose_range     = "Single 50 mg oral metoprolol tartrate tablet followed by 100 mL water after a minimum 10-h overnight fast (Eugene 2016 Methods; original cohort: Lundborg 1976).",
    regions        = "Sweden (University of Goteborg; concentration-time data digitized from Lundborg 1976)",
    notes          = "n = 7 of the original 10 study participants from Lundborg 1976 (3 females, 4 males) were retained for the popPK analysis. Eugene 2016 used MONOLIX 4.3.3 (SAEM-MCMC) to fit a one-compartment model with first-order absorption and a fixed lag time. Body weight was tested as a covariate but increased the AIC and reduced biological plausibility, so it was excluded from the final model (Discussion paragraph 1); sex on CL was the only significant covariate retained."
  )

  ini({
    # Structural parameters from Eugene 2016 Table 1 (final population pharmacokinetic
    # parameter estimates). CL is reported separately for females (59.1 L/h, the
    # structural reference of the MONOLIX model) and males (105 L/h); the published
    # log-additive male effect betaCL_Male = 0.572 reproduces CL_Male as
    # CL_Female * exp(0.572) = 59.1 * 1.772 = 104.7 L/h, matching Table 1.
    ltlag      <- log(0.469); label("Absorption lag time (hr)")                                   # Eugene 2016 Table 1: Tlag = 0.469 hr (RSE 4%)
    lka       <- log(0.235); label("First-order absorption rate constant (1/hr)")                # Eugene 2016 Table 1: Ka = 0.235 hr^-1 (RSE 8%)
    lvc       <- log(38);    label("Apparent central volume of distribution Vc/F (L)")           # Eugene 2016 Table 1: V = 38 L (RSE 95%)
    lcl       <- log(59.1);  label("Apparent clearance CL/F at female reference (L/hr)")         # Eugene 2016 Table 1: CL Females = 59.1 L/hr (RSE 12%); structural reference of the MONOLIX model
    e_sexf_cl <- 0.572;      label("Log-additive coefficient for male sex on CL (applied as (1 - SEXF))") # Eugene 2016 Table 1: betaCL_Male = 0.572 (RSE 24%, p = 2.70e-5); applied via (1 - SEXF) so SEXF = 0 (male) yields CL = CL_Female * exp(0.572) = 105 L/hr

    # Inter-individual variability (Eugene 2016 Table 1). The reported numeric
    # values are variances on the log scale (omega^2). The paper's reported CV%
    # column corresponds to the approximation CV ~= sqrt(omega^2), which matches
    # the table for the high-variance V parameter (omega^2 = 2.43, reported CV ~155%);
    # the exact log-normal CV = sqrt(exp(omega^2) - 1) is larger.
    etaltlag ~ 0.030  # Eugene 2016 Table 1: omega^2_Tlag = 0.030 (RSE 113%; reported CV ~17%)
    etalka  ~ 0.052  # Eugene 2016 Table 1: omega^2_Ka   = 0.052 (RSE 254%; reported CV ~23%)
    etalvc  ~ 2.43   # Eugene 2016 Table 1: omega^2_V    = 2.43  (RSE 27%;  reported CV ~155%; very high)
    etalcl  ~ 0.143  # Eugene 2016 Table 1: omega^2_CL   = 0.143 (RSE 39%;  reported CV ~38%)

    # Residual error: additive on the linear concentration scale (Eugene 2016
    # Table 1 footnote: 'sigma, variance of intra-individual deviation of fitted
    # metoprolol plasma concentration', reported in ng/mL as a single SD-style value).
    addSd <- 12.2; label("Additive residual error (ng/mL)")                                       # Eugene 2016 Table 1: sigma = 12.2 ng/mL (RSE 10%)
  })

  model({
    # Individual parameters with log-normal IIV; sex effect on CL is applied as
    # exp(e_sexf_cl * (1 - SEXF)) so that female (SEXF = 1) reproduces the paper's
    # CL_Female reference and male (SEXF = 0) reproduces the paper's CL_Male.
    tlag <- exp(ltlag + etaltlag)
    ka   <- exp(lka  + etalka)
    vc   <- exp(lvc  + etalvc)
    cl   <- exp(lcl  + e_sexf_cl * (1 - SEXF) + etalcl)

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    lag(depot) <- tlag

    # Dose is in mg and Vc is in L, so central / vc is in mg/L = ug/mL = 1000 ng/mL.
    # Multiply by 1000 to express Cc in ng/mL, matching the paper's reported
    # concentration units and the additive residual error magnitude (12.2 ng/mL).
    Cc <- 1000 * central / vc
    Cc ~ add(addSd)
  })
}
