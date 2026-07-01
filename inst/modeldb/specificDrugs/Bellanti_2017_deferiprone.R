Bellanti_2017_deferiprone <- function() {
  description <- "One-compartment population PK model for the oral iron chelator deferiprone in paediatric patients aged <6 years with transfusion-dependent haemoglobinopathies, with first-order absorption and first-order elimination and fixed allometric scaling of clearance and volume on body weight (Bellanti 2017)."
  reference   <- "Bellanti F, Del Vecchio GC, Putti MC, Maggio A, Filosa A, Cosmi C, Mangiarini L, Spino M, Connelly J, Ceci A, Della Pasqua O; DEEP Consortium. Population pharmacokinetics and dosing recommendations for the use of deferiprone in children younger than 6 years. Br J Clin Pharmacol. 2017 Mar;83(3):593-602. doi:10.1111/bcp.13134"
  vignette    <- "Bellanti_2017_deferiprone"
  units       <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight (kg).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Fixed allometric scaling on both CL/F (exponent 0.75) and V/F (exponent 1), centred on the median paediatric body weight of 16 kg (Bellanti 2017 Table 2 and Table 1). Reference weight of 16 kg is inferred from arithmetic consistency between Table 2 (CL/F = 8.3 L/h) and Table 3 (median AUC 340.6 umol/L*h at 25 mg/kg in a cohort with mean WT 16 kg) - see vignette Assumptions and deviations.",
      source_name        = "Weight"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age (years).",
      units       = "years",
      type        = "continuous",
      notes       = "Tested in univariate covariate screen (weight, height, age, sex) on Ka, CL/F and V/F (Bellanti 2017 Results, Covariate analysis). Not retained; body weight alone was sufficient once fixed allometric scaling was applied."
    ),
    HT = list(
      description = "Standing height (cm).",
      units       = "cm",
      type        = "continuous",
      notes       = "Tested in univariate covariate screen (Bellanti 2017 Results). Not retained in the final model."
    ),
    SEXF = list(
      description = "Biological sex indicator (1 = female, 0 = male).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested in univariate covariate screen (Bellanti 2017 Results). Not retained; the 9 male / 9 female cohort did not show a sex effect on Ka, CL/F or V/F."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 18,
    n_studies      = 1,
    age_range      = "1.2-5.9 years",
    age_median     = "3.4 years",
    weight_range   = "11-22.5 kg",
    weight_median  = "15.8 kg",
    weight_mean    = "16.08 kg (SD 3.18 kg)",
    height_median  = "99.2 cm (range 83-117 cm)",
    sex_female_pct = 9 / 18 * 100,
    disease_state  = "Transfusion-dependent haemoglobinopathies (16 beta-thalassaemia major, 2 thalasso-drepanocytosis) in children younger than 6 years.",
    dose_range     = "Single oral deferiprone dose randomised to one of three levels: 8.3, 16.7 or 33.3 mg/kg (6 subjects per dose group), administered as an 80 mg/mL solution.",
    regions        = "Multicentre paediatric haematology units in Italy (DEEP-1 study, EudraCT 2012-000658-67, NCT01740713); analysis performed at the Leiden Academic Centre for Drug Research, Netherlands.",
    notes          = "Randomised, single-blind, single-dose paediatric PK study (DEEP-1) sponsored by the DEEP Consortium. Sparse sampling with a maximum of five post-dose plasma samples per subject over 8 h following an optimised design. Bioanalytical method: HPLC-UV, LLOQ 0.238 umol/L (0.033 ug/mL). NONMEM v7.2 with informative priors on Ka (8.2/h with uncertainty 4.02) and on the CL/F and V/F IIV block (from the Bellanti 2014 healthy-adult analysis, doi:10.1111/bcp.12473); 54 degrees of freedom for the BSV prior. Bootstrap analysis performed in PsN v3.5.3. See Bellanti 2017 Methods."
  )

  ini({
    # Structural parameters - typical values at reference body weight 16 kg (paediatric
    # cohort mean; reference weight not stated explicitly in the paper, see vignette
    # Assumptions and deviations for the arithmetic derivation from Table 2 CL/F vs
    # Table 3 AUC bridging results).
    lcl <- log(8.3);  label("Apparent clearance CL/F (L/h) at reference WT 16 kg")  # Bellanti 2017 Table 2 (CL/F Estimate)
    lvc <- log(18.7); label("Apparent volume of distribution V/F (L) at reference WT 16 kg")  # Bellanti 2017 Table 2 (V/F Estimate)
    lka <- log(9.13); label("First-order absorption rate constant Ka (1/h)")  # Bellanti 2017 Table 2 (Ka Estimate; estimated with informative prior 8.2/h from Bellanti 2014)

    # Fixed allometric exponents on body weight (both marked FIX in Bellanti 2017 Table 2).
    allo_cl <- fixed(0.75); label("Fixed allometric exponent on CL/F (unitless)")  # Bellanti 2017 Table 2 (WT on CL/F, 0.75 FIX)
    allo_v  <- fixed(1);    label("Fixed allometric exponent on V/F (unitless)")   # Bellanti 2017 Table 2 (WT on V/F, 1 FIX)

    # IIV - block correlation between etalcl and etalvc; no IIV on Ka in the final model.
    # Variances (omega^2) and covariance reported in Bellanti 2017 Table 2. The BSV block
    # was estimated with an informative Normal-Inverse-Wishart prior derived from the
    # Bellanti 2014 healthy-adult population PK model (54 degrees of freedom).
    # IIV CL/F: 0.0644 (CV 25.4%), IIV V/F: 0.0392 (CV 19.8%), block CL-V covariance: 0.031
    # (implied correlation = 0.031 / sqrt(0.0644 * 0.0392) = 0.617).
    etalcl + etalvc ~ c(0.0644,
                        0.031, 0.0392)  # Bellanti 2017 Table 2 (IIV CL/F, Block CL-V, IIV V/F)

    # Residual error: proportional model reported as "Error (prop) 0.0953" in Bellanti
    # 2017 Table 2. NONMEM $SIGMA convention treats this as a variance (sigma^2), so the
    # effective proportional residual SD is sqrt(0.0953) = 0.3087 (approx. 30.9% CV).
    # The methods section describes the residual as Y_ij = F_ij + epsilon_ij * W with
    # W as a "proportional weighing factor for epsilon"; unlike the Bellanti 2014 adult
    # analysis, the 2017 paediatric Table 2 does not report a separate weighting theta,
    # so W is taken to be F_ij (standard proportional error) and 0.0953 is the sole
    # residual variance parameter. Documented in vignette Assumptions and deviations.
    propSd <- 0.3087; label("Proportional residual SD (fraction); sqrt(0.0953)")  # Bellanti 2017 Table 2 (Error (prop) 0.0953)
  })

  model({
    # Individual PK parameters with fixed allometric scaling on body weight
    # (reference weight 16 kg = paediatric cohort mean).
    cl <- exp(lcl + etalcl) * (WT / 16)^allo_cl
    vc <- exp(lvc + etalvc) * (WT / 16)^allo_v
    ka <- exp(lka)

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
