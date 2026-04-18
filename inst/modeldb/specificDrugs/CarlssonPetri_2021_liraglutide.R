CarlssonPetri_2021_liraglutide <- function() {
  description <- "Liraglutide PK model in adolescents (Carlsson Petri 2021)"
  reference <- "Carlsson Petri KC, Hale PM, Hesse D, Rathor N, Mastrandrea LD. Liraglutide pharmacokinetics and exposure-response in adolescents with obesity. Pediatric Obesity. 2021;16(10):e12799. doi:10.1111/ijpo.12799"
  units <- list(time = "hr", dosing = "mg", concentration = "TODO")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL/F and Vc/F; reference weight 100 kg per Equation 2.",
      source_name        = "WT"
    ),
    CHILD = list(
      description        = "Indicator for child age group",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (adult)",
      notes              = "Paired with ADOLESCENT; both 0 indicates adult reference. Paper's specific age cutoffs: TODO: from source paper.",
      source_name        = "CHILD"
    ),
    ADOLESCENT = list(
      description        = "Indicator for adolescent age group",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (adult)",
      notes              = "Paired with CHILD; both 0 indicates adult reference. Paper's specific age cutoffs: TODO: from source paper.",
      source_name        = "ADOLESCENT"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Source paper used SEXM (1 = male, 0 = female). Canonical SEXF inverts values; the model applies the original SEXM coefficient to (1 - SEXF) to preserve the reported effect magnitude and reference category.",
      source_name        = "SEXM"
    )
  )

  population <- list(
    n_subjects     = "TODO: from source paper",
    n_studies      = "TODO: from source paper",
    age_range      = "TODO: from source paper",
    age_median     = "TODO: from source paper",
    weight_range   = "TODO: from source paper",
    weight_median  = "TODO: from source paper",
    sex_female_pct = "TODO: from source paper",
    race_ethnicity = "TODO: from source paper",
    disease_state  = "Adolescents with obesity (pooled with child and adult reference groups).",
    dose_range     = "TODO: from source paper",
    regions        = "TODO: from source paper",
    notes          = "TODO: from source paper"
  )

  ini({
    lka <- fixed(log(0.0813)); label("Absorption rate (1/hr)")
    lcl <- log(1.01); label("Apparent clearance (L/h)")
    e_wt_cl <- 0.762; label("Body weight exponent on CL/F")
    e_sex_cl <- 1.12; label("Sex contrast (male/female) on CL/F")
    e_age_child_cl <- 1.11; label("Age contrast (child/adult) on CL/F")
    e_age_adolescent_cl <- 1.06; label("Age contrast (adolescent/adult) on CL/F")
    lvc <- fixed(log(13.8)); label("Apparent central volume of distribution (L)")
    e_wt_vc <- 0.587; label("Body weight exponent on Vc/F")

    etalcl ~ log(1.312)
    etalvc ~ log(1.317)
    propSd <- 0.433; label("Proportional residual error (fraction)")
  })
  model({
    ka <- exp(lka)
    cl_wt <- (WT / 100)^e_wt_cl # Equation 2 in the paper
    cl_sex <- (1 - SEXF)^e_sex_cl # Equation 3 in the paper; source used SEXM so (1 - SEXF) preserves the male indicator
    cl_age <- CHILD^e_age_child_cl * ADOLESCENT^e_age_adolescent_cl # Equation 4 in the paper
    cl <- exp(lcl + etalcl) * cl_wt * cl_sex * cl_age # Equation 1 in the paper
    vc_wt <- (WT / 100)^e_wt_vc # Not in the paper, based on Equation 2 in the paper
    vc <- exp(lvc + etalvc) * vc_wt # Equation 5 in the paper

    Cc <- linCmt()
    Cc ~ prop(propSd)
  })
}
