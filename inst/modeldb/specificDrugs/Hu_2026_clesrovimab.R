Hu_2026_clesrovimab <- function() {
  description <- "Two-compartment population PK model for clesrovimab in preterm and full-term infants (Hu 2026)"
  reference <- "Hu Z, Hellmann F, Zang X, et al. Population Pharmacokinetics of Clesrovimab in Preterm and Full-Term Infants. Clin Pharmacol Ther. 2026;119(4):1036-1046. doi:10.1002/cpt.70199"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying; used for allometric scaling with reference weight 5 kg.",
      source_name        = "WT"
    ),
    PNA = list(
      description        = "Postnatal age (chronological time since birth)",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying; combined with GA to produce the adjusted postnatal age used in the CL maturation function.",
      source_name        = "PNA"
    ),
    GA = list(
      description        = "Gestational age at birth",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject; used together with PNA to adjust postnatal age for prematurity via (GA - 40)/4.345.",
      source_name        = "GA"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (White/Other race)",
      notes              = "Multiplicative race effect on CL/F relative to the White/Other reference group.",
      source_name        = "ASIAN"
    ),
    RACE_BLACK = list(
      description        = "Black / African American race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (White/Other race)",
      notes              = "Multiplicative race effect on CL/F relative to the White/Other reference group.",
      source_name        = "BLACK"
    ),
    RACE_MULTI = list(
      description        = "Multiracial race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (White/Other race)",
      notes              = "Multiplicative race effect on CL/F relative to the White/Other reference group.",
      source_name        = "MULTIRACIAL"
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
    disease_state  = "Preterm and full-term infants at risk for RSV (clesrovimab population PK analysis).",
    dose_range     = "TODO: from source paper",
    regions        = "TODO: from source paper",
    notes          = "TODO: from source paper; baseline demographics per Hu 2026."
  )

  ini({
    lka  <- log(0.286);  label("Absorption rate (Ka, 1/day)")
    lcl  <- log(0.0197); label("Apparent clearance at 5 kg reference weight, full maturation, White/Other race (CL/F, L/day)")
    lvc  <- log(0.514);  label("Apparent central volume of distribution at 5 kg reference weight (Vc/F, L)")
    lvp  <- log(0.316);  label("Apparent peripheral volume of distribution at 5 kg reference weight (Vp/F, L)")
    lq   <- log(0.0406); label("Apparent intercompartmental clearance at 5 kg reference weight (Q/F, L/day)")

    allo_cl <- 0.524; label("Allometric exponent on CL/F and Q/F (unitless)")
    allo_v  <- 0.662; label("Allometric exponent on Vc/F and Vp/F (unitless)")

    # Maturation parameters for CL (Hill-type function of adjusted postnatal age)
    beta_cl <- 0.579; label("Fractional CL maturation at birth for a full-term infant (unitless)")
    t50_cl  <- 20.3;  label("Maturation half-life for CL/F (months)")

    # Race effects on CL/F relative to White/Other reference
    e_asian       <- -0.0585; label("Race effect on CL/F: Asian vs White/Other (fraction)")
    e_black       <-  0.132;  label("Race effect on CL/F: Black/African American vs White/Other (fraction)")
    e_multiracial <-  0.0872; label("Race effect on CL/F: Multiracial vs White/Other (fraction)")

    # Inter-individual variability (omega^2 = log(CV^2 + 1))
    etalka ~ 0.05376  # 23.5% CV
    etalcl ~ 0.020524 # 14.4% CV
    # Note: Vc/F IIV had high shrinkage (71.9%) in the original analysis
    etalvc ~ 0.006572 # 8.12% CV

    # Residual error
    propSd <- 0.143; label("Proportional residual error (fraction)")
    addSd  <- 0.231; label("Additive residual error (ug/mL)")
  })
  model({
    # Adjusted postnatal age for prematurity (negative values valid for preterm at birth)
    # GA - 40 converts gestational age deficit from weeks to months via 4.345 days/month
    AAGEADJ <- PNA + (GA - 40) / 4.345

    # Maturation function for CL (approaches 1 at full maturity)
    maturation_cl <- 1 - (1 - beta_cl) * exp(-AAGEADJ * log(2) / t50_cl)

    # Race effect on CL (multiplicative; White/Other is reference = 1)
    race_cl <- 1 + e_asian * RACE_ASIAN + e_black * RACE_BLACK + e_multiracial * RACE_MULTI

    ka  <- exp(lka + etalka)
    cl  <- exp(lcl + etalcl) * (WT / 5)^allo_cl * maturation_cl * race_cl
    vc  <- exp(lvc + etalvc) * (WT / 5)^allo_v
    vp  <- exp(lvp)          * (WT / 5)^allo_v
    q   <- exp(lq)           * (WT / 5)^allo_cl

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Concentration: dose in mg, volume in L -> mg/L = ug/mL
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
