PK_2cmt_mAb_Davda_2014 <- function() {
  description <- "Two compartment PK model with linear clearance for average monoclonal antibodies (Davda 2014)"
  reference <- "Davda JP, Dodds MG, Gibbs MA, Wisdom W, Gibbs JP. A model-based meta-analysis of monoclonal antibody pharmacokinetics to guide optimal first-in-human study design. MAbs. 2014;6(4):1094-1102. doi:10.4161/mabs.29095"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling with reference weight 70 kg; exponent 0.865 on CL and Q, 0.957 on Vc and Vp.",
      source_name        = "WT"
    )
  )

  population <- list(
    n_subjects     = "TODO: from source paper",
    n_studies      = "TODO: from source paper",
    age_range      = "TODO: from source paper",
    weight_range   = "TODO: from source paper",
    sex_female_pct = "TODO: from source paper",
    race_ethnicity = "TODO: from source paper",
    disease_state  = "Model-based meta-analysis across 18 therapeutic monoclonal antibodies; consensus average-mAb PK in humans.",
    dose_range     = "TODO: from source paper",
    regions        = "TODO: from source paper",
    notes          = "Meta-analysis of published population PK parameters across multiple therapeutic mAbs; detailed population-level demographics pooled from the underlying studies (see Davda 2014 Tables 1-2)."
  )

  ini({
    # Structural parameters — Davda 2014 reference weight = 70 kg
    lfdepot <- log(0.744); label("Subcutaneous bioavailability (fraction)")                            # Davda 2014 Table 3: F1 = 0.744
    lka     <- log(0.282); label("Absorption rate (Ka, 1/day)")                                        # Davda 2014 Table 3: Ka = 0.282 /day
    lcl     <- log(0.200); label("Clearance at reference weight (CL, L/day)")                          # Davda 2014 Table 3: CL = 0.200 L/day
    lvc     <- log(3.61);  label("Central volume of distribution at reference weight (Vc, L)")         # Davda 2014 Table 3: V1 = 3.61 L
    lvp     <- log(2.75);  label("Peripheral volume of distribution at reference weight (Vp, L)")      # Davda 2014 Table 3: V2 = 2.75 L
    lq      <- log(0.747); label("Intercompartmental clearance at reference weight (Q, L/day)")        # Davda 2014 Table 3: Q = 0.747 L/day

    allocl  <- 0.865; label("Allometric exponent on clearance and intercompartmental clearance (unitless)")  # Davda 2014 Table 3
    allov   <- 0.957; label("Allometric exponent on volumes of distribution (unitless)")                      # Davda 2014 Table 3

    # IIV — variances reported directly (not CV% shortcuts); correlated block for CL, Vc, Vp
    etalfdepot ~ 0
    etalka ~ 0.416
    etalcl + etalvc + etalvp ~ c(0.0987,
                                 0.0786, 0.116,
                                 0.0377, 0.0619, 0.0789)
    etalq ~ 0.699

    propSd <- 0.144; label("Proportional residual error (fraction)")
  })
  model({
    # WT is body weight in kg
    fdepot <- exp(lfdepot + etalfdepot)
    ka     <- exp(lka + etalka)
    wtnorm <- log(WT / 70)
    cl     <- exp(lcl + allocl * wtnorm + etalcl)
    q      <- exp(lq  + allocl * wtnorm + etalq)
    vc     <- exp(lvc + allov  * wtnorm + etalvc)
    vp     <- exp(lvp + allov  * wtnorm + etalvp)

    Cc <- linCmt()

    f(depot) <- fdepot
    Cc ~ prop(propSd)
  })
}
