PK_2cmt_mAb_Davda_2014 <- function() {
  description <- "Two compartment PK model with linear clearance for average monoclonal antibodies (Davda 2014)"
  reference <- "Davda JP, Dodds MG, Gibbs MA, Wisdom W, Gibbs JP. A model-based meta-analysis of monoclonal antibody pharmacokinetics to guide optimal first-in-human study design. MAbs. 2014;6(4):1094-1102. doi:10.4161/mabs.29095"
  vignette <- "PK_2cmt_mAb_Davda_2014"
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
    n_subjects     = 171,
    n_studies      = 4,
    n_observations = 2716,
    age_range      = NA,
    weight_range   = NA,
    sex_female_pct = NA,
    race_ethnicity = NA,
    disease_state  = "Healthy volunteers enrolled in first-in-human studies of four therapeutic monoclonal antibodies (mAb a/b/c: human IgG2; mAb d: IgG1). All four antibodies target soluble ligands.",
    dose_range     = "IV: 1-700 mg; SC: 2.1-700 mg",
    regions        = NA,
    notes          = "Model-based meta-analysis pooling individual-level data from 4 FIH studies (171 subjects, 2716 serum concentrations: 1153 IV, 1563 SC). Per-study demographic breakdowns (age, sex, weight, race) are not tabulated in the published paper; reference weight for allometric scaling is 70 kg."
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
