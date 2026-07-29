Bouazza_2011_lamivudine <- function() {
  description <- "Two-compartment oral popPK model for lamivudine in HIV-infected children from neonates to adolescents (Bouazza 2011)"
  reference <- "Bouazza N, Hirt D, Blanche S, Frange P, Rey E, Treluyer JM, Urien S. Developmental pharmacokinetics of lamivudine in 580 pediatric patients ranging from neonates to adolescents. Antimicrobial Agents and Chemotherapy. 2011;55(8):3498-3504. doi:10.1128/AAC.01622-10"
  vignette <- "Bouazza_2011_lamivudine"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used for allometric scaling on CL, Q, Vc, Vp with reference weight 70 kg. Paper cohort range 1-84 kg (median 23 kg).",
      source_name        = "BW"
    ),
    PAGE = list(
      description        = "Postmenstrual age (PMA = postnatal age + gestational age)",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Drives the sigmoidal CL maturation function. The paper parameterises PMA in weeks (PMA50 = 59 weeks); inside model() PAGE is multiplied by 4.345 weeks/month to recover the published equation. When gestational age is unknown the paper imputes 40 weeks (term birth) before computing PMA.",
      source_name        = "PMA"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 580,
    n_studies      = 1,
    n_observations = 2106,
    age_range      = "2 days to 18 years (median 7.41 years)",
    weight_range   = "1 to 84 kg (median 23 kg)",
    sex_female_pct = NA_real_,
    race_ethnicity = NA,
    disease_state  = "HIV-1 infection or prevention of mother-to-child transmission",
    dose_range     = "Median 7.5 (SD 3.2) mg/kg/day; tablet or oral solution; BID or OAD regimens",
    regions        = "France (Paris hospitals; retrospective therapeutic drug monitoring)",
    notes          = "Bouazza 2011 Materials and Methods 'Patients and treatment'. Race/ethnicity and sex distribution not reported in the source. The galenic form (tablet vs oral solution) was tested as a categorical covariate on bioavailability and had no significant effect."
  )

  covariatesDataExcluded <- list(
    FORMULATION = list(
      description = "Galenic form (tablet vs oral solution)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested as CA in CL = theta_CL * beta_CA but not retained: galenic form had no significant effect on bioavailability (Bouazza 2011 Results)."
    )
  )

  ini({
    # Structural parameters at the standard 70 kg reference (Bouazza 2011 Table 1)
    lka  <- log(0.432); label("Absorption rate constant (Ka, 1/h)")
    lcl  <- log(31);    label("Apparent clearance at 70 kg / fully matured (CL/F, L/h)")
    lvc  <- log(76.4);  label("Apparent central volume at 70 kg (Vc/F, L)")
    lq   <- log(5.83);  label("Apparent intercompartmental clearance at 70 kg (Q/F, L/h)")
    lvp  <- log(129);   label("Apparent peripheral volume at 70 kg (Vp/F, L)")

    # Allometric exponents (theory-based, not estimated; Bouazza 2011 Materials and Methods).
    # The paper applies the canonical 0.75 / 1 values without reporting RSE -> fixed.
    e_wt_cl_q  <- fixed(0.75); label("Allometric WT exponent shared across CL/F and Q/F (unitless, fixed)")
    e_wt_vc_vp <- fixed(1);    label("Allometric WT exponent shared across Vc/F and Vp/F (unitless, fixed)")

    # Maturation parameters for CL: F(PMA) = PMA^gamma / (PMA50^gamma + PMA^gamma)
    # Bouazza 2011 Table 1. PMA50 is in weeks; the paper's published symbols are PMA50 and gamma.
    pma50_cl <- 59;   label("PMA at which CL reaches 50% of mature value (weeks)")
    gamma_cl <- 3.02; label("Hill exponent for PMA effect on CL (unitless)")

    # IIV (exponential model -> omega is the SD of eta on log scale; variance is omega^2).
    # Bouazza 2011 Table 1 reports omega_CL = 0.32 and omega_Vc = 0.77.
    etalcl ~ 0.1024  # 0.32^2; ~32% CV
    etalvc ~ 0.5929  # 0.77^2; ~90% CV (high; eta-shrinkage 0.39 in the original fit)

    # Residual error (proportional only; the additive part was non-significant). Table 1: sigma = 0.5.
    propSd <- 0.5; label("Proportional residual error (fraction)")
  })
  model({
    # Convert canonical PAGE (months) to PMA in weeks so the published PMA50 (weeks) applies directly.
    pma_wk <- PAGE * 4.345

    # Sigmoidal maturation of CL (Bouazza 2011 final covariate model). Approaches 1 as PMA -> infinity.
    maturation_cl <- pma_wk^gamma_cl / (pma50_cl^gamma_cl + pma_wk^gamma_cl)

    # PK parameters with allometric weight scaling (reference 70 kg) and CL maturation.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q  * maturation_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    q  <- exp(lq)           * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp)          * (WT / 70)^e_wt_vc_vp

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Concentration: dose in mg, volume in L -> mg/L (= ug/mL).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
