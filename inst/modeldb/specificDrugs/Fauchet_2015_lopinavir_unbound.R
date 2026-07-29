Fauchet_2015_lopinavir_unbound <- function() {
  description <- "One-compartment first-order-absorption population PK model for lopinavir in HIV-infected pregnant and nonpregnant women parameterised on the unbound fraction, with total LPV reconstructed from a linear HSA binding term plus a saturable single-site AAG binding term (Fauchet 2015 unbound submodel)."
  reference <- "Fauchet F, Treluyer JM, Illamola SM, Pressiat C, Lui G, Valade E, Mandelbrot L, Lechedanec J, Delmas S, Blanche S, Warszawski J, Urien S, Tubiana R, Hirt D, for the ANRS 135 PRIMEVA Study Group. Population approach to analyze the pharmacokinetics of free and total lopinavir in HIV-infected pregnant women and consequences for dose adjustment. Antimicrob Agents Chemother. 2015;59(9):5727-5735. doi:10.1128/AAC.00863-15"
  vignette <- "Fauchet_2015_lopinavir"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    ALB = list(
      description        = "Human serum albumin plasma concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Drives the linear-binding contribution to total LPV in the protein-binding equation (Fauchet 2015 Table 2 / Appendix). The paper labels the column 'HSA' (human serum albumin); same analyte as the canonical ALB / serum albumin (now an ALB source alias). Time-fixed within an LPV sampling visit; median 32 g/L (range 16-48 g/L) in pregnant women (paper Results). Converted to umol/L inside model() using molecular weight 66500 g/mol.",
      source_name        = "HSA"
    ),
    AAG = list(
      description        = "Alpha-1 acid glycoprotein plasma concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Drives the saturable single-site AAG contribution to total LPV in the protein-binding equation (Fauchet 2015 Table 2 / Appendix); N_AAG fixed to 1. Time-fixed within an LPV sampling visit; median 0.55 g/L (range 0.20-1.20 g/L) in pregnant women (paper Results). Converted to umol/L inside model() using molecular weight 40000 g/mol.",
      source_name        = "AAG"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested in the covariate screen but did not significantly improve the unbound model (Fauchet 2015 Results, 'Unbound model' subsection: 'covariates were tested, but body weight, age, and gestational age had no significant effect').",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Maternal age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested in the covariate screen but not retained (Fauchet 2015 Results, 'Unbound model' subsection).",
      source_name        = "AGE"
    ),
    GA = list(
      description        = "Gestational age at sampling",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested in the covariate screen (paper Methods equation 3) but did not significantly improve the unbound model (Fauchet 2015 Results, 'Unbound model' subsection). Time-varying within a pregnancy when sampled at multiple gestational ages; this canonical (gestational age at birth) is reused for gestational-age-at-sampling here because the column has the same units and orientation; per-model notes carry the distinction.",
      source_name        = "GA"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 208L,
    n_studies      = 2L,
    n_observations = 400L,
    age_range      = "18.3-44.3 years",
    age_median     = "31.2 years",
    weight_range   = "45-122 kg",
    weight_median  = "76 kg",
    sex_female_pct = 100,
    race_ethnicity = "not reported in the source paper",
    disease_state  = "HIV-1 infection; pregnant, nonpregnant, and in-labour women treated for HIV infection or for prevention of mother-to-child transmission",
    ga_range       = "9-41 weeks (sampled gestational age across pregnant cohorts)",
    dose_range     = "LPV/r 400/100 mg twice daily oral (one subject 600 mg LPV BID)",
    regions        = "France",
    notes          = "Pooled cohort across the ANRS 135 PRIMEVA randomised trial (n=103: 69 LPV/r monotherapy, 34 LPV/r + zidovudine/lamivudine triple therapy) and the Hospital Cochin (Paris) therapeutic drug monitoring cohort (n=105: 81 pregnant, 24 nonpregnant). Total + unbound LPV pairs were measured in 79 samples of the 400 maternal samples used here; unbound model parameters were estimated against the full total-LPV set with the unbound observations linking total to free via the HSA + AAG binding equation. Cohort demographics from Table 1."
  )

  ini({
    # Structural parameters (apparent values parameterised on the unbound fraction).
    # Source: Fauchet 2015 Table 2 ("Population pharmacokinetic parameters of lopinavir
    # from the unbound fraction model").
    lka <- log(0.408); label("Absorption rate constant (1/h)")                                  # Table 2 row 'Ka' = 0.408 /h
    lcl <- log(316);   label("Apparent unbound clearance (CL_unbound/F, L/h)")                  # Table 2 row 'CL' = 316 L/h
    lvc <- log(2220);  label("Apparent unbound volume of distribution (V_unbound/F, L)")        # Table 2 row 'V' = 2,220 L

    # Protein-binding parameters (paper-mechanistic). HSA term is linear in the protein
    # concentration; AAG term is single-site saturable. Source: Fauchet 2015 Table 2 and
    # Appendix equation 2.
    lkhsa <- log(0.036);  label("Linear HSA binding constant K_HSA (L/umol; product N_HSA*K_HSA)")  # Table 2 row 'K_HSA' = 0.036 L/umol
    lkaag <- log(0.159);  label("Saturable AAG dissociation constant K_AAG (umol/L)")               # Table 2 row 'N_AAG K_AAG' with K_AAG = 0.159 umol/L
    naag  <- fixed(1);    label("Number of LPV binding sites on AAG (unitless, fixed)")             # Table 2 footnote c 'Fixed value'; Results paragraph: N_AAG was fixed to 1 in the final model

    # IIV on apparent unbound clearance only. Fauchet 2015 reports omega as the square
    # root of the between-subject variance (Table 2 footnote d), so the variance entered
    # here is omega^2.
    etalcl ~ 0.05198  # Table 2 row 'omega_Cl' = 0.228 (omega = SD on log-CL); variance = 0.228^2 = 0.05198

    # Residual error -- two separate proportional errors (Table 2 footnote e). The bare
    # propSd applies to the total-LPV output Cc; propSd_Cunbound applies to the unbound
    # output Cunbound.
    propSd          <- 0.450; label("Proportional residual error on total LPV (fraction)")       # Table 2 row 'sigma_total' = 0.450
    propSd_Cunbound <- 0.282; label("Proportional residual error on unbound LPV (fraction)")     # Table 2 row 'sigma_unbound' = 0.282
  })

  model({
    # Molecular-weight constants used to convert g/L protein and mg/L LPV to umol/L for
    # the binding equation. LPV: 628.8 g/mol; serum albumin (ALB / HSA): ~66500 g/mol
    # (66.5 kDa); AAG: ~40000 g/mol (~40 kDa). The published K_HSA and K_AAG estimates
    # carry molar units so protein concentrations and the unbound LPV concentration
    # must be in umol/L. Note the factor difference: protein concentrations are in g/L
    # (factor 1e6/MW) while LPV concentration is in mg/L (factor 1000/MW).
    MW_LPV <- 628.8
    MW_ALB <- 66500
    MW_AAG <- 40000

    # Individual PK parameters (apparent unbound values; F absorbed into CL_unbound/F
    # and V_unbound/F).
    ka  <- exp(lka)
    cl  <- exp(lcl + etalcl)
    vc  <- exp(lvc)
    kel <- cl / vc

    # Paper-mechanistic protein-binding constants.
    khsa <- exp(lkhsa)
    kaag <- exp(lkaag)

    # Convert covariate concentrations (g/L) to micromolar.
    ALB_uM <- ALB * 1e6 / MW_ALB
    AAG_uM <- AAG * 1e6 / MW_AAG

    # ODE system: standard one-compartment with first-order absorption. The central
    # state holds the amount of LPV apportioned to the unbound-equivalent space, so
    # Cunbound = central / vc directly gives unbound concentration (mg/L) and the
    # "apparent" CL/V already absorb F (Fauchet 2015 Appendix equation 1).
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Unbound concentration (mg/L = ug/mL) and its molar value.
    Cunbound    <- central / vc
    Cunbound_uM <- Cunbound * 1000 / MW_LPV

    # Protein-binding equation (Fauchet 2015 Appendix equation 2). HSA term is linear;
    # AAG term is single-site saturable (N_AAG fixed to 1).
    bound_alb_uM <- khsa * ALB_uM * Cunbound_uM
    bound_aag_uM <- naag * AAG_uM * Cunbound_uM / (kaag + Cunbound_uM)
    Ctotal_uM    <- Cunbound_uM + bound_alb_uM + bound_aag_uM

    # Total LPV output (mg/L = ug/mL).
    Cc <- Ctotal_uM * MW_LPV / 1000

    # Two-output residual error (proportional on each scale).
    Cc       ~ prop(propSd)
    Cunbound ~ prop(propSd_Cunbound)
  })
}
