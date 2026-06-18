Wahlby_2004_pefloxacin <- function() {
  description <- "One-compartment population PK model for intravenous 1-hour pefloxacin infusions in 74 critically ill adults, demonstrating Wahlby 2004's extended covariate-model formulation. Final-model clearance carries time-varying CRCL (with inter-individual variability in the CRCL effect coefficient, Eq 3), per-subject baseline total bilirubin BIL_BASE (replaces BIL), age, centre indicator, and per-subject baseline weight WT_BASE (replaces WT, with a saturating 'up to median weight' qualifier per Methods). Central volume retains the upstream Karlsson 1993 (ref [10]) WT/CRCL/BIL effects unchanged. Underlying structural PK comes from Karlsson MO, Sheiner LB. The importance of modeling interoccasion variability in population pharmacokinetic analyses. J Pharmacokin Biopharm 1993;21(6):735-750 (not on disk in this worktree)."
  reference <- paste(
    "Wahlby U, Thomson AH, Milligan PA, Karlsson MO.",
    "Models for time-varying covariates in population pharmacokinetic-pharmacodynamic analysis.",
    "Br J Clin Pharmacol 2004;58(4):367-377.",
    "doi:10.1111/j.1365-2125.2004.02170.x.",
    "Structural PK model and the reference V parameterisation carried from Karlsson MO, Sheiner LB.",
    "The importance of modeling interoccasion variability in population pharmacokinetic analyses.",
    "J Pharmacokin Biopharm 1993;21(6):735-750; the present entry encodes Wahlby 2004 Table 6 Final-Model column.",
    sep = " "
  )
  vignette <- "Wahlby_2004_time_varying_covariates"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance, time-varying within an individual (raw mL/min, NOT BSA-normalized).",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Table 2: mean 113.2, median 103.5, range 0.43-312.0 mL/min. Centered at 100 mL/min per Eq 6.",
      source_name        = "CLC"
    ),
    AGE = list(
      description        = "Subject age (years).",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline. Centered at 45 years per Eq 6.",
      source_name        = "AGE"
    ),
    CEN = list(
      description        = "Study-centre indicator (paper-specific 0/1; the source Wahlby 2004 references CEN to a single coefficient theta_CEN-CL without identifying which centre is 1 vs 0).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0",
      notes              = "Pefloxacin-specific (specific-scope) study-centre indicator carried over from Karlsson 1993. The Wahlby 2004 paper does not specify which centre corresponds to CEN = 1 versus CEN = 0; users supplying a virtual cohort should treat CEN as a sensitivity covariate.",
      source_name        = "CEN"
    ),
    WT = list(
      description        = "Body weight, time-varying within an individual.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying; enters the V model with effect e_wt_v centered at 65 kg (Eq 6 convention). Table 2: mean 69, median 67, range 42.7-125 kg.",
      source_name        = "WT"
    ),
    WT_BASE = list(
      description        = "Per-subject baseline body weight, time-fixed (BWT in the source paper).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-subject constant. Used in CL with a 'saturating up to median WT' qualifier (Methods): the effective WT in the CL model is min(WT_BASE, 65). Above 65 kg the effect plateaus. Table 2 BWT row: mean 68.8, median 67 kg.",
      source_name        = "BWT"
    ),
    TBILI = list(
      description        = "Total serum bilirubin, time-varying within an individual.",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Used in the V model (carried from Karlsson 1993). Centered at 25 umol/L per Eq 6.",
      source_name        = "BIL"
    ),
    TBILI_BASE = list(
      description        = "Per-subject baseline total serum bilirubin, time-fixed (BBIL in the source paper).",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-subject constant. Replaces BIL in the final-model CL equation because BBIL was the better predictor (RSE 11% versus BIL RSE 18%). Centered at 25 umol/L (the reference value used in Eq 6 for BIL).",
      source_name        = "BBIL"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 74L,
    n_studies      = 1L,
    age_range      = NA_character_,
    weight_range   = "Table 2 WT range 42.7-125 kg",
    sex_female_pct = NA_real_,
    disease_state  = "Critically ill adult patients receiving intravenous 1-hour pefloxacin infusions over 1-28 days (median 6 days) across 1-4 treatment periods (separated by 2.5-14 days). Wahlby 2004 re-analyses a previously-reported cohort (Karlsson MO, Sheiner LB. J Pharmacokin Biopharm 1993;21(6):735-750).",
    dose_range     = NA_character_,
    n_observations = 337L,
    follow_up      = "1-28 days per patient (median 6 days), with 1-4 treatment periods.",
    regions        = NA_character_,
    notes          = "The paper reports inter-occasion variability on CL (pi_CL = 0.32) and inter-individual variability in the CRCL effect coefficient (omega_CLC = 0.61) but no IIV on CL itself; see the vignette Errata for the IOV-vs-IIV encoding choice in this library entry."
  )

  ini({
    # Structural typical-value parameters (Wahlby 2004 Table 6 Final-Model column)
    lcl <- log(3.74); label("Population typical CL (L/h)")        # Table 6: theta_CL = 3.74 (RSE 7.6%)
    lvc <- log(61);   label("Population typical V (L)")           # Table 6: theta_V = 61 carried from Reference-Model column (Final column shows "-", interpreted as unchanged; underlying Karlsson 1993 V parameterisation)

    # Covariate effects on CL (Wahlby 2004 Eq 6 form for the Reference Model, applied to the
    # Final-Model column substitutions: BBIL replaces BIL; BWT replaces WT; SBP dropped)
    e_crcl_cl <- 0.0029;  label("CRCL effect on CL: coefficient in exp[e_crcl_cl * (CRCL - 100)]")             # Table 6: theta_CLC-CL = 0.0029 (RSE 17%)
    e_bbil_cl <- -0.0068; label("BIL_BASE effect on CL: coefficient in exp[e_bbil_cl * (BIL_BASE - 25)]")      # Table 6: theta_BBIL-CL = -0.0068 (RSE 11%)
    e_age_cl  <- -0.0086; label("AGE effect on CL: coefficient in exp[e_age_cl * (AGE - 45)]")                 # Table 6: theta_AGE-CL = -0.0086 (RSE 21%)
    e_cen_cl  <- 0.19;    label("Centre indicator effect on CL: coefficient in exp[e_cen_cl * CEN]")            # Table 6: theta_CEN-CL = 0.19 (RSE 41%)
    e_bwt_cl  <- 0.039;   label("WT_BASE effect on CL: coefficient in exp[e_bwt_cl * (min(WT_BASE,65) - 65)]") # Table 6: theta_BWT-CL = 0.039 (RSE 24%); the 'up to median WT' qualifier from Methods plateaus the effect above 65 kg

    # Covariate effects on V (Karlsson 1993 reference form carried unchanged into Wahlby's final model)
    e_wt_vc   <- 0.014;    label("WT effect on Vc: coefficient in exp[e_wt_vc * (WT - 65)]")     # Table 6: theta_WT-V = 0.014 (RSE 13%)
    e_crcl_vc <- 0.0015;   label("CRCL effect on Vc: coefficient in exp[e_crcl_vc * (CRCL - 100)]") # Table 6: theta_CLC-V = 0.0015 (RSE 23%)
    e_bil_vc  <- -0.0025;  label("BIL effect on Vc: coefficient in exp[e_bil_vc * (TBILI - 25)]")  # Table 6: theta_BIL-V = -0.0025 (RSE 28%)

    # Inter-individual variability in the CRCL effect (Wahlby 2004 Eq 3 demonstrated)
    # Per the paper this is "omega in CL ~ CLC relationship", entering multiplicatively
    # in the CRCL covariate term.
    etae_crcl_cl ~ 0.61^2  # Table 6: omega_CLC = 0.61 (RSE 36% relative to variance term)

    # IIV on CL itself (carried from the IOV term in the source paper - the paper reports
    # pi_CL = 0.32 (IOV) but no omega_CL (IIV); this library entry collapses IOV into IIV
    # to preserve the variance magnitude in single-occasion simulations. Documented in
    # vignette Errata.)
    etalcl ~ 0.32^2  # Table 6: pi_CL = 0.32 (RSE 16% relative to variance term); encoded here as IIV

    # Residual error (Wahlby 2004 Table 6 Final-Model column; proportional only)
    propSd <- 0.16; label("Proportional residual error (fraction)")  # Table 6: sigma = 0.16 (RSE 17%)
  })

  model({
    # Individual CRCL effect coefficient (Wahlby 2004 Eq 3 demonstrated on the
    # CRCL relationship)
    e_crcl_cl_i <- e_crcl_cl * exp(etae_crcl_cl)

    # Saturating WT_BASE effect: the source paper's Eq 6 qualifier "weight (WT) up to
    # median WT" plateaus the WT effect at 65 kg. In the Final-Model column this
    # applies to WT_BASE rather than time-varying WT.
    wt_base_eff <- min(WT_BASE, 65) - 65

    # Individual parameters
    cl <- exp(lcl + etalcl) * exp(
      e_crcl_cl_i * (CRCL - 100) +
      e_bbil_cl   * (TBILI_BASE - 25) +
      e_age_cl    * (AGE - 45) +
      e_cen_cl    * CEN +
      e_bwt_cl    * wt_base_eff
    )

    # Central volume (Karlsson 1993 reference V parameterisation carried unchanged into
    # Wahlby 2004 Final-Model column)
    vc <- exp(lvc) * exp(
      e_wt_vc   * (WT - 65) +
      e_crcl_vc * (CRCL - 100) +
      e_bil_vc  * (TBILI - 25)
    )

    kel <- cl / vc

    # One-compartment IV pefloxacin (administered as 1-hour infusion in the source);
    # the library model does not hard-code an infusion duration so users specify rate /
    # dur per dose in the event table.
    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
