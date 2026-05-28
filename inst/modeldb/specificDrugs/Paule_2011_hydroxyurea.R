Paule_2011_hydroxyurea <- function() {
  description <- "Two-compartment population PK + indirect-response PD models for hydroxyurea (HU) in adults with sickle cell anemia (Paule 2011): bicompartmental oral PK with first-order absorption and elimination, allometric scaling on CL/F and Vc/F; turnover models for HbF percentage and mean corpuscular volume (MCV) where HU inhibits the elimination rate of each PD response."
  reference <- "Paule I, Sassi H, Habibi A, Pham KPD, Bachir D, Galacteros F, Girard P, Hulin A, Tod M. Population pharmacokinetics and pharmacodynamics of hydroxyurea in sickle cell anemia patients, a basis for optimizing the dosing regimen. Orphanet J Rare Dis. 2011;6:30. doi:10.1186/1750-1172-6-30 (PMID 21595938)."
  vignette <- "Paule_2011_hydroxyurea"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L", HbF = "%", MCV = "fL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling reference 70 kg; exponent 0.75 on CL/F and 1.00 on Vc/F (Paule 2011 Methods, Population pharmacokinetic model section).",
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 97L,
    n_studies      = 2L,
    age_range      = "18-54 years (sparse PK/PD cohort) and 24-52 years (rich PK cohort)",
    age_median     = "30 years (sparse), 32 years (rich)",
    weight_range   = "45-163 kg (sparse) and 42-71 kg (rich)",
    weight_median  = "60 kg (sparse), 63 kg (rich)",
    sex_female_pct = 70.1,
    disease_state  = "Adult homozygous (HbSS) sickle cell anemia",
    dose_range     = "Oral hydroxyurea 500-2000 mg once daily; initial dose 20 mg/kg (10-15 mg/kg in renal insufficiency); titrated to maintain neutrophils > 3e9/L, final dose typically <= 30 mg/kg",
    regions        = "France (single center: AP-HP, GH H. Mondor, Universite Paris Est-Creteil)",
    notes          = "Two pooled datasets: a 30-month sparse-sampling observational PK/PD study (n=81, up to 9 samples/patient, baseline + every 15 days or month) and a 24-hour bioequivalence rich-sampling PK study (n=16, 10 samples/patient). Females 57/81 (sparse) and 11/16 (rich); combined 68/97 (70.1%). Baseline demographics from Paule 2011 Table 1. HbF/MCV summary statistics from Table 3."
  )

  ini({
    # ---------------------------------------------------------------------
    # PK structural parameters (reference 70 kg adult; Paule 2011 Table 2).
    # The source paper reports PK rate constants in 1/h and CL in L/h; the
    # PD parameters are in 1/day. To run PK and PD on a single time scale
    # we convert PK to per-day units (multiply rate constants by 24 h/day,
    # and CL/F by 24 h/day). The fundamental disposition is unchanged.
    # ---------------------------------------------------------------------
    lka  <- log(3.02 * 24);  label("Population absorption rate constant (Ka, 1/day)")             # Paule 2011 Table 2 (theta_ka 3.02 1/h)
    lcl  <- log(11.6 * 24);  label("Apparent oral clearance at 70 kg (CL/F, L/day)")              # Paule 2011 Table 2 (11.6 L/h)
    lvc  <- log(45.3);       label("Apparent central volume of distribution at 70 kg (Vc/F, L)")  # Paule 2011 Table 2
    lkcp <- log(0.027 * 24); label("Central-to-peripheral first-order transfer rate (kcp, 1/day)")  # Paule 2011 Table 2 (0.027 1/h)
    lkpc <- fixed(log(0.004 * 24)); label("Peripheral-to-central first-order transfer rate (kpc, 1/day) -- fixed")  # Paule 2011 Table 2 (kpc fixed at 0.004 1/h)

    # Allometric exponents (Paule 2011 Methods; held fixed at the canonical
    # values reported in the paper)
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL (unitless)")  # Paule 2011 Methods (allometric scaling exponent 0.75 on CL)
    e_wt_vc <- fixed(1.00); label("Allometric exponent of body weight on Vc (unitless)")  # Paule 2011 Methods (weight-proportional Vc)

    # ---------------------------------------------------------------------
    # HbF percentage PD parameters (Paule 2011 Table 4)
    # Turnover model: dHbF/dt = Kin_HbF - Kout_HbF * (1 - I_HbF) * HbF
    # ---------------------------------------------------------------------
    lkin_hbf      <- log(0.071); label("Typical zero-order production rate of HbF% (Kin_HbF, %/day)")  # Paule 2011 Table 4
    lkout_hbf     <- log(0.013); label("Typical first-order elimination rate of HbF% (Kout_HbF, 1/day)")  # Paule 2011 Table 4
    logitimax_hbf <- 0.276;      label("Logit-transformed Imax for HbF (LImax, unitless); back-transformed Imax = 0.569")  # Paule 2011 Table 4 (LImax)

    # ---------------------------------------------------------------------
    # MCV PD parameters (Paule 2011 Table 5)
    # Turnover model: dMCV/dt = Kin_MCV - Kout_MCV * (1 - b * C^gamma) * MCV
    # ---------------------------------------------------------------------
    lkin_mcv  <- log(3.71);  label("Typical zero-order production rate of MCV (Kin_MCV, fL/day)")  # Paule 2011 Table 5
    lkout_mcv <- log(0.042); label("Typical first-order elimination rate of MCV (Kout_MCV, 1/day)")  # Paule 2011 Table 5
    lb_mcv    <- log(0.099); label("Power-function scale constant b for MCV inhibition (units of (L/mg)^gamma)")  # Paule 2011 Table 5 (b)
    gamma_mcv <- 0.19;       label("Power-function exponent gamma for MCV inhibition (unitless)")  # Paule 2011 Table 5 (gamma)

    # ---------------------------------------------------------------------
    # Inter-individual variability (Paule 2011 Tables 2, 4, 5)
    # All omegas reported as SDs on the log scale (omega^2 = SD^2).
    # ---------------------------------------------------------------------
    # PK IIV (Paule 2011 Table 2): block on (Vc, CL, kcp); independent ka.
    # SDs: Vc 0.34, CL 0.29, kcp 0.57. Correlations: Vc-CL 0.71,
    # Vc-kcp -0.26, CL-kcp 0.37. Variances and covariances:
    #   var(Vc)  = 0.34^2 = 0.1156
    #   var(CL)  = 0.29^2 = 0.0841
    #   var(kcp) = 0.57^2 = 0.3249
    #   cov(Vc, CL)  = 0.71 * 0.34 * 0.29 = 0.0700
    #   cov(Vc, kcp) = -0.26 * 0.34 * 0.57 = -0.0504
    #   cov(CL, kcp) = 0.37 * 0.29 * 0.57 = 0.0612
    etalvc + etalcl + etalkcp ~ c(0.1156,
                                  0.0700,  0.0841,
                                 -0.0504,  0.0612, 0.3249)  # Paule 2011 Table 2
    # ka IIV (Paule 2011 Table 2): SD 1.34, no correlation listed with the
    # Vc/CL/kcp block. Paper-reported formula k_a = theta_ka * exp(eta_ka + alpha)
    # is encoded here as standard log-normal IIV on theta_ka = 3.02; the
    # reported population-marginal value k_a = 3.29 differs from theta_ka
    # but is unexplained in the trimmed source and likely reflects a small
    # population-mean shift (alpha) that the table does not quantify.
    etalka ~ 1.7956  # Paule 2011 Table 2 (SD 1.34 -> 1.34^2 = 1.7956)

    # HbF IIV (Paule 2011 Table 4): block on (Kout, LImax); independent Kin.
    # SDs: Kin 0.585, Kout 0.486, LImax 1.44. Correlation Kout-LImax = 0.892.
    etalkin_hbf ~ 0.3422  # Paule 2011 Table 4 (SD 0.585 -> 0.585^2 = 0.3422)
    etalkout_hbf + etalogitimax_hbf ~ c(0.2362,
                                        0.6242, 2.0736)  # Paule 2011 Table 4 (var, cov, var)

    # MCV IIV (Paule 2011 Table 5): full 3x3 block on (Kin, Kout, b).
    # SDs: Kin 0.191, Kout 0.186, b 0.457. Correlations: Kin-Kout 0.87,
    # Kin-b -0.98, Kout-b -0.95. Covariances:
    #   cov(Kin, Kout) = 0.87  * 0.191 * 0.186 = 0.0309
    #   cov(Kin, b)    = -0.98 * 0.191 * 0.457 = -0.0855
    #   cov(Kout, b)   = -0.95 * 0.186 * 0.457 = -0.0808
    etalkin_mcv + etalkout_mcv + etalb_mcv ~ c(0.0365,
                                                0.0309,  0.0346,
                                               -0.0855, -0.0808, 0.2088)  # Paule 2011 Table 5

    # ---------------------------------------------------------------------
    # Residual error
    # ---------------------------------------------------------------------
    # PK: mixed additive + proportional, reported separately for the sparse
    # and the rich datasets. The library default uses the sparse-dataset
    # values, which are representative of multi-occasion observational PK
    # (the typical user scenario). The dense-sampling residuals (additive
    # 0.319 mg/L, proportional 0.12) are documented in the vignette.
    addSd  <- 0.353; label("Additive residual SD on HU concentration (mg/L; sparse-dataset value)")     # Paule 2011 Table 2 (sparse-data residual error)
    propSd <- 0.435; label("Proportional residual SD on HU concentration (fraction; sparse-dataset)")    # Paule 2011 Table 2 (sparse-data residual error)

    # PD residual: proportional on each output.
    propSd_HBF <- 0.142; label("Proportional residual SD on HbF% (fraction)")  # Paule 2011 Table 4
    propSd_MCV <- 0.036; label("Proportional residual SD on MCV (fraction)")    # Paule 2011 Table 5
  })

  model({
    # -------------------------------------------------------------------
    # Individual PK parameters with allometric scaling on body weight
    # (reference 70 kg, Paule 2011 Methods). Imax interpretation of theta_ka
    # follows the paper's standard log-normal IIV: ka_i = theta_ka * exp(eta_ka).
    # -------------------------------------------------------------------
    ka  <- exp(lka  + etalka)
    cl  <- exp(lcl  + etalcl)  * (WT / 70)^e_wt_cl
    vc  <- exp(lvc  + etalvc)  * (WT / 70)^e_wt_vc
    kcp <- exp(lkcp + etalkcp)
    kpc <- exp(lkpc)
    kel <- cl / vc

    # -------------------------------------------------------------------
    # PK ODE system: 2-compartment with first-order absorption / elimination
    # -------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - kcp * central + kpc * peripheral1
    d/dt(peripheral1) <-  kcp * central - kpc * peripheral1

    Cc <- central / vc

    # -------------------------------------------------------------------
    # HbF percentage turnover (Paule 2011 Table 4)
    # Imax obtained by back-transforming the logit: Imax = 1/(1+exp(-LImax))
    # Drug-activation gating: Paule 2011 found Imax saturated across the
    # studied dose range (0.5 - 9 mg/L) and could not identify a
    # concentration-response within that range. To run forward simulation
    # we encode this as a smooth saturation Cc/(Cc + c50eps) with
    # c50eps = 0.005 mg/L (well below the assay LOQ of 0.532 mg/L); the
    # inhibition therefore equals Imax for any detectable drug
    # concentration and equals zero when no drug is on board. c50eps is
    # an engineering constant, not a fitted parameter; see vignette
    # Assumptions and deviations.
    # -------------------------------------------------------------------
    kin_hbf       <- exp(lkin_hbf  + etalkin_hbf)
    kout_hbf      <- exp(lkout_hbf + etalkout_hbf)
    logit_imax    <- logitimax_hbf + etalogitimax_hbf
    imax_hbf      <- 1 / (1 + exp(-logit_imax))
    inhib_hbf     <- imax_hbf * Cc / (Cc + 0.005)

    # Off-drug steady state: effect1(SS) = kin_hbf / kout_hbf (matches the
    # 6.3% baseline median from Paule 2011 Table 3 for the typical patient
    # at kin=0.071, kout=0.013 -> 5.46%). On-drug steady state:
    # kin_hbf / (kout_hbf * (1 - Imax)).
    d/dt(effect1) <- kin_hbf - kout_hbf * (1 - inhib_hbf) * effect1
    effect1(0)    <- kin_hbf / kout_hbf

    HBF <- effect1

    # -------------------------------------------------------------------
    # MCV turnover (Paule 2011 Table 5)
    # Inhibition is a power function of the average HU concentration. The
    # source paper used the individual posthoc average concentration Cavg
    # (Dose / (CL*tau)) as the PD input; for forward simulation we use the
    # instantaneous central concentration Cc as a proxy (documented in the
    # vignette as a deviation). With gamma = 0.19 the function is very
    # flat (b*Cc^gamma = 0.11 at 2 mg/L and 0.15 at 9 mg/L), so the
    # substitution has a small effect on long-term dynamics. A tiny
    # positive offset (1e-9) avoids the 0^positive singularity at Cc = 0.
    # -------------------------------------------------------------------
    kin_mcv  <- exp(lkin_mcv  + etalkin_mcv)
    kout_mcv <- exp(lkout_mcv + etalkout_mcv)
    b_mcv    <- exp(lb_mcv    + etalb_mcv)

    inhib_mcv <- b_mcv * (Cc + 1e-9)^gamma_mcv - b_mcv * (1e-9)^gamma_mcv

    d/dt(effect2) <- kin_mcv - kout_mcv * (1 - inhib_mcv) * effect2
    effect2(0)    <- kin_mcv / kout_mcv

    MCV <- effect2

    # -------------------------------------------------------------------
    # Observation and error
    # -------------------------------------------------------------------
    Cc  ~ add(addSd) + prop(propSd)
    HBF ~ prop(propSd_HBF)
    MCV ~ prop(propSd_MCV)
  })
}
