Diao_2016_daclizumab_treg <- function() {
  description <- "Sigmoidal Emax PK/PD model of regulatory T cell (Treg) reduction following subcutaneous daclizumab high-yield process (HYP) in adults with relapsing-remitting multiple sclerosis (Diao 2016). The PD output is the percentage of Treg (CD4+ CD127low/- Foxp3+) among all CD4+ T cells; daclizumab HYP serum concentration drives a maximum 60% reduction via a sigmoidal Emax function. The PK backbone is the two-compartment, first-order SC absorption + lag model from Othman 2014 (file inst/modeldb/specificDrugs/Othman_2014_daclizumab.R), copied verbatim with weight-based allometric scaling."
  reference <- "Diao L, Hang Y, Othman AA, Nestorov I, Tran JQ, Mehta D, Amaravadi L. Population PK/PD analyses of CD25 occupancy, CD56 bright NK cell expansion and regulatory T cell reduction by daclizumab HYP in subjects with multiple sclerosis. Br J Clin Pharmacol. 2016;82(5):1333-1342. doi:10.1111/bcp.13051 (PMID 27333593). PK backbone: Othman AA, Tran JQ, Tang MT, Dutta S. Population Pharmacokinetics of Daclizumab High-Yield Process in Healthy Volunteers. Clin Pharmacokinet. 2014;53(10):907-918. doi:10.1007/s40262-014-0159-9."
  vignette <- "Diao_2016_daclizumab_treg"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL", response = "% of CD4+ T cells")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used for allometric scaling of the inherited Othman 2014 PK parameters (CL, Q, Vc, Vp) with reference 70 kg; exponents 0.54 on CL/Q and 0.64 on Vc/Vp. Treg PD parameters do not carry weight covariates in Diao 2016.",
      source_name        = "WT"
    ),
    DOSE_50MG = list(
      description        = "Record-level indicator for the 50 mg SC dose (1 = 50 mg SC, 0 = any other SC dose or any IV dose)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (100, 150, 200, or 300 mg SC dose, or any IV dose)",
      notes              = "Inherited from the Othman 2014 PK backbone. Diao 2016 dosing is 150 or 300 mg SC every 4 weeks, so leave DOSE_50MG = 0 in clinical simulations.",
      source_name        = "(derived from AMT)"
    )
  )

  population <- list(
    n_subjects     = 1353L,
    n_records      = 8742L,
    n_studies      = 4L,
    study_names    = c("205MS201 / SELECT (Phase 2, RRMS)",
                       "205MS202 / SELECTION (Phase 2 extension with washout cohort)",
                       "205MS302 / OBSERVE (immunogenicity / PK / PD with intensive substudy)",
                       "205MS301 / DECIDE (Phase 3 vs IFN beta-1a)"),
    disease_state  = "Relapsing-remitting multiple sclerosis (RRMS)",
    dose_range     = "Daclizumab HYP 150 or 300 mg SC every 4 weeks",
    notes          = paste0(
      "Pooled PK/PD dataset of 1353 RRMS subjects with 8742 Treg records ",
      "from four daclizumab HYP clinical studies (Diao 2016 Table 2). ",
      "Treg defined as CD4+ CD127low/- Foxp3+ as a percentage of all CD4+ ",
      "T cells. Typical baseline percentage is 12.1% (Diao 2016 Table 5)."
    ),
    pd_subgroups = list(
      `205MS201/202 (SELECT/SELECTION)` = list(subjects = 545L, records = 4835L),
      `205MS302 (OBSERVE)`              = list(subjects = 106L, records =  891L),
      `205MS301 (DECIDE)`                = list(subjects = 702L, records = 3016L)
    )
  )

  ini({
    # ----------------------------------------------------------------------
    # PK backbone (copied verbatim from Othman_2014_daclizumab.R).
    # ----------------------------------------------------------------------
    lka      <- log(0.009 * 24); label("Absorption rate constant (Ka, 1/day)")                  # Othman 2014 Table 2
    lcl      <- log(0.010 * 24); label("Clearance for a 70 kg adult (CL, L/day)")               # Othman 2014 Table 2
    lvc      <- log(3.89);       label("Central volume of distribution for a 70 kg adult (Vc, L)") # Othman 2014 Table 2
    lvp      <- log(2.52);       label("Peripheral volume of distribution for a 70 kg adult (Vp, L)") # Othman 2014 Table 2
    lq       <- log(0.044 * 24); label("Inter-compartmental clearance for a 70 kg adult (Q, L/day)") # Othman 2014 Table 2
    lfdepot  <- log(0.84);       label("SC bioavailability for 100-300 mg doses (F, fraction)") # Othman 2014 Table 2
    ltlag    <- log(2 / 24);     label("Absorption lag time for SC doses (Tlag, day; 2 h)")     # Othman 2014 Table 2

    e_wt_cl_q <- 0.54; label("Allometric exponent on CL and Q (unitless)") # Othman 2014 Table 2
    e_wt_vc_vp <- 0.64; label("Allometric exponent on Vc and Vp (unitless)") # Othman 2014 Table 2

    e_dose_50mg_f <- -0.32143; label("Relative change in F for 50 mg SC vs 100-300 mg SC (fraction)") # Othman 2014 Table 2

    # Othman 2014 Table 2 (ka 58% CV, cl 27% CV, corr -0.72): ka -> omega^2 = log(1 + 0.58^2) = 0.29003;
    # cl -> omega^2 = log(1 + 0.27^2) = 0.07038; cov = -0.72 * sqrt(0.29003 * 0.07038) = -0.10290.
    etalka + etalcl ~ c(0.29003,
                        -0.10290, 0.07038)
    etalvc ~ 0.09175  # Othman 2014 Table 2 (Vc 31% CV)

    propSd <- 0.22; label("Proportional residual error on daclizumab HYP serum concentration (fraction)") # Othman 2014 Table 2
    addSd  <- 0.33; label("Additive residual error on daclizumab HYP serum concentration (ug/mL)")        # Othman 2014 Table 2

    # ----------------------------------------------------------------------
    # Treg PD parameters (Diao 2016 Table 5, sigmoidal Emax model).
    # Equation 3: Treg = E0 * (1 - Emax * Cc^gamma / (Cc^gamma + IC50^gamma))
    # The paper notes that adding a delay/effect compartment did not improve
    # fit (Discussion), so the algebraic form is used directly here.
    # ----------------------------------------------------------------------
    ltregE0    <- log(12.1);  label("Typical baseline Treg (% of CD4+ T cells)")            # Diao 2016 Table 5 (Baseline E0 = 12.1%)
    ltregIC50  <- log(3.97);  label("Treg IC50 (Cc giving 50% of Emax, mg/L)")              # Diao 2016 Table 5 (IC50 = 3.97 mg/L)
    treggamma  <- fixed(2);   label("Treg Hill coefficient (unitless)")                     # Diao 2016 Table 5 (Hill coefficient = 2; estimated, treated as fixed structural here)
    ltregEmax  <- log(0.610); label("Maximum fractional Treg reduction (unitless, 0-1)")    # Diao 2016 Table 5 (Emax = 0.610, i.e., 60% maximum reduction)

    # IIV (Diao 2016 Table 5): E0 42% CV, IC50 65% CV, Emax 12% CV.
    # omega^2 = log(1 + CV^2)
    # E0   CV 42% -> omega^2 = log(1 + 0.42^2) = 0.16127
    # IC50 CV 65% -> omega^2 = log(1 + 0.65^2) = 0.34744
    # Emax CV 12% -> omega^2 = log(1 + 0.12^2) = 0.01431
    etaltregE0   ~ 0.16127  # Diao 2016 Table 5 (Baseline E0 IIV 42% CV)
    etaltregIC50 ~ 0.34744  # Diao 2016 Table 5 (IC50 IIV 65% CV)
    etaltregEmax ~ 0.01431  # Diao 2016 Table 5 (Emax IIV 12% CV)

    # Residual error (Diao 2016 Table 5): proportional 50.1% + additive 0.416 (% units).
    propSd_treg <- 0.501; label("Proportional residual error on Treg (fraction)")                       # Diao 2016 Table 5
    addSd_treg  <- 0.416; label("Additive residual error on Treg (% of CD4+ T cells)")                  # Diao 2016 Table 5
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual PK parameters (Othman 2014 PK backbone).
    # ------------------------------------------------------------------
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    vp <- exp(lvp)          * (WT / 70)^e_wt_vc_vp
    q  <- exp(lq)           * (WT / 70)^e_wt_cl_q

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ------------------------------------------------------------------
    # 2. Two-compartment SC PK ODE system.
    # ------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    alag(depot) <- exp(ltlag)
    f(depot)    <- exp(lfdepot) * (1 + e_dose_50mg_f * DOSE_50MG)

    Cc <- central / vc

    # ------------------------------------------------------------------
    # 3. Individual Treg PD parameters.
    # ------------------------------------------------------------------
    tregE0I   <- exp(ltregE0   + etaltregE0)
    tregIC50I <- exp(ltregIC50 + etaltregIC50)
    tregEmaxI <- exp(ltregEmax + etaltregEmax)

    # ------------------------------------------------------------------
    # 4. Sigmoidal Emax Treg reduction (Diao 2016 Equation 3).
    # ------------------------------------------------------------------
    treg <- tregE0I * (1 - tregEmaxI * Cc^treggamma /
                                          (Cc^treggamma + tregIC50I^treggamma))

    # ------------------------------------------------------------------
    # 5. Observation and error model.
    # ------------------------------------------------------------------
    Cc   ~ add(addSd) + prop(propSd)
    treg ~ add(addSd_treg) + prop(propSd_treg)
  })
}
