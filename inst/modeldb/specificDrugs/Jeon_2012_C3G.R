Jeon_2012_C3G <- function() {
  description <- "One-compartment first-order absorption population PK model with an absorption lag time for cyanidin-3-glucoside (C3G) following 2-week multiple oral dosing of 1 g black bean (Phaseolus vulgaris, Cheongjakong-3-ho) seed coat extract once daily in 12 healthy adult Korean volunteers (Jeon 2012), with log-normal IIV on CL/F and V/F (with correlation rho = 0.883) and on Ka, and proportional residual error."
  reference <- paste(
    "Jeon S, Han S, Lee J, Hong T, Yim DS. (2012).",
    "The Safety and Pharmacokinetics of Cyanidin-3-Glucoside",
    "after 2-Week Administration of Black Bean Seed Coat Extract",
    "in Healthy Subjects. Korean J Physiol Pharmacol 16(4):249-253.",
    "doi:10.4196/kjpp.2012.16.4.249.",
    sep = " "
  )
  vignette <- "Jeon_2012_C3G"
  units    <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 12L,
    n_studies      = 1L,
    age_range      = "24-44 years (median 29.5)",
    age_median     = "29.5 years",
    weight_range   = "47.0-78.4 kg (median 59.6)",
    weight_median  = "59.6 kg",
    sex_female_pct = 50,
    race_ethnicity = c(Korean = 100),
    disease_state  = "Healthy adult Korean volunteers (aged 20-45 years) with no clinically relevant findings on medical history, physical examination, clinical laboratory, or ECG; subjects with any history that could alter C3G metabolism or with hypersensitivity to beans were excluded.",
    dose_range     = "Oral 1 g (200 mg x 5 capsules) of black bean (Phaseolus vulgaris, Cheongjakong-3-ho) seed coat extract once daily at approximately 9 AM for 14 days, after at least a 10 h overnight fast on the PK days (Day 1 and Day 14) and at least 4 h fasting thereafter. The C3G content per gram of extract was not reported; the C3G equivalent dose (used as the model's dose unit) is approximately 14.57 mg per 1 g of extract, back-calculated from the published CL/F and observed AUClast (see ini() comments).",
    regions        = "Republic of Korea (Clinical Trial Center, Seoul St Mary's Hospital).",
    notes          = "Demographics from Jeon 2012 Table 1; height median 166 cm (range 158-185). PK assays drew 7 mL venous blood predose and at 0.25, 0.5, 1, 1.5, 2, 4, and 6 h after dosing on Day 1 and Day 14, processed on ice within 10 min of collection and stored at -70 C until LC-MS/MS analysis (C3G LLOQ 0.2 ng/mL, linear range 0.2-50 ng/mL, intra- and inter-day accuracy 99.36-100.7%). Age, height, weight, sex, and creatinine clearance were screened as candidate covariates on CL, V, Ka, and ALAG using GAM (Xpose 4.0.4) with forward-addition/backward-elimination steps (alpha = 0.05); no covariate was retained in the final model (Jeon 2012 Results 'Final model' paragraph 2)."
  )

  ini({
    # Structural PK parameters from Jeon 2012 Table 3 (final one-compartment
    # first-order absorption model, NONMEM 6.2 ADVAN2 / TRANS2, FOCE-I). The
    # paper fits apparent oral parameters CL/F and V/F. The dose unit for
    # this packaged model is mg of C3G equivalent (not mg of raw extract):
    # 1 g of black bean seed coat extract corresponds to approximately
    # 14.57 mg of C3G equivalent dose (back-calculated as observed
    # AUClast x CL/F = 0.00426 mg*h/L x 3420 L/h = 14.57 mg, consistent with
    # the published CL/F and V/F treating bioavailability as F = 1). The
    # observed plasma C3G concentration is in ng/mL, so the model returns
    # Cc in ng/mL via Cc = 1000 * central / vc with central in mg and vc
    # in L (mg/L = ug/mL = 1000 ng/mL).
    lka    <- log(9.94);  label("Absorption rate constant Ka (1/h)")                 # Jeon 2012 Table 3: Ka = 9.94 1/h (%RSE 30.7)
    lcl    <- log(3420);  label("Apparent oral clearance CL/F (L/h)")                # Jeon 2012 Table 3: CL/F = 3420 L/h (%RSE 8.68)
    lvc    <- log(7280);  label("Apparent oral central volume of distribution V/F (L)")  # Jeon 2012 Table 3: V/F = 7280 L (%RSE 6.84)
    ltlag  <- log(0.217); label("Absorption lag time ALAG (h)")                      # Jeon 2012 Table 3: ALAG = 0.217 h (%RSE 3.67)

    # Log-normal IIV on CL/F and V/F (correlated, rho = 0.883) plus IIV on Ka.
    # The paper reports CV% per parameter; the internal variance is
    # omega^2 = log(1 + CV^2). Off-diagonal covariance for the block:
    #   var_cl = log(1 + 0.289^2) = 0.080216
    #   var_vc = log(1 + 0.213^2) = 0.044370
    #   cov    = 0.883 * sqrt(0.080216 * 0.044370) = 0.052679
    # IIV on Ka: var_ka = log(1 + 1.025^2) = 0.718145. The paper's Methods
    # section describes both BSV (eta) on CL/V/Ka/ALAG and IOV (kappa) on
    # CL/V/Ka, but the final-model Table 3 reports a single variance term on
    # Ka labelled 'Between-subject variability omega Ka'; the prose 'Final
    # model' paragraph states that inter-occasional variability between Day
    # 1 and Day 14 was retained only for Ka. This packaged model encodes the
    # 102.5% CV value as a single etalka variance, which captures the
    # absorption variability the paper reports without distinguishing BSV
    # from IOV (see vignette Assumptions and deviations). IIV on ALAG was
    # tested but not retained (Jeon 2012 Results 'Final model' paragraph 1).
    etalcl + etalvc ~ c(0.080216, 0.052679, 0.044370)  # Jeon 2012 Table 3: omega_CL = 28.9% CV, omega_V = 21.3% CV, rho_CL-V = 0.883
    etalka          ~ 0.718145                          # Jeon 2012 Table 3: omega_Ka = 102.5% CV; log(1 + 1.025^2) = 0.718145

    # Proportional residual error on the linear concentration scale. The
    # paper tested additive, proportional, and combined forms; the
    # combined-error example shown in Methods was for the exploratory phase,
    # and Table 3 reports only the retained proportional component
    # (sigma_prop = 21.3%, %RSE 9.30).
    propSd <- 0.213; label("Proportional residual error on plasma C3G concentration (fraction)")  # Jeon 2012 Table 3: sigma_prop = 21.3% (%RSE 9.30)
  })

  model({
    # Individual PK parameters (no retained covariates in the final model).
    ka    <- exp(lka + etalka)
    cl    <- exp(lcl + etalcl)
    vc    <- exp(lvc + etalvc)
    alag1 <- exp(ltlag)

    # One-compartment first-order absorption with first-order elimination
    # from the central compartment. Dose lands in the depot via the user
    # data set's cmt column and is delayed by the absorption lag alag1.
    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    alag(depot) <- alag1

    # Plasma C3G concentration in ng/mL. Dose unit is mg of C3G equivalent
    # and vc is in L, so central/vc gives mg/L = ug/mL; multiplying by 1000
    # converts to ng/mL, the bioanalytical assay unit (LLOQ 0.2 ng/mL).
    Cc <- 1000 * central / vc

    Cc ~ prop(propSd)
  })
}
