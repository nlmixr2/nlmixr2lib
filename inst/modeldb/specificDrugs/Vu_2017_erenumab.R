Vu_2017_erenumab <- function() {
  description <- "Two-compartment quasi-steady-state target-mediated drug disposition (QSS-TMDD) population PK model for erenumab (AMG 334), an anti-CGRP-receptor monoclonal antibody, with first-order SC absorption, allometric body-weight scaling of linear clearance and central volume, and a time-varying total CGRP-receptor pool; linked to a direct sigmoid inhibitory (Imax) model of the capsaicin-induced dermal blood flow (CIDBF) response in healthy subjects and migraine patients (Vu 2017)."
  reference <- "Vu T, Ma P, Chen JS, de Hoon J, Van Hecken A, Yan L, Wu LS, Hamilton L, Vargas G. Pharmacokinetic-Pharmacodynamic Relationship of Erenumab (AMG 334) and Capsaicin-Induced Dermal Blood Flow in Healthy and Migraine Subjects. Pharm Res. 2017;34(9):1784-1795. doi:10.1007/s11095-017-2183-6"
  vignette <- "Vu_2017_erenumab"

  # Time in days; doses in mg; all concentrations (erenumab, the total
  # receptor pool, Kss, ksyn, IC50) in ng/mL as printed in Tables III-IV.
  # The capsaicin-induced dermal blood flow response is in laser-Doppler
  # perfusion units (PU).
  units <- list(time = "day", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "erenumab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(
      analyte = "total (unbound + receptor-bound) erenumab",
      units = "mg",
      specimen = "serum",
      verified = TRUE
    ),
    peripheral1 = list(analyte = "unbound erenumab", units = "mg", specimen = "tissue", verified = TRUE),
    total_target = list(
      analyte = "total (free + erenumab-bound) CGRP receptor, in erenumab-mass-equivalent concentration",
      units = "ng/mL",
      specimen = "serum",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling of linear clearance (exponent fixed to 0.75) and central volume (exponent fixed to 1), centred at 70 kg (Table III footnote a). Observed range 49-104 kg (Discussion).",
      source_name = "weight"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Tested on CL and Vc; after adjusting for body weight not significant at p < 0.001 (Results, Pharmacokinetic Model).",
      source_name = "sex"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL, Vc and the PD parameters; not retained (Results; Supplementary Methods, Covariate Analyses).",
      source_name = "age"
    ),
    DIS_HEALTHY = list(
      description = "Healthy-participant indicator, 1 = healthy subject, 0 = migraine patient",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (healthy subject)",
      notes = "Population type was tested on PK and PD parameters and not retained; PK and DBF inhibition time courses were similar between healthy subjects and migraine patients (Results; Supplementary Fig. S3).",
      source_name = "type of population"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 108L,
    n_studies = 2L,
    n_pk_subjects = 78L,
    n_pd_subjects = 100L,
    n_pk_observations = 1297L,
    n_pd_observations = 676L,
    age_mean = "29.5 (SD 9.7) years",
    weight_mean = "74.9 (SD 11.9) kg",
    weight_range = "49-104 kg",
    sex_female_pct = 22,
    disease_state = "Healthy subjects (80, 74%) and episodic migraine patients (28, 26%)",
    dose_range = "SAD: SC 1, 7, 21, 70, 140 or 210 mg or IV 140 mg single dose (healthy), SC 140 mg (migraine). MAD: SC 21, 70, 140 mg Q4W x 3 or 280 mg x 1 + 210 mg x 2 Q4W (healthy); SC 21 or 140 mg Q4W x 3 (migraine).",
    regions = "Belgium (single centre, Leuven)",
    notes = "Pooled phase 1 single-ascending-dose (N = 60) and multiple-ascending-dose (N = 48) studies (Tables I-II). 30 of the 100 subjects with DBF data received placebo. 6 subjects developed binding anti-erenumab antibodies; their effect was not assessed. NONMEM 7.2.0, FOCE-I."
  )

  ini({
    # ---------------- PK (Table III; typical values at 70 kg) ----------------
    lcl <- log(0.214); label("Nonspecific linear clearance CL at 70 kg (L/day)") # Table III: CL = 0.214 L/day
    lvc <- log(4.27); label("Central volume of distribution Vc at 70 kg (L)") # Table III: Vc = 4.27 L
    lq <- log(3.34); label("Intercompartmental clearance Q (L/day)") # Table III: Q = 3.34 L/day
    lvp <- log(2.73); label("Peripheral volume of distribution Vp (L)") # Table III: Vp = 2.73 L
    lka <- log(0.432); label("First-order SC absorption rate constant ka (1/day)") # Table III: ka = 0.432 1/day
    # Eq. 1 estimates bioavailability on a logit scale (1/(1 + exp(F))); Table III
    # reports the back-transformed value 74%, carried here as logit(0.74).
    logitfdepot <- log(0.74 / (1 - 0.74)); label("SC bioavailability F (logit of fraction)") # Table III: F = 74%
    lksyn <- log(51.9); label("Zero-order receptor production rate ksyn (ng/mL/day)") # Table III: ksyn = 51.9 ng/mL/day
    lkdeg <- log(0.22); label("First-order receptor degradation rate constant kdeg (1/day)") # Table III: kdeg = 0.22 1/day
    lkss <- log(18.8); label("Quasi-steady-state binding constant Kss (ng/mL)") # Table III: Kss = 18.8 ng/mL
    lkint <- log(0.0345); label("Erenumab-receptor complex internalization rate constant kint (1/day)") # Table III: kint = 0.0345 1/day

    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL (unitless)") # Table III footnote a and Results: exponent fixed to 0.75
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on Vc (unitless)") # Table III footnote a and Results: exponent fixed to 1

    # ---------------- PD: direct sigmoid Imax on CIDBF (Table IV, Eq. 7) ----------------
    lrbase_cidbf <- log(398); label("Capsaicin-induced change in dermal blood flow in the absence of erenumab, dDBF0 (PU)") # Table IV: baseline dermal blood flow = 398 PU
    lic50 <- log(255); label("Serum erenumab concentration at half-maximal DBF inhibition IC50 (ng/mL)") # Table IV: IC50 = 255 ng/mL
    limax <- log(0.893); label("Maximum fractional inhibition of the capsaicin-induced DBF response Imax (fraction)") # Table IV: Imax = 89.3%
    lhill <- log(3.08); label("Hill coefficient of the concentration-inhibition curve (unitless)") # Table IV: Hill coefficient = 3.08

    # ---------------- Between-subject variability ----------------
    # Tables III-IV give BSV as %CV, where CV approximates sqrt(omega^2), so
    # omega^2 = CV^2. The printed covariances confirm that convention:
    #   -0.0731 / sqrt(0.251^2 * 0.595^2) = -0.490 (Results: 'correlations ... -0.49 between ka and CL')
    #    0.0682 / sqrt(0.251^2 * 0.423^2) =  0.642 (Results: '0.642 between CL and kint')
    # No ka-kint covariance is reported, so it is 0. Lower-triangle order:
    # var(CL); cov(CL,ka), var(ka); cov(CL,kint), cov(ka,kint), var(kint).
    etalcl + etalka + etalkint ~ c(
      0.063001,
      -0.0731, 0.354025,
      0.0682, 0, 0.178929
    ) # Table III: BSV CL 25.1%, ka 59.5%, kint 42.3% CV (squared); covariances CL-ka -0.0731, CL-kint 0.0682
    etalvc ~ 0.191844 # Table III: BSV in Vc = 43.8% CV, squared
    etalrbase_cidbf ~ 0.080089 # Table IV: BSV in baseline dermal blood flow = 28.3% CV, squared
    etalic50 ~ 0.885481 # Table IV: BSV in IC50 = 94.1% CV, squared
    etaRUV_cidbf ~ 0.149769 # Table IV: 'Residual error' 38.7% CV (shrinkage 33.3%), listed with the BSV rows: IIV on the DBF residual-error magnitude, squared

    # ---------------- Residual error ----------------
    # PK: fitted to ln-transformed concentrations, ln(Yobs) = ln(Ypred) + eps
    # (Supplementary Methods), i.e. log-normal error with log-scale SD sigma.
    expSd <- 0.179; label("Residual SD of ln erenumab concentration (log scale)") # Table III: residual variability = 17.9% CV
    # PD: Yobs = Ypred + W * eps, W combining additive and proportional terms
    # (Supplementary Methods; Table IV residual variability rows).
    addSd_cidbf <- 23.8; label("Additive residual error on the capsaicin-induced DBF response (PU)") # Table IV: residual variability = 23.8 PU
    propSd_cidbf <- 0.219; label("Proportional residual error on the capsaicin-induced DBF response (fraction)") # Table IV: residual variability = 21.9% CV
  })

  model({
    # ---- Individual PK parameters (Table III) ----
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q <- exp(lq)
    vp <- exp(lvp)
    ka <- exp(lka + etalka)
    fdepot <- expit(logitfdepot)
    ksyn <- exp(lksyn)
    kdeg <- exp(lkdeg)
    kss <- exp(lkss)
    kint <- exp(lkint + etalkint)

    # Micro-constants: kc0 = CL/Vc, kcp = Q/Vc, kpc = Q/Vp (Eqs. 3-4).
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- QSS algebra (Eq. 2) ----
    # Doses are mg and volumes L, so central/vc is mg/L = ug/mL; x 1000 gives
    # ng/mL, the units of Rtot, Kss and ksyn.
    ctot <- 1000 * central / vc
    rtot <- total_target
    qss <- ctot - rtot - kss
    cfree <- 0.5 * (qss + sqrt(qss * qss + 4 * kss * ctot))

    # ---- ODEs (Eqs. 1, 3, 4 and the Rtot equation under Eq. 2) ----
    # Unbound erenumab amount in central (mg) is cfree * vc / 1000; the bound
    # (complex) amount is central minus that.
    afree <- cfree * vc / 1000
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot + k21 * peripheral1 - kint * central - (k12 + kel - kint) * afree
    d/dt(peripheral1) <- k12 * afree - k21 * peripheral1
    d/dt(total_target) <- ksyn - kdeg * rtot - (kint - kdeg) * (ctot - cfree)

    # Rtot(0) = ksyn / kdeg (text under Eq. 2)
    total_target(0) <- ksyn / kdeg

    # Eq. 1: ASC(0) = F * DoseSC
    f(depot) <- fdepot

    # ---- PD: Eq. 7, driven by the serum (unbound) erenumab concentration ----
    rbase_cidbf <- exp(lrbase_cidbf + etalrbase_cidbf)
    ic50 <- exp(lic50 + etalic50)
    imax <- exp(limax)
    hill <- exp(lhill)
    inh <- imax * cfree^hill / (ic50^hill + cfree^hill)
    cidbf <- rbase_cidbf * (1 - inh)

    # ---- Observations ----
    # The measured serum erenumab concentration is the unbound concentration
    # of Eq. 2 ('Kss was estimated directly from the available unbound serum
    # erenumab concentration data').
    Cc <- cfree
    addSdi_cidbf <- addSd_cidbf * exp(etaRUV_cidbf)
    propSdi_cidbf <- propSd_cidbf * exp(etaRUV_cidbf)

    Cc ~ lnorm(expSd)
    cidbf ~ add(addSdi_cidbf) + prop(propSdi_cidbf)
  })
}
