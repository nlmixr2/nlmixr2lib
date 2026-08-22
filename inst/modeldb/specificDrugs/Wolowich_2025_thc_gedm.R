Wolowich_2025_thc_gedm <- function() {
  description <- paste(
    "Sequential population PK/PD model (Wolowich 2025 model 3GEDM) for the",
    "acute heart-rate effect of intravenous delta-9-tetrahydrocannabinol",
    "(THC) in healthy volunteers, in which the parent drug and its active",
    "metabolite 11-hydroxy-THC (11-OH-THC, written THC-OH by the authors)",
    "drive the response TOGETHER through the Gabrielsson and Weiner general",
    "empirical dynamic model (GEDM) for two interacting agonists. The PK",
    "backbone is the same joint parent-plus-metabolite model as the",
    "companion files, fitted to the first 5 h of data: THC is",
    "three-compartment with IV bolus input, all THC elimination is metabolic",
    "conversion to 11-OH-THC delayed by one transit compartment, and",
    "11-OH-THC is two-compartment. The pharmacodynamic endpoint is fHR, the",
    "increase in heart rate at time t expressed as a fraction of that",
    "individual's own maximal increase. Effect-site THC is the agonist and",
    "plasma 11-OH-THC is the interacting species; the four dimensionless",
    "hyperparameters gedm_alpha, gedm_beta, gedm_delta and gedm_gamma shape",
    "the interaction surface and Emax is fixed at 1. The fitted values",
    "satisfy the paper's criteria 4 and 5 (0 < beta < 1 with alpha > 0, and",
    "beta > delta), which the authors read as competitive plus",
    "non-competitive SYNERGISM between THC and 11-OH-THC. AIC was -2499. No",
    "covariate was retained: sex and CYP2C9 phenotype were both screened and",
    "neither was significant. Companion models from the same paper and the",
    "same PK backbone are modellib('Wolowich_2025_thc') (parent alone) and",
    "modellib('Wolowich_2025_thc_11oh') (metabolite alone)."
  )
  reference <- paste(
    "Wolowich WR, Greif R, Theiler L, Kleine-Brueggeney M.",
    "Pharmacokinetic/Pharmacodynamic Modeling of the Acute Heart Rate Effects",
    "of Delta-9 Tetrahydrocannabinol and Its Major Metabolites After",
    "Intravenous Injection in Healthy Volunteers.",
    "Eur J Drug Metab Pharmacokinet. 2025;50(3):229-241.",
    "doi:10.1007/s13318-025-00941-8.",
    "The GEDM interaction model itself is from Gabrielsson J, Weiner D.",
    "Pharmacokinetic and Pharmacodynamic Data Analysis: Concepts and",
    "Applications, cited as reference 20 of the paper."
  )
  vignette <- "Wolowich_2025_thc_heart_rate"

  # Plasma concentrations were converted to molar units by the authors
  # (THC 314.46 g/mol, 11-OH-THC 330.46 g/mol; Methods section 2.2), and
  # every PD parameter is reported in uM. Dosing therefore has to be
  # supplied in umol so that amount / volume lands in uM.
  units <- list(time = "h", dosing = "umol", concentration = "uM")

  # No covariate is retained in the final model. Methods section 2.4:
  # 'Stepwise (forward and backward, alpha = 0.05) covariate search for sex
  # and CYP2C9 phenotype were conducted.' Results section 3.2.3 closes with
  # 'Covariate analysis was conducted, but no significant covariates were
  # found.'
  covariateData <- list()

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened by stepwise forward / backward selection (alpha = 0.05) and",
        "not retained. Results section 3.2.3: 'Covariate analysis was",
        "conducted, but no significant covariates were found.' Cohort was 14",
        "of 25 female (56%). The Discussion adds that a previous study found",
        "no PK or PD difference between males and females for either THC or",
        "THC-OH."
      )
    ),
    CYP2C9_PM_IM = list(
      description = "Reduced-function CYP2C9 phenotype indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened and not retained. Three of the 25 volunteers were",
        "homozygous for CYP2C9*3. The Discussion explains the negative",
        "result: 'The CYP2C9 phenotype did not contribute to the HR effects",
        "observed as the polymorphism effect is seen in the PK of the",
        "terminal metabolite (THC-COOH) only', and THC-COOH was removed from",
        "the PK model because it had no relationship to fHR."
      )
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Used to compute the administered dose (0.1 mg/kg IV bolus) but not",
        "carried into the model: no volume or clearance is allometrically",
        "scaled, and Table 1 reports every disposition parameter as an",
        "absolute L or L/h with no per-70-kg normalisation. Cohort weight was",
        "65 (57-73) kg, median (IQR)."
      )
    )
  )

  compartmentData <- list(
    central = list(
      analyte = "delta-9-tetrahydrocannabinol", units = "umol",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "delta-9-tetrahydrocannabinol", units = "umol",
      specimen = "plasma", verified = TRUE
    ),
    peripheral2 = list(
      analyte = "delta-9-tetrahydrocannabinol", units = "umol",
      specimen = "plasma", verified = TRUE
    ),
    transit1_11oh = list(
      analyte = "11-hydroxy-delta-9-tetrahydrocannabinol", units = "umol",
      specimen = "not applicable", verified = TRUE
    ),
    central_11oh = list(
      analyte = "11-hydroxy-delta-9-tetrahydrocannabinol", units = "umol",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1_11oh = list(
      analyte = "11-hydroxy-delta-9-tetrahydrocannabinol", units = "umol",
      specimen = "plasma", verified = TRUE
    ),
    effect = list(
      analyte = "delta-9-tetrahydrocannabinol", units = "uM",
      specimen = "not applicable", verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 25L,
    n_studies      = 1L,
    n_observations = paste(
      "The full dataset held 275 THC, 275 THC-OH and 243 THC-COOH plasma",
      "concentrations plus 362 heart-rate records. PK sampling was at 0, 1,",
      "2, 5, 10, 15, 30, 45, 90, 180 and 300 min plus 24 and 48 h; heart rate",
      "was recorded before injection and at 1, 2, 5, 10, 20, 30, 45, 60, 75,",
      "90, 120, 150, 180 and 300 min. The PK/PD analysis was truncated at 5 h",
      "because the heart-rate effect did not last longer than that; Table 5",
      "reports NOBS = 817 for this model."
    ),
    age_median     = "23 years (IQR 21-25)",
    height_median  = "171 cm (IQR 167-182)",
    weight_median  = "65 kg (IQR 57-73)",
    sex_female_pct = 56,
    disease_state  = paste(
      "Healthy volunteers, cannabis naive or cannabis abstinent for at least",
      "one month. Exclusions: pregnancy, tobacco smoking within three months,",
      "suspected ischaemic heart disease, cardiac arrhythmias, BMI outside",
      "16-35 kg/m^2, hepatic P450 activity-altering medication, illicit drug",
      "use, and any lifetime history of treated or suspected psychiatric",
      "disease."
    ),
    dose_range     = paste(
      "Single 0.1 mg/kg intravenous bolus of THC (0.318 umol/kg using the",
      "paper's own MW of 314.46 g/mol)."
    ),
    genotype       = "Three of 25 volunteers were homozygous for CYP2C9*3.",
    baseline_hr    = "76 bpm (IQR 68-81); maximal change in HR 68 bpm (IQR 58-83)",
    regions        = "Switzerland (recovery room of a university anaesthesiology department, Bern)",
    notes          = paste(
      "Cantonal Ethics Committee Bern approval KEK 241-09; registered as",
      "ISRCTN53019164. Assays were LC-MS/MS with limits of quantification of",
      "0.8 ug/L (2 nM) for THC and 0.5 ug/L (1 nM) for THC-OH. Demographics",
      "are in Results section 3 paragraph 1. Fitting was done in Phoenix",
      "NLME; this model was additionally bootstrapped (n = 100) to obtain",
      "CV% and 90% confidence intervals, and the paper states the bootstrap",
      "parameters were not significantly different from the model fit",
      "parameters."
    )
  )

  ini({
    # ==================================================================
    # PK BACKBONE -- Wolowich 2025 Table 1 ('Pharmacokinetic model'),
    # the 5 h-truncated refit. Identical in all three models the paper
    # reports, because the PD layer was linked sequentially onto this
    # PK fit.
    #
    # Structural confirmation of the molar dose. Methods section 2.1
    # writes '0.1 mg/kg (3.18 uM/kg) THC intravenously', but the
    # parenthetical is wrong by a factor of 10: 0.1 mg / 314.46 g/mol =
    # 0.318 umol/kg. At 0.318 umol/kg a 65 kg volunteer receives
    # 20.7 umol and C0 = 20.7 / 5.2 = 3.97 uM, which decays to the
    # paper's own reported THC Cmax of 2.6 uM by the 1 min sample. At
    # 3.18 umol/kg C0 would be 39.7 uM, fifteen-fold above the observed
    # Cmax. The lower value is the one consistent with Table 1 and with
    # Figure 1A.
    # ==================================================================

    lvc <- log(5.2)
    label("Central volume of distribution for THC, Vc,THC (L)")                    # Table 1: VcTHC 5.2 (CV 10.0%, 95% CI 4.2-6.3)

    lvp <- log(14.9)
    label("First peripheral volume for THC, V2,THC (L)")                           # Table 1: V2THC 14.9 (CV 11.2%, 95% CI 11.6-18.2)

    lvp2 <- log(37.6)
    label("Second peripheral volume for THC, V3,THC (L)")                          # Table 1: V3THC 37.6 (CV 10.4%, 95% CI 29.9-45.3)

    lq <- log(43.7)
    label("Intercompartmental clearance to peripheral1 for THC, CLd1,THC (L/h)")   # Table 1: CLd1THC 43.7 (CV 5.9%, 95% CI 38.7-48.8)

    lq2 <- log(16.3)
    label("Intercompartmental clearance to peripheral2 for THC, CLd2,THC (L/h)")   # Table 1: CLd2THC 16.3 (CV 22.6%, 95% CI 9.0-23.5)

    lcl <- log(59.4)
    label("Clearance of THC, CL,THC (L/h)")                                        # Table 1: CLTHC 59.4 (CV 6.8%, 95% CI 51.4-67.3); Results 3.1 repeats 'the typical value of THC clearance was 59.4 l/h'

    lktr_11oh <- log(56.2)
    label("Transit rate constant from THC into 11-OH-THC, k-transit (1/h)")        # Table 1: k-transit THC-OH 56.2 (CV 25.8%, 95% CI 27.7-84.7). Table's '(l/h)' unit label is a slip -- it is a first-order rate constant, the same slip Table 5 makes when it labels Ke0 '(l/h)'. Results 3.1 lists the transit compartment as change (1) from the 48 h model.

    lvc_11oh <- log(65.3)
    label("Central volume of distribution for 11-OH-THC, Vc,THC-OH (L)")           # Table 1: VcTHC-OH 65.3 (CV 28.9%, 95% CI 28.3-102)

    lvp_11oh <- log(222)
    label("Peripheral volume for 11-OH-THC, V2,THC-OH (L)")                        # Table 1: V2THC-OH 222 (CV 19.8%, 95% CI 135-308)

    lq_11oh <- log(208)
    label("Intercompartmental clearance for 11-OH-THC, CLd,THC-OH (L/h)")          # Table 1: CLdTHCOH 208 (CV 32%, 95% CI 77.5-340)

    lcl_11oh <- log(223)
    label("Clearance of 11-OH-THC, CL,THC-OH (L/h)")                               # Table 1: CLTHC-OH 223 (CV 10.7%, 95% CI 176-270); Results 3.1 repeats 'the typical value of THC-OH clearance was 223 l/h'

    # ==================================================================
    # PD LAYER -- Wolowich 2025 Table 5, 'Model 3A: THC and THC-OH GEDM:
    # general empirical dynamic model', 'Value' column (the model fit).
    # Structure from Table 2 row 3GEDM:
    #   fHR = Emax,thc * (u1 + alpha*u2 + beta*u1*u2) /
    #                    (1  + u1 + delta*u2 + gamma*u1*u2)
    # with u1 = Ce,thc / EC50,thc and u2 = C,thc-oh / EC50,thc-oh.
    # Results 3.2.3: 'Effect site THC was the agonist and THC-OH was the
    # antagonist or synergist', so u1 uses the effect-compartment
    # concentration and u2 the plasma metabolite concentration.
    # ==================================================================

    lke0 <- log(0.26)
    label("Effect-compartment equilibration rate constant for THC, Ke0 (1/h)")     # Table 5: Ke0 0.26 (Value column; bootstrap 0.37, boot CV 60%, boot 95% CI 0.24-0.95). Table's '(l/h)' unit label is a slip -- Results 3.2.1 establishes that this row is a first-order rate constant by converting the model 1B1 value to a half-life.

    lec50 <- log(1.24)
    label("Effect-site THC concentration giving half-maximal fHR, EC50 (uM)")      # Table 5: EC50 THC 1.24 (Value; bootstrap 1.44, boot CV 18%, boot 95% CI 1.0-2.0)

    lec50_11oh <- log(0.12)
    label("Plasma 11-OH-THC concentration term of the interaction, EC50 (uM)")     # Table 5: EC50 THC-OH 0.12 (Value; bootstrap 0.13, boot CV 10%, boot 95% CI 0.1-0.15). Results 3.2.3 calls this 'EA50' in the prose describing the hyperparameters.

    lemax <- fixed(log(1))
    label("Maximum fractional heart-rate increase, Emax (unitless fraction)")      # Table 5: Emax (fHR) '1 (fixed)' in both the Value and the Boot value column

    # --- GEDM interaction hyperparameters --------------------------
    # Kept on the bare linear scale, not log-transformed: the paper's
    # Table 6 criteria compare these directly against 0, 1 and each
    # other, and 'beta = 0' (criterion 3, competitive interaction) has
    # to remain reachable.
    gedm_alpha <- 0.74
    label("GEDM efficacy of 11-OH-THC relative to THC, alpha (unitless)")          # Table 5: alpha 0.74 (Value and Boot value both 0.74; boot CV 8%, boot 95% CI 0.62-0.84)

    gedm_beta <- 0.94
    label("GEDM numerator interaction term for the doubly-occupied state, beta (unitless)")  # Table 5: beta 0.94 (Value; bootstrap 0.93, boot CV 10%, boot 95% CI 0.74-1.12)

    gedm_delta <- 0.48
    label("GEDM denominator term for 11-OH-THC occupancy, delta (unitless)")       # Table 5: delta 0.48 (Value; bootstrap 0.42, boot CV 19%, boot 95% CI 0.24-0.52)

    gedm_gamma <- 0.50
    label("GEDM denominator interaction term for the doubly-occupied state, gamma (unitless)")  # Table 5: gamma 0.50 (Value; bootstrap 0.44, boot CV 28%, boot 95% CI 0.25-0.54)

    # ------------------------------------------------------------------
    # Between-subject variability. Table 1 reports the PK IIV as
    # 'etaCV%'; Table 5 reports the PD IIV only in its 'Boot eta CV%'
    # column, i.e. derived from the n = 100 bootstrap, which is the only
    # between-subject variability the paper gives for this model.
    # omega^2 = log(1 + CV^2) for the log-normal etas.
    #
    # Eight etaCV% cells across the paper's tables are printed as the
    # string '< 1' rather than a number. Those are encoded fixed(0) --
    # the paper states only that the variability is below a 1% bound,
    # which is a statement that it is negligible rather than an estimate
    # of it. The reported '< 1' bound is recorded per line below and in
    # the vignette Errata so the information is not lost.
    # ------------------------------------------------------------------
    etalvc       ~ 3.59994e-05   # Table 1 VcTHC etaCV% 0.6;  log(1 + 0.006^2)
    etalvp       ~ fixed(0)      # Table 1 V2THC etaCV% reported as '< 1' (eta shrinkage 0.91)
    etalvp2      ~ 9.55817e-03   # Table 1 V3THC etaCV% 9.8;  log(1 + 0.098^2)
    etalq        ~ fixed(0)      # Table 1 CLd1THC etaCV% reported as '< 1' (eta shrinkage 0.51)
    etalq2       ~ 3.54637e-02   # Table 1 CLd2THC etaCV% 19;  log(1 + 0.19^2)
    etalcl       ~ 6.53957e-03   # Table 1 CLTHC etaCV% 8.1;  log(1 + 0.081^2)
    etalktr_11oh ~ fixed(0)      # Table 1 k-transit THC-OH etaCV% reported as '< 1' (eta shrinkage 0.88)
    etalvc_11oh  ~ 1.28305e-01   # Table 1 VcTHC-OH etaCV% 37;  log(1 + 0.37^2)
    etalvp_11oh  ~ fixed(0)      # Table 1 V2THC-OH etaCV% reported as '< 1' (eta shrinkage 0.91)
    etalq_11oh   ~ 6.54131e-02   # Table 1 CLdTHCOH etaCV% 26;  log(1 + 0.26^2)
    etalcl_11oh  ~ 1.42973e-02   # Table 1 CLTHC-OH etaCV% 12;  log(1 + 0.12^2)

    etalke0       ~ fixed(0)     # Table 5 Ke0 boot eta CV% reported as '< 1'
    etalec50      ~ 1.60944e+00  # Table 5 EC50 THC boot eta CV% 200;  log(1 + 2.00^2)
    etalec50_11oh ~ fixed(0)     # Table 5 EC50 THC-OH boot eta CV% reported as '< 1'
    etagedm_alpha ~ 1.20993e-04  # Table 5 alpha boot eta CV% 1.1;  log(1 + 0.011^2)
    etagedm_beta  ~ 2.49688e-03  # Table 5 beta boot eta CV% 5;  log(1 + 0.05^2)
    etagedm_delta ~ 1.59872e-03  # Table 5 delta boot eta CV% 4;  log(1 + 0.04^2)
    etagedm_gamma ~ fixed(0)     # Table 5 gamma boot eta CV% reported as '< 1'

    # ------------------------------------------------------------------
    # Residual error. Methods section 2.4 prints the PK error model as
    # 'Cobs * (1 + C,epsilon)', i.e. PROPORTIONAL, and the PD error model
    # as 'E,obs + E,epsilon', i.e. ADDITIVE. Every table labels the two
    # PK residual rows '(uM)' as though they were additive; the printed
    # equation is taken as authoritative. An additive SD of 0.22 uM would
    # exceed almost every post-distribution THC concentration (the LOQ is
    # 0.002 uM) and would put the lower limb of the Figure 1A visual
    # predictive check below zero on a log axis, which it plainly is not.
    # ------------------------------------------------------------------
    propSd <- 0.22
    label("Proportional residual error for THC (fraction)")                        # Table 5: epsilon (SD) THC 0.22 (Value and Boot value both 0.22; boot CV 7.2%, boot 95% CI 0.19-0.25)

    propSd_11oh <- 0.21
    label("Proportional residual error for 11-OH-THC (fraction)")                  # Table 5: epsilon (SD) THC-OH 0.21 (Value and Boot value both 0.21; boot CV 6.3%, boot 95% CI 0.18-0.23)

    addSd_fHR <- 0.22
    label("Additive residual error for fHR (unitless fraction)")                   # Table 5: epsilon (SD) fHR -- the 'Value' cell is BLANK in the published table; 0.22 is the 'Boot value' (boot CV 6.6%, boot 95% CI 0.20-0.26), the only number the paper gives for this row. See the vignette Errata.
  })

  model({
    # ---- 1. Individual parameters ----------------------------------
    vc  <- exp(lvc  + etalvc)
    vp  <- exp(lvp  + etalvp)
    vp2 <- exp(lvp2 + etalvp2)
    q   <- exp(lq   + etalq)
    q2  <- exp(lq2  + etalq2)
    cl  <- exp(lcl  + etalcl)

    ktr_11oh <- exp(lktr_11oh + etalktr_11oh)
    vc_11oh  <- exp(lvc_11oh  + etalvc_11oh)
    vp_11oh  <- exp(lvp_11oh  + etalvp_11oh)
    q_11oh   <- exp(lq_11oh   + etalq_11oh)
    cl_11oh  <- exp(lcl_11oh  + etalcl_11oh)

    ke0       <- exp(lke0       + etalke0)
    ec50      <- exp(lec50      + etalec50)
    ec50_11oh <- exp(lec50_11oh + etalec50_11oh)
    emax      <- exp(lemax)

    # The four hyperparameters stay on the linear scale (see ini()), but
    # Table 5 reports their between-subject variability as a bootstrap
    # CV%, which is a log-normal quantity -- so the eta enters
    # multiplicatively rather than additively.
    galpha <- gedm_alpha * exp(etagedm_alpha)
    gbeta  <- gedm_beta  * exp(etagedm_beta)
    gdelta <- gedm_delta * exp(etagedm_delta)
    ggamma <- gedm_gamma * exp(etagedm_gamma)

    # ---- 2. Micro-constants ----------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    kel_11oh <- cl_11oh / vc_11oh
    k12_11oh <- q_11oh  / vc_11oh
    k21_11oh <- q_11oh  / vp_11oh

    # ---- 3. ODE system ---------------------------------------------
    # THC: three compartments, IV bolus into central.
    d/dt(central)     <- -kel * central -
                          k12 * central + k21 * peripheral1 -
                          k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    # Metabolite formation. The paper reports no fraction metabolised, so
    # the whole of CL,THC is routed into 11-OH-THC (fm = 1) and the
    # metabolite volume and clearance are apparent (V/fm, CL/fm) -- the
    # standard identifiability convention when fm is not estimated, and
    # exact for the observed metabolite concentration. It is confirmed by
    # the paper's own numbers: Dose/CL,THC-OH = 20.7/223 = 0.093 uM*h
    # matches the metabolite AUC implied by Figure 1A (peak 0.1 uM near
    # 0.3 h, 0.005 uM at 5 h). THC and 11-OH-THC are both carried in
    # umol, so the conversion is 1:1 on a molar basis.
    d/dt(transit1_11oh)    <-  kel * central - ktr_11oh * transit1_11oh
    d/dt(central_11oh)     <-  ktr_11oh * transit1_11oh -
                               kel_11oh * central_11oh -
                               k12_11oh * central_11oh + k21_11oh * peripheral1_11oh
    d/dt(peripheral1_11oh) <-  k12_11oh * central_11oh - k21_11oh * peripheral1_11oh

    # ---- 4. Observations and error ---------------------------------
    Cc      <- central      / vc
    Cc_11oh <- central_11oh / vc_11oh

    # Effect compartment for the agonist (THC). The state carries a
    # CONCENTRATION, not an amount, so it equilibrates directly with Cc.
    d/dt(effect) <- ke0 * (Cc - effect)

    # Gabrielsson and Weiner general empirical dynamic model, Table 2
    # row 3GEDM. u1 is the reduced effect-site THC concentration and u2
    # the reduced plasma 11-OH-THC concentration.
    u1 <- effect  / ec50
    u2 <- Cc_11oh / ec50_11oh

    fHR <- emax * (u1 + galpha * u2 + gbeta  * u1 * u2) /
                  (1  + u1 + gdelta * u2 + ggamma * u1 * u2)

    Cc      ~ prop(propSd)
    Cc_11oh ~ prop(propSd_11oh)
    fHR     ~ add(addSd_fHR)
  })
}
