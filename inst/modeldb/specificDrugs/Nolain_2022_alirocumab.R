Nolain_2022_alirocumab <- function() {
  description <- "Semi-mechanistic population PK/PD model for alirocumab, total PCSK9 and LDL cholesterol in healthy adults and adults with hypercholesterolaemia (Nolain 2022). Two-compartment alirocumab disposition with first-order SC absorption (lag time, logit-normal bioavailability), quasi-steady-state target-mediated drug disposition (TMDD-QSS) binding to PCSK9 with zero-order PCSK9 synthesis, and an indirect-response (type II, inhibition of loss) LDL-C model in which free PCSK9 inhibits LDL-C degradation through a sigmoid Imax function. Statin co-administration increases Vc and Imax; baseline total PCSK9 scales the free-PCSK9 IC50."
  reference <- "Nolain P, Djebli N, Brunet A, Fabre D, Khier S. Combined Semi-mechanistic Target-Mediated Drug Disposition and Pharmacokinetic-Pharmacodynamic Models of Alirocumab, PCSK9, and Low-Density Lipoprotein Cholesterol in a Pooled Analysis of Randomized Phase I/II/III Studies. Eur J Drug Metab Pharmacokinet. 2022;47:789-802. doi:10.1007/s13318-022-00787-4"
  vignette <- "Nolain_2022_alirocumab"
  units <- list(time = "day", dosing = "mg", concentration = "nM")

  covariateData <- list(
    TPCSK9_BASE = list(
      description = "Baseline (pre-dose) total serum PCSK9 concentration, time-fixed per subject",
      units = "nM",
      type = "continuous",
      reference_category = NULL,
      notes = "Nolain 2022 uses the individual observed baseline total PCSK9 twice: (1) as the initial condition and synthesis anchor of the total-PCSK9 state (Ptot(0) = [Ptot]baseline * Vc; ksyn = kdeg * [Ptot]baseline; Table 3) and in kout(0) (Table 3), and (2) as the covariate TBSPCSK9 on the free-PCSK9 IC50, power form centred on the dataset median 6.99 nM (Eq. 6). Units are nM in this model, not the register default ng/mL. Modelling dataset mean 7.66 nM (SD 3.06), range 2.36-19.6 nM (Table 2). The paper's a priori effect of TBSPCSK9 on kdeg (Eq. 7) was removed in the covariate-reduction step (Sect. 3.2.2) and is not encoded.",
      source_name = "TBSPCSK9"
    ),
    LDLC = list(
      description = "Baseline (pre-dose) LDL cholesterol, time-fixed per subject",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Individual observed baseline LDL-C ([LDLC]baseline, Table 3). Sets the initial condition of the `ldl` turnover state and, through the pre-treatment steady state kin = kout(0) * [LDLC]baseline, the zero-order LDL-C production rate. Because the model is linear in the state, percent change from baseline is independent of this value. Modelling dataset mean 140 mg/dL (SD 33.1), range 88.5-356 (Table 2).",
      source_name = "[LDLC]baseline"
    ),
    CONMED_STATIN = list(
      description = "Concomitant statin therapy",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no statin)",
      notes = "STATIN = 1 with statin co-administration, 0 otherwise (Eqs. 3 and 5). Multiplicative 1.75-fold effect on Vc (Eq. 3) and additive +0.140 effect on the typical Imax (Eq. 5; 74.1% -> 88.1%). The a priori statin effect on kout (Eq. 4) was removed in the covariate-reduction step (Sect. 3.2.2) and is not encoded. 60.9% of the modelling dataset received a statin (Table 2).",
      source_name = "STATIN"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Additive effect of female sex on Imax (Eq. 5; SEX = 1 for women) was included a priori in the full fixed-effects model but removed in the reduction step as not clinically relevant (Sect. 3.2.2, Fig. 2c). No final-model estimate is reported.",
      source_name = "SEX"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "alirocumab", units = "nmol", specimen = "administration site", verified = TRUE),
    central = list(
      analyte = "alirocumab (total: free + PCSK9-bound)",
      units = "nmol",
      specimen = "serum",
      verified = TRUE
    ),
    peripheral1 = list(analyte = "alirocumab (free)", units = "nmol", specimen = "not applicable", verified = TRUE),
    total_target = list(
      analyte = "PCSK9 (total: free + alirocumab-bound)",
      units = "nmol",
      specimen = "serum",
      verified = TRUE
    ),
    ldl = list(analyte = "low-density lipoprotein cholesterol", units = "mg/dL", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 527L,
    n_studies = 9L,
    n_observations = 13731L,
    phases = "Phase I, II and III",
    age_range = "18-75 years",
    age_mean = "52.5 years (SD 13.0)",
    weight_range = "45.8-154 kg",
    weight_mean = "80.6 kg (SD 16.4)",
    sex_female_pct = 46.9,
    disease_state = "Healthy volunteers (28.5%) and patients with hypercholesterolaemia (71.5%), all with LDL-C >= 100 mg/dL; 30.0% heterozygous familial hypercholesterolaemia.",
    co_medication = "Statin 60.9% (low dose 38.5%, high dose 22.4%), ezetimibe 12.9%, fibrate 4.74%; 29.8% alirocumab alone.",
    dose_range = "Single IV doses 0.3-12 mg/kg (one phase I study, NCT01026597); SC 50-300 mg single dose or Q2W/Q4W for up to 104 weeks.",
    pcsk9_baseline = "Baseline total PCSK9 mean 7.66 nM (SD 3.06), range 2.36-19.6 nM; median 6.99 nM (Sect. 3.2.2).",
    ldlc_baseline = "Baseline LDL-C mean 140 mg/dL (SD 33.1), range 88.5-356 mg/dL.",
    renal_function = "MDRD creatinine clearance mean 109 mL/min/1.73 m^2 (SD 30.4), range 38.1-253.",
    notes = "Modelling dataset (Table 1 upper block, Table 2 first column). An external validation dataset of 2273 patients from four further phase II/III studies (Table 1 lower block) was used for MAP-Bayesian predictive checks only."
  )

  ini({
    # ---- Alirocumab PK (Nolain 2022 Table 4, final model; reference = no statin) ----
    lcl <- log(0.221); label("Linear clearance of free alirocumab CL (L/day)") # Table 4 'theta CL' 0.221 L/day
    lvc <- log(3.20); label("Central volume of distribution Vc without statin (L)") # Table 4 'theta Vc' 3.20 L
    lq <- log(0.557); label("Inter-compartmental clearance of free alirocumab Q (L/day)") # Table 4 'theta Q' 0.557 L/day
    lvp <- fixed(log(2.61)); label("Peripheral volume of free alirocumab Vp (L)") # Table 4 'theta VP' 2.61 L, no RSE; Sect. 3.2.1 'set at 2.61 L'
    lka <- log(0.346); label("First-order SC absorption rate constant ka (1/day)") # Table 4 'theta ka' 0.346 /day
    logitfdepot <- logit(0.681); label("Logit of absolute SC bioavailability F (F = 0.681)") # Table 4 'theta F' 68.1%; logit-normal (Sect. 3.2.1)
    ltlag <- log(0.0283); label("SC absorption lag time LAG (day)") # Table 4 'theta LAG' 0.0283 day

    # ---- TMDD-QSS binding and PCSK9 turnover (Table 3, Table 4) ----
    lkint <- log(0.127); label("Rate constant for alirocumab-PCSK9 complex clearance kclear (1/day)") # Table 4 'theta kclear' 0.127 /day
    lkdeg <- log(1.34); label("First-order degradation rate constant of free PCSK9 kdeg (1/day)") # Table 4 'theta kdeg' 1.34 /day
    lkd <- fixed(log(0.58)); label("Equilibrium dissociation constant kD = koff/kon (nM)") # Table 4 'theta kD' 0.58 nM, no RSE; Table 3 'set to 0.58 nM from in vitro experiments'
    lk1 <- fixed(log(559)); label("Association rate constant kon of the drug-target complex (1/nM/day)") # Table 4 'theta kon' 559, no RSE; Sect. 3.2.1 'set at 559/nM/day'

    # ---- LDL-C indirect response (Table 3, Table 4) ----
    lkout <- log(0.260); label("First-order LDL-C degradation rate constant during treatment kout (1/day)") # Table 4 'theta kout' 0.260 /day
    logitimax <- logit(0.741); label("Logit of maximal inhibition of LDL-C degradation by free PCSK9 without statin (Imax = 0.741)") # Table 4 'theta Imax' 74.1%; logit-normal (Sect. 3.2.1)
    lki50 <- log(6.03); label("Free PCSK9 concentration giving half of Imax at the median baseline total PCSK9, IC50 (nM)") # Table 4 'theta IC50' 6.03 nM
    lhill <- log(11.6); label("Hill coefficient gamma of the free-PCSK9 Imax function (unitless)") # Table 4 'theta gamma' 11.6

    # ---- Covariate effects (Eqs. 3, 5, 6; Table 4) ----
    e_conmed_statin_vc <- 1.75; label("Multiplicative factor on Vc with statin co-administration, Vc = theta * factor^STATIN (unitless)") # Table 4 'theta V_STATIN' 1.75; Eq. 3
    e_conmed_statin_imax <- 0.140; label("Additive increment in typical Imax with statin co-administration (fraction)") # Table 4 'theta Imax STATIN' 0.140; Eq. 5 (74.1% -> 88.1%)
    e_tpcsk9_base_ki50 <- 0.930; label("Power exponent of baseline total PCSK9 (TPCSK9_BASE / 6.99 nM) on IC50 (unitless)") # Table 4 'theta IC50 TBSPCSK' 0.930; Eq. 6

    # ---- Between-subject variability (Table 4 variances; CV = sqrt(exp(omega^2) - 1)) ----
    etalcl ~ 0.270 # Table 4 omega2 CL 0.270 (55.7%)
    etalkint ~ 0.0554 # Table 4 omega2 kclear 0.0554 (23.9%)
    etalkdeg ~ 0.124 # Table 4 omega2 kdeg 0.124 (36.4%)
    etalvc ~ 0.0648 # Table 4 omega2 Vc 0.0648 (25.9%)
    etalka ~ 0.344 # Table 4 omega2 ka 0.344 (64.1%)
    etalogitfdepot ~ 0.626 # Table 4 omega2 F 0.626 (logit scale)
    etalkout ~ 0.256 # Table 4 omega2 kout 0.256 (54.0%)
    etalogitimax ~ 0.146 # Table 4 omega2 Imax 0.146 (logit scale)
    etalki50 ~ 0.00578 # Table 4 omega2 IC50 0.00578 (7.61%)

    # ---- Residual error: combined additive + proportional per endpoint (Table 4) ----
    addSd <- 0.426; label("Additive residual error, total alirocumab (nM)") # Table 4 'ALIROCUMAB_ADD' 0.426 nM
    propSd <- 0.255; label("Proportional residual error, total alirocumab (fraction)") # Table 4 'ALIROCUMAB_PROP' 25.5%
    addSd_Ctotal_target <- 1.07; label("Additive residual error, total PCSK9 (nM)") # Table 4 'TPCSK9_ADD' 1.07 nM
    propSd_Ctotal_target <- 0.279; label("Proportional residual error, total PCSK9 (fraction)") # Table 4 'TPCSK9_PROP' 27.9%
    addSd_ldl <- 5.71; label("Additive residual error, LDL-C (mg/dL)") # Table 4 'LDLC_ADD' 5.71 mg/dL
    propSd_ldl <- 0.142; label("Proportional residual error, LDL-C (fraction)") # Table 4 'DLC_PROP' 14.2%
  })

  model({
    # Alirocumab molecular weight; not stated in Nolain 2022, which doses in nmol
    # (Table 3). 146 kDa is the approximate molecular weight on the Praluent label.
    mw_alirocumab <- 146000 # g/mol
    nmol_per_mg <- 1e6 / mw_alirocumab

    # ---- Individual parameters ----
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc) * e_conmed_statin_vc^CONMED_STATIN
    q <- exp(lq)
    vp <- exp(lvp)
    ka <- exp(lka + etalka)
    fdepot <- expit(logitfdepot + etalogitfdepot)
    tlag <- exp(ltlag)

    kint <- exp(lkint + etalkint)
    kdeg <- exp(lkdeg + etalkdeg)
    kd <- exp(lkd)
    k1 <- exp(lk1)
    # Table 3: kss = (koff + kclear) / kon = kD + kclear / kon
    kss <- kd + kint / k1

    kout <- exp(lkout + etalkout)
    # Eq. 5: statin effect is additive on the probability scale; the logit-normal
    # between-subject variability is applied around that typical value.
    imax_typ <- expit(logitimax) + e_conmed_statin_imax * CONMED_STATIN
    imax <- expit(logit(imax_typ) + etalogitimax)
    ki50 <- exp(lki50 + etalki50) * (TPCSK9_BASE / 6.99)^e_tpcsk9_base_ki50
    hill <- exp(lhill)

    # ---- Micro-constants (Table 3) ----
    kel <- cl / vc
    kcp <- q / vc
    kpc <- q / vp

    # ---- Pre-treatment steady state (Table 3) ----
    ksyn <- kdeg * TPCSK9_BASE
    kout0 <- kout * (1 - imax * TPCSK9_BASE^hill / (ki50^hill + TPCSK9_BASE^hill))
    kin <- kout0 * LDLC

    # ---- QSS free drug (Table 3) ----
    # [Afree] = (a + sqrt(a^2 + 4 kss [Atot])) / 2 with a = [Atot] - [Ptot] - kss.
    # When a < 0 (drug below target) the printed form cancels to 0 - 0; the
    # algebraically identical 2 kss [Atot] / (sqrt(.) - a) is used there.
    atot <- central / vc
    ptot <- total_target / vc
    qss_a <- atot - ptot - kss
    qss_s <- sqrt(qss_a^2 + 4 * kss * atot)
    if (qss_a < 0) {
      afree <- 2 * kss * atot / (qss_s - qss_a)
    } else {
      afree <- (qss_a + qss_s) / 2
    }
    # [Complex] = [Ptot] [Afree] / (kss + [Afree]); [Pfree] = [Ptot] - [Complex]
    fbound <- afree / (kss + afree)
    pfree <- ptot * (1 - fbound)
    cplx <- ptot * fbound

    # ---- ODEs (Table 3; amounts in nmol, LDL-C in mg/dL) ----
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - (kel + kcp) * afree * vc + kpc * peripheral1 - kint * total_target * fbound
    d/dt(peripheral1) <- kcp * afree * vc - kpc * peripheral1
    d/dt(total_target) <- ksyn * vc - kdeg * total_target - (kint - kdeg) * total_target * fbound
    d/dt(ldl) <- kin - kout * (1 - imax * pfree^hill / (ki50^hill + pfree^hill)) * ldl

    total_target(0) <- TPCSK9_BASE * vc
    ldl(0) <- LDLC

    # SC doses in mg enter the depot (bioavailability F, lag); IV doses in mg go
    # to central. Both are converted to nmol.
    f(depot) <- fdepot * nmol_per_mg
    alag(depot) <- tlag
    f(central) <- nmol_per_mg

    # ---- Observations ----
    Cc <- atot
    Ctotal_target <- ptot
    Cc ~ add(addSd) + prop(propSd)
    Ctotal_target ~ add(addSd_Ctotal_target) + prop(propSd_Ctotal_target)
    ldl ~ add(addSd_ldl) + prop(propSd_ldl)
  })
}
