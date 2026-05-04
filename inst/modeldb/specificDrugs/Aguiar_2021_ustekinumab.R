Aguiar_2021_ustekinumab <- function() {
  description <- "Population pharmacokinetic-pharmacodynamic model for ustekinumab in adults with Crohn's disease (Aguiar 2021): two-compartment quasi-equilibrium TMDD model for ustekinumab and the unbound IL-12/IL-23 p40 target, linked to fecal calprotectin via an indirect-response model with target-driven stimulation of FC production."
  reference <- "Aguiar Zdovc J, Hanzel J, Kurent T, Sever N, Smrekar N, Kozelj M, Novak G, Stabuc B, Drobne D, Grabnar I. Ustekinumab Dosing Individualization in Crohn's Disease Guided by a Population Pharmacokinetic-Pharmacodynamic Model. Pharmaceutics. 2021;13(10):1587. doi:10.3390/pharmaceutics13101587"
  vignette <- "Aguiar_2021_ustekinumab"
  units <- list(time = "day", dosing = "nmol", concentration = "nmol/L")

  covariateData <- list(
    FFM = list(
      description        = "Baseline fat-free mass",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL, Vc, Vp with reference 45 kg (Aguiar 2021 Table 2 footnotes a-c; the reference equals the cohort median FFM of 45 kg, Table 1). Derived per subject from height, weight, and sex via the Janmahasatian 2005 semi-mechanistic model (Aguiar 2021 Methods section 2.2).",
      source_name        = "FFM"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear-deviation effect on CL: multiplier (1 - 0.0165 * (ALB - 43)) (Aguiar 2021 Table 2 footnote a). Reference 43 g/L equals the cohort median (Table 1).",
      source_name        = "Serum albumin"
    ),
    CRP = list(
      description        = "Baseline C-reactive protein",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear-deviation effect on Ksyn (rate constant of target synthesis): multiplier (1 + 0.0846 * (CRP - 3)) (Aguiar 2021 Table 2 footnote d). Reference 3 mg/L equals the cohort median (Table 1). Standard CRP assay (Aguiar 2021 does not specify a high-sensitivity assay).",
      source_name        = "C-reactive protein"
    ),
    PRIOR_BIO = list(
      description        = "Prior biologic exposure indicator (any biologic; broader than anti-TNF only)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (biologic-naive). NOTE: Aguiar 2021 reports the inverted indicator 'bio-naive' (1 = biologic-naive); see notes.",
      notes              = "Aguiar 2021 Table 2 footnote a uses the inverted variable bio-naive = 1 - PRIOR_BIO. The model() block derives bio_naive <- 1 - PRIOR_BIO so the paper's reported coefficient (-0.227 on bio-naive) is preserved exactly. Effect on CL: multiplier (1 - 0.227 * bio_naive). In Aguiar 2021, prior biologic = prior anti-TNF (66.7% of cohort) and/or prior vedolizumab (17.5%); the paper's bio-naive group are subjects with no prior anti-TNF or vedolizumab.",
      source_name        = "bio-naive"
    ),
    FCGR3A_VV = list(
      description        = "FCGR3A 158 V/V homozygote indicator (rs396991)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (V/F heterozygote or F/F homozygote, combined)",
      notes              = "Aguiar 2021 Table 2: subcutaneous bioavailability F = 88.8% in V/V homozygotes vs 71.0% in V/F + F/F (combined). Modeled here on the logit scale: logit(F_pop) = log(0.71/0.29) for the reference group; e_fcgr3a_vv_fdepot adds log(0.888/0.112) - log(0.71/0.29) for V/V subjects. In the source TaqMan assay (rs396991, C_25815666_10), allele coding is A/C; the V allele is the higher-FcgammaRIIIa-affinity variant. Cohort distribution: 5/57 V/V (8.8%), 31/57 V/F (54.4%), 21/57 F/F (36.8%) (Table 1).",
      source_name        = "FCGR3A-158 V/V"
    ),
    ENDO_ULCER = list(
      description        = "Endoscopically active luminal disease at baseline (mucosal ulcerations confirmed at baseline ileocolonoscopy)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no mucosal ulcerations at baseline)",
      notes              = "Aguiar 2021 Table 3: covariate on baseline fecal calprotectin FC0. Typical FC0 = 102 mg/kg without ulcers and 213 mg/kg with ulcers; modeled here as a multiplicative effect on FC0 with reference 102 mg/kg. Cohort prevalence 44/57 (77.2%, Table 1).",
      source_name        = "Endoscopically active disease at baseline"
    )
  )

  population <- list(
    n_subjects     = 57L,
    n_studies      = 1L,
    age_range      = "median 49 years (IQR 32-56)",
    age_median     = "49 years",
    weight_range   = "median 70 kg (IQR 59-84)",
    weight_median  = "70 kg",
    ffm_range      = "median 45 kg (IQR 39-62)",
    height_range   = "median 169 cm (IQR 163-179)",
    sex_female_pct = 56.0,
    race_ethnicity = "Not reported (single-center Slovenian cohort)",
    disease_state  = "Active Crohn's disease starting ustekinumab therapy; 77.2% endoscopically active at baseline; median disease duration 14 years (IQR 7-22)",
    dose_range     = "Weight-tiered IV induction (260 mg if <=55 kg, 390 mg if 55-85 kg, 520 mg if >85 kg) at baseline, followed by 90 mg SC every 8 weeks (standard maintenance)",
    regions        = "Single-center, tertiary referral center, Slovenia",
    prior_biologic_pct = 66.7,
    prior_anti_tnf_pct = 66.7,
    prior_vedolizumab_pct = 17.5,
    fcgr3a_distribution = "V/V 8.8% (5/57), V/F 54.4% (31/57), F/F 36.8% (21/57); rs396991, Hardy-Weinberg in equilibrium",
    crp_baseline   = "median 3 mg/L (IQR 3-11)",
    albumin_baseline = "median 43 g/L (IQR 41-44)",
    fc_baseline    = "median 134 mg/kg (IQR 53-213)",
    samples        = "574 ustekinumab serum samples (5 below LOQ); 224 fecal calprotectin samples (15 below LOQ, 11 above ULOQ)",
    follow_up      = "32 weeks per subject; PK at baseline, 1 h post-IV-induction, weeks 2/4/8/9/10/12/16/20/24/32; FC at baseline and weeks 8/16/24/32",
    notes          = "Prospective observational study (Aguiar 2021). The PK assay is for unbound (free) ustekinumab via ImmunoGuide ELISA (LOQ 0.35 ug/mL); the modelled observation Cc is therefore free drug in nmol/L. Antidrug antibodies were not detected in any sample. Baseline demographics from Aguiar 2021 Table 1."
  )

  ini({
    # ----- Drug PK structure -----
    # Ustekinumab two-compartment disposition with first-order SC absorption.
    # CL, Vc, Vp typical values are for the reference subject (FFM = 45 kg,
    # serum albumin = 43 g/L, previously biologic-exposed; bio-naive = 0).
    lka  <- log(0.381); label("First-order SC absorption rate (Ka, 1/day)")              # Aguiar 2021 Table 2 final-model Ka
    lcl  <- log(0.277); label("Linear ustekinumab clearance at reference covariates (CL, L/day)") # Aguiar 2021 Table 2 final-model CL; footnote a
    lvc  <- log(3.57);  label("Central volume of distribution at reference FFM (Vc, L)") # Aguiar 2021 Table 2 final-model Vc; footnote b
    lvp  <- log(3.30);  label("Peripheral volume of distribution at reference FFM (Vp, L)") # Aguiar 2021 Table 2 final-model Vp; footnote c
    lq   <- log(1.89);  label("Drug intercompartmental clearance (Q, L/day)")            # Aguiar 2021 Table 2 final-model Q

    # Subcutaneous bioavailability on logit scale. Reference = V/F or F/F
    # (combined non-V/V); FCGR3A_VV adds an additive shift on the logit.
    logitfdepot <- log(0.710 / (1 - 0.710)); label("Logit of SC bioavailability for FCGR3A V/F or F/F (unitless)") # Aguiar 2021 Table 2 final-model F (V/F, F/F) = 71.0%
    e_fcgr3a_vv_fdepot <- log(0.888 / (1 - 0.888)) - log(0.710 / (1 - 0.710)); label("Additive logit shift on F for FCGR3A V/V vs reference (unitless)") # Aguiar 2021 Table 2 final-model F (V/V) = 88.8% vs F (V/F, F/F) = 71.0%

    # Covariate effects on CL and volumes (Aguiar 2021 Table 2 footnotes a-c).
    e_ffm_cl       <-  0.598;  label("Power exponent of FFM on CL (unitless)")                          # Aguiar 2021 Table 2 footnote a; FFM/45
    e_alb_cl       <- -0.0165; label("Linear-deviation coefficient of serum albumin on CL (per g/L, relative)") # Aguiar 2021 Table 2 footnote a; (1 - 0.0165 * (alb - 43))
    e_bionaive_cl  <- -0.227;  label("Multiplicative effect of bio-naive on CL (unitless)")             # Aguiar 2021 Table 2 footnote a; (1 - 0.227 * bio-naive); bio-naive = 1 - PRIOR_BIO
    e_ffm_vc       <-  0.590;  label("Power exponent of FFM on Vc (unitless)")                          # Aguiar 2021 Table 2 footnote b; FFM/45
    e_ffm_vp       <-  0.586;  label("Power exponent of FFM on Vp (unitless)")                          # Aguiar 2021 Table 2 footnote c; FFM/45

    # ----- Target PK structure -----
    # Two-compartment disposition for the unbound IL-12/IL-23 p40 target, with
    # binding to ustekinumab in the central compartment only. Ksyn is per-volume
    # (synthesis rate constant, nmol/L per day); Kdeg is the first-order target
    # degradation rate. Initial total-target concentration = Ksyn / Kdeg = T0.
    lksyn  <- log(9.86e-9);  label("Target synthesis rate constant at reference CRP (Ksyn, nmol/L/day)") # Aguiar 2021 Table 2 final-model Ksyn; footnote d
    lkdeg  <- log(9.26e-10); label("First-order target degradation rate constant (Kdeg, 1/day)")        # Aguiar 2021 Table 2 final-model Kdeg
    lvc_t  <- log(2.44);     label("Central volume of distribution for target (Vc-target, L)")          # Aguiar 2021 Table 2 final-model Vc-target
    lq_t   <- log(0.493);    label("Target intercompartmental clearance (Qtarget, L/day)")              # Aguiar 2021 Table 2 final-model Qtarget
    lvp_t  <- log(11.0);     label("Peripheral volume of distribution for target (Vp-target, L)")       # Aguiar 2021 Table 2 final-model Vp-target

    e_crp_ksyn <- 0.0846; label("Linear-deviation coefficient of baseline CRP on Ksyn (per mg/L, relative)") # Aguiar 2021 Table 2 footnote d; (1 + 0.0846 * (CRP - 3))

    # ----- Binding (quasi-equilibrium TMDD) -----
    lkint <- log(2.83e-6); label("Internalization rate constant of drug-target complex (Kint, 1/day)") # Aguiar 2021 Table 2 final-model Kint
    lkd   <- log(0.168);   label("Equilibrium dissociation constant for ustekinumab-target binding (Kd, nmol/L)") # Aguiar 2021 Table 2 final-model Kd

    # ----- PD: indirect-response model on fecal calprotectin -----
    # FC0 reference = patients without baseline ulcers (102 mg/kg).
    lkout       <- log(0.0581);     label("First-order FC degradation rate (Kout, 1/day)")             # Aguiar 2021 Table 3 final-model Kout
    lfc0        <- log(102);        label("Baseline fecal calprotectin in patients without baseline ulcers (FC0, mg/kg)") # Aguiar 2021 Table 3 final-model FC0 without ulcers
    e_endo_ulcer_fc0 <- log(213 / 102); label("Multiplicative log effect of baseline endoscopic ulcers on FC0 (unitless)") # Aguiar 2021 Table 3 final-model FC0 with ulcers (213 mg/kg) vs without (102 mg/kg)
    lemax_pd    <- log(2.19);       label("Maximum stimulation factor of unbound target on FC production (Emax, unitless)") # Aguiar 2021 Table 3 final-model Emax = 219%
    lc50_pd     <- log(2.46);       label("Unbound-target concentration at half-maximum FC stimulation (C50, nmol/L)") # Aguiar 2021 Table 3 final-model C50

    # ----- IIV (Aguiar 2021 Tables 2 and 3) -----
    # Log-normal (exponential model) for CL, Vc, Vp, Ksyn, FC0:
    #   omega^2 = log(CV^2 + 1) so the reported %CV maps to the log-scale variance.
    # Logit-additive for SC bioavailability (logitfdepot):
    #   IIV reported as SD on the logit scale = 0.173 -> omega^2 = 0.173^2.
    # Diagonal omega per Aguiar 2021 (no covariance reported between IIV terms).
    etalcl         ~ 0.03188   # 18.0% CV (Aguiar 2021 Table 2)  -> log(0.18^2 + 1) = 0.03188
    etalvc         ~ 0.00954   # 9.79% CV (Aguiar 2021 Table 2)  -> log(0.0979^2 + 1) = 0.00954
    etalvp         ~ 0.05645   # 24.1% CV (Aguiar 2021 Table 2)  -> log(0.241^2 + 1) = 0.05645
    etalogitfdepot ~ 0.02993   # SD = 0.173 on logit scale (Aguiar 2021 Table 2; reported as %, SD)
    etalksyn       ~ 0.68522   # 99.2% CV (Aguiar 2021 Table 2)  -> log(0.992^2 + 1) = 0.68522
    etalfc0        ~ 0.68309   # 99.0% CV (Aguiar 2021 Table 3)  -> log(0.99^2 + 1) = 0.68309

    # ----- Residual error -----
    # Combined additive + proportional on free ustekinumab in central (Cc, nmol/L).
    CcaddSd  <- 4.55;   label("Additive residual error on Cc (nmol/L)")                                # Aguiar 2021 Table 2 final-model additive RUV
    CcpropSd <- 0.0777; label("Proportional residual error on Cc (fraction)")                          # Aguiar 2021 Table 2 final-model proportional RUV (7.77%)

    # Proportional only on FC.
    propSd_fc <- 0.573;  label("Proportional residual error on FC (fraction)")                          # Aguiar 2021 Table 3 final-model proportional RUV (57.3%)
  })

  model({
    # ----- Derived covariate terms -----
    # Aguiar 2021 stores the source-paper indicator as bio-naive (1 = no prior
    # biologic). The canonical PRIOR_BIO column inverts that direction so the
    # paper's reported -0.227 coefficient applies to (1 - PRIOR_BIO).
    bio_naive <- 1 - PRIOR_BIO

    # ----- Individual ustekinumab disposition -----
    # FFM-power scaling on CL, Vc, Vp with reference 45 kg (Aguiar 2021 Table 2
    # footnotes a-c). Linear-deviation effects of serum albumin and bio-naive
    # on CL (footnote a). Reference CL/Vc/Vp values correspond to the reference
    # covariate set: FFM = 45 kg, serum albumin = 43 g/L, bio-naive = 0
    # (previously exposed to a biologic).
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) *
      (FFM / 45)^e_ffm_cl *
      (1 + e_alb_cl * (ALB - 43)) *
      (1 + e_bionaive_cl * bio_naive)
    vc <- exp(lvc + etalvc) * (FFM / 45)^e_ffm_vc
    vp <- exp(lvp + etalvp) * (FFM / 45)^e_ffm_vp
    q  <- exp(lq)

    # Subcutaneous bioavailability on the logit scale; FCGR3A_VV adds a fixed
    # shift relative to the V/F + F/F reference group.
    logit_f <- logitfdepot + etalogitfdepot + e_fcgr3a_vv_fdepot * FCGR3A_VV
    fdepot  <- 1 / (1 + exp(-logit_f))

    # ----- Individual target turnover and binding -----
    # Ksyn carries the CRP linear-deviation effect (footnote d) and a
    # log-normal eta. Initial total target concentration in central = Ksyn/Kdeg
    # (steady-state in the absence of drug). Peripheral target equilibrates at
    # the same baseline value (no drug -> only free target -> Tf = total).
    ksyn <- exp(lksyn + etalksyn) * (1 + e_crp_ksyn * (CRP - 3))
    kdeg <- exp(lkdeg)
    vc_t <- exp(lvc_t)
    vp_t <- exp(lvp_t)
    q_t  <- exp(lq_t)
    kint <- exp(lkint)
    kd   <- exp(lkd)

    t0 <- ksyn / kdeg

    # ----- PD parameters -----
    kout    <- exp(lkout)
    fc0     <- exp(lfc0 + e_endo_ulcer_fc0 * ENDO_ULCER + etalfc0)
    emax_pd <- exp(lemax_pd)
    c50_pd  <- exp(lc50_pd)

    # ----- Quasi-equilibrium TMDD relations in the central compartment -----
    # central holds the total ustekinumab amount (free + bound complex) in nmol;
    # total_target is the total target concentration (free + bound complex) in
    # nmol/L; binding occurs only in the central compartment.
    # Total drug concentration in central:
    ctot <- central / vc
    # Free drug Cf solves cf*tf = Kd*complex with ctot = cf + complex,
    # tf = total_target - complex:
    qe_disc <- ctot - total_target - kd
    cfree   <- 0.5 * (qe_disc + sqrt(qe_disc * qe_disc + 4 * kd * ctot))
    complex <- ctot - cfree
    tfree   <- total_target - complex

    # ----- ODE system -----
    # depot       : amount of ustekinumab in SC depot (nmol)
    # central     : total amount of ustekinumab in central, free + complex (nmol)
    # peripheral1 : amount of free ustekinumab in peripheral (nmol)
    # total_target    : total target concentration in central (nmol/L); only free
    #                     target degrades and distributes; only complex internalizes.
    # target_peripheral : free target concentration in peripheral (nmol/L); free
    #                     target equilibrates with central via Qtarget.
    # fc          : fecal calprotectin concentration (mg/kg).
    d/dt(depot)             <- -ka * depot
    d/dt(central)           <-  ka * depot -
                                cl * cfree -
                                q  * cfree + q * peripheral1 / vp -
                                kint * complex * vc
    d/dt(peripheral1)       <-  q * cfree - q * peripheral1 / vp
    d/dt(total_target)    <-  ksyn -
                                kdeg * tfree -
                                kint * complex -
                                (q_t / vc) * (tfree - target_peripheral)
    d/dt(target_peripheral) <-  (q_t / vp_t) * (tfree - target_peripheral)

    # FC indirect-response: unbound target stimulates FC production. Production
    # is normalized by the baseline stimulation so that FC = fc0 at SS without
    # drug (Aguiar 2021 Table 3 footnote a).
    stim_baseline <- 1 + emax_pd * t0    / (c50_pd + t0)
    stim_fc       <- 1 + emax_pd * tfree / (c50_pd + tfree)
    kin <- kout * fc0 / stim_baseline
    d/dt(fc) <- kin * stim_fc - kout * fc

    # Initial conditions: target at SS (T0), FC at baseline FC0.
    total_target(0)    <- t0
    target_peripheral(0) <- t0
    fc(0)                <- fc0

    # Bioavailability for SC depot (FCGR3A-stratified); IV induction doses go
    # directly to central (cmt = central) with f(central) = 1 by default.
    f(depot) <- fdepot

    # Observations.
    # Cc is unbound ustekinumab in central (nmol/L); the source ELISA assay
    # measures unbound drug (Aguiar 2021 Methods section 2.2, ImmunoGuide ELISA).
    Cc <- cfree
    Cc ~ add(CcaddSd) + prop(CcpropSd)

    # FC follows a proportional-only error model.
    fc ~ prop(propSd_fc)
  })
}
