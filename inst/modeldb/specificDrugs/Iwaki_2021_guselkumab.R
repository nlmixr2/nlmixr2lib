Iwaki_2021_guselkumab <- function() {
  description <- paste(
    "PK/PD model for guselkumab (anti-IL-23 p19 human IgG1-lambda mAb) and the",
    "Palmoplantar Pustulosis Area and Severity Index (PPPASI) in Japanese adults",
    "with palmoplantar pustulosis (Iwaki 2021). PK: the one-compartment model",
    "with first-order SC absorption of Yao 2018, applied unchanged. PD: the",
    "PPPASI score is a bounded outcome score analysed as an ordered categorical",
    "variable (721 categories of 0.1 points) whose probabilities are intervals",
    "of a beta distribution with mean mu and precision tau;",
    "logit(mu) = logit(baseline) - placebo effect - drug effect + linear disease",
    "progression, where the placebo effect is Pmax * (1 - exp(-kp * t)) and the",
    "drug effect is Emax * (1 - R), with R a type I indirect-response latent",
    "variable inhibited by serum guselkumab (IC50 fixed at 1e-6 ug/mL)."
  )
  reference <- paste(
    "Iwaki Y, Shibata S, Hu C. Pharmacokinetic/Pharmacodynamic Analysis of",
    "Guselkumab for Treatment of Palmoplantar Pustulosis: Clinical Implications",
    "of Guselkumab Dose, Disease Severity, and Smoking in Japanese Patients.",
    "J Clin Pharmacol. 2022;62(2):182-189. doi:10.1002/jcph.1953.",
    "PK structure and parameters taken from Yao Z et al. J Clin Pharmacol.",
    "2018;58(5):613-627. doi:10.1002/jcph.1063; see",
    "modellib('Yao_2018_guselkumab')."
  )
  vignette <- "Iwaki_2021_guselkumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL/F and V/F from Yao 2018 Table 4 footnotes b and c, reference 87.1 kg (the Yao 2018 psoriasis median). Iwaki 2021 Methods 'Population PK Modeling' applies this model to the PPP cohort, whose median baseline weight is 59.7 kg (range 37.0-98.2). Time-fixed at baseline.",
      source_name = "BWT"
    ),
    DIS_DIAB = list(
      description = "Past or current diabetes mellitus comorbidity",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no diabetes mellitus)",
      notes = "Multiplier 1.12 on CL/F from Yao 2018 Table 4 footnote b. Iwaki 2021 lists diabetic comorbidity among the Yao covariates but does not report its prevalence in the PPP cohort.",
      source_name = "DIAB"
    ),
    RACE_WHITE = list(
      description = "White race indicator (1 = White, 0 = non-White)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-White)",
      notes = "Multiplier 1.11 on CL/F for non-White subjects from Yao 2018 Table 4 footnote b, where the source RACE = 1 for non-White; RACE_WHITE = 1 - RACE, so the factor is 1.11^(1 - RACE_WHITE). Every subject in the Iwaki 2021 PPP cohort was Japanese (Results, 'Exploratory Assessments of the Observed Data'), i.e. RACE_WHITE = 0.",
      source_name = "RACE"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "guselkumab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "guselkumab", units = "mg", specimen = "serum", verified = TRUE),
    latent_pppasi = list(
      analyte = "latent drug-effect variable R(t)",
      units = "(unitless)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 207L,
    n_studies = 2L,
    weight_range = "37.0-98.2 kg",
    weight_median = 59.7,
    sex_female_pct = 77.3,
    race_ethnicity = "All Asian (Japanese)",
    smoker_pct = 60.4,
    baseline_pppasi = "median 24.6 (range 4.8-59.4) on the 0-72 scale",
    disease_state = "Adults with palmoplantar pustulosis (PPP). Pooled data from a phase 2 proof-of-concept study and a phase 3 confirmatory study of guselkumab, both randomized, double-blind and placebo-controlled at multiple sites in Japan.",
    dose_range = "Phase 2 (n = 49): guselkumab 200 mg or placebo SC at weeks 0 and 4. Phase 3 (n = 159): guselkumab 100 or 200 mg SC at weeks 0, 4 and 12 then every 8 weeks, or placebo at weeks 0, 4 and 12; the placebo group was rerandomized 1:1 to 100 or 200 mg at weeks 16 and 20 then every 8 weeks.",
    regions = "Japan",
    notes = "Iwaki 2021 Methods 'Clinical Study Designs' and Results 'Exploratory Assessments of the Observed Data': 207 patients, 160 women (77%), 125 smokers (60%), 1813 PK and 2718 PPPASI observations; efficacy data through week 24 (phase 2) and week 52 (phase 3). Age is not reported. The PK parameters come from the Yao 2018 psoriasis population (median weight 87.1 kg, 71% White, 84% men)."
  )

  ini({
    # PK: Yao 2018 Table 4, final reduced model. Iwaki 2021 Methods 'Population
    # PK Modeling' applies this model unchanged and states its structure and
    # covariates; Iwaki 2021 reports no PK parameter values of its own.
    # Reference subject: 87.1 kg, no diabetes, White.
    lcl <- log(0.516); label("Apparent clearance CL/F at reference covariates (L/day)") # Yao 2018 Table 4: CL/F = 0.516 L/day
    lvc <- log(13.5); label("Apparent volume of distribution V/F at reference covariates (L)") # Yao 2018 Table 4: V/F = 13.5 L
    lka <- log(1.11); label("First-order SC absorption rate constant ka (1/day)") # Yao 2018 Table 4: ka = 1.11 1/day

    e_wt_cl <- 0.998; label("Power exponent of body weight on CL/F (unitless)") # Yao 2018 Table 4 footnote b
    e_wt_vc <- 0.829; label("Power exponent of body weight on V/F (unitless)") # Yao 2018 Table 4 footnote c
    e_diab_cl <- 1.12; label("Multiplier on CL/F for diabetes, power form on DIS_DIAB (unitless)") # Yao 2018 Table 4 footnote b
    e_race_cl <- 1.11; label("Multiplier on CL/F for non-White race, power form on 1 - RACE_WHITE (unitless)") # Yao 2018 Table 4 footnote b

    # Yao 2018 Table 4 IIV is (variance)^(1/2) x 100%, so omega^2 = (IIV/100)^2:
    # CL/F 35.6% -> 0.126736; V/F 28.0% -> 0.078400; ka 129% -> 1.664100;
    # cov(CL/F, V/F) = 0.834 * sqrt(0.126736 * 0.078400) = 0.083133.
    etalcl + etalvc ~ c(
      0.126736,
      0.083133, 0.078400
    ) # Yao 2018 Table 4: IIV CL/F 35.6%, IIV V/F 28.0%, correlation 0.834
    etalka ~ 1.664100 # Yao 2018 Table 4: IIV ka 129%

    propSd <- 0.200; label("Proportional residual error on Cc (fraction)") # Yao 2018 Table 4: proportional residual error CV 20.0%
    addSd <- fixed(0.00289); label("Additive residual error on Cc (ug/mL)") # Yao 2018 Table 4 and Methods: additive error fixed at 0.00289 ug/mL

    # PD: Iwaki 2021 Table 1 (final PK/PD model) and the Supplemental NONMEM
    # control stream ($THETA 1-8, $OMEGA 1-5).
    ltau <- log(64.7); label("Precision of the beta distribution of the transformed PPPASI score, PhiBeta (unitless)") # Iwaki 2021 Table 1: phi = 64.7 (RSE 1.67%); THETA(1)
    logitrbase <- logit(0.356); label("Logit of the typical baseline transformed PPPASI score mu_b (logit scale; mu_b = 0.356 on the 0-1 scale)") # Iwaki 2021 Table 1: mu_beta = 0.356 (RSE 3.68%); THETA(2), entered as LOG(THETA(2)/(1-THETA(2)))
    lpmax <- log(1.55); label("Maximum placebo effect Pmax on the logit of the beta mean (logit units)") # Iwaki 2021 Table 1: Fp = 1.55 (RSE 6.71%); THETA(3) PBMX
    lkpbo <- log(0.004); label("Placebo effect onset rate constant kp (1/day)") # Iwaki 2021 Table 1: kp = 0.004 1/day (RSE 19.4%); THETA(4) PBRATE
    lkout <- log(0.00139); label("Drug effect onset rate constant kout = kin of the latent variable R (1/day)") # Iwaki 2021 Table 1: kout = 0.00139 1/day (RSE 37.5%); THETA(5) KOUT
    lic50 <- fixed(log(0.000001)); label("Guselkumab concentration giving half-maximal inhibition of the production of R, IC50 (ug/mL)") # Iwaki 2021 Table 1: IC50 = 0.000001 ug/mL fixed; THETA(6) EC50 0.000001 FIX
    lemax <- log(2.97); label("Maximum drug effect DE on the logit of the beta mean (logit units)") # Iwaki 2021 Table 1: DE = 2.97 (RSE 31.0%); THETA(7) DMX
    dp_slope <- 0.000944; label("Linear disease-progression slope kdspg on the logit of the beta mean (logit units/day)") # Iwaki 2021 Table 1: kdspg = 0.000944 1/day (RSE 33.9%); THETA(8) DSPG

    # Iwaki 2021 Table 1 'BSV Variance'; $OMEGA is diagonal. The baseline eta
    # is normal on the logit scale and the kdspg eta is proportional,
    # DSPG = THETA(8) * (1 + ETA(5)); the others are log-normal.
    etalogitrbase ~ 0.507 # Iwaki 2021 Table 1: BSV variance mu_beta 0.507 (RSE 12.0%); ETA(1)
    etalkpbo ~ 1.69 # Iwaki 2021 Table 1: BSV variance kp 1.69 (RSE 21.9%); ETA(2)
    etalkout ~ 1.04 # Iwaki 2021 Table 1: BSV variance kout 1.04 (RSE 47.1%); ETA(3)
    etalemax ~ 0.320 # Iwaki 2021 Table 1: BSV variance DE 0.320 (RSE 53.8%); ETA(4)
    etadp_slope ~ 10.3 # Iwaki 2021 Table 1: BSV variance kdspg 10.3 (RSE 71.7%); ETA(5)
  })
  model({
    # PK (Yao 2018 Table 4 footnotes b and c). Race factor uses 1 - RACE_WHITE
    # because Yao's RACE = 1 marks non-White subjects.
    cl <- exp(lcl + etalcl) * (WT / 87.1)^e_wt_cl * e_diab_cl^DIS_DIAB * e_race_cl^(1 - RACE_WHITE)
    vc <- exp(lvc + etalvc) * (WT / 87.1)^e_wt_vc
    ka <- exp(lka + etalka)
    kel <- cl / vc

    # PD individual parameters (Supplemental $PK).
    tau <- exp(ltau)
    logit_mub <- logitrbase + etalogitrbase
    pmax <- exp(lpmax)
    kpbo <- exp(lkpbo + etalkpbo)
    kout <- exp(lkout + etalkout)
    ic50 <- exp(lic50)
    emax <- exp(lemax + etalemax)
    kdspg <- dp_slope * (1 + etadp_slope)

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # dose mg / volume L = mg/L = ug/mL
    Cc <- central / vc

    # Type I indirect response with kin = kout so that R(0) = 1 (Iwaki 2021
    # Methods; Supplemental $DES: DADT(3) = KOUT * (1 - IH - A(3)), A_0(3) = 1).
    inh <- Cc / (ic50 + Cc)
    d/dt(latent_pppasi) <- kout * (1 - inh - latent_pppasi)
    latent_pppasi(0) <- 1

    # Mean of the beta distribution on the 0-1 PPPASI scale (Supplemental
    # $ERROR: MuBlgt = BSlgt - PBMX * (1 - EXP(-PBRATE * TIME)) - DMX * (1 - F)
    # + DSPG * TIME, with F = A(3) on PPPASI records). The placebo and drug
    # terms are subtracted; Iwaki 2021 prints them with plus signs (see the
    # vignette Errata). t is time since the first dose (days).
    placebo_eff <- pmax * (1 - exp(-kpbo * t))
    drug_eff <- emax * (1 - latent_pppasi)
    mu_pppasi <- expit(logit_mub - placebo_eff - drug_eff + kdspg * t)
    alpha_pppasi <- mu_pppasi * tau
    beta_pppasi <- (1 - mu_pppasi) * tau

    # Ordered-categorical likelihood of the observed score DV (0-72 in 0.1
    # steps, k = 10 * DV = 0..720): prob(k) = B((k + 1)/721) - B(k/721), where
    # B is the beta CDF; ibeta(a, b, x) is the regularized incomplete beta
    # function, equal to pbeta(x, a, b). Floored at 1e-30 as in the
    # Supplemental $ERROR.
    p_pppasi <- ibeta(alpha_pppasi, beta_pppasi, (10 * DV + 1) / 721) -
      ibeta(alpha_pppasi, beta_pppasi, 10 * DV / 721)
    ll(pppasi) ~ log(max(p_pppasi, 1e-30))

    Cc ~ add(addSd) + prop(propSd)
  })
}
