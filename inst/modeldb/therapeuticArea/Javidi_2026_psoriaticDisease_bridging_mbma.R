Javidi_2026_psoriaticDisease_bridging_mbma <- function() {
  description <- paste0(
    "MBMA. Cross-indication bridging model ('Model 1') of Javidi 2026: a ",
    "longitudinal, parametric-placebo, dose-response model-based meta-analysis ",
    "that jointly describes three study-arm responder rates over time -- ",
    "PASI75 in plaque psoriasis (prob_pasi75), PASI75 in psoriatic arthritis ",
    "(prob_pasi75_psa) and ACR20 in psoriatic arthritis (prob_acr20) -- so a ",
    "treatment studied in one indication can be projected into the other. ",
    "Fitted with nlme to aggregate data from 66 randomised controlled trials ",
    "(34 psoriasis, 32 psoriatic arthritis; 190 arms, 30,450 participants) of ",
    "11 drugs plus placebo. On the logit scale each output is a placebo term ",
    "that rises exponentially to an indication-specific plateau plus, on active ",
    "arms, a drug term that rises exponentially at a drug- and endpoint-specific ",
    "rate to a maximum effect. The maximum is an Emax function of normalised ",
    "dose for guselkumab, tofacitinib, secukinumab and deucravacitinib and a ",
    "single dose-independent shift for adalimumab, apremilast, bimekizumab, ",
    "etanercept, ixekizumab, risankizumab and ustekinumab; indication, endpoint ",
    "and study phase (phase 3 vs phase 2) enter as additive shifts. Correlated ",
    "between-trial random effects act on the placebo intercept and on the ",
    "treatment maximum of each endpoint-indication combination. Residual ",
    "variance is the binomial p(1 - p)/N of an arm proportion (sigma fixed to ",
    "1), so N_ARM is required for stochastic simulation. All three outputs are ",
    "computed on every record; read the one that matches the trial's indication ",
    "and endpoint. Time is in weeks. Simulation scope is study-arm-mean ",
    "responder trajectories, NOT individual patients. The companion ",
    "cross-endpoint models that map ACR20 to ACR50 / ACR70 are ",
    "modellib('Javidi_2026_psa_acr50_acr70_wk16_mbma') and ",
    "modellib('Javidi_2026_psa_acr50_acr70_wk24_mbma')."
  )

  reference <- paste(
    "Javidi M, Alexander GE, Kosoglou T, Vermeulen A, Beutler A, Shagroni T,",
    "Valiathan C.",
    "MBMA Bridging Models as a Tool for Exploration of Clinical Endpoints in",
    "Unstudied Indications.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15(5):726-740.",
    "doi:10.1002/psp4.70207. PMC13142780.",
    "Structural model: main-text Equations 1-4 and Supporting Information",
    "Equations S.2.1, S.3.1 and S.3.2.1-S.3.2.4, together with the authors'",
    "nlme code for Model 1 in the Supporting Information (Data S1).",
    "Parameter estimates: Supporting Information Table S1 (fixed effects) and",
    "Table S2 (random-effect variance-covariance matrix).",
    sep = " "
  )

  vignette <- "Javidi_2026_psoriaticDisease_mbma"

  units <- list(
    time = "week (weeks since randomisation; the placebo and drug onset rate constants are in week^-1)",
    dosing = "n/a (this MBMA consumes no rxode2 dose events; treatment enters through the TRT_<INN> arm indicators and, for the four dose-response drugs, the DOSE_<DRUG>_REL normalised-dose columns)",
    concentration = "probability/arm (each output is the STUDY-ARM probability that a patient meets the named responder criterion, on a 0-1 scale; it is NOT a drug concentration. The slash satisfies checkModelConventions unit parsing.)"
  )

  covariateData <- list(
    TRT_ADALIMUMAB = list(
      description = "Study-arm indicator: 1 = the arm received adalimumab (any dose), 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Adalimumab has no dose-response term (too few doses in the database): its effect is a single dose-independent shift x.drg_adalimumab on the treatment maximum. Arm-level, not patient-level.",
      source_name = "adalimumab (Javidi 2026 Data S1 Model 1 code: 'binary value which is 1 for rows of data of the treatment')"
    ),
    TRT_APREMILAST = list(
      description = "Study-arm indicator: 1 = the arm received apremilast (any dose), 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "No dose-response term; single dose-independent shift x.drg_apremilast. Arm-level.",
      source_name = "apremilast (Javidi 2026 Data S1 Model 1 code)"
    ),
    TRT_BIMEKIZUMAB = list(
      description = "Study-arm indicator: 1 = the arm received bimekizumab (any dose), 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "No dose-response term; single dose-independent shift x.drg_bimekizumab. Arm-level.",
      source_name = "bimekizumab (Javidi 2026 Data S1 Model 1 code)"
    ),
    TRT_DEUCRAVACITINIB = list(
      description = "Study-arm indicator: 1 = the arm received deucravacitinib (any dose), 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Dose-response drug: selects the deucravacitinib Emax, ED50 and onset rate; the dose itself is supplied in DOSE_DEUCRAVACITINIB_REL. The psoriatic arthritis deucravacitinib trial (IM011-084) was held out of the fit for external validation (Javidi 2026 Methods 'Model Validations'), so the PsA deucravacitinib predictions are genuine projections from psoriasis data. Arm-level.",
      source_name = "deucravacitinib (Javidi 2026 Data S1 Model 1 code)"
    ),
    TRT_ETANERCEPT = list(
      description = "Study-arm indicator: 1 = the arm received etanercept (any dose), 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "No dose-response term; single dose-independent shift x.drg_etanercept. Arm-level.",
      source_name = "etanercept (Javidi 2026 Data S1 Model 1 code)"
    ),
    TRT_GUSELKUMAB = list(
      description = "Study-arm indicator: 1 = the arm received guselkumab (any dose), 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Dose-response drug: selects the guselkumab Emax, ED50 and onset rate; the dose itself is supplied in DOSE_GUSELKUMAB_REL. Arm-level.",
      source_name = "guselkumab (Javidi 2026 Data S1 Model 1 code)"
    ),
    TRT_IXEKIZUMAB = list(
      description = "Study-arm indicator: 1 = the arm received ixekizumab (any dose), 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "No dose-response term; single dose-independent shift x.drg_ixekizumab. Arm-level.",
      source_name = "ixekizumab (Javidi 2026 Data S1 Model 1 code)"
    ),
    TRT_RISANKIZUMAB = list(
      description = "Study-arm indicator: 1 = the arm received risankizumab (any dose), 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "No dose-response term; single dose-independent shift x.drg_risankizumab. Arm-level.",
      source_name = "risankizumab (Javidi 2026 Data S1 Model 1 code)"
    ),
    TRT_SECUKINUMAB = list(
      description = "Study-arm indicator: 1 = the arm received secukinumab (any dose), 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Dose-response drug: selects the secukinumab Emax, ED50 and onset rate; the dose itself is supplied in DOSE_SECUKINUMAB_REL. Arm-level.",
      source_name = "secukinumab (Javidi 2026 Data S1 Model 1 code)"
    ),
    TRT_TOFACITINIB = list(
      description = "Study-arm indicator: 1 = the arm received tofacitinib (any dose), 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Dose-response drug: selects the tofacitinib Emax and ED50; the dose itself is supplied in DOSE_TOFACITINIB_REL. Tofacitinib is the reference level of the drug onset-rate formula (it is the one drug without a k.drug term in the Data S1 code and Table S1), so its log onset rate is 0 in psoriasis (rate 1 per week) plus the indication / endpoint shifts. Arm-level.",
      source_name = "tofacitinib (Javidi 2026 Data S1 Model 1 code)"
    ),
    TRT_USTEKINUMAB = list(
      description = "Study-arm indicator: 1 = the arm received ustekinumab (any dose), 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "No dose-response term; single dose-independent shift x.drg_ustekinumab. Arm-level.",
      source_name = "ustekinumab (Javidi 2026 Data S1 Model 1 code)"
    ),
    DOSE_GUSELKUMAB_REL = list(
      description = "Guselkumab dose intensity of the arm divided by that of the reference regimen 100 mg every 8 weeks (after induction); 0 on non-guselkumab arms.",
      units = "(ratio)",
      type = "continuous",
      reference_category = NULL,
      notes = "The paper normalises dose 'by the approved dose or the most frequent tested dose for each compound' without printing the reference regimens. The maintainers recovered them from the Model 1 typical curves of Supporting Information Figure S6 (ACR20, potential phase 3 study): 100 mg q8w reproduces the plotted plateau with 1 and 100 mg q4w with 2. A regimen's value is its average maintenance dose per unit time relative to 100 mg q8w (e.g. 50 mg q8w = 0.5).",
      source_name = "reldose (Javidi 2026 Data S1 Model 1 code: 'normalized dose')"
    ),
    DOSE_TOFACITINIB_REL = list(
      description = "Tofacitinib daily dose of the arm divided by the reference regimen 5 mg twice daily; 0 on non-tofacitinib arms.",
      units = "(ratio)",
      type = "continuous",
      reference_category = NULL,
      notes = "Reference regimen recovered by the maintainers from Supporting Information Figure S6: 5 mg BID = 1 and 10 mg BID = 2 reproduce the plotted ACR20 plateaus. See DOSE_GUSELKUMAB_REL.",
      source_name = "reldose (Javidi 2026 Data S1 Model 1 code)"
    ),
    DOSE_SECUKINUMAB_REL = list(
      description = "Secukinumab maintenance dose of the arm divided by the reference regimen 300 mg every 4 weeks; 0 on non-secukinumab arms.",
      units = "(ratio)",
      type = "continuous",
      reference_category = NULL,
      notes = "Reference regimen recovered by the maintainers from Supporting Information Figure S6: 300 mg = 1, 150 mg = 0.5 and 75 mg = 0.25 reproduce the three plotted ACR20 plateaus; the weekly loading doses and the 10 mg/kg IV loading of FUTURE 1 do not change the value (the 10 mg/kg + 150 mg curve equals the 150 mg q4w curve). See DOSE_GUSELKUMAB_REL.",
      source_name = "reldose (Javidi 2026 Data S1 Model 1 code)"
    ),
    DOSE_DEUCRAVACITINIB_REL = list(
      description = "Deucravacitinib daily dose of the arm divided by the reference regimen 6 mg once daily; 0 on non-deucravacitinib arms.",
      units = "(ratio)",
      type = "continuous",
      reference_category = NULL,
      notes = "Reference regimen 6 mg QD (the approved psoriasis dose; Javidi 2026 Results 'Applications (1)'). The value 1 for 6 mg QD and 2 for 12 mg QD reproduces the Table S4 projections for both IM011-084 deucravacitinib arms (ACR20 51.4% and 55.2%; PASI75 37.3% and 44.2%). See DOSE_GUSELKUMAB_REL.",
      source_name = "reldose (Javidi 2026 Data S1 Model 1 code)"
    ),
    STUDY_PHASE3 = list(
      description = "Study-level indicator: 1 = the arm belongs to a phase 3 trial, 0 = a phase 2 trial.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Javidi 2026 codes study phase as two arm-type-specific indicators, 'phaseP' (phase 3 AND placebo arm; coefficient p.eff.ph3 on the placebo intercept) and 'phase' (phase 3 AND treatment arm; coefficient drg.eff.ph3 on the treatment maximum, indication-specific). Both are derived inside model() from STUDY_PHASE3 and the TRT_<INN> indicators. The database also holds phase 2a / 2b, 3b and 4 trials (Table S5); the paper does not say how those were binned, and the natural reading (2a/2b -> 0; 3b/4 -> 1) is an assumption. The paper's simulations ('Phase 2' / 'Phase 3' scenarios of Table 1 and Figure 5) set this to 0 or 1.",
      source_name = "phase / phaseP (Javidi 2026 Data S1 Model 1 code)"
    ),
    N_ARM = list(
      description = "Number of participants contributing to the arm's responder proportion at the record's time point.",
      units = "participants",
      type = "count",
      reference_category = NULL,
      notes = "Study-design quantity. Sets the residual variance p(1 - p)/N_ARM of Javidi 2026 Equation 1.1 (nlme varPower(0.5, form = ~fitted(.)*(1 - fitted(.))/n, fixed = 0.5) with sigma fixed to 1). Used only by the residual: typical-value (zeroRe) simulations do not depend on it, but any stochastic simulation that keeps the residual must supply a strictly positive value.",
      source_name = "n (Javidi 2026 Data S1 Model 1 code: 'sample size'); N_ijkt (Equation 1.1)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 30450L,
    n_studies = 66L,
    n_arms = 190L,
    age_range = "not reported at arm level",
    weight_range = "not reported at arm level",
    sex_female_pct = NA_real_,
    race_ethnicity = "not reported at arm level",
    disease_state = "adults with moderate-to-severe plaque psoriasis (34 trials) or active psoriatic arthritis (32 trials) in randomised, parallel-group, placebo- or active-controlled phase 2 / phase 3 trials",
    dose_range = "approved and dose-ranging regimens of adalimumab, apremilast, bimekizumab, deucravacitinib, etanercept, guselkumab, ixekizumab, risankizumab, secukinumab, tofacitinib and ustekinumab (Supporting Information Table S5)",
    timepoints = "longitudinal arm responder rates up to 52 weeks; data after re-randomisation, a regimen change or substantial drop-out were excluded",
    regions = "international; trials published 2000-2023 from a vendor-curated database built following the Cochrane Handbook and PRISMA (PubMed, clinicaltrials.gov, company registries, conference material, FDA and EMA reviews)",
    notes = "MBMA at the study-arm level. The deucravacitinib psoriatic arthritis phase 2 trial IM011-084 was excluded from the fit and used for external validation (Javidi 2026 Methods 'Model Validations'). Model 1 uses only compounds studied in both indications (12 treatments including placebo; Results 'Data'). Random effects are between-trial (trial-indication-endpoint level), so the model must not be used to simulate individual patients. Fitted with nlme 3.1-160 in R 4.4.1."
  )

  ini({
    # ============================================================
    # Placebo component, logit scale (Javidi 2026 Equation 2 and
    # Supporting Information Equation S.2.1):
    #   pbo = p.eff + plsp * (1 - exp(-k.pbo * time)) + eta
    # with p.eff = beta1*p.eff.pso + beta2*p.eff.psa
    #            + beta3*p.eff.psa.acr20 + beta4*p.eff.ph3.
    # In the Data S1 nlme code the rate enters as exp(pk), so the
    # tabulated k.pbo values are LOG rate constants (week^-1); they
    # are carried as lkpbo_*. plsp and k.pbo are indication-level
    # (formula '~ INDI - 1'): Table S1 repeats the PsA value for the
    # PASI75.PsA and ACR20.PsA rows.
    # All values Javidi 2026 Supporting Information Table S1.
    # ============================================================
    peff_pso <- -7.08
    label("Placebo logit intercept at time 0, psoriasis (paper: p.eff.pso; unitless log-odds)") # Table S1 p.eff.pso = -7.08 (SE 0.15)
    peff_psa <- -4.09
    label("Placebo logit intercept at time 0, psoriatic arthritis (paper: p.eff.psa; unitless log-odds)") # Table S1 p.eff.psa = -4.09 (SE 0.13)
    peff_psa_acr20 <- 0.85
    label("Additional placebo logit intercept for the ACR20 endpoint in psoriatic arthritis (paper: p.eff.psa.acr20; unitless log-odds)") # Table S1 p.eff.psa.acr20 = 0.85 (SE 0.1)
    peff_ph3 <- 0.06
    label("Placebo logit shift on placebo arms of phase 3 trials (paper: p.eff.ph3; unitless log-odds)") # Table S1 p.eff.ph3 = 0.06 (SE 0.09)
    plsp_pso <- 4.78
    label("Maximum placebo-effect increment on the logit scale, psoriasis (paper: plsp PASI75,PsO; unitless log-odds)") # Table S1 plsp PASI75,PsO = 4.78 (SE 0.14)
    plsp_psa <- 2.03
    label("Maximum placebo-effect increment on the logit scale, psoriatic arthritis, shared by PASI75 and ACR20 (paper: plsp PASI75,PsA = plsp ACR20,PsA; unitless log-odds)") # Table S1 plsp PASI75,PsA = plsp ACR20,PsA = 2.03 (SE 0.1)
    lkpbo_pso <- -1.45
    label("Log rate constant of placebo-effect onset, psoriasis (paper: k.pbo PASI75,PsO; log week^-1; back-transform 0.235 /week)") # Table S1 k.pbo PASI75,PsO = -1.45 (SE 0.03)
    lkpbo_psa <- -1.29
    label("Log rate constant of placebo-effect onset, psoriatic arthritis, shared by PASI75 and ACR20 (paper: k.pbo PASI75,PsA = k.pbo ACR20,PsA; log week^-1; back-transform 0.275 /week)") # Table S1 k.pbo PASI75,PsA = k.pbo ACR20,PsA = -1.29 (SE 0.05)

    # ============================================================
    # Treatment maximum effect, logit scale (Javidi 2026 Equation 3.2
    # and Equations S.3.2.1-S.3.2.4):
    #   drg.eff = c.eff + EMAX*dose/(dose + ED50)  [dose-response drugs]
    #           + X                                 [other drugs]
    #           + gamma
    # c.eff = beta2*drg.eff.psa + beta3*drg.eff.psa.acr20
    #       + beta5*drg.eff.ph3 (phase 3 treatment arms; indication-
    #       level per the code term 'phase:INDI', so Table S1 repeats
    #       the PsA value for the PASI75.PsA and ACR20.PsA rows).
    # ============================================================
    deff_psa <- 1.97
    label("Shift in treatment maximum effect for psoriatic arthritis (paper: drg.eff.psa; unitless log-odds)") # Table S1 drg.eff.psa = 1.97 (SE 0.34)
    deff_psa_acr20 <- -0.83
    label("Additional shift in treatment maximum effect for the ACR20 endpoint in psoriatic arthritis (paper: drg.eff.psa.acr20; unitless log-odds)") # Table S1 drg.eff.psa.acr20 = -0.83 (SE 0.32)
    deff_ph3_pso <- 0.33
    label("Shift in treatment maximum effect on treatment arms of phase 3 psoriasis trials (paper: drg.eff.ph3 PASI75,PsO; unitless log-odds)") # Table S1 drg.eff.ph3 PASI75,PsO = 0.33 (SE 0.14)
    deff_ph3_psa <- -0.32
    label("Shift in treatment maximum effect on treatment arms of phase 3 psoriatic arthritis trials, shared by PASI75 and ACR20 (paper: drg.eff.ph3 PASI75,PsA = drg.eff.ph3 ACR20,PsA; unitless log-odds)") # Table S1 drg.eff.ph3 ACR20,PsA = drg.eff.ph3 PASI75,PsA = -0.32 (SE 0.24)

    # Drugs without a dose-response model: X = beta2*x.psa +
    # beta3*x.psa.acr20 + x.drg (Equation S.3.2.4). The x.psa terms
    # apply only to these drugs (code indicators PSAsingle and
    # ACR20single).
    xeff_psa <- -2
    label("Shift in the dose-independent drug effect for psoriatic arthritis, drugs without a dose-response model (paper: x.psa; unitless log-odds)") # Table S1 x.psa = -2 (SE 0.26)
    xeff_psa_acr20 <- -0.76
    label("Additional shift in the dose-independent drug effect for ACR20 in psoriatic arthritis (paper: x.psa.acr20; unitless log-odds)") # Table S1 x.psa.acr20 = -0.76 (SE 0.26)
    xeff_adalimumab <- 3.01
    label("Dose-independent adalimumab maximum effect, psoriasis scale (paper: x.drg adalimumab; unitless log-odds)") # Table S1 x.drg adalimumab = 3.01 (SE 0.16)
    xeff_apremilast <- 1.61
    label("Dose-independent apremilast maximum effect (paper: x.drg apremilast; unitless log-odds)") # Table S1 x.drg apremilast = 1.61 (SE 0.15)
    xeff_bimekizumab <- 4.51
    label("Dose-independent bimekizumab maximum effect (paper: x.drg bimekizumab; unitless log-odds)") # Table S1 x.drg bimekizumab = 4.51 (SE 0.18)
    xeff_etanercept <- 2.51
    label("Dose-independent etanercept maximum effect (paper: x.drg etanercept; unitless log-odds)") # Table S1 x.drg etanercept = 2.51 (SE 0.16)
    xeff_ixekizumab <- 3.84
    label("Dose-independent ixekizumab maximum effect (paper: x.drg ixekizumab; unitless log-odds)") # Table S1 x.drg ixekizumab = 3.84 (SE 0.16)
    xeff_risankizumab <- 4.58
    label("Dose-independent risankizumab maximum effect (paper: x.drg risankizumab; unitless log-odds)") # Table S1 x.drg risankizumab = 4.58 (SE 0.17)
    xeff_ustekinumab <- 3.06
    label("Dose-independent ustekinumab maximum effect (paper: x.drg ustekinumab; unitless log-odds)") # Table S1 x.drg ustekinumab = 3.06 (SE 0.17)

    # Dose-response drugs (set Omega = guselkumab, tofacitinib,
    # secukinumab, deucravacitinib): EMAX = beta3*emax.psa.acr20 +
    # emax.drg (S.3.2.2); ED50 = beta2*ed50.psa + beta3*ed50.psa.acr20
    # + ed50.drg (S.3.2.3). The code writes the ED50 as exp(ed50), so
    # the ED50 rows of Table S1 are LOG ED50 on the normalised-dose
    # scale and are carried as led50_*.
    emax_psa_acr20 <- -2.84
    label("Shift in Emax for ACR20 in psoriatic arthritis, dose-response drugs (paper: emax.psa.acr20; unitless log-odds)") # Table S1 emax.psa.acr20 = -2.84 (SE 0.56)
    emax_guselkumab <- 4.18
    label("Guselkumab Emax (paper: emax.drg guselkumab; unitless log-odds)") # Table S1 emax.drg guselkumab = 4.18 (SE 0.16)
    emax_tofacitinib <- 3.83
    label("Tofacitinib Emax (paper: emax.drg tofacitinib; unitless log-odds)") # Table S1 emax.drg tofacitinib = 3.83 (SE 0.2)
    emax_secukinumab <- 4.66
    label("Secukinumab Emax (paper: emax.drg secukinumab; unitless log-odds)") # Table S1 emax.drg secukinumab = 4.66 (SE 0.21)
    emax_deucravacitinib <- 4.64
    label("Deucravacitinib Emax (paper: emax.drg deucravacitinib; unitless log-odds)") # Table S1 emax.drg deucravacitinib = 4.64 (SE 0.76)
    led50_psa <- 1.76
    label("Shift in log ED50 for psoriatic arthritis, dose-response drugs (paper: ed50.psa; log normalised dose)") # Table S1 ed50.psa = 1.76 (SE 0.3)
    led50_psa_acr20 <- 0.51
    label("Additional shift in log ED50 for ACR20 in psoriatic arthritis (paper: ed50.psa.acr20; log normalised dose)") # Table S1 ed50.psa.acr20 = 0.51 (SE 0.57)
    led50_guselkumab <- -3.27
    label("Guselkumab log ED50 (paper: ed50.drg guselkumab; log normalised dose; back-transform 0.038)") # Table S1 ed50.drg guselkumab = -3.27 (SE 0.25)
    led50_tofacitinib <- 0.06
    label("Tofacitinib log ED50 (paper: ed50.drg tofacitinib; log normalised dose; back-transform 1.06)") # Table S1 ed50.drg tofacitinib = 0.06 (SE 0.18)
    led50_secukinumab <- -1.28
    label("Secukinumab log ED50 (paper: ed50.drg secukinumab; log normalised dose; back-transform 0.278)") # Table S1 ed50.drg secukinumab = -1.28 (SE 0.16)
    led50_deucravacitinib <- -0.14
    label("Deucravacitinib log ED50 (paper: ed50.drg deucravacitinib; log normalised dose; back-transform 0.869)") # Table S1 ed50.drg deucravacitinib = -0.14 (SE 0.38)

    # ============================================================
    # Treatment onset rate (Equations 3 and S.3.1):
    #   tr.eff = drg.eff * (1 - exp(-k.drg * time)),
    #   k.drg = beta2*k.drg.psa + beta3*k.drg.psa.acr20 + k.drug_j
    # entering the code as exp(k): log rate constants, week^-1.
    # Tofacitinib is the reference level of the '~ ... - 1' formula
    # and has no k.drug term (log rate 0 in psoriasis).
    # ============================================================
    lkdrug_psa <- -1.72
    label("Shift in log onset rate of the treatment effect for psoriatic arthritis (paper: k.drg.psa; log week^-1)") # Table S1 k.drg.psa = -1.72 (SE 0.1)
    lkdrug_psa_acr20 <- 1.94
    label("Additional shift in log onset rate for ACR20 in psoriatic arthritis (paper: k.drg.psa.acr20; log week^-1)") # Table S1 k.drg.psa.acr20 = 1.94 (SE 0.14)
    lkdrug_adalimumab <- -0.52
    label("Adalimumab log onset rate of the treatment effect (paper: k.drug adalimumab; log week^-1)") # Table S1 k.drug adalimumab = -0.52 (SE 0.11)
    lkdrug_apremilast <- -1.28
    label("Apremilast log onset rate of the treatment effect (paper: k.drug apremilast; log week^-1)") # Table S1 k.drug apremilast = -1.28 (SE 0.17)
    lkdrug_bimekizumab <- 0.2
    label("Bimekizumab log onset rate of the treatment effect (paper: k.drug bimekizumab; log week^-1)") # Table S1 k.drug bimekizumab = 0.2 (SE 0.11)
    lkdrug_deucravacitinib <- -1.21
    label("Deucravacitinib log onset rate of the treatment effect (paper: k.drug deucravacitinib; log week^-1)") # Table S1 k.drug deucravacitinib = -1.21 (SE 0.16)
    lkdrug_etanercept <- -1.45
    label("Etanercept log onset rate of the treatment effect (paper: k.drug etanercept; log week^-1)") # Table S1 k.drug etanercept = -1.45 (SE 0.11)
    lkdrug_guselkumab <- -1.34
    label("Guselkumab log onset rate of the treatment effect (paper: k.drug guselkumab; log week^-1)") # Table S1 k.drug guselkumab = -1.34 (SE 0.05)
    lkdrug_ixekizumab <- 0.31
    label("Ixekizumab log onset rate of the treatment effect (paper: k.drug ixekizumab; log week^-1)") # Table S1 k.drug ixekizumab = 0.31 (SE 0.14)
    lkdrug_risankizumab <- -1.38
    label("Risankizumab log onset rate of the treatment effect (paper: k.drug risankizumab; log week^-1)") # Table S1 k.drug risankizumab = -1.38 (SE 0.07)
    lkdrug_secukinumab <- -0.56
    label("Secukinumab log onset rate of the treatment effect (paper: k.drug secukinumab; log week^-1)") # Table S1 k.drug secukinumab = -0.56 (SE 0.06)
    lkdrug_ustekinumab <- -1.21
    label("Ustekinumab log onset rate of the treatment effect (paper: k.drug ustekinumab; log week^-1)") # Table S1 k.drug ustekinumab = -1.21 (SE 0.09)

    # ============================================================
    # Between-trial random effects (Equation 4; Table S2 (A.)).
    # pdBlocked: one block for the placebo-intercept etas and one for
    # the treatment-maximum gammas; no covariance across the blocks.
    # Each block has one column per endpoint-indication indicator
    # (INDIPSO, INDIC, INDIC:ENDACR20), so ACR20 in PsA carries the
    # PsA effect PLUS its own ACR20 effect. Values are variances.
    # Table S2 (A.) prints 0 for the eta PsO-PsA covariances.
    # ============================================================
    eta_study_peff_pso ~ 0.067 # Table S2 (A.) omega eta,PsO = 0.067
    eta_study_peff_psa + eta_study_peff_psa_acr20 ~ c(0.011, -0.010, 0.064) # Table S2 (A.) omega eta,PsA = 0.011; cov(eta PsA, eta ACR20.PsA) = -0.010; omega eta,ACR20.PsA = 0.064
    eta_study_deff_pso + eta_study_deff_psa + eta_study_deff_psa_acr20 ~ c(0.051, 0.003, 0.110, 0.014, -0.013, 0.628) # Table S2 (A.) omega gamma: PsO 0.051; PsO-PsA 0.003; PsA 0.110; PsO-ACR20.PsA 0.014; PsA-ACR20.PsA -0.013; ACR20.PsA 0.628

    # ============================================================
    # Residual (Equation 1.1): Var(eps) = p(1 - p)/N, implemented in
    # Data S1 as weights = varPower(0.5, form = ~fitted*(1 - fitted)/n,
    # fixed = 0.5) with nlmeControl(sigma = 1). The unit sigma is
    # carried per output as fixed(1) and scaled in model() by the
    # binomial standard error.
    # ============================================================
    addSd_prob_pasi75 <- fixed(1)
    label("Residual SD multiplier for PASI75 in psoriasis; the record SD is sqrt(p(1 - p)/N_ARM) times this (unitless; the source sets sigma = 1)") # Data S1 Model 1 code nlmeControl(sigma = 1)
    addSd_prob_pasi75_psa <- fixed(1)
    label("Residual SD multiplier for PASI75 in psoriatic arthritis; the record SD is sqrt(p(1 - p)/N_ARM) times this (unitless; the source sets sigma = 1)") # Data S1 Model 1 code nlmeControl(sigma = 1)
    addSd_prob_acr20 <- fixed(1)
    label("Residual SD multiplier for ACR20 in psoriatic arthritis; the record SD is sqrt(p(1 - p)/N_ARM) times this (unitless; the source sets sigma = 1)") # Data S1 Model 1 code nlmeControl(sigma = 1)
  })

  model({
    # No ODEs and no dose events: every output is an algebraic function
    # of time (weeks) and the arm covariates.

    # ---- Arm type -----------------------------------------------------
    # Drugs without a dose-response model (code indicator PSAsingle /
    # ACR20single when combined with PsA).
    trtSingle <- TRT_ADALIMUMAB + TRT_APREMILAST + TRT_BIMEKIZUMAB +
      TRT_ETANERCEPT + TRT_IXEKIZUMAB + TRT_RISANKIZUMAB + TRT_USTEKINUMAB
    # Dose-response drugs (code indicator PSAEMAX / ACR20EMAX).
    trtEmax <- TRT_GUSELKUMAB + TRT_TOFACITINIB + TRT_SECUKINUMAB +
      TRT_DEUCRAVACITINIB
    active <- 0
    if (trtSingle + trtEmax > 0) {
      active <- 1
    }
    # 'phaseP' (phase 3 placebo arm) and 'phase' (phase 3 treatment arm).
    ph3Pbo <- STUDY_PHASE3 * (1 - active)
    ph3Trt <- STUDY_PHASE3 * active

    # ---- Between-trial random effects ---------------------------------
    # One fixed effect + its eta per line keeps each term mu-referenced.
    peffPso <- peff_pso + eta_study_peff_pso
    peffPsa <- peff_psa + eta_study_peff_psa
    peffAcr20 <- peff_psa_acr20 + eta_study_peff_psa_acr20
    deffPso <- eta_study_deff_pso
    deffPsa <- deff_psa + eta_study_deff_psa
    deffAcr20 <- deff_psa_acr20 + eta_study_deff_psa_acr20

    # ---- Placebo component (Equation 2, S.2.1) ------------------------
    onsetPboPso <- 1 - exp(-exp(lkpbo_pso) * time)
    onsetPboPsa <- 1 - exp(-exp(lkpbo_psa) * time)
    pboPasi75Pso <- peffPso + peff_ph3 * ph3Pbo + plsp_pso * onsetPboPso
    pboPasi75Psa <- peffPsa + peff_ph3 * ph3Pbo + plsp_psa * onsetPboPsa
    pboAcr20Psa <- peffPsa + peffAcr20 + peff_ph3 * ph3Pbo +
      plsp_psa * onsetPboPsa

    # ---- Drug-specific terms ------------------------------------------
    xDrg <- xeff_adalimumab * TRT_ADALIMUMAB +
      xeff_apremilast * TRT_APREMILAST +
      xeff_bimekizumab * TRT_BIMEKIZUMAB +
      xeff_etanercept * TRT_ETANERCEPT +
      xeff_ixekizumab * TRT_IXEKIZUMAB +
      xeff_risankizumab * TRT_RISANKIZUMAB +
      xeff_ustekinumab * TRT_USTEKINUMAB
    emaxDrg <- emax_guselkumab * TRT_GUSELKUMAB +
      emax_tofacitinib * TRT_TOFACITINIB +
      emax_secukinumab * TRT_SECUKINUMAB +
      emax_deucravacitinib * TRT_DEUCRAVACITINIB
    led50Drg <- led50_guselkumab * TRT_GUSELKUMAB +
      led50_tofacitinib * TRT_TOFACITINIB +
      led50_secukinumab * TRT_SECUKINUMAB +
      led50_deucravacitinib * TRT_DEUCRAVACITINIB
    relDose <- DOSE_GUSELKUMAB_REL * TRT_GUSELKUMAB +
      DOSE_TOFACITINIB_REL * TRT_TOFACITINIB +
      DOSE_SECUKINUMAB_REL * TRT_SECUKINUMAB +
      DOSE_DEUCRAVACITINIB_REL * TRT_DEUCRAVACITINIB
    # Tofacitinib is the reference level of the onset-rate formula
    # (no k.drug term), so it contributes 0 here.
    lkDrg <- lkdrug_adalimumab * TRT_ADALIMUMAB +
      lkdrug_apremilast * TRT_APREMILAST +
      lkdrug_bimekizumab * TRT_BIMEKIZUMAB +
      lkdrug_deucravacitinib * TRT_DEUCRAVACITINIB +
      lkdrug_etanercept * TRT_ETANERCEPT +
      lkdrug_guselkumab * TRT_GUSELKUMAB +
      lkdrug_ixekizumab * TRT_IXEKIZUMAB +
      lkdrug_risankizumab * TRT_RISANKIZUMAB +
      lkdrug_secukinumab * TRT_SECUKINUMAB +
      lkdrug_ustekinumab * TRT_USTEKINUMAB

    # ---- Treatment maximum effect per endpoint-indication (3.2) -------
    dmaxPasi75Pso <- deffPso + xDrg + deff_ph3_pso * ph3Trt +
      emaxDrg * relDose / (relDose + exp(led50Drg))
    dmaxPasi75Psa <- deffPsa + xeff_psa * trtSingle + xDrg +
      deff_ph3_psa * ph3Trt +
      emaxDrg * relDose / (relDose + exp(led50_psa * trtEmax + led50Drg))
    dmaxAcr20Psa <- deffPsa + deffAcr20 +
      (xeff_psa + xeff_psa_acr20) * trtSingle + xDrg +
      deff_ph3_psa * ph3Trt +
      (emax_psa_acr20 * trtEmax + emaxDrg) * relDose /
        (relDose + exp((led50_psa + led50_psa_acr20) * trtEmax + led50Drg))

    # ---- Treatment onset (Equation 3, S.3.1) --------------------------
    onsetTrtPso <- 1 - exp(-exp(lkDrg) * time)
    onsetTrtPsa <- 1 - exp(-exp(lkdrug_psa + lkDrg) * time)
    onsetTrtAcr20 <- 1 - exp(-exp(lkdrug_psa + lkdrug_psa_acr20 + lkDrg) * time)

    # ---- Responder probabilities (Equation 1.2) -----------------------
    prob_pasi75 <- expit(pboPasi75Pso + active * dmaxPasi75Pso * onsetTrtPso)
    prob_pasi75_psa <- expit(pboPasi75Psa + active * dmaxPasi75Psa * onsetTrtPsa)
    prob_acr20 <- expit(pboAcr20Psa + active * dmaxAcr20Psa * onsetTrtAcr20)

    # ---- Residual: binomial SE of an arm proportion (Equation 1.1) ----
    sdPasi75 <- sqrt(prob_pasi75 * (1 - prob_pasi75) / N_ARM) * addSd_prob_pasi75
    sdPasi75Psa <- sqrt(prob_pasi75_psa * (1 - prob_pasi75_psa) / N_ARM) *
      addSd_prob_pasi75_psa
    sdAcr20 <- sqrt(prob_acr20 * (1 - prob_acr20) / N_ARM) * addSd_prob_acr20
    prob_pasi75 ~ add(sdPasi75)
    prob_pasi75_psa ~ add(sdPasi75Psa)
    prob_acr20 ~ add(sdAcr20)
  })
}
