AitOudhia_2024_sotatercept <- function() {
  description <- "Two-compartment population pharmacokinetic model for sotatercept in healthy post-menopausal women and patients with pulmonary arterial hypertension (Ait-Oudhia 2024). First-order subcutaneous absorption with logit-scale bioavailability (about 66%) and linear elimination from the central compartment; time-varying body weight enters as a power covariate on clearance (exponent 0.814) and central volume (exponent 1.02) with a 70 kg reference, and baseline serum albumin as a power covariate on clearance (exponent -0.849) with a 4.5 g/dL reference. Separate log-scale residual error magnitudes are used for healthy participants and for patients with PAH. Intravenous doses go to the central compartment and subcutaneous doses to the depot."
  reference <- "Ait-Oudhia S, Jaworowicz D, Hu Z, Bihorel S, Hu S, Balasubrahmanyam B, Mistry B, de Oliveira Pena J, Wenning L, Gheyas F. Population pharmacokinetic modeling of sotatercept in healthy participants and patients with pulmonary arterial hypertension. CPT Pharmacometrics Syst Pharmacol. 2024;13(8):1380-1393. doi:10.1002/psp4.13166"
  vignette <- "AitOudhia_2024_sotatercept"
  paper_specific_residual_sds <- c("expSdHealthy", "expSdPah")
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    depot = list(
      analyte = "sotatercept", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "sotatercept", units = "mg",
      specimen = "serum", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "sotatercept", units = "mg",
      specimen = "tissue", verified = FALSE
    )
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Ait-Oudhia 2024 uses the time-varying weight column (source name WTKGT) rather than the baseline weight, so a body weight recorded at each visit drives both CL and VC over the ~150-week follow-up. Enters as a power (allometric-style) term with a 70 kg reference: (WT/70)^0.814 on CL and (WT/70)^1.02 on VC (Table 2 footnotes a and c). Both exponents were estimated, not fixed at the canonical 0.75 / 1.0 allometric values, but the paper notes they came out 'close to the allometry scaling factors' (Discussion). Baseline weight in the analysis population was 71 kg mean (SD 17.6), range 39.6-136 kg (Table 1).",
      source_name        = "WTKGT"
    ),
    ALB = list(
      description        = "Baseline serum albumin concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline (time-fixed), not time-varying. Enters CL as a power term (ALB/45)^-0.849 (Table 2 footnote a). IMPORTANT UNIT NOTE: the canonical register stores ALB in SI g/L, whereas Ait-Oudhia 2024 reports albumin in US-convention g/dL and writes the reference as 4.5 g/dL. Because the covariate enters only as a reference-normalised ratio, the ratio is scale-invariant -- (ALB_gL / 45) is identically equal to (ALB_gdL / 4.5) -- so model() uses a 45 g/L reference and reproduces the published equation exactly with no conversion line. Baseline albumin in the analysis population was 4.43 g/dL mean (SD 0.33), range 2.9-5.8 g/dL, i.e. 44.3 g/L mean, range 29-58 g/L (Table 1). The negative exponent means CL decreases as albumin rises; the paper attributes this to FcRn-mediated protection from catabolism shared between albumin and IgG (Discussion).",
      source_name        = "ALB"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-participant cohort indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (patient with pulmonary arterial hypertension, from PULSAR / SPECTRA / STELLAR)",
      notes              = "1 = healthy post-menopausal woman from one of the two phase 1 trials (SAD or MAD); 0 = patient with PAH from PULSAR, SPECTRA, or STELLAR. Time-fixed per subject. Used ONLY to select the residual-error magnitude: Results 'Population PK model' states 'Two separate log error models were used to describe the RV for healthy and PAH participants', giving log-scale SD 0.239 in healthy participants and 0.189 in patients with PAH (Table 2). Disease status was also tested as a covariate on the structural PK parameters and entered during forward selection, but was NOT retained during backward elimination (Discussion), so it has no effect on CL, VC, VP, Q, KA, or F1 in the final model.",
      source_name        = "Healthy vs PAH study population"
    )
  )

  # Covariates that Ait-Oudhia 2024 screened but did not retain in the final
  # model. Documented here (not in covariateData) to preserve the provenance of
  # the paper's covariate screen without carrying "declared but not referenced"
  # convention warnings, since none of these appear in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age",
      units       = "years",
      type        = "continuous",
      notes       = "Tested (Table 1) and reported in the Figure 4 forest plot, but not a statistically significant source of IIV in sotatercept PK (Discussion). The Figure 4 GMR for age was 0.998, well inside the (0.8, 1.25) bioequivalence bounds. Population mean 50.2 years (SD 14.4), range 18-81 (Table 1)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested (Table 1) and shown in the Figure 4 forest plot; not statistically significant (Discussion). The cohort was 85.1% female (298 of 350) because both phase 1 trials enrolled only post-menopausal women, so the male stratum is small (52 subjects)."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Listed among the tested covariates (Table 1); not retained. Population mean 26.5 kg/m^2 (SD 5.48), range 15.1-50.3; not recorded in the SAD trial (Table 1 reports NA for that study). Body size entered the final model through WT rather than BMI."
    ),
    CRCL = list(
      description = "Baseline estimated glomerular filtration rate (BSA-normalized renal function)",
      units       = "mL/min/1.73 m^2",
      type        = "continuous",
      notes       = "Tested as a continuous baseline eGFR (Table 1) and, separately, as a renal-function category; neither was significant (Discussion). Population mean 84.9 mL/min/1.73 m^2 (SD 28.8), range 31.3-255. The paper notes sotatercept is a large fusion protein not expected to undergo renal filtration, and that severe impairment (eGFR < 30) was an exclusion criterion, so the tested range does not extend to severe renal impairment."
    ),
    RENALIMP_MILD = list(
      description = "Mild renal impairment indicator (eGFR 60-90 mL/min/1.73 m^2)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Level of the tested 'Renal function category' covariate (Table 1) and a Figure 4 forest-plot stratum; not significant (Discussion). 183 of 350 subjects (52.3%) were in this category."
    ),
    RENALIMP_MOD = list(
      description = "Moderate renal impairment indicator (eGFR 30-60 mL/min/1.73 m^2)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Level of the tested 'Renal function category' covariate (Table 1) and a Figure 4 forest-plot stratum; not significant (Discussion). 51 of 350 subjects (14.6%) were in this category. The reference stratum was normal renal function (eGFR 90-120), 116 subjects (33.1%). No subject had severe impairment (an exclusion criterion), so RENALIMP_SEV was not testable."
    ),
    RACE_WHITE = list(
      description = "White / Caucasian race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Level of the tested 'Racial classification' covariate (Table 1) and a Figure 4 forest-plot stratum; not significant (Discussion). 290 of 350 subjects (82.9%). Because the effect was not retained, the paper does not report which race level served as the reference category in the tested parameterisation."
    ),
    RACE_BLACK = list(
      description = "Black / African American race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Level of the tested 'Racial classification' covariate (Table 1); not significant (Discussion). 11 of 350 subjects (3.14%)."
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Level of the tested 'Racial classification' covariate (Table 1); not significant (Discussion). 25 of 350 subjects (7.14%)."
    ),
    RACE_OTHER = list(
      description = "Race-category 'Other' indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Level of the tested 'Racial classification' covariate (Table 1); not significant (Discussion). 19 of 350 subjects (5.43%) self-identified as Other. Table 1 additionally reports 1 subject (0.286%) as American Indian or Alaska Native and 4 (1.14%) as Native Hawaiian or other Pacific Islander; the paper does not state how these very small strata were grouped for covariate testing, so they are not given separate entries here."
    ),
    ADA_POS = list(
      description = "Anti-drug-antibody-positive status",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested as a TIME-VARYING covariate and not statistically significant on sotatercept PK (Discussion). The paper reports similar PK time profiles, trough concentrations, and clearance across ADA-negative, ADA-positive / neutralizing-antibody-negative, and ADA-positive / NAb-positive participants (Discussion, citing reference 34). Immunogenicity therefore has no effect in the final model."
    ),
    CONMED_ERA = list(
      description = "Concomitant endothelin-receptor antagonist alone as background PAH standard-of-care therapy",
      units       = "(binary)",
      type        = "binary",
      notes       = "The paper tested background PAH therapy as mono, double, and triple combination with standard of care and found no statistically significant effect on sotatercept PK (Discussion, Figure 4). The paper reports the covariate at the level of the NUMBER of background agents rather than by drug class, so this class-specific indicator is a documentation-only mapping onto the existing canonical family (founded by Krause 2017 selexipag, also a PAH population); Ait-Oudhia 2024 does not report a per-class coefficient."
    ),
    CONMED_PDE5I = list(
      description = "Concomitant phosphodiesterase type 5 inhibitor alone as background PAH standard-of-care therapy",
      units       = "(binary)",
      type        = "binary",
      notes       = "See the CONMED_ERA note: background PAH therapy was tested as mono / double / triple combination with standard of care and was not significant (Discussion, Figure 4). No per-class coefficient is reported."
    ),
    CONMED_ERA_PDE5I = list(
      description = "Concomitant endothelin-receptor antagonist plus phosphodiesterase type 5 inhibitor as background PAH standard-of-care therapy",
      units       = "(binary)",
      type        = "binary",
      notes       = "See the CONMED_ERA note: background PAH therapy was tested as mono / double / triple combination with standard of care and was not significant (Discussion, Figure 4). The paper concludes sotatercept is not a victim of drug-drug interaction with background PAH therapies. No per-class coefficient is reported."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 350L,
    n_studies      = 5L,
    age_range      = "18-81 years",
    age_mean       = "50.2 years (SD 14.4) (Table 1)",
    weight_range   = "39.6-136 kg (baseline)",
    weight_mean    = "71.0 kg (SD 17.6) at the baseline visit (Table 1)",
    sex_female_pct = 85.1,
    race_ethnicity = c(
      White = 82.9, Black = 3.14, Asian = 7.14,
      `American Indian or Alaska Native` = 0.286,
      `Native Hawaiian or Other Pacific Islander` = 1.14,
      Other = 5.43
    ),
    disease_state  = "Pooled healthy participants and patients with pulmonary arterial hypertension. The two phase 1 trials (a single-ascending-dose and a multiple-ascending-dose study) enrolled healthy post-menopausal women; the two phase 2 trials (PULSAR, SPECTRA) and the phase 3 trial (STELLAR) enrolled participants with PAH on background standard-of-care therapy. Disease status was tested as a PK covariate and entered during forward selection but was not retained during backward elimination; it is retained in the final model only as the switch between the two residual-error magnitudes.",
    renal_function = "Normal (eGFR 90-120 mL/min/1.73 m^2) 116 (33.1%); mild impairment (60-90) 183 (52.3%); moderate impairment (30-60) 51 (14.6%). Severe impairment (eGFR < 30) was an exclusion criterion. Baseline eGFR mean 84.9 mL/min/1.73 m^2 (SD 28.8), range 31.3-255 (Table 1).",
    albumin        = "Baseline albumin mean 4.43 g/dL (SD 0.33), range 2.9-5.8 g/dL, i.e. mean 44.3 g/L, range 29-58 g/L (Table 1).",
    bmi            = "Mean 26.5 kg/m^2 (SD 5.48), range 15.1-50.3; not collected in the SAD trial (Table 1).",
    dose_range     = "Single IV or SC doses 0.01-3.0 mg/kg (phase 1 SAD) and multiple SC doses 0.03-1.0 mg/kg Q4W (phase 1 MAD); SC 0.3 or 0.7 mg/kg Q3W in phase 2 (PULSAR, SPECTRA); SC 0.3 mg/kg initial dose followed by the 0.7 mg/kg Q3W target dose in phase 3 (STELLAR). 30 participants received sotatercept IV and 320 SC.",
    regions        = "Not reported by the paper. STELLAR, PULSAR, and SPECTRA were multinational PAH trials.",
    notes          = "Baseline demographics are from Table 1, which reports them per study (SAD n = 40, MAD n = 24, PULSAR n = 103, SPECTRA n = 21, STELLAR n = 162) and overall (n = 350). Rich PK sampling in the two phase 1 studies, sparse PK in the phase 2 and phase 3 studies (Figure 1); PK samples were collected up to a maximum of about 150 weeks. Below-limit-of-quantitation samples were under 5% of records and were excluded (M1 method). Estimation used NONMEM 7 level 3.0 FOCE with interaction. Minimum value of the objective function -6707.007; condition number 124. Qualification was by prediction-corrected VPC (500 replicates) and a 1000-replicate bootstrap in which 99.3% of runs minimized successfully and every final estimate fell inside its bootstrap 95% CI (Table 2, Results 'Model qualification')."
  )

  ini({
    # ---- Structural parameters (Ait-Oudhia 2024 Table 2) ----
    # Reference individual: 70 kg time-varying body weight, 4.5 g/dL (45 g/L)
    # baseline albumin. Time in days, amounts in mg, concentration in mg/L.
    # Every parameter below was estimated (each carries a %RSE in Table 2), so
    # none is wrapped in fixed().

    lcl     <- log(0.177)   ; label("Central clearance CL (L/day)")                          # Table 2 CL population mean = 0.177 L/day (%RSE 6.22; bootstrap 0.176, 95% CI 0.150-0.200)
    lvc     <- log(3.60)    ; label("Central volume of distribution VC (L)")                 # Table 2 VC population mean = 3.60 L (%RSE 5.50; bootstrap 3.61, 95% CI 3.25-4.11)
    lq      <- log(0.487)   ; label("Distribution clearance Q (L/day)")                      # Table 2 Q population mean = 0.487 L/day (%RSE 20.5; bootstrap 0.476, 95% CI 0.239-0.677)
    lvp     <- log(1.71)    ; label("Peripheral volume of distribution VP (L)")              # Table 2 VP population mean = 1.71 L (%RSE 15.8; bootstrap 1.66, 95% CI 0.740-2.24)
    lka     <- log(0.273)   ; label("First-order subcutaneous absorption rate constant KA (1/day)")  # Table 2 KA population mean = 0.273 1/day (%RSE 16.0; bootstrap 0.275, 95% CI 0.199-0.371)

    # Subcutaneous bioavailability F1, carried on the logit scale because the
    # paper's IIV on F1 is normal on the logit scale (Table 2 footnote d), which
    # keeps every individual F1 inside (0, 1).
    logitfdepot <- log(0.659 / (1 - 0.659)); label("Logit of subcutaneous bioavailability F1 (unitless; F1_pop = 0.659)")  # Table 2 F1 population mean = 0.659 (%RSE 6.29; bootstrap 0.657, 95% CI 0.557-0.745). logit(0.659) = 0.6588. Table 2 footnote b: F1 and CL were highly correlated (r^2 >= 0.81)

    # ---- Covariate effects (Ait-Oudhia 2024 Table 2 footnotes a and c) ----
    # Footnote a: CL = 0.177 * (WTKGT/70)^0.814 * (ALB/4.5)^-0.849, ALB in g/dL.
    # Footnote c: VC = 3.60 * (WTKGT/70)^1.02.
    e_wt_cl  <- 0.814       ; label("Exponent of the body-weight power effect on CL (unitless)")             # Table 2 CL 'Exponent of body weight power effect' = 0.814 (%RSE 11.9; bootstrap 0.809, 95% CI 0.604-1.00)
    e_alb_cl <- -0.849      ; label("Exponent of the baseline-albumin power effect on CL (unitless)")        # Table 2 CL 'Exponent of baseline albumin effect' = -0.849 (%RSE 24.2; bootstrap -0.851, 95% CI -1.27 to -0.454)
    e_wt_vc  <- 1.02        ; label("Exponent of the body-weight power effect on VC (unitless)")             # Table 2 VC 'Exponent of body weight power effect' = 1.02 (%RSE 11.8; bootstrap 0.999, 95% CI 0.732-1.25)

    # ---- Inter-individual variability ----
    # Table 2 reports IIV as %CV. For the exponential (log-normal) IIV terms the
    # paper's own back-transformation is %CV = sqrt(exp(omega^2) - 1) * 100
    # (Methods), so omega^2 = log(1 + CV^2). No correlations between etas are
    # reported, so each eta is entered on its own diagonal.
    etalcl ~ log(1 + 0.283^2)   # Table 2 CL IIV = 28.3 %CV (%RSE 16.9; bootstrap 28.2 %CV, 95% CI 23.1-32.8%) -> omega^2 = 0.07704, omega = 0.2776. Shrinkage 17.3%
    etalvc ~ log(1 + 0.247^2)   # Table 2 VC IIV = 24.7 %CV (%RSE 31.0; bootstrap 24.0 %CV, 95% CI 15.3-31.7%) -> omega^2 = 0.05922, omega = 0.2434. Shrinkage 51.9%
    etalvp ~ log(1 + 0.733^2)   # Table 2 VP IIV = 73.3 %CV (%RSE 28.2; bootstrap 81.9 %CV, 95% CI 52.9-188.2%) -> omega^2 = 0.43002, omega = 0.6558. Shrinkage 34.6%
    etalka ~ log(1 + 0.604^2)   # Table 2 KA IIV = 60.4 %CV (%RSE 24.2; bootstrap 60.4 %CV, 95% CI 41.6-77.3%) -> omega^2 = 0.31102, omega = 0.5577. Shrinkage 48.3%

    # IIV on Q is deliberately absent: Table 2 reports "NE" (not estimated) in
    # the Q variability column.

    # F1 IIV is normal on the LOGIT scale, not log-normal (Table 2 footnote d),
    # so the %CV back-transformation differs: %CV = (1 - F1_pop) * omega * 100.
    # Footnote d spells the arithmetic out as 100 * (1 - 0.659) * 0.348 = 11.9,
    # which identifies the logit-scale omega as 0.348 directly.
    etalogitfdepot ~ 0.348^2    # Table 2 F1 IIV = 11.9 %CV (%RSE 26.3; bootstrap 11.9 %CV, 95% CI 5.8-19.8%) via footnote d: omega = 0.348 on the logit scale -> omega^2 = 0.1211. Shrinkage 36.8%

    # ---- Residual variability ----
    # The paper used a log error model, log(Y_ij) = log(Yhat_ij) + eps_ij, i.e. a
    # constant variance on the natural-log concentration scale; this is
    # nlmixr2's `~ lnorm(expSd)`. Two separate magnitudes were estimated, one
    # for healthy participants and one for patients with PAH (Results
    # 'Population PK model'). Table 2's "Final parameter estimate" column holds
    # the VARIANCE and its "Magnitude of variability" column the SD:
    # sqrt(0.0570) = 0.2387 -> 0.239 and sqrt(0.0357) = 0.1889 -> 0.189, which
    # confirms which column is which.
    expSdHealthy <- sqrt(0.0570) ; label("Log-scale residual SD in healthy participants (natural-log concentration units)")  # Table 2 'Residual variability in HV (log units)' variance = 0.0570 (%RSE 20.2); reported SD 0.239 (bootstrap 0.240, 95% CI 0.194-0.289)
    expSdPah     <- sqrt(0.0357) ; label("Log-scale residual SD in patients with PAH (natural-log concentration units)")     # Table 2 'Residual variability in PAH (log units)' variance = 0.0357 (%RSE 12.6); reported SD 0.189 (bootstrap 0.189, 95% CI 0.167-0.214)
  })

  model({
    # ---- 1. Individual PK parameters ----
    # Body weight is the time-varying WTKGT column; albumin is the baseline
    # value. The albumin reference is written as 45 g/L rather than the paper's
    # 4.5 g/dL: the covariate enters only as a reference-normalised ratio, and
    # (ALB_gL / 45) is identically (ALB_gdL / 4.5), so this reproduces Table 2
    # footnote a exactly while keeping ALB in the canonical SI g/L units.
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * (ALB / 45)^e_alb_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    vp <- exp(lvp + etalvp)
    q  <- exp(lq)
    ka <- exp(lka + etalka)

    # Subcutaneous bioavailability, with IIV added on the logit scale so that
    # every individual value stays inside (0, 1) (Table 2 footnote d). The
    # logit-scale sum is kept on its own line so nlmixr2 can mu-reference
    # etalogitfdepot (same idiom as Martinez_2019_alirocumab).
    logit_f <- logitfdepot + etalogitfdepot
    fdepot  <- 1 / (1 + exp(-logit_f))

    # ---- 2. Micro-constants (Methods 'Clinical relevance of covariates') ----
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- 3. ODE system: two compartments with first-order SC absorption ----
    # Subcutaneous doses go to `depot` and carry F1; intravenous doses go
    # straight to `central` with F = 1 (the paper used data-defined infusion
    # rates for the IV arms), so one model serves both routes.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                               k12 * central - k21 * peripheral1

    # ---- 4. Bioavailability ----
    f(depot) <- fdepot

    # ---- 5. Observation and population-specific residual error ----
    # central is in mg and vc in L, so Cc is in mg/L (= ug/mL).
    Cc <- central / vc

    # Two separate log error models, selected by cohort (Results 'Population PK
    # model'). DIS_HEALTHY = 1 for the phase 1 healthy post-menopausal women,
    # 0 for patients with PAH.
    w_rv <- expSdHealthy * DIS_HEALTHY + expSdPah * (1 - DIS_HEALTHY)

    Cc ~ lnorm(w_rv)
  })
}
