Desai_2016_isavuconazole_secure <- function() {
  description <- "Two-compartment population PK model with a Weibull absorption function and first-order elimination for isavuconazole (the active moiety of the prodrug isavuconazonium sulfate), fit to pooled data from nine phase 1 studies in healthy adults and the phase 3 SECURE trial in adults with invasive aspergillosis or other filamentous-fungal infections (Desai 2016). Race (Asian vs predominantly Caucasian) on clearance; BMI and healthy-vs-patient status on peripheral volume. Clearance inter-individual variability is stratified by healthy-vs-patient cohort. Supports both p.o. and i.v. dosing."
  reference <- "Desai A, Kovanda L, Kowalski D, Lu Q, Townsend R, Bonate PL. Population Pharmacokinetics of Isavuconazole from Phase 1 and Phase 3 (SECURE) Trials in Adults and Target Attainment in Patients with Invasive Infections Due to Aspergillus and Other Filamentous Fungi. Antimicrobial Agents and Chemotherapy. 2016;60(9):5483-5491. doi:10.1128/AAC.02819-15"
  vignette <- "Desai_2016_isavuconazole_secure"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. A(1) / A(2) / A(3) of the Desai 2016 equation system on
  # p. 5484 are respectively the gut, central and peripheral amounts.
  compartmentData <- list(
    depot = list(analyte = "isavuconazole", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "isavuconazole", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "isavuconazole", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    RACE_ASIAN = list(
      description = "Asian race indicator (1 = Asian, 0 = predominantly Caucasian)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the group Desai 2016 Table 3 labels 'predominantly Caucasian'; 175/189 healthy subjects and 193/232 patients)",
      notes = "Desai 2016 Table 2 codes the source covariate as 'Race (0 for predominantly Caucasians, 1 for Asians)', which is already the canonical orientation, so no re-expression was needed. Race was the only statistically significant covariate on CL and the only one the paper considered clinically important. Asian CL is 1.51 L/h against the Caucasian 2.36 L/h (Table 5 theta_9 vs theta_1), i.e. about 36% lower, which the Discussion describes as a roughly 40% difference in exposure. The mechanism is explicitly unestablished: the Discussion rules out CYP2D6 and CYP2C19 (isavuconazole is a substrate for neither) and rules out body mass, since BMI was not significant on CL.",
      source_name = "Race"
    ),
    DIS_HEALTHY = list(
      description = "Healthy-participant cohort indicator (1 = healthy phase 1 volunteer, 0 = patient with an invasive fungal infection)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (patient with invasive aspergillosis or another filamentous-fungal infection, from the phase 3 SECURE trial)",
      notes = "Desai 2016 Table 2 codes the source covariate as 'SP (dichotomized into healthy subjects [0] and patients [1])' -- the OPPOSITE orientation to the canonical -- so the model file re-expresses it as DIS_HEALTHY = 1 - SP and anchors the structural typical values on the patient state. Two roles in this model: (a) a multiplicative fractional effect on peripheral volume, and (b) a binary stratifier selecting which of the two published clearance IIV terms applies (same dual-variance pattern as Li_2017_CC292.R and Mao_2012_vernakalant.R). Assignment of Table 5 theta_4 = 417 L to the PATIENT baseline and theta_11 = 260 L to the HEALTHY baseline is not stated in the table labels; it is pinned by back-solving each against the Discussion's typical volumes at the Table 3 median BMIs -- 417 * (1 + 0.060 * (23.6 - 24.80)) = 387 L against the stated patient value of about 390 L, and 260 * (1 + 0.060 * (25.7 - 24.80)) = 274 L against the stated healthy value of about 292 L. The opposite assignment yields 440 L and 241 L, contradicting both. The reference-BMI Vss, 49.1 + 417 = 466 L, independently matches the Discussion's 'V at steady state was approximately 460 liters'.",
      source_name = "SP"
    ),
    BMI = list(
      description = "Body mass index at baseline",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Centering value 24.80 kg/m^2 is taken verbatim from the Desai 2016 best-covariate-model equation on p. 5485. Desai 2016 prints that equation as 'V p = theta_{4,11} x (1 + theta_10) x (BMI - 24.80)', which is dimensionally incoherent as typeset -- it makes Vp zero at BMI 24.80. It is encoded here in the sensible centered-linear reading, Vp = theta_{4,11} * (1 + theta_10 * (BMI - 24.80)), which is the only form under which the Discussion's typical volumes are recovered (see the DIS_HEALTHY note). BMI was significant on Vp but explicitly NOT on CL, and the Discussion notes that although the 43 obese individuals had a larger Vp there was no corresponding difference in exposure.",
      source_name = "BMI"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 421,
    n_studies = 10,
    n_observations = 6363,
    age_range = "17-85 years (median 43 healthy, 54 patients per Table 3)",
    weight_range = "41.0-127.7 kg (median 77.8 healthy, 67.0 patients per Table 3)",
    bmi_range = "13.9-41.2 kg/m^2 (median 25.7 healthy, 23.6 patients per Table 3)",
    sex_female_pct = 35,
    race_ethnicity = "Predominantly Caucasian 175/189 (92.6%) healthy and 193/232 (83.2%) patients; Asian 14/189 (7.4%) healthy and 39/232 (16.8%) patients (Table 3).",
    disease_state = "189 healthy volunteers from nine phase 1 studies (including dedicated hepatic-impairment, renal-impairment, mass-balance, bioavailability and elderly studies) and 232 patients with invasive aspergillosis or other filamentous-fungal infections from the phase 3 SECURE trial.",
    dose_range = "Phase 1: single or multiple doses of 40 mg to 400 mg isavuconazole, p.o. or as a 1-h i.v. infusion. Phase 3 SECURE: 372 mg isavuconazonium sulfate (equivalent to 200 mg isavuconazole) i.v. every 8 h for 6 doses on days 1-2, then 372 mg once daily p.o. or i.v. from day 3.",
    notes = "Healthy subjects contributed 5,828 rich-sampling concentrations and patients contributed 535 predominantly trough concentrations. One patient was excluded as an outlier for an extremely low clearance of 0.2 L/h. Below-quantification-limit values were under 5% of the healthy-subject data and were dropped; no patient concentration was below the quantification limit. Estimation in NONMEM 7.2 (ADVAN4 TRANS4) with FOCE and no interaction, since both the data and the residual-error structure were log transformed; covariate selection by stepwise covariate modeling in PsN 3.7.6 (forward p < 0.01, backward p < 0.001). Validated by 500-replicate nonparametric bootstrap (13% of runs failed) and NPDE; condition number 41."
  )

  ini({
    # Structural typical values. The reference subject is a predominantly
    # Caucasian PATIENT with an invasive fungal infection at BMI 24.80 kg/m^2
    # (RACE_ASIAN = 0, DIS_HEALTHY = 0). All values are Desai 2016 Table 5,
    # 'Parameter estimates of the best covariate model'.
    lcl <- log(2.36);   label("Isavuconazole clearance for the Caucasian reference (L/h)")           # Desai 2016 Table 5: theta_1 (CL, Caucasian) = 2.36 L/h
    lvc <- log(49.10);  label("Isavuconazole central volume of distribution (V_1, L)")               # Desai 2016 Table 5: theta_2 (V_1) = 49.10 L
    lq  <- log(26.60);  label("Isavuconazole intercompartmental clearance (Q, L/h)")                 # Desai 2016 Table 5: theta_3 (Q) = 26.60 L/h
    lvp <- log(417.0);  label("Isavuconazole peripheral volume for patients at BMI 24.80 (V_p, L)")  # Desai 2016 Table 5: theta_4 (V_p) = 417.0 L; assigned to the patient baseline -- see covariateData DIS_HEALTHY notes

    # Weibull absorption function (Desai 2016 p. 5484, second displayed
    # equation): WB = KAMAX x [1 - exp(-(RA x TAD)^GAM1)]. This is a
    # time-varying first-order absorption rate rather than a constant ka;
    # the Weibull form was chosen over first-order absorption on a drop in
    # objective function value of more than 1,000 points.
    lkamax <- log(1.08); label("Weibull-absorption asymptotic maximum absorption rate (KAMAX, 1/h)") # Desai 2016 Table 5: theta_5 (KAMAX) = 1.08 1/h
    lra    <- log(0.72); label("Weibull-absorption rate-scaling parameter (RA, 1/h)")                # Desai 2016 Table 5: theta_6 (RA) = 0.72 1/h
    lgam1  <- log(4.88); label("Weibull-absorption shape parameter (GAM1, unitless)")                # Desai 2016 Table 5: theta_7 (GAM1) = 4.88

    # Absolute bioavailability, fixed (not estimated) at 1 on the strength of
    # prior noncompartmental analyses in healthy subjects. Applies to the p.o.
    # route only; i.v. doses are routed directly to central and bypass depot.
    lfdepot <- fixed(log(1.0)); label("Isavuconazole oral bioavailability (fraction)")               # Desai 2016 Methods, 'Structural pharmacokinetic model': 'Absolute bioavailability (F) was fixed at 1'

    # Covariate effects.
    # Race on CL. Desai 2016 reports two separate clearance THETAs rather than
    # a coefficient: theta_1 = 2.36 L/h (Caucasian) and theta_9 = 1.51 L/h
    # (Asian). Re-expressed here as a single baseline plus a multiplicative
    # fractional effect, (1.51 - 2.36) / 2.36 = -0.3602, which reproduces the
    # Discussion's 'approximately 36% lower CL value' in Asians exactly.
    e_race_asian_cl  <- -0.3602; label("Fractional effect of Asian race on CL, relative to the Caucasian reference (unitless)")  # Derived from Desai 2016 Table 5: theta_9 (1.51 L/h) vs theta_1 (2.36 L/h)

    # Healthy-vs-patient status on V_p. Same re-expression: theta_4 = 417 L
    # (patient) and theta_11 = 260 L (healthy) become a patient baseline plus
    # (260 - 417) / 417 = -0.3765. Note the paper's SP is coded 1 = patient,
    # so the canonical DIS_HEALTHY = 1 - SP and the sign is flipped relative
    # to the paper's own indicator. Same re-expression pattern as
    # Lu_2015_tacrolimus.R, Yoneyama_2017_emicizumab.R and Taubert_2018_finafloxacin.R.
    e_dis_healthy_vp <- -0.3765; label("Fractional effect of healthy-participant status on V_p, relative to the patient reference (unitless)")  # Derived from Desai 2016 Table 5: theta_11 (260 L) vs theta_4 (417 L)

    # BMI on V_p, centered-linear at 24.80 kg/m^2.
    e_bmi_vp         <-  0.060; label("Linear-deviation effect of BMI on V_p, per kg/m^2 above 24.80 (per kg/m^2)")  # Desai 2016 Table 5: theta_10 (BMI on V_p) = 0.060

    # Inter-individual variability. Desai 2016 states that 'all random effects
    # were treated as log-normally (Ln) distributed' and reports the Table 5
    # 'Variability (%)' block as percent CV, so omega^2 = log(1 + CV^2).
    # IIV was retained only on CL, V_p, Q, RA and GAM1 -- the Results section
    # states explicitly that no other parameter showed statistically
    # significant IIV, so V_1, KAMAX and F carry none.
    #
    # Clearance IIV is reported SEPARATELY for the two cohorts, 31.30% CV in
    # healthy subjects and 62.44% CV in patients. That is a categorical OMEGA
    # structure, encoded here as two mutually-exclusive etas gated on
    # DIS_HEALTHY so that exactly one contributes for any given subject; the
    # same stratified-variance pattern is used by Li_2017_CC292.R and
    # Mao_2012_vernakalant.R (on the residual error) and by
    # Pohl_2022_linzagolix_e2.R.
    etalcl_hv ~ 0.09346  # Desai 2016 Table 5 Variability: 'CL (healthy subjects)' = 31.30% CV; omega^2 = log(1 + 0.3130^2)
    etalcl_pt ~ 0.32921  # Desai 2016 Table 5 Variability: 'CL (patients)' = 62.44% CV; omega^2 = log(1 + 0.6244^2)
    etalvp    ~ 0.09622  # Desai 2016 Table 5 Variability: 'V p' = 31.78% CV; omega^2 = log(1 + 0.3178^2)
    etalq     ~ 0.21590  # Desai 2016 Table 5 Variability: 'Q' = 49.09% CV; omega^2 = log(1 + 0.4909^2)
    etalra    ~ 0.15008  # Desai 2016 Table 5 Variability: 'RA' = 40.24% CV; omega^2 = log(1 + 0.4024^2)
    etalgam1  ~ 0.18974  # Desai 2016 Table 5 Variability: 'GAM1' = 45.71% CV; omega^2 = log(1 + 0.4571^2)

    # Residual error. Desai 2016 Methods: 'The Ln-Ln transformations of both
    # the output equation in the model and the data were used to stabilize the
    # residual variance. The residual variance was modeled as additive in
    # nature.' An additive residual on the natural-log scale is exactly a
    # log-normal residual in linear space, ln(DV) = ln(IPRED) + eps, so the
    # Table 5 weighting factor W is carried as expSd rather than propSd.
    # At this magnitude the distinction matters: a proportional residual with
    # SD 0.4494 generates negative simulated concentrations, whereas lnorm
    # cannot. See the vignette's Assumptions and deviations section.
    expSd <- 0.4494; label("Isavuconazole residual error, additive on the natural-log scale (log-scale SD)")  # Desai 2016 Table 5: theta_8 (W) = 44.94%
  })

  model({
    # Race effect on clearance (multiplicative; Caucasian reference).
    # Clearance IIV is selected by cohort: healthy subjects draw etalcl_hv and
    # patients draw etalcl_pt, exactly one of which is active per subject.
    race_cl <- 1 + e_race_asian_cl * RACE_ASIAN
    cl_hv <- exp(lcl + etalcl_hv)
    cl_pt <- exp(lcl + etalcl_pt)
    cl <- (cl_hv * DIS_HEALTHY + cl_pt * (1 - DIS_HEALTHY)) * race_cl

    vc <- exp(lvc)
    q  <- exp(lq + etalq)

    # Healthy-vs-patient and BMI effects on peripheral volume (multiplicative).
    healthy_vp <- 1 + e_dis_healthy_vp * DIS_HEALTHY
    bmi_vp     <- 1 + e_bmi_vp * (BMI - 24.80)
    vp <- exp(lvp + etalvp) * healthy_vp * bmi_vp

    # Weibull absorption parameters.
    kamax <- exp(lkamax)
    ra    <- exp(lra + etalra)
    gam1  <- exp(lgam1 + etalgam1)

    # Micro-constants. Kcp and Kpc of the paper's equation system.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Time-varying absorption rate WB = KAMAX x [1 - exp(-(RA x TAD)^GAM1)].
    # tad(depot) restarts at each oral dose, so a multiple-dose regimen
    # reproduces the Weibull rise from every dose. The guard handles i.v.-only
    # subjects and the pre-first-dose window, where tad(depot) is NA and would
    # otherwise propagate NaN through the whole system; same pattern as
    # Cirincione_2017_exenatide.R and the sibling Desai_2016_isavuconazole.R.
    wb <- 0.0
    if (tad(depot) >= 0.0) wb <- kamax * (1 - exp(-(ra * tad(depot))^gam1))

    # Desai 2016 p. 5484, three-ODE structural model.
    d/dt(depot)       <- -wb * depot
    d/dt(central)     <-  wb * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                               k12 * central - k21 * peripheral1

    # Oral bioavailability anchor. i.v. doses are addressed to central in the
    # event table and so are unaffected by it.
    f(depot) <- exp(lfdepot)

    # Plasma concentration: dose in mg over volume in L gives mg/L.
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
