Desai_2016_isavuconazole <- function() {
  description <- "Two-compartment population PK model with a Weibull-function first-order absorption stage for isavuconazole administered as the prodrug isavuconazonium sulfate (p.o. or i.v.) to healthy adults and patients with invasive fungal infections (Desai 2016 SECURE pooled phase 1 / phase 3 popPK)"
  reference <- "Desai A, Kovanda L, Kowalski D, Lu Q, Townsend R, Bonate PL. Population Pharmacokinetics of Isavuconazole from Phase 1 and Phase 3 (SECURE) Trials in Adults and Target Attainment in Patients with Invasive Infections Due to Aspergillus and Other Filamentous Fungi. Antimicrob Agents Chemother. 2016;60(9):5483-5491. doi:10.1128/AAC.02819-15"
  vignette <- "Desai_2016_isavuconazole"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    BMI = list(
      description        = "Body mass index at baseline",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Centered-linear effect on V_p with reference 24.80 kg/m^2 (Desai 2016 best-model covariate paragraph and Table 5 theta_10).",
      source_name        = "BMI"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator (1 = Asian, 0 = predominantly Caucasian)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (predominantly Caucasian; per Desai 2016 Table 3 footnote, the n=368 'predominantly Caucasian' subjects also include 1 African American and 5 in other-race categories pooled into the reference)",
      notes              = "Multiplicative fractional effect on CL: Asian subjects had ~36% lower CL than the Caucasian reference. Source column RACE coded as 0 = predominantly Caucasian, 1 = Asian (Desai 2016 Table 2).",
      source_name        = "RACE (0 = predominantly Caucasian, 1 = Asian)"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-participant indicator (1 = healthy subject, 0 = patient with invasive fungal infection)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (patient with invasive aspergillosis or other IFI; SECURE trial cohort, n=232)",
      notes              = "Multiplicative fractional effect on V_p (the SP-on-Vp covariate in Desai 2016 Table 5 theta_11). Desai 2016 codes SP = 1 for patients and SP = 0 for healthy subjects; the canonical DIS_HEALTHY uses the opposite sign convention (1 = healthy), so DIS_HEALTHY = 1 - SP. Healthy subjects had ~38% lower V_p than patients at the reference BMI of 24.80 kg/m^2. Assignment of theta_4 = 417 L to patients and theta_11 = 260 L to healthy subjects was inferred from the Discussion's typical-value report (V_p ~390 L for patients vs ~292 L for healthy subjects); see vignette Assumptions for the derivation.",
      source_name        = "SP (1 = patient, 0 = healthy); DIS_HEALTHY = 1 - SP"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 421L,
    n_studies      = 10L,
    n_observations = 6363L,
    age_range      = "17-85 years",
    age_median     = "43 years (healthy subjects) / 54 years (patients)",
    weight_range   = "41.0-127.7 kg",
    weight_median  = "77.8 kg (healthy subjects) / 67.0 kg (patients)",
    sex_female_pct = 35.4,
    race_ethnicity = c(predominantly_Caucasian = 87.4, Asian = 12.6),
    disease_state  = "Pooled cohort: 189 healthy subjects from 9 phase 1 studies (single- and multiple-ascending dose, hepatic impairment, renal impairment, mass balance, bioavailability crossover, elderly) plus 232 patients with invasive aspergillosis or other filamentous-fungi IFIs enrolled in the phase 3 SECURE trial.",
    dose_range     = "Isavuconazole 40-400 mg, p.o. or i.v. (1-h infusion), single or multiple doses across phase 1; SECURE: 200 mg isavuconazole-equivalent (372 mg isavuconazonium sulfate) q8h x 6 loading doses on days 1-2, then 200 mg q24h to end of treatment.",
    regions        = "International (phase 1 sites plus the global SECURE phase 3 trial including Asian regions)",
    notes          = "Dosed as the water-soluble prodrug isavuconazonium sulfate; the model represents the active moiety isavuconazole at the molar equivalent (372 mg prodrug = 200 mg isavuconazole). Absolute bioavailability F was fixed at 1 based on noncompartmental analyses in healthy subjects. One outlier with CL = 0.2 L/h was removed prior to fitting. Demographics summarised in Desai 2016 Table 3; data sources in Table 1."
  )

  ini({
    # Structural parameters -- typical-value reference is a predominantly
    # Caucasian patient with an invasive fungal infection at BMI 24.80 kg/m^2.
    lcl <- log(2.36);  label("Clearance (CL) for the Caucasian reference (L/h)")                       # Desai 2016 Table 5 theta_1
    lvc <- log(49.10); label("Central volume of distribution (V_1, L)")                                # Desai 2016 Table 5 theta_2
    lq  <- log(26.60); label("Intercompartmental clearance (Q, L/h)")                                  # Desai 2016 Table 5 theta_3
    lvp <- log(417.0); label("Peripheral volume (V_p) for patients at BMI 24.80 (L)")                  # Desai 2016 Table 5 theta_4

    # Weibull absorption parameters (Desai 2016 WB equation):
    # WB = KAMAX * (1 - exp(-(RA * TAD)^GAM1))
    lkamax <- log(1.08); label("Maximum first-order absorption rate cap (KAMAX, 1/h)")                 # Desai 2016 Table 5 theta_5
    lra    <- log(0.72); label("Weibull rate parameter (RA, 1/h)")                                     # Desai 2016 Table 5 theta_6
    lgam1  <- log(4.88); label("Weibull shape parameter (GAM1, unitless)")                             # Desai 2016 Table 5 theta_7

    # Oral bioavailability fixed at 1 from prior NCA of healthy subjects.
    # The i.v. route bypasses depot, so this only affects p.o. dosing.
    lfdepot <- fixed(log(1)); label("Oral bioavailability F (fraction)")                               # Desai 2016 Methods (F fixed at 1)

    # Covariate effects.
    # Race on CL: CL_Asian = 1.51 L/h vs CL_Caucasian = 2.36 L/h.
    # Fractional change: (1.51 - 2.36) / 2.36 = -0.3602.
    e_race_asian_cl  <- -0.3602; label("Fractional effect of Asian race on CL (Asian vs Caucasian reference)")   # Derived from Desai 2016 Table 5 theta_9 (1.51 L/h) and theta_1 (2.36 L/h)
    # Healthy on V_p: V_p_healthy = 260 L vs V_p_patient = 417 L.
    # Fractional change: (260 - 417) / 417 = -0.3765.
    e_dis_healthy_vp <- -0.3765; label("Fractional effect of healthy status on V_p (healthy vs patient reference)") # Derived from Desai 2016 Table 5 theta_11 (260 L) and theta_4 (417 L); assignment per Discussion typical V_p ~292 L (healthy) / ~390 L (patients)
    # BMI on V_p: centered-linear, 0.060 per (BMI - 24.80) kg/m^2.
    e_bmi_vp         <-  0.060;  label("Centered-linear coefficient for BMI on V_p (per kg/m^2 above 24.80)")   # Desai 2016 Table 5 theta_10

    # IIV. Desai 2016 reports separate CL IIV for healthy (31.30% CV) and
    # patient (62.44% CV) subjects, suggesting a categorical OMEGA structure
    # in NONMEM. nlmixr2 carries a single omega per parameter; the patient
    # CV (62.44%) is used as the operational default because the model is
    # intended primarily for clinical simulations of patients with IFIs.
    # See vignette Assumptions and deviations for the rationale.
    # Convention: omega^2 = log(CV^2 + 1) for log-normal IIV.
    etalcl   ~ 0.3293                                                                                  # Desai 2016 Table 5 CL (patients) CV 62.44%; omega^2 = log(0.6244^2 + 1) = 0.3293
    etalvp   ~ 0.0962                                                                                  # Desai 2016 Table 5 V_p CV 31.78%; omega^2 = log(0.3178^2 + 1) = 0.0962
    etalq    ~ 0.2159                                                                                  # Desai 2016 Table 5 Q CV 49.09%; omega^2 = log(0.4909^2 + 1) = 0.2159
    etalra   ~ 0.1500                                                                                  # Desai 2016 Table 5 RA CV 40.24%; omega^2 = log(0.4024^2 + 1) = 0.1500
    etalgam1 ~ 0.1898                                                                                  # Desai 2016 Table 5 GAM1 CV 45.71%; omega^2 = log(0.4571^2 + 1) = 0.1898

    # Residual error. Desai 2016 used Ln-Ln transformation of observations
    # and predictions with additive residual error on the log scale (W =
    # 44.94%), which corresponds to a proportional error structure in the
    # linear concentration space (per verification-checklist.md).
    propSd <- 0.4494; label("Proportional residual error (fraction) -- Ln-Ln additive error mapped to linear-space proportional")  # Desai 2016 Table 5 theta_8 W = 44.94%
  })

  model({
    # Race effect on CL (multiplicative; reference: Caucasian, RACE_ASIAN = 0).
    race_cl <- 1 + e_race_asian_cl * RACE_ASIAN
    cl <- exp(lcl + etalcl) * race_cl

    vc <- exp(lvc)
    q  <- exp(lq + etalq)

    # SP/DIS_HEALTHY and BMI effects on V_p (multiplicative).
    healthy_vp <- 1 + e_dis_healthy_vp * DIS_HEALTHY
    bmi_vp     <- 1 + e_bmi_vp * (BMI - 24.80)
    vp <- exp(lvp + etalvp) * healthy_vp * bmi_vp

    # Weibull absorption typical values + IIV.
    kamax <- exp(lkamax)
    ra    <- exp(lra + etalra)
    gam1  <- exp(lgam1 + etalgam1)

    # Micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Time-varying first-order absorption rate per Desai 2016:
    # WB = KAMAX * (1 - exp(-(RA * TAD)^GAM1)).
    # tad() returns time since the most recent dose (hours).
    wb <- kamax * (1 - exp(-(ra * tad())^gam1))

    # ODE system (Desai 2016 three differential equations).
    d/dt(depot)       <- -wb * depot
    d/dt(central)     <-  wb * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    # Oral bioavailability (fixed at 1 from prior NCA).
    f(depot) <- exp(lfdepot)

    # Plasma concentration. dose in mg, V_c in L -> mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
