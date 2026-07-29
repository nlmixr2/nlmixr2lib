Desai_2016_isavuconazole <- function() {
  description <- "Two-compartment population PK model for isavuconazole (administered as the prodrug isavuconazonium sulfate) in healthy adults and adults with mild (Child-Pugh A) or moderate (Child-Pugh B) hepatic impairment, following single 100 mg oral or 2-h intravenous doses (Desai 2016). Weibull absorption for the oral route; hepatic-impairment-group-specific typical CL and Q; linear BMI effect on peripheral volume."
  reference <- "Desai A, Schmitt-Hoffmann A-H, Mujais S, Townsend R. Population Pharmacokinetics of Isavuconazole in Subjects with Mild or Moderate Hepatic Impairment. Antimicrobial Agents and Chemotherapy. 2016;60(5):3025-3031. doi:10.1128/AAC.02942-15"
  vignette <- "Desai_2016_isavuconazole"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    HEPIMP_MILD = list(
      description        = "Mild hepatic impairment indicator (1 = Child-Pugh Class A, 0 = healthy or moderate)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy or moderate; paired with HEPIMP_MOD so all-zero on both indicators selects the healthy reference)",
      notes              = "Classification scheme: Child-Pugh Class A (composite score 5-6 across bilirubin, albumin, prothrombin time). Desai 2016 used Child-Pugh A for the mild stratum (Methods 'Subjects'; Table 3 reports mild Child-Pugh scores 5.18-5.75 across the two source studies). Renamed from the paper's liver-function index HEP1 to canonical HEPIMP_MILD per inst/references/covariate-columns.md.",
      source_name        = "HEP1"
    ),
    HEPIMP_MOD = list(
      description        = "Moderate hepatic impairment indicator (1 = Child-Pugh Class B, 0 = healthy or mild)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy or mild; paired with HEPIMP_MILD so all-zero on both indicators selects the healthy reference)",
      notes              = "Classification scheme: Child-Pugh Class B (composite score 7-9 across bilirubin, albumin, prothrombin time). Desai 2016 used Child-Pugh B for the moderate stratum (Methods 'Subjects'; Table 3 reports moderate Child-Pugh scores 7.43-8.31 across the two source studies). Renamed from the paper's liver-function index HEP2 to canonical HEPIMP_MOD per inst/references/covariate-columns.md.",
      source_name        = "HEP2"
    ),
    BMI = list(
      description        = "Body mass index at baseline",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference BMI 27.00 kg/m^2 (Desai 2016 best-covariate model equation on p. 3028; corresponds to the median of the pooled cohort). Linear-deviation effect on Vp: Vp = Vp_pop * (1 + e_bmi_vp * (BMI - 27)).",
      source_name        = "BMI"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 96,
    n_studies       = 2,
    age_range       = "37-64 years (median 50 healthy, 54 mild, 54 moderate per Table 2)",
    weight_range    = "53-107 kg (median 78 healthy, 76 mild, 76 moderate per Table 2)",
    bmi_range       = "21-34 kg/m^2 (median 27 healthy, 28 mild, 26 moderate per Table 2)",
    sex_female_pct  = 34,
    disease_state   = "Healthy controls and adults with mild (Child-Pugh A) or moderate (Child-Pugh B) hepatic impairment due to either alcoholic cirrhosis (Study 1, 48 subjects) or chronic hepatitis B/C (Study 2, 48 subjects). Equal n=32 in each of three hepatic-impairment strata.",
    hepatic_function = "Child-Pugh A (mild; mean composite score 5.18 alcoholic cirrhosis cohort / 5.75 hepatitis B/C cohort) and Child-Pugh B (moderate; mean 7.43 / 8.31). Healthy controls matched to impaired subjects on age (within +-7 years), sex, body weight (within +-8 kg), and BMI (within +-4 kg/m^2). See Table 3.",
    dose_range      = "Single 100 mg isavuconazole equivalent (186 mg isavuconazonium sulfate prodrug) administered either orally (p.o.) or as a 2-h continuous intravenous infusion.",
    smoking_status  = "Healthy 59% smokers; mild 75%; moderate 75% (Table 2).",
    notes           = "Data pooled across two Phase 1 single-dose hepatic-impairment studies for combined popPK analysis. Plasma sampling 0.5 to 480 h post-dose. 2,016 total isavuconazole plasma concentrations. Estimation in NONMEM 7.2 with FOCE method. Severe hepatic impairment (Child-Pugh C) was not studied."
  )

  ini({
    # Structural typical values for the reference (healthy) subject -- Desai 2016
    # Table 4 'Parameter estimates of the best covariate model'. Linear-scale
    # values in the table are in ml/h (CL, Q) and ml (V2, V3); converted to L/h
    # and L for nlmixr2lib unit consistency (time = h, concentration = mg/L).
    lcl <- log(2.54);    label("Isavuconazole clearance for the healthy reference subject (CL, L/h)")  # Desai 2016 Table 4: theta_8 (CL, healthy) = 2540 mL/h
    lvc <- log(51.4);    label("Isavuconazole central volume of distribution (V2 = Vc, L)")  # Desai 2016 Table 4: theta_2 (V2) = 51400 mL
    lq  <- log(33.678);  label("Isavuconazole intercompartmental clearance for the healthy reference subject (Q, L/h)")  # Desai 2016 Table 4: theta_10 (Q, healthy) = 33678 mL/h
    lvp <- log(410);     label("Isavuconazole peripheral volume of distribution at reference BMI 27 (V3 = Vp, L)")  # Desai 2016 Table 4: theta_4 (V3) = 410000 mL (printed in the table as '41,0000' with the comma misplaced; bootstrap mean 410661 mL confirms 410 L)

    # Weibull absorption parameters (Piotrovskij saturating-ka form, operator-
    # approved interpretation -- sidecar request-001 q3 = A). The paper text
    # names the three parameters but does not write the equation; the canonical
    # form for these three named parameters is
    #   ka(t) = kamax * (1 - exp(-(ra * tad)^gam1))
    # See Assumptions and deviations in the vignette.
    lra    <- log(0.653);  label("Weibull-absorption rate-scaling parameter (RA, 1/h)")  # Desai 2016 Table 4: theta_5 (RA) = 0.653 1/h
    lgam1  <- log(4.57);   label("Weibull-absorption shape / sigmoidicity parameter (GAM1, unitless)")  # Desai 2016 Table 4: theta_6 (GAM1) = 4.57
    lkamax <- log(0.86);   label("Weibull-absorption asymptotic maximum absorption rate constant (KAMAX, 1/h)")  # Desai 2016 Table 4: theta_7 (KAMAX) = 0.86 1/h

    # Bioavailability anchor. Desai 2016 Table 4 reports F1 = 1.00 (fixed). The
    # 100 mg dose is the isavuconazole equivalent of the 186 mg prodrug; the
    # paper assumes complete conversion of isavuconazonium sulfate to active
    # isavuconazole in vivo.
    lfdepot <- fixed(log(1.0));  label("Isavuconazole oral bioavailability (F1, fraction; fixed)")  # Desai 2016 Table 4: F1 = 1.00 (fixed)

    # Categorical covariate effects on CL and Q for hepatic-impairment status.
    # Desai 2016 reports three group-specific typical values (theta_{1,8,9} for
    # CL and theta_{3,10,11} for Q); the log-additive shift encoding below is
    # mathematically equivalent to the paper's three-typical-value form:
    #   exp(lcl + e_hepimp_mild_cl) = 2.54 * exp(log(1.55/2.54)) = 1.55 L/h (mild)
    #   exp(lcl + e_hepimp_mod_cl)  = 2.54 * exp(log(1.326/2.54)) = 1.326 L/h (moderate)
    # Same for Q. The healthy reference (HEPIMP_MILD=HEPIMP_MOD=0) recovers
    # lcl / lq directly.
    e_hepimp_mild_cl <- log(1.55  / 2.54);    label("Log-additive shift on CL for mild Child-Pugh A vs healthy reference (unitless)")  # Desai 2016 Table 4: theta_1 (CL, mild) = 1550 mL/h relative to theta_8 (CL, healthy) = 2540 mL/h
    e_hepimp_mod_cl  <- log(1.326 / 2.54);    label("Log-additive shift on CL for moderate Child-Pugh B vs healthy reference (unitless)")  # Desai 2016 Table 4: theta_9 (CL, moderate) = 1326 mL/h relative to theta_8 (CL, healthy) = 2540 mL/h
    e_hepimp_mild_q  <- log(38.8  / 33.678);  label("Log-additive shift on Q for mild Child-Pugh A vs healthy reference (unitless)")    # Desai 2016 Table 4: theta_3 (Q, mild) = 38800 mL/h relative to theta_10 (Q, healthy) = 33678 mL/h
    e_hepimp_mod_q   <- log(63.554/ 33.678);  label("Log-additive shift on Q for moderate Child-Pugh B vs healthy reference (unitless)") # Desai 2016 Table 4: theta_11 (Q, moderate) = 63554 mL/h relative to theta_10 (Q, healthy) = 33678 mL/h

    # Linear-deviation BMI effect on V3 (Vp). Desai 2016 best-covariate model on
    # p. 3028: V3 = theta_4 * [1 + theta_12 * (BMI - 27.00)].
    e_bmi_vp <- 0.058; label("Linear-deviation effect of BMI on Vp, Vp = Vp_pop * (1 + e_bmi_vp * (BMI - 27)) (per kg/m^2)")  # Desai 2016 Table 4: theta_12 (BMI on V3) = 0.058 per kg/m^2

    # Inter-individual variability. Desai 2016 Table 4 reports the IIVs as
    # percent CV (linear-space); the log-normal omega^2 is computed as
    # log(1 + (CV/100)^2). The shared eta on CL applies to all three hepatic-
    # impairment strata (paper's exp(eta_jCL) form); same for Q.
    etalcl    ~ 0.17315  # 43.47% CV; omega^2 = log(1 + 0.4347^2). Desai 2016 Table 4 Variability section
    etalvc    ~ 0.04408  # 21.23% CV; omega^2 = log(1 + 0.2123^2). Desai 2016 Table 4 Variability section
    etalra    ~ 0.10070  # 32.55% CV; omega^2 = log(1 + 0.3255^2). Desai 2016 Table 4 Variability section
    etalkamax ~ 0.09621  # 31.78% CV; omega^2 = log(1 + 0.3178^2). Desai 2016 Table 4 Variability section
    etalgam1  ~ 0.13534  # 38.07% CV; omega^2 = log(1 + 0.3807^2). Desai 2016 Table 4 Variability section
    etalvp    ~ 0.07061  # 27.05% CV; omega^2 = log(1 + 0.2705^2). Desai 2016 Table 4 Variability section
    etalq     ~ 0.12480  # 36.46% CV; omega^2 = log(1 + 0.3646^2). Desai 2016 Table 4 Variability section

    # Residual error. Desai 2016 Methods 'Structural pharmacokinetic model'
    # states that ln-ln transformations of both the model and the data were
    # used and the residual variance was modeled as additive on the log scale,
    # which corresponds to proportional residual error in linear space for
    # nlmixr2. The Table 4 value reported as 'sigma^2 = 17.88' is the proportional
    # SD on the natural scale (17.88% CV), confirmed by the SE/RSE pair
    # 0.0029 / 9% which gives RSE on sigma^2 itself (0.0029 / 0.1788^2 ~= 9%).
    propSd <- 0.1788; label("Isavuconazole proportional residual error (fraction)")  # Desai 2016 Table 4: Residual error sigma^2 = 17.88
  })

  model({
    # Individual structural parameters. Hepatic-impairment indicators apply a
    # log-additive shift to the reference (healthy) typical CL and Q. The IIV
    # is applied to the group-specific log typical value, mirroring the paper's
    # CL = theta_g * exp(eta_jCL) form.
    cl <- exp(lcl + e_hepimp_mild_cl * HEPIMP_MILD + e_hepimp_mod_cl * HEPIMP_MOD + etalcl)
    vc <- exp(lvc + etalvc)
    q  <- exp(lq  + e_hepimp_mild_q  * HEPIMP_MILD + e_hepimp_mod_q  * HEPIMP_MOD + etalq)
    vp <- exp(lvp + etalvp) * (1 + e_bmi_vp * (BMI - 27))

    # Weibull absorption (Piotrovskij saturating-ka form, operator-approved
    # interpretation -- sidecar request-001 q3 = A). Time after most recent
    # dose to the depot compartment drives the time-varying absorption rate;
    # tad(depot) restarts at each new oral dose so multi-dose regimens
    # (e.g., 200 mg q8h loading followed by 200 mg q24h maintenance per the
    # paper's clinical-dose simulation) reproduce the Weibull rise from each
    # individual dose. The if-gate handles IV-only subjects (no PO dose) where
    # tad(depot) is NaN before any dose and would otherwise propagate as NaN
    # through the entire ODE; same pattern as Cirincione_2017_exenatide.R and
    # Horita_2018_rifampicin.R.
    ra    <- exp(lra    + etalra)
    gam1  <- exp(lgam1  + etalgam1)
    kamax <- exp(lkamax + etalkamax)
    ka <- 0.0
    if (tad(depot) >= 0.0) ka <- kamax * (1 - exp(-(ra * tad(depot))^gam1))

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Oral bioavailability anchor (F1 = 1.00 fixed). IV doses bypass the depot
    # entirely by routing the amt event to the central compartment via the cmt
    # column in the event table.
    f(depot) <- exp(lfdepot)

    # Concentration in mg/L (= ug/mL; Desai 2016 reports concentrations in
    # ng/mL in the trough-concentration summary table, requiring a 1e3 scale
    # factor for any direct comparison). Dose mg / volume L gives mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
