Wang_2015_rucaparib <- function() {
  description <- "Three-compartment IV population PK model coupled to a direct-effect Emax PK/PD model for inhibition of poly(ADP-ribose) polymerase (PARP-1) activity in peripheral blood lymphocytes (PBL) by rucaparib (AG-014699 / PF-01367338) in adult cancer patients (Wang 2015 Phase 1 study A4991002), with a power covariate effect of baseline PBL PARP activity on the residual maximum-inhibition parameter Emin."
  reference <- paste(
    "Wang DD, Li C, Sun W, Zhang S, Shalinsky DR, Kern KA, Curtin NJ,",
    "Sam W-J, Kirkpatrick TR, Plummer R.",
    "PARP activity in peripheral blood lymphocytes as a predictive",
    "biomarker for PARP inhibition in tumor tissues -- a population",
    "pharmacokinetic/pharmacodynamic analysis of rucaparib.",
    "Clin Pharmacol Drug Dev. 2015;4(2):89-98. doi:10.1002/cpdd.176."
  )
  vignette <- "Wang_2015_rucaparib"
  units <- list(time = "hr", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    BLPARP = list(
      description        = "Subject-specific pre-dose baseline PARP-1 activity in peripheral blood lymphocytes",
      units              = "pmol/10^6 PBL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Per-subject baseline biomarker; appears in the residual maximum-inhibition equation",
        "Emin = TV(Emin) * (BLPARP / BLB_median)^alpha (Wang 2015 equation 2). The numeric BLB_median",
        "is not published; the model uses 90.8 pmol/10^6 PBL (= TV(E0) from Wang 2015 Table 2)",
        "as a defensible stand-in for the cohort median because Wang 2015 reports only the typical",
        "E0 and the exponent value but not the actual centring constant. For typical-value",
        "simulation, set BLPARP to 90.8 so Emin reduces to TV(Emin)."
      ),
      source_name        = "BLB"
    )
  )

  population <- list(
    n_subjects     = 32L,                                        # Wang 2015 Methods (Subjects and Study Design): "A total of 32 patients participated in the study"
    n_studies      = 1L,                                         # Single Phase 1 first-in-patient study A4991002
    age_range      = "not reported in main text",                # Demographic table not in main text; Methods lists demographics tested as covariates but no summary table
    weight_range   = "not reported in main text",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = "Adults with advanced solid tumors (Part 1: any tumor type meeting eligibility criteria; Part 2: metastatic malignant melanoma) participating in a Phase 1 first-in-patient dose-escalation study A4991002.",
    dose_range     = "Rucaparib 1, 2, 4, 8, 12, or 18 mg/m^2 IV over 30 min. Part 1 used escalating single doses up to the PARP Inhibitory Dose (PID, 12 mg/m^2). Part 2 administered the PID and 18 mg/m^2 once daily for 5 days each 28-day cycle alongside escalating temozolomide (TMZ) doses (135, 170, 200 mg/m^2/day).",
    regions        = "Multicentre clinical trial; specific regions not reported in main text.",
    co_medication  = "Temozolomide co-administration in combination cohorts; TMZ had no evident effect on rucaparib PK so PK data with and without TMZ were pooled (Methods: Population PK and PK/PD Analyses).",
    n_pk_samples   = 1022L,                                      # Wang 2015 Results: "A total of 1022 rucaparib PK samples obtained from the 26 patients who had valid PK data"
    n_subjects_pk  = 26L,                                        # 26 of 32 patients with valid PK data (Results)
    n_pd_samples_pbl = 348L,                                     # Wang 2015 Results: "A total of 348 PARP activity data values in PBL from all 32 patients participated in the study"
    n_subjects_pd_pbl = 32L,
    n_pd_samples_tumor = 30L,                                    # Wang 2015 Results: "only 30 data points were available for analysis" (tumor PD)
    n_subjects_pd_tumor = 15L,                                   # 14 patients in Part 2 plus one Part 1 patient at 4 mg/m^2 (Results)
    pk_sampling    = "Plasma sampled pre-infusion, 0.25 and 0.5 h after start of infusion, and 0.25, 0.5, 1, 2, 4, 6, 8, 24 h after end of infusion on Day -7 (single-agent), Day 1, and Day 4 of the first cycle.",
    pd_sampling_pbl = "PBL PARP activity sampled pre-dose, end-of-infusion, 4-6 h and 24 h after end of infusion on Days -7, 1, 4 of cycle 1 plus an additional Day 8 sample (3 days after last dose) for the duration-of-inhibition window.",
    pd_sampling_tumor = "Tumor biopsies were collected at baseline and 4-6 h or 24 h after treatment with rucaparib (Day 1) in Part 2 patients only.",
    notes          = "S-PLUS 7.0 / NONMEM 7.1.2 (FOCE) used for analysis. M6 BQL handling per Ahn 2008 (12% of PK samples below LLOQ 2 ng/mL). Demographic / physiological covariates tested (age, gender, weight, body surface area, serum creatinine, AST, ALT, disease stage, PAR baseline in PBL, PAR baseline in tumor) -- only PAR baseline in PBL (BLPARP) and PAR baseline in tumor were retained. None of the demographic covariates entered the final model."
  )

  ini({
    # Three-compartment IV PK model -- Wang 2015 Table 1 (CL, V1) and Discussion
    # text page 95 (V2, V3, Q2, Q3 derived in the comparison-of-distribution
    # paragraph: "the typical volume of distribution of rucaparib in peripheral
    # compartment (311 L for V2 and 48 L for V3) was much larger than that in
    # the central compartment (15.5 L)... the high intercompartmental
    # clearances (21.7 L/hr for Q2, and 52.9 L/hr for Q3)").
    lcl  <- log(17.5);  label("Total clearance CL (L/h)")                          # Wang 2015 Table 1: CL = 17.5 L/hr (RSE 6.57%)
    lvc  <- log(15.5);  label("Central volume of distribution V1 (L)")             # Wang 2015 Table 1: V1 = 15.5 L (RSE 11.6%)
    lq   <- log(21.7);  label("Distributional clearance to peripheral1 Q2 (L/h)")  # Wang 2015 Discussion p.95: Q2 = 21.7 L/hr
    lvp  <- log(311);   label("Peripheral volume 1 V2 (L)")                        # Wang 2015 Discussion p.95: V2 = 311 L
    lq2  <- log(52.9);  label("Distributional clearance to peripheral2 Q3 (L/h)")  # Wang 2015 Discussion p.95: Q3 = 52.9 L/hr
    lvp2 <- log(48);    label("Peripheral volume 2 V3 (L)")                        # Wang 2015 Discussion p.95: V3 = 48 L

    # PD parameters -- direct-effect Emax model for PARP activity in PBL
    # (Wang 2015 equation 1: E = E0 - (E0 - Emin) * C / (IC50 + C))
    le0   <- log(90.8); label("Baseline PARP activity in PBL E0 (pmol/10^6 PBL)") # Wang 2015 Table 2: E0 = 90.8 pmol/10^6 PBL (RSE 19.7%)
    lemin <- log(8.24); label("Typical residual PARP activity at maximum inhibition Emin (pmol/10^6 PBL)") # Wang 2015 Table 2: TV(Emin) = 8.24 pmol/10^6 PBL (RSE 10.6%)
    lic50 <- log(1.05); label("Concentration producing 50% maximum inhibition IC50 (ng/mL)") # Wang 2015 Table 2: IC50 = 1.05 ng/mL (RSE 24.0%)

    # Covariate effect -- power exponent on (BLPARP / BLB_median) for Emin
    # (Wang 2015 equation 2: Emin = TV(Emin) * (BLB / BLB_median)^alpha)
    e_blparp_emin <- 0.620; label("Power exponent on (BLPARP / 90.8) ratio for Emin (unitless)") # Wang 2015 Table 2: alpha1 = 0.620 (RSE 14.2%)

    # IIV -- Wang 2015 Table 1 reports CV% on the log-normal scale; omega^2 = log(1 + CV^2)
    # The paper text states "IIV on all PK parameters" but Table 1 reports IIV
    # CV% only for CL and V1 (the supplementary Table 1 with IIV for V2/V3/Q2/Q3
    # was not on disk for this extraction). IIV on the additional PK parameters
    # is therefore omitted -- see vignette Assumptions and deviations.
    etalcl ~ 0.2329  # log(1 + 0.512^2); Wang 2015 Table 1: IIV CL = 51.2 %CV
    etalvc ~ 0.2840  # log(1 + 0.573^2); Wang 2015 Table 1: IIV V1 = 57.3 %CV

    # IIV on PD parameters -- Wang 2015 Table 2; assumed diagonal (paper text
    # states covariance was estimated only when strong correlation was observed
    # on diagnostic plots, and Table 2 reports no covariance terms).
    etale0   ~ 0.8528  # log(1 + 1.16^2);  Wang 2015 Table 2: IIV E0   = 116  %CV
    etalemin ~ 0.1990  # log(1 + 0.469^2); Wang 2015 Table 2: IIV Emin = 46.9 %CV
    etalic50 ~ 0.3229  # log(1 + 0.617^2); Wang 2015 Table 2: IIV IC50 = 61.7 %CV

    # Residual error -- Wang 2015 Methods: "residual error for the
    # log-transformed PK and PD data was assumed to follow a normal
    # distribution and was modeled using additive structure". Additive on
    # log-scale equals proportional in linear space.
    propSd         <- 0.208; label("Proportional residual error on plasma rucaparib (fraction)") # Wang 2015 Table 1: residual SD = 0.208 (RSE 18.1%) on log-scale
    propSd_parpPbl <- 0.529; label("Proportional residual error on PBL PARP activity (fraction)") # Wang 2015 Table 2: residual SD = 0.529 (RSE 12.5%) on log-scale
  })

  model({
    # Individual structural PK parameters
    cl  <- exp(lcl + etalcl)
    vc  <- exp(lvc + etalvc)
    q   <- exp(lq)
    vp  <- exp(lvp)
    q2  <- exp(lq2)
    vp2 <- exp(lvp2)

    # Three-compartment IV ODE system. Dose enters `central`. Concentration in
    # the central compartment in mg/L is multiplied by 1000 to express as
    # ng/mL (= ug/L), the bioanalytical assay output and publication unit.
    kel <- cl  / vc
    k12 <- q   / vc
    k21 <- q   / vp
    k13 <- q2  / vc
    k31 <- q2  / vp2

    d/dt(central)     <- -kel * central -
                          k12 * central + k21 * peripheral1 -
                          k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    Cc <- 1000 * central / vc

    # Individual PD parameters
    e0      <- exp(le0   + etale0)
    emin_tv <- exp(lemin + etalemin)
    ic50    <- exp(lic50 + etalic50)

    # Power covariate effect of baseline PBL PARP activity on Emin. The paper's
    # cohort median BLB_median is not published; 90.8 pmol/10^6 PBL is used as a
    # defensible stand-in (TV(E0) from Wang 2015 Table 2); for typical-value
    # simulation set BLPARP = 90.8 so emin = emin_tv.
    emin <- emin_tv * (BLPARP / 90.8)^e_blparp_emin

    # Direct-effect Emax inhibition (Wang 2015 equation 1)
    parpPbl <- e0 - (e0 - emin) * Cc / (ic50 + Cc)

    Cc      ~ prop(propSd)
    parpPbl ~ prop(propSd_parpPbl)
  })
}
