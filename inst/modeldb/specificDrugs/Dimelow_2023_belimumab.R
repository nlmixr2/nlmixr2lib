Dimelow_2023_belimumab <- function() {
  description <- "Belimumab two-compartment population PK (route-specific IV and SC typical-value parameter sets with body weight and BMI effects) driving a circulating / extra-circulatory memory B-cell model with saturable belimumab stimulation of memory B-cell transit into the circulation and of peripheral-pool depletion, in adult and pediatric patients with systemic lupus erythematosus (Dimelow 2023)"
  reference <- paste(
    "Dimelow R, Gillespie WR, van Maurik A. Population model-based analysis",
    "of the memory B-cell response following belimumab therapy in the",
    "treatment of systemic lupus erythematosus. CPT Pharmacometrics Syst",
    "Pharmacol. 2023;12:462-473. doi:10.1002/psp4.12919.",
    "The PK typical values (Dimelow 2023 Table S1) come from Struemper H,",
    "Chen C, Cai W. J Clin Pharmacol. 2013;53:711-720 (IV) and Struemper H,",
    "Thapar M, Roth D. Clin Pharmacokinet. 2018;57:717-728 (SC).",
    sep = " "
  )
  vignette <- "Dimelow_2023_belimumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")
  dosing <- c("depot", "central")

  # Model structure (Dimelow 2023 Figure 2a, the stimulus-based model; the
  # final Model D of Table 1):
  #   dM/dt = KM * MEMBL * ST(C) * P - KM * M,       M(0) = MEMBL
  #   dP/dt = KP - KP * STP(C) * P,                   P(0) = 1
  #   ST(C)  = 1 + Emax  * C / (EC50 + C)
  #   STP(C) = 1 + EmaxP * C / (EC50 + C)
  #   KP = theta * KM  (Dimelow 2023 Results, 'Parameter identifiability')
  # M is the circulating memory B-cell concentration (cells/uL) and P the
  # unobserved extra-circulatory memory B-cell pool, relative to its baseline.
  # C is the belimumab serum concentration (ug/mL) from the PK sub-model.

  compartmentData <- list(
    depot = list(analyte = "belimumab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "belimumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "belimumab", units = "mg", specimen = "tissue", verified = TRUE),
    bcell = list(
      analyte = "CD19+ CD20+ CD27+ memory B cells",
      units = "cells/uL",
      specimen = "blood cell",
      verified = TRUE
    ),
    bcell_periph = list(
      analyte = "memory B cells outside the circulation (lymphoid organs and/or inflamed tissues)",
      units = "fraction of baseline",
      specimen = "tissue",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effects on the belimumab PK only (Dimelow 2023 Table S1): IV",
        "set CL (WT/66.3)^0.506 and V1 (WT/66.3)^1.14; SC set CL and Q",
        "(WT/67)^0.75 and V1 and V2 (WT/67)^1.",
        sep = " "
      ),
      source_name = "WT"
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on the belimumab central volume only (Dimelow 2023 Table S1): IV set (BMI/25.1)^-0.616, SC set (BMI/24.7)^-0.610.",
      source_name = "BMI"
    ),
    ROUTE_IV = list(
      description = "Belimumab administration route indicator: 1 = intravenous, 0 = subcutaneous",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (SC)",
      notes = paste(
        "Selects which typical-value PK parameter set of Dimelow 2023 Table S1",
        "applies: the IV set (Struemper 2013) for ROUTE_IV = 1, the SC set",
        "(Struemper 2018) for ROUTE_IV = 0. Dose into central when",
        "ROUTE_IV = 1 and into depot when ROUTE_IV = 0. Dimelow 2023 derived",
        "each patient's PK profile from the model for the route that patient",
        "received; the memory B-cell model itself does not depend on route.",
        sep = " "
      ),
      source_name = "route (IV / SC)"
    ),
    SBLYS = list(
      description = "Baseline serum B-lymphocyte stimulator (BLyS) concentration",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on the memory B-cell removal rate KM, (BLyS/1.13)^-0.125, reference 1.13 ng/mL = cohort median (Dimelow 2023 Table 2 and Table S2).",
      source_name = "BLyS-BL"
    ),
    ANTIDSDNA_BL = list(
      description = "Baseline serum anti-double-stranded DNA antibody concentration",
      units = "IU/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effects on Emax (+0.0535), EmaxP (+0.192) and the memory B-cell removal rate KM (-0.0756), each as (anti-dsDNA/79)^theta, reference 79 IU/mL = cohort median (Dimelow 2023 Table 2 and Table S2). Must be > 0 (cohort range 1-15001 IU/mL).",
      source_name = "anti-dsDNA-BL"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened against the week-8 simulated memory B-cell response (Dimelow 2023 Table S3, p = 0.028, adjusted R^2 0.002); not retained (p >= 0.01)."
    ),
    C3_BL = list(
      description = "Baseline serum complement C3",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened (Dimelow 2023 Table S3, p = 0.018); not retained. Reported in mg/dL (cohort median 99 mg/dL = 0.99 g/L)."
    ),
    SCORE_SLEDAI = list(
      description = "Baseline SELENA-SLEDAI score",
      units = "(score)",
      type = "continuous",
      notes = "Screened (Dimelow 2023 Table S3, p = 0.097); not retained. Baseline complement C4 (p = 0.189) was also screened and not retained; it has no canonical column."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1880,
    n_studies = 7,
    age_range = "6-77 years (pooled belimumab + placebo dataset)",
    age_median = "39 years",
    weight_range = "17.0-169.5 kg",
    weight_median = "68.5 kg",
    bmi_median = "25.68 kg/m^2",
    disease_state = "systemic lupus erythematosus; all patients starting belimumab for the first time",
    dose_range = "belimumab 1, 4 or 10 mg/kg IV on days 0, 14 and 28 then every 4 weeks; 200 mg SC weekly; 100 mg SC regimens (weekly-to-every-2-weeks, or 3 times a week after a 200 mg SC load) in the open-label BEL112232 study (Dimelow 2023 Figure 1)",
    regions = "multinational, including Northeast Asia (BEL113750) and patients of Black African ancestry (EMBRACE); studies LBSL02, BEL112232, PLUTO (pediatric), BLISS-SC, BEL113750, BLISS-76 and EMBRACE",
    notes = paste(
      "1880 belimumab-treated patients contributed memory B-cell data: 1216",
      "adults IV, 611 adults SC and 53 children aged 5-17 years IV (Dimelow",
      "2023 Results). 873 placebo patients were excluded from the fit.",
      "Baseline (pooled belimumab + placebo, N = 2753; Table S2): BLyS median",
      "1.13 ng/mL (IQR 0.70-1.77), anti-dsDNA median 79 IU/mL (IQR 29-201),",
      "memory B cells median 15 cells/uL (IQR 7-30). Memory B cells were",
      "CD19+ CD20+ CD27+ by flow cytometry.",
      sep = " "
    )
  )

  ini({
    # --- Belimumab PK (Dimelow 2023 Table S1) ---------------------------------
    # Population typical values only: Dimelow 2023 'Pharmacokinetics' states
    # the profiles were derived with the population median PK parameters,
    # accounting for between-patient variability explained by body weight and
    # BMI. The PK was not estimated in this analysis, so every PK value is
    # fixed. Volumes and clearances converted from mL to L.
    lcl_iv <- fixed(log(0.215)); label("IV clearance CL at 66.3 kg (L/day)")                      # Table S1 IV: 215 mL/day
    lq_iv <- fixed(log(0.459)); label("IV inter-compartmental clearance Q (L/day)")               # Table S1 IV: 459 mL/day
    lvc_iv <- fixed(log(2.56)); label("IV central volume V1 at 66.3 kg and BMI 25.1 (L)")       # Table S1 IV: 2560 mL
    lvp_iv <- fixed(log(2.73)); label("IV peripheral volume V2 (L)")                             # Table S1 IV: 2730 mL
    e_wt_cl_iv <- fixed(0.506); label("IV power exponent of body weight on CL (unitless)")       # Table S1 IV: (WT/66.3)^0.506
    e_wt_vc_iv <- fixed(1.14); label("IV power exponent of body weight on V1 (unitless)")        # Table S1 IV: (WT/66.3)^1.14
    e_bmi_vc_iv <- fixed(-0.616); label("IV power exponent of BMI on V1 (unitless)")             # Table S1 IV: (BMI/25.1)^-0.616

    lcl_sc <- fixed(log(0.204)); label("SC clearance CL at 67 kg (L/day)")                        # Table S1 SC: 204 mL/day
    lq_sc <- fixed(log(0.698)); label("SC inter-compartmental clearance Q at 67 kg (L/day)")      # Table S1 SC: 698 mL/day
    lvc_sc <- fixed(log(2.30)); label("SC central volume V1 at 67 kg and BMI 24.7 (L)")         # Table S1 SC: 2300 mL
    lvp_sc <- fixed(log(2.65)); label("SC peripheral volume V2 at 67 kg (L)")                    # Table S1 SC: 2650 mL
    e_wt_cl_sc <- fixed(0.75); label("SC power exponent of body weight on CL (unitless)")        # Table S1 SC: (WT/67)^0.75
    e_wt_q_sc <- fixed(0.75); label("SC power exponent of body weight on Q (unitless)")          # Table S1 SC: (WT/67)^0.75
    e_wt_vc_sc <- fixed(1); label("SC power exponent of body weight on V1 (unitless)")           # Table S1 SC: printed 'WT/67', exponent 1
    e_wt_vp_sc <- fixed(1); label("SC power exponent of body weight on V2 (unitless)")           # Table S1 SC: printed 'WT/67', exponent 1
    e_bmi_vc_sc <- fixed(-0.610); label("SC power exponent of BMI on V1 (unitless)")             # Table S1 SC: (BMI/24.7)^-0.610
    lka <- fixed(log(0.235)); label("SC absorption rate constant (1/day)")                       # Table S1 SC: 0.235 1/day
    lalag <- fixed(log(0.179)); label("SC dose-absorption lag time (day)")                       # Table S1 SC: 0.179 day
    lfdepot <- fixed(log(0.742)); label("SC bioavailability (fraction)")                         # Table S1 SC: 74.2%

    # --- Memory B-cell model (Dimelow 2023 Table 2, Bayesian posterior medians)
    lec50 <- log(2.44); label("Belimumab concentration for half-maximal stimulation EC50 (ug/mL)")              # Table 2: EC50 2.44 ug/ml
    lemax <- log(2.11); label("Maximal stimulation of memory B-cell input to the circulation Emax (unitless)") # Table 2: Emax 2.11
    lemax_bcell_periph <- log(2.91); label("Maximal stimulation of peripheral memory B-cell loss EmaxP (unitless)") # Table 2: EmaxP 2.91
    lkout <- log(0.0631); label("Circulating memory B-cell removal rate KM (1/day)")                           # Table 2: KM 0.0631 day^-1
    lkout_bcell_periph_ratio <- log(0.0237); label("Ratio of peripheral to circulating memory B-cell removal rate KP/KM (unitless)") # Table 2: Ratio KP/KM 0.0237
    lrbase <- log(14.7); label("Baseline circulating memory B-cell concentration MEMBL (cells/uL)")            # Table 2: MEMBL 14.7 cells/ul

    e_antidsdna_emax <- 0.0535; label("Power exponent of baseline anti-dsDNA (/79 IU/mL) on Emax (unitless)")              # Table 2: x (Anti-dsDNA/79)^theta, theta 0.0535
    e_antidsdna_emax_bcell_periph <- 0.192; label("Power exponent of baseline anti-dsDNA (/79 IU/mL) on EmaxP (unitless)") # Table 2: x (Anti-dsDNA/79)^theta, theta 0.192
    e_sblys_kout <- -0.125; label("Power exponent of baseline BLyS (/1.13 ng/mL) on KM (unitless)")                        # Table 2: x (BLyS/1.13)^theta, theta -0.125
    e_antidsdna_kout <- -0.0756; label("Power exponent of baseline anti-dsDNA (/79 IU/mL) on KM (unitless)")               # Table 2: x (Anti-dsDNA/79)^theta, theta -0.0756

    # Between-patient variability: Table 2 reports log-scale SDs (posterior
    # medians); variance = SD^2. No IIV on EC50 or on the KP/KM ratio. The
    # correlation matrix was estimated (LKJ(2) prior) but its values are not
    # reported, so the etas are entered as independent.
    etalemax ~ 0.484416              # Table 2: Emax log-scale SD 0.696; 0.696^2
    etalemax_bcell_periph ~ 0.898704 # Table 2: EmaxP log-scale SD 0.948; 0.948^2
    etalkout ~ 1.3456                # Table 2: KM log-scale SD 1.16; 1.16^2
    etalrbase ~ 0.966289             # Table 2: MEMBL log-scale SD 0.983; 0.983^2

    # Residual error y = f * exp(eps) (Dimelow 2023 Methods 'Random effects'),
    # i.e. additive on the log scale.
    expSd <- 0.348; label("Exponential residual error SD on memory B cells (log scale)") # Table 2: y = f x exp(eps), SD 0.348
  })

  model({
    # 1. Belimumab PK (Dimelow 2023 Table S1); route-specific typical values
    # selected per subject by ROUTE_IV (1 = IV, 0 = SC).
    cl_iv <- exp(lcl_iv) * (WT / 66.3)^e_wt_cl_iv
    q_iv <- exp(lq_iv)
    vc_iv <- exp(lvc_iv) * (WT / 66.3)^e_wt_vc_iv * (BMI / 25.1)^e_bmi_vc_iv
    vp_iv <- exp(lvp_iv)

    cl_sc <- exp(lcl_sc) * (WT / 67)^e_wt_cl_sc
    q_sc <- exp(lq_sc) * (WT / 67)^e_wt_q_sc
    vc_sc <- exp(lvc_sc) * (WT / 67)^e_wt_vc_sc * (BMI / 24.7)^e_bmi_vc_sc
    vp_sc <- exp(lvp_sc) * (WT / 67)^e_wt_vp_sc

    cl <- ROUTE_IV * cl_iv + (1 - ROUTE_IV) * cl_sc
    q <- ROUTE_IV * q_iv + (1 - ROUTE_IV) * q_sc
    vc <- ROUTE_IV * vc_iv + (1 - ROUTE_IV) * vc_sc
    vp <- ROUTE_IV * vp_iv + (1 - ROUTE_IV) * vp_sc
    ka <- exp(lka)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    f(depot) <- exp(lfdepot)
    alag(depot) <- exp(lalag)

    Cc <- central / vc # mg/L = ug/mL

    # 2. Memory B-cell model (Dimelow 2023 Figure 2a and Table 2, Model D)
    ec50 <- exp(lec50)
    emax <- exp(lemax + etalemax) * (ANTIDSDNA_BL / 79)^e_antidsdna_emax
    emax_bcell_periph <- exp(lemax_bcell_periph + etalemax_bcell_periph) *
      (ANTIDSDNA_BL / 79)^e_antidsdna_emax_bcell_periph
    kout <- exp(lkout + etalkout) * (SBLYS / 1.13)^e_sblys_kout *
      (ANTIDSDNA_BL / 79)^e_antidsdna_kout
    kout_bcell_periph <- exp(lkout_bcell_periph_ratio) * kout
    rbase <- exp(lrbase + etalrbase)

    stim <- 1 + emax * Cc / (ec50 + Cc)
    stim_bcell_periph <- 1 + emax_bcell_periph * Cc / (ec50 + Cc)

    d/dt(bcell) <- kout * rbase * stim * bcell_periph - kout * bcell
    d/dt(bcell_periph) <- kout_bcell_periph - kout_bcell_periph * stim_bcell_periph * bcell_periph
    bcell(0) <- rbase
    bcell_periph(0) <- 1

    # 3. Observation: circulating memory B cells (cells/uL), y = f * exp(eps)
    bcell ~ lnorm(expSd)
  })
}
