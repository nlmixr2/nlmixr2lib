Kir_2025_atenolol_rat_pbpk <- function() {
  description <- "Preclinical (rat, Sprague-Dawley). mPBPK (minimal physiologically-based, Monolix 2023R1 population fit). Atenolol (ATN) oral absorption in non-malnourished (control) and malnourished rats. Four mass-balance states: a well-mixed blood pool, two lumped tissue groups (Tissue 1 = rapidly perfused liver + lung, high Kp; Tissue 2 = slowly perfused remainder) and the kidney, which is the eliminating organ (ATN is hydrophilic, logP 0.16, and is cleared almost entirely renally). Distribution is PERMEABILITY-limited: the fraction of cardiac output actually equilibrating with each tissue group is small and estimated (fd1 = fd2 = 0.134, so only ~27% of cardiac output is distributive), which is what encodes atenolol's BCS Class III low permeability. Oral input is NOT first-order: the paper's point-area deconvolution showed a near-constant absorption rate over successive time windows, so absorption is modelled as two (control) or three (malnourished) SEQUENTIAL ZERO-ORDER processes with author-fixed window boundaries (0-120, 120-300, and 300-420 min). Malnutrition acts on ABSORPTION only -- disposition (fd1, Kp1, CL) is shared across groups, and the malnourished group gains a third absorption phase (k03 = 1.16 mg/min/kg, fixed to 0 in controls by goodness of fit) that raises apparent bioavailability from 0.43 to 0.67. Blood-to-plasma ratio differs by group (1.000 vs 1.014) because malnutrition changes haematocrit and albumin binding. The model as shipped simulates the ORAL arms; the literature IV reference profile is reproduced by setting k01/k02/k03 to 0 and giving a bolus into a_blood (see the vignette). IIV is on clearance only. NOTE: the absorption-rate scale in Tables 2/3 required reconciliation against the deposited Monolix code -- see the vignette Errata."
  reference <- paste(
    "Kir F, Sahin S, Jusko WJ. (2025).",
    "Minimal Physiologically-Based Pharmacokinetic Modeling of Atenolol and",
    "Metoprolol Absorption in Malnourished Rats.",
    "Eur J Drug Metab Pharmacokinet 50:243-255.",
    "doi:10.1007/s13318-025-00943-6. PMCID PMC12081501.",
    "Population parameter estimates: Table 2 (Monolix). Naive-pooled (ADAPT 5)",
    "estimates for the same model: Supplementary Table S3.",
    "Rat physiology (tissue volumes and blood flows): Supplementary Table S2.",
    "Model equations: Article Equations 1-5 (ATN), 9 (tissue-volume closure),",
    "10 (permeability-limited fd constraint) and 12 (blood-to-plasma ratio).",
    "ODEs, fixed system constants and the absorption-window logic were taken from the",
    "author-deposited Monolix (MLXTRAN) and ADAPT-5 source listings in the",
    "Supplementary Materials, which agree with each other exactly.",
    sep = " "
  )
  vignette <- "Kir_2025_atenolol_metoprolol_malnutrition"

  units <- list(
    time          = "min",
    dosing        = "ug",
    concentration = "ug/mL",
    amount        = "ug",
    weight        = "kg"
  )

  # Every volume, flow, amount and clearance in this model is normalised PER KG
  # of body weight (mL/kg, mL/min/kg, ug/kg), exactly as the deposited code is
  # written. A simulated "amount" is therefore ug per kg body weight.

  covariateData <- list(
    MAL_NOURISH = list(
      description        = "Malnutrition status at study entry",
      units              = "unitless",
      type               = "binary",
      reference_category = "0 (control, non-malnourished)",
      notes              = paste(
        "1 = malnourished, produced experimentally by feeding a 5% protein isocaloric diet",
        "for 17-20 days; 0 = control, fed the 20% protein isocaloric diet (Methods 2.2.1).",
        "Malnutrition was confirmed biochemically in Table 1: body weight 216 g vs 300 g",
        "(p < 0.05), serum albumin 3.75 vs 4.40 g/dL and total cholesterol 50.3 vs",
        "78.4 mg/dL (both p < 0.001).",
        "In this model MAL_NOURISH switches (a) the zero-order absorption rate of every",
        "phase, (b) the presence of the third absorption phase and its end time, (c) the",
        "blood-to-plasma ratio, and (d) the proportional residual error magnitude.",
        "It deliberately does NOT change any disposition parameter: fd1, Kp1 and CL were",
        "shared across groups in the published fit (Results 3.1)."
      ),
      source_name        = "Group"
    )
  )

  # Screened but NOT retained. Results 3: "Assessment of these measures as
  # covariates did not improve the population modeling." Body weight, serum
  # albumin and total cholesterol are reported per group in Table 1 and all
  # differ significantly between groups, but none entered the final model --
  # the malnutrition effect is carried by MAL_NOURISH alone.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Table 1: 300 (SD 18.27) g control vs 216 (SD 14.99) g malnourished; screened, not retained."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/dL",
      type        = "continuous",
      notes       = paste(
        "Table 1: 4.40 (SD 0.07) g/dL control vs 3.75 (SD 0.17) g/dL malnourished.",
        "Not retained as a covariate, but albumin IS used indirectly: the fraction unbound",
        "was rescaled from the control value by the albumin ratio (fu 0.970 -> 1.00),",
        "which then enters the blood-to-plasma ratio via Equation 12."
      )
    ),
    TCHOL = list(
      description = "Total cholesterol",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Table 1: 78.4 (SD 15.01) mg/dL control vs 50.3 (SD 8.34) mg/dL malnourished; screened, not retained."
    )
  )

  compartmentData <- list(
    a_blood            = list(analyte = "atenolol", units = "ug", specimen = "whole blood", verified = TRUE),
    a_rapidly_perfused = list(analyte = "atenolol", units = "ug", specimen = "tissue", verified = TRUE),
    a_slowly_perfused  = list(analyte = "atenolol", units = "ug", specimen = "tissue", verified = TRUE),
    a_kidney           = list(analyte = "atenolol", units = "ug", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species        = "rat (Sprague-Dawley)",
    n_subjects     = 8L,
    n_studies      = 2L,
    age_range      = NA_character_,
    weight_range   = "216 (SD 14.99) g malnourished to 300 (SD 18.27) g control",
    sex_female_pct = 0,
    disease_state  = "Experimental protein-calorie malnutrition (5% protein isocaloric diet, 17-20 days) vs control (20% protein isocaloric diet)",
    dose_range     = "250 mg/kg single oral dose (suspension, feeding tube); literature IV reference 1 mg/kg",
    regions        = "Turkey (Kobay Experimental Animals Laboratory, Ankara)",
    notes          = paste(
      "n = 4 control + 4 malnourished male Sprague-Dawley rats (Methods 2.2.1); the oral",
      "arms are the only rats actually dosed in this study. Blood was sampled from the tail",
      "vein at 0, 30, 60, 90, 120, 180, 240, 300, 360, 420, 450 and 480 min.",
      "The IV reference profile that informs disposition was NOT generated here: it was",
      "digitised from the literature (1 mg/kg IV in male SD rats, reference 18 of the",
      "paper), and the per-timepoint standard deviations were digitised and treated as",
      "3 additional rats (Methods 2.2.6). External validation used a further literature",
      "IV dataset at 1.67 mg/kg (reference 19).",
      "Rats were fasted overnight before dosing and fed ad libitum otherwise."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Rat physiological system parameters, per kg body weight.
    # Volumes: Supplementary Table S2 (lung 3.27, liver 33.81, kidney 6.16,
    # blood 58.82 mL/kg). Tissue 1 = liver + lung = 33.81 + 3.27 = 37.08.
    # Tissue 2 closes the body volume by Equation 9 (Body weight = Vb + V1 +
    # V2 + V3): 1000 - 58.82 - 37.08 - 6.16 = 897.94. All four values are
    # written literally in the deposited Monolix and ADAPT listings.
    # ------------------------------------------------------------------
    v_blood            <- fixed(58.82)  ; label("Blood volume (mL/kg)")                                          # Supp Table S2 (Shah & Betts 2012); Monolix/ADAPT Vb = 58.82
    v_rapidly_perfused <- fixed(37.08)  ; label("Tissue 1 volume, rapidly perfused = liver + lung (mL/kg)")      # Supp Table S2 liver 33.81 + lung 3.27; Monolix/ADAPT V1 = 37.08
    v_slowly_perfused  <- fixed(897.94) ; label("Tissue 2 volume, slowly perfused remainder (mL/kg)")            # Equation 9 closure; Monolix/ADAPT V2 = 897.94
    v_kidney           <- fixed(6.16)   ; label("Kidney volume (mL/kg)")                                         # Supp Table S2; Monolix/ADAPT Vkid = 6.16
    q_co               <- fixed(179.91) ; label("Cardiac output (mL/min/kg)")                                    # Supp Table S2; Monolix/ADAPT Qco = 179.91 (beta-blocker-reduced, Methods 2.2.3)
    q_kidney           <- fixed(29.00)  ; label("Renal blood flow (mL/min/kg)")                                  # Supp Table S2; Monolix/ADAPT Qkid = 29.00

    # Kidney partition coefficient was PREDICTED with GastroPlus 9.9, not
    # fitted (Methods 2.2.2 quotes liver 2.87, kidney 3, lung 2.63).
    kp_kidney          <- fixed(3.0)    ; label("Kidney-to-plasma partition coefficient (unitless)")             # Methods 2.2.2 (GastroPlus 9.9 prediction); Monolix/ADAPT Kpkid = 3.0

    # Blood-to-plasma ratio, Equation 12: rho = (HCT - 1 + Rb)/(HCT * fu),
    # i.e. Rb = 1 - HCT + rho*HCT*fu. HCT 47.5% control / 45.5% malnourished
    # (Methods 2.2.3, reference 30); fu 0.970 control / 1.00 malnourished.
    # Both reproduce the paper's quoted Rb of 1.000 and 1.014 exactly.
    bpr_control        <- fixed(1.000)  ; label("Blood-to-plasma ratio, control (unitless)")                     # Methods 2.2.3 "Rb values of ATN for the control and malnourished groups were 1.00 and 1.014"; Monolix Rb_H
    bpr_malnourished   <- fixed(1.014)  ; label("Blood-to-plasma ratio, malnourished (unitless)")                # Methods 2.2.3; ADAPT Rb_M = 1.014

    # Study design constant. The oral input rate scales with this value
    # because the deposited code writes input = k0 * Dose (see model()).
    dose_po            <- fixed(250)    ; label("Oral atenolol dose (mg/kg)")                                    # Methods 2.2.1: "ATN was suspended in purified water at a dose of 250 mg/kg"; Monolix input_H = k01*250000 (ug/kg)

    # ------------------------------------------------------------------
    # Absorption-window boundaries. All three were FIXED by the authors, not
    # estimated: t1 "determined based on deconvolution and trial and error,
    # and then fixed to 120 min for both groups"; tmax "fixed at 300 min for
    # both control and malnourished groups according to visual inspection"
    # (Results 3.1). The third window ends at 360 min in controls and 420 min
    # in the malnourished group (Table 2 phase-3 row headings; ADAPT tmal=420).
    # ------------------------------------------------------------------
    t_abs1             <- fixed(120)    ; label("End of absorption phase 1 (min)")                               # Results 3.1 "fixed to 120 min for both groups"; Monolix/ADAPT tfirst = 120
    t_abs2             <- fixed(300)    ; label("End of absorption phase 2 (min)")                               # Results 3.1 "t max was fixed at 300 min"; Monolix/ADAPT tmax = 300
    t_abs3_control     <- fixed(360)    ; label("End of absorption phase 3, control (min)")                      # Table 2 k03 C row: "t = 300-360 min"
    t_abs3_malnourished <- fixed(420)   ; label("End of absorption phase 3, malnourished (min)")                 # Table 2 k03 M row: "t = 300-420 min"; ADAPT tmal = 420

    # ------------------------------------------------------------------
    # Estimated disposition parameters (Table 2). Shared across groups:
    # "The disposition-related parameters (fd1, Kp1, and CL) were shared in
    # all groups" (Results 3.1). Tissue 2 is constrained to Tissue 1 by the
    # Table 2 footnotes (a: fd1 = fd2; b: Kp1 = Kp2) -- applied in model().
    # ------------------------------------------------------------------
    fd_rapidly_perfused <- 0.134        ; label("Fractional distribution parameter for Tissue 1 (unitless)")     # Table 2: 0.134 (SE 0.016, RSE 11.9%); footnote a: fd1 = fd2
    kp_rapidly_perfused <- 1.43         ; label("Tissue 1-to-plasma partition coefficient (unitless)")           # Table 2: 1.43 (SE 0.0973, RSE 6.81%); footnote b: Kp1 = Kp2
    lcl                 <- log(16.04)   ; label("Total (renal) clearance (mL/min/kg)")                           # Table 2: CL = 16.04 (SE 1.28, RSE 8.00%); Monolix CL = 16.04

    # ------------------------------------------------------------------
    # Apparent zero-order absorption rate constants (Table 2), in the units
    # the paper tabulates them (mg/min/kg). k01 and k02 "were shared in the
    # joint fitting" because control and malnourished did not differ
    # significantly, but Table 2 still reports the separate population
    # estimates used here (Results 3.1).
    # ------------------------------------------------------------------
    k01_control         <- 1.19         ; label("Zero-order absorption rate, phase 1, control (mg/min/kg)")      # Table 2: 1.19 (SE 0.0758, RSE 6.38%), t = 0-120 min
    k01_malnourished    <- 1.012        ; label("Zero-order absorption rate, phase 1, malnourished (mg/min/kg)") # Table 2: 1.012 (SE 0.0633, RSE 6.26%), significantly lower than k02
    k02_control         <- 1.86         ; label("Zero-order absorption rate, phase 2, control (mg/min/kg)")      # Table 2: 1.86 (SE 0.138, RSE 7.41%), t = 120-300 min
    k02_malnourished    <- 1.68         ; label("Zero-order absorption rate, phase 2, malnourished (mg/min/kg)") # Table 2: 1.68 (SE 0.134, RSE 7.99%)
    k03_control         <- fixed(0)     ; label("Zero-order absorption rate, phase 3, control (mg/min/kg)")      # Table 2 footnote c: "Assumed to be 0 according to GoF" -- no SE reported
    k03_malnourished    <- 1.16         ; label("Zero-order absorption rate, phase 3, malnourished (mg/min/kg)") # Table 2: 1.16 (SE 0.147, RSE 12.6%), t = 300-420 min

    # ------------------------------------------------------------------
    # Between-subject variability. "According to the evaluations, the random
    # effect of CL was included in the model for both drugs" (Methods 2.2.7)
    # -- CL is the ONLY parameter carrying an eta. Table 2 reports the
    # standard deviation of the random effect (the block is headed "Standard
    # deviation of the random effects"), so the variance is 0.124^2 =
    # 0.015376. Cross-check: sqrt(exp(0.124^2) - 1) = 12.45%, matching the
    # 12.5 CV% Table 2 prints alongside it.
    # ------------------------------------------------------------------
    etalcl ~ 0.015376                                                                                            # Table 2: omega CL = 0.124 (CV 12.5%), SE 0.0589, RSE 47.4%

    # ------------------------------------------------------------------
    # Residual error. "The proportional error model was used for all groups"
    # (Methods 2.2.7). Table 2 reports three separate proportional terms --
    # one per dataset. The two ORAL terms are selected by MAL_NOURISH in
    # model() and are the ones this model uses. The third, for the
    # literature IV reference arm, is b IV = 0.154 (SE 0.0328, RSE 21.3%);
    # it is recorded here rather than as a parameter because this model
    # simulates the oral arms, and an unused ini() entry does not parse.
    # ------------------------------------------------------------------
    propSd_control      <- 0.256        ; label("Proportional residual error, oral control (fraction)")          # Table 2 b C: 0.256 (SE 0.03204, RSE 12.5%)
    propSd_malnourished <- 0.234        ; label("Proportional residual error, oral malnourished (fraction)")     # Table 2 b M: 0.234 (SE 0.0272, RSE 11.6%)
  })

  model({
    # ================================================================
    # 1. Malnutrition switches. Disposition is deliberately NOT switched --
    #    fd1, Kp1 and CL are shared across groups (Results 3.1).
    # ================================================================
    bpr <- bpr_control + (bpr_malnourished - bpr_control) * MAL_NOURISH
    k01 <- k01_control + (k01_malnourished - k01_control) * MAL_NOURISH
    k02 <- k02_control + (k02_malnourished - k02_control) * MAL_NOURISH
    k03 <- k03_control + (k03_malnourished - k03_control) * MAL_NOURISH
    t_abs3 <- t_abs3_control + (t_abs3_malnourished - t_abs3_control) * MAL_NOURISH
    propSd <- propSd_control + (propSd_malnourished - propSd_control) * MAL_NOURISH

    # ================================================================
    # 2. Individual parameters. Tissue 2 is constrained to Tissue 1 by the
    #    Table 2 footnotes (fd1 = fd2, Kp1 = Kp2). fd1 + fd2 = 0.268 <= 1,
    #    satisfying the permeability-limited constraint of Equation 10.
    # ================================================================
    cl <- exp(lcl + etalcl)
    fd_slowly_perfused <- fd_rapidly_perfused
    kp_slowly_perfused <- kp_rapidly_perfused

    # ================================================================
    # 3. Sequential zero-order oral input (ug/min/kg).
    #
    #    The deposited code writes the input as a fraction of the dose per
    #    minute times the dose in ug/kg (Monolix: input_H = k01*250000;
    #    ADAPT: input_H = k01*Bolus(1,1)). The value Table 2 tabulates is
    #    that Monolix parameter scaled by 1000, so the mass input rate is
    #    k0_table * dose_po ug/min/kg. See the vignette Errata: taking the
    #    tabulated numbers directly as mg/min/kg overpredicts the paper's
    #    own fitted peak 4-fold and implies a bioavailability above 1.
    # ================================================================
    inp <- 0
    if (t <= t_abs1) {
      inp <- k01 * dose_po
    } else if (t <= t_abs2) {
      inp <- k02 * dose_po
    } else if (t <= t_abs3) {
      inp <- k03 * dose_po
    }

    # ================================================================
    # 4. Concentrations (ug/mL). Cc is the PLASMA concentration; the flow
    #    terms below carry BLOOD, hence the bpr factors.
    # ================================================================
    Cc                 <- a_blood            / v_blood
    c_rapidly_perfused <- a_rapidly_perfused / v_rapidly_perfused
    c_slowly_perfused  <- a_slowly_perfused  / v_slowly_perfused
    c_kidney           <- a_kidney           / v_kidney

    q_tissue <- q_co - q_kidney

    # ================================================================
    # 5. mPBPK mass balance -- Article Equations 1-5, transcribed from the
    #    deposited Monolix/ADAPT listings (which agree exactly). The states
    #    here are amounts (ug/kg) whose derivatives are the numerators of
    #    the deposited concentration ODEs; a_blood is scaled so that
    #    Cc = a_blood/v_blood, matching Monolix depot(target=B, p=1/58.82).
    #    Elimination is renal and is driven by the concentration ENTERING
    #    the kidney (CL * Cc), as written in the deposited code.
    # ================================================================
    d/dt(a_blood) <- (inp
      + q_tissue * fd_rapidly_perfused * c_rapidly_perfused * bpr / kp_rapidly_perfused
      + q_tissue * fd_slowly_perfused  * c_slowly_perfused  * bpr / kp_slowly_perfused
      + q_kidney * c_kidney * bpr / kp_kidney
      - (q_tissue * fd_rapidly_perfused * bpr
         + q_tissue * fd_slowly_perfused * bpr
         + q_kidney * bpr) * Cc) / bpr

    d/dt(a_rapidly_perfused) <-
      q_tissue * fd_rapidly_perfused * bpr * Cc -
      q_tissue * fd_rapidly_perfused * c_rapidly_perfused * bpr / kp_rapidly_perfused

    d/dt(a_slowly_perfused) <-
      q_tissue * fd_slowly_perfused * bpr * Cc -
      q_tissue * fd_slowly_perfused * c_slowly_perfused * bpr / kp_slowly_perfused

    d/dt(a_kidney) <-
      q_kidney * bpr * Cc -
      q_kidney * c_kidney * bpr / kp_kidney -
      cl * Cc

    # ================================================================
    # 6. Observation and error (Methods 2.2.7: proportional for ATN).
    # ================================================================
    Cc ~ prop(propSd)
  })
}
