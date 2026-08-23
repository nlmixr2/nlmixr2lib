Kir_2025_metoprolol_rat_pbpk <- function() {
  description <- "Preclinical (rat, Sprague-Dawley). mPBPK (minimal physiologically-based, Monolix 2023R1 population fit). Metoprolol (MET) oral absorption in non-malnourished (control) and malnourished rats. Sister model to modellib('Kir_2025_atenolol_rat_pbpk') from the same paper, and deliberately contrasted with it: metoprolol is lipophilic (logP 1.76, BCS Class I) and hepatically extracted, so the eliminating organ is the LIVER rather than the kidney, distribution is PERFUSION-limited (fd1 + fd2 = 1 exactly, versus 0.268 for atenolol), and Tissue 1 is kidney + lung rather than liver + lung. Four mass-balance states: a well-mixed blood pool, two lumped tissue groups (Tissue 1 rapidly perfused, high Kp; Tissue 2 slowly perfused remainder) and the liver. Crucially the oral zero-order input is delivered INTO THE LIVER, so hepatic first-pass extraction is structural rather than a fitted F -- roughly 50% of an oral dose is lost on first pass. Absorption is two SEQUENTIAL ZERO-ORDER processes with author-fixed windows (0-60 min, then 60-135 min in controls and 60-110 min in the malnourished group). Malnutrition raises both absorption rates and shortens the second window, taking apparent bioavailability from 0.42 to 0.84. The oral arms also carry an estimated FRACTION of the IV intrinsic clearance (fr = 0.336 control, 0.256 malnourished) rather than the full CLint, which is how the published fit reconciles the oral decline phases with the literature IV data. Blood-to-plasma ratio differs by group (1.508 vs 1.607). The model as shipped simulates the ORAL arms; the literature IV reference profile is reproduced by setting k01/k02 to 0, fr to 1 and bpr to 1.70, and giving a bolus into a_blood (see the vignette). IIV is on intrinsic clearance only. NOTE: the absorption-rate scale in Tables 2/3 required reconciliation against the deposited Monolix code -- see the vignette Errata."
  reference <- paste(
    "Kir F, Sahin S, Jusko WJ. (2025).",
    "Minimal Physiologically-Based Pharmacokinetic Modeling of Atenolol and",
    "Metoprolol Absorption in Malnourished Rats.",
    "Eur J Drug Metab Pharmacokinet 50:243-255.",
    "doi:10.1007/s13318-025-00943-6. PMCID PMC12081501.",
    "Population parameter estimates: Table 3 (Monolix). Naive-pooled (ADAPT 5)",
    "estimates for the same model: Supplementary Table S4.",
    "Rat physiology (tissue volumes and blood flows): Supplementary Table S2.",
    "Model equations: Article Equations 6-8 (MET), 9 (tissue-volume closure),",
    "11 (perfusion-limited fd constraint) and 12 (blood-to-plasma ratio).",
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
        "Malnutrition was confirmed biochemically in Table 1: body weight 233 g vs 302 g",
        "(p < 0.001), serum albumin 4.15 vs 4.68 g/dL and total cholesterol 49.5 vs",
        "66.5 mg/dL (both p < 0.001).",
        "In this model MAL_NOURISH switches (a) the zero-order absorption rate of both",
        "phases, (b) the end time of the second absorption phase, (c) the blood-to-plasma",
        "ratio, (d) the fraction of intrinsic clearance operating in the oral arm, and",
        "(e) the additive residual error magnitude.",
        "It deliberately does NOT change fd1, Kp1 or the underlying CLint, which were",
        "shared across groups in the published fit (Results 3.1)."
      ),
      source_name        = "Group"
    )
  )

  # Screened but NOT retained. Results 3: "Assessment of these measures as
  # covariates did not improve the population modeling."
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Table 1: 302 (SD 16.47) g control vs 233 (SD 14.25) g malnourished; screened, not retained."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/dL",
      type        = "continuous",
      notes       = paste(
        "Table 1: 4.68 (SD 0.13) g/dL control vs 4.15 (SD 0.15) g/dL malnourished.",
        "Not retained as a covariate, but albumin IS used indirectly: the fraction unbound",
        "was rescaled from the control value by the albumin ratio (fu 0.805 -> 0.908),",
        "which then enters the blood-to-plasma ratio via Equation 12."
      )
    ),
    TCHOL = list(
      description = "Total cholesterol",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Table 1: 66.5 (SD 10.01) mg/dL control vs 49.5 (SD 2.69) mg/dL malnourished; screened, not retained."
    )
  )

  compartmentData <- list(
    a_blood            = list(analyte = "metoprolol", units = "ug", specimen = "whole blood", verified = TRUE),
    a_rapidly_perfused = list(analyte = "metoprolol", units = "ug", specimen = "tissue", verified = TRUE),
    a_slowly_perfused  = list(analyte = "metoprolol", units = "ug", specimen = "tissue", verified = TRUE),
    a_liver            = list(analyte = "metoprolol", units = "ug", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species        = "rat (Sprague-Dawley)",
    n_subjects     = 8L,
    n_studies      = 2L,
    age_range      = NA_character_,
    weight_range   = "233 (SD 14.25) g malnourished to 302 (SD 16.47) g control",
    sex_female_pct = 0,
    disease_state  = "Experimental protein-calorie malnutrition (5% protein isocaloric diet, 17-20 days) vs control (20% protein isocaloric diet)",
    dose_range     = "312 mg/kg metoprolol base single oral dose (400 mg/kg metoprolol tartrate); literature IV reference 0.5, 1 and 2 mg/kg",
    regions        = "Turkey (Kobay Experimental Animals Laboratory, Ankara)",
    notes          = paste(
      "n = 4 control + 4 malnourished male Sprague-Dawley rats (Methods 2.2.1); the oral",
      "arms are the only rats actually dosed in this study. Metoprolol tartrate was",
      "dissolved in purified water at 400 mg/kg, equal to 312 mg/kg metoprolol base, and",
      "given by feeding tube. Blood was sampled from the tail vein at 0, 5, 15, 30, 45, 60,",
      "90, 105, 120, 150, 180, 240, 300, 360 and 420 min.",
      "The IV reference profile that informs disposition was NOT generated here: it was",
      "digitised from the literature (0.5, 1 and 2 mg/kg IV in male SD rats, reference 10",
      "of the paper, shown to be dose-linear), and the per-timepoint standard deviations",
      "were digitised and treated as 3 additional rats per dose (Methods 2.2.6).",
      "External validation used a further literature dataset dosed with R-metoprolol at",
      "1.5 mg/kg IV (reference 19); because that is a single enantiomer rather than the",
      "racemate, Kp_liver had to be refixed to 39.1 for that simulation only (Results 3.3).",
      "Rats were fasted overnight before dosing and fed ad libitum otherwise."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Rat physiological system parameters, per kg body weight.
    # Volumes: Supplementary Table S2 (lung 3.27, liver 33.81, kidney 6.16,
    # blood 58.82 mL/kg). For metoprolol the LIVER is the eliminating organ,
    # so Tissue 1 = kidney + lung = 6.16 + 3.27 = 9.43 (Methods 2.2.3: "The
    # kidney, and lung were defined as Tissue 1 with high Kp values").
    # Tissue 2 closes the body volume by Equation 9 (Body weight = Vb + V1 +
    # V2 + V3): 1000 - 58.82 - 9.43 - 33.81 = 897.94.
    # ------------------------------------------------------------------
    v_blood            <- fixed(58.82)  ; label("Blood volume (mL/kg)")                                          # Supp Table S2 (Shah & Betts 2012); Monolix/ADAPT Vb = 58.82
    v_rapidly_perfused <- fixed(9.43)   ; label("Tissue 1 volume, rapidly perfused = kidney + lung (mL/kg)")     # Supp Table S2 kidney 6.16 + lung 3.27; Monolix/ADAPT V1 = 9.43
    v_slowly_perfused  <- fixed(897.94) ; label("Tissue 2 volume, slowly perfused remainder (mL/kg)")            # Equation 9 closure; Monolix/ADAPT V2 = 897.94
    v_liver            <- fixed(33.81)  ; label("Liver volume (mL/kg)")                                          # Supp Table S2; Monolix/ADAPT Vliver = 33.81
    q_co               <- fixed(179.91) ; label("Cardiac output (mL/min/kg)")                                    # Supp Table S2; Monolix/ADAPT Qco = 179.91 (beta-blocker-reduced, Methods 2.2.3)
    q_liver            <- fixed(75.34)  ; label("Hepatic blood flow (mL/min/kg)")                                # Supp Table S2 75.348; Monolix/ADAPT Qliver = 75.34

    # Liver partition coefficient was PREDICTED, not fitted (Methods 2.2.3
    # quotes kidney 14 and lung 21.15 for Tissue 1, and Kp_hep = 62.25).
    kp_liver           <- fixed(62.25)  ; label("Liver-to-plasma partition coefficient (unitless)")              # Methods 2.2.3 "Kp_hep is the tissue-to-plasma partition coefficient for the liver (62.25)"; Monolix/ADAPT Kpliver = 62.25

    # Blood-to-plasma ratio, Equation 12: rho = (HCT - 1 + Rb)/(HCT * fu),
    # i.e. Rb = 1 - HCT + rho*HCT*fu, with rho fixed at 2.57 for MET.
    # HCT 47.5% control / 45.5% malnourished; fu 0.805 control / 0.908
    # malnourished (Methods 2.2.3). Both reproduce the paper's quoted Rb of
    # 1.508 and 1.607 exactly. The IV reference used Rb = 1.70 (the average
    # of R- and S-metoprolol) -- override bpr_control to simulate that arm.
    bpr_control        <- fixed(1.508)  ; label("Blood-to-plasma ratio, control (unitless)")                     # Methods 2.2.3 "the Rb values of MET were 1.508 and 1.607"; Monolix Rb_H = 1.508
    bpr_malnourished   <- fixed(1.607)  ; label("Blood-to-plasma ratio, malnourished (unitless)")                # Methods 2.2.3; Monolix (MET M) Rb = 1.607

    # Study design constant. The oral input rate scales with this value
    # because the deposited code writes input = k0 * Dose (see model()).
    dose_po            <- fixed(312.3274) ; label("Oral metoprolol base dose (mg/kg)")                           # Methods 2.2.1: 400 mg/kg tartrate = 312 mg/kg base; Monolix input_PO = k01_H*312327.4 (ug/kg)

    # ------------------------------------------------------------------
    # Absorption-window boundaries, all FIXED by the authors: "The t max was
    # fixed at 135 and 110 min for the control and malnourished groups,
    # respectively, and t 1 was fixed to 60 min for both groups by
    # optimization" (Results 3.1).
    # ------------------------------------------------------------------
    t_abs1             <- fixed(60)     ; label("End of absorption phase 1 (min)")                               # Results 3.1 "t 1 was fixed to 60 min for both groups"; Monolix/ADAPT tfirst = 60
    t_abs2_control     <- fixed(135)    ; label("End of absorption phase 2, control (min)")                      # Results 3.1; Monolix tmax_H = 135; Table 3 k02 C row "t = 60-135 min"
    t_abs2_malnourished <- fixed(110)   ; label("End of absorption phase 2, malnourished (min)")                 # Results 3.1; Monolix tmax_M = 110; Table 3 k02 M row "t = 60-110 min"

    # ------------------------------------------------------------------
    # Estimated disposition parameters (Table 3). Shared across groups:
    # "Parameters related to disposition (fd1, Kp2 and CL) were shared in all
    # groups" (Results 3.1). Distribution is PERFUSION-limited, so Table 3
    # footnote a constrains fd2 = 1 - fd1 (Equation 11) -- applied in model().
    # Kp1 = Kp2 as for atenolol.
    # ------------------------------------------------------------------
    fd_rapidly_perfused <- 0.464        ; label("Fractional distribution parameter for Tissue 1 (unitless)")     # Table 3: 0.464 (SE 0.104, RSE 22.4%); footnote a: fd2 = 1 - fd1, giving 0.536
    kp_rapidly_perfused <- 4.033        ; label("Tissue 1-to-plasma partition coefficient (unitless)")           # Table 3: 4.033 (SE 0.54, RSE 13.4%); Monolix Kp1 = Kp2 = 4.033
    lclint              <- log(148.63)  ; label("Hepatic intrinsic clearance (mL/min/kg)")                       # Table 3: CL int = 148 (SE 26.7, RSE 18.0%); deposited Monolix carries the unrounded CL = 148.63

    # ------------------------------------------------------------------
    # Fraction of the IV intrinsic clearance operating in each ORAL arm.
    # "The population models were optimized by using 0.34 and 0.27 fractions
    # of the CL obtained from the IV fitting for control and malnourished
    # groups according to goodness of fit" (Results 3.1); Table 3 gives the
    # unrounded estimates. Enters the liver ODE as CLint * fr (Monolix
    # writes this factor as `b`). Set to 1 to simulate the IV reference arm.
    # ------------------------------------------------------------------
    fr_control          <- 0.336        ; label("Fraction of intrinsic clearance, oral control (unitless)")      # Table 3 fr C: 0.336 (SE 0.0611, RSE 18.2%)
    fr_malnourished     <- 0.256        ; label("Fraction of intrinsic clearance, oral malnourished (unitless)") # Table 3 fr M: 0.256 (SE 0.0387, RSE 15.1%)

    # ------------------------------------------------------------------
    # Apparent zero-order absorption rate constants (Table 3), in the units
    # the paper tabulates them (mg/min/kg). Unlike atenolol these were NOT
    # shared: "The absorption rate constants were not shared in the joint
    # fitting and were evaluated as separate parameters to compare the two
    # groups" (Results 3.1).
    # ------------------------------------------------------------------
    k01_control         <- 4.5          ; label("Zero-order absorption rate, phase 1, control (mg/min/kg)")      # Table 3: 4.5 (SE 0.463, RSE 10.3%), t = 0-60 min; significantly lower than k01 M
    k01_malnourished    <- 7.56         ; label("Zero-order absorption rate, phase 1, malnourished (mg/min/kg)") # Table 3: 7.56 (SE 0.611, RSE 8.076%)
    k02_control         <- 5.57         ; label("Zero-order absorption rate, phase 2, control (mg/min/kg)")      # Table 3: 5.57 (SE 1.02, RSE 18.4%), t = 60-135 min
    k02_malnourished    <- 7.12         ; label("Zero-order absorption rate, phase 2, malnourished (mg/min/kg)") # Table 3: 7.12 (SE 1.12, RSE 15.8%), t = 60-110 min

    # ------------------------------------------------------------------
    # Between-subject variability on intrinsic clearance only (Methods
    # 2.2.7). Table 3 reports the standard deviation of the random effect,
    # so the variance is 0.258^2 = 0.066564. Cross-check:
    # sqrt(exp(0.258^2) - 1) = 26.2%, matching the CV% Table 3 prints.
    # ------------------------------------------------------------------
    etalclint ~ 0.066564                                                                                         # Table 3: omega CL = 0.258 (CV 26.2%), SE 0.137, RSE 53.2%

    # ------------------------------------------------------------------
    # Residual error. "The proportional error model was used for all groups,
    # while the constant error model was used to define MET oral profiles"
    # (Methods 2.2.7) -- so the two ORAL arms this model simulates use an
    # ADDITIVE error, and only the literature IV arm is proportional.
    # That IV term is b IV = 0.365 (SE 0.033, RSE 9.11%); it is recorded
    # here rather than as a parameter because this model simulates the oral
    # arms, and an unused ini() entry does not parse.
    # ------------------------------------------------------------------
    addSd_control       <- 1.76         ; label("Additive residual error, oral control (ug/mL)")                 # Table 3 a C: 1.76 (SE 0.194, RSE 11.04%)
    addSd_malnourished  <- 2.54         ; label("Additive residual error, oral malnourished (ug/mL)")            # Table 3 a M: 2.54 (SE 0.274, RSE 10.78%)
  })

  model({
    # ================================================================
    # 1. Malnutrition switches. fd1, Kp1 and the underlying CLint are
    #    deliberately NOT switched -- they were shared (Results 3.1).
    # ================================================================
    bpr    <- bpr_control + (bpr_malnourished - bpr_control) * MAL_NOURISH
    k01    <- k01_control + (k01_malnourished - k01_control) * MAL_NOURISH
    k02    <- k02_control + (k02_malnourished - k02_control) * MAL_NOURISH
    fr     <- fr_control  + (fr_malnourished  - fr_control)  * MAL_NOURISH
    t_abs2 <- t_abs2_control + (t_abs2_malnourished - t_abs2_control) * MAL_NOURISH
    addSd  <- addSd_control  + (addSd_malnourished  - addSd_control)  * MAL_NOURISH

    # ================================================================
    # 2. Individual parameters. Distribution is PERFUSION-limited, so
    #    Equation 11 / Table 3 footnote a set fd2 = 1 - fd1 (0.536); the two
    #    fractions sum to exactly 1, unlike the permeability-limited
    #    atenolol model where they sum to 0.268.
    # ================================================================
    clint <- exp(lclint + etalclint)
    fd_slowly_perfused <- 1 - fd_rapidly_perfused
    kp_slowly_perfused <- kp_rapidly_perfused

    # ================================================================
    # 3. Sequential zero-order oral input (ug/min/kg), delivered INTO THE
    #    LIVER so that hepatic first-pass extraction is structural.
    #
    #    The deposited code writes the input as a fraction of the dose per
    #    minute times the dose in ug/kg (Monolix: input_PO = k01_H*312327.4;
    #    ADAPT: k01*Bolus). The value Table 3 tabulates is that Monolix
    #    parameter scaled by 1000, so the mass input rate is
    #    k0_table * dose_po ug/min/kg. See the vignette Errata.
    # ================================================================
    inp <- 0
    if (t <= t_abs1) {
      inp <- k01 * dose_po
    } else if (t <= t_abs2) {
      inp <- k02 * dose_po
    }

    # ================================================================
    # 4. Concentrations (ug/mL). Cc is the PLASMA concentration; the flow
    #    terms below carry BLOOD, hence the bpr factors.
    # ================================================================
    Cc                 <- a_blood            / v_blood
    c_rapidly_perfused <- a_rapidly_perfused / v_rapidly_perfused
    c_slowly_perfused  <- a_slowly_perfused  / v_slowly_perfused
    c_liver            <- a_liver            / v_liver

    q_tissue <- q_co - q_liver

    # ================================================================
    # 5. mPBPK mass balance -- Article Equations 6-8, transcribed from the
    #    deposited Monolix/ADAPT listings (which agree exactly). The states
    #    are amounts (ug/kg) whose derivatives are the numerators of the
    #    deposited concentration ODEs; a_blood is scaled so that
    #    Cc = a_blood/v_blood, matching Monolix depot(target=B, p=1/58.82).
    #    Hepatic elimination is driven by the concentration LEAVING the
    #    liver (c_liver * bpr / kp_liver), i.e. a well-stirred liver, in
    #    contrast to the atenolol model where renal elimination is driven by
    #    the concentration entering the kidney.
    # ================================================================
    d/dt(a_blood) <- (
        q_tissue * fd_rapidly_perfused * c_rapidly_perfused * bpr / kp_rapidly_perfused
      + q_tissue * fd_slowly_perfused  * c_slowly_perfused  * bpr / kp_slowly_perfused
      + q_liver  * c_liver * bpr / kp_liver
      - (q_tissue * fd_rapidly_perfused * bpr
         + q_tissue * fd_slowly_perfused * bpr
         + q_liver * bpr) * Cc) / bpr

    d/dt(a_rapidly_perfused) <-
      q_tissue * fd_rapidly_perfused * bpr * Cc -
      q_tissue * fd_rapidly_perfused * c_rapidly_perfused * bpr / kp_rapidly_perfused

    d/dt(a_slowly_perfused) <-
      q_tissue * fd_slowly_perfused * bpr * Cc -
      q_tissue * fd_slowly_perfused * c_slowly_perfused * bpr / kp_slowly_perfused

    d/dt(a_liver) <- inp +
      q_liver * bpr * Cc -
      q_liver * c_liver * bpr / kp_liver -
      clint * fr * c_liver * bpr / kp_liver

    # ================================================================
    # 6. Observation and error (Methods 2.2.7: constant/additive for the
    #    MET oral profiles).
    # ================================================================
    Cc ~ add(addSd)
  })
}
