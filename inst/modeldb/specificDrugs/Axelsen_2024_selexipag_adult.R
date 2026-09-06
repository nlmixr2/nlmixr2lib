Axelsen_2024_selexipag_adult <- function() {
  description <- "Joint two-compartment parent + two-compartment metabolite population PK model for oral selexipag and its active metabolite JNJ-68006861 (ACT-333679) in adults with pulmonary arterial hypertension (Axelsen 2024, adult column of Table 1; the GRIPHON adult model of Krause 2017 re-estimated in NONMEM ahead of the pediatric analysis). First-order absorption into a two-compartment selexipag disposition; selexipag leaves the central compartment by a linear apparent clearance CL/F AND, in parallel, by the first-order metabolite-formation rate constant kmet, so total apparent selexipag clearance is CL/F + Vp/F * kmet. The metabolite has its own two-compartment disposition with first-order elimination km. The absorption lag time is logit-bounded on (0, 2) h with a fixed typical value of 0.668 h. Body weight (power on CL/F, Vp/F and Vm/F), total bilirubin (power on CL/F), male sex (exponential on km) and a four-level PAH-comedication categorical (naive / ERA only / PDE5 inhibitor only / ERA + PDE5; exponential on km) are the retained covariates."
  reference   <- "Axelsen LN, Kummel A, Perez Ruixo JJ, Russu A. Population pharmacokinetics of selexipag for dose selection and confirmation in pediatric patients with pulmonary arterial hypertension. CPT Pharmacometrics Syst Pharmacol. 2024;13(12):2185-2195. doi:10.1002/psp4.13231"
  vignette    <- "Axelsen_2024_selexipag"
  units       <- list(time = "h", dosing = "ug", concentration = "ng/mL")

  # Unit convention is stated verbatim at the head of the Supplementary
  # Table S3 NONMEM control stream: "Dose: ug / Concentration: ng/mL /
  # Time: hours". Table 1 labels CL "(1/h)" and the metabolite residual
  # error "(ug/mL)"; both are typographical slips -- see the vignette
  # Errata section.
  compartmentData <- list(
    depot           = list(analyte = "selexipag",     units = "ug", specimen = "administration site", verified = TRUE),
    central         = list(analyte = "selexipag",     units = "ug", specimen = "plasma",              verified = TRUE),
    peripheral1     = list(analyte = "selexipag",     units = "ug", specimen = "plasma",              verified = TRUE),
    central_act     = list(analyte = "JNJ-68006861",  units = "ug", specimen = "plasma",              verified = TRUE),
    peripheral1_act = list(analyte = "JNJ-68006861",  units = "ug", specimen = "plasma",              verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight at baseline.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effects on selexipag CL/F (exponent 0.546), selexipag Vp/F (exponent 1.04) and JNJ-68006861 Vm/F (exponent 0.803). Reference body weight 70 kg, hard-coded as log(WEIGHTBL/70) in the Supplementary Table S3 control stream MU_3 / MU_4 / MU_10 blocks and stated as 'centered around 70 kg' in Axelsen 2024 Table 1. Baseline (not time-varying) weight was used.",
      source_name        = "WEIGHTBL"
    ),
    TBILI = list(
      description        = "Total serum bilirubin at baseline.",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form negative effect on selexipag CL/F (exponent -0.442); reference 10 umol/L, hard-coded as log(BILIBL/10) in the Supplementary Table S3 control stream MU_3 block and stated as 'centered around 10 umol/L' in Axelsen 2024 Table 1. Hepatic-function marker; higher bilirubin lowers selexipag CL/F. Already in SI units in the source, so no mg/dL conversion is applied.",
      source_name        = "BILIBL"
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (female). The Supplementary Table S3 control stream builds SEXN_0 (' ; reference: 1') as the indicator for SEXN == 0, and Axelsen 2024 Table 1 labels the coefficient 'Gender male on km' -- so females are the reference and the coefficient is carried on the male indicator.",
      notes              = "Exponential effect on the JNJ-68006861 elimination rate constant km: km *= exp(0.148 * (1 - SEXF)). To preserve the paper's female-reference parameterisation while using the canonical SEXF column (1 = female), the model() block applies the coefficient to the male indicator (1 - SEXF). Males therefore have ~16% higher km (faster metabolite elimination) and correspondingly lower metabolite exposure than females.",
      source_name        = "SEXN (0 = male, 1 = female), entering the control stream as the derived indicator SEXN_0 = as.integer(SEXN == 0)"
    ),
    CONMED_ERA = list(
      description        = "Concomitant endothelin-receptor-antagonist (ERA) monotherapy indicator (1 = on an ERA but not on a PDE5 inhibitor, 0 = otherwise). One of three orthogonal mutually-exclusive indicators decomposing a four-level PAH-comedication categorical whose reference level is PAH-comedication-naive (all three indicators = 0).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the PAH-comedication-naive stratum; mutually exclusive with CONMED_PDE5I and CONMED_ERA_PDE5I).",
      notes              = "Exponential effect on the JNJ-68006861 elimination rate constant km: km *= exp(0.164 * CONMED_ERA), i.e. +17.8% relative to PAH-comedication-naive (Axelsen 2024 Table 1, beta_km(COPAH1_1) = 0.164, RSE 41.5%). The source data set already carries the decomposed binary columns COPAH1/COPAH2/COPAH3 (the control stream data path is named 'S12_data_binCOPAH'), so no re-coding of a multi-level factor is involved.",
      source_name        = "COPAH1"
    ),
    CONMED_PDE5I = list(
      description        = "Concomitant phosphodiesterase type 5 inhibitor (PDE5I) monotherapy indicator (1 = on a PDE5I but not on an ERA, 0 = otherwise). One of three orthogonal mutually-exclusive indicators decomposing a four-level PAH-comedication categorical whose reference level is PAH-comedication-naive.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the PAH-comedication-naive stratum; mutually exclusive with CONMED_ERA and CONMED_ERA_PDE5I).",
      notes              = "Exponential effect on the JNJ-68006861 elimination rate constant km: km *= exp(0.0633 * CONMED_PDE5I), i.e. +6.5% relative to PAH-comedication-naive (Axelsen 2024 Table 1, beta_km(COPAH2_1) = 0.0633). The coefficient is imprecise (RSE 96.1%) but is retained because the four-level PAH-comedication categorical was kept intact end-to-end.",
      source_name        = "COPAH2"
    ),
    CONMED_ERA_PDE5I = list(
      description        = "Concomitant ERA + PDE5-inhibitor combination indicator (1 = on both an ERA and a PDE5 inhibitor, 0 = otherwise). One of three orthogonal mutually-exclusive indicators decomposing a four-level PAH-comedication categorical whose reference level is PAH-comedication-naive.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the PAH-comedication-naive stratum; mutually exclusive with CONMED_ERA and CONMED_PDE5I).",
      notes              = "Exponential effect on the JNJ-68006861 elimination rate constant km: km *= exp(0.358 * CONMED_ERA_PDE5I), i.e. +43.0% relative to PAH-comedication-naive (Axelsen 2024 Table 1, beta_km(COPAH3_1) = 0.358, RSE 17.2%). The combined stratum carries its own coefficient rather than the sum of the ERA-only and PDE5I-only coefficients, because the control stream adds one THETA per non-reference level of the four-level categorical (Supplementary Table S3, MU7WRAP_1 / MU7WRAP_2).",
      source_name        = "COPAH3"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 510,
    n_studies      = 1,
    age_range      = "Adult PAH patients enrolled in GRIPHON (NCT01106014). Per-subject ages are not tabulated in Axelsen 2024; the adult cohort is characterised in the GRIPHON popPK source (Krause 2017) and the GRIPHON main paper (Sitbon 2015 NEJM).",
    weight_range   = "40-148 kg -- the body-weight span over which Axelsen 2024 Results ('Selexipag pediatric dose selection based on adult data') states the continuous body-weight-exposure relationship was identified. Reference weight 70 kg.",
    sex_female_pct = "Not reported in Axelsen 2024; sex was retained as a covariate on the metabolite elimination rate constant km with females as the reference category.",
    disease_state  = "Pulmonary arterial hypertension (WHO Group I), on stable background PAH comedication (endothelin-receptor antagonist and/or phosphodiesterase type 5 inhibitor) or PAH-comedication-naive.",
    dose_range     = "Selexipag 200 ug twice daily, up-titrated weekly in 200 ug increments to the individual maximum tolerated dose, maximum 1600 ug twice daily.",
    regions        = "Multinational phase III GRIPHON study (NCT01106014).",
    co_medication  = "PAH-specific comedication as a four-level categorical (naive / ERA only / PDE5 inhibitor only / ERA + PDE5), reference = naive. Approximately 20% of the adult GRIPHON participants received no other PAH comedication (Axelsen 2024 Results, 'PK data collected in pediatric participants').",
    notes          = "The N = 510 count is the number of adult participants contributing model-based AUCtau,ss,combined values in Axelsen 2024 Table 3, which is the adult reference distribution the pediatric exposures were compared against. The parameter values here are the adult column of Axelsen 2024 Table 1: the Krause 2017 GRIPHON adult model, originally fitted in Monolix, re-estimated in NONMEM because NONMEM was planned for the pediatric analysis. They are reproduced at full precision in the Supplementary Table S3 control stream $THETAP / $OMEGAP / $SIGMAP blocks, where they serve as the Bayesian prior for the pediatric update (see Axelsen_2024_selexipag_pediatric)."
  )

  ini({
    # ---------------- Parent (selexipag) structural parameters ----------------
    # All volumes and clearances are apparent (implicit "/F"); the control
    # stream sets F1 = 1 (Supplementary Table S3 $PK, "Dosing compartments
    # info"). Reference covariate values are 70 kg body weight and
    # 10 umol/L total bilirubin, hard-coded in the control stream as
    # log(WEIGHTBL/70) and log(BILIBL/10).
    lka  <- log(0.688);  label("First-order absorption rate constant ka (1/h)")                                     # Table 1 adult column: ka = 0.688 1/h, RSE 4.27% (control stream $THETA 1 = -0.373883 on the log scale)
    lcl  <- log(18.0);   label("Apparent selexipag clearance CL/F at 70 kg / 10 umol/L bilirubin (L/h)")            # Table 1 adult column: CL = 18 L/h, RSE 7.56% (control stream $THETA 3 = 2.89048). Excludes the metabolite-formation route; total apparent clearance is CL + Vp*kmet
    lvc  <- log(12.2);   label("Apparent selexipag central volume Vp/F at 70 kg (L)")                               # Table 1 adult column: Vp = 12.2 L, RSE 8.57% (control stream $THETA 4 = 2.50310)
    lk12 <- log(0.101);  label("Selexipag central-to-peripheral rate constant kx12 (1/h)")                          # Table 1 adult column: kx12 = 0.101 1/h, RSE 10.1% (control stream $THETA 5 = -2.29026)
    lk21 <- log(0.0521); label("Selexipag peripheral-to-central rate constant kx21 (1/h)")                          # Table 1 adult column: kx21 = 0.0521 1/h, RSE 21.6% (control stream $THETA 6 = -2.95489)

    # Absorption lag time, carried on the logit scale. The source
    # parameterises the FIXED quantity Tlag1Half = expit(THETA(11)), which is
    # bounded on (0, 1), and sets the NONMEM lag as ALAG1 = 2 * Tlag1Half --
    # i.e. the lag time itself is bounded on (0, 2) h. THETA(11) = -0.69 FIX
    # back-transforms to Tlag1Half = 0.334 h (Table 1) and hence a typical lag
    # of 0.668 h, matching the value Krause 2017 fixed from the healthy-subject
    # popPK model. The logit-normal IIV is what makes the bound load-bearing.
    logittlag <- fixed(-0.69); label("Logit-scale absorption lag time; tlag = 2 * expit(logittlag), bounded on (0, 2) h")  # Supplementary Table S3 control stream $THETA 11: "-0.69 FIX ; 11 log(Tlag1Half/(1-Tlag1Half)) (0.334)"; Table 1 adult column: Tlag1Half = 0.334 h (FIX)

    # ---------------- Metabolite (JNJ-68006861 / ACT-333679) ------------------
    # kmet is the first-order rate constant carrying selexipag out of the
    # parent central compartment INTO the metabolite central compartment
    # (control stream ADVAN5 transfer k2T4 = kmet), in addition to the CL/Vp
    # elimination (k2T0 = CL/Vp). The paper's printed AUCtau,ss,combined
    # formula uses the same convention: AUCtau,ss,selexipag = Dose / (Vp*kmet + CL/F).
    lkmet    <- log(0.887);  label("Selexipag-to-JNJ-68006861 metabolite-formation rate constant kmet (1/h)")       # Table 1 adult column: kmet = 0.887 1/h, RSE 4.06% (control stream $THETA 2 = -0.120236)
    lvc_act  <- log(5.88);   label("Apparent JNJ-68006861 central volume Vm/F at 70 kg (L)")                        # Table 1 adult column: Vm = 5.88 L, RSE 4.37% (control stream $THETA 10 = 1.77121)
    lkm_act  <- log(0.468);  label("JNJ-68006861 elimination rate constant km at the female / PAH-naive reference (1/h)")  # Table 1 adult column: km = 0.468 1/h, RSE 7.52% (control stream $THETA 7 = -0.759585)
    lk34_act <- log(0.898);  label("JNJ-68006861 central-to-peripheral rate constant kx34 (1/h)")                   # Table 1 adult column: kx34 = 0.898 1/h, RSE 10.8% (control stream $THETA 8 = -0.107333)
    lk43_act <- log(0.183);  label("JNJ-68006861 peripheral-to-central rate constant kx43 (1/h)")                   # Table 1 adult column: kx43 = 0.183 1/h, RSE 15% (control stream $THETA 9 = -1.69993)

    # ---------------- Covariate effects ---------------------------------------
    # Every covariate effect is EXPONENTIAL: the control stream builds each
    # MU_i as a sum on the log scale and then exponentiates
    # (e.g. MU_3 = THETA(3) + THETA(12)*log(WEIGHTBL/70) + THETA(13)*log(BILIBL/10);
    # CL = EXP(T_CL)). For the continuous covariates that is the usual power
    # model; for the categorical covariates it is exp(beta * indicator), NOT
    # the (1 + beta * indicator) form used by the sibling Krause_2017_selexipag
    # extraction.
    e_wt_cl        <-  0.546;  label("Power exponent of body weight on selexipag CL/F (unitless)")                  # Table 1 adult column: beta_CL(WEIGHTBL) = 0.546, RSE 36.4% (control stream $THETA 12 = 0.545841)
    e_tbili_cl     <- -0.442;  label("Power exponent of total bilirubin on selexipag CL/F (unitless)")              # Table 1 adult column: beta_CL(BILIBL) = -0.442, RSE 21.6% (control stream $THETA 13 = -0.441952)
    e_wt_vc        <-  1.04;   label("Power exponent of body weight on selexipag Vp/F (unitless)")                  # Table 1 adult column: beta_Vp(WEIGHTBL) = 1.04, RSE 17.7% (control stream $THETA 15 = 1.03797)
    e_wt_vc_act    <-  0.803;  label("Power exponent of body weight on JNJ-68006861 Vm/F (unitless)")               # Table 1 adult column: beta_Vm(WEIGHTBL) = 0.803, RSE 12.7% (control stream $THETA 14 = 0.80287)
    e_sexf_km_act  <-  0.148;  label("Log-scale effect of male sex on JNJ-68006861 km (unitless)")                  # Table 1 adult column: beta_km(SEXN_0) = 0.148, RSE 31.6%; carried on the male indicator (1 - SEXF) because females are the reference
    e_era_km_act   <-  0.164;  label("Log-scale effect of ERA-only PAH comedication on JNJ-68006861 km")            # Table 1 adult column: beta_km(COPAH1_1) = 0.164, RSE 41.5%
    e_pde5_km_act  <-  0.0633; label("Log-scale effect of PDE5I-only PAH comedication on JNJ-68006861 km")          # Table 1 adult column: beta_km(COPAH2_1) = 0.0633, RSE 96.1% (imprecise but retained to keep the four-level categorical intact)
    e_combo_km_act <-  0.358;  label("Log-scale effect of ERA + PDE5I PAH comedication on JNJ-68006861 km")         # Table 1 adult column: beta_km(COPAH3_1) = 0.358, RSE 17.2%

    # ---------------- IIV (variances on the log / logit scale) ----------------
    # Table 1's note states "omega values reported as standard deviation", and
    # the control stream $OMEGA / $OMEGAP blocks carry the STANDARD keyword,
    # confirming the tabulated omegas are SDs. ini() takes variances, so each
    # entry below is omega^2 with the source omega quoted in the comment.
    etalka        ~ 0.155236    # 0.394^2;  Table 1 adult column: omega(ka) = 0.394, RSE 8.36%
    etalkmet      ~ 0.00695556  # 0.0834^2; Table 1 adult column: omega(kmet) = 0.0834, RSE 36.8% (86% shrinkage)
    etalcl        ~ 0.641601    # 0.801^2;  Table 1 adult column: omega(CL) = 0.801, RSE 6.67%
    etalvc        ~ 0.126736    # 0.356^2;  Table 1 adult column: omega(Vp) = 0.356, RSE 17.1%
    etalk12       ~ 0.185761    # 0.431^2;  Table 1 adult column: omega(kx12) = 0.431, RSE 23.7%
    etalk21       ~ 1.21        # 1.1^2;    Table 1 adult column: omega(kx21) = 1.1, RSE 12%
    etalkm_act    ~ 0.073984    # 0.272^2;  Table 1 adult column: omega(km) = 0.272, RSE 8.92%
    etalk34_act   ~ 0.096721    # 0.311^2;  Table 1 adult column: omega(kx34) = 0.311, RSE 31.5%
    etalk43_act   ~ 0.651249    # 0.807^2;  Table 1 adult column: omega(kx43) = 0.807, RSE 13.6%
    etalvc_act    ~ 0.00848241  # 0.0921^2; Table 1 adult column: omega(Vm) = 0.0921, RSE 41.2% (85% shrinkage)
    etalogittlag  ~ 3.0976      # 1.76^2;   Table 1 adult column: omega(Tlag1Half) = 1.76, RSE 12%, distribution LogitNormal -- this SD is on the logit scale

    # ---------------- Residual error (proportional, by output) ----------------
    # The control stream $ERROR uses Y = IPRED + EPS*IPRED for both outputs
    # (pure proportional), and $SIGMA carries the variances 0.567954 =
    # 0.753627^2 and 0.239527 = 0.489415^2 -- confirming that the error_PROP
    # values printed in Table 1 are proportional standard deviations.
    propSd     <- 0.754; label("Selexipag proportional residual SD (fraction)")                                     # Table 1 adult column: error_PROP1 = 0.754, RSE 3.21%
    propSd_act <- 0.489; label("JNJ-68006861 proportional residual SD (fraction)")                                  # Table 1 adult column: error_PROP2 = 0.489, RSE 1.9%
  })

  model({
    # Reference covariate values, hard-coded in the Supplementary Table S3
    # control stream as log(WEIGHTBL/70) and log(BILIBL/10).
    ref_wt    <- 70    # kg
    ref_tbili <- 10    # umol/L

    # ----- Individual parameters: selexipag (parent) --------------------------
    ka   <- exp(lka + etalka)
    kmet <- exp(lkmet + etalkmet)
    cl   <- exp(lcl + etalcl) * (WT / ref_wt)^e_wt_cl * (TBILI / ref_tbili)^e_tbili_cl
    vc   <- exp(lvc + etalvc) * (WT / ref_wt)^e_wt_vc
    k12  <- exp(lk12 + etalk12)
    k21  <- exp(lk21 + etalk21)

    # Lag time bounded on (0, 2) h: ALAG1 = 2 * Tlag1Half with
    # Tlag1Half = expit(THETA(11) + ETA(11)) (control stream $PK, "Parameter
    # transformations" and "Dosing compartments info").
    tlag <- 2 * expit(logittlag + etalogittlag)

    # ----- Individual parameters: JNJ-68006861 (metabolite) -------------------
    # The sex and PAH-comedication coefficients sit inside exp() because the
    # control stream adds them to MU_7 on the log scale before exponentiating.
    # The male indicator is (1 - SEXF): the control stream's SEXN_0 flags
    # SEXN == 0 with females (SEXN == 1) as the reference.
    vc_act  <- exp(lvc_act + etalvc_act) * (WT / ref_wt)^e_wt_vc_act
    km_act  <- exp(lkm_act + etalkm_act +
                     e_sexf_km_act  * (1 - SEXF) +
                     e_era_km_act   * CONMED_ERA +
                     e_pde5_km_act  * CONMED_PDE5I +
                     e_combo_km_act * CONMED_ERA_PDE5I)
    k34_act <- exp(lk34_act + etalk34_act)
    k43_act <- exp(lk43_act + etalk43_act)

    # Selexipag elimination via routes other than metabolite formation.
    kel <- cl / vc

    # ----- ODE system ---------------------------------------------------------
    # Mirrors the ADVAN5 transfer-rate block of the Supplementary Table S3
    # control stream, whose non-zero rates simplify to:
    #   k1T2 = ka          depot -> selexipag central
    #   k2T0 = CL/Vp       selexipag central -> eliminated
    #   k2T3 = kx12        selexipag central -> selexipag peripheral
    #   k3T2 = kx21        selexipag peripheral -> selexipag central
    #   k2T4 = kmet        selexipag central -> metabolite central
    #   k4T0 = km          metabolite central -> eliminated
    #   k4T5 = kx34        metabolite central -> metabolite peripheral
    #   k5T4 = kx43        metabolite peripheral -> metabolite central
    # Selexipag therefore leaves the central compartment at kel + kmet, and
    # the metabolite-formation flux kmet * central is mass-conserving.
    d/dt(depot)           <- -ka * depot
    d/dt(central)         <-  ka * depot - kel * central - kmet * central -
                              k12 * central + k21 * peripheral1
    d/dt(peripheral1)     <-  k12 * central - k21 * peripheral1
    d/dt(central_act)     <-  kmet * central - km_act * central_act -
                              k34_act * central_act + k43_act * peripheral1_act
    d/dt(peripheral1_act) <-  k34_act * central_act - k43_act * peripheral1_act

    alag(depot) <- tlag

    # ----- Observations -------------------------------------------------------
    Cc     <- central     / vc
    Cc_act <- central_act / vc_act

    Cc     ~ prop(propSd)
    Cc_act ~ prop(propSd_act)
  })
}
