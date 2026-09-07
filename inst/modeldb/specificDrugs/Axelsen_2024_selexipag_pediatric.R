Axelsen_2024_selexipag_pediatric <- function() {
  description <- "Joint two-compartment parent + two-compartment metabolite population PK model for oral selexipag and its active metabolite JNJ-68006861 (ACT-333679) in pediatric patients aged 2 to <18 years with pulmonary arterial hypertension (Axelsen 2024, pediatric column of Table 1; NCT03492177). Structurally identical to the companion adult model, updated by NONMEM BAYES estimation using the adult parameter estimates as priors. First-order absorption into a two-compartment selexipag disposition; selexipag leaves the central compartment by a linear apparent clearance CL/F AND, in parallel, by the first-order metabolite-formation rate constant kmet, so total apparent selexipag clearance is CL/F + Vp/F * kmet. The metabolite has its own two-compartment disposition with first-order elimination km. The absorption lag time is logit-bounded on (0, 2) h with a fixed typical value of 0.668 h. Body weight (power on CL/F, Vp/F and Vm/F), total bilirubin (power on CL/F), male sex (exponential on km) and a four-level PAH-comedication categorical (naive / ERA only / PDE5 inhibitor only / ERA + PDE5; exponential on km) are the retained covariates. The body-weight exponent on CL/F updated from 0.546 to 0.828, close to the standard allometric 0.75."
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
      notes              = "Power-form effects on selexipag CL/F (exponent 0.828), selexipag Vp/F (exponent 0.805) and JNJ-68006861 Vm/F (exponent 0.542). Reference body weight 70 kg, hard-coded as log(WEIGHTBL/70) in the Supplementary Table S3 control stream MU_3 / MU_4 / MU_10 blocks and stated as 'centered around 70 kg' in Axelsen 2024 Table 1. Baseline (not time-varying) weight was used deliberately: Axelsen 2024 Methods states weight change over the 12-week PK sampling window was expected to be limited. Observed pediatric range 9.9-93.5 kg. Note the 70 kg reference lies well above the pediatric cohort, so the weight term is an extrapolation downwards for every subject.",
      source_name        = "WEIGHTBL"
    ),
    TBILI = list(
      description        = "Total serum bilirubin at baseline.",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form negative effect on selexipag CL/F (exponent -0.363); reference 10 umol/L, hard-coded as log(BILIBL/10) in the Supplementary Table S3 control stream MU_3 block. Observed pediatric means (SD) [range] by age cohort were 15.1 (16.8) [3-84], 10.4 (6.77) [3-32] and 7.9 (5.57) [3-23] umol/L for the 12-17, 6-11 and 2-5 year cohorts respectively (Table S3). Already in SI units in the source, so no mg/dL conversion is applied.",
      source_name        = "BILIBL"
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (female). The Supplementary Table S3 control stream builds SEXN_0 (' ; reference: 1') as the indicator for SEXN == 0, and Axelsen 2024 Table 1 labels the coefficient 'Gender male on km' -- so females are the reference and the coefficient is carried on the male indicator.",
      notes              = "Exponential effect on the JNJ-68006861 elimination rate constant km: km *= exp(0.145 * (1 - SEXF)). To preserve the paper's female-reference parameterisation while using the canonical SEXF column (1 = female), the model() block applies the coefficient to the male indicator (1 - SEXF). The pediatric cohort was 57.1% female overall (Table S1).",
      source_name        = "SEXN (0 = male, 1 = female), entering the control stream as the derived indicator SEXN_0 = as.integer(SEXN == 0)"
    ),
    CONMED_ERA = list(
      description        = "Concomitant endothelin-receptor-antagonist (ERA) monotherapy indicator (1 = on an ERA but not on a PDE5 inhibitor, 0 = otherwise). One of three orthogonal mutually-exclusive indicators decomposing a four-level PAH-comedication categorical whose reference level is PAH-comedication-naive (all three indicators = 0).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the PAH-comedication-naive stratum; mutually exclusive with CONMED_PDE5I and CONMED_ERA_PDE5I).",
      notes              = "Exponential effect on the JNJ-68006861 elimination rate constant km: km *= exp(0.186 * CONMED_ERA), i.e. +20.4% relative to PAH-comedication-naive (Axelsen 2024 Table 1, beta_km(COPAH1_1) = 0.186, RSE 34.8%). ERA-only prevalence in the pediatric cohort was 14.3% / 9.52% / 10% across the 12-17, 6-11 and 2-5 year cohorts (Table S3).",
      source_name        = "COPAH1"
    ),
    CONMED_PDE5I = list(
      description        = "Concomitant phosphodiesterase type 5 inhibitor (PDE5I) monotherapy indicator (1 = on a PDE5I but not on an ERA, 0 = otherwise). One of three orthogonal mutually-exclusive indicators decomposing a four-level PAH-comedication categorical whose reference level is PAH-comedication-naive.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the PAH-comedication-naive stratum; mutually exclusive with CONMED_ERA and CONMED_ERA_PDE5I).",
      notes              = "Exponential effect on the JNJ-68006861 elimination rate constant km: km *= exp(0.0495 * CONMED_PDE5I), i.e. +5.1% relative to PAH-comedication-naive (Axelsen 2024 Table 1, beta_km(COPAH2_1) = 0.0495). The coefficient is imprecise (RSE 116%) but is retained because the four-level PAH-comedication categorical was kept intact end-to-end. PDE5I-only prevalence was 23.8% / 57.1% / 15% across the three age cohorts (Table S3).",
      source_name        = "COPAH2"
    ),
    CONMED_ERA_PDE5I = list(
      description        = "Concomitant ERA + PDE5-inhibitor combination indicator (1 = on both an ERA and a PDE5 inhibitor, 0 = otherwise). One of three orthogonal mutually-exclusive indicators decomposing a four-level PAH-comedication categorical whose reference level is PAH-comedication-naive.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the PAH-comedication-naive stratum; mutually exclusive with CONMED_ERA and CONMED_PDE5I).",
      notes              = "Exponential effect on the JNJ-68006861 elimination rate constant km: km *= exp(0.368 * CONMED_ERA_PDE5I), i.e. +44.5% relative to PAH-comedication-naive (Axelsen 2024 Table 1, beta_km(COPAH3_1) = 0.368, RSE 15.3%). This is by far the most common stratum in the pediatric cohort (52.4% / 28.6% / 75% across the three age cohorts, Table S3); only 3 of 59 pediatric participants with full PK profiles received no PAH comedication at all.",
      source_name        = "COPAH3"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 59,
    n_studies      = 1,
    age_range      = "2-17 years (inclusion >=2 to <18 years), enrolled in three age cohorts: >=12 to <18 years (N = 21), >=6 to <12 years (N = 21), >=2 to <6 years (N = 17 with full PK profiles, +3 with trough concentrations only). Cohort mean (SD) ages 14.1 (1.71), 8.57 (1.43) and 3.8 (1.28) years (Table S3).",
    weight_range   = "9.9-93.5 kg. Cohort mean (SD) [range] baseline weights 51.1 (14.3) [30-93.5], 25.0 (6.49) [16.8-36.4] and 16.4 (6.95) [9.9-41.5] kg for the 12-17, 6-11 and 2-5 year cohorts (Table S3).",
    sex_female_pct = "57.1% female overall (36/63 enrolled; Table S1). By age cohort: 66.7%, 52.4% and 50.0% female for the 12-17, 6-11 and 2-5 year cohorts (Table S3).",
    race_ethnicity = "White 61.9%, Asian 25.4%, Other 3.2%, Unknown 9.5%; Hispanic or Latino 1.6% (Table S1, all 63 enrolled participants). Race was not tested as a covariate.",
    disease_state  = "Pediatric pulmonary arterial hypertension: idiopathic PAH, heritable PAH, PAH associated with congenital heart disease, and PAH associated with HIV, connective tissue disease, or drug/toxin exposure; WHO Functional Class II or III.",
    dose_range     = "Body-weight-banded starting doses twice daily -- 100 ug for >=9 to <25 kg, 150 ug for >=25 to <50 kg, 200 ug for >=50 kg -- up-titrated weekly in increments equal to the starting dose over a 12-week titration period to the individual maximum tolerated dose, capped at 8-fold the starting dose (800 / 1200 / 1600 ug twice daily). Observed doses ranged 50-1600 ug twice daily; one participant assigned 100 ug initially received 50 ug in error.",
    regions        = "Prospective, multicenter, open-label, single-arm phase II study NCT03492177.",
    co_medication  = "PAH-specific comedication as a four-level categorical (naive / ERA only / PDE5 inhibitor only / ERA + PDE5), reference = naive. Only 3 of 59 participants (5.1%) with full PK profiles received no other PAH comedication, versus roughly 20% in the adult GRIPHON study.",
    n_observations = "1167 selexipag and JNJ-68006861 observations analysed, of 1198 in the dataset; 31 (2.6%) were excluded, 27 of them (2.3%) because selexipag concentrations were below the 0.01 ng/mL limit of quantification.",
    notes          = "Parameters were obtained by NONMEM BAYES estimation (500 burn-in, 2000 iterations) using the adult model estimates as priors -- normal on the natural logarithm of the fixed effects and inverse-Wishart for the variances (Supplementary Table S3 $PRIOR NWPRI / $THETAP / $OMEGAP / $SIGMAP blocks). No additional covariate search was performed; the covariate set was inherited from the adult model. Serial PK profiles (predose and 1, 2, 4, 6, 8 and 12 h postdose) were taken at steady state either at Week 1 (first five participants per cohort) or Week 12, plus three steady-state troughs per participant during titration. Selexipag was given after a light meal on PK days. Shrinkage is high for omega(kmet) (76.2%) and omega(Vm) (79.6%), so those two IIV terms are weakly informed by the pediatric data and are effectively carried from the adult prior."
  )

  ini({
    # ---------------- Parent (selexipag) structural parameters ----------------
    # All volumes and clearances are apparent (implicit "/F"); the control
    # stream sets F1 = 1 (Supplementary Table S3 $PK, "Dosing compartments
    # info"). Reference covariate values are 70 kg body weight and
    # 10 umol/L total bilirubin, hard-coded in the control stream as
    # log(WEIGHTBL/70) and log(BILIBL/10).
    lka  <- log(0.663);  label("First-order absorption rate constant ka (1/h)")                                     # Table 1 pediatric column: ka = 0.663 1/h, RSE 3.67%
    lcl  <- log(18.9);   label("Apparent selexipag clearance CL/F at 70 kg / 10 umol/L bilirubin (L/h)")            # Table 1 pediatric column: CL = 18.9 L/h, RSE 6.71%. Excludes the metabolite-formation route; total apparent clearance is CL + Vp*kmet
    lvc  <- log(12.4);   label("Apparent selexipag central volume Vp/F at 70 kg (L)")                               # Table 1 pediatric column: Vp = 12.4 L, RSE 7.32%
    lk12 <- log(0.115);  label("Selexipag central-to-peripheral rate constant kx12 (1/h)")                          # Table 1 pediatric column: kx12 = 0.115 1/h, RSE 8.48%
    lk21 <- log(0.0563); label("Selexipag peripheral-to-central rate constant kx21 (1/h)")                          # Table 1 pediatric column: kx21 = 0.0563 1/h, RSE 16%

    # Absorption lag time, carried on the logit scale. The source
    # parameterises the FIXED quantity Tlag1Half = expit(THETA(11)), which is
    # bounded on (0, 1), and sets the NONMEM lag as ALAG1 = 2 * Tlag1Half --
    # i.e. the lag time itself is bounded on (0, 2) h. THETA(11) = -0.69 FIX
    # back-transforms to Tlag1Half = 0.334 h (Table 1) and hence a typical lag
    # of 0.668 h. The value was not updated by the pediatric fit; it stays
    # fixed at the adult value.
    logittlag <- fixed(-0.69); label("Logit-scale absorption lag time; tlag = 2 * expit(logittlag), bounded on (0, 2) h")  # Supplementary Table S3 control stream $THETA 11: "-0.69 FIX ; 11 log(Tlag1Half/(1-Tlag1Half)) (0.334)"; Table 1 pediatric column: Tlag1Half = 0.334 h (FIX)

    # ---------------- Metabolite (JNJ-68006861 / ACT-333679) ------------------
    # kmet is the first-order rate constant carrying selexipag out of the
    # parent central compartment INTO the metabolite central compartment
    # (control stream ADVAN5 transfer k2T4 = kmet), in addition to the CL/Vp
    # elimination (k2T0 = CL/Vp). The paper's printed AUCtau,ss,combined
    # formula uses the same convention: AUCtau,ss,selexipag = Dose / (Vp*kmet + CL/F).
    lkmet    <- log(0.868);  label("Selexipag-to-JNJ-68006861 metabolite-formation rate constant kmet (1/h)")       # Table 1 pediatric column: kmet = 0.868 1/h, RSE 3.77%
    lvc_act  <- log(6.05);   label("Apparent JNJ-68006861 central volume Vm/F at 70 kg (L)")                        # Table 1 pediatric column: Vm = 6.05 L, RSE 4.04%
    lkm_act  <- log(0.445);  label("JNJ-68006861 elimination rate constant km at the female / PAH-naive reference (1/h)")  # Table 1 pediatric column: km = 0.445 1/h, RSE 6.52%
    lk34_act <- log(0.83);   label("JNJ-68006861 central-to-peripheral rate constant kx34 (1/h)")                   # Table 1 pediatric column: kx34 = 0.83 1/h, RSE 8.35%
    lk43_act <- log(0.168);  label("JNJ-68006861 peripheral-to-central rate constant kx43 (1/h)")                   # Table 1 pediatric column: kx43 = 0.168 1/h, RSE 12.5%

    # ---------------- Covariate effects ---------------------------------------
    # Every covariate effect is EXPONENTIAL: the control stream builds each
    # MU_i as a sum on the log scale and then exponentiates
    # (e.g. MU_3 = THETA(3) + THETA(12)*log(WEIGHTBL/70) + THETA(13)*log(BILIBL/10);
    # CL = EXP(T_CL)). For the continuous covariates that is the usual power
    # model; for the categorical covariates it is exp(beta * indicator), NOT
    # the (1 + beta * indicator) form used by the sibling Krause_2017_selexipag
    # extraction.
    e_wt_cl        <-  0.828;  label("Power exponent of body weight on selexipag CL/F (unitless)")                  # Table 1 pediatric column: beta_CL(WEIGHTBL) = 0.828, RSE 12.1%; Results note this is close to the standard allometric value of 0.75
    e_tbili_cl     <- -0.363;  label("Power exponent of total bilirubin on selexipag CL/F (unitless)")              # Table 1 pediatric column: beta_CL(BILIBL) = -0.363, RSE 22.5%
    e_wt_vc        <-  0.805;  label("Power exponent of body weight on selexipag Vp/F (unitless)")                  # Table 1 pediatric column: beta_Vp(WEIGHTBL) = 0.805, RSE 11.3%
    e_wt_vc_act    <-  0.542;  label("Power exponent of body weight on JNJ-68006861 Vm/F (unitless)")               # Table 1 pediatric column: beta_Vm(WEIGHTBL) = 0.542, RSE 11.4%; Results note this is ~30% smaller than the adult value of 0.803
    e_sexf_km_act  <-  0.145;  label("Log-scale effect of male sex on JNJ-68006861 km (unitless)")                  # Table 1 pediatric column: beta_km(SEXN_0) = 0.145, RSE 29.1%; carried on the male indicator (1 - SEXF) because females are the reference
    e_era_km_act   <-  0.186;  label("Log-scale effect of ERA-only PAH comedication on JNJ-68006861 km")            # Table 1 pediatric column: beta_km(COPAH1_1) = 0.186, RSE 34.8%
    e_pde5_km_act  <-  0.0495; label("Log-scale effect of PDE5I-only PAH comedication on JNJ-68006861 km")          # Table 1 pediatric column: beta_km(COPAH2_1) = 0.0495, RSE 116% (imprecise but retained to keep the four-level categorical intact)
    e_combo_km_act <-  0.368;  label("Log-scale effect of ERA + PDE5I PAH comedication on JNJ-68006861 km")         # Table 1 pediatric column: beta_km(COPAH3_1) = 0.368, RSE 15.3%

    # ---------------- IIV (variances on the log / logit scale) ----------------
    # Table 1's note states "omega values reported as standard deviation", and
    # the control stream $OMEGA / $OMEGAP blocks carry the STANDARD keyword,
    # confirming the tabulated omegas are SDs. ini() takes variances, so each
    # entry below is omega^2 with the source omega quoted in the comment.
    etalka        ~ 0.179776    # 0.424^2;  Table 1 pediatric column: omega(ka) = 0.424, RSE 3.27%
    etalkmet      ~ 0.00698896  # 0.0836^2; Table 1 pediatric column: omega(kmet) = 0.0836, RSE 3.22% (76.2% shrinkage)
    etalcl        ~ 0.646416    # 0.804^2;  Table 1 pediatric column: omega(CL) = 0.804, RSE 3.02%
    etalvc        ~ 0.125316    # 0.354^2;  Table 1 pediatric column: omega(Vp) = 0.354, RSE 3.05%
    etalk12       ~ 0.1936      # 0.44^2;   Table 1 pediatric column: omega(kx12) = 0.44, RSE 3.22%
    etalk21       ~ 1.2321      # 1.11^2;   Table 1 pediatric column: omega(kx21) = 1.11, RSE 3.15%
    etalkm_act    ~ 0.072361    # 0.269^2;  Table 1 pediatric column: omega(km) = 0.269, RSE 3.06%
    etalk34_act   ~ 0.097344    # 0.312^2;  Table 1 pediatric column: omega(kx34) = 0.312, RSE 3.13%
    etalk43_act   ~ 0.6561      # 0.81^2;   Table 1 pediatric column: omega(kx43) = 0.81, RSE 3.2%
    etalvc_act    ~ 0.00848241  # 0.0921^2; Table 1 pediatric column: omega(Vm) = 0.0921, RSE 3.06% (79.6% shrinkage)
    etalogittlag  ~ 3.1329      # 1.77^2;   Table 1 pediatric column: omega(Tlag1Half) = 1.77, RSE 3.17%, distribution LogitNormal -- this SD is on the logit scale

    # ---------------- Residual error (proportional, by output) ----------------
    # The control stream $ERROR uses Y = IPRED + EPS*IPRED for both outputs
    # (pure proportional), and $SIGMA carries variances that are the squares of
    # the printed error_PROP values -- confirming that the error_PROP values in
    # Table 1 are proportional standard deviations.
    propSd     <- 0.694; label("Selexipag proportional residual SD (fraction)")                                     # Table 1 pediatric column: error_PROP1 = 0.694, RSE 4.28% (epsilon shrinkage 6.53%)
    propSd_act <- 0.455; label("JNJ-68006861 proportional residual SD (fraction)")                                  # Table 1 pediatric column: error_PROP2 = 0.455, RSE 3.62% (epsilon shrinkage 5.99%)
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
