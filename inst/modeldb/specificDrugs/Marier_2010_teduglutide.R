Marier_2010_teduglutide <- function() {
  description <- "One-compartment population PK model with first-order subcutaneous absorption and a fixed absorption lag time for the GLP-2 analog teduglutide, fit to pooled data from eight phase I/II/III studies in healthy participants and in patients with short bowel syndrome, Crohn's disease, or moderate hepatic impairment (Marier 2010). Site-specific typical absorption rate constant (abdomen vs pooled arm/thigh), a power effect of body weight on both apparent clearance and apparent central volume, and a categorical effect of sex on apparent clearance."
  reference   <- "Marier JF, Mouksassi MS, Gosselin NH, Beliveau M, Cyran J, Wallens J. Population pharmacokinetics of teduglutide following repeated subcutaneous administrations in healthy participants and in patients with short bowel syndrome and Crohn's disease. J Clin Pharmacol. 2010;50(1):36-49. doi:10.1177/0091270009342252"
  vignette    <- "Marier_2010_teduglutide"
  units       <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Both entries verified against Marier 2010: the dose is
  # administered subcutaneously into the abdomen, arm or thigh (Table I "Dose
  # and Route"), and teduglutide was quantified in plasma by LC/MS/MS with an
  # analytical range of 1.00-120 ng/mL (Methods, "Blood Samples and Analytical
  # Assay").
  compartmentData <- list(
    depot   = list(analyte = "teduglutide", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "teduglutide", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight (baseline)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline. Power exponents on Vc/F (2.37) and CL/F (0.323); reference 70.9 kg = the cohort MEDIAN weight of Table II (mean is 71.5 kg -- Table VI writes the denominator as 70.9, and reproducing the Table VI half-lives and the Table VII AUC/Cmax grid confirms the median is the normalising value). The 2.37 exponent is far above any allometric expectation and the authors flag it as such (Results, 'Covariates Analysis': the power function of weight on Vc/F 'suggest[s] a very important effect of body weight on Vc/F'); it is an empirical descriptor over the observed 40.1-127 kg range and should not be extrapolated. The source paper labels this column WT.",
      source_name        = "WT"
    ),
    SEXF = list(
      description        = "Biological sex indicator: 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Time-fixed. Marier 2010 uses male as the reference: Table VI gives CL/F = 12.4 * (WT/70.9)^0.323 for men and 10.5 * (WT/70.9)^0.323 for women, and Table V reports the 'Sex on CL/F' coefficient as the multiplicative ratio 0.843 (12.4 * 0.843 = 10.45, printed as 10.5). The canonical 1 = female orientation therefore needs no sign inversion: the effect is applied as e_sexf_cl^SEXF, which is 1 for men and 0.843 for women. The paper's own framing is that men have ~18% HIGHER CL/F than women (Abstract; Results, 'Covariates Analysis').",
      source_name        = "SEX"
    ),
    INJSITE_ARM = list(
      description        = "SC injection-site indicator: 1 = arm, 0 = abdomen",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (abdomen)",
      notes              = "Per-dose-record covariate. Marier 2010 tested site of administration on ka, on the lag time and on relative bioavailability (Table III mod 4/5/6) and retained it on ka only; the arm and thigh rate constants were then POOLED into a single 'other' category (Table III mod 16, dMOF = +1.552 vs the three-category model, i.e. no significant loss of fit). This model therefore reads INJSITE_ARM and INJSITE_THIGH as a pooled non-abdomen indicator selecting lka_armthigh in place of lka_abdomen; the two indicators share one estimated typical value because the paper estimated one. Study CL0600-015 is the only study that dosed outside the abdomen for PK (10 mg SC in abdomen, thigh and arm, 3-way crossover), so a subject can switch sites between occasions -- supply the indicator on each dose record.",
      source_name        = "Site of administration (Table III narrative; abdomen vs arm vs thigh, arm and thigh pooled in the final model)"
    ),
    INJSITE_THIGH = list(
      description        = "SC injection-site indicator: 1 = thigh, 0 = abdomen",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (abdomen)",
      notes              = "Per-dose-record covariate, pooled with INJSITE_ARM into the single non-abdomen absorption stratum -- see the INJSITE_ARM notes. Both indicators 0 denotes the abdomen reference. Note that the phase II study CL0600-008 and the phase III study CL0600-004 allowed self-administration 'into the abdomen or the thigh' (Methods, 'Clinical Studies'), so the thigh arm of the pooled stratum carries most of the patient (as opposed to healthy-volunteer) data.",
      source_name        = "Site of administration (Table III narrative; abdomen vs arm vs thigh, arm and thigh pooled in the final model)"
    )
  )

  # Covariates that Marier 2010 SCREENED but did not retain in the final model.
  # These are documentation only -- none is referenced in model(). They are
  # recorded because the paper's headline clinical conclusion ("no dose
  # adjustment was deemed necessary in patients with altered renal or liver
  # function", Abstract) rests entirely on this screen coming back negative.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age", units = "y", type = "continuous",
      notes = "Table II: mean 39.1 y, median 38 y, range 18-79 y. Table IV step 1 mod1 tested age on CL/F (dMOF = -10.074, statistically significant at the 1% forward-inclusion threshold) but it was not selected over weight on Vc/F and dropped out by step 4; no age effect survives in Table VI."
    ),
    BMI = list(
      description = "Body mass index", units = "kg/m^2", type = "continuous",
      notes = "Explored graphically only (Figure 2, relationships between PK parameters and demographic covariates); never entered the formal stepwise covariate model of Table IV and no point estimate is reported."
    ),
    CRCL_BASE = list(
      description = "Baseline creatinine clearance (NOT BSA-normalized)", units = "mL/min", type = "continuous",
      notes = "Table II: mean 109 mL/min, median 107, range 26.0-236, n = 264. Entered Table IV as 'CLCRT', creatinine clearance TRUNCATED to 150 mL/min (Table III footnote). Significant in isolation on both CL/F (mod2, dMOF = -25.22) and Vc/F (mod3, dMOF = -63.893) and still significant on CL/F at step 3 (mod18, dMOF = -10.267), but not retained: Results, 'Covariates Analysis' states 'Creatinine clearance and markers of hepatic impairment were not detected as significant components explaining the variability of PK parameters of teduglutide' and 'No covariates were removed during the backward elimination step'. Units are raw mL/min here, not the mL/min/1.73 m^2 of the general CRCL canonical."
    ),
    CREAT = list(
      description = "Serum creatinine", units = "mg/dL", type = "continuous",
      notes = "Table II: mean 0.923 mg/dL, range 0.50-2.0, n = 264. Reported as a demographic descriptor and used to derive the renal-impairment strata (60/264 participants with some degree of renal impairment; 45 mild, 14 moderate, 1 severe); never tested directly as a PK covariate."
    ),
    TBILI = list(
      description = "Total bilirubin", units = "mg/dL", type = "continuous",
      notes = "Table II: mean 0.561 mg/dL, range 0.10-7.8, n = 255; 5 participants (2%) above the 1.4 mg/dL normal limit. Screened graphically in Figure 3 as a marker of hepatic impairment. Results, 'Covariates Analysis': 'Because no relationship was observed with PK parameters of teduglutide, markers of hepatic impairment were not tested in the population PK model' -- so no coefficient exists for it."
    ),
    AST = list(
      description = "Aspartate aminotransferase", units = "U/L", type = "continuous",
      notes = "Table II: mean 28.5 U/L, range 10-297, n = 264; 37 participants (14.0%) above the 5-40 U/L normal range. Screened graphically in Figure 3 only; never entered the formal covariate model (see the TBILI note)."
    ),
    ALT = list(
      description = "Alanine aminotransferase", units = "U/L", type = "continuous",
      notes = "Table II: mean 31.6 U/L, range 6.0-412, n = 264; 24 participants (9.1%) above the 5-60 U/L normal range. Screened graphically in Figure 3 only; never entered the formal covariate model (see the TBILI note)."
    ),
    GGT = list(
      description = "Gamma-glutamyltransferase", units = "U/L", type = "continuous",
      notes = "Table II: mean 44.5 U/L, range 7.0-382, n = 235; 36 participants (15.3%) above 78 U/L. Screened graphically in Figure 3 only; never entered the formal covariate model (see the TBILI note)."
    ),
    HEPIMP_MOD = list(
      description = "Moderate hepatic impairment indicator: 1 = moderate, 0 = normal", units = "(binary)", type = "binary",
      notes = "Child-Pugh class B (total score 7-9) per Methods, 'Covariate Analysis'. Carried by the 24 participants of the dedicated phase I hepatic-impairment study CL0600-017. Tested as a dichotomous covariate on CL/F (Table IV mod8 dMOF = -5.382 NS; mod16 dMOF = -2.105 NS; mod20 dMOF = -2.631 NS) and on Vc/F (mod17 dMOF = -0.444 NS; mod21 dMOF = -0.465 NS) and rejected at every step."
    ),
    RACE_BLACK = list(
      description = "Black / African American race indicator", units = "(binary)", type = "binary",
      notes = "Results, 'Summary of Data': 12.8% African American. Race was tested on Vc/F only, as a single grouped term rather than as per-group indicators, and was non-significant at every step (Table IV mod4 dMOF = -0.824; mod12 dMOF = 0; mod22 dMOF = -0.002). The paper reports no per-race coefficient, so no indicator can be reconstructed."
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator", units = "(binary)", type = "binary",
      notes = "Results, 'Summary of Data': 0.75% Asian (i.e. 2 of 266 participants) -- the cohort cannot support an Asian-specific estimate. Pooled into the same non-significant grouped race term as RACE_BLACK; see that note."
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 266L,
    n_studies        = 8L,
    age_range        = "18 - 79 years",
    age_median       = "38 years",
    age_mean         = "39.1 years",
    weight_range     = "40.1 - 127 kg",
    weight_median    = "70.9 kg",
    weight_mean      = "71.5 kg",
    sex_female_pct   = 38.3,
    race_ethnicity   = c(White = 83.8, `African American` = 12.8, Other = 2.63, Asian = 0.75),
    disease_state    = "Mixed: healthy participants (5 phase I studies), patients with moderately active Crohn's disease of at least 6 months' history (phase II CL0600-008), patients with parenteral-nutrition-dependent short bowel syndrome secondary to vascular ischemic disease, malrotation or volvulus (phase II ALX-0600-92001 and phase III CL0600-004), and participants with moderate hepatic impairment (phase I CL0600-017).",
    dose_range       = "Daily subcutaneous teduglutide, 2.5 to 80 mg as fixed doses (2.5, 5, 7, 10, 15, 20, 25, 30, 50, 80 mg) and 0.03 to 0.20 mg/kg/day as weight-based doses, injected into the abdomen, arm or thigh. Treatment durations ranged from a single dose to 24 weeks.",
    renal_function   = "Creatinine clearance 26.0 - 236 mL/min (n = 264). 60 participants (22.7%) had creatinine clearance in an impaired range: 45 mild (17.0%), 14 moderate (5.30%), 1 severe (0.38%), classified per regulatory guidance on PK studies in renal impairment.",
    hepatic_function = "Includes 24 participants with moderate (Child-Pugh class B, score 7-9) hepatic impairment from the dedicated study CL0600-017. Across the whole cohort, 37 participants (14.0%) had AST and 24 (9.1%) had ALT above the normal limit, 36 (15.3%) had GGT > 78 U/L, and 5 (2%) had bilirubin > 1.4 mg/dL.",
    regions          = "Not reported. Sponsor NPS Pharmaceuticals (Parsippany, NJ, USA); bioanalysis at Allelix Biopharmaceuticals (Mississauga, ON, Canada) and the NPS DMPK Department (Salt Lake City, UT, USA).",
    studies          = "ALX-0600 1621/1613 (phase I SAD, n = 32), CL0600-006 (phase I IV/SC bioavailability crossover, n = 14), CL0600-015 (phase I 3-way injection-site crossover, n = 18), CL0600-017 (phase I moderate hepatic impairment, n = 24), CL0600-022 (phase I MAD, n = 95), CL0600-008 (phase II Crohn's disease, n = 100), ALX-0600-92001 (phase II short bowel syndrome, n = 17), CL0600-004 (phase III short bowel syndrome, n = 83). Table I.",
    notes            = "Demographics from Table II (n = 266; 164 men / 102 women). Estimation by NONMEM VI with FOCE-with-interaction. Plasma teduglutide by LC/MS/MS, analytical range 1.00-120 ng/mL, LOQ 1.00 ng/mL. Most concentrations fell below the limit of quantitation by 24 h postdose and only 60 of 183 predose samples (32.8%) after repeated dosing were measurable, so accumulation on daily dosing is minimal. The final model was confirmed by a nonparametric bootstrap (899 successfully minimised runs); every median bootstrap estimate in Table V is within 2% of the corresponding original estimate. The Abstract's participant count of 256 disagrees with the 266 of the Results and of Table II -- see the vignette Errata."
  )

  ini({
    # ---- Structural PK parameters (Marier 2010 Table V / Table VI) ----
    # Reference subject: 70.9 kg (cohort median weight, Table II), male
    # (SEXF = 0). Absorption is stratified into two parallel typical values
    # rather than a reference plus an offset, because Table V reports Ka_abd
    # and Ka_other as two separate typical values each with its own bootstrap
    # 95% percentile interval, under "Population PK Parameters" rather than
    # under "Covariate model". Both strata share the single BSV on Ka that
    # Table V reports (39.5%), so etalka keeps its bare canonical name.
    lka_abdomen  <- log(0.299); label("First-order SC absorption rate constant for abdominal injection (Ka_abd, 1/h)")        # Marier 2010 Table V: Ka_abd = 0.299 1/h (bootstrap 0.280-0.320)
    lka_armthigh <- log(0.206); label("First-order SC absorption rate constant for arm or thigh injection (Ka_other, 1/h)")   # Marier 2010 Table V: Ka_other = 0.206 1/h (bootstrap 0.171-0.242)
    lcl          <- log(12.4);  label("Apparent clearance at reference weight for men (CL/F, L/h)")                           # Marier 2010 Table V: CL/F = 12.4 L/h (bootstrap 11.9-13.0); Table VI
    lvc          <- log(32.8);  label("Apparent central volume of distribution at reference weight (Vc/F, L)")                # Marier 2010 Table V: Vc/F = 32.8 L (bootstrap 30.6-34.9); Table VI

    # Absorption lag time. Marier 2010 fixed both the typical value and its
    # BSV: Results, "Structural Model Buildup" -- "The absorption lag time
    # (ALAG) as well as BSV were fixed to population values previously
    # estimated with the model to stabilize the model and facilitate
    # convergence" -- and Table VI prints "0.208 (Fixed)" / "31.6 (fixed)".
    # Consistent with this, ALAG carries no row in the Table V bootstrap.
    ltlag <- fixed(log(0.208)); label("Absorption lag time (ALAG, h)")                                                        # Marier 2010 Table VI: ALAG = 0.208 h (Fixed)

    # ---- Covariate effects (Marier 2010 Table V "Covariate model"; Table VI) ----
    # Both weight effects are power functions on WT / 70.9 (Methods, "Covariate
    # Analysis": "Continuous covariates were included in the structural model
    # with a power function"). Sex enters CL/F as a multiplicative ratio.
    e_wt_vc   <- 2.37;  label("Power exponent of WT on Vc/F (unitless)")                              # Marier 2010 Table V: body weight on Vc/F = 2.37 (bootstrap 2.08-2.67); Table VI
    e_wt_cl   <- 0.323; label("Power exponent of WT on CL/F (unitless)")                              # Marier 2010 Table V: body weight on CL/F = 0.323 (bootstrap 0.102-0.521); Table VI
    e_sexf_cl <- 0.843; label("Multiplicative ratio of CL/F for women relative to men (unitless)")    # Marier 2010 Table V: sex on CL/F = 0.843 (bootstrap 0.764-0.930); Table VI 12.4 -> 10.5 L/h

    # ---- IIV (Marier 2010 Table V "Between-subject variability") ----
    # Reported as "BSV, %". Methods, "Population PK Model Development": "The
    # between-subject variability (BSV) was modeled as an exponential random
    # effect model ... assumed to follow the lognormal distribution", so the
    # percentages are coefficients of variation on the lognormal and the
    # internal variance is omega^2 = log(CV^2 + 1). The vignette Errata records
    # the arbitration that rules out reading these percentages as raw omega^2
    # values, and the size of the residual exact-vs-approximate-CV ambiguity.
    etalka  ~ 0.14499   # Marier 2010 Table V BSV on Ka   = 39.5% -> log(0.395^2 + 1)
    etalcl  ~ 0.07600   # Marier 2010 Table V BSV on CL/F = 28.1% -> log(0.281^2 + 1)
    etalvc  ~ 0.14430   # Marier 2010 Table V BSV on Vc/F = 39.4% -> log(0.394^2 + 1)
    # The BSV on the lag time was held constant alongside its typical value
    # (Table VI prints it as "31.6 (fixed)"), and it likewise has no Table V
    # bootstrap row. rxode2 adopts an eta's trailing comment as its label, so
    # that comment must not restate what the fixed() wrapper already says.
    etaltlag ~ fixed(0.09518)   # Marier 2010 Table VI BSV on ALAG = 31.6% -> log(0.316^2 + 1)

    # ---- Residual error (Marier 2010 Table V "Error model") ----
    # Combined additive + proportional, per Results, "Structural Model
    # Buildup": "residual variability in plasma concentrations of teduglutide
    # was fitted using a proportional and additive error model". Both rows are
    # printed in SD-scale units ("ng/mL" and "%"), not as variances.
    addSd  <- 6.80;  label("Additive residual error SD on plasma teduglutide Cc (ng/mL)")   # Marier 2010 Table V: additive error = 6.80 ng/mL (bootstrap 5.47-8.06)
    propSd <- 0.243; label("Proportional residual error on plasma teduglutide Cc (fraction)") # Marier 2010 Table V: proportional error = 24.3% (bootstrap 21.9-26.5)
  })

  model({
    # ---- 1. Derived covariate term ----
    # Marier 2010 pooled arm and thigh into a single non-abdomen absorption
    # stratum (Table III mod 16). INJSITE_ARM and INJSITE_THIGH are mutually
    # exclusive by the register definition of the INJSITE_<site> family, with
    # both zero denoting the abdomen reference, so their sum is the 0/1
    # non-abdomen indicator.
    injsite_nonabd <- INJSITE_ARM + INJSITE_THIGH

    # ---- 2. Individual PK parameters ----
    # The absorption rate constant switches between the two parallel typical
    # values on the log scale; the shared etalka preserves a subject's relative
    # ka deviation across a change of injection site.
    ka <- exp(lka_abdomen * (1 - injsite_nonabd) + lka_armthigh * injsite_nonabd + etalka)

    cl <- exp(lcl + etalcl) * (WT / 70.9)^e_wt_cl * e_sexf_cl^SEXF
    vc <- exp(lvc + etalvc) * (WT / 70.9)^e_wt_vc

    tlag <- exp(ltlag + etaltlag)

    # ---- 3. Micro-constants ----
    kel <- cl / vc

    # ---- 4. ODE system ----
    # depot and central hold teduglutide amounts in mg; the SC dose enters
    # depot in mg. CL/F and Vc/F are apparent, so no separate bioavailability
    # parameter is carried (Marier 2010 tested relative bioavailability by
    # injection site as Table III mod 6 and did not retain it).
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ---- 5. Lag time ----
    alag(depot) <- tlag

    # ---- 6. Observation and error ----
    # central / vc is mg/L = ug/mL; multiply by 1000 to report ng/mL, the units
    # of the assay (analytical range 1.00-120 ng/mL) and of addSd.
    Cc <- central / vc * 1000

    Cc ~ add(addSd) + prop(propSd)
  })
}
