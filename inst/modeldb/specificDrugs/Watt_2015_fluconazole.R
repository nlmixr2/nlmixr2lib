Watt_2015_fluconazole <- function() {
  description <- paste(
    "One-compartment population PK model for intravenous fluconazole in",
    "critically ill children (1 day to 17 years; n=40) supported with",
    "extracorporeal membrane oxygenation (ECMO) or matched non-ECMO",
    "controls (Watt 2015). Clearance scales linearly with body weight and",
    "is modulated by serum creatinine via a power function (CREAT/0.4)^-0.29",
    "centered at the cohort median initial SCR of 0.4 mg/dL; allometric 3/4-",
    "power scaling on CL was tested and rejected (delta-OFV +9.7) so a",
    "linear weight scaling was retained. Central volume scales linearly with",
    "body weight and is increased 1.39-fold in subjects on ECMO support via",
    "a multiplicative power factor 1.39^ECMO_STATUS. Proportional residual",
    "error (15.3% CV); diagonal Omega with IIV on CL (33.2% CV) and V (22.2%",
    "CV) only -- the paper retained the proportional-only error model after",
    "showing the proportional-plus-additive form could not precisely",
    "estimate the additive component. Used by the authors to derive dosing",
    "recommendations for invasive candidiasis prevention (12 mg/kg loading",
    "then 6 mg/kg daily) and treatment (35 mg/kg loading then 12 mg/kg",
    "daily) in children on ECMO."
  )
  reference <- paste(
    "Watt KM, Gonzalez D, Benjamin DK Jr, Brouwer KLR, Wade KC,",
    "Capparelli E, Barrett J, Cohen-Wolkowiez M. 2015. Fluconazole",
    "population pharmacokinetics and dosing for prevention and treatment",
    "of invasive candidiasis in children supported with extracorporeal",
    "membrane oxygenation. Antimicrob Agents Chemother 59(7):3935-3943.",
    "doi:10.1128/AAC.00102-15"
  )
  vignette <- "Watt_2015_fluconazole"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at the value on the day of the first study dose",
        "(Watt 2015 Methods, 'Population PK analysis': 'PNA and weight",
        "were calculated on the day of the first dose of study drug, and",
        "those values were imputed forward'). Linear weight scaling on CL",
        "and V (no allometric exponent); the 3/4-power exponent on CL was",
        "tested and rejected in model development with delta-OFV +9.7",
        "(Watt 2015 Results, 'Population PK model development': 'Allometric",
        "scaling of weight (3/4 power) on CL did not improve model fit and",
        "increased the objective function value by 9.7 points. Similarly,",
        "the use of a sigmoidal maximum effect (Emax) maturation",
        "relationship between postmenstrual age and CL resulted in an",
        "increase in the objective function value by 4.8 points.",
        "Consequently, weight was scaled to the power of 1 for both CL and",
        "V.'). Cohort median 3.4 kg, range 1.9-77 kg (Table 1)."
      ),
      source_name        = "wt"
    ),
    CREAT = list(
      description        = "Serum creatinine",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying within subject (Watt 2015 Methods, 'Population PK",
        "analysis': 'For children with multiple measurements of SCR,",
        "albumin, AST, or ALT, values were allowed to change with time.').",
        "Enters CL via power scaling (CREAT/0.4)^-0.29 with reference 0.4",
        "mg/dL (cohort median initial SCR; Watt 2015 Methods: 'All",
        "continuous variables were centered using the median value'; Table",
        "1 'Initial' SCR median = 0.4 mg/dL; Table 2 univariable analysis",
        "footnote 'CL = theta_CL * wt * (creatinine/0.4)^-0.29'). Note: the",
        "Watt 2015 Results section main text reports the final-model CL",
        "equation with a centering value of 0.5 mg/dL ('CL (in liters per",
        "hour) was calculated as 0.019 * weight * (SCR/0.5)^-0.29 *",
        "exp(eta_CL)'), but the abstract, the Table 2 footnote on the",
        "univariable creatinine model, the Methods median-centering rule,",
        "and the Table 1 initial-SCR median all consistently report 0.4;",
        "the 0.5 in the Results section main text is therefore a",
        "transcription typo and the model is encoded with the consistent",
        "0.4 mg/dL reference. See vignette Errata."
      ),
      source_name        = "SCR"
    ),
    ECMO_STATUS = list(
      description        = "Extracorporeal membrane oxygenation (ECMO) support indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "1 = the subject was receiving ECMO support during the modelled",
        "period; 0 = no ECMO support. Per-subject in the source data;",
        "time-fixed at the subject level (Watt 2015 Methods 'Population PK",
        "analysis' lists ECMO support as a screened covariate; Table 1",
        "reports 21 of 40 subjects on ECMO support and 19 of 40 not on",
        "ECMO). Enters V via a multiplicative power factor",
        "e_ecmo_status_vc^ECMO_STATUS so that V increases 39% for ECMO",
        "subjects (V = 0.93 * WT * 1.39^ECMO; Watt 2015 Table 3,",
        "'Coefficient for ECMO on V' = 1.39 with %RSE 7.8). Hemofiltration",
        "(univariable delta-OFV -10.61 on V) and CVVHD were screened in",
        "the univariable analysis on V but were collinear with ECMO",
        "support (all 5 children on hemofiltration were on ECMO; the 2 on",
        "CVVHD were among the hemofiltration subjects); ECMO support",
        "subsumed both effects and was retained alone in the final",
        "multivariable model. ECMO prime volume and the ratio of prime",
        "volume to the patient's native blood volume were also tested as",
        "alternatives to the binary indicator but did not improve fit",
        "more than the binary ECMO_STATUS."
      ),
      source_name        = "ECMO"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 40L,
    n_studies      = 3L,
    age_range      = "1 day - 17 years",
    age_median     = "22 days",
    weight_range   = "1.9-77 kg",
    weight_median  = "3.4 kg",
    sex_female_pct = 35.0,
    race_ethnicity = c(White = 43, Black = 45, Other = 12),
    disease_state  = paste(
      "Critically ill children at high risk for invasive candidiasis, with",
      "or without ECMO support. 23 of 40 (57%) received fluconazole for",
      "prophylaxis and 17 of 40 (43%) for treatment of suspected fungal",
      "infection. 21 of 40 (53%) were on ECMO; of those, 5 had concomitant",
      "hemofiltration during PK sample collection and 2 of those 5 also",
      "subsequently required CVVHD. None of the children developed",
      "culture-confirmed invasive candidiasis during the study. Cohort",
      "median initial serum creatinine 0.4 mg/dL (range 0.1-1.3),",
      "max-during-study 0.6 mg/dL (range 0.1-3.2); ECMO subjects had",
      "higher max SCR than non-ECMO (0.7 vs 0.5 mg/dL, p = 0.03) but the",
      "initial SCR difference was not significant (0.5 vs 0.3, p = 0.13)."
    ),
    dose_range     = paste(
      "Intravenous fluconazole; first dose 2.7-26.5 mg/kg (cohort median",
      "25 mg/kg). Three constituent prospective trials: study 1 (n=20",
      "ECMO subjects who received 25 mg/kg once weekly for prophylaxis or",
      "standard-of-care dosing for presumed fungal infection); study 2",
      "(n=12 critically ill infants <1 year of age, of whom 1 was on",
      "ECMO, who received a 25 mg/kg loading dose IV once followed by",
      "daily 12 mg/kg maintenance); study 3 (n=8 infants 23-42 weeks",
      "gestational age at birth and <120 days of age; only the >=36 weeks",
      "GA subset was included to limit PK variability from prematurity)."
    ),
    pma_range      = "35-76 weeks postmenstrual age (median 41; reported for n=33 infants <1 year)",
    ga_range       = "24-41 weeks gestational age (reported for n=33 infants <1 year)",
    regions        = paste(
      "United States. Constituent studies were enrolled at Duke",
      "University Medical Center (study 1, single-center ECMO study), at",
      "a single-center pediatric ICU (study 2, infant fluconazole loading",
      "dose), and at a multicenter neonatal network (study 3, infant",
      "antifungal PK)."
    ),
    notes          = paste(
      "Pooled cohort from three prospective fluconazole PK trials. PK",
      "samples: 360 plasma concentrations included in the population PK",
      "analysis (Watt 2015 Results, 'Study infants and PK specimens';",
      "median 8 per child, range 1-22). 55 of 360 (15%) were scavenge",
      "samples (study 1 n=32, study 2 n=6, study 3 n=17). Median 10 SCR",
      "samples per child during the PK collection period (range 1-23);",
      "all children had an SCR sample at the time of the first dose.",
      "Baseline demographics from Watt 2015 Table 1; Bayesian estimates",
      "of CL and V stratified by age and ECMO status are reported in",
      "Table 4. NONMEM 7.2 with PsN 3.6.2 and Pirana 2.8.0; FOCEI;",
      "bootstrap n=1000 (999 successful runs)."
    )
  )

  ini({
    # -----------------------------------------------------------------------
    # STRUCTURAL PK PARAMETERS -- Watt 2015 Table 3 (Final population PK
    # model parameter estimates). The final-model equations are reported in
    # both the abstract and in the Results section, 'Population PK model
    # development':
    #     V  = 0.93 * weight * 1.39^ECMO  * exp(eta_V)
    #     CL = 0.019 * weight * (SCR/0.4)^-0.29 * exp(eta_CL)
    # Weight enters linearly (allometric 3/4-power on CL was tested and
    # rejected with delta-OFV +9.7; Methods 'allometric scaling of weight
    # (3/4 power) on CL did not improve model fit'). Linear weight scaling
    # implies CL and V are reported in L/h/kg and L/kg respectively (Watt
    # 2015 Table 3 'V (liters/kg)' and 'CL (liters/h/kg)').
    # -----------------------------------------------------------------------
    lcl <- log(0.019); label("Typical CL per kg body weight (L/h/kg)")  # Watt 2015 Table 3: CL = 0.019 L/h/kg (%RSE 5.6)
    lvc <- log(0.93);  label("Typical V  per kg body weight (L/kg)")    # Watt 2015 Table 3: V  = 0.93  L/kg   (%RSE 5.8)

    # -----------------------------------------------------------------------
    # COVARIATE EFFECTS -- Watt 2015 Table 3.
    #
    # e_ecmo_status_vc is the multiplicative power factor on V when on
    # ECMO; the model body applies e_ecmo_status_vc^ECMO_STATUS so the
    # factor is 1 for non-ECMO subjects (1^0) and 1.39 for ECMO subjects
    # (1.39^1). Goti 2018 vancomycin (HEMODIAL on CL and Vc) and Patel
    # 2011 fluconazole (FILT_AGE_HI on CL_CVVHDF) are the closest
    # precedents for binary power-of-indicator encodings.
    #
    # e_creat_cl is the power exponent on the centered serum-creatinine
    # ratio (CREAT/0.4); the model body applies (CREAT/0.4)^e_creat_cl.
    # The negative sign of the exponent encodes the physiologic
    # expectation that CL decreases as SCR rises (Watt 2015 Discussion:
    # 'CL decreased with increasing SCR level. This relationship also was
    # expected, as fluconazole is primarily excreted by the kidneys').
    # -----------------------------------------------------------------------
    e_ecmo_status_vc <- 1.39;  label("Multiplicative power factor on V for ECMO_STATUS = 1 (vs 0)")  # Watt 2015 Table 3: Coefficient for ECMO on V = 1.39 (%RSE 7.8)
    e_creat_cl       <- -0.29; label("Power exponent on (CREAT/0.4) for CL")                          # Watt 2015 Table 3: Exponent for creatinine on CL = -0.29 (%RSE 9.9)

    # -----------------------------------------------------------------------
    # INTER-INDIVIDUAL VARIABILITY -- Watt 2015 Table 3 (Random effects, CV%
    # column). Methods state: 'An exponential model for interindividual
    # variance was used.' Under the lognormal exponential interpretation,
    # %CV maps to log-scale variance as omega^2 = log(1 + CV^2). See
    # references/verification-checklist.md item A 'CV% vs. variance'.
    #
    #     IIV V  = 22.2 % CV -> omega^2 = log(1 + 0.222^2)
    #     IIV CL = 33.2 % CV -> omega^2 = log(1 + 0.332^2)
    #
    # The eta shrinkage was reported as low (4.6% for CL and 7.0% for V),
    # supporting the precision of the IIV estimates (Watt 2015 Results,
    # 'Population PK model development'). CL and V were modelled with a
    # diagonal Omega matrix: 'CL and V were not correlated, and use of a
    # covariance term between CL and V did not improve the model fit.'
    # -----------------------------------------------------------------------
    etalcl ~ log(1 + 0.332^2)  # Watt 2015 Table 3: CL IIV = 33.2 % CV (%RSE 21.3)
    etalvc ~ log(1 + 0.222^2)  # Watt 2015 Table 3: V  IIV = 22.2 % CV (%RSE 28.6)

    # -----------------------------------------------------------------------
    # RESIDUAL ERROR -- Watt 2015 Table 3. Proportional only: 'Residual
    # variability was best described by a proportional error model. While
    # a proportional-plus-additive error model resulted in a significant
    # drop in the objective function, we were unable to precisely estimate
    # the additive error component. Because goodness-of-fit plots and
    # estimates were virtually identical between the two error models, we
    # used the proportional error model in the final model.' The eps
    # shrinkage was low (9.6%), supporting model identifiability.
    # 15.3 % CV maps to a proportional SD of 0.153 on the fraction scale.
    # -----------------------------------------------------------------------
    propSd <- 0.153; label("Proportional residual error (fraction)")  # Watt 2015 Table 3: Residual proportional error = 15.3 % CV (%RSE 13.1)
  })

  model({
    # ------------------------------------------------------------------
    # Individual PK parameters.
    #
    # CL scales linearly with body weight and is modulated by serum
    # creatinine via a centered power function. V scales linearly with
    # body weight and is multiplied by 1.39 for subjects on ECMO support
    # (via the power-of-indicator encoding e_ecmo_status_vc^ECMO_STATUS,
    # which evaluates to 1 when ECMO_STATUS = 0 and to 1.39 when
    # ECMO_STATUS = 1).
    # ------------------------------------------------------------------
    cl <- exp(lcl + etalcl) * WT * (CREAT / 0.4)^e_creat_cl
    vc <- exp(lvc + etalvc) * WT * e_ecmo_status_vc^ECMO_STATUS

    # Micro-constant.
    kel <- cl / vc

    # ------------------------------------------------------------------
    # ODE system.
    #
    # One-compartment IV-infusion disposition: a single central
    # compartment with first-order elimination. Dose records target
    # central with amt = dose (mg) and rate or dur for the infusion
    # duration (typical clinical fluconazole IV infusions are over
    # 1-2 hours; Watt 2015 sampling windows quote 'after the end of the
    # infusion').
    # ------------------------------------------------------------------
    d/dt(central) <- -kel * central

    # ------------------------------------------------------------------
    # Observation.
    # Cc = plasma concentration in mg/L (dose mg, vc L -> mg/L).
    # ------------------------------------------------------------------
    Cc <- central / vc

    Cc ~ prop(propSd)
  })
}
