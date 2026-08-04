Walsh_2024_buprenorphine_drugLiking <- function() {
  description <- paste(
    "Direct maximum-inhibition (Imax) exposure-response model of opioid",
    "drug-liking blockade by subcutaneous long-acting buprenorphine",
    "(CAM2038 Q1W) in non-treatment-seeking adults with moderate-to-severe",
    "opioid use disorder. The endpoint is the period-corrected drug liking",
    "maximum-effect (Emax) bipolar visual analog scale (VAS) score measured",
    "after an intramuscular hydromorphone 18 mg challenge; the driver is the",
    "time-matched buprenorphine plasma concentration CP_BPN_NGML sampled",
    "60 min before each challenge. Structure (paper Methods equation):",
    "VAS(Cp) = BASE * (1 - Imax * Cp^gamma / (IC50^gamma + Cp^gamma)),",
    "with Imax and gamma both fixed to 1, so the relationship reduces to",
    "VAS = BASE * IC50 / (IC50 + Cp) and IC90 = 9 * IC50 = 0.675 ng/mL",
    "(reproducing the paper's independently reported IC90).",
    "Predictions are logit-transformed onto the range -1 to 52 with an",
    "additive residual error on that logit scale, encoded as",
    "logitNorm(addSd, -1, 52) (Supplementary Methods, RUV equations).",
    "The model has NO PK ODE and no compartments: it is a static algebraic",
    "exposure-response relationship, and CP_BPN_NGML is a required",
    "time-varying input covariate supplied per observation record. Users",
    "who need to drive it from a CAM2038 dose regimen must supply",
    "buprenorphine concentrations from an external PK source; the CAM2038",
    "popPK model is not reported in this paper and is not packaged in",
    "nlmixr2lib at extraction time.",
    sep = " "
  )
  reference <- paste(
    "Walsh SL, Comer SD, Aguiar Zdovc J, Sarr C, Bjornsson M,",
    "Strandgarden K, Hjelmstrom P, Tiberg F.",
    "Pharmacokinetic-pharmacodynamic analysis of drug liking blockade by",
    "buprenorphine subcutaneous depot (CAM2038) in participants with opioid",
    "use disorder.",
    "Neuropsychopharmacology. 2024;49(7):1050-1057.",
    "doi:10.1038/s41386-023-01793-z.",
    "Parameter estimates from Table 1; structural equation and the logit",
    "residual-error transformation from the Supplementary Appendix",
    "(41386_2023_1793_MOESM1_ESM.pdf).",
    "Underlying phase 2 study NCT02611752, reported in",
    "Walsh SL, Comer SD, Lofwall MR, et al. JAMA Psychiatry.",
    "2017;74(9):894-902. doi:10.1001/jamapsychiatry.2017.1874.",
    "Sister endpoint models from the same paper:",
    "modellib('Walsh_2024_buprenorphine_desireToUse').",
    sep = " "
  )
  vignette <- "Walsh_2024_buprenorphine_opioidBlockade"

  depends <- c("CP_BPN_NGML")

  units <- list(
    time          = "day",
    dosing        = "N/A (PD-only; buprenorphine plasma concentration is a required input covariate)",
    concentration = "VAS units (observation: period-corrected drug liking Emax bipolar VAS score); ng/mL (CP_BPN_NGML input covariate)"
  )

  covariateData <- list(
    CP_BPN_NGML = list(
      description        = "Time-matched buprenorphine plasma concentration driving the Imax inhibition of the drug liking Emax VAS score. Time-varying; supplied per observation record in the event table rather than computed from a coupled PK model.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Paper Methods: 'Plasma samples to quantify BPN concentration were collected approximately 60 min before each hydromorphone administration'; the analysis used one BPN concentration per challenge day, time-matched to the drug liking Emax VAS score. Observed BPN plasma concentration range across the phase 2 study was 0.636-12.3 ng/mL (Supplementary Methods, 'Model application'). Member of the canonical CP_<DRUG>_<UNITS> plasma-PD-driver family (siblings CP_MORPH_NGML, CP_OXY_NGML, CP_FBX_NGML). The upstream CAM2038 popPK model is not reported in this paper and is not packaged in nlmixr2lib; concentrations must be supplied externally.",
      source_name        = "Cp (Walsh 2024 Methods, structural Imax equation)"
    )
  )

  # ---------------------------------------------------------------------
  # Covariates screened on IC50 but NOT retained in the final model. The
  # paper reports only the negative result ("No statistically significant
  # effect of age, sex, ethnicity, race, body weight, height, or BMI was
  # found for IC50 [p > 0.2]", Results / Final PK/PD models) with no point
  # estimates, so there is nothing to encode in model(). Documented here
  # to preserve the provenance of the covariate screen.
  # ---------------------------------------------------------------------
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age in years.",
      units       = "years",
      type        = "continuous",
      notes       = "Screened on IC50 by single forward step in PsN stepwise covariate modelling (Methods, 'Covariate model building'); not retained (p > 0.2). Cohort median 36.0 (18.0-54.0) years, mean 35.8 (SD 9.1) (Table S1)."
    ),
    WT = list(
      description = "Body weight.",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened on IC50 as a power model; not retained (p > 0.2). Cohort median 73.4 (53.0-110.0) kg, mean 75.9 (SD 14.0) (Table S1)."
    ),
    HT = list(
      description = "Height.",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened on IC50 as a power model; not retained (p > 0.2). Cohort median 176 (154-192) cm (Table S1)."
    ),
    BMI = list(
      description = "Body mass index.",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened on IC50 as a power model; not retained (p > 0.2). Cohort median 24.3 (16.8-34.0) kg/m^2 (Table S1)."
    ),
    SEXF = list(
      description = "Sex, 1 = female, 0 = male.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on IC50 as a fractional difference to the most common category (male); not retained (p > 0.2). 12/47 (26%) female (Table S1)."
    ),
    RACE_BLACK = list(
      description = "Race indicator, 1 = Black or African American.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on IC50 as a fractional difference to the most common category; not retained (p > 0.2). 24/47 (51%) Black or African American, 22/47 (47%) White, 1/47 (2%) Other (Table S1)."
    ),
    ETHNIC_HISPANIC = list(
      description = "Ethnicity indicator, 1 = Hispanic or Latino.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on IC50 as a fractional difference to the most common category (not Hispanic or Latino); not retained (p > 0.2). 1/47 (2%) Hispanic or Latino (Table S1)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 47L,
    n_studies      = 1L,
    n_observations = 231L,
    age_range      = "18-54 years (median 36.0; mean 35.8, SD 9.1)",
    age_median     = "36.0 years",
    weight_range   = "53.0-110.0 kg (median 73.4; mean 75.9, SD 14.0)",
    weight_median  = "73.4 kg",
    sex_female_pct = 25.5,
    race_ethnicity = c(Black = 51, White = 47, Other = 2),
    disease_state  = "Non-treatment-seeking adults with moderate to severe opioid use disorder (DSM-5), physically dependent on opioids.",
    dose_range     = "CAM2038 (subcutaneous long-acting buprenorphine) 24 mg (n = 22) or 32 mg (n = 25) once weekly on days 0 and 7. Intramuscular hydromorphone challenges of 0, 6 or 18 mg once daily over each of the day -3 to -1, 1-3, 4-6, 8-10 and 11-13 periods; only the 18 mg challenge sessions contributed to this model.",
    regions        = "USA (multisite: University of Kentucky; Columbia University / New York State Psychiatric Institute)",
    biomarkers     = "Period-corrected drug liking Emax bipolar VAS score (0-100, 50 = neutral), recorded 30 min before and at 15 timepoints up to 300 min after each hydromorphone challenge; the Emax (maximum post-challenge value) was period-corrected against the pre-challenge value.",
    notes          = "Baseline demographics from Supplementary Table S1; observation counts from Supplementary Table S2 (drug liking: 231 PK and 231 PD records, 109 + 122 across the 24 mg and 32 mg arms). Study NCT02611752. Table S1 reports arm sizes 22 / 25; the Fig. 1 caption's '32 mg: n = 24' refers to the PK-evaluable subset. Data from the hydromorphone 0 mg and 6 mg challenges were excluded from PK/PD modelling because they produced no or lower drug liking responses."
  )

  ini({
    # =====================================================================
    # Structural parameters. Walsh 2024 Table 1 ("Parameter estimates of
    # the final PK/PD models", drug liking Emax VAS score block).
    #
    # Paper structural equation (Methods, "Model development"):
    #
    #                            Imax * Cp^gamma
    #   VAS(Cp) = BASE * ( 1 - ---------------------- )
    #                          IC50^gamma + Cp^gamma
    #
    # with Imax and gamma both FIXED to 1 (Results, "Model development":
    # "Imax and gamma were fixed to 1, but all other parameters were
    # estimated"). The relationship therefore reduces to the simple
    # hyperbolic form VAS = BASE * IC50 / (IC50 + Cp).
    #
    # Back-check on the fixed Hill coefficient: for a sigmoid Imax model
    # IC90 = IC50 * 9^(1/gamma). With gamma = 1 and IC50 = 0.075 ng/mL,
    # IC90 = 9 * 0.075 = 0.675 ng/mL, which reproduces exactly the IC90 of
    # 0.675 ng/mL the paper reports independently in Table 1 and in the
    # Abstract / Results. This confirms gamma = 1 and Imax = 1 as encoded.
    # =====================================================================

    # Baseline period-corrected drug liking Emax VAS score before any
    # CAM2038 administration, i.e. the full unblocked response to an
    # intramuscular hydromorphone 18 mg challenge. Log-transformed here for
    # a positively-constrained parameter carrying the paper's exponential
    # (log-normal) IIV.
    lbase <- log(45.3)
    label("BASE: baseline period-corrected drug liking Emax VAS score before CAM2038 (VAS units)")
    # Walsh 2024 Table 1, drug liking block: Baseline = 45.3 VAS units (RSE 2.40%).

    # Buprenorphine plasma concentration producing half of the maximum
    # inhibition of drug liking.
    lic50 <- log(0.075)
    label("IC50: BPN plasma concentration giving half-maximal drug-liking inhibition (ng/mL)")
    # Walsh 2024 Table 1, drug liking block: IC50 = 0.075 ng/mL (RSE 30.4%).

    # Maximum fractional inhibition of the drug liking response. Fixed to 1
    # (complete inhibition) by the authors, flagged "(FIX)" in Table 1.
    imax <- fixed(1)
    label("Imax: maximum fractional inhibition of drug liking (unitless)")
    # Walsh 2024 Table 1, drug liking block: Imax = 1.00 (FIX).

    # Sigmoidicity (Hill) coefficient. Fixed to 1 by the authors; not given
    # a row in Table 1 but stated in Results, "Model development".
    hill <- fixed(1)
    label("gamma: sigmoidicity (Hill) coefficient of the Imax relationship (unitless)")
    # Walsh 2024 Results, "Model development": "Imax and gamma were fixed to 1".

    # =====================================================================
    # Inter-individual variability. Table 1 reports the IIV magnitudes in a
    # column headed "SD"; the table footnote confirms "The relative standard
    # error (RSE) for interindividual variability (IIV) is reported on the
    # approximate standard deviation (SD) scale". The Supplementary Methods
    # ("Parameter Estimation") state IIV was carried in exponential form,
    # P_i = TVP * exp(eta_i) with eta_i ~ N(0, omega^2). The tabulated
    # values are therefore omega (SD on the log scale) and are squared here
    # to give the variances rxode2 expects.
    #
    # IIV on Imax is reported as SD = 0 (FIX) in Table 1, so no eta is
    # placed on imax.
    # =====================================================================
    etalbase ~ 0.106^2
    # Walsh 2024 Table 1: IIV Baseline SD = 0.106 (RSE 19.8%, shrinkage 11.6%).
    etalic50 ~ 1.90^2
    # Walsh 2024 Table 1: IIV IC50 SD = 1.90 (RSE 10.3%, shrinkage 13.3%).

    # =====================================================================
    # Residual unexplained variability. Additive on the LOGIT-transformed
    # prediction scale. Supplementary Methods ("Residual Unexplained
    # Variability") give the transformation explicitly:
    #
    #   TI_yhat  = (yhat + 1) / 53
    #   PHI_yhat = ln( TI_yhat / (1 - TI_yhat) )
    #   y        = expit(PHI_yhat + eps_add) * 53 - 1
    #
    # i.e. a logit-normal residual bounded on (-1, 52) -- exactly
    # rxode2's logitNorm(sd, low = -1, hi = 52). The bounds were chosen so
    # "predictions could fit in the range -1 to 52"; 52 is just above the
    # largest attainable increase on the bipolar 0-100 VAS from a neutral
    # score of 50, and only 3/231 observations were below 0.
    # =====================================================================
    addSd <- 0.744
    label("Additive residual SD on the logit-transformed drug liking VAS scale (unitless)")
    # Walsh 2024 Table 1: Additive RUV (logit scale) = 0.744 (RSE 8.62%, shrinkage 9.39%).
  })

  model({
    # -------------------------------------------------------------------
    # Individual parameters. Exponential (log-normal) IIV per the
    # Supplementary Methods parameter-estimation section.
    # -------------------------------------------------------------------
    base <- exp(lbase + etalbase)
    ic50 <- exp(lic50 + etalic50)

    # -------------------------------------------------------------------
    # Direct (no-delay) Imax exposure-response. The paper explicitly checked
    # for and found no hysteresis ("No hysteresis in the relationship
    # between BPN plasma concentration and drug liking Emax VAS scores was
    # observed", Results), which is why a direct rather than effect-
    # compartment model is appropriate.
    # -------------------------------------------------------------------
    vas <- base * (1 - imax * CP_BPN_NGML^hill / (ic50^hill + CP_BPN_NGML^hill))

    # -------------------------------------------------------------------
    # Logit-normal residual error on the bounded range -1 to 52.
    # -------------------------------------------------------------------
    vas ~ logitNorm(addSd, -1, 52)
  })
}
