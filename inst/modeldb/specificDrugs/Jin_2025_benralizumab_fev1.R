Jin_2025_benralizumab_fev1 <- function() {
  description <- "Longitudinal pre-bronchodilator FEV1 model for benralizumab in patients with severe eosinophilic asthma (Jin 2025), with an additive baseline, an exponential-onset placebo effect and a time-driven (exposure-independent) Emax treatment effect; sex and theophylline co-medication on baseline FEV1 and baseline blood eosinophil count on the treatment Emax"
  reference <- paste(
    "Jin Y, Guiastrennec B, Stuke M, Yao Y, Zhang Y, Barker P, Jison M,",
    "Penland RC, Ding J, Lukka PB.",
    "Population pharmacokinetics and exposure-response analysis of benralizumab in",
    "Chinese adults, adolescents, and pediatric participants with severe eosinophilic asthma.",
    "Clin Pharmacokinet. 2025;64:1233-1245. doi:10.1007/s40262-025-01538-9.",
    "Model equation and parameter values from Resource 4 of the electronic supplementary material",
    "('Global legacy model (full data) - base/final model' column).",
    "Updates the legacy exposure-response model of Chia YL, Yan L, Yu B, et al.",
    "Clin Pharmacol Ther. 2019;106:383-90; doi:10.1002/cpt.1371.",
    sep = " "
  )
  vignette <- "Jin_2025_benralizumab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Log-scale multiplicative effect on baseline FEV1: FEV1baseline * exp(-0.09 * SEXF), i.e. an",
        "8.6% lower baseline pre-bronchodilator FEV1 in women (Jin 2025 Resource 4,",
        "'Beta_FEVB (SEXF_1)' = -0.09, RSE 16%; the Resource 4 footnote defines",
        "'FEVB(SEXF_1) female on FEVB'). The source column SEXF is already 1 for female",
        "(Jin 2025 Resource 1: 'SEXF Gender (0: male, 1: female)'), so it matches the canonical",
        "orientation and needs no value inversion. Time-fixed per subject."
      ),
      source_name        = "SEXF"
    ),
    CONMED_THEOPHYLLINE = list(
      description        = "Concomitant theophylline / aminophylline use, 1 = using, 0 = not using",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no theophylline or aminophylline co-medication)",
      notes              = paste(
        "Log-scale multiplicative effect on baseline FEV1: FEV1baseline * exp(-0.089 *",
        "CONMED_THEOPHYLLINE), i.e. an 8.5% lower baseline pre-bronchodilator FEV1 in users",
        "(Jin 2025 Resource 4, 'Beta_FEVB (CTHEO_1)' = -0.089, RSE 18%; the Resource 4 footnote",
        "defines 'FEVB(CTHEO_1) theophylline/aminophylline use on FEVB'). The effect is a severity",
        "marker rather than a drug interaction -- theophylline is a late-line asthma controller, so",
        "its users have worse baseline lung function. Jin 2025 Resource 1 defines the column as",
        "'CTHEO Theophylline/aminophylline use (0: no, 1: yes)'. Time-fixed per subject.",
        "Source column CTHEO."
      ),
      source_name        = "CTHEO"
    ),
    EOS = list(
      description        = "Baseline blood eosinophil count",
      units              = "cells/uL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power effect on the benralizumab treatment Emax: Emax * (EOS/380)^0.699",
        "(Jin 2025 Resource 4, 'Beta_EMAX (BEOSL)' = 0.699, RSE 15%; the Resource 4 footnote defines",
        "'EMAX(BEOSL) continuous baseline eosinophil count in cell/uL on EMAX (centered around",
        "380 cell/uL)', which supplies both the covariate units and the 380 cells/uL reference).",
        "The power (rather than log-linear) form follows the continuous-covariate equation printed in",
        "Jin 2025 Resource 2, p_i = p_pop * (COV_i/REF)^beta. Higher baseline eosinophil counts give a",
        "larger FEV1 response, consistent with the legacy analysis (Chia 2019: 'Patients with greater",
        "baseline eosinophil counts were associated with superior FEV1 response'). Baseline",
        "(pre-first-dose) value, time-fixed per subject. Source column BEOSL."
      ),
      source_name        = "BEOSL"
    ),
    TRT_BENRALIZUMAB = list(
      description        = "Benralizumab treatment-arm indicator, 1 = randomised to benralizumab, 0 = randomised to placebo",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (placebo)",
      notes              = paste(
        "Gates the treatment-effect term on and off. Unlike a concentration-driven exposure-response",
        "model, the Jin 2025 FEV1 treatment effect is a function of TIME ALONE",
        "(Emax * t / (T50 + t), Resource 4) and so does not vanish of its own accord in the placebo",
        "arm; an explicit arm indicator is therefore required to reproduce the placebo profile. This",
        "matches how the legacy model was fitted -- Chia 2019 Methods: 'When modeling the drug effect",
        "(Emax), only FEV1 data from benralizumab-treated patients were analyzed.' The placebo arm is",
        "described by the baseline plus the exponential-onset placebo term only. Time-fixed per",
        "subject. Derived from the randomised treatment assignment (source column TRTAN, coded",
        "0 = placebo, 1 = 30 mg Q4W, 2 = 30 mg Q8W in Jin 2025 Resource 3: set",
        "TRT_BENRALIZUMAB = 1 when TRTAN > 0)."
      ),
      source_name        = "TRTAN"
    )
  )

  # Covariates the source model retains but whose functional form is NOT
  # recoverable from any document on disk, plus covariates that were screened
  # and dropped. Deliberately absent from model(); see the vignette Errata.
  covariatesDataExcluded <- list(
    HT = list(
      description = "Baseline standing height",
      units       = "cm",
      type        = "continuous",
      notes       = paste(
        "RETAINED IN THE SOURCE MODEL BUT NOT ENCODABLE. Jin 2025 Resource 4 states that 'the model",
        "included covariate effects of age and height on Pmax, age, height, sex, theophylline",
        "comedication on FEV1baseline', and tabulates a single shared coefficient",
        "'RHGHT' = 2.4 (RSE 5.1%), which its footnote defines only as 'effect of height on PMAX and",
        "FEV1'. Neither the functional form (linear in cm? power? normalised to a reference?) nor a",
        "reference height is stated, and height is not reported in ANY demographics table in Jin 2025",
        "or its supplement, so no median is available to centre on. The upstream legacy publication",
        "(Chia 2019, acquired and checked) prints the same base equation and names height as",
        "significant on both Pmax and baseline FEV1, but likewise gives no functional form, no",
        "reference value and no parameter table; its Figures S2 and S3 are individual-prediction",
        "diagnostic scatter plots with LOESS trends by age group, which cannot identify a covariate",
        "form because they carry the (large) inter-individual variability. Guessing a form and a",
        "reference height would move baseline FEV1 by hundreds of mL, so the effect is omitted rather",
        "than invented. The tabulated FEVB (1710 mL) and PmaxB (172 mL) are consistent with being the",
        "typical values AT the reference covariates, so omitting the centred term leaves a coherent",
        "typical-subject model."
      )
    ),
    AGE = list(
      description = "Baseline subject age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "RETAINED IN THE SOURCE MODEL BUT NOT ENCODABLE, for the same reason as HT. Jin 2025",
        "Resource 4 tabulates two age coefficients -- 'RPAGE' = -5.92 ('effect of age on PMAX') and",
        "'RFAGE' = -14.6 ('effect of age on FEV1') -- with no functional form and no reference age.",
        "The upstream Chia 2019 text further indicates the age effect on the placebo response is NOT",
        "a simple linear slope: 'adolescents had the strongest placebo effect compared with other age",
        "groups (with the same standing height). No placebo effect was observed for patients aged",
        "> 65 years', and its Figure S2 stratifies the LOESS trends into four age bands",
        "(12-<18, 18-<50, 50-<65, >=65-75 years). Whether the published coefficients are linear",
        "slopes in years, band contrasts, or transformed percentages cannot be determined from any",
        "document on disk, so the age effects are omitted rather than invented."
      )
    ),
    FEV1_BL = list(
      description = "Observed baseline pre-bronchodilator FEV1 supplied as a covariate column",
      units       = "mL",
      type        = "continuous",
      notes       = paste(
        "Screened on the placebo effect and on Emax in the legacy analysis (Chia 2019 Methods) and",
        "not retained. Baseline FEV1 is a MODEL PARAMETER here (rbase with IIV etalrbase), not a",
        "covariate input."
      )
    ),
    WT = list(
      description = "Baseline body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened on Emax in the legacy analysis (Chia 2019 Methods) and not retained."
    ),
    BMI = list(
      description = "Baseline body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened on the placebo effect and on Emax in the legacy analysis (Chia 2019 Methods) and not retained."
    ),
    CRCL = list(
      description = "Baseline creatinine clearance",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Screened on the placebo effect in the legacy analysis (Chia 2019 Methods) and not retained."
    ),
    RACE_CHINESE = list(
      description = "Chinese-heritage race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Tested in Jin 2025 and NOT retained -- this is the paper's headline negative result for the",
        "FEV1 endpoint. Results 3.4.1: 'The global FEV1 model was predictive of the MIRACLE data, and",
        "no statistically significant effect on participants from China was found. Consequently, no",
        "additional simulations were conducted.' Conclusions: 'There was no statistically significant",
        "effect of Chinese race upon FEV1.' The base model was therefore retained as the final model,",
        "which is why this file has no Chinese term while its sibling",
        "modellib('Jin_2025_benralizumab_aaer') does."
      )
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Tested in Jin 2025 and not retained; Resource 6 reports that 'no statistically significant",
        "effect on participants from Asia was found in the FEV1 model'. Resource 5 shows the",
        "corresponding visual predictive checks stratified by both Chinese and Asian race."
      )
    ),
    REGION_EUROPE = list(
      description = "Geographic enrollment region",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Screened on the placebo effect and on Emax in the legacy analysis (Chia 2019 Methods) and not retained on either."
    ),
    CONMED_STEROID = list(
      description = "Maintenance oral corticosteroid use at baseline",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened on Emax in the legacy analysis (Chia 2019 Methods, 'OCS use') and not retained on",
        "the FEV1 endpoint. It IS retained on the exacerbation-rate baseline -- see",
        "modellib('Jin_2025_benralizumab_aaer')."
      )
    ),
    NEXAC12M = list(
      description = "Number of asthma exacerbations in the 12 months before study entry",
      units       = "count (events in the prior 12 months)",
      type        = "count",
      notes       = paste(
        "Screened on Emax in the legacy analysis as a '>= 3 exacerbations' indicator (Chia 2019",
        "Methods) and not retained on the FEV1 endpoint. It IS retained on the exacerbation-rate",
        "baseline -- see modellib('Jin_2025_benralizumab_aaer')."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 3205L,
    n_studies      = 3L,
    age_range      = "12-75 years",
    age_median     = "not published for the exposure-response subset",
    weight_range   = "not published for the exposure-response subset",
    weight_median  = "not published for the exposure-response subset",
    sex_female_pct = NA_real_,
    race_ethnicity = "Chinese, Asian and non-Asian strata were tested on this endpoint and none was significant; the exposure-response dataset is the SIROCCO + CALIMA + MIRACLE pool.",
    disease_state  = "Severe, uncontrolled eosinophilic asthma on medium-to-high-dose inhaled corticosteroid plus a long-acting beta2-agonist, with benralizumab or placebo as add-on maintenance therapy.",
    dose_range     = "Placebo, benralizumab 30 mg subcutaneously every 4 weeks, or benralizumab 30 mg subcutaneously every 8 weeks (first three doses every 4 weeks).",
    regions        = "Multi-regional (SIROCCO, CALIMA and MIRACLE).",
    notes          = paste(
      "Jin 2025 Methods 2.5: 'Data from three phase III studies (SIROCCO, CALIMA, and MIRACLE) were",
      "used in the ER analysis'. THE PAPER DOES NOT PRINT AN ANALYSIS-SET N for the FEV1 dataset, so",
      "n_subjects is the combined RANDOMISED total of those three studies (SIROCCO 1204 +",
      "CALIMA 1306 + MIRACLE 695 = 3205) and is an UPPER BOUND rather than an exact analysis-set",
      "size. Note additionally that only benralizumab-treated patients contributed to estimating the",
      "treatment Emax (Chia 2019 Methods), while all arms contributed to the baseline and placebo",
      "terms, so no single N describes the whole fit. See the vignette Errata.",
      "The endpoint is PRE-BRONCHODILATOR FEV1 in mL, and every parameter of this model is in mL, as",
      "Jin 2025 Resource 4 reports them; no unit conversion has been applied.",
      "This model is EXPOSURE-INDEPENDENT. Because the exposure-response relationship for FEV1 was",
      "flat over the studied dose range, the treatment effect was modelled as a function of TIME",
      "rather than of concentration, so the model carries no PK layer and no concentration term and",
      "is valid only for the 30 mg subcutaneous regimens studied. Chia 2019 (the legacy analysis)",
      "gives the reason: the estimated FEV1 EC50 fell below the assay LLOQ, 'therefore, accurate",
      "estimation of EC50 is not possible ... Given the flat exposure-response relationship, the",
      "prebronchodilator FEV1 data were modeled longitudinally.'",
      "The placebo effect is large relative to the drug effect (Pmax 172 mL versus Emax 104 mL) and",
      "carries very large inter-individual variability (SD 292 mL, additive).",
      "Onset timing is only weakly identified: Chia 2019 notes that estimation of the onset of both",
      "the placebo and the treatment effect 'is limited by lack of measurements prior to week 4'.",
      "TWO RETAINED COVARIATE EFFECTS ARE OMITTED because their functional forms are unreported in",
      "both Jin 2025 and the upstream Chia 2019 -- see covariatesDataExcluded$HT and",
      "covariatesDataExcluded$AGE and the vignette Errata."
    )
  )

  ini({
    # ------------------------------------------------------------------------
    # Longitudinal pre-bronchodilator FEV1 model. Jin 2025 Resource 4 prints
    # the structure as
    #
    #   FEV1(t) = FEV1baseline
    #             + Pmax * (1 - exp(-kpbo * t))
    #             + Emax * t / (T50 + t)
    #
    # where Pmax is the maximum placebo effect as t -> infinity, kpbo the
    # placebo onset rate, Emax the maximum benralizumab effect and T50 the time
    # to half of it. All values are the "Global legacy model (full data) -
    # base/final model" column of Resource 4, which Jin 2025 Results 3.4.1
    # confirms is also the FINAL model ("No statistically significant covariate
    # effect was found, and the base model was retained as the final model").
    #
    # EVERY parameter here is in mL or 1/day, exactly as Resource 4 reports.
    # ------------------------------------------------------------------------

    # Baseline FEV1. Log-transformed because the two retained covariate effects
    # on it (sex, theophylline) are log-scale multipliers and its IIV is
    # lognormal.
    lrbase <- log(1710);  label("Typical baseline pre-bronchodilator FEV1 in a male non-theophylline user (mL)")  # Jin 2025 Resource 4, "FEVB" = 1710 (RSE 1.1%)

    # Placebo effect. Pmax carries an ADDITIVE (mL-scale) random effect: its
    # typical value row is "PmaxB" = 172 mL while a separate row, "PmaxIIV",
    # has a typical value of 0 FIXED and an omega of 292 -- i.e. the authors
    # carry the placebo-effect IIV on a mean-zero additive placeholder
    # parameter. Pmax is therefore kept on the natural scale, not logged.
    plbmax <- 172;         label("Maximum placebo effect on FEV1 (mL)")                              # Jin 2025 Resource 4, "PmaxB" = 172 (RSE 5.3%)
    lkplb  <- log(0.0356); label("First-order rate constant of placebo-effect onset (1/day)")        # Jin 2025 Resource 4, "KPBL" = 0.0356 (RSE 6%)

    # Benralizumab treatment effect. Time-driven, not concentration-driven.
    lemax <- log(104);   label("Maximum benralizumab treatment effect on FEV1 (mL)")                 # Jin 2025 Resource 4, "Emax" = 104 (RSE 11%)
    lt50  <- log(9.21);  label("Time to half of the maximum benralizumab FEV1 effect T50 (day)")     # Jin 2025 Resource 4, "T50" = 9.21 (RSE 14%)

    # Covariate effects.
    e_sexf_rbase <- -0.09;   label("Log-scale effect of female sex on baseline FEV1 (unitless; -8.6%)")           # Jin 2025 Resource 4, "Beta_FEVB (SEXF_1)" = -0.09 (RSE 16%)
    e_theo_rbase <- -0.089;  label("Log-scale effect of theophylline co-medication on baseline FEV1 (unitless; -8.5%)")  # Jin 2025 Resource 4, "Beta_FEVB (CTHEO_1)" = -0.089 (RSE 18%)
    e_eos_emax   <-  0.699;  label("Power exponent of EOS/380 on the benralizumab Emax (unitless)")               # Jin 2025 Resource 4, "Beta_EMAX (BEOSL)" = 0.699 (RSE 15%)

    # ------------------------------------------------------------------------
    # Inter-individual variability. Resource 4 carries no explicit SD-versus-
    # variance footnote, unlike Resources 8, 10 and 11 which all state that
    # "Omega values in this table are presented as standard deviations". Two
    # independent checks confirm the SD reading here as well:
    #   (1) Omega(PMAXIIV) = 292 is an ADDITIVE effect in mL. Read as an SD it
    #       gives a +/-2 SD span of about +/-580 mL, which matches the spread of
    #       individual Pmax values in Chia 2019 Figure S2 (roughly -900 to
    #       +1400 mL). Read as a variance it would give an SD of 17 mL, which
    #       that figure flatly contradicts.
    #   (2) Omega(KPBL) = 1.32 read as a lognormal SD is 132% CV, essentially
    #       identical to the 128% CV that the sibling asthma FEV1 model
    #       modellib('Zhang_2025_dupilumab_fev1') reports for the same placebo
    #       onset-rate parameter.
    # nlmixr2's `~` takes a VARIANCE, so each SD is squared below.
    # ------------------------------------------------------------------------
    etalrbase ~ 0.071289  # Jin 2025 Resource 4: Omega(FEVB) SD 0.267 (RSE 1.6%) -> 0.267^2; lognormal on baseline FEV1
    etaplbmax ~ 85264     # Jin 2025 Resource 4: Omega(PMAXIIV) SD 292 mL (RSE 1.9%) -> 292^2; ADDITIVE eta, Pmax = PmaxB + eta
    etalkplb  ~ 1.7424    # Jin 2025 Resource 4: Omega(KPBL) SD 1.32 (RSE 6.5%) -> 1.32^2; lognormal
    etalemax  ~ 1.3456    # Jin 2025 Resource 4: Omega(EMAX) SD 1.16 (RSE 6.7%) -> 1.16^2; lognormal

    # ------------------------------------------------------------------------
    # Residual error: combined additive (mL) plus proportional (fraction) on
    # predicted FEV1, per the Resource 4 footnote ("ADD1 additive error (mL) -
    # predicted FEV1 (mL)", "PROP1 proportional error (fraction) - predicted
    # FEV1 (mL)").
    # ------------------------------------------------------------------------
    addSd_FEV1  <- 137;     label("Additive residual error on FEV1 (mL)")             # Jin 2025 Resource 4, "Error_ADD1" = 137 (RSE 0.98%)
    propSd_FEV1 <- 0.0883;  label("Proportional residual error on FEV1 (fraction)")   # Jin 2025 Resource 4, "Error_PROP1" = 0.0883 (RSE 0.75%)
  })
  model({
    # --- Baseline FEV1 ------------------------------------------------------
    # Reference subject: male, not using theophylline / aminophylline.
    rbase <- exp(lrbase + etalrbase) *
             exp(e_sexf_rbase * SEXF) *
             exp(e_theo_rbase * CONMED_THEOPHYLLINE)

    # --- Placebo effect -----------------------------------------------------
    # Pmax uses an additive (mL) random effect; see the ini() note.
    plbmaxi <- plbmax + etaplbmax
    kplb    <- exp(lkplb + etalkplb)
    plbeff  <- plbmaxi * (1 - exp(-kplb * t))

    # --- Benralizumab treatment effect --------------------------------------
    # Driven by TIME, not by concentration, so it must be gated by the
    # treatment-arm indicator to reproduce the placebo profile.
    emaxi   <- exp(lemax + etalemax) * (EOS / 380)^e_eos_emax
    t50     <- exp(lt50)
    drugeff <- emaxi * t / (t50 + t) * TRT_BENRALIZUMAB

    # --- Observation --------------------------------------------------------
    FEV1 <- rbase + plbeff + drugeff
    FEV1 ~ add(addSd_FEV1) + prop(propSd_FEV1)
  })
}
