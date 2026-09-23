Sood_2026_lu177dotatate_dosimetry_adult <- function() {
  description <- paste(
    "Empirical exposure-dosimetry model for the somatostatin-receptor-",
    "targeted radiopharmaceutical [177Lu]Lu-DOTATATE (177Lu-DOTATATE,",
    "lutetium Lu 177 dotatate) in ADULTS with gastroenteropancreatic",
    "neuroendocrine tumors, fitted to 47 patients pooled from NETTER-1",
    "(n = 20) and ERASMUS (n = 27). It predicts the PER-CYCLE ABSORBED",
    "RADIATION DOSE in gray delivered to the two organs at risk -- the",
    "kidneys and the red bone marrow -- from the administered",
    "radioactivity per cycle and baseline creatinine clearance, as a",
    "product of power terms normalized to the 7.4 GBq label activity and",
    "the 99 mL/min cohort median (Sood 2026 Eqs. 1-4). THERE IS NO PK",
    "LAYER AND NO ODE: only one dosimetry value per subject was",
    "available, so each organ was fitted by its own nonlinear regression",
    "with a proportional residual, and the two organ regressions share no",
    "parameter and no random effect. Radioactivity, not model-predicted",
    "AUC, is the exposure metric -- popPK-predicted exposure was",
    "available only for the NETTER-1 subset and did not correlate with",
    "kidney dosimetry (R = -0.027, P = 0.91), so activity was used as the",
    "surrogate that also admits the ERASMUS patients. Renal function is",
    "the dominant covariate and acts in the protective direction for both",
    "organs (exponents -0.552 and -1.11): faster clearance means less",
    "absorbed dose. Cumulative absorbed dose over the approved 4-cycle",
    "course is four times the per-cycle prediction. This model drove the",
    "NETTER-P sample-size determination. Companion models from the same",
    "paper: Sood_2026_lu177dotatate_dosimetry_pooled.R (the adult +",
    "adolescent refit that adds a study effect on the kidney renal-",
    "function term), Sood_2026_lu177dotatate_adult.R and",
    "Sood_2026_lu177dotatate_adolescent.R (the plasma popPK models).",
    sep = " "
  )
  reference <- paste(
    "Sood M, Lachi Silva L, Ho YY, Blumenstein L, Cherfi A, Xu L,",
    "Khanshan F. [177Lu]Lu-DOTATATE population pharmacokinetics and",
    "dosimetry modeling for adolescent and adult patients with",
    "somatostatin receptor-positive gastroenteropancreatic neuroendocrine",
    "tumors. J Nucl Med. 2026;67(6):887-894.",
    "doi:10.2967/jnumed.125.270202",
    sep = " "
  )
  vignette <- "Sood_2026_lu177dotatate"
  units <- list(
    time = "n/a (static per-cycle absorbed-dose regression; the model has no time dimension and no ODE state)",
    dosing = "n/a (no dose events; the administered activity enters as the covariate DOSE_LU177DOTATATE_GBQ in GBq)",
    concentration = "n/a (the outputs absDoseKidney and absDoseBoneMarrow are absorbed radiation doses in Gy per treatment cycle, not concentrations)"
  )

  covariateData <- list(
    DOSE_LU177DOTATATE_GBQ = list(
      description = "Administered [177Lu]Lu-DOTATATE radioactivity for the treatment cycle. The approved adult activity is 7.4 GBq per cycle for 4 cycles 8 weeks apart.",
      units = "GBq",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters both organ regressions as the power term",
        "(DOSE_LU177DOTATATE_GBQ / 7.4)^exponent, with the 7.4 GBq label",
        "activity as the normalizing reference printed inside Eqs. 3 and",
        "4. Used as a SURROGATE FOR EXPOSURE rather than as a dose",
        "amount: popPK-predicted AUC from time zero to infinity was",
        "available only for the 20 NETTER-1 patients and, for the kidney,",
        "did not correlate with dosimetry at all (R = -0.027, P = 0.91),",
        "whereas activity is recorded for all 47 adults including the 27",
        "ERASMUS patients for whom no popPK exposure metric was ever",
        "derived. The paper's own simulations exercise the term over",
        "1-8 GBq per cycle (Fig. 3), so that is the range over which the",
        "exponent is supported.",
        sep = " "
      ),
      source_name = "activity (GBq)"
    ),
    CRCL = list(
      description = paste(
        "Baseline creatinine clearance, raw (NOT body-surface-area",
        "normalized) mL/min. Table 1 gives the adult dosimetry set",
        "(n = 47) as median 98.82 mL/min, range 46.97-189.77; the",
        "equations round that median to the 99 mL/min normalizing",
        "reference. The paper does not name the estimating equation for",
        "the adult cohort.",
        sep = " "
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters both organ regressions as the power term (CRCL / 99)^",
        "exponent. Both exponents are NEGATIVE, so lower renal function",
        "raises the absorbed dose -- the mechanism the supplement",
        "attributes to Svensson 2015, 'patients with inferior renal",
        "function were exposed to higher renal absorbed dose and",
        "developed hematological toxicity'. Baseline value, held constant",
        "per subject. The paper's simulations exercise the term over",
        "35-180 mL/min (Methods) and its dosing conclusions are stated",
        "only down to 55 mL/min ('median kidney and bone marrow radiation",
        "exposures would remain within safe limits if CrCL remained above",
        "55 mL/min'). Note the SIGN REVERSAL for the kidney in",
        "adolescents, which is the whole reason the companion pooled",
        "model carries a study effect on this exponent.",
        sep = " "
      ),
      source_name = "CrCL (mL/min)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 47,
    n_studies = 2,
    age_median = "56 years",
    age_range = "29-83 years",
    weight_median = "75 kg",
    weight_range = "48-145 kg",
    sex_female_pct = 48.9,
    disease_state = "Gastroenteropancreatic neuroendocrine tumors (47 of 47, 100%).",
    renal_function = "Creatinine clearance median 98.82 mL/min, range 46.97-189.77. Kidney mass median 339 g, range 201-575.",
    dose_range = "Intravenous [177Lu]Lu-DOTATATE 7.4 GBq per cycle, 4 cycles 8 weeks apart (cumulative 29.6 GBq).",
    regions = "NETTER-1 (NCT01578239, n = 20, 42.55%) and ERASMUS (n = 27, 57.45%)",
    notes = paste(
      "One organ absorbed-dose value per subject, so the fit is an",
      "ordinary nonlinear regression rather than a mixed-effects model:",
      "there is no between-subject random effect to estimate and the",
      "proportional residual carries all of the unexplained variability.",
      "Candidate covariates screened by inspecting correlation plots of",
      "log organ dosimetry against log covariate (age, weight, body mass",
      "index, body surface area, height, creatinine clearance, activity)",
      "and visually for sex; at a p-value criterion of 0.01 'No strong",
      "correlation ... was observed between kidney or bone marrow",
      "dosimetry and covariates such as age, weight, BMI, BSA, height or",
      "sex', leaving creatinine clearance and activity. Model selection",
      "used the Bayesian information criterion plus visual predictive",
      "checks (Fig. 2) and goodness-of-fit plots (Supplemental Fig. 4).",
      sep = " "
    )
  )

  ini({
    # ==================================================================
    # Kidney absorbed dose per cycle -- Table 3 ('Final Dosimetry Model
    # Parameters Estimated for Adult Population (n = 47)'), written out
    # as Eq. 3:
    #
    #   kidney dosimetry = 4.3 * (activity(GBq) / 7.4(GBq))^0.66
    #                          * (CrCL(mL/min) / 99(mL/min))^-0.552
    #
    # which is the generic Eq. 1 form
    #   A_pop * (COV1 / median COV1)^B_pop * (COV2 / median COV2)^C_pop
    # with COV1 = activity and COV2 = creatinine clearance.
    #
    # A_pop is log-transformed here because it is a strictly positive
    # multiplicative reference value (the typical per-cycle absorbed dose
    # at the 7.4 GBq label activity and the 99 mL/min median creatinine
    # clearance), so the 4.3 recovered by exp() is exactly the printed
    # estimate.
    # ==================================================================
    lrbase_absDoseKidney     <- log(4.3)   ; label("Typical kidney absorbed dose per cycle at 7.4 GBq and 99 mL/min creatinine clearance (Gy)")  # Table 3 'A pop, kidney baseline' = 4.3 (%RSE 8.37); printed inside Eq. 3 as the leading multiplier
    e_activity_absDoseKidney <- 0.66       ; label("Power exponent of administered activity on the kidney absorbed dose per cycle (unitless)")   # Table 3 'B pop, activity effect on kidney' = 0.66 (%RSE 23); Eq. 3 exponent on (activity / 7.4)
    e_crcl_absDoseKidney     <- -0.552     ; label("Power exponent of creatinine clearance on the kidney absorbed dose per cycle (unitless)")    # Table 3 'C pop, CrCL effect on kidney' = -0.552 (%RSE 37.3); Eq. 3 exponent on (CrCL / 99)

    # ==================================================================
    # Bone marrow absorbed dose per cycle -- Table 3, written out as
    # Eq. 4:
    #
    #   bone marrow dosimetry = 0.246 * (activity(GBq) / 7.4(GBq))^0.597
    #                                 * (CrCL(mL/min) / 99(mL/min))^-1.11
    #
    # i.e. the generic Eq. 2 form
    #   D_pop * (COV3 / median COV3)^E_pop * (COV4 / median COV4)^F_pop.
    # ==================================================================
    lrbase_absDoseBoneMarrow     <- log(0.246) ; label("Typical bone marrow absorbed dose per cycle at 7.4 GBq and 99 mL/min creatinine clearance (Gy)")  # Table 3 'D pop, bone marrow baseline' = 0.246 (%RSE 11); printed inside Eq. 4 as the leading multiplier
    e_activity_absDoseBoneMarrow <- 0.597      ; label("Power exponent of administered activity on the bone marrow absorbed dose per cycle (unitless)")   # Table 3 'E pop, activity effect on bone marrow' = 0.597 (%RSE 35.6); Eq. 4 exponent on (activity / 7.4)
    e_crcl_absDoseBoneMarrow     <- -1.11      ; label("Power exponent of creatinine clearance on the bone marrow absorbed dose per cycle (unitless)")    # Table 3 'F pop, CrCL effect on bone marrow' = -1.11 (%RSE 24.4); Eq. 4 exponent on (CrCL / 99)

    # ==================================================================
    # Residual error. The supplement is explicit that each organ was
    # fitted 'assuming proportional error'; Table 3's footnote names b1
    # and b2 as 'proportional residual error for kidney model' and 'for
    # bone marrow model'. There is no between-subject random effect
    # because each subject contributed a single absorbed-dose value.
    # ==================================================================
    propSd_absDoseKidney     <- 0.515  ; label("Proportional residual SD for the kidney absorbed dose (fraction)")       # Table 3 'b 1' = 0.515 (%RSE 12.8)
    propSd_absDoseBoneMarrow <- 0.675  ; label("Proportional residual SD for the bone marrow absorbed dose (fraction)")  # Table 3 'b 2' = 0.675 (%RSE 14.3)
  })

  model({
    # Both normalizing constants are printed inside Eqs. 3 and 4: the
    # 7.4 GBq approved per-cycle activity and the 99 mL/min rounding of
    # the Table 1 cohort median creatinine clearance (98.82 mL/min).
    activityRatio <- DOSE_LU177DOTATATE_GBQ / 7.4
    crclRatio     <- CRCL / 99

    # Eq. 3 -- per-cycle kidney absorbed dose (Gy). Multiply by the
    # number of identical cycles to obtain the cumulative absorbed dose
    # the 23 Gy / 29 Gy external-beam thresholds are compared against.
    absDoseKidney <- exp(lrbase_absDoseKidney) *
      activityRatio^e_activity_absDoseKidney *
      crclRatio^e_crcl_absDoseKidney

    # Eq. 4 -- per-cycle red bone marrow absorbed dose (Gy), compared
    # against the conservative 2 Gy threshold after 4 cycles.
    absDoseBoneMarrow <- exp(lrbase_absDoseBoneMarrow) *
      activityRatio^e_activity_absDoseBoneMarrow *
      crclRatio^e_crcl_absDoseBoneMarrow

    absDoseKidney ~ prop(propSd_absDoseKidney)
    absDoseBoneMarrow ~ prop(propSd_absDoseBoneMarrow)
  })
}
