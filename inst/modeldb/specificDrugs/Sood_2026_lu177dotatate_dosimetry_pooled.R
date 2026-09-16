Sood_2026_lu177dotatate_dosimetry_pooled <- function() {
  description <- paste(
    "Empirical exposure-dosimetry model for the somatostatin-receptor-",
    "targeted radiopharmaceutical [177Lu]Lu-DOTATATE (177Lu-DOTATATE,",
    "lutetium Lu 177 dotatate) fitted to the POOLED adult and adolescent",
    "cohort (n = 57: 47 adults from NETTER-1 and ERASMUS plus 10",
    "adolescents from NETTER-P). It predicts the PER-CYCLE ABSORBED",
    "RADIATION DOSE in gray to the kidneys and the red bone marrow from",
    "the administered radioactivity per cycle and baseline creatinine",
    "clearance, as a product of power terms normalized to the 7.4 GBq",
    "label activity and 99 mL/min (Sood 2026 Eqs. 5-6). THERE IS NO PK",
    "LAYER AND NO ODE: one dosimetry value per subject, so each organ is",
    "its own nonlinear regression with a proportional residual and the",
    "two share no parameter. The one structural difference from the",
    "adult-only companion model is a STUDY EFFECT ON THE KIDNEY",
    "RENAL-FUNCTION EXPONENT: the creatinine-clearance exponent is -0.55",
    "in the adult studies but -0.55 + 1.58 = +1.03 in NETTER-P, so the",
    "relationship REVERSES SIGN in adolescents and kidney absorbed dose",
    "rises rather than falls with renal function. The paper flags this as",
    "provisional ('this result should be interpreted with caution because",
    "of the small sample size of the adolescent population (n = 10)') and",
    "notes that the adolescent cohort contains no renal impairment. Bone",
    "marrow needed no such term: 'no improvement in parameter estimates",
    "was observed when additional covariates were added'. Cumulative",
    "absorbed dose over the approved 4-cycle course is four times the",
    "per-cycle prediction. Companion models from the same paper:",
    "Sood_2026_lu177dotatate_dosimetry_adult.R (the adult-only fit that",
    "drove the NETTER-P sample size), Sood_2026_lu177dotatate_adult.R and",
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
      description = "Administered [177Lu]Lu-DOTATATE radioactivity for the treatment cycle. Both cohorts received the same flat 7.4 GBq per cycle for 4 cycles 8 weeks apart.",
      units = "GBq",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters both organ regressions as the power term",
        "(DOSE_LU177DOTATATE_GBQ / 7.4)^exponent, with the 7.4 GBq label",
        "activity as the normalizing reference printed inside Eqs. 5 and",
        "6. Used as a surrogate for exposure because popPK-predicted AUC",
        "was derived only for the NETTER-1 subset of the adults.",
        sep = " "
      ),
      source_name = "activity (GBq)"
    ),
    CRCL = list(
      description = paste(
        "Baseline creatinine clearance, raw (NOT body-surface-area",
        "normalized) mL/min. Table 1 medians are 98.82 mL/min (adults,",
        "range 46.97-189.77) and 122.1 mL/min (adolescents, range",
        "86-160); the equations use 99 mL/min as the single normalizing",
        "reference for the pooled fit. The adolescent scale is anchored",
        "to Piepsz 2008, whose subject is escaping the body-surface-area",
        "correction, so the column is absolute mL/min in both cohorts.",
        sep = " "
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters the bone marrow regression with a single negative",
        "exponent (-1.18) and the kidney regression with a",
        "STUDY-DEPENDENT exponent: -0.55 in the adult studies and",
        "-0.55 + 1.58 = +1.03 in NETTER-P. Baseline value, held constant",
        "per subject. Both cohorts skew to preserved or supranormal renal",
        "function (NETTER-P required at least 70 mL/min at entry), and",
        "the Discussion states the resulting limitation directly:",
        "'Interpretation is limited by small sample size and absence of",
        "renal impairment.'",
        sep = " "
      ),
      source_name = "CrCL (mL/min)"
    ),
    STUDY_NETTER_P = list(
      description = "Indicator for the NETTER-P adolescent study cohort: 1 for the 10 NETTER-P adolescents, 0 for the 47 NETTER-1 and ERASMUS adults.",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = adult studies (NETTER-1 and ERASMUS pooled)",
      notes = paste(
        "Printed inside Eq. 5 as the conditional exponent modifier",
        "'(-0.55 [+1.58 if study is NETTER-P])', so the column is a",
        "STUDY-LEVEL, not a subject-level, attribute and it modifies the",
        "creatinine-clearance EXPONENT rather than the kidney absorbed",
        "dose directly. The paper is explicit that the column stands in",
        "for age: 'The study effect, which can be assumed as adult versus",
        "adolescent population, was used instead of age.' It was added",
        "only for the kidney; the bone marrow model gained nothing from",
        "additional covariates. Two cautions the paper itself raises:",
        "the term rests on 10 adolescents, and it makes the",
        "renal-function relationship change sign, which is",
        "physiologically surprising and is flagged as needing cautious",
        "interpretation.",
        sep = " "
      ),
      source_name = "study (NETTER-P)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 57,
    n_studies = 3,
    age_median = "56 years in the 47 adults; 15 years in the 10 adolescents",
    age_range = "29-83 years in the adults; 13-17 years in the adolescents",
    weight_median = "75 kg in the adults; 55 kg in the adolescents",
    weight_range = "48-145 kg in the adults; 39.5-71 kg in the adolescents",
    disease_state = "Gastroenteropancreatic neuroendocrine tumors in all 47 adults; in the adolescent cohort, gastroenteropancreatic neuroendocrine tumors (4 of 11) or pheochromocytomas and paragangliomas (7 of 11).",
    renal_function = "Creatinine clearance median 98.82 mL/min (range 46.97-189.77) in the adults and 122.1 mL/min (range 86-160) in the adolescents; no renally impaired adolescent was enrolled.",
    dose_range = "Intravenous [177Lu]Lu-DOTATATE 7.4 GBq per cycle, 4 cycles 8 weeks apart (cumulative 29.6 GBq), flat-dosed in both cohorts.",
    regions = "NETTER-1 (NCT01578239, n = 20), ERASMUS (n = 27) and NETTER-P (NCT04711135, n = 10)",
    notes = paste(
      "The pooled dosimetry set carries 10 of the 11 NETTER-P",
      "adolescents: 'One patient with PPGL was excluded from the exposure",
      "dosimetry analysis, due to issues with whole-body planar and",
      "single-photon emission computed tomography (SPECT)/computed",
      "tomography (CT) images, wherein these images could not be utilized",
      "for dosimetry analysis.' Adolescent organ doses were computed from",
      "whole-body conjugate planar and abdominal SPECT/CT images acquired",
      "about 1-2, 18-26, 36-48 and 156-168 h after injection of cycle 1,",
      "reduced with the RADAR/MIRD method in OLINDA by CDE Dosimetry",
      "Services. One absorbed-dose value per subject, so the fit is an",
      "ordinary nonlinear regression with no between-subject random",
      "effect.",
      sep = " "
    )
  )

  ini({
    # ==================================================================
    # Kidney absorbed dose per cycle -- Table 5 ('Final Dosimetry Model
    # Parameters Estimated for Pooled Adults and Adolescents (n = 57)'),
    # written out as Eq. 5:
    #
    #   kidney dosimetry = 4.37 * (activity(GBq) / 7.4)^0.65
    #                           * (CrCL / 99)^(-0.55 [+1.58 if study is NETTER-P])
    #
    # A_pop is log-transformed here because it is a strictly positive
    # multiplicative reference value -- the typical per-cycle absorbed
    # dose at the 7.4 GBq label activity and 99 mL/min -- so exp()
    # recovers exactly the printed 4.37. The paper confirms that reading
    # in its own words: 'For kidney dosimetry, the predicted cumulative
    # absorbed dose based on the final model for adults and adolescents
    # for cycle 1 was calculated as 4.37 and 5.39 Gy, respectively.'
    # 4.37 is the adult reference patient exactly, and the adolescent
    # 5.39 is recovered by putting the NETTER-P exponent +1.03 on the
    # adolescent median creatinine clearance.
    # ==================================================================
    lrbase_absDoseKidney     <- log(4.37) ; label("Typical kidney absorbed dose per cycle at 7.4 GBq and 99 mL/min creatinine clearance in the adult studies (Gy)")  # Table 5 'A pop, kidney baseline' = 4.37 (%RSE 7.22); printed inside Eq. 5 as the leading multiplier
    e_activity_absDoseKidney <- 0.65      ; label("Power exponent of administered activity on the kidney absorbed dose per cycle (unitless)")                        # Table 5 'B pop, activity effect on kidney' = 0.65 (%RSE 22.0); Eq. 5 exponent on (activity / 7.4)
    e_crcl_absDoseKidney     <- -0.55     ; label("Power exponent of creatinine clearance on the kidney absorbed dose per cycle in the adult studies (unitless)")    # Table 5 'C pop, CrCL effect on kidney' = -0.55 (%RSE 35.1); Eq. 5 exponent on (CrCL / 99) for the reference (adult) studies

    # Study effect: an ADDITIVE shift of the creatinine-clearance
    # exponent, not a multiplier on the absorbed dose. It takes the
    # kidney exponent from -0.55 in the adult studies to +1.03 in
    # NETTER-P, which is the sign reversal the Discussion describes:
    # 'For kidney, "study" was added as a covariate to address differing
    # CrCL effects, which showed an opposite trend in adolescents
    # (dosimetry increased with CrCL).'
    e_study_netter_p_crcl_absDoseKidney <- 1.58 ; label("Additive shift of the creatinine-clearance exponent on kidney absorbed dose for the NETTER-P adolescent cohort versus the adult studies (unitless)")  # Table 5 'beta_C_STUDY_NETTER_P, CrCL effect on kidney dosimetry based on adult or adolescent populations' = 1.58 (%RSE 37.5); printed inside Eq. 5 as '[+1.58 if study is NETTER-P]'

    # ==================================================================
    # Bone marrow absorbed dose per cycle -- Table 5, written out as
    # Eq. 6:
    #
    #   bone marrow dosimetry = 0.24 * (activity(GBq) / 7.4)^0.52
    #                                * (CrCL / 99)^-1.18
    #
    # No study effect: 'For bone marrow dosimetry, no improvement in
    # parameter estimates was observed when additional covariates were
    # added.'
    # ==================================================================
    lrbase_absDoseBoneMarrow     <- log(0.24) ; label("Typical bone marrow absorbed dose per cycle at 7.4 GBq and 99 mL/min creatinine clearance (Gy)")  # Table 5 'D pop, bone marrow baseline' = 0.24 (%RSE 8.90); printed inside Eq. 6 as the leading multiplier
    e_activity_absDoseBoneMarrow <- 0.52      ; label("Power exponent of administered activity on the bone marrow absorbed dose per cycle (unitless)")   # Table 5 'E pop, activity effect on bone marrow' = 0.52 (%RSE 37.8); Eq. 6 exponent on (activity / 7.4)
    e_crcl_absDoseBoneMarrow     <- -1.18     ; label("Power exponent of creatinine clearance on the bone marrow absorbed dose per cycle (unitless)")    # Table 5 'F pop, CrCL effect on bone marrow' = -1.18 (%RSE 19.7); Eq. 6 exponent on (CrCL / 99)

    # ==================================================================
    # Residual error. Each organ was fitted by nonlinear regression
    # 'assuming proportional error' (supplement); Table 5's footnote
    # names b1 and b2 as the proportional error for the kidney and bone
    # marrow models. There is no between-subject random effect because
    # each subject contributed a single absorbed-dose value.
    # ==================================================================
    propSd_absDoseKidney     <- 0.48  ; label("Proportional residual SD for the kidney absorbed dose (fraction)")       # Table 5 'b 1' = 0.48 (%RSE 11.3)
    propSd_absDoseBoneMarrow <- 0.63  ; label("Proportional residual SD for the bone marrow absorbed dose (fraction)")  # Table 5 'b 2' = 0.63 (%RSE 12.5)
  })

  model({
    # Normalizing constants printed inside Eqs. 5 and 6: the 7.4 GBq
    # approved per-cycle activity and 99 mL/min.
    activityRatio <- DOSE_LU177DOTATATE_GBQ / 7.4
    crclRatio     <- CRCL / 99

    # Eq. 5 -- the creatinine-clearance exponent is shifted for the
    # NETTER-P adolescent cohort, reversing its sign.
    crclExpKidney <- e_crcl_absDoseKidney +
      e_study_netter_p_crcl_absDoseKidney * STUDY_NETTER_P

    # Eq. 5 -- per-cycle kidney absorbed dose (Gy). Multiply by the
    # number of identical cycles to obtain the cumulative absorbed dose
    # the 23 Gy / 29 Gy external-beam thresholds are compared against.
    absDoseKidney <- exp(lrbase_absDoseKidney) *
      activityRatio^e_activity_absDoseKidney *
      crclRatio^crclExpKidney

    # Eq. 6 -- per-cycle red bone marrow absorbed dose (Gy), compared
    # against the conservative 2 Gy threshold after 4 cycles.
    absDoseBoneMarrow <- exp(lrbase_absDoseBoneMarrow) *
      activityRatio^e_activity_absDoseBoneMarrow *
      crclRatio^e_crcl_absDoseBoneMarrow

    absDoseKidney ~ prop(propSd_absDoseKidney)
    absDoseBoneMarrow ~ prop(propSd_absDoseBoneMarrow)
  })
}
