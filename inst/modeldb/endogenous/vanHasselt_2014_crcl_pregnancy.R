vanHasselt_2014_crcl_pregnancy <- function() {
  description <- paste(
    "Nonlinear mixed-effect model for the gestation-induced rise in maternal",
    "creatinine clearance (CrCL, mL/min) across pregnancy. CrCL follows an",
    "additive hyperbolic (Anderson-Holford, hill = 1) trajectory in maternal",
    "gestational age, rising from a nonpregnant baseline to a plateau of",
    "baseline + span, with half of the span reached at crcl_ega50 weeks.",
    "The three structural values were FIXED from an upstream meta-analysis of",
    "gestational CrCL dynamics; only the three between-subject variances and",
    "the proportional residual error were estimated here. The model has no",
    "drug compartments and no dosing: it predicts CrCL at each user-supplied",
    "gestational week. It is the covariate-generating layer of the sibling",
    "model vanHasselt_2014_cefazolin_semiphysiological, and is reusable as a",
    "gestational renal-function prior for other renally cleared drugs.",
    sep = " "
  )
  reference <- paste(
    "van Hasselt JGC, Allegaert K, van Calsteren K, Beijnen JH, Schellens JHM, Huitema ADR.",
    "Semiphysiological versus empirical modelling of the population pharmacokinetics of free and total cefazolin during pregnancy.",
    "Biomed Res Int. 2014;2014:897216. doi:10.1155/2014/897216.",
    "Corrigendum: Biomed Res Int. 2015;2015:124035. doi:10.1155/2015/124035.",
    "The structural CrCL trajectory values were fixed from the upstream meta-analysis",
    "cited as reference 13 of the source paper.",
    sep = " "
  )
  vignette <- "vanHasselt_2014_cefazolin"
  units <- list(
    time = "week (maternal estimated gestational age, EGA; 0 = prepregnancy)",
    dosing = "n/a (no exogenous dosing; endogenous renal-function model)",
    concentration = "mL/min (creatinine clearance, Cockcroft-Gault with body weight)"
  )

  population <- list(
    species = "human",
    n_subjects = 94L,
    n_studies = 3L,
    age_range = "20-42 years",
    age_median = "31 years",
    weight_range = "54-99 kg",
    weight_median = "72 kg",
    sex_female_pct = 100,
    race_ethnicity = "not reported in the source paper",
    disease_state = "pregnant women undergoing in utero surgical intervention, elective caesarean delivery, or fetal intervention",
    dose_range = "n/a (no exogenous drug)",
    ga_range = "17-40 weeks (median 33)",
    regions = "Belgium (University Hospitals Leuven) plus two previously published cohorts",
    renal_function = "serum creatinine median 0.64 mg/dL (range 0.33-0.88); creatinine clearance computed by Cockcroft-Gault using body weight",
    notes = paste(
      "TWO-LEVEL PROVENANCE. The three STRUCTURAL parameters were not estimated from",
      "the cefazolin cohort at all: Methods 2.4.1 states they 'were fixed to the",
      "previously estimated values (Table 3), which were based on a meta-analysis of",
      "the literature reporting gestational changes in CrCL', itself a sample-size",
      "weighted regression over two studies of CrCL dynamics during pregnancy",
      "(source-paper references 13-15). Only the three RANDOM effects and the",
      "proportional residual error were estimated, against the observed CrCL values",
      "of the 94-woman pooled cefazolin cohort described in Table 1.",
      "Serum creatinine was missing for 58% of patients (Methods 2.5); those values",
      "were imputed either by the pooled median or, in the semiphysiological arm, by",
      "the typical CrCL change predicted by this very model. The paper is explicit",
      "(Results 3.3.1) that Figure 3 'is merely intended as an illustration of how",
      "predictions of CrCL were generated and specifically not as a goodness-of-fit",
      "plot'."
    )
  )

  ini({
    # ---- Structural trajectory (van Hasselt 2014 Table 3) ----
    # Methods eq (4), verbatim from the JATS MathML of the source article:
    #     CrCL(t) = CrCL_0 + (CrCL_MAX * t) / (CrCL_50 + t)
    # with t = maternal gestational age in weeks, and (Methods 2.4.1)
    # 'CrCL0 represented baseline (nonpregnant CrCL), CrCLMAX represented
    # maximum typical increase in CrCL, and CrCL50 represented the time of
    # half-maximum change in CrCL'.
    #
    # This is the additive Anderson-Holford maturation shape with hill = 1,
    # but on a GESTATIONAL-AGE axis rather than the registered family's
    # postnatal-age axis, so every token is kept inside the crcl_ namespace
    # (operator ruling 2026-09-21; see inst/references/parameter-names.md).
    # crcl_matspan is the SPAN, not the plateau: the plateau CrCL reaches at
    # full term is crcl_ega0 + crcl_matspan * 40/(crcl_ega50 + 40).
    #
    # All three are FIXED. Table 3 marks each with an asterisk footnoted
    # '*These values were fixed during estimation of the mixed effect model;
    # that is, only random effects were estimated.' The RSE% printed beside
    # each is the precision carried over from the upstream meta-analysis, not
    # a precision estimated here.
    lcrcl_ega0 <- fixed(log(97.83))
    label("Baseline nonpregnant creatinine clearance at EGA 0 (mL/min)") # Table 3, row 'Baseline CrCL' CrCL0 = 97.83 mL/min (upstream RSE 3.91%), FIXED
    lcrcl_matspan <- fixed(log(83.83))
    label("Maximum gestational increase in creatinine clearance, the span above baseline (mL/min)") # Table 3, row 'Maximum CrCL' CrCLMAX = 83.83 (upstream RSE 12.48%), FIXED
    lcrcl_ega50 <- fixed(log(13.3))
    label("Gestational age at half of the creatinine-clearance span (weeks)") # Table 3, row 'Time of half-maximum CrCL' CrCL50 = 13.3 weeks (upstream RSE 37.59%), FIXED

    # ---- Between-subject variability (van Hasselt 2014 Table 3) ----
    # Methods eq (1) is exponential: P_i = P * exp(eta_i). Table 3 heads this
    # block 'Between subject variability (CV%)', so the internal variance is
    # taken as (CV/100)^2, the usual NONMEM CV% = 100*sqrt(omega^2) reading.
    # See the vignette Errata for the exact-lognormal alternative
    # log(1 + CV^2), which matters most for the 111.8 CV% term.
    # These three ARE estimated here (Results 3.3.1).
    etalcrcl_ega0 ~ 0.098596 # Table 3, row 'Baseline CrCL' omega_CrCL0 = 31.4 CV% (RSE 16%) -> 0.314^2
    etalcrcl_matspan ~ 0.123201 # Table 3, row 'Maximum CrCL' omega_CrCLMAX = 35.1 CV% (RSE 17%) -> 0.351^2
    etalcrcl_ega50 ~ 1.249924 # Table 3, row 'Time of half-maximum CrCL' omega_CrCL50 = 111.8 CV% (RSE 121%) -> 1.118^2

    # ---- Residual error (van Hasselt 2014 Table 3) ----
    # Methods eq (5) is purely proportional:
    #     CrCL_obs,ij(t) = CrCL_pred,ij(t) * (1 + eps_CrCL,ij)
    # Table 3 heads this block 'Residual unexplained variability variance', so
    # the tabulated 0.0291 is a VARIANCE and the SD entered here is its root.
    propSd <- 0.170587
    label("Proportional residual error on creatinine clearance (fraction)") # Table 3, row 'Proportional error' sigma_CrCL = 0.0291 (variance, RSE 43%) -> sqrt = 0.170587
  })

  model({
    # The rxode2 time variable IS maternal gestational age in weeks, matching
    # the paper's own CrCL(t) notation in eq (4) and the x-axis of Figure 3
    # ('time in gestational weeks'). Same convention as the shipped sibling
    # ontogeny models, which drive an age axis off `time` directly
    # (Wu_2024_gfr_maturation.R, Codaccioni_2024_cyp3a4_hepatic_ontogeny.R).
    # A user holding the canonical EGA covariate column instead can simply
    # solve this model on EGA as the time grid.
    ega <- time

    crcl_ega0 <- exp(lcrcl_ega0 + etalcrcl_ega0)
    crcl_matspan <- exp(lcrcl_matspan + etalcrcl_matspan)
    crcl_ega50 <- exp(lcrcl_ega50 + etalcrcl_ega50)

    # Methods eq (4). At ega = 0 this collapses exactly to the nonpregnant
    # baseline crcl_ega0, which is the anchor the sibling semiphysiological
    # PK model divides by to form its CrCL(t)/CrCL(0) ratio.
    CrCL <- crcl_ega0 + crcl_matspan * ega / (crcl_ega50 + ega)

    CrCL ~ prop(propSd)
  })
}
