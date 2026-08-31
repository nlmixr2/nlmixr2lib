Codaccioni_2024_cyp3a4_hepatic_ontogeny <- function() {
  description <- paste(
    "PBPK system-parameter (hepatic enzyme ontogeny) model. The two in vivo",
    "derived hepatic CYP3A4 ontogeny profiles compared by Codaccioni 2024",
    "within the Simcyp v22.1 paediatric P-PBPK population: the Salem",
    "profile (Salem 2014, deconvoluted from paediatric intravenous",
    "midazolam data) and the MODIFIED Upreti profile (Upreti and Wahlstrom",
    "2016, deconvoluted from paediatric sufentanil data, with three",
    "modifications introduced by Codaccioni 2024). Both are emitted side by",
    "side, together with their ratio, so the head-to-head comparison that",
    "is the subject of the paper can be run directly. This is a SYSTEM",
    "layer intended to be coupled to a PBPK drug model -- it contains no",
    "drug, no dosing, no compartments and no ODEs. Each output is an",
    "algebraic function of the rxode2 time variable, which the model",
    "interprets as postnatal age in YEARS, so a solve over time 0 to 25",
    "traces both profiles from birth to the adult plateau. The output is",
    "the fraction of adult hepatic CYP3A4 activity, which enters the Simcyp",
    "hepatic clearance cascade as the 'ontogeny fraction' multiplier (see",
    "the vignette; the remaining terms of that cascade are NOT encoded here",
    "because their age functions and the drug-specific intrinsic clearances",
    "live in Simcyp library files that the paper does not reproduce). The",
    "three Codaccioni modifications to Upreti are: fractional adult",
    "expression at birth raised from 0.10 to 0.15; the C3 term removed from",
    "the second (post-Age-Cap) equation; and the profile clamped to the",
    "adult value of 1 above 12.5 years, because the unclamped equation",
    "declines indefinitely below 1. Deterministic: Codaccioni 2024 report",
    "no variability, standard errors or confidence intervals for any",
    "ontogeny coefficient, so every parameter is fixed and the model",
    "carries no between-subject variability and no residual error.",
    sep = " "
  )
  reference <- paste(
    "Codaccioni M, Southall RL, Dinh J, Johnson TN.",
    "Prediction of pediatric pharmacokinetics for CYP3A4 metabolized drugs:",
    "comparison of the performance of two hepatic ontogeny within a",
    "physiologically based pharmacokinetic model.",
    "J Clin Pharmacol. 2024;64(9):1083-1094. doi:10.1002/jcph.2452.",
    "Equation (1) and Equation (2) are printed in Methods, section",
    "'Modified Upreti and Salem Hepatic CYP3A4 Ontogeny Functions';",
    "all coefficient values are from Supplementary Table S2.",
    "The two profiles originate in Salem F, Johnson TN, Abduljalil K,",
    "Tucker GT, Rostami-Hodjegan A. Clin Pharmacokinet.",
    "2014;53(7):625-636 (doi:10.1007/s40262-014-0140-7) and",
    "Upreti VV, Wahlstrom JL. J Clin Pharmacol. 2016;56(3):266-283",
    "(doi:10.1002/jcph.585); both are restated in full, with the",
    "Codaccioni modifications applied, in Supplementary Table S2 of the",
    "present paper, which is the source used here.",
    sep = " "
  )
  vignette <- "Codaccioni_2024_cyp3a4_hepatic_ontogeny"
  units <- list(
    time = paste(
      "year (postnatal age; the Salem profile is defined to its Age Cap of",
      "25 years and the modified Upreti profile is clamped to the adult",
      "value above 12.5 years, so a solve over 0 to 25 covers both)"
    ),
    dosing = "n/a (no exogenous dosing; hepatic enzyme ontogeny system model)",
    concentration = paste(
      "n/a (no drug concentration). Both outputs are the fraction of adult",
      "hepatic CYP3A4 activity and are unitless, as is their ratio."
    )
  )

  # No covariate columns: both profiles are functions of postnatal age
  # alone. Codaccioni 2024 Methods 'Comparative Dataset Building' record
  # that ethnicity and disease were considered when the two source
  # ontogeny profiles were derived, but neither profile carries an
  # ethnicity, sex or disease term -- the published functions depend on
  # age only.

  population <- list(
    species    = "human",
    n_subjects = 9L,
    n_studies  = 13L,
    age_range  = paste(
      "Birth to 15.71 years across the paediatric verification dataset",
      "(Codaccioni 2024 Table 1); the ontogeny functions themselves are",
      "defined from birth to the adult plateau."
    ),
    sex_female_pct = NA_real_,
    race_ethnicity = paste(
      "Not carried by either ontogeny function. The paediatric",
      "verification simulations used the Simcyp Sim-Pediatric population",
      "except for the Hamano 2019 midazolam scenario, which used",
      "Sim-Japanese Pediatric (Codaccioni 2024 Methods, 'Simulation Trial",
      "Design')."
    ),
    disease_state = paste(
      "Relatively healthy children undergoing minor or reshaping surgery,",
      "day-case surgery, or elective non-cardiac surgery, plus one",
      "status-epilepticus cohort, one malaria cohort, one single-ventricle",
      "cardiac-catheterisation cohort and one paediatric intensive care",
      "cohort (Codaccioni 2024 Table 1). Critically ill children,",
      "cardiopulmonary-bypass patients, ventilated children receiving",
      "midazolam, children on CYP3A4 comedication, children given",
      "halothane or isoflurane with fentanyl or sufentanil, and preterm",
      "infants below 36 gestational weeks were all excluded."
    ),
    dose_range = "n/a (no drug in this model)",
    notes = paste(
      "n_subjects = 9 is the number of full-term neonates less than 1 week",
      "old, pooled from two studies, whose individual intravenous",
      "midazolam clearances Codaccioni 2024 deconvoluted to raise the",
      "modified Upreti fractional adult expression at birth from 0.10 to",
      "0.15 (Methods, 'Modified Upreti and Salem Hepatic CYP3A4 Ontogeny",
      "Functions'). It is the only subject count the paper reports for a",
      "value that is newly derived here; every other coefficient in both",
      "profiles is inherited from Salem 2014 or Upreti and Wahlstrom 2016",
      "and restated in Supplementary Table S2. n_studies = 13 is the",
      "independent paediatric verification dataset, contributing 17",
      "single-dose clinical scenarios across four intravenous CYP3A4",
      "substrates (alfentanil, fentanyl, midazolam, sildenafil), each",
      "selected for fmCYP3A4 at or above 80% in adults and a",
      "low-to-intermediate hepatic extraction ratio (Codaccioni 2024",
      "Methods and Table 1). Those 13 studies verify the profiles; they do",
      "not contribute parameters to them."
    )
  )

  ini({
    # -----------------------------------------------------------------
    # Every coefficient below is a published point estimate read from
    # Codaccioni 2024 Supplementary Table S2. The paper reports no
    # standard error, RSE or confidence interval for any of them, so all
    # are encoded as fixed().
    # -----------------------------------------------------------------

    # ---- Salem hepatic CYP3A4 ontogeny (Eq. 1 only) ----
    # Fraction of adult = (FMax - FBirth) * Age^n / (Age50^n + Age^n) + FBirth
    salem_birth <- fixed(0.11)
    label("Salem profile: fraction of adult hepatic CYP3A4 activity at birth (unitless)")  # Table S2, F_Birth, Salem column
    salem_max <- fixed(1.06)
    label("Salem profile: maximal fraction of adult hepatic CYP3A4 activity (unitless)")  # Table S2, F_Max, Salem column
    salem_t50 <- fixed(0.64)
    label("Salem profile: postnatal age at half-maximal CYP3A4 activity (years)")  # Table S2, Age_50, Salem column
    salem_hill <- fixed(1.91)
    label("Salem profile: Hill exponent on postnatal age (unitless)")  # Table S2, n, Salem column
    salem_agecap <- fixed(25)
    label("Salem profile: postnatal age above which Eq. 1 is held at its capped value (years)")  # Table S2, Age Cap, Salem column

    # ---- Modified Upreti hepatic CYP3A4 ontogeny (Eq. 1 then Eq. 2) ----
    # Eq. 1, to the Age Cap:
    #   Fraction of adult = (FMax - FBirth) * Age^n / (Age50^n + Age^n) + FBirth
    # Eq. 2, above the Age Cap:
    #   Fraction of adult = C0 + C1 * exp(C2 * (Age - AgeCap))
    upreti_birth <- fixed(0.15)
    label("Modified Upreti profile: fraction of adult hepatic CYP3A4 activity at birth (unitless)")  # Table S2, F_Birth, Modified Upreti column; raised from the published Upreti value of 0.10 by deconvoluting IV midazolam clearance in 9 full-term neonates (Methods)
    upreti_max <- fixed(1.7)
    label("Modified Upreti profile: maximal fraction of adult hepatic CYP3A4 activity (unitless)")  # Table S2, F_Max, Modified Upreti column
    upreti_t50 <- fixed(0.1)
    label("Modified Upreti profile: postnatal age at half-maximal CYP3A4 activity (years)")  # Table S2, Age_50, Modified Upreti column
    upreti_hill <- fixed(1.3)
    label("Modified Upreti profile: Hill exponent on postnatal age (unitless)")  # Table S2, n, Modified Upreti column
    upreti_agecap <- fixed(2.5)
    label("Modified Upreti profile: postnatal age at which Eq. 1 hands over to Eq. 2 (years)")  # Table S2, Age Cap, Modified Upreti column; Methods states Eq. 1 applies "up to 2.5 years"
    upreti_c0 <- fixed(0.7)
    label("Modified Upreti profile: additive constant C0 of the post-Age-Cap decline (unitless)")  # Table S2, C_0, Modified Upreti column
    upreti_c1 <- fixed(1)
    label("Modified Upreti profile: amplitude C1 of the post-Age-Cap exponential decline (unitless)")  # Table S2, C_1, Modified Upreti column
    upreti_c2 <- fixed(-0.1)
    label("Modified Upreti profile: rate constant C2 of the post-Age-Cap exponential decline (1/year)")  # Table S2, C_2, Modified Upreti column
    upreti_agecap2 <- fixed(12.5)
    label("Modified Upreti profile: postnatal age above which the profile is clamped to the adult value (years)")  # Table S2, "Age Cap (second equation)", Modified Upreti column
    upreti_adult <- fixed(1)
    label("Modified Upreti profile: clamped fraction of adult CYP3A4 activity above 12.5 years (unitless)")  # Methods: "after age 12.5, the fractional activity of CYP3A4 in children was set to be comparable to adults", i.e. a fraction of adult of exactly 1
  })

  model({
    # The rxode2 time variable IS postnatal age in years. Solve over 0 to
    # 25 to trace both profiles from birth through to the adult plateau.
    age <- time

    # ---- Salem hepatic CYP3A4 ontogeny (Eq. 1) ----
    # Table S2 gives Salem a single Age Cap of 25 years and no C0 / C1 /
    # C2, so Salem is Eq. 1 throughout, evaluated at the capped age. The
    # cap barely bites: Eq. 1 at 25 years is 1.0591, against the F_Max
    # asymptote of 1.06, so holding the value from 25 years onward and
    # letting it run to F_Max differ by 0.08%.
    age_salem <- (age < salem_agecap) * age +
      (age >= salem_agecap) * salem_agecap
    fcyp3a4_salem <- (salem_max - salem_birth) * age_salem^salem_hill /
      (salem_t50^salem_hill + age_salem^salem_hill) + salem_birth

    # ---- Modified Upreti hepatic CYP3A4 ontogeny (Eq. 1, then Eq. 2, then clamp) ----
    # Eq. 1 to the 2.5-year Age Cap. The capped age keeps the branch
    # finite above the cap, where it is multiplied out by zero anyway.
    age_upreti <- (age < upreti_agecap) * age +
      (age >= upreti_agecap) * upreti_agecap
    f_upreti_eq1 <- (upreti_max - upreti_birth) * age_upreti^upreti_hill /
      (upreti_t50^upreti_hill + age_upreti^upreti_hill) + upreti_birth

    # Eq. 2 from the Age Cap to 12.5 years. The published Upreti form
    # carries a further C3 term; Codaccioni 2024 removed it ("removing
    # the C3 term from the published Upreti model", Methods), which is
    # equivalent to C3 = 0, so it does not appear here.
    f_upreti_eq2 <- upreti_c0 +
      upreti_c1 * exp(upreti_c2 * (age - upreti_agecap))

    # Above 12.5 years the profile is clamped to the adult value. Without
    # the clamp Eq. 2 keeps falling, crossing 1 at about 14.5 years and
    # tending to C0 = 0.7, which is why Codaccioni 2024 impose it.
    fcyp3a4_upreti <-
      (age <= upreti_agecap) * f_upreti_eq1 +
      (age > upreti_agecap) * (age <= upreti_agecap2) * f_upreti_eq2 +
      (age > upreti_agecap2) * upreti_adult

    # ---- Salem-to-modified-Upreti ratio ----
    # Every other term of the Simcyp hepatic clearance cascade (CLint,
    # CYP3A4 abundance, MPPGL, liver weight, Q_H, fu_B) is identical
    # between the two simulations of a given scenario, because only the
    # ontogeny profile was swapped and the same fixed random seed was
    # retained (Methods, 'Simulation Trial Design'). This ratio is
    # therefore the whole of the difference in unbound intrinsic
    # clearance per liver between the two arms at a given age, and in the
    # low-extraction-ratio limit it is also the ratio of their predicted
    # hepatic clearances.
    fcyp3a4_ratio_salem_upreti <- fcyp3a4_salem / fcyp3a4_upreti
  })
}
