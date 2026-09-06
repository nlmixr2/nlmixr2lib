Miyagi_2011_ugt1a1_ugt1a6_hepatic_ontogeny <- function() {
  description <- paste(
    "Hepatic enzyme ontogeny + IVIVE system model. Postnatal development of",
    "UDP-glucuronosyltransferase (UGT) 1A1 and 1A6 catalytic activity in 50",
    "individual pediatric human liver microsome preparations spanning 13 days",
    "to 20 years (Miyagi and Collier 2011), together with the in vitro to in",
    "vivo extrapolation (IVIVE) cascade the authors used to turn those",
    "activities into normalized hepatic clearances. Each isoform's activity is",
    "a one-phase exponential association in postnatal age (Eq. 1), rising from",
    "a non-zero activity at birth to an adult plateau: UGT1A1 (bilirubin",
    "glucuronidation) matures fast, reaching 90% of plateau by 3.8 months,",
    "whereas UGT1A6 (serotonin glucuronidation) starts above half of adult",
    "activity and reaches 90% of plateau only at 14 months, which is the",
    "paper's central finding of independent regulation of the two isoforms.",
    "Activity is converted to an unbound intrinsic clearance by",
    "Michaelis-Menten back-calculation at the assay substrate concentration",
    "with literature Km values, scaled to the whole liver through microsomal",
    "protein per gram liver and liver weight, and passed through BOTH the well",
    "stirred (Eq. 2) and parallel-tube (Eq. 3) hepatic extraction models, then",
    "allometrically scaled to the child by weight (Eq. 4) and normalized to",
    "body weight. The Simcyp Pediatric physiology functions the paper also",
    "used are emitted alongside: liver weight from body surface area (Eq. 5),",
    "microsomal protein per gram liver versus age (Eq. 7), plasma albumin",
    "versus age (Eq. 8), the pediatric unbound fraction that follows from it",
    "(Eq. 9) and body surface area from height and weight (Eq. 10). This is a",
    "SYSTEM layer, not a drug model: it has no dosing, no compartments and no",
    "ODEs, and every output is an algebraic function of the rxode2 time",
    "variable, which the model interprets as postnatal age in MONTHS.",
    "Deterministic: the paper reports no between-subject variance component",
    "for any quantity, so every parameter is fixed and there is no residual",
    "error. IMPORTANT: the paper's Eq. 6 (Simcyp hepatic blood flow versus",
    "age) is NOT encoded, because it is non-physical as printed -- see the",
    "vignette Errata and the note in model() below -- so only the allometric",
    "clearance branch is emitted.",
    sep = " "
  )
  reference <- paste(
    "Miyagi SJ, Collier AC.",
    "The development of UDP-glucuronosyltransferases 1A1 and 1A6 in the",
    "pediatric liver.",
    "Drug Metab Dispos. 2011;39(5):912-919. doi:10.1124/dmd.110.037192.",
    "PMCID: PMC3082376.",
    "Equations 1 to 10 are printed in Methods, section 'Pharmacokinetics,",
    "Scaling, and Statistical Analyses'; the fitted activity parameters are in",
    "Results, section 'UGT1A1 and 1A6 Activity' and in the Fig. 1A and Fig. 2A",
    "captions; the modeled clearances are in Table 1.",
    sep = " "
  )
  vignette <- "Miyagi_2011_ugt1a1_ugt1a6_hepatic_ontogeny"
  units <- list(
    time = paste(
      "month (postnatal age). The donor livers span 13 days (0.427 months) to",
      "20 years (240 months); Eq. 8 contains LN(age) and is singular at age 0,",
      "so solve from the youngest donor age rather than from 0."
    ),
    dosing = "n/a (no exogenous dosing; hepatic enzyme ontogeny system model)",
    concentration = paste(
      "n/a (no drug concentration). Enzyme activities are nmol substrate",
      "metabolized per min per mg microsomal protein; normalized hepatic",
      "clearances are l/h/kg body weight."
    )
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight of the liver donor",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Used twice: as the allometric scalar in Eq. 4, (WT / 70)^0.75, and as",
        "the normalization denominator that turns hepatic clearance in l/h",
        "into the normalized l/h/kg the paper reports. Also feeds Eq. 10 (body",
        "surface area). Miyagi 2011 Methods state that individual children's",
        "weights were used, with the 50th percentile for age and gender from",
        "the National Center for Health Statistics (2000) growth charts",
        "substituted for the 8 donors whose weight was missing."
      ),
      source_name        = "Wi"
    ),
    HT = list(
      description        = "Height (body length) of the liver donor",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Feeds Eq. 10 (body surface area) only, which in turn drives the",
        "Simcyp liver-weight function Eq. 5. Miyagi 2011 Methods state that",
        "the 50th percentile for age and gender from the National Center for",
        "Health Statistics (2000) growth charts was substituted for the 15",
        "donors whose height and/or weight was missing."
      ),
      source_name        = "Height"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "categorical",
      notes       = paste(
        "Screened but not retained. Miyagi 2011 Results: 'Neither UGT1A1 nor",
        "1A6 activities differed significantly with gender or ethnicity', and",
        "no sex coefficient is reported for any equation. Sex enters the",
        "paper only through the growth-chart percentile used to impute a",
        "missing weight or height, which is a data-preparation step rather",
        "than a model term."
      )
    ),
    RACE_WHITE = list(
      description = "White race indicator",
      units       = "(binary)",
      type        = "categorical",
      notes       = "Screened but not retained; see the SEXF note above. 64% of the 50 pediatric donors."
    ),
    RACE_BLACK = list(
      description = "Black or African American race indicator",
      units       = "(binary)",
      type        = "categorical",
      notes       = "Screened but not retained; see the SEXF note above. 16% of the 50 pediatric donors."
    ),
    RACE_HISPANIC = list(
      description = "Hispanic ethnicity indicator",
      units       = "(binary)",
      type        = "categorical",
      notes       = "Screened but not retained; see the SEXF note above. 14% of the 50 pediatric donors."
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units       = "(binary)",
      type        = "categorical",
      notes       = "Screened but not retained; see the SEXF note above. 3% of the 50 pediatric donors."
    ),
    RACE_INDIAN = list(
      description = "American Indian race indicator",
      units       = "(binary)",
      type        = "categorical",
      notes       = "Screened but not retained; see the SEXF note above. 3% of the 50 pediatric donors."
    )
  )

  population <- list(
    species    = "human",
    n_subjects = 50L,
    n_studies  = 1L,
    age_range  = "13 days to 20 years (mean 8.0 years)",
    weight_range = paste(
      "Not reported. Miyagi 2011 state that individual donor weights were",
      "used where available and that National Center for Health Statistics",
      "(2000) 50th-percentile-for-age-and-gender values were substituted for",
      "8 donors (weight) and 15 donors (height and/or weight), but no",
      "weight summary statistic or range is published."
    ),
    sex_female_pct = 30,
    race_ethnicity = paste(
      "Pediatric donors (n = 50): white 64%, African American 16%, Hispanic",
      "14%, Asian 3%, American Indian 3%. Neither activity differed",
      "significantly with ethnicity."
    ),
    disease_state = paste(
      "None. Postmortem donors with healthy livers obtained through the",
      "United Network for Organ Sharing, free of infectious disease, with",
      "causes of death of anoxia, cerebrovascular aneurysm, head trauma or",
      "motor vehicle accident, processed within 8 h postmortem."
    ),
    dose_range = "n/a (no drug in this model)",
    notes = paste(
      "15 female and 35 male pediatric donors. The 'adult' reference values",
      "come from a separate pooled adult human liver microsome preparation",
      "(XTreme 200 Pool, XenoTech) of 200 donors, 100 men and 100 women,",
      "mean age 50 years (range 17 to 78), white 84%, African American 6%,",
      "Hispanic 6%, Asian 4%. That pool measured 0.58 +/- 0.056 nmol/min/mg",
      "for bilirubin (UGT1A1) and 4.3 +/- 0.49 nmol/min/mg for serotonin",
      "(UGT1A6), against the fitted pediatric plateaus of 0.7690 and 4.737",
      "encoded below. UGT1A1 and 1A6 PROTEIN expression was also measured by",
      "Western blot and is deliberately NOT modelled here: it showed no age",
      "dependence over any age window and did not correlate with activity",
      "for either isoform, which is a null result with no fitted equation."
    )
  )

  ini({
    # -----------------------------------------------------------------
    # Every value below is a published point estimate from Miyagi 2011.
    # No variance component, standard error or confidence interval is
    # reported for any quantity that enters the model as a coefficient,
    # so all parameters are fixed and the model carries no IIV and no
    # residual error. The two maturation rate constants are the only
    # derived values; their derivation is stated on the line.
    # -----------------------------------------------------------------

    # ---- UGT1A1 activity ontogeny (Eq. 1, bilirubin glucuronidation) ----
    # Y = Y0 + (Ymax - Y0) * (1 - exp(-K * age))
    act0_ugt1a1 <- fixed(0.0007848)
    label("UGT1A1: bilirubin glucuronidation activity at birth (nmol/min/mg microsomal protein)")  # Results 'UGT1A1 and 1A6 Activity' and Fig. 1A caption: 0.7848 pmol/min/mg protein, i.e. 0.0007848 nmol/min/mg in the units used for the plateau
    actmax_ugt1a1 <- fixed(0.7690)
    label("UGT1A1: maximal (adult plateau) bilirubin glucuronidation activity (nmol/min/mg microsomal protein)")  # Results and Fig. 1A caption: 0.7690 +/- 0.081 nmol/min/mg protein [95% CI 0.6028-0.9351]
    kmat_ugt1a1 <- fixed(0.605675)
    label("UGT1A1: first-order maturation rate constant of activity (1/month)")  # DERIVED, not printed: Miyagi 2011 report the Eq. 1 fit only through Y0, Ymax and the '90th percentile' age, defined in Methods as "the youngest age at which activity reached within 90% of the plateau level". Solving Y(t90) = 0.9 * Ymax for K with Y0 = 0.0007848, Ymax = 0.7690 and t90 = 3.8 months (Results, Fig. 1A caption) gives K = 0.605675 /month exactly. See vignette 'Derived quantities'.

    # ---- UGT1A6 activity ontogeny (Eq. 1, serotonin glucuronidation) ----
    act0_ugt1a6 <- fixed(2.560)
    label("UGT1A6: serotonin glucuronidation activity at birth (nmol/min/mg microsomal protein)")  # Results 'UGT1A1 and 1A6 Activity' and Fig. 2A caption: 2.560 nmol/min/mg protein at birth
    actmax_ugt1a6 <- fixed(4.737)
    label("UGT1A6: maximal (adult plateau) serotonin glucuronidation activity (nmol/min/mg microsomal protein)")  # Results and Fig. 2A caption: 4.737 +/- 0.33 nmol/min/mg protein (95% CI 4.041-5.434)
    kmat_ugt1a6 <- fixed(0.108938)
    label("UGT1A6: first-order maturation rate constant of activity (1/month)")  # DERIVED, not printed: same construction as kmat_ugt1a1, with Y0 = 2.560, Ymax = 4.737 and t90 = 14 months (Results, Fig. 2A caption), giving K = 0.108938 /month exactly. See vignette 'Derived quantities'.

    # ---- In vitro assay conditions and enzyme kinetics ----
    km_bilirubin <- fixed(5.0)
    label("Michaelis constant of bilirubin for UGT1A1 (micromolar)")  # Methods 'Pharmacokinetics, Scaling...': "The literature Km values, which were used to calculate intrinsic clearances, were 5.0 uM for bilirubin (UGT1A1) (Ciotti et al., 1998)"
    sub_bilirubin <- fixed(125)
    label("Bilirubin concentration in the UGT1A1 activity assay (micromolar)")  # Methods: "Michaelis-Menten kinetics (125 uM bilirubin/1A1 ...)"; matches the 0.125 mM bilirubin of the UGT1A1 assay
    fu_bilirubin <- fixed(0.001)
    label("Adult plasma unbound fraction of bilirubin (unitless)")  # Methods: "Plasma unbound fraction for bilirubin and serotonin were 0.001 (Ostrow et al., 2003) and 0.17 ... respectively"
    km_serotonin <- fixed(5200)
    label("Michaelis constant of serotonin for UGT1A6 (micromolar)")  # Methods: "5.2 mM for serotonin (UGT1A6) (Krishnaswamy et al., 2003)", i.e. 5200 uM
    sub_serotonin <- fixed(100)
    label("Serotonin concentration in the UGT1A6 activity assay (micromolar)")  # Methods: "100 uM serotonin/1A6"; chosen from the linear portion of the Michaelis-Menten curve because the Km exceeds physiological serotonin concentrations
    fu_serotonin <- fixed(0.17)
    label("Adult plasma unbound fraction of serotonin (unitless)")  # Methods: 0.17 (Breyer-Pfaff et al., 1989)

    # ---- Adult physiology used for the allometric branch ----
    mppgl_adult <- fixed(45)
    label("Adult microsomal protein per gram of liver (mg/g)")  # Methods: "Microsomal protein per gram liver (MPPGL) was unknown, so the standard variable of 45 mg/g was used (Houston, 1994)"
    liverwt_adult <- fixed(1500)
    label("Adult liver weight (g)")  # Methods: "assuming a liver size of 1500 g ... for adults"
    qh_adult <- fixed(90)
    label("Adult hepatic blood flow (l/h)")  # Methods: "a hepatic flow rate of 1.5 l/min for adults", i.e. 90 l/h
    wt_adult <- fixed(70)
    label("Reference adult body weight for allometric scaling, Wstd (kg)")  # Methods define Wstd as "the weight of an average adult (20 years of age)" taken as "the 50th percentile at 20 years for each gender" but never print a number. 70 kg is the rounded standard; it is corroborated by reproducing the paper's own pooled adult UGT1A1 clearance (see vignette).
    allom_exp <- fixed(0.75)
    label("Allometric exponent on body weight in Eq. 4 (unitless)")  # Eq. 4 prints the exponent as the fraction 3/4; Discussion: "the allometric model to three-forths power scaling was chosen"

    # ---- Simcyp Pediatric physiology functions (Eqs. 5, 7, 8, 9, 10) ----
    liverwt_coef <- fixed(0.722)
    label("Eq. 5 coefficient relating liver weight to body surface area (kg/m^2.352)")  # Eq. 5: Liver Size = (Body Surface Area)^1.176 x 0.722
    liverwt_exp <- fixed(1.176)
    label("Eq. 5 exponent on body surface area (unitless)")  # Eq. 5
    mppgl_c0 <- fixed(1.407)
    label("Eq. 7 intercept of log10 microsomal protein per gram liver (unitless)")  # Eq. 7: MPPGL (mg/g) = 10^(1.407 + 0.0158*Age - 0.000382*Age^2 + 0.0000024*Age^3), Age in years
    mppgl_c1 <- fixed(0.0158)
    label("Eq. 7 linear age coefficient of log10 MPPGL (1/year)")  # Eq. 7
    mppgl_c2 <- fixed(-0.000382)
    label("Eq. 7 quadratic age coefficient of log10 MPPGL (1/year^2)")  # Eq. 7
    mppgl_c3 <- fixed(0.0000024)
    label("Eq. 7 cubic age coefficient of log10 MPPGL (1/year^3)")  # Eq. 7
    alb_slope <- fixed(1.1287)
    label("Eq. 8 slope of pediatric plasma albumin on natural-log age (g/l per ln-year)")  # Eq. 8: [P]Pediatric (g/l) = 1.1287 x LN(Age) + 33.746, Age in years
    alb_intercept <- fixed(33.746)
    label("Eq. 8 intercept of pediatric plasma albumin (g/l)")  # Eq. 8
    alb_adult <- fixed(44)
    label("Adult plasma albumin concentration, [P]adult (g/l)")  # Methods: "where [P]adult is 44 g/l (McNamara and Alcorn, 2002)"
    bsa_coef <- fixed(0.007184)
    label("Eq. 10 body surface area coefficient (m^2 / (cm^0.725 kg^0.425))")  # Eq. 10: BSA (m^2) = 0.007184 x Height^0.725 x Weight^0.425, height in cm and weight in kg
    bsa_exp_ht <- fixed(0.725)
    label("Eq. 10 exponent on height (unitless)")  # Eq. 10
    bsa_exp_wt <- fixed(0.425)
    label("Eq. 10 exponent on weight (unitless)")  # Eq. 10
  })

  model({
    # The rxode2 time variable IS postnatal age in MONTHS. The donor
    # livers span 13 days (0.427 months) to 20 years (240 months).
    age_month <- time
    # Eqs. 7 and 8 are written in years; Eq. 1 is fitted in months.
    age_year <- age_month / 12

    # =================================================================
    # 1. Enzyme activity ontogeny -- Eq. 1, the paper's central result
    # =================================================================
    # One-phase exponential association: an enzyme starts at some
    # non-zero activity at birth and rises at a constant fractional
    # rate K to a plateau. Miyagi 2011 tested this against a
    # zero-intercept exponential, a biphasic rise-then-fall and a
    # sigmoidal form, and F tests selected this one for both isoforms.
    act_ugt1a1 <- act0_ugt1a1 +
      (actmax_ugt1a1 - act0_ugt1a1) * (1 - exp(-kmat_ugt1a1 * age_month))
    act_ugt1a6 <- act0_ugt1a6 +
      (actmax_ugt1a6 - act0_ugt1a6) * (1 - exp(-kmat_ugt1a6 * age_month))

    # Fraction of the adult plateau, the form in which an ontogeny
    # profile is normally consumed by a PBPK clearance cascade.
    fmat_ugt1a1 <- act_ugt1a1 / actmax_ugt1a1
    fmat_ugt1a6 <- act_ugt1a6 / actmax_ugt1a6

    # =================================================================
    # 2. Activity -> unbound intrinsic clearance (Michaelis-Menten)
    # =================================================================
    # Methods: "Vmax was calculated for each sample within the study by
    # using experimental enzyme activity, substrate concentrations, and
    # Michaelis-Menten kinetics". Inverting v = Vmax*S/(Km+S) gives
    # Vmax = v*(Km+S)/S, and CLint = Vmax/Km. Because a micromolar
    # concentration is nmol/ml, CLint comes out in ml/min/mg protein.
    vmax_ugt1a1 <- act_ugt1a1 * (km_bilirubin + sub_bilirubin) / sub_bilirubin
    vmax_ugt1a6 <- act_ugt1a6 * (km_serotonin + sub_serotonin) / sub_serotonin
    clint_ugt1a1 <- vmax_ugt1a1 / km_bilirubin
    clint_ugt1a6 <- vmax_ugt1a6 / km_serotonin

    # Scale to the whole adult liver: ml/min/mg x mg/g x g -> ml/min,
    # then /1000 to l/min. This is the adult-physiology scaling that
    # the allometric branch uses; the Simcyp branch would instead use
    # the individualized mppgl_simcyp and liverwt_simcyp emitted below.
    clint_liver_ugt1a1 <- clint_ugt1a1 * mppgl_adult * liverwt_adult / 1000
    clint_liver_ugt1a6 <- clint_ugt1a6 * mppgl_adult * liverwt_adult / 1000

    # =================================================================
    # 3. Hepatic extraction -- Eq. 2 (well stirred) and Eq. 3 (parallel tube)
    # =================================================================
    # Both are written with hepatic blood flow in l/min to match the
    # paper's stated adult flow of 1.5 l/min; qh_adult is held in l/h,
    # so divide by 60 here and multiply the result back by 60.
    qh_min <- qh_adult / 60
    clh_ugt1a1_ws <- qh_min * fu_bilirubin * clint_liver_ugt1a1 /
      (qh_min + fu_bilirubin * clint_liver_ugt1a1)
    clh_ugt1a6_ws <- qh_min * fu_serotonin * clint_liver_ugt1a6 /
      (qh_min + fu_serotonin * clint_liver_ugt1a6)
    clh_ugt1a1_pt <- qh_min *
      (1 - exp(-clint_liver_ugt1a1 * fu_bilirubin / qh_min))
    clh_ugt1a6_pt <- qh_min *
      (1 - exp(-clint_liver_ugt1a6 * fu_serotonin / qh_min))

    # =================================================================
    # 4. Allometric pediatric scaling -- Eq. 4, then normalize by weight
    # =================================================================
    # Eq. 4: CLpediatric = CLhepatic * (Wi / Wstd)^(3/4). Miyagi 2011
    # then report every clearance "normalized to individual body
    # weights", which is the division by WT that follows. Units go
    # l/min -> l/h (x 60) -> l/h/kg.
    allom <- (WT / wt_adult)^allom_exp
    clnorm_ugt1a1_ws <- clh_ugt1a1_ws * allom * 60 / WT
    clnorm_ugt1a1_pt <- clh_ugt1a1_pt * allom * 60 / WT
    clnorm_ugt1a6_ws <- clh_ugt1a6_ws * allom * 60 / WT
    clnorm_ugt1a6_pt <- clh_ugt1a6_pt * allom * 60 / WT

    # =================================================================
    # 5. Simcyp Pediatric physiology -- Eqs. 10, 5, 7, 8, 9
    # =================================================================
    # Emitted as system quantities. They are the individualized
    # physiology the paper's second scaling branch used in place of the
    # fixed adult constants above.

    # Eq. 10, body surface area from height (cm) and weight (kg).
    bsa <- bsa_coef * HT^bsa_exp_ht * WT^bsa_exp_wt

    # Eq. 5, liver weight from body surface area. The paper prints the
    # left-hand side as "Liver Size (g)", but the equation returns
    # KILOGRAMS: at the adult anchor the paper itself uses elsewhere
    # (1.85 m^2 for a 175 cm, 70 kg adult) it gives 1.487, against the
    # "liver size of 1500 g" the same Methods paragraph assumes. A
    # gram-valued reading would make an adult liver 1.5 g. The printed
    # unit is therefore a typo; the value is converted to grams here so
    # it is interchangeable with liverwt_adult.
    liverwt_simcyp <- liverwt_coef * bsa^liverwt_exp * 1000

    # Eq. 7, microsomal protein per gram liver versus age in years.
    mppgl_simcyp <- 10^(mppgl_c0 + mppgl_c1 * age_year +
                          mppgl_c2 * age_year^2 + mppgl_c3 * age_year^3)

    # Eq. 8, pediatric plasma albumin versus age in years. Singular at
    # age 0; the youngest donor is 13 days (0.0356 years).
    alb_simcyp <- alb_slope * log(age_year) + alb_intercept

    # Eq. 9, pediatric unbound fraction implied by that albumin. The
    # form is self-consistent: setting alb_simcyp = alb_adult returns
    # the adult fu exactly.
    fu_ugt1a1_simcyp <- 1 /
      (1 + (1 - fu_bilirubin) * alb_simcyp / (alb_adult * fu_bilirubin))
    fu_ugt1a6_simcyp <- 1 /
      (1 + (1 - fu_serotonin) * alb_simcyp / (alb_adult * fu_serotonin))

    # -----------------------------------------------------------------
    # Eq. 6 (Simcyp hepatic blood flow versus age) is deliberately NOT
    # encoded, and therefore neither is the Simcyp clearance branch
    # that consumes it. As printed,
    #   Qhepatic (l/h) = 0.265 x 10^(-0.6492 + 1.943*Age - 0.8118*Age^2
    #                                 + 0.08891*Age^3)
    # with Age in years is non-physical across the whole donor age
    # range: 0.059 l/h at birth and 0.99 l/h at 1 year against a
    # physiologic 12 to 30 l/h, a spurious maximum of 1.33 l/h at 2
    # years, a fall to 0.20 l/h by 5 years, then divergence to 9e25 l/h
    # at 10 years and to infinity by 20 years. It is non-monotone and
    # unbounded, so it cannot be a hepatic blood flow. No reading of
    # "Age" recovers a physiologic profile (log10 and natural log, in
    # years, months or weeks, and base e in place of base 10 were all
    # checked), and the paper attributes Eq. 6 to the commercial Simcyp
    # Pediatric platform rather than to any citable reference, so the
    # correct coefficients are not recoverable from any source on disk.
    # Note that Eq. 7, printed immediately below it in the same style,
    # IS correct and is independently corroborated by
    # Han_2025_alfentanil_pbpk.R, which encodes the same polynomial.
    # See the vignette Errata.
    # -----------------------------------------------------------------
  })
}
