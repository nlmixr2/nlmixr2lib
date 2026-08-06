Wahlquist_2024_propofol <- function() {
  description <- "Three-compartment population PK model for propofol with a symbolic-regression (neural-network) derived covariate model on all four inter-compartmental rate constants, the elimination rate constant and the central volume; trained on the pooled Eleveld 30-study dataset spanning neonates to the elderly" # nolint: line_length_linter.
  reference <- "Wahlquist Y, Sundell J, Soltesz K. Learning pharmacometric covariate model structures with symbolic regression networks. J Pharmacokinet Pharmacodyn. 2024;51(2):155-167. doi:10.1007/s10928-023-09887-3. PMCID: PMC11416364. Code and data deposit (cited as reference 27 of the paper): https://github.com/wahlquisty/learning-pharmacometric-covariate-structures" # nolint: line_length_linter.
  vignette <- "Wahlquist_2024_propofol"
  units <- list(time = "s", dosing = "mg", concentration = "mg/L")
  # Time is seconds because every rate constant in Eq. 10 is published in
  # s^-1; the values are kept exactly as printed rather than rescaled.
  # Concentration mg/L is numerically identical to the ug/mL used throughout
  # the paper (central holds mg, vc is in L).

  covariateData <- list(
    AGE = list(
      description = "Chronological age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters every expression normalised as AGE / AGEMAX with AGEMAX = 88 years",
        "(Results section following Eq. 10). Dataset range 27 weeks postmenstrual",
        "to 88 years (Table 1 gives 0-88 years)."
      ),
      source_name = "AGE"
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters normalised as WGT / WGTMAX with WGTMAX = 160 kg (Results section",
        "following Eq. 10). Dataset range 0.68-160 kg (Table 1). Baseline, not",
        "time-varying, in the source dataset."
      ),
      source_name = "WGT"
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters normalised as BMI / BMIMAX. The paper prints BMIMAX = 52.8 kg/m^2;",
        "the authors' code (get_results.jl in the reference-27 deposit) uses the",
        "full-precision dataset maximum 52.84713965, which is what this model uses.",
        "Dataset range 6.2-52.8 kg/m^2 (Table 1)."
      ),
      source_name = "BMI"
    ),
    SEXF = list(
      description = "Sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "The paper prints separate male and female expressions for k13 (Eqs. 10c,",
        "10d) and k31 (Eqs. 10f, 10g). The two forms differ only in the sign of the",
        "AGE-linear terms, so they are encoded here with a single sign switch",
        "sexsign = 1 - 2 * SEXF (+1 male, -1 female). This reproduces both printed",
        "expressions exactly. Internally the network coded sex as -0.5 male /",
        "+0.5 female (Covariate model section; get_results.jl)."
      ),
      source_name = "M1F2 (1 = male, 2 = female); SEXF = M1F2 - 1"
    )
  )

  covariatesDataExcluded <- list(
    SAMPLE_ARTERIAL = list(
      description = "Blood sampling site, 1 = arterial, 0 = venous",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (venous)",
      notes = paste(
        "Offered to the symbolic-regression network as one of the five candidate",
        "covariates (Table 1; 727 arterial / 306 venous subjects) but automatically",
        "pruned out of the final model: 'The blood sampling site (arterial or venous)",
        "was available as a modeling covariate, but was automatically pruned by the",
        "symbolic regression algorithm' (Results, after Eq. 10). Documentation only -",
        "not referenced in model(), and deliberately NOT proposed as a canonical",
        "covariate column, because no retained model uses it. Distinct from the",
        "registered SAMPLE_CAPILLARY canonical, which is capillary-vs-venous."
      ),
      source_name = "A1V2 (1 = arterial, 2 = venous)"
    )
  )

  compartmentData <- list(
    central = list(analyte = "propofol", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "propofol", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "propofol", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 1031,
    n_studies = 30,
    age_range = "27 weeks postmenstrual to 88 years",
    weight_range = "0.68-160 kg",
    bmi_range = "6.2-52.8 kg/m^2",
    sex_female_pct = 35,
    disease_state = "surgical and volunteer cohorts receiving propofol for sedation or general anaesthesia",
    dose_range = "intravenous bolus and infusion regimens as administered in the 30 pooled studies",
    notes = paste(
      "The pooled dataset published by Eleveld et al. (2018)",
      "doi:10.1016/j.bja.2018.01.018, comprising 15,433 plasma concentration",
      "observations (11,530 arterial, 3,903 venous) from 1,031 individuals",
      "(670 male, 361 female) in 30 clinical studies. Eleveld's original release",
      "holds 1,033 individuals; two were excluded here for having no observations",
      "(Data set section, footnote 1). Table 1 lists the covariate ranges."
    )
  )

  ini({
    # --------------------------------------------------------------------
    # All coefficients are point estimates from a deterministic (non-mixed-
    # effects) gradient-based fit, so every one is fixed(). The paper reports
    # no standard errors, no inter-individual variability and no residual-
    # error model; see the vignette's Assumptions and deviations section.
    #
    # Notation used below (Results section following Eq. 10):
    #   agen = AGE / 88          wtn = WT / 160          bmin = BMI / 52.84713965
    # --------------------------------------------------------------------

    # -- k10, elimination rate constant (Eq. 10a) ------------------------
    lkel <- fixed(log(0.00342))
    label("k10 intercept of the covariate expression (1/s)") # Eq. 10a
    e_wt_kel <- fixed(0.00441)
    label("k10 slope on normalised weight (1/s)") # Eq. 10a

    # -- k12, central -> peripheral1 (Eq. 10b) ---------------------------
    # Eq. 10b is printed WITHOUT the network output normalisation constant
    # that Eqs. 10a and 10c-10h all carry; k12_norm restores it. See the
    # vignette Errata: the factor is the k12 entry of
    # `pkparams_normalization` in checkresults/volumesclearances.jl of the
    # paper's own reference-27 code deposit, and applying it reproduces the
    # authors' deposited per-observation predictions (median ratio 1.006).
    k12_norm <- fixed(0.013328808)
    label("k12 network-output normalisation constant (1/s)") # reference-27 deposit, not printed in the paper
    e_age_k12_num <- fixed(0.158)
    label("k12 numerator coefficient on normalised age (unitless)") # Eq. 10b
    e_bmi_k12_num <- fixed(-0.00431)
    label("k12 numerator coefficient on normalised BMI (unitless)") # Eq. 10b
    k12_num_int <- fixed(-0.188)
    label("k12 numerator intercept (unitless)") # Eq. 10b
    e_age_k12_den <- fixed(0.64)
    label("k12 denominator coefficient on normalised age (unitless)") # Eq. 10b
    e_bmi_k12_den <- fixed(-0.0174)
    label("k12 denominator coefficient on normalised BMI (unitless)") # Eq. 10b
    k12_den_int <- fixed(-0.743)
    label("k12 denominator intercept (unitless)") # Eq. 10b

    # -- k13, central -> peripheral2 (Eqs. 10c male, 10d female) ---------
    # The male and female forms differ only in the sign of the two
    # age-linear terms, applied below through sexsign.
    e_age_k13_num_quad <- fixed(0.0058)
    label("k13 numerator coefficient on squared normalised age (1/s)") # Eqs. 10c, 10d
    e_age_k13_num <- fixed(0.00208)
    label("k13 numerator coefficient on normalised age, sign-switched by sex (1/s)") # Eqs. 10c, 10d
    k13_num_int <- fixed(0.0026)
    label("k13 numerator intercept (1/s)") # Eqs. 10c, 10d
    e_age_k13_den_quad <- fixed(2.75)
    label("k13 denominator coefficient on squared normalised age (unitless)") # Eqs. 10c, 10d
    e_age_k13_den <- fixed(0.985)
    label("k13 denominator coefficient on normalised age, sign-switched by sex (unitless)") # Eqs. 10c, 10d
    k13_den_int <- fixed(0.601)
    label("k13 denominator intercept (unitless)") # Eqs. 10c, 10d

    # -- k21, peripheral1 -> central (Eq. 10e) ---------------------------
    lk21 <- fixed(log(0.00218))
    label("k21 intercept of the covariate expression (1/s)") # Eq. 10e
    e_bmi_k21_quad <- fixed(0.00408)
    label("k21 coefficient on squared normalised BMI (1/s)") # Eq. 10e
    e_bmi_k21 <- fixed(-0.000816)
    label("k21 coefficient on normalised BMI (1/s)") # Eq. 10e
    e_bmi_wt_k21 <- fixed(-0.0057)
    label("k21 coefficient on the normalised BMI by normalised weight interaction (1/s)") # Eq. 10e

    # -- k31, peripheral2 -> central (Eqs. 10f male, 10g female) ---------
    lk31 <- fixed(log(4.52e-5))
    label("k31 intercept of the covariate expression (1/s)") # Eqs. 10f, 10g
    e_age_k31 <- fixed(1.92e-5)
    label("k31 coefficient on normalised age, sign-switched by sex (1/s)") # Eqs. 10f, 10g

    # -- V1, central volume (Eq. 10h) ------------------------------------
    e_age_vc <- fixed(0.0596)
    label("Central volume coefficient on normalised age (L)") # Eq. 10h
    e_wt_vc <- fixed(18.7)
    label("Central volume coefficient on normalised weight (L)") # Eq. 10h
    e_wt_vc_quad <- fixed(-13.7)
    label("Central volume coefficient on squared normalised weight (L)") # Eq. 10h
    e_age_wt_vc <- fixed(-3.5)
    label("Central volume coefficient on the normalised age by normalised weight interaction (L)") # Eq. 10h
    vc_int <- fixed(-0.0557)
    label("Central volume intercept (L)") # Eq. 10h

    # -- Residual error --------------------------------------------------
    # The paper fits a deterministic covariate model by minimising mean
    # MdALE and reports no residual-error model, so no magnitude can be
    # sourced. Held at zero rather than invented; see vignette Errata.
    propSd <- fixed(0)
    label("Proportional residual error (fraction; not published)")
  })

  model({
    # 1. Covariate normalisation.
    #    AGEMAX / WGTMAX / BMIMAX are the input-normalisation constants
    #    stated in the Results section immediately after Eq. 10. BMIMAX is
    #    printed as 52.8 kg/m^2 in the paper and as the full-precision
    #    dataset maximum 52.84713965 in the authors' code deposit; the
    #    latter is used here.
    agen <- AGE / 88
    wtn <- WT / 160
    bmin <- BMI / 52.84713965
    # +1 for male, -1 for female: selects between the printed male (Eqs.
    # 10c, 10f) and female (Eqs. 10d, 10g) forms.
    sexsign <- 1 - 2 * SEXF

    # 2. Covariate model. abs() on every output is the final-layer base
    #    expression g3(z3) = |z3| of Eq. 8, which guarantees a positive PK
    #    parameter; the paper prints the bars explicitly only on k12 and
    #    k21, where the argument can change sign over the covariate range.
    kel <- abs(exp(lkel) + e_wt_kel * wtn)
    k12 <- k12_norm *
      abs((e_age_k12_num * agen + e_bmi_k12_num * bmin + k12_num_int) /
        (e_age_k12_den * agen + e_bmi_k12_den * bmin + k12_den_int))
    k13 <- abs((e_age_k13_num_quad * agen^2 + sexsign * e_age_k13_num * agen + k13_num_int) /
      (e_age_k13_den_quad * agen^2 + sexsign * e_age_k13_den * agen + k13_den_int))
    k21 <- abs(exp(lk21) + e_bmi_k21_quad * bmin^2 + e_bmi_k21 * bmin +
      e_bmi_wt_k21 * bmin * wtn)
    k31 <- abs(exp(lk31) + sexsign * e_age_k31 * agen)
    vc <- abs(e_age_vc * agen + e_wt_vc * wtn + e_wt_vc_quad * wtn^2 +
      e_age_wt_vc * agen * wtn + vc_int)

    # 3. Mammillary three-compartment system, Eq. 1, written on amounts.
    #    Eq. 1 is stated for xi = Ai / V1; multiplying through by V1 gives
    #    the amount form below, so peripheral1 and peripheral2 hold true
    #    amounts and V2 = (k12 / k21) * V1, V3 = (k13 / k31) * V1 per Eq. 2.
    d/dt(central) <- -(kel + k12 + k13) * central + k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # 4. Observation. mg / L is the ug/mL used throughout the paper.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
