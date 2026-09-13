Tobin_2026_atomoxetine <- function() {
  description <- paste0(
    "One-compartment population pharmacokinetic model for oral atomoxetine in ",
    "86 children and adolescents (6-17 years) with ADHD across three studies ",
    "(159 participant-occasions, 1946 plasma concentrations), pooling ",
    "single-dose and steady-state occasions. Absorption is sequential ",
    "zero-order input into a depot over D1 = 0.75 h followed by first-order ",
    "transfer (ka = 27.35/h) into the central compartment, with linear ",
    "elimination. Apparent volume (Vc/F = 175.94 L) and apparent clearance ",
    "(CL/F = 39.03 L/h) are allometrically scaled on actual body weight ",
    "centred at 70 kg with fixed exponents of 1.0 and 0.75. Both apparent ",
    "parameters are divided by a relative-bioavailability factor Frel built ",
    "from CYP2D6 poor (3.02) and intermediate (1.33) metabolizer status and ",
    "CYP2C19 poor metabolizer status (2.32), all relative to normal ",
    "metabolizers; apparent clearance additionally carries an 81% reduction ",
    "in CYP2D6 poor metabolizers. CYP2D6 normal and ultrarapid metabolizers ",
    "share the reference category because ultrarapid effects did not explain ",
    "variability. Between-subject variability is large on the absorption ",
    "parameters (127% CV on D1, 183% CV on ka) and moderate on disposition ",
    "(64.5% CV on Vc/F, 75.5% CV on CL/F); residual error is combined ",
    "additive (0.81 ng/mL) plus proportional (21.4% CV). Parameters come from ",
    "a two-stage analysis in which each participant-occasion was treated as a ",
    "separate individual, so the reported between-subject variability pools ",
    "between-subject and between-occasion variability."
  )
  reference <- paste(
    "Tobin KVT, Gobburu J, Leeder JS, Pritchett A, Dunn A. (2026).",
    "Understanding Atomoxetine Exposure Variability in Children and",
    "Adolescents With ADHD Through Population Pharmacokinetics.",
    "The Journal of Clinical Pharmacology 66(4):e70168.",
    "doi:10.1002/jcph.70168.",
    sep = " "
  )
  vignette <- "Tobin_2026_atomoxetine"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ng/mL"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Tobin 2026 Results 'Structural Model'
  # paragraph ("a one-compartment distribution model with zero-order transit
  # into a depot compartment and first-order absorption into the central
  # compartment") and Methods 'Study Populations' (atomoxetine plasma
  # concentrations).
  compartmentData <- list(
    depot   = list(analyte = "atomoxetine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "atomoxetine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Actual body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Allometric scaling on both Vc/F (exponent 1.0) and CL/F (exponent ",
        "0.75), centred at a 70 kg typical adult. Tobin 2026 Equation (3) ",
        "defines the term as (wt_i / 70)^b_theta with '70 is the typical ",
        "adult body weight in kg'; the Discussion 'Impact of Covariates' ",
        "paragraph confirms 'centered around a 70-kg typical adult'. NOTE: ",
        "the printed Equation (5) and the Table 3 covariate-equation footer ",
        "both show (wt_i / 10)^0.75 for CL/F. That is a typographical error ",
        "-- a 10 kg centring is falsified by the paper's own NCA half-lives ",
        "(it predicts 3.7 h for CYP2D6 poor metabolizers against an observed ",
        "17.3 +/- 4.9 h, and 0.67 h for normal metabolizers against 2.5 +/- ",
        "0.6 h), whereas 70 kg reproduces all four phenotype half-lives. See ",
        "the vignette's Assumptions and deviations section. Cohort mean ",
        "54.5 +/- 26.4 kg (Tobin 2026 Table 1)."
      ),
      source_name        = "wt"
    ),
    CYP2D6_PM = list(
      description        = paste0(
        "1 = CYP2D6 poor-metabolizer phenotype (CPIC activity score 0), ",
        "0 = otherwise. Together with CYP2D6_IM = 0 the reference category ",
        "is the pooled CYP2D6 normal-plus-ultrarapid metabolizer group."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "Normal or ultrarapid metabolizer (CYP2D6_PM = 0 and CYP2D6_IM = 0)",
      notes              = paste0(
        "Time-fixed (germline genotype-derived phenotype). Tobin 2026 ",
        "Methods 'Genotyping and Phenotyping': phenotype assigned from CPIC ",
        "activity score, poor metabolizer = activity score 0. Carries two ",
        "distinct effects: a relative bioavailability of 3.02 versus normal ",
        "metabolizers, which divides BOTH Vc/F and CL/F (Tobin 2026 Table 3 ",
        "covariate equations), and a separate -81% effect on CL/F only ",
        "(E_CYP2D6,PM = -0.81). Cohort prevalence 6 of 86 participants ",
        "(7.0%), all male (Tobin 2026 Table 1). The CYP2D6_PM + CYP2D6_IM ",
        "pairing with a pooled normal-plus-ultrarapid reference follows the ",
        "Nguyen 2025 valbenazine precedent recorded in the CYP2D6_IM ",
        "register entry; ultrarapid metabolizers sit in the reference ",
        "because Tobin 2026 Results 'Covariate Model' reports that ",
        "'relative bioavailability or clearance effects of CYP2D6 UMs did ",
        "not explain variability in Vc/F nor CL/F'."
      ),
      source_name        = "CYP2D6 phenotype (poor metabolizer level)"
    ),
    CYP2D6_IM = list(
      description        = paste0(
        "1 = CYP2D6 intermediate-metabolizer phenotype (CPIC activity score ",
        "0.25-1), 0 = otherwise. Together with CYP2D6_PM = 0 the reference ",
        "category is the pooled CYP2D6 normal-plus-ultrarapid metabolizer ",
        "group."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "Normal or ultrarapid metabolizer (CYP2D6_PM = 0 and CYP2D6_IM = 0)",
      notes              = paste0(
        "Time-fixed (germline genotype-derived phenotype). Tobin 2026 ",
        "Methods 'Genotyping and Phenotyping': intermediate metabolizer = ",
        "CPIC activity score 0.25-1. Carries a relative bioavailability of ",
        "1.33 versus normal metabolizers, dividing BOTH Vc/F and CL/F. It ",
        "carries NO clearance effect: Tobin 2026 Results 'Covariate Model' ",
        "reports 'the effect of CYP2D6 IM on CL/F was investigated, but did ",
        "not explain the variability'. Cohort prevalence 41 of 86 ",
        "participants (47.7%), the largest phenotype group (Tobin 2026 ",
        "Table 1)."
      ),
      source_name        = "CYP2D6 phenotype (intermediate metabolizer level)"
    ),
    CYP2C19_PM = list(
      description        = paste0(
        "1 = CYP2C19 poor-metabolizer phenotype, 0 = any other CYP2C19 ",
        "phenotype (intermediate, normal, rapid, or ultrarapid)."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "Any non-poor CYP2C19 phenotype (CYP2C19_PM = 0)",
      notes              = paste0(
        "Time-fixed (germline genotype-derived phenotype). Carries a ",
        "relative bioavailability of 2.32 versus CYP2C19 normal ",
        "metabolizers, dividing BOTH Vc/F and CL/F, multiplicatively with ",
        "the CYP2D6 factor (Tobin 2026 Table 3 covariate equations). The ",
        "reference pools every non-poor CYP2C19 phenotype because only the ",
        "poor-metabolizer level entered the final model, so CYP2C19_IM is ",
        "deliberately not carried -- the same construction as ",
        "Marathe_2023_belzutifan.R. Cohort prevalence 2 of 86 participants ",
        "(2.3%, 4 dosing occasions), both of whom were also CYP2D6 ",
        "intermediate metabolizers, so the two relative-bioavailability ",
        "factors multiply for those subjects. Tobin 2026 Discussion ",
        "'Impact of Covariates' cautions that this effect rests on two ",
        "individuals and 'should be interpreted with caution'."
      ),
      source_name        = "CYP2C19 phenotype (poor metabolizer level)"
    )
  )

  # Covariates that Tobin 2026 screened but did not retain in the final model.
  # Documentation only -- none is referenced in model().
  covariatesDataExcluded <- list(
    CYP2D6_UM = list(
      description = "1 = CYP2D6 ultrarapid-metabolizer phenotype (CPIC activity score > 2), 0 = otherwise",
      units       = "(binary)",
      type        = "binary",
      notes       = paste0(
        "Screened and rejected. Tobin 2026 Results 'Covariate Model': ",
        "'Relative bioavailability or clearance effects of CYP2D6 UMs did ",
        "not explain variability in Vc/F nor CL/F.' Ultrarapid metabolizers ",
        "(4 of 86 participants) therefore sit in the model's reference ",
        "category alongside normal metabolizers."
      )
    ),
    CYP2C19_IM = list(
      description = "1 = CYP2C19 intermediate-metabolizer phenotype, 0 = otherwise",
      units       = "(binary)",
      type        = "binary",
      notes       = paste0(
        "Screened as part of the CYP2C19 activity covariate exploration ",
        "(Tobin 2026 Methods 'Covariate Model Development': 'the effect of ",
        "CYP2C19 enzyme activity on clearance compared to CYP2C19 NMs was ",
        "explored similarly to Equation (5)'). Only the poor-metabolizer ",
        "level was retained; intermediate metabolizers (21 of 86, 24%) sit ",
        "in the reference category."
      )
    ),
    BMI = list(
      description = "Body-mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = paste0(
        "Obesity (BMI at or above the 95th percentile for age) was screened ",
        "against the empirical Bayes estimates of Vc/F and CL/F (Tobin 2026 ",
        "Results 'Covariate Model', Figure S4) but not retained as a ",
        "separate term: allometric scaling on actual body weight 'adequately ",
        "explained the BSV in apparent volume and clearance by body weight ",
        "and obesity' (Discussion 'Impact of Covariates'). Cohort mean ",
        "21.8 +/- 6.4 kg/m^2, 22 of 86 participants (26%) obese."
      )
    ),
    RACE_ASIAN = list(
      description = "1 = self-identified Asian race, 0 = otherwise",
      units       = "(binary)",
      type        = "binary",
      notes       = paste0(
        "Screened and rejected for lack of data. Tobin 2026 Results ",
        "'Covariate Model': 'When stratifying by ethnicity, Asians showed ",
        "relatively lower Ka, Vc/F, and CL/F but were not included as a ",
        "covariate due to having only one Asian participant (with two PK ",
        "visits) in the dataset.'"
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 86L,
    n_studies      = 3L,
    age_range      = "6-17 years",
    age_median     = "12.6 years (mean; SD 3.2)",
    weight_range   = "mean 54.5 kg (SD 26.4); range not reported",
    weight_median  = "54.5 kg (mean; SD 26.4)",
    sex_female_pct = 18.6,
    race_ethnicity = c(
      White = 51, Black = 33, `Mixed race` = 13, Asian = 1.2,
      `Hispanic/Latino` = 1.2, `Native Hawaiian/Pacific Islander` = 1.2
    ),
    disease_state  = "attention-deficit/hyperactivity disorder (ADHD), genotyped for CYP2D6 and CYP2C19",
    dose_range     = paste0(
      "Study A: single oral 0.5 mg/kg weight-based dose. Studies B and C: ",
      "single oral doses individualised by CYP2D6 phenotype to target a ",
      "steady-state Cmax of 400 ng/mL, followed in Study C by consistent ",
      "once-daily oral dosing (twice daily for two participants) with ",
      "steady-state PK sessions at weeks 6 +/- 2 and 18 +/- 2"
    ),
    regions        = "United States (Children's Mercy Kansas City)",
    n_observations = "1946 atomoxetine plasma concentrations over 159 participant-occasions",
    cyp2d6_phenotype = c(Intermediate = 41L, Normal = 35L, Poor = 6L, `Ultra rapid` = 4L),
    cyp2c19_phenotype = c(Normal = 31L, Rapid = 26L, Intermediate = 21L, `Ultra rapid` = 6L, Poor = 2L),
    notes          = paste0(
      "Baseline demographics are Tobin 2026 Table 1. Three pooled studies, ",
      "all approved by the Children's Mercy Kansas City IRB: Study A ",
      "(single 0.5 mg/kg dose, reported by Brown et al), Study B (single ",
      "dose model-validation with a repeat visit 15-386 days later to ",
      "assess day-to-day variability), and Study C (the GOLDILOKs ",
      "exposure-escalation study, NCT03154359, grant P50HD090258: first ",
      "dose plus up to two steady-state sessions). 33 participants had one ",
      "single-dose visit, 21 had two single-dose visits, 12 had one ",
      "single-dose and one steady-state visit, and 20 had one single-dose ",
      "and two steady-state visits. Sampling ran to 72 h post-dose for poor ",
      "metabolizers, 24 h for intermediate and normal metabolizers, and ",
      "12 h for ultrarapid metabolizers. Modelling used Pumas v2.5.1, not ",
      "NONMEM. Because absorption was unpredictably rapid or delayed both ",
      "between and within subjects, a one-stage analysis failed to ",
      "converge and each participant-occasion was treated as a separate ",
      "individual in a two-stage approach; the reported between-subject ",
      "variability therefore pools true between-subject with ",
      "between-occasion variability, and the authors note that two-stage ",
      "approaches tend to overestimate it."
    )
  )

  ini({
    # ========================================================================
    # Structural PK. Tobin 2026 Table 3 ('Covariate Model Parameter and
    # Uncertainty Estimates'), final covariate model. Typical values are for
    # a 70 kg CYP2D6 normal metabolizer who is not a CYP2C19 poor metabolizer
    # (Frel = 1, E_CYP2D6,PM = 0).
    # ========================================================================
    ld1 <- log(0.75)
    label("Zero-order duration of input into the depot, D1 (h)")
    # Tobin 2026 Table 3: Dur = 0.75 h (95% CI 0.64-0.88; RSE 7.76%). Also
    # quoted in Results 'Structural Model' ('the estimated zero-order
    # duration into the depot was 0.75 h').

    lka <- log(27.35)
    label("First-order absorption rate constant from depot to central, ka (1/h)")
    # Tobin 2026 Table 3: ka = 27.35 1/h (95% CI 22.44-32.77; RSE 10.1%).
    # Discussion 'Structural Model' notes this parameter was bounded during
    # estimation because profiles with near-instantaneous absorption pushed
    # it arbitrarily high, and that 'first-order rate constants above 50/h
    # had minimal visual differences in atomoxetine concentrations'.

    lvc <- log(175.94)
    label("Apparent central volume of distribution, Vc/F (L) at 70 kg for the reference phenotype")
    # Tobin 2026 Table 3: Vc/F = 175.94 L (95% CI 148.48-217.50; RSE 9.46%).
    # The base structural model (no covariates) gave 98.96 L, which the
    # Discussion compares with the label-reported 110.25 L apparent volume.

    lcl <- log(39.03)
    label("Apparent clearance, CL/F (L/h) at 70 kg for the reference phenotype")
    # Tobin 2026 Table 3: CL/F = 39.03 L/h (95% CI 33.90-46.33; RSE 8.14%).
    # The base structural model gave 23.45 L/h.

    # ========================================================================
    # Allometric scaling on actual body weight, centred at 70 kg.
    # Tobin 2026 Equation (3) and Table 3 covariate equations. Both exponents
    # are conventional values carried without uncertainty, so both are fixed:
    # Equation (3) states b_theta 'was 1.0 when applied to volume and 0.75
    # when applied to clearance', and the Discussion calls them 'conventional
    # allometric scaling exponents'. Neither appears in Table 3's estimate
    # column.
    # ========================================================================
    e_wt_vc <- fixed(1.0)
    label("Allometric exponent on Vc/F for body weight centred at 70 kg (unitless)")
    # Tobin 2026 Equation (3) and Equation (4); Results 'Covariate Model'
    # ('a typical weight of 70 kg on Vc/F (exponent: 1.0)').

    e_wt_cl <- fixed(0.75)
    label("Allometric exponent on CL/F for body weight centred at 70 kg (unitless)")
    # Tobin 2026 Equation (3); Results 'Covariate Model' ('and CL/F
    # (exponent: 0.75)'). The centring weight is 70 kg, not the 10 kg shown
    # in the printed Equation (5) and Table 3 footer -- see the WT entry in
    # covariateData and the vignette's Assumptions and deviations section.

    # ========================================================================
    # Relative bioavailability (Frel) by metabolizer phenotype, relative to
    # CYP2D6 / CYP2C19 normal metabolizers (Frel,NM = 1). Tobin 2026
    # Equation (4), Equation (5), and the Table 3 covariate equations apply
    # Frel as a DIVISOR of both Vc/F and CL/F, and the CYP2D6 and CYP2C19
    # factors multiply. Held on the log scale so the factors compose
    # additively inside model(); the linear-scale value is in each comment.
    # ========================================================================
    e_cyp2d6_pm_vc_cl <- log(3.02)
    label("Log relative bioavailability of CYP2D6 poor metabolizers, dividing Vc/F and CL/F (unitless)")
    # Tobin 2026 Table 3: Frel CYP2D6,PM = 3.02 (95% CI 2.38-3.77; RSE
    # 11.67%). Results 'Covariate Model': 'The relative bioavailabilities of
    # CYP2D6 PMs and IMs were estimated as 3.02 and 1.33, respectively,
    # compared to NMs'. log(3.02) = 1.1053.

    e_cyp2d6_im_vc_cl <- log(1.33)
    label("Log relative bioavailability of CYP2D6 intermediate metabolizers, dividing Vc/F and CL/F (unitless)")
    # Tobin 2026 Table 3: relative bioavailability of CYP2D6 IMs = 1.33
    # (95% CI 0.95-2.01; RSE 20.97%). The Table 3 parameter label is
    # misprinted as 'Frel CYP2C6,IM'; its Description column reads
    # 'Relative Bioavailability of CYP2D6 IMs', and the Results text
    # confirms 1.33 for CYP2D6 IMs. log(1.33) = 0.2852. The 95% CI includes
    # 1, so this effect is imprecisely estimated; the point estimate is
    # retained as published.

    e_cyp2c19_pm_vc_cl <- log(2.32)
    label("Log relative bioavailability of CYP2C19 poor metabolizers, dividing Vc/F and CL/F (unitless)")
    # Tobin 2026 Table 3: Frel CYP2C19,PM = 2.32 (95% CI 1.65-3.24; RSE
    # 18.27%). Results 'Covariate Model': 'the relative bioavailability of
    # CYP2C19 PMs was estimated as 2.32'. log(2.32) = 0.8416.

    # ========================================================================
    # CYP2D6 poor-metabolizer effect on apparent clearance only. Applied as
    # the multiplicative factor (1 + E_CYP2D6) in Tobin 2026 Equation (5),
    # with E_CYP2D6,NM = 0.
    # ========================================================================
    e_cyp2d6_pm_cl <- -0.81
    label("Proportional effect of CYP2D6 poor-metabolizer status on CL/F, applied as (1 + effect) (unitless)")
    # Tobin 2026 Table 3: E_CYP2D6,PM = -0.81 (95% CI -0.84 to -0.78; RSE
    # 1.95%). Discussion 'Impact of Covariates': 'This effect was estimated
    # as an 81% reduction in clearance for CYP2D6 PMs compared to CYP2D6
    # NMs.' The factor is 1 - 0.81 = 0.19. Note the Results text renders this
    # as '-0.81%'; the Discussion, the Table 3 value and the CI all agree it
    # is a proportion (-81%), not -0.81%.

    # ========================================================================
    # Between-subject variability. Tobin 2026 Table 3 reports these as CV%
    # (the column header reads 'CV%' and Results 'Covariate Model' quotes
    # the same numbers as percentages: 'reduced the Vc/F BSV from 81.6% to
    # 64.5% and reduced the CL/F BSV from 92.4% to 75.5%'). Converted to the
    # log-normal variance scale by omega^2 = log(1 + CV^2) to match
    # Equation (2), theta_i = theta_tv * exp(eta).
    #
    # Because the analysis treated each participant-occasion as a separate
    # individual (Methods 'Data Analysis'), these variances pool
    # between-subject with between-occasion variability.
    # ========================================================================
    etald1 ~ log(1 + 1.2742^2)
    # Tobin 2026 Table 3, row 'BSV on Zero-order Duration into Depot (CV%)'
    # = 127.42 (95% CI 100.64-150.30; RSE 13.79%; shrinkage 0.00%).
    # omega^2 = log(1 + 1.2742^2) = 0.9646.

    etalka ~ log(1 + 1.827^2)
    # Tobin 2026 Table 3, row 'BSV on Absorption Rate Constant (CV%)' =
    # 182.7 (95% CI 146.09-244.82; RSE 14.03%; shrinkage 0.00%).
    # omega^2 = log(1 + 1.827^2) = 1.4674.

    etalvc ~ log(1 + 0.6452^2)
    # Tobin 2026 Table 3, row 'BSV on Apparent Central Volume (CV%)' =
    # 64.52 (95% CI 45.08-76.48; RSE 21.32%; shrinkage 15.88%).
    # omega^2 = log(1 + 0.6452^2) = 0.3480.

    etalcl ~ log(1 + 0.7548^2)
    # Tobin 2026 Table 3, row 'BSV on Apparent Central Clearance (CV%)' =
    # 75.48 (95% CI 45.82-93.55; RSE 26.29%; shrinkage 15.52%).
    # omega^2 = log(1 + 0.7548^2) = 0.4509.

    # ========================================================================
    # Residual unexplained variability. Combined additive plus proportional,
    # Tobin 2026 Equation (1): OBS = PRED + PRED * eps_prop + eps_add.
    # ========================================================================
    addSd <- 0.81
    label("Additive residual standard deviation (ng/mL)")
    # Tobin 2026 Table 3, row 'Additive RUV (ng/mL)' = 0.81 (95% CI
    # 0.38-1.21; RSE 50.78%). Results 'Structural Model' quotes '0.81 ng/mL'
    # directly, confirming the standard-deviation scale despite the sigma^2
    # symbol in the row label. The Discussion attributes the high RSE to
    # 'one subject with an estimated 150 ng/mL additive error'.

    propSd <- 0.2139
    label("Proportional residual standard deviation (fraction)")
    # Tobin 2026 Table 3, row 'Proportional RUV (CV%)' = 21.39 (95% CI
    # 19.36-23.63; RSE 10.17%). Results 'Structural Model' quotes '21 CV%',
    # confirming the CV / standard-deviation scale.
  })

  model({
    # ----------------------------------------------------------------------
    # 1. Relative bioavailability factor. Tobin 2026 Table 3 covariate
    #    equations carry the CYP2D6 and CYP2C19 factors as separate
    #    multiplicative terms, (1 / Frel_CYP2D6) * (1 / Frel_CYP2C19), so on
    #    the log scale they add. CYP2D6 poor and intermediate are mutually
    #    exclusive; the two CYP2C19 poor metabolizers in the cohort were both
    #    CYP2D6 intermediate metabolizers, so the product term is exercised.
    #    Normal and ultrarapid CYP2D6 metabolizers give frel = 1.
    # ----------------------------------------------------------------------
    frel <- exp(e_cyp2d6_pm_vc_cl * CYP2D6_PM +
                  e_cyp2d6_im_vc_cl * CYP2D6_IM +
                  e_cyp2c19_pm_vc_cl * CYP2C19_PM)

    # ----------------------------------------------------------------------
    # 2. Individual parameters.
    #    Equation (4): Vc/F_i = Vc/F_tv * (wt/70)^1.0 * (1/Frel) * exp(eta).
    #    Equation (5): CL/F_i = CL/F_tv * (wt/70)^0.75 * (1/Frel)
    #                           * (1 + E_CYP2D6) * exp(eta).
    #    Both are APPARENT parameters (oral dosing only, so F is not
    #    identifiable); Frel divides them rather than multiplying an explicit
    #    bioavailability term, exactly as printed.
    # ----------------------------------------------------------------------
    ka <- exp(lka + etalka)
    d1 <- exp(ld1 + etald1)
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc / frel
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl / frel *
      (1 + e_cyp2d6_pm_cl * CYP2D6_PM)

    # ----------------------------------------------------------------------
    # 3. Micro-constant.
    # ----------------------------------------------------------------------
    kel <- cl / vc

    # ----------------------------------------------------------------------
    # 4. ODE system. One-compartment disposition with sequential zero-order
    #    input into the depot followed by first-order absorption into the
    #    central compartment and linear elimination (Tobin 2026 Results
    #    'Structural Model').
    # ----------------------------------------------------------------------
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ----------------------------------------------------------------------
    # 5. Zero-order input duration. The dose is delivered uniformly into the
    #    depot over [0, d1] and then drains by ka. Dosing records must set
    #    rate = -2 (modelled duration); without it rxode2 silently treats the
    #    dose as a bolus and this line has no effect.
    # ----------------------------------------------------------------------
    dur(depot) <- d1

    # ----------------------------------------------------------------------
    # 6. Observation. Doses are in mg and vc is in L, so central / vc is in
    #    mg/L = ug/mL; multiply by 1000 to report ng/mL, the units used
    #    throughout Tobin 2026.
    # ----------------------------------------------------------------------
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
