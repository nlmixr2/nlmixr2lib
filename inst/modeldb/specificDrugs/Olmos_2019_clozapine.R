# Simultaneous one-compartment parent-plus-metabolite population PK model for
# oral clozapine and its active metabolite norclozapine, with a relative
# bioavailability term comparing two Uruguayan clozapine brands (Olmos 2019,
# BioMed Res Int 2019:3163502; doi:10.1155/2019/3163502).

Olmos_2019_clozapine <- function() {
  description <- paste(
    "Simultaneous one-compartment parent-plus-metabolite population PK model",
    "for oral clozapine (CZP) and its active metabolite norclozapine (NCZP)",
    "in 98 Uruguayan adult inpatients (76 male, 22 female) with DSM-IV",
    "schizophrenia, fit to 171 steady-state morning trough observations per",
    "analyte (Olmos 2019). First-order absorption (ka fixed at 1.24 1/h from",
    "Jerling 1996) into a clozapine central compartment with first-order",
    "elimination; complete (f = 1) conversion of clozapine to norclozapine is",
    "assumed, so the whole clozapine elimination flux feeds a second",
    "one-compartment metabolite compartment after a molecular-weight",
    "correction. Both apparent volumes of distribution are fixed from Golden",
    "and Honigfeld (750 L clozapine, 1860 L norclozapine at 70 kg) and scale",
    "linearly with body weight; both apparent clearances scale with body",
    "weight to the fixed allometric 0.75 power. Smoking status was the only",
    "covariate retained in the final model: clozapine apparent clearance is",
    "estimated separately in nonsmokers (28.1 L/h) and smokers (36.5 L/h).",
    "The study switched patients from the brand-name product (Leponex) to a",
    "similar product (Luverina), so a relative bioavailability of 0.892 for",
    "Luverina versus the Leponex reference (whose F is the fixed 1 anchor) is",
    "estimated together with its own between-subject variability. Clozapine",
    "and norclozapine apparent clearances carry correlated between-subject",
    "variability; residual error is proportional and separate per analyte.",
    sep = " "
  )
  reference <- paste(
    "Olmos I, Ibarra M, Vazquez M, Maldonado C, Fagiolino P, Giachetto G",
    "(2019). Population Pharmacokinetics of Clozapine and Norclozapine and",
    "Switchability Assessment between Brands in Uruguayan Patients with",
    "Schizophrenia. BioMed Research International 2019:3163502.",
    "doi:10.1155/2019/3163502.",
    sep = " "
  )
  vignette <- "Olmos_2019_clozapine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Units follow the units block (mg amounts, ng/mL
  # concentrations after the 1000x conversion in model()).
  compartmentData <- list(
    depot = list(analyte = "clozapine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "clozapine", units = "mg", specimen = "plasma", verified = TRUE),
    central_norcloz = list(analyte = "norclozapine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Applied to both apparent volumes of distribution linearly and to",
        "both apparent clearances as a fixed 0.75 power, each centred on a",
        "70 kg reference subject. Methods Eq. (1) 'Vi = V * (BWi / 70)'",
        "where 'Vi is the apparent volume of distribution for the i-th",
        "subject, BWi its body weight in kilograms, and V the apparent",
        "volume of distribution for a 70-kg body weight subject'; Methods",
        "Eq. (2) 'CLapi = CLap * (BWi / 70)^0.75'. Eq. (2) is written for a",
        "generic CLap and the paper reports a single allometric statement",
        "('the effect of body weight on clearance'), so the same fixed 0.75",
        "exponent is applied to the clozapine and the norclozapine apparent",
        "clearance. Cohort body weight (Table 1) was a median 78 kg, range",
        "48-137 kg."
      ),
      source_name = "BW (body weight, kg)"
    ),
    SMOKE = list(
      description = "Current-smoker status indicator",
      units = "(binary)",
      type = "binary",
      reference_category = paste(
        "0 (nonsmoker). Olmos 2019 estimates a separate typical clozapine",
        "apparent clearance in each stratum rather than a reference value",
        "plus a fractional offset, so neither level is a reference in the",
        "usual covariate-coefficient sense; the nonsmoker estimate is the",
        "one the paper labels the basal value in the Discussion."
      ),
      notes = paste(
        "Self-reported, dichotomised into smokers and nonsmokers with no",
        "assessment of the magnitude of smoking; the authors list this as a",
        "study limitation (Discussion: 'smoking status was assessed using",
        "patient self-reporting. We dichotomized patients into smokers and",
        "nonsmokers but did not assess the magnitude of smoking'). Cohort",
        "split (Table 1): 46 smokers, 52 nonsmokers, and smokers were 45%",
        "within each sex. The only covariate retained in the final model,",
        "and only on clozapine apparent clearance -- norclozapine apparent",
        "clearance was not affected (Discussion: 'in our study only CLap CZP",
        "seemed to be affected')."
      ),
      source_name = "smoking status"
    ),
    FORM_CZP_LUVERINA = list(
      description = paste(
        "Clozapine drug-product indicator: 1 = Luverina (Celsius",
        "Laboratories, the similar product), 0 = Leponex (Novartis",
        "Laboratories, the brand-name reference product). Both are 100 mg",
        "oral clozapine tablets."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = paste(
        "0 (Leponex, the brand-name reference product, whose oral",
        "bioavailability F is fixed to 1 and is what makes the Luverina",
        "relative bioavailability identifiable)."
      ),
      notes = paste(
        "Time-varying within subject: the study is a sequential switch, so",
        "73 of 98 patients contributed one trough on Leponex and one on",
        "Luverina after two months on the new brand, and the indicator must",
        "be set per dose record. Methods: 'F was fixed to 1 for Leponex and",
        "evaluated for Luverina as an estimate of the relative CZP",
        "bioavailability between both drug products.' Because clozapine and",
        "norclozapine apparent clearances are both defined with F in the",
        "denominator (CLapCZP = CLCZP/F, CLapNCZP = CLNCZP/F/f), the",
        "bioavailability term shifts both analytes together while each",
        "apparent clearance shifts only its own analyte; that is what makes",
        "the bioavailability random effect identifiable alongside the two",
        "clearance random effects."
      ),
      source_name = "CZP formulation (Leponex / Luverina)"
    )
  )

  # Covariates that Olmos 2019 screened but did not retain in the final
  # model. Documented here so the covariate screen is preserved without
  # declaring a covariate that model() never references. Methods: 'Covariate
  # search was performed for CLapCZP, CLapNCZP, and F, evaluating the effect
  # of sex, smoking status, CZP formulation, beginning of treatment, caffeine
  # consumption, and concomitant treatments: valproic acid, benzodiazepines,
  # antidepressants, antipsychotics, antidiabetics, and oral hypoglycemic
  # drugs. The effect of daily CZP dose on both apparent clearances was also
  # assessed.' Forward inclusion used dOFV = 3.84 (p < 0.05, 1 df) and
  # backward elimination the stricter dOFV = 10.83 (p < 0.001, 1 df).
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Biological sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Screened on both apparent clearances and not retained. Results:",
        "'Sex did not produce a significant impact in these data: an",
        "increase in the AIC was observed after including it as a covariate",
        "on CZP and NCZP apparent clearance, and no significant differences",
        "between male and female estimates were obtained. This covariate was",
        "reevaluated after smoking factor was included to discard a masking",
        "effect, obtaining similar results.' Cohort: 22 women, 76 men",
        "(Table 1). No point estimate is published, so the effect cannot be",
        "encoded even optionally."
      )
    ),
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      notes = paste(
        "Listed in the study objective as a covariate of interest ('with a",
        "focus on covariates such as cigarette smoking, age, sex, caffeine",
        "consumption, brands available of CZP, and comedication') and",
        "tabulated for the cohort (Table 1: median 39 years, range 20-68),",
        "but no age effect is reported in the covariate-search Methods list",
        "or the Results, and none is retained."
      )
    ),
    CAFFEINE_USE = list(
      description = "Habitual caffeine (coffee) consumption indicator",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Screened on compound elimination and discarded. Results: 'The",
        "inclusion of caffeine intake as a covariate on compound elimination",
        "was also discarded, after estimating a very small (less than 1%",
        "increase) impact on CZP and NCZP CLap.' 75 of 98 patients were",
        "consumers (Table 1), similarly distributed across smoking strata.",
        "The sub-1% figure is a magnitude statement, not a published point",
        "estimate, so nothing is encoded."
      )
    ),
    CONMED_VALPROIC_ACID = list(
      description = "Concomitant valproic acid indicator",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Screened and not retained. Results: 'Regarding drug-drug",
        "interactions, the most influential coadministered drug was valproic",
        "acid (VPA), increasing NCZP CLap by 10%, a relationship that was",
        "not retained for the final model because it was not found to be",
        "statistically significant.' 37 of 98 patients received VPA at",
        "approximately 400 mg/day (Table 1, Discussion). The 10% figure",
        "failed the backward-elimination criterion and is therefore not",
        "encoded."
      )
    ),
    CONMED_BENZODIAZEPINE = list(
      description = "Concomitant benzodiazepine indicator",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Named in the Methods covariate-search list and tabulated (41 of 98",
        "patients, Table 1); no effect is reported in Results and none is",
        "retained."
      )
    ),
    CONMED_ANTIDEPRESSANT = list(
      description = "Concomitant antidepressant (sertraline or escitalopram) indicator",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Screened and not retained. Results: 'Among all subjects in the",
        "study, 28 also received antidepressants, sertraline, or",
        "escitalopram. However, no impact on CZP and/or NCZP CLap was",
        "observed with this medication.' The Discussion explains the",
        "mechanism: neither agent inhibits CYP1A2 or CYP3A4, unlike",
        "fluvoxamine."
      )
    ),
    DOSE_CZP_MGD = list(
      description = "Total daily oral clozapine dose",
      units = "mg/day",
      type = "continuous",
      notes = paste(
        "Screened on both apparent clearances as a test for dose-dependent",
        "(nonlinear) kinetics and not retained. Results: 'No correlation was",
        "observed between CZP and NCZP apparent clearances with the daily",
        "dose.' Cohort daily dose was a median 350 mg, range 150-700 mg",
        "(Table 1). The dose itself is of course carried on the dose",
        "records; this entry documents only the screened covariate effect."
      )
    ),
    TIME_TRT_START = list(
      description = "Beginning of treatment (time since clozapine therapy was started)",
      units = "not specified",
      type = "continuous",
      notes = paste(
        "Named in the Methods covariate-search list as 'beginning of",
        "treatment'. All patients had been on Leponex for more than a year",
        "at enrolment, so the cohort carries little contrast on this axis.",
        "No effect is reported in Results and none is retained."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 98L,
    n_observations = 171L,
    n_studies = 1L,
    age_range = "20-68 years (median 39; Table 1)",
    age_median = "39 years",
    weight_range = "48-137 kg (median 78; Table 1)",
    weight_median = "78 kg",
    bmi_range = "15-43 kg/m^2 (median 26; Table 1)",
    sex_female_pct = 22.4,
    race_ethnicity = c(White = 100),
    disease_state = paste(
      "DSM-IV-diagnosed schizophrenia, inpatients of Hospital Vilardebo,",
      "Montevideo, Uruguay. All patients had been treated with brand-name",
      "clozapine (Leponex, Novartis) for more than one year before the",
      "hospital's purchasing switched to the similar product (Luverina,",
      "Celsius); patients were on Luverina for two months before the second",
      "blood sampling. The final dataset describes the cohort as Caucasian."
    ),
    dose_range = paste(
      "Oral clozapine 150-700 mg/day (median 350 mg/day; Table 1),",
      "administered twice a day with each brand. 68 of the 73 patients who",
      "completed both periods (93%) kept the same dosage regimen after the",
      "brand change."
    ),
    smoke_strata = "46 smokers (37 male), 52 nonsmokers (39 male) (Table 1)",
    regions = "Uruguay (single centre, Hospital Vilardebo, Montevideo)",
    notes = paste(
      "Very sparse therapeutic-drug-monitoring design: a single morning",
      "predose (trough) sample per subject per treatment period, taken at",
      "steady state under unchanged comedication. 171 trough observations",
      "were recorded for each of clozapine and norclozapine, of which 146",
      "came from the 73 patients who completed both periods; 25 patients",
      "contributed one period only (17 Luverina, 8 Leponex). Because only",
      "one observation per subject per period was available, interoccasion",
      "variability was not identifiable and was not included, and Cmax,ss /",
      "Tmax,ss could not be estimated (Discussion limitation).",
      "Concentrations were measured by HPLC with UV detection at 230 nm",
      "using medazepam as internal standard; the assay was linear over",
      "54.8-1086 ng/mL (clozapine) and 72.3-1085 ng/mL (norclozapine).",
      "All clozapine observations were above the LLOQ; left-censored",
      "norclozapine observations were under 4% of the total and were",
      "included as such. Mean (SD) concentrations were 421 (262) ng/mL for",
      "clozapine and 275 (180) ng/mL for norclozapine (Table 1).",
      "Estimation used NONMEM 7.4 with Pirana-PsN-Xpose; the final model was",
      "evaluated by numerical predictive check and NPDE, and parameter",
      "precision came from a 200-sample nonparametric bootstrap."
    )
  )

  ini({
    # ---- Absorption ----------------------------------------------------
    # Methods 2.3: 'CZP first-order constant rate for absorption (ka) was
    # fixed to a value of 1.24 h-1 as estimated by Jerling et al. [31].'
    # Not estimated here, so no uncertainty is reported for it.
    lka <- fixed(log(1.24))
    label("Clozapine absorption rate constant ka (1/h), taken from Jerling 1996")

    # ---- Clozapine apparent clearance, by smoking stratum --------------
    # Table 2 reports two separate typical values, each with its own RSE
    # and bootstrap interval, rather than a reference value plus a
    # fractional covariate coefficient, so both strata carry an explicit
    # suffix (symmetric stratum-suffix convention).
    lcl_nonsmoke <- log(28.1)
    label("Clozapine apparent elimination clearance CLap CZP in nonsmokers at 70 kg (L/h)")
    lcl_smoke <- log(36.5)
    label("Clozapine apparent elimination clearance CLap CZP in smokers at 70 kg (L/h)")

    # ---- Norclozapine apparent clearance -------------------------------
    # Table 2 row 3 is printed with the parameter name 'CLap CZP' but the
    # Description column reads 'Norclozapine apparent elimination
    # clearance'; this is a typographical slip for CLap NCZP (see the
    # vignette Errata). Smoking was NOT retained on this parameter.
    lcl_norcloz <- log(53.6)
    label("Norclozapine apparent elimination clearance CLap NCZP at 70 kg (L/h)")

    # ---- Apparent volumes of distribution (both fixed) -----------------
    # Methods 2.3: 'fixing the apparent volumes of distribution (V/F) of
    # 750 L and 1860 L for CZP and NCZP, respectively. These mean values
    # were estimated by Golden and Honigfeld after conducting a multiple
    # dosing bioequivalence study with extensive sampling in 30 patients
    # with schizophrenia and were attributed in this work to a 70-kg body
    # weight subject under a proportional centered model.'
    lvc <- fixed(log(750))
    label("Clozapine apparent volume of distribution V/F at 70 kg (L), taken from Golden and Honigfeld")
    lvc_norcloz <- fixed(log(1860))
    label("Norclozapine apparent volume of distribution V/F at 70 kg (L), taken from Golden and Honigfeld")

    # ---- Body-weight scaling (all fixed) -------------------------------
    # Methods Eq. (2): CLapi = CLap * (BWi/70)^0.75, with the text 'a power
    # model was included evaluating different coefficients and finally
    # fixing this value to the allometric standard of 0.75'.
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent on (WT/70) for clozapine apparent clearance (unitless)")
    e_wt_cl_norcloz <- fixed(0.75)
    label("Allometric exponent on (WT/70) for norclozapine apparent clearance (unitless)")

    # Methods Eq. (1): Vi = V * (BWi/70), i.e. an exponent of exactly 1
    # ('a proportional centered model'), not estimated.
    e_wt_vc <- fixed(1)
    label("Exponent on (WT/70) for clozapine apparent volume of distribution (unitless)")
    e_wt_vc_norcloz <- fixed(1)
    label("Exponent on (WT/70) for norclozapine apparent volume of distribution (unitless)")

    # ---- Relative bioavailability of the two clozapine products --------
    # Methods 2.3: 'F was fixed to 1 for Leponex and evaluated for Luverina
    # as an estimate of the relative CZP bioavailability between both drug
    # products.' Both strata carry an explicit suffix; the Leponex value is
    # the fixed identifiability anchor.
    lfdepot_leponex <- fixed(log(1))
    label("Oral bioavailability of the Leponex reference product (unitless), the identifiability anchor")
    lfdepot_luverina <- log(0.892)
    label("Bioavailability of Luverina relative to Leponex (unitless)")

    # ---- Between-subject variability -----------------------------------
    # Table 2 'Between-subject CV' section reports CV percentages and the
    # clozapine/norclozapine clearance association as a percentage. For the
    # exponential (log-normal) model of Methods Eq. (3), omega^2 =
    # log(CV^2 + 1):
    #   CLap CZP  CV = 43.3% -> omega^2 = log(1 + 0.433^2) = 0.171841
    #   CLap NCZP CV = 49.9% -> omega^2 = log(1 + 0.499^2) = 0.222344
    #   F         CV = 43.6% -> omega^2 = log(1 + 0.436^2) = 0.174034
    # The 'cov CLap CZP - CLap NCZP (%)' row of 55.7% must be read as a
    # CORRELATION coefficient of 0.557, not as a covariance: every
    # covariance reading of 55.7% implies a correlation above 1 and is
    # therefore inadmissible (0.557 / sqrt(0.171841 * 0.222344) = 2.85;
    # log(1 + 0.557^2) / sqrt(...) = 1.38). At r = 0.557 the covariance is
    #   0.557 * sqrt(0.171841 * 0.222344) = 0.108876
    # and the 2x2 block is positive definite (determinant 0.0264).
    etalcl + etalcl_norcloz ~ c(
      0.171841,
      0.108876, 0.222344
    )
    etalfdepot_luverina ~ 0.174034

    # ---- Residual unexplained variability -------------------------------
    # Methods Eq. (4): Cik = Cpred * (1 + eps_ik), a pure proportional
    # error with a separate magnitude per analyte.
    propSd <- 0.0954
    label("Proportional residual error for clozapine (fraction)")
    propSd_norcloz <- 0.153
    label("Proportional residual error for norclozapine (fraction)")
  })

  model({
    # Molecular weights used to convert the clozapine elimination flux
    # (mg of clozapine per hour) into the norclozapine formation flux
    # (mg of norclozapine per hour). Methods 2.3: "Complete conversion of
    # CZP into NCZP was assumed and a factor was included in NCZP formation
    # to account for the molecular weight differences." Olmos 2019 does not
    # print the two molecular weights, so they are taken from the compound
    # formulae (clozapine C18H19ClN4, norclozapine C17H17ClN4) and are the
    # same values already used by Li_2012_clozapine.R in this package.
    # NOT paper-derived; see the vignette Errata.
    mw_cloz <- 326.83 # g/mol, clozapine
    mw_norcloz <- 312.80 # g/mol, norclozapine (N-desmethylclozapine)

    # 1. Individual parameters.
    ka <- exp(lka)

    # Clozapine apparent clearance: the smoking stratum selects which of
    # the two typical values applies, then the fixed allometric weight
    # term and the log-normal between-subject variability are applied.
    cl_typ <- exp(lcl_nonsmoke) * (1 - SMOKE) + exp(lcl_smoke) * SMOKE
    cl <- cl_typ * (WT / 70)^e_wt_cl * exp(etalcl)

    # Norclozapine apparent clearance. Smoking was screened here and not
    # retained, so a single typical value covers both strata.
    cl_norcloz <- exp(lcl_norcloz + etalcl_norcloz) * (WT / 70)^e_wt_cl_norcloz

    # Both apparent volumes are fixed literature values scaling linearly
    # with body weight; neither carries between-subject variability.
    vc <- exp(lvc) * (WT / 70)^e_wt_vc
    vc_norcloz <- exp(lvc_norcloz) * (WT / 70)^e_wt_vc_norcloz

    # 2. Micro-constants.
    kel <- cl / vc
    kel_norcloz <- cl_norcloz / vc_norcloz

    # 3. Relative bioavailability. Leponex is the fixed F = 1 reference;
    #    Luverina carries the estimated relative bioavailability and the
    #    only bioavailability random effect. FORM_CZP_LUVERINA is set per
    #    dose record, so a subject who crosses over changes branch between
    #    the two study periods.
    #    The Luverina branch is built on its own line so the eta stays
    #    mu-referenced for estimation.
    f_luverina <- exp(lfdepot_luverina + etalfdepot_luverina)
    frel <- exp(lfdepot_leponex) * (1 - FORM_CZP_LUVERINA) +
      f_luverina * FORM_CZP_LUVERINA
    f(depot) <- frel

    # 4. ODE system. One-compartment disposition for each substance, with
    #    the whole clozapine elimination flux feeding norclozapine
    #    (complete conversion, f fixed to 1) after the molecular-weight
    #    correction.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    d/dt(central_norcloz) <- (mw_norcloz / mw_cloz) * kel * central -
      kel_norcloz * central_norcloz

    # 5. Observations. Doses are in mg and the volumes in L, giving mg/L;
    #    the factor of 1000 converts to the ng/mL in which Olmos 2019
    #    reports both analytes (Table 1; assay range 54.8-1086 ng/mL for
    #    clozapine and 72.3-1085 ng/mL for norclozapine).
    Cc <- 1000 * central / vc
    Cc_norcloz <- 1000 * central_norcloz / vc_norcloz

    Cc ~ prop(propSd)
    Cc_norcloz ~ prop(propSd_norcloz)
  })
}
