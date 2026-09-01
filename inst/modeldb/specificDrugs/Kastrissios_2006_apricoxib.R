Kastrissios_2006_apricoxib <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order absorption",
    "and an absorption lag time for the selective cyclooxygenase-2",
    "(COX-2) inhibitor apricoxib (development code CS-706) in 104",
    "healthy adult volunteers across three phase 1 studies (2 to 800",
    "mg single or twice-daily doses for up to 14 days). The model",
    "carries three nonlinearities the authors identified and retained.",
    "(1) Relative bioavailability saturates with dose according to a",
    "standard saturation (Emax) model, Frel = D50 / (Dose + D50) with",
    "D50 = 221 mg, so Frel halves at 221 mg. (2) The same Frel is",
    "increased exp(0.351) = 1.42-fold for an evening dose relative to",
    "a morning dose, reproducing the apparent diurnal variation in",
    "exposure seen on both day 1 and day 14 of the twice-daily study.",
    "(3) Apparent oral clearance takes a separate typical value at the",
    "two supratherapeutic dose levels (400 and 800 mg, DOSE_HIGH = 1),",
    "19.5 L/h versus 34.1 L/h. Covariate effects on CL/F are sex",
    "(SEXF), the pooled poor-or-intermediate CYP2D6 phenotype",
    "(CYP2D6_PM_IM) and the CYP2C9 reduced-hydroxylator phenotype",
    "(CYP2C9_RH); central volume scales as a power of body weight",
    "normalised to the 73.3 kg cohort median. Exponential IIV on all",
    "six disposition and absorption parameters, proportional residual",
    "error. Absolute bioavailability F was not identifiable and was",
    "set to 1 by the authors, so all clearance and volume terms are",
    "apparent (CL/F, Vc/F, Vp/F, Q/F)."
  )
  reference <- paste(
    "Kastrissios H, Rohatagi S, Moberly J, Truitt K, Gao Y, Wada R,",
    "Takahashi M, Kawabata K, Salazar D. (2006).",
    "Development of a Predictive Pharmacokinetic Model for a Novel",
    "Cyclooxygenase-2 Inhibitor.",
    "Journal of Clinical Pharmacology 46(5):537-548.",
    "doi:10.1177/0091270006287122.",
    "The paper names the compound only by its Sankyo development code",
    "CS-706; the INN subsequently assigned to that molecule",
    "(2-(4-ethoxyphenyl)-4-methyl-1-(4-sulfamoylphenyl)-pyrrole,",
    "CAS 197904-84-0) is apricoxib, which this file uses per the",
    "library's generic-name-over-development-code convention.",
    sep = " "
  )
  vignette <- "Kastrissios_2006_apricoxib"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ng/mL"
  )

  compartmentData <- list(
    depot       = list(analyte = "apricoxib", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "apricoxib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "apricoxib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = paste(
        "Body weight. Enters the apparent central volume as a power",
        "term normalised to 73.3 kg, the development-cohort median",
        "weight, which Kastrissios 2006 prints literally inside its",
        "equation 8: Vc/F = (Vc/F)_Typ * (Weight / 73.3)^K_Vc/F-WT."
      ),
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Reference (normalising) value 73.3 kg is both the printed",
        "denominator of equation 8 and the Table II development-set",
        "median weight (range 49.5 to 104 kg), so the two agree and",
        "no assumption is needed. Treated as time-fixed per subject:",
        "the studies are 1 to 14 days long in healthy volunteers and",
        "the paper states only that 'baseline values were included in",
        "the NONMEM data set'. Body weight was ALSO tested on CL/F as",
        "a substitute for sex and was rejected (Table III run 18,",
        "'Effect of body weight on CL/F not significant, suggesting",
        "that gender is an independent covariate that cannot be",
        "explained by differences in body weight'), so WT acts on Vc/F",
        "only."
      ),
      source_name        = "Weight"
    ),
    SEXF = list(
      description        = paste(
        "Biological sex, 1 = female, 0 = male. Kastrissios 2006 uses",
        "the opposite orientation: its equation 7 indicator 'Gender'",
        "is 1 for MALE and 0 for female, so the canonical column is",
        "the value inversion SEXF = 1 - Gender."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "VALUE INVERSION, sign confirmed by the operator (sidecar",
        "request-001 q4, option A). The published coefficient",
        "K_CL/F-SEX = +0.325 multiplies the paper's MALE indicator, so",
        "it is applied here to (1 - SEXF) and the typical value",
        "lcl = log(34.1) is the FEMALE apparent oral clearance. Males",
        "are exp(0.325) = 1.384-fold higher, 47.2 L/h, which is",
        "exactly the abstract's headline 'Typical apparent clearance",
        "(CL/F) was 47.2 L/h' - an independent numeric anchor that",
        "pins the orientation. NOTE the paper's own prose mis-states",
        "the direction of this effect ('Apparent clearance was",
        "decreased by 38% in female subjects'; 'Female subjects were",
        "observed to have 38.4% lower clearance than were male",
        "subjects'): 38.4% is exp(0.325) - 1, the amount by which",
        "MALES exceed females, whereas the model-implied female-vs-",
        "male reduction is 1 - exp(-0.325) = 27.7%. The authors",
        "applied the correct '1 - exp(K)' arithmetic to their two",
        "negative coefficients (63.6% and 15.0%) and then reused the",
        "same phrasing on a POSITIVE coefficient, inverting the",
        "reference. This file encodes the estimate, not the prose;",
        "see the vignette Assumptions and deviations section."
      ),
      source_name        = "Gender"
    ),
    CYP2D6_PM_IM = list(
      description        = paste(
        "Pooled CYP2D6 poor-or-intermediate-metabolizer phenotype",
        "indicator. 1 = poor OR intermediate metabolizer, 0 =",
        "extensive OR ultrafast metabolizer."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "Kastrissios 2006 equation 7 indicator 'CYPD', defined in the",
        "paragraph following equation 8: 'CYPD is an indicator",
        "variable taking the value 0 for subjects with the",
        "extensive/ultrafast CYP 2D6 phenotype and the value 1 for",
        "subjects with the poor/intermediate phenotype.' Both poolings",
        "were tested and could not be rejected: Table III run 12",
        "combines poor with intermediate ('Objective function not",
        "significantly different ... indicating that poor and",
        "intermediate CYP 2D6 status can be combined') and run 17",
        "combines extensive with ultrafast. The Discussion gives the",
        "reason: intermediate and ultrafast metabolizers were fewer",
        "than 4% of treated subjects, so their separate estimates",
        "'may have been unreliable'. Mapping from the paper's Table II",
        "phenotype codes (1 = extensive, 2 = intermediate, 3 = poor,",
        "4 = ultrafast; counts 92/1/8/3 of 104): CYP2D6_PM_IM = 1 for",
        "codes 2 and 3 (n = 9, 8.7%), 0 for codes 1 and 4 (n = 95).",
        "Table I independently corroborates the split, quoting the",
        "Western proportions of the '(1&4, 2&3)' groups as 0.91 and",
        "0.09. Because the paper estimated ONE theta for the pooled",
        "group, this single pooled column is used rather than the",
        "CYP2D6_PM + CYP2D6_IM pair sharing a coefficient."
      ),
      source_name        = "CYPD"
    ),
    CYP2C9_RH = list(
      description        = paste(
        "CYP2C9 reduced-hydroxylator phenotype indicator. 1 = reduced",
        "hydroxylator, 0 = intermediate metabolizer OR normal",
        "hydroxylator."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "Kastrissios 2006 equation 7 indicator 'CYPC', defined in the",
        "paragraph following equation 8: 'CYPC is an indicator",
        "variable taking the value 0 for subjects with the",
        "intermediate/normal hydroxylator CYP 2C9 phenotype and the",
        "value 1 for those with the reduced hydroxylator phenotype.'",
        "The 1-level is phenotype code 5 in the paper's shared 1-6",
        "code legend; Table II reports 'CYP 2C9, 2/5/6 = 45/16/43' of",
        "104 subjects, so CYP2C9_RH = 1 for n = 16 (15.4%) and 0 for",
        "the 88 intermediate-metabolizer and normal-hydroxylator",
        "subjects pooled. Table I corroborates the split, quoting the",
        "Western proportions of the '(1&6, 5)' groups as 0.85 and",
        "0.15. Note that the paper's REFERENCE group contains its",
        "'intermediate metabolizer' subjects, which is why the",
        "registered CYP2C9_PM_IM canonical must NOT be substituted",
        "here: doing so would move those 45 subjects into the 1-level",
        "and apply the -0.163 clearance effect to them. A minor",
        "internal inconsistency in the paper (the footnote legend and",
        "Table II give CYP2C9 levels 2/5/6 while Table I groups the",
        "same data as '(1&6, 5)') does not affect the encoding: the",
        "1-level is unambiguously the level-5 reduced-hydroxylator",
        "group either way. The effect was only marginally significant",
        "under FOCE (P < .1; Table III run 19 notes the '95%",
        "confidence interval includes 0, but not 90% CI') and was",
        "retained by the authors 'because of its pharmacologic",
        "significance'."
      ),
      source_name        = "CYPC"
    ),
    DOSE_HIGH = list(
      description        = paste(
        "Supratherapeutic-dose-cohort indicator. 1 = the subject",
        "received a 400 mg or 800 mg dose of apricoxib, 0 = the",
        "subject received a dose in the 2 to 200 mg range."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "Kastrissios 2006 equation 7 indicator 'D_hi': 'Dhi is an",
        "indicator variable that took the value 0 for doses in the",
        "range 2 to 200 mg and 1 for doses of 400 or 800 mg.' Gates",
        "the bracketed two-typical-value apparent oral clearance",
        "split of equation 7, CL/F = [(CL/F)_Typ * (1 - Dhi) +",
        "(CL/F)_Typ,hi * Dhi] * (covariate factors): DOSE_HIGH = 0",
        "selects 34.1 L/h and DOSE_HIGH = 1 selects 19.5 L/h",
        "(Table IV), a 42.8% reduction matching the abstract's",
        "'reduced by 43% at doses greater than 200 mg'. Time-fixed",
        "per subject in this cohort: the 400 and 800 mg levels were",
        "studied only as single doses in study 2, so a subject who",
        "carries DOSE_HIGH = 1 received no other dose level.",
        "Deliberately carried as its own indicator column rather than",
        "derived as an inequality on the dose amount, per the",
        "DOSE_LOW register entry's rule - the authors defined this",
        "group by the studied dose levels (400, 800) and not by a",
        "threshold, so an inequality would silently classify an",
        "intermediate dose that was never studied. Distinct from the",
        "saturable relative-bioavailability nonlinearity, which is a",
        "smooth function of the actual dose amount and is read from",
        "the event record via podo(depot); the paper notes that",
        "'because both apparent oral clearance and Frel were reduced",
        "at doses greater than 200 mg, it is difficult to delineate",
        "these 2 effects' but retained both because each improved the",
        "fit."
      ),
      source_name        = "Dhi"
    )
  )

  # Screened by the authors and NOT retained in the final model, so these are
  # documented for provenance only and are deliberately absent from model().
  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Age.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Tested as a candidate covariate (Methods, 'Age, gender, race,",
        "body weight, height, body mass index, and phenotypes of CYP",
        "2D6, CYP 2C9, and CYP 2C19 were considered as possible",
        "covariates') and not retained. Discussion: 'Although age did",
        "not influence CS-706 pharmacokinetics, the study population",
        "consisted of young, healthy adults, aged 19 to 46 years,",
        "therefore, variability in CS-706 pharmacokinetics due to age",
        "may not be adequately addressed by this study population.'",
        "No point estimate is reported, so no effect can be encoded."
      ),
      source_name        = "Age"
    ),
    HT = list(
      description        = "Height.",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened per Methods and not retained; no point estimate",
        "reported. Table II development-set median 173 cm (150-191)."
      ),
      source_name        = "Height"
    ),
    BMI = list(
      description        = "Body mass index.",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened per Methods and not retained; no point estimate",
        "reported. Table II development-set median 25.3 kg/m^2",
        "(19.3-29.3)."
      ),
      source_name        = "BMI"
    ),
    RACE_WHITE = list(
      description        = "White race indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "Race was screened as a categorical covariate on CL/F",
        "(Methods: 'indicator variables were used to evaluate the",
        "effects of gender, race, and cytochrome P450 phenotype on",
        "apparent oral clearance') and was not retained - the Results",
        "state 'gender, CYP 2D6- and CYP 2C9-predicted phenotypes",
        "correlated with CL/F, whereas Vc/F was correlated with body",
        "weight. No other parameter-covariate relationships were found",
        "to be significant.' Table II development-set race counts",
        "W/B/A/H/O = 92/3/2/4/3 of 104. No point estimate is reported.",
        "Race nonetheless matters to the paper's conclusion",
        "indirectly, through the CYP phenotype frequencies used in the",
        "Japanese-versus-Western simulation (Table I), not through a",
        "direct effect on any PK parameter."
      ),
      source_name        = "Race"
    ),
    CYP2C19_PM = list(
      description        = "CYP2C19 poor-metabolizer phenotype indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "CYP2C19 phenotype was genotyped and screened (Methods) but",
        "no CYP2C19 effect was retained. Table II reports 'CYP 2C19,",
        "1/3 = 103/1' - a single poor metabolizer in 104 subjects,",
        "which cannot support an estimate. No point estimate reported."
      ),
      source_name        = "CYP 2C19"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 104L,
    n_studies      = 3L,
    age_range      = "19-46 years",
    age_median     = "27.0 years",
    weight_range   = "49.5-104 kg",
    weight_median  = "73.3 kg",
    sex_female_pct = 50,
    race_ethnicity = c(White = 88.5, Black = 2.9, Asian = 1.9, Hispanic = 3.8, Other = 2.9),
    disease_state  = "healthy adult volunteers",
    dose_range     = paste(
      "Oral apricoxib (CS-706). Study 1: single doses of 2, 5, 10,",
      "25, 50, 100 or 200 mg (n = 56 active), serial sampling to 72 h.",
      "Study 2: single doses of 400 or 800 mg (n = 16 active), serial",
      "sampling to 72 h (400 mg) or 144 h (800 mg). Study 3: 25, 100",
      "or 200 mg twice daily for 14 days (n = 32 active), serial",
      "sampling over 12 h after the first two and final two doses plus",
      "troughs on days 2-13 and to 72 h after the last dose. All doses",
      "administered under fasted conditions."
    ),
    regions        = "United States (single phase 1 site: MDS Pharma Services, Lincoln, Nebraska)",
    notes          = paste(
      "Baseline characteristics from Table II, development data set",
      "column. n_subjects = 104 counts only apricoxib-treated",
      "subjects (8 per dose level, except 16 subjects at 100 mg twice",
      "daily in study 3; 56/16/32 across studies 1/2/3). The",
      "abstract's '130 subjects' counts everyone randomised including",
      "placebo, because subjects were randomised to active or placebo",
      "at a 4:1 ratio (8 active and 2 placebo per dose level) and",
      "'placebo and missing CS-706 plasma concentration data were",
      "excluded from the NONMEM data sets' (Methods). Sex was balanced",
      "at 52 male / 52 female. Race percentages are the Table II",
      "counts 92/3/2/4/3 expressed as percentages of 104. Estimation",
      "was FOCE in NONMEM version 5 (final model, Table III run 19,",
      "objective function 15000). A SEPARATE 80-subject phase 1",
      "endoscopic-surveillance study supplied the external validation",
      "data set (one trough concentration 24 h after the seventh daily",
      "100 or 200 mg dose); its subjects are NOT part of the 104 and",
      "its CYP phenotypes were not collected, so the validation",
      "simulation drew phenotype distributions from the development",
      "set."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters. All values are Kastrissios 2006 Table IV
    # ("Final Pharmacokinetic Parameter Estimates"), a FOCE fit in
    # NONMEM v5 (Table III run 19). Nothing in Table IV is flagged
    # FIX and every row carries a CV% and a 95% confidence interval,
    # so every parameter below is estimated, not fixed. The single
    # genuinely fixed quantity in the paper is absolute oral
    # bioavailability F, which "was not uniquely estimable and
    # therefore was assumed to be equal to the product of 1 and Frel"
    # - i.e. the leading constant 1 of equation 6. That constant is
    # not a free parameter and so appears in model() as the implicit
    # unit multiplier of frel rather than as an ini() entry; every
    # clearance and volume below is consequently APPARENT (divided by
    # the unknown true F).
    # ------------------------------------------------------------------
    lvc <- log(166)
    label("Apparent central volume Vc/F at the 73.3 kg reference weight (L)")
    # Table IV: Vc/F = 166 L, CV 8%, 95% CI (141, 191).
    # Abstract: "Typical apparent volume of the central compartment was 166 L".

    lvp <- log(483)
    label("Apparent peripheral volume Vp/F (L)")
    # Table IV: Vp/F = 483 L, CV 7%, 95% CI (417, 549).

    lq <- log(75)
    label("Apparent inter-compartmental clearance Q/F (L/h)")
    # Table IV: Q/F = 75 L/h, CV 7%, 95% CI (64, 86).

    lcl <- log(34.1)
    label("Apparent oral clearance CL/F, 2-200 mg doses, female reference subject (L/h)")
    # Table IV row "CL/F, L/h, extensive metabolizer phenotype / 2- to
    # 200-mg doses" = 34.1 L/h, CV 6%, 95% CI (29.8, 38.4).
    # This is the value of the bracketed typical-value term of equation 7
    # at Dhi = 0, and it is the FEMALE value: equation 7's Gender
    # indicator is 1 for MALE, so the male typical value is
    # 34.1 * exp(0.325) = 47.2 L/h, matching the abstract exactly
    # ("Typical apparent clearance (CL/F) was 47.2 L/h"). See the
    # SEXF covariateData notes.

    lcl_highdose <- log(19.5)
    label("Apparent oral clearance CL/F, 400-800 mg doses, female reference subject (L/h)")
    # Table IV row "CL/F ... / 400- to 800-mg doses" = 19.5 L/h,
    # CV 17%, 95% CI (12.9, 26.1). Selected inside model() by
    # DOSE_HIGH = 1. 1 - 19.5/34.1 = 42.8%, the abstract's "reduced by
    # 43% at doses greater than 200 mg".

    lka <- log(0.542)
    label("First-order absorption rate constant KA (1/h)")
    # Table IV: KA = 0.542 1/h, CV 5%, 95% CI (0.489, 0.595).
    # Discussion: "a first-order KA with a median absorption half-life
    # of 1.3 hours (23% intersubject variability)";
    # log(2)/0.542 = 1.279 h and sqrt(0.052) = 22.8%.

    ltlag <- log(0.236)
    label("Absorption lag time TLAG (h)")
    # Table IV: TLAG = 0.236 h, CV 1%, 95% CI (0.232, 0.240).
    # Discussion: "a population median TLAG of 0.23 hours (14 minutes)";
    # 0.236 h * 60 = 14.2 min.

    led50 <- log(221)
    label("Dose at which relative bioavailability is half maximal, D50 (mg)")
    # Table IV row "D50 for F, mg" = 221 mg, CV 14%, 95% CI (161, 281).
    # Definition from the text preceding equation 6: "the dose at which
    # Frel is half maximal". Abstract: bioavailability "decreased
    # saturably with increasing dose (50% reduction at 221 mg)";
    # 221/(221 + 221) = 0.500 exactly.

    # ------------------------------------------------------------------
    # Covariate effects. Kastrissios 2006 equation 5 gives the
    # categorical form theta_i = theta_T * exp(Cov_i * Kcov), so every
    # coefficient below is an exponential (log-scale) shift; equation 4
    # gives the continuous form theta_i = theta_T * (Cov_i/Cov_med)^Kcov,
    # so the body-weight term is a power.
    # ------------------------------------------------------------------
    e_evening_fdepot <- 0.351
    label("Log-scale increase in relative bioavailability for an evening dose (unitless)")
    # Table IV row "KFrel-PM Dose" = 0.351, CV 13%, 95% CI (0.261, 0.441).
    # Equation 6: Frel = 1 * exp(KFrel-PM * PM) * D50/(DOSE + D50), where
    # "PM is an indicator variable taking on the value of 0 for daytime
    # doses and 1 for evening doses". exp(0.351) - 1 = 42.0%, the
    # abstract's "Bioavailability increased by 42% after nighttime doses".

    e_sexf_cl <- 0.325
    label("Log-scale effect of MALE sex on apparent oral clearance (unitless)")
    # Table IV row "KCL/F-SEX" = 0.325, CV 24%, 95% CI (0.171, 0.479).
    # Equation 7 applies this to the paper's Gender indicator (1 = male),
    # so model() applies it to (1 - SEXF). exp(0.325) = 1.384.

    e_cyp2d6_pmim_cl <- -1.01
    label("Log-scale effect of pooled poor/intermediate CYP2D6 phenotype on CL/F (unitless)")
    # Table IV row "KCL/F-CYP2D6 (poor and intermediate metabolizer)"
    # = -1.01, CV 16%, 95% CI (-1.33, -0.69).
    # 1 - exp(-1.01) = 63.6%, the Discussion's "63.6% lower apparent
    # oral clearance compared to extensive and ultrafast metabolizers".

    e_cyp2c9_rh_cl <- -0.163
    label("Log-scale effect of reduced-hydroxylator CYP2C9 phenotype on CL/F (unitless)")
    # Table IV row "KCL/F-CYP2C9 (reduced hydroxylator)". Table IV prints
    # the estimate as "-00.163", a typesetting artifact of "-0.163";
    # the value is confirmed by its own printed 95% CI (-0.350, 0.024)
    # and CV 59%, since -0.163 +/- 1.96 * 0.59 * 0.163 = (-0.351, 0.025).
    # 1 - exp(-0.163) = 15.0%, the Discussion's "Apparent oral clearance
    # was decreased by 15% in individuals expressing the reduced
    # hydroxylator CYP 2C9 phenotype".

    e_wt_vc <- 0.831
    label("Power exponent for body weight on apparent central volume (unitless)")
    # Table IV row "KVc/F-WT" = 0.831, CV 33%, 95% CI (0.288, 1.374).
    # Equation 8: Vc/F = (Vc/F)_Typ * (Weight/73.3)^KVc/F-WT.
    # Discussion cross-check: "an increase of 10 kg greater than the
    # median body weight resulted in an approximately 10% increase in
    # Vc/F"; (83.3/73.3)^0.831 = 1.112, i.e. +11.2%.

    # ------------------------------------------------------------------
    # Inter-individual variability. Equation 1 is theta_i = theta_T *
    # exp(eta_i), i.e. exponential IIV, and Table IV reports the
    # omega^2 VARIANCES directly (the "Estimated Variability (%)"
    # column is footnote b's "intersubject variability = sqrt(omega^2)",
    # which reproduces each row: sqrt(0.127) = 35.6%, sqrt(0.234) =
    # 48.4%, sqrt(0.229) = 47.9%, sqrt(0.125) = 35.4%, sqrt(0.052) =
    # 22.7%, sqrt(0.002) = 4.5% vs the printed 4.36%). The values below
    # are therefore entered as variances without transformation. The
    # paper reports no off-diagonal covariances, so the omega matrix is
    # diagonal here; see the vignette Assumptions and deviations
    # section. IIV is carried on all six disposition and absorption
    # parameters but on none of the covariate effects or on D50.
    # ------------------------------------------------------------------
    etalvc   ~ 0.127   # Table IV: omega^2 Vc/F  = 0.127, CV 28%, 95% CI (0.057, 0.197); 35.6% CV
    etalvp   ~ 0.234   # Table IV: omega^2 Vp/F  = 0.234, CV 15%, 95% CI (0.166, 0.302); 48.4% CV
    etalq    ~ 0.229   # Table IV: omega^2 Q/F   = 0.229, CV 20%, 95% CI (0.140, 0.318); 47.9% CV
    etalcl   ~ 0.125   # Table IV: omega^2 CL/F  = 0.125, CV 16%, 95% CI (0.086, 0.164); 35.4% CV
    etalka   ~ 0.052   # Table IV: omega^2 Ka    = 0.052, CV 47%, 95% CI (0.004, 0.099); 22.7% CV
    etaltlag ~ 0.002   # Table IV: omega^2 Tlag  = 0.002, CV 64%, 95% CI (0.000, 0.004); 4.36% CV

    # ------------------------------------------------------------------
    # Residual error. Equation 3 is the FINAL model, y_ij = yhat_ij *
    # (1 + eps_ij), i.e. proportional in linear space. (Equation 2, the
    # log-normal alternative, was tested as Table III run 5 and rejected
    # - "No improvement in diagnostic plots compared to #4".) Table IV
    # reports sigma^2 = 0.069, a VARIANCE, so the nlmixr2 SD is its
    # square root; footnote b's "residual variability = sqrt(sigma^2)"
    # gives the printed 26.3%.
    # ------------------------------------------------------------------
    propSd <- sqrt(0.069)
    label("Proportional residual error (fraction)")
    # Table IV: sigma^2 = 0.069, CV 6%, 95% CI (0.061, 0.077);
    # sqrt(0.069) = 0.2627, matching the printed 26.3%.
  })

  model({
    # ------------------------------------------------------------------
    # 1. Evening-dose indicator.
    #
    #    Kastrissios 2006 equation 6 carries a per-dose indicator PM,
    #    "0 for daytime doses and 1 for evening doses". Derived here
    #    from the model clock rather than from a covariate column,
    #    following the in-repo precedent of
    #    Zhang_2013_lopinavir_ritonavir.R, which encodes the same
    #    concept (its Table 2 "Evening effect on bioavailability")
    #    the same way. rxode2's `t` is the integration time in hours
    #    from the first event, so `t mod 24` is the hour of the dosing
    #    day and a dose in the second half of it is the evening dose.
    #    With the source studies' regimens this is exact: the
    #    twice-daily study doses ~12 h apart (morning PM = 0, evening
    #    PM = 1) and the single-dose and once-daily studies dose in the
    #    morning (t mod 24 = 0, PM = 0). A user simulating a different
    #    clock should shift the event times so that hour 0 is the
    #    morning dose.
    # ------------------------------------------------------------------
    evening <- (t - 24 * floor(t / 24)) >= 12

    # ------------------------------------------------------------------
    # 2. Individual parameters.
    #
    #    Apparent oral clearance, equation 7:
    #      CL/F = [(CL/F)_Typ * (1 - Dhi) + (CL/F)_Typ,hi * Dhi]
    #             * exp(Sex  * K_CL/F-SEX)
    #             * exp(CYPD * K_CL/F-CYP2D6)
    #             * exp(CYPC * K_CL/F-CYP2C9)
    #    The bracketed two-value switch is written on the log scale as
    #    lcl + (lcl_highdose - lcl) * DOSE_HIGH, which returns exactly
    #    exp(lcl) = 34.1 at DOSE_HIGH = 0 and exp(lcl_highdose) = 19.5
    #    at DOSE_HIGH = 1 (the Viberg_2012_AZD6088_rat.R precedent).
    #    The paper's Gender indicator is 1 for MALE, hence (1 - SEXF).
    # ------------------------------------------------------------------
    cl <- exp(lcl + (lcl_highdose - lcl) * DOSE_HIGH +
                e_sexf_cl * (1 - SEXF) +
                e_cyp2d6_pmim_cl * CYP2D6_PM_IM +
                e_cyp2c9_rh_cl * CYP2C9_RH +
                etalcl)

    # Apparent central volume, equation 8: power of weight normalised
    # to the 73.3 kg cohort median printed in the equation itself.
    vc <- exp(lvc + etalvc) * (WT / 73.3)^e_wt_vc

    # Vp/F, Q/F, KA and TLAG carry no covariate effects: "No other
    # parameter-covariate relationships were found to be significant."
    vp   <- exp(lvp + etalvp)
    q    <- exp(lq + etalq)
    ka   <- exp(lka + etalka)
    tlag <- exp(ltlag + etaltlag)

    # ------------------------------------------------------------------
    # 3. Relative bioavailability, equation 6:
    #      Frel = 1 * exp(K_Frel-PM * PM) * D50 / (DOSE + D50)
    #    DOSE is the amount of the dose being administered, read
    #    directly off the event record with podo(depot) so it can never
    #    drift out of sync with `amt`. The leading 1 is the authors'
    #    non-identifiable absolute bioavailability F ("assumed to be
    #    equal to the product of 1 and Frel"), so Frel -> 1 as the dose
    #    approaches 0 and Frel = 0.5 at DOSE = D50 = 221 mg.
    # ------------------------------------------------------------------
    ed50 <- exp(led50)
    frel <- exp(e_evening_fdepot * evening) * ed50 / (podo(depot) + ed50)

    # ------------------------------------------------------------------
    # 4. Micro-constants and the two-compartment system with
    #    first-order absorption. "A 2-compartment linear model captured
    #    the apparent biphasic disposition of CS-706"; the Discussion
    #    states "The model assumes that drug distribution and
    #    elimination occur from a central compartment, with a linked
    #    peripheral compartment." A third compartment (Table III run 3)
    #    and a parallel second absorption process (run 9) were both
    #    tested and rejected.
    # ------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 5. Dose-level and time-of-day dependent bioavailability, plus the
    #    absorption lag time. "Absorption of CS-706 after oral dose
    #    administration was characterized by first-order absorption and
    #    an absorption lag time (TLAG)."
    f(depot)    <- frel
    alag(depot) <- tlag

    # ------------------------------------------------------------------
    # 6. Observation. The assay reports apricoxib in ng/mL over a
    #    2-1000 ng/mL linear range (LLOQ 2 ng/mL), and doses are in mg
    #    with volumes in L, so central/vc is in mg/L = ug/mL and the
    #    1000 converts to ng/mL. Residual error is proportional
    #    (equation 3), so the unit scaling cancels out of the fit.
    # ------------------------------------------------------------------
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
