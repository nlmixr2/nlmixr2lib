Nakashima_2015_valproic_acid_ge19y <- function() {
  description <- paste0(
    "Logistic-regression exposure-response model for the probability of an ",
    "over-50% reduction in seizure frequency in Japanese patients with ",
    "epilepsy aged 19 YEARS OR OLDER on sustained-release valproic acid ",
    "(VPA), driven by the model-predicted steady-state trough VPA ",
    "concentration (Nakashima 2015 Eq 8; n = 21, mean age 26.16 +/- 5.34 ",
    "years, 132 concentration measurement points). The smaller of the two ",
    "age-subgroup analyses run alongside the all-age final model, and the ",
    "most structurally unusual of the three: only the SCN1A A/A term ",
    "multiplies the trough concentration, so G/G and G/A adults get a ",
    "CONCENTRATION-INDEPENDENT response probability. That is the equation ",
    "exactly as printed -- Eq 8 carries no parentheses, unlike Eqs 6 and 7 ",
    "-- and is encoded verbatim by operator ruling; the ambiguity is ",
    "documented in the ini() block and in the vignette Errata, and no ",
    "source anchor exists to test it. A one-compartment first-order ",
    "absorption model with lag time supplies the trough; its parameters ",
    "were FIXED into this PD fit from the group's upstream population PK ",
    "analysis (Ogusu 2014 Eqs 5-8 and Table 2). The paper carries a PLOS ",
    "ONE Expression of Concern on funding-disclosure grounds only (no ",
    "concern about data, methods or results). Companions: ",
    "Nakashima_2015_valproic_acid (Eq 6, all ages) and ",
    "Nakashima_2015_valproic_acid_le18y (Eq 7, <= 18 y)."
  )
  reference <- paste(
    "Nakashima H, Oniki K, Nishimura M, Ogusu N, Shimomasuda M, Ono T,",
    "Matsuda K, Yasui-Furukori N, Nakagawa K, Ishitsu T, Saruwatari J.",
    "Determination of the Optimal Concentration of Valproic Acid in Patients",
    "with Epilepsy: A Population Pharmacokinetic-Pharmacodynamic Analysis.",
    "PLoS ONE. 2015;10(10):e0141266. doi:10.1371/journal.pone.0141266",
    "(Eq 8, the subgroup aged 19 years or older).",
    "SUBJECT OF AN EXPRESSION OF CONCERN:",
    "The PLOS ONE Editors. Expression of Concern: Determination of the",
    "Optimal Concentration of Valproic Acid in Patients with Epilepsy: A",
    "Population Pharmacokinetic-Pharmacodynamic Analysis.",
    "PLoS ONE. 2023;18(1):e0279487. doi:10.1371/journal.pone.0279487",
    "(grounds: the study was part-funded by the Smoking Research Foundation,",
    "which received tobacco-industry support, contrary to the journal's 2010",
    "tobacco-funding policy; the notice raises no concern about the data,",
    "methods or results and is not a retraction).",
    "The fixed PK layer is taken from the upstream population PK analysis:",
    "Ogusu N, Saruwatari J, Nakashima H, Noai M, Nishimura M, Eshima N,",
    "Oniki K, Ogata Y, Nakagawa K, Ishitsu T. Impact of the superoxide",
    "dismutase 2 Val16Ala polymorphism on the relationship between valproic",
    "acid exposure and elevation of gamma-glutamyltransferase in patients",
    "with epilepsy: a population pharmacokinetic-pharmacodynamic analysis.",
    "PLoS One. 2014;9(11):e111066. doi:10.1371/journal.pone.0111066.",
    sep = " "
  )
  vignette <- "Nakashima_2015_valproic_acid_seizure_control"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ug/mL"
  )

  covariateData <- list(
    SEXF = list(
      description        = "Sex indicator; 1 = female, 0 = male. Acts on the PK layer only (apparent oral clearance).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste0(
        "Ogusu 2014 Eq 7 / Table 2 'Gender on CL/F' = 0.917, applied as a ",
        "power form 0.917^SEXF (females have 8.3% lower CL/F). Sex carries ",
        "no PD coefficient in any of the three Nakashima 2015 models. The ",
        "subgroup's own sex split is not reported; the full cohort was ",
        "37.7% female (Table 1)."
      ),
      source_name        = "female"
    ),
    DOSE_VPA_MGD = list(
      description        = "Patient's own TOTAL daily valproic acid dose, not normalised by body weight. Drives both apparent clearance and apparent volume in the fixed PK layer.",
      units              = "mg/d",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Ogusu 2014 Eqs 6-7 scale both Vd/F and CL/F by (Dose/1000)^theta ",
        "with Dose the daily VPA dose in mg/day, so 1000 mg/day is the ",
        "reference. Must be strictly positive. This covariate is the DAILY ",
        "total, which equals the per-administration event-table amt only ",
        "for once-daily dosing. The subgroup's own dose range is not ",
        "reported; the full cohort was 1120.0 +/- 592.5 [50-3200] mg/day ",
        "(Table 1). Adults in this cohort were started on 400-1200 mg/day ",
        "and escalated by 200 mg/day steps (Methods)."
      ),
      source_name        = "Dose"
    ),
    CONMED_CBZ = list(
      description        = "Concomitant carbamazepine indicator; 1 = co-administered, 0 = not. Acts on the PK layer only in this model.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant carbamazepine)",
      notes              = paste0(
        "Ogusu 2014 Table 2 'CBZ on CL/F' = 1.19, power form ",
        "1.19^CONMED_CBZ (+19% CL/F). Carbamazepine carries a PD ",
        "coefficient in the all-age Eq 6 (-1.75 intercept) and a large one ",
        "in the <= 18 y Eq 7 (-4.88 intercept, +4.73 slope) but is ABSENT ",
        "from Eq 8, so it has no PD effect here. Full-cohort prevalence 47 ",
        "of 77 (61.0%; Table 1)."
      ),
      source_name        = "CBZ"
    ),
    CONMED_PB = list(
      description        = "Concomitant phenobarbital indicator; 1 = co-administered, 0 = not. Acts on the PK layer only in this model.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant phenobarbital)",
      notes              = paste0(
        "Ogusu 2014 Table 2 'PB on CL/F' = 1.12, power form ",
        "1.12^CONMED_PB (+12% CL/F). Phenobarbital carries a PD coefficient ",
        "in the <= 18 y Eq 7 (-1.93 intercept) but is absent from Eq 8. ",
        "Full-cohort prevalence 31 of 77 (40.3%; Table 1)."
      ),
      source_name        = "PB"
    ),
    CONMED_PHT = list(
      description        = "Concomitant phenytoin indicator; 1 = co-administered, 0 = not. Acts on the PK layer only in this model.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant phenytoin)",
      notes              = paste0(
        "Ogusu 2014 Table 2 'PHT on CL/F' = 1.43, power form ",
        "1.43^CONMED_PHT (+43% CL/F). Phenytoin carries a PD slope ",
        "coefficient in the all-age Eq 6 (-3.62) and in the <= 18 y Eq 7 ",
        "(-3.86) but is absent from Eq 8. Full-cohort prevalence 22 of 77 ",
        "(28.6%; Table 1)."
      ),
      source_name        = "PHT"
    ),
    CONMED_CLB = list(
      description        = "Concomitant clobazam indicator; 1 = co-administered, 0 = not. Acts on the PK layer only.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant clobazam)",
      notes              = paste0(
        "Ogusu 2014 Table 2 'CLB on CL/F' = 0.906, power form ",
        "0.906^CONMED_CLB (-9.4% CL/F). Clobazam carries no PD coefficient ",
        "in any of the three Nakashima 2015 models. Full-cohort prevalence ",
        "33 of 77 (42.9%; Table 1)."
      ),
      source_name        = "CLB"
    ),
    CONMED_CZP = list(
      description        = "Concomitant clonazepam indicator; 1 = co-administered, 0 = not. Acts on the PD logit intercept only.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant clonazepam)",
      notes              = paste0(
        "Nakashima 2015 Eq 8 gives -2.56 on the logit intercept, more than ",
        "twice the all-age Eq 6 value of -1.18, and clonazepam is the only ",
        "concomitant AED retained at all in this subgroup. Not a PK ",
        "covariate: clonazepam was not retained on VPA clearance in the ",
        "upstream Ogusu 2014 model. The authors state the mechanism of the ",
        "CZP association is unclear (Discussion). Full-cohort prevalence 21 ",
        "of 77 (27.3%; Table 1)."
      ),
      source_name        = "CZP"
    ),
    SNP_SCN1A_RS3812718_GA = list(
      description        = "SCN1A rs3812718 (IVS5-91 G>A) heterozygous G/A genotype indicator; 1 = G/A, 0 = otherwise. Acts on the PD logit INTERCEPT only in this model.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 in combination with SNP_SCN1A_RS3812718_AA = 0, i.e. the G/G genotype group",
      notes              = paste0(
        "Nakashima 2015 Eq 8 gives -9.88 on the intercept, the largest ",
        "genotype effect anywhere in the paper (compare -5.87 in the ",
        "all-age Eq 6 and -4.75 in the <= 18 y Eq 7). Unlike the other two ",
        "models, this term does NOT modify the exposure slope: as printed, ",
        "Eq 8's only concentration-dependent term is the A/A one -- see the ",
        "ini() block note on e_snp_scn1a_rs3812718_aa_slope. Full-cohort ",
        "genotype frequencies G/G 11.7%, G/A 53.2%, A/A 35.1% (Table 1); ",
        "with only 21 patients in this subgroup the G/G reference group is ",
        "very small and no subgroup genotype counts are reported."
      ),
      source_name        = "SCN1A G/A genotype"
    ),
    SNP_SCN1A_RS3812718_AA = list(
      description        = "SCN1A rs3812718 (IVS5-91 G>A) homozygous A/A genotype indicator; 1 = A/A, 0 = otherwise. As Eq 8 is printed, this is the ONLY term that multiplies the trough concentration.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 in combination with SNP_SCN1A_RS3812718_GA = 0, i.e. the G/G genotype group",
      notes              = paste0(
        "Nakashima 2015 Eq 8 gives -14.3, and in the printed equation this ",
        "coefficient is the one multiplied by the predicted trough ",
        "concentration. Consequence, which a user of this model must ",
        "understand: G/G and G/A adults have NO concentration-response at ",
        "all under Eq 8 as printed, so their predicted probability is ",
        "constant in dose. See the extended note on ",
        "e_snp_scn1a_rs3812718_aa_slope in ini() for why the literal ",
        "reading is also the only parenthesisation consistent with the ",
        "printed signs, and for the absence of any source anchor to test ",
        "it. Full-cohort prevalence 27 of 77 (35.1%; Table 1)."
      ),
      source_name        = "SCN1A A/A genotype"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Patient age.",
      units       = "years",
      type        = "continuous",
      notes       = paste0(
        "Age DEFINES this subgroup (19 years or older) rather than acting ",
        "as a covariate within it, so it carries no coefficient in Eq 8 -- ",
        "unlike the all-age Eq 6, where age enters the logit intercept at ",
        "0.98 per year. Subgroup mean age 26.16 +/- 5.34 years; the cohort ",
        "maximum was 36.9 years (Results and Table 1). Listed here so the ",
        "subgroup definition is machine-visible without implying an ",
        "unpublished coefficient."
      )
    ),
    CONMED_TPM = list(
      description = "Concomitant topiramate indicator; 1 = co-administered, 0 = not.",
      units       = "(binary)",
      type        = "binary",
      notes       = paste0(
        "Retained on the exposure slope in the all-age Eq 6 (-1.73) but ",
        "absent from Eq 8, so no coefficient exists for this subgroup. ",
        "Full-cohort prevalence 10 of 77 (13.0%; Table 1)."
      )
    ),
    SEIZURE_LOCUS_PARTIAL = list(
      description = "Seizure-locus indicator; 1 = partial (focal), 0 = generalized.",
      units       = "(binary)",
      type        = "binary",
      notes       = paste0(
        "Retained on the exposure slope in the all-age Eq 6 (+2.41) but ",
        "absent from Eq 8. Note that the gloss printed immediately AFTER ",
        "Eq 8 still defines 'partial seizure = 1 if a partial seizure was ",
        "present, and 0 if partial seizures were absent' even though ",
        "neither Eq 7 nor Eq 8 contains such a term; that gloss is shared ",
        "boilerplate carried over from Eq 6 (it also lists 'CBZ, PB, CZP or ",
        "PHT' when Eq 8 uses only CZP) and must NOT be read as evidence ",
        "that a term was dropped in typesetting. Full-cohort split partial ",
        "87.0% / generalized 13.0% (Table 1)."
      )
    ),
    MENT_DISABLED = list(
      description = "Complication with intellectual disability; 1 = present, 0 = absent.",
      units       = "(binary)",
      type        = "binary",
      notes       = paste0(
        "Screened on the all-age model and dropped in backward elimination ",
        "(Nakashima 2015 Table 2 and Results); never appears in Eq 8. No ",
        "point estimate is published. Full-cohort prevalence 59 of 77 ",
        "(76.6%; Table 1)."
      )
    )
  )

  compartmentData <- list(
    depot   = list(analyte = "valproic acid", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "valproic acid", units = "mg", specimen = "serum", verified = TRUE)
  )

  # See Nakashima_2015_valproic_acid.R: `prob_seizure50` is an algebraic
  # probability output, registered as a canonical member of the
  # `prob_<endpoint>` PD-output family in compartment-names.md.

  population <- list(
    species        = "human",
    n_subjects     = 21L,
    n_studies      = 1L,
    n_observations = "132 VPA concentration measurement points; one binary efficacy outcome per patient",
    age_range      = "19-36.9 years (subgroup lower bound is the inclusion cut-off; cohort maximum 36.9 years)",
    age_median     = "mean 26.16 +/- 5.34 years (median not reported)",
    sex_female_pct = NA_real_,
    race_ethnicity = c(Asian = 100),
    disease_state  = paste0(
      "Epilepsy on maintenance sustained-release valproic acid, adult ",
      "subgroup (aged 19 years or older). Idiopathic epilepsy and severe ",
      "retardation were exclusion criteria."
    ),
    dose_range     = "not reported for this subgroup; full cohort 1120.0 +/- 592.5 [50-3200] mg/day sustained-release valproic acid (adults initiated at 400-1200 mg/day)",
    regions        = "Japan (Kumamoto Saishunso National Hospital), June 1989 to April 2011",
    notes          = paste0(
      "Smallest of the paper's three models: 21 patients and 132 ",
      "concentration measurement points (Results, Effect of the patient age ",
      "on the PK-PD modeling). Nakashima 2015 reports NO parameter table ",
      "for Eq 8 -- only the printed equation, plus the between-subject ",
      "variability (9.71) in the following paragraph -- so every value in ",
      "this file traces to the equation as typeset, at the 3 significant ",
      "figures it carries. Baseline demographics, concomitant-AED ",
      "prevalences and genotype frequencies are reported for the FULL ",
      "cohort only (Table 1) and are not broken out by age subgroup. The ",
      "authors explicitly caution that this subgroup is the most ",
      "underpowered part of the analysis: 'the patient population of the ",
      "present study was small, especially that of the patients who were 19 ",
      "years of age or older, thus further studies with larger numbers of ",
      "patients are needed' (Discussion, limitation 4). Treat this model as ",
      "the paper's most exploratory result."
    )
  )

  ini({
    # ==================================================================
    # PK LAYER -- ALL FIXED, ALL FROM THE UPSTREAM SOURCE (Ogusu 2014)
    #
    # Identical to the all-age companion; see
    # Nakashima_2015_valproic_acid.R for the full nine-constant table
    # documenting why the upstream Ogusu 2014 values are used in
    # preference to Nakashima 2015's own restatement in its Eqs 1-2
    # (operator ruling, sidecar oare_PMC9833507 q2 = A), and for the
    # physiological falsifier on the restated Vd/F.
    # ==================================================================

    lka   <- fixed(log(0.109))  ; label("Absorption rate constant, from Ogusu 2014 (Ka, 1/h)")                                # Ogusu 2014 Eq 5 and Table 2 (Ka 0.109 1/h)
    lvc   <- fixed(log(21.4))   ; label("Apparent central volume at a 1000 mg/day dose, from Ogusu 2014 (Vd/F, L)")           # Ogusu 2014 Eq 6 and Table 2 (Vd/F 21.4 L)
    lcl   <- fixed(log(0.559))  ; label("Apparent oral clearance at a 1000 mg/day dose in a male on no other AED, from Ogusu 2014 (CL/F, L/h)")  # Ogusu 2014 Eq 7 and Table 2 (CL/F 0.559 L/h)
    ltlag <- fixed(log(3.00))   ; label("Absorption lag time, from Ogusu 2014 (ALAG, h)")                                     # Ogusu 2014 Eq 8 and Table 2 ("3.00 (Fixed)")

    e_dose_vc <- fixed(1.52)    ; label("Power exponent on (DOSE_VPA_MGD/1000) for apparent volume, from Ogusu 2014 (unitless)")    # Ogusu 2014 Eq 6 and Table 2 "Dose on Vd/F"
    e_dose_cl <- fixed(0.596)   ; label("Power exponent on (DOSE_VPA_MGD/1000) for apparent clearance, from Ogusu 2014 (unitless)") # Ogusu 2014 Eq 7 and Table 2 "Dose on CL/F"

    e_sexf_cl <- fixed(0.917)   ; label("Fold-change in apparent clearance for female vs male, from Ogusu 2014 (unitless)")               # Ogusu 2014 Eq 7 and Table 2 "Gender on CL/F"
    e_cbz_cl  <- fixed(1.19)    ; label("Fold-change in apparent clearance with concomitant carbamazepine, from Ogusu 2014 (unitless)")   # Ogusu 2014 Eq 7 and Table 2 "CBZ on CL/F"
    e_pb_cl   <- fixed(1.12)    ; label("Fold-change in apparent clearance with concomitant phenobarbital, from Ogusu 2014 (unitless)")   # Ogusu 2014 Eq 7 and Table 2 "PB on CL/F"
    e_pht_cl  <- fixed(1.43)    ; label("Fold-change in apparent clearance with concomitant phenytoin, from Ogusu 2014 (unitless)")       # Ogusu 2014 Eq 7 and Table 2 "PHT on CL/F"
    e_clb_cl  <- fixed(0.906)   ; label("Fold-change in apparent clearance with concomitant clobazam, from Ogusu 2014 (unitless)")        # Ogusu 2014 Eq 7 and Table 2 "CLB on CL/F"

    etalka   ~ fixed(7.77e-7)   # Ogusu 2014 Table 2 'omega^2 on Ka'
    etalvc   ~ fixed(1.83e-7)   # Ogusu 2014 Table 2 'omega^2 on Vd/F'
    etalcl   ~ fixed(0.0587)    # Ogusu 2014 Table 2 'omega^2 on CL/F'
    etaltlag ~ fixed(4.48e-9)   # Ogusu 2014 Table 2 'omega^2 on ALAG'

    # ==================================================================
    # PD LAYER -- Nakashima 2015 Eq 8 (subjects aged 19 years or older)
    #
    # Eq 8 EXACTLY as printed, verified against the rendered page image
    # and the PDF text layer:
    #
    #   Logit(Pr) = 10.3 - 2.56^CZP - 9.88^{SCN1A GA} - 14.3^{SCN1A AA}
    #             x Predicted trough concentration of VPA           (8)
    #
    # THERE ARE NO PARENTHESES. Eqs 6 and 7 both read
    #   "... - ( base_slope + covariate terms ) x Ctrough"
    # but Eq 8 has no grouping at all, so read literally only the A/A
    # term multiplies the concentration and the G/G and G/A groups get a
    # concentration-INDEPENDENT probability.
    #
    # Why the literal reading is encoded (operator ruling, sidecar
    # oare_PMC9833507 q4 = A -- full fidelity plus errata):
    #
    #  1. The literal reading and the "restore a dropped opening paren"
    #     reading COINCIDE. An opening parenthesis before 9.88 would
    #     require the inner operator to be "+", because the group as a
    #     whole is preceded by "-" (that is exactly how Eqs 6 and 7 are
    #     typeset: the group's interior signs are inverted relative to
    #     the expanded form). The printed operator before 9.88 is "-",
    #     so no consistent parenthesisation puts 9.88 inside the slope
    #     group.
    #  2. There is NO anchor to test Eq 8 against. Table 4's six
    #     simulated patients are aged 5 and 10, so they exercise Eq 6
    #     only; Fig 1's visual predictive check is the all-age model; and
    #     Eq 8 has no parameter table.
    #  3. Eq 8 also has no base slope term of its own -- there is no
    #     e_ctrough_logit here, unlike Eqs 6 and 7 -- which is itself a
    #     consequence of the missing grouping.
    #
    # A user who wants the alternative reading (all three genotype terms
    # inside a slope group, with a base slope that Eq 8 does not print)
    # must supply that base slope themselves; it is not recoverable from
    # any source on disk. See the vignette Errata.
    #
    # The concentration scaling (trough enters as Cc/100) is inherited
    # from Eq 6, where it is established six times over against Table 4;
    # see Nakashima_2015_valproic_acid.R.
    # ==================================================================

    logit_ref     <- 10.3       ; label("Logit of the probability of an over-50% seizure-frequency reduction for an adult with SCN1A G/G and no concomitant clonazepam (unitless logit)")  # Nakashima 2015 Eq 8, intercept 10.3
    e_czp_logit   <- -2.56      ; label("Log-odds shift on the logit intercept with concomitant clonazepam (unitless logit)")   # Nakashima 2015 Eq 8, CZP term -2.56
    e_snp_scn1a_rs3812718_ga_logit <- -9.88 ; label("Log-odds shift on the logit intercept for SCN1A rs3812718 G/A vs G/G (unitless logit)")  # Nakashima 2015 Eq 8, SCN1A G/A term -9.88; NOT a slope term as printed (see the block comment above)

    # The ONLY concentration-dependent term in Eq 8 as printed.
    e_snp_scn1a_rs3812718_aa_slope <- -14.3 ; label("Log-odds per 100 ug/mL of predicted steady-state trough VPA, applying to SCN1A rs3812718 A/A adults only (unitless logit per 100 ug/mL)")  # Nakashima 2015 Eq 8, "- 14.3^{SCN1A AA} x Predicted trough concentration of VPA"; the trailing concentration factor attaches to this term alone because Eq 8 carries no parentheses

    # Between-subject variability on the logit, read as a variance on the
    # logit scale, consistent with the all-age model (see
    # Nakashima_2015_valproic_acid.R for the 12.9 -> 11.3 = -12.4%
    # arithmetic that fixes the scale).
    etalogit_ref ~ 9.71         # Nakashima 2015 Results: 'The variations among individuals of the PK-PD models of the patients aged 18 years or younger and those aged 19 years or older were 11.5 and 9.71, respectively.'

    # No residual error is reported or estimable: the source likelihood is
    # Bernoulli on a binary per-patient outcome. See
    # Nakashima_2015_valproic_acid.R for the full rationale.
    addSd_prob_seizure50 <- fixed(0.001) ; label("Additive residual SD on the typical-value response probability (placeholder; the source likelihood is Bernoulli, so no residual is published)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # ---- PK layer: one compartment, first-order absorption, lag time ----
    # Ogusu 2014 Eqs 5-8, referenced to a 1000 mg/day daily dose.
    dosescale <- DOSE_VPA_MGD / 1000

    ka   <- exp(lka + etalka)
    vc   <- exp(lvc + etalvc) * dosescale^e_dose_vc
    cl   <- exp(lcl + etalcl) * dosescale^e_dose_cl *
      e_sexf_cl^SEXF *
      e_cbz_cl^CONMED_CBZ *
      e_pb_cl^CONMED_PB *
      e_pht_cl^CONMED_PHT *
      e_clb_cl^CONMED_CLB
    tlag <- exp(ltlag + etaltlag)

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    alag(depot) <- tlag

    Cc <- central / vc

    # ---- PD layer: Nakashima 2015 Eq 8, encoded exactly as printed ----
    # Concentration-INDEPENDENT part.
    logit_int <- logit_ref + etalogit_ref +
      e_czp_logit * CONMED_CZP +
      e_snp_scn1a_rs3812718_ga_logit * SNP_SCN1A_RS3812718_GA

    # Concentration-DEPENDENT part: the A/A term only, because Eq 8 is
    # printed without parentheses (see the ini() block comment).
    logit_slope <- e_snp_scn1a_rs3812718_aa_slope * SNP_SCN1A_RS3812718_AA

    logit_seizure50 <- logit_int + logit_slope * (Cc / 100)

    prob_seizure50 <- expit(logit_seizure50)

    prob_seizure50 ~ add(addSd_prob_seizure50)
  })
}
