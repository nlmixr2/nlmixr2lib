Nakashima_2015_valproic_acid_le18y <- function() {
  description <- paste0(
    "Logistic-regression exposure-response model for the probability of an ",
    "over-50% reduction in seizure frequency in Japanese patients with ",
    "epilepsy aged 18 YEARS OR YOUNGER on sustained-release valproic acid ",
    "(VPA), driven by the model-predicted steady-state trough VPA ",
    "concentration (Nakashima 2015 Eq 7; n = 56, mean age 11.7 +/- 4.41 ",
    "years, 597 concentration measurement points). One of two age-subgroup ",
    "analyses the authors ran alongside their all-age final model: relative ",
    "to that model (Nakashima_2015_valproic_acid, Eq 6) this subgroup fit ",
    "drops the age, seizure-locus, clonazepam and topiramate terms, adds a ",
    "phenobarbital effect on the logit intercept and a carbamazepine effect ",
    "on the exposure slope, and carries no parameter table -- every value ",
    "here comes from the printed Eq 7 itself. The slope coefficients are ",
    "per 100 ug/mL, so the trough enters as Cc/100. A one-compartment ",
    "first-order absorption model with lag time supplies the trough; its ",
    "parameters were FIXED into this PD fit from the group's upstream ",
    "population PK analysis (Ogusu 2014 Eqs 5-8 and Table 2), so every PK ",
    "value here is fixed rather than estimated. The paper carries a PLOS ",
    "ONE Expression of Concern on funding-disclosure grounds only (no ",
    "concern about data, methods or results). Companions: ",
    "Nakashima_2015_valproic_acid (Eq 6, all ages) and ",
    "Nakashima_2015_valproic_acid_ge19y (Eq 8, >= 19 y)."
  )
  reference <- paste(
    "Nakashima H, Oniki K, Nishimura M, Ogusu N, Shimomasuda M, Ono T,",
    "Matsuda K, Yasui-Furukori N, Nakagawa K, Ishitsu T, Saruwatari J.",
    "Determination of the Optimal Concentration of Valproic Acid in Patients",
    "with Epilepsy: A Population Pharmacokinetic-Pharmacodynamic Analysis.",
    "PLoS ONE. 2015;10(10):e0141266. doi:10.1371/journal.pone.0141266",
    "(Eq 7, the subgroup aged 18 years or younger).",
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
        "no PD coefficient: Nakashima 2015 screened it on the logit and ",
        "dropped it in backward elimination (Table 2). The subgroup's own ",
        "sex split is not reported separately; the full cohort was 37.7% ",
        "female (Table 1)."
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
        "reported separately; the full cohort was 1120.0 +/- 592.5 ",
        "[50-3200] mg/day (Table 1)."
      ),
      source_name        = "Dose"
    ),
    CONMED_CBZ = list(
      description        = "Concomitant carbamazepine indicator; 1 = co-administered, 0 = not. Acts on the PK layer AND on BOTH the PD logit intercept and the PD exposure slope.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant carbamazepine)",
      notes              = paste0(
        "PK: Ogusu 2014 Table 2 'CBZ on CL/F' = 1.19, power form ",
        "1.19^CONMED_CBZ (+19% CL/F). PD: Nakashima 2015 Eq 7 gives -4.88 ",
        "on the intercept and +4.73 on the slope. Both are much larger in ",
        "magnitude than the all-age Eq 6 values (-1.75 intercept, no slope ",
        "term at all), and carbamazepine is the only covariate that acts on ",
        "both PD blocks in this subgroup. Full-cohort prevalence 47 of 77 ",
        "(61.0%; Table 1); the subgroup's own prevalence is not reported."
      ),
      source_name        = "CBZ"
    ),
    CONMED_PB = list(
      description        = "Concomitant phenobarbital indicator; 1 = co-administered, 0 = not. Acts on the PK layer AND on the PD logit intercept.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant phenobarbital)",
      notes              = paste0(
        "PK: Ogusu 2014 Table 2 'PB on CL/F' = 1.12, power form ",
        "1.12^CONMED_PB (+12% CL/F). PD: Nakashima 2015 Eq 7 gives -1.93 on ",
        "the intercept. This is a term the all-age Eq 6 does NOT have: ",
        "phenobarbital was screened and rejected on the all-age logit ",
        "(Table 2, difference of objective function 0.38 on the intercept ",
        "and 0.00 on the slope) yet is retained here, which is one of the ",
        "structural differences between the age subgroups the authors ",
        "highlight in the Discussion. Full-cohort prevalence 31 of 77 ",
        "(40.3%; Table 1)."
      ),
      source_name        = "PB"
    ),
    CONMED_PHT = list(
      description        = "Concomitant phenytoin indicator; 1 = co-administered, 0 = not. Acts on the PK layer AND on the PD exposure slope.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant phenytoin)",
      notes              = paste0(
        "PK: Ogusu 2014 Table 2 'PHT on CL/F' = 1.43, power form ",
        "1.43^CONMED_PHT (+43% CL/F). PD: Nakashima 2015 Eq 7 gives -3.86 ",
        "on the slope, close to the all-age Eq 6 value of -3.62. ",
        "Full-cohort prevalence 22 of 77 (28.6%; Table 1)."
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
    SNP_SCN1A_RS3812718_GA = list(
      description        = "SCN1A rs3812718 (IVS5-91 G>A) heterozygous G/A genotype indicator; 1 = G/A, 0 = otherwise. Acts on BOTH the PD logit intercept and the PD exposure slope.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 in combination with SNP_SCN1A_RS3812718_AA = 0, i.e. the G/G genotype group",
      notes              = paste0(
        "Nakashima 2015 Eq 7: -4.75 on the intercept, +7.62 on the slope. ",
        "Paired with SNP_SCN1A_RS3812718_AA; both are 0 for the G/G ",
        "reference. Both coefficients are smaller in magnitude than the ",
        "all-age Eq 6 values (-5.87 and +10.1) and much smaller than the ",
        ">= 19 y Eq 8 value (-9.88), which is the age-dependence of the ",
        "genotype effect the authors flag in the Discussion. Full-cohort ",
        "genotype frequencies G/G 11.7%, G/A 53.2%, A/A 35.1% (Table 1); ",
        "the subgroup's own frequencies are not reported."
      ),
      source_name        = "SCN1A G/A genotype"
    ),
    SNP_SCN1A_RS3812718_AA = list(
      description        = "SCN1A rs3812718 (IVS5-91 G>A) homozygous A/A genotype indicator; 1 = A/A, 0 = otherwise. Acts on BOTH the PD logit intercept and the PD exposure slope.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 in combination with SNP_SCN1A_RS3812718_GA = 0, i.e. the G/G genotype group",
      notes              = paste0(
        "Nakashima 2015 Eq 7: -4.30 on the intercept, +7.60 on the slope. ",
        "Paired with SNP_SCN1A_RS3812718_GA; both are 0 for the G/G ",
        "reference. Note that in this subgroup the G/A and A/A slope ",
        "coefficients are nearly identical (+7.62 and +7.60), i.e. the ",
        "effect behaves as A-allele CARRIAGE rather than as an allele-dose ",
        "gradient -- unlike the all-age Eq 6, where they differ (+10.1 vs ",
        "+9.48), and unlike Eq 8. Full-cohort prevalence 27 of 77 (35.1%; ",
        "Table 1)."
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
        "Age DEFINES this subgroup (18 years or younger) rather than acting ",
        "as a covariate within it, so it carries no coefficient in Eq 7 -- ",
        "unlike the all-age Eq 6, where age enters the logit intercept at ",
        "0.98 per year. Subgroup mean age 11.7 +/- 4.41 years (Results, ",
        "Effect of the patient age on the PK-PD modeling). Listed here so ",
        "the subgroup definition is machine-visible without implying an ",
        "unpublished coefficient."
      )
    ),
    CONMED_CZP = list(
      description = "Concomitant clonazepam indicator; 1 = co-administered, 0 = not.",
      units       = "(binary)",
      type        = "binary",
      notes       = paste0(
        "Retained on the logit intercept in the all-age Eq 6 (-1.18) and in ",
        "the >= 19 y Eq 8 (-2.56) but ABSENT from Eq 7, so no coefficient ",
        "exists for this subgroup. Full-cohort prevalence 21 of 77 (27.3%; ",
        "Table 1)."
      )
    ),
    CONMED_TPM = list(
      description = "Concomitant topiramate indicator; 1 = co-administered, 0 = not.",
      units       = "(binary)",
      type        = "binary",
      notes       = paste0(
        "Retained on the exposure slope in the all-age Eq 6 (-1.73) but ",
        "absent from Eq 7, so no coefficient exists for this subgroup. ",
        "Full-cohort prevalence 10 of 77 (13.0%; Table 1)."
      )
    ),
    SEIZURE_LOCUS_PARTIAL = list(
      description = "Seizure-locus indicator; 1 = partial (focal), 0 = generalized.",
      units       = "(binary)",
      type        = "binary",
      notes       = paste0(
        "Retained on the exposure slope in the all-age Eq 6 (+2.41) but ",
        "absent from Eq 7, so no coefficient exists for this subgroup. Note ",
        "that the gloss printed after Eq 8 still defines 'partial seizure' ",
        "even though neither Eq 7 nor Eq 8 uses the term -- that gloss is ",
        "shared boilerplate carried over from Eq 6 and must not be read as ",
        "evidence of an omitted term. Full-cohort split partial 87.0% / ",
        "generalized 13.0% (Table 1)."
      )
    ),
    MENT_DISABLED = list(
      description = "Complication with intellectual disability; 1 = present, 0 = absent.",
      units       = "(binary)",
      type        = "binary",
      notes       = paste0(
        "Screened on the all-age model and dropped in backward elimination ",
        "(Nakashima 2015 Table 2 and Results); never appears in Eq 7. No ",
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
    n_subjects     = 56L,
    n_studies      = 1L,
    n_observations = "597 VPA concentration measurement points; one binary efficacy outcome per patient",
    age_range      = "0.8-18 years (subgroup upper bound is the 18-year inclusion cut-off; cohort minimum 0.8 years)",
    age_median     = "mean 11.7 +/- 4.41 years (median not reported)",
    sex_female_pct = NA_real_,
    race_ethnicity = c(Asian = 100),
    disease_state  = paste0(
      "Epilepsy on maintenance sustained-release valproic acid, paediatric ",
      "and adolescent subgroup (aged 18 years or younger). Idiopathic ",
      "epilepsy and severe retardation were exclusion criteria."
    ),
    dose_range     = "not reported for this subgroup; full cohort 1120.0 +/- 592.5 [50-3200] mg/day sustained-release valproic acid",
    regions        = "Japan (Kumamoto Saishunso National Hospital), June 1989 to April 2011",
    notes          = paste0(
      "Subgroup of the 77-patient cohort behind Nakashima 2015 Eq 6, split ",
      "at the 18/19-year boundary (Results, Effect of the patient age on ",
      "the PK-PD modeling). Nakashima 2015 reports NO parameter table for ",
      "Eq 7 -- only the printed equation, plus the between-subject ",
      "variability (11.5) in the following paragraph -- so every value in ",
      "this file traces to the equation as typeset, at the 3 significant ",
      "figures it carries. Baseline demographics, concomitant-AED ",
      "prevalences and genotype frequencies are reported for the FULL ",
      "cohort only (Table 1) and are not broken out by age subgroup. The ",
      "authors caution that the subgroup analyses are underpowered and that ",
      "the differences between Eqs 7 and 8 need confirmation in a larger ",
      "cohort (Discussion, limitation 4)."
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
    # physiological falsifier on the restated Vd/F. Nakashima 2015
    # Methods: "The individual PK parameters, determined from our
    # previously reported population PK model [8], were fixed in the
    # PK-PD analysis."
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
    # PD LAYER -- Nakashima 2015 Eq 7 (subjects aged 18 years or younger)
    #
    # Eq 7 as printed:
    #   logit(Pr) = 7.73 - 4.88^CBZ - 1.93^PB - 4.75^{SCN1A GA}
    #                    - 4.30^{SCN1A AA}
    #             - (10.9 - 4.73^CBZ + 3.86^PHT - 7.62^{SCN1A GA}
    #                - 7.60^{SCN1A AA}) x Predicted trough concentration
    #
    # (the paper typesets the indicator variables as superscripts; the
    # Eq 6 / Eq 8 gloss "otherwise the term was assigned a value of 0"
    # confirms these are additive terms, not exponents). Distributing the
    # leading minus over the parenthesised group gives the signed slope
    # coefficients used below: -10.9, +4.73, -3.86, +7.62, +7.60. Note
    # the sign flip on CBZ -- it is the one covariate printed with a
    # MINUS inside the group, so it is the one whose slope coefficient
    # comes out POSITIVE.
    #
    # Unlike Eq 6, Eq 7 has no parameter table, so there is no
    # higher-precision column to prefer and no independent anchor to
    # test it against (Table 4's six simulated patients are aged 5 and
    # 10 but were computed from Eq 6). The concentration scaling is
    # inherited from Eq 6, where it is established six times over
    # against Table 4: the trough enters as Cc/100. The three base
    # slopes are mutually consistent with a single shared scaling
    # (Eq 6 -13.5, Eq 7 -10.9, Eq 8 -14.3 per 100 ug/mL).
    # ==================================================================

    logit_ref     <- 7.73       ; label("Logit of the probability of an over-50% seizure-frequency reduction at zero trough, for a patient aged 18 y or younger with SCN1A G/G and no concomitant AED (unitless logit)")  # Nakashima 2015 Eq 7, intercept 7.73
    e_cbz_logit   <- -4.88      ; label("Log-odds shift on the logit intercept with concomitant carbamazepine (unitless logit)")   # Nakashima 2015 Eq 7, intercept block, CBZ term -4.88
    e_pb_logit    <- -1.93      ; label("Log-odds shift on the logit intercept with concomitant phenobarbital (unitless logit)")   # Nakashima 2015 Eq 7, intercept block, PB term -1.93
    e_snp_scn1a_rs3812718_ga_logit <- -4.75 ; label("Log-odds shift on the logit intercept for SCN1A rs3812718 G/A vs G/G (unitless logit)")  # Nakashima 2015 Eq 7, intercept block, SCN1A G/A term -4.75
    e_snp_scn1a_rs3812718_aa_logit <- -4.30 ; label("Log-odds shift on the logit intercept for SCN1A rs3812718 A/A vs G/G (unitless logit)")  # Nakashima 2015 Eq 7, intercept block, SCN1A A/A term -4.30

    e_ctrough_logit <- -10.9    ; label("Log-odds of an over-50% seizure-frequency reduction per 100 ug/mL of predicted steady-state trough VPA, in the reference patient (unitless logit)")  # Nakashima 2015 Eq 7, slope group base term 10.9 (negated by the leading minus)
    e_cbz_slope     <- 4.73     ; label("Shift in the trough exposure slope with concomitant carbamazepine (unitless logit per 100 ug/mL)")  # Nakashima 2015 Eq 7, slope group "- 4.73^CBZ" inside a group preceded by a minus, hence +4.73
    e_pht_slope     <- -3.86    ; label("Shift in the trough exposure slope with concomitant phenytoin (unitless logit per 100 ug/mL)")      # Nakashima 2015 Eq 7, slope group "+ 3.86^PHT" inside a group preceded by a minus, hence -3.86
    e_snp_scn1a_rs3812718_ga_slope <- 7.62 ; label("Shift in the trough exposure slope for SCN1A rs3812718 G/A vs G/G (unitless logit per 100 ug/mL)")  # Nakashima 2015 Eq 7, slope group "- 7.62^{SCN1A GA}", hence +7.62
    e_snp_scn1a_rs3812718_aa_slope <- 7.60 ; label("Shift in the trough exposure slope for SCN1A rs3812718 A/A vs G/G (unitless logit per 100 ug/mL)")  # Nakashima 2015 Eq 7, slope group "- 7.60^{SCN1A AA}", hence +7.60

    # Between-subject variability on the logit, read as a variance on the
    # logit scale, consistent with the all-age model (see
    # Nakashima_2015_valproic_acid.R for the 12.9 -> 11.3 = -12.4%
    # arithmetic that fixes the scale).
    etalogit_ref ~ 11.5         # Nakashima 2015 Results: 'The variations among individuals of the PK-PD models of the patients aged 18 years or younger and those aged 19 years or older were 11.5 and 9.71, respectively.'

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

    # ---- PD layer: Nakashima 2015 Eq 7 ----
    logit_int <- logit_ref + etalogit_ref +
      e_cbz_logit * CONMED_CBZ +
      e_pb_logit  * CONMED_PB +
      e_snp_scn1a_rs3812718_ga_logit * SNP_SCN1A_RS3812718_GA +
      e_snp_scn1a_rs3812718_aa_logit * SNP_SCN1A_RS3812718_AA

    logit_slope <- e_ctrough_logit +
      e_cbz_slope * CONMED_CBZ +
      e_pht_slope * CONMED_PHT +
      e_snp_scn1a_rs3812718_ga_slope * SNP_SCN1A_RS3812718_GA +
      e_snp_scn1a_rs3812718_aa_slope * SNP_SCN1A_RS3812718_AA

    logit_seizure50 <- logit_int + logit_slope * (Cc / 100)

    prob_seizure50 <- expit(logit_seizure50)

    prob_seizure50 ~ add(addSd_prob_seizure50)
  })
}
