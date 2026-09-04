Nakashima_2015_valproic_acid <- function() {
  description <- paste0(
    "Logistic-regression exposure-response model for the probability of an ",
    "over-50% reduction in seizure frequency in Japanese patients with ",
    "epilepsy on sustained-release valproic acid (VPA), driven by the ",
    "model-predicted steady-state trough VPA concentration (Nakashima 2015 ",
    "Eq 6, the paper's FINAL all-age model; n = 77, age 0.8-36.9 years, 729 ",
    "concentration measurement points). logit(Pr) is linear in the ",
    "predicted trough: an intercept carrying age, carbamazepine, ",
    "clonazepam and SCN1A rs3812718 genotype effects, plus a slope carrying ",
    "seizure locus, phenytoin, topiramate and SCN1A genotype effects ",
    "(Nakashima 2015 Table 3). The published slope coefficients are per 100 ",
    "ug/mL, so the trough enters as Cc/100. A one-compartment first-order ",
    "absorption model with lag time supplies the trough; its parameters ",
    "were FIXED into this PD fit from the group's upstream population PK ",
    "analysis (Ogusu 2014 Eqs 5-8 and Table 2), so every PK value here is ",
    "fixed rather than estimated. NOTE: Nakashima 2015 restates that PK ",
    "model in its own Eqs 1-2 and disagrees with Ogusu 2014 on all nine ",
    "constants; the upstream source values are used here per operator ",
    "ruling -- see the vignette Errata. The paper carries a PLOS ONE ",
    "Expression of Concern on funding-disclosure grounds only (no concern ",
    "about data, methods or results). Age-subgroup companions in ",
    "Nakashima_2015_valproic_acid_le18y (Eq 7) and ",
    "Nakashima_2015_valproic_acid_ge19y (Eq 8)."
  )
  reference <- paste(
    "Nakashima H, Oniki K, Nishimura M, Ogusu N, Shimomasuda M, Ono T,",
    "Matsuda K, Yasui-Furukori N, Nakagawa K, Ishitsu T, Saruwatari J.",
    "Determination of the Optimal Concentration of Valproic Acid in Patients",
    "with Epilepsy: A Population Pharmacokinetic-Pharmacodynamic Analysis.",
    "PLoS ONE. 2015;10(10):e0141266. doi:10.1371/journal.pone.0141266.",
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
    AGE = list(
      description        = "Patient age. Enters the logit intercept linearly, at 0.98 per YEAR.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Nakashima 2015 Eq 6 prints the age term as (Age/10) x 1.0, but ",
        "Table 3 labels the row 'Age (years)' with estimate 0.98. The two ",
        "readings differ 10-fold and Table 4 settles it: inverting Eq 6 at ",
        "the ROC cut-off logit(Pr) = 0.1 for the six simulated patients in ",
        "Table 4 reproduces all six printed optimal troughs to <= 0.30% ",
        "with 0.98 PER YEAR, and misses by 88.9% with the Age/10 reading. ",
        "The printed /10 is therefore a typesetting artefact. Cohort ",
        "15.2 +/- 8.2 [0.8-36.9] years (Table 1)."
      ),
      source_name        = "Age"
    ),
    SEXF = list(
      description        = "Sex indicator; 1 = female, 0 = male. Acts on the PK layer only (apparent oral clearance).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste0(
        "Ogusu 2014 Eq 7 / Table 2 'Gender on CL/F' = 0.917, applied as a ",
        "power form 0.917^SEXF (i.e. females have 8.3% lower CL/F). Sex was ",
        "screened by Nakashima 2015 on the PD logit and dropped in backward ",
        "elimination (Table 2), so it carries no PD coefficient here. ",
        "Cohort 29 of 77 female (37.7%; Table 1)."
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
        "reference. Must be strictly positive. Distinct from the per-",
        "administration dose amount in the event table: this covariate is ",
        "the DAILY total, which equals the per-administration amt only for ",
        "once-daily dosing. Cohort 1120.0 +/- 592.5 [50-3200] mg/day of ",
        "sustained-release VPA (Table 1)."
      ),
      source_name        = "Dose"
    ),
    CONMED_CBZ = list(
      description        = "Concomitant carbamazepine indicator; 1 = co-administered, 0 = not. Acts on BOTH the PK layer (induction of VPA clearance) and the PD logit intercept.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant carbamazepine)",
      notes              = paste0(
        "PK: Ogusu 2014 Table 2 'CBZ on CL/F' = 1.19, power form ",
        "1.19^CONMED_CBZ (+19% CL/F). PD: Nakashima 2015 Table 3 intercept ",
        "row 'CBZ' = -1.75 log-odds. Cohort 47 of 77 (61.0%; Table 1)."
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
        "1.12^CONMED_PB (+12% CL/F). Phenobarbital was screened on the PD ",
        "logit by Nakashima 2015 and was not significant (Table 2, DOBF ",
        "0.38 intercept / 0.00 slope), so it carries no PD coefficient in ",
        "Eq 6; it DOES carry one in the <= 18 y companion (Eq 7). Cohort 31 ",
        "of 77 (40.3%; Table 1)."
      ),
      source_name        = "PB"
    ),
    CONMED_PHT = list(
      description        = "Concomitant phenytoin indicator; 1 = co-administered, 0 = not. Acts on BOTH the PK layer (induction of VPA clearance) and the PD exposure slope.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant phenytoin)",
      notes              = paste0(
        "PK: Ogusu 2014 Table 2 'PHT on CL/F' = 1.43, power form ",
        "1.43^CONMED_PHT (+43% CL/F, the largest induction coefficient in ",
        "that model). PD: Nakashima 2015 Table 3 slope row 'PHT' = -3.62, ",
        "i.e. phenytoin makes the exposure-response slope more negative. ",
        "Cohort 22 of 77 (28.6%; Table 1)."
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
        "0.906^CONMED_CLB (-9.4% CL/F; the only AED in that model that ",
        "LOWERS VPA clearance). Clobazam was screened on the PD logit by ",
        "Nakashima 2015 and was not significant (Table 2, DOBF 0.27 / ",
        "0.48). Cohort 33 of 77 (42.9%; Table 1)."
      ),
      source_name        = "CLB"
    ),
    CONMED_CZP = list(
      description        = "Concomitant clonazepam indicator; 1 = co-administered, 0 = not. Acts on the PD logit intercept only.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant clonazepam)",
      notes              = paste0(
        "Nakashima 2015 Table 3 intercept row 'CZP' = -1.18 log-odds. Not ",
        "a PK covariate: clonazepam was not retained in the Ogusu 2014 ",
        "clearance model. The authors state the mechanism for the CZP ",
        "association is unclear (Discussion). Cohort 21 of 77 (27.3%; ",
        "Table 1)."
      ),
      source_name        = "CZP"
    ),
    CONMED_TPM = list(
      description        = "Concomitant topiramate indicator; 1 = co-administered, 0 = not. Acts on the PD exposure slope only.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant topiramate)",
      notes              = paste0(
        "Nakashima 2015 Table 3 slope row 'TPM' = -1.73. Not a PK ",
        "covariate: topiramate was not retained in the Ogusu 2014 clearance ",
        "model. Cohort 10 of 77 (13.0%; Table 1)."
      ),
      source_name        = "TPM"
    ),
    SEIZURE_LOCUS_PARTIAL = list(
      description        = "Seizure-locus indicator from the ILAE classification; 1 = partial (focal) seizures present, 0 = generalized. Acts on the PD exposure slope.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (generalized seizure locus)",
      notes              = paste0(
        "Nakashima 2015 Table 3 slope row 'Partial seizure' = +2.41, i.e. ",
        "partial-onset patients have a LESS negative exposure-response ",
        "slope than generalized patients. Eq 6's gloss defines it: 'partial ",
        "seizure = 1 if a partial seizure was present, and 0 if partial ",
        "seizures were absent'. Seizures were classified per the ILAE ",
        "guidelines (Methods, Definitions). Cohort: partial 67 of 77 ",
        "(87.0%), generalized 10 (13.0%) (Table 1). Note the cohort is ",
        "heavily partial-dominant, so the generalized reference group is ",
        "only 10 patients -- the +2.41 estimate has a wide bootstrap CI ",
        "(-0.53 to 7.82, Table 3) that spans zero."
      ),
      source_name        = "partial seizure"
    ),
    SNP_SCN1A_RS3812718_GA = list(
      description        = "SCN1A rs3812718 (IVS5-91 G>A) heterozygous G/A genotype indicator; 1 = G/A, 0 = otherwise. Acts on BOTH the PD logit intercept and the PD exposure slope.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 in combination with SNP_SCN1A_RS3812718_AA = 0, i.e. the G/G genotype group (9 of 77)",
      notes              = paste0(
        "Nakashima 2015 Table 3: intercept -5.87, slope +10.1. Paired with ",
        "SNP_SCN1A_RS3812718_AA; both are 0 for the G/G reference. ",
        "Genotyped by custom TaqMan allelic discrimination (Methods, ",
        "Genotyping). rs3812718 sits in the splice-donor site of SCN1A exon ",
        "5 and alters exon-5N / exon-5A alternative splicing. Cohort G/G ",
        "11.7%, G/A 53.2%, A/A 35.1%; A-allele frequency 61.7%; consistent ",
        "with Hardy-Weinberg equilibrium (Results, Table 1)."
      ),
      source_name        = "SCN1A G/A genotype"
    ),
    SNP_SCN1A_RS3812718_AA = list(
      description        = "SCN1A rs3812718 (IVS5-91 G>A) homozygous A/A genotype indicator; 1 = A/A, 0 = otherwise. Acts on BOTH the PD logit intercept and the PD exposure slope.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 in combination with SNP_SCN1A_RS3812718_GA = 0, i.e. the G/G genotype group (9 of 77)",
      notes              = paste0(
        "Nakashima 2015 Table 3: intercept -4.88, slope +9.48. Paired with ",
        "SNP_SCN1A_RS3812718_GA; both are 0 for the G/G reference. The ",
        "authors read the A allele as associated with reduced VPA efficacy ",
        "(Discussion), consistent with their earlier finding that A/A was ",
        "associated with carbamazepine-resistant epilepsy. Cohort 27 of 77 ",
        "(35.1%; Table 1)."
      ),
      source_name        = "SCN1A A/A genotype"
    )
  )

  covariatesDataExcluded <- list(
    MENT_DISABLED = list(
      description = "Complication with intellectual disability; 1 = present, 0 = absent.",
      units       = "(binary)",
      type        = "binary",
      notes       = paste0(
        "Screened by Nakashima 2015 on both the logit intercept and the ",
        "slope and significant in FORWARD inclusion (Table 2, DOBF 11.84 ",
        "and 12.45) but dropped in BACKWARD elimination (DOBF 0.01 and ",
        "1.27, both NS); the Results state explicitly that 'the patient ",
        "gender and complication with intellectual disability were removed ",
        "from the full covariate model'. No point estimate is published, so ",
        "no coefficient exists to encode. Cohort 59 of 77 (76.6%; Table 1). ",
        "Documented here to preserve the covariate screen without carrying ",
        "an unused-covariate convention warning."
      )
    ),
    SEIZURE_TYPE_SYMPTOMATIC = list(
      description = "Seizure aetiology; 1 = symptomatic, 0 = cryptogenic. NOT the same axis as SEIZURE_LOCUS_PARTIAL.",
      units       = "(binary)",
      type        = "binary",
      notes       = paste0(
        "Screened by Nakashima 2015 and not significant at either step ",
        "(Table 2 'Seizure type', DOBF 3.05 intercept and 0.67 slope, both ",
        "NS), so it never entered the model and no point estimate exists. ",
        "Cohort symptomatic 38 of 77 (49.4%), cryptogenic 39 (50.6%) ",
        "(Table 1). Listed here because the paper tabulates and screens ",
        "seizure TYPE separately from seizure LOCUS, and conflating the two ",
        "would mis-encode SEIZURE_LOCUS_PARTIAL. Name is provisional: it is ",
        "documentation only and is not registered as a canonical, because ",
        "no model uses it."
      )
    )
  )

  compartmentData <- list(
    depot   = list(analyte = "valproic acid", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "valproic acid", units = "mg", specimen = "serum", verified = TRUE)
  )

  # The observation variable `prob_seizure50` (probability of an over-50%
  # reduction in seizure frequency) is an algebraic output, not an ODE state.
  # It is registered as a canonical member of the established `prob_<endpoint>`
  # PD-output family in `inst/references/compartment-names.md`, alongside
  # `prob_roc`, `prob_scc` and the `Fukae_2024_valemetostat_*` outputs. The
  # `50` suffix is definitional -- a 75%-responder or seizure-freedom endpoint
  # would be a different canonical.

  population <- list(
    species        = "human",
    n_subjects     = 77L,
    n_studies      = 1L,
    n_observations = "729 VPA concentration measurement points; one binary efficacy outcome per patient",
    age_range      = "0.8-36.9 years",
    age_median     = "mean 15.2 +/- 8.2 years (median not reported)",
    sex_female_pct = 37.7,
    race_ethnicity = c(Asian = 100),
    disease_state  = paste0(
      "Epilepsy on maintenance sustained-release valproic acid; seizure ",
      "locus partial 87.0% / generalized 13.0%; seizure type symptomatic ",
      "49.4% / cryptogenic 50.6%; intellectual disability 76.6%. Idiopathic ",
      "epilepsy and severe retardation were exclusion criteria."
    ),
    dose_range     = "1120.0 +/- 592.5 [50-3200] mg/day sustained-release valproic acid",
    regions        = "Japan (Kumamoto Saishunso National Hospital), June 1989 to April 2011",
    notes          = paste0(
      "Retrospective single-centre cohort; the same patients as the ",
      "upstream Ogusu 2014 population PK analysis (Nakashima 2015 Methods: ",
      "'conducted among the same epileptic patients reported in our ",
      "previous study'). Mean follow-up 5.1 +/- 4.4 [0-14.5] years. ",
      "Predicted steady-state trough VPA 69.3 +/- 19.9 [11.8-130.1] ug/mL ",
      "(Table 1; note the Table 1 row mislabels the units as 'mg/day' -- ",
      "Table 4 and the Fig 1 x-axis both give ug/ml). Concomitant AEDs: ",
      "CBZ 61.0%, CZP 27.3%, CLB 42.9%, GBP 10.4%, LTG 3.9%, PB 40.3%, ",
      "PHT 28.6%, TPM 13.0%, ZNS 33.8%. Efficacy endpoint: over-50% ",
      "reduction in seizure frequency versus the pre-VPA baseline, assessed ",
      "over at least five times each patient's seizure interval. The ROC ",
      "curve for logit(Pr) had AUC 0.823 (95% CI 0.793-0.853); the optimal ",
      "cut-off was logit(Pr) = 0.1 (sensitivity 71.8%, specificity 80.4%). ",
      "Upstream PK cohort (Ogusu 2014): 237 patients, 827 steady-state VPA ",
      "concentrations."
    )
  )

  ini({
    # ==================================================================
    # PK LAYER -- ALL FIXED, ALL FROM THE UPSTREAM SOURCE (Ogusu 2014)
    #
    # Nakashima 2015 Methods, PK-PD modeling: "The individual PK
    # parameters, determined from our previously reported population PK
    # model [8], were fixed in the PK-PD analysis." Reference [8] is
    # Ogusu 2014. Every PK value below is therefore fixed(), not
    # estimated, in this model.
    #
    # WHY OGUSU 2014 AND NOT NAKASHIMA 2015's OWN Eqs 1-2:
    # Nakashima 2015 restates the inherited PK model in its Eqs 1-2 and
    # disagrees with Ogusu 2014 on ALL NINE constants:
    #
    #   constant           Ogusu 2014 (source)   Nakashima 2015 (restated)
    #   Vd/F coefficient   21.4                  110
    #   Vd/F dose exponent 1.52                  1.51
    #   CL/F coefficient   0.559                 0.577
    #   CL/F dose exponent 0.596                 0.535
    #   female on CL/F     0.917                 0.875
    #   CBZ on CL/F        1.19                  1.22
    #   PB on CL/F         1.12                  1.10
    #   PHT on CL/F        1.43                  1.40
    #   CLB on CL/F        0.906                 0.915
    #
    # These are not rounding, and are not Ogusu's bootstrap medians.
    # Ogusu 2014's printed Eqs 5-8 and its Table 2 agree with each other
    # EXACTLY, so the source is internally consistent and the
    # restatement agrees with neither. Physiological falsifier: the
    # restated Vd/F = 110 * (Dose/1000)^1.51 is 110-130 L at this
    # cohort's doses, i.e. about 2.8 L/kg, impossible for valproic acid
    # (V/F about 0.1-0.4 L/kg); Ogusu's 21.4 L is about 0.5 L/kg.
    # Operator ruling (sidecar oare_PMC9833507 q2 = A): use the upstream
    # source values and document the misquotation. See vignette Errata.
    # ==================================================================

    lka   <- fixed(log(0.109))  ; label("Absorption rate constant, from Ogusu 2014 (Ka, 1/h)")                                # Ogusu 2014 Eq 5 and Table 2 (Ka 0.109 1/h)
    lvc   <- fixed(log(21.4))   ; label("Apparent central volume at a 1000 mg/day dose, from Ogusu 2014 (Vd/F, L)")           # Ogusu 2014 Eq 6 and Table 2 (Vd/F 21.4 L)
    lcl   <- fixed(log(0.559))  ; label("Apparent oral clearance at a 1000 mg/day dose in a male on no other AED, from Ogusu 2014 (CL/F, L/h)")  # Ogusu 2014 Eq 7 and Table 2 (CL/F 0.559 L/h)
    ltlag <- fixed(log(3.00))   ; label("Absorption lag time, from Ogusu 2014 (ALAG, h)")                                     # Ogusu 2014 Eq 8 and Table 2 ("3.00 (Fixed)" -- fixed in the upstream fit too)

    e_dose_vc <- fixed(1.52)    ; label("Power exponent on (DOSE_VPA_MGD/1000) for apparent volume, from Ogusu 2014 (unitless)")    # Ogusu 2014 Eq 6 and Table 2 "Dose on Vd/F"
    e_dose_cl <- fixed(0.596)   ; label("Power exponent on (DOSE_VPA_MGD/1000) for apparent clearance, from Ogusu 2014 (unitless)") # Ogusu 2014 Eq 7 and Table 2 "Dose on CL/F"

    # Power-form categorical effects on CL/F: applied as <ratio>^INDICATOR,
    # exactly the form Ogusu 2014 Eq 7 prints (0.917^female etc.), so each
    # value below is the published fold-change on CL/F.
    e_sexf_cl <- fixed(0.917)   ; label("Fold-change in apparent clearance for female vs male, from Ogusu 2014 (unitless)")               # Ogusu 2014 Eq 7 and Table 2 "Gender on CL/F"
    e_cbz_cl  <- fixed(1.19)    ; label("Fold-change in apparent clearance with concomitant carbamazepine, from Ogusu 2014 (unitless)")   # Ogusu 2014 Eq 7 and Table 2 "CBZ on CL/F"
    e_pb_cl   <- fixed(1.12)    ; label("Fold-change in apparent clearance with concomitant phenobarbital, from Ogusu 2014 (unitless)")   # Ogusu 2014 Eq 7 and Table 2 "PB on CL/F"
    e_pht_cl  <- fixed(1.43)    ; label("Fold-change in apparent clearance with concomitant phenytoin, from Ogusu 2014 (unitless)")       # Ogusu 2014 Eq 7 and Table 2 "PHT on CL/F"
    e_clb_cl  <- fixed(0.906)   ; label("Fold-change in apparent clearance with concomitant clobazam, from Ogusu 2014 (unitless)")        # Ogusu 2014 Eq 7 and Table 2 "CLB on CL/F"

    # Upstream PK IIV variances (Ogusu 2014 Table 2, NONMEM Final
    # Estimates column; the rows are labelled omega^2). The upstream fit
    # collapsed IIV on Ka, Vd/F and ALAG to numerically negligible values
    # (1e-7 to 1e-9) and retained meaningful IIV only on CL/F. These are
    # reproduced as published rather than rounded to zero, so that
    # zeroRe() and the published variance structure are both available;
    # note the bootstrap medians in Table 2 are much larger (8.82e-4,
    # 4.28e-4, 6.69e-5), which is the usual signature of a variance
    # pinned against its lower bound.
    etalka   ~ fixed(7.77e-7)   # Ogusu 2014 Table 2 "omega^2 on Ka"
    etalvc   ~ fixed(1.83e-7)   # Ogusu 2014 Table 2 "omega^2 on Vd/F"
    etalcl   ~ fixed(0.0587)    # Ogusu 2014 Table 2 "omega^2 on CL/F"
    etaltlag ~ fixed(4.48e-9)   # Ogusu 2014 Table 2 "omega^2 on ALAG"

    # ==================================================================
    # PD LAYER -- Nakashima 2015 Eq 6, values from Table 3
    #
    # Eq 6 as printed groups the exposure terms:
    #   logit(Pr) = 6.1 + (Age/10)*1.0 - 1.8^CBZ - 1.2^CZP - 5.9^GA - 4.9^AA
    #             - (13.3 + 3.6^PHT + 1.7^TPM - 2.4^partial - 10.1^GA
    #                - 9.5^AA) x Predicted trough concentration of VPA
    # (the paper typesets the indicator variables as superscripts; the Eq 6
    # gloss "otherwise the term was assigned a value of 0" confirms these
    # are additive terms, not exponents).
    #
    # Distributing the leading minus over the parenthesised group gives
    # the signed slope coefficients -13.3, -3.6, -1.7, +2.4, +10.1, +9.5,
    # which is exactly the sign pattern of Table 3's Slope block. Table 3
    # is used here because it carries more significant figures, and
    # because its base slope -13.5 (not Eq 6's -13.3) is the value that
    # reproduces Table 4 -- see the source-trace note on e_ctrough_logit.
    # ==================================================================

    logit_ref     <- 6.09       ; label("Logit of the probability of an over-50% seizure-frequency reduction at zero trough, for a newborn male with generalized seizures, SCN1A G/G and no concomitant AED (unitless logit)")  # Nakashima 2015 Table 3, Intercept, NONMEM Estimate 6.09 (SE 2.3)
    e_age_logit   <- 0.98       ; label("Log-odds shift on the logit intercept per YEAR of age (unitless logit)")                     # Nakashima 2015 Table 3, Intercept block, "Age (years)" 0.98 (SE 0.41); the /10 printed in Eq 6 is falsified by Table 4 (see covariateData$AGE notes)
    e_cbz_logit   <- -1.75      ; label("Log-odds shift on the logit intercept with concomitant carbamazepine (unitless logit)")      # Nakashima 2015 Table 3, Intercept block, "CBZ" -1.75 (SE 0.50)
    e_czp_logit   <- -1.18      ; label("Log-odds shift on the logit intercept with concomitant clonazepam (unitless logit)")         # Nakashima 2015 Table 3, Intercept block, "CZP" -1.18 (SE 0.66)
    e_snp_scn1a_rs3812718_ga_logit <- -5.87 ; label("Log-odds shift on the logit intercept for SCN1A rs3812718 G/A vs G/G (unitless logit)")  # Nakashima 2015 Table 3, Intercept block, "SCN1A G/A genotype" -5.87 (SE 2.53)
    e_snp_scn1a_rs3812718_aa_logit <- -4.88 ; label("Log-odds shift on the logit intercept for SCN1A rs3812718 A/A vs G/G (unitless logit)")  # Nakashima 2015 Table 3, Intercept block, "SCN1A A/A genotype" -4.88 (SE 2.56)

    # The exposure slope. All slope coefficients are PER 100 ug/mL: the
    # trough enters model() as Cc/100. Back-solving the concentration
    # scale from each of Table 4's six printed optimal troughs returns
    # 0.010017, 0.010029, 0.010021, 0.010025, 0.010029, 0.010028 -- i.e.
    # 1/100 to within the printed rounding, six times independently.
    # Without the /100 the logit would swing by about 1600 over the
    # cohort's observed 11.8-130.1 ug/mL trough range.
    e_ctrough_logit <- -13.5    ; label("Log-odds of an over-50% seizure-frequency reduction per 100 ug/mL of predicted steady-state trough VPA, in the reference patient (unitless logit)")  # Nakashima 2015 Table 3, Slope, NONMEM Estimate -13.5 (SE 2.3). Eq 6 prints 13.3; Table 3's -13.5 reproduces all six Table 4 troughs to <= 0.29% where -13.3 gives 3.24%
    e_partial_slope <- 2.41     ; label("Shift in the trough exposure slope for a partial vs generalized seizure locus (unitless logit per 100 ug/mL)")   # Nakashima 2015 Table 3, Slope block, "Partial seizure" 2.41 (SE 1.12)
    e_pht_slope     <- -3.62    ; label("Shift in the trough exposure slope with concomitant phenytoin (unitless logit per 100 ug/mL)")                   # Nakashima 2015 Table 3, Slope block, "PHT" -3.62 (SE 1.12)
    e_tpm_slope     <- -1.73    ; label("Shift in the trough exposure slope with concomitant topiramate (unitless logit per 100 ug/mL)")                  # Nakashima 2015 Table 3, Slope block, "TPM" -1.73 (SE 1.12)
    e_snp_scn1a_rs3812718_ga_slope <- 10.1  ; label("Shift in the trough exposure slope for SCN1A rs3812718 G/A vs G/G (unitless logit per 100 ug/mL)")   # Nakashima 2015 Table 3, Slope block, "SCN1A G/A genotype" 10.1 (SE 3.17)
    e_snp_scn1a_rs3812718_aa_slope <- 9.48  ; label("Shift in the trough exposure slope for SCN1A rs3812718 A/A vs G/G (unitless logit per 100 ug/mL)")   # Nakashima 2015 Table 3, Slope block, "SCN1A A/A genotype" 9.48 (SE 3.03)

    # Between-subject variability on the logit. Eq 4 places the random
    # effect additively on the logit scale (logit(Pr) = Intercept +
    # trough x Slope + eta), and Table 3's "Individual random effect"
    # row is 11.3 (SE 1.37). Read as a variance on the logit scale
    # (SD about 3.4 logit units): the Results state the base model was
    # 12.9 and the final model 11.3, "a 12.4% decrease", and
    # (12.9 - 11.3)/12.9 = 12.4% confirms the two numbers are the same
    # quantity on the same scale.
    etalogit_ref ~ 11.3         # Nakashima 2015 Table 3, "Individual random effect" 11.3 (SE 1.37); Results narrative cross-check 12.9 -> 11.3 = -12.4%

    # No residual error is reported or estimable: the source likelihood is
    # Bernoulli on a binary per-patient outcome (NONMEM LAPLACIAN), so
    # there is no sigma to transcribe. A small additive placeholder is
    # attached to the probability so the model has a valid observation
    # equation for simulation. The upstream PK residual (Ogusu 2014
    # Table 2, sigma^2 = 0.0617 proportional, i.e. 24.8% CV) is NOT a
    # parameter of this model -- the PK parameters were fixed into this
    # fit, so the PK residual does not enter the PD likelihood.
    addSd_prob_seizure50 <- fixed(0.001) ; label("Additive residual SD on the typical-value response probability (placeholder; the source likelihood is Bernoulli, so no residual is published)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # ---- PK layer: one compartment, first-order absorption, lag time ----
    # Ogusu 2014 Eqs 5-8. Dose scaling is on the DAILY dose in mg/day,
    # referenced to 1000 mg/day.
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

    # Serum VPA concentration in ug/mL (mg / L). Read at the end of a
    # steady-state dosing interval this is the paper's "predicted trough
    # concentration of VPA", which is the PD driver.
    Cc <- central / vc

    # ---- PD layer: Nakashima 2015 Eq 6 ----
    # Intercept block (Table 3, Intercept rows).
    logit_int <- logit_ref + etalogit_ref +
      e_age_logit * AGE +
      e_cbz_logit * CONMED_CBZ +
      e_czp_logit * CONMED_CZP +
      e_snp_scn1a_rs3812718_ga_logit * SNP_SCN1A_RS3812718_GA +
      e_snp_scn1a_rs3812718_aa_logit * SNP_SCN1A_RS3812718_AA

    # Slope block (Table 3, Slope rows), per 100 ug/mL.
    logit_slope <- e_ctrough_logit +
      e_partial_slope * SEIZURE_LOCUS_PARTIAL +
      e_pht_slope     * CONMED_PHT +
      e_tpm_slope     * CONMED_TPM +
      e_snp_scn1a_rs3812718_ga_slope * SNP_SCN1A_RS3812718_GA +
      e_snp_scn1a_rs3812718_aa_slope * SNP_SCN1A_RS3812718_AA

    # The published slope coefficients are per 100 ug/mL (see the
    # e_ctrough_logit source-trace comment).
    logit_seizure50 <- logit_int + logit_slope * (Cc / 100)

    # Eq 3: Pr = exp(logit)/(1 + exp(logit)).
    prob_seizure50 <- expit(logit_seizure50)

    prob_seizure50 ~ add(addSd_prob_seizure50)
  })
}
