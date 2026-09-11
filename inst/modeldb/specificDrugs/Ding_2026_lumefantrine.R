# Joint parent-metabolite population PK model for oral lumefantrine and its
# active metabolite desbutyl-lumefantrine in patients with uncomplicated
# Plasmodium falciparum malaria, pooled from the TRACII (NCT02453308) and
# TACT-CV (NCT03355664) trials (Ding 2026, Br J Clin Pharmacol
# 92(2):589-605; doi:10.1002/bcp.70301).

Ding_2026_lumefantrine <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model for oral lumefantrine and",
    "its active metabolite desbutyl-lumefantrine in adults and children",
    "with acute uncomplicated Plasmodium falciparum malaria, given the",
    "standard six-dose artemether-lumefantrine regimen alone or together",
    "with amodiaquine as a triple artemisinin-based combination therapy",
    "(Ding 2026, pooled TRACII + TACT-CV, n = 885).",
    "Five-transit-compartment absorption feeds a two-compartment",
    "lumefantrine disposition model, with complete molar-corrected",
    "bioconversion to a two-compartment desbutyl-lumefantrine disposition",
    "model. Allometric body-weight scaling on all apparent clearances",
    "(fixed exponent 0.75) and apparent volumes (fixed exponent 1.0) at a",
    "reference weight of 45 kg. Relative bioavailability is anchored at 1",
    "and falls with admission parasitaemia, with the milligram-per-kilogram",
    "dose and with admission body temperature; the apparent central volume",
    "is lower in the TACT-CV trial; and desbutyl-lumefantrine clearance",
    "matures with age. Coadministered amodiaquine was not a significant",
    "covariate on any parameter. Predictions are plasma lumefantrine and",
    "desbutyl-lumefantrine concentrations in ng/mL.",
    sep = " "
  )
  reference <- paste(
    "Ding J, Hoglund RM, van der Pluijm RW, Callery JJ, Peto TJ, Tripura R,",
    "Das S, Nguyen HC, Promnarate C, Mukaka M, Dysoley L, Fanello C,",
    "Onyamboko MA, Anvikar AR, Mayxay M, Smithuis F, von Seidlein L,",
    "Dhorda M, Amaratunga C, Faiz MA, Ho DTN, White NJ, Day NPJ,",
    "Dondorp AM, Tarning J (2026). Population pharmacokinetics of",
    "artemether-lumefantrine plus amodiaquine in patients with",
    "uncomplicated Plasmodium falciparum malaria. British Journal of",
    "Clinical Pharmacology 92(2):589-605. doi:10.1002/bcp.70301.",
    sep = " "
  )
  vignette <- "Ding_2026_artemether_lumefantrine_amodiaquine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Ding 2026 Methods ('PK sampling
  # scheme', 'Drug quantification': venous plasma lumefantrine and
  # desbutyllumefantrine by LC-MS/MS) and Figure S8.
  compartmentData <- list(
    depot           = list(analyte = "lumefantrine",           units = "mg", specimen = "administration site", verified = TRUE),
    transit1        = list(analyte = "lumefantrine",           units = "mg", specimen = "administration site", verified = TRUE),
    transit2        = list(analyte = "lumefantrine",           units = "mg", specimen = "administration site", verified = TRUE),
    transit3        = list(analyte = "lumefantrine",           units = "mg", specimen = "administration site", verified = TRUE),
    transit4        = list(analyte = "lumefantrine",           units = "mg", specimen = "administration site", verified = TRUE),
    transit5        = list(analyte = "lumefantrine",           units = "mg", specimen = "administration site", verified = TRUE),
    central         = list(analyte = "lumefantrine",           units = "mg", specimen = "plasma",              verified = TRUE),
    peripheral1     = list(analyte = "lumefantrine",           units = "mg", specimen = "plasma",              verified = TRUE),
    central_desbutlum     = list(analyte = "desbutyl-lumefantrine",  units = "mg", specimen = "plasma",              verified = TRUE),
    peripheral1_desbutlum = list(analyte = "desbutyl-lumefantrine",  units = "mg", specimen = "plasma",              verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at baseline. Ding 2026 Methods ('Covariates model'):",
        "'bodyweight was included on all clearance and volume parameters",
        "using a conventional allometric function with fixed exponents of",
        "0.75 and 1.0, respectively'. The Table 4 footnote fixes the",
        "reference: 'Population estimates are given for a typical adult",
        "patient weighing 45 kg', so 45 kg is the normalising constant",
        "encoded here. Strongly supported for both analytes",
        "(Results 3.1.3: delta-OFV = -87.52 for lumefantrine and -730.44",
        "for desbutyl-lumefantrine). WT also enters the bioavailability",
        "dose term as DOSE / WT. Cohort median 41.5 kg (TRACII) and 52.2 kg",
        "(TACT-CV); range 9.0-98.8 kg (Table 1).",
        sep = " "
      ),
      source_name        = "BW"
    ),
    DOSE = list(
      description        = "Lumefantrine amount administered at the current dose event (mg)",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Use case (a) of the canonical DOSE entry, supplied per dose",
        "record. Required so the model can form the milligram-per-kilogram",
        "dose DOSE / WT that drives relative bioavailability. The Table 4",
        "footnote gives the form verbatim: 'dose (mg/kg) was added on F as",
        "[1 + (theta x (Dose-9.8)]', and the same footnote names the",
        "reference patient as 'receiving 9.8 mg/kg dose of lumefantrine',",
        "so 9.8 mg/kg is a PER-DOSE value, not a daily one (Table 1 reports",
        "18.3-20.9 mg/kg/day for the twice-daily regimen). Standard",
        "fixed-dose artemether-lumefantrine tablets contain 120 mg",
        "lumefantrine, and the per-dose tablet count is 1 (5-14.9 kg),",
        "2 (15-24.9 kg), 3 (25-34.9 kg) or 4 (>35 kg) (Table S1), so DOSE",
        "is 120, 240, 360 or 480 mg. Set DOSE at each dose event in the",
        "rxode2 event table alongside amt (mg); the two are numerically",
        "equal here because lumefantrine is dosed into the depot directly.",
        "The relationship is EMPIRICAL and bounded: Results 3.1.3 states",
        "'this linear dose-covariate relationship on relative",
        "bioavailability is not suitable for extrapolation beyond the dose",
        "range studied here' (per-dose 4.0-16.0 mg/kg over the weight",
        "bands).",
        sep = " "
      ),
      source_name        = "Dose"
    ),
    PARA = list(
      description        = "Admission asexual Plasmodium falciparum parasite count",
      units              = "parasites/uL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Admission-only / time-fixed. The Table 4 footnote gives the form",
        "verbatim: 'Baseline parasitaemia (PBM) was implemented on F as",
        "[1 + (theta x (log(PBM)-4.56)]'. The logarithm is base 10, pinned",
        "by the reference value itself: the centring constant 4.56",
        "corresponds to 10^4.56 = 36,300 parasites/uL, which sits squarely",
        "inside the pooled cohort's per-trial medians of 47,500-52,500",
        "(TRACII) and 14,390-21,500 (TACT-CV) parasites/uL (Table 1),",
        "whereas a natural-log reading would imply exp(4.56) = 96",
        "parasites/uL, far below every reported median. The log10 transform",
        "is applied inside model() with the max(PARA, 1) gating convention",
        "shared by the sibling Mahidol-Oxford malaria models",
        "(Kloprogge_2014_quinine.R, Kloprogge_2018_lumefantrine.R,",
        "Tarning_2012_dihydroartemisinin.R, Ding_2024_piperaquine.R); the",
        "Table 1 cohort minimum is 1 parasite/uL. The negative coefficient",
        "matches the sibling Ding_2024_piperaquine.R and Hoglund 2017",
        "(higher admission parasitaemia, lower relative bioavailability).",
        sep = " "
      ),
      source_name        = "PBM"
    ),
    BODYTEMP = list(
      description        = "Admission body temperature",
      units              = "degC",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Admission-only / time-fixed. The Table 4 footnote gives the form",
        "verbatim: 'baseline temperature (TEMP) was included on F using",
        "[1 + (theta x (TEMP-37.5)]', and the same footnote names the",
        "reference patient as having a 'temperature of 37.5 C', so 37.5",
        "degC is the centring constant encoded here. Cohort medians 37.5",
        "and 37.7 degC, range 35.0-40.9 degC (Table 1). Discussion:",
        "'Baseline temperature was a significant covariate on",
        "bioavailability, with reduced bioavailability in patients with",
        "high baseline body temperature. This could be a result of reduced",
        "absorption due to malaria illness.' This is a different reference",
        "value and a different affected parameter from the sibling",
        "Kloprogge_2013_lumefantrine.R, which centres at 36.9 degC and puts",
        "the effect on mean transit time.",
        sep = " "
      ),
      source_name        = "TEMP"
    ),
    STUDY_TACTCV = list(
      description        = "TACT-CV trial indicator (1 = TACT-CV, 0 = TRACII)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = TRACII (NCT02453308)",
      notes              = paste(
        "Binary study-of-origin indicator for the pooled two-trial",
        "analysis. The Table 4 footnote gives the form verbatim: 'Study",
        "effect (TACT) was implemented on VC/F LF as a proportional",
        "equation [1 + theta x TACT]', with theta = -28.1%. Results 3.1.3",
        "describes it as 'a study effect on the central volume of",
        "distribution', so the indicator distinguishes the two TRIALS, not",
        "the two treatment ARMS -- coadministered amodiaquine was tested",
        "separately and found not to affect lumefantrine PK. TACT-CV is",
        "coded 1 because it is the trial whose acronym the footnote",
        "abbreviates, leaving TRACII as the reference (0). The effect is",
        "small in consequence: Discussion, 'a study effect on central",
        "volume had minimal impact on terminal elimination half-life",
        "(182 vs. 191 h), Day 7 concentration (444 vs. 454 ng/mL) and no",
        "change in total drug exposure. As a result, this study effect was",
        "not considered to be clinically meaningful.' The paper attributes",
        "it to the different fatty food given with the dose (a fatty snack",
        "in TRACII, 80 mL milk in TACT-CV).",
        sep = " "
      ),
      source_name        = "TACT"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at baseline. Enters desbutyl-lumefantrine clearance",
        "only, as the maturation function the Table 4 footnote gives",
        "verbatim: 'Age was added on CL/F DLF as [age/(theta + age)]', a",
        "plain hyperbolic maturation whose theta is the age at half of full",
        "maturation. Discussion: 'age substantially influenced the",
        "clearance of desbutyllumefantrine in an age-dependent maturation",
        "manner ... This finding could be interpreted by the age-relevant",
        "maturation of UGT enzyme (the main metabolism enzyme involved in",
        "the metabolism of desbutyl-lumefantrine) in children.' Because the",
        "function is bounded above by 1, the tabulated CL/F DLF of 744 L/h",
        "is the fully-matured adult value and the typical patient's",
        "clearance is always below it. Cohort median 17.0 years (TRACII)",
        "and 25.0 years (TACT-CV); range 1.6-65.0 years (Table 1). No",
        "age-maturation effect was evaluated for amodiaquine or",
        "desethylamodiaquine (Discussion, limitation 4).",
        sep = " "
      ),
      source_name        = "age"
    )
  )

  covariatesDataExcluded <- list(
    CONMED_AMODIAQUINE = list(
      description = "Coadministration of amodiaquine (triple ACT versus artemether-lumefantrine alone)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened as a binary drug-drug-interaction covariate on every PK",
        "parameter of lumefantrine and desbutyl-lumefantrine, and",
        "additionally assessed by a 500-bootstrap full covariate model, but",
        "NOT retained: Ding 2026 Results 3.1.3, 'Coadministration of",
        "amodiaquine did not affect the PK properties of lumefantrine. This",
        "was further confirmed by the full covariate model, which showed no",
        "significant DDI effect of amodiaquine on the main primary PK",
        "parameters of lumefantrine (Figure S9)', and the same for",
        "desbutyl-lumefantrine. The final model therefore has no",
        "amodiaquine term, and this single model serves both the",
        "artemether-lumefantrine and the artemether-lumefantrine-amodiaquine",
        "arms.",
        sep = " "
      )
    ),
    PARA_GAMETOCYTE = list(
      description = "Admission Plasmodium falciparum gametocyte density",
      units       = "parasites/uL",
      type        = "continuous",
      notes       = paste(
        "Identified as statistically significant on desbutyl-lumefantrine",
        "inter-compartmental clearance by the stepwise covariate search but",
        "deliberately discarded: Ding 2026 Results 3.1.3, 'One additional",
        "covariate was identified in the stepwise covariate process;",
        "gametocyte density on intercompartment clearance. However,",
        "considering its biologically implausible nature, this covariate",
        "was not retained in the final model.' Methods ('Covariates model')",
        "confirms the general policy: 'In the final covariate model,",
        "covariates that were statistically significant but biologically",
        "implausible were not retained in the final model.' No point",
        "estimate is published for the discarded coefficient, which is why",
        "no canonical entry is proposed for it in",
        "inst/references/covariate-columns.md; the name is written to sit",
        "beside the registered asexual-parasitaemia canonical PARA.",
        "Gametocytes were present in only 5.2-10.3% of patients (Table 1).",
        sep = " "
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 885L,
    n_studies      = 2L,
    n_observations = paste(
      "1448 lumefantrine (96 below the lower limit of quantification, 93 of",
      "them in the absorption phase) and 1448 desbutyl-lumefantrine (319",
      "below, 283 of them in the absorption phase) plasma concentrations;",
      "censored values were retained by the Beal M3 likelihood method",
      "(Results 3.1.3)"
    ),
    age_range      = "17.0 years (TRACII) and 25.0 years (TACT-CV), medians (range 1.6-65.0 years) (Table 1)",
    weight_range   = "41.5 kg (TRACII) and 52.2 kg (TACT-CV), medians (range 9.0-98.8 kg) (Table 1)",
    sex_female_pct = 24.6,
    disease_state  = paste(
      "Acute uncomplicated Plasmodium falciparum malaria. Median admission",
      "asexual parasite count 47,500-52,500 parasites/uL (TRACII) and",
      "14,390-21,500 parasites/uL (TACT-CV); median admission body",
      "temperature 37.5-37.7 degC, range 35.0-40.9 degC (Table 1)."
    ),
    dose_range     = paste(
      "Standard fixed-dose artemether-lumefantrine (20 mg artemether +",
      "120 mg lumefantrine per tablet) given orally as six doses over 3 days",
      "at 0, 8, 24, 36, 48 and 60 h, directly observed, with a fatty snack",
      "(TRACII) or 80 mL milk (TACT-CV) to improve lumefantrine absorption.",
      "Tablets per dose by weight band: 1 (5-14.9 kg), 2 (15-24.9 kg),",
      "3 (25-34.9 kg), 4 (>35 kg) (Table S1). Cohort median lumefantrine",
      "dose 18.3-20.9 mg/kg/day (range 9.5-32.0) (Table 1)."
    ),
    regions        = paste(
      "TRACII (NCT02453308, n = 575): Bangladesh, India, Myanmar,",
      "Democratic Republic of Congo and Lao PDR. TACT-CV (NCT03355664,",
      "n = 310): western and eastern Cambodia and Vietnam."
    ),
    notes          = paste(
      "All 885 randomised patients contribute, combining dense PK sampling",
      "in 79 patients (1, 2, 4, 6, 8, 12, 24, 64 h and Days 4, 7, 14, 28",
      "after the first dose; additionally 52 h in TRACII) with sparse",
      "sampling at baseline, Day 7 and at any recurrent infection detected",
      "during 42-day follow-up in the remainder. Lumefantrine and",
      "desbutyl-lumefantrine were fitted SEQUENTIALLY rather than",
      "simultaneously (Methods, 'Population PK analysis'): 'lumefantrine is",
      "the major contributor to the malaria parasite killing effect",
      "compared to its active metabolite, and to avoid the risk of biasing",
      "the modelling of the parent drug ... we therefore used a sequential",
      "modelling approach where lumefantrine data were modelled first, then",
      "desbutyllumefantrine was modelled with fixed individual parameters",
      "for lumefantrine.' The assay lower limits of quantification were",
      "9.71 ng/mL for lumefantrine and 1.01 ng/mL for",
      "desbutyl-lumefantrine."
    )
  )

  ini({
    # ---- Absorption --------------------------------------------------
    # Ding 2026 Results 3.1.3: "Lumefantrine plasma concentration-time data
    # were best described by a model with two disposition compartments and
    # five transit compartments characterizing the absorption process."
    # Table 4 reports a mean transit time and a fixed transit-compartment
    # count but NO separate absorption rate constant, so a single rate
    # governs all six transfers (depot -> transit1 -> ... -> transit5 ->
    # central) and the mean transit time spans all six:
    # ktr = (NN + 1) / MTT = 6 / MTT. Same reading as the sibling
    # Ding_2026_artemether.R.
    lmtt <- log(5.43)
    label("Mean transit time of the lumefantrine absorption chain (h)")
    # Ding 2026 Table 4: Mean transit time = 5.43 h (%RSE 5.8; SIR median
    # 5.56, 95% CI 5.00-6.25; eta shrinkage 69.9%)

    # ---- Lumefantrine disposition (two compartments) -----------------
    # Ding 2026 Table 4, "NONMEM estimates" column. Values are apparent
    # (relative to F = 1) and reported on the linear scale for a typical
    # adult patient weighing 45 kg with a baseline parasite density of
    # 10^4.56 parasites/uL and a temperature of 37.5 degC, receiving
    # 9.8 mg/kg of lumefantrine; log() is applied here for the nlmixr2
    # internal log scale.
    lcl <- log(4.35)
    label("Apparent lumefantrine elimination clearance CL/F at WT = 45 kg (L/h)")
    # Ding 2026 Table 4: CL/F LF = 4.35 L/h (%RSE 3.4; SIR median 4.35,
    # 95% CI 4.04-4.61). Discussion cross-checks this against a previously
    # reported pooled value of 5.3 L/h.

    lvc <- log(101)
    label("Apparent lumefantrine central volume of distribution Vc/F at WT = 45 kg in the TRACII trial (L)")
    # Ding 2026 Table 4: Vc/F LF = 101 L (%RSE 8.5; SIR median 101,
    # 95% CI 88.0-121; eta shrinkage 64.1%)

    lq <- log(1.65)
    label("Apparent lumefantrine inter-compartmental clearance Q/F at WT = 45 kg (L/h)")
    # Ding 2026 Table 4: Q/F LF = 1.65 L/h (%RSE 6.0; SIR median 1.66,
    # 95% CI 1.45-1.83)

    lvp <- log(311)
    label("Apparent lumefantrine peripheral volume of distribution Vp/F at WT = 45 kg (L)")
    # Ding 2026 Table 4: Vp/F LF = 311 L (%RSE 5.3; SIR median 311,
    # 95% CI 275-343)

    # ---- Desbutyl-lumefantrine disposition (two compartments) --------
    lcl_desbutlum <- log(744)
    label("Apparent fully-matured desbutyl-lumefantrine elimination clearance CL/F at WT = 45 kg (L/h)")
    # Ding 2026 Table 4: CL/F DLF = 744 L/h (%RSE 3.5; SIR median 721,
    # 95% CI 677-775; eta shrinkage 71.5%). This is the asymptotic adult
    # value: the age maturation term age/(theta + age) is bounded above by
    # 1, so a typical patient's clearance is always below 744 L/h.
    # Discussion notes the estimate is higher than a previously reported
    # 298 L/h while the apparent half-life is comparable (139 vs 148 h).

    lvc_desbutlum <- log(8530)
    label("Apparent desbutyl-lumefantrine central volume of distribution Vc/F at WT = 45 kg (L)")
    # Ding 2026 Table 4: Vc/F DLF = 8530 L (%RSE 10.0; SIR median 8580,
    # 95% CI 7170-10,500; eta shrinkage 73.6%)

    lq_desbutlum <- log(1100)
    label("Apparent desbutyl-lumefantrine inter-compartmental clearance Q/F at WT = 45 kg (L/h)")
    # Ding 2026 Table 4: Q/F DLF = 1100 L/h (%RSE 7.5; SIR median 1050,
    # 95% CI 879-1200)

    lvp_desbutlum <- log(62600)
    label("Apparent desbutyl-lumefantrine peripheral volume of distribution Vp/F at WT = 45 kg (L)")
    # Ding 2026 Table 4: Vp/F DLF = 62,600 L (%RSE 3.3; SIR median 60,900,
    # 95% CI 56,500-64,900)

    # ---- Relative bioavailability ------------------------------------
    lfdepot <- fixed(log(1))
    label("Relative oral bioavailability of lumefantrine at the reference covariate values (unitless)")
    # Ding 2026 Table 4: F = 1 Fix. Methods ('Population PK analysis'):
    # "Relative bioavailability (F) was fixed to unity in the population,
    # allowing for quantification of the IIV in the absorption process."

    # ---- Allometric exponents ----------------------------------------
    e_wt_cl <- fixed(0.75)
    label("Allometric WT exponent on all apparent clearance parameters (CL/F and Q/F of both analytes)")
    # Ding 2026 Methods, 'Covariates model': clearance exponent fixed 0.75

    e_wt_vc <- fixed(1.00)
    label("Allometric WT exponent on all apparent volume parameters (Vc/F and Vp/F of both analytes)")
    # Ding 2026 Methods, 'Covariates model': volume exponent fixed 1.0

    # ---- Covariate effects on relative bioavailability ---------------
    # All three are linear-deviation terms multiplying F, given verbatim in
    # the Table 4 footnote. All three coefficients are NEGATIVE: the
    # printed NONMEM point estimates carry a leading unicode minus that
    # several text extractors drop, but the accompanying SIR confidence
    # intervals are unambiguous (e.g. -11.6 with 95% CI -18.0 to -5.6), and
    # the Discussion states the temperature effect's direction in words
    # ("reduced bioavailability in patients with high baseline body
    # temperature").
    e_para_f <- -0.116
    label("Linear fractional change in relative bioavailability per log10 unit of admission parasitaemia above 10^4.56 parasites/uL (per log10 unit)")
    # Ding 2026 Table 4: Baseline parasite density on F (%) = -11.6
    # (%RSE 26.4; SIR median -11.6, 95% CI -18.0 to -5.6)

    e_dose_f <- -0.0647
    label("Linear fractional change in relative bioavailability per mg/kg of lumefantrine dose above 9.8 mg/kg (per mg/kg)")
    # Ding 2026 Table 4: Dose (mg/kg) on F = -6.47 (%RSE 14.3; SIR median
    # -6.54, 95% CI -8.50 to -4.82). Read as a PERCENTAGE per mg/kg, like
    # its three table-mates: the Dose row is the only one of the four
    # covariate rows whose label omits the '(%)' suffix, but a raw-fraction
    # reading is arithmetically impossible -- it would give
    # F = 1 - 6.47 * (10.67 - 9.8) = -4.6 for a 45 kg adult on the standard
    # four-tablet dose. At -6.47% per mg/kg the same patient gets F = 0.944
    # and the model reproduces the published typical lumefantrine AUC
    # (625 vs 600 h*ug/mL) and Day 7 concentration (465 vs 452 ng/mL).
    # Dose-dependent (saturable) lumefantrine absorption is independently
    # established: the sibling Kloprogge_2018_lumefantrine.R, which is
    # reference 34 of this paper, fits it as a saturable Emax term.

    e_bodytemp_f <- -0.122
    label("Linear fractional change in relative bioavailability per degC of admission body temperature above 37.5 degC (per degC)")
    # Ding 2026 Table 4: Baseline temperature on F (%) = -12.2 (%RSE 16.4;
    # SIR median -12.4, 95% CI -16.2 to -8.43)

    # ---- Covariate effect on the lumefantrine central volume ---------
    e_study_tactcv_vc <- -0.281
    label("Proportional change in apparent lumefantrine central volume in the TACT-CV trial relative to TRACII (unitless)")
    # Ding 2026 Table 4: Study effect (TACT) on VC/F LF (%) = -28.1
    # (%RSE 14.7; SIR median -27.7, 95% CI -36.3 to -20.0)

    # ---- Covariate effect on desbutyl-lumefantrine clearance ---------
    e_age_cl_desbutlum <- 10.1
    label("Age at half of full maturation of apparent desbutyl-lumefantrine clearance (years)")
    # Ding 2026 Table 4: Age on CL/F DLF (year) = 10.1 (%RSE 9.2; SIR
    # median 9.8, 95% CI 8.1-11.7). The Results and Discussion narrative
    # instead quotes an Age50 of 10.6 years in two places; the Table 4
    # value is used here because it is the final parameter estimate and is
    # the value the accompanying SIR interval brackets. The difference
    # changes a typical 20-year-old's desbutyl-lumefantrine clearance by
    # 1.6%.

    # ---- Inter-individual and inter-site variability -----------------
    # Ding 2026 Table 4 footnote: "Coefficients of variation for
    # interindividual variability and intersite variability (IIV and ISV)
    # were calculated as 100 x (e^variance - 1)^(1/2)", so the internal
    # log-scale variance is recovered as omega^2 = log(CV^2 + 1).
    #
    #   MTT       IIV 62.7% -> omega^2 = log(0.627^2 + 1) = 0.3315523
    #   F         IIV 57.0% -> omega^2 = log(0.570^2 + 1) = 0.2813370
    #   F         ISV 13.3% -> omega^2 = log(0.133^2 + 1) = 0.0175344
    #   Vc/F LF   IIV 89.3% -> omega^2 = log(0.893^2 + 1) = 0.5863684
    #   CL/F DLF  IIV 14.9% -> omega^2 = log(0.149^2 + 1) = 0.0219581
    #   Vc/F DLF  IIV 101%  -> omega^2 = log(1.010^2 + 1) = 0.7031470
    #
    # Table 4 reports no random effect on CL/F LF, Q/F LF, Vp/F LF, Q/F DLF
    # or Vp/F DLF, so no eta slots are created for those parameters. Unlike
    # the artemether and amodiaquine models, Results 3.1.3 reports that
    # "Inclusion of interoccasion variability on the relative
    # bioavailability and MTT did not improve model fit, and interoccasion
    # variability was therefore not retained in the final model", so there
    # is no IOV here.
    etalmtt ~ 0.3315523
    # Ding 2026 Table 4: IIV on mean transit time = 62.7% CV (%RSE 15.7;
    # SIR median 58.5, 95% CI 49.5-72.0; eta shrinkage 69.9%)

    etalfdepot ~ 0.2988714
    # Ding 2026 Table 4: IIV on F = 57.0% CV (%RSE 12.3; SIR median 56.9,
    # 95% CI 49.5-64.9; eta shrinkage 44.6%) PLUS ISV on F = 13.3% CV
    # (%RSE 19.1; SIR median 13.3, 95% CI 10.3-15.5; eta shrinkage 24.3%),
    # added on the final model (Results 3.1.3: delta-OFV = -34.69).
    # nlmixr2 supports a single (subject) level of random effects, so the
    # published two-level structure is collapsed to its exact SUBJECT-LEVEL
    # MARGINAL: each patient belongs to exactly one site, and the two
    # log-scale random effects are independent, so the marginal variance of
    # log(F) across patients is the sum,
    # 0.2813370 + 0.0175344 = 0.2988714 (equivalently 58.5% CV rather than
    # 57.0%). This is the encoding that reproduces the published
    # across-cohort exposure percentiles; what it cannot reproduce is the
    # correlation between patients treated at the same site. See the
    # vignette's Assumptions and deviations section.

    etalvc ~ 0.5863684
    # Ding 2026 Table 4: IIV on Vc/F LF = 89.3% CV (%RSE 14.3; SIR median
    # 89.3, 95% CI 71.2-103; eta shrinkage 64.1%)

    etalcl_desbutlum ~ 0.0219581
    # Ding 2026 Table 4: IIV on CL/F DLF = 14.9% CV (%RSE 12.5; SIR median
    # 15.0, 95% CI 11.5-18.7; eta shrinkage 71.5%)

    etalvc_desbutlum ~ 0.7031470
    # Ding 2026 Table 4: IIV on Vc/F DLF = 101% CV (%RSE 7.8; SIR median
    # 97.5, 95% CI 77.1-122; eta shrinkage 73.6%)

    # ---- Residual unexplained variability ----------------------------
    # As for the sibling models, the paper's additive-on-log-scale residual
    # maps to a proportional residual in linear concentration space, and
    # the Table 4 footnote states "RUV is the residual error variance", so
    # the tabulated number is a variance and the SD is its square root.
    propSd <- sqrt(0.297)
    label("Proportional residual SD for lumefantrine plasma concentration (SD on log scale)")
    # Ding 2026 Table 4: RUV = 0.297 (variance; %RSE 6.5; SIR median 0.304,
    # 95% CI 0.270-0.344; epsilon shrinkage 24.3%)

    propSd_desbutlum <- sqrt(0.178)
    label("Proportional residual SD for desbutyl-lumefantrine plasma concentration (SD on log scale)")
    # Ding 2026 Table 4: RUV = 0.178 (variance; %RSE 5.5; SIR median 0.177,
    # 95% CI 0.160-0.199; epsilon shrinkage 0.1%)
  })

  model({
    # Molecular weights (g/mol). Ding 2026 Methods ('Population PK
    # analysis'): "Parent drugs were assumed to be completely metabolized
    # to their metabolites due to identifiability issues with other model
    # structures." The conversion factor is not printed, so the mass flux
    # leaving lumefantrine central is molar-corrected before it enters
    # desbutyl-lumefantrine central, matching the sibling
    # Ding_2024_amodiaquine.R from the same group. Desbutyl-lumefantrine is
    # lumefantrine less a butyl group (C4H8, 56.11 g/mol). The paper's own
    # secondary parameters support the correction: for a typical 20-year-old
    # 45 kg patient on the standard 2880 mg total dose, the published
    # desbutyl-lumefantrine AUC of 4.59 h*ug/mL is reproduced as 4.78 with
    # the molar correction against 5.35 without it.
    mwLF        <- 528.94
    mwDLF       <- 472.83
    molarFactor <- mwDLF / mwLF

    # Absorption chain rate constant. NN = 5 transit compartments and the
    # single rate governs all NN + 1 = 6 transfers, so ktr = 6 / MTT (see
    # the lmtt annotation in ini()).
    mtt <- exp(lmtt + etalmtt)
    ktr <- 6 / mtt

    # Covariate terms on relative bioavailability, in the linear-deviation
    # forms of the Ding 2026 Table 4 footnote. The log10 transform of
    # parasitaemia is applied here (rather than at dataset assembly) with
    # the max(PARA, 1) gating convention of the sibling Mahidol-Oxford
    # malaria models. The dose term uses the milligram-per-kilogram dose of
    # the current dose event.
    doseMgKg <- DOSE / WT
    f_para   <- 1 + e_para_f     * (log10(max(PARA, 1)) - 4.56)
    f_dose   <- 1 + e_dose_f     * (doseMgKg - 9.8)
    f_temp   <- 1 + e_bodytemp_f * (BODYTEMP - 37.5)

    # Individual PK parameters. Allometric weight scaling on all apparent
    # clearances (exponent 0.75) and apparent volumes (exponent 1) centred
    # on the 45 kg reference of the Table 4 footnote. The lumefantrine
    # central volume additionally carries the proportional study effect,
    # and desbutyl-lumefantrine clearance carries the age maturation.
    cl <- exp(lcl)         * (WT / 45)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 45)^e_wt_vc *
          (1 + e_study_tactcv_vc * STUDY_TACTCV)
    q  <- exp(lq)  * (WT / 45)^e_wt_cl
    vp <- exp(lvp) * (WT / 45)^e_wt_vc

    cl_desbutlum <- exp(lcl_desbutlum + etalcl_desbutlum) * (WT / 45)^e_wt_cl *
              (AGE / (e_age_cl_desbutlum + AGE))
    vc_desbutlum <- exp(lvc_desbutlum + etalvc_desbutlum) * (WT / 45)^e_wt_vc
    q_desbutlum  <- exp(lq_desbutlum)  * (WT / 45)^e_wt_cl
    vp_desbutlum <- exp(lvp_desbutlum) * (WT / 45)^e_wt_vc

    # Micro-rate constants (1/h).
    kel     <- cl / vc
    k12     <- q  / vc
    k21     <- q  / vp
    kel_desbutlum <- cl_desbutlum / vc_desbutlum
    k12_desbutlum <- q_desbutlum  / vc_desbutlum
    k21_desbutlum <- q_desbutlum  / vp_desbutlum

    # ODE system (Ding 2026 Figure S8). Compartment amounts are in mg of
    # analyte and volumes are in L, so amount/volume is mg/L and is scaled
    # to ng/mL below.
    #
    # Absorption: depot -> transit1 -> ... -> transit5 -> lumefantrine
    # central, all six transfers at the same rate ktr.
    d/dt(depot)    <- -ktr * depot
    d/dt(transit1) <-  ktr * depot    - ktr * transit1
    d/dt(transit2) <-  ktr * transit1 - ktr * transit2
    d/dt(transit3) <-  ktr * transit2 - ktr * transit3
    d/dt(transit4) <-  ktr * transit3 - ktr * transit4
    d/dt(transit5) <-  ktr * transit4 - ktr * transit5

    # Lumefantrine central plus one peripheral compartment. The entire mass
    # flux leaving lumefantrine central by elimination (kel * central) is
    # routed to desbutyl-lumefantrine central under the complete-conversion
    # assumption, after the molar correction.
    d/dt(central)     <- ktr * transit5 - kel * central -
                         k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Desbutyl-lumefantrine central plus one peripheral compartment.
    d/dt(central_desbutlum)     <- molarFactor * kel * central -
                             kel_desbutlum * central_desbutlum -
                             k12_desbutlum * central_desbutlum + k21_desbutlum * peripheral1_desbutlum
    d/dt(peripheral1_desbutlum) <- k12_desbutlum * central_desbutlum - k21_desbutlum * peripheral1_desbutlum

    # Relative oral bioavailability on the lumefantrine dose: the fixed
    # unity anchor, the subject-level random effect, and the three
    # multiplicative covariate terms.
    f(depot) <- exp(lfdepot + etalfdepot) * f_para * f_dose * f_temp

    # Plasma concentrations in ng/mL: amount (mg) / volume (L) is mg/L,
    # multiplied by 1000 to give ng/mL, the units used for the Day 7
    # concentration in Ding 2026 Table 4 and throughout Figure 4.
    Cc     <- 1000 * central     / vc
    Cc_desbutlum <- 1000 * central_desbutlum / vc_desbutlum

    # Proportional residual error on the linear-concentration scale, the
    # linear-space equivalent of the paper's additive-on-log-scale error.
    Cc     ~ prop(propSd)
    Cc_desbutlum ~ prop(propSd_desbutlum)
  })
}
