Du_2025_repotrectinib <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order absorption, an",
    "absorption lag time, and empirical trough-concentration-driven",
    "autoinduction of clearance for repotrectinib, a next-generation ROS1 /",
    "TRK / ALK tyrosine kinase inhibitor, in 644 pooled subjects (118 healthy",
    "volunteers and 526 patients with advanced solid tumors harboring ALK,",
    "ROS1 or NTRK1-3 rearrangements) across eight studies including 24",
    "pediatric patients. Repotrectinib is cleared by CYP3A4 and induces its",
    "own metabolism. Rather than the semi-mechanistic enzyme-turnover ODE or",
    "a discrete dose-driven step function, the authors drive the induction",
    "with the model-predicted trough concentration captured at each dosing",
    "event, multiplied by a hyperbolic function of time since the first dose.",
    "That formulation avoids the abrupt clearance jumps a dose-driven model",
    "produces when the regimen steps from 160 mg once daily to twice daily on",
    "day 15. Maximum induction raises clearance to exp(1.59) = 4.9 times",
    "baseline. Body weight scales clearances and volumes with estimated (not",
    "fixed) allometric exponents of 0.477 and 0.962 about a 70 kg reference,",
    "age below 18 years increases the maximum induction, central volume is",
    "much smaller in healthy volunteers than in patients, and prandial state",
    "at each dose selects among four absorption-rate and four bioavailability",
    "typical values. Between-subject variability on clearance and the",
    "combined proportional plus additive residual error are both stratified",
    "between healthy volunteers and patients."
  )
  reference <- paste(
    "Du S, Hu Z, Shen J, Hamuro L, Lam J, Lu M, Zhu L, Roy A, Kondic A. A",
    "Novel Empirical Autoinduction Model to Characterize the Population",
    "Pharmacokinetics and Recommend Dose for Repotrectinib in Adult and",
    "Adolescents With Advanced Solid Tumors Harboring ALK, ROS1, or NTRK1-3",
    "Rearrangements. CPT Pharmacometrics Syst Pharmacol. 2025;14(7):1179-1190.",
    "doi:10.1002/psp4.70036"
  )
  vignette <- "Du_2025_repotrectinib"
  units    <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # `ctroughHold` is a bookkeeping state, not a physiological compartment: it
  # is the rxode2 realization of the source control stream's sample-and-hold
  # of the pre-dose plasma concentration. See the extended comment in model().
  paper_specific_compartments <- c("ctroughHold")

  # Clearance IIV is stratified by cohort (two separate NONMEM ETAs, not one
  # eta with a covariate effect), and the bioavailability IIV is carried on
  # the logit scale, so none of the three matches the bare canonical names.
  paper_specific_etas <- c("etalclHv", "etalclPt", "etalogitfdepot")

  # Residual error is stratified by cohort; the canonical matcher recognises
  # only the bare names propSd / addSd / expSd (same pattern as
  # Ozdin_2025_dexamethasone.R and Shoji_2011_pregabalin.R).
  paper_specific_residual_sds <- c(
    "propSdHv", "addSdHv",
    "propSdPt", "addSdPt"
  )

  compartmentData <- list(
    depot       = list(analyte = "repotrectinib", units = "mg",    specimen = "administration site", verified = TRUE),
    central     = list(analyte = "repotrectinib", units = "mg",    specimen = "plasma",              verified = TRUE),
    peripheral1 = list(analyte = "repotrectinib", units = "mg",    specimen = "plasma",              verified = TRUE),
    # Holds a CONCENTRATION (ng/mL), not an amount: it latches the value of Cc
    # at each dosing event and holds it until the next one.
    ctroughHold = list(analyte = "repotrectinib", units = "ng/mL", specimen = "not applicable",      verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline (time-fixed) body weight, source column WTB. Allometric",
        "power scaling about a 70 kg reference, with BOTH exponents ESTIMATED",
        "rather than fixed at the canonical 0.75 / 1 values -- this is one of",
        "the four stated Stage II enhancements over the Stage I model, whose",
        "control stream (Data S2) carries `(0.75) FIX` and `(1) FIX`. Du 2025",
        "Table 2 reports CLQWT = 0.477 (RSE 14.3%) shared by CL and Q and",
        "VCVPWT = 0.962 (RSE 5.4%) shared by Vc and Vp; supplement Data S3",
        "writes them as ALLOM1 = (WTB/70)**THETA(18) and",
        "ALLOM2 = (WTB/70)**THETA(19). Cohort range 5.90-169 kg, median 70.6",
        "kg (Table 1). Du 2025 Results section 3.2 reports the resulting",
        "effect sizes at the 95th weight percentile (102.05 kg) relative to",
        "the 70.55 kg reference: CL 19% higher and Vc / Vp both 43% greater,",
        "which reproduce as (102.05/70.55)^0.477 = 1.193 and",
        "(102.05/70.55)^0.962 = 1.427."
      ),
      source_name        = "WTB"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column AGE2. Enters ONLY through the maximum induction",
        "exponent cl_time_max, and only below 18 years: Du 2025 Equation (4)",
        "and supplement Data S3 write",
        "`IF (AGE2.GE.18) CLMAX = TVCLMAX` /",
        "`IF (AGE2.LT.18) CLMAX = TVCLMAX*(AGE2/18)**THETA(25)`, so 18 years",
        "is both the reference age and a hard knee -- the effect is clamped",
        "for every adult rather than extrapolated. The model file encodes the",
        "clamp as `(min(AGE, 18) / 18)^e_age_cl_time_max`, which is",
        "algebraically identical to the two-branch IF and avoids an indicator",
        "product. Estimate -0.292 (RSE 9.2%, Table 2 CLMAXAGE), so younger",
        "subjects autoinduce MORE. Du 2025 Results section 3.2 and the",
        "Discussion give four checks that all reproduce from -0.292:",
        "CLMAX 21.7% lower at the 95th pediatric percentile (15.7 y) than at",
        "the 6.8 y reference ((15.7/6.8)^-0.292 = 0.783); 77% higher at the",
        "5th percentile (0.965 y) ((0.965/6.8)^-0.292 = 1.769); about 13%",
        "higher for a typical 12-year-old than for an adult; and about 40%",
        "higher for a typical 6-year-old. Note the Discussion prints the",
        "exponent as -0.291 while Table 2 reports -0.292; the Table 2 final",
        "estimate is used here. Cohort age range 0.800-93.0 years, median",
        "51.5 (Table 1), with 16 patients under 12 years and 8 adolescents",
        "aged 12 to under 18 years."
      ),
      source_name        = "AGE2"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-volunteer cohort indicator; 0 = patient with an advanced solid tumor",
      units              = "(binary)",
      type               = "binary",
      reference_category = paste(
        "0 (patient with an advanced or metastatic solid tumor harboring an",
        "ALK, ROS1 or NTRK1-3 rearrangement, enrolled in TRIDENT-1",
        "(NCT03093116) or CARE (NCT04094610))"
      ),
      notes              = paste(
        "Source column HV, same orientation (1 = healthy volunteer), so no",
        "value transformation is needed. Time-fixed per subject. THREE",
        "distinct roles in this model.",
        "(1) Multiplicative LINEAR-deviation effect on central volume:",
        "`vc * (1 + e_dis_healthy_vc * DIS_HEALTHY)` with",
        "e_dis_healthy_vc = -0.854, i.e. healthy volunteers have a central",
        "volume 14.6% of the patient value (19.8 -> 2.89 L). The form is",
        "linear, NOT the exponential `exp(theta * X)` of Du 2025 Equation (2):",
        "supplement Data S3 writes `IF(HV.EQ.1) V2HV = (1 + THETA(20))` and",
        "`TVV2 = TVV2*V2HV`, and Table 2 names the parameter VCPOP with",
        "estimate -0.854 (RSE 2.6%). The patient cohort is the reference",
        "(V2HV = 1), which is why lvc carries the patient typical value. The",
        "large split is credible because the healthy-volunteer studies are",
        "richly sampled single-dose phase 1 trials that resolve the early",
        "distribution phase, whereas the sparse patient sampling cannot.",
        "(2) Binary stratifier selecting which of two clearance IIV terms",
        "applies (etalclHv versus etalclPt). Du 2025 Results section 3.1",
        "gives the reason: healthy volunteers received at most 3 doses, so",
        "autoinduction of CL is far better characterised in patients, and a",
        "single shared eta gave unacceptable shrinkage. Supplement Data S3",
        "writes this as an explicit IF/ELSE over ETA(1) and ETA(7).",
        "(3) Binary stratifier selecting which of the two combined",
        "proportional-plus-additive residual error magnitudes applies",
        "(propSdHv / addSdHv versus propSdPt / addSdPt; Data S3 $ERROR",
        "block, Table 2 'Residual error' rows).",
        "Cohort split: 118 healthy volunteers (18.3%) across six phase 1",
        "studies, 526 patients (81.7%) across TRIDENT-1 and CARE (Methods",
        "section 2.1, Table 1)."
      ),
      source_name        = "HV"
    ),
    FED = list(
      description        = "Fed-state indicator for the dose record",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not fed; the record is fasted, modified fasted, or of unknown prandial state)",
      notes              = paste(
        "Per dose record, NOT per subject: Du 2025 Table 1 footnote a states",
        "that subjects who switched between fasted and fed states were",
        "counted twice, so the indicator varies within subject across days.",
        "One of three indicator columns that jointly select among the four",
        "prandial levels of the source column FED2 (0 = fasted, 1 = fed,",
        "2 = modified fasted, -99 = unknown). FED = 1 selects the fed level.",
        "Fasted is the model's prandial reference and is derived inside",
        "model() as the complement, following the FED_HIGHFAT / FED_LOWFAT",
        "convention of carrying only the non-reference indicators. The four",
        "levels each carry their OWN absorption rate constant and their own",
        "bioavailability -- they are not deviations from the fasted value --",
        "so all eight typical values appear in ini(). Du 2025 Results section",
        "3.2: relative to fasted, ka is 129% higher and F1 46% higher in the",
        "fed state (0.124 vs 0.0541 /h; 0.76 vs 0.52). Cohort: 96 records",
        "fed (13.4%) of 715 prandial records (Table 1)."
      ),
      source_name        = "FED2 (level 1)"
    ),
    FASTED_STRICT = list(
      description        = "Strict-versus-relaxed fasting indicator for the dose record",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (relaxed fast; Du 2025's 'modified fasted' level)",
      notes              = paste(
        "Per dose record. READ ONLY WHEN FED = 0 AND FED_MISSING = 0 -- it",
        "distinguishes the two food-free levels of the source column FED2 and",
        "is not consulted on fed or unknown-prandial-state records, where the",
        "model() selection terms zero it out. Du 2025 defines its 'modified",
        "fasted' level in the Equations (4)-(12) abbreviation list as 'no food",
        "or beverages 1 h before and 2 h after dosing'; its plain 'fasted'",
        "level is the conventional phase 1 overnight fast. Both levels are",
        "food-free and they differ in how tightly food was excluded around the",
        "dose, which is exactly the axis this canonical carries, so",
        "FASTED_STRICT = 1 is Du 2025's 'fasted' and FASTED_STRICT = 0 is its",
        "'modified fasted'. Per the register's standing instruction to encode",
        "the protocol and never the phrase, the assignment was made from the",
        "printed 1 h / 2 h window rather than from the words 'modified",
        "fasted', which a sibling paper (Mauro 2025 nilotinib) uses for the",
        "unrelated meal-then-dose protocol carried by MEAL_PREDOSE_2H. The",
        "orientation is corroborated by the estimates themselves: modified",
        "fasted F1 = 0.639 sits between fasted 0.52 and fed 0.76, the",
        "ordering expected when a relaxed fast admits a partial food effect.",
        "Cohort: 115 records fasted (16.1%) and 71 modified fasted (9.9%) of",
        "715 prandial records (Table 1)."
      ),
      source_name        = "FED2 (levels 0 and 2)"
    ),
    FED_MISSING = list(
      description        = "Prandial-state-not-recorded indicator for the dose record",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (prandial state was recorded, i.e. fasted, fed or modified fasted)",
      notes              = paste(
        "Per dose record. Selects the FED2 = -99 level of the source column.",
        "This is a genuine estimated NUISANCE level, not a dropped-data flag:",
        "Du 2025 Discussion states that 'the unknown food effect was also",
        "estimated as a nuisance parameter' precisely so that the 434 records",
        "(60.6% of the 715 prandial records, Table 1) with no recorded",
        "prandial state could be retained rather than excluded. It therefore",
        "carries its own absorption rate constant (0.193 /h, the FASTEST of",
        "the four) and its own bioavailability (0.533, essentially the fasted",
        "value), and it must be supplied to reproduce the paper's pooled fit.",
        "Simulations of a defined prandial state should set it to 0.",
        "Du 2025 Results section 3.4 notes that subjects of unknown prandial",
        "state showed exposures comparable to the fed state on day 1 but 25%",
        "lower steady-state Cmin than the modified fasted state."
      ),
      source_name        = "FED2 (level -99)"
    )
  )

  covariatesDataExcluded <- list(
    FORM_CAPSULE = list(
      description = "Capsule versus oral-suspension formulation indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened but NOT retained. Du 2025 Methods section 2.2 enhancement",
        "(4) states that a formulation effect (capsules versus oral",
        "suspension) was assessed on ka alone, on F1 alone, and on both",
        "simultaneously; supplement Table S3 lists all three as candidate",
        "models. None improved BIC enough to be carried forward, and the",
        "final control stream (Data S3) selects ka and F1 from the prandial",
        "column FED2 only -- FORMN appears in $INPUT solely to set F1 = 1 for",
        "the intravenous formulation. No point estimate is reported."
      )
    ),
    RENIMP = list(
      description = "Renal impairment category (normal / mild / moderate, by eGFR)",
      units       = "(categorical)",
      type        = "categorical",
      notes       = paste(
        "Screened but NOT retained. Du 2025 Discussion reports that renal",
        "impairment was identified as a significant covariate on peripheral",
        "volume during the Stage I univariate search, but was dropped from",
        "the final model because of instability and convergence problems and",
        "because few subjects had severe impairment. Repotrectinib shows",
        "minimal renal elimination and the EBE-based exposure distributions",
        "overlap (Figure S5). Cohort: 448 normal (69.6%), 157 mild (24.4%),",
        "33 moderate (5.1%), 6 missing (Table 1). No point estimate is",
        "reported."
      )
    ),
    HEPIMP = list(
      description = "Liver dysfunction category (normal / mild / moderate)",
      units       = "(categorical)",
      type        = "categorical",
      notes       = paste(
        "Screened but NOT retained. Du 2025 Discussion reports that hepatic",
        "impairment was identified for CL or Vc during the Stage I univariate",
        "search but was not retained in the final model for the same",
        "instability and sample-size reasons as renal impairment; EBE-based",
        "exposures for mild impairment overlap those of normal hepatic",
        "function (Figure S6). Cohort: 582 normal (90.4%), 59 mild (9.2%),",
        "1 moderate (0.2%), 2 missing (Table 1). No point estimate is",
        "reported."
      )
    ),
    PRVTRT = list(
      description = "Prior tyrosine-kinase-inhibitor treatment status (naive / pretreated)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened but NOT retained. Du 2025 Discussion and Figure S7 report",
        "that TKI-naive and TKI-pretreated subjects showed similar",
        "steady-state exposures (differences below 15%). Cohort: 186 naive",
        "(28.9%), 340 pretreated (52.8%) (Table 1). No point estimate is",
        "reported."
      )
    ),
    GENMT = list(
      description = "Driver-mutation type (ALK / ROS1 / NTRK1 / NTRK2 / NTRK3)",
      units       = "(categorical)",
      type        = "categorical",
      notes       = paste(
        "Screened but NOT retained. Du 2025 Results section 3.4 and Figure S4",
        "report steady-state exposure differences below 20% across mutation",
        "groups. Cohort: 35 ALK, 363 ROS1, 54 NTRK1, 9 NTRK2, 65 NTRK3",
        "(Table 1). No point estimate is reported."
      )
    ),
    RACE = list(
      description = "Race group (White / Black or African American / Asian / Other)",
      units       = "(categorical)",
      type        = "categorical",
      notes       = paste(
        "Screened but NOT retained as a structural covariate. Du 2025 Results",
        "section 3.4 and Figure S3 report that Asian patients had higher",
        "geometric-mean steady-state exposures than White or other races, but",
        "the difference was below 38% with substantially overlapping",
        "distributions and was judged not clinically relevant; the authors",
        "attribute much of it to the 18% lower median body weight of the",
        "Asian subgroup, which the retained WT covariate already captures.",
        "Cohort: 319 White (49.5%), 43 Black or African American (6.7%), 245",
        "Asian (38.0%), 8 Other (1.2%), 29 unknown (4.5%) (Table 1). No point",
        "estimate is reported."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 644L,
    n_studies      = 8L,
    n_observations = 9220L,
    age_range      = "0.800-93.0 years",
    age_median     = "51.5 years",
    weight_range   = "5.90-169 kg",
    weight_median  = "70.6 kg",
    sex_female_pct = 46.6,
    race_ethnicity = c(White = 49.5, Black = 6.7, Asian = 38.0, Other = 1.2, Unknown = 4.5),
    disease_state  = paste(
      "Advanced or metastatic solid tumors harboring ALK, ROS1 or NTRK1-3",
      "rearrangements (363 ROS1, 65 NTRK3, 54 NTRK1, 35 ALK, 9 NTRK2), plus",
      "118 healthy volunteers"
    ),
    dose_range     = paste(
      "Oral capsules; the approved and simulated regimen is 160 mg once daily",
      "for 14 days followed by 160 mg twice daily"
    ),
    pediatric      = paste(
      "24 pediatric patients from the CARE study (NCT04094610): 16 under 12",
      "years and 8 adolescents aged 12 to under 18 years"
    ),
    renal_function = "448 normal (69.6%), 157 mild (24.4%), 33 moderate (5.1%) by eGFR",
    hepatic_function = "582 normal (90.4%), 59 mild (9.2%), 1 moderate (0.2%)",
    notes          = paste(
      "Baseline demographics are Du 2025 Table 1. Eight studies: six phase 1",
      "trials in healthy volunteers plus TRIDENT-1 (NCT03093116, adults) and",
      "CARE (NCT04094610, pediatric). This is the Stage II analysis, which",
      "supported the June 2024 US accelerated approval for adult and",
      "adolescent patients with NTRK-fusion-positive solid tumors; the Stage",
      "I model that supported the November 2023 ROS1-positive NSCLC approval",
      "is a different, dose-driven autoinduction model and is NOT packaged",
      "here (its control stream is Du 2025 supplement Data S2 and its",
      "estimates are Table S2)."
    )
  )

  ini({
    # -----------------------------------------------------------------
    # All values are the Stage II FINAL estimates from Du 2025 Table 2.
    # The supplement control stream (Data S3) is the authority for model
    # STRUCTURE only -- its $THETA / $OMEGA blocks carry INITIAL estimates
    # (e.g. CL 6.933 vs final 7.1, CLMAX 1.468 vs final 1.59, age effect
    # -0.336 vs final -0.292) and must not be read as results.
    # -----------------------------------------------------------------

    # Disposition. Typical values are for the reference subject: a 70 kg
    # PATIENT (the healthy-volunteer central volume is a covariate effect).
    lcl <- log(7.1);   label("Baseline (uninduced) clearance (L/h)")               # Table 2: CL = 7.1 (RSE 4.5%)
    lvc <- log(19.8);  label("Central volume of distribution, patients (L)")       # Table 2: VC = 19.8 (RSE 9.5%)
    lvp <- log(221);   label("Peripheral volume of distribution (L)")              # Table 2: VP = 221 (RSE 4.0%)
    lq  <- log(4.98);  label("Intercompartmental clearance (L/h)")                 # Table 2: Q = 4.98 (RSE 6.5%)

    # Absorption rate constant: one typical value per prandial level.
    lkaFasted    <- log(0.0541); label("Absorption rate constant, fasted (1/h)")           # Table 2: KAFASTED = 0.0541 (RSE 3.8%)
    lkaFed       <- log(0.124);  label("Absorption rate constant, fed (1/h)")              # Table 2: KAFED = 0.124 (RSE 3.8%)
    lkaModFasted <- log(0.141);  label("Absorption rate constant, modified fasted (1/h)")  # Table 2: KAMODIFIED = 0.141 (RSE 5.8%)
    lkaUnknown   <- log(0.193);  label("Absorption rate constant, unknown prandial state (1/h)")  # Table 2: KAUnknown = 0.193 (RSE 3.6%)

    # Bioavailability: one typical value per prandial level, carried on the
    # LOGIT scale because Data S3 adds the IIV there
    # (LGT_TVF1 = LOG(TVF1/(1-TVF1)); F1 = 1/(1+EXP(-(LGT_TVF1 + ETA(6))))).
    # Table 2 prints the natural-scale values reproduced in each comment.
    logitfdepotFasted    <- logit(0.52);  label("Bioavailability, fasted (logit scale)")           # Table 2: F1FASTED = 0.52 (RSE 5.3%)
    logitfdepotFed       <- logit(0.76);  label("Bioavailability, fed (logit scale)")              # Table 2: F1FED = 0.76 (RSE 4.8%)
    logitfdepotModFasted <- logit(0.639); label("Bioavailability, modified fasted (logit scale)")  # Table 2: F1MODIFIED = 0.639 (RSE 6.1%)
    logitfdepotUnknown   <- logit(0.533); label("Bioavailability, unknown prandial state (logit scale)")  # Table 2: F1Unknown = 0.533 (RSE 5.6%)

    ltlag <- log(0.421); label("Absorption lag time (h)")                          # Table 2: ALAG1 = 0.421 (RSE 1.7%)

    # Allometric scaling about a 70 kg reference. ESTIMATED in Stage II
    # (they were FIXED at 0.75 / 1 in the Stage I model, Data S2).
    e_wt_cl <- 0.477; label("Allometric exponent of body weight on CL and Q (unitless)")    # Table 2: CLQWT = 0.477 (RSE 14.3%)
    e_wt_vc <- 0.962; label("Allometric exponent of body weight on Vc and Vp (unitless)")   # Table 2: VCVPWT = 0.962 (RSE 5.4%)

    # Cohort effect on central volume. LINEAR deviation (1 + theta), per
    # Data S3 `IF(HV.EQ.1) V2HV = (1 + THETA(20))`; patients are the
    # reference, so healthy volunteers get vc * (1 - 0.854) = 0.146 * vc.
    e_dis_healthy_vc <- -0.854; label("Linear-deviation effect of healthy-volunteer status on Vc (unitless)")  # Table 2: VCPOP = -0.854 (RSE 2.6%)

    # Empirical autoinduction of clearance. The induction is an Emax function
    # of the held pre-dose (trough) concentration MULTIPLIED by a hyperbolic
    # function of time since the first dose, and it acts on the LOG scale:
    #   cl = cl_baseline * exp(cl_time_max * Emax(Ctrough) * Hyperbola(t))
    # so cl_time_max is a log-scale multiplier, not a clearance. Du 2025
    # Discussion: "The CLMAX was estimated to be 4.9 times the baseline CL",
    # which reproduces as exp(1.59) = 4.90.
    # Carried on the log scale in ini() and back-transformed in model(),
    # following the cl_time_ family convention of Kuchimanchi_2024_dostarlimab.R
    # and Masters_2022_avelumab.R.
    lcl_time_max  <- log(1.59);     label("log maximum induction exponent; cl reaches exp(cl_time_max) = 4.9-fold baseline (unitless)")  # Table 2: CLMAX = 1.59 (RSE 4.0%)
    lcl_ec50      <- log(77);       label("log trough concentration giving half-maximal induction (ng/mL)")   # Table 2: EC50 = 77 (RSE 9.3%)
    lcl_conc_hill <- fixed(log(1)); label("log Hill coefficient on trough concentration in the induction term (unitless)")  # Table 2: GAMMA = 1 FIX
    lcl_t50       <- log(47.2);     label("log time since first dose giving half-maximal induction (h)")      # Table 2: TC50 = 47.2 (RSE 22.6%)

    e_age_cl_time_max <- -0.292; label("Power exponent of age (clamped at 18 y) on the maximum induction exponent (unitless)")  # Table 2: CLMAXAGE = -0.292 (RSE 9.2%)

    # Between-subject variability. Table 2 reports these as VARIANCES: the
    # column is headed "omega^2" and the values match the $OMEGA diagonal of
    # Data S3 (e.g. 0.0392 initial vs 0.0391 final for the healthy-volunteer
    # clearance term).
    etalclHv       ~ 0.0391;  label("IIV on clearance, healthy volunteers (variance, log scale)")  # Table 2: omega^2 CLH = 0.0391 (RSE 20.6%, shrinkage 17.3%)
    etalclPt       ~ 0.289;   label("IIV on clearance, patients (variance, log scale)")            # Table 2: omega^2 CLP = 0.289 (RSE 8.8%, shrinkage 14.6%)
    etalvc         ~ 0.768;   label("IIV on central volume (variance, log scale)")                 # Table 2: omega^2 VC = 0.768 (RSE 10.1%, shrinkage 32.6%)
    etalq          ~ 0.458;   label("IIV on intercompartmental clearance (variance, log scale)")   # Table 2: omega^2 Q = 0.458 (RSE 13.0%, shrinkage 31.8%)
    etalka         ~ 0.0961;  label("IIV on absorption rate constant (variance, log scale)")       # Table 2: omega^2 KA = 0.0961 (RSE 13.2%, shrinkage 36.2%)
    etalogitfdepot ~ 0.355;   label("IIV on bioavailability (variance, logit scale)")              # Table 2: omega^2 F1 = 0.355 (RSE 16.3%, shrinkage 40.2%)

    # Residual error. These are STANDARD DEVIATIONS, not variances: Data S3
    # squares each one, `W = SQRT(THETA(14)**2*IPRED**2 + THETA(15)**2)` with
    # $SIGMA fixed to 1, so the proportional terms are fractional CVs and the
    # additive terms are in ng/mL. (Table 2's footnote a, which says random
    # effects and residual error are shown "as variance", is correct for the
    # omega^2 rows above but not for these four; the control stream settles
    # it. Reading 0.424 as a variance would state a 65% proportional error
    # instead of the actual 42.4%.)
    propSdHv <- 0.335;          label("Proportional residual SD, healthy volunteers (fraction)")  # Table 2: Proportional error in healthy volunteers = 0.335 (RSE 1.3%)
    addSdHv  <- fixed(0.00001); label("Additive residual SD, healthy volunteers (ng/mL)")         # Table 2: Additive error in healthy volunteers = 0.00001 FIX
    propSdPt <- 0.424;          label("Proportional residual SD, patients (fraction)")            # Table 2: Proportional error in patients = 0.424 (RSE 1.3%)
    addSdPt  <- 10.6;           label("Additive residual SD, patients (ng/mL)")                   # Table 2: Additive error in patients = 10.6 (RSE 5.3%)
  })

  model({
    # ---------------------------------------------------------------
    # 1. Prandial-level selection.
    #
    # The source column FED2 has four mutually exclusive levels (0 fasted,
    # 1 fed, 2 modified fasted, -99 unknown), each with its OWN ka and its
    # own F1. They are reconstructed here from three canonical binary
    # columns, with fasted as the complement. FASTED_STRICT is consulted
    # only on records that are neither fed nor unknown.
    # ---------------------------------------------------------------
    isFed       <- FED
    isUnknown   <- FED_MISSING
    isModFasted <- (1 - FED) * (1 - FED_MISSING) * (1 - FASTED_STRICT)
    isFasted    <- (1 - FED) * (1 - FED_MISSING) * FASTED_STRICT

    lkaSel <-
      lkaFasted * isFasted +
      lkaFed * isFed +
      lkaModFasted * isModFasted +
      lkaUnknown * isUnknown

    logitfdepotSel <-
      logitfdepotFasted * isFasted +
      logitfdepotFed * isFed +
      logitfdepotModFasted * isModFasted +
      logitfdepotUnknown * isUnknown

    # ---------------------------------------------------------------
    # 2. Empirical trough-driven autoinduction of clearance
    #    (Du 2025 Equations (4) and (5)).
    #
    #   CLMAX_i   = cl_time_max                              , age >= 18
    #             = cl_time_max * (AGE/18)^e_age_cl_time_max , age <  18
    #   induction = CLMAX_i * CONC^g/(EC50^g + CONC^g) * TIME/(TC50 + TIME)
    #   TVCL      = CL_baseline * exp(induction)
    #
    # `min(AGE, 18)` reproduces the two-branch IF exactly: it returns 18 for
    # every adult, so the power term collapses to 1.
    # TIME is time since the FIRST dose, i.e. tafd().
    # ---------------------------------------------------------------
    cl_time_max  <- exp(lcl_time_max)
    cl_ec50      <- exp(lcl_ec50)
    cl_conc_hill <- exp(lcl_conc_hill)
    cl_t50       <- exp(lcl_t50)

    clTimeMaxInd <- cl_time_max * (min(AGE, 18) / 18)^e_age_cl_time_max

    indConc <- ctroughHold^cl_conc_hill /
      (cl_ec50^cl_conc_hill + ctroughHold^cl_conc_hill)
    # tafd() is assigned once and then reused: calling an rxode2 time function
    # twice inside a single expression fails to parse in the ini()/model()
    # interface (`z <- tafd()/(t50 + tafd())` is rejected).
    #
    # tafd() returns NA until the first dose record, and an NA here would
    # propagate through cl into every state derivative and blank the whole
    # solve. The guard is also exact rather than merely defensive: NONMEM
    # evaluates TIME/(TC50 + TIME) with the dataset clock, which is 0 at and
    # before the first dose, so the pre-dose value of this term is 0.
    tsfd <- tafd()
    if (tsfd > 0) {
      indTime <- tsfd / (cl_t50 + tsfd)
    } else {
      indTime <- 0
    }
    induction <- clTimeMaxInd * indConc * indTime

    # ---------------------------------------------------------------
    # 3. Individual parameters.
    #
    # Clearance IIV is drawn from one of two distributions depending on
    # cohort; Data S3 writes this as IF(HV.EQ.0) ... ETA(7) ELSE ETA(1).
    # Exactly one of the two indicator products is 1 for any subject.
    # ---------------------------------------------------------------
    etalclSel <- etalclHv * DIS_HEALTHY + etalclPt * (1 - DIS_HEALTHY)

    wtCl <- (WT / 70)^e_wt_cl
    wtV  <- (WT / 70)^e_wt_vc

    cl <- exp(lcl + etalclSel + induction) * wtCl
    vc <- exp(lvc + etalvc) * (1 + e_dis_healthy_vc * DIS_HEALTHY) * wtV
    vp <- exp(lvp) * wtV
    q  <- exp(lq + etalq) * wtCl
    ka <- exp(lkaSel + etalka)

    fdepot <- expit(logitfdepotSel + etalogitfdepot)
    tlag   <- exp(ltlag)

    # ---------------------------------------------------------------
    # 4. Disposition (ADVAN4 TRANS4: two compartments, first-order input).
    # ---------------------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot)    <- fdepot
    alag(depot) <- tlag

    # Amounts are in mg and vc in L, so central/vc is ug/mL; the factor of
    # 1000 converts to the ng/mL of the bioanalytical assay. Data S3 does the
    # same with S2 = V2/1000.
    Cc <- 1000 * central / vc

    # ---------------------------------------------------------------
    # 5. Sample-and-hold of the pre-dose (trough) concentration.
    #
    # Data S3 captures the induction driver with
    #     IF(NEWIND.NE.2) CONC = 0
    #     IF(EVID.EQ.1) THEN
    #       CONC = A(2)/S2
    #     ELSE
    #       CONC = CONC
    #     ENDIF
    # i.e. CONC starts at 0, is refreshed to the current central
    # concentration at every dosing event, and is held constant in between.
    # NONMEM can do this because $PK executes once per data record; rxode2
    # evaluates the model continuously along the solution, so the latch is
    # realized as a state that is charged only inside a narrow window
    # immediately after each dose and is inert otherwise.
    #
    # The window (0.01 h = 36 s) sits entirely inside the 0.421 h absorption
    # lag, so none of the dose just administered has reached the central
    # compartment yet and the value captured IS the pre-dose trough. The rate
    # (2000 /h) gives 20 window time constants, so the state reaches the
    # target to within exp(-20). Both are numerical implementation constants
    # chosen for this purpose, NOT parameters of the published model, which
    # is why they are literals here rather than entries in ini().
    #
    # Accuracy was measured against the true pre-dose Cc across a five-day
    # 160 mg QD course and an eight-day 160 mg BID course, and is invariant
    # to solver tolerance from atol/rtol 1e-8/1e-6 through 1e-12/1e-10:
    # the held value is high by at most 0.53% (QD) and 1.16% (BID). Because
    # the driver enters through a saturating Emax term, that propagates to
    # under 0.3% on clearance. The vignette re-runs this check as a gate.
    # ---------------------------------------------------------------
    # tad() is NA before the first dose, so the gate is closed until one has
    # been given (the latch correctly stays at its initial value of 0, which
    # is what NONMEM's `IF(NEWIND.NE.2) CONC = 0` establishes).
    tsld <- tad()
    if (tsld >= 0) {
      holdOn <- (tsld < 0.01)
    } else {
      holdOn <- 0
    }
    d/dt(ctroughHold) <- 2000 * holdOn * (Cc - ctroughHold)

    # ---------------------------------------------------------------
    # 6. Cohort-stratified combined proportional plus additive residual
    #    error (Data S3 $ERROR).
    # ---------------------------------------------------------------
    propSd <- propSdHv * DIS_HEALTHY + propSdPt * (1 - DIS_HEALTHY)
    addSd  <- addSdHv  * DIS_HEALTHY + addSdPt  * (1 - DIS_HEALTHY)

    Cc ~ prop(propSd) + add(addSd)
  })
}
