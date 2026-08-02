Ganesan_2023_tebipenem <- function() {
  description <- paste(
    "Two-compartment population PK model with two absorption transit",
    "compartments for tebipenem, the active moiety of the oral carbapenem",
    "pro-drug tebipenem pivoxil hydrobromide (TBP-PI-HBr), pooled across",
    "three phase 1 studies and the phase 3 ADAPT-PO trial in adults with",
    "complicated urinary tract infection / acute pyelonephritis (Ganesan",
    "2023; 746 subjects, 3448 plasma concentrations). Apparent oral",
    "clearance is split into an additive non-renal arm (power function of",
    "creatinine clearance) and a renal arm (sigmoidal Hill function of",
    "creatinine clearance) that drives a cumulative urine compartment;",
    "the summed CL/F carries a linear body-surface-area effect. Central",
    "volume scales with height and peripheral volume with body surface",
    "area, both shifted by infection status. The absorption rate constant",
    "switches on fed status, is shifted by infection status, and carries a",
    "dose effect confined to the crossover thorough-QT study (study 104).",
    "Interindividual variability on CL/F is cohort-specific (healthy",
    "subjects vs infected patients) and Ka carries both IIV and",
    "two-occasion interoccasion variability."
  )
  reference <- paste(
    "Ganesan H, Gupta VK, Safir MC, Bhavnani SM, Talley AK, Melnick D,",
    "Rubino CM. Population Pharmacokinetic Analyses for Tebipenem after",
    "Oral Administration of Pro-Drug Tebipenem Pivoxil Hydrobromide.",
    "Antimicrob Agents Chemother. 2023 Jun 15;67(6):e0145122.",
    "doi:10.1128/aac.01451-22."
  )
  vignette <- "Ganesan_2023_tebipenem"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description        = "Baseline creatinine clearance, Cockcroft-Gault, normalized to body surface area",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Ganesan 2023 Table 1 footnote b: calculated from baseline serum",
        "creatinine, age, sex, and body weight using the Cockcroft-Gault",
        "equation and then normalized to body surface area. Pooled median",
        "74.2 mL/min/1.73 m^2 (range 6.90-192). Enters the model twice:",
        "(a) as a power function on the non-renal clearance arm,",
        "(CRCL / 79.32)^0.722 (Eq. 1) and (b) as a sigmoidal Hill",
        "function on the renal clearance arm with CRCL50 = 44.7 and Hill",
        "= 2.13 (Eq. 2). The 79.32 mL/min/1.73 m^2 normalizer is printed",
        "in Eq. 1; it is not the pooled median (74.2) and is presumably",
        "the analysis-set mean.  The paper describes this two-arm",
        "renal-function structure as 'somewhat unconventional' but",
        "necessary to avoid bias in subjects with severe renal impairment",
        "(Ganesan 2023 Discussion)."
      ),
      source_name        = "CLcr"
    ),
    BSA = list(
      description        = "Baseline body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Pooled median 1.86 m^2 (range 1.32-2.70; Ganesan 2023 Table 1),",
        "which is the centering value used in both Eq. 3 and Eq. 5.",
        "Enters as a linear centered slope on total CL/F,",
        "1 + 0.479 * (BSA - 1.86) (Eq. 3), and on Vp/F,",
        "1 + 0.491 * (BSA - 1.86) (Eq. 5). Note that the slopes are per",
        "m^2 and are NOT power exponents."
      ),
      source_name        = "BSA"
    ),
    HT = list(
      description        = "Baseline body height",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Pooled median 169 cm (range 110-202; Ganesan 2023 Table 1),",
        "which is the reference value in Eq. 4. Enters Vc/F as a power",
        "function (HT / 169)^2.09. Height, not weight, was the body-size",
        "descriptor retained on Vc/F."
      ),
      source_name        = "HTCM"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-participant cohort indicator (1 = healthy phase 1 subject, 0 = infected phase 3 patient)",
      units              = "(binary)",
      type               = "binary",
      reference_category = paste(
        "0 (infected patient with cUTI / AP). Ganesan 2023 Eq. 4-6 write",
        "the contrast as (1 - Infected); the canonical register codes the",
        "same indicator as DIS_HEALTHY = 1 - Infected, so the structural",
        "typical values lvc, lvp, and lka are the infected-patient values",
        "and DIS_HEALTHY shifts them to the healthy-subject state."
      ),
      notes              = paste(
        "Time-fixed per subject; healthy = the 99 phase 1 subjects",
        "(studies 101, 102, 104), infected = the 647 ADAPT-PO patients.",
        "Three multiplicative effects, all from Ganesan 2023 Table 2 and",
        "Eq. 4-6: Vc/F * (1 - 0.290 * DIS_HEALTHY), Vp/F *",
        "(1 - 0.245 * DIS_HEALTHY), and Ka * (1 + 0.368 * DIS_HEALTHY).",
        "The same indicator also selects which of the two CL/F IIV",
        "variances applies (24.8 %CV healthy vs 57.2 %CV infected;",
        "Ganesan 2023 Table 2). Note that the Discussion states",
        "absorption was 'faster in infected patients' while Eq. 6 as",
        "printed makes Ka 36.8% FASTER in healthy subjects; the equation",
        "is implemented here (see the vignette Errata)."
      ),
      source_name        = "Infected (reverse-coded)"
    ),
    FED = list(
      description        = "Fed-vs-fasted state at the dose record (1 = fed, 0 = fasted)",
      units              = "(binary)",
      type               = "binary",
      reference_category = paste(
        "0 (fasted). Ganesan 2023 Eq. 6 uses the complementary FAST flag;",
        "FED = 1 - FAST, so the structural lka is the fasted typical value",
        "3.04 1/h and the fed state is reached via e_fed_ka."
      ),
      notes              = paste(
        "Only the food-effect cohorts of study 101 contribute fed",
        "records. Food slowed but did not reduce absorption: Ka is",
        "3.04 1/h fasted and 1.23 1/h fed (Ganesan 2023 Eq. 6), a 60%",
        "reduction in rate with no effect on extent. Ganesan 2023 Table 2",
        "labels these two rows the other way round ('Ka (fasted) 1.23',",
        "'Ka (fed) 3.04'); Eq. 6, the Results narrative, and Figure 1",
        "(fasted profiles peak near 0.75 h, fed near 2 h) all agree that",
        "the fasted value is the larger one, so the Table 2 row labels are",
        "treated as transposed (see the vignette Errata)."
      ),
      source_name        = "FAST (reverse-coded)"
    ),
    DOSE_TBPPI_MG = list(
      description        = "Administered TBP-PI-HBr dose level for the record",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Used only inside the study-104 absorption term of Eq. 6,",
        "1 - 0.478 * (DOSE_TBPPI_MG / 1200) * STUDY_SPR994_104, so the",
        "reference",
        "value is the 1200 mg arm of the thorough-QT crossover study.",
        "Ganesan 2023 Discussion: the dose effect on absorption rate was",
        "detectable only in study 104, where the same 24 subjects",
        "received 600 mg and 1200 mg in crossover fashion; it was not",
        "apparent across the 100-900 mg single doses of study 101. Set",
        "DOSE_TBPPI_MG to the administered milligram amount on every",
        "record; it is",
        "inert unless STUDY_SPR994_104 = 1."
      ),
      source_name        = "DOSEMG"
    ),
    STUDY_SPR994_104 = list(
      description        = "Study SPR994-104 (thorough-QT crossover) cohort indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (studies 101, 102, and ADAPT-PO / SPR994-301)",
      notes              = paste(
        "Ganesan 2023 Eq. 6 flag variable S104. Gates the dose effect on",
        "Ka so that it applies only to the 24 subjects of the four-way",
        "crossover thorough-QT study (NCT04238195). Time-fixed per",
        "subject."
      ),
      source_name        = "S104"
    ),
    OCC = list(
      description        = "Occasion index for interoccasion variability on Ka",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Ganesan 2023 Table 2 reports two IOV variances on Ka, 'IOV Ka",
        "(Occ. 1)' and 'IOV Ka (Occ. 2)', both 0.201 (44.8 %CV); the",
        "second is reported without a standard error, i.e. constrained",
        "equal to the first in the NONMEM '$OMEGA BLOCK(1) ... SAME'",
        "idiom. IOV on Ka was required to fit the multiple-ascending-dose",
        "cohorts of study 101 (Ganesan 2023 Results). Decomposed inside",
        "model() into the binary indicators oc1 and oc2; supply OCC = 1",
        "for single-occasion simulations."
      ),
      source_name        = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Baseline age",
      units              = "years",
      type               = "continuous",
      notes              = paste(
        "Pooled median 60.0 years (range 18-91; Ganesan 2023 Table 1).",
        "Age on Vp/F survived forward selection and backward elimination",
        "but was removed from the final model after the SIR analysis",
        "because the effect magnitude was trivial (0.00224 unit change in",
        "Vp/F per year of age); no other age effect was retained. The",
        "paper reports no dose adjustment is warranted on the basis of",
        "age beyond the renal-function adjustment."
      )
    ),
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      notes              = paste(
        "Pooled median 76.0 kg (range 42.0-142; Ganesan 2023 Table 1).",
        "Screened in the covariate analysis but not retained: body size",
        "entered the final model through BSA (CL/F, Vp/F) and height",
        "(Vc/F) instead."
      )
    ),
    BMI = list(
      description        = "Baseline body mass index",
      units              = "kg/m^2",
      type               = "continuous",
      notes              = paste(
        "Pooled median 26.3 kg/m^2 (range 15.3-57.9; Ganesan 2023",
        "Table 1). BMI on Ka survived forward selection but was removed",
        "after the SIR analysis showed the 'true' value of the",
        "corresponding parameter was not significantly different from 0",
        "(Ganesan 2023 Results)."
      )
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      notes              = paste(
        "51.3% female overall (Ganesan 2023 Table 1). Tested but not a",
        "statistically significant predictor of interindividual",
        "variability in tebipenem PK; the small apparent difference in",
        "Vss/F between sexes is attributed by the authors to the lower",
        "height and BSA of the female subjects, which the retained",
        "body-size covariates already capture."
      )
    ),
    RACE_BLACK = list(
      description        = "Black race indicator",
      units              = "(binary)",
      type               = "binary",
      notes              = paste(
        "2.80% of the pooled population (Ganesan 2023 Table 1). Race was",
        "tested in the covariate analysis but not retained; the study",
        "population was 95.9% White overall and 98.6% White in ADAPT-PO,",
        "so the analysis had little power to detect race effects."
      )
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator",
      units              = "(binary)",
      type               = "binary",
      notes              = paste(
        "1.20% of the pooled population (Ganesan 2023 Table 1). Tested",
        "but not retained; see the RACE_BLACK note on the limited racial",
        "diversity of the analysis population."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 746L,
    n_studies      = 4L,
    age_range      = "18-91 years",
    age_median     = "60.0 years",
    weight_range   = "42.0-142 kg",
    weight_median  = "76.0 kg",
    height_range   = "110-202 cm",
    height_median  = "169 cm",
    bsa_range      = "1.32-2.70 m^2",
    bsa_median     = "1.86 m^2",
    sex_female_pct = 51.3,
    race_ethnicity = c(White = 95.9, Black = 2.80, Asian = 1.20, Other = 0.13),
    renal_function = paste(
      "Creatinine clearance (Cockcroft-Gault, BSA-normalized) median 74.2",
      "mL/min/1.73 m^2, range 6.90-192. Study 102 enrolled by renal",
      "function stratum (normal, mild, moderate, severe, and end-stage",
      "renal disease on hemodialysis; median CLcr 42.2). ADAPT-PO",
      "enrolled patients with CLcr as low as 30 mL/min."
    ),
    disease_state  = paste(
      "Pooled: 99 healthy adults or adults with varying degrees of renal",
      "impairment from three phase 1 studies (SPR994-101 single and",
      "multiple ascending dose, SPR994-102 renal impairment, SPR994-104",
      "thorough QT) and 647 patients with complicated urinary tract",
      "infection or acute pyelonephritis from the phase 3 ADAPT-PO trial",
      "(SPR994-301)."
    ),
    dose_range     = paste(
      "TBP-PI-HBr immediate-release tablet, oral. Study 101: 100, 300,",
      "600, 900 mg single dose and 300 or 600 mg q8h for 14 days. Study",
      "102: 600 mg single dose. Study 104: 600 or 1200 mg single dose in",
      "four-way crossover. ADAPT-PO: 600 mg q8h for 7-10 days, reduced to",
      "300 mg q8h for baseline CLcr 30-50 mL/min/1.73 m^2."
    ),
    n_observations = "3448 plasma concentrations (1985 from ADAPT-PO patients) plus urine concentrations from 67 phase 1 subjects and 37 phase 3 patients",
    notes          = paste(
      "NONMEM 7.4, first-order conditional estimation with interaction.",
      "Demographics in Ganesan 2023 Table 1; study designs in supplemental",
      "Table S1; final parameter estimates and sampling-importance-",
      "resampling statistics in Table 2; structural schematic in",
      "supplemental Figure S3. Doses are recorded as milligrams of",
      "TBP-PI (pro-drug), so all apparent parameters are conditioned on",
      "the pro-drug-to-tebipenem mass ratio and on the fraction converted",
      "(expected to be near 100%, since pro-drug was not measurable",
      "systemically). Whole-blood concentrations were converted to plasma",
      "with a factor of 3.6."
    )
  )

  ini({
    # ---------------------------------------------------------------
    # Structural parameters. All values from Ganesan 2023 Table 2
    # ('Final model / Estimate' column), cross-checked against the
    # covariate equations 1-6 printed in the Results section.
    # ---------------------------------------------------------------

    # Absorption. Ka is a fasted/fed switch (Eq. 6). lka carries the
    # FASTED typical value; e_fed_ka moves it to the fed value. See the
    # FED covariateData note: Table 2's two Ka row labels are transposed
    # relative to Eq. 6, the Results narrative ("the presence of food
    # slowed the rate of tebipenem absorption"), and Figure 1.
    lka        <- log(3.04);      label("Absorption / transit rate constant Ka in the fasted state (1/h)")          # Ganesan 2023 Eq. 6 coefficient on FAST = 3.04 1/h (Table 2 row mislabeled as Ka (fed); %SEM 4.21)

    # Clearance. CL/F is the sum of an additive non-renal and renal arm
    # (Eq. 1-3), unlike the more common single-CL parameterization.
    lcl_nonren <- log(15.6);      label("Non-renal clearance CL_NR at the reference CRCL of 79.32 mL/min/1.73 m^2 (L/h)")  # Ganesan 2023 Table 2 CL_NR = 15.6 L/h (%SEM 5.87; SIR mean 15.7, 90% CI [14.4, 17.1]). Eq. 1 prints 10.2 for this coefficient -- see the vignette Errata; 15.6 is used because it is the value in the parameter table, is corroborated three ways by the SIR resample statistics, and is the only one of the two that reproduces the published Table 3 CL/F and AUC0-24 medians for the 600 mg q8h cohort (n = 595)
    lcl_renal  <- log(21.1);      label("Maximum renal clearance CL_R,MAX (L/h)")                                   # Ganesan 2023 Table 2 CL_R,MAX = 21.1 L/h (%SEM 5.96); Eq. 2 numerator

    # Sigmoidal (Hill) relationship between creatinine clearance and the
    # renal clearance arm, Eq. 2:
    #   CL_R = CL_R,MAX * CRCL^hill / (CRCL^hill + CRCL50^hill)
    crcl50_cl_renal <- 44.7;      label("Creatinine clearance giving half-maximal renal clearance CRCL50 (mL/min/1.73 m^2)")  # Ganesan 2023 Table 2 CL_R,CLcr50 = 44.7 (%SEM 7.06)
    hill_cl_renal   <- 2.13;      label("Hill coefficient of the CRCL - CL_R relationship (unitless)")              # Ganesan 2023 Table 2 CL_R,Hill = 2.13 (%SEM 5.13)

    # Distribution.
    lvc        <- log(38.5);      label("Apparent central volume Vc/F in infected patients at HT = 169 cm (L)")     # Ganesan 2023 Table 2 Vc/F = 38.5 L (%SEM 3.81); Eq. 4 leading coefficient
    lvp        <- log(4.84);      label("Apparent peripheral volume Vp/F in infected patients at BSA = 1.86 m^2 (L)")  # Ganesan 2023 Table 2 Vp/F = 4.84 L (%SEM 6.60); Eq. 5 leading coefficient
    lq         <- log(2.23);      label("Apparent distributional clearance CLd/F (L/h)")                            # Ganesan 2023 Table 2 CLd/F = 2.23 L/h (%SEM 11.0)

    # ---------------------------------------------------------------
    # Covariate effects, Ganesan 2023 Eq. 1-6 with the estimates from
    # Table 2. Sign convention: the equations print the infection-status
    # terms as [1 - |theta| * (1 - Infected)]; because the canonical
    # register codes DIS_HEALTHY = 1 - Infected, the coefficients below
    # keep the Table 2 signs and multiply DIS_HEALTHY directly.
    # ---------------------------------------------------------------
    e_crcl_cl_nonren <- 0.722;    label("Power exponent: CRCL on the non-renal clearance arm (unitless)")           # Ganesan 2023 Table 2 CL_NR, CLcr power = 0.722 (%SEM 5.92); Eq. 1
    e_bsa_cl         <- 0.479;    label("Linear BSA slope on total CL/F, centered at 1.86 m^2 (per m^2)")           # Ganesan 2023 Table 2 CL/F:BSA (slope) = 0.479 (%SEM 21.0); Eq. 3
    e_ht_vc          <- 2.09;     label("Power exponent: height on Vc/F, referenced to 169 cm (unitless)")          # Ganesan 2023 Table 2 Vc/F:HTCM (power) = 2.09 (%SEM 24.7); Eq. 4
    e_healthy_vc     <- -0.290;   label("Fractional shift in Vc/F for healthy subjects vs infected patients")       # Ganesan 2023 Table 2 Vc/F:Infection status = -0.29 (%SEM 20.4); Eq. 4 term [1 - 0.290 * (1 - Infected)]
    e_bsa_vp         <- 0.491;    label("Linear BSA slope on Vp/F, centered at 1.86 m^2 (per m^2)")                 # Ganesan 2023 Table 2 Vp/F:BSA (slope) = 0.491 (%SEM 21.3); Eq. 5
    e_healthy_vp     <- -0.245;   label("Fractional shift in Vp/F for healthy subjects vs infected patients")       # Ganesan 2023 Table 2 Vp/F:Infection status = -0.245 (%SEM 19.7); Eq. 5 term [1 - 0.245 * (1 - Infected)]
    e_fed_ka         <- 1.23 / 3.04 - 1; label("Fractional shift in Ka in the fed state vs fasted")                 # Ganesan 2023 Eq. 6: Ka = 3.04 fasted, 1.23 fed; 1.23 / 3.04 - 1 = -0.5954 reproduces both printed values exactly
    e_healthy_ka     <- 0.368;    label("Fractional shift in Ka for healthy subjects vs infected patients")         # Ganesan 2023 Table 2 Ka:Infection status = 0.368 (%SEM 39.4); Eq. 6 term [1 + 0.368 * (1 - Infected)]
    e_dose_ka        <- -0.478;   label("Fractional shift in Ka per 1200 mg of dose, study 104 only")               # Ganesan 2023 Table 2 Ka:Dose effect = -0.478 (%SEM 3.45); Eq. 6 term [1 - 0.478 * (DOSEMG / 1200) * S104]

    # ---------------------------------------------------------------
    # Interindividual variability. Ganesan 2023 Table 2 reports the
    # variance omega^2 with the corresponding %CV in parentheses using
    # the sqrt(omega^2) approximation (e.g. 0.0614 -> 24.8 %CV), so the
    # tabulated numbers are used directly as variances without the
    # log(CV^2 + 1) back-transform.
    #
    # CL/F carries two IIV terms, invoked separately for healthy
    # subjects and infected patients (Ganesan 2023 Results: "separate
    # IIV terms for CL/F were invoked for healthy subjects and infected
    # patients, which resulted in a substantial improvement in the
    # PC-VPC plots"). DIS_HEALTHY selects which one applies.
    #
    # Ganesan 2023 reports NO infection-status effect on the typical
    # value of CL/F -- only on its variance -- so the two cohort
    # anchors below are log-scale multipliers fixed at 0 (i.e. a factor
    # of exactly 1). They exist so that each cohort-specific eta pairs
    # with a fixed-effect parameter of the same name, following the
    # paired-structural-mean idiom used by Klunder_2017_upadacitinib
    # (whose Vc/F means are likewise identical across cohorts because
    # the source reports a cohort split on the variance only).
    # ---------------------------------------------------------------
    lcl_healthy    <- fixed(0);   label("Log-scale CL/F anchor for healthy subjects (hosts the phase 1 IIV)")   # Ganesan 2023 Table 2: no infection-status covariate on typical CL/F; the cohort split affects only the CL/F IIV variance
    lcl_patient    <- fixed(0);   label("Log-scale CL/F anchor for infected patients (hosts the phase 3 IIV)")  # Ganesan 2023 Table 2: no infection-status covariate on typical CL/F; the cohort split affects only the CL/F IIV variance
    etalcl_healthy ~ 0.0614       # Ganesan 2023 Table 2 IIV CL (phase 1) = 0.0614 (24.8 %CV; %SEM 14.8; shrinkage 66.3%)
    etalcl_patient ~ 0.328        # Ganesan 2023 Table 2 IIV CL (phase 3) = 0.328 (57.2 %CV; %SEM 5.24; shrinkage 13.3%)
    etalvc         ~ 0.197        # Ganesan 2023 Table 2 IIV Vc/F = 0.197 (44.4 %CV; %SEM 10.6; shrinkage 41.5%)
    etalvp         ~ 0.0115       # Ganesan 2023 Table 2 IIV Vp/F = 0.0115 (10.7 %CV; %SEM 46.4; shrinkage 77.3%). The Results text describing the full-multivariable-model simplification says the Vp/F IIV term was removed due to excessive shrinkage, but Table 2 (final model) reports it and both the Results and Discussion quote 10.7% for Vp/F as the low end of the final-model IIV range, so it is retained here
    etalka         ~ 0.518        # Ganesan 2023 Table 2 IIV Ka = 0.518 (71.9 %CV; %SEM 8.46; shrinkage 24.0%)

    # Interoccasion variability on Ka. Table 2 lists occasion 1 and
    # occasion 2 with the identical variance and reports a standard
    # error only for occasion 1, i.e. the NONMEM '$OMEGA BLOCK(1) ...
    # SAME' idiom; occasion 2 is therefore fixed equal to occasion 1.
    etaiov_ka_1    ~ 0.201        # Ganesan 2023 Table 2 IOV Ka (Occ. 1) = 0.201 (44.8 %CV; %SEM 37.8; shrinkage 90.2%)
    etaiov_ka_2    ~ fix(0.201)   # Ganesan 2023 Table 2 IOV Ka (Occ. 2) = 0.201 (44.8 %CV; shrinkage 92.0%); no %SEM reported

    # ---------------------------------------------------------------
    # Residual variability. Table 2 reports variances with the SD in
    # parentheses; the SDs are used here (sqrt(0.209) = 0.457,
    # sqrt(0.298) = 0.546, sqrt(482) = 21.9). The additive component of
    # the plasma error was dropped during model development ("removal of
    # the additive portion of the residual variability model for
    # plasma"), so plasma error is proportional only.
    # ---------------------------------------------------------------
    propSd            <- 0.457;   label("Proportional residual error on plasma concentration (fraction)")           # Ganesan 2023 Table 2 RV prop, plasma = 0.209 variance -> 45.7 %CV (%SEM 2.79)
    propSd_urineAmt   <- 0.546;   label("Proportional residual error on cumulative urinary amount (fraction)")      # Ganesan 2023 Table 2 RV prop, urine = 0.298 variance -> 54.6 %CV (%SEM 11.3)
  })

  model({
    # 1. Occasion indicators and the per-occasion IOV eta on Ka.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2

    # 2. Cohort-specific interindividual variability on CL/F. Healthy
    # phase 1 subjects draw the 24.8 %CV eta; infected phase 3 patients
    # draw the 57.2 %CV eta.
    eta_cl <- DIS_HEALTHY * (lcl_healthy + etalcl_healthy) +
      (1 - DIS_HEALTHY) * (lcl_patient + etalcl_patient)

    # 3. Apparent oral clearance, Ganesan 2023 Eq. 1-3.
    #      CL_NR = 15.6 * (CRCL / 79.32)^0.722
    #      CL_R  = 21.1 * CRCL^2.13 / (CRCL^2.13 + 44.7^2.13)
    #      CL/F  = (CL_NR + CL_R) * [1 + 0.479 * (BSA - 1.86)]
    # The BSA factor and the log-normal IIV act on the summed clearance,
    # so they are applied identically to both arms.
    cl_bsa    <- 1 + e_bsa_cl * (BSA - 1.86)
    cl_nonren <- exp(lcl_nonren) * (CRCL / 79.32)^e_crcl_cl_nonren * cl_bsa * exp(eta_cl)
    cl_renal  <- exp(lcl_renal) *
      CRCL^hill_cl_renal / (CRCL^hill_cl_renal + crcl50_cl_renal^hill_cl_renal) *
      cl_bsa * exp(eta_cl)

    # 4. Volumes, Ganesan 2023 Eq. 4-5, and distributional clearance.
    vc <- exp(lvc + etalvc) * (HT / 169)^e_ht_vc * (1 + e_healthy_vc * DIS_HEALTHY)
    vp <- exp(lvp + etalvp) * (1 + e_bsa_vp * (BSA - 1.86)) * (1 + e_healthy_vp * DIS_HEALTHY)
    q  <- exp(lq)

    # 5. Absorption rate constant, Ganesan 2023 Eq. 6. The same Ka
    # governs every transfer in the depot -> transit -> transit ->
    # central chain (supplemental Figure S3).
    ka <- exp(lka + etalka + iov_ka) *
      (1 + e_fed_ka * FED) *
      (1 + e_healthy_ka * DIS_HEALTHY) *
      (1 + e_dose_ka * (DOSE_TBPPI_MG / 1200) * STUDY_SPR994_104)

    # 6. Micro-constants. The renal and non-renal arms are parallel
    # elimination pathways out of the central compartment; only the
    # renal arm feeds the urine compartment.
    kel_nonren <- cl_nonren / vc
    kel_renal  <- cl_renal / vc
    k12        <- q / vc
    k21        <- q / vp

    # 7. ODE system: two absorption transit compartments feeding a
    # two-compartment disposition model, plus a cumulative urine
    # compartment (supplemental Figure S3).
    d/dt(depot)       <- -ka * depot
    d/dt(transit1)    <-  ka * depot    - ka * transit1
    d/dt(transit2)    <-  ka * transit1 - ka * transit2
    d/dt(central)     <-  ka * transit2 - (kel_nonren + kel_renal) * central -
      k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(urine)       <-  kel_renal * central

    # 8. Observations. Cc is plasma tebipenem in mg/L (equivalently
    # ug/mL) with dose in mg and vc in L; urineAmt is the cumulative
    # amount of tebipenem excreted unchanged in urine, in mg.
    Cc <- central / vc
    urineAmt <- urine

    Cc ~ prop(propSd)
    urineAmt ~ prop(propSd_urineAmt)
  })
}
