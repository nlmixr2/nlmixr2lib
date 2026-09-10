# Joint parent-metabolite one-compartment population PK model for oral
# fluoxetine and its active metabolite norfluoxetine in Chinese psychiatric
# patients (Han 2025, Pharmaceutics 17:1516;
# doi:10.3390/pharmaceutics17121516).

Han_2025_fluoxetine <- function() {
  description <- paste(
    "Joint parent-metabolite population pharmacokinetic model for oral",
    "fluoxetine and its active metabolite norfluoxetine in 198 Chinese",
    "psychiatric adolescent and adult patients (146 female, 52 male; median",
    "age 17 years, range 12-56) receiving 20-60 mg once daily, developed",
    "from routine steady-state trough therapeutic-drug-monitoring records",
    "(Han 2025). Structural model: two connected one-compartment models --",
    "first-order absorption (Ka fixed at 0.3 1/h from the Panchaud 2011",
    "literature model, because only trough samples were available) into a",
    "fluoxetine central compartment with first-order elimination, the whole",
    "of which (fraction metabolised fixed at 1) feeds a one-compartment",
    "norfluoxetine disposition with its own apparent clearance and apparent",
    "volume. Sex was the sole covariate retained by forward inclusion",
    "(p < 0.01) / backward elimination (p < 0.001): apparent fluoxetine",
    "clearance is 16.5% higher in males than in females, and the Table 4",
    "typical values are the female-reference estimates. Between-subject",
    "variability is exponential on the apparent clearance of each analyte",
    "only; residual error is combined additive-plus-proportional, fitted",
    "separately for parent and metabolite. Because the fit used trough-only",
    "steady-state data, the apparent volumes (particularly the",
    "norfluoxetine volume, RSE 57%) are weakly identified and should be",
    "read as trough-reproducing rather than physiologic. The same paper",
    "externally evaluated the two previously published fluoxetine popPK",
    "models (Panchaud 2011, Wilens 2002) against this dataset and found",
    "both to underpredict at the population level.",
    sep = " "
  )
  reference <- paste(
    "Han B, Xu N, Ma C, Ju G, Xi X, Qian C, Guo N, Liu X, Zhu X, Li C,",
    "Liu L (2025). Bridging Literature and Real-World Evidence: External",
    "Evaluation and Development of Fluoxetine Population Pharmacokinetics",
    "Model. Pharmaceutics 17(12):1516.",
    "doi:10.3390/pharmaceutics17121516.",
    sep = " "
  )
  vignette <- "Han_2025_fluoxetine"
  units    <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Amount units follow the mg dosing unit declared
  # above; the observation variables convert mg/L to the ng/mL reporting
  # scale used throughout Han 2025 (Table 1 LOQ 1 ng/mL; Table 4
  # Add.err.sd in ng/mL; Sect. 3.6 target range 120-500 ng/mL).
  # verified = TRUE: Han 2025 reports CL/F in L/h and V/F in L for both
  # analytes (Table 4) and doses in mg/day (Table 1), which fixes the
  # amount unit of each state.
  compartmentData <- list(
    depot            = list(analyte = "fluoxetine",    units = "mg", specimen = "administration site", verified = TRUE),
    central          = list(analyte = "fluoxetine",    units = "mg", specimen = "plasma",              verified = TRUE),
    central_norfluox = list(analyte = "norfluoxetine", units = "mg", specimen = "plasma",              verified = TRUE)
  )

  covariateData <- list(
    SEXF = list(
      description        = "Biological sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male) in the canonical column. Han 2025's own reference category is FEMALE: the Table 4 footnote gives CL/F,males = CL/F,females x (1 + 0.165) with Females = 0 and Males = 1, so the printed CL/F = 2.91 L/h is the female typical value.",
      notes              = paste(
        "Han 2025 encodes sex as a male-indicator (Table 4 footnote 1:",
        "'Females = 0; Males = 1') with female as the reference category.",
        "To store the covariate under the canonical SEXF (1 = female,",
        "0 = male) while preserving Han's female-reference CL/F of",
        "2.91 L/h verbatim, the effect is applied in model() as",
        "(1 + e_sex_cl * (1 - SEXF)): SEXF = 1 (female) yields factor 1,",
        "SEXF = 0 (male) yields the paper's +16.5% increment. This is the",
        "same construction used by Bajaj_2017_nivolumab.R,",
        "Wada_2023_sparsentan.R and Li_2012_clozapine.R. Sex was the sole",
        "covariate retained by the stepwise search (Sect. 2.5 forward",
        "inclusion p < 0.01, backward elimination p < 0.001). The cohort",
        "was 73.7% female (146 of 198; Table 1), and the effect acts on",
        "fluoxetine clearance only -- no sex effect was retained on the",
        "fluoxetine volume or on either norfluoxetine parameter."
      ),
      source_name        = "Sex (1 = male, 0 = female)"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Screened as a candidate covariate (Sect. 2.5: 'Candidate",
        "covariates were derived from patients' demographic data,",
        "including body weight, sex, and other clinically relevant",
        "parameters') and NOT retained. The Discussion attributes the",
        "null result to the narrow weight distribution of the cohort",
        "(IQR 48.15-67.75 kg), explicitly contrasting it with Wilens",
        "2002, which did find weight significant on both CL/F and V/F.",
        "Cohort weight: median 59 kg, range 35.9-115 (Table 1).",
        "Documented here to preserve the covariate screen without",
        "carrying a declared-but-unused convention warning."
      )
    ),
    AGE = list(
      description = "Subject age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Collected for every subject and used to stratify the cohort",
        "into pediatric (< 18 years) and adult (>= 18 years) groups",
        "(Sect. 3.3), but not retained as a covariate in the final",
        "model. The Discussion argues that CYP2C19 and CYP2D6 are fully",
        "mature beyond ~6 years of age, so no maturation term is",
        "expected in a 12-56 year cohort. Listed among the covariates",
        "the Limitations section says the homogeneous population",
        "prevented identifying ('the relatively homogeneous study",
        "population may have restricted the identification of important",
        "covariates, such as age and body weight')."
      )
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 198L,
    n_observations   = 482L,
    n_studies        = 1L,
    age_range        = "12-56 years (median 17; Table 1)",
    age_median       = "17 years",
    weight_range     = "35.9-115 kg (median 59; Table 1)",
    weight_median    = "59 kg",
    sex_female_pct   = 73.7,
    race_ethnicity   = c(Asian = 100),
    disease_state    = paste(
      "Chinese psychiatric inpatients and outpatients treated with",
      "fluoxetine (predominantly depression; fluoxetine is first-line for",
      "major depressive disorder in adults and adolescents). Retrospective",
      "therapeutic-drug-monitoring records from Hunan Brain Hospital (the",
      "Second People's Hospital of Hunan Province) collected 2021-2024.",
      "Subjects were required to have complete dosing histories,",
      "demographic data, and TDM records."
    ),
    dose_range       = paste(
      "20-60 mg once daily (Table 1 'External dataset' row). Note the",
      "paper is internally inconsistent about the upper bound: Sect. 3.3",
      "states 'The daily dosage in this cohort ranged from 20 to 40 mg',",
      "while Table 1 and the Discussion both give 20-60 mg/day with a",
      "median of 40 mg/day."
    ),
    regions          = "China (Changsha, Hunan Province)",
    age_strata       = paste(
      "Pediatric (< 18 years) n = 102, median age 15 (range 12-17),",
      "median weight 55 kg (range 35.9-96); adult (>= 18 years) n = 92,",
      "median age 21 (range 18-56), median weight 59 kg (range",
      "35.9-115) (Sect. 3.3). The two strata sum to 194 rather than the",
      "198 subjects reported for the whole cohort; the paper does not",
      "reconcile the difference."
    ),
    sampling         = paste(
      "Trough-only: 'All plasma samples were collected at steady state,",
      "prior to the next scheduled dose' (Sect. 2.3); Table 1 records",
      "'All concentration was trough concentration'. This is why Ka",
      "could not be estimated and was fixed."
    ),
    notes            = paste(
      "241 fluoxetine and 241 norfluoxetine plasma concentrations (482",
      "observations total) from 198 subjects. Assay: validated LC-MS/MS",
      "(HPLC-MS/MS) with a lower limit of quantification of 1 ng/mL.",
      "Median observed fluoxetine trough 140.97 ng/mL (range",
      "7.43-979.93); median observed norfluoxetine trough 129.2 ng/mL",
      "(range 9.6-490.41) (Sect. 3.3). Fitted in NONMEM 7.5 with PsN",
      "5.2.6; internal validation by goodness-of-fit plots,",
      "prediction-corrected VPC, and a 1000-replicate stratified",
      "bootstrap (Table 4 reports bootstrap medians and 95% CIs). The",
      "same dataset was used both for the external evaluation of the two",
      "published literature models and for developing this model."
    )
  )

  ini({
    # ---- Structural typical values (Han 2025 Table 4) ----
    # Every fluoxetine typical value in Table 4 is the estimate for the
    # FEMALE reference (Table 4 footnote 1: CL/F,males = CL/F,females x
    # (1 + 0.165); Females = 0, Males = 1). The male value is recovered
    # multiplicatively in model().

    # Absorption rate constant, fixed. Sect. 2.5: 'Given the clinical
    # data limitation (only steady-state trough concentrations were
    # available), the absorption process could not be reliably
    # estimated. Therefore, the absorption rate constant (ka) was fixed
    # to literature-derived values (0.3 h-1) from the model repository
    # to ensure model identifiability.' The 0.3 1/h value is the fixed
    # Ka of the Panchaud 2011 model tabulated in Table 2.
    lka <- fixed(log(0.3)); label("Fluoxetine absorption rate constant Ka (1/h); from the Panchaud 2011 literature model, see Sect. 2.5")  # Han 2025 Table 4: Ka (fixed) = 0.3 1/h

    # Fluoxetine apparent oral clearance and apparent volume of
    # distribution at the female reference.
    lcl <- log(2.91); label("Apparent oral clearance of fluoxetine CLP/F at the female reference (L/h)")  # Han 2025 Table 4: CL/F = 2.91 L/h (RSE 23%); bootstrap median 2.76 [1.53-4.13]
    lvc <- log(24.9); label("Apparent volume of distribution of fluoxetine VP/F (L)")                     # Han 2025 Table 4: V/F = 24.9 L (RSE 38%); bootstrap median 22.76 [9.06-38.48]

    # Fraction of fluoxetine elimination routed to norfluoxetine, fixed
    # at 1 (Table 4 'FM (fixed)' = 1; the abbreviations line defines FM
    # as 'the fraction of metabolism from fluoxetine to norfluoxetine').
    # With FM = 1 the whole of the parent's apparent elimination flux
    # feeds the metabolite compartment, so any true fraction < 1 and any
    # molar-mass correction (fluoxetine 309.33 vs norfluoxetine 295.31
    # g/mol) are absorbed into the apparent CLM/F and VM/F below.
    fm <- fixed(1); label("Fraction of fluoxetine elimination forming norfluoxetine (unitless, in (0, 1])")  # Han 2025 Table 4: FM (fixed) = 1

    # Norfluoxetine apparent clearance and apparent volume of
    # distribution. Both are 'apparent' in the parent-bioavailability
    # sense because norfluoxetine is never dosed directly.
    lcl_norfluox <- log(3.24); label("Apparent clearance of norfluoxetine CLM/F (L/h)")                # Han 2025 Table 4 (Norfluoxetine block): CL/F = 3.24 L/h (RSE 20%); bootstrap median 3.06 [1.77-4.53]
    lvc_norfluox <- log(1.52); label("Apparent volume of distribution of norfluoxetine VM/F (L)")      # Han 2025 Table 4 (Norfluoxetine block): V/F = 1.52 L (RSE 57%); bootstrap median 1.17 [0.67-1.98]

    # ---- Covariate effect (Han 2025 Table 4) ----
    # Linear multiplicative fractional form, taken verbatim from the
    # Table 4 footnote: CL/F,males = CL/F,females x (1 + 0.165), with
    # Males = 1 and Females = 0. Stored against the canonical SEXF
    # (1 = female) column by applying it as (1 - SEXF); see
    # covariateData[[SEXF]]$notes. Sect. 3.5: 'clearance in males was
    # estimated to be 16.5% higher than in females'.
    e_sex_cl <- 0.165; label("Male-sex effect on fluoxetine CLP/F (fraction; applied as (1 + e_sex_cl * (1 - SEXF)))")  # Han 2025 Table 4: Sex effects on CL/F = 16.50% (RSE 44%); bootstrap median 16.43 [6.64-27.21]

    # ---- Inter-individual variability (Han 2025 Table 4) ----
    # Sect. 2.5: 'Inter-individual variability (IIV) was modeled using
    # exponential error structures.' Table 4 reports IIV as a percent
    # for the apparent clearance of each analyte only -- there is no
    # IIV row for either volume, for Ka (fixed) or for FM (fixed).
    # Converted to the NONMEM OMEGA variance with the log-normal
    # identity omega^2 = log(CV^2 + 1):
    #   fluoxetine    CL/F: CV = 31.6% -> omega^2 = log(0.316^2 + 1) = 0.0951793
    #   norfluoxetine CL/F: CV = 20.9% -> omega^2 = log(0.209^2 + 1) = 0.0427539
    etalcl          ~ 0.0951793  # Han 2025 Table 4: 'IIV CL/F, %' fluoxetine = 31.6 (RSE 27%) [eta-shrinkage 10%] -> omega^2 = log(0.316^2 + 1)
    etalcl_norfluox ~ 0.0427539  # Han 2025 Table 4: 'IIV CL/F, %' norfluoxetine = 20.9 (RSE 41%) [eta-shrinkage 48%] -> omega^2 = log(0.209^2 + 1)

    # ---- Residual error (Han 2025 Table 4) ----
    # Sect. 2.5: 'residual unexplained variability (RUV) was described
    # by a combined additive-proportional error model', fitted
    # separately for the parent and the metabolite. Table 4 names the
    # rows 'Prop.err.sd, %' and 'Add.err.sd, ng/mL', i.e. both are
    # already on the SD scale and need no variance-to-SD conversion.
    propSd          <- 0.341; label("Proportional residual error for fluoxetine (fraction)")        # Han 2025 Table 4: Prop.err.sd = 34.1% (RSE 12%) [eps-shrinkage 20.8%]
    addSd           <- 14.9;  label("Additive residual error for fluoxetine (ng/mL)")               # Han 2025 Table 4: Add.err.sd = 14.9 ng/mL (RSE 49%) [eps-shrinkage 20.8%]
    propSd_norfluox <- 0.305; label("Proportional residual error for norfluoxetine (fraction)")     # Han 2025 Table 4: Prop.err.sd = 30.5% (RSE 16%) [eps-shrinkage 37.3%]
    addSd_norfluox  <- 22.9;  label("Additive residual error for norfluoxetine (ng/mL)")            # Han 2025 Table 4: Add.err.sd = 22.9 ng/mL (RSE 31%) [eps-shrinkage 37.3%]
  })

  model({
    # Han 2025 encodes sex as a male-indicator with FEMALE as the
    # reference category (Table 4 footnote 1); (1 - SEXF) reproduces the
    # paper's Males = 1 column while keeping the canonical SEXF
    # (1 = female) storage convention.
    sex_male <- 1 - SEXF

    # Fluoxetine disposition parameters. IIV is exponential on apparent
    # clearance only; the sex effect multiplies the female-reference
    # typical value.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (1 + e_sex_cl * sex_male)
    vc <- exp(lvc)

    # Norfluoxetine disposition parameters. No covariate effect was
    # retained on either.
    cl_norfluox <- exp(lcl_norfluox + etalcl_norfluox)
    vc_norfluox <- exp(lvc_norfluox)

    # Two connected one-compartment models (Sect. 3.5: 'The final joint
    # parent-metabolite model for fluoxetine and norfluoxetine was
    # parameterized as two connected one-compartment models').
    # First-order absorption into the fluoxetine central compartment;
    # a fraction FM (fixed at 1) of the parent's elimination flux is
    # delivered to the norfluoxetine central compartment, which is
    # eliminated with its own apparent clearance.
    d/dt(depot)            <- -ka * depot
    d/dt(central)          <-  ka * depot - (cl / vc) * central
    d/dt(central_norfluox) <-  fm * (cl / vc) * central -
      (cl_norfluox / vc_norfluox) * central_norfluox

    # Plasma concentrations. The compartments hold mg and the volumes
    # are in L, so central / vc is mg/L = ug/mL; the factor 1000
    # converts to the ng/mL scale Han 2025 reports throughout (Table 1
    # LOQ 1 ng/mL, Table 4 Add.err.sd in ng/mL, Sect. 3.6 active-moiety
    # target 120-500 ng/mL). The clinical "active moiety" monitored in
    # TDM is the arithmetic sum Cc + Cc_norfluox.
    Cc          <- 1000 * central          / vc
    Cc_norfluox <- 1000 * central_norfluox / vc_norfluox

    # Combined additive-plus-proportional residual error, fitted
    # separately for parent and metabolite (Sect. 2.5, Table 4).
    Cc          ~ add(addSd)          + prop(propSd)
    Cc_norfluox ~ add(addSd_norfluox) + prop(propSd_norfluox)
  })
}
