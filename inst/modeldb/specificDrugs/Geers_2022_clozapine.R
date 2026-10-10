Geers_2022_clozapine <- function() {
  description <- "One-compartment population pharmacokinetic model with first-order absorption for oral clozapine in adults with schizophrenia on a stable dose, developed from plasma and dried-blood-spot concentrations by an iterative two-stage Bayesian procedure in MwPharm for Bayesian-forecasting limited-sampling AUC estimation. The volume of distribution scales linearly with body weight; the absorption rate constant is fixed at its literature prior and bioavailability at 1."
  reference <- paste(
    "Geers LM, Cohen D, Wehkamp LM, van Wattum HJ, Kosterink JGW,",
    "Loonen AJM, Touw DJ. Population pharmacokinetic model and limited",
    "sampling strategy for clozapine using plasma and dried blood spot",
    "samples. Ther Adv Psychopharmacol. 2022;12:20451253211065857.",
    "doi:10.1177/20451253211065857. PMCID: PMC9066631.",
    sep = " "
  )
  vignette <- "Geers_2022_clozapine"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Multiplies the per-kg volume of distribution: V = Vd * WT, a",
        "linear (exponent 1) scaling with no reference weight. Methods,",
        "Population PK model: 'for simplicity clozapine volume of",
        "distribution is expressed in L/kg bodyweight'; Results: 'distribution",
        "volume per kg body weight (Vd)'. Table 3 prints the unit as",
        "'L/kgLBMc', the MwPharm fat-corrected lean body mass",
        "LBMc = LBM + (WT - LBM) * fd. The paper tested and rejected",
        "estimating the fat-distribution fraction ('percentage of distribution",
        "into fatty tissue', Methods; 'into fat tissue also did not improve",
        "the model', Results), and the prose states the per-kg-bodyweight",
        "simplification twice, which is LBMc with fd = 1 (LBMc = WT). See the",
        "vignette Assumptions section. Elimination is a rate constant, so",
        "clearance CL = Kel * V inherits the same linear weight scaling.",
        "Cohort mean 90 kg, SD 13 kg (Table 1, printed as '90 (77-103)')."
      ),
      source_name = "Bodyweight"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "clozapine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "clozapine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 15,
    n_studies = 1,
    n_observations = 121,
    age_range = "18-55 years (inclusion); mean 44, SD 11 years",
    weight_range = "mean 90 kg, SD 13 kg",
    height_range = "mean 1.81 m, SD 0.065 m",
    sex_female_pct = 20,
    race_ethnicity = "Western European descent only (inclusion criterion, to minimise CYP genetic variation)",
    disease_state = "Schizophrenia, outpatients on a stable clozapine dose for at least 2 weeks",
    dose_range = "Clozapine 75-800 mg/day, mean 287 mg/day; once-daily users split the dose into two halves (evening before and morning of sampling)",
    regions = "Netherlands (Mental Health Organization Noord-Holland Noord)",
    renal_function = "Serum creatinine mean 80 umol/L, SD 15 umol/L; renal elimination fixed at zero",
    co_medication = "Smokers 10 of 15; plasma-level-affecting co-medication: fluvoxamine 3, omeprazole 3, insulin 1, sodium valproate 1 (none modelled as covariates)",
    notes = paste(
      "Baseline characteristics from Table 1. Table 1 labels every",
      "continuous row 'mean +/- SD' and prints the bounds of the mean +/- SD",
      "interval in parentheses (e.g. age '44 (33-55)'); the dose row",
      "'287 (75-800)' is not symmetric about its mean and is read as the",
      "range. 121 concentrations: 49 venous plasma and 72 capillary",
      "dried-blood-spot (DBS) samples, drawn pre-dose and 2, 4, 6 and 8 h",
      "after a morning dose (Table 2). DBS whole-blood results were converted",
      "to plasma-equivalent concentrations with a conversion factor fixed",
      "from the authors' earlier clinical validation study of the same",
      "patients, so the model predicts PLASMA concentration for both sample",
      "types. Patients with infection or inflammation on the sampling day",
      "were excluded."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Table 3, column 'Clozapine model Final (mean +/- SD)', the selected
    # one-compartment model without lag time (AIC -119 vs 267 for the
    # literature default model, Jerling 1994, ref 21). The values are
    # population means with the interindividual SD on the natural scale of
    # the parameter, the MwPharm KinPop iterative two-stage Bayesian
    # output. Parameters were 'assumed to be log-normally distributed'
    # (Methods), so the mean is carried as the typical value and the SD is
    # converted to a log-scale variance by omega^2 = log(1 + (SD/mean)^2).
    # Same convention as the MwPharm-built Daskapan_2019_darunavir model.
    # ---------------------------------------------------------------------

    # Absorption. Held at the default-model (literature) value because 'a
    # Bayesian-estimated value did not improve the model, probably due to
    # lack of data in the absorption phase' (Results).
    lka <- fixed(log(1.37))
    label("Absorption rate constant (1/h; literature prior, not estimated)") # Table 3, Ka = 1.37 +/- 0.68 (fixed), identical to the default-model column

    # Metabolic elimination rate constant. The published elimination model is
    # Kel = Kelm + Kelr * CLcr with Kelr 'set to zero' because renal
    # elimination of unmetabolised clozapine is negligible (Methods), so
    # Kel = Kelm and no creatinine-clearance column enters the model.
    lkel <- log(0.0641)
    label("Elimination rate constant Kelm (1/h)") # Table 3, Kelm = 0.0641 +/- 0.0365 (bootstrap median 0.0635)

    # Volume of distribution per kg body weight, apparent (F fixed at 1, so
    # Vd is Vd/F). Paper-named per-kg coefficient, so 'vd' rather than
    # 'vc'; the absolute volume vc = vd * WT is derived in model().
    lvd <- log(5.21)
    label("Apparent volume of distribution per kg body weight Vd/F (L/kg)") # Table 3, Vd/F = 5.21 +/- 3.82 L/kg (bootstrap median 5.20); Methods 'expressed in L/kg bodyweight'

    # Bioavailability. 'The bioavailabity (F) of clozapine was fixed at 1, so
    # volume of distribution (Vd) is expressed as Vd/F' (Methods).
    lfdepot <- fixed(log(1))
    label("Oral bioavailability (fraction; F = 1, so V is apparent V/F)") # Methods; Table 3 footnote 'F was fixed at 1'

    # IIV, log-scale variances from the Table 3 natural-scale SDs:
    #   Kelm 0.0365 / 0.0641 -> CV 56.9 % -> log(1 + 0.569423^2) = 0.280840
    #   Vd   3.82 / 5.21     -> CV 73.3 % -> log(1 + 0.733205^2) = 0.430216
    #   Ka   0.68 / 1.37     -> CV 49.6 % -> log(1 + 0.496350^2) = 0.220230
    # No covariance is reported, so the block is diagonal. The Ka SD is the
    # literature prior SD the paper kept with the fixed Ka; it is carried as
    # a fixed variance so the full MwPharm prior is preserved.
    etalka ~ fixed(0.220230) # Table 3, Ka 1.37 +/- 0.68, held at the literature prior -> 49.6 % CV
    etalkel ~ 0.280840 # Table 3, Kelm 0.0641 +/- 0.0365 -> 56.9 % CV; shrinkage 0.079
    etalvd ~ 0.430216 # Table 3, Vd/F 5.21 +/- 3.82 -> 73.3 % CV; shrinkage 0.054

    # Residual error. 'The assay error was concentration dependent and set to
    # SD = 10 + 0.1 x C where C is the observed or calculated (DBS) clozapine
    # plasma concentration in ug/L' (Methods). A linear sum of additive and
    # proportional terms is nlmixr2's combined1() form; both are assay
    # settings, not estimates, so both are fixed.
    addSd <- fixed(10)
    label("Additive residual SD (ug/L; assay SD intercept, not estimated)") # Methods 'SD = 10 + 0.1 x C'
    propSd <- fixed(0.1)
    label("Proportional residual SD (fraction; assay SD slope, not estimated)") # Methods 'SD = 10 + 0.1 x C'
  })

  model({
    # Individual parameters
    ka <- exp(lka + etalka)
    kel <- exp(lkel + etalkel)
    vd <- exp(lvd + etalvd)

    # Absolute apparent volume (L): per-kg Vd times total body weight
    # (Methods, 'expressed in L/kg bodyweight').
    vc <- vd * WT
    cl <- kel * vc

    # One compartment with first-order absorption and elimination; lag time,
    # a peripheral compartment and fat-tissue distribution were each tested
    # and rejected on AIC (Results, Population PK model).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    f(depot) <- exp(lfdepot)

    # Dose in mg and volume in L give mg/L; x 1000 gives ug/L, the unit of the
    # paper's assay-error model and figures.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
