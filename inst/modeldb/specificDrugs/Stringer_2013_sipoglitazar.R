Stringer_2013_sipoglitazar <- function() {
  description <- "Two-compartment population PK model for sipoglitazar with UGT2B15 genotype-stratified apparent clearance in healthy volunteers and adults with type 2 diabetes"
  reference <- paste(
    "Stringer F, Ploeger BA, DeJongh J, Scott G, Urquhart R, Karim A,",
    "Danhof M. Evaluation of the Impact of UGT Polymorphism on the",
    "Pharmacokinetics and Pharmacodynamics of the Novel PPAR Agonist",
    "Sipoglitazar. J Clin Pharmacol. 2013;53(3):256-263.",
    "doi:10.1177/0091270012447121.",
    sep = " "
  )
  vignette <- "Stringer_2013_sipoglitazar"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The residual variability is stratified by study population (phase I
  # healthy volunteers vs phase II patients), so the canonical propSd
  # consumed by the error model is derived inside model() from two
  # study-specific ini() magnitudes. Same construction as
  # Xie_2025_aztreonam_avibactam and Cammarata_2024_sulbactam_durlobactam.
  paper_specific_residual_sds <- c("propSdPhase1", "propSdPhase2")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Doses are oral tablets in mg; plasma concentrations
  # were reported in ng/mL (Stringer 2013 Bioanalysis) and are carried here
  # in mg/L, which is the same number divided by 1000.
  compartmentData <- list(
    depot       = list(analyte = "sipoglitazar", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "sipoglitazar", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "sipoglitazar", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    UGT2B15_STAR2_HET = list(
      description        = "UGT2B15 *1/*2 heterozygote at the D85Y (rs1902023) coding polymorphism",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 with UGT2B15_STAR2_HOM = 0, i.e. the *1/*1 wild-type homozygote stratum",
      notes              = paste(
        "Germline genotype, time-fixed. Paired with UGT2B15_STAR2_HOM to select one",
        "of three stratum-specific typical clearances. Genotype frequencies in the",
        "pooled analysis population (Stringer 2013 Table 1, n = 1151): *1/*1 22%,",
        "*1/*2 51%, *2/*2 27%."
      ),
      source_name        = "UGT2B15*1/*2"
    ),
    UGT2B15_STAR2_HOM = list(
      description        = "UGT2B15 *2/*2 homozygous-variant (poor metaboliser) genotype at D85Y (rs1902023)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 with UGT2B15_STAR2_HET = 0, i.e. the *1/*1 wild-type homozygote stratum",
      notes              = paste(
        "Germline genotype, time-fixed. Paired with UGT2B15_STAR2_HET. The *2/*2",
        "stratum has the lowest apparent clearance (1.53 vs 4.46 L/h in *1/*1)."
      ),
      source_name        = "UGT2B15*2/*2"
    ),
    FFM = list(
      description        = "Fat-free mass",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Linear centered effect on the central volume per Stringer 2013 equation 4,",
        "P(mean) = theta(i) * (1 + theta(f) * (COV - COV(median))). Stringer 2013",
        "cites Anderson & Holford (Annu Rev Pharmacol Toxicol 2008;48:303-332) for",
        "the FFM derivation but does not report the cohort median FFM, so the",
        "centering value is an assumption of this extraction (50 kg, the register's",
        "reference-adult fat-free mass). A user simulating a specific cohort should",
        "substitute that cohort's median FFM; the effect is centered, so the choice",
        "only rescales where the typical V = 9.03 L is anchored."
      ),
      source_name        = "FFM"
    ),
    STUDY_SIPO_PHASE2 = list(
      description        = "Observation belongs to the phase II type-2-diabetes patient studies (EC201, EC202)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (phase I healthy-volunteer study 006)",
      notes              = paste(
        "Selects the phase II proportional residual magnitude and switches on the",
        "inter-individual variability Stringer 2013 places on that residual. The",
        "phase II samples were all intended as troughs but the actual time after",
        "dose was never recorded (supplementary methods), so the phase II residual",
        "absorbs an unrecorded-sampling-time component in addition to assay error."
      ),
      source_name        = "study"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1151,
    n_studies      = 3,
    age_range      = "18-75 years",
    age_median     = "26 years (phase I, study 006); 56 years (EC201); 57 years (EC202)",
    weight_range   = "48-149 kg",
    weight_median  = "71 kg (phase I, study 006); 90 kg (EC201); 87 kg (EC202)",
    sex_female_pct = 54.1,
    disease_state  = paste(
      "Pooled: 524 healthy volunteers (phase I study 006) and 627 adults with type",
      "2 diabetes mellitus and no prior exposure to antidiabetic medication",
      "(phase II studies EC201 and EC202)"
    ),
    dose_range     = paste(
      "Single oral 64 mg (study 006); once- or twice-daily oral 8, 16, 32 and 64 mg",
      "total daily dose for 12 weeks (EC201, EC202)"
    ),
    genotype_frequency = c(`*1/*1` = 22, `*1/*2` = 51, `*2/*2` = 27),
    notes          = paste(
      "Demographics from Stringer 2013 Table 1 (median and range). Female percentage",
      "is the pooled (304 + 155 + 164) / 1151. The phase I study enrolled a range of",
      "ethnicities but Table 1 does not tabulate race, so race_ethnicity is omitted.",
      "Sampling was dense in the phase I study (0, 1, 2, 3, 4, 6, 8, 12, 16 and 24 h",
      "postdose) and sparse in the phase II studies (3 trough samples at weeks 4, 6",
      "and 8). A predetermined chronic-exposure safety limit of AUC > 73 mg*h/L is",
      "quoted in Methods and drawn on Figure 2."
    )
  )

  ini({
    # Apparent clearance is stratified by UGT2B15 genotype: Stringer 2013
    # optimised a separate median clearance for each of the three genotype
    # groups rather than a reference value plus offsets, so every stratum
    # carries a suffix (parameter-names.md, "Stratum-suffixed parameters").
    # CL and V are apparent (CL/F, V/F) - bioavailability is unknown
    # (Table 2 footnote a).
    lcl_s1s1 <- log(4.46); label("Apparent clearance, UGT2B15 *1/*1 (CL/F, L/h)")   # Table 2: CL (*1/*1) = 4.46 L/h (RSE 2.5%)
    lcl_s1s2 <- log(3.25); label("Apparent clearance, UGT2B15 *1/*2 (CL/F, L/h)")   # Table 2: CL (*1/*2) = 3.25 L/h (RSE 2.2%)
    lcl_s2s2 <- log(1.53); label("Apparent clearance, UGT2B15 *2/*2 (CL/F, L/h)")   # Table 2: CL (*2/*2) = 1.53 L/h (RSE 2.2%)

    lvc <- log(9.03);  label("Apparent central volume of distribution (V/F, L)")    # Table 2: V = 9.03 L (RSE 2.4%)
    lvp <- log(1.71);  label("Apparent peripheral volume of distribution (V2/F, L)") # Table 2: V2 = 1.71 L (RSE 4.9%); Table 2 footnote b - implemented as a fraction of the central volume
    lq  <- log(0.313); label("Intercompartmental clearance (Q/F, L/h)")             # Table 2: Q = 0.313 L/h (RSE 6.6%)

    # Absorption. Table 2 reports only ka and D1 with no fraction splitting
    # the dose between the two routes, and the NONMEM duration variable is
    # named D1 (the depot slot). The parameter set is therefore complete only
    # for a zero-order release of duration D1 into the depot followed by
    # first-order transfer at ka; see the vignette Errata for the check
    # against Supplementary Figure S2a.
    lka <- log(2.07);  label("First-order absorption rate constant from the depot (1/h)") # Table 2: ka = 2.07 1/h (RSE 4.8%)
    ld1 <- log(0.568); label("Duration of the zero-order input into the depot (h)")       # Table 2: D1 = 0.568 h (RSE 6.81%)

    # Fat-free mass on the central volume, linear and centered on the cohort
    # median per Stringer 2013 equation 4. The paper does not report the
    # median FFM; see covariateData$FFM$notes and the vignette Errata.
    e_ffm_vc <- 0.00349; label("Linear effect of (FFM - 50 kg) on the central volume (1/kg)") # Table 2: theta(f) = 0.00349 (RSE 27.2%)

    # IIV, log-normal (equation 1). omega^2 = log(CV^2 + 1) from the CV%
    # reported in the Table 2 "IIV (%, CV%)" column.
    etalcl ~ 0.150086  # Table 2: IIV on CL = 40.25% CV (RSE 7.72%);  log(1 + 0.4025^2)
    etalvc ~ 0.110665  # Table 2: IIV on V  = 34.21% CV (RSE 13.0%);  log(1 + 0.3421^2)
    etald1 ~ 0.478053  # Table 2: IIV on D1 = 78.29% CV (RSE 14.8%);  log(1 + 0.7829^2)

    # Residual error: proportional, with a separate magnitude for the phase I
    # healthy-volunteer study and the phase II patient studies (equations 2
    # and 3). The phase II residual additionally carries IIV.
    propSdPhase1 <- 0.234947; label("Proportional residual SD, phase I healthy volunteers (fraction)") # Table 2: sigma^2 phase I = 0.0552 (RSE 8.8%); sqrt(0.0552)
    propSdPhase2 <- 0.408656; label("Proportional residual SD, phase II patients (fraction)")          # Table 2: sigma^2 phase II = 0.167 (RSE 10.2%); sqrt(0.167)
    etapropSdPhase2 ~ 0.464430 # Table 2: IIV on the phase II residual = 76.88% CV (RSE 14.9%); log(1 + 0.7688^2)
  })

  model({
    # 1. Derived covariate terms.
    #    UGT2B15 genotype selects one of the three stratum-specific typical
    #    clearances; both indicators 0 selects the *1/*1 reference stratum.
    lclGeno <- lcl_s1s1 * (1 - UGT2B15_STAR2_HET - UGT2B15_STAR2_HOM) +
      lcl_s1s2 * UGT2B15_STAR2_HET +
      lcl_s2s2 * UGT2B15_STAR2_HOM
    #    Equation 4: P(mean) = theta(i) * (1 + theta(f) * (COV - COV(median))).
    #    50 kg is the assumed cohort-median FFM (not reported by the paper).
    ffmVc <- 1 + e_ffm_vc * (FFM - 50)

    # 2. Individual parameters.
    cl <- exp(lclGeno + etalcl)
    vc <- exp(lvc + etalvc) * ffmVc
    #    Table 2 footnote b: the peripheral volume was implemented as a
    #    fraction of the central compartment, so it inherits the central
    #    volume's IIV and FFM effect. exp(lvp - lvc) = 1.71 / 9.03 = 0.189.
    vp <- exp(lvp - lvc) * vc
    q  <- exp(lq)
    ka <- exp(lka)
    d1 <- exp(ld1 + etald1)

    # 3. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system. The dose fills the depot at a zero-order rate over
    #    [0, d1] and the depot then drains into the central compartment at ka.
    d/dt(depot) <- -ka * depot
    dur(depot) <- d1
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 5. Observation and error.
    Cc <- central / vc
    propSd <- propSdPhase1 * (1 - STUDY_SIPO_PHASE2) +
      propSdPhase2 * exp(etapropSdPhase2) * STUDY_SIPO_PHASE2
    Cc ~ prop(propSd)
  })
}
