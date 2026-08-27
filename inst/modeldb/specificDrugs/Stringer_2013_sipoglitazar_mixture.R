Stringer_2013_sipoglitazar_mixture <- function() {
  description <- "Two-compartment population PK model for sipoglitazar with apparent clearance assigned to one of three latent metabolizer subpopulations by a NONMEM mixture model, fitted without using UGT2B15 genotype"
  reference <- paste(
    "Stringer F, Ploeger BA, DeJongh J, Scott G, Urquhart R, Karim A,",
    "Danhof M. Evaluation of the Impact of UGT Polymorphism on the",
    "Pharmacokinetics and Pharmacodynamics of the Novel PPAR Agonist",
    "Sipoglitazar. J Clin Pharmacol. 2013;53(3):256-263.",
    "doi:10.1177/0091270012447121. Mixture-model parameters from the",
    "supplementary online material, Table A1.",
    sep = " "
  )
  vignette <- "Stringer_2013_sipoglitazar"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # As in the companion genotype model, residual variability is stratified by
  # study population, so the canonical propSd consumed by the error model is
  # derived inside model() from two study-specific ini() magnitudes.
  paper_specific_residual_sds <- c("propSdPhase1", "propSdPhase2")

  compartmentData <- list(
    depot       = list(analyte = "sipoglitazar", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "sipoglitazar", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "sipoglitazar", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    UGT2B15_IM = list(
      description        = "Subject assigned to the UGT2B15 intermediate-metaboliser subpopulation (POP2)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 with UGT2B15_PM = 0, i.e. the extensive-metaboliser subpopulation (POP1)",
      notes              = paste(
        "Latent subpopulation membership. Stringer 2013 estimated this with the",
        "NONMEM $MIX subroutine from apparent clearance alone, deliberately without",
        "using the UGT2B15 genotype. rxode2 / nlmixr2 has no $MIX analogue, so the",
        "assignment enters here as a covariate indicator; the mixture prior weights",
        "the paper estimated are recorded in population$mixture_probabilities and",
        "are what a simulation should sample the indicator from."
      ),
      source_name        = "POP2 (CL2 IM)"
    ),
    UGT2B15_PM = list(
      description        = "Subject assigned to the UGT2B15 poor-metaboliser subpopulation (POP3)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 with UGT2B15_IM = 0, i.e. the extensive-metaboliser subpopulation (POP1)",
      notes              = paste(
        "Latent subpopulation membership; see UGT2B15_IM notes. Stringer 2013",
        "reports that 61 of the 744 subjects genotyped *1/*1 or *1/*2 (8%) were",
        "assigned to this subpopulation, i.e. their phenotype did not match their",
        "genotype (supplementary Results, Figure S3)."
      ),
      source_name        = "POP3 (CL3 PM)"
    ),
    FFM = list(
      description        = "Fat-free mass",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Linear centered effect on the central volume per Stringer 2013 equation 4.",
        "The paper does not report the cohort median FFM, so the centering value is",
        "an assumption of this extraction (50 kg); see the vignette Errata."
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
        "inter-individual variability on that residual."
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
    mixture_probabilities = c(EM = 0.18, IM = 0.522, PM = 0.30),
    notes          = paste(
      "Same pooled analysis population as Stringer_2013_sipoglitazar; the two models",
      "differ only in how a subject's clearance stratum is assigned. The mixture",
      "prior weights are supplementary Table A1 rows 'Probability of belonging to",
      "POP 1/2/3' (0.18 with RSE 18.6%, 0.522, 0.30). They are metadata rather than",
      "ini() parameters because the subpopulation enters this model as a covariate",
      "indicator. Table A1 also reports PROB = 0.367 with the footnote",
      "POP2 = (1 - POP1) * PROB and POP3 = (1 - POP1) * (1 - PROB); evaluating those",
      "expressions gives 0.301 and 0.519, i.e. the footnote's two assignments are",
      "swapped relative to the tabulated 0.522 / 0.30. The tabulated values are used",
      "here because they are the ones that track the observed genotype frequencies",
      "(22 / 51 / 27%); see the vignette Errata.",
      "Note that the companion PD paper Stringer 2014 (J Clin Pharmacol",
      "54(4):453-461, doi:10.1002/jcph.227) fixes its per-genotype clearances at",
      "5.04 / 3.35 / 1.53 L/h, i.e. it uses THIS model's latent-subpopulation",
      "estimates and labels them by UGT2B15 genotype rather than using the",
      "genotype-stratified 4.46 / 3.25 / 1.53 L/h of Stringer 2013 Table 2."
    )
  )

  ini({
    # Apparent clearance per latent metabolizer subpopulation. As in the
    # companion genotype model every stratum carries a suffix; CL and V are
    # apparent (CL/F, V/F) because bioavailability is unknown (Table A1
    # footnote a).
    lcl_em <- log(5.04); label("Apparent clearance, extensive-metaboliser subpopulation POP1 (CL/F, L/h)")     # Suppl. Table A1: CL POP1 (EM) = 5.04 L/h (RSE 3.85%)
    lcl_im <- log(3.35); label("Apparent clearance, intermediate-metaboliser subpopulation POP2 (CL/F, L/h)")  # Suppl. Table A1: CL POP2 (IM) = 3.35 L/h (RSE 2.38%)
    lcl_pm <- log(1.53); label("Apparent clearance, poor-metaboliser subpopulation POP3 (CL/F, L/h)")          # Suppl. Table A1: CL POP3 (PM) = 1.53 L/h (RSE 2.64%)

    lvc <- log(9.06);   label("Apparent central volume of distribution (V/F, L)")     # Suppl. Table A1: V = 9.06 L (RSE 2.41%)
    lvp <- log(1.703);  label("Apparent peripheral volume of distribution (V2/F, L)") # Suppl. Table A1: V2 = 0.188 (RSE 4.93%) as a fraction of the central volume (footnote b); 0.188 * 9.06 L = 1.703 L
    lq  <- log(0.311);  label("Intercompartmental clearance (Q/F, L/h)")              # Suppl. Table A1: Q = 0.311 L/h (RSE 6.72%)

    # Absorption: zero-order release of duration D1 into the depot followed by
    # first-order transfer at ka. See the companion genotype model and the
    # vignette Errata for why the reported parameter set fixes this structure.
    lka <- log(2.15);  label("First-order absorption rate constant from the depot (1/h)") # Suppl. Table A1: ka1 = 2.15 1/h (RSE 6.19%)
    ld1 <- log(0.637); label("Duration of the zero-order input into the depot (h)")       # Suppl. Table A1: D1 = 0.637 h (RSE 3.69%)

    e_ffm_vc <- 0.00556; label("Linear effect of (FFM - 50 kg) on the central volume (1/kg)") # Suppl. Table A1: FFM on central volume of distribution = 0.00556 (RSE 16.5%)

    # IIV, log-normal. omega^2 = log(CV^2 + 1) from the Table A1
    # "IIV (%, CV%)" column.
    etalcl ~ 0.136258  # Suppl. Table A1: IIV on CL = 38.21% CV (RSE 10.7%); log(1 + 0.3821^2)
    etalvc ~ 0.112447  # Suppl. Table A1: IIV on V  = 34.50% CV (RSE 13.4%); log(1 + 0.3450^2)
    etald1 ~ 0.467514  # Suppl. Table A1: IIV on D1 = 77.20% CV (RSE 16.9%); log(1 + 0.7720^2)

    propSdPhase1 <- 0.234243; label("Proportional residual SD, phase I healthy volunteers (fraction)") # Suppl. Table A1: sigma^2 phase I = 0.05487 (RSE 9.23%); sqrt(0.05487)
    propSdPhase2 <- 0.408656; label("Proportional residual SD, phase II patients (fraction)")          # Suppl. Table A1: sigma^2 phase II = 0.167 (RSE 10.2%); sqrt(0.167)
    etapropSdPhase2 ~ 0.419432 # Suppl. Table A1: IIV on the phase II residual = 72.18% CV (RSE 15.6%); log(1 + 0.7218^2)
  })

  model({
    # 1. Derived covariate terms. Both indicators 0 selects the
    #    extensive-metaboliser subpopulation POP1.
    lclPop <- lcl_em * (1 - UGT2B15_IM - UGT2B15_PM) +
      lcl_im * UGT2B15_IM +
      lcl_pm * UGT2B15_PM
    ffmVc <- 1 + e_ffm_vc * (FFM - 50)

    # 2. Individual parameters.
    cl <- exp(lclPop + etalcl)
    vc <- exp(lvc + etalvc) * ffmVc
    #    Table A1 footnote b: the peripheral volume was implemented as a
    #    fraction of the central compartment. exp(lvp - lvc) = 0.188.
    vp <- exp(lvp - lvc) * vc
    q  <- exp(lq)
    ka <- exp(lka)
    d1 <- exp(ld1 + etald1)

    # 3. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system.
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
