Eechoute_2012_imatinib <- function() {
  description <- paste0(
    "Two-compartment population PK model with a five-compartment transit ",
    "absorption chain for oral imatinib in adults with gastrointestinal ",
    "stromal tumor, from a long-term prospective study (Eechoute 2012). ",
    "This is the only two-compartment model and the only transit-absorption ",
    "model among the 15 evaluated by Yang 2025, and the only one whose ",
    "clearance and volumes are true (not apparent) because it carries an ",
    "explicit relative bioavailability term. Clearance falls linearly with ",
    "the volume of liver metastases. Both the absorption rate constant and ",
    "relative bioavailability decay exponentially with time since the first ",
    "dose, with a shared rate constant of 0.0256 per day (a half-life of ",
    "about 27 days), capturing the well-described fall in imatinib exposure ",
    "over the first months of continuous therapy. Residual error is ",
    "log-scale additive (equivalent to proportional in linear space). ",
    "TRANSCRIBED FROM A SECONDARY SOURCE: the parameter values come from ",
    "Table 1 of Yang 2025, an external evaluation of 15 published imatinib ",
    "population PK models, not from the primary publication. Re-extract ",
    "from Eechoute 2012 when that paper is obtained; the topology of the ",
    "transit chain in particular is an inference, flagged below."
  )
  reference <- paste0(
    "Eechoute K, Fransson MN, Reyners AK, de Jong FA, Sparreboom A, van ",
    "der Graaf WTA, Verweij J, Wiemer EAC, Steeghs N, Mathijssen RHJ, ",
    "Friberg LE. A long-term prospective population pharmacokinetic study ",
    "on imatinib plasma concentrations in GIST patients. Clin Cancer Res. ",
    "2012;18(20):5780-5787. doi:10.1158/1078-0432.CCR-12-0490. ",
    "PARAMETER SOURCE (secondary): Yang T, Rasmussen ASB, Weimann A, ",
    "Thastrup M, Rank CU, Als-Nielsen B, Malmros J, Wik HS, Lohi O, ",
    "Overgaard U, Johannsdottir IMR, Vaitkeviciene G, Dalhoff K, ",
    "Schmiegelow K, Lund TM. Published population pharmacokinetic models ",
    "of imatinib perform poorly on TDM data from pediatric patients. ",
    "Target Oncol. 2025;20(5):871-886. Table 1, row 'Eechoute et al. ",
    "(2012)'. doi:10.1007/s11523-025-01172-2."
  )
  vignette <- "Yang_2025_imatinib_external_evaluation"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot       = list(analyte = "imatinib", units = "mg", specimen = "administration site", verified = FALSE),
    transit1    = list(analyte = "imatinib", units = "mg", specimen = "administration site", verified = FALSE),
    transit2    = list(analyte = "imatinib", units = "mg", specimen = "administration site", verified = FALSE),
    transit3    = list(analyte = "imatinib", units = "mg", specimen = "administration site", verified = FALSE),
    transit4    = list(analyte = "imatinib", units = "mg", specimen = "administration site", verified = FALSE),
    transit5    = list(analyte = "imatinib", units = "mg", specimen = "administration site", verified = FALSE),
    central     = list(analyte = "imatinib", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "imatinib", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    TUM_VOL = list(
      description        = "Volume of liver metastases",
      units              = "mm^3",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "UNIT CONVERSION IS REQUIRED AT THE CALL SITE. Yang 2025 Table 1 ",
        "abbreviation list defines the source column as 'LIVM volume of ",
        "liver metastasis (cm^3)' and prints the effect as ",
        "CL = 9.12 x (1 - 0.000381 x LIVM) with LIVM in cm^3. The ",
        "canonical TUM_VOL column is carried in mm^3 (see ",
        "inst/references/covariate-columns.md, which instructs models to ",
        "convert on ingestion and rescale the per-model reference so the ",
        "covariate term stays numerically invariant), and 1 cm^3 = 1000 ",
        "mm^3. The model therefore writes the term as ",
        "(1 - e_tumvol_cl * TUM_VOL / 1000) so that the PUBLISHED ",
        "coefficient 0.000381 remains visible verbatim in ini() and the ",
        "conversion is explicit rather than folded into a rescaled ",
        "constant. A data assembler holding cm^3 must multiply by 1000 ",
        "before using this model.",
        " ",
        "SCOPE: TUM_VOL is the register's canonical for any volumetric ",
        "tumor-burden measurement (scope promoted 2026-07-30 to pool ",
        "measurement modalities). Eechoute 2012 measured specifically the ",
        "hepatic metastatic burden, which in a GIST cohort is the dominant ",
        "measurable tumor volume; the specificity is recorded here rather ",
        "than by registering a separate liver-only canonical. The effect ",
        "is LINEAR, not power-form, and is negative: clearance falls with ",
        "increasing hepatic tumor burden, consistent with reduced ",
        "functional hepatic mass. Note that the linear form has a root -- ",
        "the term reaches zero at LIVM = 2625 cm^3 and turns NEGATIVE ",
        "beyond it -- so a simulation must not extrapolate past that ",
        "burden. Set TUM_VOL = 0 to recover the metastasis-free typical ",
        "clearance of 9.12 L/h."
      ),
      source_name        = "LIVM"
    ),
    T_FIRSTDOSE = list(
      description        = "Time elapsed since the first imatinib dose of the treatment course",
      units              = "h",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Yang 2025 Table 1 abbreviation list: 'TAF time from first dose ",
        "administration (h)'. This is a treatment-duration clock that runs ",
        "monotonically across the whole record and does NOT reset at each ",
        "dose; it must not be confused with time-after-dose. It drives two ",
        "terms that share a single rate constant of 0.0256 per day: ",
        "ka = 0.699 x (1 + 1.18 x exp(-0.0256 x TAF/24)) and ",
        "F = 1 + 0.482 x exp(-0.0256 x TAF/24). The division by 24 is ",
        "printed inside the exponent in the source and converts the ",
        "canonical hours to days so the published per-day rate constant is ",
        "used unchanged. Both quantities START HIGH and DECAY to an ",
        "asymptote with a half-life of ln(2)/0.0256 = 27 days: ka falls ",
        "from 1.524 to 0.699 per hour and relative bioavailability falls ",
        "from 1.482 to 1. The combined effect is a roughly one-third fall ",
        "in imatinib exposure over the first months of continuous therapy, ",
        "which is the phenomenon this long-term prospective study was ",
        "designed to characterise. For a simulation whose first dose is at ",
        "time zero, T_FIRSTDOSE equals the model time."
      ),
      source_name        = "TAF"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 50L,
    n_studies      = 1L,
    n_observations = "1743 imatinib plasma concentrations (Yang 2025 Table 1)",
    age_range      = "39-82 years",
    disease_state  = "Adults with gastrointestinal stromal tumor (GIST)",
    dose_range     = "Oral imatinib 400-800 mg total daily dose",
    regions        = "Europe",
    bioanalytical  = "LC-MS/MS; limit of quantification not reported in Yang 2025 Table 1",
    notes          = paste0(
      "The longest-followed cohort among the 15 models evaluated by Yang ",
      "2025 (1743 samples from 50 patients), which is what makes the ",
      "time-since-first-dose effect estimable. Demographic detail beyond ",
      "the row above (weight range, sex split, race) is not reported by ",
      "the secondary source and must be read from the primary publication."
    )
  )

  ini({
    # ----- Disposition (Yang 2025 Table 1, Eechoute row) -----
    # These are TRUE clearances and volumes, not apparent ones: this is the
    # only model among the 15 that carries an explicit bioavailability
    # term, so F is not folded into CL and V. Yang 2025 Table 1 writes them
    # as 'CL', 'Vc', 'Q', 'Vp' (no '/F'), in contrast to the 'CL/F', 'Vc/F'
    # of every other row.
    lcl <- log(9.12); label("Clearance CL in a patient with no liver metastases (L/h)")  # Yang 2025 Table 1: CL = 9.12 x (1 - 0.000381 x LIVM)
    lvc <- log(128); label("Central volume of distribution Vc (L)")  # Yang 2025 Table 1: Vc = 128
    lq <- log(24.9); label("Intercompartmental clearance Q (L/h)")  # Yang 2025 Table 1: Q = 24.9
    lvp <- log(197); label("Peripheral volume of distribution Vp (L)")  # Yang 2025 Table 1: Vp = 197

    # ----- Liver-metastasis effect on clearance (LINEAR, not power) -----
    # Coefficient is per cm^3 as published; the model divides the canonical
    # mm^3 column by 1000 at the call site.
    e_tumvol_cl <- 0.000381; label("Linear coefficient of liver-metastasis volume on CL, per cm^3 (unitless per cm^3)")  # Yang 2025 Table 1: (1 - 0.000381 x LIVM)

    # ----- Absorption -----
    # Five-compartment transit chain plus a first-order step into the
    # central compartment. Yang 2025 Table 1 gives the structural model as
    # '2-CMT, T5' and its abbreviation list defines 'T5 5 transit
    # compartment absorption model' and 'Ktr transit rate constant between
    # the transit compartments (1/h)'.
    lktr <- log(15.8); label("Transit rate constant between transit compartments Ktr (1/h)")  # Yang 2025 Table 1: Ktr = 15.8
    lka <- log(0.699); label("Asymptotic first-order absorption rate constant ka at long treatment duration (1/h)")  # Yang 2025 Table 1: Ka = 0.699 x (1 + 1.18 x e^(-0.0256 x TAF/24))

    # ----- Time-since-first-dose decay on ka and on bioavailability -----
    # Both terms share the SAME rate constant (0.0256 per day) but have
    # different amplitudes, so they are carried as three parameters.
    e_tfirstdose_ka <- 1.18; label("Amplitude of the time-since-first-dose decay on ka (unitless)")  # Yang 2025 Table 1: (1 + 1.18 x e^(...))
    e_tfirstdose_fdepot <- 0.482; label("Amplitude of the time-since-first-dose decay on relative bioavailability (unitless)")  # Yang 2025 Table 1: F = 1 + 0.482 x e^(...)
    k_tfirstdose <- 0.0256; label("Rate constant of the time-since-first-dose decay, shared by ka and F (1/day)")  # Yang 2025 Table 1: e^(-0.0256 x TAF/24), shared by the Ka and F terms

    # ----- Inter-individual variability -----
    # Yang 2025 Table 1 reports 'CV%(CL): 49.5%', 'CV%(Vc-Ka): 70.9%',
    # 'CV%(Ktr): 160%' and 'CV%(Vp): 65.9%'. The tabulated CV% is taken as
    # omega (the log-scale standard deviation), so the variance is
    # (CV/100)^2 -- the convention used throughout this transcription.
    #
    # The 'Vc-Ka' entry denotes a SINGLE random effect SHARED by the
    # central volume and the absorption rate constant: the secondary source
    # prints one omega under a hyphenated two-parameter label, in contrast
    # to the separate per-parameter entries on the same row. It is
    # implemented below by multiplying both vc and ka by exp(etalvc). This
    # is the reading that the notation supports; if the primary instead
    # reports two separate omegas of equal magnitude, or a full block, the
    # correction is confined to this one eta.
    etalcl ~ 0.245025  # Yang 2025 Table 1: CV%(CL) 49.5% -> omega^2 = 0.495^2
    etalvc ~ 0.502681  # Yang 2025 Table 1: CV%(Vc-Ka) 70.9% -> omega^2 = 0.709^2; SHARED by vc and ka
    etalktr ~ 2.56  # Yang 2025 Table 1: CV%(Ktr) 160% -> omega^2 = 1.60^2
    etalvp ~ 0.434281  # Yang 2025 Table 1: CV%(Vp) 65.9% -> omega^2 = 0.659^2

    # ----- Residual unexplained variability -----
    # Yang 2025 Table 1 reports 'Log-scale add: 35%'. An additive residual
    # on the log-transformed observation is equivalent to a proportional
    # residual in linear space, and the secondary source reports the
    # magnitude as a PERCENT, i.e. already on the coefficient-of-variation
    # scale, so it is encoded as prop().
    propSd <- 0.35; label("Proportional residual error (fraction; log-scale additive in the source)")  # Yang 2025 Table 1: log-scale add 35%
  })

  model({
    # ----- 1. Time-since-first-dose decay factor -----
    # T_FIRSTDOSE is carried in hours; the published exponent divides by 24
    # so that the rate constant keeps its per-day units.
    tfd_decay <- exp(-k_tfirstdose * T_FIRSTDOSE / 24)

    # ----- 2. Individual parameters -----
    # The liver-metastasis term is linear and the canonical TUM_VOL column
    # is mm^3, so it is divided by 1000 to reach the cm^3 in which the
    # published coefficient is expressed.
    cl <- exp(lcl + etalcl) * (1 - e_tumvol_cl * TUM_VOL / 1000)

    # etalvc is SHARED between vc and ka (see the ini() comment).
    vc <- exp(lvc + etalvc)
    ka <- exp(lka + etalvc) * (1 + e_tfirstdose_ka * tfd_decay)

    q <- exp(lq)
    vp <- exp(lvp + etalvp)
    ktr <- exp(lktr + etalktr)

    # Relative bioavailability decays from 1.482 at the first dose towards
    # 1 at long treatment duration.
    fdepot <- 1 + e_tfirstdose_fdepot * tfd_decay

    # ----- 3. ODE system -----
    # Absorption: the dose enters `depot` and passes through five transit
    # compartments at rate ktr before entering the central compartment at
    # rate ka.
    #
    # STRUCTURAL INFERENCE: Yang 2025 Table 1 gives the absorption model
    # only as 'T5' with the two rate constants Ktr and Ka, and does not
    # print the chain topology. The layout used here -- a dosing depot
    # followed by five named transit compartments, all inter-compartment
    # steps at ktr, and the final step into central at ka -- is the
    # standard Savic-type chain and is the reading that yields exactly five
    # transit compartments. An alternative reading (dosing directly into
    # transit1, so that there is no separate depot) would shorten the mean
    # absorption time by one 1/ktr step = 0.063 h, which is under 4% of the
    # total mean absorption time of about 1.75 h; the choice is therefore
    # numerically almost immaterial, but it is recorded here as an
    # inference to confirm against the primary.
    d/dt(depot)    <- -ktr * depot
    d/dt(transit1) <-  ktr * depot    - ktr * transit1
    d/dt(transit2) <-  ktr * transit1 - ktr * transit2
    d/dt(transit3) <-  ktr * transit2 - ktr * transit3
    d/dt(transit4) <-  ktr * transit3 - ktr * transit4
    d/dt(transit5) <-  ktr * transit4 - ka  * transit5

    d/dt(central) <-
      ka * transit5 -
      cl * central / vc -
      q * central / vc +
      q * peripheral1 / vp
    d/dt(peripheral1) <-
      q * central / vc - q * peripheral1 / vp

    f(depot) <- fdepot

    # ----- 4. Observation and error -----
    # central is mg and vc is L, so central/vc is mg/L; x 1000 gives ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
