Li_2017_cediranib <- function() {
  description <- "Two-compartment population PK model for oral cediranib (AZD2171) in adult cancer patients (Li 2017), with sequential zero- and first-order absorption (zero-order release into depot followed by first-order absorption to central), bioavailability fixed to 1, allometric power scaling on apparent clearance ((WT/73 kg)^0.517 and (Age/59 y)^-0.409) and on apparent central volume ((WT/73 kg)^0.65), correlated inter-individual variability between CL/F and Vc/F (correlation 0.839), independent IIV on Ka, and proportional residual error (rich-sampling estimate)."
  reference <- paste(
    "Li J, Al-Huniti N, Henningsson A, Tang W, Masson E. (2017).",
    "Population pharmacokinetic and exposure simulation analysis for cediranib (AZD2171)",
    "in pooled Phase I/II studies in patients with cancer.",
    "Br J Clin Pharmacol 83(8):1723-1733. doi:10.1111/bcp.13266"
  )
  vignette <- "Li_2017_cediranib"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline in Li 2017. Allometric power scaling with reference 73 kg (study-population median; Table 1) on CL/F (exponent 0.517) and Vc/F (exponent 0.65). Study median 73 kg (range 35-150).",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline in Li 2017. Power scaling (AGE/59)^-0.409 on CL/F. Reference 59 years is the study-population median (Table 1; range 19-89).",
      source_name        = "AGE"
    )
  )

  population <- list(
    n_subjects     = 625L,
    n_studies      = 19L,
    n_observations = 7011L,
    age_range      = "19-89 years",
    age_median     = "59 years",
    weight_range   = "35-150 kg",
    weight_median  = "73 kg",
    sex_female_pct = 42.0,
    race_ethnicity = c(Caucasian = 85.1, Asian = 12.6, Black = 1.9, Other = 0.3),
    disease_state  = "Adult patients with various advanced solid tumours, including ovarian, prostate, non-small-cell lung, gastric, colorectal, renal cell, gastrointestinal stromal, and soft-tissue sarcoma cancers (pooled across 19 Phase I and II monotherapy and combination-chemotherapy studies). Includes 17 ovarian cancer patients and 232 patients on antihypertensive co-medication.",
    dose_range     = "Oral cediranib (cediranib maleate tablets) 0.5-90 mg once daily; majority of subjects on 20 mg (24%), 30 mg (31%), or 45 mg (42%) starting doses.",
    regions        = "Multinational (USA, Canada, Europe, Japan).",
    platinum_chemo_pct = 7.0,
    notes          = "Baseline demographics from Li 2017 Table 1. Race was reduced to Asian vs. non-Asian for the covariate analysis because African American (n=12) and Other (n=2) subjects were too few to support separate effect estimation. Race, sex, and platinum-containing chemotherapy were tested but not retained in the final covariate model; only WT and AGE survived backward elimination."
  )

  ini({
    # Structural PK -- Li 2017 Table 2 final-model estimates. Time in hours;
    # apparent clearances (CL/F, Q/F) in L/h; apparent volumes (Vc/F, Vp/F)
    # in L; first-order absorption rate Ka in 1/h; zero-order absorption
    # duration D1 in h. The Table 2 reference subject is a 73-kg, 59-year-old
    # adult; covariate-equation factors below recreate individual-specific
    # CL/F and Vc/F inside model(). Concentration units are ng/mL after the
    # 1000x conversion from internal mg/L (dose mg, vc L) inside model().
    lka <- log(2.70)  ; label("First-order absorption rate constant Ka (1/h)")                                # Li 2017 Table 2 final Ka = 2.70 1/h
    lcl <- log(26.3)  ; label("Apparent clearance CL/F at WT 73 kg, Age 59 y (L/h)")                          # Li 2017 Table 2 final CL/F = 26.3 L/h
    lvc <- log(489)   ; label("Apparent central volume Vc/F at WT 73 kg (L)")                                  # Li 2017 Table 2 final Vc/F = 489 L
    lq  <- log(11.8)  ; label("Apparent inter-compartmental clearance Q/F (L/h)")                              # Li 2017 Table 2 final Q/F = 11.8 L/h
    lvp <- log(213)   ; label("Apparent peripheral volume Vp/F (L)")                                           # Li 2017 Table 2 final Vp/F = 213 L
    ld1 <- log(1.68)  ; label("Zero-order absorption duration D1 (h)")                                         # Li 2017 Table 2 final D1 = 1.68 h

    # Bioavailability anchor -- Li 2017 Methods (page 1727): "F1 ... fixed to
    # the reference value (F = 1)" for sparse-sampling occasions; rich-sampling
    # F1 carried IOV (44.6% CV) that is dropped from this single-occasion
    # simulation library model (see vignette Assumptions and deviations).
    lfdepot <- fixed(log(1.0)) ; label("Bioavailability into depot (F1, fixed)")                               # Li 2017 Table 2 F1 = 1 (fixed)

    # Covariate effects -- Li 2017 page 1727 final covariate equations:
    #   CL/F = 26.3 * (Age/59)^-0.409 * (WT/73)^0.517
    #   Vc/F = 489  * (WT/73)^0.65
    e_age_cl <- -0.409 ; label("Age power exponent on CL/F (Age/59 y; unitless)")                              # Li 2017 Table 2 AGECL = -0.409
    e_wt_cl  <-  0.517 ; label("Body-weight power exponent on CL/F (WT/73 kg; unitless)")                      # Li 2017 Table 2 WTCL = 0.517
    e_wt_vc  <-  0.65  ; label("Body-weight power exponent on Vc/F (WT/73 kg; unitless)")                      # Li 2017 Table 2 WTVc = 0.65

    # Inter-individual variability -- Li 2017 Table 2 reports IIV as %CV.
    # Convert to log-scale variance via omega^2 = log(CV^2 + 1):
    #   CL/F : 53.7% CV  -> log(0.537^2 + 1) = 0.2532
    #   Vc/F : 61.5% CV  -> log(0.615^2 + 1) = 0.3205
    #   Ka   : 151%  CV  -> log(1.51^2  + 1) = 1.1879
    # CL/F and Vc/F are correlated (r = 0.839); Ka is independent. Off-diagonal
    # covariance: cov(CL/F, Vc/F) = 0.839 * sqrt(0.2532) * sqrt(0.3205) = 0.2390.
    etalcl + etalvc ~ c(0.2532,
                        0.2390, 0.3205)                                                                        # Li 2017 Table 2 IIV CL/F 53.7%, IIV Vc/F 61.5%, correlation 0.839
    etalka          ~ 1.1879                                                                                   # Li 2017 Table 2 IIV Ka 151%

    # Residual unexplained variability -- Li 2017 reports an additive error on
    # the log-transformed concentration (proportional in nlmixr2 linear space)
    # with separate magnitudes for rich (26.5% CV) and sparse (47.3% CV)
    # sampling occasions. The library uses the rich-sampling value as the
    # representative analytical-method residual error; the sparse-sampling
    # excess (which the paper attributes to imputed dosing times and
    # uncaptured IOV; see Discussion page 1731) is documented as a deviation
    # in the validation vignette rather than embedded as an occasion-flagged
    # mixture that this single-occasion simulation library cannot exercise.
    propSd <- 0.265 ; label("Proportional residual error (fraction; rich-sampling estimate)")                  # Li 2017 Table 2 residual variability for rich profile = 26.5% CV
  })

  model({
    # Individual PK parameters with Li 2017 covariate equations. Reference
    # subject: 73-kg, 59-year-old adult.
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (AGE / 59) ^ e_age_cl * (WT / 73) ^ e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 73) ^ e_wt_vc
    q  <- exp(lq)
    vp <- exp(lvp)
    d1 <- exp(ld1)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment oral PK with sequential zero- and first-order absorption:
    # the dose enters `depot` over a zero-order window of duration D1
    # (representing tablet dissolution / gastric emptying), then is absorbed
    # first-order at rate Ka into `central`. Bioavailability F1 is fixed at 1.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    dur(depot) <- d1
    f(depot)   <- exp(lfdepot)

    # Cediranib plasma concentration: dose mg, vc L -> mg/L; multiply by 1000
    # to report in ng/mL (the unit used in the source paper's exposure
    # comparisons against in vitro VEGFR IC50 values).
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
