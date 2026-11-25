Felton_2014_piperacillin <- function() {
  description <- "Two-compartment IV-infusion population PK model for piperacillin in critically ill adults with nosocomial infection, fitted nonparametrically (Pmetrics NPAG) and built as the Bayesian prior of the BestDose dose-optimisation software. Clearance is an additive non-renal intercept plus a slope on Cockcroft-Gault creatinine clearance, and central volume is an additive intercept plus a slope on body weight; distribution uses explicit k12 / k21 micro-constants (Felton 2014, 'linear with covariates' model)"
  reference <- "Felton TW, Roberts JA, Lodise TP, Van Guilder M, Boselli E, Neely MN, Hope WW. Individualization of piperacillin dosing for critically ill patients: dosing software to optimize antimicrobial therapy. Antimicrob Agents Chemother. 2014;58(7):4094-4102. doi:10.1128/AAC.02664-14"
  vignette <- "Felton_2014_piperacillin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Felton 2014 Fig. 2 legend: X1 and X2 are the amounts of piperacillin (mg)
  # in the central and peripheral compartments; concentrations were measured
  # in plasma by HPLC (Materials and Methods, 'Pharmacokinetic studies').
  compartmentData <- list(
    central = list(analyte = "piperacillin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "piperacillin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Estimated creatinine clearance by the Cockcroft-Gault equation",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Raw Cockcroft-Gault creatinine clearance in mL/min, NOT BSA-normalized (Felton 2014 Materials and Methods cites Cockcroft and Gault; Table 1 reports it in ml/min). Enters clearance as an additive linear term, CL = Cli + CRCL * Cls (Fig. 2, equation 1d), with no reference value; the intercept Cli is the clearance extrapolated to zero creatinine clearance. The mL/min scale is what makes the slope dimensionally sensible: at the cohort mean of 85.9 mL/min the typical CL is 13.3 L/h, whereas reading CRCL in L/h would give 4.4 L/h and predict concentrations several-fold above the Fig. 4 observations. Combined-cohort mean 85.9, median 81.0, range 14.0-245.7 mL/min (Table 1).",
      source_name = "Cr_cl"
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters central volume as an additive linear term, Vc = Vi + WT * Vs (Felton 2014 Fig. 2, equation 1d), not an allometric power and not centred on a reference weight. Combined-cohort mean 71.7, median 69.0, range 38.1-122.5 kg (Table 1).",
      source_name = "Wt"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 146L,
    n_studies = 3L,
    n_observations = 803L,
    age_range = "18-88 years",
    age_median = "56 years (mean 55.8)",
    weight_range = "38.1-122.5 kg",
    weight_median = "69.0 kg (mean 71.7)",
    sex_female_pct = 40.4,
    disease_state = "Critically ill or hospitalised adults with nosocomial infection: patients undergoing abdominal, thoracic or colorectal surgery (Lodise et al.), hospitalised patients (Felton et al.), and ventilated patients with ventilator-associated pneumonia (Boselli et al.).",
    dose_range = "Piperacillin 2-4 g by 30-min infusion every 6-8 h, 3 g by 4-h extended infusion every 8 h, or 12-16 g/day by continuous infusion.",
    renal_function = "Cockcroft-Gault creatinine clearance mean 85.9 mL/min, median 81.0, range 14.0-245.7",
    regions = "Pooled from three previously published studies (Felton et al., Boselli et al., Lodise et al.)",
    notes = "Felton 2014 Table 1 (Combined column) and Materials and Methods. Lodise et al.: 76 surgical patients given 2 g over 30 min (first dosing interval sampled) plus 18 colorectal-surgery patients given 4 g over 30 min every 6 h (first and second intervals). Felton et al.: 12 patients given 3 g over 4 h every 8 h, sampled mainly on day 3. Boselli et al.: 40 ventilated patients given 12 or 16 g/day by continuous infusion, sampled during the second 24 h. 803 concentrations, 2-10 per patient; sex 87 male : 59 female. Fitted with the nonparametric adaptive grid (NPAG) algorithm of Pmetrics 1.1.1. External validation used 8 critically ill patients from Roberts et al. (4 g over 20 min every 6 or 8 h; creatinine clearance mean 165.5 mL/min, weight mean 86.8 kg), who are not part of the estimation cohort."
  )

  ini({
    # Structural parameters: Felton 2014 Table 3 ('Mean [median] +/- SD' of the
    # NPAG marginal distributions for the linear-clearance-with-covariates
    # model, the model carried into the dose-optimisation software). The MEANS
    # are used as typical values: at the validation cohort's covariates the
    # means reproduce the observed first-dose decline of Fig. 4, whereas the
    # medians (much smaller Vi, Vs, kcp and kpc) give a half-life of about
    # 0.7 h and concentrations several-fold below the observations from 2 h
    # on. See the vignette 'Assumptions and deviations'.
    lcl_nonren <- log(3.83)
    label("Non-renal clearance, the intercept of the linear CL-versus-CrCl relation (L/h)") # Table 3, 'Cli' 3.83 [2.79] +/- 3.35
    e_crcl_cl_renal <- 0.11
    label("Creatinine-clearance slope of the linear CL-versus-CrCl relation (L/h per mL/min)") # Table 3, 'Cls' 0.11 [0.10] +/- 0.07
    lvc <- log(4.54)
    label("Central volume intercept of the linear Vc-versus-weight relation (L)") # Table 3, 'Vi' 4.54 [1.83] +/- 4.43
    e_wt_vc <- 0.12
    label("Body-weight slope of the linear Vc-versus-weight relation (L per kg)") # Table 3, 'Vs' 0.12 [0.06] +/- 0.12
    lk12 <- log(6.74)
    label("Central-to-peripheral first-order rate constant kcp (1/h)") # Table 3, 'kcp' 6.74 [0.85] +/- 11.38
    lk21 <- log(9.14)
    label("Peripheral-to-central first-order rate constant kpc (1/h)") # Table 3, 'kpc' 9.14 [1.65] +/- 13.13

    # Inter-individual variability. NPAG estimates a discrete joint density
    # and Table 3 reports only each marginal's mean, median and SD; each
    # marginal is approximated by a log-normal with the same coefficient of
    # variation, omega^2 = log(1 + (SD/mean)^2), centred on the mean. No
    # covariances are reported, so the etas are independent.
    etalcl_nonren ~ 0.568181 # Table 3, 'Cli' SD 3.35 on mean 3.83 -> log(1 + (3.35/3.83)^2)
    etae_crcl_cl_renal ~ 0.340008 # Table 3, 'Cls' SD 0.07 on mean 0.11 -> log(1 + (0.07/0.11)^2)
    etalvc ~ 0.668921 # Table 3, 'Vi' SD 4.43 on mean 4.54 -> log(1 + (4.43/4.54)^2)
    etae_wt_vc ~ 0.693147 # Table 3, 'Vs' SD 0.12 on mean 0.12 -> log(1 + (0.12/0.12)^2)
    etalk12 ~ 1.348277 # Table 3, 'kcp' SD 11.38 on mean 6.74 -> log(1 + (11.38/6.74)^2)
    etalk21 ~ 1.119609 # Table 3, 'kpc' SD 13.13 on mean 9.14 -> log(1 + (13.13/9.14)^2)

    # Residual error: the Pmetrics assay error polynomial SD = 1.04 + 0.14 * C
    # (Materials and Methods, 'Development of a population pharmacokinetic
    # model'), set by the authors from a regression on the Boselli
    # continuous-infusion data before the fit. Pmetrics sums the two terms
    # linearly, hence combined1(). The variance was also multiplied by an
    # estimated scalar gamma whose value is not reported, so these SDs
    # correspond to gamma = 1.
    addSd <- fixed(1.04)
    label("Additive term of the assay error polynomial (mg/L)") # Methods: SD = 1.04 + (0.14 x C)
    propSd <- fixed(0.14)
    label("Proportional term of the assay error polynomial (fraction)") # Methods: SD = 1.04 + (0.14 x C)
  })
  model({
    # Felton 2014 Fig. 2, equation 1d. In a nonparametric fit each intercept
    # and slope is its own subject-level parameter, so each carries an eta.
    cl_nonren <- exp(lcl_nonren + etalcl_nonren)
    cl_renal <- e_crcl_cl_renal * exp(etae_crcl_cl_renal) * CRCL
    cl <- cl_nonren + cl_renal
    wt_vc <- e_wt_vc * exp(etae_wt_vc)
    vc <- exp(lvc + etalvc) + wt_vc * WT
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)
    kel <- cl / vc

    # Fig. 2, equations 1d and 2; R(1), the intravenous infusion into the
    # central compartment, is supplied on the dosing rows.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # mg / L = mg/L, the units of the reported concentrations.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
