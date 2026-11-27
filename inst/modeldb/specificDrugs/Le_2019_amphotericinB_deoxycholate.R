Le_2019_amphotericinB_deoxycholate <- function() {
  description <- "Two-compartment IV-infusion population PK model for amphotericin B deoxycholate in adults with HIV-associated disseminated talaromycosis (Talaromyces marneffei), with clearance an additive linear function of body weight, linked to a pharmacodynamic model of monophasic first-order killing of T. marneffei in blood driven by a sigmoid Emax function of the total plasma amphotericin B concentration (no fungal growth term) (Le 2019)"
  reference <- "Le T, Ly VT, Thu NTM, Nguyen A, Thanh NT, Chau NVV, Thwaites G, Perfect J, Kolamunnage-Dona R, Hope W. Population pharmacodynamics of amphotericin B deoxycholate for disseminated infection caused by Talaromyces marneffei. Antimicrob Agents Chemother. 2019;63(2):e01739-18. doi:10.1128/AAC.01739-18"
  vignette <- "Le_2019_amphotericinB_deoxycholate"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Le 2019 Materials and Methods, 'Population pharmacokinetic modeling': X(1)
  # and X(2) are the masses (mg) of amphotericin B in the central and
  # peripheral compartments, and Y(1) = X(1)/Vc is the plasma concentration
  # measured by HPLC. 'Pharmacokinetic and pharmacodynamic modeling': N is the
  # T. marneffei density in blood; Table 4 gives its initial condition in
  # CFU/mL, from quantitative blood cultures.
  compartmentData <- list(
    central = list(analyte = "amphotericin B", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "amphotericin B", units = "mg", specimen = "plasma", verified = TRUE),
    cfu = list(analyte = "Talaromyces marneffei", units = "CFU/mL", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Le 2019 final PK model: CL = int + slope * weight, an ADDITIVE LINEAR term on the natural scale with a separately estimated intercept (L/h) and slope (L/h per kg), not an allometric power and not centred on a reference weight. Weight was the only clinically relevant covariate available. Time-fixed at the subject level. Cohort median 45 kg (IQR 39-45) in the 78 PK patients and 44 kg (IQR 38-50) in the 55 PK-PD patients (Le 2019 Table 1).",
      source_name = "weight"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 78L,
    n_studies = 1L,
    age_range = "adults > 18 years; median 32 years (IQR 27-37)",
    age_median = "32 years",
    weight_range = "IQR 39-45 kg",
    weight_median = "45 kg",
    sex_female_pct = 32.1,
    race_ethnicity = c(Vietnamese = 100),
    disease_state = "Advanced HIV infection (median CD4 count 7 cells/uL, IQR 4-14) with culture-confirmed disseminated talaromycosis (Talaromyces marneffei isolated from blood, skin lesions, lymph nodes and/or serous fluid). Median serum creatinine 62 umol/L (IQR 52-80). 57 of 78 (73.1%) were blood-culture positive, median baseline fungal count 2.53 log10 CFU/mL (IQR 1.40-3.40). Patients with CNS infection, transaminases > 400 U/L, neutrophils < 500 cells/mm^3, creatinine clearance < 30 mL/min, interacting medication, or more than 48 h of prior antifungal therapy were excluded.",
    dose_range = "Amphotericin B deoxycholate 0.7 mg/kg once daily by intravenous infusion for induction therapy (14 days)",
    regions = "Vietnam (Hospital for Tropical Diseases, Ho Chi Minh City)",
    notes = "Pharmacokinetic-pharmacodynamic substudy of the IVAP randomised trial (itraconazole versus amphotericin B for talaromycosis, ISRCTN59144167), patients recruited October 2012 to December 2015. The PK model was fitted to 78 patients; the PD model to the 55 of them with serial quantitative blood cultures (Le 2019 Table 1). Amphotericin B assayed in plasma by HPLC-UV, calibration range 0.05-8.0 mg/L, LLOQ 0.05 mg/L, CV < 9.3%. Quantitative blood cultures daily for the first 4 days then every other day until negative over the first 14 days (100 uL of blood plated). MICs (CLSI M27-A3) ranged 0.25-1 mg/L, modal 0.5 mg/L. Fitting used the nonparametric adaptive grid (NPAG) algorithm of Pmetrics in two steps: PK first, then the PD model with each patient's Bayesian posterior PK parameters fixed."
  )

  ini({
    # Pharmacokinetics: Le 2019 Table 3 (final PK model). Pmetrics NPAG
    # reports the mean, median and SD of each parameter's nonparametric
    # marginal distribution. The MEANS are the typical values, following the
    # convention of the same group's sibling Pmetrics amphotericin B models
    # (Stott 2018, Stott 2022); the medians are recorded in each comment and
    # compared in the vignette.
    lcl <- log(0.320)
    label("Clearance intercept of the linear CL-versus-weight relation (L/h)") # Le 2019 Table 3, 'Intercept (liters/h)' mean = 0.320 (median 0.069, SD 0.430)
    e_wt_cl <- 0.048
    label("Body-weight slope of the linear CL-versus-weight relation (L/h per kg)") # Le 2019 Table 3, 'Slope' mean = 0.048 (median 0.049, SD 0.018); units L/h/kg per the Table 3 footnote
    lvc <- log(23.907)
    label("Central volume of distribution (L)") # Le 2019 Table 3, 'Vc (liters)' mean = 23.907 (median 22.650, SD 13.836)
    lk12 <- log(8.996)
    label("Central-to-peripheral first-order rate constant Kcp (1/h)") # Le 2019 Table 3, 'Kcp (h-1)' mean = 8.996 (median 1.001, SD 12.812)
    lk21 <- log(4.237)
    label("Peripheral-to-central first-order rate constant Kpc (1/h)") # Le 2019 Table 3, 'Kpc (h-1)' mean = 4.237 (median 0.184, SD 8.003)

    # Pharmacodynamics: Le 2019 Table 4 (final PD model), means as above.
    lemax <- log(0.133)
    label("Maximal first-order rate of drug-induced fungal killing Kkillmax (1/h)") # Le 2019 Table 4, 'Kkillmax' mean = 0.133 (median 0.118, SD 0.096); Table 4 prints the unit 'no. of log10 CFU/ml/h' but the PD equation makes it a natural-log first-order rate, see vignette
    lhill <- log(13.815)
    label("Hill coefficient of the killing function H (unitless)") # Le 2019 Table 4, 'H' mean = 13.815 (median 15.447, SD 5.964)
    lec50 <- log(0.016)
    label("Plasma amphotericin B concentration giving half-maximal killing C50k (mg/L)") # Le 2019 Table 4, 'C50k (mg/liter)' mean = 0.016 (median 0.010, SD 0.016)
    lrbase <- log(4318.564)
    label("Fungal density in blood at the start of treatment IC (CFU/mL)") # Le 2019 Table 4, 'IC (no. of CFU/ml)' mean = 4,318.564 (median 733.770, SD 7,560.737)

    # Inter-individual variability. Pmetrics NPAG estimates a discrete joint
    # distribution over support points and Le 2019 reports only each
    # marginal's mean and SD, so each marginal is approximated by an
    # independent log-normal with the same coefficient of variation,
    # omega^2 = log(1 + (SD / mean)^2), carrying the mean as the log-normal
    # median. No correlations are published.
    etalcl ~ 1.03164 # Table 3, 'Intercept' SD 0.430 on mean 0.320 -> log(1 + (0.430/0.320)^2)
    etae_wt_cl ~ 0.131576 # Table 3, 'Slope' SD 0.018 on mean 0.048 -> log(1 + (0.018/0.048)^2)
    etalvc ~ 0.288889 # Table 3, 'Vc' SD 13.836 on mean 23.907 -> log(1 + (13.836/23.907)^2)
    etalk12 ~ 1.10801 # Table 3, 'Kcp' SD 12.812 on mean 8.996 -> log(1 + (12.812/8.996)^2)
    etalk21 ~ 1.51901 # Table 3, 'Kpc' SD 8.003 on mean 4.237 -> log(1 + (8.003/4.237)^2)
    etalemax ~ 0.419369 # Table 4, 'Kkillmax' SD 0.096 on mean 0.133 -> log(1 + (0.096/0.133)^2)
    etalhill ~ 0.170897 # Table 4, 'H' SD 5.964 on mean 13.815 -> log(1 + (5.964/13.815)^2)
    etalec50 ~ 0.693147 # Table 4, 'C50k' SD 0.016 on mean 0.016 -> log(1 + 1^2)
    etalrbase ~ 1.40245 # Table 4, 'IC' SD 7560.737 on mean 4318.564 -> log(1 + (7560.737/4318.564)^2)

    # Residual error. Pmetrics weighted the PK observations by an assay-error
    # polynomial, and the PD fit used a per-patient variance model
    # (intercept + slope * Y)^2 estimated in ADAPT 5; Le 2019 publishes the
    # values of neither. The PK value below is the HPLC assay imprecision
    # reported in Materials and Methods (the same stand-in the sibling Stott
    # 2018 model uses), so that stochastic simulation is possible; it is not a
    # fitted residual error. The PD error keeps the published form (SD linear
    # in the prediction, nlmixr2 combined1) with both terms fixed at zero
    # because no magnitude is reported.
    propSd <- fixed(0.093)
    label("Proportional residual SD for plasma amphotericin B, taken from the HPLC assay imprecision (fraction)") # Le 2019 Materials and Methods, 'Bioanalytical methodology': 'The coefficient of variation was <9.3% over the concentration range of 0.05 to 8 mg/liter'
    addSd_cfu <- fixed(0)
    label("Additive residual SD for blood fungal density (CFU/mL); value not published") # Le 2019 Materials and Methods: variance = [intercept + slope x Y(1)]^2, values not reported
    propSd_cfu <- fixed(0)
    label("Proportional residual SD for blood fungal density (fraction); value not published") # Le 2019 Materials and Methods: variance = [intercept + slope x Y(1)]^2, values not reported
  })
  model({
    # Final PK model: CL = int + slope * weight. In a nonparametric fit both
    # the intercept and the slope are subject-level parameters, so each
    # carries its own log-normal eta.
    wt_cl <- e_wt_cl * exp(etae_wt_cl)
    cl <- exp(lcl + etalcl) + wt_cl * WT
    vc <- exp(lvc + etalvc)
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)
    kel <- cl / vc

    emax <- exp(lemax + etalemax)
    hill <- exp(lhill + etalhill)
    ec50 <- exp(lec50 + etalec50)
    rbase <- exp(lrbase + etalrbase)

    # Standard two-compartment mammillary model with zero-order infusion R(1)
    # into central (encoded on the dose rows). The printed equations omit the
    # -Kcp*X(1) loss from central and write -Kcp*X(2) for the return from
    # peripheral; the prose describes Kcp and Kpc moving mass back and forth,
    # which only this form does. See the vignette.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Y(1) = X(1)/Vc, mg/L.
    Cc <- central / vc

    # PD: dN/dt = -Kkillmax * C^H / (C^H + C50k^H) * N, with C = X(1)/Vc
    # (the printed X(1)^H / Vc is a typesetting slip; C50k is a
    # concentration). No growth term: untreated growth could not be observed.
    # N(0) is the estimated fungal density at the start of treatment.
    ckill <- Cc^hill / (Cc^hill + ec50^hill)
    d/dt(cfu) <- -emax * ckill * cfu
    cfu(0) <- rbase

    Cc ~ prop(propSd)
    cfu ~ add(addSd_cfu) + prop(propSd_cfu) + combined1()
  })
}
