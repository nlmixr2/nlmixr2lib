Horvatits_2014_ganciclovir <- function() {
  description <- paste(
    "Two-compartment population PK model for IV ganciclovir in critically ill",
    "patients receiving continuous venovenous hemodiafiltration (Horvatits 2014).",
    "No covariates were retained because of the small cohort (n = 9), so every",
    "disposition parameter is a covariate-free typical value with between-subject",
    "variability.",
    "Parameters transcribed from the Yang 2023 ganciclovir / valganciclovir",
    "population-PK model repository (Table 3), not from the primary publication;",
    "re-verify against Horvatits 2014 when the primary is obtained.",
    sep = " "
  )
  reference <- paste(
    "Horvatits T, Kitzberger R, Drolz A, Zauner C, Jager W, Bohmdorfer M, Kraff S,",
    "Fritsch A, Thalhammer F, Fuhrmann V, et al. Pharmacokinetics of ganciclovir",
    "during continuous venovenous hemodiafiltration in critically ill patients.",
    "Antimicrob Agents Chemother. 2014;58(1):94-101. doi:10.1128/AAC.00892-13.",
    "Parameters transcribed from Yang W, Mak W, Gwee A, Gu M, Wu Y, Shi Y, He Q,",
    "Xiang X, Han B, Zhu X. Establishment and Evaluation of a Parametric Population",
    "Pharmacokinetic Model Repository for Ganciclovir and Valganciclovir.",
    "Pharmaceutics. 2023;15(7):1801. doi:10.3390/pharmaceutics15071801 (Table 3).",
    sep = " "
  )
  vignette <- "Yang_2023_ganciclovir_model_repository"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # No covariateData: Yang 2023 Table 4 records "NR" for every covariate column of
  # this study, with the footnote "Covariates were not included in the model due to
  # the small number of patients (9 patients)." Continuous venovenous
  # hemodiafiltration characterises the entire cohort rather than acting as a
  # covariate, and is recorded in `population` below.

  population <- list(
    species        = "human",
    n_subjects     = 9L,
    n_studies      = 1L,
    age_median     = "mean 56 (SD 9) years",
    weight_median  = "mean 86 (SD 25) kg",
    sex_female_pct = 11.1,
    race_ethnicity = "Not reported.",
    disease_state  = paste(
      "Critically ill patients with suspected or proven CMV infection, all",
      "receiving continuous venovenous hemodiafiltration (CVVHDF) as",
      "renal-replacement therapy.",
      sep = " "
    ),
    renal_function = paste(
      "All subjects were on CVVHDF. The final model controlled for the",
      "extracorporeal therapy given to the patients (Yang 2023 Section 3.3.1), so",
      "the reported clearance of 2.2 L/h is the total clearance of a",
      "CVVHDF-supported patient rather than a native renal clearance.",
      sep = " "
    ),
    dose_range     = "IV ganciclovir 5 mg/kg infused over 0.5 h via a central line.",
    regions        = paste(
      "Reported as 'Australia (Prospective)' in Yang 2023 Table 2. Flagged for",
      "verification against the primary publication: the Horvatits 2014",
      "author group is based in Vienna, so the country field in the review's",
      "table may be a transcription error for Austria.",
      sep = " "
    ),
    bioassay       = "HPLC, LLOQ 5 ng/mL.",
    notes          = paste(
      "Demographics and dosing from Yang 2023 Table 2. Intensive sampling at",
      "0 (pre-dose), 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12 and 24 h post dose; the number",
      "of observations was not reported. No covariates were investigated because of",
      "the limited cohort size (Yang 2023 Table 4 footnote). Yang 2023 Section 3.3.1",
      "notes that the simulated half-life from this model was longer than that of",
      "the other repository models because the underlying data came from critically",
      "ill patients on CVVHDF; Section 4.1 additionally notes that this model was",
      "limited by the small cohort and that its between-subject variability was not",
      "clearly explained. Simulation targets in the primary study were an AUC0-24h",
      "of 50 mg*h/L and a trough concentration above 2 mg/L.",
      sep = " "
    )
  )

  ini({
    # Structural PK -- Yang 2023 Table 3, Horvatits et al. (2014) row. No
    # covariates; these are the covariate-free typical values. Clearances in L/h,
    # volumes in L. IV-only model, so no absorption parameters.
    lcl <- log(2.2) ; label("Total clearance during CVVHDF (CL, L/h)")            # Yang 2023 Table 3 (Horvatits 2014): CL = 2.2
    lvc <- log(32.4); label("Central volume of distribution (Vc, L)")             # Yang 2023 Table 3 (Horvatits 2014): Vc = 32.4
    lq  <- log(16.8); label("Inter-compartmental clearance (Q, L/h)")             # Yang 2023 Table 3 (Horvatits 2014): Q = 16.8
    lvp <- log(33.5); label("Peripheral volume of distribution (Vp, L)")          # Yang 2023 Table 3 (Horvatits 2014): Vp = 33.5

    # Between-subject variability. Yang 2023 Methods: %CV = sqrt(omega^2) * 100%,
    # so variance = (BSV% / 100)^2. This is the only model in the repository that
    # reports BSV on all four disposition parameters.
    etalcl ~ 0.378225  # Yang 2023 Table 3 (Horvatits 2014): BSV CL = 61.5% -> 0.615^2
    etalvc ~ 0.112896  # Yang 2023 Table 3 (Horvatits 2014): BSV Vc = 33.6% -> 0.336^2
    etalq  ~ 0.120409  # Yang 2023 Table 3 (Horvatits 2014): BSV Q  = 34.7% -> 0.347^2
    etalvp ~ 0.367236  # Yang 2023 Table 3 (Horvatits 2014): BSV Vp = 60.6% -> 0.606^2

    # Residual unexplained variability: proportional.
    propSd <- 0.0722; label("Proportional residual error (fraction)")  # Yang 2023 Table 3 (Horvatits 2014): 7.22% proportional error
  })

  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q  <- exp(lq  + etalq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Dose in mg, volume in L -> concentration in mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
