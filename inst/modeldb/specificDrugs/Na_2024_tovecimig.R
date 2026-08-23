Na_2024_tovecimig <- function() {
  description <- "Two-compartment population PK model for tovecimig (ABL001/CTX-009), a bispecific antibody targeting DLL4 and VEGF-A, in adult patients with relapsed or refractory solid tumors, with parallel linear and Michaelis-Menten elimination from the central compartment and a power body-weight effect on the central volume."
  reference   <- "Na JY, Jeon J, Huh KY, Eom J, Ahn J, You WK, Oh J. Population pharmacokinetic model of ABL001/CTX-009 (anti-VEGF/DLL4) in adult cancer patients with solid tumor. Cancer Sci. 2024;115(12):3943-3951. doi:10.1111/cas.16363"
  vignette    <- "Na_2024_tovecimig"
  units       <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline body weight. Power effect on the central volume V1 only: V1 = theta_V1 * (WT/70)^0.598. Reference 70 kg is the 'typically accepted value' the authors centered on (Na 2024 Methods Section 2.4, Equation 1; Table 2 footnote b). The only covariate retained in the final model; it explained approximately 7% of the between-subject variability in V1 (Na 2024 Discussion). Cohort mean 64.2 kg (SD 16.4), range 35.6-110.2 kg.",
      source_name        = "Weight"
    )
  )

  compartmentData <- list(
    central     = list(analyte = "tovecimig", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "tovecimig", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 31L,
    n_observations = 712L,
    n_studies      = 1L,
    phases         = "Phase 1 (first-in-human dose escalation), NCT03292783",
    age_range      = "25-81 years",
    age_mean       = "53.5 years (SD 12.4)",
    weight_range   = "35.6-110.2 kg",
    weight_mean    = "64.2 kg (SD 16.4)",
    sex_female_pct = 38.7,
    race_ethnicity = "Not reported; all three study sites were in the Republic of Korea",
    disease_state  = "Adults (aged >= 19 years) with progressive locally advanced or metastatic solid tumors relapsed or refractory to standard therapy; ECOG 0-2 (87% ECOG 1). Primary tumor locations: stomach 38.7%, colon/rectal 35.5%, other 16.1%, biliary tract / liver / ovary 3.2% each. Median 5 prior lines of treatment.",
    dose_range     = "0.3, 1, 2.5, 5, 7.5, 10, 12.5, 15, and 17.5 mg/kg as a 1-h IV infusion in 2-week cycles (the second dose used a 3-week interval)",
    regions        = "Republic of Korea (Seoul National University Hospital Bundang, Samsung Medical Center, Asan Medical Center)",
    notes          = "Demographics from Na 2024 Table 1 (dose-escalation cohort, N = 31; the dose-expansion cohort of 14 patients was NOT used for model development). Only 6 of 712 samples (< 1%) were below the 250 ng/mL LLOQ and were treated as missing. Covariates screened by stepwise selection included dose group (categorical, exponential form) and body weight (continuous, power form centered on 70 kg); only body weight on V1 was retained. Immunogenicity (anti-drug antibody) and tumor burden were not evaluated as covariates (Na 2024 Discussion, limitations)."
  )

  ini({
    # Structural PK parameters - Na 2024 Table 2 final-model estimates.
    # Reference covariate: WT = 70 kg. Time in hours, amounts in mg, concentrations in mg/L.
    lcl   <- log(0.0184);        label("Linear clearance CL (L/h)")                                      # Na 2024 Table 2 (RSE 18%)
    lvc   <- log(3.87);          label("Central volume of distribution V1 at WT = 70 kg (L)")            # Na 2024 Table 2 (RSE 5%)
    lq    <- log(0.0207);        label("Intercompartmental clearance Q (L/h)")                           # Na 2024 Table 2 (RSE 16%)
    lvp   <- log(1.31);          label("Peripheral volume of distribution V2 (L)")                       # Na 2024 Table 2 (RSE 28%)
    lvmax <- log(0.0968);        label("Maximum Michaelis-Menten elimination rate Vmax (mg/h)")          # Na 2024 Table 2 (RSE 26%; printed unit 'ug/h' is a mislabel - see vignette Errata)
    lkm   <- fixed(log(4.15));   label("Michaelis-Menten constant Km (mg/L)")                            # Na 2024 Table 2, fixed at the value estimated in exploration model 4 (Table S1); printed unit 'ug/L' is a mislabel - see vignette Errata

    # Covariate effect - power exponent of body weight on the central volume.
    e_wt_vc <- 0.598;            label("Power exponent of WT/70 on central volume V1 (unitless)")        # Na 2024 Table 2 (theta_weight,V1; RSE 32%) and Table 2 footnote b

    # Inter-individual variability. Na 2024 Table 2 reports IIV as %CV for a log-normal
    # (exponential) random-effect model, so omega^2 = log(CV^2 + 1).
    # IIV on Vmax and Q did not improve the model and was not estimated (Na 2024 Section 3.2).
    etalcl ~ 0.258394                                                                                    # Na 2024 Table 2 (IIV CL 54.3 %CV, shrinkage 1%); log(0.543^2 + 1)
    etalvc ~ 0.024967                                                                                    # Na 2024 Table 2 (IIV V1 15.9 %CV, shrinkage 9%); log(0.159^2 + 1)
    etalvp ~ 0.699147                                                                                    # Na 2024 Table 2 (IIV V2 100.6 %CV, shrinkage 19%); log(1.006^2 + 1)

    # Residual error - proportional only (Na 2024 Section 3.2).
    propSd <- 0.245;             label("Proportional residual error (SD, fraction)")                      # Na 2024 Table 2 (RSE 1%); read as an SD, i.e. 24.5% CV - see vignette Errata
  })

  model({
    # Individual PK parameters. Body weight enters as a power term on V1 only,
    # centered on 70 kg (Na 2024 Table 2 footnote b).
    cl   <- exp(lcl + etalcl)
    vc   <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q    <- exp(lq)
    vp   <- exp(lvp + etalvp)
    vmax <- exp(lvmax)
    km   <- exp(lkm)

    # Plasma concentration in the central compartment.
    Cc <- central / vc

    # Two-compartment disposition with parallel linear (CL) and saturable
    # Michaelis-Menten elimination from the central compartment (Na 2024 Figure 1).
    # The MM term vmax * Cc / (km + Cc) is in mg/h when vmax is in mg/h.
    d/dt(central)     <- -(cl / vc) * central -
                          vmax * Cc / (km + Cc) -
                          (q / vc) * central +
                          (q / vp) * peripheral1
    d/dt(peripheral1) <-  (q / vc) * central - (q / vp) * peripheral1

    Cc ~ prop(propSd)
  })
}
