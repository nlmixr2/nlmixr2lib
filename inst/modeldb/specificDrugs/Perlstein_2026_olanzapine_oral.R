Perlstein_2026_olanzapine_oral <- function() {
  description <- "Two-compartment population PK model for oral olanzapine (ZYPREXA tablets) with first-order absorption and an absorption lag time, in healthy participants and adults with schizophrenia or schizoaffective disorder (Perlstein 2026). Companion oral comparator model to the TV-44749 subcutaneous long-acting injectable model Perlstein_2026_olanzapine_lai; the two were fitted to separate datasets in the same Phase 1 study. No covariates were retained. Carries the published Emax dopamine D2 receptor occupancy layer (Mamo 2008) as an algebraic observable."
  reference <- paste(
    "Perlstein I, Cherniakov I, Elgart A, Gomeni R, Gutman D, Merenlender Wagner A, Singh R.",
    "Population Pharmacokinetic Model-Based Dose Selection of Extended-Release Injectable Olanzapine",
    "(TV-44749) for Subcutaneous Use in Phase 3 Clinical Trial in Adults with Schizophrenia.",
    "The Journal of Clinical Pharmacology. 2026;66(1):e70144. doi:10.1002/jcph.70144.",
    "Parameter values from Supplemental Table 2 of the Supplementary Online Material.",
    "D2 receptor occupancy layer from Mamo D, Kapur S, Keshavan M, et al.",
    "D2 receptor occupancy of olanzapine pamoate depot using positron emission tomography:",
    "an open-label study in patients with schizophrenia.",
    "Neuropsychopharmacology. 2008;33(2):298-304.",
    sep = " "
  )
  vignette <- "Perlstein_2026_olanzapine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot       = list(analyte = "olanzapine", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "olanzapine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "olanzapine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight at baseline",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened in the stepwise covariate search (Methods, Covariate Analysis) and retained in the TV-44749 model, but Supplemental Table 2 reports NO covariate effect for the oral olanzapine model; the oral model is therefore covariate-free."
    ),
    AGE = list(
      description = "Age at baseline",
      units       = "years",
      type        = "continuous",
      notes       = "Screened in the stepwise covariate search but not present in the reported oral olanzapine model (Supplemental Table 2)."
    ),
    BMI = list(
      description = "Body mass index at baseline",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened in the stepwise covariate search but not present in the reported oral olanzapine model (Supplemental Table 2)."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened in the stepwise covariate search but not present in the reported oral olanzapine model (Supplemental Table 2)."
    ),
    SMOKER = list(
      description = "Current smoker indicator (1 = smoker, 0 = non-smoker)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened because smoking induces CYP1A2-mediated olanzapine clearance (Discussion) but not present in the reported oral olanzapine model (Supplemental Table 2)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 90,
    n_studies      = 1,
    age_range      = "18-60 years (healthy participants) and 18-65 years (patients); cohort mean 45.6 years (SD 10.3)",
    weight_mean    = "85.8 kg (SD 15.7)",
    bmi_mean       = "28.3 kg/m^2 (SD 4.7)",
    sex_female_pct = 27,
    race_ethnicity = c(Black = 79, White = 29),
    disease_state  = "Healthy participants, and adults with a DSM-5-confirmed diagnosis of schizophrenia or schizoaffective disorder on a clinically stable oral olanzapine regimen.",
    dose_range     = "Daily oral olanzapine (ZYPREXA) 2.5 or 5 mg in healthy participants; 10, 15 or 20 mg in patients.",
    regions        = "United States",
    notes          = paste(
      "Phase 1 study TV-44749-SAD-10154. The oral olanzapine lead-in period preceded the",
      "TV-44749 treatment period; 1196 PK samples from 90 oral olanzapine recipients entered",
      "this model (Results, Data). The model was used to quantify the carryover of residual",
      "oral olanzapine into the post-injection samples, and as the exposure comparator for",
      "TV-44749 dose selection. Race percentages are quoted as printed in the paper and sum",
      "to more than 100. Estimation in NONMEM 7.4."
    )
  )

  ini({
    # Structural parameters -- Supplemental Table 2, 'Oral olanzapine
    # population parameter estimates'. The main text never cites this table
    # by number; it is in the Supplementary Online Material PDF.
    lka   <- log(0.0292);  label("First-order oral absorption rate constant (ka, 1/h)")     # Supplemental Table 2: ka = 0.0292 1/h (RSE 3.5%)
    lvc   <- log(9.61);    label("Apparent central volume of distribution V/F (L)")         # Supplemental Table 2: V/F = 9.61 L (RSE 15.8%)
    lcl   <- log(22.6);    label("Apparent clearance CL/F (L/h)")                           # Supplemental Table 2: CL/F = 22.6 L/h (RSE 6.2%)
    lk12  <- log(0.254);   label("First-order transfer rate constant central to peripheral1 (k12, 1/h)")  # Supplemental Table 2: k12 = 0.254 1/h (RSE < 0.01%)
    lk21  <- log(0.015);   label("First-order transfer rate constant peripheral1 to central (k21, 1/h)")  # Supplemental Table 2: k21 = 0.015 1/h (RSE < 0.01%)
    ltlag <- log(0.937);   label("Oral absorption lag time (h)")                            # Supplemental Table 2: lag = 0.937 h (RSE 1%)

    # Published Emax dopamine D2 receptor occupancy layer, inherited from
    # Mamo 2008 (the paper's reference 13) rather than estimated here, so both
    # parameters are fixed. Perlstein 2026 Methods, PK/D2RO Model:
    #   D2RO = ROmax * Cp / (EC50 + Cp)
    emax  <- fixed(100);      label("Maximal attainable dopamine D2 receptor occupancy ROmax (%)")            # Methods PK/D2RO Model: ROmax fixed to 100% (Mamo 2008)
    lec50 <- fixed(log(11));  label("Plasma olanzapine concentration giving 50% D2 receptor occupancy (EC50, ng/mL)")  # Methods PK/D2RO Model: EC50 = 11 ng/mL, estimated by Mamo 2008

    # Inter-individual variability -- Supplemental Table 2, Random-effect
    # block. The rows are VARIANCES: at N = 90 subjects the RSE of a variance
    # cannot fall much below sqrt(2/N) = 15.0%, and the three rows report
    # 22.4%, 21.6% and 15.1% -- the CL/F row sitting exactly on the bound.
    # Log-normal (exponential) IIV assumed; the paper does not state the
    # transform. See the vignette Errata.
    etalka ~ 0.0444  # Supplemental Table 2 Random-effect ka = 0.0444 (RSE 22.4%, shrinkage 23%); 21.2% CV
    etalvc ~ 0.897   # Supplemental Table 2 Random-effect V/F = 0.897 (RSE 21.6%, shrinkage 20%); 121% CV
    etalcl ~ 0.227   # Supplemental Table 2 Random-effect CL/F = 0.227 (RSE 15.1%, shrinkage 0.43%); 50.5% CV

    # No IIV on k12 or k21. The Supplemental Table 2 footnote states 'The IIV
    # of k12 and k21 were fixed to zero.' The etas are omitted here rather
    # than written as `~ fixed(0)`, because a zero-variance diagonal makes
    # OMEGA singular and breaks the Cholesky sampler used by rxSolve; the two
    # encodings are otherwise identical. Same treatment as
    # Wattanakul_2024_primaquine_motherinfant.R.

    # Residual error -- Supplemental Table 2. These rows are SDs, not
    # variances: at N = 1196 observations the RSE of a variance cannot fall
    # much below sqrt(2/N) = 4.09%, and the proportional row reports 1%.
    # The additive term carries the table's `*` fixed footnote and is
    # effectively switched off at 1e-4 ng/mL.
    addSd  <- fixed(0.0001);  label("Additive residual error (ng/mL)")         # Supplemental Table 2: Additive RSE model = 0.0001, marked * (fixed)
    propSd <- 0.147;          label("Proportional residual error (fraction)")  # Supplemental Table 2: Proportional RSE model = 0.147 (RSE 1%)
  })

  model({
    # 1. Individual parameters
    ka   <- exp(lka + etalka)
    vc   <- exp(lvc + etalvc)
    cl   <- exp(lcl + etalcl)
    k12  <- exp(lk12)
    k21  <- exp(lk21)
    tlag <- exp(ltlag)

    # 2. Micro-constants
    kel <- cl / vc

    # 3. ODE system
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 4. Absorption lag
    alag(depot) <- tlag

    # 5. Observation. Dose in mg over volume in L gives mg/L; the factor of
    #    1000 converts to the ng/mL the paper reports.
    Cc <- 1000 * central / vc

    # Published Emax D2 receptor occupancy (percent), Mamo 2008 via
    # Perlstein 2026 Methods. Reported without residual error -- it is a
    # deterministic transform of the plasma concentration, not a measured
    # endpoint of this study.
    ec50 <- exp(lec50)
    D2RO <- emax * Cc / (ec50 + Cc)

    Cc ~ add(addSd) + prop(propSd)
  })
}
