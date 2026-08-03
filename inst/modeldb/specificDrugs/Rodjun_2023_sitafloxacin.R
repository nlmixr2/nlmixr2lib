Rodjun_2023_sitafloxacin <- function() {
  description <- "Population PK model for oral sitafloxacin as encoded by Rodjun 2023 for Monte Carlo probability-of-target-attainment simulation against carbapenem-, multidrug- and colistin-resistant Acinetobacter baumannii. One compartment with first-order absorption in the fasted state. Apparent clearance is entirely proportional to creatinine clearance (CL/F = 2.58 x CrCL, with CrCL in L/h) and apparent volume is proportional to body weight (V/F = 1.72 L/kg). Unbound concentration is returned as Ccu using the reported unbound fraction of 0.388, because the paper's PK/PD target is the unbound AUC ratio fAUC/MIC > 30."
  reference <- paste(
    "Rodjun V, Montakantikul P, Houngsaitong J, Jitaree K, Nosoongnoen W.",
    "Pharmacokinetic/pharmacodynamic (PK/PD) simulation for dosage optimization of",
    "colistin and sitafloxacin, alone and in combination, against carbapenem-,",
    "multidrug-, and colistin-resistant Acinetobacter baumannii.",
    "Front Microbiol. 2023;14:1275909. doi:10.3389/fmicb.2023.1275909.",
    "Parameter values (Table 2) and the unbound fraction are reproduced by Rodjun 2023",
    "from Tanigawara Y, Kaku M, Totsuka K, Tsuge H, Saito A. Population pharmacokinetics",
    "and pharmacodynamics of sitafloxacin in patients with community-acquired respiratory",
    "tract infections. J Infect Chemother. 2013;19(5):858-866.",
    "doi:10.1007/s10156-013-0580-2.",
    "See also modellib('Rodjun_2023_colistin') for the companion agent.",
    sep = " "
  )
  vignette <- "Rodjun_2023_colistin_sitafloxacin"

  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    depot   = list(analyte = "sitafloxacin", units = "mg", specimen = "administration site", verified = FALSE),
    central = list(analyte = "sitafloxacin", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance (raw, NOT BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Apparent clearance is entirely proportional to creatinine clearance: Rodjun 2023 Table 2 gives CL_t/F = 2.58 x CrCL with no intercept. The Table 2 unit column prints 'L' for this row, which cannot be right for a slope; the slope is dimensionless and requires CrCL in L/h, so model() converts the data-side mL/min to L/h via CRCL * 60 / 1000. This is confirmed numerically in the validation vignette: with CrCL in L/h the closed-form PTA reproduces all 14 evaluable cells of Rodjun 2023 Table 5 at CrCL 90 mL/min to within 0.6 percentage points, whereas taking CrCL in mL/min inflates CL/F 60-fold and drives every PTA to 0%. Rodjun 2023 simulated the four discrete values 90, 50, 30 and 10 mL/min (Materials and methods, Simulated dosage regimens). Stored under canonical CRCL per inst/references/covariate-columns.md, which accepts raw mL/min when the source paper does not apply BSA normalization.",
      source_name        = "CrCL"
    ),
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Rodjun 2023 Table 2 reports the apparent volume of distribution per kilogram (V/F = 1.72 L/kg), so WT enters as a linear (exponent 1) structural multiplier on vc and does NOT scale clearance. Rodjun 2023 simulated a single 60 kg patient (Materials and methods, Simulated dosage regimens / Sitafloxacin: 'administered orally to an inpatient in a fasted state who weighed 60 kg and was under 65 years of age'), giving V/F = 103.2 L.",
      source_name        = "weight"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = NA_integer_,
    n_simulated    = 10000L,
    n_studies      = 1L,
    age_range      = "Simulated as under 65 years. Rodjun 2023 Table 2 footnote flags the volume of distribution as the 'age < 65 years' value, so the packaged parameters apply to the under-65 stratum only. The underlying Tanigawara 2013 estimation data pooled healthy volunteers, elderly volunteers, and renally impaired patients with patients enrolled in a respiratory-tract-infection PK/PD study; Rodjun 2023 does not restate that study's demographics.",
    weight_range   = "Simulated at a single body weight of 60 kg (Materials and methods, Simulated dosage regimens / Sitafloxacin).",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not reported on disk; the Tanigawara 2013 estimation studies were conducted in Japan.",
    disease_state  = "Simulation target population is inpatients with carbapenem-resistant (CRAB), multidrug-resistant (MDR-AB) or colistin-resistant (CoR-AB) Acinetobacter baumannii infection. The underlying Tanigawara 2013 estimation data came from clinical pharmacology studies in healthy, elderly and renally impaired subjects plus a PK/PD study in patients with community-acquired respiratory tract infections.",
    dose_range     = "Oral doses from 50 mg q48h to 1500 mg q12h according to creatinine clearance, administered in the fasted state (Materials and methods, Simulated dosage regimens / Sitafloxacin). The manufacturer-recommended regimens evaluated as comparators were 50 mg q12h, 100 mg q24h and 100 mg q12h.",
    regions        = "Thailand (simulation study, Mahidol University, Bangkok); the underlying Tanigawara 2013 PK cohort was Japanese.",
    renal_function = "Simulated at creatinine clearance 90, 50, 30 and 10 mL/min.",
    notes          = "Monte Carlo simulation of 10,000 virtual subjects per dosage regimen (Crystal Ball version 2017, Decisioneering Inc.) with log-normal between-patient variability on every PK parameter (Materials and methods, Monte Carlo simulation). Rodjun 2023 does not report the size of the Tanigawara 2013 estimation cohort, so n_subjects is NA; n_simulated records the virtual cohort size instead. The PK/PD target is fAUC/MIC > 30 (Tanigawara 2013, 91 respiratory-tract-infection isolates), applied to the unbound 24-h sitafloxacin AUC. Sitafloxacin absorption is reduced by antacids, ferrous sulfate and other metal-cation-containing compounds; the model represents the fasted state without such co-medication (Discussion, limitations)."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters - Rodjun 2023 Table 2 ("Population
    # pharmacokinetic parameters of sitafloxacin (Tanigawara et al.,
    # 2013)"). The model is the one-compartment first-order-absorption
    # system of Rodjun 2023 Eq. 4, dX/dt = Ka * Xa - K * X.
    # ------------------------------------------------------------------

    lka <- log(1.67); label("First-order absorption rate constant ka in the fasted state (1/h)")  # Rodjun 2023 Table 2: ka = 1.67 h^-1 (omega^2 = 4.57)

    # Apparent clearance. Table 2 gives CL_t/F = 2.58 x CrCL with no
    # intercept term, so the whole of CL/F is renal-function driven and
    # lcl is the log of a dimensionless SLOPE rather than a clearance.
    # model() multiplies it by CrCL expressed in L/h. See the CRCL entry
    # in covariateData for the numeric evidence that the Table 2 unit
    # column ('L') is a misprint and that CrCL must be in L/h.
    lcl <- log(2.58); label("Apparent clearance slope on creatinine clearance, CL_t/F = 2.58 x CrCL (unitless; CrCL in L/h)")  # Rodjun 2023 Table 2: CL t/F = 2.58 x CrCL (omega^2 = 0.0757)

    # Apparent volume of distribution, reported per kilogram of body
    # weight; model() multiplies it by WT.
    lvc <- log(1.72); label("Apparent volume of distribution V/F (L/kg; age < 65 years)")  # Rodjun 2023 Table 2: V/F = 1.72 L/kg (omega^2 = 0.087)

    # Unbound fraction in plasma, used to convert total sitafloxacin
    # concentration to the unbound concentration that drives the paper's
    # fAUC/MIC target. Reported as a single value with no distribution,
    # unlike the colistin unbound fraction.
    fu <- fixed(0.388); label("Sitafloxacin plasma unbound fraction (unitless)")  # Rodjun 2023 Materials and methods, Pharmacokinetic model / Sitafloxacin: "An unbound fraction of sitafloxacin of 0.388 was used (Tanigawara et al., 2013)"

    # ------------------------------------------------------------------
    # Inter-individual variability. Rodjun 2023 Table 2 reports omega^2
    # and its footnote states "SD were calculated from the square root of
    # omega^2 x 100% x estimate" - i.e. sqrt(omega^2) is an ARITHMETIC CV
    # which Crystal Ball then used as the spread of a log-normal draw.
    # The equivalent log-scale variance is omega^2_log = log(CV^2 + 1) =
    # log(omega^2_reported + 1).
    # ------------------------------------------------------------------
    etalcl ~ 0.072959  # log(0.0757 + 1); Rodjun 2023 Table 2 CL t/F omega^2 = 0.0757 (arithmetic CV 27.5%)
    etalvc ~ 0.083410  # log(0.0870 + 1); Rodjun 2023 Table 2 V/F omega^2 = 0.0870 (arithmetic CV 29.5%)
    etalka ~ 1.717400  # log(4.5700 + 1); Rodjun 2023 Table 2 ka omega^2 = 4.57 (arithmetic CV 214%) - see the validation vignette; this variance is implausibly large for an absorption rate constant but is encoded exactly as printed, and it does not affect the AUC-driven PK/PD target

    # ------------------------------------------------------------------
    # Residual error. Rodjun 2023 reports none: the Monte Carlo
    # simulation is deterministic once each virtual subject's PK
    # parameters have been drawn (Materials and methods, Monte Carlo
    # simulation). propSd is therefore fixed at 0 so the packaged model
    # reproduces the paper's noise-free simulation; a user fitting real
    # data should re-estimate it.
    # ------------------------------------------------------------------
    propSd <- fixed(0); label("Proportional residual error (fraction; 0 - not reported)")  # Rodjun 2023 Materials and methods (no residual error model reported)
  })

  model({
    # 1. Derived covariate term. Table 2's slope of 2.58 is dimensionless
    #    and therefore requires creatinine clearance in L/h, while the
    #    canonical CRCL data column is in mL/min.
    crcl_lph <- CRCL * 60 / 1000

    # 2. Individual PK parameters.
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * crcl_lph   # CL_t/F = 2.58 x CrCL (L/h)
    vc <- exp(lvc + etalvc) * WT         # V/F = 1.72 L/kg x WT

    # 3. Micro-constant.
    kel <- cl / vc

    # 4. ODE system - Rodjun 2023 Eq. 4. `depot` holds Xa (the mass of
    #    absorbable drug at the extravascular site) and `central` holds X
    #    (the mass of drug in plasma). Bioavailability F is folded into
    #    the apparent CL/F and V/F, so no f(depot) term is applied.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 5. Observations. Cc is the total plasma sitafloxacin concentration
    #    and Ccu the unbound concentration that drives the paper's
    #    fAUC/MIC > 30 target.
    Cc  <- central / vc
    Ccu <- Cc * fu

    Cc ~ prop(propSd)
  })
}
