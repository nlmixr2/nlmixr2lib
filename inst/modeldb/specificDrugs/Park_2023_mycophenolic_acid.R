Park_2023_mycophenolic_acid <- function() {
  description <- "Population PK model for total mycophenolic acid (MPA, the active moiety of mycophenolate mofetil MMF) in pediatric allogeneic haematopoietic stem cell transplantation (HSCT) recipients receiving oral MMF for acute graft-versus-host disease prophylaxis or treatment (Park 2023). One-compartment disposition with first-order absorption and no lag time (NONMEM ADVAN2 TRANS2). Body surface area is the only retained covariate, entering the apparent volume of distribution as a centered-linear term Vd/F = 89.83 * (1 + 0.854 * (BSA - 1.11)). Inter-individual variability is log-normal on CL/F and Vd/F; residual error is combined proportional plus additive. Doses are MMF mass (mg) with no MMF-to-MPA molecular-weight conversion -- CL/F and Vd/F are apparent parameters that absorb both the molecular-weight ratio and oral bioavailability."
  reference <- paste(
    "Park HJ, Hong KT, Han N, Kim I-W, Oh JM, Kang HJ.",
    "Body Surface Area-Based Dosing of Mycophenolate Mofetil in Pediatric",
    "Hematopoietic Stem Cell Transplant Recipients: A Prospective Population",
    "Pharmacokinetic Study.",
    "Pharmaceutics. 2023;15(12):2741.",
    "doi:10.3390/pharmaceutics15122741.",
    sep = " "
  )
  vignette <- "Park_2023_mycophenolic_acid"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    depot   = list(analyte = "mycophenolic acid", units = "mg", specimen = "administration site", verified = FALSE),
    central = list(analyte = "mycophenolic acid", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    BSA = list(
      description        = "Body surface area on the day of blood sampling.",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Centered-linear effect on the apparent volume of distribution: Vd/F = 89.83 * (1 + 0.854 * (BSA - 1.11)), per Park 2023 Equation 2. The centering value 1.11 m^2 is described in the text immediately below Equation 2 as 'the median value of BSA for patients in this study'; Park 2023 Table 1 independently reports the cohort median BSA as 1.12 m^2. The 1.11 m^2 value from the equation is used here because the covariate model was estimated against it (see vignette Errata). Park 2023 does not state which BSA formula (DuBois, Mosteller, Haycock) was used to derive the covariate; the register entry requires this to be recorded as 'unspecified' when the source is silent. Cohort median 1.12 m^2, range 0.49-1.60 m^2 (Table 1). The simulation in Park 2023 Figure 5 spans BSA 0.5, 1.0 and 1.5 m^2. Note this is a centered-linear parameterization rather than the power scaling (BSA/ref)^exponent documented as the common form in the BSA register entry; the linear form can go negative for BSA below 1.11 - 1/0.854 = -0.06 m^2, which is outside any physiological range, so no clamping is required.",
      source_name        = "BSA"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 20L,
    n_studies      = 1L,
    n_observations = 80L,
    age_range      = "1.7-15.6 years",
    age_median     = "9.7 years",
    age_group      = c(age_1_to_11_pct = 70.0, age_12_to_17_pct = 30.0),
    weight_range   = "9.9-51.0 kg",
    weight_median  = "31.2 kg",
    height_range   = "83.6-176.9 cm",
    height_median  = "136.0 cm",
    bsa_range      = "0.49-1.60 m^2",
    bsa_median     = "1.12 m^2 (Table 1); the covariate model in Equation 2 is centered on 1.11 m^2",
    sex_female_pct = 40.0,
    race_ethnicity = "Not reported. Single-center Korean study (Seoul National University Hospital); the cohort is presumed predominantly Korean but Park 2023 does not tabulate race or ethnicity.",
    disease_state  = "Pediatric recipients of allogeneic haematopoietic stem cell transplantation started on mycophenolate mofetil for prophylaxis (16 patients, 80.0%) or treatment (4 patients, 20.0%) of acute graft-versus-host disease. Underlying diagnoses: acute lymphoblastic leukemia (7, 35.0%), acute myeloid leukemia (5, 25.0%), aplastic anemia (2, 10.0%), non-Hodgkin lymphoma (2, 10.0%), congenital neutropenia (1, 5.0%), haemophagocytic lymphohistiocytosis (1, 5.0%), Krabbe disease (1, 5.0%), therapy-related myelodysplastic syndrome (1, 5.0%). Conditioning: BuFluCy (15, 75.0%), TBI-FluCyATG (2, 10.0%), BuFluATG / BuFludaVPATG / FluMelATG (1 each, 5.0%). Donor source haploidentical family (16, 80.0%) or matched unrelated (4, 20.0%). ABO compatible (16, 80.0%) or incompatible (4, 20.0%). Renal and hepatic function were essentially normal: SCr median 0.42 mg/dL (0.32-0.72), eGFR median 112.6 mL/min/1.73 m^2 (86.3-165.3), total bilirubin median 0.5 mg/dL (0.4-2.0), albumin median 3.8 g/dL (3.3-4.2). Acute GVHD developed in 11 patients (55.0%), 4 with grade III-IV.",
    dose_range     = "Oral mycophenolate mofetil 15-20 mg/kg twice daily as capsule (6 patients, 30.0%) or oral suspension (14 patients, 70.0%). Median total daily dose 1100 mg/day (range 380-2000); median 17.9 mg/kg/dose (range 16.1-19.8). Doses in this model file are MMF mass in mg with no MMF-to-MPA molecular-weight conversion applied.",
    regions        = "Single center: Seoul National University Hospital, Seoul, Republic of Korea.",
    co_medication  = "All 20 patients (100.0%) had previous (4, 20.0%) or concomitant (16, 80.0%) tacrolimus on the day of blood sampling. Ciprofloxacin 15 (75.0%), esomeprazole 4 (20.0%), famotidine 4 (20.0%), lansoprazole 3 (15.0%), itraconazole 2 (10.0%), voriconazole 1 (5.0%). Proton-pump inhibitor use (35%) is noted in the Discussion as a plausible contributor to low MPA exposure but was not retained as a covariate.",
    notes          = "Prospective single-center study conducted 1 September 2020 - 30 June 2022 (IRB No. 2006-120-1133). Sampling was a limited-sampling strategy at pre-dose (0 h) and 1, 2 and 6 h post-dose, drawn at least 3 days after MMF initiation so that MPA had reached steady state; median time post-HSCT 31 days (20-181), median MMF duration 23 days (16-123). Serum total MPA was assayed by particle-enhanced turbidimetric inhibition immunoassay (PETINIA) on a Dimension EXL 200 with a lower limit of quantification of 0.1 mg/L. Estimation was FOCE-I in NONMEM 7.5.0. Enterohepatic recirculation of MPA-7-O-glucuronide back to MPA was investigated during model building but minimization terminated due to model instability, so it is absent from the final model (Park 2023 Discussion, limitations). Genetic polymorphisms (e.g. UGT2B7) were not screened. Baseline demographics per Park 2023 Table 1; final-model parameter estimates per Park 2023 Table 4 and Equation 2."
  )

  ini({
    # Structural parameters. Park 2023 Table 4, "Final Model" column.
    # Absorption is first order with no lag time (ADVAN2 TRANS2; the
    # Methods state lag-time, Erlang and transit-compartment absorption
    # models were all evaluated and the plain first-order model was
    # selected). No IIV was retained on ka.
    lka <- log(5.18);   label("Absorption rate constant (1/h)")                    # Park 2023 Table 4 final model, ka = 5.18 h^-1 (RSE 21%); bootstrap median 5.22 (95% CI 1.83-7.04)

    # Vd/F typical value is the value at the covariate reference
    # BSA = 1.11 m^2, at which the centered-linear factor equals 1.
    # Equation 2 prints 89.83 L; Table 4 prints the same estimate rounded
    # to 89.8 L. The 4-significant-figure value from the equation is used.
    lvc <- log(89.83);  label("Apparent volume of distribution at BSA = 1.11 m^2 (L)")  # Park 2023 Equation 2 (89.83 L); Table 4 final model rounds to 89.8 L (RSE 16%); bootstrap median 85.13 (95% CI 60.65-121.45)

    # No covariate was retained on CL/F -- BSA was the only significant
    # covariate in the stepwise search and it entered Vd/F alone.
    lcl <- log(16.6);   label("Apparent clearance (L/h)")                          # Park 2023 Table 4 final model, CL/F = 16.6 L/h (RSE 17%); bootstrap median 17.14 (95% CI 12.03-23.60)

    # Covariate effect: centered-linear BSA on Vd/F.
    # Vd/F = 89.83 * (1 + 0.854 * (BSA - 1.11)). Equation 2 prints the
    # slope rounded to 0.85; Table 4 gives 0.854 with RSE 24%. The
    # 3-significant-figure table value is used.
    e_bsa_vc <- 0.854;  label("Linear coefficient for centered BSA on Vd/F (per m^2)")  # Park 2023 Table 4 final model, "BSA on Vd/F" = 0.854 (RSE 24%); bootstrap median 0.839 (95% CI 0.32-1.20); form per Equation 2

    # Inter-individual variability. Park 2023 reports IIV as a percent CV
    # and gives the conversion explicitly in the Table 4 footnote:
    #   %CV = sqrt(exp(OMEGA) - 1) * 100
    # so the tabulated %CV is the exact log-normal CV and OMEGA is the
    # variance of eta. Inverting: omega^2 = log((%CV/100)^2 + 1).
    #   IIV Vd/F 37.71% -> log(0.3771^2 + 1) = 0.13296
    #   IIV CL/F 89.96% -> log(0.8996^2 + 1) = 0.59293
    # This footnote removes the usual CV%-vs-omega ambiguity: no
    # approximation is needed. No IIV was reported on ka.
    etalvc ~ 0.13296    # Park 2023 Table 4 final model, IIV Vd/F = 37.71 %CV (RSE 69%); omega^2 = log(0.3771^2 + 1) per the Table 4 footnote conversion
    etalcl ~ 0.59293    # Park 2023 Table 4 final model, IIV CL/F = 89.96 %CV (RSE 33%); omega^2 = log(0.8996^2 + 1) per the Table 4 footnote conversion

    # Residual error: combined proportional plus additive (Park 2023
    # Results 3.2.2, "A combined residual error model was selected").
    #
    # SCALE AMBIGUITY: Park 2023 Table 4 lists "Proportional error 0.660"
    # and "Additive error 0.111" as bare decimals with no units and no
    # conversion footnote, unlike the IIV rows which are explicitly
    # labelled %CV. They are therefore either (A) standard deviations --
    # 66.0% proportional CV and 0.111 mg/L additive SD -- or (B) raw
    # NONMEM $SIGMA variances -- 81.2% proportional CV and 0.333 mg/L
    # additive SD. Reading (A) is encoded here on four grounds, all
    # pointing the same way:
    #   1. An additive SD of 0.111 mg/L sits essentially on the assay
    #      lower limit of quantification of 0.1 mg/L (Methods 2.2), the
    #      textbook magnitude for an additive term; under reading (B) it
    #      would be 3.3x the LLOQ.
    #   2. Simulating the study design and applying the paper's own
    #      limited-sampling NCA (0, 1, 2, 6 h) gives a mean Cmax of about
    #      9.0 mg/L under (A) and about 9.7 mg/L under (B), against the
    #      reported 8.50 mg/L (Table 2) -- roughly +6% versus +14%.
    #   3. The same simulation gives a mean observed Tmax of about 1.70 h
    #      under (A) and about 1.77 h under (B), against the 1.37 h implied
    #      by Table 3 (1.29 h in 7 aGVHD and 1.44 h in 9 non-aGVHD
    #      patients) -- roughly +24% versus +29%. Observed Tmax is a direct
    #      probe of residual-error
    #      magnitude, because residual noise is what moves a subject's
    #      observed maximum off the 1 h sample.
    #   4. Under reading (B) the additive SD of 0.333 mg/L drives the
    #      predicted 5th-percentile concentration band below zero at 6 h,
    #      which is not what Park 2023 Figure 4 (VPC) shows.
    # No single test is decisive -- both readings overshoot Cmax and Tmax,
    # so the data would favour an even smaller residual error than (A) --
    # but (A) is closer on every test. The choice is flagged in the
    # vignette Errata. Under reading (B) replace with
    # propSd <- sqrt(0.660) and addSd <- sqrt(0.111).
    propSd <- 0.660;    label("Proportional residual error (fraction)")            # Park 2023 Table 4 final model, "Proportional error" 0.660 (RSE 9%); bootstrap median 0.656 (95% CI 0.547-0.763); read as an SD -- see comment above
    addSd  <- 0.111;    label("Additive residual error (mg/L)")                    # Park 2023 Table 4 final model, "Additive error" 0.111 (RSE 17%); bootstrap median 0.109 (95% CI 0.083-0.327); read as an SD -- see comment above
  })

  model({
    # Absorption rate constant -- no IIV retained (Park 2023 Table 4
    # reports IIV on Vd/F and CL/F only).
    ka <- exp(lka)

    # Apparent clearance, log-normal IIV, no covariate.
    cl <- exp(lcl + etalcl)

    # Centered-linear BSA effect on the apparent volume of distribution
    # (Park 2023 Equation 2). The reference BSA of 1.11 m^2 is the study
    # median quoted in the text below the equation.
    bsa_factor <- 1 + e_bsa_vc * (BSA - 1.11)
    vc <- exp(lvc + etalvc) * bsa_factor

    kel <- cl / vc

    # ODE system: depot -> central, first-order in and out. NONMEM
    # ADVAN2 TRANS2.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Bioavailability is not separately identifiable from an oral-only
    # dataset, so F is absorbed into the apparent parameters CL/F and
    # Vd/F and no f(depot) term is applied. Doses are MMF mass in mg;
    # the MMF-to-MPA molecular-weight ratio is likewise absorbed into the
    # apparent parameters (Park 2023 applies no molecular-weight
    # conversion -- see vignette Errata).
    #
    # Observation: dose in mg, vc in L -> central/vc in mg/L, matching the
    # serum total MPA concentration units used throughout Park 2023.
    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
