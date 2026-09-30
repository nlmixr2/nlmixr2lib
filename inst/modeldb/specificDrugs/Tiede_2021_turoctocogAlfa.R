Tiede_2021_turoctocogAlfa <- function() {
  description <- "One-compartment population PK model for factor VIII activity (FVIII:C, IU/dL = %) after intravenous turoctocog alfa (NovoEight, B-domain-truncated recombinant FVIII) in previously treated children, adolescents and adults with severe hemophilia A without inhibitors (guardian 1, 2 and 3 trials; Tiede 2021, parameters from Jimenez-Yuste 2015). Clearance and volume scale allometrically with body weight (reference 70 kg; estimated exponents 0.95 on CL and 0.86 on V), and clearance decreases linearly with age (-1% per year, reference 20 years). Log-normal uncorrelated IIV on CL and V. The combined proportional + additive residual error magnitudes are not reported and are fixed to 0. The paper's exposure-response analysis (negative-binomial annualized bleeding rate per predicted FVIII:C category) is a statistical regression and is reproduced in the vignette, not in the model."
  reference <- "Tiede A, Abdul Karim F, Jimenez-Yuste V, Klamroth R, Lejniece S, Suzuki T, Groth A, Santagostino E. Factor VIII activity and bleeding risk during prophylaxis for severe hemophilia A: a population pharmacokinetic model. Haematologica. 2021;106(7):1902-1909. doi:10.3324/haematol.2019.241554. Population PK parameters (Table 2) cite: Jimenez-Yuste V, Lejniece S, Klamroth R, et al. The pharmacokinetics of a B-domain truncated recombinant factor VIII, turoctocog alfa (NovoEight), in patients with hemophilia A. J Thromb Haemost. 2015;13(3):370-379. doi:10.1111/jth.12816."
  vignette <- "Tiede_2021_turoctocogAlfa"
  units <- list(time = "h", dosing = "IU", concentration = "IU/dL")

  compartmentData <- list(
    central = list(analyte = "turoctocog alfa", units = "IU", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric power scaling of CL (exponent 0.95) and V (exponent 0.86) with reference weight Wref = 70 kg (Tiede 2021 Results equations for V(W) and CL(W, A); Table 2). Analysis-population mean (SD) 73.5 (18.13) kg in adults/adolescents and 24.6 (10.03) kg in children (Table 1).",
      source_name = "W"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear effect on CL, (1 + k_age * (AGE - Aref)) with Aref = 20 years and k_age = -0.01 per year (Tiede 2021 Results equation for CL(W, A); Table 2). Table 1 reports age on entering the trial programme: mean (SD) 28.98 (12.15) years in adults/adolescents and 6.08 (2.91) years in children.",
      source_name = "A"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 231L,
    n_studies = 3L,
    age_range = "0-<12 years (children) and >=12 years (adults/adolescents); range not reported",
    age_mean = "adults/adolescents 28.98 (SD 12.15) years; children 6.08 (SD 2.91) years (Table 1, all patients)",
    weight_mean = "adults/adolescents 73.5 (SD 18.13) kg; children 24.6 (SD 10.03) kg (Table 1, all patients)",
    sex_female_pct = 0,
    disease_state = "Severe hemophilia A (FVIII <= 1%) without inhibitors, previously treated (>= 150 exposure days for adults/adolescents, >= 50 for children).",
    dose_range = "Intravenous turoctocog alfa prophylaxis 20-50 IU/kg every second day or 20-60 IU/kg three times weekly depending on age; PK sessions used a single 50 IU/kg dose with sampling to 48 h.",
    regions = "Multinational (guardian 1: NCT00840086; guardian 3: NCT01138501; guardian 2 extension: NCT00984126).",
    notes = "Tiede 2021 Table 1: 168 adults/adolescents (22 in the PK subgroup) and 63 children (28 in the PK subgroup). The PK data pool was the rich PK profiles of the pivotal-trial PK subgroups plus post-dose FVIII:C from routine visits for all 231 patients (Results). The abstract reports n = 187 patients for the exposure-response analysis, which conflicts with the 231 in Results and Table 1. FVIII:C by one-stage clot assay at a central laboratory. Sex is not reported; hemophilia A patients in these trials were male (0% female assumed)."
  )

  ini({
    # Reference individual: 70 kg, 20 years. Volume encoded in dL (10x the L
    # value in Table 2) and clearance in dL/h (mL/h / 100) so that
    # central / vc is FVIII:C in IU/dL (= %).
    lcl <- log(3.02) ; label("Clearance CL for a 70-kg, 20-year-old patient (dL/h)")   # Table 2: CL 70 kg, 20 y = 302 mL/h (= 3.02 dL/h)
    lvc <- log(34.6) ; label("Volume of distribution V for a 70-kg patient (dL)")      # Table 2: V 70 kg = 3.46 L (= 34.6 dL)

    e_wt_cl <- 0.95 ; label("Allometric exponent of (WT/70) on CL (unitless)")          # Table 2: eCL = 0.95
    e_wt_vc <- 0.86 ; label("Allometric exponent of (WT/70) on V (unitless)")           # Table 2: eV = 0.86
    e_age_cl <- -0.01 ; label("Linear age effect on CL, fraction per year from 20 years (1/year)") # Table 2: Age effect on CL = -0.01 1/year

    # Table 2 reports inter-individual variability as a CV (0.32 CL, 0.22 V);
    # log-normal, uncorrelated (Methods). omega^2 = log(1 + CV^2).
    etalcl ~ 0.09749   # Table 2: CL CV 0.32; omega^2 = log(1 + 0.32^2)
    etalvc ~ 0.04727   # Table 2: V CV 0.22; omega^2 = log(1 + 0.22^2)

    # Methods: 'combined proportional and additive residual error model';
    # magnitudes are not reported in Tiede 2021, so both are fixed to 0.
    propSd <- fixed(0) ; label("Proportional residual error (fraction; magnitude not reported)")   # Methods: combined proportional + additive; value not reported
    addSd <- fixed(0) ; label("Additive residual error (IU/dL; magnitude not reported)")          # Methods: combined proportional + additive; value not reported
  })

  model({
    # Results equations: V(W) = (W/Wref)^eV and
    # CL(W, A) = (W/Wref)^eCL * (1 + k_age * (A - Aref)), multiplying the
    # Table 2 reference values; Wref = 70 kg, Aref = 20 years.
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * (1 + e_age_cl * (AGE - 20))
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc

    kel <- cl / vc

    # One-compartment IV model, Cp(t) = D/V * exp(-CL/V * t) (Results).
    d/dt(central) <- -kel * central

    Cc <- central / vc

    Cc ~ prop(propSd) + add(addSd)
  })
}
