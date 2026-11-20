HourcadePotelleret_2014_trastuzumab <- function() {
  description <- "Two-compartment population PK model with first-order subcutaneous absorption and parallel linear and Michaelis-Menten elimination from the central compartment for intravenous and subcutaneous (rHuPH20 co-formulated) trastuzumab in healthy male volunteers and women with HER2-positive early breast cancer (Hourcade-Potelleret 2014)"
  reference <- "Hourcade-Potelleret F, Lemenuel-Diot A, McIntyre C, Brewster M, Lum B, Bittner B. Use of a population pharmacokinetic approach for the clinical development of a fixed-dose subcutaneous formulation of trastuzumab. CPT Pharmacometrics Syst Pharmacol. 2014;3(1):e87. doi:10.1038/psp.2013.63"
  vignette <- "HourcadePotelleret_2014_trastuzumab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline body weight. Power effects on linear CL (exponent 0.591) and central volume V2 (exponent 0.48), both normalised to 68 kg and both exponents fixed to the values of the earlier intravenous trastuzumab popPK model (Hourcade-Potelleret 2014 Results 'PopPK modeling' and the supplementary NONMEM control stream: CLWT = (WT/68)**THETA(7), V2WT = (WT/68)**THETA(8)).",
      source_name = "WT"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "trastuzumab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "trastuzumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "trastuzumab", units = "mg", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 66L,
    n_studies = 1L,
    n_observations = 786L,
    sex_female_pct = 63.6,
    disease_state = "Healthy male volunteers (n = 24, part 1) and women previously treated for HER2-positive early breast cancer (n = 42, parts 1 and 2) of the phase I/Ib dose-finding / dose-confirmation trial BP22023 (NCT00800436).",
    dose_range = "Single dose of trastuzumab: 6 mg/kg IV (6 healthy men, 6 women with EBC), or 6, 8, or 10 mg/kg SC in healthy men (6 per dose) and 8 or 12 mg/kg SC in women with EBC (20 per dose). Four women received both 6 mg/kg IV and 8 mg/kg SC, giving 70 doses in 66 individuals.",
    regions = "Not stated in the source (single phase I/Ib trial; design reported in Wynne 2013, J Clin Pharmacol 53:192-201).",
    notes = "Hourcade-Potelleret 2014 Results 'Pharmacokinetic data set' and Methods 'Pharmacokinetic data sets'. Demographics are not tabulated in this paper (they are in Wynne 2013). 786 of 789 serum concentrations were analysed (three pre-dose quantifiable samples excluded); LLOQ 0.156 mg/L. Female percentage = 42/66 individuals. Simulations in the paper used a virtual population with body weight mean 68 kg, SD 11 kg."
  )

  ini({
    # Final NONMEM estimates (Hourcade-Potelleret 2014 Table 1; identical to the
    # FINAL PARAMETER ESTIMATE block of the supplementary NONMEM output, which
    # gives THETA(10) = KM as 6.12E+03 ug/L because the control stream scales
    # concentrations by S2 = V2/1000). Typical values are for a 68 kg subject.
    lcl <- log(0.187); label("Linear clearance CL for a 68 kg subject (L/day)") # Table 1 CL = 0.187 L/day; NONMEM TH1 1.87E-01
    lvc <- log(3.43); label("Central volume V2 for a 68 kg subject (L)") # Table 1 V2 = 3.43 L; NONMEM TH2 3.43E+00
    lq <- log(0.24); label("Intercompartmental clearance Q (L/day)") # Table 1 Q = 0.24 L/day; NONMEM TH3 2.40E-01
    lvp <- log(1.71); label("Peripheral volume V3 (L)") # Table 1 V3 = 1.71 L; NONMEM TH4 1.71E+00
    lfdepot <- log(0.873); label("Absolute subcutaneous bioavailability F1 (fraction)") # Table 1 F1 = 0.873; NONMEM TH5 8.73E-01
    lka <- log(0.28); label("First-order subcutaneous absorption rate constant Ka (1/day)") # Table 1 Ka = 0.28 /day; NONMEM TH6 2.80E-01
    lvmax <- log(3.49); label("Maximum Michaelis-Menten elimination rate Vmax (mg/day)") # Table 1 Vm = 3.49 mg/day; NONMEM TH9 3.49E+00
    lkm <- log(6.12); label("Michaelis-Menten constant Km (mg/L)") # Table 1 Km = 6.12 mg/L; NONMEM TH10 6.12E+03 ug/L

    # Body-weight exponents fixed to the earlier IV trastuzumab popPK model
    # (Results 'PopPK modeling'; Table 1 CVSE column reads 'Fix'; control
    # stream $THETA 0.591 FIX and 0.48 FIX).
    e_wt_cl <- fixed(0.591); label("Power exponent of body weight on CL (unitless)") # Table 1 'Influence of body weight on CL' = 0.591 Fix
    e_wt_vc <- fixed(0.48); label("Power exponent of body weight on V2 (unitless)") # Table 1 'Influence of body weight on V2' = 0.48 Fix

    # IIV: exponential (log-normal) etas; variances are the OMEGA diagonal of
    # the supplementary NONMEM output. Table 1 reports their square roots as
    # apparent CVs (22.6%, 21.2%, 16.9%, 38.6%).
    etalcl ~ 0.0512 # NONMEM OMEGA(1,1) 5.12E-02; Table 1 IIV CL 22.6 pct
    etalvc ~ 0.0449 # NONMEM OMEGA(2,2) 4.49E-02; Table 1 IIV V2 21.2 pct
    etalfdepot ~ 0.0285 # NONMEM OMEGA(4,4) 2.85E-02; Table 1 IIV F1 16.9 pct
    etalka ~ 0.149 # NONMEM OMEGA(5,5) 1.49E-01; Table 1 IIV Ka 38.6 pct

    # Proportional residual error: Y = F + F*EPS(1), SIGMA = 7.11E-02 in the
    # NONMEM output; SD = sqrt(0.0711) = 0.2667. Table 1 prints the square root
    # as 0.26 (truncated).
    propSd <- 0.2667; label("Proportional residual error (fraction)") # NONMEM SIGMA(1,1) 7.11E-02, sqrt = 0.2667; Table 1 0.26
  })

  model({
    # Individual parameters (supplementary NONMEM $PK):
    #   CL = THETA(1) * (WT/68)**THETA(7) * EXP(ETA(1))
    #   V2 = THETA(2) * (WT/68)**THETA(8) * EXP(ETA(2))
    #   F1 = THETA(5) * EXP(ETA(4));  KA = THETA(6) * EXP(ETA(5))
    cl <- exp(lcl + etalcl) * (WT / 68)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 68)^e_wt_vc
    q <- exp(lq)
    vp <- exp(lvp)
    ka <- exp(lka + etalka)
    vmax <- exp(lvmax)
    km <- exp(lkm)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    Cc <- central / vc

    # Supplementary NONMEM $DES (ADVAN6): DADT(2) = KA*A(1) - K12*A(2) +
    # K21*A(3) - K*A(2) - C2*VM/(KM+C2). Amounts in mg, volumes in L, so Cc is
    # mg/L and both elimination fluxes are mg/day.
    d / dt(depot) <- -ka * depot
    d / dt(central) <- ka * depot - kel * central - vmax * Cc / (km + Cc) - k12 * central + k21 * peripheral1
    d / dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Bioavailability applies to the SC depot only (NONMEM F1 on compartment 1);
    # IV doses go directly into central. F1 is log-normal as in the control
    # stream, so individual values above 1 are possible.
    f(depot) <- exp(lfdepot + etalfdepot)

    Cc ~ prop(propSd)
  })
}
