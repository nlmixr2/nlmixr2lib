KhanAsa_2020_voriconazole <- function() {
  description <- "One-compartment population pharmacokinetic model with first-order absorption and linear elimination for oral voriconazole at steady state in Thai adults with hematologic diseases (Khan-asa 2020); apparent clearance falls linearly with serum albumin below 3.2 g/dL and by 30.6% with concomitant omeprazole 40 mg/day or more, and the absorption rate constant is fixed at 1.1 per hour"
  reference <- "Khan-asa B, Punyawudho B, Singkham N, Chaivichacharn P, Karoopongse E, Montakantikul P, Chayakulkeeree M. Impact of Albumin and Omeprazole on Steady-State Population Pharmacokinetics of Voriconazole and Development of a Voriconazole Dosing Optimization Model in Thai Patients with Hematologic Diseases. Antibiotics (Basel). 2020;9(9):574. doi:10.3390/antibiotics9090574"
  vignette <- "KhanAsa_2020_voriconazole"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # All 65 patients took voriconazole orally ('A PK study of oral voriconazole
  # was conducted', Section 4.1), so the model carries a depot and CL and V are
  # apparent oral values (CL/F, V/F); bioavailability is not identifiable.
  compartmentData <- list(
    depot = list(analyte = "voriconazole", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "voriconazole", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    ALB = list(
      description = "Plasma (serum) albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "The paper reports albumin in g/dL (Table 1: mean 3.04 +/- 0.61, range 1.60-4.3 g/dL) and centres the linear clearance term on the cohort median of 3.2 g/dL (Section 2.4: 'median: 3.2 g/dL, IQR: 2.78-3.6 g/dL'; Section 4.5: continuous covariates 'were centered by their median'). The canonical column is g/L, so model() converts with alb_gdL <- ALB * 0.1 before applying the published coefficient. The term is LINEAR, CL/F x [1 + 0.249 x (ALB - 3.2)], so clearance falls by 24.9% of the typical value per 1 g/dL drop in albumin and would reach zero at about -0.8 g/dL; it is therefore well defined over the whole observed range (1.6-4.5 g/dL).",
      source_name = "ALB"
    ),
    DOSE_OMEPRAZOLE_MGD = list(
      description = "Total daily dose of concomitant omeprazole (0 when not co-administered)",
      units = "mg/day",
      type = "continuous",
      reference_category = "0 mg/day (no omeprazole); any dose below 40 mg/day is equivalent to the reference",
      notes = "The paper's indicator OME is '0 for patients receiving omeprazole 20 mg/day or less, and 1 for patients receiving a dose of omeprazole >= 40 mg/day' (Section 2.2). model() recomputes it as ome40 <- DOSE_OMEPRAZOLE_MGD >= 40. A plain presence indicator (CONMED_OMEPRAZOLE) would be WRONG here: omeprazole 20 mg/day (44.6% of patients) was explicitly tested and not retained, so only the 40 mg/day patients (26.2%) take the effect. Esomeprazole (2 patients) and rabeprazole (1 patient) were not part of OME; enter 0 for them.",
      source_name = "OME"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Allometric scaling of CL/F (exponent 0.75) and V/F (exponent 1) on median-centred weight was tested and did not reduce the OFV significantly (Section 2.2), so no weight term is in the final model."
    ),
    CYP2C19_PHENOTYPE = list(
      description = "CYP2C19 metabolizer phenotype / genotype",
      units = "categorical",
      type = "categorical",
      notes = "Genotyped for *2, *3 and *17 (1 UM, 33 EM, 24 IM, 7 PM; Table 1). CYP2C19 *2/*3 was significant only in forward selection and was removed on backward deletion; no genotype or phenotype is retained (Section 2.2)."
    ),
    CONMED_TMPSMX = list(
      description = "Concomitant sulfamethoxazole/trimethoprim",
      units = "(binary)",
      type = "binary",
      notes = "Significant only in forward selection; not retained (Section 2.2)."
    )
  )

  population <- list(
    n_subjects = 65,
    n_studies = 1,
    n_observations = 237,
    species = "human",
    age_range = "20-78 years (mean 47.65, SD 15.22)",
    weight_range = "27.2-105 kg (mean 58.57, SD 13.83)",
    sex_female_pct = 36.92,
    race_ethnicity = "Thai (100%)",
    disease_state = "Hematologic diseases (AML 56.9%, lymphoma 21.5%, ALL 13.9%, other 7.7%) receiving voriconazole for prophylaxis or treatment of invasive aspergillosis",
    dose_range = "Oral loading dose 400 mg every 12 h (2 doses), then 200 mg every 12 h (weight >= 40 kg); mean 7.18 mg/kg/day (range 3.81-14.71)",
    regions = "Thailand (Siriraj Hospital, Bangkok)",
    route = "oral",
    albumin = "3.04 +/- 0.61 g/dL (range 1.60-4.3); median 3.2 g/dL",
    co_medication = "Omeprazole 20 mg/day 44.62%, omeprazole 40 mg/day 26.15%, esomeprazole 3.08%, rabeprazole 1.54%, sulfamethoxazole/trimethoprim 29.23%",
    genotype = "CYP2C19: 1 UM, 33 EM, 24 IM, 7 PM (not retained)",
    notes = "18 patients had intensive sampling on day 7 (0, 1, 1.5, 2, 4, 8, 12 h); 47 contributed trough-type samples between day 7 and day 40. All samples at steady state, 0.75-12 h post dose. HPLC assay, LLOQ 0.2 mg/L. NONMEM VII, FOCE-I. Baseline characteristics from Table 1."
  )

  ini({
    # Structural parameters: Table 2, 'NONMEM Results / Estimated Value'.
    lka <- fixed(log(1.1)) # Table 2, row 'Ka (/h)' = 'FIX 1.100'; fixed to the value of Pascual 2012 (ref 25) because the estimate was unreliable (Section 2.2). Table S1 shows CL/F moved only 2.99-3.15 L/h across Ka 0.163-1.38 /h.
    label("Absorption rate constant (1/h)")
    lcl <- log(3.43) # Table 2, row 'CL/F (L/h)' = 3.430 (SE 0.287; bootstrap median 3.430); also the leading constant of the Section 2.2 final-model equation.
    label("Apparent oral clearance at ALB 3.2 g/dL without omeprazole >= 40 mg/day (L/h)")
    lvc <- log(47.6) # Table 2, row 'V/F (L)' = 47.600 (SE 6.600; bootstrap median 47.811).
    label("Apparent volume of distribution (L)")

    # Covariate effects: Section 2.2 final-model equation
    #   CL/F (L/h) = 3.43 x [1 + 0.249 x (ALB - 3.2)] x [1 + (-0.306 x OME)]
    e_alb_cl <- 0.249 # Table 2, row 'CL-albumin' = 0.249 (SE 0.0925; bootstrap median 0.250). Linear slope per g/dL on the median-centred albumin (3.2 g/dL).
    label("Linear effect of albumin on CL/F, per g/dL from 3.2 g/dL (unitless)")
    e_ome40_cl <- -0.306 # Table 2, row 'CL-omeprazole >= 40 mg/day' = -0.306 (SE 0.084; 95% CI -0.471 to -0.141; bootstrap median -0.300); matches the Section 2.2 equation '[1 + (-0.306 x OME)]'.
    label("Fractional change in CL/F with omeprazole >= 40 mg/day (unitless)")

    etalcl ~ 0.226 # Table 2, row 'IIV-CL' = 0.226 with '(%CV) 50.40%'; sqrt(exp(0.226) - 1) = 0.504 confirms 0.226 is the log-scale variance. No IIV on V/F or Ka (Section 2.2).

    # Additive residual error. Table 2 prints 'RUV (mg/L) = 2.67'. Read as
    # the NONMEM $SIGMA VARIANCE (mg/L)^2, so the SD is sqrt(2.67) = 1.634
    # mg/L. Reasons: (i) the same table reports IIV-CL as a variance (0.226,
    # 50.4% CV), i.e. raw $OMEGA/$SIGMA output; (ii) the VPC in Figure S1
    # discriminates: the simulated 5th-percentile band is centred near 2.7
    # mg/L at 1.5 h, ~2 mg/L at 5 h and ~0 at 12 h, which an SD of 1.63
    # reproduces (2.7, 1.8, -0.1 mg/L) while an SD of 2.67 puts it at 1.4,
    # 0.6 and -1.4 mg/L - on or below the lower edge of the band.
    addSd <- 1.634 # sqrt(2.67); Table 2, row 'RUV (mg/L)' = 2.67 (SE 0.706; bootstrap median 2.576), interpreted as the $SIGMA variance - see above.
    label("Additive residual error (mg/L)")
  })

  model({
    # Albumin: canonical column is g/L; the published slope is per g/dL.
    alb_gdL <- ALB * 0.1
    # Paper's OME indicator: 1 for omeprazole >= 40 mg/day, 0 otherwise.
    ome40 <- 0
    if (DOSE_OMEPRAZOLE_MGD >= 40) ome40 <- 1

    ka <- exp(lka)
    cl <- exp(lcl + etalcl) *
      (1 + e_alb_cl * (alb_gdL - 3.2)) *
      (1 + e_ome40_cl * ome40)
    vc <- exp(lvc)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
