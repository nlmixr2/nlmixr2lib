Daud_2011_peginterferon_alfa_2b <- function() {
  description <- paste(
    "Integrated PK/PD model for subcutaneous high-dose peginterferon alfa-2b",
    "(PEG-Intron) and the absolute neutrophil count in adults with resected",
    "high-risk (stage IIB-IIIC) melanoma (Daud 2011). PK: one-compartment",
    "model with first-order absorption and first-order elimination, volume",
    "proportional to body weight. PD: indirect-response model in which serum",
    "peginterferon alfa-2b fully inhibits the zero-order production (kin) of",
    "the ANC response, expressed as percent of baseline, with first-order",
    "loss (kout). Parameters are geometric means of individual WinNonlin",
    "fits; IIV is the between-individual CV of those individual estimates."
  )
  reference <- paste(
    "Daud AI, Xu C, Hwu WJ, Urbas P, Andrews S, Papadopoulos NE, Floren LC,",
    "Yver A, DeConti RC, Sondak VK.",
    "Pharmacokinetic/pharmacodynamic analysis of adjuvant pegylated",
    "interferon alpha-2b in patients with resected high-risk melanoma.",
    "Cancer Chemother Pharmacol. 2011;67(3):657-666.",
    "doi:10.1007/s00280-010-1326-9"
  )
  vignette <- "Daud_2011_peginterferon_alfa_2b"
  units <- list(
    time = "h",
    dosing = "ug",
    concentration = "pg/mL serum peginterferon alfa-2b (Cc); ANC as percent of baseline (ANC)"
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Daud 2011 dosed per kilogram (6 ug/kg/week induction, 3 ug/kg/week maintenance) and reported V/F per kilogram (1.08 L/kg, supplementary Table 1). The model therefore scales V/F linearly with WT (exponent 1, implied by the per-kg normalisation, not estimated); ka and kel are not weight-scaled, so CL/F = kel * V/F is also proportional to WT. Supply the dose as the total ug (ug/kg x WT). Table 1 median 82 kg (range 46-137).",
      source_name = "Weight"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "peginterferon alfa-2b", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "peginterferon alfa-2b", units = "ug", specimen = "serum", verified = TRUE),
    anc = list(
      analyte = "absolute neutrophil count, percent of baseline",
      units = "% of baseline",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 32L,
    n_studies = 1L,
    age_range = "20-77 years (mean 49, median 50; Daud 2011 Table 1)",
    weight_range = "46-137 kg (mean 84, median 82; one missing; Daud 2011 Table 1)",
    sex_female_pct = 31,
    race_ethnicity = c(White = 97, Black = 3),
    disease_state = "Adults with histologically documented, completely resected AJCC stage IIB-IIIC melanoma (IIB 6%, IIIA 16%, IIIB 59%, IIIC 19%), ECOG 0-1, lymphadenectomy within 90 days, no prior interferon, chemotherapy or immunotherapy.",
    dose_range = "Peginterferon alfa-2b 6 ug/kg SC once weekly for 8 weeks (induction) then 3 ug/kg SC once weekly (maintenance); PK phase = first 12 weeks.",
    regions = "Four sites (phase II, open-label, single-arm; NCT00457418).",
    notes = "Compartmental PK parameters are geometric means of individual one-compartment WinNonlin fits (supplementary Table 1, N = 32); ANC PD parameters are geometric means of individual fits in 31 patients (one discontinued without sufficient PD data). The noncompartmental analysis (Table 3) used the 20 patients with full doses and complete sampling. PK sampled predose and at 24, 48, 96 and 168 h in weeks 1, 8 and 12, plus predose troughs at weeks 10 and 11; ANC approximately weekly."
  )

  ini({
    # PK -- one-compartment, first-order absorption and elimination
    # (Daud 2011 supplementary Methods Scheme 1; supplementary Table 1,
    # 'Estimate of geometric mean', N = 32)
    lka <- log(0.161); label("First-order SC absorption rate constant K01 (1/h)") # Supp Table 1: K01 = 0.161 h-1 (CV 83%)
    lkel <- log(0.0137); label("First-order elimination rate constant K10 (1/h)") # Supp Table 1: K10 = 0.0137 h-1 (CV 35%)
    lvc <- log(1.08); label("Apparent volume of distribution V/F per kg body weight (L/kg)") # Supp Table 1: V/F = 1.08 l/kg (CV 38%)

    # PD -- indirect response, inhibition of production (Daud 2011
    # supplementary Methods Scheme 2 and supplementary Table 1 equation
    # dR/dt = Kin*[1 - Cp/(Cp + IC50)] - Kout*R; 'PD parameters of the ANC
    # PK/PD model', N = 31). Rate constants are reported per day; model()
    # converts them to per hour.
    lkin <- log(57.9); label("Zero-order production rate of the ANC response kin (% of baseline per day)") # Supp Table 1 ANC PK/PD: Kin = 57.9 day-1 (CV 145%)
    lkout <- log(0.526); label("First-order loss rate constant of the ANC response kout (1/day)") # Supp Table 1 ANC PK/PD: Kout = 0.526 day-1 (CV 150%)
    lic50 <- log(453); label("Serum concentration giving half-maximal inhibition of ANC production IC50 (pg/mL)") # Supp Table 1 ANC PK/PD: IC50 = 453 pg/mL (CV 84%)

    # IIV -- the %CV columns of supplementary Table 1 are the between-patient
    # variability of the individual WinNonlin estimates around their
    # geometric means; encoded log-normally as omega^2 = log(1 + CV^2).
    # The covariance between individual estimates was not reported, so all
    # etas are independent (see the vignette's Assumptions and deviations).
    etalka ~ 0.5241 # Supp Table 1: K01 CV 83% -> log(1 + 0.83^2)
    etalkel ~ 0.1156 # Supp Table 1: K10 CV 35% -> log(1 + 0.35^2)
    etalvc ~ 0.1349 # Supp Table 1: V/F CV 38% -> log(1 + 0.38^2)
    etalkin ~ 1.1322 # Supp Table 1: Kin CV 145% -> log(1 + 1.45^2)
    etalkout ~ 1.1787 # Supp Table 1: Kout CV 150% -> log(1 + 1.50^2)
    etalic50 ~ 0.5339 # Supp Table 1: IC50 CV 84% -> log(1 + 0.84^2)

    # Residual error -- not reported (individual WinNonlin fits); fixed to
    # zero so the model simulates individual predictions.
    propSd <- fixed(0); label("Proportional residual error on Cc (fraction); not reported") # not reported in Daud 2011
    addSd_ANC <- fixed(0); label("Additive residual error on ANC (% of baseline); not reported") # not reported in Daud 2011
  })

  model({
    # PK (time in hours)
    ka <- exp(lka + etalka)
    kel <- exp(lkel + etalkel)
    vc <- exp(lvc + etalvc) * WT

    # PD rate constants: per day in ini(), per hour here
    kin <- exp(lkin + etalkin) / 24
    kout <- exp(lkout + etalkout) / 24
    ic50 <- exp(lic50 + etalic50)

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # central in ug, vc in L -> ug/L = ng/mL; x 1000 for pg/mL
    Cc <- central / vc * 1000

    # Indirect response: drug inhibits production of the ANC response
    # completely at saturating concentration (no Imax parameter in the
    # supplementary equation). R(0) is the pre-dose steady state kin/kout.
    d/dt(anc) <- kin * (1 - Cc / (Cc + ic50)) - kout * anc
    anc(0) <- kin / kout

    ANC <- anc

    Cc ~ prop(propSd)
    ANC ~ add(addSd_ANC)
  })
}
