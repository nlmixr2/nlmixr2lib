dePillis_2025_covid19Vaccine <- function() {
  description <- paste(
    "Nonlinear mixed-effects ODE model of SARS-CoV-2 neutralizing antibody",
    "(NAb) activity after a two-dose mRNA COVID-19 vaccine primary series",
    "(BNT162b2 or mRNA-1273) in COVID-19-naive and convalesced adults.",
    "Two states: NAb activity (antibody, fraction neutralization 0-1) and a",
    "proxy for vaccine-transfected cells (depot_antigen, mL). Each vaccine",
    "dose is a unit pulse (amt = 1 over 1 day) into depot_antigen, which is",
    "eliminated by Michaelis-Menten kinetics (vmax = 10 mL/day, km = 50 mL,",
    "both fixed). NAb rises by vaccine-driven initiation (r1 * V) and",
    "antibody-dependent boosting (r2 * A * V) and otherwise follows logistic",
    "dynamics A * (r3 - r4 * A). The four rate constants and the initial NAb",
    "level are logit-normal with correlated IIV on r1-r4; prior COVID-19",
    "infection (convalesced) shifts all five and the BNT162b2 (Pfizer) vaccine",
    "shifts r1-r3. Logit-normal residual error with Monolix bounds",
    "(-0.00954434, 1.01970) and a combined1 error on the logit scale.",
    sep = " "
  )
  reference <- paste(
    "de Pillis L, Caffrey R, Chen G, Dela M, Eldevik L, Liu H, McConnell JP,",
    "Shabahang S, Varvel SA. (2025). Comparing neutralizing antibody activity",
    "over time between naive and convalesced COVID-19 vaccinated individuals.",
    "Sci Rep 15:34800. doi:10.1038/s41598-025-18673-9.",
    "Correction (author-name only): doi:10.1038/s41598-025-26404-3.",
    "Monolix project and data deposited at",
    "https://github.com/depillis/Covid_AntibodyResponse_Analysis.",
    sep = " "
  )
  vignette <- "dePillis_2025_covid19Vaccine"
  units <- list(
    time = "day",
    dosing = "dimensionless (unit vaccine pulse)",
    concentration = "fraction neutralization (0-1)"
  )

  paper_specific_compartments <- c("depot_antigen", "antibody")

  compartmentData <- list(
    depot_antigen = list(
      analyte = "vaccine-transfected cells (proxy for mRNA vaccine antigen exposure)",
      units = "mL",
      specimen = "administration site",
      verified = TRUE
    ),
    antibody = list(
      analyte = "SARS-CoV-2 neutralizing antibody activity",
      units = "fraction neutralization (0-1)",
      specimen = "serum",
      verified = TRUE
    )
  )

  covariateData <- list(
    SARS_SEROPOS = list(
      description = "Convalesced status: 1 = history of COVID-19 infection before the first vaccine dose (convalesced), 0 = COVID-19-naive",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "The paper's indicator C (data column isConvalesced = Yes / No).",
        "Convalesced status was taken from the subject's COVID-19 history;",
        "subjects without a known infection but with nucleocapsid antibodies",
        "in their pre-vaccine sample were excluded, so in this vaccine-naive",
        "population 1 means seropositive from prior infection and 0 means",
        "nucleocapsid-seronegative at baseline."
      ),
      source_name = "isConvalesced (C in Eq. 8)"
    ),
    TRT_TOZINAMERAN = list(
      description = "Primary-series vaccine product: 1 = BNT162b2 (tozinameran, Pfizer-BioNTech), 0 = mRNA-1273 (elasomeran, Moderna)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "The paper's indicator D (data column dose1Type = Pfizer / Moderna). Reference is mRNA-1273; every subject received one of the two mRNA vaccines.",
      source_name = "dose1Type (D in Eq. 8)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 92L,
    n_studies = 1L,
    age_range = "adults; mean (SD) age 53.6 (13.3) years naive, 53.7 (10.5) years convalesced",
    weight_range = "mean (SD) BMI 27.4 (5.8) kg/m^2 naive, 29.5 (4.2) kg/m^2 convalesced",
    sex_female_pct = 47.8,
    race_ethnicity = NA_character_,
    disease_state = "Adults vaccinated against COVID-19: 68 COVID-19-naive and 24 convalesced (prior COVID-19 infection; mean 110.7 days from symptom onset to first vaccine dose)",
    dose_range = "Two-dose mRNA primary series: mRNA-1273 (Moderna; 27 naive, 14 convalesced) or BNT162b2 (Pfizer-BioNTech; 41 naive, 10 convalesced)",
    regions = "United States",
    notes = paste(
      "Aditxt prospective cohort NCT05379478, samples November 2020 to April",
      "2022; 324 NAb observations (244 naive, 80 convalesced), from the first",
      "vaccine dose until the first booster or breakthrough infection.",
      "Table 1 of de Pillis 2025. Percent female is pooled from Table 1",
      "(52.9% of 68 naive and 33.3% of 24 convalesced)."
    )
  )

  ini({
    # Fixed vaccine-cell kinetics. The paper states k1 and k2 are fixed for
    # the whole population; the values are in the deposited Monolix model
    # (4ParamModel.txt: 'K1 = 10; K2 = 50', carried from the earlier
    # de Pillis 2023 J Theor Biol model).
    lvmax <- fixed(log(10))
    label("Maximum elimination rate of vaccine-transfected cells, paper k1 (mL/day)") # Eq. 6 k1; value from deposited 4ParamModel.txt
    lkm <- fixed(log(50))
    label("Transfected-cell amount at half-maximal elimination, paper k2 (mL)") # Eq. 6 k2; value from deposited 4ParamModel.txt

    # Typical values for a COVID-19-naive subject given mRNA-1273, held on the
    # logit scale (Monolix logitNormal individual parameters, Eq. 8).
    logitr1 <- logit(0.089)
    label("logit of r1, vaccine-driven NAb initiation rate (1/(mL*day))") # Table 3 r10 = 0.089 (SE 0.015)
    logitr2 <- logit(0.028)
    label("logit of r2, antibody-dependent vaccine boost rate (1/(mL*day))") # Table 3 r20 = 0.028 (SE 0.041)
    logitr3 <- logit(0.002)
    label("logit of r3, intrinsic logistic NAb growth rate (1/day)") # Table 3 r30 = 0.002 (SE 0.001)
    logitr4 <- logit(0.010)
    label("logit of r4, logistic NAb self-limitation rate (1/day)") # Table 3 r40 = 0.010 (SE 0.002)
    logitnab0 <- logit(0.00003)
    label("logit of A0, NAb activity at the first vaccine dose (fraction)") # Table 3 A00 = 0.00003 (SE 0.002)

    # Covariate effects, additive on the logit scale (Eq. 8)
    e_sars_seropos_r1 <- -3.670
    label("Convalesced effect on logit(r1) (unitless)") # Table 3 beta11 = -3.670 (SE 1.100)
    e_trt_tozinameran_r1 <- -0.905
    label("BNT162b2 vs mRNA-1273 effect on logit(r1) (unitless)") # Table 3 beta12 = -0.905 (SE 0.236)
    e_sars_seropos_r2 <- -5.860
    label("Convalesced effect on logit(r2) (unitless)") # Table 3 beta21 = -5.860 (SE 8.070)
    e_trt_tozinameran_r2 <- 1.760
    label("BNT162b2 vs mRNA-1273 effect on logit(r2) (unitless)") # Table 3 beta22 = 1.760 (SE 1.490)
    e_sars_seropos_r3 <- -1.440
    label("Convalesced effect on logit(r3) (unitless)") # Table 3 beta31 = -1.440 (SE 1.390)
    e_trt_tozinameran_r3 <- -0.862
    label("BNT162b2 vs mRNA-1273 effect on logit(r3) (unitless)") # Table 3 beta32 = -0.862 (SE 1.110)
    e_sars_seropos_r4 <- -2.100
    label("Convalesced effect on logit(r4) (unitless)") # Table 3 beta41 = -2.100 (SE 0.431)
    e_sars_seropos_nab0 <- 16.500
    label("Convalesced effect on logit(A0) (unitless)") # Table 3 beta01 = 16.500 (SE 72.500)

    # IIV on the logit scale. Table 3 reports SDs omega1-omega4 and the six
    # pairwise correlations; covariance = rho * omega_i * omega_j.
    etalogitr1 + etalogitr2 + etalogitr3 + etalogitr4 ~ c(
      0.345^2,
      -0.260 * 0.345 * 0.611,
      0.611^2,
      0.662 * 0.345 * 0.323,
      0.151 * 0.611 * 0.323,
      0.323^2,
      -0.145 * 0.345 * 0.675,
      0.841 * 0.611 * 0.675,
      0.346 * 0.323 * 0.675,
      0.675^2
    ) # Table 3 omega1-omega4 = 0.345, 0.611, 0.323, 0.675; rho(r2,r1) = -0.260, rho(r3,r1) = 0.662, rho(r4,r1) = -0.145, rho(r3,r2) = 0.151, rho(r4,r2) = 0.841, rho(r4,r3) = 0.346
    etalogitnab0 ~ 3.630^2 # Table 3 omega0 = 3.630 (SE 1.320); uncorrelated with r1-r4 per Results

    # Residual error on the logit scale (Monolix combined1, Eq. 7)
    addSd <- 0.348
    label("Additive residual error on the logit-transformed NAb scale, Monolix a (unitless)") # Table 3 a = 0.348 (SE 0.0266)
    propSd <- 0.206
    label("Proportional residual error on the logit-transformed NAb scale, Monolix b (fraction)") # Table 3 b = 0.206 (SE 0.0263)
  })

  model({
    vmax <- exp(lvmax)
    km <- exp(lkm)

    # Individual parameters: logit-normal (Eq. 8)
    r1 <- expit(logitr1 + e_sars_seropos_r1 * SARS_SEROPOS + e_trt_tozinameran_r1 * TRT_TOZINAMERAN + etalogitr1)
    r2 <- expit(logitr2 + e_sars_seropos_r2 * SARS_SEROPOS + e_trt_tozinameran_r2 * TRT_TOZINAMERAN + etalogitr2)
    r3 <- expit(logitr3 + e_sars_seropos_r3 * SARS_SEROPOS + e_trt_tozinameran_r3 * TRT_TOZINAMERAN + etalogitr3)
    r4 <- expit(logitr4 + e_sars_seropos_r4 * SARS_SEROPOS + etalogitr4)
    nab0 <- expit(logitnab0 + e_sars_seropos_nab0 * SARS_SEROPOS + etalogitnab0)

    # Vaccine kinetics (Eq. 6): u(t) enters as a dose of amt = 1 infused over
    # 1 day into depot_antigen; V(0) = 0.
    d/dt(depot_antigen) <- -vmax * depot_antigen / (km + depot_antigen)

    # NAb dynamics (Eq. 5)
    d/dt(antibody) <- r1 * depot_antigen + r2 * antibody * depot_antigen + antibody * (r3 - r4 * antibody)
    antibody(0) <- nab0

    # Observation: NAb activity as a 0-1 fraction. Eq. 7 applies the logit
    # with Monolix's default bounds m = -0.00954434 and M = 1.01970 to both
    # the observation and the prediction, then a combined1 error a + b * f on
    # that transformed scale.
    nab_activity <- antibody
    nab_activity ~ add(addSd) + propT(propSd) + combined1() + logitNorm(NA, -0.00954434, 1.01970)
  })
}
