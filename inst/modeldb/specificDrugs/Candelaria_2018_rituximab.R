Candelaria_2018_rituximab <- function() {
  description <- "Two-compartment population PK model of rituximab (and its biosimilar RTXM83) with linear distribution and linear elimination from the central compartment in patients with diffuse large B-cell lymphoma (DLBCL) treated with rituximab-CHOP or RTXM83-CHOP (Candelaria 2018; pooled-arm fit, all 5341 concentrations from both treatment arms)"
  reference <- "Candelaria M, Gonzalez D, Fernandez Gomez FJ, Paravisini A, Del Campo Garcia A, Perez L, Miguel-Lillo B, Millan S. Comparative assessment of pharmacokinetics, and pharmacodynamics between RTXM83, a rituximab biosimilar, and rituximab in diffuse large B-cell lymphoma patients: a population PK model approach. Cancer Chemother Pharmacol. 2018 Mar;81(3):515-527. doi:10.1007/s00280-018-3524-9"
  vignette <- "Candelaria_2018_rituximab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    BSA = list(
      description        = "Body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. Power effect on V1 (central volume) with exponent 1.11; reference 1.72 m^2 (Candelaria 2018 Table 1 median). The paper uses median-centered power equations for continuous covariates (Methods 'Covariate analysis'). BSA computation formula is not specified in the paper (assume unspecified).",
      source_name        = "BSA"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 251L,
    n_studies      = 1L,
    age_range      = "18-66 years",
    age_median     = "51 years",
    weight_range   = "33.0-137.3 kg",
    weight_median  = "66.0 kg",
    sex_female_pct = 45.4,
    bsa_range      = "1.14-2.54 m^2",
    bsa_median     = "1.72 m^2",
    disease_state  = "Diffuse large B-cell lymphoma (DLBCL), Ann Arbor stage II-IV with bulky disease, ECOG performance status 0-2, IPI score 0-1; first-line treatment.",
    dose_range     = "375 mg/m^2 IV every 3 weeks for 1-6 cycles, co-administered with CHOP chemotherapy (cyclophosphamide, doxorubicin, vincristine, prednisone). Loading and steady-state PK sampled at cycle 1 and cycle 6.",
    regions        = "58 research sites in 12 countries (Argentina, Brazil, Colombia, India, Indonesia, Iran, Malaysia, Philippines, Russian Federation, South Africa).",
    n_observations = 5341L,
    treatment_arms = "RTXM83-CHOP n=127 (50.6%); rituximab reference-CHOP n=124 (49.4%). PK pooled across both arms (similarity demonstrated).",
    notes          = "Final pooled-arm population PK model in Candelaria 2018 Table 2. Both RTXM83 (biosimilar) and rituximab reference product were fitted with a single structural model demonstrating bioequivalence (AUC and Cmax 90% CI within 0.80-1.25 at both cycle 1 and cycle 6). Covariates evaluated and not retained: age, weight, sex, CD19+ count, bone marrow involvement, ECOG PS, bulky disease, extranodal lesions, IPI score, ADA status. Reference (median) covariate values: BSA 1.72 m^2 (Table 1). NCT02268045 (RTXM83-AC-01-11). Source NONMEM 7.3.0 with PsN 4.4.8; FOCE with INTERACTION."
  )

  ini({
    # Structural parameters (typical values for the reference subject;
    # Candelaria 2018 Table 2 "Fixed-effect parameters"). The paper reports
    # CL and Q in mL/h, V1 and V2 in mL. Internal time unit here is days, so
    # clearances are converted: L/h * 24 = L/day. Volumes converted from mL
    # to L by /1000.
    lcl <- log(12.5 * 24 / 1000); label("Typical clearance (CL, L/day) at reference BSA 1.72 m^2")               # Candelaria 2018 Table 2 CL = 12.5 mL/h (RSE 2.35%)
    lvc <- log(3191 / 1000);      label("Typical central volume of distribution (V1, L) at reference BSA 1.72 m^2") # Candelaria 2018 Table 2 V1 = 3191 mL (RSE 1.36%)
    lq  <- log(18.6 * 24 / 1000); label("Typical inter-compartmental clearance (Q, L/day)")                       # Candelaria 2018 Table 2 Q = 18.6 mL/h (RSE 6.40%)
    lvp <- log(4154 / 1000);      label("Typical peripheral volume of distribution (V2, L)")                      # Candelaria 2018 Table 2 V2 = 4154 mL (RSE 2.77%)

    # Covariate effect: BSA on V1 only (median-centered power form per Methods
    # 'Covariate analysis'). The paper also evaluated BSA on CL but discarded
    # it (Results: "the contribution of the effect of BSA on CL was not
    # finally included as the contribution to the effect was limited with an
    # exponent estimated to be 1"). BSA inclusion on V1 reduced IIV from 39%
    # to 14% (Results paragraph following Table 1).
    e_bsa_vc <- 1.11; label("Power exponent of BSA on V1 (unitless; V1 = V1_typ * (BSA/1.72)^e_bsa_vc)")           # Candelaria 2018 Table 2 V1-BSA = 1.11 (RSE 9.91%)

    # Inter-individual variability (Candelaria 2018 Table 2 "Random-effect
    # parameters (CV, %)"). The paper's footnote (b) explicitly gives the
    # log-normal mapping: CV(%) = 100 * sqrt(exp(omega^2) - 1), so the
    # stored variance is omega^2 = log(1 + (CV/100)^2).
    etalcl ~ 0.05923   # 24.7% CV -> log(1 + 0.247^2) = 0.05923; Candelaria 2018 Table 2 IIV CL
    etalvc ~ 0.01996   # 14.2% CV -> log(1 + 0.142^2) = 0.01996; Candelaria 2018 Table 2 IIV V1
    etalq  ~ 0.07549   # 28.0% CV -> log(1 + 0.280^2) = 0.07549; Candelaria 2018 Table 2 IIV Q
    etalvp ~ 0.07037   # 27.0% CV -> log(1 + 0.270^2) = 0.07037; Candelaria 2018 Table 2 IIV V2

    # Inter-occasion variability (IOV) on CL (35.9% CV) and V1 (16.8% CV) was
    # reported in Candelaria 2018 Table 2 "Random-effect parameters", with
    # occasions defined as cycle 1 vs cycles 2-6 (Methods 'Statistical
    # model'). IOV is not encoded structurally here because nlmixr2lib model
    # files target a single subject-level eta per parameter and there is no
    # standardized OCC indicator in the package's event-data schema.
    # Downstream users who want to simulate IOV can add an OCC indicator
    # column to the event dataset and a per-occasion eta in rxode2. See the
    # vignette 'Assumptions and deviations' section for this point.

    # Residual error (Candelaria 2018 Table 2 "Residual variability").
    # Combined proportional + additive on the linear scale: the proportional
    # value 27% is interpreted as the SD as a fraction of the prediction
    # (CV) and the additive 278 ng/mL is interpreted as an SD in ng/mL
    # which converts to 0.278 ug/mL to match this model's concentration
    # unit declaration.
    propSd <- 0.27;  label("Proportional residual error (SD, fraction)")  # Candelaria 2018 Table 2 proportional = 27% (RSE 3.29%)
    addSd  <- 0.278; label("Additive residual error (ug/mL)")             # Candelaria 2018 Table 2 additive = 278 ng/mL = 0.278 ug/mL (RSE 26.5%)
  })

  model({
    # Individual PK parameters. Reference subject: BSA = 1.72 m^2 (median
    # of the pooled DLBCL population, Candelaria 2018 Table 1). The BSA
    # power effect is applied only to V1 per the final covariate model.
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc) * (BSA / 1.72)^e_bsa_vc
    q  <- exp(lq  + etalq)
    vp <- exp(lvp + etalvp)

    # Two-compartment model with IV (central) dosing and linear elimination
    # from the central compartment (Candelaria 2018 Methods "Structural PK
    # model").
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Concentration: dose in mg, volume in L -> mg/L = ug/mL.
    Cc <- central / vc

    Cc ~ add(addSd) + prop(propSd)
  })
}
