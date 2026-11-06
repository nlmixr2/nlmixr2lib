Lindemalm_2005_cladribine <- function() {
  description <- "Three-compartment population PK model for cladribine (2-chlorodeoxyadenosine) in adults with indolent lymphoid and myeloid malignancies (Lindemalm 2005), pooled across intravenous infusion, oral solution and subcutaneous injection: first-order absorption from separate oral (F = 0.353) and subcutaneous (F = 1 fixed) depots, exponential IIV with one variance shared by Q2, Q3 and V3, and route-specific proportional residual error with between-subject variability on the residual magnitude. No covariates retained."
  reference <- "Lindemalm S, Savic RM, Karlsson MO, Juliusson G, Liliemark J, Albertioni F. Application of population pharmacokinetics to cladribine. BMC Pharmacol. 2005;5:4. doi:10.1186/1471-2210-5-4. PMID: 15757511; PMCID: PMC1079880."
  vignette <- "Lindemalm_2005_cladribine"
  units <- list(time = "h", dosing = "mg", concentration = "nmol/L")

  compartmentData <- list(
    depot_oral = list(analyte = "cladribine", units = "mg", specimen = "administration site", verified = TRUE),
    depot_sc = list(analyte = "cladribine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "cladribine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cladribine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "cladribine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    ROUTE_ORAL = list(
      description = "Observation follows an oral dose (1) or a parenteral dose (0)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (intravenous infusion or subcutaneous injection; the two are separated by ROUTE_IV)",
      notes = "Per-record indicator read by the residual-error model only: Lindemalm 2005 Table 3 estimates a separate proportional residual error for infusion (0.191), oral (0.232) and subcutaneous (0.162) data. The structural route branch (oral ka and F vs subcutaneous ka and F) is carried by the dose record's target compartment (depot_oral vs depot_sc), so ROUTE_ORAL does not enter any structural equation. Set it on observation rows according to the route of the dose the observation follows.",
      source_name = "route of administration (Table 3 column 'Oral')"
    ),
    ROUTE_IV = list(
      description = "Observation follows an intravenous dose (1) or a subcutaneous dose (0); ignored when ROUTE_ORAL = 1",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (subcutaneous injection)",
      notes = "Per-record indicator read by the residual-error model only, selecting the infusion proportional residual error (0.191) against the subcutaneous one (0.162) on non-oral records (Lindemalm 2005 Table 3). Intravenous doses (2-h or continuous infusion) go straight into 'central' via the dose record's cmt, so the indicator does not enter any structural equation.",
      source_name = "route of administration (Table 3 columns 'Infusion' / 'Subcutaneous')"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Screened (gender) by the generalised additive model and stepwise NONMEM search; not retained (Lindemalm 2005 Results: 'None of the covariate relations ... were found to be significant')."
    ),
    AGE = list(
      description = "Age at treatment",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened; not retained (Lindemalm 2005 Results)."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened; not retained. The Discussion reports that body weight and BSA did not reduce the interindividual variability."
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened; not retained (Lindemalm 2005 Results)."
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened; not retained (Lindemalm 2005 Discussion)."
    ),
    CRCL = list(
      description = "Cockcroft-Gault creatinine clearance",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Table 2: mean 81, SD 27, range 5-162 mL/min); not retained. No patient with severe renal impairment was included."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 161L,
    n_studies = 5L,
    n_observations = 1102L,
    age_range = "22-89 years (mean 60, SD 13)",
    weight_range = "48-118 kg (mean 76, SD 14)",
    height_range = "152-198 cm (mean 174, SD 9)",
    bsa_range = "1.5-2.4 m^2 (mean 1.9, SD 0.2)",
    sex_female_pct = 20.9,
    disease_state = "Indolent B- and T-cell lymphoid and myeloid malignancies: hairy cell leukaemia (84), chronic lymphocytic leukaemia (63), non-Hodgkin lymphoma (7), chronic myeloid leukaemia (4), acute myeloid leukaemia (3), Langerhans cell histiocytosis (2) among the 163 patients initially included",
    dose_range = "2-h IV infusion 5 mg/m^2 or 0.12 mg/kg (93 doses / 63 courses); oral solution 10 mg/m^2 or 0.24 mg/kg after overnight fast (24 doses / 22 courses); SC injection 5 mg/m^2 into abdominal-wall adipose tissue (83 doses / 83 courses); six patients received a 4-7 day continuous IV infusion of 4.0-5.6 mg/m^2/24 h",
    regions = "Sweden (Karolinska, Huddinge and Linkoping hospitals) and the United Kingdom (Taunton)",
    renal_function = "Creatinine clearance (Cockcroft-Gault) mean 81 mL/min (SD 27, range 5-162); no patient with severe renal or hepatic impairment",
    notes = "Retrospective pooled analysis of five previously published PK studies (1990-1996): 168 courses in 161 patients analysed. Demographics (Lindemalm 2005 Table 1) describe the 163 patients initially included (129 male / 34 female). Seven patients contributed both intravenous and oral data."
  )

  ini({
    # Structural parameters: Lindemalm 2005 Table 3 (final model, typical individual)
    lcl <- log(39.3); label("Clearance (L/h)") # Table 3 Clearance 39.3 L/h (RSE 4.9%)
    lvc <- log(71.7); label("Central volume of distribution V1 (L)") # Table 3 V1 71.7 L (RSE 13%)
    lq <- log(51.1); label("Intercompartmental clearance Q2 to the large peripheral compartment (L/h)") # Table 3 Q2 51.1 L/h (RSE 6.8%)
    lvp <- log(475); label("Peripheral volume of distribution V2 (L)") # Table 3 V2 475 L (RSE 1.8%)
    lq2 <- log(105); label("Intercompartmental clearance Q3 to the small peripheral compartment (L/h)") # Table 3 Q3 105 L/h (RSE 21%)
    lvp2 <- log(73.6); label("Peripheral volume of distribution V3 (L)") # Table 3 V3 73.6 L (RSE 13%)
    lka_oral <- log(1.31); label("Oral first-order absorption rate constant (1/h)") # Table 3 Oral Ka 1.31 1/h (RSE 14%)
    lfdepot_oral <- log(0.353); label("Oral bioavailability (fraction)") # Table 3 Oral F 0.353 (RSE 7.9%)
    lka_sc <- log(2.48); label("Subcutaneous first-order absorption rate constant (1/h)") # Table 3 Subcutaneous Ka 2.48 1/h (RSE 9.6%)
    lfdepot_sc <- fixed(log(1)); label("Subcutaneous bioavailability (fraction)") # Results: subcutaneous F 'was fixed to this value' (1)

    # IIV: exponential models (Methods 'Statistical model'); Table 3 reports CV%,
    # converted as omega^2 = log(1 + CV^2). Q2, Q3 and V3 carry one identical
    # estimate and RSE (61%, RSE 17%), encoded as a single eta shared by the three.
    etalcl ~ log(1 + 0.54^2) # Table 3 IIV clearance 54% (RSE 17%)
    etalvc ~ log(1 + 0.34^2) # Table 3 IIV V1 34% (RSE 62%)
    etalq ~ log(1 + 0.61^2) # Table 3 IIV Q2 = Q3 = V3 61% (RSE 17%); shared eta on lq, lq2 and lvp2
    etalvp ~ log(1 + 0.70^2) # Table 3 IIV V2 70% (RSE 31%)
    etalka_oral ~ log(1 + 0.75^2) # Table 3 IIV oral Ka 75% (RSE 50%)
    etalfdepot_oral ~ log(1 + 0.041^2) # Table 3 IIV oral F 4.1% (RSE 63%)

    # Residual error: proportional, separate per route (Table 3, footnote 2),
    # with IIV on the residual magnitude (Table 3 footnote 3; Methods ref 19)
    propSdIv <- 0.191; label("Proportional residual error after IV infusion (fraction)") # Table 3 Residual error, Infusion 0.191 (RSE 8.2%)
    propSdOral <- 0.232; label("Proportional residual error after oral dosing (fraction)") # Table 3 Residual error, Oral 0.232 (RSE 9.8%)
    propSdSc <- 0.162; label("Proportional residual error after SC injection (fraction)") # Table 3 Residual error, Subcutaneous 0.162 (RSE 8.8%)
    etapropSd ~ log(1 + 0.24^2) # Table 3 IIV in residual variability 24% (RSE 28%)
  })

  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)
    q2 <- exp(lq2 + etalq)
    vp2 <- exp(lvp2 + etalq)
    ka_oral <- exp(lka_oral + etalka_oral)
    ka_sc <- exp(lka_sc)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # Route is selected by the dose record's cmt: 'central' for IV infusion,
    # 'depot_oral' for the oral solution, 'depot_sc' for subcutaneous injection.
    d/dt(depot_oral) <- -ka_oral * depot_oral
    d/dt(depot_sc) <- -ka_sc * depot_sc
    d/dt(central) <- ka_oral * depot_oral + ka_sc * depot_sc - kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    f(depot_oral) <- exp(lfdepot_oral + etalfdepot_oral)
    f(depot_sc) <- exp(lfdepot_sc)

    # Dose in mg, volume in L -> mg/L; x 1e6 / 285.69 (cladribine molar mass,
    # g/mol) -> nmol/L, the unit of the paper's assay and figures.
    Cc <- central / vc * 1e6 / 285.69

    propSdRoute <- ROUTE_ORAL * propSdOral + (1 - ROUTE_ORAL) * (ROUTE_IV * propSdIv + (1 - ROUTE_IV) * propSdSc)
    propSdi <- propSdRoute * exp(etapropSd)
    Cc ~ prop(propSdi)
  })
}
