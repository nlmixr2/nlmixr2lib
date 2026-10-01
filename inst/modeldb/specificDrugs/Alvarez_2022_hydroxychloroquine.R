# Joint population pharmacokinetic model of oral hydroxychloroquine and its
# three metabolites desethylhydroxychloroquine, desethylchloroquine and
# bisdesethylchloroquine (didesethylchloroquine) in whole blood of adults
# hospitalised with COVID-19 (Alvarez 2022, Pharmaceuticals 15:256;
# doi:10.3390/ph15020256).

Alvarez_2022_hydroxychloroquine <- function() {
  description <- paste(
    "Joint parent + three-metabolite population PK model for oral",
    "hydroxychloroquine in whole blood of 100 adults hospitalised with",
    "COVID-19 (medicine wards and ICU) (Alvarez 2022). One-compartment",
    "hydroxychloroquine disposition with first-order absorption and a lag",
    "time (both fixed to published values), a non-metabolic apparent",
    "clearance, and three parallel first-order formation clearances into",
    "one-compartment desethylhydroxychloroquine, desethylchloroquine and",
    "bisdesethylchloroquine (didesethylchloroquine) compartments whose",
    "volumes equal the individual hydroxychloroquine volume. Formation",
    "fluxes carry a molar correction so each metabolite state holds mg of",
    "that metabolite. No covariates were retained (age, weight, height,",
    "BMI, sex, azithromycin co-treatment and ICU stay were screened).",
    "Proportional plus additive residual error on every analyte.",
    sep = " "
  )
  reference <- paste(
    "Alvarez JC, Davido B, Moine P, Etting I, Annane D, Larabi IA,",
    "Simon N (2022). Population Pharmacokinetics of Hydroxychloroquine",
    "and 3 Metabolites in COVID-19 Patients and",
    "Pharmacokinetic/Pharmacodynamic Application. Pharmaceuticals",
    "15(2):256. doi:10.3390/ph15020256.",
    sep = " "
  )
  vignette <- "Alvarez_2022_hydroxychloroquine"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. The dose is the labelled tablet strength
  # (Plaquenil 200 mg) in mg; the metabolite states hold mg of the named
  # metabolite after the molar correction in model(). Whole blood is the
  # assayed matrix (Section 4.2: blood collected on lithium heparin and
  # quantified by turbulent-flow LC-MS/MS).
  compartmentData <- list(
    depot = list(analyte = "hydroxychloroquine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "hydroxychloroquine", units = "mg", specimen = "whole blood", verified = TRUE),
    central_dhcq = list(
      analyte = "desethylhydroxychloroquine",
      units = "mg",
      specimen = "whole blood",
      verified = TRUE
    ),
    central_dcq = list(analyte = "desethylchloroquine", units = "mg", specimen = "whole blood", verified = TRUE),
    central_bdcq = list(
      analyte = "bisdesethylchloroquine (didesethylchloroquine)",
      units = "mg",
      specimen = "whole blood",
      verified = TRUE
    )
  )

  # The final model is covariate-free (Section 2.1.1: "Among all the
  # covariates tested on the PK parameters none had a significant effect";
  # the same for the metabolite parameters).
  covariateData <- list()

  # Screened with a power model (Section 4.3: "age, body weight, height,
  # body mass index (BMI), gender, AZT combination, and clinical unit")
  # and not retained. Documentation only.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Table 1: 60.7 +/- 15.9 years (median 62.5, range 20-94). Screened; not retained."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Table 1: 83.6 +/- 20.1 kg (median 82, range 37.5-190). Screened; not retained."
    ),
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      notes = "Table 1 reports metres: 1.71 +/- 0.094 m (median 1.73, range 1.52-1.93). Screened; not retained."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Table 1: 28.9 +/- 5.41 kg/m^2 (median 27.5, range 18.5-52.4). Screened; not retained."
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Table 1: 34 female, 66 male. Screened as 'gender'; not retained."
    ),
    CONMED_AZITHROMYCIN = list(
      description = "Concomitant azithromycin (250 mg bid on day 1 then 250 mg qd)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no azithromycin)",
      notes = "Table 1: 78 of 100 patients. Screened as 'AZT combination'; not retained."
    ),
    DIS_CRITILL = list(
      description = "Intensive care unit stay (1) versus medicine ward (0)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (medicine ward)",
      notes = "Table 1: 25 ICU, 75 medicine. Screened as 'clinical unit'; not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 100,
    n_studies = 1,
    n_observations = "333 whole-blood samples, each assayed for hydroxychloroquine and the three metabolites; 1 to 9 samples per patient (Section 2.1.1)",
    age_range = "20-94 years (mean 60.7, SD 15.9; median 62.5)",
    weight_range = "37.5-190 kg (mean 83.6, SD 20.1; median 82)",
    sex_female_pct = 34,
    disease_state = "Hospitalised COVID-19 confirmed by SARS-CoV-2 RT-PCR and/or compatible chest CT; 75 on medicine wards, 25 in the ICU",
    dose_range = "Plaquenil 200 mg bid or tid orally, preceded in 42 patients by a 400 mg bid loading dose on day 1",
    regions = "France (Raymond Poincare Hospital, Garches)",
    co_medication = "Azithromycin in 78 patients",
    notes = "Retrospective therapeutic drug monitoring cohort; samples roughly every two days over the first two weeks of treatment (Figure 2)."
  )

  ini({
    # Hydroxychloroquine absorption: fixed to Lim 2009 (reference 16) after
    # a sensitivity analysis (Section 2.1.1).
    ltlag <- fixed(log(0.389))
    label("Absorption lag time (h)") # Table 2 'Lag (fixed)' = 0.389 h
    lka <- fixed(log(1.15))
    label("First-order absorption rate constant (1/h)") # Table 2 'KA (fixed)' = 1.15 1/h

    # Hydroxychloroquine disposition.
    lcl <- log(5.60)
    label("Apparent non-metabolic hydroxychloroquine clearance CL/F (L/h)") # Table 2 'CL/F HCQ' = 5.60 L/h
    lvc <- log(1850)
    label("Apparent hydroxychloroquine volume VP/F, shared by the metabolites (L)") # Table 2 'VP/F HCQ' = 1850 L

    # Formation clearances from hydroxychloroquine to each metabolite.
    lcl_form_dhcq <- log(9.63)
    label("Hydroxychloroquine to desethylhydroxychloroquine formation clearance (L/h)") # Table 2 'CL HCQ_DesHCQ' = 9.63 L/h
    lcl_form_dcq <- log(4.99)
    label("Hydroxychloroquine to desethylchloroquine formation clearance (L/h)") # Table 2 'CL HCQ_DesCQ' = 4.99 L/h
    lcl_form_bdcq <- log(1.84)
    label("Hydroxychloroquine to bisdesethylchloroquine formation clearance (L/h)") # Table 2 'CL HCQ_DiDesCQ' = 1.84 L/h

    # Metabolite elimination clearances.
    lcl_dhcq <- log(8.89)
    label("Desethylhydroxychloroquine elimination clearance (L/h)") # Table 2 'CL DesHCQ' = 8.89 L/h
    lcl_dcq <- log(49.8)
    label("Desethylchloroquine elimination clearance (L/h)") # Table 2 'CL DesCQ' = 49.8 L/h
    lcl_bdcq <- log(11.6)
    label("Bisdesethylchloroquine elimination clearance (L/h)") # Table 2 'CL DiDesCQ' = 11.6 L/h

    # Between-subject variability, exponential model (Section 4.3). Table 2
    # prints these under 'Inter Individual Variability (omega)' with no
    # scale. They are entered as variances: the same group's companion
    # lopinavir paper (Alvarez 2021, same table layout and labelling) was
    # shown to report variances, and in this paper the Figure 8 trough
    # quantiles cannot separate the two readings (see the vignette).
    etalcl ~ 1.327 # Table 2 IIV 'CLHCQ' = 1.327
    etalvc ~ 0.889 # Table 2 IIV 'VPHCQ' = 0.889
    etalcl_dhcq ~ 0.860 # Table 2 IIV 'CL DesHCQ' = 0.860
    etalcl_dcq ~ 0.362 # Table 2 IIV 'CL DesCQ' = 0.362
    etalcl_bdcq ~ 0.953 # Table 2 IIV 'CL DiDesCQ' = 0.953

    # Residual error: proportional plus additive per analyte (Table 2,
    # 'Residual Unexplained Variability (sigma)'); additive terms are in
    # ug/L as printed.
    propSd <- 0.448
    label("Proportional residual SD, hydroxychloroquine (fraction)") # Table 2 'Proportional HCQ' = 0.448
    addSd <- 86.9
    label("Additive residual SD, hydroxychloroquine (ug/L)") # Table 2 'Additive HCQ' = 86.9 ug/L
    propSd_dhcq <- 0.428
    label("Proportional residual SD, desethylhydroxychloroquine (fraction)") # Table 2 'Proportional DesHCQ' = 0.428
    addSd_dhcq <- 6.69
    label("Additive residual SD, desethylhydroxychloroquine (ug/L)") # Table 2 'Additive DesHCQ' = 6.69 ug/L
    propSd_dcq <- 0.322
    label("Proportional residual SD, desethylchloroquine (fraction)") # Table 2 'Proportional DesCQ' = 0.322
    addSd_dcq <- 5.78
    label("Additive residual SD, desethylchloroquine (ug/L)") # Table 2 'Additive DesCQ' = 5.78 ug/L
    propSd_bdcq <- 0.0574
    label("Proportional residual SD, bisdesethylchloroquine (fraction)") # Table 2 'Proportional DiDescCQ' = 0.0574
    addSd_bdcq <- 2.49
    label("Additive residual SD, bisdesethylchloroquine (ug/L)") # Table 2 'Additive DiDescCQ' = 2.49 ug/L
  })

  model({
    # Molecular weights (g/mol) of the free bases. The source fitted molar
    # concentrations (Section 4.3: "All concentrations were expressed as
    # umol/L"), so each metabolite is formed 1:1 in moles; the ratio below
    # converts mg of hydroxychloroquine leaving central into mg of the
    # metabolite formed. HCQ C18H26ClN3O; DHCQ C16H22ClN3O; DCQ
    # C16H22ClN3; BDCQ C14H18ClN3.
    mw_hcq <- 335.87
    mw_dhcq <- 307.82
    mw_dcq <- 291.82
    mw_bdcq <- 263.77

    ka <- exp(lka)
    tlag <- exp(ltlag)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    cl_form_dhcq <- exp(lcl_form_dhcq)
    cl_form_dcq <- exp(lcl_form_dcq)
    cl_form_bdcq <- exp(lcl_form_bdcq)
    cl_dhcq <- exp(lcl_dhcq + etalcl_dhcq)
    cl_dcq <- exp(lcl_dcq + etalcl_dcq)
    cl_bdcq <- exp(lcl_bdcq + etalcl_bdcq)

    # Metabolite volumes are fixed to the individual parent volume
    # (Section 2.1.1: "a metabolite volume (VM/F) fixed to the volume of
    # the parent HCQ (VP/F)").
    vc_dhcq <- vc
    vc_dcq <- vc
    vc_bdcq <- vc

    # Figure 3 / NONMEM ADVAN5: hydroxychloroquine leaves central by its
    # own clearance CL/F and by the three formation clearances; the
    # fraction metabolised is sum(CL_form) / (CL/F + sum(CL_form)) = 0.75
    # (Section 2.1.1).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot -
      (cl + cl_form_dhcq + cl_form_dcq + cl_form_bdcq) / vc * central
    d/dt(central_dhcq) <- (mw_dhcq / mw_hcq) * cl_form_dhcq / vc * central -
      cl_dhcq / vc_dhcq * central_dhcq
    d/dt(central_dcq) <- (mw_dcq / mw_hcq) * cl_form_dcq / vc * central -
      cl_dcq / vc_dcq * central_dcq
    d/dt(central_bdcq) <- (mw_bdcq / mw_hcq) * cl_form_bdcq / vc * central -
      cl_bdcq / vc_bdcq * central_bdcq

    alag(depot) <- tlag

    # Whole-blood concentrations: mg / L x 1000 = ug/L.
    Cc <- 1000 * central / vc
    Cc_dhcq <- 1000 * central_dhcq / vc_dhcq
    Cc_dcq <- 1000 * central_dcq / vc_dcq
    Cc_bdcq <- 1000 * central_bdcq / vc_bdcq

    Cc ~ add(addSd) + prop(propSd)
    Cc_dhcq ~ add(addSd_dhcq) + prop(propSd_dhcq)
    Cc_dcq ~ add(addSd_dcq) + prop(propSd_dcq)
    Cc_bdcq ~ add(addSd_bdcq) + prop(propSd_bdcq)
  })
}
