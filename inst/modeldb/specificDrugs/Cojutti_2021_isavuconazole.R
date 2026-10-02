Cojutti_2021_isavuconazole <- function() {
  description <- paste(
    "Two-compartment population PK model for isavuconazole (dosed as the prodrug",
    "isavuconazonium sulfate, doses expressed as isavuconazole) in hospitalized adults",
    "treated for invasive fungal disease, mostly invasive pulmonary aspergillosis, with",
    "first-order oral absorption, oral bioavailability, 1-h intravenous infusion into",
    "the central compartment, and linear elimination (Cojutti 2021). No covariates were",
    "retained. Fitted non-parametrically with the NPAG algorithm in Pmetrics; the",
    "Table 3 medians are encoded as lognormal medians and the tabulated CV percentages",
    "as independent lognormal marginal variances.",
    sep = " "
  )
  reference <- paste(
    "Cojutti PG, Carnelutti A, Lazzarotto D, Sozio E, Candoni A, Fanin R, Tascini C, Pea F.",
    "Population Pharmacokinetics and Pharmacodynamic Target Attainment of Isavuconazole",
    "against Aspergillus fumigatus and Aspergillus flavus in Adult Patients with Invasive",
    "Fungal Diseases: Should Therapeutic Drug Monitoring for Isavuconazole Be Considered as",
    "Mandatory as for the Other Mold-Active Azoles? Pharmaceutics. 2021;13(12):2099.",
    "doi:10.3390/pharmaceutics13122099. PMCID PMC8708495.",
    sep = " "
  )
  vignette <- "Cojutti_2021_isavuconazole"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Checked against Cojutti 2021 Methods 2.1-2.2 and
  # Results 3.2.
  compartmentData <- list(
    depot = list(
      analyte = "isavuconazole",
      units = "mg",
      specimen = "administration site",
      verified = TRUE,
      notes = "Oral doses only (Results 3.2: 'first-order input (for orally administered doses)'). Intravenous doses go directly into central as a 1-h infusion (Methods 2.1)."
    ),
    central = list(
      analyte = "isavuconazole",
      units = "mg",
      specimen = "serum",
      verified = TRUE,
      notes = "Methods 2.1: blood samples 'were centrifuged to obtain serum. Isavuconazole serum concentrations were estimated by means of a validated liquid chromatography-tandem mass spectrometry analytic method'. The paper also calls these 'plasma' concentrations."
    ),
    peripheral1 = list(
      analyte = "isavuconazole",
      units = "mg",
      specimen = "serum",
      verified = TRUE,
      notes = "Kinetic peripheral compartment of the two-compartment model; amounts expressed in the same serum-concentration reference as central."
    )
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Tested on isavuconazole CL in the population model but 'did not improve model fit' (Results 3.2). Associated with a 3.7 percent lower Ctrough per year in the separate mixed-effect regression of Table 2. Cohort median 61.5 years (IQR 51.3-72.0), Table 1."
    ),
    CONMED_CYP3A4_INH = list(
      description = "Concomitant mild or moderate CYP3A4 inhibitor (1 = yes, 0 = no)",
      units = "(binary)",
      type = "binary",
      notes = "Tested on isavuconazole CL in the population model but 'did not improve model fit' (Results 3.2), although it was associated with a 2.154 mg/L higher Ctrough in the Table 2 multivariate regression. 10 of 50 patients were cotreated with mild or moderate CYP3A4 inhibitors (loperamide, haloperidol, venetoclax, cyclosporine, letermovir, sorafenib); none received strong CYP3A4 inhibitors or inducers (Results)."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened for Ctrough in the Table 2 regression (not retained); the population model retained no covariate. Cohort median 65.0 kg (IQR 55.5-71.5), Table 1."
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "Screened for Ctrough in the Table 2 regression (not retained). 19 of 50 patients female, Table 1."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened for Ctrough in the Table 2 regression (not retained). Cohort median 35.0 g/L (IQR 28.4-40.0), Table 1."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Screened for Ctrough in the Table 2 regression (not retained). Cohort median 0.28 mg/dL (IQR 0.2-0.4) in Table 1, about 4.8 umol/L."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened for Ctrough in the Table 2 regression (not retained). Cohort median 21.0 IU/L (IQR 15.0-38.0), Table 1."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened for Ctrough in the Table 2 regression (not retained). Cohort median 20.0 IU/L (IQR 15.0-31.0), Table 1."
    ),
    GGT = list(
      description = "Gamma-glutamyltransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened for Ctrough in the Table 2 regression (retained in the multivariate regression but not significant, p = 0.751); not in the population model. Cohort median 70.0 IU/L (IQR 42.0-173.0), Table 1."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 50L,
    n_studies = 1L,
    n_observations = 199L,
    age_median = "61.5 years (IQR 51.3-72.0)",
    weight_median = "65.0 kg (IQR 55.5-71.5)",
    sex_female_pct = 38,
    race_ethnicity = "not reported in the source paper",
    disease_state = paste(
      "Hospitalized adults treated with isavuconazole for invasive fungal disease:",
      "invasive pulmonary aspergillosis 80 percent, invasive fusariosis 4 percent,",
      "cerebral mucormycosis, Scedosporium osteomyelitis and Aspergillus brain abscess",
      "2 percent each, and unspecified invasive fungal disease 10 percent. Underlying",
      "disease: oncohaematological malignancy 50 percent, nosocomial pneumonia 22 percent,",
      "immunosuppression (solid organ transplant, solid malignancy, rheumatological",
      "disease) 18 percent, other 10 percent. 45 of 50 received isavuconazole first line.",
      sep = " "
    ),
    hepatic_function = "Albumin median 35.0 g/L; total bilirubin median 0.28 mg/dL; ALT median 21.0 IU/L; AST median 20.0 IU/L; gamma-GT median 70.0 IU/L (Table 1)",
    dose_range = "Loading 200 mg every 8 h for 2 days, then 200 mg once daily orally (76 percent of patients) or as a 1-h intravenous infusion; median treatment 48 days (IQR 19-91)",
    regions = "Italy (Azienda Sanitaria Universitaria Friuli Centrale, Udine)",
    notes = paste(
      "Monocentric retrospective observational TDM study, September 2017 to November 2020",
      "(Methods 2.1). 199 serum concentrations: 175 troughs drawn about 5 min before the",
      "daily dose and 24 peaks drawn 2 h after an oral dose or 0.5 h after the end of a",
      "1-h infusion (Table 1, Results 3.1). LC-MS/MS assay LOQ 0.11 mg/L, linear",
      "0.1-10 mg/L. Median 2 TDM instances per patient. Baseline characteristics from",
      "Table 1.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------------
    # Structural parameters: the MEDIAN row of Cojutti 2021 Table 3
    # ('Parameter estimates for the final population pharmacokinetic model of
    # isavuconazole'), which gives mean, SD, CV (%) and median of the NPAG
    # non-parametric distribution for each parameter. The median column is
    # encoded because it reproduces the paper's own Monte Carlo Ctrough
    # probabilities (Table 4) and the mean column does not (the mean column
    # puts 18.5 percent of troughs below 1 mg/L at the end of loading against
    # the published 1.7 percent); see the vignette 'Choice of typical values'.
    # ------------------------------------------------------------------------
    lka <- log(22.64); label("First-order absorption rate constant (1/h)")
    # Table 3, 'Ka (h-1)' median = 22.64 (mean 22.64, SD 3.54, CV 15.66 percent)
    lcl <- log(1.33); label("Clearance (L/h)")
    # Table 3, 'CL (L/h)' median = 1.33 (mean 1.52, SD 0.97, CV 64.03 percent)
    lvc <- log(102.58); label("Central volume of distribution (L)")
    # Table 3, 'V (L)' median = 102.58 (mean 89.50, SD 42.38, CV 47.35 percent)
    lq <- log(5.08); label("Intercompartmental clearance (L/h)")
    # Table 3, 'Q (L/h)' median = 5.08 (mean 16.78, SD 18.35, CV 109.37 percent)
    lvp <- log(385.93); label("Peripheral volume of distribution (L)")
    # Table 3, 'Vp (L)' median = 385.93 (mean 735.24, SD 633.89, CV 86.22 percent)
    lfdepot <- log(1.00); label("Oral bioavailability (fraction)")
    # Table 3, 'Fos (%)' median = 1.00 (mean 0.95, SD 0.07, CV 7.42 percent); the
    # column header says percent but the values are fractions

    # ------------------------------------------------------------------------
    # Inter-individual variability. NPAG estimates a discrete joint density,
    # not a parametric omega. Table 3's CV column is SD / mean (0.97 / 1.52 =
    # 64 percent, and likewise for every row), carried here as a lognormal
    # marginal with omega^2 = log(CV^2 + 1). No covariances are reported, so
    # the etas are independent. The 7.42 percent CV on Fos is NOT encoded: the
    # median sits at the upper bound of 1, so any lognormal or logit-normal
    # marginal around it is either physically impossible (F > 1) or undefined.
    # ------------------------------------------------------------------------
    etalka ~ 0.0242277 # Table 3 Ka CV 15.66 percent -> log(0.1566^2 + 1)
    etalcl ~ 0.343578 # Table 3 CL CV 64.03 percent -> log(0.6403^2 + 1)
    etalvc ~ 0.202289 # Table 3 V CV 47.35 percent -> log(0.4735^2 + 1)
    etalq ~ 0.786719 # Table 3 Q CV 109.37 percent -> log(1.0937^2 + 1)
    etalvp ~ 0.555831 # Table 3 Vp CV 86.22 percent -> log(0.8622^2 + 1)

    # ------------------------------------------------------------------------
    # Residual error. Methods 2.2: assay error SD = C0 + C1 * C with
    # C0 = 0.006 and C1 = 0.189, and 'extra process noise was captured with a
    # gamma (G) model (G = 2)'. Pmetrics multiplies the assay SD by gamma, so
    # the residual SD is 2 * (0.006 + 0.189 * C) = 0.012 + 0.378 * C. Pmetrics
    # sums the two parts linearly, which is nlmixr2's combined1(). Not fixed:
    # gamma is the estimated part of the Pmetrics error model.
    # ------------------------------------------------------------------------
    addSd <- 0.012; label("Additive residual SD (mg/L); gamma 2 x C0 0.006")
    # Methods 2.2: C0 = 0.006, G = 2
    propSd <- 0.378; label("Proportional residual SD (fraction); gamma 2 x C1 0.189")
    # Methods 2.2: C1 = 0.189, G = 2
  })
  model({
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- exp(lfdepot)

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
