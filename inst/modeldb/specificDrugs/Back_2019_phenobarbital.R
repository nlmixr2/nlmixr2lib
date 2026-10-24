Back_2019_phenobarbital <- function() {
  description <- "One-compartment population PK model for oral and intravenous phenobarbital in Korean preterm and term neonates and infants (Back 2019), with first-order absorption, allometric size scaling on CL (exponent 0.75) and Vd (exponent 1) referenced to 70 kg, and an estimated sigmoidal postmenstrual-age maturation function on CL. The printed Vd (5.51 L at 70 kg, about 0.08 L/kg) gives an 11-h half-life in the typical neonate, which does not reproduce the flat trough profile of the paper's own VPC (Figure 3b); the Table 3 values are encoded verbatim (see vignette Errata)."
  reference <- "Back HM, Lee JB, Han N, Goo S, Jung E, Kim J, Song B, An SH, Kim JT, Rhie SJ, Ree YS, Chae JW, Kim JW, Yun HY. Application of Size and Maturation Functions to Population Pharmacokinetic Modeling of Pediatric Patients. Pharmaceutics. 2019;11(6):259. doi:10.3390/pharmaceutics11060259"
  vignette <- "Back_2019_pediatric_size_maturation"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric size function referenced to WTnormal = 70 kg (Back 2019 Methods 2.4, Eq 4); exponent 0.75 on CL (Eq 8) and 1 on Vd (Eq 7).",
      source_name = "WT"
    ),
    PAGE = list(
      description = "Postmenstrual (post-conceptional) age",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = "Back 2019 calls this post-conceptional age (PCA); Table 1 values equal gestational age plus postnatal age (36.7 weeks + 32.4 days = 41.3 weeks), i.e. the standard postmenstrual age. Drives the sigmoidal maturation function on CL (Eq 5 / Eq 8), with TM50 in weeks.",
      source_name = "PCA"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "phenobarbital", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "phenobarbital", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 28L,
    n_studies = 1L,
    age_range = "Postnatal age 3-150 days (mean 32.4, SD 30.7); postmenstrual age 31-51.1 weeks (mean 41.3, SD 3.9)",
    weight_range = "1-6.9 kg (mean 3.3, SD 1)",
    sex_female_pct = 60.7,
    race_ethnicity = "Korean (single-centre NICU, Kyunghee University Hospital at Gang-dong, Seoul)",
    disease_state = "Preterm and term neonates and infants receiving phenobarbital for seizure control in the neonatal intensive care unit",
    dose_range = "Loading 15-20 mg/kg, then maintenance 3-5 mg/kg/day starting 12-24 h after the loading dose; oral or intravenous",
    regions = "South Korea",
    ga_range = "23.6-41.7 weeks (mean 36.7, SD 4.4)",
    birth_weight_range = "0.4-3.81 kg (mean 2.64, SD 0.87)",
    notes = "Retrospective therapeutic-drug-monitoring data; samples taken between 5 min and 3 h before the next dose and assayed by turbidimetric immunoassay (Cobas 6000). Demographics from Back 2019 Table 1."
  )

  ini({
    # Back 2019 Table 3 'Population Mean' for PHB: values at the 70-kg,
    # fully mature reference subject of Eq 7 and Eq 8.
    lcl <- log(0.569); label("Clearance at 70 kg, fully mature (L/h)") # Back 2019 Table 3 PHB: CL = 0.569 L/hr (RSE 5.0%)
    lvc <- log(5.51); label("Volume of distribution at 70 kg (L)") # Back 2019 Table 3 PHB: Vd = 5.51 L (RSE 2.1%); see vignette Errata
    lka <- fixed(log(50)); label("First-order absorption rate constant (1/h)") # Back 2019 Table 3 PHB: ka = 50 hr-1, footnote '*: Fixed parameter'
    lfdepot <- log(0.724); label("Oral bioavailability (fraction)") # Back 2019 Table 3 PHB: Bioavailability = 0.724 (RSE 7.2%)

    ltmat50 <- log(48.2); label("Postmenstrual age at half-maximal CL maturation (weeks)") # Back 2019 Table 3 PHB: TM50 = 48.2 week (RSE 2.1%)
    lhill_mat <- log(5.99); label("Hill coefficient of the CL maturation function (unitless)") # Back 2019 Table 3 PHB: Hill coefficient = 5.99 (RSE 1.2%)

    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL (unitless)") # Back 2019 Eq 8 (PHB, VAN branch): exponent printed as 0.75
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on Vd (unitless)") # Back 2019 Eq 7: exponent printed as 1

    # IIV reported as CV%; omega^2 = log(CV^2 + 1)
    etalcl ~ 0.154 # Back 2019 Table 3 PHB: IIV CL 40.8% CV
    etalvc ~ 0.482 # Back 2019 Table 3 PHB: IIV Vd 78.7% CV

    propSd <- 0.356; label("Proportional residual error (fraction)") # Back 2019 Table 3 PHB: proportional error 35.6% (Eq 6, CCV model)
  })

  model({
    tmat50 <- exp(ltmat50)
    hill_mat <- exp(lhill_mat)
    # Back 2019 Eq 5 / Eq 8: sigmoidal maturation on CL only
    fmat <- PAGE^hill_mat / (PAGE^hill_mat + tmat50^hill_mat)

    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * fmat
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    f(depot) <- exp(lfdepot)

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
