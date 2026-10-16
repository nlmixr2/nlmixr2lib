Back_2019_vancomycin <- function() {
  description <- "One-compartment population PK model for intravenous vancomycin in Korean preterm neonates to infants (Back 2019), with allometric size scaling on CL (exponent 0.75) and Vd (exponent 1) referenced to 70 kg and a sigmoidal postmenstrual-age maturation function on CL with TM50 and Hill held at reference values. The printed CL (69.4 L/h) and Vd (3.23 L) at 70 kg give a 1.3-min half-life in the typical neonate and do not reproduce the paper's own VPC (Figure 3c); the Table 3 values are encoded verbatim (see vignette Errata)."
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
      notes = "Back 2019 calls this post-conceptional age (PCA), i.e. gestational age plus postnatal age, the standard postmenstrual age. Drives the sigmoidal maturation function on CL (Eq 5 / Eq 8), with TM50 in weeks.",
      source_name = "PCA"
    )
  )

  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 93L,
    n_studies = 1L,
    age_range = "Postnatal age 0.1-80.4 weeks (mean 9.3, SD 12.4); postmenstrual age 25.6-110 weeks (mean 41.2, SD 14.2)",
    weight_range = "0.4-14.9 kg (mean 3.2, SD 2.6)",
    sex_female_pct = 38.7,
    race_ethnicity = "Korean (single-centre, Chungnam National University Hospital, Daejeon)",
    disease_state = "Preterm neonates to infants receiving intravenous vancomycin",
    dose_range = "Weight- and symptom-based intravenous doses (amounts not reported) at dosing intervals of 6 to 24 h",
    regions = "South Korea",
    ga_range = "22.9-40.3 weeks (mean 31.9, SD 4.7)",
    notes = "Retrospective therapeutic-drug-monitoring data; samples collected 0.5-23.5 h after a dose and assayed by turbidimetric immunoassay (TDx). Body surface area 0.2 +/- 0.1 m^2, serum albumin 2.7 +/- 0.6 g/dL. Demographics from Back 2019 Table 1."
  )

  ini({
    # Back 2019 Table 3 'Population Mean' for VAN: values at the 70-kg,
    # fully mature reference subject of Eq 7 and Eq 8.
    lcl <- log(69.4); label("Clearance at 70 kg, fully mature (L/h)") # Back 2019 Table 3 VAN: CL = 69.4 L/hr (RSE 13.7%); see vignette Errata
    lvc <- log(3.23); label("Volume of distribution at 70 kg (L)") # Back 2019 Table 3 VAN: Vd = 3.23 L (RSE 6.1%); see vignette Errata

    tmat50 <- fixed(33.3); label("Postmenstrual age at half-maximal CL maturation (weeks)") # Back 2019 Table 3 VAN: TM50 = 33.3 week, footnote '*: Fixed parameter'; reference values per Results 3.3 (ref 2, Holford 2013)
    hill_mat <- fixed(3.68); label("Hill coefficient of the CL maturation function (unitless)") # Back 2019 Table 3 VAN: Hill coefficient = 3.68, footnote '*: Fixed parameter'

    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL (unitless)") # Back 2019 Eq 8 (PHB, VAN branch): exponent printed as 0.75
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on Vd (unitless)") # Back 2019 Eq 7: exponent printed as 1

    # IIV reported as CV%; omega^2 = log(CV^2 + 1)
    etalcl ~ 0.01076 # Back 2019 Table 3 VAN: IIV CL 10.4% CV
    etalvc ~ 0.2459 # Back 2019 Table 3 VAN: IIV Vd 52.8% CV

    propSd <- 0.408; label("Proportional residual error (fraction)") # Back 2019 Table 3 VAN: proportional error 40.8% (Eq 6, CCV model)
  })

  model({
    # Back 2019 Eq 5 / Eq 8: sigmoidal maturation on CL only
    fmat <- PAGE^hill_mat / (PAGE^hill_mat + tmat50^hill_mat)

    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * fmat
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    kel <- cl / vc

    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
