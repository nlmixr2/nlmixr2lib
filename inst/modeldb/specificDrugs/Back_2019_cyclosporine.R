Back_2019_cyclosporine <- function() {
  description <- "One-compartment population PK model for intravenous high-dose cyclosporine A in Korean infants and children with retinoblastoma (Back 2019), with allometric size scaling on CL (exponent 0.75) and Vd (exponent 1) referenced to 70 kg and no maturation function (not significant for this cohort)."
  reference <- "Back HM, Lee JB, Han N, Goo S, Jung E, Kim J, Song B, An SH, Kim JT, Rhie SJ, Ree YS, Chae JW, Kim JW, Yun HY. Application of Size and Maturation Functions to Population Pharmacokinetic Modeling of Pediatric Patients. Pharmaceutics. 2019;11(6):259. doi:10.3390/pharmaceutics11060259"
  vignette <- "Back_2019_pediatric_size_maturation"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric size function referenced to WTnormal = 70 kg (Back 2019 Methods 2.4, Eq 4); exponent 0.75 on CL (Eq 8) and 1 on Vd (Eq 7).",
      source_name = "WT"
    )
  )

  compartmentData <- list(
    central = list(analyte = "cyclosporine A", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 34L,
    n_studies = 1L,
    age_range = "Postnatal age 1-79 months (mean 26.8, SD 17.8)",
    weight_range = "5-24 kg (mean 12.9, SD 3.8)",
    sex_female_pct = 41.2,
    race_ethnicity = "Korean (single-centre, Yonsei University Hospital, Seoul)",
    disease_state = "Infants and children with retinoblastoma receiving high-dose cyclosporine A as a P-glycoprotein inhibitor to enhance chemotherapy",
    dose_range = "High-dose cyclosporine A infused over 24 h, starting 3 h before chemotherapy on days 1 and 2; the dose amount is not reported",
    regions = "South Korea",
    notes = "Retrospective therapeutic-drug-monitoring data; samples drawn 20 h after the start of each daily infusion and assayed by LC/MS/MS. Chemotherapy cycles 5 +/- 2.8 (1-12). Demographics from Back 2019 Table 1."
  )

  ini({
    # Back 2019 Table 3 'Population Mean' for CsA: values at the 70-kg
    # reference subject of Eq 7 and Eq 8.
    lcl <- log(21.3); label("Clearance at 70 kg (L/h)") # Back 2019 Table 3 CsA: CL = 21.3 L/hr (RSE 4.4%)
    lvc <- log(218); label("Volume of distribution at 70 kg (L)") # Back 2019 Table 3 CsA: Vd = 218 L (RSE 25.5%)

    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL (unitless)") # Back 2019 Eq 8 (CsA branch): exponent printed as 0.75
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on Vd (unitless)") # Back 2019 Eq 7: exponent printed as 1

    # IIV reported as CV%; omega^2 = log(CV^2 + 1)
    etalcl ~ 0.02783 # Back 2019 Table 3 CsA: IIV CL 16.8% CV
    etalvc ~ 0.01502 # Back 2019 Table 3 CsA: IIV Vd 12.3% CV

    propSd <- 0.468; label("Proportional residual error (fraction)") # Back 2019 Table 3 CsA: proportional error 46.8% (Eq 6, CCV model)
  })

  model({
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    kel <- cl / vc

    d/dt(central) <- -kel * central

    # dose in mg, volume in L -> mg/L; x 1000 for ng/mL
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
