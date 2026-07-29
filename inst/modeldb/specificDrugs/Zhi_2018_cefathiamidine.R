Zhi_2018_cefathiamidine <- function() {
  description <- "Two-compartment population PK model for intravenous cefathiamidine (a first-generation cephalosporin) in 54 children (age 2.0-11.8 years; weight 8.0-36.0 kg) with hematologic disease, developed in NONMEM v7.2 (FOCE-I) from 120 sparse plasma samples. Structural model: first-order elimination from a central compartment, with allometric body-weight scaling on CL, Q (exponent 0.75) and V1, V2 (exponent 1), reference weight 17.75 kg (the cohort median current weight). Inter-individual variability (exponential) is estimated for CL and V2 only; residual variability is exponential (lognormal on the linear scale). Bodyweight was the only retained covariate; age and creatinine clearance were not significant in the limited cohort (CrCL range 130-462 mL/min)."
  reference <- paste(
    "Zhi LJ, Wang L, Chen XK, Zhai XY, Wen L, Dong L, Jacqz-Aigrain E,",
    "Shi ZR, Zhao W. (2018).",
    "Population pharmacokinetics and dosing optimization of cefathiamidine",
    "in children with hematologic infection.",
    "Drug Des Devel Ther 12:1845-1853.",
    "doi:10.2147/DDDT.S160329"
  )
  vignette <- "Zhi_2018_cefathiamidine"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Current body weight (day of study)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling on CL, Q (exponent 0.75) and V1, V2 (exponent 1) with reference weight 17.75 kg (Zhi 2018 Table 2 formula and Note: median current weight in the cohort was 17.78 kg, rounded to 17.75 kg in the model).",
      source_name        = "CW"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 54L,
    n_studies      = 1L,
    age_range      = "2.0-11.8 years",
    age_median     = "4.9 years",
    weight_range   = "8.0-36.0 kg",
    weight_median  = "17.8 kg",
    sex_female_pct = 31.5,
    disease_state  = "Children with hematologic disease (leukemia n=23, immune thrombocytopenia n=10, anemia n=9, infectious mononucleosis syndrome n=2, neuroblastoma n=2, hemophagocytic syndrome n=2, others n=6) with confirmed or suspected bacterial infection; serum creatinine 10-47 micromol/L; Cockcroft-Gault-style creatinine clearance 130-462 mL/min.",
    dose_range     = "Cefathiamidine 100 mg/kg/day IV given as q12h infusion over 3-5 minutes (median 950 mg/dose, range 500-1800 mg/dose); 28-112 mg/kg/dose (median 50)",
    regions        = "Single-centre, Children's Hospital of Hebei Province, Shijiazhuang, China; enrolment April 2015 - March 2016",
    notes          = "Baseline demographics from Zhi 2018 Table 1 (54 patients, 120 samples). Mean age 5.2 (SD 2.4) years; mean current weight 18.4 (SD 6.5) kg. Concentrations were quantified by UPLC-MS/MS with LLOQ 30 ng/mL; observed range in the dataset was 65-245000 ng/mL (= 0.065-245 mg/L). All subjects received the same dosing regimen (100 mg/kg/day q12h)."
  )

  ini({
    # Structural parameters -- typical values for a child at the reference
    # weight of 17.75 kg. All four estimates come from Zhi 2018 Table 2
    # (full-data final estimates). The paper used NONMEM v7.2 with FOCE-I.
    lcl <- log(3.93); label("Clearance for a 17.75 kg child (CL, L/h)")                          # Zhi 2018 Table 2 (theta_1)
    lvc <- log(4.18); label("Central volume of distribution for a 17.75 kg child (V1, L)")       # Zhi 2018 Table 2 (theta_2)
    lq  <- log(0.50); label("Inter-compartmental clearance for a 17.75 kg child (Q, L/h)")       # Zhi 2018 Table 2 (theta_3)
    lvp <- log(1.30); label("Peripheral volume of distribution for a 17.75 kg child (V2, L)")    # Zhi 2018 Table 2 (theta_4)

    # Allometric exponents on body weight, fixed a priori per Zhi 2018
    # Methods / Results paragraph 2 ("allometric coefficients of 0.75 for
    # CL and Q, and 1 for V1 and V2") -- the standard Anderson-Holford
    # paediatric size-scaling choice; no RSE/CI reported, consistent with
    # values held fixed during estimation.
    e_wt_cl_q  <- fixed(0.75); label("Allometric (WT) exponent on CL and Q (unitless)")  # Zhi 2018 Results, covariate analysis paragraph
    e_wt_vc_vp <- fixed(1.00); label("Allometric (WT) exponent on V1 and V2 (unitless)")  # Zhi 2018 Results, covariate analysis paragraph

    # Inter-individual variability. Zhi 2018 reports IIV only for CL and
    # V2 ("interindividual variability ... was then estimated for CL and
    # V2"). Table 2 reports IIV as a percentage CV; consistent with the
    # individual-PK formula theta_i = theta_mean * exp(eta_i) where
    # eta ~ N(0, omega^2) and omega is reported as a percentage (i.e.,
    # CV% = omega * 100 on the log scale, NOT log(1 + CV^2)), the
    # back-transformations are:
    #   CL : 41.83% -> omega^2 = 0.4183^2 = 0.174975
    #   V2 : 74.23% -> omega^2 = 0.7423^2 = 0.551009
    etalcl ~ 0.174975   # Zhi 2018 Table 2 (IIV CL, 41.83% CV)
    etalvp ~ 0.551009   # Zhi 2018 Table 2 (IIV V2, 74.23% CV)

    # Residual error -- exponential model in the paper ("This exponential
    # model best described the residual variability"). NONMEM's exponential
    # error Y = IPRED * exp(eps), eps ~ N(0, sigma^2) maps to lnorm(expSd)
    # in nlmixr2. Reported as 35.78% -> expSd = 0.3578 (SD on the log
    # scale).
    expSd <- 0.3578; label("Lognormal residual error (SD on the log scale)")  # Zhi 2018 Table 2 (err(1), 35.78%)
  })
  model({
    # Individual PK parameters: allometric scaling on body weight with
    # reference 17.75 kg (Zhi 2018 Table 2 formulae and Note).
    cl <- exp(lcl + etalcl) * (WT / 17.75)^e_wt_cl_q
    vc <- exp(lvc)          * (WT / 17.75)^e_wt_vc_vp
    q  <- exp(lq)           * (WT / 17.75)^e_wt_cl_q
    vp <- exp(lvp + etalvp) * (WT / 17.75)^e_wt_vc_vp

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Plasma concentration; dose in mg, vc in L -> mg/L (= ug/mL = 1000 ng/mL).
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
