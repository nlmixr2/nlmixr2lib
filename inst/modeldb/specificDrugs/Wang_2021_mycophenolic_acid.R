Wang_2021_mycophenolic_acid <- function() {
  description <- paste0(
    "Two-compartment population PK model for mycophenolic acid (MPA) after ",
    "oral mycophenolate mofetil (MMF) dispersible tablets in 91 adult Chinese ",
    "heart transplant recipients on tacrolimus and corticosteroids (Wang ",
    "2021). First-order absorption after a lag time, first-order elimination ",
    "from the central compartment, and bioavailability fixed at 0.95. ",
    "Concomitant proton-pump inhibitor use lowers F by a factor of 0.724; ",
    "eGFR (MDRD, modified for Chinese patients) acts linearly on CL/F around ",
    "57 mL/min/1.73 m^2; serum albumin acts as a power function on V2/F ",
    "around 40 g/L (exponent -7.31). The source fitted doses and ",
    "concentrations in molar units, so the MMF-to-MPA molecular-weight ratio ",
    "320.3/433.5 is carried inside the bioavailability: doses are given as mg ",
    "of MMF and Cc is MPA in mg/L. IIV is exponential on ka, CL/F, V2/F, Q/F, ",
    "V3/F, lag time and F; residual error is combined proportional (26.1%) ",
    "plus additive (0.144 mg/L)."
  )
  reference <- paste(
    "Wang X, Wu Y, Huang J, Shan S, Mai M, Zhu J, Yang M, Shang D, Wu Z,",
    "Lan J, Zhong S, Wu M. (2021). Estimation of Mycophenolic Acid Exposure",
    "in Heart Transplant Recipients by Population Pharmacokinetic and",
    "Limited Sampling Strategies. Front Pharmacol 12:748609.",
    "doi:10.3389/fphar.2021.748609",
    sep = " "
  )
  vignette <- "Wang_2021_mycophenolic_acid"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The MMF -> MPA molecular-weight ratio is applied in f(depot), so every
  # state holds an MPA-equivalent amount in mg (Wang 2021 Methods,
  # 'Pharmacokinetic Modeling': "The units of MMF doses and MPA
  # concentrations were unified as moles").
  compartmentData <- list(
    depot = list(analyte = "mycophenolic acid", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "mycophenolic acid", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "mycophenolic acid", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CONMED_PPI = list(
      description = "Concomitant proton-pump inhibitor (omeprazole or pantoprazole, intravenous or oral) during the PK sampling cycle; 1 = yes, 0 = no.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no PPI co-medication)",
      notes = "Multiplicative effect on F: 0.724^CONMED_PPI (Table 2 footnote 'F = TVF x theta_PPI-F, where PPI equals one if co-medication of PPI during PK sampling'). Recorded per PK sampling cycle, so it may change between occasions in the same patient. 45 of 91 patients (49.5%) received a PPI (Table 1). Route of PPI administration had no detectable effect (Discussion).",
      source_name = "PPI"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate by the MDRD equation modified for Chinese patients (Ma 2006): eGFR = 175 x Scr^-1.234 x Age^-0.179 x 0.79 if female, Scr in mg/dL.",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear effect on CL/F centred at 57 mL/min/1.73 m^2 (cohort median 57.2, range 6.3-197.1; Table 1): CL/F = TVCL x (1 + (eGFR - 57) x 0.00791) (Table 2 footnote). BSA-normalised, as the MDRD equation returns it.",
      source_name = "eGFR"
    ),
    ALB = list(
      description = "Serum albumin.",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on V2/F centred at 40 g/L (cohort median 40.50, range 28.56-57.90; Table 1): V2/F = TVV2 x (ALB/40)^-7.31 (Table 2 footnote; Equation 3).",
      source_name = "ALB"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight.",
      units = "kg",
      type = "continuous",
      notes = "Significant alone on CL/F (Supplementary Table S1 model 4), V2/F (model 12) and Q/F (model 15), but not retained once eGFR, albumin and PPI were in the model."
    ),
    CONMED_DIURETIC = list(
      description = "Concomitant diuretic use (1 = yes, 0 = no).",
      units = "(binary)",
      type = "binary",
      notes = "Significant alone on CL/F (Supplementary Table S1 model 18, dOFV -18.4) but dropped when PPI on F was added (model 19, dOFV +0.7); diuretic and PPI use overlapped in most patients (Discussion)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 91L,
    n_studies = 1L,
    n_observations = "508 plasma MPA samples (507 above the LLOQ of 0.1 mg/L) from 105 PK sampling cycles; 77 patients contributed one cycle, 13 two and one three. 14 patients had intensive 12-h profiles (pre-dose, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12 h); the rest were sparse (0.5, 1.5, 4, 9 h).",
    age_range = "21-74 years",
    age_median = "50 years",
    weight_range = "33.4-95.0 kg",
    weight_median = "60.0 kg",
    height_median = "1.65 m (range 1.51-1.78)",
    sex_female_pct = 7.7,
    race_ethnicity = c(Asian = 100),
    disease_state = "First heart transplantation, at least 7 days post-operation (post-operative time 7-1067 days, median 37).",
    dose_range = "Oral MMF dispersible tablets (Cycopin) 250-750 mg twice daily, median 500 mg.",
    co_medication = "Tacrolimus 1-8 mg/day (trough 2.2-30.0 ng/mL), methylprednisolone 8-20 mg/day, diuretics in 50.5%, proton-pump inhibitors in 49.5%.",
    renal_function = "MDRD (Chinese) eGFR median 57.2 mL/min/1.73 m^2 (range 6.3-197.1); serum creatinine 1.21 mg/dL (0.35-6.55).",
    albumin = "Serum albumin median 40.50 g/L (range 28.56-57.90).",
    assay = "EMIT (Siemens Viva-E), linear range 0.1-15 mg/L; EMIT reads higher than HPLC for MPA (Discussion).",
    regions = "China (Guangdong Provincial People's Hospital, Guangzhou)",
    notes = "Single-centre, prospective, open-label observational study (ChiCTR2000030903). Demographics from Table 1."
  )

  ini({
    lka <- log(0.781); label("Absorption rate constant ka (1/h)") # Table 2 'Ka, 1/h' 0.781
    lcl <- log(7.36); label("Apparent clearance CL/F at eGFR 57 mL/min/1.73 m^2 (L/h)") # Table 2 'CL/F, L/h' 7.36
    lvc <- log(5.69); label("Apparent central volume V2/F at albumin 40 g/L (L)") # Table 2 'V2/F, L' 5.69
    lq <- log(17.0); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 2 'Q/F, L/h' 17.0
    lvp <- log(560); label("Apparent peripheral volume V3/F (L)") # Table 2 'V3/F, L' 560
    ltlag <- log(0.408); label("Absorption lag time (h)") # Table 2 'Tlag, h' 0.408
    lfdepot <- fixed(log(0.95)); label("Oral bioavailability of MMF as MPA, before the molecular-weight ratio (unitless)") # Table 2 'F 0.95 FIX'; Methods cites Bullingham 1996 and Armstrong 2005

    e_conmed_ppi_f <- 0.724; label("Multiplicative factor on F with PPI co-medication (applied as factor^CONMED_PPI, unitless)") # Table 2 'theta PPI-F' 0.724
    e_crcl_cl <- 0.00791; label("Linear slope of eGFR on CL/F (per mL/min/1.73 m^2)") # Table 2 'theta eGFR-CL' 0.00791
    e_alb_vc <- -7.31; label("Power exponent of albumin on V2/F (unitless)") # Table 2 'theta ALB-V2' -7.31

    # IIV: Table 2 prints each as a percentage; converted with
    # omega^2 = log(1 + CV^2). No covariances are reported, so the matrix is
    # diagonal.
    etalka ~ 0.158904 # Table 2 'IIV of Ka' 41.5%
    etalcl ~ 0.156785 # Table 2 'IIV of CL/F' 41.2%
    etalvc ~ 1.499227 # Table 2 'IIV of V2/F' 186.5%
    etalq ~ 0.105161 # Table 2 'IIV of Q/F' 33.3%
    etalvp ~ 1.524103 # Table 2 'IIV of V3/F' 189.5%
    etaltlag ~ 0.017015 # Table 2 'IIV of Tlag' 13.1%
    etalfdepot ~ 0.047686 # Table 2 'IIV of F' 22.1%

    propSd <- 0.261; label("Proportional residual error (fraction)") # Table 2 'Prop. res. error, %' 26.1
    addSd <- 0.144; label("Additive residual error (mg/L)") # Table 2 'Add. res. error, mg/L' 0.144
  })

  model({
    # Individual parameters (Equation 1 exponential IIV; covariate forms from
    # the Table 2 footnote and Equations 2-4).
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (1 + (CRCL - 57) * e_crcl_cl)
    vc <- exp(lvc + etalvc) * (ALB / 40)^e_alb_vc
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)
    tlag <- exp(ltlag + etaltlag)

    # The source fitted MMF doses and MPA concentrations in molar units
    # (Methods; the Supplementary Figure S1 pcVPC axis is in mol/L). Doses here
    # are mg of MMF, so the molecular-weight ratio MPA/MMF = 320.3/433.5
    # converts them to mg of MPA and Cc comes out in mg/L.
    mw_ratio <- 320.3 / 433.5
    fdepot <- exp(lfdepot + etalfdepot) * e_conmed_ppi_f^CONMED_PPI * mw_ratio

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fdepot
    alag(depot) <- tlag

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
