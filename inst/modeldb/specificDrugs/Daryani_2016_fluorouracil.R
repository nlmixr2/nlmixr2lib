Daryani_2016_fluorouracil <- function() {
  description <- paste(
    "Population PK model of intravenous bolus 5-fluorouracil (5-FU) in",
    "children and young adults with recurrent or refractory ependymoma",
    "(Daryani 2016, phase I trial). Two-compartment model with first-order",
    "elimination; every clearance and volume is proportional to body surface",
    "area (the authors fitted BSA-normalised doses), with IIV on all four",
    "parameters and inter-occasion variability on clearance. The paper",
    "linked this model to its mouse tumor extracellular fluid (tECF) model to",
    "predict tumor exposure; that link is carried here as the simulated-only",
    "output Cecf, with the mouse tECF rate constants and a human plasma",
    "fraction unbound of 0.9. Dose in mg; concentrations in uM."
  )
  reference <- paste(
    "Daryani VM, Patel YT, Tagen M, Turner DC, Carcaboso AM, Atkinson JM,",
    "Gajjar A, Gilbertson RJ, Wright KD, Stewart CF. Translational",
    "Pharmacokinetic-Pharmacodynamic Modeling and Simulation: Optimizing",
    "5-Fluorouracil Dosing in Children With Pediatric Ependymoma.",
    "CPT Pharmacometrics Syst Pharmacol. 2016;5(4):211-221.",
    "doi:10.1002/psp4.12075.",
    "Structure, occasion handling and residual-error form are from the",
    "authors' NONMEM control stream (Supplementary Material 9)."
  )
  vignette <- "Daryani_2016_fluorouracil"
  units <- list(time = "h", dosing = "mg", concentration = "umol/L")

  compartmentData <- list(
    central = list(analyte = "5-fluorouracil", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "5-fluorouracil", units = "mg", specimen = "tissue", verified = TRUE),
    ecf = list(
      analyte = "5-fluorouracil (unbound, tumor extracellular fluid; predicted from the mouse model)",
      units = "mg",
      specimen = "tumor",
      verified = TRUE
    ),
    int_tumor = list(
      analyte = "5-fluorouracil (intracellular tumor; predicted from the mouse model)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    )
  )

  covariateData <- list(
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NA,
      notes = paste(
        "Every clearance and volume is linear in BSA (exponent 1, no",
        "centring): the authors fitted doses normalised to BSA, so the Table 2",
        "values are per m^2. The BSA formula is not stated."
      ),
      source_name = "BSA"
    ),
    OCC = list(
      description = "Treatment occasion: one bolus dose and the concentrations that follow it",
      units = "(count)",
      type = "categorical",
      reference_category = NA,
      notes = paste(
        "Occasions 1, 2 and 3 each have their own inter-occasion eta on CL;",
        "every occasion after the third shares a fourth eta",
        "(control stream: IF (OCC.GT.3) OCL = DEXP(ETA(8)))."
      ),
      source_name = "OCC"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 26L,
    n_studies = 1L,
    age_range = "3.8-21.1 years (median 7.1)",
    weight_range = "10.8-83.0 kg (median 25.5)",
    sex_female_pct = 69.2,
    race_ethnicity = c(White = 73.1, Black = 19.2, Other = 7.7),
    disease_state = "Recurrent or refractory pediatric ependymoma",
    dose_range = "400, 500 or 650 mg/m^2 5-FU intravenous bolus",
    regions = "USA (St. Jude Children's Research Hospital phase I trial)",
    bsa_range = "0.6-2.0 m^2 (median 0.9)",
    height_range = "96.5-173.0 cm (median 118.0)",
    n_observations = 570L,
    notes = paste(
      "Demographics from Supplementary Table 2. About 40% of the 570",
      "concentrations were below the 15 uM limit of quantitation and were",
      "handled with the M3 likelihood method (Laplacian + Monte Carlo",
      "importance sampling, NONMEM 7.3). No covariate beyond the a priori BSA",
      "normalisation was retained."
    )
  )

  ini({
    # Plasma PK (Table 2), per m^2 of BSA.
    lcl <- log(16.6); label("Clearance per m^2 BSA (L/h/m^2)") # Table 2 CLp = 16.6 +/- 1.2
    lvc <- log(4.2); label("Central volume per m^2 BSA (L/m^2)") # Table 2 V1,p = 4.2 +/- 0.3
    lq <- log(1.6); label("Intercompartmental clearance per m^2 BSA (L/h/m^2)") # Table 2 Qp = 1.6 +/- 0.7
    lvp <- log(20.1); label("Peripheral volume per m^2 BSA (L/m^2)") # Table 2 V2,p = 20.1 +/- 8.6

    # Tumor ECF link: the mouse tECF parameters of Table 1 (as held in the
    # PK-PD control stream, Supplementary Material 9) with the human fu.
    fu <- fixed(0.9); label("Fraction unbound in human plasma (fraction)") # Results 'Application of scaled-adult model': 5-FU plasma fraction unbound 0.9 (ref 23)
    lvecf <- fixed(log(0.001)); label("Tumor ECF volume per m^2 BSA (L/m^2)") # Table 1 V3,m = 0.001, used on the same per-m^2 scale as V1,p
    lkin_ecf <- fixed(log(0.00426)); label("Rate constant from plasma (unbound) into tumor ECF (1/h)") # Supplementary Material 9 K13 = 0.00426; Table 1 K13,m = 0.0043
    lkout_ecf <- fixed(log(3.24)); label("Rate constant from tumor ECF back to plasma (1/h)") # Table 1 K31,m = 3.24
    lkin_int_tumor <- fixed(log(1.82)); label("Rate constant from tumor ECF into intracellular tumor (1/h)") # Table 1 K34,m = 1.82
    lkout_int_tumor <- fixed(log(0.334)); label("Rate constant from intracellular tumor back to tumor ECF (1/h)") # Table 1 K43,m = 0.334
    lvmax_ecf <- fixed(log(0.00633)); label("Maximum elimination rate from tumor ECF per m^2 BSA (umol/h/m^2)") # Supplementary Material 9 K30 = 0.00633; Table 1 Vmax,tumor,m = 0.0063
    lkm_ecf <- fixed(log(0.0107)); label("Michaelis-Menten constant for tumor ECF elimination (umol/L)") # Supplementary Material 9 KM2 = 0.0107; Table 1 prints 0.012

    # IIV and IOV: exponential, Table 2 % CV converted as omega^2 = log(1 + CV^2).
    etalcl ~ 0.07968 # Table 2 IIV CLp 28.8% CV
    etalvc ~ 0.06836 # Table 2 IIV V1,p 26.6% CV
    etalq ~ 0.4957 # Table 2 IIV Qp 80.1% CV
    etalvp ~ 0.6322 # Table 2 IIV V2,p 93.9% CV
    etaiov_cl_1 ~ 0.01676 # Table 2 IOV CLp 13.0% CV, occasion 1
    etaiov_cl_2 ~ fixed(0.01676) # same variance as occasion 1 (control stream OMEGA BLOCK(1) SAME)
    etaiov_cl_3 ~ fixed(0.01676) # same variance as occasion 1 (control stream OMEGA BLOCK(1) SAME)
    etaiov_cl_4 ~ fixed(0.01676) # same variance as occasion 1 (control stream OMEGA BLOCK(1) SAME), occasions 4 and later

    # Residual error: SD = sqrt(add^2 + (prop * F)^2) (control stream $ERROR).
    addSd <- 1.98; label("Additive residual error (umol/L)") # Table 2 additive 1.98 +/- 1.5
    propSd <- 0.51; label("Proportional residual error (fraction)") # Table 2 proportional 0.51 +/- 0.04
  })

  model({
    # 5-FU molecular weight (g/mol), for the mg -> umol conversion.
    mw <- 130.08

    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC >= 4)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 + oc4 * etaiov_cl_4

    cl <- exp(lcl + etalcl + iov_cl) * BSA
    vc <- exp(lvc + etalvc) * BSA
    q <- exp(lq + etalq) * BSA
    vp <- exp(lvp + etalvp) * BSA

    vecf <- exp(lvecf) * BSA
    kin_ecf <- exp(lkin_ecf)
    kout_ecf <- exp(lkout_ecf)
    kin_int_tumor <- exp(lkin_int_tumor)
    kout_int_tumor <- exp(lkout_int_tumor)
    vmax_ecf <- exp(lvmax_ecf) * BSA
    km_ecf <- exp(lkm_ecf)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    Cc <- central / vc * 1000 / mw
    Cecf <- ecf / vecf * 1000 / mw

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(ecf) <- kin_ecf * fu * central - kout_ecf * ecf - kin_int_tumor * ecf + kout_int_tumor * int_tumor -
      vmax_ecf * mw / 1000 * Cecf / (km_ecf + Cecf)
    d/dt(int_tumor) <- kin_int_tumor * ecf - kout_int_tumor * int_tumor

    Cc ~ add(addSd) + prop(propSd) + combined2()
  })
}
