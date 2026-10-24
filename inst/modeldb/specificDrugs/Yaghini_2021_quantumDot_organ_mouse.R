Yaghini_2021_quantumDot_organ_mouse <- function() {
  description <- paste(
    "Preclinical (mouse, female Balb/c). Three-compartment population PK model of",
    "blood, liver and spleen disposition of indium-based, cadmium-free PEGylated",
    "quantum dots (bio CFQD, InP core) after a single IV bolus, fitted to indium",
    "concentrations over 90 days. Blood is the central compartment with first-order",
    "elimination; liver and spleen exchange with blood through intercompartmental",
    "clearances. No between-animal variability (destructive sampling). The printed",
    "parameters do not reproduce the paper's own Figure 8 blood profile and the",
    "residual-error magnitudes are not reported; see the vignette Errata. Companion",
    "0-24 h blood model: modellib('Yaghini_2021_quantumDot_blood_mouse')."
  )
  reference <- paste(
    "Yaghini E, Tacconi E, Pilling A, Rahman P, Broughton J, Naasani I, Keshtgar MRS,",
    "MacRobert A, Della Pasqua O. (2021). Population pharmacokinetic modelling of",
    "indium-based quantum dot nanoparticles: preclinical in vivo studies.",
    "Eur J Pharm Sci 157:105639. doi:10.1016/j.ejps.2020.105639."
  )
  vignette <- "Yaghini_2021_quantumDot"
  units <- list(time = "h", dosing = "ug", concentration = "ug/L")

  covariateData <- list()

  compartmentData <- list(
    central = list(analyte = "indium (quantum dot)", units = "ug", specimen = "serum", verified = TRUE),
    liver = list(analyte = "indium (quantum dot)", units = "ug", specimen = "tissue", verified = TRUE),
    spleen = list(analyte = "indium (quantum dot)", units = "ug", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "mouse (female Balb/c)",
    n_subjects = 30L,
    n_studies = 1L,
    sex_female_pct = 100,
    disease_state = "healthy",
    dose_range = "single IV bolus of 20 mg/kg PEGylated bio CFQD quantum dots (200 uL in PBS via the tail vein)",
    notes = paste(
      "The analysis dataset is not tabulated. Nine biodistribution groups of 3 mice",
      "were sacrificed at 5 min, 1 h, 4 h and 1, 3, 10, 30, 60 and 90 days; Figure 8",
      "also carries the 20- and 40-min serum samples of a separate set of 3 serially",
      "sampled mice, so 30 mice is inferred. Indium was quantified by ICP-MS in serum",
      "and in nitric-acid digests of 0.1 g tissue. The indium mass fraction of the",
      "quantum dot is not reported, so the dose units are not stated (see vignette",
      "Errata). Body weight is not reported."
    )
  )

  ini({
    lcl <- log(0.0021); label("Clearance from blood (L/h)") # Table 2: CL = 0.0021 L/h (29% CV)
    lvc <- log(3.41); label("Blood volume of distribution Vb (L)") # Table 2: Vb = 3.41 L (33% CV)
    lq <- log(0.0018); label("Blood-liver intercompartmental clearance QT1 (L/h)") # Table 2: QT1 = 0.0018 L/h (35% CV)
    lvp <- log(0.13); label("Liver volume of distribution VT1 (L)") # Table 2: VT1 = 0.13 L (34% CV)
    lq2 <- log(0.1606); label("Blood-spleen intercompartmental clearance QT2 (L/h)") # Table 2: QT2 = 0.1606 L/h (23% CV)
    lvp2 <- log(0.22); label("Spleen volume of distribution VT2 (L)") # Table 2: VT2 = 0.22 L (34% CV)

    # Results: residual variability 'was characterised by a proportional and
    # additive error model'; Table 2 prints no residual-error values, so every
    # magnitude is a zero placeholder.
    propSd <- fixed(0); label("Proportional residual error on blood (fraction; not reported)") # not reported in Table 2
    addSd <- fixed(0); label("Additive residual error on blood (ug/L; not reported)") # not reported in Table 2
    propSd_Cliver <- fixed(0); label("Proportional residual error on liver (fraction; not reported)") # not reported in Table 2
    addSd_Cliver <- fixed(0); label("Additive residual error on liver (ug/L; not reported)") # not reported in Table 2
    propSd_Cspleen <- fixed(0); label("Proportional residual error on spleen (fraction; not reported)") # not reported in Table 2
    addSd_Cspleen <- fixed(0); label("Additive residual error on spleen (ug/L; not reported)") # not reported in Table 2
  })
  model({
    cl <- exp(lcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)
    q2 <- exp(lq2)
    vp2 <- exp(lvp2)

    # Eqs (6)-(10): rate constants from flows and volumes. Table 2 prints the
    # resulting K12 0.0005, K21 0.014, K13 0.047, K31 0.73 and ke 0.0006 1/h.
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2
    kel <- cl / vc

    # Eqs (3)-(5) are printed on concentrations C1-C3 but with k = Q/V rate
    # constants, which is the mass-balanced system only when the states are
    # amounts (NONMEM ADVAN6 integrates amounts). Compartment 2 is the liver
    # (VT1) and compartment 3 the spleen (VT2), per eqs (4)-(5) and Figure 1B.
    d/dt(central) <- k21 * liver + k31 * spleen - (k12 + k13 + kel) * central
    d/dt(liver) <- k12 * central - k21 * liver
    d/dt(spleen) <- k13 * central - k31 * spleen

    Cc <- central / vc
    Cliver <- liver / vp
    Cspleen <- spleen / vp2

    Cc ~ add(addSd) + prop(propSd)
    Cliver ~ add(addSd_Cliver) + prop(propSd_Cliver)
    Cspleen ~ add(addSd_Cspleen) + prop(propSd_Cspleen)
  })
}
