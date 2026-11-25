Sou_2019_SEN032_rat <- function() {
  description <- paste(
    "Preclinical (rat, Sprague-Dawley male). Population PK model for the",
    "PqsR-antagonist quorum-sensing inhibitor SEN032 (a structural analogue",
    "of SEN001) after intratracheal administration of a solubilised dose,",
    "fit simultaneously to plasma and lung-homogenate concentrations. The",
    "dose enters a central lung compartment (volume fixed to the 2 g rat",
    "lung) that exchanges with a deep lung compartment and, through a",
    "bidirectional lung-plasma clearance, with a two-compartment systemic",
    "disposition eliminated from plasma (the full four-compartment model).",
    "Parameters are for a typical 0.3 kg rat.",
    sep = " "
  )
  reference <- paste(
    "Sou T, Kukavica-Ibrulj I, Soukarieh F, Halliday N, Levesque RC,",
    "Williams P, Stocks M, Camara M, Friberg LE, Bergstrom CAS.",
    "Model-Based Drug Development in Pulmonary Delivery: Pharmacokinetic",
    "Analysis of Novel Drug Candidates for Treatment of Pseudomonas",
    "aeruginosa Lung Infection. J Pharm Sci. 2019;108(1):630-640.",
    "doi:10.1016/j.xphs.2018.09.017. PMCID: PMC6307981.",
    "Structural ODEs: Equations 1-4 and Figure 2. Parameter estimates:",
    "Table 5 (SEN032 column).",
    sep = " "
  )
  vignette <- "Sou_2019_quorumSensingInhibitors"
  units <- list(time = "h", dosing = "ug", concentration = "ug/L")
  # Intratracheal dose into the central lung compartment
  dosing <- "lung"

  covariateData <- list()

  compartmentData <- list(
    lung = list(analyte = "SEN032", units = "ug", specimen = "tissue", verified = TRUE),
    lung_deep = list(analyte = "SEN032", units = "ug", specimen = "tissue", verified = TRUE),
    central = list(analyte = "SEN032", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "SEN032", units = "ug", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "rat (Sprague-Dawley)",
    n_subjects = 18L,
    n_studies = 1L,
    sex_female_pct = 0,
    weight_range = "not reported; parameters scaled to a typical 0.3 kg rat",
    disease_state = "healthy (no infection)",
    dose_range = "0.28 mg/kg single intratracheal dose of a solution in PEG400 2% + polysorbate 80 2% w/w",
    regions = "United Kingdom (XenoGesis / Saretius, Nottingham)",
    notes = paste(
      "Table 2: 18 rats, terminal sampling of plasma and lung (no BAL) at",
      "0.167, 0.5, 1, 3, 6 and 24 h. Data below the LLOQ handled with the",
      "M3 method. IIV not estimated because most samples were terminal.",
      sep = " "
    )
  )

  ini({
    lq_lung <- log(0.150); label("Bidirectional lung-plasma clearance CLLP (L/h)") # Table 5: CLLP = 0.150 L/h (RSE 10%)
    lv_lung <- fixed(log(0.002)); label("Central lung volume VL (L)") # Table 5: VL = 0.002 L, no RSE; Methods: fixed to the 2 g rat lung
    lcl <- log(3.74); label("Central (systemic) clearance CLC (L/h)") # Table 5: CLC = 3.74 L/h (RSE 6%)
    lvc <- log(1.51); label("Central plasma volume VC (L)") # Table 5: VC = 1.51 L (RSE 12%)
    lq <- log(2.35); label("Plasma-peripheral intercompartmental clearance QC (L/h)") # Table 5: QC = 2.35 L/h (RSE 25%)
    lvp <- log(1.36); label("Peripheral tissue volume VP (L)") # Table 5: VP = 1.36 L (RSE 35%)
    lq_lung_lung_deep <- log(0.0107); label("Lung to deep-lung intercompartmental clearance QL (L/h)") # Table 5: QL = 0.0107 L/h (RSE 17%)
    lv_lung_deep <- log(0.00360); label("Deep lung volume VLP (L)") # Table 5: VLP = 0.00360 L (RSE 15%)

    propSd <- 0.253; label("Proportional residual error, plasma (fraction)") # Table 5: PROPPL = 25.3% (RSE 15%)
    addSd <- 0.711; label("Additive residual error, plasma (ug/L)") # Table 5: ADDPL = 0.505 reported as a variance (footnote a); SD = sqrt(0.505)
    propSd_Clung <- 0.461; label("Proportional residual error, lung homogenate (fraction)") # Table 5: PROPLG = 46.1% (RSE 30%)
  })

  model({
    q_lung <- exp(lq_lung)
    v_lung <- exp(lv_lung)
    cl <- exp(lcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)
    q_lung_lung_deep <- exp(lq_lung_lung_deep)
    v_lung_deep <- exp(lv_lung_deep)

    # Equations 1-4; dose is given into `lung`.
    d/dt(lung) <- -q_lung / v_lung * lung + q_lung / vc * central -
      q_lung_lung_deep / v_lung * lung + q_lung_lung_deep / v_lung_deep * lung_deep
    d/dt(central) <- q_lung / v_lung * lung - q_lung / vc * central - cl / vc * central -
      q / vc * central + q / vp * peripheral1
    d/dt(peripheral1) <- q / vc * central - q / vp * peripheral1
    d/dt(lung_deep) <- q_lung_lung_deep / v_lung * lung - q_lung_lung_deep / v_lung_deep * lung_deep

    # ug / L (numerically ng/mL)
    Cc <- central / vc
    Clung <- lung / v_lung

    Cc ~ add(addSd) + prop(propSd)
    Clung ~ prop(propSd_Clung)
  })
}
