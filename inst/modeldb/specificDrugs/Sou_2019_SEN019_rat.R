Sou_2019_SEN019_rat <- function() {
  description <- paste(
    "Preclinical (rat, Sprague-Dawley male). Population PK model for the",
    "PqsR-antagonist quorum-sensing inhibitor SEN019 (a structural analogue",
    "of SEN001) after intratracheal administration of a solubilised dose,",
    "fit simultaneously to plasma and lung-homogenate concentrations. The",
    "dose enters a central lung compartment (volume fixed to the 2 g rat",
    "lung) that exchanges through a bidirectional lung-plasma clearance with",
    "a one-compartment systemic disposition eliminated from plasma; neither",
    "a deep lung nor a peripheral tissue compartment was supported. Parameters",
    "are for a typical 0.3 kg rat.",
    sep = " "
  )
  reference <- paste(
    "Sou T, Kukavica-Ibrulj I, Soukarieh F, Halliday N, Levesque RC,",
    "Williams P, Stocks M, Camara M, Friberg LE, Bergstrom CAS.",
    "Model-Based Drug Development in Pulmonary Delivery: Pharmacokinetic",
    "Analysis of Novel Drug Candidates for Treatment of Pseudomonas",
    "aeruginosa Lung Infection. J Pharm Sci. 2019;108(1):630-640.",
    "doi:10.1016/j.xphs.2018.09.017. PMCID: PMC6307981.",
    "Structural ODEs: Equations 1-2 (QL = QC = 0) and Figure 2. Parameter",
    "estimates: Table 5 (SEN019 column).",
    sep = " "
  )
  vignette <- "Sou_2019_quorumSensingInhibitors"
  units <- list(time = "h", dosing = "ug", concentration = "ug/L")
  # Intratracheal dose into the central lung compartment
  dosing <- "lung"

  covariateData <- list()

  compartmentData <- list(
    lung = list(analyte = "SEN019", units = "ug", specimen = "tissue", verified = TRUE),
    central = list(analyte = "SEN019", units = "ug", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "rat (Sprague-Dawley)",
    n_subjects = 18L,
    n_studies = 1L,
    sex_female_pct = 0,
    weight_range = "not reported; parameters scaled to a typical 0.3 kg rat",
    disease_state = "healthy (no infection)",
    dose_range = "0.12 mg/kg single intratracheal dose of a solution in PEG400 2% + polysorbate 80 2% w/w",
    regions = "United Kingdom (XenoGesis / Saretius, Nottingham)",
    notes = paste(
      "Table 2: 18 rats, terminal sampling of plasma and lung (no BAL) at",
      "0.167, 0.5, 1, 3, 6 and 24 h. Data below the LLOQ handled with the",
      "M3 method. IIV not estimated because most samples were terminal.",
      sep = " "
    )
  )

  ini({
    lq_lung <- log(0.0724); label("Bidirectional lung-plasma clearance CLLP (L/h)") # Table 5: CLLP = 0.0724 L/h (RSE 1%)
    lv_lung <- fixed(log(0.002)); label("Central lung volume VL (L)") # Table 5: VL = 0.002 L, no RSE; Methods: fixed to the 2 g rat lung
    lcl <- log(1.11); label("Central (systemic) clearance CLC (L/h)") # Table 5: CLC = 1.11 L/h (RSE 5%)
    lvc <- log(0.583); label("Central plasma volume VC (L)") # Table 5: VC = 0.583 L (RSE 5%)

    propSd <- 0.248; label("Proportional residual error, plasma (fraction)") # Table 5: PROPPL = 24.8% (RSE 7%)
    addSd <- 0.533; label("Additive residual error, plasma (ug/L)") # Table 5: ADDPL = 0.284 reported as a variance (footnote a); SD = sqrt(0.284)
    propSd_Clung <- 0.181; label("Proportional residual error, lung homogenate (fraction)") # Table 5: PROPLG = 18.1% (RSE 25%)
  })

  model({
    q_lung <- exp(lq_lung)
    v_lung <- exp(lv_lung)
    cl <- exp(lcl)
    vc <- exp(lvc)

    # Equations 1 and 2 without the deep-lung and peripheral terms; dose is
    # given into `lung`.
    d/dt(lung) <- -q_lung / v_lung * lung + q_lung / vc * central
    d/dt(central) <- q_lung / v_lung * lung - q_lung / vc * central - cl / vc * central

    # ug / L (numerically ng/mL)
    Cc <- central / vc
    Clung <- lung / v_lung

    Cc ~ add(addSd) + prop(propSd)
    Clung ~ prop(propSd_Clung)
  })
}
