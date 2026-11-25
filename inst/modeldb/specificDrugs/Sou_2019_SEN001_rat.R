Sou_2019_SEN001_rat <- function() {
  description <- paste(
    "Preclinical (rat, Sprague-Dawley male). Population PK model for the",
    "PqsR-antagonist quorum-sensing inhibitor SEN001 after intratracheal",
    "administration of a solubilised dose, fit simultaneously to plasma,",
    "lung-homogenate and bronchoalveolar-lavage (BAL) concentrations. The",
    "dose enters a central lung compartment (volume fixed to the 2 g rat",
    "lung) that exchanges with a deep lung compartment and, through a",
    "bidirectional lung-plasma clearance, with a one-compartment systemic",
    "disposition eliminated from plasma. BAL recovers a fraction FBAL of",
    "the central-lung amount, and the lung homogenate (taken after lavage)",
    "holds the remainder. Parameters are for a typical 0.3 kg rat.",
    sep = " "
  )
  reference <- paste(
    "Sou T, Kukavica-Ibrulj I, Soukarieh F, Halliday N, Levesque RC,",
    "Williams P, Stocks M, Camara M, Friberg LE, Bergstrom CAS.",
    "Model-Based Drug Development in Pulmonary Delivery: Pharmacokinetic",
    "Analysis of Novel Drug Candidates for Treatment of Pseudomonas",
    "aeruginosa Lung Infection. J Pharm Sci. 2019;108(1):630-640.",
    "doi:10.1016/j.xphs.2018.09.017. PMCID: PMC6307981.",
    "Structural ODEs: Equations 1-4 and Figure 2. BAL observation:",
    "Equation 5. Parameter estimates: Table 5 (SEN001 column).",
    sep = " "
  )
  vignette <- "Sou_2019_quorumSensingInhibitors"
  units <- list(time = "h", dosing = "ug", concentration = "ug/L")
  # Intratracheal dose into the central lung compartment
  dosing <- "lung"

  covariateData <- list()

  compartmentData <- list(
    lung = list(analyte = "SEN001", units = "ug", specimen = "tissue", verified = TRUE),
    lung_deep = list(analyte = "SEN001", units = "ug", specimen = "tissue", verified = TRUE),
    central = list(analyte = "SEN001", units = "ug", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "rat (Sprague-Dawley)",
    n_subjects = 21L,
    n_studies = 1L,
    sex_female_pct = 0,
    weight_range = "not reported; parameters scaled to a typical 0.3 kg rat",
    disease_state = "healthy (no infection)",
    dose_range = "0.067 mg/kg single intratracheal dose, 250 uL of solution in PEG400 2% + polysorbate 80 2% w/w",
    regions = "Canada (Laval University, Quebec)",
    notes = paste(
      "Table 2: 21 rats, terminal sampling of plasma, lung and BAL at",
      "0.5, 2, 4, 8, 12, 24 and 48 h. Data below the LLOQ handled with",
      "the M3 method. IIV not estimated because most samples were",
      "terminal (one sample per animal per matrix).",
      sep = " "
    )
  )

  ini({
    lq_lung <- log(0.0175); label("Bidirectional lung-plasma clearance CLLP (L/h)") # Table 5: CLLP = 0.0175 L/h (RSE 7%)
    lv_lung <- fixed(log(0.002)); label("Central lung volume VL (L)") # Table 5: VL = 0.002 L, no RSE; Methods: fixed to the 2 g rat lung
    lcl <- log(0.0745); label("Central (systemic) clearance CLC (L/h)") # Table 5: CLC = 0.0745 L/h (RSE 10%)
    lvc <- log(0.0855); label("Central plasma volume VC (L)") # Table 5: VC = 0.0855 L (RSE 31%)
    lq_lung_lung_deep <- log(0.00596); label("Lung to deep-lung intercompartmental clearance QL (L/h)") # Table 5: QL = 0.00596 L/h (RSE 37%)
    lv_lung_deep <- log(0.0582); label("Deep lung volume VLP (L)") # Table 5: VLP = 0.0582 L (RSE 19%)
    logitfbal <- log(0.980 / (1 - 0.980)); label("Fraction of the central-lung amount recovered in BAL, FBAL (logit)") # Table 5: FBAL = 0.980 (RSE 0.3%)
    lv_bal <- fixed(log(0.008)); label("Volume of lavage fluid used to collect BAL, VBAL (L)") # Equation 5; value not printed -- back-solved from Figure 3 (BAL/lung median ratio), see vignette

    propSd <- 0.459; label("Proportional residual error, plasma (fraction)") # Table 5: PROPPL = 45.9% (RSE 12%)
    propSd_Clung_postbal <- 0.435; label("Proportional residual error, lung homogenate (fraction)") # Table 5: PROPLG = 43.5% (RSE 13%)
    propSd_Cbal <- 0.380; label("Proportional residual error, BAL (fraction)") # Table 5: PROPBAL = 38.0% (RSE 22%)
  })

  model({
    q_lung <- exp(lq_lung)
    v_lung <- exp(lv_lung)
    cl <- exp(lcl)
    vc <- exp(lvc)
    q_lung_lung_deep <- exp(lq_lung_lung_deep)
    v_lung_deep <- exp(lv_lung_deep)
    fbal <- expit(logitfbal)
    v_bal <- exp(lv_bal)

    # Equations 1, 2 and 4; dose is given into `lung`.
    d/dt(lung) <- -q_lung / v_lung * lung + q_lung / vc * central -
      q_lung_lung_deep / v_lung * lung + q_lung_lung_deep / v_lung_deep * lung_deep
    d/dt(central) <- q_lung / v_lung * lung - q_lung / vc * central - cl / vc * central
    d/dt(lung_deep) <- q_lung_lung_deep / v_lung * lung - q_lung_lung_deep / v_lung_deep * lung_deep

    # ug / L (numerically ng/mL)
    Cc <- central / vc
    # Total central-lung concentration before lavage (the quantity of Figure 4)
    Clung <- lung / v_lung
    # Lung homogenate collected after lavage holds the amount BAL did not recover
    Clung_postbal <- (1 - fbal) * lung / v_lung
    # Equation 5
    Cbal <- fbal / v_bal * lung

    Cc ~ prop(propSd)
    Clung_postbal ~ prop(propSd_Clung_postbal)
    Cbal ~ prop(propSd_Cbal)
  })
}
