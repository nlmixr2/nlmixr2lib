Verscheijden_2021_morphine_pbpk <- function() {
  description <- paste(
    "PBPK/PD (whole-body, permeability-limited 4-compartment brain, R/deSolve).",
    "Morphine and its active metabolite morphine-6-glucuronide (M6G) in",
    "adults, children and neonates (Verscheijden et al. 2021, PLoS Comput",
    "Biol). Thirteen perfusion-limited body compartments per analyte",
    "(venous and arterial blood, lung, adipose, bone, heart, kidney, muscle,",
    "skin, spleen, gut, liver, rest of body) with Rodgers and Rowland tissue",
    "partition coefficients computed from age-dependent tissue composition,",
    "plus the Gaohua / Verscheijden brain model (brain blood, brain mass,",
    "cranial CSF, spinal CSF) with passive blood-brain and blood-CSF barrier",
    "permeability and active P-glycoprotein efflux of morphine from brain",
    "mass, scaled from MDCKII-Pgp transwell data by in vitro-in vivo",
    "extrapolation. Organ volumes, blood flows, haematocrit and CSF flows",
    "are age-, sex-, weight- and height-dependent (paediatric equations",
    "below 18 y, adult equations of Verscheijden 2019 from 18 y). Morphine",
    "is cleared from venous blood; a fixed fraction of the cleared mass",
    "(corrected for molecular weight) is formed as M6G. The PD part converts",
    "unbound brain-mass morphine and M6G into competitive mu-opioid receptor",
    "occupancy and a sigmoid relative analgesic response. Forward-simulation",
    "model: all parameters are fixed; IIV reproduces the virtual-population",
    "variability of the deposited code; no residual error was reported."
  )
  reference <- paste(
    "Verscheijden LFM, Litjens CHC, Koenderink JB, Mathijssen RHJ,",
    "Verbeek MM, de Wildt SN, Russel FGM. Physiologically based",
    "pharmacokinetic/pharmacodynamic model for the prediction of morphine",
    "brain disposition and analgesia in adults and children. PLoS Comput",
    "Biol. 2021;17(3):e1008786. doi:10.1371/journal.pcbi.1008786.",
    "Adult and paediatric physiology (S1 Table) and brain framework:",
    "Verscheijden LFM, Koenderink JB, de Wildt SN, Russel FGM. Development",
    "of a physiologically-based pharmacokinetic pediatric brain model for",
    "prediction of cerebrospinal fluid drug concentrations and the",
    "influence of meningitis. PLoS Comput Biol. 2019;15(6):e1007117.",
    "doi:10.1371/journal.pcbi.1007117."
  )
  vignette <- "Verscheijden_2021_morphine"

  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "mg/L"
  )

  compartmentData <- list(
    venous = list(analyte = "morphine", units = "mg", specimen = "whole blood", verified = TRUE),
    arterial = list(analyte = "morphine", units = "mg", specimen = "whole blood", verified = TRUE),
    lung = list(analyte = "morphine", units = "mg", specimen = "tissue", verified = TRUE),
    adipose = list(analyte = "morphine", units = "mg", specimen = "tissue", verified = TRUE),
    bone = list(analyte = "morphine", units = "mg", specimen = "tissue", verified = TRUE),
    heart = list(analyte = "morphine", units = "mg", specimen = "tissue", verified = TRUE),
    kidney = list(analyte = "morphine", units = "mg", specimen = "tissue", verified = TRUE),
    muscle = list(analyte = "morphine", units = "mg", specimen = "tissue", verified = TRUE),
    skin = list(analyte = "morphine", units = "mg", specimen = "tissue", verified = TRUE),
    spleen = list(analyte = "morphine", units = "mg", specimen = "tissue", verified = TRUE),
    other = list(analyte = "morphine", units = "mg", specimen = "tissue", verified = TRUE),
    gut = list(analyte = "morphine", units = "mg", specimen = "tissue", verified = TRUE),
    liver = list(analyte = "morphine", units = "mg", specimen = "tissue", verified = TRUE),
    brain_vascular = list(analyte = "morphine", units = "mg", specimen = "whole blood", verified = TRUE),
    brain = list(analyte = "morphine", units = "mg", specimen = "tissue", verified = TRUE),
    brain_csf_sas_cranial = list(analyte = "morphine", units = "mg", specimen = "CSF", verified = TRUE),
    brain_csf_sas_spinal = list(analyte = "morphine", units = "mg", specimen = "CSF", verified = TRUE),
    venous_m6g = list(analyte = "morphine-6-glucuronide", units = "mg", specimen = "whole blood", verified = TRUE),
    arterial_m6g = list(analyte = "morphine-6-glucuronide", units = "mg", specimen = "whole blood", verified = TRUE),
    lung_m6g = list(analyte = "morphine-6-glucuronide", units = "mg", specimen = "tissue", verified = TRUE),
    adipose_m6g = list(analyte = "morphine-6-glucuronide", units = "mg", specimen = "tissue", verified = TRUE),
    bone_m6g = list(analyte = "morphine-6-glucuronide", units = "mg", specimen = "tissue", verified = TRUE),
    heart_m6g = list(analyte = "morphine-6-glucuronide", units = "mg", specimen = "tissue", verified = TRUE),
    kidney_m6g = list(analyte = "morphine-6-glucuronide", units = "mg", specimen = "tissue", verified = TRUE),
    muscle_m6g = list(analyte = "morphine-6-glucuronide", units = "mg", specimen = "tissue", verified = TRUE),
    skin_m6g = list(analyte = "morphine-6-glucuronide", units = "mg", specimen = "tissue", verified = TRUE),
    spleen_m6g = list(analyte = "morphine-6-glucuronide", units = "mg", specimen = "tissue", verified = TRUE),
    other_m6g = list(analyte = "morphine-6-glucuronide", units = "mg", specimen = "tissue", verified = TRUE),
    gut_m6g = list(analyte = "morphine-6-glucuronide", units = "mg", specimen = "tissue", verified = TRUE),
    liver_m6g = list(analyte = "morphine-6-glucuronide", units = "mg", specimen = "tissue", verified = TRUE),
    brain_vascular_m6g = list(
      analyte = "morphine-6-glucuronide",
      units = "mg",
      specimen = "whole blood",
      verified = TRUE
    ),
    brain_m6g = list(analyte = "morphine-6-glucuronide", units = "mg", specimen = "tissue", verified = TRUE),
    brain_csf_sas_cranial_m6g = list(
      analyte = "morphine-6-glucuronide",
      units = "mg",
      specimen = "CSF",
      verified = TRUE
    ),
    brain_csf_sas_spinal_m6g = list(analyte = "morphine-6-glucuronide", units = "mg", specimen = "CSF", verified = TRUE)
  )

  covariateData <- list(
    AGE = list(
      description = paste(
        "Age. Selects the neonatal (< 0.25 y), paediatric (0.25 to < 18 y)",
        "or adult (>= 18 y) parameter sets and enters the paediatric organ",
        "volume, blood flow, haematocrit and tissue-composition equations",
        "directly. Must be > 0 (several tissue-composition terms use",
        "log10(AGE))."
      ),
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Verscheijden 2021 simulated adults 18-69 y, children 1-18 y and",
        "neonates 1-29 days postnatal age (Table 2, Fig 2B). The < 0.25 y",
        "neonatal cut-off is an assumption of this implementation (the",
        "paper does not define its neonatal age boundary)."
      ),
      source_name = "age"
    ),
    WT = list(
      description = paste(
        "Total body weight. Drives organ volumes, the rest-of-body volume,",
        "morphine and M6G clearance and the paediatric spinal CSF volume."
      ),
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The deposited R code draws weight from sex-specific age-to-height",
        "and height-to-weight equations with log-normal variability",
        "(Verscheijden 2019 S1 Table); supply the individual weight here."
      ),
      source_name = "Weight"
    ),
    HT = list(
      description = paste(
        "Body height. Drives body surface area (Du Bois) and hence cardiac",
        "output and skin volume, and enters the lung, heart, spleen, liver,",
        "adipose, bone, muscle, kidney and blood volume equations."
      ),
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = "Body surface area is computed internally as 0.007184 * HT^0.725 * WT^0.425.",
      source_name = "Height"
    ),
    SEXF = list(
      description = paste(
        "Female sex indicator. Selects the sex-specific organ volume, blood",
        "flow and haematocrit equations and the sex-specific haematocrit",
        "variability."
      ),
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = male",
      notes = paste(
        "1 = female, 0 = male. The deposited code simulates males and",
        "females with separate ODE functions; the 358.5 * gender term of",
        "the infant adipose equation uses gender = 1 for female."
      ),
      source_name = "gender"
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 11L,
    age_range = "10 days to 69 years (neonates, children 1-18 y, adults 18-69 y)",
    weight_range = "not reported; generated from age and height in the virtual populations",
    sex_female_pct = NA_real_,
    disease_state = paste(
      "Mixed. PK verification: neurological / neurosurgical patients",
      "(Meineke 2002), traumatic brain injury (Bouw 2001, Ederoth 2003,",
      "Ketharanathan 2019), paediatric acute leukaemia (Hain 1999),",
      "neonates (Pokela 1993; Radboudumc CSF biobank). PD verification:",
      "healthy adult volunteers (Sarton 2000, Dahan 2004, Skarke 2003),",
      "paediatric cancer patients (Mashayekhi 2009) and preterm neonates",
      "on CPAP (Enders 2008)."
    ),
    dose_range = paste(
      "Single IV morphine 0.025-0.38 mg/kg or 10-28 mg; continuous IV",
      "infusion 0.03 mg/kg/h (children) and 0.25 mg/kg/h for 1 h."
    ),
    regions = "Europe (Netherlands, Germany, Sweden, UK, Austria)",
    notes = paste(
      "Forward-simulation PBPK/PD: no parameter was estimated from",
      "individual data in this paper. Drug-specific inputs come from",
      "literature and in-house MDCKII-Pgp transwell experiments (Table 1);",
      "clinical data from 11 published studies plus 19 neonatal CSF",
      "samples (Table 2) were used only for verification. The paper ran",
      "500 virtual individuals per scenario, matched to each study's age",
      "range, dose and sex ratio. The deposited S1 File is the paediatric",
      "(2.6-16.4 y) version of the code; the adult physiology is that of",
      "Verscheijden 2019 S1 Table (adult column)."
    )
  )

  ini({
    # ---- Morphine clearance from venous blood (Table 1). The paediatric
    # form is the Wang 2013 bodyweight-dependent-exponent model: 1.62 L/min
    # for a 70 kg individual, times 60 min/h.
    lcl_adult <- fixed(log(1.962)); label("Morphine blood clearance per kg body weight, adults (L/h/kg)") # Table 1, CLiv (adult) morphine = 1.962*BW
    lcl_child <- fixed(log(1.62)); label("Morphine clearance of a 70 kg individual, children and neonates (L/min)") # Table 1, CLiv (2-18y) and (neonates) = 60*1.62*(weight/70)^(...); S1 File line 706
    kmax_bde <- fixed(1.47); label("Baseline allometric exponent of the bodyweight-dependent-exponent clearance model (unitless)") # Table 1, CLiv (2-18y): exponent 1.47 - ...
    kdec <- fixed(0.59); label("Maximum decrease of the allometric exponent (unitless)") # Table 1, CLiv (2-18y): 0.59*weight^4.62/(...)
    khal <- fixed(4.01); label("Body weight at half-maximal exponent decrease (kg)") # Table 1, CLiv (2-18y): 4.01^4.62
    hill_bde <- fixed(4.62); label("Hill coefficient of the exponent decrease (unitless)") # Table 1, CLiv (2-18y): exponent 4.62

    # ---- M6G clearance and formation (Table 1).
    lcl_m6g_adult <- fixed(log(0.131)); label("M6G clearance per kg body weight, adults (L/h/kg)") # Table 1, CLiv (adult) M6G = 0.131*BW
    lcl_m6g_child <- fixed(log(0.114)); label("M6G clearance per kg body weight, children (L/h/kg)") # Table 1, CLiv (2-18y) M6G = 0.114*BW
    lcl_m6g_neo <- fixed(log(0.017)); label("M6G clearance per kg body weight, neonates (L/h/kg)") # Table 1, CLiv (neonates) M6G = 0.017*BW
    fm_m6g <- fixed(0.1); label("Fraction of morphine clearance forming M6G, from 3 months of age (fraction)") # Table 1, Fractional M6G formation 0.1 (2 years-adult)
    fmneo_m6g <- fixed(0.044); label("Fraction of morphine clearance forming M6G, neonates (fraction)") # Table 1, Fractional M6G formation 0.044 (neonates)
    mw <- fixed(285.343); label("Molecular weight of morphine (g/mol)") # Table 1, MW morphine = 285.343
    mw_m6g <- fixed(461.467); label("Molecular weight of M6G (g/mol)") # Table 1, MW M6G = 461.467

    # ---- Physicochemical inputs to the Rodgers and Rowland partition
    # coefficients (Table 1; olive oil:water values from the S1 File).
    logp <- fixed(0.89); label("Morphine log10 octanol:water partition coefficient (unitless)") # Table 1, LogP morphine = 0.89
    logpvo <- fixed(-0.35765); label("Morphine log10 olive oil:water partition coefficient (unitless)") # S1 File line 572: Povo = 10^(-0.35765), calculated using Simcyp
    pka <- fixed(8.21); label("Morphine basic pKa (unitless)") # Table 1, pKa morphine = 8.21 (base)
    ep <- fixed(1.34); label("Morphine erythrocyte:plasma concentration ratio (unitless)") # Table 1, EP morphine = 1.34
    logp_m6g <- fixed(-2.9); label("M6G log10 octanol:water partition coefficient (unitless)") # Table 1, LogP M6G = -2.9
    logpvo_m6g <- fixed(-4.5835); label("M6G log10 olive oil:water partition coefficient (unitless)") # S1 File line 591: PovoM = 10^(-4.5835), calculated using Simcyp
    pka_m6g <- fixed(2.87); label("M6G acidic pKa (unitless)") # Table 1, pKa M6G = 2.87 (acid)
    pkb_m6g <- fixed(9.12); label("M6G basic pKa (unitless)") # Table 1, pKa M6G = 9.12 (base)
    ep_m6g <- fixed(0.15); label("M6G erythrocyte:plasma concentration ratio (unitless)") # Table 1, EP M6G = 0.15
    kpscalar_m6g <- fixed(0.5); label("Scalar applied to all M6G tissue partition coefficients (unitless)") # Methods: 'M6G Kp values were optimized using a kp-scalar of 0.5'; S1 File line 599

    # ---- Binding (Table 1; blood unbound fractions from the S1 File).
    fu <- fixed(0.64); label("Morphine fraction unbound in plasma (fraction)") # Table 1, Fupl morphine = 0.64
    bp_fubb <- fixed(1.14); label("Morphine blood:plasma ratio used to convert fu plasma to fu blood (unitless)") # S1 File line 475: Fubb = 0.64/1.14 (Fupl/BP)
    fu_br <- fixed(0.5); label("Morphine fraction unbound in brain mass (fraction)") # Table 1, Fubm morphine = 0.5
    fu_csf <- fixed(1); label("Morphine fraction unbound in CSF (fraction)") # Table 1, FuCSF morphine = 1
    fu_m6g <- fixed(0.83); label("M6G fraction unbound in plasma, also used for brain blood (fraction)") # Table 1, Fupl M6G = 0.83; S1 File line 480: FubbM = 0.83 (=Fupl)
    fu_br_m6g <- fixed(0.99); label("M6G fraction unbound in brain mass (fraction)") # Table 1, Fubm M6G = 0.99
    fu_csf_m6g <- fixed(1); label("M6G fraction unbound in CSF (fraction)") # Table 1, FuCSF M6G = 1

    # ---- Brain barrier permeability-surface-area products (Table 1),
    # expressed per kg of brain (Vbrain * 1.04 kg/L).
    lpsb <- fixed(log(0.2112)); label("Morphine blood-brain barrier PS product per kg brain (L/h/kg)") # Table 1, PSb morphine = 0.2112*(Vbrain*1.04)
    lpsb_m6g <- fixed(log(0.0072)); label("M6G blood-brain barrier PS product per kg brain (L/h/kg)") # Table 1, PSb M6G = 0.0072*(Vbrain*1.04)
    f_psc <- fixed(0.5); label("Blood-CSF barrier PS product as a fraction of the blood-brain barrier PS product (fraction)") # Table 1, PSc = PSb*0.5, 'Surface area ~50% of BBB'
    lpse <- fixed(log(300)); label("Brain mass to cranial CSF PS product, both analytes (L/h)") # Table 1, PSe = 300, 'Assumed to be no barrier'

    # ---- Active P-gp efflux of morphine at the blood-brain barrier
    # (Methods Eqs 2-3). CLefflux,vivo = 2*(ER-1)*Papp*SA/Procell *
    # abundance(ex vivo)/abundance(in vitro) * BMvPGB * BW = 0.14 L/h.
    er_pgp <- fixed(1.30); label("Net in vitro efflux ratio in MDCKII-Pgp transwells (unitless)") # Table 1, Net in vitro ER = 1.30
    papp_ab <- fixed(2.12e-6); label("In vitro apical-to-basolateral permeability with Pgp inhibited (cm/s)") # Table 1, Papp,AB (inhibited) = 2.12e-6
    sa_filter <- fixed(0.33); label("Transwell filter surface area (cm^2)") # Methods Eq 2: SA = 0.33 cm2
    prot_cell <- fixed(81.4); label("Cellular protein per transwell filter (ug)") # Methods Eq 2: Procell = 81,4 ug
    abund_vivo <- fixed(4.21); label("Pgp abundance in adult isolated brain microvessels (pmol/mg total protein)") # Methods Eq 3: abundance (ex vivo) = 4.21
    abund_vitro <- fixed(0.19); label("Pgp abundance in MDCKII-Pgp cells (pmol/mg total protein)") # Methods Eq 3 and Results: maintained at 0.19 pmol/mg
    bmvpgb <- fixed(0.244); label("Brain microvessel protein per gram brain (mg/g)") # Methods Eq 3: BMvPGB = 0.244 mg/g
    brain_wt_pgp <- fixed(1400); label("Brain weight used in the Pgp in vitro-in vivo extrapolation (g)") # Methods Eq 3: BW = 1400 g
    f_pgp_term <- fixed(0.41); label("BBB Pgp expression at term birth relative to adults (fraction)") # Methods: 'at term 41% of adult expression'
    age_pgp_mature <- fixed(0.5); label("Postnatal age at which BBB Pgp expression is fully mature (years)") # Methods: 'fully matured at 6 months of age'

    # ---- CSF production rate (L/h).
    lq_csf_prod_adult <- fixed(log(0.021)); label("CSF production rate, adults (L/h)") # Verscheijden 2019 S1 Table, Qproductionrate adult = 0.021
    lq_csf_prod_child <- fixed(log(0.024)); label("CSF production rate, 3 months to 18 years (L/h)") # Verscheijden 2019 S1 Table, Qproductionrate 3m-18y = 0.024; S1 File line 399
    lq_csf_prod_neo <- fixed(log(0.010)); label("CSF production rate, neonates (L/h)") # Methods: 'CSF production rate in neonates was assumed to be 10 mL/h'
    lf_bulk <- fixed(log(0.25)); label("Bulk flow brain mass to cranial CSF as a fraction of CSF production (fraction)") # Verscheijden 2019 S1 Table, Qbulk = 0.25*Qproductionrate; S1 File line 403
    lf_ssink <- fixed(log(0.38)); label("Spinal CSF to blood flow as a fraction of total CSF outflow (fraction)") # Verscheijden 2019 S1 Table, Qssink = 0.38*(0.75*Qproductionrate+Qbulk); S1 File line 407
    lf_sout <- fixed(log(0.9)); label("Spinal to cranial CSF flow as a fraction of the spinal sink flow (fraction)") # Verscheijden 2019 S1 Table, Qsout = 0.9*Qssink; S1 File line 411

    # ---- Adult haematocrit (children: age equations in model()).
    lhct_adult_male <- fixed(log(0.43)); label("Haematocrit, adult males (fraction)") # Verscheijden 2019 S1 Table, Hematocrit adult male = 0.43
    lhct_adult_female <- fixed(log(0.38)); label("Haematocrit, adult females (fraction)") # Verscheijden 2019 S1 Table, Hematocrit adult female = 0.38

    # ---- PD (Methods Eqs 4-6).
    km_nm <- fixed(11.8); label("Morphine mu-opioid receptor equilibrium dissociation constant (nM)") # Methods: 'average values of 11.8 nM for morphine'
    km_m6g_nm <- fixed(42.5); label("M6G mu-opioid receptor equilibrium dissociation constant (nM)") # Methods: '42.5 nM for M6G'
    br50 <- fixed(59.26); label("Percentage of receptors bound by morphine giving 50% relative response (percent)") # Methods: BR50 = 59.3 (95% CI 54.1-64.4); S1 File line 965 uses 59.26
    hill <- fixed(4.217); label("Hill slope of the morphine binding-response relationship (unitless)") # Methods: Hill slope 4.2 (95% CI 2.6-5.8); S1 File line 965 uses 4.217
    br50_m6g <- fixed(17); label("Percentage of receptors bound by M6G giving 50% relative response (percent)") # Methods: M6G 'BR50: 17'
    hill_m6g <- fixed(2.2); label("Hill slope of the M6G binding-response relationship (unitless)") # Methods: M6G 'Hill slope: 2.2'

    # ---- Between-subject variability of the deposited virtual-population
    # code (S1 File): parameter * exp(rnorm(0, SD)), so variance = SD^2.
    etalcl ~ fixed(0.16) # S1 File line 707: CLivCV = rnorm(N, 0, 0.4)
    etalcl_m6g ~ fixed(0.0324) # S1 File line 711: CLMivCV = rnorm(N, 0, 0.18)
    etalq_csf_prod ~ fixed(0.01) # S1 File line 400: QproductionrateCV = rnorm(N, 0, 0.1); Verscheijden 2019 S1 Table 10%
    etalf_bulk ~ fixed(0.0064) # S1 File line 404: QbulkCV = rnorm(N, 0, 0.08); Verscheijden 2019 S1 Table 8%
    etalf_ssink ~ fixed(0.09) # S1 File line 408: QssinkCV = rnorm(N, 0, 0.3); Verscheijden 2019 S1 Table 30%
    etalf_sout ~ fixed(1) # S1 File line 412: QsoutCV = rnorm(N, 0, 1); Verscheijden 2019 S1 Table 100%
    etalhct ~ fixed(0.004225) # S1 File line 487: HematocritmaleCV = rnorm(N, 0, 0.065); females use SD 0.071 (line 492) via the SEXF rescaling in model()

    # ---- Residual error: not reported (forward simulation only).
    propSd <- fixed(0); label("Proportional residual error on plasma morphine (fraction; not reported)") # not reported: forward-simulation model
  })

  model({
    # =====================================================================
    # Age / sex switches. adult: >= 18 y (Verscheijden 2019 S1 Table adult
    # column). neo: < 0.25 y (Table 1 neonatal clearance / fm values and the
    # 10 mL/h neonatal CSF production; the boundary is assumed).
    # =====================================================================
    adult <- 0
    if (AGE >= 18) adult <- 1
    neo <- 0
    if (AGE < 0.25) neo <- 1
    fem <- 0
    if (SEXF > 0.5) fem <- 1

    bsa <- 0.007184 * HT^0.725 * WT^0.425 # S1 File line 27 (Du Bois)
    htm <- HT / 100

    # =====================================================================
    # Organ volumes (L). Paediatric: S1 File lines 266-390. Adult:
    # Verscheijden 2019 S1 Table adult column.
    # =====================================================================
    # Brain
    v_brain_p <- (10 * (AGE + 0.315) / (9 + 6.92 * AGE)) / 1.04
    v_brain_a <- (1.449 - 3.62 / WT) / 1.04
    v_brain <- v_brain_p * (1 - adult) + v_brain_a * adult
    v_bb <- 0.05 * v_brain
    v_endo <- 0.005 * v_brain
    csf_frac <- 0.105 * (1 - fem) + 0.092 * fem
    v_scsf_p <- (1.94 * WT + 0.13) / 1000
    if (v_scsf_p > (0.143 / 0.8) * 0.2) v_scsf_p <- (0.143 / 0.8) * 0.2
    v_ccsf <- 0.143 * (1 - adult) + v_brain * csf_frac * 0.8 * adult
    v_scsf <- v_scsf_p * (1 - adult) + v_brain * csf_frac * 0.2 * adult
    v_bm <- v_brain - v_endo - v_bb - v_ccsf - v_scsf

    # Lung (same equations for children and adults)
    v_lung_m <- ((29.08 * htm * WT^0.5 + 11.06 + 35.47 * htm * WT^0.5 + 5.53) / 1000) / 1.05
    v_lung_f <- ((31.46 * htm * WT^0.5 + 1.43 + 35.3 * htm * WT^0.5 + 1.53) / 1000) / 1.05
    v_lung <- v_lung_m * (1 - fem) + v_lung_f * fem

    # Adipose
    v_ad01 <- ((908.4 + 0.706 * (WT * 1000) - 53 * HT + 358.5 * fem - 3.057 * (AGE * 365)) / 1000) / 0.92
    if (v_ad01 < 0) v_ad01 <- 0.01
    v_ad13 <- ((908.4 + 0.706 * (WT * 1000) - 53 * HT + 358.5 * fem - 3.057 * (1.1 * 365)) / 1000) / 0.92
    if (v_ad13 < 0) v_ad13 <- 0.01
    v_ad312 <- ((0.534 * WT - 1.59 * AGE + 3.03) / 0.92) * (1 - fem) +
      ((0.642 * WT - 0.12 * HT - 0.606 * AGE + 8.98) / 0.92) * fem
    if (v_ad312 < 0) v_ad312 <- 0.01
    v_ad12 <- ((1.36 * WT) / htm - 42) * (1 - fem) + ((1.61 * WT) / htm - 38.3) * fem
    if (v_ad12 < 0) v_ad12 <- 0.01
    v_adipose <- v_ad12
    if (AGE < 12.3) v_adipose <- v_ad312
    if (AGE < 3) v_adipose <- v_ad13
    if (AGE < 1) v_adipose <- v_ad01

    # Bone
    v_bo01 <- (((77.24 + 24.94 * WT + 0.21 * (AGE * 365) - 1.889 * HT) / 0.33) / 1000) / 1.3
    v_bo13 <- (((77.24 + 24.94 * WT + 0.21 * (1.1 * 365) - 1.889 * HT) / 0.33) / 1000) / 1.3
    if (v_bo13 < 0) v_bo13 <- 0.01
    v_bo318 <- (((7.4 * HT + 28.9 * WT - 789.6) / 0.4) / 1000) / 1.3
    v_bo_a <- ((WT - v_adipose * 0.92) / 1.3 * 0.058) * (1 - fem) + ((WT - v_adipose * 0.92) / 1.3 * 0.051) * fem
    v_bone <- v_bo_a
    if (AGE < 18) v_bone <- v_bo318
    if (AGE < 3) v_bone <- v_bo13
    if (AGE < 1) v_bone <- v_bo01

    # Heart
    v_heart_p <- (((22.81 * htm * WT^0.5 - 4.15) / 1000) / 1.05) * (1 - fem) +
      (((19.99 * htm * WT^0.5 - 1.53) / 1000) / 1.05) * fem
    v_heart_a <- ((155.18 * bsa^1.29) / 1000 / 1.05) * (1 - fem) + ((124.13 * bsa^1.242) / 1000 / 1.05) * fem
    v_heart <- v_heart_p * (1 - adult) + v_heart_a * adult

    # Kidney
    v_kidney_p <- (4.214 * WT^0.823 + 4.456 * WT^0.795) / 1000
    v_kidney_a <- ((15.4 + 2.04 * WT + 51.8 * htm^2) / 1000) / 1.05
    v_kidney <- v_kidney_p * (1 - adult) + v_kidney_a * adult

    # Muscle
    v_muscle_p <- ((0.3 + ((0.54 - 0.3) / 18) * AGE) * (1 - fem) + (0.3 + ((0.489 - 0.3) / 18) * AGE) * fem) *
      (WT - v_adipose * 0.92) / 1.04
    v_muscle_a <- ((0.244 * WT + 7.8 * htm - 0.098 * AGE + 3.3) / 1.04) * (1 - fem) +
      ((0.244 * WT + 7.8 * htm - 0.098 * AGE - 3.3) / 1.04) * fem
    v_muscle <- v_muscle_p * (1 - adult) + v_muscle_a * adult

    # Skin (same equation for children and adults)
    v_skin <- (bsa / 1000) * 45.655 + (bsa / 1000) * 1240

    # Spleen (the paediatric male equation multiplies by 1.06, as printed)
    v_spleen_p <- (((8.74 * htm * WT^0.5 + 11.06) / 1000) * 1.06) * (1 - fem) +
      (((9.36 * htm * WT^0.5 + 7.98) / 1000) / 1.06) * fem
    v_spleen_a <- (6.516 * WT^0.797) / 1000
    v_spleen <- v_spleen_p * (1 - adult) + v_spleen_a * adult

    # Gut
    gut_frac <- 0.021 * (1 - adult) + (0.021 * (1 - fem) + 0.027 * fem) * adult
    v_gut <- (gut_frac * (WT - v_adipose * 0.92)) / 1.05

    # Liver
    v_liver_p <- (((576.9 * htm + 8.9 * WT - 159.7) / 1000) / 1.05) * (1 - fem) +
      (((674.3 * htm + 6.5 * WT - 214.4) / 1000) / 1.05) * fem
    v_liver_a <- (1072.8 * bsa - 345.7) / 1000
    v_liver <- v_liver_p * (1 - adult) + v_liver_a * adult

    # Arterial and venous blood (each half of total blood volume)
    v_bl01 <- ((10^(0.7891 * (log10(WT) + 0.004132 * HT + 1.8117))) / 1000) / 2
    v_bl_c1 <- ((10^(0.6459 * log10(WT) + 0.002743 * HT + 2.0324)) / 1000) / 2
    v_bl_c2 <- ((10^(0.6412 * log10(WT) + 0.00127 * HT + 2.2169)) / 1000) / 2
    v_bl_a <- ((((13.1 * HT + 18.05 * WT - 480) / 0.5723) / 1000) / 2) * (1 - fem) +
      ((((35.5 * HT + 2.27 * WT - 3382) / 0.6178) / 1000) / 2) * fem
    v_blood_half <- v_bl_a
    if (AGE < 18) v_blood_half <- v_bl_c1
    if (AGE < 18 && AGE >= 6 && fem == 1) v_blood_half <- v_bl_c2
    if (AGE < 1) v_blood_half <- v_bl01
    v_venous <- v_blood_half
    v_arterial <- v_blood_half

    # Rest of body (body density 1 kg/L)
    v_other <- WT - v_venous - v_arterial - v_liver - v_gut - v_spleen - v_skin - v_muscle -
      v_kidney - v_heart - v_bone - v_adipose - v_lung - v_brain
    if (v_other < 0) v_other <- 0.01

    # =====================================================================
    # Blood and CSF flows (L/h). Paediatric: S1 File lines 393-470. Adult:
    # Verscheijden 2019 S1 Table adult column.
    # =====================================================================
    q_co_p <- bsa * (110 + (184.974 * (exp(-0.0378 * AGE) - exp(-0.2477 * AGE))))
    q_co_a <- bsa * 60 * (3 - 0.01 * (AGE - 20))
    q_co <- q_co_p * (1 - adult) + q_co_a * adult
    q_brain <- (q_co * ((10 + 2290 * (exp(-0.608 * AGE) - exp(-0.639 * AGE))) / 100)) * (1 - adult) + q_co * 0.12 * adult

    q_adipose_p <- 0.05 * (1 - fem) + ((5 + (3.59 * AGE^5 / (14.49^5 + AGE^5))) / 100) * fem
    q_adipose <- q_co * (q_adipose_p * (1 - adult) + (0.05 * (1 - fem) + 0.085 * fem) * adult)
    q_bone <- q_co * 0.05
    q_heart <- q_co * (0.04 * (1 - fem) + 0.05 * fem)
    q_kidney_p <- ((4.53 + (14.63 * AGE / (0.188 + AGE))) / 100) * (1 - fem) +
      ((4.53 + (13 * AGE^1.15 / (0.188^1.15 + AGE^1.15))) / 100) * fem
    q_kidney <- q_co * (q_kidney_p * (1 - adult) + (0.19 * (1 - fem) + 0.17 * fem) * adult)
    q_muscle_p <- ((6.03 + (12 * AGE^2.5 / (11^2.5 + AGE^2.5))) / 100) * (1 - fem) +
      ((6.03 + (7 * AGE^2.5 / (12^2.5 + AGE^2.5))) / 100) * fem
    q_muscle <- q_co * (q_muscle_p * (1 - adult) + (0.17 * (1 - fem) + 0.12 * fem) * adult)
    q_skin_p <- (1.0335 + (4 * AGE^5 / (5.16^5 + AGE^5))) / 100
    q_skin <- q_co * (q_skin_p * (1 - adult) + 0.05 * adult)
    q_spleen <- q_co * (0.02 * (1 - fem) + 0.03 * fem)
    q_gut <- q_co * (0.15 * (1 - fem) + 0.17 * fem)
    q_ha <- q_co * 0.065
    q_liver <- q_co * (0.235 * (1 - fem) + 0.265 * fem)

    # Rest-of-body flow closes the cardiac output; if it would be negative
    # the lung flow absorbs the deficit and the rest flow is set to 0.1 L/h
    # (S1 File lines 428-431).
    q_rest_raw <- q_co - q_brain - q_adipose - q_bone - q_heart - q_kidney - q_muscle - q_skin -
      q_spleen - q_gut - q_ha
    q_lung <- q_co
    q_rest <- q_rest_raw
    if (q_rest_raw < 0) {
      q_lung <- q_co - q_rest_raw
      q_rest <- 0.1
    }

    # CSF flows (S1 File lines 399-416). The spinal-sink and spinal-outflow
    # typical values are built from the TYPICAL production and bulk flows;
    # the individual values then receive their own variability.
    q_csf_prod_typ <- exp(lq_csf_prod_child) * (1 - adult) + exp(lq_csf_prod_adult) * adult
    if (neo == 1) q_csf_prod_typ <- exp(lq_csf_prod_neo)
    q_csf_prod <- q_csf_prod_typ * exp(etalq_csf_prod)
    q_bulk_typ <- exp(lf_bulk) * q_csf_prod_typ
    q_bulk <- q_bulk_typ * exp(etalf_bulk)
    q_ssink_typ <- exp(lf_ssink) * (0.75 * q_csf_prod_typ + q_bulk_typ)
    q_ssink <- q_ssink_typ * exp(etalf_ssink)
    q_sout <- exp(lf_sout) * q_ssink_typ * exp(etalf_sout)
    q_sin <- q_ssink + q_sout
    q_csink <- 0.75 * q_csf_prod + q_bulk - q_sin + q_sout

    # =====================================================================
    # Haematocrit and blood:plasma ratios (S1 File lines 484-501).
    # =====================================================================
    hct_p <- ((53 - ((43.0 * AGE^1.12 / (0.05^1.12 + AGE^1.12)) * (1 + (-0.93 * AGE^0.25 / (0.10^0.25 + AGE^0.25))))) / 100) * (1 - fem) +
      ((53 - ((37.4 * AGE^1.12 / (0.05^1.12 + AGE^1.12)) * (1 + (-0.80 * AGE^0.25 / (0.10^0.25 + AGE^0.25))))) / 100) * fem
    hct_a <- exp(lhct_adult_male) * (1 - fem) + exp(lhct_adult_female) * fem
    hct <- (hct_p * (1 - adult) + hct_a * adult) * exp(etalhct * (1 + (0.071 / 0.065 - 1) * fem))
    bp <- 1 - hct + ep * hct
    bp_m6g <- 1 - hct + ep_m6g * hct
    fu_bb <- fu / bp_fubb
    fu_bb_m6g <- fu_m6g

    # =====================================================================
    # Tissue composition (fractions of tissue volume): EW extracellular
    # water, IW intracellular water, NL neutral lipid, NP neutral
    # phospholipid. Paediatric: S1 File lines 512-569. Adult: Verscheijden
    # 2019 S1 Table adult column.
    # =====================================================================
    lage <- log10(AGE)
    few_ad <- ((32.154 - 2.7863 * lage) * (14.1 / 18 / 100)) * (1 - adult) + 0.141 * adult
    fiw_ad <- ((32.154 - 2.7863 * lage) * (3.9 / 18 / 100)) * (1 - adult) + 0.039 * adult
    fnl_ad <- ((35.5 + 43.46 * AGE / (1.5 + AGE)) * (79 / 79.2 / 100)) * (1 - adult) + 0.79 * adult
    fnp_ad <- ((35.5 + 43.46 * AGE / (1.5 + AGE)) * (0.2 / 79.2 / 100)) * (1 - adult) + 0.002 * adult

    few_bo <- ((64.179 - 1.2697 * AGE) * (9.8 / 43.9 / 100)) * (1 - adult) + 0.098 * adult
    fiw_bo <- ((64.179 - 1.2697 * AGE) * (34.1 / 43.9 / 100)) * (1 - adult) + 0.341 * adult
    fnl_bo <- ((0.2 + 0.3655 * AGE) * (7.4 / 7.51 / 100)) * (1 - adult) + 0.074 * adult
    fnp_bo <- ((0.2 + 0.3655 * AGE) * (0.11 / 7.51 / 100)) * (1 - adult) + 0.0011 * adult

    few_gu <- ((75.378 - 0.3932 * lage) * (26.7 / 71.8 / 100)) * (1 - adult) + 0.267 * adult
    fiw_gu <- ((75.378 - 0.3932 * lage) * (45.1 / 71.8 / 100)) * (1 - adult) + 0.451 * adult
    fnl_gu <- ((2.5 + 0.185 * AGE) * (4.87 / 6.5 / 100)) * (1 - adult) + 0.0487 * adult
    fnp_gu <- ((2.5 + 0.185 * AGE) * (1.63 / 6.5 / 100)) * (1 - adult) + 0.0163 * adult

    few_he <- ((84.523 - 0.4249 * AGE) * (31.3 / 75.8 / 100)) * (1 - adult) + 0.313 * adult
    fiw_he <- ((84.523 - 0.4249 * AGE) * (44.5 / 75.8 / 100)) * (1 - adult) + 0.445 * adult
    fnl_he <- ((2.3159 + 0.0797 * AGE) * (1.15 / 2.81 / 100)) * (1 - adult) + 0.0115 * adult
    fnp_he <- ((2.3159 + 0.0797 * AGE) * (1.66 / 2.81 / 100)) * (1 - adult) + 0.0166 * adult

    few_ki <- ((83.278 - 0.2162 * AGE) * (28.3 / 78.3 / 100)) * (1 - adult) + 0.283 * adult
    fiw_ki <- ((83.278 - 0.2162 * AGE) * (50 / 78.3 / 100)) * (1 - adult) + 0.50 * adult
    fnl_ki <- ((2.73 + 1.995 * AGE / (2.59 + AGE)) * (2.07 / 3.69 / 100)) * (1 - adult) + 0.0207 * adult
    fnp_ki <- ((2.73 + 1.995 * AGE / (2.59 + AGE)) * (1.62 / 3.69 / 100)) * (1 - adult) + 0.0162 * adult

    few_li <- ((75.69 - 0.573 * lage) * (16.5 / 75.1 / 100)) * (1 - adult) + 0.165 * adult
    fiw_li <- ((75.69 - 0.573 * lage) * (58.6 / 75.1 / 100)) * (1 - adult) + 0.586 * adult
    fnl_li <- ((3 + 3.089 * AGE / (1.8 + AGE)) * (3.48 / 6 / 100)) * (1 - adult) + 0.0348 * adult
    fnp_li <- ((3 + 3.089 * AGE / (1.8 + AGE)) * (2.52 / 6 / 100)) * (1 - adult) + 0.0252 * adult

    few_lu <- ((80.973 - 0.4916 * lage) * (34.8 / 81.1 / 100)) * (1 - adult) + 0.348 * adult
    fiw_lu <- ((80.973 - 0.4916 * lage) * (46.3 / 81.1 / 100)) * (1 - adult) + 0.463 * adult
    fnl_lu <- ((1.857 - 0.211 * lage) * (0.3 / 1.2 / 100)) * (1 - adult) + 0.003 * adult
    fnp_lu <- ((1.857 - 0.211 * lage) * (0.9 / 1.2 / 100)) * (1 - adult) + 0.009 * adult

    few_mu <- ((77.211 - 0.4321 * lage) * (9.1 / 76 / 100)) * (1 - adult) + 0.091 * adult
    fiw_mu <- ((77.211 - 0.4321 * lage) * (66.9 / 76 / 100)) * (1 - adult) + 0.667 * adult
    fnl_mu <- ((1.9852 + 0.0649 * AGE) * (2.38 / 3.1 / 100)) * (1 - adult) + 0.0238 * adult
    fnp_mu <- ((1.9852 + 0.0649 * AGE) * (0.72 / 3.1 / 100)) * (1 - adult) + 0.0072 * adult

    few_sk <- ((72.395 - 1.1462 * lage) * (62.3 / 71.77 / 100)) * (1 - adult) + 0.623 * adult
    fiw_sk <- ((72.395 - 1.1462 * lage) * (9.47 / 71.77 / 100)) * (1 - adult) + 0.0947 * adult
    fnl_sk <- (3.95 * (2.84 / 3.95 / 100)) * (1 - adult) + 0.0284 * adult
    fnp_sk <- (3.95 * (1.11 / 3.95 / 100)) * (1 - adult) + 0.0111 * adult

    few_sp <- ((79.952 - 0.4178 * lage) * (20.8 / 78.7 / 100)) * (1 - adult) + 0.208 * adult
    fiw_sp <- ((79.952 - 0.4178 * lage) * (57.9 / 78.7 / 100)) * (1 - adult) + 0.579 * adult
    fnl_sp <- ((1.5 + 0.015 * AGE) * (2.01 / 3.99 / 100)) * (1 - adult) + 0.0201 * adult
    fnp_sp <- ((1.5 + 0.015 * AGE) * (1.98 / 3.99 / 100)) * (1 - adult) + 0.0198 * adult

    fiw_rbc <- 66 / 100
    fnl_rbc <- 0.3 * (0.17 / 0.46 / 100)
    fnp_rbc <- 0.3 * (0.29 / 0.46 / 100)

    # =====================================================================
    # Rodgers and Rowland partition coefficients, morphine (monoprotic
    # base; S1 File lines 571-588). Tissue acidic-phospholipid
    # concentrations (mg/g) are the constants multiplying KaAP.
    # =====================================================================
    pow <- 10^logp
    povo <- 10^logpvo
    xb7 <- 1 + 10^(pka - 7)
    xb74 <- 1 + 10^(pka - 7.4)
    kubc <- (hct - 1 + bp) / (fu_bb * hct)
    kaap <- (kubc - ((1 + 10^(pka - 7.22)) / xb74) * fiw_rbc - ((pow * fnl_rbc + (0.3 * pow + 0.7) * fnp_rbc) / xb74)) *
      (xb74 / (0.44 * 10^(pka - 7.22)))
    kp_adipose <- (few_ad + (xb7 * fiw_ad) / xb74 + (povo * fnl_ad + (0.3 * povo + 0.7) * fnp_ad) / xb74 + (kaap * 0.4 * 10^(pka - 7)) / xb74) * fu
    kp_lung <- (few_lu + (xb7 * fiw_lu) / xb74 + (pow * fnl_lu + (0.3 * pow + 0.7) * fnp_lu) / xb74 + (kaap * 0.5 * 10^(pka - 7)) / xb74) * fu
    kp_bone <- (few_bo + (xb7 * fiw_bo) / xb74 + (pow * fnl_bo + (0.3 * pow + 0.7) * fnp_bo) / xb74 + (kaap * 0.67 * 10^(pka - 7)) / xb74) * fu
    kp_heart <- (few_he + (xb7 * fiw_he) / xb74 + (pow * fnl_he + (0.3 * pow + 0.7) * fnp_he) / xb74 + (kaap * 3.07 * 10^(pka - 7)) / xb74) * fu
    kp_kidney <- (few_ki + (xb7 * fiw_ki) / xb74 + (pow * fnl_ki + (0.3 * pow + 0.7) * fnp_ki) / xb74 + (kaap * 2.48 * 10^(pka - 7)) / xb74) * fu
    kp_muscle <- (few_mu + (xb7 * fiw_mu) / xb74 + (pow * fnl_mu + (0.3 * pow + 0.7) * fnp_mu) / xb74 + (kaap * 2.49 * 10^(pka - 7)) / xb74) * fu
    kp_skin <- (few_sk + (xb7 * fiw_sk) / xb74 + (pow * fnl_sk + (0.3 * pow + 0.7) * fnp_sk) / xb74 + (kaap * 1.32 * 10^(pka - 7)) / xb74) * fu
    kp_spleen <- (few_sp + (xb7 * fiw_sp) / xb74 + (pow * fnl_sp + (0.3 * pow + 0.7) * fnp_sp) / xb74 + (kaap * 2.81 * 10^(pka - 7)) / xb74) * fu
    kp_gut <- (few_gu + (xb7 * fiw_gu) / xb74 + (pow * fnl_gu + (0.3 * pow + 0.7) * fnp_gu) / xb74 + (kaap * 2.84 * 10^(pka - 7)) / xb74) * fu
    kp_liver <- (few_li + (xb7 * fiw_li) / xb74 + (pow * fnl_li + (0.3 * pow + 0.7) * fnp_li) / xb74 + (kaap * 5.09 * 10^(pka - 7)) / xb74) * fu
    # Rest of body: Kp = 1 (S1 File line 588), written inline below.

    # =====================================================================
    # Rodgers and Rowland partition coefficients, M6G (zwitterion; S1 File
    # lines 590-610). Lung intracellular pH 6.7 and kidney 7.2 as coded;
    # the KaAP estimate is floored at zero; all Kp are scaled by 0.5.
    # =====================================================================
    powm <- 10^logp_m6g
    povom <- 10^logpvo_m6g
    zi <- 10^(pkb_m6g - pka_m6g)
    z74 <- 1 + 10^(pkb_m6g - 7.4) + 10^(7.4 - pka_m6g) + zi
    z7 <- 1 + 10^(pkb_m6g - 7) + 10^(7 - pka_m6g) + zi
    z67 <- 1 + 10^(pkb_m6g - 6.7) + 10^(6.7 - pka_m6g) + zi
    z72 <- 1 + 10^(pkb_m6g - 7.2) + 10^(7.2 - pka_m6g) + zi
    a7 <- 10^(pkb_m6g - 7) + 10^(7 - pka_m6g)
    a67 <- 10^(pkb_m6g - 6.7) + 10^(6.7 - pka_m6g)
    a72 <- 10^(pkb_m6g - 7.2) + 10^(7.2 - pka_m6g)
    y74 <- 1 + 10^(pkb_m6g - 7.4) + 10^(7.4 - pka_m6g)
    y722 <- 10^(pkb_m6g - 7.22) + 10^(7.22 - pka_m6g)
    kubc_m6g <- (hct - 1 + bp_m6g) / (fu_bb_m6g * hct)
    kaap_m6g <- (kubc_m6g - ((1 + y722) / y74) * fiw_rbc - (powm * fnl_rbc + (0.3 * powm + 0.7) * fnp_rbc) / y74) *
      (y74 / (0.44 * y722))
    if (kaap_m6g < 0) kaap_m6g <- 0
    kp_adipose_m6g <- (few_ad + (z7 * fiw_ad) / z74 + (povom * fnl_ad + (0.3 * povom + 0.7) * fnp_ad) / z74 + (kaap_m6g * 0.4 * a7) / z74) * fu_m6g * kpscalar_m6g
    kp_lung_m6g <- (few_lu + (z67 * fiw_lu) / z74 + (powm * fnl_lu + (0.3 * powm + 0.7) * fnp_lu) / z74 + (kaap_m6g * 0.5 * a67) / z74) * fu_m6g * kpscalar_m6g
    kp_bone_m6g <- (few_bo + (z7 * fiw_bo) / z74 + (powm * fnl_bo + (0.3 * powm + 0.7) * fnp_bo) / z74 + (kaap_m6g * 0.67 * a7) / z74) * fu_m6g * kpscalar_m6g
    kp_heart_m6g <- (few_he + (z7 * fiw_he) / z74 + (powm * fnl_he + (0.3 * powm + 0.7) * fnp_he) / z74 + (kaap_m6g * 3.07 * a7) / z74) * fu_m6g * kpscalar_m6g
    kp_kidney_m6g <- (few_ki + (z72 * fiw_ki) / z74 + (powm * fnl_ki + (0.3 * powm + 0.7) * fnp_ki) / z74 + (kaap_m6g * 2.48 * a72) / z74) * fu_m6g * kpscalar_m6g
    kp_muscle_m6g <- (few_mu + (z7 * fiw_mu) / z74 + (powm * fnl_mu + (0.3 * powm + 0.7) * fnp_mu) / z74 + (kaap_m6g * 2.49 * a7) / z74) * fu_m6g * kpscalar_m6g
    kp_skin_m6g <- (few_sk + (z7 * fiw_sk) / z74 + (powm * fnl_sk + (0.3 * powm + 0.7) * fnp_sk) / z74 + (kaap_m6g * 1.32 * a7) / z74) * fu_m6g * kpscalar_m6g
    kp_spleen_m6g <- (few_sp + (z7 * fiw_sp) / z74 + (powm * fnl_sp + (0.3 * powm + 0.7) * fnp_sp) / z74 + (kaap_m6g * 2.81 * a7) / z74) * fu_m6g * kpscalar_m6g
    kp_gut_m6g <- (few_gu + (z7 * fiw_gu) / z74 + (powm * fnl_gu + (0.3 * powm + 0.7) * fnp_gu) / z74 + (kaap_m6g * 2.84 * a7) / z74) * fu_m6g * kpscalar_m6g
    kp_liver_m6g <- (few_li + (z7 * fiw_li) / z74 + (powm * fnl_li + (0.3 * powm + 0.7) * fnp_li) / z74 + (kaap_m6g * 5.09 * a7) / z74) * fu_m6g * kpscalar_m6g
    kp_other_m6g <- kpscalar_m6g # S1 File line 610: KpMrest = 1 * KpscalarM

    # =====================================================================
    # Clearance, metabolite formation, brain transfer.
    # =====================================================================
    e_cl <- kmax_bde - kdec * WT^hill_bde / (khal^hill_bde + WT^hill_bde)
    cl_p <- 60 * exp(lcl_child) * (WT / 70)^e_cl
    cl_a <- exp(lcl_adult) * WT
    cl <- (cl_p * (1 - adult) + cl_a * adult) * exp(etalcl)

    cl_m6g_typ <- exp(lcl_m6g_child) * WT * (1 - adult) + exp(lcl_m6g_adult) * WT * adult
    if (neo == 1) cl_m6g_typ <- exp(lcl_m6g_neo) * WT
    cl_m6g <- cl_m6g_typ * exp(etalcl_m6g)

    frac_m6g <- fm_m6g
    if (neo == 1) frac_m6g <- fmneo_m6g
    # mass of M6G formed per mass of morphine cleared (S1 File line 503:
    # FMformed = 0.1 * 1.62, the M6G:morphine molecular-weight ratio)
    f_form <- frac_m6g * mw_m6g / mw

    brain_kg <- v_brain * 1.04
    psb <- exp(lpsb) * brain_kg
    psc <- f_psc * psb
    psb_m6g <- exp(lpsb_m6g) * brain_kg
    psc_m6g <- f_psc * psb_m6g
    pse <- exp(lpse)

    # P-gp efflux clearance (L/h), Methods Eqs 2-3: uL-scale unit chain
    # cm3/s/ug -> x1000 ug/mg -> x(ex vivo/in vitro abundance) x BMvPGB x
    # brain weight -> mL/s -> x3.6 -> L/h.
    cl_pgp_vitro <- 2 * (er_pgp - 1) * papp_ab * sa_filter / prot_cell
    cl_pgp_adult <- cl_pgp_vitro * 1000 * (abund_vivo / abund_vitro) * bmvpgb * brain_wt_pgp * 3.6
    # Linear postnatal maturation from 41% at birth to 100% at 6 months
    # (the paper gives the anchors, not the shape).
    f_pgp <- min(f_pgp_term + (1 - f_pgp_term) * AGE / age_pgp_mature, 1)
    q_pgp <- cl_pgp_adult * f_pgp

    # =====================================================================
    # Concentrations (mg/L). Blood compartments and brain blood hold whole-
    # blood concentrations; tissues hold total tissue concentrations.
    # =====================================================================
    c_ven <- venous / v_venous
    c_art <- arterial / v_arterial
    c_lung <- lung / v_lung
    c_ad <- adipose / v_adipose
    c_bo <- bone / v_bone
    c_he <- heart / v_heart
    c_ki <- kidney / v_kidney
    c_mu <- muscle / v_muscle
    c_sk <- skin / v_skin
    c_sp <- spleen / v_spleen
    c_ot <- other / v_other
    c_gu <- gut / v_gut
    c_li <- liver / v_liver
    c_bb <- brain_vascular / v_bb
    c_bm <- brain / v_bm
    c_ccsf <- brain_csf_sas_cranial / v_ccsf
    c_scsf <- brain_csf_sas_spinal / v_scsf

    m_ven <- venous_m6g / v_venous
    m_art <- arterial_m6g / v_arterial
    m_lung <- lung_m6g / v_lung
    m_ad <- adipose_m6g / v_adipose
    m_bo <- bone_m6g / v_bone
    m_he <- heart_m6g / v_heart
    m_ki <- kidney_m6g / v_kidney
    m_mu <- muscle_m6g / v_muscle
    m_sk <- skin_m6g / v_skin
    m_sp <- spleen_m6g / v_spleen
    m_ot <- other_m6g / v_other
    m_gu <- gut_m6g / v_gut
    m_li <- liver_m6g / v_liver
    m_bb <- brain_vascular_m6g / v_bb
    m_bm <- brain_m6g / v_bm
    m_ccsf <- brain_csf_sas_cranial_m6g / v_ccsf
    m_scsf <- brain_csf_sas_spinal_m6g / v_scsf

    # =====================================================================
    # Morphine ODEs (S1 File lines 91-109, multiplied through by the
    # compartment volume so the states are amounts in mg). IV doses enter
    # venous blood. Clearance acts on the venous concentration and removes
    # mass from venous blood, as in the Verscheijden 2019 framework code
    # (2019 S1 File line 81). The 2021 S1 File instead subtracts CL x
    # C_venous from the ARTERIAL compartment (S1 File line 97); that form
    # falls below the paper's own Fig 2D adult ECF profile by 41-68% and
    # its Fig 2A paediatric ECF:plasma ratio by about 20%, and turns
    # arterial blood negative when CL nears cardiac output. The venous form
    # reproduces both figures and is used here (see the vignette).
    # =====================================================================
    d/dt(venous) <- q_adipose * c_ad / kp_adipose * bp + q_bone * c_bo / kp_bone * bp +
      q_heart * c_he / kp_heart * bp + q_kidney * c_ki / kp_kidney * bp +
      q_muscle * c_mu / kp_muscle * bp + q_skin * c_sk / kp_skin * bp +
      q_liver * c_li / kp_liver * bp + q_brain * c_bb + q_rest * c_ot * bp -
      q_lung * c_ven - cl * c_ven
    d/dt(arterial) <- q_lung * c_lung / kp_lung * bp -
      (q_rest + q_brain + q_adipose + q_bone + q_heart + q_kidney + q_muscle + q_skin + q_spleen + q_gut + q_ha) * c_art
    d/dt(lung) <- q_lung * (c_ven - c_lung / kp_lung * bp)
    d/dt(adipose) <- q_adipose * (c_art - c_ad / kp_adipose * bp)
    d/dt(bone) <- q_bone * (c_art - c_bo / kp_bone * bp)
    d/dt(heart) <- q_heart * (c_art - c_he / kp_heart * bp)
    d/dt(kidney) <- q_kidney * (c_art - c_ki / kp_kidney * bp)
    d/dt(muscle) <- q_muscle * (c_art - c_mu / kp_muscle * bp)
    d/dt(skin) <- q_skin * (c_art - c_sk / kp_skin * bp)
    d/dt(spleen) <- q_spleen * (c_art - c_sp / kp_spleen * bp)
    d/dt(other) <- q_rest * (c_art - c_ot * bp)
    d/dt(gut) <- q_gut * (c_art - c_gu / kp_gut * bp)
    d/dt(liver) <- q_ha * c_art + q_gut * c_gu / kp_gut * bp + q_spleen * c_sp / kp_spleen * bp -
      q_liver * c_li / kp_liver * bp
    # Brain model (Verscheijden 2019 Eqs 2-5) with Pgp efflux from brain
    # mass to brain blood (S1 File Qtransport term).
    d/dt(brain_vascular) <- q_brain * (c_art - c_bb) + psb * (fu_br * c_bm - fu_bb * c_bb) +
      psc * (fu_csf * c_ccsf - fu_bb * c_bb) + q_ssink * c_scsf + q_csink * c_ccsf +
      q_pgp * c_bm * fu_br
    d/dt(brain) <- psb * (fu_bb * c_bb - fu_br * c_bm) + pse * (fu_csf * c_ccsf - fu_br * c_bm) -
      q_bulk * c_bm - q_pgp * c_bm * fu_br
    d/dt(brain_csf_sas_cranial) <- pse * (fu_br * c_bm - fu_csf * c_ccsf) + psc * (fu_bb * c_bb - fu_csf * c_ccsf) +
      q_bulk * c_bm + q_sout * c_scsf - q_sin * c_ccsf - q_csink * c_ccsf
    d/dt(brain_csf_sas_spinal) <- q_sin * c_ccsf - q_sout * c_scsf - q_ssink * c_scsf

    # =====================================================================
    # M6G ODEs (S1 File lines 111-128). Formed in venous blood from the
    # cleared morphine mass and cleared from venous blood (same placement
    # as morphine, see above); no Pgp transport.
    # =====================================================================
    d/dt(venous_m6g) <- cl * c_ven * f_form +
      q_adipose * m_ad / kp_adipose_m6g * bp_m6g + q_bone * m_bo / kp_bone_m6g * bp_m6g +
      q_heart * m_he / kp_heart_m6g * bp_m6g + q_kidney * m_ki / kp_kidney_m6g * bp_m6g +
      q_muscle * m_mu / kp_muscle_m6g * bp_m6g + q_skin * m_sk / kp_skin_m6g * bp_m6g +
      q_liver * m_li / kp_liver_m6g * bp_m6g + q_brain * m_bb + q_rest * m_ot / kp_other_m6g * bp_m6g -
      q_lung * m_ven - cl_m6g * m_ven
    d/dt(arterial_m6g) <- q_lung * m_lung / kp_lung_m6g * bp_m6g -
      (q_rest + q_brain + q_adipose + q_bone + q_heart + q_kidney + q_muscle + q_skin + q_spleen + q_gut + q_ha) * m_art
    d/dt(lung_m6g) <- q_lung * (m_ven - m_lung / kp_lung_m6g * bp_m6g)
    d/dt(adipose_m6g) <- q_adipose * (m_art - m_ad / kp_adipose_m6g * bp_m6g)
    d/dt(bone_m6g) <- q_bone * (m_art - m_bo / kp_bone_m6g * bp_m6g)
    d/dt(heart_m6g) <- q_heart * (m_art - m_he / kp_heart_m6g * bp_m6g)
    d/dt(kidney_m6g) <- q_kidney * (m_art - m_ki / kp_kidney_m6g * bp_m6g)
    d/dt(muscle_m6g) <- q_muscle * (m_art - m_mu / kp_muscle_m6g * bp_m6g)
    d/dt(skin_m6g) <- q_skin * (m_art - m_sk / kp_skin_m6g * bp_m6g)
    d/dt(spleen_m6g) <- q_spleen * (m_art - m_sp / kp_spleen_m6g * bp_m6g)
    d/dt(other_m6g) <- q_rest * (m_art - m_ot / kp_other_m6g * bp_m6g)
    d/dt(gut_m6g) <- q_gut * (m_art - m_gu / kp_gut_m6g * bp_m6g)
    d/dt(liver_m6g) <- q_ha * m_art + q_gut * m_gu / kp_gut_m6g * bp_m6g + q_spleen * m_sp / kp_spleen_m6g * bp_m6g -
      q_liver * m_li / kp_liver_m6g * bp_m6g
    d/dt(brain_vascular_m6g) <- q_brain * (m_art - m_bb) + psb_m6g * (fu_br_m6g * m_bm - fu_bb_m6g * m_bb) +
      psc_m6g * (fu_csf_m6g * m_ccsf - fu_bb_m6g * m_bb) + q_ssink * m_scsf + q_csink * m_ccsf
    d/dt(brain_m6g) <- psb_m6g * (fu_bb_m6g * m_bb - fu_br_m6g * m_bm) + pse * (fu_csf_m6g * m_ccsf - fu_br_m6g * m_bm) -
      q_bulk * m_bm
    d/dt(brain_csf_sas_cranial_m6g) <- pse * (fu_br_m6g * m_bm - fu_csf_m6g * m_ccsf) +
      psc_m6g * (fu_bb_m6g * m_bb - fu_csf_m6g * m_ccsf) +
      q_bulk * m_bm + q_sout * m_scsf - q_sin * m_ccsf - q_csink * m_ccsf
    d/dt(brain_csf_sas_spinal_m6g) <- q_sin * m_ccsf - q_sout * m_scsf - q_ssink * m_scsf

    # =====================================================================
    # PD (Methods Eqs 4-6; S1 File lines 964-983). Unbound brain-mass
    # concentrations drive competitive receptor binding; KM in mg/L is the
    # nM constant times the molecular weight.
    # =====================================================================
    cu_brain <- c_bm * fu_br
    cu_brain_m6g <- m_bm * fu_br_m6g
    km <- km_nm * mw / 1e6
    km_m6g <- km_m6g_nm * mw_m6g / 1e6
    br_pct <- 100 * cu_brain / (cu_brain + km * (1 + cu_brain_m6g / km_m6g))
    br_pct_m6g <- 100 * cu_brain_m6g / (cu_brain_m6g + km_m6g * (1 + cu_brain / km))
    rr <- br_pct^hill / (br_pct^hill + br50^hill)
    rr_m6g <- br_pct_m6g^hill_m6g / (br_pct_m6g^hill_m6g + br50_m6g^hill_m6g)
    effect_rel <- rr + rr_m6g

    # =====================================================================
    # Outputs. Plasma = venous blood / blood:plasma ratio (S1 File lines
    # 837 and 896). Cbrain_u is the unbound brain-mass concentration the
    # paper compares with microdialysis ECF data.
    # =====================================================================
    Cbrain_u <- cu_brain
    Ccsf <- c_ccsf
    Ccsf_spinal <- c_scsf
    Cbrain_u_m6g <- cu_brain_m6g
    Ccsf_m6g <- m_ccsf
    Ccsf_spinal_m6g <- m_scsf
    Cc <- c_ven / bp
    Cc_m6g <- m_ven / bp_m6g
    Cc ~ prop(propSd)
  })
}
