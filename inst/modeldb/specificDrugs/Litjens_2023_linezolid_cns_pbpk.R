Litjens_2023_linezolid_cns_pbpk <- function() {
  description <- paste(
    "PBPK (permeability-limited 4-compartment CNS model, adapted from the",
    "Simcyp Simulator V19R1 brain module). Linezolid disposition in plasma,",
    "brain blood, brain mass, cranial cerebrospinal fluid (CSF) and spinal",
    "CSF in adults and children, developed to predict cranial CSF exposure",
    "and AUC0-24:MIC target attainment in tuberculous meningitis (Litjens",
    "et al. 2023, Antibiotics). The four CNS ODEs, the CNS physiological",
    "volumes and the CSF/brain fluid flows are those of the upstream",
    "framework paper Verscheijden et al. 2019 (PLoS Comput Biol; Eqs 2-5,",
    "S1 Table), which Litjens 2023 adapted with permission; the linezolid",
    "drug-specific inputs (Vss, plasma clearance, ka, fa, B:P, fu, fu-brain,",
    "PSB, PSC, PSE, and the BCRP / P-gp efflux clearances measured in-house",
    "in MDCKII monolayers) are from Litjens 2023 Supplementary Table S1 and",
    "S2. Adult and paediatric brain physiology branch on AGE (< 18 y uses",
    "the paediatric equations of Verscheijden 2019 S1 Table). IMPORTANT",
    "DEVIATION: the systemic side of the published model is the proprietary",
    "Simcyp full-PBPK whole-body distribution model (Method 2, Rodgers and",
    "Rowland) whose tissue:plasma partition coefficients are not reported",
    "anywhere in the paper or its supplement; here it is replaced by a",
    "single well-stirred plasma compartment parameterised with the exact",
    "aggregate values the paper does report - Vss (L/kg) and the per-study",
    "total plasma clearance carried in the Simcyp 'additional clearance'",
    "slot. Plasma AUC (= Dose * F / CL, the quantity every PK-PD conclusion",
    "of the paper rests on) is preserved exactly; the distribution phase of",
    "the plasma curve is approximate. See the vignette Assumptions and",
    "deviations section. Typical-value forward-simulation model: the paper",
    "reports no IIV or residual-error estimates."
  )
  reference <- paste(
    "Litjens CHC, Verscheijden LFM, Svensson EM, van den Broek PHH,",
    "van Hove H, Koenderink JB, Russel FGM, Aarnoutse RE, te Brake LHM.",
    "Physiologically-Based Pharmacokinetic Modelling to Predict the",
    "Pharmacokinetics and Pharmacodynamics of Linezolid in Adults and",
    "Children with Tuberculous Meningitis. Antibiotics. 2023;12(4):702.",
    "doi:10.3390/antibiotics12040702.",
    "Structural CNS framework and physiological parameters:",
    "Verscheijden LFM, Koenderink JB, de Wildt SN, Russel FGM. Development",
    "of a physiologically-based pharmacokinetic pediatric brain model for",
    "prediction of cerebrospinal fluid drug concentrations and the",
    "influence of meningitis. PLoS Comput Biol. 2019;15(6):e1007117.",
    "doi:10.1371/journal.pcbi.1007117.",
    sep = " "
  )
  vignette <- "Litjens_2023_linezolid_tuberculous_meningitis"

  # The four CNS states are the paper-mechanistic compartments of the
  # Gaohua 2016 / Verscheijden 2019 permeability-limited brain model. None
  # of them exists as a canonical brain_<region> entry in
  # inst/references/compartment-names.md: brain_blood and brain_mass are
  # the vascular and parenchymal halves of the blood-brain barrier, and the
  # two CSF states are anatomically distinct (cranial vs spinal), following
  # the csf_<site> idiom already used by Westerhout_2012_acetaminophen_rat_pbpk.
  paper_specific_compartments <- c(
    "brain_blood", "brain_mass", "csf_cranial", "csf_spinal"
  )

  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "mg/L"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    depot       = list(analyte = "Linezolid", units = "mg", specimen = "administration site", verified = FALSE),
    central     = list(analyte = "Linezolid", units = "mg", specimen = "plasma", verified = FALSE),
    brain_blood = list(analyte = "Linezolid", units = "mg", specimen = "tissue", verified = FALSE),
    brain_mass  = list(analyte = "Linezolid", units = "mg", specimen = "tissue", verified = FALSE),
    csf_cranial = list(analyte = "Linezolid", units = "mg", specimen = "CSF", verified = FALSE),
    csf_spinal  = list(analyte = "Linezolid", units = "mg", specimen = "CSF", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = paste(
        "Total body weight. Scales the steady-state volume of distribution",
        "(V = Vss * WT), the adult total brain volume, and the paediatric",
        "spinal CSF volume."
      ),
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Litjens 2023 Supplementary File S1: adult TB patients had a mean",
        "reported weight of ~60 kg but the Simcyp virtual population mean",
        "is ~75 kg, which is why Vss was re-expressed per kg (0.47 L/kg at",
        "60 kg -> 0.4 L/kg at 75 kg). Supply the weight of the simulated",
        "individual; the model does the scaling."
      ),
      source_name        = "body weight"
    ),
    AGE = list(
      description        = paste(
        "Age. Selects the adult (>= 18 y) or paediatric (< 18 y) brain",
        "physiology equation set of Verscheijden 2019 S1 Table and enters",
        "the cardiac-output, brain-blood-flow, brain-volume and CSF-",
        "production-rate equations directly."
      ),
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The paediatric CSF production rate has a further branch below",
        "3 months (0.25 y). Litjens 2023 simulated children aged 0.25-21 y",
        "(critically ill validation) and 0.6-9.4 y (paediatric MDR-TB",
        "target attainment); adults 18-78 y."
      ),
      source_name        = "age"
    ),
    SEXF = list(
      description        = paste(
        "Female sex indicator. Selects the sex-specific adult total-CSF",
        "fraction of brain volume (10.5 percent male, 9.2 percent female;",
        "Verscheijden 2019 S1 Table adult column)."
      ),
      units              = "(binary)",
      type               = "categorical",
      reference_category = "0 = male",
      notes              = paste(
        "1 = female, 0 = male. Only used in the adult branch; the",
        "paediatric branch uses a fixed cranial CSF volume of 0.143 L and a",
        "weight-based spinal CSF volume that are not sex-specific."
      ),
      source_name        = "sex"
    ),
    BSA = list(
      description        = paste(
        "Body surface area. Drives cardiac output, and hence brain blood",
        "flow (Qbrain = 0.12 * cardiac output in adults)."
      ),
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Verscheijden 2019 S1 File computes BSA with the Du Bois formula",
        "BSA = 0.007184 * HT^0.725 * WT^0.425 (HT in cm, WT in kg); supply",
        "the same value here."
      ),
      source_name        = "body surface area"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 131L,
    n_studies      = 10L,
    age_range      = "0.25-78 years (adults 18-78 y; children 0.25-21 y)",
    weight_range   = paste(
      "not tabulated per subject; adult TB patients mean ~60 kg,",
      "Simcyp virtual adult population mean ~75 kg"
    ),
    race_ethnicity = c(White = 100),
    disease_state  = paste(
      "Two distinct simulated populations. (1) Critically ill adults and",
      "children with non-tuberculous CNS conditions (subarachnoid /",
      "intraventricular / intracerebral haemorrhage, traumatic brain",
      "injury, hydrocephalus drainage, brain tumour excision, proven or",
      "suspected neurosurgical CNS infection), used to verify predicted",
      "cranial CSF exposure against measured plasma + CSF data",
      "(Supplementary Table S2, 9 study arms). (2) Adults and children with",
      "drug-sensitive or multidrug-resistant pulmonary tuberculosis, whose",
      "published plasma PK was used to predict cranial CSF exposure and",
      "AUC0-24:MIC target attainment for tuberculous meningitis",
      "(Supplementary Table S2, 5 study arms). Plasma clearance in",
      "tuberculous meningitis patients was assumed equal to that in",
      "pulmonary TB patients (Discussion, limitation 3)."
    ),
    dose_range     = paste(
      "Verification: 600 mg IV single dose and 600 mg IV BID in adults;",
      "10 mg/kg IV single dose and BID in children. Target attainment:",
      "300 mg oral BID, 600 mg oral BID and 1200 mg oral QD in adults;",
      "~9.24 mg/kg oral BID in children."
    ),
    regions        = paste(
      "Simulated as Simcyp healthy North European Caucasian volunteers",
      "(adult) and the Simcyp paediatric population; source clinical data",
      "from Austria, Italy, Greece, Spain, the Netherlands, South Africa",
      "and the USA."
    ),
    n_virtual      = paste(
      "100 virtual individuals per simulation (Supplementary File S1);",
      "age range, proportion female and dose matched to the clinical",
      "data set being reproduced. Ages above 65 y were capped at 65 y in",
      "the simulation (Supplementary Table S2 footnote)."
    ),
    notes          = paste(
      "No individual-level data were fitted in this paper: every simulation",
      "is a forward prediction. n_subjects = 131 is the number of patients",
      "in the 10 published studies (13 dosing arms) whose plasma and CSF",
      "concentration data were digitised (WebPlotDigitizer 4.1) for model",
      "verification: 65 critically ill (Beer 2007 n=5, Viaggi 2011 n=7,",
      "Tsona 2010 n=18, Yogev 2010 n=10, Luque 2014 n=11, Myrianthefs 2006",
      "n=14) and 66 with tuberculosis (Alffenaar 2010 Ther Drug Monit n=14,",
      "Alffenaar 2010 Clin Pharmacokinet n=8, Garcia-Prats 2019 n=13,",
      "Diacon 2020 n=15 + n=16). Male proportions per arm range 39-93",
      "percent (Supplementary Table S2); no pooled female percentage is",
      "reported.",
      "Litjens 2023 Supplementary Table S2 lists the 13 clinical study arms",
      "(Beer 2007, Viaggi 2011, Tsona 2010, Yogev 2010, Luque 2014,",
      "Myrianthefs 2006, Alffenaar 2010 x2, Garcia-Prats 2019, Diacon 2020",
      "x2) with the plasma clearance applied to each. Clearance could not",
      "be modelled mechanistically (linezolid elimination is ~35 percent",
      "renal and ~65 percent non-enzymatic chemical oxidation), so the",
      "reported total plasma clearance of each study was entered into the",
      "Simcyp 'enzyme kinetics - additional clearance' slot; the same",
      "approach is used here via the lcl parameter. Reported clearances:",
      "3.95 L/h (adult MDR-TB, Alffenaar 2010), 4.3 L/h (paediatric MDR-TB,",
      "Garcia-Prats 2019), 4.9 L/h (adult DS-TB 1200 mg QD, Diacon 2020),",
      "5.46 L/h (adult DS-TB 300 mg BID, Diacon 2020), 7.3-16.6 L/h",
      "(critically ill adults), 11-13 L/h (critically ill children).",
      "NOTE on the paediatric clearance: 4.3 L/h applied literally to a",
      "10-30 kg child gives a plasma AUC0-24 of ~40-60 mg*h/L, not the 202",
      "mg*h/L reported in Table 1. Supplementary File S1 states paediatric",
      "clearance 'was optimized by visual inspection'. Interpreting 4.3 L/h",
      "as a 70-kg-normalised typical value scaled allometrically,",
      "CL_i = 4.3 * (WT/70)^0.75, reproduces Table 1 to within a few",
      "percent across the 0.6-9.4 y range; this interpretation is used in",
      "the vignette and is flagged there as an inference, not a statement",
      "of the paper. Supply lcl as the INDIVIDUAL clearance in L/h."
    )
  )

  ini({
    # =========================================================================
    # All parameters are FIXED. Litjens 2023 is a forward-simulation PBPK
    # study: no parameter was estimated from individual-level data in this
    # paper, and neither IIV nor a residual-error model is reported.
    #
    # Default scenario = the adult pulmonary-TB oral simulation that produces
    # the paper's headline result (600 mg BID -> plasma AUC0-24:MIC 281,
    # cranial CSF AUC0-24:MIC 181; Table 1 row 2), i.e. the Alffenaar 2010
    # clearance with the TB volume of distribution and first-order
    # absorption. See the vignette for the parameter swaps that reproduce
    # the critically-ill IV arms and the paediatric arm.
    # =========================================================================

    # ---- Absorption (Supplementary Table S1, "Absorption: First-Order
    # Absorption Model"). The first-order model with ka = 0.9 /h was applied
    # to the Alffenaar 2010 and Garcia-Prats 2019 simulations; the Simcyp
    # ADAM model (ka = 1.20 /h, fa = 0.96) was applied to Diacon 2020.
    lka   <- fixed(log(0.9))   ; label("First-order absorption rate constant (1/h)")               # Suppl. Table S1, Absorption first-order: ka = 0.9 (Alffenaar 2010)
    fa    <- fixed(1)          ; label("Fraction available for absorption from dosage form")       # Suppl. Table S1: fa = 1 (Simcyp default); Suppl. Table S2 footnote f: bioavailability considered 100%

    # ---- Systemic disposition. The paper's Simcyp full-PBPK distribution
    # model is collapsed here to a single well-stirred plasma compartment of
    # volume Vss * WT (see description). Vss for the TB simulations was set
    # to 0.4 L/kg; the critically-ill simulations used 0.66 L/kg.
    lvss  <- fixed(log(0.4))   ; label("Steady-state volume of distribution (L/kg)")               # Suppl. Table S1, Distribution: "Optimized for TB studies to 0.4"; 0.66 [16] for critically ill
    lcl   <- fixed(log(3.95))  ; label("Total plasma clearance (L/h)")                             # Suppl. Table S2, Alffenaar 2010 adult MDR-TB arms: Cl 3.95; entered as Simcyp 'additional clearance'

    # ---- Blood binding (Supplementary Table S1, "Blood binding properties")
    bp    <- fixed(0.603)      ; label("Blood:plasma concentration ratio")                         # Suppl. Table S1: Blood-to-plasma ratio 0.603 [14, 15]
    fu    <- fixed(0.69)       ; label("Fraction unbound in plasma")                               # Suppl. Table S1: Fraction unbound in plasma 0.69 [14, 15]

    # ---- Brain model unbound fractions (Supplementary Table S1, "Brain model")
    fu_br  <- fixed(0.81)      ; label("Fraction unbound in brain mass")                           # Suppl. Table S1, BBB: fu,br = 0.81 (predicted by Simcyp)
    fu_csf <- fixed(1)         ; label("Fraction unbound in cerebrospinal fluid")                  # Suppl. Table S1, Blood-cranial CSF barrier: fuCSF = 1 [23]

    # ---- Permeability-surface-area products (Supplementary Table S1,
    # "Brain model"). PSB was optimised on the Viaggi 2011 single-dose study
    # and not changed for any multiple-dose simulation; PSC was assumed to be
    # half of PSB; PSE is large enough to represent "no barrier" between
    # brain mass and cranial CSF.
    lpsb  <- fixed(log(1.0))   ; label("Permeability-surface-area product, blood-brain barrier (L/h)")            # Suppl. Table S1, BBB: PSB = 1.0, optimized [6]
    lpsc  <- fixed(log(0.5))   ; label("Permeability-surface-area product, blood-cranial-CSF barrier (L/h)")      # Suppl. Table S1, Blood-cranial CSF barrier: PSC = 0.5, assumed half of PSB [7]
    lpse  <- fixed(log(300))   ; label("Permeability-surface-area product, brain mass to cranial CSF (L/h)")      # Suppl. Table S1, Brain-CSF barrier: PSE = 300 [7]

    # ---- Blood-brain-barrier efflux transporters. In-vitro efflux
    # clearances were measured in-house in MDCKII-BCRP and MDCKII-P-gp
    # monolayers (Supplementary File S2, Eq. S2) and scaled to whole-brain
    # in-vivo clearance with the relative activity factors of Eq. S3.
    # Efflux at the blood-CSF barrier was NOT incorporated (Suppl. File S2).
    cl_bcrp_vitro <- fixed(16)  ; label("In vitro BCRP (ABCG2) efflux clearance (uL/min/mg protein)")   # Suppl. Table S1: CL_ABCG2,vitro = 16; MDCK exp Eq 2. Results 3.1: net efflux ratio 3.59
    cl_pgp_vitro  <- fixed(2.1) ; label("In vitro P-gp (ABCB1) efflux clearance (uL/min/mg protein)")   # Suppl. Table S1: CL_ABCB1,vitro = 2.1; MDCK exp Eq 2. Results 3.1: net efflux ratio 1.43
    raf_bcrp      <- fixed(143) ; label("Relative activity factor, BCRP (mg microvessel protein)")      # Suppl. Table S1: RAF for ABCG2 = 143; Suppl. File S2 Eq. S3 = (5.5/13.18) * 0.244 * 1400
    raf_pgp       <- fixed(140) ; label("Relative activity factor, P-gp (mg microvessel protein)")      # Suppl. Table S1: RAF for ABCB1 = 140; Suppl. File S2 Eq. S3 = (4.21/10.3) * 0.244 * 1400

    # ---- Placeholder residual error. Litjens 2023 fitted no error model
    # (all simulations are forward predictions from fixed inputs); this entry
    # exists so rxode2 / nlmixr2 stochastic simulation is syntactically valid.
    propSd <- fixed(0.10) ; label("Proportional residual error placeholder (fraction)")  # not reported by the paper; placeholder for forward simulation
  })

  model({
    # =======================================================================
    # CNS physiology. Verscheijden 2019 S1 Table ("Physiological
    # parameters"), which gives a paediatric and an adult column for every
    # quantity; Litjens 2023 Section 2.1 states the paediatric population
    # "incorporates age-related differences in brain volume, brain blood
    # flow, CSF production rate, CSF volumes and blood-brain barrier and
    # blood-cerebrospinal fluid barrier surface area".
    # =======================================================================
    paed <- (AGE < 18)

    # Total brain volume (L). Adult: (1.449 - 3.62/WT)/1.04. Paediatric:
    # (10*(AGE+0.315)/(9+6.92*AGE))/1.04. 1.04 kg/L is the brain density
    # used to convert the published brain-weight equations to volumes.
    v_brain <- ifelse(paed,
                      (10 * (AGE + 0.315) / (9 + 6.92 * AGE)) / 1.04,   # Verscheijden 2019 S1 Table, Vbrain total, paediatric column
                      (1.449 - 3.62 / WT) / 1.04)                       # Verscheijden 2019 S1 Table, Vbrain total, adult column

    v_bb   <- 0.05 * v_brain     # Verscheijden 2019 S1 Table: Vbrain blood = 0.05 * Vbrain total
    v_endo <- 0.005 * v_brain    # Verscheijden 2019 S1 Table: Vendothelial mass = Vbrain * 0.005

    # Total CSF as a fraction of brain volume differs by sex in the adult
    # column (10.5% male, 9.2% female); 80% of it is cranial and 20% spinal.
    csf_frac  <- ifelse(SEXF > 0.5, 0.092, 0.105)   # Verscheijden 2019 S1 Table, Vccsf / Vscsf adult column
    v_ccsf_ad <- v_brain * csf_frac * 0.8
    v_scsf_ad <- v_brain * csf_frac * 0.2

    # Paediatric: cranial CSF is a constant 0.143 L; spinal CSF is weight-
    # based and capped at 20% of total CSF (the same 20% split as the adult).
    v_scsf_pd_raw <- (1.94 * WT + 0.13) / 1000                          # Verscheijden 2019 S1 Table, Vscsf paediatric column
    v_scsf_pd_cap <- (0.143 / 0.8) * 0.2                                # Verscheijden 2019 S1 Table: Limit Vscsf <= (0.143/0.8)*0.2
    v_scsf_pd     <- ifelse(v_scsf_pd_raw > v_scsf_pd_cap, v_scsf_pd_cap, v_scsf_pd_raw)

    v_ccsf <- ifelse(paed, 0.143, v_ccsf_ad)   # Verscheijden 2019 S1 Table, Vccsf
    v_scsf <- ifelse(paed, v_scsf_pd, v_scsf_ad)

    # Brain mass is the residual of total brain volume.
    v_bm <- v_brain - v_endo - v_bb - v_ccsf - v_scsf   # Verscheijden 2019 S1 Table: Vbrain mass = Vbrain - Vbrainendo - Vbb - Vccsf - Vscsf

    # ---- CSF fluid flows (L/h). Verscheijden 2019 S1 Table "Fluid flow
    # rates"; the relative CSF flows (as fractions of the CSF production
    # rate) are assumed identical for adults and children.
    # log10(AGE) is written as log(AGE)/log(10) for rxode2 portability.
    q_csf_prod <- ifelse(paed,
                         ifelse(AGE < 0.25,
                                (4.007 * log(AGE) / log(10) + 7.088) / 1000,  # Verscheijden 2019 S1 Table: Qproductionrate < 3 months
                                0.024),                                       # Verscheijden 2019 S1 Table: Qproductionrate 3m-18y = 0.024
                         0.021)                                               # Verscheijden 2019 S1 Table: Qproductionrate adult = 0.021

    q_bulk  <- 0.25 * q_csf_prod                          # Verscheijden 2019 S1 Table: Qbulk = 0.25 * Qproductionrate (brain mass -> cranial CSF)
    q_ssink <- 0.38 * (0.75 * q_csf_prod + q_bulk)        # Verscheijden 2019 S1 Table: Qssink = 0.38*(0.75*Qproductionrate + Qbulk) (spinal CSF -> blood)
    q_sout  <- 0.9 * q_ssink                              # Verscheijden 2019 S1 Table: Qsout = 0.9 * Qssink (spinal CSF -> cranial CSF)
    q_sin   <- q_ssink + q_sout                           # Verscheijden 2019 S1 Table: Qsin = Qssink + Qsout (cranial CSF -> spinal CSF)
    q_csink <- 0.75 * q_csf_prod + q_bulk - q_sin + q_sout # Verscheijden 2019 S1 Table: Qcsink (cranial CSF -> brain blood)

    # ---- Cardiac output and brain blood flow (L/h).
    q_co <- ifelse(paed,
                   BSA * (110 + 184.974 * (exp(-0.0378 * AGE) - exp(-0.2477 * AGE))),  # Verscheijden 2019 S1 Table, Qcarout paediatric column
                   BSA * 60 * (3 - 0.01 * (AGE - 20)))                                 # Verscheijden 2019 S1 Table, Qcarout adult column
    q_brain <- ifelse(paed,
                      q_co * (10 + 2290 * (exp(-0.608 * AGE) - exp(-0.639 * AGE))) / 100,  # Verscheijden 2019 S1 Table, Qbrain paediatric column
                      q_co * 0.12)                                                          # Verscheijden 2019 S1 Table, Qbrain adult column = Qcarout * 0.12

    # =======================================================================
    # Drug-specific transfer terms.
    # =======================================================================
    ka  <- exp(lka)
    cl  <- exp(lcl)
    v_c <- exp(lvss) * WT     # plasma-referenced volume of distribution (L)

    # Unbound fraction referenced to WHOLE BLOOD, because the brain-blood
    # state carries a blood (not plasma) concentration. Verscheijden 2019
    # S1 File line "Fubb = Fupl/1.09  # Fraction unbound blood (Fupl/BP)".
    # For linezolid B:P < 1 (drug excluded from erythrocytes), so fu_bb > fu.
    fu_bb <- fu / bp

    # Paediatric permeability-surface-area products are scaled by the ratio
    # of brain volume to the adult reference brain (1.36 kg / 1.04 kg/L =
    # 1.3077 L), which is how Verscheijden 2019 S2 Table propagates the
    # blood-brain-barrier surface area to children.
    v_brain_ref <- 1.36 / 1.04
    ps_scale <- ifelse(paed, v_brain / v_brain_ref, 1)   # Verscheijden 2019 S2 Table, paediatric column: PSb = Vbrain/(1.36/1.04) * PSb_adult
    psb <- exp(lpsb) * ps_scale
    psc <- exp(lpsc) * ps_scale
    pse <- exp(lpse)                                     # Verscheijden 2019 S2 Table: PSe = 300 in both columns (no barrier assumed)

    # Whole-brain efflux-transporter clearance at the blood-brain barrier.
    # Litjens Suppl. File S2 Eq. S3: CLefflux,vivo = CLefflux,vitro * RAF,
    # with CLefflux,vitro in uL/min/mg and RAF in mg, giving uL/min; the
    # 60/1e6 factor converts uL/min to L/h. BCRP: 16 * 143 = 2288 uL/min;
    # P-gp: 2.1 * 140 = 294 uL/min; total 2582 uL/min = 0.15492 L/h. The RAF
    # carries an explicit brain-weight factor (1400 g), so the paediatric
    # value is scaled with the same brain-volume ratio as PSB / PSC.
    cl_bout <- (cl_bcrp_vitro * raf_bcrp + cl_pgp_vitro * raf_pgp) * 60 / 1e6 * ps_scale

    # =======================================================================
    # Concentrations. All states hold AMOUNT of linezolid (mg).
    # =======================================================================
    Cp    <- central / v_c        # plasma concentration (mg/L)
    Cbla  <- Cp * bp              # arterial blood concentration entering the brain (mg/L)
    Cbb   <- brain_blood / v_bb   # brain blood concentration (mg/L)
    Cbm   <- brain_mass / v_bm    # brain mass concentration (mg/L)
    Cccsf <- csf_cranial / v_ccsf # cranial CSF concentration (mg/L)
    Cscsf <- csf_spinal / v_scsf  # spinal CSF concentration (mg/L)

    # =======================================================================
    # ODE system. The four CNS equations are Verscheijden 2019 Eqs 2-5
    # verbatim (multiplied through by the compartment volume so the states
    # are amounts), plus the BBB efflux-transporter term CLbout that Litjens
    # 2023 added (Figure 1: "CLbout represents overall clearance of efflux
    # transporters P-gp and BCRP expressing at BBB"), which moves drug from
    # brain mass back to brain blood and therefore appears with opposite
    # signs in the two equations.
    # =======================================================================
    d/dt(depot) <- -ka * depot

    # Systemic plasma compartment. Loses drug to elimination and to the
    # brain-blood exchange, and receives the brain venous return so that the
    # CNS sub-model is mass-balanced against the body.
    d/dt(central) <- ka * depot * fa - cl * Cp - q_brain * Cbla + q_brain * Cbb

    # Verscheijden 2019 Eq. 2 (brain blood) + CLbout influx.
    d/dt(brain_blood) <- q_brain * (Cbla - Cbb) +
      psb * (fu_br * Cbm - fu_bb * Cbb) +
      psc * (fu_csf * Cccsf - fu_bb * Cbb) +
      q_ssink * Cscsf + q_csink * Cccsf +
      cl_bout * fu_br * Cbm

    # Verscheijden 2019 Eq. 3 (brain mass) - CLbout efflux.
    d/dt(brain_mass) <- psb * (fu_bb * Cbb - fu_br * Cbm) +
      pse * (fu_csf * Cccsf - fu_br * Cbm) -
      q_bulk * Cbm -
      cl_bout * fu_br * Cbm

    # Verscheijden 2019 Eq. 4 (cranial CSF).
    d/dt(csf_cranial) <- pse * (fu_br * Cbm - fu_csf * Cccsf) +
      psc * (fu_bb * Cbb - fu_csf * Cccsf) +
      q_bulk * Cbm + q_sout * Cscsf -
      q_sin * Cccsf - q_csink * Cccsf

    # Verscheijden 2019 Eq. 5 (spinal CSF).
    d/dt(csf_spinal) <- q_sin * Cccsf - q_sout * Cscsf - q_ssink * Cscsf

    # =======================================================================
    # Outputs. Cc is the total (bound + unbound) plasma concentration that
    # the paper compares against measured plasma / serum data and from which
    # the plasma AUC0-24:MIC index is computed. Ccsf is the cranial CSF
    # concentration driving the paper's central conclusion; Cspinalcsf and
    # Cbrain reproduce Supplementary Figures S2-S3.
    # =======================================================================
    Ccsf       <- Cccsf
    Cspinalcsf <- Cscsf
    Cbrain     <- Cbm
    Cc         <- Cp
    Cc ~ prop(propSd)
  })
}
