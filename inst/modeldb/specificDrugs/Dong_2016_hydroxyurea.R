Dong_2016_hydroxyurea <- function() {
  description <- paste(
    "One-compartment population PK model for oral hydroxyurea in pediatric",
    "patients with sickle cell anaemia (Dong 2016, HUSTLE trial",
    "NCT00305175; n = 96 children aged 1.2-16.6 years on a 20 mg/kg",
    "starting dose, all African-American). Saturable Michaelis-Menten",
    "elimination from the central compartment (Vmax 490 mg/h per 70 kg,",
    "Km fixed at 25 mg/L based on a prior report) with allometric scaling",
    "of Vmax (exponent 0.75 fixed) and apparent central volume V/F",
    "(exponent 1.0 fixed) on total body weight (reference 70 kg).",
    "Cystatin C is a power-model covariate on Vmax (exponent -0.509,",
    "reference 0.74 mg/L), and was the only covariate retained over",
    "serum creatinine, eGFR, and direct 99mTc-DTPA-measured GFR. Oral",
    "absorption is described by a Savic 2007 transit-compartment chain",
    "(NN = 12.4 transit compartments fitted, MTT = 0.158 h) feeding a",
    "depot with first-order absorption Ka = 8.19 /h. Residual error is",
    "combined additive (0.117 mg/L) and proportional (39.7% CV)."
  )
  reference <- paste(
    "Dong M, McGann PT, Mizuno T, Ware RE, Vinks AA (2016). Development",
    "of a pharmacokinetic-guided dose individualization strategy for",
    "hydroxyurea treatment in children with sickle cell anaemia.",
    "Br J Clin Pharmacol 81(4):742-752.",
    "doi:10.1111/bcp.12851.",
    sep = " "
  )
  vignette <- "Dong_2016_hydroxyurea"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline body weight (kg). Drives the allometric scaling of",
        "Vmax (exponent 0.75 fixed) and the apparent central volume",
        "V/F (exponent 1.0 fixed); both scalings use 70 kg as the",
        "reference weight per Dong 2016 Results 'Population pharmacokinetic",
        "modelling' paragraph 4 ('We used fixed theoretical values to",
        "scale allometrically the influence of body size on hydroxyurea",
        "PK') and Table 2 'Vmax = theta1 * (WT/70)^0.75' and 'V/F = theta3",
        "* (WT/70)'. Cohort median 26.6 kg, range 8.8-88.3 kg on day 1",
        "(Dong 2016 Table 1)."
      ),
      source_name        = "WT"
    ),
    CYSC = list(
      description        = "Serum cystatin C concentration",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline serum cystatin C (turbidimetric assay with",
        "anti-cystatin-C-coated latex particles, Roche Diagnostics).",
        "Used as a power-model covariate on Vmax: Vmax_typ * (CYSC /",
        "0.74)^(-0.509). Reference 0.74 mg/L is the cohort median per",
        "Dong 2016 Table 1 ('Cystatin C, mg/L: 0.74 (0.57-1.25)') and",
        "Table 2 footnote ('typical patient with cystatin C level of",
        "0.74 mg/L'). Cystatin C was the only renal-function covariate",
        "retained in the final model after the backward elimination step;",
        "serum creatinine and Cockcroft-Gault-style creatinine clearance",
        "(CLcr) reduced OFV by only 5.05 and 4.31 respectively (p > 0.01)",
        "vs the cystatin C drop of 24.5 OFV points (Dong 2016 Results",
        "'Population pharmacokinetic modelling' paragraph 4)."
      ),
      source_name        = "CysC"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened in the Dong 2016 stepwise covariate analysis but not",
        "retained in the final model (only body weight and cystatin C",
        "survived backward elimination). Cohort median 8.8 years, range",
        "1.2-16.6 years on day 1 (Dong 2016 Table 1). Documented here for",
        "provenance of the covariate screen; not referenced inside model()."
      ),
      source_name        = "AGE"
    ),
    BSA = list(
      description        = "Body surface area",
      units              = "m2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened in the Dong 2016 stepwise covariate analysis but not",
        "retained (Dong 2016 Methods 'Covariate analysis' first paragraph",
        "lists the tested covariates including BSA; only body weight and",
        "cystatin C were retained). Documented for screen provenance;",
        "not referenced inside model()."
      ),
      source_name        = "BSA"
    ),
    HT = list(
      description        = "Body height",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened in the Dong 2016 stepwise covariate analysis but not",
        "retained. Documented for screen provenance; not referenced",
        "inside model()."
      ),
      source_name        = "HT"
    ),
    SEXF = list(
      description        = "Sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Screened in the Dong 2016 stepwise covariate analysis but not",
        "retained. Cohort 61 male / 35 female on day 1 (Dong 2016 Table 1)."
      ),
      source_name        = "Gender"
    ),
    CREAT = list(
      description        = "Serum creatinine",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened but not retained; cystatin C outperformed serum",
        "creatinine as a clearance covariate. Serum creatinine reduced",
        "OFV by only 5.05 points (p > 0.01) vs 24.5 for cystatin C",
        "(Dong 2016 Results 'Population pharmacokinetic modelling'",
        "paragraph 4). Cohort median 0.3 mg/dL on day 1 (Table 1)."
      ),
      source_name        = "SCr"
    ),
    CRCL = list(
      description        = "Estimated creatinine clearance (Schwartz formula)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened but not retained. Schwartz-formula estimated CrCl",
        "reduced OFV by only 4.31 points (p > 0.01) vs 24.5 for cystatin",
        "C (Dong 2016 Results 'Population pharmacokinetic modelling'",
        "paragraph 4)."
      ),
      source_name        = "CLcr"
    ),
    AST = list(
      description        = "Aspartate aminotransferase",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened but not retained.",
      source_name        = "AST"
    ),
    ALT = list(
      description        = "Alanine aminotransferase",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened but not retained.",
      source_name        = "ALT"
    ),
    TBILI = list(
      description        = "Total bilirubin at baseline",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened but not retained. Canonical units standardized to SI",
        "umol/L per the 2026-06-19 canonical-register audit; the source",
        "paper reports total bilirubin in mg/dL (1 mg/dL = 17.1 umol/L).",
        "No inline conversion is needed because TBILI is an excluded",
        "covariate and is not referenced in model()/ini()."
      ),
      source_name        = "Bilirubin"
    ),
    BUN = list(
      description        = "Blood urea nitrogen",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened but not retained. Cohort median 7 mg/dL, range 3-21",
        "mg/dL on day 1 (Dong 2016 Table 1)."
      ),
      source_name        = "BUN"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 96L,
    n_studies       = 1L,
    age_range       = "1.2-16.6 years",
    age_median      = "8.8 years",
    weight_range    = "8.8-88.3 kg",
    weight_median   = "26.6 kg",
    sex_female_pct  = 100 * 35 / 96,
    race_ethnicity  = "African American (96 of 96 day-1 subjects, 100%; per Dong 2016 Table 1 'Race: African American / Others = 96 / 0')",
    disease_state   = paste(
      "Children and adolescents with sickle cell anaemia enrolled in",
      "the HUSTLE clinical trial (NCT00305175) at St. Jude Children's",
      "Research Hospital. All subjects received a single oral 20 mg/kg",
      "hydroxyurea dose on day 1 of the dose-escalation protocol; 63 of",
      "the 96 subjects contributed paired pharmacokinetic profiles",
      "after reaching individual maximum tolerated dose (MTD) over a",
      "minimum 3-month titration window."
    ),
    dose_range      = paste(
      "Day 1: a single 20 mg/kg oral hydroxyurea dose (range 18-30",
      "mg/kg/day per Table 1). At MTD: 14.2-35.5 mg/kg/day (median 23.4",
      "mg/kg/day, n = 63). Most subjects received a liquid hydroxyurea",
      "formulation; 10 patients early in the study received capsules",
      "for their day-1 dosing."
    ),
    regions         = "United States (single-centre: St. Jude Children's Research Hospital, Memphis, TN)",
    notes           = paste(
      "712 hydroxyurea plasma concentrations from 96 children with",
      "sickle cell anaemia. Plasma samples were collected predose and",
      "at 20 min, 40 min, 1, 2, 4, 6, and 8 h after oral hydroxyurea",
      "administration on both day 1 and at MTD. Hydroxyurea was assayed",
      "by a colorimetric technique. Baseline demographics from Dong 2016",
      "Table 1 day-1 column. The race distribution of the cohort is",
      "entirely African American (96/0) by self-report, consistent with",
      "the typical pediatric sickle cell anaemia population in the",
      "United States. Renal function was assessed by direct 99mTc-DTPA",
      "plasma clearance (cohort median direct GFR 153 mL/min, range",
      "77-308 mL/min), Schwartz-formula estimated creatinine clearance,",
      "and serum cystatin C; cystatin C was the only renal-function",
      "covariate retained in the final model."
    )
  )

  ini({
    # =========================================================================
    # Structural disposition (Dong 2016 Table 2 'Estimates' column).
    # Reference subject: WT 70 kg, CYSC 0.74 mg/L. Allometric exponents
    # 0.75 on Vmax and 1.0 on V/F are FIXED at theoretical values per
    # Dong 2016 Results 'Population pharmacokinetic modelling' paragraph 4
    # ('We used fixed theoretical values to scale allometrically the
    # influence of body size on hydroxyurea PK ... it is rational to apply
    # the fixed exponent of 0.75 and 1 in our model.').
    # =========================================================================
    lvmax <- log(490);   label("Maximum elimination rate Vmax (mg/h, at WT 70 kg and CYSC 0.74 mg/L reference)")  # Dong 2016 Table 2: theta1 = 490.0 mg/h (RSE 2%)
    lkm   <- fixed(log(25));  label("Michaelis-Menten constant Km (mg/L)")                                        # Dong 2016 Table 2: Km = 25 mg/L (Fixed); Results paragraph 3 ('Km was fixed to 25 mg/L based on a previous report'); literature reference [24]
    lvc   <- log(49.6);  label("Apparent central volume of distribution V/F (L, at WT 70 kg reference)")          # Dong 2016 Table 2: theta3 = 49.6 L (RSE 2%)

    # =========================================================================
    # Absorption (Dong 2016 Table 2; Results paragraph 3 'a transit
    # absorption model provided better flexibility'). Savic 2007 transit-
    # compartment chain with Ka, MTT, and NN reported separately. The
    # transit() function in rxode2 implements the analytical gamma-PDF
    # input rate for non-integer NN; the drug then absorbs from the depot
    # into central with first-order rate Ka.
    # =========================================================================
    lka  <- log(8.19);   label("First-order absorption rate constant Ka (1/h)")                                   # Dong 2016 Table 2: Ka = 8.19 /h (RSE 21%)
    lmtt <- log(0.158);  label("Mean absorption transit time MTT (h)")                                            # Dong 2016 Table 2: MTT = 0.158 h (RSE 15%)
    lnn  <- log(12.4);   label("Number of absorption transit compartments NN (continuous, dimensionless)")        # Dong 2016 Table 2: N = 12.4 (RSE 15%)

    lfdepot <- fixed(log(1));  label("Oral bioavailability F (anchored at 1; apparent CL/F and V/F parameterisation)")  # Dong 2016 reports apparent CL/F (=Vmax/Km in steady-state limit) and V/F; absolute F not estimable from oral-only data

    # =========================================================================
    # Allometric exponents (Dong 2016 Results 'Population pharmacokinetic
    # modelling' paragraph 4: 'fixed theoretical values to scale allometrically
    # the influence of body size on hydroxyurea PK' and 'it is rational to
    # apply the fixed exponent of 0.75 and 1 in our model').
    # =========================================================================
    allo_vmax <- fixed(0.75);  label("Allometric exponent on Vmax (unitless)")                                    # Dong 2016 Table 2: 'Vmax = theta1 * (WT/70)^0.75'; fixed at theoretical 0.75 per Anderson-Holford
    allo_vc   <- fixed(1.00);  label("Allometric exponent on V/F (unitless)")                                     # Dong 2016 Table 2: 'V/F = theta3 * (WT/70)'; fixed at theoretical 1.0 per Anderson-Holford

    # =========================================================================
    # Cystatin C effect on Vmax (Dong 2016 Table 2: power-model exponent
    # theta2 = -0.509, RSE 20%, on cystatin C normalised to the cohort
    # median 0.74 mg/L).
    # =========================================================================
    e_cysc_vmax <- -0.509;  label("Power exponent of (CYSC / 0.74) on Vmax (unitless)")                           # Dong 2016 Table 2: theta2 = -0.509 (RSE 20%); 'Vmax = theta1 * (WT/70)^0.75 * (cysC/0.74)^theta2'

    # =========================================================================
    # Inter-individual variability (Dong 2016 Table 2 'IIV (%CV)' column;
    # converted to log-normal variance via omega^2 = log(1 + CV^2)).
    # IIV on Km and on N are FIXED to 0 per Table 2 ('omega(Km) = 0, Fixed'
    # and 'omega(N) = 0, Fixed' / Results paragraph 3 'high correlation
    # between Km and Vmax values suggested that the data did not support
    # estimation of all the parameters'). Single etas -- no correlation
    # block reported in Table 2.
    # =========================================================================
    etalvmax ~ 0.020807   # Dong 2016 Table 2 'omega(Vmax) = 14.5 (%CV)'; omega^2 = log(1 + 0.145^2) = 0.020807
    etalvc   ~ 0.018061   # Dong 2016 Table 2 'omega(V/F) = 13.5 (%CV)'; omega^2 = log(1 + 0.135^2) = 0.018061
    etalka   ~ 0.764975   # Dong 2016 Table 2 'omega(Ka) = 107.2 (%CV)'; omega^2 = log(1 + 1.072^2) = 0.764975
    etalmtt  ~ 0.528000   # Dong 2016 Table 2 'omega(MTT) = 83.4 (%CV)'; omega^2 = log(1 + 0.834^2) = 0.528000

    # =========================================================================
    # Residual error (Dong 2016 Table 2 'RV' column). Combined additive +
    # proportional model on linear-scale plasma concentration.
    # =========================================================================
    addSd  <- 0.117;  label("Additive residual error (mg/L)")                                                     # Dong 2016 Table 2: sigma_add = 0.117 mg/L (RSE 5.4%)
    propSd <- 0.397;  label("Proportional residual error (fraction, CV)")                                         # Dong 2016 Table 2: sigma_prop = 39.7 %CV (RSE 5.7%)
  })

  model({
    # --- 1. Individual PK parameters ----------------------------------------
    # Vmax: typical-value 490 mg/h scaled allometrically to body weight
    # (exponent 0.75 fixed, reference 70 kg) and modified by a power-model
    # cystatin C effect (exponent -0.509, reference 0.74 mg/L).
    vmax <- exp(lvmax + etalvmax) * (WT / 70)^allo_vmax * (CYSC / 0.74)^e_cysc_vmax

    # Km: fixed at 25 mg/L (no IIV by paper design).
    km   <- exp(lkm)

    # Apparent central volume of distribution: allometric on body weight
    # (exponent 1.0 fixed, reference 70 kg).
    vc   <- exp(lvc + etalvc) * (WT / 70)^allo_vc

    # Absorption rate constant: log-normal IIV.
    ka   <- exp(lka + etalka)

    # Mean transit time: log-normal IIV.
    mtt  <- exp(lmtt + etalmtt)

    # Number of transit compartments: no IIV (fixed by paper to a single
    # typical-value estimate of 12.4).
    nn   <- exp(lnn)

    # Bioavailability: anchored at 1.0 (apparent CL/F and V/F parameterisation).
    fdepot <- exp(lfdepot)

    # --- 2. Observation variable (concentration) ----------------------------
    Cc <- central / vc

    # --- 3. ODE system ------------------------------------------------------
    # Savic 2007 transit-compartment chain feeds the depot via the rxode2
    # transit() function (analytical gamma-PDF input rate for non-integer
    # NN; bioavailability gating via the bio argument). Hydroxyurea then
    # absorbs from the depot into the central compartment with first-order
    # rate Ka. Elimination from central follows saturable Michaelis-Menten
    # kinetics: rate = Vmax * Cc / (Km + Cc); multiplied by Vc gives an
    # amount/time elimination rate from the central amount compartment.
    d/dt(depot)   <- transit(nn, mtt, fdepot) - ka * depot
    d/dt(central) <- ka * depot - vmax * Cc / (km + Cc)

    f(depot) <- 0  # suppress dose-event bolus into depot; transit() reads podo(depot) and delivers the input

    # --- 4. Residual error -------------------------------------------------
    # Dose mg, V/F L -> Cc in mg/L (matches Dong 2016 Table 1 / Figure 2
    # plot axes and Table 2 residual-error units).
    Cc ~ add(addSd) + prop(propSd)
  })
}
