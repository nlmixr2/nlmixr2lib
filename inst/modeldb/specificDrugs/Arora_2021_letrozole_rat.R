Arora_2021_letrozole_rat <- function() {
  description <- paste(
    "Preclinical (rat). PBPK-derived reduced model (Simcyp Animal Simulator V17)",
    "for letrozole in plasma and brain extracellular fluid (ECF) of male and",
    "female Sprague-Dawley rats after single and once-daily 4 mg/kg",
    "extravascular doses. The source paper built a bottom-up whole-body PBPK",
    "model with the Simcyp multi-compartment brain model (brain blood, brain",
    "mass, CSF); the rat organ volumes, blood flows, tissue compositions and the",
    "Rodgers-Rowland tissue partition coefficients are Simcyp library values",
    "that are not printed, so the whole-body structure is not reproduced. What",
    "the paper does print (Table 1) is a complete systemic card -- first-order",
    "absorption rate and fraction absorbed, the in vivo oral clearance entered",
    "as the clearance input, and the predicted steady-state volume of",
    "distribution -- plus the brain permeability-surface area product across the",
    "blood-brain barrier and the unbound fractions in plasma and brain. This",
    "file encodes a one-compartment first-order-absorption plasma model with",
    "V = Vss and CL/F = CLpo, and a brain-mass compartment exchanging unbound",
    "drug with plasma across the blood-brain barrier at the printed PSB; brain",
    "ECF concentration is the unbound brain concentration fu,brain times total",
    "brain concentration, exactly as the paper derives it. The CSF compartment",
    "(printed PSC and PSE, but unprinted CSF volume, bulk flow and CSF sink",
    "flow) and the brain blood compartment (unprinted volume and cerebral blood",
    "flow) are omitted. Clearance and absorption rate differ by sex (females",
    "clear letrozole about 3.7-fold more slowly); volume and brain parameters",
    "are shared. With zero fitted parameters, the reduction reproduces all 18",
    "PBPK-predicted plasma and brain ECF exposure metrics of Tables 4 and 5",
    "within 12% (14 within 5%). Typical-value model: the paper reports no",
    "interindividual or residual variability, so no etas are declared and all",
    "residual error terms are fixed at zero.",
    sep = " "
  )
  reference <- paste(
    "Arora P, Gudelsky G, Desai PB. Gender-based differences in brain and plasma",
    "pharmacokinetics of letrozole in sprague-dawley rats: Application of",
    "physiologically-based pharmacokinetic modeling to gain quantitative",
    "insights. PLoS ONE. 2021;16(4):e0248579. doi:10.1371/journal.pone.0248579.",
    "PMCID PMC8018653. Drug-specific PBPK inputs are Table 1 (with footnotes a",
    "and b giving the per-kg clearance and per-gram PSB they were derived from)",
    "and Eq 3 (Crone-Renkin PSB). PBPK-predicted versus observed NCA metrics are",
    "Table 4 (single dose) and Table 5 (steady state); observed NCA is Tables 2",
    "and 3; individual observed concentrations are Supporting Information Tables",
    "S1 and S2 (figshare doi:10.6084/m9.figshare.14184959.v1).",
    sep = " "
  )
  vignette <- "Arora_2021_letrozole_rat"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    SEXF = list(
      description = "Biological sex, female indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "1 = female, 0 = male. Selects the sex-specific first-order absorption",
        "rate and in vivo oral clearance of Table 1 (males 0.29 1/h and",
        "0.77 mL/min; females 0.49 1/h and 0.21 mL/min). The paper simulated",
        "and studied each sex separately.",
        sep = " "
      ),
      source_name = "males / females"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Scales clearance and volume linearly (per kg). Table 1 reports Vss in",
        "L/kg and footnote a derives the entered CLpo (mL/min) from a per-kg",
        "clearance (0.185 L/h/kg males, 0.051 L/h/kg females) at an assumed",
        "250 g body weight, so the per-kg scaling is the paper's own. Because",
        "doses are 4 mg/kg, plasma concentrations do not depend on WT. The",
        "brain parameters (PSB for an assumed 1.8 g brain) are absolute and are",
        "not scaled. Study rats weighed 201-225 g (females) and 301-325 g",
        "(males) at purchase.",
        sep = " "
      ),
      source_name = "body weight"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "letrozole", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "letrozole", units = "mg", specimen = "plasma", verified = TRUE),
    # Total drug in brain tissue; brain ECF (microdialysate) is fu_brain times
    # its concentration.
    brain = list(analyte = "letrozole", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "rat (Sprague-Dawley)",
    n_subjects = 20,
    n_studies = 1,
    age_range = "7-9 weeks (age-matched adults)",
    weight_range = "Females 201-225 g, males 301-325 g at purchase",
    sex_female_pct = 50,
    disease_state = paste(
      "Healthy, jugular-vein-cannulated rats with a striatal microdialysis",
      "probe for brain ECF sampling.",
      sep = " "
    ),
    dose_range = paste(
      "Letrozole 4 mg/kg intraperitoneally, single dose, and once daily for 5",
      "days (males) or 11 days (females) to steady state. The PBPK simulations",
      "used extravascular first-order absorption of the same doses.",
      sep = " "
    ),
    regions = "United States (University of Cincinnati)",
    notes = paste(
      "The PBPK model was bottom-up, not fitted: absorption and clearance were",
      "taken from an earlier oral-gavage rat study (Liu 2000, reference 16),",
      "Vss was predicted by the Rodgers-Rowland method, and the brain",
      "parameters from in situ perfusion and in vivo fu,brain measurements. The",
      "observed study (N = 10 per sex; 3-6 rats per group with data) was used",
      "only to evaluate the predictions (Tables 4 and 5). Virtual rats came",
      "from the Simcyp rodent library.",
      sep = " "
    )
  )

  ini({
    # Absorption. Table 1 'Ka (1/h) 0.29 (males); 0.49 (females)' (from
    # reference 16, first-order absorption model).
    lka <- fixed(log(0.29)); label("First-order absorption rate, males (1/h)") # Table 1 Ka 0.29 (males)
    e_sexf_ka <- fixed(log(0.49) - log(0.29)); label("Log-scale effect of female sex on ka (unitless)") # Table 1 Ka 0.49 (females) vs 0.29 (males)

    # Fraction absorbed. Table 1 'Fa 0.99'. Gut and hepatic first-pass
    # fractions are Simcyp-internal and unprinted; the reduction takes
    # F = Fa (see the vignette assumptions).
    lfdepot <- fixed(log(0.99)); label("Bioavailability, taken as the fraction absorbed Fa (fraction)") # Table 1 Fa 0.99

    # In vivo oral clearance CLpo, entered in Simcyp as the clearance input.
    # Table 1 'CLpo (mL/min) 0.77 (males); 0.21 (females)', footnote a: based
    # on CL = 0.185 (males) and 0.051 (females) L/h/kg assuming a 250 g rat.
    # Stored per kg: 0.77 mL/min * 0.06 / 0.25 kg = 0.1848 L/h/kg.
    lcl <- fixed(log(0.77 * 0.06 / 0.25)); label("Apparent oral clearance CL/F per kg body weight, males (L/h/kg)") # Table 1 CLpo 0.77 mL/min at 250 g (footnote a)
    e_sexf_cl <- fixed(log(0.21) - log(0.77)); label("Log-scale effect of female sex on CL/F (unitless)") # Table 1 CLpo 0.21 (females) vs 0.77 (males) mL/min

    # Volume. Table 1 'Vss (L/kg) 3.29', predicted by Rodgers and Rowland
    # (Method 2), Kp scalar 1. Used as the single volume of the reduction.
    lvc <- fixed(log(3.29)); label("Volume of distribution Vss per kg body weight (L/kg)") # Table 1 Vss 3.29 L/kg

    # Brain. Table 1 multi-compartment brain model.
    # PSB 0.84 mL/min = 0.469 mL/min/g (Eq 3: Fpf 0.036 mL/s/g, Kin 0.422
    # mL/min/g) times an assumed 1800 mg brain weight (footnote b).
    lps_bbb <- fixed(log(0.84 * 0.06)); label("Passive permeability-surface area product across the blood-brain barrier PSB (L/h)") # Table 1 PSB 0.84 mL/min
    lvbrain <- fixed(log(1.8 / 1000)); label("Brain volume (L)") # Table 1 footnote b, 1800 mg brain weight
    fu_plasma <- fixed(0.4); label("Unbound fraction in plasma (fraction)") # Table 1 fraction unbound in plasma 0.4
    fu_brain <- fixed(0.58); label("Unbound fraction in brain (fraction)") # Table 1 fraction unbound in brain 0.58

    # Residual error. The paper reports no residual error model (the PBPK
    # predictions are compared with observed means); fixed at zero.
    propSd <- fixed(0); label("Proportional residual error, plasma (fraction; not reported)")
    addSd <- fixed(0); label("Additive residual error, plasma (ng/mL; not reported)")
    propSd_Cecf <- fixed(0); label("Proportional residual error, brain ECF (fraction; not reported)")
    addSd_Cecf <- fixed(0); label("Additive residual error, brain ECF (ng/mL; not reported)")
  })

  model({
    # Sex-specific absorption and clearance (Table 1).
    ka <- exp(lka + e_sexf_ka * SEXF)
    clpo <- exp(lcl + e_sexf_cl * SEXF) * WT
    vc <- exp(lvc) * WT
    fdepot <- exp(lfdepot)

    # CLpo is the in vivo oral clearance (CL/F); systemic clearance is
    # F * CLpo, so AUC = F * Dose / (F * CLpo) = Dose / CLpo.
    kel <- clpo * fdepot / vc

    ps_bbb <- exp(lps_bbb)
    vbrain <- exp(lvbrain)

    # Unbound plasma and brain concentrations (mg/L); passive exchange
    # across the blood-brain barrier is driven by their difference.
    Cp_u <- fu_plasma * central / vc
    Cbrain_u <- fu_brain * brain / vbrain

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - ps_bbb * (Cp_u - Cbrain_u)
    d/dt(brain) <- ps_bbb * (Cp_u - Cbrain_u)

    f(depot) <- fdepot

    # Plasma (total) and brain ECF (= fu,brain x total brain) in ng/mL.
    Cc <- central / vc * 1000
    Cecf <- Cbrain_u * 1000

    Cc ~ add(addSd) + prop(propSd)
    Cecf ~ add(addSd_Cecf) + prop(propSd_Cecf)
  })
}
