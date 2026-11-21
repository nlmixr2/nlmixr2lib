Sakares_2025_deferasirox <- function() {
  description <- paste(
    "PBPK/PD (iron-chelation PD layer only). Two-compartment iron-homeostasis",
    "model (chelatable iron in the blood volume, storage iron in the liver)",
    "with deferasirox (DFX) binding chelatable iron at 2:1 stoichiometry and",
    "the DFX-iron complex eliminated first-order, in adults with",
    "transfusion-dependent thalassaemia (TDT) and healthy adults (Sakares",
    "2025). Uptake from chelatable to storage iron is Michaelis-Menten plus a",
    "linear chelatable-iron-driven term; efflux from storage iron is",
    "Michaelis-Menten scaled by a hepcidin-like factor that rises with liver",
    "iron concentration (LIC) and falls with chelatable iron. TDT lowers",
    "dietary iron absorption and adds a transfusional iron input. The paper's",
    "DFX PK layer is a Simcyp V22 whole-body PBPK model (ADAM absorption,",
    "retrograde enzyme kinetics, platform physiology) whose organ volumes,",
    "blood flows and partition coefficients are not printed, so it is not",
    "reproduced: total plasma DFX is supplied as the time-varying covariate",
    "CP_DEFERASIROX_UM. The output liver_iron is liver iron concentration (mg Fe/g",
    "dry weight). No inter-individual variability or residual error is",
    "reported; the model is deterministic.",
    sep = " "
  )
  reference <- paste(
    "Sakares W, Udomnilobol U, Pan X, Yampayon K, Ekwattanakit S, Viprakasit V,",
    "Limwongse C, Srichairatanakool S, Teerawonganan P, Kunsa-ngiem S,",
    "Khaowroongrueng V, Techatanawat I, Prueksaritanont T, Yodsurang V.",
    "Bridging Pharmacokinetics and Pharmacodynamics: A PBPK/PD Model-Based",
    "Approach for Deferasirox Dosing in Transfusion-Dependent Thalassemia.",
    "Clin Transl Sci. 2025;18(9):e70355. doi:10.1111/cts.70355",
    "(PMID 40965322; PMCID PMC12445128).",
    sep = " "
  )
  vignette <- "Sakares_2025_deferasirox"

  units <- list(
    time = "h",
    dosing = "(not applicable; deferasirox enters as the CP_DEFERASIROX_UM covariate, not as a dose record)",
    concentration = "umol/L (iron and DFX-iron complex states; the output liver_iron is in mg/g dry weight)"
  )

  covariateData <- list(
    CP_DEFERASIROX_UM = list(
      description = "Total (bound + unbound) plasma deferasirox concentration driving iron chelation",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying. Sakares 2025 Eq. 3 and Eq. 6 use [DFXsys], the systemic",
        "DFX concentration after the iron interaction, multiplied by the unbound",
        "fraction fu inside the binding term, so the column carries TOTAL plasma",
        "DFX. In the source this profile is produced by the authors' Simcyp V22",
        "PBPK model (Table 1: MW 373.4 g/mol, fu 0.004, B:P 0.61, Vss 0.181",
        "L/kg, CLiv 3.53 L/h, 65% enterohepatic reabsorption of biliary DFX),",
        "which is not reproducible from the publication. Supply observed or",
        "independently simulated DFX concentrations; divide ng/mL by 373.4 to",
        "get umol/L. Set to 0 for untreated periods."
      ),
      source_name = "[DFXsys]"
    ),
    TRANSF_IRON_MGKGD = list(
      description = "Daily transfusional iron loading (mg iron per kg body weight per day)",
      units = "mg/kg/day",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed per subject. Converted to the zero-order transfusional iron",
        "input kblood (umol/L/h) into the chelatable-iron compartment as",
        "TRANSF_IRON_MGKGD * WT / 55.845 * 1000 / 24 / vol_cfe (iron atomic",
        "mass 55.845 g/mol). Table S1 gives 0.31-0.62 mg/kg/day for the TDT",
        "datasets and 0 for healthy subjects; Table S6 maps transfusion",
        "regimens to 0.16 (1 unit every 4 weeks), 0.22 (1 unit/3 wk), 0.32",
        "(2 units/4 wk), 0.45 (2 units/3 wk) and 0.67 (3 units/3 wk) mg/kg/day."
      ),
      source_name = "k_blood (Table 1); 'Iron loading from blood transfusion (mg/kg/day)' (Table S1)"
    ),
    LIVER_IRON_BL = list(
      description = "Baseline liver iron concentration (mg iron per g liver dry weight)",
      units = "mg/g",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed per subject; sets the storage-iron initial condition",
        "sfe(0) = LIVER_IRON_BL * 2886.2 umol/L. The constant 2886.2 is printed in",
        "Eq. 1-4 and equals the healthy initial storage iron 2828.48 umol/L",
        "(Table 1) divided by the healthy LIC 0.98 mg/g (Table S1 footnote).",
        "Healthy adults: 0.98 mg/g. TDT datasets: 2.5-31.44 mg/g (Table S1)."
      ),
      source_name = "Baseline LIC (Table S1); 'Initial storage iron concentration' (Table 1)"
    ),
    DIS_THALASSEMIA_TD = list(
      description = "Transfusion-dependent thalassaemia indicator (1 = TDT, 0 = healthy)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy adult)",
      notes = paste(
        "Selects the iron-absorption rate kin: 0.28 umol/L/h in healthy adults",
        "and 0.07 umol/L/h in TDT (Table 1, PD parameters, Healthy and TDT",
        "columns). The other population difference in Table 1, the",
        "transfusional input kblood, enters through TRANSF_IRON_MGKGD."
      ),
      source_name = "Population column (Healthy / TDT) of Table 1"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Used only to convert the per-kg transfusional iron loading into an",
        "absolute iron input. The paper does not print the weight it used.",
        "The Table S6 regimen-to-loading mapping implies about 43-45 kg at",
        "200 mg iron per red-cell unit, consistent with the Thai TDT",
        "population (Table S2 weight-height equations)."
      ),
      source_name = "Weight (Simcyp population model, Table S2)"
    )
  )

  # cfe and sfe are the paper's own chelatable-iron and storage-iron pools.
  # The registered `iron` state is reserved for serum iron as a PD biomarker.
  paper_specific_compartments <- c("cfe", "sfe")

  compartmentData <- list(
    cfe = list(
      analyte = "iron (chelatable: non-transferrin-bound and labile iron)",
      units = "umol/L",
      specimen = "whole blood",
      verified = TRUE
    ),
    sfe = list(
      analyte = "iron (storage: ferritin and haemosiderin)",
      units = "umol/L",
      specimen = "tissue",
      verified = TRUE
    ),
    complex = list(analyte = "deferasirox-iron complex", units = "umol/L", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = "not applicable (literature-calibrated PD layer; Simcyp virtual populations of 10 trials per scenario)",
    n_studies = "21 datasets from 12 sources (Table S1): 8 healthy Caucasian, 9 TDT Caucasian, 1 healthy Thai, 3 TDT Thai",
    age_range = "18-66 years across datasets (Table S1)",
    sex_female_pct = "0-75% across datasets (Table S1)",
    disease_state = "Transfusion-dependent thalassaemia (TDT) and healthy adults",
    dose_range = "Deferasirox 5-40 mg/kg/day oral (TDT, multiple dose); 65-130 mg IV and 250 mg-35 mg/kg oral single dose",
    regions = "Caucasian (Europe/North America literature) and Thailand (Siriraj Thalassemia Center; Government Pharmaceutical Organization)",
    notes = paste(
      "The PD layer was calibrated to LIC progression rates and daily",
      "transfusional iron in healthy and TDT subjects (COPASI 4.44), and",
      "verified against LIC 6 or 12 months after DFX in four Caucasian TDT",
      "datasets (Piga 2006, Cappellini 2006) and two Thai TDT datasets",
      "(Table S5). Baseline LIC 2.5-31.44 mg/g and transfusional iron",
      "0.31-0.62 mg/kg/day across the TDT datasets (Table S1)."
    )
  )

  ini({
    # Iron turnover (Table 1, PD parameters). kin and kout are zero-order
    # rates into / out of the chelatable-iron compartment (umol/L/h).
    kin_fe <- fixed(0.28)
    label("Iron absorption rate into chelatable iron, healthy (umol/L/h)") # Table 1 'k in' Healthy column = 0.28 (Ganz)
    kin_fe_tdt <- fixed(0.07)
    label("Iron absorption rate into chelatable iron, TDT (umol/L/h)") # Table 1 'k in' TDT column = 0.07 (Porter)
    kout_fe <- fixed(0.28)
    label("Iron loss rate from chelatable iron (umol/L/h)") # Table 1 'k out' = 0.28 in both columns (Ganz)

    # Compartment volumes (Table 1)
    vol_sfe <- fixed(1.69)
    label("Storage-iron compartment volume, liver (L)") # Table 1 'Storage iron (sFe) compartment volume' = 1.69 (Peters)
    vol_cfe <- fixed(5.40)
    label("Chelatable-iron compartment volume, total blood (L)") # Table 1 'Chelatable iron (cFe) compartment volume' = 5.40 (Peters)

    # Uptake transporter, chelatable -> storage iron (Table 1)
    vmax_uptake <- 0.306687
    label("Maximum uptake velocity, chelatable to storage iron (umol/L/h)") # Table 1 'V max,uptake' = 0.306687 (estimated from Parkes)
    km_uptake <- fixed(15.2)
    label("Michaelis constant of the uptake transporter (umol/L)") # Table 1 'Kmuptake' = 15.2 (Parkes)
    k_cfe_uptake <- 1e-6
    label("Linear chelatable-iron-driven uptake rate constant (1/h)") # Table 1 'k a' = 1e-6 (estimated from Parkes)

    # Efflux transporter, storage -> chelatable iron (Table 1)
    vmax_efflux <- 43
    label("Maximum efflux velocity, storage to chelatable iron (umol/L/h)") # Table 1 'V max,efflux' = 43 (estimated)
    km_efflux <- 1e6
    label("Michaelis constant of the efflux transporter (umol/L)") # Table 1 'Kmefflux' = 1e6 (estimated)

    # Hepcidin-like modulation of efflux (Table 1; Eq. 1-4)
    hep_a <- fixed(1.66)
    label("Effect of liver iron concentration on hepcidin (g/mg)") # Table 1 'a' = 1.66 (Aoki)
    hep_b <- fixed(0.036)
    label("Effect of chelatable iron on hepcidin (L/umol)") # Table 1 'b' = 0.036 (Gehrke)
    hep_c <- fixed(0.135)
    label("Exponent of the hepcidin factor on efflux (unitless)") # Table 1 'c' = 0.135 (Chung)
    sfe_per_lic <- fixed(2886.2)
    label("Storage iron per unit liver iron concentration (umol/L per mg/g)") # Eq. 1-4 constant 2886.2 = 2828.48 / 0.98

    # Deferasirox-iron binding (Table 1; Eq. 3, 6)
    kon <- 200
    label("DFX-iron complex formation rate constant (L^2/umol^2/h)") # Table 1 'k on' = 200 (manually optimised)
    koff <- fixed(0.96)
    label("DFX-iron complex dissociation rate constant (1/h)") # Table 1 'k off' = 0.96 (Suzuki and Iki)
    keli <- 0.4
    label("DFX-iron complex elimination rate constant (1/h)") # Table 1 'k eli' = 0.4 (manually optimised)
    fu <- fixed(0.004)
    label("Unbound fraction of deferasirox in plasma (fraction)") # Table 1 'Fraction of unbound in plasma' = 0.004 (Weiss)

    # Chelatable-iron initial condition: Table 1 gives it only as depending
    # on initial storage iron. 2.2853 umol/L is the healthy steady state of
    # Eq. 1-2 at the Table 1 values with sfe = 2828.48 umol/L.
    cfe0 <- fixed(2.2853)
    label("Initial chelatable iron concentration (umol/L)") # derived: steady state of Eq. 1-2, Table 1 healthy inputs

    addSd <- fixed(0)
    label("Additive residual error on LIC (mg/g; zero, not reported)") # not reported
  })

  model({
    # Population switch on dietary iron absorption (Table 1 Healthy / TDT)
    kin <- kin_fe * (1 - DIS_THALASSEMIA_TD) + kin_fe_tdt * DIS_THALASSEMIA_TD

    # Transfusional iron input kblood (umol/L/h) from mg/kg/day; 55.845 is
    # the atomic mass of iron (mg/mmol).
    kblood <- TRANSF_IRON_MGKGD * WT / 55.845 * 1000 / 24 / vol_cfe

    # Transporter fluxes (Eq. 1-4), each in the concentration units of its
    # own source compartment.
    uptake <- vmax_uptake * cfe / (km_uptake + cfe) + k_cfe_uptake * cfe
    hepcidin <- (hep_a * sfe / sfe_per_lic * exp(-hep_b * cfe))^hep_c
    efflux <- vmax_efflux * sfe / (km_efflux + sfe) * hepcidin

    # Net DFX-iron complex formation (Eq. 3 and 6); total plasma DFX times
    # fu is the unbound DFX, two molecules per iron.
    binding <- kon * cfe * (CP_DEFERASIROX_UM * fu)^2 - koff * complex

    d/dt(cfe) <- kin + kblood - kout_fe - uptake + efflux * vol_sfe / vol_cfe - binding
    d/dt(sfe) <- uptake * vol_cfe / vol_sfe - efflux
    d/dt(complex) <- binding - keli * complex

    cfe(0) <- cfe0
    sfe(0) <- LIVER_IRON_BL * sfe_per_lic

    liver_iron <- sfe / sfe_per_lic
    liver_iron ~ add(addSd)
  })
}
