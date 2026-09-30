Jeong_2021_fexuprazan_pbpk <- function() {
  description <- paste0(
    "PBPK (whole-body, 13 compartments, coded in Berkeley Madonna). ",
    "Fexuprazan (DWP14012, a potassium-competitive acid blocker) after oral ",
    "dosing in healthy adults (Jeong et al. 2021, Pharmaceutics). Arterial ",
    "and venous blood pools, lung, and ten perfusion-limited tissues ",
    "(adipose, adrenal gland, brain, heart, kidney, liver, spleen, stomach, ",
    "small and large intestine). Stomach, spleen and both intestines drain ",
    "into the liver; all other tissues and a residual arterio-venous flow ",
    "drain into venous blood. First-order absorption of the fraction ",
    "absorbed (Fa) directly into the liver via the portal vein. Hepatic ",
    "elimination only: two CYP3A4 Michaelis-Menten pathways (formation of ",
    "M14 and M11) driven by the unbound liver concentration, plus a linear ",
    "additional unbound intrinsic clearance. Human tissue partition ",
    "coefficients are rat Kp,SS values scaled by a single Kp scalar to the ",
    "allometric human Vss, with the liver Kp further corrected for hepatic ",
    "extraction. Deterministic: the paper reports no inter-individual ",
    "variance and no residual-error model, so the model is for ",
    "typical-value simulation."
  )
  reference <- paste0(
    "Jeong YS, Kim MS, Lee N, Lee A, Chae YJ, Chung SJ, Lee KR. ",
    "Development of Physiologically Based Pharmacokinetic Model for Orally ",
    "Administered Fexuprazan in Humans. Pharmaceutics. 2021;13(6):813. ",
    "doi:10.3390/pharmaceutics13060813"
  )
  vignette <- "Jeong_2021_fexuprazan"
  units <- list(
    time = "min",
    dosing = "mg",
    concentration = "ng/mL"
  )

  # The two intestinal segments are perfused tissue compartments (not gut
  # lumen segments); the registered small/large-intestine names are luminal
  # (a_small_intestine) or membrane-limited sub-compartments (vp_ / is_).
  paper_specific_compartments <- c("small_intestine", "large_intestine")

  # Every state is an AMOUNT of fexuprazan in mg; model() divides by the
  # Table 1 volume (mL) to obtain the concentration in mg/mL.
  compartmentData <- list(
    depot = list(analyte = "fexuprazan", units = "mg", specimen = "administration site", verified = TRUE),
    venous = list(analyte = "fexuprazan", units = "mg", specimen = "whole blood", verified = TRUE),
    lung = list(analyte = "fexuprazan", units = "mg", specimen = "tissue", verified = TRUE),
    arterial = list(analyte = "fexuprazan", units = "mg", specimen = "whole blood", verified = TRUE),
    adipose = list(analyte = "fexuprazan", units = "mg", specimen = "tissue", verified = TRUE),
    adrenal = list(analyte = "fexuprazan", units = "mg", specimen = "tissue", verified = TRUE),
    brain = list(analyte = "fexuprazan", units = "mg", specimen = "tissue", verified = TRUE),
    heart = list(analyte = "fexuprazan", units = "mg", specimen = "tissue", verified = TRUE),
    kidney = list(analyte = "fexuprazan", units = "mg", specimen = "tissue", verified = TRUE),
    stomach = list(analyte = "fexuprazan", units = "mg", specimen = "tissue", verified = TRUE),
    spleen = list(analyte = "fexuprazan", units = "mg", specimen = "tissue", verified = TRUE),
    small_intestine = list(analyte = "fexuprazan", units = "mg", specimen = "tissue", verified = TRUE),
    large_intestine = list(analyte = "fexuprazan", units = "mg", specimen = "tissue", verified = TRUE),
    liver = list(analyte = "fexuprazan", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 24L,
    n_studies = 2L,
    age_range = NULL,
    weight_range = "70 kg reference adult (Table 1 physiology)",
    sex_female_pct = NA_real_,
    race_ethnicity = "Training set: healthy volunteers (NCT02757144, Korean). Validation set 2: Korean, Caucasian and Japanese healthy volunteers (NCT03574415).",
    disease_state = "Healthy volunteers",
    dose_range = "20, 40 and 80 mg once daily orally for 7 days (training / validation set 1); 40 and 80 mg once daily (validation set 2).",
    regions = "Korea (NCT02757144); Korean, Caucasian and Japanese subjects (NCT03574415)",
    notes = paste0(
      "Fa and CLu,add were optimized by fitting the day-1 plasma profiles of ",
      "24 volunteers (8 per 20, 40, 80 mg dose group) of the multiple ",
      "ascending dose study NCT02757144 (Results 3.1); the model was ",
      "validated on day 7 of that study and on the first and eighth doses ",
      "of NCT03574415 (24 Korean, Caucasian and Japanese volunteers per dose). ",
      "Physiological volumes and flows are for a 70 kg adult. The paper ",
      "does not tabulate subject demographics."
    )
  )

  ini({
    # ----------------------------------------------------------------
    # Organ volumes (mL) -- Table 1 'Volume (mL)' column. Physiological
    # constants taken from the literature (refs 15, 24); none estimated.
    # ----------------------------------------------------------------
    v_adipose <- fixed(15000)
    label("Volume of adipose tissue (mL)") # Table 1 row 'Adipose'
    v_adrenal <- fixed(14)
    label("Volume of adrenal gland (mL)") # Table 1 row 'Adrenal gland'
    v_brain <- fixed(1400)
    label("Volume of brain (mL)") # Table 1 row 'Brain'
    v_heart <- fixed(329)
    label("Volume of heart (mL)") # Table 1 row 'Heart'
    v_kidney <- fixed(308)
    label("Volume of kidney (mL)") # Table 1 row 'Kidney'
    v_large_intestine <- fixed(371)
    label("Volume of large intestine (mL)") # Table 1 row 'Large Intestine'
    v_liver <- fixed(1800)
    label("Volume of liver (mL)") # Table 1 row 'Liver'
    v_lung <- fixed(532)
    label("Volume of lung (mL)") # Table 1 row 'Lung'
    v_small_intestine <- fixed(520)
    label("Volume of small intestine (mL)") # Table 1 row 'Small Intestine'
    v_spleen <- fixed(182)
    label("Volume of spleen (mL)") # Table 1 row 'Spleen'
    v_stomach <- fixed(147)
    label("Volume of stomach (mL)") # Table 1 row 'Stomach'
    v_venous <- fixed(3470)
    label("Volume of venous blood (mL)") # Table 1 row 'Venous blood'
    v_arterial <- fixed(1730)
    label("Volume of arterial blood (mL)") # Table 1 row 'Arterial blood'

    # ----------------------------------------------------------------
    # Organ blood flows (mL/min) -- Table 1 'Blood Flow (mL/min)' column.
    # q_lung equals the cardiac output (Table 1 caption: 5200 mL/min for a
    # 70 kg human). q_liver is the TOTAL hepatic flow, including the portal
    # inflow from stomach, spleen and intestines (Equation 10 subtracts
    # those four flows to obtain the hepatic-artery inflow).
    # ----------------------------------------------------------------
    q_adipose <- fixed(270)
    label("Blood flow to adipose tissue (mL/min)") # Table 1 row 'Adipose'
    q_adrenal <- fixed(15.6)
    label("Blood flow to adrenal gland (mL/min)") # Table 1 row 'Adrenal gland'
    q_brain <- fixed(593)
    label("Blood flow to brain (mL/min)") # Table 1 row 'Brain'
    q_heart <- fixed(208)
    label("Blood flow to heart (mL/min)") # Table 1 row 'Heart'
    q_kidney <- fixed(910)
    label("Blood flow to kidney (mL/min)") # Table 1 row 'Kidney'
    q_large_intestine <- fixed(208)
    label("Blood flow to large intestine (mL/min)") # Table 1 row 'Large Intestine'
    q_liver <- fixed(1326)
    label("Total blood flow to liver (mL/min)") # Table 1 row 'Liver'
    q_lung <- fixed(5200)
    label("Blood flow to lung = cardiac output (mL/min)") # Table 1 row 'Lung' and caption
    q_small_intestine <- fixed(520)
    label("Blood flow to small intestine (mL/min)") # Table 1 row 'Small Intestine'
    q_spleen <- fixed(104)
    label("Blood flow to spleen (mL/min)") # Table 1 row 'Spleen'
    q_stomach <- fixed(52)
    label("Blood flow to stomach (mL/min)") # Table 1 row 'Stomach'

    # ----------------------------------------------------------------
    # Human tissue-to-plasma partition coefficients -- Table 3
    # 'Distribution (Kp)' rows, 'Corrected by Kp,scalar'. Each is the rat
    # Kp,SS of Table 2 times Kp,scalar = 0.371 (Section 2.3); the liver value
    # is additionally corrected for hepatic extraction.
    # ----------------------------------------------------------------
    lkp_adrenal <- fixed(log(20.8))
    label("Log adrenal gland-to-plasma partition coefficient (unitless)") # Table 3 Kp 'Adrenal gland' = 20.8
    lkp_adipose <- fixed(log(4.32))
    label("Log adipose-to-plasma partition coefficient (unitless)") # Table 3 Kp 'Adipose' = 4.32
    lkp_brain <- fixed(log(1.32))
    label("Log brain-to-plasma partition coefficient (unitless)") # Table 3 Kp 'Brain' = 1.32
    lkp_heart <- fixed(log(4.60))
    label("Log heart-to-plasma partition coefficient (unitless)") # Table 3 Kp 'Heart' = 4.60
    lkp_kidney <- fixed(log(16.4))
    label("Log kidney-to-plasma partition coefficient (unitless)") # Table 3 Kp 'Kidney' = 16.4
    lkp_liver <- fixed(log(303))
    label("Log liver-to-plasma partition coefficient (unitless)") # Table 3 Kp 'Liver' = 303
    lkp_lung <- fixed(log(87.6))
    label("Log lung-to-plasma partition coefficient (unitless)") # Table 3 Kp 'Lung' = 87.6
    lkp_large_intestine <- fixed(log(40.8))
    label("Log large intestine-to-plasma partition coefficient (unitless)") # Table 3 Kp 'Large Intestine' = 40.8
    lkp_small_intestine <- fixed(log(124))
    label("Log small intestine-to-plasma partition coefficient (unitless)") # Table 3 Kp 'Small Intestine' = 124
    lkp_spleen <- fixed(log(17.8))
    label("Log spleen-to-plasma partition coefficient (unitless)") # Table 3 Kp 'Spleen' = 17.8
    lkp_stomach <- fixed(log(193))
    label("Log stomach-to-plasma partition coefficient (unitless)") # Table 3 Kp 'Stomach' = 193

    # ----------------------------------------------------------------
    # Drug-specific constants -- Table 3.
    # ----------------------------------------------------------------
    fu <- fixed(0.0645)
    label("Fraction of fexuprazan unbound in plasma, fup (unitless)") # Table 3 'fup' = 0.0645 (Determined, ref 16)
    bpr <- fixed(0.8)
    label("Blood-to-plasma concentration ratio R (unitless)") # Table 3 'B/P ratio (R)' = 0.8; Section 2.3 'R = 0.8'
    fumic <- fixed(0.904)
    label("Fraction of fexuprazan unbound in human liver microsomes, fu,mic (unitless)") # Table 3 'fu,mic' = 0.904 (Predicted in silico, ref 21)
    mw <- fixed(410.4)
    label("Fexuprazan molecular weight (g/mol)")
    # NOT printed in Jeong 2021. Chemical constant for the fexuprazan molecular formula C19H17F3N2O3S (410.4 g/mol; the structure is the M11 hydroxylamine of Section 2.4 minus its N-hydroxy oxygen). Needed only as the unit bridge between the mg/mL concentrations of the ODEs and the uM scale of the Table 3 Km values; not a fitted model parameter.

    # ----------------------------------------------------------------
    # Absorption -- Table 3 'Absorption' rows.
    # ----------------------------------------------------------------
    lka <- fixed(log(0.0606))
    label("Log first-order absorption rate constant Ka (1/min)")
    # Table 3 'Ka (min-1)' = 0.0606, 'Predicted' from Ka = 2 * Peff / radius (Section 2.2) with Peff from the Caco-2 correlation of Equation 1; not estimated.
    lfdepot <- log(0.761)
    label("Log fraction absorbed Fa (unitless)")
    # Table 3 'Fa' = 0.761 'Optimized'; Results 3.1: average of the per-subject fits of 24 volunteers (0.627, 0.767, 0.890 at 20, 40, 80 mg).

    # ----------------------------------------------------------------
    # Hepatic elimination -- Table 3 'Elimination' and 'M14 / M11 Formation
    # by CYP3A4' rows. Vmax is the whole-liver value (MPPGL 39.79 mg/g and
    # 1800 g liver, Section 2.4).
    # ----------------------------------------------------------------
    lclint_add <- log(12900)
    label("Log additional unbound intrinsic clearance CLu,add (mL/min)")
    # Table 3 'CLu,add (L/min)' = 12.9 'Optimized', stored in mL/min to match the Table 1 flows; Results 3.1: average of per-subject fits (15.8, 13.9, 8.95 L/min at 20, 40, 80 mg).
    lvmax_m14 <- fixed(log(248))
    label("Log whole-liver Vmax of CYP3A4 M14 formation (nmol/min)") # Table 3 'M14 Formation by CYP3A4' Vmax = 248
    lkm_m14 <- fixed(log(0.093))
    label("Log Km of CYP3A4 M14 formation in liver microsomes (uM)") # Table 3 'M14 Formation by CYP3A4' Km = 0.093
    lvmax_m11 <- fixed(log(800))
    label("Log whole-liver Vmax of CYP3A4 M11 formation (nmol/min)") # Table 3 'M11 Formation by CYP3A4' Vmax = 800
    lkm_m11 <- fixed(log(15.95))
    label("Log Km of CYP3A4 M11 formation in liver microsomes (uM)") # Table 3 'M11 Formation by CYP3A4' Km = 15.95

    # ----------------------------------------------------------------
    # Residual error. The paper fitted Fa and CLu,add per subject in
    # WinNonlin but reports no residual-error model or magnitude, and no
    # inter-individual variances (only the between-subject SD of the
    # per-subject estimates, Results 3.1). Encoded as zero rather than
    # invented; see the vignette Assumptions and deviations.
    # ----------------------------------------------------------------
    propSd <- fixed(0)
    label("Proportional residual error (fraction; ZERO - not reported in source)")
  })

  model({
    # Absorption and bioavailability (Section 2.5): the depot holds Fa * dose
    # and empties first-order into the liver (Equation 6 and the Ka * Xa
    # term of Equation 10).
    ka <- exp(lka)
    fdepot <- exp(lfdepot)
    clint_add <- exp(lclint_add)
    vmax_m14 <- exp(lvmax_m14)
    km_m14 <- exp(lkm_m14)
    vmax_m11 <- exp(lvmax_m11)
    km_m11 <- exp(lkm_m11)

    kp_adrenal <- exp(lkp_adrenal)
    kp_adipose <- exp(lkp_adipose)
    kp_brain <- exp(lkp_brain)
    kp_heart <- exp(lkp_heart)
    kp_kidney <- exp(lkp_kidney)
    kp_liver <- exp(lkp_liver)
    kp_lung <- exp(lkp_lung)
    kp_large_intestine <- exp(lkp_large_intestine)
    kp_small_intestine <- exp(lkp_small_intestine)
    kp_spleen <- exp(lkp_spleen)
    kp_stomach <- exp(lkp_stomach)

    # Cardiac output QCO is the lung flow (Table 1 caption: 5200 mL/min).
    qco <- q_lung
    # Hepatic-artery inflow, Equation 10: (QLI - QST - QSP - QSm,IN - QLa,IN).
    q_hepatic_artery <- q_liver - q_stomach - q_spleen - q_small_intestine -
      q_large_intestine
    # Residual flow QRE of Equation 11: the part of the cardiac output that
    # returns to venous blood without passing a modelled tissue, so that
    # venous inflow balances QCO.
    q_residual <- qco - q_adrenal - q_adipose - q_brain - q_heart - q_kidney -
      q_liver

    # Concentrations (mg/mL) = amount (mg) / volume (mL). Blood pools hold
    # blood concentrations (Cart, Cven); tissues hold tissue concentrations.
    c_venous <- venous / v_venous
    c_arterial <- arterial / v_arterial
    c_lung <- lung / v_lung
    c_adipose <- adipose / v_adipose
    c_adrenal <- adrenal / v_adrenal
    c_brain <- brain / v_brain
    c_heart <- heart / v_heart
    c_kidney <- kidney / v_kidney
    c_stomach <- stomach / v_stomach
    c_spleen <- spleen / v_spleen
    c_small_intestine <- small_intestine / v_small_intestine
    c_large_intestine <- large_intestine / v_large_intestine
    c_liver <- liver / v_liver

    # Hepatic intrinsic clearance, Equation 5:
    #   CLu,int = sum_j Vmax_j / (Km_j * fu,mic + CLI * fu,LI) + CLu,add
    # with fu,LI = fup / Kp,LI (Section 2.4), so CLI * fu,LI is the unbound
    # liver concentration. That concentration is converted to uM for the
    # Km terms: mg/mL * 1e6 = ng/mL, and ng/mL / MW (g/mol) = nmol/mL = uM.
    # Vmax (nmol/min) / uM (nmol/mL) gives mL/min, the unit of CLu,add.
    cu_liver <- c_liver * fu / kp_liver
    cu_liver_um <- cu_liver * 1e6 / mw
    clint_u <- vmax_m14 / (km_m14 * fumic + cu_liver_um) +
      vmax_m11 / (km_m11 * fumic + cu_liver_um) +
      clint_add

    d/dt(depot) <- -ka * depot
    f(depot) <- fdepot

    # Perfusion-limited non-eliminating tissues, Equation 9:
    #   VT dCT/dt = QT * (Cart - CT * R / Kp)
    d/dt(adipose) <- q_adipose * (c_arterial - c_adipose * bpr / kp_adipose)
    d/dt(adrenal) <- q_adrenal * (c_arterial - c_adrenal * bpr / kp_adrenal)
    d/dt(brain) <- q_brain * (c_arterial - c_brain * bpr / kp_brain)
    d/dt(heart) <- q_heart * (c_arterial - c_heart * bpr / kp_heart)
    d/dt(kidney) <- q_kidney * (c_arterial - c_kidney * bpr / kp_kidney)
    d/dt(stomach) <- q_stomach * (c_arterial - c_stomach * bpr / kp_stomach)
    d/dt(spleen) <- q_spleen * (c_arterial - c_spleen * bpr / kp_spleen)
    d/dt(small_intestine) <- q_small_intestine *
      (c_arterial - c_small_intestine * bpr / kp_small_intestine)
    d/dt(large_intestine) <- q_large_intestine *
      (c_arterial - c_large_intestine * bpr / kp_large_intestine)

    # Liver, Equation 10: absorbed drug, hepatic-artery inflow and the
    # portal outflows of stomach, spleen and intestines in; total hepatic
    # outflow and unbound intrinsic clearance out.
    d/dt(liver) <- ka * depot +
      q_hepatic_artery * c_arterial +
      q_stomach * c_stomach * bpr / kp_stomach +
      q_spleen * c_spleen * bpr / kp_spleen +
      q_small_intestine * c_small_intestine * bpr / kp_small_intestine +
      q_large_intestine * c_large_intestine * bpr / kp_large_intestine -
      q_liver * c_liver * bpr / kp_liver -
      clint_u * cu_liver

    # Venous blood, Equation 11.
    d/dt(venous) <- q_adrenal * c_adrenal * bpr / kp_adrenal +
      q_adipose * c_adipose * bpr / kp_adipose +
      q_brain * c_brain * bpr / kp_brain +
      q_heart * c_heart * bpr / kp_heart +
      q_kidney * c_kidney * bpr / kp_kidney +
      q_liver * c_liver * bpr / kp_liver +
      q_residual * c_arterial -
      qco * c_venous

    # Lung, Equation 12.
    d/dt(lung) <- qco * (c_venous - c_lung * bpr / kp_lung)

    # Arterial blood, Equation 13.
    d/dt(arterial) <- qco * (c_lung * bpr / kp_lung - c_arterial)

    # Observation: the paper reports PLASMA concentrations (ng/mL). The
    # blood pools hold blood concentrations, so plasma = venous blood / R;
    # mg/mL * 1e6 = ng/mL.
    Cc <- c_venous / bpr * 1e6
    Cc ~ prop(propSd)
  })
}
