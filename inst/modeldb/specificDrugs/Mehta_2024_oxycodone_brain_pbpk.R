Mehta_2024_oxycodone_brain_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, permeability-limited 4-compartment brain, Pumas 2.2.0).",
    "Oxycodone in adults after IV dosing (Mehta et al. 2024, the human",
    "oxycodone instantiation of an interspecies brain PBPK platform).",
    "Thirteen perfusion-limited body compartments (venous and arterial",
    "blood, lung, adipose, bone, heart, kidney, muscle, skin, spleen, gut,",
    "liver, rest of body) with Poulin and Theil tissue partition",
    "coefficients computed from adult tissue composition, plus the",
    "Verscheijden 2019 brain model (brain blood, brain mass, cranial CSF,",
    "spinal CSF) with passive blood-brain and blood-CSF barrier",
    "permeability. The blood-brain barrier permeability was optimised on",
    "rat neuropharmacokinetic data and scaled to humans by BBB surface",
    "area. Organ volumes and blood flows are age-, sex-, weight- and",
    "height-dependent adult equations. Total plasma clearance acts on",
    "venous blood. Forward-simulation model: all parameters are fixed;",
    "IIV is on clearance and the CSF flows; no residual error was",
    "reported. With the published Poulin and Theil inputs the terminal",
    "half-life is about 1.3 h, shorter than the about 3 h of the paper's",
    "Figure 5C (see the vignette)."
  )
  reference <- paste(
    "Mehta P, Soliman A, Rodriguez-Vera L, Schmidt S, Muniz P, Rodriguez M,",
    "Forcadell M, Gonzalez-Perez E, Vozmediano V. Interspecies brain PBPK",
    "modeling platform to predict passive transport through the blood-brain",
    "barrier and assess target site disposition. Pharmaceutics.",
    "2024;16(2):226. doi:10.3390/pharmaceutics16020226."
  )
  vignette <- "Mehta_2024_brain_pbpk"

  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "mg/L"
  )
  dosing <- "venous"

  compartmentData <- list(
    venous = list(analyte = "oxycodone", units = "mg", specimen = "whole blood", verified = TRUE),
    arterial = list(analyte = "oxycodone", units = "mg", specimen = "whole blood", verified = TRUE),
    lung = list(analyte = "oxycodone", units = "mg", specimen = "tissue", verified = TRUE),
    adipose = list(analyte = "oxycodone", units = "mg", specimen = "tissue", verified = TRUE),
    bone = list(analyte = "oxycodone", units = "mg", specimen = "tissue", verified = TRUE),
    heart = list(analyte = "oxycodone", units = "mg", specimen = "tissue", verified = TRUE),
    kidney = list(analyte = "oxycodone", units = "mg", specimen = "tissue", verified = TRUE),
    muscle = list(analyte = "oxycodone", units = "mg", specimen = "tissue", verified = TRUE),
    skin = list(analyte = "oxycodone", units = "mg", specimen = "tissue", verified = TRUE),
    spleen = list(analyte = "oxycodone", units = "mg", specimen = "tissue", verified = TRUE),
    other = list(analyte = "oxycodone", units = "mg", specimen = "tissue", verified = TRUE),
    gut = list(analyte = "oxycodone", units = "mg", specimen = "tissue", verified = TRUE),
    liver = list(analyte = "oxycodone", units = "mg", specimen = "tissue", verified = TRUE),
    brain_vascular = list(analyte = "oxycodone", units = "mg", specimen = "whole blood", verified = TRUE),
    brain = list(analyte = "oxycodone", units = "mg", specimen = "tissue", verified = TRUE),
    brain_csf_sas_cranial = list(analyte = "oxycodone", units = "mg", specimen = "CSF", verified = TRUE),
    brain_csf_sas_spinal = list(analyte = "oxycodone", units = "mg", specimen = "CSF", verified = TRUE)
  )

  covariateData <- list(
    AGE = list(
      description = paste(
        "Age. Enters cardiac output and the adult muscle-volume equation."
      ),
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Adult equations only (the platform's human physiology is the",
        "adult branch of Verscheijden 2019). The oxycodone verification",
        "study (Kokki et al.) enrolled 26-60 y."
      ),
      source_name = "age"
    ),
    WT = list(
      description = paste(
        "Total body weight. Drives the brain, lung, adipose, bone, kidney,",
        "muscle, spleen, gut, blood and rest-of-body volumes and, through",
        "body surface area, cardiac output."
      ),
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The supplementary code draws weight from a sex-specific",
        "height-to-weight equation with log-normal variability (Table S2);",
        "supply the individual weight here."
      ),
      source_name = "Weight"
    ),
    HT = list(
      description = paste(
        "Body height. Drives body surface area (Du Bois) and hence cardiac",
        "output, heart, skin and liver volume, and enters the lung, adipose,",
        "kidney, muscle and blood volume equations."
      ),
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = "Body surface area is computed internally as 0.007184 * HT^0.725 * WT^0.425.",
      source_name = "Height"
    ),
    SEXF = list(
      description = paste(
        "Female sex indicator. Selects the sex-specific organ volume, CSF",
        "volume and blood flow equations."
      ),
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = male",
      notes = paste(
        "1 = female, 0 = male. The supplementary code builds male and",
        "female virtual subjects from separate equation sets."
      ),
      source_name = "male / female virtual subjects"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 11L,
    n_studies = 1L,
    age_range = "26-60 years (Kokki et al., Table 1)",
    weight_range = "not reported; generated from age and height in the virtual population",
    sex_female_pct = NA_real_,
    disease_state = paste(
      "Adults with plasma and CSF sampling after a single IV oxycodone",
      "dose (Kokki et al., used for verification only)."
    ),
    dose_range = "0.092 mg/kg IV, single dose",
    regions = NA_character_,
    notes = paste(
      "Forward-simulation PBPK: no human parameter was estimated. Clearance",
      "and binding come from the literature (Table S3); the blood-brain",
      "barrier permeability was optimised on rat blood and brain data",
      "(Ball et al., 0.3 mg/kg 1-h IV infusion, n = 11) and scaled to humans",
      "by the BBB surface area (Table 3). Human data were used only for",
      "verification (Table 4, Fig 5A-B). The virtual-population size and",
      "sex ratio were not reported (placeholders in the supplementary code)."
    )
  )

  ini({
    # ---- Clearance (Table S3; Methods 2.4). Total plasma clearance, applied
    # in the code as tvCL / BP to venous whole blood. Table S3 lists a
    # 'CL-scalar' of 0.8 for oxycodone, optimised on the human data; it is
    # applied to the literature CL of 49.2 L/h (0.8 x 49.2 = 39.36 L/h).
    lcl <- fixed(log(39.36)); label("Total plasma clearance (L/h)") # Table S3, 'CL - Human (L/hr)' oxycodone = 49.2 x 'CL-scalar' 0.8

    # ---- Physicochemical and binding inputs (Table S3).
    logp <- fixed(1.2); label("Oxycodone log10 octanol:water partition coefficient (unitless)") # Table S3, Log P oxycodone = 1.2
    pka <- fixed(8.53); label("Oxycodone basic pKa (unitless)") # Table S3, Pka oxycodone = 8.53 (monoprotic base)
    bpr <- fixed(1.3); label("Blood:plasma concentration ratio (unitless)") # Table S3, BP Humans oxycodone = 1.3
    fu <- fixed(0.59); label("Fraction unbound in plasma (fraction)") # Table S3, Fupl Humans oxycodone = 0.59
    fu_br <- fixed(0.45); label("Fraction unbound in brain mass (fraction)") # Table S3, Fubm Humans oxycodone = 0.45
    fu_csf <- fixed(1); label("Fraction unbound in cranial CSF (fraction)") # Table S3, Fuccsf oxycodone = 1 (assumed)

    # ---- Brain barrier permeability-surface-area products (Table 3, Table S3).
    lpsb <- fixed(log(14.47)); label("Blood-brain barrier PS product, human (L/h)") # Table 3, PSB in Humans oxycodone = 14.47
    f_psc <- fixed(0.1); label("Blood-CSF barrier PS product as a fraction of the BBB PS product (fraction)") # Table S3 footnote beta 'Assumed to be 1/10th of PSB'; Methods 2.3.2
    lpse <- fixed(log(300)); label("Brain mass to cranial CSF PS product (L/h)") # Table S3, PSE (Brain-CCSF) Humans = 300, 'Assumed to be no barrier'

    # ---- CSF flows (Table S2).
    lq_csf_prod <- fixed(log(0.021)); label("CSF production rate (L/h)") # Table S2, Qproductionrate = 0.021, CV 10%
    lf_bulk <- fixed(log(0.25)); label("Bulk flow brain mass to cranial CSF as a fraction of CSF production (fraction)") # Table S2, Qbulk = 0.25*Qproductionrate, CV 8%
    lf_ssink <- fixed(log(0.38)); label("Spinal CSF to blood flow as a fraction of total CSF outflow (fraction)") # Table S2, Qssink = 0.38*(0.75*Qproductionrate + Qbulk), CV 30%
    lf_sout <- fixed(log(0.9)); label("Spinal to cranial CSF flow as a fraction of the spinal sink flow (fraction)") # Table S2, Qsout = 0.9*Qssink, CV 100%

    # ---- Between-subject variability. Clearance: log-normal CV 30%
    # (Table S3), variance log(1 + 0.30^2). CSF flows: the supplementary
    # code multiplies each flow by exp(rand(Normal(0, SD))), so variance =
    # SD^2.
    etalcl ~ fixed(log(1 + 0.30^2)) # Table S3, CL - Human oxycodone CV (%) 30
    etalq_csf_prod ~ fixed(0.1^2) # Table S2, Qproductionrate CV 10%; code rand(Normal(0,0.1))
    etalf_bulk ~ fixed(0.08^2) # Table S2, Qbulk CV 8%; code rand(Normal(0,0.08))
    etalf_ssink ~ fixed(0.3^2) # Table S2, Qssink CV 30%; code rand(Normal(0,0.3))
    etalf_sout ~ fixed(1^2) # Table S2, Qsout CV 100%; code rand(Normal(0,1))

    # ---- Residual error: not reported (forward simulation only).
    propSd <- fixed(0); label("Proportional residual error on plasma oxycodone (fraction; not reported)") # not reported: forward-simulation model
  })

  model({
    # =====================================================================
    # Adult physiology (Table S2; supplementary model code, human section).
    # =====================================================================
    fem <- 0
    if (SEXF > 0.5) fem <- 1
    bsa <- 0.007184 * HT^0.725 * WT^0.425
    htm <- HT / 100

    # Brain (cranial CSF fraction 0.105 male, 0.092 female)
    v_brain <- (1.449 - 3.62 / WT) / 1.04
    v_bb <- 0.05 * v_brain
    csf_frac <- 0.105 * (1 - fem) + 0.092 * fem
    v_ccsf <- v_brain * csf_frac * 0.8
    v_scsf <- v_brain * csf_frac * 0.2
    v_endo <- 0.005 * v_brain
    v_bm <- v_brain - v_endo - v_bb - v_ccsf - v_scsf
    # Brain extracellular fluid volume, used only for the ECF output
    v_ecf <- 0.29 # code: Vecf = 0.29 L

    v_lung <- (((29.08 * htm * WT^0.5 + 11.06 + 35.47 * htm * WT^0.5 + 5.53) / 1000) / 1.05) * (1 - fem) +
      (((31.46 * htm * WT^0.5 + 1.43 + 35.3 * htm * WT^0.5 + 1.53) / 1000) / 1.05) * fem
    v_adipose <- ((1.36 * WT) / htm - 42) * (1 - fem) + ((1.61 * WT) / htm - 38.3) * fem
    if (v_adipose < 0) v_adipose <- 0.01
    v_bone <- ((WT - v_adipose * 0.92) / 1.3 * 0.058) * (1 - fem) + ((WT - v_adipose * 0.92) / 1.3 * 0.051) * fem
    if (v_bone < 0) v_bone <- 0.01
    v_heart <- ((155.18 * bsa^1.29) / 1000 / 1.05) * (1 - fem) + ((124.13 * bsa^1.242) / 1000 / 1.05) * fem
    v_kidney <- ((15.4 + 2.04 * WT + 51.8 * htm^2) / 1000) / 1.05
    v_muscle <- ((0.244 * WT + 7.8 * htm - 0.098 * AGE + 3.3) / 1.04) * (1 - fem) +
      ((0.244 * WT + 7.8 * htm - 0.098 * AGE - 3.3) / 1.04) * fem
    v_skin <- (bsa / 1000) * 45.655 + (bsa / 1000) * 1240
    v_spleen <- (6.516 * WT^0.797) / 1000
    v_gut <- ((0.021 * (1 - fem) + 0.027 * fem) * (WT - v_adipose * 0.92)) / 1.05
    v_liver <- (1072.8 * bsa - 345.7) / 1000
    v_blood_half <- ((((13.1 * HT + 18.05 * WT - 480) / 0.5723) / 1000) / 2) * (1 - fem) +
      ((((35.5 * HT + 2.27 * WT - 3382) / 0.6178) / 1000) / 2) * fem
    v_venous <- v_blood_half
    v_arterial <- v_blood_half
    # Rest of body (body density 1 kg/L)
    v_other <- WT - v_venous - v_arterial - v_liver - v_gut - v_spleen - v_skin - v_muscle -
      v_kidney - v_heart - v_bone - v_adipose - v_lung - v_brain
    if (v_other < 0) v_other <- 0.01

    # Blood flows (L/h): cardiac output times fixed fractional tissue flows
    q_co <- bsa * 60 * (3 - 0.01 * (AGE - 20))
    q_brain <- q_co * 0.12
    q_adipose <- q_co * (0.05 * (1 - fem) + 0.085 * fem)
    q_bone <- q_co * 0.05
    q_heart <- q_co * (0.04 * (1 - fem) + 0.05 * fem)
    q_kidney <- q_co * (0.19 * (1 - fem) + 0.17 * fem)
    q_muscle <- q_co * (0.17 * (1 - fem) + 0.12 * fem)
    q_skin <- q_co * 0.05
    q_spleen <- q_co * (0.02 * (1 - fem) + 0.03 * fem)
    q_gut <- q_co * (0.15 * (1 - fem) + 0.17 * fem)
    q_ha <- q_co * 0.065
    q_liver <- q_co * (0.235 * (1 - fem) + 0.265 * fem)
    q_rest <- q_co * (0.095 * (1 - fem) + 0.09 * fem)
    q_lung <- q_co

    # CSF flows (L/h). In the code the ODE bulk flow is 0.25 x the
    # individual production rate WITHOUT the bulk-flow variability; the
    # varied bulk flow enters only the spinal-sink and cranial-sink flows.
    q_csf_prod <- exp(lq_csf_prod) * exp(etalq_csf_prod)
    q_bulk <- exp(lf_bulk) * q_csf_prod
    q_bulk_var <- q_bulk * exp(etalf_bulk)
    q_ssink <- exp(lf_ssink) * (0.75 * q_csf_prod + q_bulk_var) * exp(etalf_ssink)
    q_sout <- exp(lf_sout) * q_ssink * exp(etalf_sout)
    q_sin <- q_ssink + q_sout
    q_csink <- 0.75 * q_csf_prod + q_bulk_var - q_sin + q_sout

    # =====================================================================
    # Poulin and Theil partition coefficients (code calcKp_PT, type 3 =
    # monoprotic base). Tissue composition from Table S2 (Poulin and Theil
    # column): fractional water FW, neutral lipid NL and neutral
    # phospholipid NP per tissue; plasma FW 0.945, NL 0.0035, NP 0.0023.
    # As coded: adipose uses the ionisation-corrected olive-oil:water
    # coefficient D* at the adipose pH 7 (logD = 1.115 logP - 1.35) and
    # is multiplied by fu; the other tissues use the octanol:water P and
    # are multiplied by fu / fut, fut = 1 / (1 + 0.5 (1 - fu) / fu).
    # =====================================================================
    pow <- 10^logp
    dstar <- 10^((1.115 * logp - 1.35) - log10(1 + 10^(pka - 7)))
    fut <- 1 / (1 + ((1 - fu) / fu) * 0.5)
    den_pl <- pow * (0.0035 + 0.3 * 0.0023) + (0.945 + 0.7 * 0.0023)
    kp_adipose <- ((dstar * (0.79 + 0.3 * 0.002) + (0.18 + 0.7 * 0.002)) /
      (dstar * (0.0035 + 0.3 * 0.0023) + (0.945 + 0.7 * 0.0023))) * fu
    kp_bone <- ((pow * (0.074 + 0.3 * 0.0011) + (0.439 + 0.7 * 0.0011)) / den_pl) * (fu / fut)
    kp_gut <- ((pow * (0.0487 + 0.3 * 0.0163) + (0.718 + 0.7 * 0.0163)) / den_pl) * (fu / fut)
    kp_heart <- ((pow * (0.0115 + 0.3 * 0.0166) + (0.758 + 0.7 * 0.0166)) / den_pl) * (fu / fut)
    kp_kidney <- ((pow * (0.0207 + 0.3 * 0.0162) + (0.783 + 0.7 * 0.0162)) / den_pl) * (fu / fut)
    kp_liver <- ((pow * (0.0348 + 0.3 * 0.0252) + (0.751 + 0.7 * 0.0252)) / den_pl) * (fu / fut)
    kp_lung <- ((pow * (0.003 + 0.3 * 0.009) + (0.811 + 0.7 * 0.009)) / den_pl) * (fu / fut)
    kp_muscle <- ((pow * (0.0238 + 0.3 * 0.0072) + (0.76 + 0.7 * 0.0072)) / den_pl) * (fu / fut)
    kp_skin <- ((pow * (0.0284 + 0.3 * 0.0111) + (0.7177 + 0.7 * 0.0111)) / den_pl) * (fu / fut)
    kp_spleen <- ((pow * (0.0201 + 0.3 * 0.0198) + (0.787 + 0.7 * 0.0198)) / den_pl) * (fu / fut)
    # Rest of body: arithmetic mean of the ten tissue Kp values (code Kprest_human)
    kp_other <- (kp_adipose + kp_lung + kp_bone + kp_heart + kp_kidney + kp_muscle + kp_skin + kp_spleen + kp_gut + kp_liver) / 10

    # =====================================================================
    # Clearance and brain transfer.
    # =====================================================================
    cl <- exp(lcl + etalcl)
    fu_bb <- fu / bpr
    psb <- exp(lpsb)
    psc <- f_psc * psb
    pse <- exp(lpse)

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

    # =====================================================================
    # ODEs (code @dynamics, written in concentrations there; multiplied
    # through by the compartment volume so the states are amounts in mg).
    # IV doses enter venous blood (Methods 2.3.1). The code's elimination
    # term (tvCL / BP) * Cblv removes CL x plasma concentration.
    # =====================================================================
    d/dt(venous) <- q_adipose * c_ad / kp_adipose * bpr + q_bone * c_bo / kp_bone * bpr +
      q_heart * c_he / kp_heart * bpr + q_kidney * c_ki / kp_kidney * bpr +
      q_muscle * c_mu / kp_muscle * bpr + q_skin * c_sk / kp_skin * bpr +
      q_liver * c_li / kp_liver * bpr + q_brain * c_bb + q_rest * c_ot / kp_other * bpr -
      q_lung * c_ven - (cl / bpr) * c_ven
    d/dt(arterial) <- q_lung * c_lung / kp_lung * bpr -
      (q_rest + q_brain + q_adipose + q_bone + q_heart + q_kidney + q_muscle + q_skin + q_spleen + q_gut + q_ha) * c_art
    d/dt(lung) <- q_lung * (c_ven - c_lung / kp_lung * bpr)
    d/dt(adipose) <- q_adipose * (c_art - c_ad / kp_adipose * bpr)
    d/dt(bone) <- q_bone * (c_art - c_bo / kp_bone * bpr)
    d/dt(heart) <- q_heart * (c_art - c_he / kp_heart * bpr)
    d/dt(kidney) <- q_kidney * (c_art - c_ki / kp_kidney * bpr)
    d/dt(muscle) <- q_muscle * (c_art - c_mu / kp_muscle * bpr)
    d/dt(skin) <- q_skin * (c_art - c_sk / kp_skin * bpr)
    d/dt(spleen) <- q_spleen * (c_art - c_sp / kp_spleen * bpr)
    d/dt(other) <- q_rest * (c_art - c_ot / kp_other * bpr)
    d/dt(gut) <- q_gut * (c_art - c_gu / kp_gut * bpr)
    d/dt(liver) <- q_ha * c_art + q_gut * c_gu / kp_gut * bpr + q_spleen * c_sp / kp_spleen * bpr -
      q_liver * c_li / kp_liver * bpr
    # 4Brain model (Fig 3): brain blood, brain mass, cranial and spinal CSF
    d/dt(brain_vascular) <- q_brain * (c_art - c_bb) + psb * (fu_br * c_bm - fu_bb * c_bb) +
      psc * (fu_csf * c_ccsf - fu_bb * c_bb) + q_ssink * c_scsf + q_csink * c_ccsf
    d/dt(brain) <- psb * (fu_bb * c_bb - fu_br * c_bm) + pse * (fu_csf * c_ccsf - fu_br * c_bm) -
      q_bulk * c_bm
    d/dt(brain_csf_sas_cranial) <- pse * (fu_br * c_bm - fu_csf * c_ccsf) + psc * (fu_bb * c_bb - fu_csf * c_ccsf) +
      q_bulk * c_bm + q_sout * c_scsf - q_sin * c_ccsf - q_csink * c_ccsf
    d/dt(brain_csf_sas_spinal) <- q_sin * c_ccsf - q_sout * c_scsf - q_ssink * c_scsf

    # =====================================================================
    # Outputs (code @derived). Cecf scales the total brain-mass
    # concentration by Vecf / Vbm, as printed (Methods 2.3.2: ECF
    # concentrations approximated from the ECF volume).
    # =====================================================================
    Cecf <- c_bm * (v_ecf / v_bm)
    Ccsf <- c_ccsf
    Ccsf_spinal <- c_scsf
    Cc <- c_ven / bpr
    Cc ~ prop(propSd)
  })
}
