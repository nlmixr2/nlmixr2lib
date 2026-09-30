Cristea_2021_oat13_renal_ontogeny_pbpk <- function() {
  description <- paste(
    "PBPK (renal clearance, popPBPK; NONMEM 7.3). Paediatric renal",
    "clearance (CLR) of the OAT1,3 probe pair clavulanic acid and",
    "amoxicillin in critically ill children aged 1 month to 15 years",
    "(Cristea 2021), plus the paper's PBPK predictions for two further",
    "OAT1,3 substrates, piperacillin and cefazolin. CLR is glomerular",
    "filtration (GF) plus active tubular secretion (ATS) in series:",
    "CLR = fu * GFR + (QR - GFR) * fu * CLsec / (QR + fu * CLsec / BP).",
    "Clavulanic acid is cleared by GF only; amoxicillin by GF and ATS.",
    "OAT1,3 secretion is CLsec = CLint * ont * KW, and the OAT1,3",
    "ontogeny is a sigmoid (Hill) function of postnatal age with",
    "TM50 = 27.3 weeks and Hill 1.17. GFR, renal blood flow, kidney",
    "weight, albumin-driven fu and hematocrit-driven BP all follow",
    "published age functions (Table S1). IIV is on a GF correction",
    "factor and on CLint. The model has no compartments, ODEs or",
    "dosing. Each record gives CLR (L/h) from its covariates (PNA,",
    "GA, WT, BSA); the time column is not used.",
    sep = " "
  )
  reference <- paste(
    "Cristea S, Krekels EHJ, Allegaert K, De Paepe P, de Jaeger A,",
    "De Cock P, Knibbe CAJ. Estimation of Ontogeny Functions for Renal",
    "Transporters Using a Combined Population Pharmacokinetic and",
    "Physiology-Based Pharmacokinetic Approach: Application to OAT1,3.",
    "AAPS J. 2021;23(3):65. doi:10.1208/s12248-021-00595-9.",
    "Model equations 1-6 and the Results estimates are from the main",
    "article; the system-parameter age functions (Table S1) and the",
    "retrospective IVIVE for piperacillin and cefazolin (equation S1 and",
    "Table S2) are from the Supplemental Material (ESM 1). The individual",
    "CLR values used as dependent variables came from the popPK model of",
    "De Cock PAJG et al. Antimicrob Agents Chemother. 2015;59(11):7027-7035",
    "(doi:10.1128/AAC.01368-15).",
    sep = " "
  )
  vignette <- "Cristea_2021_oat13_renal_ontogeny_pbpk"
  units <- list(
    time = paste(
      "n/a -- the model has no dynamics. Each record gives the renal",
      "clearance implied by its covariates, so the time column is ignored."
    ),
    dosing = "n/a (renal clearance model; no dosing)",
    concentration = paste(
      "n/a (no drug concentration). Outputs are renal clearances in L/h",
      "(clr_clav, clr_amox, clr_amox_gf, clr_amox_ats, clr_pip, clr_cef),",
      "the dimensionless OAT1,3 ontogeny fraction ont_oat13, and the",
      "system quantities gfr, qr (mL/min) and kw (g)."
    )
  )

  covariateData <- list(
    PNA = list(
      description = "Postnatal age (chronological age since birth)",
      units = "months",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Drives the OAT1,3 ontogeny (in weeks; Eq. 6), the albumin function",
        "(in days; Table S1), and the cardiac output, renal fraction of",
        "cardiac output and hematocrit functions (in years; Table S1).",
        "model() converts months to those units using 1 month = 30.4375",
        "days. PNA must be > 0: the albumin function is ln(AGE in days).",
        "Cohort range 1 month to 15 years (median 2.6 years)."
      ),
      source_name = "PNA"
    ),
    GA = list(
      description = "Gestational age at birth",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Only used to form postmenstrual age (PMA = GA + PNA in weeks) for",
        "the GFR maturation function in Table S1. The Table S1 legend",
        "assumes GA = 40 weeks when the individual value is unknown."
      ),
      source_name = "GA"
    ),
    WT = list(
      description = "Current body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric size term of the GFR function, (WT / 70)^0.63, and the",
        "kidney-weight function (Table S1)."
      ),
      source_name = "WT"
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Scales cardiac output, and so renal blood flow (Table S1). The",
        "paper does not say how BSA was computed. It is taken as a",
        "covariate so the user's own BSA formula is used."
      ),
      source_name = "BSA"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 50L,
    n_studies = 1L,
    age_range = "1 month to 15 years (median 2.6 years)",
    weight_range = "Not reported in Cristea 2021 (see De Cock 2015 for the source cohort)",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not reported",
    disease_state = paste(
      "Critically ill children in paediatric intensive care, without renal",
      "dysfunction, who received intravenous amoxicillin-clavulanic acid",
      "at a fixed 10:1 dose ratio (De Cock 2015 cohort, Ghent University",
      "Hospital)."
    ),
    dose_range = "Clinical amoxicillin-clavulanic acid dosing (10:1 amoxicillin:clavulanic acid)",
    regions = "Belgium",
    notes = paste(
      "The dependent variables were the individual post hoc CLR values of",
      "clavulanic acid and amoxicillin from the De Cock 2015 popPK model,",
      "one pair per child (Cristea 2021 Methods). The piperacillin and",
      "cefazolin outputs are external predictions, not fitted. Piperacillin:",
      "47 critically ill children aged 2.5 months to 15 years (median 2.83",
      "years; De Cock 2017). Cefazolin: 26 near-term neonates with GA over",
      "35 weeks and PNA 1 to 30 days (median 8 days; De Cock 2014)."
    )
  )

  ini({
    # ---- OAT1,3 ontogeny and intrinsic clearance (Results) ----
    # CLint is the Results value '15.8 ml/h/g kidney'. It is per gram of
    # kidney, so it already includes the proximal-tubule-cells-per-gram
    # factor PTCPGK = 60 x 10^6 cells/g of Eq. 5 (Table S1). Checks:
    # Figure 1 plots CLsec/KW in ml/h/g kidney and its typical curve
    # levels off at about 15.8 at 15 years. 15.8 ml/h/g / 60 = 4.39
    # uL/min per 10^6 cells, the '4.4' amoxicillin CLint quoted in the
    # Discussion.
    lclint_oat13 <- log(15.8)
    label("Adult OAT1,3-mediated in vivo intrinsic clearance of amoxicillin per gram kidney, CLint x PTCPGK (mL/h/g)") # Results: 'CLint, OAT1,3, in vivo was estimated to be 15.8 ml/h/g kidney (RSE% of 5%) at 15 years'
    lpna50 <- log(27.3)
    label("Postnatal age at half of adult OAT1,3 activity, TM50 (weeks)") # Results: 'half of the adult capacity at a PNA of 27.3 weeks (RSE of 28%)'
    lhill <- log(1.17)
    label("Hill coefficient of the OAT1,3 ontogeny function (unitless)") # Results: 'hill exponent of 1.17 (%RSE of 36%)'

    # ---- Glomerular filtration correction for critical illness (Eq. 3) ----
    lfcorr_gfr <- log(1.83)
    label("Typical GF correction factor, theta_corr, for critically ill children (unitless)") # Results: 'GF correction factor ... estimated at 1.83 (RSE of 4%)'

    # ---- Drug-specific inputs, adult values (Methods; Table S2) ----
    fu_clav <- fixed(0.75)
    label("Adult fraction unbound in plasma, clavulanic acid (unitless)") # Methods: 'fu clav.acid = 0.75'; Table S1
    fu_amox <- fixed(0.82)
    label("Adult fraction unbound in plasma, amoxicillin (unitless)") # Methods: 'fu amox = 0.82'; Table S1
    bpr_amox <- fixed(0.55)
    label("Adult blood-to-plasma ratio, amoxicillin (unitless)") # Methods: 'BPamox.= 0.55'
    fu_pip <- fixed(0.8)
    label("Adult fraction unbound in plasma, piperacillin (unitless)") # Methods: 'fu, adult of 0.8 (27)'; Table S2
    fu_cef <- fixed(0.31)
    label("Adult fraction unbound in plasma, cefazolin (unitless)") # Methods: '0.31 (25) for piperacillin and cefazolin, respectively'; Table S2
    bpr_pip <- fixed(0.55)
    label("Adult blood-to-plasma ratio, piperacillin (unitless)") # Methods: 'BP adult of 0.55 for both drugs'
    bpr_cef <- fixed(0.55)
    label("Adult blood-to-plasma ratio, cefazolin (unitless)") # Methods: 'BP adult of 0.55 for both drugs'

    # ---- Retrospective IVIVE for piperacillin and cefazolin (Eq. S1; Table S2) ----
    clint_vitro_pip <- fixed(1.95)
    label("In vitro OAT3-mediated intrinsic clearance, piperacillin (uL/min/mg protein)") # Table S2: 1.95 uL/min/mg protein (Wen 2018)
    clint_vitro_cef <- fixed(7.1)
    label("In vitro OAT3-mediated intrinsic clearance, cefazolin (uL/min/mg protein)") # Table S2: 7.1 uL/min/mg protein (Mathialagan 2017)
    prot_hek <- fixed(0.25)
    label("Protein content of OAT-transfected HEK293 cells (mg protein per 10^6 cells)") # Supplement IVIVE: '0.25 mg protein per 10^6 ... HEK293 ... cells'
    raf_oat3 <- fixed(4.6)
    label("Relative activity factor for OAT3 (unitless)") # Supplement IVIVE: 'For OAT3, the reported value was 4.6'
    aaf_pip <- fixed(11.6)
    label("Activity adjustment factor, piperacillin (unitless)") # Table S2: AAF = 11.6
    aaf_cef <- fixed(0.65)
    label("Activity adjustment factor, cefazolin (unitless)") # Table S2: AAF = 0.65
    ptcpgk <- fixed(60)
    label("Proximal tubule cells per gram kidney (10^6 cells/g)") # Table S1: 'PTCPKG = 60 (adult value)'

    # ---- Adult reference for the albumin and hematocrit scaling ----
    # Not printed. The fu scaling of Table S1 needs an adult albumin
    # concentration, and the BP scaling needs a partition coefficient kp,
    # which is back-solved from the adult BP at an adult hematocrit. Both
    # adult values come from the printed Table S1 age functions at this
    # reference age (albumin 44.2 g/L, hematocrit 41.0%). See the vignette
    # section 'Assumptions and deviations'.
    age_adult_ref <- fixed(30)
    label("Adult reference age for the adult albumin and hematocrit values (years)") # Not printed; the maintainers' choice (vignette, 'Assumptions and deviations')

    # ---- Inter-individual variability (Results) ----
    # Printed as CV%. omega^2 = log(1 + CV^2).
    etalfcorr_gfr ~ 0.05785 # Results: 'GF correction factor ... with an IIV of 24.4%'; log(1 + 0.244^2) = 0.05785
    etalclint_oat13 ~ 0.47999 # Results: 'CLint ... with an IIV of 78.5%'; log(1 + 0.785^2) = 0.47999
  })

  model({
    # ---- Age axes (PNA in months -> weeks, days, years) ----
    pna_wk <- PNA * 30.4375 / 7
    pna_d <- PNA * 30.4375
    age_yr <- PNA * 30.4375 / 365.25
    pma_wk <- GA + pna_wk

    # ---- Individual parameters ----
    fcorr <- exp(lfcorr_gfr + etalfcorr_gfr)
    clint <- exp(lclint_oat13 + etalclint_oat13)
    clint_typ <- exp(lclint_oat13)
    pna50 <- exp(lpna50)
    hill <- exp(lhill)

    # ---- OAT1,3 ontogeny, Eq. 6 with COV = PNA in weeks ----
    ont_oat13 <- pna_wk^hill / (pna_wk^hill + pna50^hill)

    # ---- Glomerular filtration rate (mL/min), Table S1 ----
    gfr <- 112 * (WT / 70)^0.63 * pma_wk^3.3 / (pma_wk^3.3 + 55.4^3.3)

    # ---- Kidney weight (g), Table S1 ----
    kw <- 1050 * (4.214 * WT^0.823 + 4.456 * WT^0.795) / 1000

    # ---- Renal blood flow (mL/min), Table S1 ----
    # Cardiac output as printed, AGE in years. The equation is in L/h per
    # m^2, not the 'ml/min' of the Table S1 header: in mL/min, QR in a
    # 15-year-old would be about 60 mL/min, below GFR, which makes
    # QR - GFR in Eq. 1 negative. The printed form has no bracket around
    # the two exponentials and uses 184 rather than Simcyp's 184.974; the
    # vignette quantifies the effect of that difference.
    co <- BSA * (110 + 184 * exp(-0.0378 * age_yr) - exp(-0.24477 * age_yr))
    # Fraction of cardiac output to the kidneys, in percent, averaged
    # over the male and female functions.
    fr_m <- 4.53 + 14.63 * age_yr / (0.1888 + age_yr)
    fr_f <- 4.53 + 13 * age_yr^1.15 / (0.188^1.15 + age_yr^1.15)
    fr <- (fr_m + fr_f) / 2 / 100
    qr <- co * fr * 1000 / 60

    # ---- Plasma albumin (g/L) and fraction unbound, Table S1 ----
    # Albumin uses AGE in days (Table S1 legend). The adult albumin is the
    # same function at age_adult_ref.
    hsa <- 1.1287 * log(pna_d) + 33.746
    hsa_adult <- 1.1287 * log(age_adult_ref * 365.25) + 33.746
    fu_clav_i <- 1 / (1 + (1 - fu_clav) * hsa / (hsa_adult * fu_clav))
    fu_amox_i <- 1 / (1 + (1 - fu_amox) * hsa / (hsa_adult * fu_amox))
    fu_pip_i <- 1 / (1 + (1 - fu_pip) * hsa / (hsa_adult * fu_pip))
    fu_cef_i <- 1 / (1 + (1 - fu_cef) * hsa / (hsa_adult * fu_cef))

    # ---- Hematocrit (fraction) and blood-to-plasma ratio, Table S1 ----
    # Male and female functions, in percent, averaged; AGE in years.
    hct_m <- 53 - ((43 * age_yr^1.12 / (0.05^1.12 + age_yr^1.12)) *
      (1 - 0.93 * age_yr^0.25 / (0.10^0.25 + age_yr^0.25)))
    hct_f <- 53 - ((37.4 * age_yr^1.12 / (0.05^1.12 + age_yr^1.12)) *
      (1 - 0.80 * age_yr^0.25 / (0.10^0.25 + age_yr^0.25)))
    hct <- (hct_m + hct_f) / 2 / 100
    hct_m_ad <- 53 - ((43 * age_adult_ref^1.12 / (0.05^1.12 + age_adult_ref^1.12)) *
      (1 - 0.93 * age_adult_ref^0.25 / (0.10^0.25 + age_adult_ref^0.25)))
    hct_f_ad <- 53 - ((37.4 * age_adult_ref^1.12 / (0.05^1.12 + age_adult_ref^1.12)) *
      (1 - 0.80 * age_adult_ref^0.25 / (0.10^0.25 + age_adult_ref^0.25)))
    hct_adult <- (hct_m_ad + hct_f_ad) / 2 / 100
    # BP = 1 + hct * (fu * kp - 1), with kp chosen so that the adult
    # fu and hematocrit give back the printed adult BP.
    kp_amox <- (bpr_amox - 1 + hct_adult) / (hct_adult * fu_amox)
    kp_pip <- (bpr_pip - 1 + hct_adult) / (hct_adult * fu_pip)
    kp_cef <- (bpr_cef - 1 + hct_adult) / (hct_adult * fu_cef)
    bp_amox_i <- 1 + hct * (fu_amox_i * kp_amox - 1)
    bp_pip_i <- 1 + hct * (fu_pip_i * kp_pip - 1)
    bp_cef_i <- 1 + hct * (fu_cef_i * kp_cef - 1)

    # ---- Clavulanic acid: GF only, Eq. 3 (mL/min -> L/h) ----
    clr_clav <- gfr * fu_clav_i * fcorr * 60 / 1000

    # ---- Amoxicillin: GF + OAT1,3 secretion, Eqs. 4 and 5 ----
    # CLint is per gram kidney (PTCPGK already included), mL/h -> mL/min.
    # As printed in Eq. 4, the secretion term uses the uncorrected GFR.
    clsec_amox <- clint * ont_oat13 * kw / 60
    clr_amox_gf <- gfr * fu_amox_i * fcorr * 60 / 1000
    clr_amox_ats <- (qr - gfr) * fu_amox_i * clsec_amox /
      (qr + fu_amox_i * clsec_amox / bp_amox_i) * 60 / 1000
    clr_amox <- clr_amox_gf + clr_amox_ats

    # ---- Piperacillin and cefazolin: PBPK predictions, Eqs. 1-2 and S1 ----
    # Typical-value predictions, so no theta_corr and no IIV (Methods,
    # 'Predictive Properties of the OAT1,3 Ontogeny Function').
    # CLint in vivo (uL/min per 10^6 cells) = CLint in vitro x protein x
    # RAF x AAF (Eq. S1). CLsec (mL/min) = CLint x ont x PTCPGK x KW / 1000.
    clsec_pip <- clint_vitro_pip * prot_hek * raf_oat3 * aaf_pip *
      ont_oat13 * ptcpgk * kw / 1000
    clsec_cef <- clint_vitro_cef * prot_hek * raf_oat3 * aaf_cef *
      ont_oat13 * ptcpgk * kw / 1000
    clr_pip <- (fu_pip_i * gfr + (qr - gfr) * fu_pip_i * clsec_pip /
      (qr + fu_pip_i * clsec_pip / bp_pip_i)) * 60 / 1000
    clr_cef <- (fu_cef_i * gfr + (qr - gfr) * fu_cef_i * clsec_cef /
      (qr + fu_cef_i * clsec_cef / bp_cef_i)) * 60 / 1000

    # Typical OAT1,3 secretion per gram kidney (mL/h/g), the Figure 1 curve.
    clsec_per_kw_typ <- clint_typ * ont_oat13
  })
}
