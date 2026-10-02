Scotcher_2021_gadoxetate_rat_pbpk <- function() {
  description <- paste(
    "Preclinical (rat). Reduced PBPK model (seven compartments,",
    "permeability-limited liver) for the hepatobiliary MRI contrast agent",
    "gadoxetate in the 250 g male Wistar-Han rat after a 25 umol/kg",
    "intravenous dose, given alone or 1 h after a single 10 mg/kg intravenous",
    "rifampicin dose. Blood, spleen and splanchnic extracellular spaces are",
    "perfusion-limited; the rest of the body is split into a vascular and an",
    "interstitial space linked by a permeability-surface product; the liver",
    "has an extracellular space and a hepatocyte space linked by saturable",
    "OATP-mediated active uptake (linearised, CLactive) and bidirectional",
    "passive diffusion, with biliary (Mrp2) excretion from the hepatocyte and",
    "renal clearance from blood. Active uptake, biliary clearance and PS were",
    "fitted (naive pooled) to dynamic contrast-enhanced MRI Delta R1 profiles",
    "of blood, spleen and liver at 4.7 T and 7 T; the defaults are the",
    "simultaneous control plus rifampicin fit (Table 4), in which rifampicin",
    "inhibits active uptake by 96 percent. Delta R1 outputs at both field",
    "strengths are derived with the ex vivo relaxivities of Table 1.",
    "Deterministic: no between-animal variability or residual error was",
    "estimated."
  )
  reference <- paste(
    "Scotcher D, Melillo N, Tadimalla S, Darwich AS, Ziemian S, Ogungbenro K,",
    "Schutz G, Sourbron S, Galetin A. Physiologically Based Pharmacokinetic",
    "Modeling of Transporter-Mediated Hepatic Disposition of Imaging Biomarker",
    "Gadoxetate in Rats. Mol Pharm. 2021;18(8):2997-3009.",
    "doi:10.1021/acs.molpharmaceut.1c00206 (PMC8397403).",
    "Model equations (S1-S3) and physiological parameters (Table S2) are in",
    "the paper's Supporting Information (mp1c00206_si_001.pdf)."
  )
  vignette <- "Scotcher_2021_gadoxetate_rat_pbpk"

  # State names: `blood` = systemic blood; `vp_remainder` / `is_remainder` =
  # rest-of-body (ROB) vascular / interstitial space; `is_spleen` = spleen
  # extracellular space (organ blood + interstitium); `is_liver` / `int_liver`
  # = liver extracellular space / hepatocytes; `urine` and `a_bile` are
  # cumulative excretion records. The splanchnic lump (stomach + gut +
  # pancreas extracellular space) has no registered organ token, so it is
  # declared paper-specific.
  paper_specific_compartments <- c("is_splanchnic")

  # Every state is an AMOUNT in umol (equation system S1 is written in
  # amounts); concentrations are umol/L (= uM). Delta R1 is in 1/s.
  units <- list(
    time = "h",
    dosing = "umol",
    concentration = "umol/L"
  )

  covariateData <- list(
    CONMED_RIFAMPICIN_SD = list(
      description = paste(
        "1 = gadoxetate given 1 h after a single 10 mg/kg intravenous",
        "rifampicin dose (inhibitory phase), 0 = gadoxetate alone (control",
        "phase)"
      ),
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Selects the rifampicin-phase estimates of active hepatic uptake",
        "(CLactive,inh) and biliary clearance (CLbiliary,inh) of Table 4 in",
        "place of the control estimates; PS is shared by both phases. The",
        "paper estimated the two phases as separate parameter values rather",
        "than as an inhibition-constant model of rifampicin concentration, so",
        "the indicator applies only to the studied design (single 10 mg/kg IV",
        "rifampicin 1 h before gadoxetate)."
      ),
      source_name = "rifampicin phase"
    )
  )

  compartmentData <- list(
    blood = list(analyte = "gadoxetate", units = "umol", specimen = "whole blood", verified = TRUE),
    vp_remainder = list(analyte = "gadoxetate", units = "umol", specimen = "whole blood", verified = TRUE),
    is_remainder = list(analyte = "gadoxetate", units = "umol", specimen = "tissue", verified = TRUE),
    is_spleen = list(analyte = "gadoxetate", units = "umol", specimen = "tissue", verified = TRUE),
    is_splanchnic = list(analyte = "gadoxetate", units = "umol", specimen = "tissue", verified = TRUE),
    is_liver = list(analyte = "gadoxetate", units = "umol", specimen = "tissue", verified = TRUE),
    int_liver = list(analyte = "gadoxetate", units = "umol", specimen = "tissue", verified = TRUE),
    urine = list(analyte = "gadoxetate", units = "umol", specimen = "urine", verified = TRUE),
    a_bile = list(analyte = "gadoxetate", units = "umol", specimen = "bile", verified = TRUE)
  )

  population <- list(
    species = "rat (Wistar-Han)",
    n_subjects = 76L,
    n_studies = 1L,
    weight_median = "250 g (model reference rat; Table 3 footnote e and Supporting Information Section 5)",
    sex_female_pct = 0,
    disease_state = "healthy",
    dose_range = paste(
      "gadoxetate 25 umol/kg IV, alone (control) or 1 h after rifampicin",
      "10 mg/kg IV"
    ),
    regions = "multicentre preclinical imaging study (two sites per field strength)",
    notes = paste(
      "DCE-MRI Delta R1 profiles of blood, spleen and liver: control arm 43",
      "profiles at 4.7 T (n = 33 animals) and 52 at 7 T (n = 43 animals; some",
      "animals scanned twice); rifampicin arm 7 profiles at 4.7 T and 6 at",
      "7 T (Experimental Section 'DCE-MRI Dataset'; Figures 4-6 captions).",
      "In vitro uptake kinetics were measured separately in plated",
      "hepatocytes from male Sprague-Dawley rats (250-300 g, n = 4; Table 2)",
      "and supply CLpassive and fu,liv,cell. Fits were naive pooled."
    )
  )

  ini({
    # Gadoxetate transporter parameters from the simultaneous fit of the
    # control and rifampicin phases (Table 4; mean of 1000 case-bootstrap
    # samples, CV in parentheses). The control-only top-down fit (Table 3:
    # CLactive 2.17, CLbiliary 0.07, PS 0.62 L/h) and the bottom-up IVIVE
    # values (Table 3: CLactive 0.23, CLbiliary 0.014, PS 0.014 L/h) are
    # reproduced in the vignette by overriding these three parameters.
    lps_act_inf <- log(2.38)
    label("Active (OATP) hepatic uptake clearance CLactive, control phase (L/h)") # Table 4 CLactive control 2.38 (CV 13.6%)
    lps_act_inf_inh <- log(0.095)
    label("Active hepatic uptake clearance CLactive,inh, rifampicin phase (L/h)") # Table 4 CLactive with rifampicin 0.095 (CV 16.1%)
    lcl_bile <- log(0.07)
    label("Biliary clearance from hepatocytes CLbiliary, control phase (L/h)") # Table 4 CLbiliary control 0.07 (CV 3.1%)
    lcl_bile_inh <- log(0.08)
    label("Biliary clearance from hepatocytes CLbiliary,inh, rifampicin phase (L/h)") # Table 4 CLbiliary with rifampicin 0.08 (CV 16.7%)
    lps_remainder <- log(0.71)
    label("Rest-of-body vascular-interstitial permeability-surface product PS (L/h)") # Table 4 PS 0.71 (CV 5.7%), shared by both phases

    # Parameters held at in vitro / literature values in every fit
    # (Experimental Section 'PBPK Analysis Overview'; Table 3 footnotes d, e).
    lps_dif <- fixed(log(0.014))
    label("Bidirectional passive diffusion clearance across the hepatocyte membrane CLpassive (L/h)") # Table 3 CLpassive 0.014 (Table 2 mean 0.193 uL/min/10^6 cells x 120e6 cells/g x 40 g/kg x 0.25 kg)
    fu_liver_cell <- fixed(0.648)
    label("Gadoxetate fraction unbound in hepatocytes fu,liv,cell (fraction)") # Table 3 fu,liv,cell 0.648 = Table 2 mean of four animals
    lcl_renal <- fixed(log(0.17))
    label("Renal clearance from blood CLr (L/h)") # Table 3 CLr 0.17; footnote e: 36.7 mL/min/kg x fe 0.305 x 0.25 kg

    # System parameters for the 250 g reference rat (Supporting Information
    # Section 5 and Table S2).
    qc <- fixed(6.62)
    label("Cardiac output QCO (L/h)") # Table S2 blood / lungs blood flow 6.62 L/h
    hct <- fixed(0.4183)
    label("Haematocrit (fraction)") # Supporting Information Section 5: haematocrit 0.4183
  })

  model({
    # ------------------------------------------------------------------
    # 1. Physiology of the 250 g rat (Table S2). Organ volume (L) = mass
    #    (g) / density (kg/L) / 1000; vascular and interstitial volumes are
    #    the organ volume times the Table S2 vascular / interstitial
    #    fractions (Kawai 1998).
    # ------------------------------------------------------------------
    v_blood <- 15.77 / 1 / 1000 # Table S2 blood 15.77 g, density 1 kg/L

    # Rest of the body = lungs, brain, heart, kidneys, bone, muscle, skin, fat
    # (Supporting Information Section 2). Lumped volumes are sums (eq S2).
    # Table S2 rows: vascular fraction * mass / density.
    v_vas_remainder <- (0.262 * 1.25 / 1.0505 + # lungs
      0.037 * 1.44 / 1.0355 + # brain
      0.262 * 0.84 / 1.03 + # heart
      0.105 * 1.84 / 1.05 + # kidneys
      0.041 * 15 / 1.4303 + # bone
      0.026 * 101 / 1.041 + # muscle
      0.019 * 47.58 / 1.183 + # skin
      0.01 * 15.63 / 0.916) / 1000 # fat
    # Table S2 rows: interstitial fraction * mass / density.
    v_is_remainder <- (0.188 * 1.25 / 1.0505 + # lungs
      0.004 * 1.44 / 1.0355 + # brain
      0.1 * 0.84 / 1.03 + # heart
      0.2 * 1.84 / 1.05 + # kidneys
      0.1 * 15 / 1.4303 + # bone
      0.12 * 101 / 1.041 + # muscle
      0.302 * 47.58 / 1.183 + # skin
      0.135 * 15.63 / 0.916) / 1000 # fat

    # Spleen: extracellular space = organ blood + interstitium (Table S2
    # spleen 0.5 g, density 1.054, vascular 0.282, interstitial 0.15).
    v_spleen <- 0.5 / 1.054 / 1000
    vb_spleen <- 0.282 * v_spleen
    vi_spleen <- 0.15 * v_spleen
    v_is_spleen <- vb_spleen + vi_spleen

    # Splanchnic organs = stomach + gut + pancreas (Supporting Information
    # Section 2); extracellular space only.
    vb_splanchnic <- (0.032 * 1.15 / 1.05 + # stomach
      0.024 * 5.6 / 1.043 + # gut
      0.18 * 0.8 / 1.045) / 1000 # pancreas
    vi_splanchnic <- (0.1 * 1.15 / 1.05 + # stomach
      0.094 * 5.6 / 1.043 + # gut
      0.12 * 0.8 / 1.045) / 1000 # pancreas
    v_is_splanchnic <- vb_splanchnic + vi_splanchnic

    # Liver: extracellular = organ blood + interstitium; hepatocytes = whole
    # liver minus extracellular (Supporting Information Section 2). Table S2
    # liver 9.15 g, density 1.08, vascular 0.115, interstitial 0.163.
    v_liver <- 9.15 / 1.08 / 1000
    vb_liver <- 0.115 * v_liver
    vi_liver <- 0.163 * v_liver
    v_is_liver <- vb_liver + vi_liver
    v_int_liver <- v_liver - v_is_liver

    # Tissue(extracellular)-to-blood partition coefficients, eq S3: gadoxetate
    # does not enter red cells and plasma equilibrates with the interstitium.
    # For the splanchnic lump the volume-weighted mean of eq S2 over the
    # organs' extracellular volumes equals eq S3 applied to the summed volumes.
    # The ROB vascular-to-blood coefficient is 1 (Supporting Information
    # Section 2), so it does not appear below.
    kp_spleen <- (vi_spleen + vb_spleen * (1 - hct)) / ((vi_spleen + vb_spleen) * (1 - hct))
    kp_splanchnic <- (vi_splanchnic + vb_splanchnic * (1 - hct)) /
      ((vi_splanchnic + vb_splanchnic) * (1 - hct))
    kp_liver <- (vi_liver + vb_liver * (1 - hct)) / ((vi_liver + vb_liver) * (1 - hct))

    # Blood flows (L/h, Table S2). Hepatic blood flow = hepatic artery 0.14
    # plus the portal tributaries (spleen 0.053; stomach 0.08 + gut 0.85 +
    # pancreas 0.03). The ROB flow is the remainder of cardiac output so
    # that blood leaving (QCO) equals blood returning (Qrob + Qh).
    q_spleen <- 0.053
    q_splanchnic <- 0.08 + 0.85 + 0.03
    q_hepatic_artery <- 0.14
    qh <- q_hepatic_artery + q_spleen + q_splanchnic
    q_remainder <- qc - qh

    # ------------------------------------------------------------------
    # 2. Gadoxetate parameters; the rifampicin phase swaps in its own
    #    active-uptake and biliary estimates (Table 4).
    # ------------------------------------------------------------------
    ps_act_inf <- exp(lps_act_inf * (1 - CONMED_RIFAMPICIN_SD) + lps_act_inf_inh * CONMED_RIFAMPICIN_SD)
    cl_bile <- exp(lcl_bile * (1 - CONMED_RIFAMPICIN_SD) + lcl_bile_inh * CONMED_RIFAMPICIN_SD)
    ps_remainder <- exp(lps_remainder)
    ps_dif <- exp(lps_dif)
    cl_renal <- exp(lcl_renal)

    # ------------------------------------------------------------------
    # 3. Concentrations (umol/L).
    # ------------------------------------------------------------------
    c_blood <- blood / v_blood
    c_vas_remainder <- vp_remainder / v_vas_remainder
    c_is_remainder <- is_remainder / v_is_remainder
    c_is_spleen <- is_spleen / v_is_spleen
    c_is_splanchnic <- is_splanchnic / v_is_splanchnic
    c_is_liver <- is_liver / v_is_liver
    c_int_liver <- int_liver / v_int_liver

    # ------------------------------------------------------------------
    # 4. ODEs: equation system S1 (amounts). The blood equation's hepatic
    #    return term is printed with subscripts 'liv,int' and 'Q_liv'; it
    #    is the hepatic venous outflow of the liver extracellular space,
    #    Qh * a_liv,extr / (V_liv,extr * K_liv,extr-b), as in the liver
    #    extracellular equation of S1 and eq 6.
    # ------------------------------------------------------------------
    d/dt(blood) <- -qc * c_blood - cl_renal * c_blood + q_remainder * c_vas_remainder +
      qh * c_is_liver / kp_liver
    d/dt(vp_remainder) <- q_remainder * (c_blood - c_vas_remainder) -
      ps_remainder * (c_vas_remainder - c_is_remainder)
    d/dt(is_remainder) <- ps_remainder * (c_vas_remainder - c_is_remainder)
    d/dt(is_spleen) <- q_spleen * (c_blood - c_is_spleen / kp_spleen)
    d/dt(is_splanchnic) <- q_splanchnic * (c_blood - c_is_splanchnic / kp_splanchnic)
    d/dt(is_liver) <- (qh - q_spleen - q_splanchnic) * c_blood +
      q_spleen * c_is_spleen / kp_spleen +
      q_splanchnic * c_is_splanchnic / kp_splanchnic -
      qh * c_is_liver / kp_liver -
      ps_act_inf * c_is_liver -
      ps_dif * (c_is_liver - fu_liver_cell * c_int_liver)
    d/dt(int_liver) <- ps_act_inf * c_is_liver +
      ps_dif * (c_is_liver - fu_liver_cell * c_int_liver) -
      cl_bile * fu_liver_cell * c_int_liver
    d/dt(urine) <- cl_renal * c_blood
    d/dt(a_bile) <- cl_bile * fu_liver_cell * c_int_liver

    # ------------------------------------------------------------------
    # 5. Observations. Cc is the gadoxetate blood concentration (umol/L).
    #    Whole-organ concentrations are volume-weighted over the organ's
    #    compartments (eq 7 without the relaxivities).
    # ------------------------------------------------------------------
    Cc <- c_blood
    C_spleen <- c_is_spleen * v_is_spleen / v_spleen
    C_liver <- (is_liver + int_liver) / v_liver

    # Delta R1 (1/s), eq 7, with concentrations converted from umol/L to
    # mmol/L. Ex vivo relaxivities (1/s per mM), Table 1: blood 6.4 at 4.7 T
    # and 6.2 at 7 T; hepatocytes 7.6 at 4.7 T and 6 at 7 T. Spleen and liver
    # extracellular relaxivities equal the blood value (Table 1 footnote a;
    # text below eq 7).
    r1_blood_4p7t <- 6.4
    r1_hep_4p7t <- 7.6
    r1_blood_7t <- 6.2
    r1_hep_7t <- 6
    dR1_blood_4p7t <- c_blood / 1000 * r1_blood_4p7t
    dR1_spleen_4p7t <- c_is_spleen / 1000 * v_is_spleen * r1_blood_4p7t / v_spleen
    dR1_liver_4p7t <- (c_is_liver * v_is_liver * r1_blood_4p7t + c_int_liver * v_int_liver * r1_hep_4p7t) /
      (1000 * v_liver)
    dR1_blood_7t <- c_blood / 1000 * r1_blood_7t
    dR1_spleen_7t <- c_is_spleen / 1000 * v_is_spleen * r1_blood_7t / v_spleen
    dR1_liver_7t <- (c_is_liver * v_is_liver * r1_blood_7t + c_int_liver * v_int_liver * r1_hep_7t) /
      (1000 * v_liver)
  })
}
