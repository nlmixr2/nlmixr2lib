Zhu_2023_irinotecan_mouse <- function() {
  description <- paste(
    "Preclinical (tumour-bearing mouse). QSP / minimal-PBPK. Joint irinotecan",
    "and SN-38 minimal-PBPK model: two-compartment plasma disposition for each",
    "of irinotecan and its active metabolite SN-38, plus an explicit tumour",
    "compartment. Irinotecan enters the tumour interstitial space by tumour",
    "plasma flow and permeates into tumour cells at a permeability-surface-area",
    "product; SN-38 enters a single lumped tumour-tissue space governed by a",
    "tumour/plasma partition coefficient. This is the preclinical fit that",
    "supplied the tumour-distribution parameters (KP_IRI, KP_SN, PS_IRI) for",
    "the human model, see modellib('Zhu_2023_irinotecan'); it carries no PD",
    "component because Zhu 2023 took the drug effect from patient-derived",
    "tumour organoids rather than from mouse tumour growth."
  )
  reference <- paste(
    "Zhu J, Zhang Y, Zhao Y, Zhang J, Hao K, He H. Translational",
    "Pharmacokinetic/Pharmacodynamic Modeling and Simulation of Oxaliplatin and",
    "Irinotecan in Colorectal Cancer. Pharmaceutics. 2023;15(9):2274.",
    "doi:10.3390/pharmaceutics15092274.",
    sep = " "
  )
  vignette <- "Zhu_2023_oxaliplatin_irinotecan_colorectal_cancer"

  paper_specific_compartments <- c("tumor_is", "tumor_cell", "tumor_sn38")

  units <- list(time = "h", dosing = "umol", concentration = "umol/L")

  compartmentData <- list(
    central = list(
      analyte = "irinotecan", units = "umol", specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "irinotecan", units = "umol", specimen = "plasma", verified = TRUE
    ),
    tumor_is = list(
      analyte = "irinotecan", units = "umol", specimen = "tumor", verified = TRUE
    ),
    tumor_cell = list(
      analyte = "irinotecan", units = "umol", specimen = "tumor", verified = TRUE
    ),
    central_sn38 = list(
      analyte = "SN-38", units = "umol", specimen = "plasma", verified = TRUE
    ),
    peripheral1_sn38 = list(
      analyte = "SN-38", units = "umol", specimen = "plasma", verified = TRUE
    ),
    tumor_sn38 = list(
      analyte = "SN-38", units = "umol", specimen = "tumor", verified = TRUE
    )
  )

  covariateData <- list()

  population <- list(
    species = "mouse (tumour-bearing)",
    n_subjects = NA_integer_,
    n_studies = NA_integer_,
    disease_state = "subcutaneous tumour xenograft",
    dose_range = "intravenous bolus or infusion of irinotecan solution (doses not reported)",
    notes = paste(
      "Plasma and tumour concentration-time profiles of irinotecan and SN-38 in",
      "tumour-bearing mice were digitised from published figures with GetData",
      "Graph Digitizer (Zhu 2023 Section 2.1); the individual source",
      "publications, strain, tumour model, doses and animal counts are not",
      "reported. Studies with co-administered antitumour drugs and whole-blood",
      "only studies were excluded. Tumour volumes and tumour plasma flow were",
      "fixed at species-specific physiologic values (Zhu 2023 Table 3)."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Irinotecan plasma disposition, Zhu 2023 Equations (6)-(7).
    # Zhu 2023 Table 3, "Mouse Estimates (RSE%)" column. Table 3 labels the
    # second irinotecan and SN-38 volumes VP_IRI / VP_SN ("peripheral")
    # while Equations (7) and (11) call the same volumes VO_IRI / VO_SN
    # (the lumped "other" tissue compartment); they are the same parameter.
    # ------------------------------------------------------------------
    lvc <- log(0.0349);  label("Central volume of irinotecan VC_IRI (L)")                           # Zhu 2023 Table 3, VC_IRI mouse = 0.0349 L (RSE 32%)
    lvp <- log(0.0493);  label("Volume of the lumped other-tissue compartment for irinotecan VO_IRI (L)")  # Zhu 2023 Table 3, VP_IRI mouse = 0.0493 L (RSE 25.8%)
    lcl <- log(0.0527);  label("Systemic clearance of irinotecan CL_IRI (L/h)")                     # Zhu 2023 Table 3, CL_IRI mouse = 0.0527 L/h (RSE 19.7%)
    lq  <- log(0.0156);  label("Plasma flow to the other-tissue compartment for irinotecan Q_IRI (L/h)")   # Zhu 2023 Table 3, Q_IRI mouse = 0.0156 L/h (RSE 40.8%)

    # ------------------------------------------------------------------
    # SN-38 formation and plasma disposition, Zhu 2023 Equations (10)-(11).
    # ------------------------------------------------------------------
    lcl_met  <- log(1.65e-4); label("Formation clearance of SN-38 from irinotecan CLM_SN (L/h)")    # Zhu 2023 Table 3, CLM_SN mouse = 1.65e-4 L/h (RSE 93.8%)
    lvc_sn38 <- log(0.00122); label("Central volume of SN-38 VC_SN (L)")                            # Zhu 2023 Table 3, VC_SN mouse = 0.00122 L (RSE 15.9%)
    lvp_sn38 <- log(0.108);   label("Volume of the lumped other-tissue compartment for SN-38 VO_SN (L)")   # Zhu 2023 Table 3, VP_SN mouse = 0.108 L (RSE 33.1%)
    lcl_sn38 <- log(0.0402);  label("Systemic clearance of SN-38 CL_SN (L/h)")                      # Zhu 2023 Table 3, CL_SN mouse = 0.0402 L/h (RSE 19.8%)
    lq_sn38  <- log(0.0369);  label("Plasma flow to the other-tissue compartment for SN-38 Q_SN (L/h)")    # Zhu 2023 Table 3, Q_SN mouse = 0.0369 L/h (RSE 38.2%)

    # ------------------------------------------------------------------
    # Tumour physiology and distribution, Zhu 2023 Equations (8), (9), (12).
    # Volumes converted from mL to L and the permeability-surface-area
    # product from cm3/h to L/h (1 cm3 = 1 mL = 0.001 L).
    # ------------------------------------------------------------------
    v_tumor_is   <- fixed(1e-4);    label("Tumour interstitial-space volume V_TIS (L)")             # Zhu 2023 Table 3, V_TIS mouse = 0.1 mL (assumed)
    v_tumor_cell <- fixed(4e-4);    label("Tumour cell volume V_TC (L)")                            # Zhu 2023 Table 3, V_TC mouse = 0.4 mL (assumed)
    v_tumor      <- fixed(5e-4);    label("Total tumour volume V_T (L)")                            # Zhu 2023 Table 3, V_T mouse = 0.5 mL (reference [26])
    q_tumor      <- fixed(3.38e-3); label("Plasma flow between central and tumour compartments Q_T (L/h)")  # Zhu 2023 Table 3, Q_T mouse = 3.38e-3 L/h (references [26,27])

    lps_tumor      <- log(4.48e-4); label("Permeability-surface-area product of irinotecan into tumour cells PS_IRI (L/h)")  # Zhu 2023 Table 3, PS_IRI mouse = 0.448 cm3/h = 4.48e-4 L/h (RSE >100%)
    lkp_tumor      <- log(3.43);    label("Tumour/plasma partition coefficient of irinotecan KP_IRI (unitless)")             # Zhu 2023 Table 3, KP_IRI mouse = 3.43 (RSE 74.4%)
    lkp_tumor_sn38 <- log(7.32);    label("Tumour/plasma partition coefficient of SN-38 KP_SN (unitless)")                   # Zhu 2023 Table 3, KP_SN mouse = 7.32 (RSE 71.6%)

    fu      <- fixed(0.35); label("Fraction unbound of irinotecan in plasma fu_IRI (unitless); literature value")  # Zhu 2023 Table 3, fu_IRI = 0.35 (reference [28])
    fu_sn38 <- fixed(0.05); label("Fraction unbound of SN-38 in plasma fu_SN (unitless); literature value")        # Zhu 2023 Table 3, fu_SN = 0.05 (reference [28])

    # ------------------------------------------------------------------
    # Interindividual variability. Monolix reports omega as the standard
    # deviation of the log-scale random effect, so the variance entered here
    # is the square of the tabulated IIV. VP_IRI and VP_SN carry no IIV in
    # Zhu 2023 Table 3.
    # ------------------------------------------------------------------
    etalvc      ~ 0.762129;  # Zhu 2023 Table 3, VC_IRI mouse IIV = 0.873 (RSE 20.4%); variance = 0.873^2
    etalcl      ~ 0.393129;  # Zhu 2023 Table 3, CL_IRI mouse IIV = 0.627 (RSE 20.4%); variance = 0.627^2
    etalq       ~ 0.535824;  # Zhu 2023 Table 3, Q_IRI mouse IIV = 0.732 (RSE 26%); variance = 0.732^2
    etalcl_met  ~ 0.724201;  # Zhu 2023 Table 3, CLM_SN mouse IIV = 0.851 (RSE 35.7%); variance = 0.851^2
    etalvc_sn38 ~ 0.244036;  # Zhu 2023 Table 3, VC_SN mouse IIV = 0.494 (RSE 59.3%); variance = 0.494^2
    etalcl_sn38 ~ 0.054756;  # Zhu 2023 Table 3, CL_SN mouse IIV = 0.234 (RSE 78.9%); variance = 0.234^2
    etalq_sn38  ~ 0.367236;  # Zhu 2023 Table 3, Q_SN mouse IIV = 0.606 (RSE 34.3%); variance = 0.606^2

    etalps_tumor      ~ 3.61;    # Zhu 2023 Table 3, PS_IRI mouse IIV = 1.9 (RSE 36.5%); variance = 1.9^2
    etalkp_tumor      ~ 1.7689;  # Zhu 2023 Table 3, KP_IRI mouse IIV = 1.33 (RSE 37.1%); variance = 1.33^2
    etalkp_tumor_sn38 ~ 2.6896;  # Zhu 2023 Table 3, KP_SN mouse IIV = 1.64 (RSE 30%); variance = 1.64^2

    # Zhu 2023 reports no residual-error estimates for any output, so all
    # residual SDs are held at zero. See the vignette Errata.
    propSd <- fixed(0);              label("Proportional residual error on plasma irinotecan (fraction; not reported in the source)")  # Zhu 2023 reports no residual-error model
    propSd_sn38 <- fixed(0);      label("Proportional residual error on plasma SN-38 (fraction; not reported in the source)")       # Zhu 2023 reports no residual-error model
    propSd_Ctumor <- fixed(0);       label("Proportional residual error on tumour irinotecan (fraction; not reported in the source)")  # Zhu 2023 reports no residual-error model
    propSd_Ctumor_sn38 <- fixed(0);  label("Proportional residual error on tumour SN-38 (fraction; not reported in the source)")       # Zhu 2023 reports no residual-error model
  })

  model({
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp)
    cl <- exp(lcl + etalcl)
    q  <- exp(lq + etalq)

    cl_met  <- exp(lcl_met + etalcl_met)
    vc_sn38 <- exp(lvc_sn38 + etalvc_sn38)
    vp_sn38 <- exp(lvp_sn38)
    cl_sn38 <- exp(lcl_sn38 + etalcl_sn38)
    q_sn38  <- exp(lq_sn38 + etalq_sn38)

    ps_tumor      <- exp(lps_tumor + etalps_tumor)
    kp_tumor      <- exp(lkp_tumor + etalkp_tumor)
    kp_tumor_sn38 <- exp(lkp_tumor_sn38 + etalkp_tumor_sn38)

    # States are amounts (umol); the paper writes its equations as V * dC/dt,
    # so each amount derivative is exactly the published right-hand side.
    Cc          <- central / vc
    Cp          <- peripheral1 / vp
    Ctumor_is   <- tumor_is / v_tumor_is
    Ctumor_cell <- tumor_cell / v_tumor_cell
    Cc_sn38     <- central_sn38 / vc_sn38
    Cp_sn38     <- peripheral1_sn38 / vp_sn38
    Ctumor_sn38 <- tumor_sn38 / v_tumor

    # Zhu 2023 Equation (6)
    d/dt(central) <- -(cl + q + q_tumor) * Cc + q * Cp + q_tumor * Ctumor_is
    # Zhu 2023 Equation (7)
    d/dt(peripheral1) <- q * (Cc - Cp)
    # Zhu 2023 Equation (8)
    d/dt(tumor_is) <- q_tumor * (Cc * fu - Ctumor_is) -
      ps_tumor * (Ctumor_is - Ctumor_cell / kp_tumor)
    # Zhu 2023 Equation (9)
    d/dt(tumor_cell) <- ps_tumor * (Ctumor_is - Ctumor_cell / kp_tumor)
    # Zhu 2023 Equation (10). The formation term carries the published
    # VC_IRI / VC_SN volume ratio; see the vignette Errata.
    d/dt(central_sn38) <- cl_met * Cc * (vc / vc_sn38) +
      q_sn38 * Cp_sn38 +
      q_tumor * Ctumor_sn38 / kp_tumor_sn38 -
      (cl_sn38 + q_sn38 + q_tumor) * Cc_sn38
    # Zhu 2023 Equation (11)
    d/dt(peripheral1_sn38) <- q_sn38 * (Cc_sn38 - Cp_sn38)
    # Zhu 2023 Equation (12)
    d/dt(tumor_sn38) <- q_tumor * (Cc_sn38 * fu_sn38 - Ctumor_sn38 / kp_tumor_sn38)

    # Whole-tumour irinotecan concentration, the quantity measured in a tumour
    # homogenate (interstitial space plus cells over the total tumour volume);
    # this is what Zhu 2023 Figure 3B plots as "Tumor irinotecan in mice".
    Ctumor <- (tumor_is + tumor_cell) / v_tumor

    Cc ~ prop(propSd)
    Cc_sn38 ~ prop(propSd_sn38)
    Ctumor ~ prop(propSd_Ctumor)
    Ctumor_sn38 ~ prop(propSd_Ctumor_sn38)
  })
}
