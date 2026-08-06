Zhu_2023_irinotecan <- function() {
  description <- paste(
    "Minimal-PBPK. Joint irinotecan and SN-38 population PK model in adults with",
    "metastatic colorectal cancer: two-compartment plasma disposition for each of",
    "irinotecan and its active metabolite SN-38, plus an explicit tumour",
    "compartment. Irinotecan enters the tumour interstitial space by tumour",
    "plasma flow and permeates into tumour cells at a permeability-surface-area",
    "product; SN-38 enters a single lumped tumour-tissue space governed by a",
    "tumour/plasma partition coefficient. Zhu 2023 uses the tumour SN-38",
    "concentration as the driver of an in vitro-derived tumour-killing model;",
    "that in vivo tumour-growth layer is NOT included here because its growth",
    "equation and limiting parameter are not reported anywhere in the paper -",
    "see the vignette Errata. The tumour-distribution parameters are transferred",
    "from the tumour-bearing-mouse fit, see",
    "modellib('Zhu_2023_irinotecan_mouse'); the in vitro drug-effect half is",
    "modellib('Zhu_2023_sn38_organoid')."
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
    species = "human",
    n_subjects = NA_integer_,
    n_studies = NA_integer_,
    disease_state = "metastatic colorectal cancer",
    dose_range = "350 mg/m2 IV every 3 weeks or 125 mg/m2 IV weekly x4 then 2 weeks off (simulated regimens, Zhu 2023 Table 1)",
    notes = paste(
      "The human plasma PK parameters for irinotecan and SN-38 were estimated",
      "from published human concentration-time profiles digitised from the",
      "literature with GetData Graph Digitizer (Zhu 2023 Section 2.1); the",
      "individual source publications, subject counts and demographics are not",
      "reported. Because no human tumour concentrations exist, the tumour",
      "distribution parameters (KP_IRI, KP_SN, PS_IRI) were obtained in",
      "tumour-bearing mice and translated: KP_IRI and KP_SN were held constant",
      "across species and PS_IRI was scaled by tumour volume. Tumour volumes and",
      "tumour plasma flow were fixed at species-specific physiologic values. The",
      "reported interindividual variability therefore mixes between-subject with",
      "between-study variability."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Irinotecan plasma disposition, Zhu 2023 Equations (6)-(7).
    # Zhu 2023 Table 3, "Human Estimates (RSE%)" column. Table 3 labels the
    # second irinotecan and SN-38 volumes VP_IRI / VP_SN ("peripheral")
    # while Equations (7) and (11) call the same volumes VO_IRI / VO_SN
    # (the lumped "other" tissue compartment); they are the same parameter.
    # ------------------------------------------------------------------
    lvc <- log(72.1);   label("Central volume of irinotecan VC_IRI (L)")                            # Zhu 2023 Table 3, VC_IRI human = 72.1 L (RSE 6.78%)
    lvp <- log(93.4);   label("Volume of the lumped other-tissue compartment for irinotecan VO_IRI (L)")  # Zhu 2023 Table 3, VP_IRI human = 93.4 L (RSE 15.2%)
    lcl <- log(22.8);   label("Systemic clearance of irinotecan CL_IRI (L/h)")                      # Zhu 2023 Table 3, CL_IRI human = 22.8 L/h (RSE 5.69%)
    lq  <- log(24.6);   label("Plasma flow to the other-tissue compartment for irinotecan Q_IRI (L/h)")   # Zhu 2023 Table 3, Q_IRI human = 24.6 L/h (RSE 28.8%)

    # ------------------------------------------------------------------
    # SN-38 formation and plasma disposition, Zhu 2023 Equations (10)-(11).
    # ------------------------------------------------------------------
    lcl_met  <- log(0.216); label("Formation clearance of SN-38 from irinotecan CLM_SN (L/h)")      # Zhu 2023 Table 3, CLM_SN human = 0.216 L/h (RSE 51.9%)
    lvc_sn38 <- log(11.2);  label("Central volume of SN-38 VC_SN (L)")                              # Zhu 2023 Table 3, VC_SN human = 11.2 L (RSE 34.5%)
    lvp_sn38 <- log(706);   label("Volume of the lumped other-tissue compartment for SN-38 VO_SN (L)")    # Zhu 2023 Table 3, VP_SN human = 706 L (RSE 52.7%)
    lcl_sn38 <- log(42.8);  label("Systemic clearance of SN-38 CL_SN (L/h)")                        # Zhu 2023 Table 3, CL_SN human = 42.8 L/h (RSE 32.5%)
    lq_sn38  <- log(43.5);  label("Plasma flow to the other-tissue compartment for SN-38 Q_SN (L/h)")     # Zhu 2023 Table 3, Q_SN human = 43.5 L/h (RSE 30.7%)

    # ------------------------------------------------------------------
    # Tumour physiology and distribution, Zhu 2023 Equations (8), (9), (12).
    # Volumes converted from mL to L and the permeability-surface-area
    # product from cm3/h to L/h (1 cm3 = 1 mL = 0.001 L).
    # ------------------------------------------------------------------
    v_tumor_is   <- fixed(0.002);  label("Tumour interstitial-space volume V_TIS (L)")              # Zhu 2023 Table 3, V_TIS human = 2 mL (assumed)
    v_tumor_cell <- fixed(0.008);  label("Tumour cell volume V_TC (L)")                             # Zhu 2023 Table 3, V_TC human = 8 mL (assumed)
    v_tumor      <- fixed(0.010);  label("Total tumour volume V_T (L)")                             # Zhu 2023 Table 3, V_T human = 10 mL (reference [26])
    q_tumor      <- fixed(0.06);   label("Plasma flow between central and tumour compartments Q_T (L/h)")  # Zhu 2023 Table 3, Q_T human = 0.06 L/h (references [26,27])

    lps_tumor <- fixed(log(0.052));      label("Permeability-surface-area product of irinotecan into tumour cells PS_IRI (L/h); scaled to human tumour volume from the mouse estimate")  # Zhu 2023 Table 3, PS_IRI human = 52 cm3/h = 0.052 L/h
    lkp_tumor <- fixed(log(3.43));       label("Tumour/plasma partition coefficient of irinotecan KP_IRI (unitless); taken from the mouse fit, held constant across species")            # Zhu 2023 Table 3, KP_IRI human = 3.43
    lkp_tumor_sn38 <- fixed(log(7.32));  label("Tumour/plasma partition coefficient of SN-38 KP_SN (unitless); taken from the mouse fit, held constant across species")                  # Zhu 2023 Table 3, KP_SN human = 7.32

    fu      <- fixed(0.35); label("Fraction unbound of irinotecan in plasma fu_IRI (unitless); literature value")  # Zhu 2023 Table 3, fu_IRI = 0.35 (reference [28])
    fu_sn38 <- fixed(0.05); label("Fraction unbound of SN-38 in plasma fu_SN (unitless); literature value")        # Zhu 2023 Table 3, fu_SN = 0.05 (reference [28])

    # ------------------------------------------------------------------
    # Interindividual variability. Monolix reports omega as the standard
    # deviation of the log-scale random effect, so the variance entered here
    # is the square of the tabulated IIV. VP_IRI and VP_SN carry no IIV in
    # Zhu 2023 Table 3, and the human tumour-distribution parameters are all
    # fixed.
    # ------------------------------------------------------------------
    etalvc      ~ 2.6244;    # Zhu 2023 Table 3, VC_IRI human IIV = 1.62 (RSE 39.3%); variance = 1.62^2
    etalcl      ~ 0.022201;  # Zhu 2023 Table 3, CL_IRI human IIV = 0.149 (RSE 27.4%); variance = 0.149^2
    etalq       ~ 0.463761;  # Zhu 2023 Table 3, Q_IRI human IIV = 0.681 (RSE 35.8%); variance = 0.681^2
    etalcl_met  ~ 0.443556;  # Zhu 2023 Table 3, CLM_SN human IIV = 0.666 (RSE 29%); variance = 0.666^2
    etalvc_sn38 ~ 0.019321;  # Zhu 2023 Table 3, VC_SN human IIV = 0.139 (RSE 71.9%); variance = 0.139^2
    etalcl_sn38 ~ 0.25;      # Zhu 2023 Table 3, CL_SN human IIV = 0.5 (RSE 50.6%); variance = 0.5^2
    etalq_sn38  ~ 0.228484;  # Zhu 2023 Table 3, Q_SN human IIV = 0.478 (RSE 32.4%); variance = 0.478^2

    # Zhu 2023 reports no residual-error estimates, so both residual SDs are
    # held at zero. See the vignette Errata.
    propSd <- fixed(0);         label("Proportional residual error on plasma irinotecan (fraction; not reported in the source)")  # Zhu 2023 reports no residual-error model
    propSd_sn38 <- fixed(0);    label("Proportional residual error on plasma SN-38 (fraction; not reported in the source)")       # Zhu 2023 reports no residual-error model
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

    ps_tumor       <- exp(lps_tumor)
    kp_tumor       <- exp(lkp_tumor)
    kp_tumor_sn38  <- exp(lkp_tumor_sn38)

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
    # homogenate (interstitial space plus cells over the total tumour volume).
    # Zhu 2023 Figure 3D simulates this in humans but has no data to fit it, so
    # it is a derived output rather than a declared endpoint.
    Ctumor <- (tumor_is + tumor_cell) / v_tumor

    Cc ~ prop(propSd)
    Cc_sn38 ~ prop(propSd_sn38)
  })
}
