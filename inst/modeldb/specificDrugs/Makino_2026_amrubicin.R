Makino_2026_amrubicin <- function() {
  description <- "Simplified population PK/PD model of amrubicin, its active metabolite amrubicinol and neutropenia in Japanese adults with small cell lung cancer during the first treatment cycle (Makino 2026). PK: three-compartment parent amrubicin (5-min IV infusion on days 1-3) whose only elimination is first-order metabolic clearance CLp into a one-compartment amrubicinol model with clearance CLm, all in molar units; body surface area enters CLp, CLm and the amrubicinol volume Vm as power effects centred on the cohort mean of 1.688 m^2. PD: Friberg semi-mechanistic myelosuppression model with four maturation transit compartments, a linear amrubicinol effect Slope * C_AMROH on proliferation, a (Circ0 / Circ)^gamma feedback on proliferation and a second (Circ0 / Circ)^gamma_m feedback that shortens the mean maturation time MMT when neutrophils fall below baseline; ECOG performance status >= 1 shortens MMT. The enterohepatic delay compartment and SLC28A3 genotype covariate of the Makino 2019 predecessor model are removed."
  reference <- paste(
    "Makino Y, Ogawa T, Makihara-Ando R, Sakiyama N, Yamazaki M, Naito C, Takayama H, Ueda G, Kohyama S,",
    "Todo M, Kuwata Y, Hirasaki M, Hamaguchi T. (2026).",
    "Development of a simplified population pharmacokinetic-pharmacodynamic framework for exposure-informed",
    "risk stratification of neutropenia during amrubicin treatment.",
    "Cancer Chemother Pharmacol 96(1):81. doi:10.1007/s00280-026-04931-8.",
    "Residual-error structure, covariate centring and PD inter-individual variability from the predecessor model:",
    "Makino Y, Makihara-Ando R, Ogawa T, Sato H, Goto Y, Kanda S, Horinouchi H, Fujiwara Y, Ohe Y, Yamamoto N.",
    "(2019). Individual optimal dose of amrubicin to prevent severe",
    "neutropenia in Japanese patients with lung cancer. Cancer Sci 110:3573-3583. doi:10.1111/cas.14194.",
    sep = " "
  )
  vignette <- "Makino_2026_amrubicin"
  # Box-Cox shape parameter of the ANC residual error (same pattern as
  # Wilbaux_2014_ovarianCancer_ca125); not a residual SD.
  paper_specific_residual_sds <- c("lambda_ANC")
  units <- list(time = "h", dosing = "umol", concentration = "umol/L", anc = "10^9/L")
  # Unit notes.
  # 1. PK. Makino 2019 Section 2.5 states that amrubicin and amrubicinol
  #    concentrations were converted to molar concentrations for modelling, and
  #    Makino 2026 Figure 2 plots both analytes in nmol/mL (= umol/L). Doses
  #    therefore enter in umol of amrubicin; with volumes in L and clearances in
  #    L/h, central / vc is in umol/L. Formation of amrubicinol from amrubicin is
  #    1:1 on a molar basis, so the metabolite needs no molecular-weight factor.
  #    To dose a label dose expressed in mg of amrubicin hydrochloride, divide by
  #    519.93 g/mol and multiply by 1000 (umol).
  # 2. PD. Slope multiplies the amrubicinol concentration in umol/L (L/umol).
  #    MMT is in hours (Makino 2019 Discussion: 148.3 vs 172.9 h). Circ0 and ANC
  #    are in 10^3 cells/mm^3 = 10^9 cells/L.
  covariateData <- list(
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effects (BSA / 1.688)^e on CLp, CLm and Vm. Makino 2019 Section 2.5: continuous covariates were centred at the population mean and entered as power models, P = tvP * (cov / cov_mean)^cov_theta. The mean, 1.688 m^2, is the Makino 2026 Results value (printed as 1.688 +/- 0.522 m^2; the SD is not reproducible from the 1.21-1.97 m^2 range in Supplementary Table 1). Calculated by the DuBois formula in the Makino 2026 application (Supplementary Fig. 2, Patient info tab).",
      source_name = "BSA"
    ),
    ECOG_GE1 = list(
      description = "ECOG performance status >= 1 indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ECOG PS 0)",
      notes = "Categorical exponential effect exp(covPS * PS) on the mean maturation time MMT (Makino 2026 Table 1 equation image; Makino 2019 Section 2.5 categorical-covariate template). The analysis cohort contained only PS 0 (n = 20) and PS 1 (n = 30) patients (Supplementary Table 1), so PS enters as a 0/1 indicator; extrapolation to PS >= 2 is not supported by the data.",
      source_name = "PS"
    )
  )
  covariatesDataExcluded <- list(
    SNP_SLC28A3_RS7853758 = list(
      description = "SLC28A3 rs7853758 minor-allele carrier",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-carrier)",
      notes = "Covariate on the enterohepatic delay rate constant Kdc in the Makino 2019 predecessor model; removed together with the enterohepatic delay compartment in Makino 2026 (Results 'Model development and structural modification').",
      source_name = "SLC28A3"
    )
  )
  compartmentData <- list(
    central = list(analyte = "amrubicin", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "amrubicin", units = "umol", specimen = "tissue", verified = TRUE),
    peripheral2 = list(analyte = "amrubicin", units = "umol", specimen = "tissue", verified = TRUE),
    central_amroh = list(analyte = "amrubicinol", units = "umol", specimen = "plasma", verified = TRUE),
    prol = list(analyte = "neutrophil progenitor cells", units = "10^9/L", specimen = "tissue", verified = TRUE),
    transit1 = list(analyte = "maturing neutrophils", units = "10^9/L", specimen = "tissue", verified = TRUE),
    transit2 = list(analyte = "maturing neutrophils", units = "10^9/L", specimen = "tissue", verified = TRUE),
    transit3 = list(analyte = "maturing neutrophils", units = "10^9/L", specimen = "tissue", verified = TRUE),
    transit4 = list(analyte = "maturing neutrophils", units = "10^9/L", specimen = "tissue", verified = TRUE),
    circ = list(analyte = "neutrophils", units = "10^9/L", specimen = "whole blood", verified = TRUE)
  )
  population <- list(
    species = "human",
    n_subjects = 50L,
    n_studies = 1L,
    age_range = "39-81 years (median 63.5)",
    weight_range = "34.9-80.3 kg (median 59.2)",
    sex_female_pct = 16,
    race_ethnicity = c(Asian = 100),
    disease_state = "Lung cancer treated with second-line or later amrubicin monotherapy: small cell lung cancer (42/50; limited disease 13, extensive disease 29), large-cell neuroendocrine carcinoma (6), squamous cell carcinoma (1), other (1). ECOG PS 0 (20) or 1 (30); 35/50 had one prior chemotherapy regimen.",
    dose_range = "Amrubicin hydrochloride 40 mg/m^2 (n = 46) or 35 mg/m^2 (n = 4) as a 5-minute IV infusion on days 1-3 of the first cycle.",
    regions = "Japan (National Cancer Center Hospital, Tokyo)",
    bsa_range = "1.21-1.97 m^2 (median 1.65; mean 1.688)",
    observations = "388 plasma concentrations of amrubicin and amrubicinol and 357 absolute neutrophil counts from the first treatment cycle.",
    notes = "Demographics from Makino 2026 Supplementary Table 1 (identical to Makino 2019 Table 2). The Makino 2026 Results text quotes 80% male and 65% PS 0, which disagree with Supplementary Table 1 (84% male, 40% PS 0); the table values are used here. Same 50-patient dataset as Makino 2019. Fitted in Phoenix NLME 1.3, PK first and PD second with the PK parameters fixed at their final estimates."
  )
  ini({
    # ---- Amrubicin PK: Makino 2026 Table 1A ----
    lcl <- log(19.1); label("CLp: metabolic clearance of amrubicin to amrubicinol at BSA 1.688 m^2 (L/h)") # Table 1A 'tvCLp, L/h' = 19.1 (bootstrap CV 3.0%)
    lvc <- log(8.3); label("Vp: central volume of distribution of amrubicin (L)") # Table 1A 'tvVp, L' = 8.3 (bootstrap CV 11.0%)
    lq <- log(10.6); label("CL2: intercompartmental clearance between central and peripheral-1 (L/h)") # Table 1A 'tvCL2, L/h' = 10.6 (bootstrap CV 10.4%)
    lvp <- log(31.7); label("V2: peripheral-1 volume of distribution of amrubicin (L)") # Table 1A 'tvV2, L' = 31.7 (bootstrap CV 6.2%)
    lq2 <- log(44.0); label("CL3: intercompartmental clearance between central and peripheral-2 (L/h)") # Table 1A 'tvCL3, L/h' = 44.0 (bootstrap CV 9.8%)
    lvp2 <- log(28.1); label("V3: peripheral-2 volume of distribution of amrubicin (L)") # Table 1A 'tvV3, L' = 28.1 (bootstrap CV 9.9%)
    # ---- Amrubicinol PK: Makino 2026 Table 1A ----
    lcl_amroh <- log(140.3); label("CLm: clearance of amrubicinol at BSA 1.688 m^2 (L/h)") # Table 1A 'tvCLm, L/h' = 140.3 (bootstrap CV 4.3%)
    lvc_amroh <- log(1224.2); label("Vm: central volume of distribution of amrubicinol at BSA 1.688 m^2 (L)") # Table 1A 'tvVm, L' = 1224.2 (bootstrap CV 4.2%)
    # Power effects (BSA / 1.688)^e; centring and power form per Makino 2019 Section 2.5
    e_bsa_cl <- 1.0; label("Exponent of the (BSA / 1.688 m^2) power effect on CLp (unitless)") # Table 1A 'CovBSA (CLp)' = 1.0 (bootstrap CV 20.4%)
    e_bsa_cl_amroh <- 1.9; label("Exponent of the (BSA / 1.688 m^2) power effect on CLm (unitless)") # Table 1A 'CovBSA (CLm)' = 1.9 (bootstrap CV 27.7%)
    e_bsa_vc_amroh <- 1.6; label("Exponent of the (BSA / 1.688 m^2) power effect on Vm (unitless)") # Table 1A 'CovBSA (Vm)' = 1.6 (bootstrap CV 30.3%)
    # ---- Neutrophil PD: Makino 2026 Table 1B ----
    lcirc0 <- log(3.7); label("Circ0: baseline absolute neutrophil count (10^9/L)") # Table 1B 'tvCirc0' = 3.7 (bootstrap CV 5.9%)
    lmtt <- log(180.5); label("MMT: mean maturation time at ECOG PS 0 (h)") # Table 1B 'tvMMT' = 180.5 (bootstrap CV 4.6%)
    lgamma <- log(0.4); label("gamma: feedback exponent of (Circ0 / Circ) on proliferation (unitless)") # Table 1B 'tvGamma, gamma' = 0.4 (bootstrap CV 17.9%)
    lgamma_mtt <- log(0.2); label("gamma_m: feedback exponent of (Circ0 / Circ) on the maturation rate 1 / MMT (unitless)") # Table 1B 'tvGamma-m, gamma_m' = 0.2 (bootstrap CV 28.2%)
    lslope <- log(45.9); label("Slope: linear amrubicinol effect on proliferation (L/umol)") # Table 1B 'tvSlope' = 45.9 (bootstrap CV 13.1%)
    e_ecog_ge1_mtt <- -0.2; label("Exponential effect of ECOG PS >= 1 on MMT (log-scale shift, unitless)") # Table 1B 'CovPS (MMT)' = -0.2 (bootstrap CV 42.0%)
    # ---- Inter-individual variability ----
    # Log-normal IIV; Table 1C gives CV% = 100 * sqrt(exp(omega^2) - 1), so
    # omega^2 = log(CV^2 + 1). No covariances are reported for the 2026 model.
    etalcl ~ 0.03731 # Table 1C 'IIV(CLp)' = 19.5%: log(1 + 0.195^2)
    etalvc ~ 0.14741 # Table 1C 'IIV(Vp)' = 39.8%: log(1 + 0.398^2)
    etalq ~ 0.00475 # Table 1C 'IIV(CL2)' = 6.9%: log(1 + 0.069^2)
    etalvp ~ 0.02374 # Table 1C 'IIV(V2)' = 15.5%: log(1 + 0.155^2)
    etalq2 ~ 0.11683 # Table 1C 'IIV(CL3)' = 35.2%: log(1 + 0.352^2)
    etalvp2 ~ 0.11366 # Table 1C 'IIV(V3)' = 34.7%: log(1 + 0.347^2)
    etalcl_amroh ~ 0.04768 # Table 1C 'IIV(CLm)' = 22.1%: log(1 + 0.221^2)
    etalvc_amroh ~ 0.06345 # Table 1C 'IIV(Vm)' = 25.6%: log(1 + 0.256^2)
    # PD IIV is not tabulated in Makino 2026 (its MMT equation carries eta_MMT);
    # carried forward from the Makino 2019 final PD model fitted to the same
    # 50 patients, whose PD structure Makino 2026 states is unchanged.
    etalcirc0 ~ 0.02138 # Makino 2019 Table 1 final model 'omega Circ0, %' = 14.7: log(1 + 0.147^2)
    etalmtt ~ 0.09405 # Makino 2019 Table 1 final model 'omega MMT, %' = 31.4: log(1 + 0.314^2)
    etalslope ~ 0.24213 # Makino 2019 Table 1 final model 'omega Slope, %' = 52.4: log(1 + 0.524^2)
    # ---- Residual error ----
    # Error model forms from Makino 2019 Section 2.5 (proportional for PK;
    # additive on the Box-Cox scale with lambda 0.2 for ANC). Neither paper
    # reports the residual magnitudes, so they are held at zero.
    propSd <- fixed(0); label("Proportional residual error on amrubicin (fraction); magnitude not reported") # Makino 2019 Section 2.5 proportional error model; value not reported
    propSd_amroh <- fixed(0); label("Proportional residual error on amrubicinol (fraction); magnitude not reported") # Makino 2019 Section 2.5 proportional error model; value not reported
    addSd_ANC <- fixed(0); label("Additive residual error on Box-Cox-transformed ANC (transformed units); magnitude not reported") # Makino 2019 Section 2.5 additive error on the Box-Cox scale; value not reported
    lambda_ANC <- fixed(0.2); label("Box-Cox transformation factor for ANC (unitless)") # Makino 2019 Section 2.5 and Makino 2026 Figure 2 caption: Box-Cox factor 0.2
  })
  model({
    # ---- 1. Individual PK parameters ----
    bsa_ratio <- BSA / 1.688
    cl <- exp(lcl + etalcl) * bsa_ratio^e_bsa_cl
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)
    q2 <- exp(lq2 + etalq2)
    vp2 <- exp(lvp2 + etalvp2)
    cl_amroh <- exp(lcl_amroh + etalcl_amroh) * bsa_ratio^e_bsa_cl_amroh
    vc_amroh <- exp(lvc_amroh + etalvc_amroh) * bsa_ratio^e_bsa_vc_amroh

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2
    kel_amroh <- cl_amroh / vc_amroh

    # ---- 2. PK ODEs (Makino 2026 Figure 1) ----
    # CLp is the only route out of the amrubicin central compartment and feeds
    # amrubicinol 1:1 on a molar basis.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2
    d/dt(central_amroh) <- kel * central - kel_amroh * central_amroh

    Cc <- central / vc
    Cc_amroh <- central_amroh / vc_amroh

    # ---- 3. Individual PD parameters ----
    circ0 <- exp(lcirc0 + etalcirc0)
    mtt <- exp(lmtt + etalmtt + e_ecog_ge1_mtt * ECOG_GE1)
    gamma <- exp(lgamma)
    gamma_mtt <- exp(lgamma_mtt)
    slope <- exp(lslope + etalslope)
    # Figure 1: four transit compartments, MMT = (n + 1) / Ktr, Kprol = Kcirc = Ktr.
    ktr <- 5 / mtt

    # ---- 4. Myelosuppression chain (Makino 2026 Figure 1) ----
    # Drug effect Slope * C_AMROH on proliferation. Feedback (Circ0 / Circ)^gamma
    # acts on proliferation; the second feedback (Circ0 / Circ)^gamma_m shortens
    # the maturation time, i.e. scales the transit rate out of the progenitor
    # pool and through the four maturation compartments (Makino 2019 Section
    # 2.5, after Quartino 2014). Kprol and Kcirc keep the baseline Ktr.
    edrug <- slope * Cc_amroh
    fb_prol <- (circ0 / circ)^gamma
    ktr_fb <- ktr * (circ0 / circ)^gamma_mtt
    d/dt(prol) <- ktr * prol * (1 - edrug) * fb_prol - ktr_fb * prol
    d/dt(transit1) <- ktr_fb * prol - ktr_fb * transit1
    d/dt(transit2) <- ktr_fb * transit1 - ktr_fb * transit2
    d/dt(transit3) <- ktr_fb * transit2 - ktr_fb * transit3
    d/dt(transit4) <- ktr_fb * transit3 - ktr_fb * transit4
    d/dt(circ) <- ktr_fb * transit4 - ktr * circ
    prol(0) <- circ0
    transit1(0) <- circ0
    transit2(0) <- circ0
    transit3(0) <- circ0
    transit4(0) <- circ0
    circ(0) <- circ0

    # ---- 5. Observations and error ----
    ANC <- circ
    Cc ~ prop(propSd)
    Cc_amroh ~ prop(propSd_amroh)
    ANC ~ add(addSd_ANC) + boxCox(lambda_ANC)
  })
}
