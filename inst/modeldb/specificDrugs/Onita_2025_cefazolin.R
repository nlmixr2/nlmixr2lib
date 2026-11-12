Onita_2025_cefazolin <- function() {
  description <- paste(
    "Hybrid plasma-plus-cerebrospinal-fluid (CSF) population PK model for",
    "intravenous cefazolin in adults (Onita 2025). The plasma layer is the",
    "two-compartment Lanoiselee 2021 total-hip-arthroplasty model carried",
    "with every value fixed (clearance a power function of CKD-EPI creatinine",
    "clearance centred at 80 mL/min/1.73 m^2); a third CSF compartment is fed",
    "from the central concentration by the physiological CSF flow (0.021 L/h,",
    "fixed) and drains by the same flow divided by a CSF-to-plasma partition",
    "coefficient (0.0525, fixed from paired serum/CSF data in 8 adults with",
    "suspected meningitis), so steady-state CSF concentration is 0.0525 x",
    "plasma. Only the interindividual variances of the three CSF parameters",
    "were estimated, against the paired serum/CSF data. Used for a 100% T > MIC",
    "probability-of-target-attainment analysis of high-dose regimens for",
    "methicillin-susceptible Staphylococcus aureus meningitis."
  )
  reference <- paste(
    "Onita T, Ikawa K, Ishihara N, Tamaki H, Yano T. Cerebrospinal",
    "pharmacokinetic modeling and pharmacodynamic simulation of high-dose",
    "cefazolin for meningitis caused by methicillin-susceptible Staphylococcus",
    "aureus. Antibiotics (Basel). 2025;14(10):1008.",
    "doi:10.3390/antibiotics14101008. PMCID PMC12561440.",
    "All parameter values are from Table 2; the ODE system is the display",
    "equation set in Methods section 4.1. The plasma layer is Lanoiselee J et",
    "al. Sci Rep. 2021;11:19763 (doi:10.1038/s41598-021-99162-7), also",
    "available in this package as Lanoiselee_2021_cefazolin.",
    sep = " "
  )
  vignette <- "Onita_2025_cefazolin"
  # Dose in mg and volumes in L, so concentrations are mg/L, numerically equal
  # to the ug/mL the paper uses for concentrations and MICs.
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Plasma concentrations behind the Lanoiselee 2021 layer were measured by
  # liquid chromatography; the CSF (and paired serum) concentrations of Ikuno
  # et al. used to estimate KPCSF were measured by bioassay (Bacillus subtilis
  # ATCC 6633) (Methods 4.1; Supplementary Table S3). Both are total cefazolin:
  # 'The total CSF concentration was not adjusted for the free fraction'
  # (Methods 4.4).
  compartmentData <- list(
    central = list(analyte = "cefazolin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cefazolin", units = "mg", specimen = "plasma", verified = TRUE),
    csf = list(analyte = "cefazolin", units = "mg", specimen = "CSF", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance estimated by the CKD-EPI equation (BSA-normalized)",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Onita 2025 Table 2: CL (L/h) = thetaCL x (CLcr/80)^thetaCLcr on CL,",
        "with CLcr 'calculated by the CKD-EPI formula' (Methods 4.1), the sole",
        "covariate retained by the upstream Lanoiselee 2021 model. Onita prints",
        "the unit as mL/min and simulates at CLcr = 90, 60 and 30 mL/min",
        "(Table 1); Lanoiselee 2021 Table 1 and Onita Supplementary Table S3",
        "report the CKD-EPI value in its native mL/min/1.73 m^2, which is the",
        "unit used here (same reading as Lanoiselee_2021_cefazolin)."
      ),
      source_name = "CLcr"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 8L,
    n_studies = 3L,
    age_range = "38-60 years (mean 47.9) in the 8 paired serum/CSF patients (Supplementary Table S3)",
    sex_female_pct = 50,
    disease_state = paste(
      "Adults with suspected meningitis (Ikuno et al.; 19 serum samples from",
      "7 patients and 24 CSF samples from 8 patients after 1-3 g cefazolin,",
      "sampled 1-9 h post-dose). Two patients had CSF protein > 300 mg/dL",
      "(inflamed; mean CSF/serum ratio 0.081) and six < 100 mg/dL (mean ratio",
      "0.038). The plasma layer comes from 100 adults undergoing total hip",
      "arthroplasty (Lanoiselee 2021) and the CSF flow and volume from 12",
      "children (Cutler et al.)."
    ),
    dose_range = "1-3 g cefazolin in the paired serum/CSF data; simulations 4-12 g/day",
    regions = "Japan (Shimane University Hospital analysis); source data from Japan (Ikuno), France (Lanoiselee) and the USA (Cutler)",
    renal_function = "Not reported for the CSF cohort; the Lanoiselee 2021 plasma cohort had CKD-EPI CrCL mean 83 (range 17-129) mL/min/1.73 m^2",
    notes = paste(
      "Hybrid model assembled from three literature sources (Supplementary",
      "Table S3): Lanoiselee 2021 (plasma two-compartment PK; N = 100, 51 male,",
      "age 24-91 years, weight 48-123 kg), Cutler et al. (CSF flow 0.021 L/h",
      "and CSF volume 0.092 L; N = 12 children aged 4-13 years with subacute",
      "sclerosing panencephalitis or pontine glioma) and Ikuno et al. (paired",
      "serum/CSF; N = 8 adults, 4 male, 4 female). Only the three CSF IIV",
      "variances were estimated (NONMEM 7.4, no priors). The paper notes the",
      "model is not applicable to children."
    )
  )

  ini({
    # ---- Plasma layer (Onita 2025 Table 2, footnote a: 'parameters derived
    # from [25]', i.e. Lanoiselee 2021; every value printed 'Fixed') ----
    lcl <- fixed(log(2.86))
    label("Clearance at CRCL = 80 mL/min/1.73 m^2 (L/h)") # Table 2: thetaCL(L/h) 2.86 Fixed
    lvc <- fixed(log(5.2))
    label("Central volume of distribution (L)") # Table 2: Vcentral (L) = thetaVc 5.2 Fixed
    lq <- fixed(log(10.9))
    label("Intercompartmental clearance central-peripheral (L/h)") # Table 2: Q (L/h) = thetaQ 10.9 Fixed
    lvp <- fixed(log(4.56))
    label("Peripheral volume of distribution (L)") # Table 2: Vperipheral (L) = thetaVp 4.56 Fixed
    e_crcl_cl <- fixed(0.79)
    label("Power exponent on (CRCL/80) for clearance (unitless)") # Table 2: theta CLcr on CL 0.79 Fixed; CL = thetaCL x (CLcr/80)^thetaCLcr

    # ---- CSF layer (Table 2) ----
    lkp_csf <- fixed(log(0.0525))
    label("CSF-to-plasma partition coefficient (unitless)") # Table 2: KPCSF = thetaKPCSF 0.0525 Fixed; Methods 4.2 mean CSF/serum ratio 0.0525
    lqcsf <- fixed(log(0.021))
    label("CSF flow (L/h)") # Table 2: QCSF (L/h) = thetaQCSF 0.021 Fixed (footnote c, physiological value)
    lvcsf <- fixed(log(0.092))
    label("CSF volume (L)") # Table 2: VCSF (L) = thetaVCSF 0.092 Fixed (footnote c, physiological value)

    # ---- Interindividual variability (Table 2, 'exponential error model').
    # Entries are VARIANCES (footnote: eta 'normally distributed with a mean
    # of zero and variance'); the four plasma-layer values are the squares of
    # Lanoiselee 2021's log-scale SDs (0.32, 0.57, 0.66, 0.10). Table 2 lists
    # no covariances, so the Lanoiselee 2021 CL-Vc correlation (0.83) is not
    # carried; see the vignette. ----
    etalcl ~ fixed(0.102) # Table 2: etaCL 0.102, held constant (from Lanoiselee 2021)
    etalvc ~ fixed(0.325) # Table 2: etaVcentral 0.325, held constant (from Lanoiselee 2021)
    etalq ~ fixed(0.437) # Table 2: etaQ 0.437, held constant (from Lanoiselee 2021)
    etalvp ~ fixed(0.01) # Table 2: etaVperipheral 0.01, held constant (from Lanoiselee 2021)
    etalkp_csf ~ 1.39 # Table 2: etaKPCSF 1.39 (RSE 203.0%), SIR 95% CI 0.0188-43.6
    etalqcsf ~ 1.22 # Table 2: etaQCSF 1.22 (RSE 0.75%), SIR 95% CI 0.00245-172.7
    etalvcsf ~ 1.22 # Table 2: etaVCSF 1.22 (RSE 0.38%), SIR 95% CI 0.00443-178.4

    # ---- Residual variability (Table 2, 'proportion error model'). The
    # single epsilon 0.0144 is a variance, so SD = sqrt(0.0144) = 0.12 (the
    # Lanoiselee 2021 12% proportional error). Table 2 lists one epsilon for
    # both matrices; nlmixr2 needs one residual parameter per endpoint, so the
    # same fixed value is carried on each. ----
    propSd <- fixed(0.12)
    label("Proportional residual SD, plasma cefazolin Cc (fraction)") # Table 2: epsilon proportional 0.0144 Fixed; sqrt(0.0144) = 0.12
    propSd_Ccsf <- fixed(0.12)
    label("Proportional residual SD, CSF cefazolin Ccsf (fraction)") # Table 2: same single epsilon proportional 0.0144 Fixed
  })

  model({
    # Individual parameters (Table 2 parameter equations, exponential IIV).
    cl <- exp(lcl + etalcl) * (CRCL / 80)^e_crcl_cl
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)
    kp_csf <- exp(lkp_csf + etalkp_csf)
    qcsf <- exp(lqcsf + etalqcsf)
    vcsf <- exp(lvcsf + etalvcsf)

    # Methods 4.1 display equations, in amounts (mg); the infusion rate Rinf
    # enters central through the dosing records:
    #   dX(central)/dt    = Rinf - (CL/Vc + Q/Vc) X(central) + Q X(peripheral)/Vp
    #   dX(peripheral)/dt = Q X(central)/Vc - Q X(peripheral)/Vp
    #   dX(CSF)/dt        = QCSF X(central)/Vc - QCSF X(CSF)/VCSF/KPCSF
    # The CSF exchange does not appear in the central equation: the CSF
    # compartment is a satellite of the plasma model, so the plasma profile is
    # exactly that of Lanoiselee 2021. At equilibrium Ccsf = KPCSF x Cc.
    d/dt(central) <- -(cl / vc + q / vc) * central + q * peripheral1 / vp
    d/dt(peripheral1) <- q * central / vc - q * peripheral1 / vp
    d/dt(csf) <- qcsf * central / vc - qcsf * csf / vcsf / kp_csf

    Cc <- central / vc
    Ccsf <- csf / vcsf

    Cc ~ prop(propSd)
    Ccsf ~ prop(propSd_Ccsf)
  })
}
