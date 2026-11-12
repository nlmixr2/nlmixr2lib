Park_2025_SNB101 <- function() {
  description <- paste(
    "Eleven-compartment population PK model for SNB-101, a polymeric",
    "nanoparticle co-encapsulating irinotecan and SN-38 (about 1:1 molar),",
    "in adults with advanced solid tumours (Park 2025, phase I, N = 21).",
    "Each analyte is split into a nanoparticle (NP, `_np` states) and a",
    "dissolved (S, unsuffixed states) form fitted to total concentrations:",
    "NP-irinotecan (two compartments) releases into S-irinotecan (three",
    "compartments, parameters fixed from Minichmayr 2021) at krel;",
    "NP-SN-38 shares the NP-irinotecan volumes and distribution clearance",
    "and releases into S-SN-38 (two compartments) at krel_sn38; 3% of",
    "S-irinotecan clearance forms S-SN-38; all S-SN-38 clearance forms",
    "SN-38 glucuronide (central compartment with first-order clearance and",
    "a second, one-way first-order loss at Q6 into a peripheral compartment",
    "that does not return, as in the published control stream). Release is",
    "the only loss of each NP form. Amounts are in",
    "umol and concentrations in nmol/L. The NP:S split of each dose",
    "(irinotecan 0.15:0.85, SN-38 0.98:0.02) is applied in the dosing",
    "records, not in the model.",
    sep = " "
  )
  reference <- paste(
    "Park J, Bae SH, Jeon S, Park YH, Lee DC, Han S. Model-based",
    "prediction of nanoparticle and dissolved",
    "form ratios using total concentration data: a case study of SNB-101.",
    "Front Pharmacol. 2025;16:1556618. doi:10.3389/fphar.2025.1556618.",
    "PMCID: PMC12408558. Parameter values from Table 2; model structure and",
    "the fixed S-irinotecan variances from the NONMEM control stream in",
    "Supplementary Material 1. The fixed S-irinotecan parameters are from",
    "Minichmayr IK, Karlsson MO, Jonsson S. Pharm Res. 2021;38:593-605",
    "(doi:10.1007/s11095-021-03024-w), as cited by the source.",
    sep = " "
  )
  vignette <- "Park_2025_SNB101"
  units <- list(
    time = "h",
    dosing = paste(
      "umol. Give each SNB-101 dose as four records: 0.15 x the irinotecan",
      "amount into central_np, 0.85 x into central, 0.98 x the SN-38 amount",
      "into central_sn38_np and 0.02 x into central_sn38 (Methods 2.2)."
    ),
    concentration = "nmol/L (Cc = total irinotecan, Cc_sn38 = total SN-38, Cc_sn38g = SN-38 glucuronide)"
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      notes = paste(
        "Present in the NONMEM dataset ($INPUT BSA) and used only to convert",
        "the mg/m^2 dose into an amount; no parameter depends on it. The",
        "covariate analysis 'did not identify any statistically significant",
        "influencing factors' (Results 3.1); the screened covariates are not",
        "listed."
      )
    )
  )

  compartmentData <- list(
    central_np = list(analyte = "irinotecan (nanoparticle form)", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1_np = list(
      analyte = "irinotecan (nanoparticle form)",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    central = list(analyte = "irinotecan (dissolved form)", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "irinotecan (dissolved form)", units = "umol", specimen = "tissue", verified = TRUE),
    peripheral2 = list(analyte = "irinotecan (dissolved form)", units = "umol", specimen = "tissue", verified = TRUE),
    central_sn38_np = list(analyte = "SN-38 (nanoparticle form)", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1_sn38_np = list(
      analyte = "SN-38 (nanoparticle form)",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    central_sn38 = list(analyte = "SN-38 (dissolved form)", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1_sn38 = list(analyte = "SN-38 (dissolved form)", units = "umol", specimen = "tissue", verified = TRUE),
    central_sn38g = list(analyte = "SN-38 glucuronide", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1_sn38g = list(analyte = "SN-38 glucuronide", units = "umol", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 21L,
    n_studies = 1L,
    age_range = "about 58-65 years (cohort means 58-65; Table 1)",
    weight_range = "about 68-75 kg (cohort means 68-75; Table 1)",
    sex_female_pct = 33.3,
    disease_state = paste(
      "Advanced solid tumours refractory to, or without, standard therapy;",
      "ECOG 0-1; no prior topoisomerase I inhibitor"
    ),
    dose_range = paste(
      "SNB-101 5-50 mg/m^2 as SN-38 with 8-80 mg/m^2 irinotecan, single",
      "intravenous dose in Cycle 1, seven escalation cohorts (Table 1)"
    ),
    regions = "Not stated in the source",
    notes = paste(
      "Phase I trial NCT04640480. Cohort sizes 1, 1, 1, 5, 7, 3 and 3",
      "(14 men, 7 women; Table 1). Cycle 1 sampling at 0, 0.5, 1, 2, 4, 6, 8,",
      "10, 12, 24, 36, 48, 72, 96, 120, 144 and 169.5 h. Total (NP + S)",
      "irinotecan, total SN-38 and SN-38 glucuronide were measured; the NP",
      "and S forms were never measured separately. The infusion duration is",
      "not stated in the source."
    )
  )

  ini({
    # --- S-irinotecan clearance (estimated) and NP disposition
    lcl <- log(29.7)
    label("Clearance of dissolved (S) irinotecan, CL1 (L/h)") # Table 2, CL1 = 29.7 L/h; $THETA(1)
    lvc_np <- log(54.4)
    label("Central volume of the nanoparticle forms, V1 = V6 (L)") # Table 2, 'V1 (=V6)' = 54.4 L; $THETA(2), $PK V6 = V1
    lvp_np <- log(287)
    label("Peripheral volume of the nanoparticle forms, V2 = V7 (L)") # Table 2, 'V2 (=V7)' = 287 L; $THETA(3), $PK V7 = V2
    lq_np <- log(6.45)
    label("Intercompartmental clearance of the nanoparticle forms, Q1 = Q4 (L/h)") # Table 2, 'Q1 (=Q4)' = 6.45 L/h; $THETA(4), $PK Q4 = Q1
    lkrel <- log(64.4)
    label("Release rate constant, NP-irinotecan to S-irinotecan, k13 (1/h)") # Table 2, k13 = 64.4 1/h; $THETA(5)

    # --- S-irinotecan distribution, fixed from Minichmayr 2021
    lvc <- fixed(log(68.6))
    label("Central volume of dissolved irinotecan, V3 (L)") # Table 2, V3 = '68.6 FIX'; $THETA(6)
    lvp <- fixed(log(67.2))
    label("Rapid peripheral volume of dissolved irinotecan, V4 (L)") # Table 2, V4 = '67.2 FIX'; $THETA(7)
    lvp2 <- fixed(log(127))
    label("Slow peripheral volume of dissolved irinotecan, V5 (L)") # Table 2, V5 = '127 FIX'; $THETA(8)
    lq <- fixed(log(114))
    label("Intercompartmental clearance, S-irinotecan central to rapid peripheral, Q2 (L/h)") # Table 2, Q2 = '114 FIX'; $THETA(9)
    lq2 <- fixed(log(9.89))
    label("Intercompartmental clearance, S-irinotecan central to slow peripheral, Q3 (L/h)") # Table 2, Q3 = '9.89 FIX'; $THETA(10)
    fm_sn38 <- fixed(0.03)
    label("Molar fraction of S-irinotecan clearance forming S-SN-38 (fraction)") # Methods 2.2, 'fixed at 0.03 (3%)'; $PK K38 = CL1*0.03/V3

    # --- SN-38 and SN-38 glucuronide
    lkrel_sn38 <- log(0.712)
    label("Release rate constant, NP-SN-38 to S-SN-38, k68 (1/h)") # Table 2, k68 = 0.712 1/h; $THETA(11)
    lvc_sn38 <- log(40.3)
    label("Central volume of dissolved SN-38, V8 (L)") # Table 2, V8 = 40.3 L; $THETA(12)
    lvp_sn38 <- log(491)
    label("Peripheral volume of dissolved SN-38 and of SN-38 glucuronide, V9 = V11 (L)") # Table 2, 'V9 (=V11)' = 491 L; $THETA(13), $PK V11 = V9
    lq_sn38 <- log(121)
    label("Intercompartmental clearance of dissolved SN-38, Q5 (L/h)") # Table 2, Q5 = 121 L/h; $THETA(14)
    lcl_sn38 <- log(167)
    label("Clearance of dissolved SN-38 to SN-38 glucuronide, CL2 (L/h)") # Table 2, CL2 = 167 L/h; $THETA(15)
    lcl_sn38g <- log(16.9)
    label("Clearance of SN-38 glucuronide, CL3 (L/h)") # Table 2, CL3 = 16.9 L/h; $THETA(17)
    lq_sn38g <- log(11.3)
    label("Clearance of SN-38 glucuronide into its non-returning peripheral compartment, Q6 (L/h)") # Table 2, Q6 = 11.3 L/h; $THETA(18)
    lvc_sn38g <- log(7.84)
    label("Central volume of SN-38 glucuronide, V10 (L)") # Table 2, V10 = 7.84 L; $THETA(19)

    # --- Between-subject variability, log-normal and diagonal. Table 2
    # prints BSV as sqrt(exp(omega^2) - 1): the five fixed S-irinotecan
    # variances in the control stream $OMEGA reproduce the printed 15.4%,
    # 16.2%, 29.1%, 86.1% and 49.5% exactly, so the estimated ones are
    # converted back the same way. Etas the control stream fixes to 0 (Q1,
    # k13, V8, V9, Q5, Q6) are omitted.
    etalcl ~ log(1 + 0.347^2) # Table 2, CL1 BSV 34.7%
    etalvc_np ~ log(1 + 0.409^2) # Table 2, V1 BSV 40.9%
    etalvp_np ~ log(1 + 0.945^2) # Table 2, V2 BSV 94.5%
    etalvc ~ fixed(0.023302744) # $OMEGA 6 (V3 FIX); Table 2 V3 BSV 15.4% FIX
    etalvp ~ fixed(0.025809644) # $OMEGA 7 (V4 FIX); Table 2 V4 BSV 16.2% FIX
    etalvp2 ~ fixed(0.08134926) # $OMEGA 8 (V5 FIX); Table 2 V5 BSV 29.1% FIX
    etalq ~ fixed(0.554644021) # $OMEGA 9 (Q2 FIX); Table 2 Q2 BSV 86.1% FIX
    etalq2 ~ fixed(0.21915561) # $OMEGA 10 (Q3 FIX); Table 2 Q3 BSV 49.5% FIX
    etalkrel_sn38 ~ log(1 + 0.396^2) # Table 2, k68 BSV 39.6%
    etalcl_sn38 ~ log(1 + 0.449^2) # Table 2, CL2 BSV 44.9%
    etalcl_sn38g ~ log(1 + 1.159^2) # Table 2, CL3 BSV 115.9%
    etalvc_sn38g ~ log(1 + 1.221^2) # Table 2, V10 BSV 122.1%

    # --- Residual error, one proportional SD per analyte ($ERROR
    # W = SQRT(add^2 + prop^2 * IPRED^2) with $SIGMA 1 FIX)
    propSd <- 0.21
    label("Proportional residual error, total irinotecan (fraction)") # Table 2, PROP1 = 0.21; $THETA(21)
    propSd_sn38 <- 0.416
    label("Proportional residual error, total SN-38 (fraction)") # Table 2, PROP2 = 0.416; $THETA(23)
    propSd_sn38g <- 0.308
    label("Proportional residual error, SN-38 glucuronide (fraction)") # Table 2, PROP3 = 0.308; $THETA(25)
  })

  model({
    # Individual parameters ($PK)
    cl <- exp(lcl + etalcl)
    vc_np <- exp(lvc_np + etalvc_np)
    vp_np <- exp(lvp_np + etalvp_np)
    q_np <- exp(lq_np)
    krel <- exp(lkrel)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    vp2 <- exp(lvp2 + etalvp2)
    q <- exp(lq + etalq)
    q2 <- exp(lq2 + etalq2)
    krel_sn38 <- exp(lkrel_sn38 + etalkrel_sn38)
    vc_sn38 <- exp(lvc_sn38)
    vp_sn38 <- exp(lvp_sn38)
    q_sn38 <- exp(lq_sn38)
    cl_sn38 <- exp(lcl_sn38 + etalcl_sn38)
    cl_sn38g <- exp(lcl_sn38g + etalcl_sn38g)
    q_sn38g <- exp(lq_sn38g)
    vc_sn38g <- exp(lvc_sn38g + etalvc_sn38g)

    # The two nanoparticle forms share one set of disposition parameters
    # (V6 = V1, V7 = V2, Q4 = Q1); the SN-38G peripheral compartment shares
    # the S-SN-38 peripheral volume (V11 = V9).
    vp_sn38g <- vp_sn38

    # Micro-constants, named after the control stream's K-block
    k12 <- q_np / vc_np # K12 = K67
    k21 <- q_np / vp_np # K21 = K76
    k34 <- q / vc
    k43 <- q / vp
    k35 <- q2 / vc
    k53 <- q2 / vp2
    k30 <- cl * (1 - fm_sn38) / vc
    k38 <- cl * fm_sn38 / vc
    k89 <- q_sn38 / vc_sn38
    k98 <- q_sn38 / vp_sn38
    # FMS = 1 FIX: all S-SN-38 clearance forms SN-38 glucuronide (K80 = 0).
    k810 <- cl_sn38 / vc_sn38
    k100 <- cl_sn38g / vc_sn38g
    k1011 <- q_sn38g / vc_sn38g
    k1110 <- q_sn38g / vp_sn38g

    # $DES, compartments 1-11 in control-stream order. Amounts are molar, so
    # the irinotecan -> SN-38 -> SN-38G conversions carry no weight ratio.
    d/dt(central_np) <- -(krel + k12) * central_np + k21 * peripheral1_np
    d/dt(peripheral1_np) <- k12 * central_np - k21 * peripheral1_np
    d/dt(central) <- krel * central_np - (k34 + k35 + k38 + k30) * central + k43 * peripheral1 + k53 * peripheral2
    d/dt(peripheral1) <- k34 * central - k43 * peripheral1
    d/dt(peripheral2) <- k35 * central - k53 * peripheral2
    d/dt(central_sn38_np) <- -(k12 + krel_sn38) * central_sn38_np + k21 * peripheral1_sn38_np
    d/dt(peripheral1_sn38_np) <- k12 * central_sn38_np - k21 * peripheral1_sn38_np
    d/dt(central_sn38) <- k38 * central + krel_sn38 * central_sn38_np - (k89 + k810) * central_sn38 + k98 * peripheral1_sn38
    d/dt(peripheral1_sn38) <- k89 * central_sn38 - k98 * peripheral1_sn38
    # DADT(10) of the published control stream has no return term from
    # compartment 11, so Q6 acts as a second elimination route from the
    # SN-38G central compartment and peripheral1_sn38g (and V11) never feed
    # back into any observation. The estimates were obtained with this code,
    # and it reproduces the observed SN-38G tail (Figure 1F) where a
    # reversible exchange over-predicts it about five-fold; see the vignette.
    d/dt(central_sn38g) <- k810 * central_sn38 - (k100 + k1011) * central_sn38g
    d/dt(peripheral1_sn38g) <- k1011 * central_sn38g - k1110 * peripheral1_sn38g

    # Form-specific concentrations (nmol/L = 1000 x umol/L; $PK S1 = V1/1000).
    # Never measured; the source predicts them from the fit (Figure 5, Table 3).
    Cc_np <- 1000 * central_np / vc_np
    Cc_dissolved <- 1000 * central / vc
    Cc_sn38_np <- 1000 * central_sn38_np / vc_np
    Cc_sn38_dissolved <- 1000 * central_sn38 / vc_sn38

    # Observed totals ($ERROR IND = 1, 2, 3)
    Cc <- Cc_np + Cc_dissolved
    Cc_sn38 <- Cc_sn38_np + Cc_sn38_dissolved
    Cc_sn38g <- 1000 * central_sn38g / vc_sn38g

    Cc ~ prop(propSd)
    Cc_sn38 ~ prop(propSd_sn38)
    Cc_sn38g ~ prop(propSd_sn38g)
  })
}
