vanDyk_2026_atazanavir <- function() {
  description <- paste(
    "Two-compartment population PK model for oral ritonavir-boosted",
    "atazanavir (ATV/r) in African children and adolescents with HIV on",
    "subsequent-line antiretroviral therapy (CHAPAS-4 trial and VirTUAL",
    "study; van Dyk 2026). Absorption is an absorption lag (1 h, fixed)",
    "followed by sequential zero-order input into a depot over D1 (1.04 h)",
    "and first-order absorption (ka 2.23 1/h); elimination is first order.",
    "Clearance, intercompartmental clearance and both volumes are",
    "allometrically scaled on total body weight with a 27 kg reference and",
    "fixed 0.75 / 1 exponents (CL 4.8 L/h, Vc 41.2 L, Q 0.418 L/h, Vp",
    "47.5 L at 27 kg). Atazanavir and ritonavir given as separate tablets",
    "(non-fixed-dose combination) have 26.7% lower bioavailability than the",
    "300/100 mg fixed-dose combination tablet. Rifampicin effects are NOT",
    "estimated from these data: they are the adult-derived simulation",
    "assumptions the paper fixed (clearance x3.12 and bioavailability -53%",
    "with once-daily ATV/r; clearance x2.04 with twice-daily ATV/r).",
    "Random effects are between-subject variability on clearance and",
    "between-occasion variability on lag time, D1, ka and bioavailability,",
    "the last inflated 1.43-fold for doses that were not directly observed.",
    "Residual error is combined proportional (11.8%) plus additive, the",
    "additive term fixed to 20% of the study assay's lower limit of",
    "quantification (0.021 mg/L CHAPAS-4, 0.006 mg/L VirTUAL)."
  )
  reference <- paste(
    "van Dyk J, Waitt C, Mugerwa H, Wiesner L, Orrell C, Bamford A,",
    "Gibb DM, McIlleron HM, Colbers A, Burger DM, Denti P, Wasmann RE (2026).",
    "Population pharmacokinetics of ritonavir-boosted atazanavir in",
    "subsequent-line treatment in African children with HIV.",
    "Antimicrob Agents Chemother 70(8):e00592-26. doi:10.1128/aac.00592-26",
    sep = " "
  )
  vignette <- "vanDyk_2026_atazanavir"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(
      analyte = "atazanavir",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "atazanavir",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "atazanavir",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric size descriptor with a 27 kg reference: exponent 0.75 on",
        "CL and Q, exponent 1 on Vc and Vp (van Dyk 2026 Methods, fixed a",
        "priori; Table 2 footnote a 'Allometrically scaled for a 27-kg",
        "individual'). Total body weight was chosen over fat-free mass,",
        "which fitted similarly (dOFV = 1.885)."
      ),
      source_name = "WT"
    ),
    FORM_FDC = list(
      description = "Fixed-dose-combination tablet indicator (1 = atazanavir/ritonavir 300/100 mg fixed-dose combination tablet, 0 = separate atazanavir and ritonavir tablets)",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (300/100 mg ATV/r fixed-dose combination; the typical-value reference of Table 2)",
      notes = paste(
        "The FDC arm is the bioavailability reference (F = 1 fixed); separate",
        "tablets enter as a deviation, F x (1 - 0.267) when FORM_FDC = 0.",
        "Table 2 footnote c: non-FDC covers children < 25 kg receiving",
        "2 x 100 mg atazanavir with 3 x 25 mg or 1 x 100 mg ritonavir tablets",
        "and children >= 25 kg receiving 150 or 200 mg atazanavir with a",
        "100 mg ritonavir tablet (VirTUAL). The paper tested ritonavir dose",
        "(75 vs 100 mg) against formulation and retained formulation; in",
        "CHAPAS-4 formulation is confounded with the < 25 kg weight bands."
      ),
      source_name = "formulation (non-FDC)"
    ),
    CONMED_RIF = list(
      description = "Concomitant rifampicin co-administration indicator (1 = on rifampicin at full induction, 0 = no rifampicin)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no rifampicin; all observed data)",
      notes = paste(
        "Simulation-only covariate. No atazanavir concentrations were",
        "collected with rifampicin (CHAPAS-4 stopped ATV/r during",
        "tuberculosis treatment), so the effects are fixed adult-derived",
        "assumptions, not estimates: van Dyk 2026 Methods 'Simulations' and",
        "Table 2 footnote d ('Rifampicin assumptions extracted from the Kengo",
        "et al. model'). Chronic-induction semantics: the indicator",
        "represents steady-state induction. Combined with REGI_BID it selects",
        "clearance x3.12 and bioavailability -53% with once-daily ATV/r, or",
        "clearance x2.04 and unchanged bioavailability with twice-daily ATV/r."
      ),
      source_name = "RIF (simulation scenario)"
    ),
    REGI_BID = list(
      description = "ATV/r dosing-regimen indicator (1 = twice daily, 0 = once daily)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (once-daily ATV/r)",
      notes = paste(
        "Interaction-only use, as in Kengo_2025_atazanavir: REGI_BID has no",
        "main effect and only switches which rifampicin induction factor",
        "applies (x2.04 on CL and no bioavailability penalty when twice",
        "daily, versus x3.12 on CL and -53% on F when once daily). All",
        "observed data were once daily; REGI_BID = 1 with CONMED_RIF = 0",
        "selects the reference parameters."
      ),
      source_name = "OD / BD (simulation scenario)"
    ),
    OCC = list(
      description = "Integer dosing-occasion index (1-6) multiplexing the between-occasion random effects",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "van Dyk 2026 Methods: 'each dose is considered a separate",
        "occasion'. The paper does not state the number of occasions; six",
        "are carried here, enough for the VirTUAL design (three visits, each",
        "with the self-administered dose before the visit and the observed",
        "clinic dose) and for CHAPAS-4 (the home dose behind the pre-dose",
        "sample and the observed PK-day dose). OCC values outside 1-6 get",
        "no between-occasion deviation. For simulation, assign OCC per",
        "dosing interval (cycling 1-6), or OCC = 1 throughout to make the",
        "occasion deviations act as subject-level variability."
      ),
      source_name = "OCC"
    ),
    SELFADMIN = list(
      description = "Self-administered (not directly observed) dosing-occasion indicator (1 = dose taken unobserved, 0 = directly observed)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (directly observed dose)",
      notes = paste(
        "Variance-model role: multiplies the between-occasion eta on",
        "bioavailability by the fixed scaling factor 1.43 for concentrations",
        "measured after an unobserved dose (van Dyk 2026 Results 'To account",
        "for high variability in pre-dose concentrations, a scaling factor",
        "was applied ... to inflate the between-occasion variability on",
        "bioavailability'; Table 2). Unlike Kengo 2025, only the",
        "bioavailability eta is scaled. Set SELFADMIN = 0 to simulate",
        "supervised dosing."
      ),
      source_name = "unobserved dose"
    ),
    STUDY_VIRTUAL = list(
      description = "VirTUAL study indicator (1 = VirTUAL, 0 = CHAPAS-4)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CHAPAS-4)",
      notes = paste(
        "Selects only the fixed additive residual error, which the paper",
        "fixed to 20% of each laboratory's lower limit of quantification",
        "(Table 2 footnote e: CHAPAS-4 LLOQ 0.105 mg/L, VirTUAL LLOQ",
        "0.0300 mg/L). No study effect on any structural parameter or on the",
        "proportional error was found. Has no effect on typical-value or",
        "IPRED simulations."
      ),
      source_name = "study"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Tested after allometric scaling; not significant (van Dyk 2026 Results)."
    ),
    FFM = list(
      description = "Fat-free mass",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Tested as the alternative allometric size descriptor; performed",
        "similarly to total body weight (dOFV = 1.885), and total body weight",
        "was retained for practicality."
      )
    ),
    CONMED_RTV_AUC = list(
      description = "Ritonavir AUC over the 0-24 h dosing interval",
      units = "mg*h/L",
      type = "continuous",
      notes = paste(
        "Power model on atazanavir clearance entered in forward inclusion",
        "(dOFV = -6.268, inverse relationship) but removed at backward",
        "elimination. Ritonavir dose (75 vs 100 mg) was also tested and lost",
        "to the formulation effect."
      )
    ),
    CONMED_TAF = list(
      description = "Concomitant emtricitabine/tenofovir alafenamide backbone",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "No significant effect of the TAF/emtricitabine, abacavir/lamivudine",
        "or zidovudine/lamivudine NRTI backbones on atazanavir PK (van Dyk",
        "2026 Results)."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 70L,
    n_studies = 2L,
    n_observations = 518L,
    age_range = "3.2-17.7 years (median 10.9)",
    age_median = "10.9 years",
    weight_range = "15-85 kg (median 29)",
    weight_median = "29 kg",
    height_range = "102-168 cm (median 133)",
    sex_female_pct = 44,
    race_ethnicity = "African (Zambia, Uganda, Zimbabwe, South Africa); race not tabulated",
    disease_state = paste(
      "Children and adolescents with HIV on subsequent-line antiretroviral",
      "therapy with ritonavir-boosted atazanavir. NRTI backbones: TAF/FTC 49%,",
      "ABC/3TC 31%, ZDV/3TC 13%, TDF/3TC 7.1%."
    ),
    dose_range = paste(
      "Once-daily ATV/r 150/100 mg (2.9%), 200/75 mg (34%), 200/100 mg",
      "(7.1%) as separate tablets, or 300/100 mg fixed-dose combination",
      "(56%), by WHO weight band (Table S1)."
    ),
    regions = "Zambia, Uganda, Zimbabwe (CHAPAS-4); South Africa, Uganda (VirTUAL)",
    notes = paste(
      "CHAPAS-4 (ISRCTN22964075): 60 children, intensive sampling at week 6",
      "(pre-dose, 1, 2, 4, 6, 8, 12, 24 h; plus 0.5 h on the TAF arm), 488",
      "concentrations, dosed with food under supervision. VirTUAL",
      "(NCT03923231): 10 children, sparse opportunistic sampling, 30",
      "concentrations. Two post-dose samples were below the LLOQ (M7+",
      "handling). Baseline characteristics are van Dyk 2026 Table 1."
    )
  )

  ini({
    # Final joint (CHAPAS-4 + VirTUAL) estimates, van Dyk 2026 Table 2, with
    # 95% CIs from sampling importance resampling. Disposition values are for
    # a 27 kg child on the 300/100 mg fixed-dose combination without
    # rifampicin. F is fixed to 1, so CL, Vc, Q and Vp are apparent oral values.
    lcl <- log(4.8)
    label("Apparent clearance CL/F at 27 kg (L/h)") # Table 2 CL 4.8 (4.37-5.41)
    lvc <- log(41.2)
    label("Apparent central volume Vc/F at 27 kg (L)") # Table 2 Vc 41.2 (37.8-46.0)
    lq <- log(0.418)
    label("Apparent intercompartmental clearance Q/F at 27 kg (L/h)") # Table 2 Q/F 0.418 (0.336-0.680)
    lvp <- log(47.5)
    label("Apparent peripheral volume Vp/F at 27 kg (L)") # Table 2 Vp/F 47.5 (18.9-113)
    ltlag <- fixed(log(1))
    label("Absorption lag time (h)") # Table 2 lag time 1 FIX (footnote g: estimated, then fixed for stability)
    ld1 <- log(1.04)
    label("Duration of zero-order input into the depot, D1 (h)") # Table 2 D1 1.04 (0.829-1.28)
    lka <- log(2.23)
    label("First-order absorption rate constant (1/h)") # Table 2 ka 2.23 (1.62-2.73)
    lfdepot <- fixed(log(1))
    label("Oral bioavailability of the fixed-dose combination (fraction)") # Table 2 F 1 FIX

    e_fdc0_fdepot <- -0.267
    label("Change in F with separate (non-FDC) ATV and RTV tablets (fraction)") # Table 2 non-FDC change in F -26.7% (-35.0 to -13.1)

    # Rifampicin induction, fixed adult-derived simulation assumptions
    # (Table 2 footnote d, Methods 'Simulations'); not estimated here.
    e_rif_qd_cl <- fixed(3.12)
    label("Fold-change in CL with rifampicin and once-daily ATV/r (-fold)") # Table 2 change in CL due to rifampicin (OD) 3.12 FIX
    e_rif_bid_cl <- fixed(2.04)
    label("Fold-change in CL with rifampicin and twice-daily ATV/r (-fold)") # Table 2 change in CL due to rifampicin (BD) 2.04 FIX
    e_rif_qd_fdepot <- fixed(-0.53)
    label("Change in F with rifampicin and once-daily ATV/r (fraction)") # Table 2 change in F due to rifampicin (OD) -53.0% FIX

    e_selfadmin_iov <- fixed(1.43)
    label("Multiplier on the between-occasion eta of F after an unobserved dose (-fold)") # Table 2 scaling factor on BOV for the unobserved dose 1.43 FIX

    e_wt_cl <- fixed(0.75)
    label("Allometric exponent of body weight on CL and Q (unitless)") # Methods: fixed exponent 0.75 on clearances
    e_wt_vc <- fixed(1)
    label("Allometric exponent of body weight on Vc and Vp (unitless)") # Methods: fixed exponent 1 on volumes

    # Table 2 footnote f: %CV = sqrt(exp(omega^2) - 1) * 100, so
    # omega^2 = log((CV/100)^2 + 1).
    etalcl ~ 0.0350978 # Table 2 BSV CL 18.9% CV; log(0.189^2 + 1)

    # Between-occasion variability: one shared variance per parameter across
    # the six occasion slots; slots 2-6 are pinned to slot 1.
    etaiov_tlag_1 ~ fixed(0.311905) # Table 2 BOV lag time 60.5% CV FIX; log(0.605^2 + 1)
    etaiov_tlag_2 ~ fixed(0.311905)
    etaiov_tlag_3 ~ fixed(0.311905)
    etaiov_tlag_4 ~ fixed(0.311905)
    etaiov_tlag_5 ~ fixed(0.311905)
    etaiov_tlag_6 ~ fixed(0.311905)
    etaiov_d1_1 ~ 0.802943 # Table 2 BOV D1 111% CV; log(1.11^2 + 1)
    etaiov_d1_2 ~ fixed(0.802943)
    etaiov_d1_3 ~ fixed(0.802943)
    etaiov_d1_4 ~ fixed(0.802943)
    etaiov_d1_5 ~ fixed(0.802943)
    etaiov_d1_6 ~ fixed(0.802943)
    etaiov_ka_1 ~ 1.28769 # Table 2 BOV ka 162% CV; log(1.62^2 + 1)
    etaiov_ka_2 ~ fixed(1.28769)
    etaiov_ka_3 ~ fixed(1.28769)
    etaiov_ka_4 ~ fixed(1.28769)
    etaiov_ka_5 ~ fixed(1.28769)
    etaiov_ka_6 ~ fixed(1.28769)
    etaiov_fdepot_1 ~ 0.113694 # Table 2 BOV F 34.7% CV; log(0.347^2 + 1)
    etaiov_fdepot_2 ~ fixed(0.113694)
    etaiov_fdepot_3 ~ fixed(0.113694)
    etaiov_fdepot_4 ~ fixed(0.113694)
    etaiov_fdepot_5 ~ fixed(0.113694)
    etaiov_fdepot_6 ~ fixed(0.113694)

    propSd <- 0.118
    label("Proportional residual error (fraction)") # Table 2 proportional error 11.8% (10.9-13.3)
    addSdChapas4 <- fixed(0.021)
    label("Additive residual error, CHAPAS-4 assay (mg/L)") # Table 2 footnote e: 20% of CHAPAS-4 LLOQ 0.105 mg/L
    addSdVirtual <- fixed(0.006)
    label("Additive residual error, VirTUAL assay (mg/L)") # Table 2 footnote e: 20% of VirTUAL LLOQ 0.0300 mg/L
  })

  model({
    # Occasion multiplexing of the between-occasion etas.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    oc6 <- (OCC == 6)
    iov_tlag <- oc1 * etaiov_tlag_1 + oc2 * etaiov_tlag_2 + oc3 * etaiov_tlag_3 +
      oc4 * etaiov_tlag_4 + oc5 * etaiov_tlag_5 + oc6 * etaiov_tlag_6
    iov_d1 <- oc1 * etaiov_d1_1 + oc2 * etaiov_d1_2 + oc3 * etaiov_d1_3 +
      oc4 * etaiov_d1_4 + oc5 * etaiov_d1_5 + oc6 * etaiov_d1_6
    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2 + oc3 * etaiov_ka_3 +
      oc4 * etaiov_ka_4 + oc5 * etaiov_ka_5 + oc6 * etaiov_ka_6
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2 +
      oc3 * etaiov_fdepot_3 + oc4 * etaiov_fdepot_4 +
      oc5 * etaiov_fdepot_5 + oc6 * etaiov_fdepot_6

    # An unobserved dose inflates only the bioavailability occasion eta.
    iov_scale_f <- 1 + (e_selfadmin_iov - 1) * SELFADMIN

    # Rifampicin regimen factors; all collapse to 1 when CONMED_RIF = 0.
    rif_qd <- CONMED_RIF * (1 - REGI_BID)
    rif_bid <- CONMED_RIF * REGI_BID
    cl_rif <- 1 + (e_rif_qd_cl - 1) * rif_qd + (e_rif_bid_cl - 1) * rif_bid
    fdepot_rif <- 1 + e_rif_qd_fdepot * rif_qd

    # Separate tablets lower bioavailability relative to the FDC reference.
    fdepot_form <- 1 + e_fdc0_fdepot * (1 - FORM_FDC)

    cl <- exp(lcl + etalcl) * (WT / 27)^e_wt_cl * cl_rif
    vc <- exp(lvc) * (WT / 27)^e_wt_vc
    q <- exp(lq) * (WT / 27)^e_wt_cl
    vp <- exp(lvp) * (WT / 27)^e_wt_vc
    tlag <- exp(ltlag + iov_tlag)
    d1 <- exp(ld1 + iov_d1)
    ka <- exp(lka + iov_ka)
    fdepot <- exp(lfdepot + iov_fdepot * iov_scale_f) * fdepot_form * fdepot_rif

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Doses enter the depot after the lag as a zero-order input lasting D1
    # (dose records need rate = -2 so rxode2 uses the modelled duration).
    alag(depot) <- tlag
    dur(depot) <- d1
    f(depot) <- fdepot

    Cc <- central / vc

    addSdLLOQ <- addSdChapas4 * (1 - STUDY_VIRTUAL) + addSdVirtual * STUDY_VIRTUAL
    Cc ~ add(addSdLLOQ) + prop(propSd)
  })
}
