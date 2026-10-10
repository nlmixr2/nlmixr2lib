Simons_2022_s_ketamine <- function() {
  description <- "Joint population PK model of S-ketamine (three-compartment), S-norketamine (two-compartment) and S-hydroxynorketamine (two-compartment) after a sublingual or buccal S-ketamine oral thin film (50 or 100 mg) followed by a 20 mg intravenous S-ketamine infusion in healthy adult volunteers (Simons 2022). The film dose is absorbed by two parallel routes: an oral-mucosal route (bioavailability F1 = 26.3%, zero-order input over D1 = 13.1 min into a depot drained first order by KA1 into the S-ketamine central compartment) and a swallowed route (F2 = 116% of the film dose, zero-order input over D2 = 29.9 min into a gut depot drained by KA2, then a gut-to-liver delay compartment with mean transit time MTTG) that delivers S-ketamine directly into hepatic metabolism without reaching the systemic S-ketamine pool. Both the swallowed S-ketamine and systemic S-ketamine cleared by CLK1 pass through a two-compartment metabolism delay chain (mean transit time 20.1 min) of which a fixed 80% forms S-norketamine; S-norketamine cleared by CLN1 passes through a second two-compartment chain (1.12 min) of which a fixed 70% forms S-hydroxynorketamine. The S-norketamine central volume equals the S-ketamine central volume. Clearances and volumes are referenced to 70 kg. Every OTF dose needs TWO dose records of the full film dose with rate = -2, one into `depot` and one into `depot2`; concentrations are in nmol/mL."
  reference <- paste(
    "Simons P, Olofsen E, van Velzen M, van Lemmen M, Mooren R, van Dasselaar T,",
    "Mohr P, Hammes F, van der Schrier R, Niesters M, Dahan A.",
    "S-Ketamine Oral Thin Film-Part 1: Population Pharmacokinetics of S-Ketamine,",
    "S-Norketamine and S-Hydroxynorketamine.",
    "Front Pain Res. 2022;3:946486. doi:10.3389/fpain.2022.946486",
    sep = " "
  )
  vignette <- "Simons_2022_s_ketamine"
  units <- list(time = "min", dosing = "mg", concentration = "nmol/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight (kg).",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Simons 2022 Table 4 reports every clearance and volume '@ 70 kg' but does not print the scaling exponents. The standard allometric exponents (0.75 for clearances and intercompartmental clearances, 1 for volumes) are used here, fixed; see the vignette Assumptions section.",
      source_name = "WT"
    ),
    OCC = list(
      description = "Study occasion (visit) index for inter-occasion variability: 1 or 2.",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = "Simons 2022 is a two-visit crossover (one visit with a 50 mg film, one with a 100 mg film, in random order, at least 7 days apart); the occasion is the visit. Inter-occasion variability is applied to F1, D1, KA1, F2, D2, KA2, MTTG, MTT K->NK, VH2 and CLH1. Set OCC = 1 or 2 per visit; a value outside 1-2 switches every IOV term off.",
      source_name = "OCC"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "S-ketamine", units = "mg", specimen = "administration site", verified = TRUE),
    depot2 = list(analyte = "S-ketamine", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "S-ketamine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "S-ketamine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "S-ketamine", units = "mg", specimen = "tissue", verified = TRUE),
    peripheral2 = list(analyte = "S-ketamine", units = "mg", specimen = "tissue", verified = TRUE),
    transit1_snk = list(analyte = "S-ketamine", units = "mg", specimen = "not applicable", verified = TRUE),
    transit2_snk = list(analyte = "S-ketamine", units = "mg", specimen = "not applicable", verified = TRUE),
    central_snk = list(
      analyte = "S-norketamine",
      units = "mg S-ketamine equivalents",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1_snk = list(
      analyte = "S-norketamine",
      units = "mg S-ketamine equivalents",
      specimen = "tissue",
      verified = TRUE
    ),
    transit1_shnk = list(
      analyte = "S-norketamine",
      units = "mg S-ketamine equivalents",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit2_shnk = list(
      analyte = "S-norketamine",
      units = "mg S-ketamine equivalents",
      specimen = "not applicable",
      verified = TRUE
    ),
    central_shnk = list(
      analyte = "S-hydroxynorketamine",
      units = "mg S-ketamine equivalents",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1_shnk = list(
      analyte = "S-hydroxynorketamine",
      units = "mg S-ketamine equivalents",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 20L,
    n_studies = 1L,
    age_range = "19-32 years",
    age_mean = "24 years (SD 3)",
    weight_range = "53-93 kg",
    weight_mean = "73 kg (SD 12)",
    height_mean = "179 cm (SD 10; range 161-197)",
    bmi_mean = "23 kg/m^2 (SD 2; range 19-27)",
    sex_female_pct = 50,
    disease_state = "Healthy volunteers.",
    dose_range = "S-ketamine oral thin film 50 mg (one film) and 100 mg (two films) on separate visits in random order, placed sublingually (n = 15) or buccally (n = 5) and held for 10 min without swallowing; 6 h after film placement every subject received 20 mg S-ketamine intravenously over 20 min.",
    regions = "Netherlands (Leiden University Medical Center).",
    n_observations = "Arterial plasma S-ketamine, S-norketamine and S-hydroxynorketamine at 0, 5, 10, 20, 40, 60, 90, 120, 180, 240, 300 and 360 min after film placement and at 2, 4, 10, 15, 20, 30, 40, 60, 75, 90 and 120 min after the start of the intravenous infusion, on each visit.",
    notes = "Demographics from Simons 2022 Table 1. 19 of 20 subjects completed both visits; one subject withdrew after the first (100 mg) visit. Sublingual and buccal placement were pooled because no PK difference was observed. Concentrations were converted from ng/mL to nmol/mL before modelling."
  )

  ini({
    # Oral-mucosal absorption from the film (Simons 2022 Table 4)
    lfdepot <- log(0.263); label("Oral-mucosal bioavailability of the film dose, F1 (fraction)") # Table 4: F1 = 26.3%, SEE 1.2
    ld1 <- log(13.1); label("Duration of zero-order mucosal input, D1 (min)") # Table 4: D1 = 13.1 min, SEE 1.0
    lka <- log(0.04); label("Mucosal absorption rate constant, KA1 (1/min)") # Table 4: KA1 = 0.04 1/min, SEE 0.002

    # Swallowed (gastrointestinal) route (Simons 2022 Table 4)
    lfdepot2 <- log(1.16); label("Gastrointestinal bioavailability of the film dose into hepatic metabolism, F2 (fraction)") # Table 4: F2 = 116%, SEE 6
    ld2 <- log(29.9); label("Duration of zero-order gastrointestinal input, D2 (min)") # Table 4: D2 = 29.9 min, SEE 3.5
    lka2 <- log(0.049); label("Gastrointestinal absorption rate constant, KA2 (1/min)") # Table 4: KA2 = 0.049 1/min, SEE 0.007
    lmtt <- log(10.7); label("Gut-to-liver mean transit time, MTTG (min)") # Table 4: MTT GUT = 10.7 min, SEE 1.7

    # S-ketamine disposition, 70 kg (Simons 2022 Table 4)
    lvc <- log(11.6); label("S-ketamine central volume VK1 at 70 kg (L)") # Table 4: VK1 = 11.6 L, SEE 0.9
    lvp <- log(39.0); label("S-ketamine first peripheral volume VK2 at 70 kg (L)") # Table 4: VK2 = 39.0 L, SEE 2.9
    lvp2 <- log(174); label("S-ketamine second peripheral volume VK3 at 70 kg (L)") # Table 4: VK3 = 174 L, SEE 11
    lcl <- log(1.48); label("S-ketamine clearance into the metabolism chain CLK1 at 70 kg (L/min)") # Table 4: CLK1 = 1.48 L/min, SEE 0.06
    lq <- log(2.43); label("S-ketamine intercompartmental clearance VK1-VK2, CLK2, at 70 kg (L/min)") # Table 4: CLK2 = 2.43 L/min, SEE 0.24
    lq2 <- log(1.21); label("S-ketamine intercompartmental clearance VK1-VK3, CLK3, at 70 kg (L/min)") # Table 4: CLK3 = 1.21 L/min, SEE 0.08

    # S-ketamine -> S-norketamine metabolism (Simons 2022 Table 4, Methods, Figure 2)
    lmtt_snk <- log(20.1); label("Mean transit time of the two-compartment S-ketamine metabolism chain, MTT K->NK (min)") # Table 4: MTT K->NK = 20.1 min, SEE 1.0
    fm_snk <- fixed(0.8); label("Fraction of metabolised S-ketamine forming S-norketamine (fraction)") # Methods 'Population Pharmacokinetic Analysis': 80% assumed; Figure 2 '20% loss'

    # S-norketamine disposition, 70 kg (VN1 = VK1; Simons 2022 Table 4, Figure 2)
    lvp_snk <- log(221); label("S-norketamine peripheral volume VN2 at 70 kg (L)") # Table 4: VN2 = 221 L, SEE 13
    lcl_snk <- log(1.00); label("S-norketamine clearance into the metabolism chain CLN1 at 70 kg (L/min)") # Table 4: CLN1 = 1.00 L/min, SEE 0.04
    lq_snk <- log(2.63); label("S-norketamine intercompartmental clearance CLN2 at 70 kg (L/min)") # Table 4: CLN2 = 2.63 L/min, SEE 0.15

    # S-norketamine -> S-hydroxynorketamine metabolism (Simons 2022 Table 4, Methods, Figure 2)
    lmtt_shnk <- log(1.12); label("Mean transit time of the two-compartment S-norketamine metabolism chain, MTT NK->HNK (min)") # Table 4: MTT NK->HNK = 1.12 min, SEE 0.51
    fm_shnk <- fixed(0.7); label("Fraction of metabolised S-norketamine forming S-hydroxynorketamine (fraction)") # Methods 'Population Pharmacokinetic Analysis': 70% assumed; Figure 2 '30% loss'

    # S-hydroxynorketamine disposition, 70 kg (Simons 2022 Table 4)
    lvc_shnk <- log(4.4); label("S-hydroxynorketamine central volume VH1 at 70 kg (L)") # Table 4: VH1 = 4.4 L, SEE 2.0
    lvp_shnk <- log(87.5); label("S-hydroxynorketamine peripheral volume VH2 at 70 kg (L)") # Table 4: VH2 = 87.5 L, SEE 6.5
    lcl_shnk <- log(0.933); label("S-hydroxynorketamine terminal clearance CLH1 at 70 kg (L/min)") # Table 4: CLH1 = 0.933 L/min, SEE 0.068
    lq_shnk <- log(1.70); label("S-hydroxynorketamine intercompartmental clearance CLH2 at 70 kg (L/min)") # Table 4: CLH2 = 1.70 L/min, SEE 0.25

    # Allometric exponents: Table 4 values are '@ 70 kg'; the exponents are not printed
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent of (WT/70) on all clearances (unitless)") # not printed; standard allometric value (vignette Assumptions)
    e_wt_vc_vp <- fixed(1); label("Allometric exponent of (WT/70) on all volumes (unitless)") # not printed; standard allometric value (vignette Assumptions)

    # Inter-subject variability, omega^2 (Simons 2022 Table 4)
    etalvc ~ 0.057 # Table 4 VK1 (and VN1, shared) omega^2 = 0.057, SEE 0.019
    etalcl ~ 0.029 # Table 4 CLK1 omega^2 = 0.029, SEE 0.012
    etalq2 ~ 0.026 # Table 4 CLK3 omega^2 = 0.026, SEE 0.014
    etalcl_snk ~ 0.050 # Table 4 CLN1 omega^2 = 0.050, SEE 0.012
    etalmtt_snk ~ 0.021 # Table 4 MTT K->NK omega^2 = 0.021, SEE 0.122
    etalvc_shnk ~ 1.22 # Table 4 VH1 omega^2 = 1.22, SEE 0.91
    etalcl_shnk ~ 0.103 # Table 4 CLH1 omega^2 = 0.103, SEE 0.042
    etalq_shnk ~ 0.287 # Table 4 CLH2 omega^2 = 0.287, SEE 0.124

    # Inter-occasion variability, nu^2, two occasions (Simons 2022 Table 4)
    etaiov_fdepot_1 ~ 0.060 # Table 4 F1 nu^2 = 0.060, SEE 0.019
    etaiov_fdepot_2 ~ fixed(0.060) # same variance on occasion 2
    etaiov_d1_1 ~ 0.154 # Table 4 D1 nu^2 = 0.154, SEE 0.033
    etaiov_d1_2 ~ fixed(0.154) # same variance on occasion 2
    etaiov_ka_1 ~ 0.062 # Table 4 KA1 nu^2 = 0.062, SEE 0.014
    etaiov_ka_2 ~ fixed(0.062) # same variance on occasion 2
    etaiov_fdepot2_1 ~ 0.057 # Table 4 F2 nu^2 = 0.057, SEE 0.031
    etaiov_fdepot2_2 ~ fixed(0.057) # same variance on occasion 2
    etaiov_d2_1 ~ 0.611 # Table 4 D2 nu^2 = 0.611, SEE 0.120
    etaiov_d2_2 ~ fixed(0.611) # same variance on occasion 2
    etaiov_ka2_1 ~ 0.376 # Table 4 KA2 nu^2 = 0.376, SEE 0.150
    etaiov_ka2_2 ~ fixed(0.376) # same variance on occasion 2
    etaiov_mtt_1 ~ 0.937 # Table 4 MTT GUT nu^2 = 0.937, SEE 0.312
    etaiov_mtt_2 ~ fixed(0.937) # same variance on occasion 2
    etaiov_mtt_snk_1 ~ 0.751 # Table 4 nu^2 = 0.751, SEE 0.349, printed on the S-norketamine additive-error row; assigned to MTT K->NK (vignette Assumptions)
    etaiov_mtt_snk_2 ~ fixed(0.751) # same variance on occasion 2
    etaiov_vp_shnk_1 ~ 0.152 # Table 4 VH2 nu^2 = 0.152, SEE 0.031
    etaiov_vp_shnk_2 ~ fixed(0.152) # same variance on occasion 2
    etaiov_cl_shnk_1 ~ 0.008 # Table 4 CLH1 nu^2 = 0.008, SEE 0.004
    etaiov_cl_shnk_2 ~ fixed(0.008) # same variance on occasion 2

    # Residual error, SD scale (Simons 2022 Table 4; concentrations in nmol/mL)
    propSd <- 0.1095; label("S-ketamine proportional residual SD (fraction)") # Table 4: sigma Relative = 0.012, SEE 0.0004, read as a variance (sqrt = 0.1095); see vignette Assumptions
    propSd_snk <- 0.102; label("S-norketamine proportional residual SD (fraction)") # Table 4: sigma Relative = 0.102, SEE 0.007
    addSd_snk <- 0.058; label("S-norketamine additive residual SD (nmol/mL)") # Table 4: sigma Additive = 0.058, SEE 0.018
    propSd_shnk <- 0.079; label("S-hydroxynorketamine proportional residual SD (fraction)") # Table 4: sigma Relative = 0.079, SEE 0.005
    addSd_shnk <- 0.020; label("S-hydroxynorketamine additive residual SD (nmol/mL)") # Table 4: sigma Additive = 0.020, SEE 0.003
  })

  model({
    # Molecular weight of S-ketamine free base (g/mol). Every state is carried
    # in mg S-ketamine (equivalents) so that metabolite formation conserves
    # moles; dividing by this weight gives umol/L = nmol/mL for all analytes.
    mw_ketamine <- 237.73

    # Occasion indicators for inter-occasion variability
    oc1 <- OCC == 1
    oc2 <- OCC == 2
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2
    iov_d1 <- oc1 * etaiov_d1_1 + oc2 * etaiov_d1_2
    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2
    iov_fdepot2 <- oc1 * etaiov_fdepot2_1 + oc2 * etaiov_fdepot2_2
    iov_d2 <- oc1 * etaiov_d2_1 + oc2 * etaiov_d2_2
    iov_ka2 <- oc1 * etaiov_ka2_1 + oc2 * etaiov_ka2_2
    iov_mtt <- oc1 * etaiov_mtt_1 + oc2 * etaiov_mtt_2
    iov_mtt_snk <- oc1 * etaiov_mtt_snk_1 + oc2 * etaiov_mtt_snk_2
    iov_vp_shnk <- oc1 * etaiov_vp_shnk_1 + oc2 * etaiov_vp_shnk_2
    iov_cl_shnk <- oc1 * etaiov_cl_shnk_1 + oc2 * etaiov_cl_shnk_2

    # Allometric multipliers (reference 70 kg)
    wt_cl <- (WT / 70)^e_wt_cl_q
    wt_v <- (WT / 70)^e_wt_vc_vp

    # Absorption
    fdepot <- exp(lfdepot + iov_fdepot)
    d1 <- exp(ld1 + iov_d1)
    ka <- exp(lka + iov_ka)
    fdepot2 <- exp(lfdepot2 + iov_fdepot2)
    d2 <- exp(ld2 + iov_d2)
    ka2 <- exp(lka2 + iov_ka2)
    mtt <- exp(lmtt + iov_mtt)

    # S-ketamine
    vc <- exp(lvc + etalvc) * wt_v
    vp <- exp(lvp) * wt_v
    vp2 <- exp(lvp2) * wt_v
    cl <- exp(lcl + etalcl) * wt_cl
    q <- exp(lq) * wt_cl
    q2 <- exp(lq2 + etalq2) * wt_cl

    # S-norketamine (central volume shared with S-ketamine, Figure 2 'VN1 = VK1')
    mtt_snk <- exp(lmtt_snk + etalmtt_snk + iov_mtt_snk)
    vc_snk <- vc
    vp_snk <- exp(lvp_snk) * wt_v
    cl_snk <- exp(lcl_snk + etalcl_snk) * wt_cl
    q_snk <- exp(lq_snk) * wt_cl

    # S-hydroxynorketamine
    mtt_shnk <- exp(lmtt_shnk)
    vc_shnk <- exp(lvc_shnk + etalvc_shnk) * wt_v
    vp_shnk <- exp(lvp_shnk + iov_vp_shnk) * wt_v
    cl_shnk <- exp(lcl_shnk + etalcl_shnk + iov_cl_shnk) * wt_cl
    q_shnk <- exp(lq_shnk + etalq_shnk) * wt_cl

    # Rate constants. Each two-compartment metabolism chain has a total mean
    # transit time MTT, so each of its compartments drains at 2 / MTT.
    ktr <- 1 / mtt
    ktr_snk <- 2 / mtt_snk
    ktr_shnk <- 2 / mtt_shnk
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2
    kel_snk <- cl_snk / vc_snk
    k12_snk <- q_snk / vc_snk
    k21_snk <- q_snk / vp_snk
    kel_shnk <- cl_shnk / vc_shnk
    k12_shnk <- q_shnk / vc_shnk
    k21_shnk <- q_shnk / vp_shnk

    # Oral-mucosal route into systemic S-ketamine
    d/dt(depot) <- -ka * depot
    # Swallowed route: gut depot -> gut-to-liver delay -> hepatic metabolism
    d/dt(depot2) <- -ka2 * depot2
    d/dt(transit1) <- ka2 * depot2 - ktr * transit1

    d/dt(central) <- ka * depot - (kel + k12 + k13) * central + k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # S-ketamine metabolism chain, fed by systemic clearance and the swallowed route
    d/dt(transit1_snk) <- kel * central + ktr * transit1 - ktr_snk * transit1_snk
    d/dt(transit2_snk) <- ktr_snk * transit1_snk - ktr_snk * transit2_snk

    d/dt(central_snk) <- fm_snk * ktr_snk * transit2_snk - (kel_snk + k12_snk) * central_snk + k21_snk * peripheral1_snk
    d/dt(peripheral1_snk) <- k12_snk * central_snk - k21_snk * peripheral1_snk

    # S-norketamine metabolism chain
    d/dt(transit1_shnk) <- kel_snk * central_snk - ktr_shnk * transit1_shnk
    d/dt(transit2_shnk) <- ktr_shnk * transit1_shnk - ktr_shnk * transit2_shnk

    d/dt(central_shnk) <- fm_shnk * ktr_shnk * transit2_shnk - (kel_shnk + k12_shnk) * central_shnk + k21_shnk * peripheral1_shnk
    d/dt(peripheral1_shnk) <- k12_shnk * central_shnk - k21_shnk * peripheral1_shnk

    # Bioavailability and zero-order input durations (dose records need rate = -2)
    f(depot) <- fdepot
    dur(depot) <- d1
    f(depot2) <- fdepot2
    dur(depot2) <- d2

    # Concentrations in nmol/mL (mg/L divided by g/mol gives mmol/L; x 1000 = umol/L = nmol/mL)
    Cc <- 1000 * central / vc / mw_ketamine
    Cc_snk <- 1000 * central_snk / vc_snk / mw_ketamine
    Cc_shnk <- 1000 * central_shnk / vc_shnk / mw_ketamine

    Cc ~ prop(propSd)
    Cc_snk ~ add(addSd_snk) + prop(propSd_snk)
    Cc_shnk ~ add(addSd_shnk) + prop(propSd_shnk)
  })
}
