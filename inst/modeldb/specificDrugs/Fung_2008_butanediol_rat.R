Fung_2008_butanediol_rat <- function() {
  description <- paste(
    "Preclinical (rat).",
    "Population PK model for 1,4-butanediol (BD) and its bioactivation",
    "pathway in adult male Sprague-Dawley rats after intravenous and",
    "oral dosing, jointly with the unmeasured semialdehyde intermediate",
    "(ALD), the gamma-hydroxybutyric acid metabolite (GHB), and",
    "co-administered ethanol (ETOH). Each of the four substances follows",
    "two-compartment disposition with Michaelis-Menten (mixed-order)",
    "elimination; GHB additionally has a parallel first-order elimination.",
    "The metabolic flux is BD -> ALD -> GHB; ETOH has no metabolic",
    "connection to the BD/ALD/GHB chain. Mutual competitive inhibition",
    "is encoded for BD elimination inhibited by GHB and by ETOH; ALD",
    "elimination inhibited by BD; GHB elimination inhibited by BD; ETOH",
    "elimination inhibited by BD. Oral BD absorption is first-order",
    "with a 7.5 min lag-time; the absorbed fraction (F = 0.93,",
    "dose-independent) is split between BD and ALD central compartments",
    "to represent pre-systemic conversion of BD to ALD. The fraction",
    "entering as BD was dose-dependent in the source paper (30 percent",
    "at 1.58 mmol/kg and 55 percent at 6.34 mmol/kg); the model defaults",
    "to the low-dose 1.58 mmol/kg fractions (30 percent as BD, 70",
    "percent as ALD); see the vignette Assumptions and deviations for",
    "the high-dose alternative. Volume of distribution of ALD is",
    "mathematically non-identifiable and was constrained so that Vss",
    "of ALD equals Vss of BD (source paper footnote c). Volume of the",
    "central compartment for GHB was fixed at the average rat plasma",
    "volume of 0.010 L (source paper footnote b). Fit by NONMEM VI",
    "ADVAN9 with FOCE-I."
  )
  reference <- paste(
    "Fung HL, Tsou PS, Bulitta JB, Tran DC, Page NA, Soda D, Fung SM.",
    "Pharmacokinetics of 1,4-butanediol in rats: bioactivation to",
    "gamma-hydroxybutyric acid, interaction with ethanol, and oral",
    "bioavailability.",
    "The AAPS Journal. 2008;10(1):56-69.",
    "doi:10.1208/s12248-007-9006-3.",
    sep = " "
  )
  vignette <- "Fung_2008_butanediol_rat"
  paper_specific_compartments <- c(
    "central_ald", "peripheral1_ald",
    "central_ghb", "peripheral1_ghb",
    "central_etoh", "peripheral1_etoh"
  )
  paper_specific_etas <- c(
    "etalkm_ghb", "etalvc_ghb",
    "etalvmax_etoh", "etalkm_etoh", "etalvc_etoh",
    "etalki_bd_ald", "etalki_bd_ghb", "etalki_bd_etoh", "etalki_ghb_bd"
  )
  paper_specific_residual_sds <- c(
    "propSd_ghb", "addSd_ghb",
    "propSd_etoh", "addSd_etoh"
  )

  units <- list(time = "min", dosing = "mmol", concentration = "mmol/L")

  covariateData <- list()

  population <- list(
    species        = "rat (male Sprague-Dawley)",
    n_subjects     = 42L,
    n_studies      = 13L,
    age_range      = "adult (exact age not reported)",
    weight_range   = "~ 300 g (uniform within Harlan-supplied cohort)",
    sex_female_pct = 0,
    disease_state  = paste(
      "Healthy male Sprague-Dawley rats (Harlan, Indianapolis, IN) of",
      "about 300 g body weight. Jugular vein cannulated for blood",
      "withdrawal and femoral vein cannulated for IV dosing."
    ),
    dose_range     = paste(
      "Thirteen studies (Table I): 11 IV studies covering BD 1.58 or",
      "6.34 mmol/kg, GHB 1.58, 1.79, or 6.34 mmol/kg, and ETOH 6.34 or",
      "12.7 mmol/kg in mono- and pairwise co-administration; plus 2 oral",
      "BD studies (1.58 or 6.34 mmol/kg by gastric gavage). Each IV dose",
      "was a 5-min constant-rate infusion. Groups of 3 or 4 rats per study."
    ),
    regions        = "USA (University at Buffalo)",
    notes          = paste(
      "Animal protocols reviewed and approved by the University at",
      "Buffalo Institutional Animal Care and Use Committee. Plasma",
      "samples assayed for BD and GHB by LCMS (LOQ 60 uM) and for",
      "ETOH by GC-FID. The semialdehyde intermediate ALD was not",
      "measured experimentally; its parameters were inferred indirectly,",
      "with V_ALD constrained so that Vss_ALD equals Vss_BD (source",
      "paper footnote c). Model fit by NONMEM VI ADVAN9 with the FOCE",
      "interaction option. See Fung 2008 Methods (Animal experiments;",
      "Population PK modeling) and Tables I and II."
    )
  )

  ini({
    # ============================================================
    # 1,4-Butanediol (BD) -- parent drug; canonical (unsuffixed)
    # parameter names. Source: Fung 2008 Table II, BD column.
    # Vmax_BD is "calculated" (= CLic_BD * Km_BD, Table II
    # footnote a), so the estimated primary parameters were
    # CLic_BD (with IIV 11 percent) and Km_BD (no IIV). Because
    # Km_BD has no IIV, IIV on CLic and IIV on Vmax are
    # mathematically equivalent (Vmax = CLic * Km_typ is a
    # constant scaling of CLic); the 11 percent IIV is therefore
    # carried on lvmax here without information loss.
    # ============================================================
    lvmax <- log(0.0362)  ; label("BD Michaelis-Menten Vmax (mmol/min)")        # Table II BD: Vmax = 0.0362 (calc as CLic 0.0232 * Km 1.56)
    lkm   <- log(1.56)    ; label("BD Michaelis-Menten Km (mmol/L)")            # Table II BD: Km = 1.56
    lvc   <- log(0.248)   ; label("BD central volume of distribution V (L)")     # Table II BD: V = 0.248
    lk12  <- log(0.00231) ; label("BD inter-compartmental rate constant k12 (1/min)")  # Table II BD: K12 = 0.00231
    lk21  <- log(0.0150)  ; label("BD inter-compartmental rate constant k21 (1/min)")  # Table II BD: K21 = 0.0150

    # ============================================================
    # Semialdehyde intermediate (ALD) -- unmeasured species.
    # Vmax_ALD was the estimated primary parameter; CLic_ALD is
    # calculated (Table II footnote a). No IIV on any ALD
    # parameter. The central volume V_ALD is structurally
    # constrained so that Vss_ALD = V_ALD * (1 + K12/K21) equals
    # Vss_BD = 0.286 L (Table II footnote c); the numeric
    # K12_ALD / K21_ALD / V_ALD set below satisfies this
    # constraint exactly (0.0476 * (1 + 0.0526 / 0.0105) =
    # 0.286 L).
    # ============================================================
    lvmax_ald <- log(0.0168) ; label("ALD Michaelis-Menten Vmax (mmol/min)")           # Table II ALD: Vmax = 0.0168 (estimated)
    lkm_ald   <- log(0.446)  ; label("ALD Michaelis-Menten Km (mmol/L)")               # Table II ALD: Km = 0.446
    lvc_ald   <- log(0.0476) ; label("ALD central volume of distribution V (L)")        # Table II ALD: V = 0.0476 (constrained s.t. Vss_ALD = Vss_BD; footnote c)
    lk12_ald  <- log(0.0526) ; label("ALD inter-compartmental rate constant k12 (1/min)")  # Table II ALD: K12 = 0.0526
    lk21_ald  <- log(0.0105) ; label("ALD inter-compartmental rate constant k21 (1/min)")  # Table II ALD: K21 = 0.0105

    # ============================================================
    # gamma-Hydroxybutyric acid (GHB) -- downstream metabolite.
    # Estimated primary parameters were CL (linear, no IIV),
    # CLic (no IIV listed), Km (IIV 35 percent), V (fixed to
    # rat plasma volume 0.010 L; Table II footnote b), V (IIV
    # 16 percent on log of fixed value -> NONMEM IIV applied
    # to an effectively constant typical value; included here
    # as reported), K12 (no IIV), K21 (no IIV). Vmax_GHB is
    # calculated (Table II footnote a) and shown without IIV.
    # ============================================================
    lcl_ghb   <- log(0.00046) ; label("GHB first-order linear clearance CL (L/min)")        # Table II GHB: CL = 0.00046
    lvmax_ghb <- log(0.00361) ; label("GHB Michaelis-Menten Vmax (mmol/min)")               # Table II GHB: Vmax = 0.00361 (calc as CLic 0.0398 * Km 0.0906)
    lkm_ghb   <- log(0.0906)  ; label("GHB Michaelis-Menten Km (mmol/L)")                   # Table II GHB: Km = 0.0906
    lvc_ghb   <- fixed(log(0.010)) ; label("GHB central volume V fixed at rat plasma volume (L)") # Table II GHB footnote b: V FIXED to 0.010 L
    lk12_ghb  <- log(0.551)   ; label("GHB inter-compartmental rate constant k12 (1/min)")  # Table II GHB: K12 = 0.551
    lk21_ghb  <- log(0.0554)  ; label("GHB inter-compartmental rate constant k21 (1/min)")  # Table II GHB: K21 = 0.0554

    # ============================================================
    # Ethanol (ETOH) -- co-administered substance with its own
    # two-compartment MM disposition. CLic and Km both estimated
    # (CLic IIV 1.9 percent, Km IIV 14 percent); V estimated
    # (IIV 20 percent); K12 and K21 estimated without IIV.
    # ============================================================
    lvmax_etoh <- log(0.0410) ; label("ETOH Michaelis-Menten Vmax (mmol/min)")              # Table II ETOH: Vmax = 0.0410 (calc as CLic 0.0550 * Km 0.746)
    lkm_etoh   <- log(0.746)  ; label("ETOH Michaelis-Menten Km (mmol/L)")                  # Table II ETOH: Km = 0.746
    lvc_etoh   <- log(0.205)  ; label("ETOH central volume of distribution V (L)")           # Table II ETOH: V = 0.205
    lk12_etoh  <- log(0.0338) ; label("ETOH inter-compartmental rate constant k12 (1/min)")  # Table II ETOH: K12 = 0.0338
    lk21_etoh  <- log(0.0589) ; label("ETOH inter-compartmental rate constant k21 (1/min)")  # Table II ETOH: K21 = 0.0589

    # ============================================================
    # Mutual competitive inhibition constants. Each Ki appears in
    # the denominator of the inhibited substance's Michaelis-
    # Menten rate as Km_target * (1 + C_inhibitor / Ki_inhib/target).
    # Source: Fung 2008 Table II inhibition-constants section
    # and Eqs. 1-4.
    # ============================================================
    lki_bd_ald  <- log(3.67)  ; label("Ki for BD inhibiting ALD elimination (mmol/L)")  # Table II: Ki_BD/ALD = 3.67
    lki_bd_ghb  <- log(3.56)  ; label("Ki for BD inhibiting GHB elimination (mmol/L)")  # Table II: Ki_BD/GHB = 3.56
    lki_bd_etoh <- log(2.24)  ; label("Ki for BD inhibiting ETOH elimination (mmol/L)") # Table II: Ki_BD/ETOH = 2.24
    lki_ghb_bd  <- log(15.2)  ; label("Ki for GHB inhibiting BD elimination (mmol/L)")  # Table II: Ki_GHB/BD = 15.2
    lki_etoh_bd <- log(0.615) ; label("Ki for ETOH inhibiting BD elimination (mmol/L)") # Table II: Ki_ETOH/BD = 0.615

    # ============================================================
    # Oral BD absorption (Fung 2008 Results "BD Oral
    # Bioavailability"). First-order absorption with a 7.5 min
    # lag-time. Total absorbed fraction (BD plus ALD entering
    # systemic circulation) was 93 percent and dose-independent.
    # The split between systemic BD and systemic ALD was
    # dose-dependent: 30 percent BD / 70 percent ALD at the
    # 1.58 mmol/kg dose; 55 percent BD / 45 percent ALD at the
    # 6.34 mmol/kg dose. The model defaults to the low-dose
    # 1.58 mmol/kg fraction (30 percent as BD). To simulate the
    # high-dose split, edit fr_bd in the loaded model or override
    # logitfr_bd directly. See vignette Assumptions and deviations.
    # ============================================================
    lka        <- log(log(2) / 0.98) ; label("BD oral absorption rate constant ka (1/min); ka = ln(2) / t_half_abs with t_half_abs = 0.98 min")  # Results modeling: half-life of absorption 0.98 min
    ltlag      <- log(7.5)           ; label("BD oral absorption lag time t_lag (min)")  # Results modeling: lag-time 7.5 min
    lfdepot    <- log(0.93)          ; label("BD oral total bioavailability F (sum of systemic BD plus ALD)")  # Results modeling: F_total = 93 percent, dose-independent
    logitfr_bd <- qlogis(0.30)       ; label("Logit of BD oral fraction entering systemic circulation as BD (defaults to low-dose 30 percent)")  # Results modeling: 30 percent at 1.58 mmol/kg dose

    # ============================================================
    # Between-subject variability (exponential model;
    # omega^2 = log(CV^2 + 1)). Fung 2008 Table II reports
    # percent CV for IIV.
    # ============================================================
    etalvmax       ~ 0.012027     # Table II BD CLic IIV = 11 percent -> log(1 + 0.11^2) = 0.012027 (carried on Vmax; equivalent under fixed Km)
    etalvc         ~ 0.012027     # Table II BD V IIV = 11 percent
    etalkm_ghb     ~ 0.115558     # Table II GHB Km IIV = 35 percent -> log(1 + 0.35^2) = 0.115558
    etalvc_ghb     ~ 0.025278     # Table II GHB V IIV = 16 percent -> log(1 + 0.16^2) = 0.025278
    etalvmax_etoh  ~ 0.000361     # Table II ETOH CLic IIV = 1.9 percent -> log(1 + 0.019^2) = 0.000361 (carried on Vmax; equivalent under fixed Km)
    etalkm_etoh    ~ 0.019410     # Table II ETOH Km IIV = 14 percent -> log(1 + 0.14^2) = 0.019410
    etalvc_etoh    ~ 0.039221     # Table II ETOH V IIV = 20 percent -> log(1 + 0.20^2) = 0.039221
    etalki_bd_ald  ~ 0.191942     # Table II Ki_BD/ALD IIV = 46 percent -> log(1 + 0.46^2) = 0.191942
    etalki_bd_ghb  ~ 0.389403     # Table II Ki_BD/GHB IIV = 69 percent -> log(1 + 0.69^2) = 0.389403
    etalki_bd_etoh ~ 0.128305     # Table II Ki_BD/ETOH IIV = 37 percent -> log(1 + 0.37^2) = 0.128305
    etalki_ghb_bd  ~ 0.060625     # Table II Ki_GHB/BD IIV = 25 percent -> log(1 + 0.25^2) = 0.060625
    etaltlag       ~ 0.009950     # Results modeling: 10 percent BSV on absorption lag -> log(1 + 0.10^2) = 0.009950

    # ============================================================
    # Residual unidentified variability -- combined proportional
    # plus additive model (Fung 2008 "Observation model"): Y = C *
    # (1 + eps_CV) + eps_SD. Reported separately for each
    # measured analyte (BD, GHB, ETOH); ALD was not measured.
    # ============================================================
    propSd      <- 0.175  ; label("BD proportional residual error (fraction)")     # Table II BD: CV_CP = 17.5 percent
    addSd       <- 0.0030 ; label("BD additive residual error (mmol/L)")           # Table II BD: SD_CP = 0.0030 mmol/L
    propSd_ghb  <- 0.128  ; label("GHB proportional residual error (fraction)")    # Table II GHB: CV_CP = 12.8 percent
    addSd_ghb   <- 0.0186 ; label("GHB additive residual error (mmol/L)")          # Table II GHB: SD_CP = 0.0186 mmol/L
    propSd_etoh <- 0.0232 ; label("ETOH proportional residual error (fraction)")   # Table II ETOH: CV_CP = 2.32 percent
    addSd_etoh  <- 0.233  ; label("ETOH additive residual error (mmol/L)")         # Table II ETOH: SD_CP = 0.233 mmol/L
  })

  model({
    # ============================================================
    # 1. Individual parameters
    # ============================================================
    vmax    <- exp(lvmax + etalvmax)
    km      <- exp(lkm)
    vc      <- exp(lvc + etalvc)
    k12     <- exp(lk12)
    k21     <- exp(lk21)

    vmax_ald <- exp(lvmax_ald)
    km_ald   <- exp(lkm_ald)
    vc_ald   <- exp(lvc_ald)
    k12_ald  <- exp(lk12_ald)
    k21_ald  <- exp(lk21_ald)

    cl_ghb   <- exp(lcl_ghb)
    vmax_ghb <- exp(lvmax_ghb)
    km_ghb   <- exp(lkm_ghb + etalkm_ghb)
    vc_ghb   <- exp(lvc_ghb + etalvc_ghb)
    k12_ghb  <- exp(lk12_ghb)
    k21_ghb  <- exp(lk21_ghb)

    vmax_etoh <- exp(lvmax_etoh + etalvmax_etoh)
    km_etoh   <- exp(lkm_etoh   + etalkm_etoh)
    vc_etoh   <- exp(lvc_etoh   + etalvc_etoh)
    k12_etoh  <- exp(lk12_etoh)
    k21_etoh  <- exp(lk21_etoh)

    ki_bd_ald   <- exp(lki_bd_ald   + etalki_bd_ald)
    ki_bd_ghb   <- exp(lki_bd_ghb   + etalki_bd_ghb)
    ki_bd_etoh  <- exp(lki_bd_etoh  + etalki_bd_etoh)
    ki_ghb_bd   <- exp(lki_ghb_bd   + etalki_ghb_bd)
    ki_etoh_bd  <- exp(lki_etoh_bd)

    ka    <- exp(lka)
    tlag  <- exp(ltlag + etaltlag)
    fr_bd <- 1 / (1 + exp(-logitfr_bd))

    # ============================================================
    # 2. Concentrations (molar; mmol/L)
    # ============================================================
    Cc      <- central       / vc
    Cc_ald  <- central_ald   / vc_ald
    Cc_ghb  <- central_ghb   / vc_ghb
    Cc_etoh <- central_etoh  / vc_etoh

    # ============================================================
    # 3. Michaelis-Menten elimination rates with competitive
    # inhibition. Source: Fung 2008 Eqs. 1-4. BD elimination is
    # inhibited competitively by GHB and ETOH; ALD elimination
    # by BD; GHB elimination by BD; ETOH elimination by BD. The
    # bioactivation flux BD -> ALD enters the ALD compartment;
    # the subsequent ALD -> GHB flux enters the GHB compartment.
    # ============================================================
    rate_bd_to_ald <- (vmax     * Cc)     / (km     * (1 + Cc_ghb / ki_ghb_bd + Cc_etoh / ki_etoh_bd) + Cc)
    rate_ald_to_ghb <- (vmax_ald * Cc_ald) / (km_ald * (1 + Cc / ki_bd_ald) + Cc_ald)
    rate_ghb_mm    <- (vmax_ghb * Cc_ghb) / (km_ghb * (1 + Cc / ki_bd_ghb) + Cc_ghb)
    rate_etoh_out  <- (vmax_etoh * Cc_etoh) / (km_etoh * (1 + Cc / ki_bd_etoh) + Cc_etoh)

    # ============================================================
    # 4. ODE system (Fung 2008 Eqs. 1-8). The oral depot empties
    # at first-order rate ka; the absorbed mass is split between
    # the BD central compartment (fraction fr_bd) and the ALD
    # central compartment (fraction 1 - fr_bd) to represent
    # pre-systemic conversion of BD to ALD. IV doses enter
    # central / central_ghb / central_etoh directly.
    # ============================================================
    d/dt(depot) <- -ka * depot

    d/dt(central)     <-  ka * depot * fr_bd - rate_bd_to_ald - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    d/dt(central_ald)     <-  ka * depot * (1 - fr_bd) + rate_bd_to_ald - rate_ald_to_ghb -
                              k12_ald * central_ald + k21_ald * peripheral1_ald
    d/dt(peripheral1_ald) <-  k12_ald * central_ald - k21_ald * peripheral1_ald

    d/dt(central_ghb)     <-  rate_ald_to_ghb - rate_ghb_mm - cl_ghb * Cc_ghb -
                              k12_ghb * central_ghb + k21_ghb * peripheral1_ghb
    d/dt(peripheral1_ghb) <-  k12_ghb * central_ghb - k21_ghb * peripheral1_ghb

    d/dt(central_etoh)     <- -rate_etoh_out - k12_etoh * central_etoh + k21_etoh * peripheral1_etoh
    d/dt(peripheral1_etoh) <-  k12_etoh * central_etoh - k21_etoh * peripheral1_etoh

    # ============================================================
    # 5. Bioavailability and absorption lag for the oral BD depot
    # ============================================================
    f(depot)    <- exp(lfdepot)
    alag(depot) <- tlag

    # ============================================================
    # 6. Observations and combined proportional + additive error
    # ============================================================
    Cc      ~ add(addSd)      + prop(propSd)
    Cc_ghb  ~ add(addSd_ghb)  + prop(propSd_ghb)
    Cc_etoh ~ add(addSd_etoh) + prop(propSd_etoh)
  })
}
