Eaton_2025_dabigatran_sheep <- function() {
  description <- "Preclinical (sheep). Integrated two-compartment IV pharmacokinetic plus effect-compartment sigmoid Emax pharmacodynamic model for dabigatran and its reversal agent idarucizumab in anaesthetised sheep, describing thromboelastographic reaction time (R-time). Parameters are allometrically standardised to a 70 kg body weight for cross-species comparison (Eaton 2025 Tables 2-3 and the supplementary NM-TRAN control stream)."
  reference <- "Eaton MP, Nadtochiy SM, Stefanos T, Anderson BJ. Dabigatran pharmacokinetic-pharmacodynamic in sheep: Informing dose for anticoagulation during cardiopulmonary bypass. Perfusion. 2025;40(1):183-191. doi:10.1177/02676591231226291"
  vignette <- "Eaton_2025_dabigatran_sheep"
  units <- list(time = "min", dosing = "mg", concentration = "mg/L")
  # Dabigatran doses into `central` (IV); the idarucizumab unit dose goes into
  # `depot_kpd`. Declared explicitly because buildModelDb()'s heuristic only
  # recognises `depot` / `central`.
  dosing <- c("central", "depot_kpd")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against the supplementary NM-TRAN $MODEL /
  # $DES / $ERROR blocks (COMP(DIGAB), COMP(PERIPH), COMP(EFFECT), COMP(IDA);
  # S1 = V1 and S2 = V2, so states 1-2 are amounts, while state 3 is scaled by
  # the NONMEM default S3 = 1 and is therefore already a concentration).
  compartmentData <- list(
    central     = list(analyte = "dabigatran",   units = "mg",   specimen = "plasma",              verified = TRUE),
    peripheral1 = list(analyte = "dabigatran",   units = "mg",   specimen = "tissue",              verified = TRUE),
    effect      = list(analyte = "dabigatran",   units = "mg/L", specimen = "not applicable",      verified = TRUE),
    depot_kpd   = list(analyte = "idarucizumab", units = "unit", specimen = "administration site", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric size descriptor, reference weight 70 kg. Exponent 0.75 on CL and Q, 1 on Vc and Vp (Eaton 2025 Methods, 'Pharmacokinetics'), and -0.25 on ke0 -- the NM-TRAN control stream scales the equilibration HALF-TIME T1/2keo by (WT/70)^0.25 (FSZT), which is the same as scaling ke0 by (WT/70)^-0.25. Sheep body weights were 28.7-41.8 kg (Table 1); the 70 kg standard is a reporting convention, not a studied weight.",
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "sheep",
    n_subjects     = 5L,                                 # Eaton 2025 Methods, 'Animals and materials'; Table 1
    n_studies      = 1L,                                 # Single-centre study, University of Rochester
    age_range      = "6 months",                         # Eaton 2025 Table 1: all five sheep 6 months
    age_median     = "6 months",                         # Eaton 2025 Table 1
    weight_range   = "28.7-41.8 kg",                     # Eaton 2025 Table 1
    weight_median  = "33.9 kg",                          # Eaton 2025 Table 1: 28.7, 33.9, 33.8, 33.9, 41.8 kg
    sex_female_pct = 100,                                # Eaton 2025 Table 1: all five sheep female
    race_ethnicity = NULL,                               # Not applicable (animal study)
    disease_state  = "Healthy sheep under general anaesthesia (ketamine 4 mg/kg plus midazolam 0.4 mg/kg induction, isoflurane 1-4% maintenance); no cardiopulmonary bypass was run during the PKPD study itself.",
    dose_range     = "Dabigatran 4 mg/kg IV over 1 min at time 0; idarucizumab 15 mg/kg IV over 30 s at 120 min.",
    regions        = "USA (University of Rochester Medical Center).",
    n_observations = "Plasma dabigatran and thromboelastographic R-time at baseline and 5, 15, 30, 60, 90, 120 min after dabigatran, then 5, 15, 30, 60, 120, 240, 480 min and 24 h after idarucizumab (Eaton 2025 Methods).",
    notes          = "Demographics from Eaton 2025 Table 1. Parameter estimates were obtained in NONMEM 7.5 (ADVAN6, FOCE with INTERACTION) using a sequential PPP&D approach for the pharmacodynamic parameters; the full control stream is reproduced in the article's Supplementary NM-TRAN Code."
  )

  ini({
    # -----------------------------------------------------------------
    # Dabigatran pharmacokinetics -- standardised to a 70 kg body weight
    # (Eaton 2025 Table 2, 'Standardized dabigatran population
    # pharmacokinetic parameter estimates').
    # -----------------------------------------------------------------
    lcl <- log(0.0453); label("Clearance standardised to 70 kg (L/min)")                        # Eaton 2025 Table 2: CL_STD = 0.0453 L/min/70 kg (95% CI 0.02-0.053)
    lvc <- log(2.94);   label("Central volume of distribution standardised to 70 kg (L)")       # Eaton 2025 Table 2: V1_STD = 2.94 L/70 kg (95% CI 2.27-3.14)
    lq  <- log(0.268);  label("Intercompartmental clearance standardised to 70 kg (L/min)")     # Eaton 2025 Table 2: Q_STD = 0.268 L/min/70 kg (95% CI 0.175-0.400)
    lvp <- log(9.51);   label("Peripheral volume of distribution standardised to 70 kg (L)")    # Eaton 2025 Table 2: V2_STD = 9.51 L/70 kg (95% CI 7.07-14.7)

    # Allometric exponents, reference weight 70 kg. Eaton 2025 Methods
    # ('Pharmacokinetics', Eq. 1): P_i = P_STD * (W_i/W_STD)^EXP with
    # "The EXP exponent was 0.75 for clearance and 1 for distribution
    # volumes". Supplementary NM-TRAN: FSZV = WT/70, FSZCL = FSZV**0.75,
    # applied as CL = FSZCL*CL, Q = FSZCL*Q, V1 = FSZV*V1, V2 = FSZV*V2.
    e_wt_cl_q  <- fixed(0.75); label("Shared allometric exponent on CL and Q (unitless)")   # Eaton 2025 Methods Eq. 1; supplement FSZCL = (WT/70)**0.75
    e_wt_vc_vp <- fixed(1);    label("Shared allometric exponent on Vc and Vp (unitless)")  # Eaton 2025 Methods Eq. 1; supplement FSZV = WT/70

    # -----------------------------------------------------------------
    # Dabigatran pharmacodynamics -- sigmoid Emax on R-time driven by the
    # effect-site concentration (Eaton 2025 Table 3 and Eq. 2).
    # -----------------------------------------------------------------
    le0   <- log(0.681);      label("Baseline R-time in the absence of drug (min)")                    # Eaton 2025 Table 3: E0 = 0.681 min (95% CI 0.547-0.796)
    lemax <- fixed(log(180)); label("Maximum drug-attributable increase in R-time (min)")              # Eaton 2025 Table 3: Emax = 180 min FIX; supplement $THETA (100,180.,300) FIX
    lec50 <- log(64.2);       label("Effect-site concentration producing half of Emax (mg/L)")         # Eaton 2025 Table 3: Ce50 = 64.2 mg/L (95% CI 56.6-71.92)
    lhill <- fixed(log(1));   label("Hill coefficient of the concentration-response curve (unitless)") # Eaton 2025 Table 3: N = 1 FIX

    # Plasma / effect-site equilibration. The paper reports the half-time
    # T1/2keo = 1.04 min (Table 3); ke0 = ln(2)/T1/2keo = 0.6665 1/min. The
    # supplementary control stream computes exactly this (TEQ then
    # KEQ = LN2/TEQ) and scales the HALF-TIME allometrically by
    # FSZT = (WT/70)**0.25, which is an exponent of -0.25 on ke0 itself.
    lke0     <- log(log(2) / 1.04); label("Effect-compartment equilibration rate constant ke0 (1/min)")  # Eaton 2025 Table 3: T1/2keo = 1.04 min (95% CI 0.936-2.14) -> ke0 = ln(2)/1.04
    e_wt_ke0 <- fixed(-0.25);       label("Allometric exponent on ke0 (unitless)")                       # supplement: FSZT = (WT/70)**0.25 applied to T1/2keo, i.e. -0.25 on ke0

    # -----------------------------------------------------------------
    # Idarucizumab reversal arm. Eaton 2025 Methods ('Pharmacodynamics'):
    # "The impact of idarucizumab was modelled by a unit bolus input into a
    # fourth compartment with first order elimination described by a rate
    # constant (KIDA). Input duration (DURIDA) was estimated as a parameter.
    # The [idarucizumab] concentration (CIDA) was assumed to directly
    # decrease R time and this relationship was described using a slope
    # constant (SLOPEIDA). Observed R time was the sum of Dabigatran
    # (EFFECTDABI) and idarucizumab effects (EFFECTIDA)."
    # A dose of amt = 1 into depot_kpd represents one 15 mg/kg idarucizumab
    # administration; SLOPEIDA therefore carries units of min per unit dose.
    # -----------------------------------------------------------------
    lkel_ida   <- log(0.0218); label("Idarucizumab K-PD elimination rate constant (1/min)")              # Eaton 2025 Table 3: K_IDA = 0.0218 (95% CI 0.009-0.075)
    ldur_ida   <- log(0.0751); label("Zero-order input duration of the idarucizumab unit dose (min)")    # Eaton 2025 Table 3: DUR_IDA = 0.0751
    lslope_ida <- log(8.24);   label("Idarucizumab effect slope on R-time (min per unit dose)")          # Eaton 2025 Table 3: SLOPE_IDA = 8.24 (95% CI 6.26-12.45)

    # -----------------------------------------------------------------
    # Between-subject variability. Eaton 2025 reports %BSV in Tables 2-3;
    # the supplementary $OMEGA initial estimates confirm that %BSV is the
    # approximate CV sqrt(omega^2)*100 rather than the exact log-normal CV
    # (e.g. Q %BSV 51.4 -> 0.514^2 = 0.264 vs the supplement's $OMEGA
    # initial 0.265; V1 7.4 -> 0.005476 vs 0.00591; C50 19.9 -> 0.0396 vs
    # 0.0395; E0 15.1 -> 0.0228 vs 0.0227). omega^2 = (%BSV/100)^2 is
    # therefore used throughout. The authors fit OMEGA BLOCK(4) on the PK
    # parameters and OMEGA BLOCK(2) on C50/E0, but no final off-diagonal
    # estimates are published, so only the diagonals are encoded here --
    # see the vignette's "Assumptions and deviations" section.
    # -----------------------------------------------------------------
    etalcl      ~ 0.045369  # Eaton 2025 Table 2: %BSV CL  = 21.3  -> 0.213^2
    etalvc      ~ 0.005476  # Eaton 2025 Table 2: %BSV V1  = 7.4   -> 0.074^2
    etalq       ~ 0.264196  # Eaton 2025 Table 2: %BSV Q   = 51.4  -> 0.514^2
    etalvp      ~ fixed(0.156816)  # Eaton 2025 Table 2: %BSV V2 = 39.6 -> 0.396^2; supplement $OMEGA BLOCK(4) carries PPVV2 as FIX
    etale0      ~ 0.022801  # Eaton 2025 Table 3: %BSV E0   = 15.1  -> 0.151^2
    etalec50    ~ 0.039601  # Eaton 2025 Table 3: %BSV Ce50 = 19.9  -> 0.199^2
    etalkel_ida ~ 1.110916  # Eaton 2025 Table 3: %BSV KIDA = 105.4 -> 1.054^2

    # -----------------------------------------------------------------
    # Residual unidentified variability. The supplementary $ERROR block
    # builds a combined error SD as SQRT((PROP*PROP)+(ADD*ADD)) against a
    # $SIGMA fixed to 1, which is exactly nlmixr2's default combined2
    # parameterisation for `add() + prop()`.
    #
    # Scale of the proportional terms: both tables print the proportional
    # RUV under a "(%)" header, so the tabulated number is a percentage and
    # the fraction is that number / 100. The PD row settles the convention
    # -- RUV PROP 37.0 can only be the fraction 0.370 (the supplement's
    # $THETA initial for RUV_CVPD is 0.412), so the PK row's 0.47 is the
    # fraction 0.0047 (supplement $THETA initial for RUV_CVCP is 0.0005).
    # -----------------------------------------------------------------
    addSd  <- 1.32;   label("Additive residual error on dabigatran concentration (mg/L)")     # Eaton 2025 Table 2: RUV_ADD = 1.32 mg/L (95% CI 0.79-1.49)
    propSd <- 0.0047; label("Proportional residual error on dabigatran concentration (fraction)")  # Eaton 2025 Table 2: RUV_PROP = 0.47% (95% CI 0.42-0.49%)

    addSd_Rtime  <- 0.149; label("Additive residual error on R-time (min)")             # Eaton 2025 Table 3: RUV_ADD = 0.149 min (95% CI 0.002-0.231)
    propSd_Rtime <- 0.370; label("Proportional residual error on R-time (fraction)")    # Eaton 2025 Table 3: RUV_PROP = 37.0% (95% CI 34.3-46.9%)
  })
  model({
    # 1. Individual parameters, allometrically scaled to the 70 kg standard
    #    (supplement $PK: CL = FSZCL*CL*EXP(PPVCL), V1 = FSZV*V1*EXP(PPVV1),
    #    TEQ = FSZT*TEQ, KEQ = LN2/TEQ).
    cl  <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q
    q   <- exp(lq  + etalq)  * (WT / 70)^e_wt_cl_q
    vc  <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    vp  <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc_vp
    ke0 <- exp(lke0)         * (WT / 70)^e_wt_ke0

    e0   <- exp(le0 + etale0)
    emax <- exp(lemax)
    ec50 <- exp(lec50 + etalec50)
    hill <- exp(lhill)

    kel_ida   <- exp(lkel_ida + etalkel_ida)
    dur_ida   <- exp(ldur_ida)
    slope_ida <- exp(lslope_ida)

    # 2. Two-compartment dabigatran disposition. Supplement $DES:
    #    DADT(1) = -DCP*CL - DCP*Q + DC2*Q, DADT(2) = Q*DCP - DC2*Q,
    #    with DCP = A(1)/V1 and DC2 = A(2)/V2.
    d/dt(central)     <- q * (peripheral1 / vp) - (cl + q) * (central / vc)
    d/dt(peripheral1) <- q * (central / vc)     -  q       * (peripheral1 / vp)
    Cc                <- central / vc

    # 3. Effect compartment. Supplement $DES: DADT(3) = KEQ*(DCP - DCE)
    #    with DCE = A(3); NONMEM's default S3 = 1 makes state 3 a
    #    concentration, so no volume scaling is applied here either.
    #    The lower guard mirrors "IF (CE.LE.0) CE = 1D-10" in $ERROR.
    d/dt(effect) <- ke0 * (Cc - effect)
    Ce           <- max(effect, 1e-10)

    # 4. Idarucizumab K-PD arm. Supplement $DES: DADT(4) = -DIDA*KIDA with
    #    DIDA = A(4), and $PK: D4 = DURIDA. Dose amt = 1 into depot_kpd.
    d/dt(depot_kpd) <- -kel_ida * depot_kpd
    dur(depot_kpd)  <- dur_ida

    # 5. Observed R-time is the dabigatran sigmoid Emax effect (Eaton 2025
    #    Eq. 2: Effect_DABI = E0 + Emax*Ce^N/(Ce50^N + Ce^N)) minus the
    #    idarucizumab effect (supplement $ERROR: FX1 = E0 + EMAX*CEN/(C50N+CEN),
    #    FX2 = SLOPE*CIDA, FX = FX1 - FX2).
    Rtime <- e0 + emax * Ce^hill / (ec50^hill + Ce^hill) - slope_ida * depot_kpd

    # 6. Residual error (combined additive + proportional, combined2 form)
    Cc    ~ add(addSd)       + prop(propSd)
    Rtime ~ add(addSd_Rtime) + prop(propSd_Rtime)
  })
}
