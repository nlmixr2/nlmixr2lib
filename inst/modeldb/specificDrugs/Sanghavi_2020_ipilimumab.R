Sanghavi_2020_ipilimumab <- function() {
  description <- "Two-compartment population PK model for intravenous ipilimumab (anti-CTLA-4 IgG1) with time-varying clearance via a sigmoid Emax function in patients with advanced solid tumors receiving ipilimumab alone or in combination with nivolumab (Sanghavi 2020)"
  reference <- "Sanghavi K, Zhang J, Zhao X, et al. Population Pharmacokinetics of Ipilimumab in Combination With Nivolumab in Patients With Advanced Solid Tumors. CPT Pharmacometrics Syst Pharmacol. 2020;9(1):29-39. doi:10.1002/psp4.12477"
  vignette <- "Sanghavi_2020_ipilimumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL and Q (exponent CL_BBWT) and on VC and VP (exponent V_BBWT) with reference weight 80 kg (Sanghavi 2020 Figure 1 reference patient).",
      source_name        = "BBWT"
    ),
    LDH = list(
      description        = "Baseline serum lactate dehydrogenase",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Effect on CL via the literal Sanghavi 2020 form (log(LDH)/log(217))^CL_logBLDH — a power of a ratio of logs, not the conventional (LDH/ref)^exponent. Reference 217 U/L (Figure 1 reference patient). The unusual functional form is what makes the BLDH effect span only ~10-15% across the 5th/95th LDH percentiles, matching the paper's narrative that the BLDH effect was <20% and not clinically meaningful; the conventional power-of-ratio form would inflate the effect to >100% at the same percentile range.",
      source_name        = "BLDH"
    ),
    TUMTP_SCLC = list(
      description        = "Tumor-type indicator for small cell lung cancer",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (all other tumor types: melanoma, NSCLC, CRC, HCC, RCC)",
      notes              = "Exponential effect on CL (Sanghavi 2020 Table 2: CL_SCLC). Derived from a categorical tumor-type column TUMTP as TUMTP_SCLC = as.integer(TUMTP == 'SCLC'). The reference group is melanoma; effects of NSCLC, CRC, HCC, and RCC vs. melanoma were tested in the full model but did not survive backward elimination.",
      source_name        = "TUMTP"
    ),
    LINE_1L = list(
      description        = "Line-of-therapy indicator: 1 = first-line, 0 = second-line or greater",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (2L+, second-line or greater)",
      notes              = "Exponential effect on CL (Sanghavi 2020 Table 2: CL_LINE = -0.0949). 1L treatment is associated with ~9% lower CL relative to 2L+.",
      source_name        = "LINE"
    ),
    NIVO_1Q3W = list(
      description        = "Nivolumab dose-regimen indicator: 1 = co-administered nivolumab 1 mg/kg every 3 weeks, 0 = otherwise",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no nivolumab or any other nivolumab regimen)",
      notes              = "Exponential effect on ipilimumab CL (Sanghavi 2020 Table 2: CL_N1Q3W = 0.0950). Derived from the source NIVO_REGIMEN column as NIVO_1Q3W = as.integer(NIVO_REGIMEN == '1 mg/kg Q3W'). Other nivolumab regimens (0.3 mg/kg Q3W, 1 mg/kg Q2W, 3 mg/kg Q3W) were tested but only the 1 mg/kg Q3W and 3 mg/kg Q2W indicators were retained in the final model.",
      source_name        = "NIVO_REGIMEN"
    ),
    NIVO_3Q2W = list(
      description        = "Nivolumab dose-regimen indicator: 1 = co-administered nivolumab 3 mg/kg every 2 weeks, 0 = otherwise",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no nivolumab or any other nivolumab regimen)",
      notes              = "Exponential effect on ipilimumab CL (Sanghavi 2020 Table 2: CL_N3Q2W = 0.191). Derived from the source NIVO_REGIMEN column as NIVO_3Q2W = as.integer(NIVO_REGIMEN == '3 mg/kg Q2W'). Paired with NIVO_1Q3W; both indicators are 0 for ipilimumab monotherapy and for nivolumab regimens not retained in the final model.",
      source_name        = "NIVO_REGIMEN"
    ),
    COMBO_NIVO = list(
      description        = "Combination-therapy indicator: 1 = ipilimumab co-administered with any nivolumab regimen, 0 = ipilimumab monotherapy",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ipilimumab monotherapy)",
      notes              = "Additive effect on the Emax parameter of the time-varying CL function (Sanghavi 2020 Table 2: Emax_COMBO = -0.202). Captures the greater decrease in ipilimumab CL over time observed when ipilimumab is given with nivolumab (any regimen) compared with monotherapy. Distinct from the per-regimen NIVO_1Q3W / NIVO_3Q2W indicators on baseline CL.",
      source_name        = "COMBO"
    )
  )

  population <- list(
    n_subjects     = 3411L,
    n_studies      = 16L,
    age_range      = "not tabulated in main-text Table 1 (pooled solid-tumor adult oncology population)",
    weight_range   = "36.8-181 kg",
    weight_median  = "76.8 kg",
    disease_state  = "Advanced / metastatic solid tumors: melanoma 50.4%, non-small cell lung cancer 17.2%, renal cell carcinoma 13.1%, small cell lung cancer 5.2%, hepatocellular carcinoma 3.8%, colorectal cancer 3.6% (Sanghavi 2020 Table 1).",
    dose_range     = "Ipilimumab 0.3-10 mg/kg IV; regimens included 1 mg/kg Q3W, 3 mg/kg Q3W, 1 mg/kg every 6 weeks, and 1 mg/kg every 12 weeks. Nivolumab co-administration regimens: 0.3 mg/kg Q3W, 1 mg/kg Q2W, 1 mg/kg Q3W, 3 mg/kg Q2W, or 3 mg/kg Q3W. Approved combination evaluated: ipilimumab 3 mg/kg Q3W + nivolumab 1 mg/kg Q3W for 4 doses.",
    regions        = "Multinational (16 trials spanning two phase I, two phase I/II, eight phase II, three phase III, and one phase IIIb/IV).",
    notes          = "Pooled data from 16 clinical trials, 12,545 ipilimumab serum concentrations. Monotherapy N = 893; combination with nivolumab N = 2,518. Baseline LDH median 217 U/L (range 74-6,245); baseline ALB median 4.1 g/dL (range 1.8-5.3); baseline tumor size median 6.29 cm. Performance status 0 in 57.3%, 1 in 41.3%, >=2 in 1.4%. Baseline demographics from Sanghavi 2020 Table 1."
  )

  ini({
    # Structural parameters at the reference patient: 80 kg BBWT, 217 U/L
    # BLDH, melanoma tumor type, ipilimumab monotherapy, 2L+ line of
    # therapy (Sanghavi 2020 Figure 1 reference-patient definition).
    # CL0 and Q are reported in the source as mL/hour and converted to
    # L/day (mL/h * 24 / 1000 = mL/h * 0.024) so the model carries time
    # in days; T50 below is also converted from hours to days.
    lcl    <- log(14.1 * 0.024); label("Baseline CL0 at reference covariates (L/day)")              # Sanghavi 2020 Table 2: CL0_REF = 14.1 mL/hour
    lvc    <- log(3.95);         label("Central volume VC at reference (L)")                        # Sanghavi 2020 Table 2: VC_REF = 3.95 L
    lq     <- log(27.9 * 0.024); label("Intercompartmental clearance Q at reference (L/day)")       # Sanghavi 2020 Table 2: Q_REF = 27.9 mL/hour
    lvp    <- log(3.18);         label("Peripheral volume VP at reference (L)")                     # Sanghavi 2020 Table 2: VP_REF = 3.18 L

    # Time-varying CL parameters: CL(t) = CL0 * exp(Emax_i * t^HILL / (T50^HILL + t^HILL)),
    # where Emax_i = Emax + e_combo_emax * COMBO_NIVO + etaEmax.
    # The fixed-effect parameter is named Emax (the monotherapy reference
    # value) so that the IIV pairs with it as 'etaEmax' per the
    # eta<param> naming convention. Emax is signed and dimensionless and
    # enters additively, so it is kept on the linear scale rather than
    # log-transformed. T50 and HILL are positive and log-transformed;
    # T50 converted from 2,540 h to days.
    Emax     <- -0.0644;          label("Reference (monotherapy) Emax of time-varying CL (unitless)") # Sanghavi 2020 Table 2: Emax_REF
    lt50     <- log(2540 / 24);   label("log T50 - time at which CL change reaches half of Emax (day)") # Sanghavi 2020 Table 2: T50 = 2,540 hour
    lhill    <- log(7.43);        label("log HILL - sigmoid shape coefficient of the time-on-CL function (unitless)") # Sanghavi 2020 Table 2: HILL = 7.43

    # Allometric exponents on baseline body weight (reference 80 kg).
    # CL_BBWT applies to both CL and Q; V_BBWT applies to both VC and VP.
    e_wt_cl <- 0.694;             label("Power exponent of WT on CL and Q (unitless)")              # Sanghavi 2020 Table 2: CL_BBWT
    e_wt_v  <- 0.600;             label("Power exponent of WT on VC and VP (unitless)")             # Sanghavi 2020 Table 2: V_BBWT

    # Power exponent on the ratio of logs: (log(LDH)/log(217))^e_logldh_cl.
    # This is a literal reading of the Sanghavi 2020 final-model equation
    # (Results "Final model"), where BLDH is log-transformed and the result
    # is then raised to a power. It is NOT the conventional (LDH/ref)^theta
    # form (verified by checking that this form yields the ~10-15% effect
    # at the 5th/95th BLDH percentiles described in Results, whereas the
    # conventional form would yield >100%).
    e_logldh_cl <- 0.703;         label("Power exponent on log(LDH)/log(217) for CL (unitless)")     # Sanghavi 2020 Table 2: CL_log-BLDH

    # Categorical covariate effects on CL applied as exp(coef * indicator).
    e_sclc_cl   <- -0.124;        label("Exponential coefficient of SCLC tumor type on CL")         # Sanghavi 2020 Table 2: CL_SCLC
    e_line_cl   <- -0.0949;       label("Exponential coefficient of 1L vs 2L+ line of therapy on CL") # Sanghavi 2020 Table 2: CL_LINE
    e_n1q3w_cl  <-  0.0950;       label("Exponential coefficient of nivolumab 1 mg/kg Q3W on CL")   # Sanghavi 2020 Table 2: CL_N1Q3W
    e_n3q2w_cl  <-  0.191;        label("Exponential coefficient of nivolumab 3 mg/kg Q2W on CL")   # Sanghavi 2020 Table 2: CL_N3Q2W

    # Additive effect of any-regimen nivolumab combination therapy on
    # Emax (linear, not exponential). Applied on the same linear scale
    # as Emax_ref and the etaEmax random effect.
    e_combo_emax <- -0.202;       label("Additive effect of nivolumab combination therapy on Emax") # Sanghavi 2020 Table 2: Emax_COMBO

    # IIV: log-normal on CL and VC (correlated 2x2 block); additive
    # normal on Emax (independent). Variance / covariance entered in
    # lower-triangular row-major order:
    #   row 1: omega^2_CL                  = 0.112    (33.4% CV)
    #   row 2: cov(CL,VC), omega^2_VC      = 0.0404, 0.0884 (29.7% CV; r 0.406)
    # Source: Sanghavi 2020 Table 2 random-effects rows.
    etalcl + etalvc ~ c(
      0.112,
      0.0404, 0.0884
    )
    # Independent normal eta on Emax. The "(0.126)" parenthetical in
    # Table 2 is sqrt(0.0158), reported as a standard deviation rather
    # than CV% because Emax itself is on a linear (not log) scale.
    etaEmax ~ 0.0158                                                                                # Sanghavi 2020 Table 2: omega^2_Emax

    # Combined proportional + additive residual error model. Treated as
    # standard deviations on the linear scale, matching the convention
    # used by Masters 2022 avelumab and Budha 2023 tislelizumab in this
    # package (both papers report sigma in the same column form).
    propSd <- 0.223;              label("Proportional residual error SD (fraction)")                # Sanghavi 2020 Table 2: Proportional
    addSd  <- 0.607;              label("Additive residual error SD (ug/mL)")                       # Sanghavi 2020 Table 2: Additive (ug/mL)
  })
  model({
    # Baseline (time-zero) CL with covariates per the Sanghavi 2020 final
    # model (Results "Final model" equation):
    #   CL0_i = CL0_REF * (BBWT/80)^CL_BBWT * (log(BLDH)/log(217))^CL_logBLDH
    #         * exp(CL_SCLC*I_SCLC) * exp(CL_LINE*I_1L)
    #         * exp(CL_N1Q3W*I_N1Q3W) * exp(CL_N3Q2W*I_N3Q2W) * exp(eta_CL).
    # Note the BLDH term raises a ratio of logs to a power (literal source
    # form), not the conventional (BLDH/ref)^theta.
    cl0 <- exp(lcl + etalcl) *
      (WT / 80)^e_wt_cl *
      (log(LDH) / log(217))^e_logldh_cl *
      exp(e_sclc_cl  * TUMTP_SCLC) *
      exp(e_line_cl  * LINE_1L) *
      exp(e_n1q3w_cl * NIVO_1Q3W) *
      exp(e_n3q2w_cl * NIVO_3Q2W)

    # Volumes and Q with allometric scaling (BBWT exponents shared:
    # CL and Q use CL_BBWT; VC and VP use V_BBWT).
    vc <- exp(lvc + etalvc) * (WT / 80)^e_wt_v
    q  <- exp(lq)           * (WT / 80)^e_wt_cl
    vp <- exp(lvp)          * (WT / 80)^e_wt_v

    # Individual Emax: monotherapy reference (Emax) + additive combo
    # offset for any-regimen nivolumab + additive normal eta.
    Emax_i <- Emax + e_combo_emax * COMBO_NIVO + etaEmax

    # Hill exponent and T50 in linear space.
    hill <- exp(lhill)
    t50  <- exp(lt50)

    # Time-varying CL: t is the within-subject simulation time (days
    # since first dose). At t = 0, CL = cl0; as t -> infinity,
    # CL -> cl0 * exp(Emax_i).
    cl <- cl0 * exp(Emax_i * t^hill / (t50^hill + t^hill))

    # Two-compartment micro-constants and ODEs.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Dose in mg, V in L -> central/V has units mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
