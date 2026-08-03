Lee_2023_tripegfilgrastim <- function() {
  description <- paste(
    "Semi-mechanistic population PK/PD model for tripegfilgrastim (a",
    "PEGylated long-acting recombinant human G-CSF) and absolute",
    "neutrophil count (ANC) in healthy Korean adults and Korean pediatric",
    "patients with solid tumors receiving chemotherapy. Subcutaneous drug",
    "enters a depot compartment and is absorbed at first-order rate KSC",
    "into a total-drug compartment. Pharmacodynamics-mediated drug",
    "disposition (PDMDD) uses a quasi-equilibrium quadratic between total",
    "drug and the circulating G-CSF-receptor pool, giving free drug (FDC)",
    "and bound complex (RDC). Free drug is cleared linearly (CLD/VD) and",
    "bound drug is internalised (KINT). Granulopoiesis is a five-state",
    "receptor chain (stem -> mitotic -> post-mitotic I -> post-mitotic II",
    "-> circulating blood receptors) with baselines KP/KTR and KP/KC;",
    "ANC = 1000 * RB / SR. Drug binding stimulates receptor production",
    "(ST1 = 1 + STM1*driver) and bone-marrow transit (ST2 = 1 +",
    "STM2*driver), where the driver is the fraction of receptors bound by",
    "exogenous drug plus endogenous G-CSF. Endogenous G-CSF is carried as",
    "its own turnover compartment (fixed KIN, KEL, GCSF0 from Quartino",
    "2014) and a negative-feedback term FB = (RB0/RB)^GAM modulates",
    "receptor production. Chemotherapy is a KPD virtual compartment with",
    "lag LAG; its output (KCHM*CHM) scaled by CHMSL adds to the",
    "mitotic-cell elimination rate. Study population (DIS_HEALTHY) shifts",
    "VD and KD, age scales KSC (exponent -0.97, reference 18.5 y), body",
    "weight scales KINT (exponent 1.7, reference 55.1 kg), and baseline",
    "ANC scales KP (exponential, reference 2106 cells/uL)."
  )
  reference <- paste(
    "Lee S, Hong KT, Jang I-J, Yu K-S, Kang HJ, Oh J. (2023).",
    "Semimechanistic pharmacokinetic-pharmacodynamic model of",
    "tripegfilgrastim for pediatric patients after chemotherapy.",
    "CPT Pharmacometrics Syst Pharmacol 12(9):1319-1334.",
    "doi:10.1002/psp4.13012.",
    sep = " "
  )
  vignette <- "Lee_2023_tripegfilgrastim"

  paper_specific_compartments <- c("depot_kpd_chemotherapy", "endogenous_gcsf")

  units <- list(
    time          = "hour",
    dosing        = "ug (tripegfilgrastim); mg (chemotherapy KPD input)",
    concentration = "ug/L (serum G-CSF, free exogenous drug plus free endogenous G-CSF)",
    ANC           = "cells/uL",
    notes         = paste(
      "Tripegfilgrastim dose enters compartment 'depot' in ug (the paper's",
      "weight-based doses are 30, 60, 100 and 300 ug/kg; the proposed",
      "fixed doses are 1.5, 2.5, 4 and 6 mg = 1500, 2500, 4000 and 6000",
      "ug). Chemotherapy dose in mg enters compartment",
      "'depot_kpd_chemotherapy'; its effect on the mitotic pool is delayed",
      "by LAG = 171 h and then decays first-order at KCHM. Receptor-pool",
      "states (precursor1-4, circ) are concentrations in ug/L; ANC is",
      "recovered as 1000 * circ / SR in cells/uL. NOTE: Lee 2023 does not",
      "report the numeric chemotherapy mass entered into the KPD",
      "compartment (16 different regimens were pooled), so that dose is a",
      "user-supplied simulation input -- see the vignette Errata."
    )
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power-form covariate on KINT (the G-CSF-receptor-mediated",
        "internalisation rate) with exponent 1.7 and reference 55.1 kg, the",
        "population-weighted mean weight of the pooled analysis set",
        "(Lee 2023 Equation 6 and Table 2). Doubling body weight raises KINT",
        "about 2^1.7 = 3.25-fold, which the paper describes as a three-fold",
        "increase (Discussion). Observed range 18-75 kg (Table 1)."
      ),
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power-form covariate on KSC (the first-order subcutaneous",
        "absorption rate) with exponent -0.97 and reference 18.5 years, the",
        "population-weighted mean age of the pooled analysis set",
        "(Lee 2023 Equation 7 and Table 2). Absorption is faster in younger",
        "subjects: relative to age 18, KSC is about nine-fold higher at age",
        "2 and three-fold higher at age 6, and about 40 percent lower at age",
        "30 (Results). The paper attributes this to increasing subcutaneous",
        "tissue thickness with age (Discussion). Observed range 6-38 years",
        "(Table 1); the paper extrapolates to ages 2-6 in its simulations."
      ),
      source_name        = "AGE"
    ),
    NEUT = list(
      description        = "Baseline absolute neutrophil count before tripegfilgrastim treatment",
      units              = "cells/uL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Exponential-form covariate on KP (the G-CSF-receptor production",
        "rate): KP_i = KP_pop * exp(0.56 * NEUT / 2106), where 2106 cells/uL",
        "is the population-weighted mean baseline ANC of the pooled analysis",
        "set (Lee 2023 Equation 8 and Table 2). The paper's source column is",
        "BSLD. Note this is a ratio-in-the-exponent form, NOT the more usual",
        "centred-deviation form, so KP is already multiplied by exp(0.56) at",
        "the reference value. Worked values from Results / Discussion that",
        "confirm the encoding: KP = 0.095 ug/L/h at NEUT = 500 (severe",
        "neutropenia), 0.12 ug/L/h at NEUT = 1500, and roughly double that at",
        "NEUT = 4000. Observed median (min-max) 2411 (408-5602) cells/uL",
        "overall; 2607 (1430-5602) in healthy adults and 1274 (408-4611) in",
        "pediatric patients (Table 1)."
      ),
      source_name        = "BSLD"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-participant indicator (1 = healthy Korean adult volunteer, 0 = Korean pediatric patient with a solid tumor receiving chemotherapy)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (pediatric patient with a solid tumor on chemotherapy)",
      notes              = paste(
        "The paper's 'study population' covariate, retained on VD and KD",
        "(Lee 2023 Equations 4 and 5, Table 2). Reference category is the",
        "pediatric chemotherapy cohort, matching the orientation used by the",
        "upstream Melhem_2018_g_csf.R extraction. Table 2 reports both",
        "typical values directly rather than the exponential coefficients:",
        "VD/F = 4.7 L (healthy) vs 12.7 L (pediatric patients) and",
        "KD = 42.2 ug/L (healthy) vs 16.2 ug/L (pediatric patients), i.e. a",
        "62 percent lower KD in patients (Abstract, Results). IMPORTANT",
        "CONFOUND: in this analysis every healthy participant is an adult",
        "(age 20-38) and every patient is pediatric (age 6-17) and receiving",
        "chemotherapy, so the indicator conflates health status, age stratum",
        "and chemotherapy exposure. The paper acknowledges this as its",
        "principal limitation ('the PK-PD model was constructed using limited",
        "numbers of healthy subjects and pediatric patients without adult",
        "patients')."
      ),
      source_name        = "POP"
    )
  )

  covariatesDataExcluded <- list(
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened in the stepwise covariate analysis (Lee 2023 Methods 'Covariate analysis') but not retained in the final model."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened in the stepwise covariate analysis (Lee 2023 Methods 'Covariate analysis') but not retained in the final model."
    ),
    BSA = list(
      description = "Body surface area",
      units       = "m^2",
      type        = "continuous",
      notes       = "Screened in the stepwise covariate analysis (Lee 2023 Methods 'Covariate analysis') but not retained in the final model."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened in the stepwise covariate analysis (Lee 2023 Methods 'Covariate analysis') but not retained in the final model. The cohort was 88.1 percent male (all 40 healthy adults were male), so the sex effect was poorly informed."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 67L,
    n_studies      = 2L,
    age_range      = "6-38 years overall; healthy adults 24 (20-38) years, pediatric patients 12 (6-17) years (median (min-max), Table 1)",
    weight_range   = "18-75 kg overall; healthy adults 69.7 (60-75) kg, pediatric patients 43 (18-67) kg (median (min-max), Table 1)",
    sex_female_pct = 11.9,
    race_ethnicity = "Korean",
    disease_state  = paste(
      "Two pooled populations: (i) 40 healthy Korean adult male volunteers",
      "receiving single subcutaneous tripegfilgrastim (8 of whom received",
      "placebo and contribute endogenous ANC profiles only); (ii) 27 Korean",
      "pediatric patients with solid tumors receiving a single subcutaneous",
      "dose 24 h after the end of chemotherapy, across 16 different",
      "chemotherapy regimen combinations."
    ),
    dose_range     = paste(
      "Healthy adults: single SC tripegfilgrastim 1.8, 3.6, 6 and 18 mg",
      "(equivalent to 30, 60, 100 and 300 ug/kg), plus placebo.",
      "Pediatric patients: single SC 60 or 100 ug/kg given 24 h after the",
      "end of chemotherapy."
    ),
    regions        = "Republic of Korea (Seoul National University Hospital)",
    notes          = paste(
      "Trial registrations NCT00959777 (healthy adults, Ahn 2013) and",
      "NCT02963389 (pediatric patients, Lee 2022). 876 tripegfilgrastim",
      "concentration samples (842 above LLOQ) and 811 ANC samples (104 from",
      "placebo-treated healthy adults). Assay LLOQ 0.7 ug/L in healthy",
      "adults and 0.017 ug/L in pediatric patients (different ELISA kits).",
      "BLQ handled by the M3 method except predose BLQ set to zero.",
      "Estimation software: Monolix 2023R1 (SAEM); bootstrap via Rsmlx 4.0",
      "(200 runs, 100 percent convergence). The structural framework, and",
      "the fixed chemotherapy parameters KCHM and CHMSL, are inherited from",
      "Melhem 2018 (see Melhem_2018_g_csf.R); the fixed endogenous G-CSF",
      "turnover parameters are from Quartino 2014."
    )
  )

  ini({
    # -----------------------------------------------------------------
    # PK structural parameters. The DIS_HEALTHY reference category is the
    # pediatric chemotherapy cohort, so lvd / lkd carry the PATIENT
    # typical values and e_hv_* shift them to the healthy-adult values.
    # -----------------------------------------------------------------
    lfsc      <- fixed(log(1))                ; label("Log FSC (subcutaneous bioavailability)")                                          # Lee 2023 Table 2 (Fsc = 1, no RSE reported) and Results 'Population PK-PD model' ('The bioavailability after subcutaneous administration of tripegfilgrastim was fixed to 1')

    lksc      <- log(0.027)                   ; label("Log KSC (first-order SC absorption rate) at AGE 18.5 y (1/h)")                                  # Lee 2023 Table 2 (Ksc = 0.027, RSE 6.1%)
    e_age_ksc <- -0.97                        ; label("Exponent of the power relationship between AGE / 18.5 y and KSC")                               # Lee 2023 Table 2 (Age effect - beta_Ksc = -0.97, RSE 13.4%) and Equation 7

    lvd       <- log(12.7)                    ; label("Log VD/F (apparent volume of distribution) -- pediatric-patient reference (L)")                 # Lee 2023 Table 2 (Vd/F pediatric patients = 12.7 L, RSE 9.2%)
    e_hv_vd   <- log(4.7 / 12.7)              ; label("Log shift on VD/F for healthy adults (VD healthy / VD pediatric)")                              # Lee 2023 Table 2 (Vd/F healthy adults = 4.7 L, RSE 9.2%) and Equation 4

    lcld      <- log(0.18)                    ; label("Log CL/F (ANC-independent linear apparent clearance) (L/h)")                                    # Lee 2023 Table 2 (CL/F = 0.18, RSE 10.8%)

    lkd       <- log(16.2)                    ; label("Log KD (tripegfilgrastim-G-CSF receptor dissociation constant) -- pediatric-patient reference (ug/L)") # Lee 2023 Table 2 (Kd pediatric patients = 16.2 ug/L, RSE 26.7%)
    e_hv_kd   <- log(42.2 / 16.2)             ; label("Log shift on KD for healthy adults (KD healthy / KD pediatric)")                                # Lee 2023 Table 2 (Kd healthy adults = 42.2 ug/L, RSE 26.7%) and Equation 5

    lkint     <- log(0.41)                    ; label("Log KINT (receptor-complex internalisation rate) at WT 55.1 kg (1/h)")                          # Lee 2023 Table 2 (Kint = 0.41, RSE 15%)
    e_wt_kint <- 1.7                          ; label("Exponent of the power relationship between WT / 55.1 kg and KINT")                              # Lee 2023 Table 2 (Weight effect - beta_Kint = 1.7, RSE 20.5%) and Equation 6

    # -----------------------------------------------------------------
    # Granulopoiesis parameters.
    # -----------------------------------------------------------------
    lkp       <- log(0.083)                   ; label("Log KP (G-CSF receptor production rate) at NEUT 0 (ug/L/h)")                                    # Lee 2023 Table 2 (Kp = 0.083, RSE 14.4%)
    e_neut_kp <- 0.56                         ; label("Coefficient of the exponential relationship between NEUT / 2106 cells/uL and KP")               # Lee 2023 Table 2 (Baseline ANC effect - beta_Kp = 0.56, RSE 13.1%) and Equation 8

    lktr      <- fixed(log(0.033))            ; label("Log KTR (bone-marrow receptor transit rate) (1/h) -- the literature value")            # Lee 2023 Table 2 (Ktr = 0.033, no RSE reported) and Methods 'Population PK-PD model development' ('The transit rate between the receptor compartments in the bone marrow was fixed to a literature value of 0.033 h-1, as 5 days have been reported for the maturation and migration of neutrophils'); note 4 transits / 120 h = 0.0333
    lkc       <- fixed(log(0.1155))           ; label("Log KC (elimination rate of neutrophils from blood into tissues) (1/h)")               # Lee 2023 Appendix S1 Monolix code comment ('KC ... (fixed to 0.1155 h-1)'); Methods reports 0.116 h-1 and Table 2 rounds to 0.12; 0.1155 = ln(2)/6 h for the stated 6-hour blood half-life
    lsr       <- log(0.54)                    ; label("Log SR (scaling factor between receptor concentration and ANC) (g per 10^9 cells)")             # Lee 2023 Table 2 (Scale = 0.54, RSE 8.5%)

    lstm1     <- log(14.4)                    ; label("Log STM1 (stimulation of the G-CSF receptor production rate)")                                  # Lee 2023 Table 2 (STM1 = 14.4, RSE 13.5%)
    lstm2     <- log(13.8)                    ; label("Log STM2 (stimulation of the bone-marrow receptor transit rate)")                               # Lee 2023 Table 2 (STM2 = 13.8, RSE 7.4%)

    # -----------------------------------------------------------------
    # Endogenous G-CSF turnover and ANC negative feedback. All fixed:
    # Lee 2023 Results, 'Parameters were not adequately estimated when the
    # self-regulating endogenous component was applied ... Therefore, we
    # fixed those parameters to the literature value.'
    # -----------------------------------------------------------------
    lgam      <- fixed(log(0.145))            ; label("Log GAM (exponent of the baseline-ANC-to-ANC negative-feedback term)")                 # Lee 2023 Results ('The negative feedback function is governed by gamma parameter, which was fixed to 0.145', ref 32); Table 2 rounds to 0.15
    lgcsf0    <- fixed(log(0.0243))           ; label("Log GCSF0 (baseline endogenous G-CSF concentration) (ug/L)")                           # Lee 2023 Results ('baseline endogenous G-CSF concentration (0.0243 ug/L)', ref 27 Quartino 2014); Table 2 rounds to 0.024
    lkel_endo <- fixed(log(0.592))            ; label("Log KEL (nonspecific linear elimination rate of endogenous G-CSF) (1/h)")              # Lee 2023 Results ('nonspecific linear elimination rate constant (0.592 h-1)', ref 27 Quartino 2014); Table 2 rounds to 0.59
    lkin_endo <- fixed(log(0.498))            ; label("Log KIN (zero-order production rate of endogenous G-CSF) (ug/L/h)")                    # Lee 2023 Results ('The endogenous G-CSF is explained with zero-order production rate (0.498 ug/L/h)', ref 27 Quartino 2014); Table 2 rounds to 0.5

    # -----------------------------------------------------------------
    # Chemotherapy KPD sub-model. Lee 2023 Results: 'Because of our
    # limited data sets, we fixed the parameters involving chemotherapy
    # effect with the literature values except for the lag time' -- so
    # KCHM and CHMSL (and their IIV) are fixed to Melhem 2018 values and
    # only LAG is estimated.
    # -----------------------------------------------------------------
    ltlag_chem <- log(171)                    ; label("Log LAG (lag between chemotherapy administration and its effect on mitotic cells) (h)")         # Lee 2023 Table 2 (Lag = 171 h, RSE 12.4%)
    lkel_chem  <- fixed(log(0.072))           ; label("Log KCHM (chemotherapy KPD elimination rate) (1/h) -- the literature value")           # Lee 2023 Table 2 (Kchemo = 0.072, no RSE reported); fixed to Melhem 2018 (ref 23), whose Table 2 reports 0.0724
    lchmsl     <- fixed(log(668))             ; label("Log CHMSL (slope relating chemotherapy KPD output to mitotic-cell loss rate) (1/mg)")           # Lee 2023 Table 2 (CHMslope = 668, no RSE reported); fixed to Melhem 2018 (ref 23), whose Table 2 also reports 668

    # -----------------------------------------------------------------
    # Inter-individual variability. The Table 2 'Random effects' rows are
    # log-scale SDs (omega), not variances -- confirmed against the
    # Results text, which reports the LAG IIV as 62.6% CV, and
    # sqrt(exp(0.57^2) - 1) = 0.626. ini() takes variances, so each entry
    # is omega^2.
    # -----------------------------------------------------------------
    etalksc       ~ 0.39^2                                                                                                                            # Lee 2023 Table 2 (Omega Ksc = 0.39, RSE 12.8%)
    etalvd        ~ 0.39^2                                                                                                                            # Lee 2023 Table 2 (Omega Vd/F = 0.39, RSE 15.1%)
    etalcld       ~ 0.5^2                                                                                                                             # Lee 2023 Table 2 (Omega CL/F = 0.5, RSE 15%)
    etalkp        ~ 0.17^2                                                                                                                            # Lee 2023 Table 2 (Omega Kp = 0.17, RSE 41.8%)
    etalkd        ~ 1.3^2                                                                                                                             # Lee 2023 Table 2 (Omega Kd = 1.3, RSE 12.5%)
    etalstm1      ~ 0.67^2                                                                                                                            # Lee 2023 Table 2 (Omega STM1 = 0.67, RSE 18.3%)
    etalstm2      ~ 0.11^2                                                                                                                            # Lee 2023 Table 2 (Omega STM2 = 0.11, RSE 34.6%)
    etalkint      ~ 0.56^2                                                                                                                            # Lee 2023 Table 2 (Omega Kint = 0.56, RSE 26%)
    etaltlag_chem ~ 0.57^2                                                                                                                            # Lee 2023 Table 2 (Omega Lag = 0.57, RSE 22.1%); Results reports the corresponding CV as 62.6%
    etalkel_chem  ~ fixed(0.26^2)                                                                                                                     # Lee 2023 Table 2 (Omega Kchemo = 0.26, no RSE); Results: 'Random effects were all estimated except for Kchemo and CHMslope, which were with the reported value' (Melhem 2018 reports 0.259)
    etalchmsl     ~ fixed(2.3^2)                                                                                                                      # Lee 2023 Table 2 (Omega CHMslope = 2.3, no RSE); per the same Results sentence (Melhem 2018 reports 2.28)

    # -----------------------------------------------------------------
    # Residual error. PK: combined additive + proportional (Results, 'a
    # combined error model was used in the PK model'). PD: additive in the
    # log domain (Methods, 'The additive residual error model with
    # log-transform on both sides in ANC was used for the PD model'),
    # which maps to lnorm() in nlmixr2.
    # -----------------------------------------------------------------
    propSd    <- 0.37                         ; label("Proportional residual SD on tripegfilgrastim concentration Cc")                                 # Lee 2023 Table 2 (b1 = 0.37, RSE 3.8%)
    addSd     <- 0.005                        ; label("Additive residual SD on tripegfilgrastim concentration Cc (ug/L)")                              # Lee 2023 Table 2 (a1 = 0.005 ug/L, RSE 27%)
    expSd_ANC <- 0.48                         ; label("Log-normal residual SD on absolute neutrophil count ANC")                                       # Lee 2023 Table 2 (a2 = 0.48, RSE 3.1%; 'Additive error for predicted ANC in log domain')
  })

  model({
    # -----------------------------------------------------------------
    # Individual parameters (Lee 2023 Equations 4-8 and Appendix S1).
    # KSC carries the AGE power effect, KINT the WT power effect, KP the
    # NEUT exponential effect, and VD / KD the DIS_HEALTHY shift.
    # -----------------------------------------------------------------
    fsc      <- exp(lfsc)
    ksc      <- exp(lksc + etalksc) * (AGE / 18.5) ^ e_age_ksc
    vd       <- exp(lvd + e_hv_vd * DIS_HEALTHY + etalvd)
    cld      <- exp(lcld + etalcld)
    kd       <- exp(lkd + e_hv_kd * DIS_HEALTHY + etalkd)
    kint     <- exp(lkint + etalkint) * (WT / 55.1) ^ e_wt_kint

    # Equation 8 normalises baseline ANC as a ratio inside the exponent
    # (NEUT / 2106), not as a centred deviation.
    kp       <- exp(lkp + etalkp) * exp(e_neut_kp * (NEUT / 2106))
    ktr      <- exp(lktr)
    kc       <- exp(lkc)
    sr       <- exp(lsr)
    stm1     <- exp(lstm1 + etalstm1)
    stm2     <- exp(lstm2 + etalstm2)

    gam      <- exp(lgam)
    gcsf0    <- exp(lgcsf0)
    kel_endo <- exp(lkel_endo)
    kin_endo <- exp(lkin_endo)

    tlag_chem <- exp(ltlag_chem + etaltlag_chem)
    kel_chem  <- exp(lkel_chem + etalkel_chem)
    chmsl     <- exp(lchmsl + etalchmsl)

    # -----------------------------------------------------------------
    # Drug elimination micro-constant. CL/F acts on FREE drug only; bound
    # drug is removed by receptor internalisation (KINT).
    # -----------------------------------------------------------------
    ke_drug  <- cld / vd

    # -----------------------------------------------------------------
    # Quasi-equilibrium binding (Lee 2023 Appendix S1 EQUATION block).
    # trc = total G-CSF receptor concentration = the circulating blood
    # receptor pool. Free drug is the positive root of the binding
    # quadratic; the same quadratic is applied independently to
    # endogenous G-CSF against the same receptor pool, exactly as
    # published (see the vignette Errata -- this is an approximation to,
    # not a solution of, true competitive binding).
    # -----------------------------------------------------------------
    trc      <- circ
    tdc      <- central / vd
    qroot_d  <- tdc - trc - kd
    fdc      <- 0.5 * (qroot_d + sqrt(qroot_d * qroot_d + 4 * kd * tdc))
    rdc      <- tdc - fdc

    tec      <- endogenous_gcsf / vd
    qroot_e  <- tec - trc - kd
    fec      <- 0.5 * (qroot_e + sqrt(qroot_e * qroot_e + 4 * kd * tec))
    rec      <- tec - fec

    # -----------------------------------------------------------------
    # Linear stimulation of receptor production (ST1) and bone-marrow
    # transit (ST2) by the bound fraction of the receptor pool. Both
    # exogenous drug (RDC) and endogenous G-CSF (REC) contribute.
    # -----------------------------------------------------------------
    st1      <- 1 + stm1 * ((rdc + rec) / trc)
    st2      <- 1 + stm2 * ((rdc + rec) / trc)

    # -----------------------------------------------------------------
    # Chemotherapy KPD: the zero-order output rate DCHM = KCHM * CHM is
    # scaled by CHMSL into a first-order mitotic-cell elimination rate.
    # -----------------------------------------------------------------
    dchm     <- kel_chem * depot_kpd_chemotherapy
    stm_chem <- chmsl * dchm

    # -----------------------------------------------------------------
    # ANC negative feedback FB = (RB0 / RB)^GAM, where RB0 = KP / KC is
    # the individual's own circulating-receptor baseline.
    # -----------------------------------------------------------------
    circ0    <- kp / kc
    fb       <- (circ0 / circ) ^ gam

    # -----------------------------------------------------------------
    # ODE system (Lee 2023 Appendix S1 EQUATION block).
    #   depot                  : SC tripegfilgrastim depot (ug)
    #   central                : total drug DT in the body (ug)
    #   endogenous_gcsf        : endogenous G-CSF EN (ug)
    #   precursor1             : SM, stem-cell receptors (ug/L)
    #   precursor2             : MT, mitotic-cell receptors (ug/L; chemo acts here)
    #   precursor3             : PM1, first post-mitotic receptors (ug/L)
    #   precursor4             : PM2, second post-mitotic receptors (ug/L)
    #   circ                   : RB, circulating blood receptors (ug/L)
    #   depot_kpd_chemotherapy : CHM, KPD chemotherapy amount (mg)
    # -----------------------------------------------------------------
    d/dt(depot)                  <- -ksc * depot
    d/dt(central)                <-  ksc * depot - ke_drug * (fdc * vd) - kint * (rdc * vd)
    d/dt(endogenous_gcsf)        <-  kin_endo - kel_endo * (fec * vd) - kint * (rec * vd)
    d/dt(precursor1)             <-  kp * st1 * fb            - ktr * st2 * precursor1
    d/dt(precursor2)             <-  ktr * st2 * precursor1   - (ktr * st2 + stm_chem) * precursor2
    d/dt(precursor3)             <-  ktr * st2 * precursor2   - ktr * st2 * precursor3
    d/dt(precursor4)             <-  ktr * st2 * precursor3   - ktr * st2 * precursor4
    d/dt(circ)                   <-  ktr * st2 * precursor4   - kc * circ
    d/dt(depot_kpd_chemotherapy) <- -kel_chem * depot_kpd_chemotherapy

    # -----------------------------------------------------------------
    # Initial conditions (Lee 2023 Appendix S1). The granulopoiesis chain
    # starts at its unstimulated steady state; the endogenous G-CSF
    # compartment starts at the fixed literature baseline GCSF0.
    # -----------------------------------------------------------------
    precursor1(0)       <- kp / ktr
    precursor2(0)       <- kp / ktr
    precursor3(0)       <- kp / ktr
    precursor4(0)       <- kp / ktr
    circ(0)             <- kp / kc
    endogenous_gcsf(0)  <- gcsf0

    # -----------------------------------------------------------------
    # Bioavailability and the chemotherapy effect lag.
    # -----------------------------------------------------------------
    f(depot)                     <- fsc
    alag(depot_kpd_chemotherapy) <- tlag_chem

    # -----------------------------------------------------------------
    # Observations (Lee 2023 Appendix S1 OUTPUT block):
    #   CONC = FDC + FEC   (free exogenous drug plus free endogenous G-CSF)
    #   ANC  = 1000 * RB / SR
    # -----------------------------------------------------------------
    Cc  <- fdc + fec
    Cc  ~ prop(propSd) + add(addSd)
    ANC <- 1000 * (circ / sr)
    ANC ~ lnorm(expSd_ANC)
  })
}
