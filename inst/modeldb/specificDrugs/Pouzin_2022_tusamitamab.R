Pouzin_2022_tusamitamab <- function() {
  description <- "Integrated multi-analyte semi-mechanistic population PK model of tusamitamab ravtansine (SAR408701, anti-CEACAM5 IgG1-SPDB-DM4 ADC) in adults with advanced solid tumors (Pouzin 2022): explicit two-compartment disposition for DAR1-DAR8 ADC species and a separate naked-antibody (NAB) chain sharing Vc/Vp/Q, irreversible first-order DAR_n -> DAR_(n-1) deconjugation feeding a one-compartment DM4 catabolite that converts to MeDM4."
  reference <- "Pouzin C, Gibiansky L, Fagniez N, Tod M, Chadjaa M, Nguyen L. Integrated multiple analytes and semi-mechanistic population pharmacokinetic model of tusamitamab ravtansine, a DM4 anti-CEACAM5 antibody-drug conjugate. J Pharmacokinet Pharmacodyn. 2022;49(4):381-394. doi:10.1007/s10928-021-09799-0"
  vignette <- "Pouzin_2022_tusamitamab"
  units <- list(
    time = "day",
    dosing = "umol of antibody (= mg of ADC / 150)",
    concentration = "uM (umol/L) in antibody-equivalent for ADC and NAB; uM for free DM4 and MeDM4"
  )

  covariateData <- list()

  population <- list(
    n_subjects     = 254,
    n_studies      = 1,
    age_range      = "adults (specific range not tabulated in the modelling paper)",
    weight_range   = "not tabulated in the modelling paper",
    sex_female_pct = NA_real_,
    race_ethnicity = NA,
    disease_state  = "Adults with advanced solid tumors expressing CEACAM5 (colorectal, gastric, NSCLC high/low CEACAM5, SCLC).",
    dose_range     = "5-190 mg/m^2 IV in escalation cohorts (Q2W or Q3W) and 100 mg/m^2 IV Q2W in expansion cohorts; 1- to 3-h infusions.",
    regions        = "TED13751 first-in-human Phase 1 study (NCT02187848); regions not specified in the modelling paper.",
    study          = "TED13751 (NCT02187848) — open-label, non-randomized, multi-cohort first-in-human Phase 1 of SAR408701.",
    analytes       = "Conjugated ADC (DAR>=1) by Gyrolab xP immunoassay (LLOQ 0.5 ug/mL); naked antibody NAB by competitive immunoenzymatic assay after anti-DM4 immunodepletion (LLOQ 1-9.6 ug/mL); DM4 and MeDM4 by LC-MS/MS in acidified plasma (LLOQ 0.2 ng/mL); per-DAR proportions by LC-HRMS in 13 patients.",
    n_observations = "3746 ADC, 3740 DM4, 3734 MeDM4, 3734 NAB plasma concentrations across <=58 cycles (median 4 cycles per patient).",
    notes          = "Pouzin 2022 reports the structural integrated model; baseline demographics (age/weight/sex/race) for the 254-patient analysis are not tabulated in the model paper itself. The companion covariate paper (Pouzin et al., CPT Pharmacometrics Syst Pharmacol 2022, PMID 35191618) tabulates demographics and is the basis for any subsequent covariate-aware refit; this entry packages the structural base model only.",
    dose_distribution = list(
      note     = "DAR distribution in the administered drug solution. Pouzin 2022 Table 4 'F_DARi (%)' typical fixed-effect values; CV from Table 3 (administered batches). FNAB was estimated (the others were fixed at the median of 4 measured batches).",
      F_DAR8   = 0.009,  # 0.9% (fixed; IIV back-calculated to keep sum=1)
      F_DAR7   = 0.028,  # 2.8% (fixed; CV 27.2%)
      F_DAR6   = 0.071,  # 7.1% (fixed; CV 13.2%)
      F_DAR5   = 0.142,  # 14.2% (fixed; CV 10.2%)
      F_DAR4   = 0.199,  # 19.9% (fixed; CV 5.9%)
      F_DAR3   = 0.218,  # 21.8% (fixed; CV 6.6%)
      F_DAR2   = 0.175,  # 17.5% (fixed; CV 9.1%)
      F_DAR1   = 0.085,  # 8.5%  (fixed; CV 9.1%)
      F_NAB    = 0.071   # 7.1%  (estimated; IIV CV 41.8%)
    ),
    molecular_weights = list(
      note    = "Pouzin 2022 Table 2. SAR408701 average MW corresponds to an ADC carrying three [linker-DM4] moieties; the model converts DM4/MeDM4/NAB observations to ADC molar equivalent by normalising with the SAR408701 MW.",
      ADC     = 150000,  # g/mol  (average MW for DAR ~ 3 product)
      NAB     = 144522,  # g/mol  (naked antibody)
      DM4     = 780,     # g/mol
      MeDM4   = 794      # g/mol
    )
  )

  ini({
    # ----- Distribution (shared across all DAR0..DAR8 ADC species) -----
    lvc       <- log(3.37);  label("Central volume of distribution Vc, shared across DAR0-DAR8 (L)")            # Pouzin 2022 Table 4: Vc = 3.37 L
    lvp       <- log(2.54);  label("Peripheral volume of distribution Vp, shared across DAR0-DAR8 (L)")          # Pouzin 2022 Table 4: Vp = 2.54 L
    lq        <- log(0.543); label("Inter-compartmental clearance Q, shared across DAR0-DAR8 (L/day)")           # Pouzin 2022 Table 4: Q = 0.543 L/day

    # ----- Proteolytic / catabolic clearances -----
    lcladc    <- log(0.392); label("Proteolytic clearance of conjugated ADC species (DAR1-DAR8) CL_ADC (L/day)") # Pouzin 2022 Table 4: CL_ADC = 0.392 L/day
    lclnab    <- log(0.408); label("Proteolytic clearance of naked antibody (DAR0/NAB) CL_NAB (L/day)")          # Pouzin 2022 Table 4: CL_NAB = 0.408 L/day

    # ----- DAR_n -> DAR_(n-1) deconjugation rates (1/day) -----
    # All deconjugation occurs in the central compartment as a first-order
    # process. Pouzin 2022 estimates kdec1..kdec6 individually and sets
    # kdec7 = kdec8 = kdec6 (DAR7/DAR8 plasma levels too low to identify).
    lkdec1    <- log(0.0565); label("Deconjugation rate DAR1 -> DAR0 / NAB, kdec1 (1/day)")                       # Pouzin 2022 Table 4: kdec1 = 0.0565/day
    lkdec2    <- log(0.181);  label("Deconjugation rate DAR2 -> DAR1, kdec2 (1/day)")                             # Pouzin 2022 Table 4: kdec2 = 0.181/day
    lkdec3    <- log(0.340);  label("Deconjugation rate DAR3 -> DAR2, kdec3 (1/day)")                             # Pouzin 2022 Table 4: kdec3 = 0.340/day
    lkdec4    <- log(0.525);  label("Deconjugation rate DAR4 -> DAR3, kdec4 (1/day)")                             # Pouzin 2022 Table 4: kdec4 = 0.525/day
    lkdec5    <- log(0.751);  label("Deconjugation rate DAR5 -> DAR4, kdec5 (1/day)")                             # Pouzin 2022 Table 4: kdec5 = 0.751/day
    lkdec6    <- log(0.938);  label("Deconjugation rate DAR6 -> DAR5, kdec6 (1/day); also reused for kdec7, kdec8") # Pouzin 2022 Table 4: kdec6 = 0.938/day; kdec7 = kdec8 = kdec6 per Methods (DAR7/DAR8 not separately identifiable)

    # ----- DM4 and MeDM4 disposition -----
    # Pouzin 2022: V_DM4 and V_MeDM4 are each fixed at 1 L because formation-
    # limited kinetics prevented simultaneous identification of V and FR_MeDM4.
    # CL_DM4 and CL_MeDM4 are therefore APPARENT clearances (CL/V = kel) with
    # V fixed to 1 L.
    lcldm4    <- log(240);    label("Apparent DM4 clearance CL_DM4 (L/day; V_DM4 fixed to 1 L)")                 # Pouzin 2022 Table 4: CL_DM4 = 240 L/day
    lclmedm4  <- log(0.256);  label("Apparent MeDM4 clearance CL_MeDM4 (L/day; V_MeDM4 fixed to 1 L)")           # Pouzin 2022 Table 4: CL_MeDM4 = 0.256 L/day
    lvdm4     <- fixed(log(1));   label("Apparent DM4 volume V_DM4 (L; FIXED)")                                  # Pouzin 2022 Methods: V_DM4 fixed to 1 L
    lvmedm4   <- fixed(log(1));   label("Apparent MeDM4 volume V_MeDM4 (L; FIXED)")                              # Pouzin 2022 Methods: V_MeDM4 fixed to 1 L
    lfrmedm4  <- log(0.0107); label("Apparent fraction of DM4 elimination producing MeDM4, FR_MeDM4 (unitless)") # Pouzin 2022 Table 4: FR_MeDM4 = 0.0107

    # ----- Administered DAR fractions (Pouzin 2022 Tables 3 & 4) -----
    # The Monolix supplement parameterises the dose composition as ratios FR_i
    # over an implicit FR_DAR8 = 1 reference, with SUM = 1 + sum(FR_i) and
    # F_DARi = FR_i / SUM, F_DAR8 = 1/SUM. Typical values below are derived
    # from the published F_DARi (%) by dividing by F_DAR8 = 0.9%; this
    # exactly reproduces Table 4 typical fractions through the SUM transform.
    lfrdar1   <- fixed(log(0.085 / 0.009)); label("log(FR_DAR1) ratio (typical FR_DAR1 = F_DAR1/F_DAR8)")        # Pouzin 2022 Table 4: F_DAR1 = 8.5% (fixed)
    lfrdar2   <- fixed(log(0.175 / 0.009)); label("log(FR_DAR2) ratio (typical FR_DAR2 = F_DAR2/F_DAR8)")        # Pouzin 2022 Table 4: F_DAR2 = 17.5% (fixed)
    lfrdar3   <- fixed(log(0.218 / 0.009)); label("log(FR_DAR3) ratio (typical FR_DAR3 = F_DAR3/F_DAR8)")        # Pouzin 2022 Table 4: F_DAR3 = 21.8% (fixed)
    lfrdar4   <- fixed(log(0.199 / 0.009)); label("log(FR_DAR4) ratio (typical FR_DAR4 = F_DAR4/F_DAR8)")        # Pouzin 2022 Table 4: F_DAR4 = 19.9% (fixed)
    lfrdar5   <- fixed(log(0.142 / 0.009)); label("log(FR_DAR5) ratio (typical FR_DAR5 = F_DAR5/F_DAR8)")        # Pouzin 2022 Table 4: F_DAR5 = 14.2% (fixed)
    lfrdar6   <- fixed(log(0.071 / 0.009)); label("log(FR_DAR6) ratio (typical FR_DAR6 = F_DAR6/F_DAR8)")        # Pouzin 2022 Table 4: F_DAR6 = 7.1% (fixed)
    lfrdar7   <- fixed(log(0.028 / 0.009)); label("log(FR_DAR7) ratio (typical FR_DAR7 = F_DAR7/F_DAR8)")        # Pouzin 2022 Table 4: F_DAR7 = 2.8% (fixed)
    lfrnab    <- log(0.071 / 0.009);        label("log(FR_NAB) ratio (typical FR_NAB = F_NAB/F_DAR8); estimated") # Pouzin 2022 Table 4: F_NAB = 7.1% (estimated)

    # ----- Inter-individual variability (omega = SD on log scale; Pouzin 2022 Table 4) -----
    # Variance entered as omega^2; omega is the SD reported in the "standard
    # deviation of the random effect" column of Table 4. For fixed CVs the
    # corresponding etas are wrapped in fixed(); FR_DAR8 carries no eta — its
    # variability is implicit through the FR/SUM renormalisation.
    etalvc      ~ 0.060025                  # 0.245^2  (Vc CV 24.5%)
    etalvp      ~ 0.366025                  # 0.605^2  (Vp CV 60.5%)
    etalq       ~ 0.279841                  # 0.529^2  (Q  CV 52.9%)
    etalcladc   ~ 0.219961                  # 0.469^2  (CL_ADC  CV 46.9%)
    etalclnab   ~ 0.119025                  # 0.345^2  (CL_NAB  CV 34.5%)
    etalkdec1   ~ 0.040804                  # 0.202^2  shared across all kdec_i (per Pouzin 2022 Methods)
    etalcldm4   ~ 0.133225                  # 0.365^2  (CL_DM4   CV 36.5%)
    etalclmedm4 ~ 0.427716                  # 0.654^2  (CL_MeDM4 CV 65.4%)
    etalfrmedm4 ~ 0.522729                  # 0.723^2  (FR_MeDM4 CV 72.3%)
    etalfrdar1  ~ fixed(0.008281)           # 0.091^2  (FR_DAR1; CV from Table 3 administered batch)
    etalfrdar2  ~ fixed(0.008281)           # 0.091^2  (FR_DAR2)
    etalfrdar3  ~ fixed(0.004356)           # 0.066^2  (FR_DAR3)
    etalfrdar4  ~ fixed(0.003481)           # 0.059^2  (FR_DAR4)
    etalfrdar5  ~ fixed(0.010404)           # 0.102^2  (FR_DAR5)
    etalfrdar6  ~ fixed(0.017424)           # 0.132^2  (FR_DAR6)
    etalfrdar7  ~ fixed(0.073984)           # 0.272^2  (FR_DAR7)
    etalfrnab   ~ 0.174724                  # 0.418^2  (FR_NAB; estimated, not fixed)

    # ----- Residual error (Pouzin 2022 Table 4) -----
    # ADC: combined error model. The published additive SD of 1.03 ug/mL is
    # converted to uM (umol/L) by dividing by the SAR408701 MW = 150,000 g/mol
    # and multiplying by 1000 mL/L: 1.03 / 150,000 * 1000 = 0.006867 uM.
    CcaddSd      <- 0.006867; label("Additive residual error on ADC concentration (uM)")                          # Pouzin 2022 Table 4: a_ADC = 1.03 ug/mL (= 0.006867 uM)
    CcpropSd     <- 0.089;    label("Proportional residual error on ADC concentration (fraction)")                # Pouzin 2022 Table 4: b_ADC   = 8.9%
    CnabpropSd   <- 0.260;    label("Proportional residual error on NAB concentration (fraction)")                # Pouzin 2022 Table 4: b_NAB   = 26.0%
    Cdm4propSd   <- 0.335;    label("Proportional residual error on DM4 concentration (fraction)")                # Pouzin 2022 Table 4: b_DM4   = 33.5%
    Cmedm4propSd <- 0.500;    label("Proportional residual error on MeDM4 concentration (fraction)")              # Pouzin 2022 Table 4: b_MeDM4 = 50.0%
    DARavgaddSd  <- 0.219;    label("Additive residual error on average DAR (DAR units)")                         # Pouzin 2022 Table 4: a_DARaverage = 0.219
  })

  model({
    # ----- Individual structural parameters -----
    vc       <- exp(lvc + etalvc)
    vp       <- exp(lvp + etalvp)
    q        <- exp(lq  + etalq)
    cladc    <- exp(lcladc   + etalcladc)
    clnab    <- exp(lclnab   + etalclnab)
    cldm4    <- exp(lcldm4   + etalcldm4)
    clmedm4  <- exp(lclmedm4 + etalclmedm4)
    vdm4     <- exp(lvdm4)
    vmedm4   <- exp(lvmedm4)
    frmedm4  <- exp(lfrmedm4 + etalfrmedm4)

    # Shared deconjugation eta (etalkdec1) applied to every kdec_i — same
    # IIV across the chain per Pouzin 2022 Methods. kdec7 and kdec8 are
    # constrained equal to kdec6.
    kdec1 <- exp(lkdec1 + etalkdec1)
    kdec2 <- exp(lkdec2 + etalkdec1)
    kdec3 <- exp(lkdec3 + etalkdec1)
    kdec4 <- exp(lkdec4 + etalkdec1)
    kdec5 <- exp(lkdec5 + etalkdec1)
    kdec6 <- exp(lkdec6 + etalkdec1)
    kdec7 <- kdec6
    kdec8 <- kdec6

    # ----- Administered DAR fractions: FR/SUM transform (Pouzin 2022 supplement) -----
    frdar1 <- exp(lfrdar1 + etalfrdar1)
    frdar2 <- exp(lfrdar2 + etalfrdar2)
    frdar3 <- exp(lfrdar3 + etalfrdar3)
    frdar4 <- exp(lfrdar4 + etalfrdar4)
    frdar5 <- exp(lfrdar5 + etalfrdar5)
    frdar6 <- exp(lfrdar6 + etalfrdar6)
    frdar7 <- exp(lfrdar7 + etalfrdar7)
    frnab  <- exp(lfrnab  + etalfrnab)
    sumfr  <- 1 + frdar1 + frdar2 + frdar3 + frdar4 + frdar5 + frdar6 + frdar7 + frnab
    fdar8  <- 1      / sumfr
    fdar7  <- frdar7 / sumfr
    fdar6  <- frdar6 / sumfr
    fdar5  <- frdar5 / sumfr
    fdar4  <- frdar4 / sumfr
    fdar3  <- frdar3 / sumfr
    fdar2  <- frdar2 / sumfr
    fdar1  <- frdar1 / sumfr
    fnab   <- frnab  / sumfr

    # ----- Micro-constants (shared distribution kinetics) -----
    keladc   <- cladc   / vc
    kelnab   <- clnab   / vc
    keldm4   <- cldm4   / vdm4
    kelmedm4 <- clmedm4 / vmedm4
    k12      <- q / vc
    k21      <- q / vp

    # ----- ODE system: 9 parallel two-compartment chains + DM4 + MeDM4 -----
    # All ADC and NAB amounts are tracked in umol of antibody equivalent.
    # DAR_n -> DAR_(n-1) is an irreversible first-order central process.

    # DAR8 chain: receives no upstream input (top of chain).
    d/dt(dar8_central)     <- -keladc * dar8_central - kdec8 * dar8_central -
                               k12 * dar8_central + k21 * dar8_peripheral1
    d/dt(dar8_peripheral1) <-  k12 * dar8_central - k21 * dar8_peripheral1

    # DAR7..DAR1 chains: each receives kdec_(n+1) * upstream_central and loses
    # via proteolysis (keladc) and deconjugation (kdec_n) plus 2-compartment
    # exchange.
    d/dt(dar7_central)     <-  kdec8 * dar8_central - keladc * dar7_central -
                               kdec7 * dar7_central -
                               k12 * dar7_central + k21 * dar7_peripheral1
    d/dt(dar7_peripheral1) <-  k12 * dar7_central - k21 * dar7_peripheral1

    d/dt(dar6_central)     <-  kdec7 * dar7_central - keladc * dar6_central -
                               kdec6 * dar6_central -
                               k12 * dar6_central + k21 * dar6_peripheral1
    d/dt(dar6_peripheral1) <-  k12 * dar6_central - k21 * dar6_peripheral1

    d/dt(dar5_central)     <-  kdec6 * dar6_central - keladc * dar5_central -
                               kdec5 * dar5_central -
                               k12 * dar5_central + k21 * dar5_peripheral1
    d/dt(dar5_peripheral1) <-  k12 * dar5_central - k21 * dar5_peripheral1

    d/dt(dar4_central)     <-  kdec5 * dar5_central - keladc * dar4_central -
                               kdec4 * dar4_central -
                               k12 * dar4_central + k21 * dar4_peripheral1
    d/dt(dar4_peripheral1) <-  k12 * dar4_central - k21 * dar4_peripheral1

    d/dt(dar3_central)     <-  kdec4 * dar4_central - keladc * dar3_central -
                               kdec3 * dar3_central -
                               k12 * dar3_central + k21 * dar3_peripheral1
    d/dt(dar3_peripheral1) <-  k12 * dar3_central - k21 * dar3_peripheral1

    d/dt(dar2_central)     <-  kdec3 * dar3_central - keladc * dar2_central -
                               kdec2 * dar2_central -
                               k12 * dar2_central + k21 * dar2_peripheral1
    d/dt(dar2_peripheral1) <-  k12 * dar2_central - k21 * dar2_peripheral1

    d/dt(dar1_central)     <-  kdec2 * dar2_central - keladc * dar1_central -
                               kdec1 * dar1_central -
                               k12 * dar1_central + k21 * dar1_peripheral1
    d/dt(dar1_peripheral1) <-  k12 * dar1_central - k21 * dar1_peripheral1

    # NAB (DAR0) chain: receives from DAR1 deconjugation; NAB clearance only,
    # no further deconjugation.
    d/dt(nab_central)     <-  kdec1 * dar1_central - kelnab * nab_central -
                              k12 * nab_central + k21 * nab_peripheral1
    d/dt(nab_peripheral1) <-  k12 * nab_central - k21 * nab_peripheral1

    # DM4 catabolite: every DAR>=1 deconjugation event releases exactly one
    # DM4 molecule (1:1 molar). Pouzin 2022 supplement Online Resource 1
    # input_1 expression.
    dm4_input <- kdec1 * dar1_central + kdec2 * dar2_central +
                 kdec3 * dar3_central + kdec4 * dar4_central +
                 kdec5 * dar5_central + kdec6 * dar6_central +
                 kdec7 * dar7_central + kdec8 * dar8_central
    d/dt(dm4)   <- dm4_input - keldm4 * dm4

    # MeDM4 catabolite: formed sequentially from DM4 elimination, with
    # apparent fraction frmedm4 routed into MeDM4.
    d/dt(medm4) <- frmedm4 * keldm4 * dm4 - kelmedm4 * medm4

    # ----- Bioavailability: split a single dose across 9 ADC chains -----
    # rxode2 / nlmixr2 do not natively split a single dose row across
    # multiple compartments; the vignette therefore issues 9 dose events at
    # each administration (one per chain). Setting f(<chain>) = f<chain>
    # encodes the per-chain administered fraction so the user enters the
    # TOTAL antibody amount on each row and rxode2 scales by f<chain>.
    f(dar1_central) <- fdar1
    f(dar2_central) <- fdar2
    f(dar3_central) <- fdar3
    f(dar4_central) <- fdar4
    f(dar5_central) <- fdar5
    f(dar6_central) <- fdar6
    f(dar7_central) <- fdar7
    f(dar8_central) <- fdar8
    f(nab_central)  <- fnab

    # ----- Observed concentrations -----
    Cc     <- (dar1_central + dar2_central + dar3_central + dar4_central +
               dar5_central + dar6_central + dar7_central + dar8_central) / vc  # SAR408701 (DAR>=1) by Gyrolab
    Cnab   <- nab_central / vc                                                  # NAB (DAR0) by depletion immunoassay
    Cdm4   <- dm4   / vdm4                                                      # Free DM4 by LC-MS/MS
    Cmedm4 <- medm4 / vmedm4                                                    # Free MeDM4 by LC-MS/MS

    # Average DAR: total conjugated DM4 / total antibody (ADC + NAB).
    Ctab    <- Cc + Cnab
    DARavg  <- (1 * dar1_central + 2 * dar2_central + 3 * dar3_central +
                4 * dar4_central + 5 * dar5_central + 6 * dar6_central +
                7 * dar7_central + 8 * dar8_central) / vc / Ctab

    Cc     ~ add(CcaddSd) + prop(CcpropSd)
    Cnab   ~ prop(CnabpropSd)
    Cdm4   ~ prop(Cdm4propSd)
    Cmedm4 ~ prop(Cmedm4propSd)
    DARavg ~ add(DARavgaddSd)
  })
}
