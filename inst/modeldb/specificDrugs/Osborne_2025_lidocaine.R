Osborne_2025_lidocaine <- function() {
  description <- "Two-compartment population PK model for unbound lidocaine given by intraperitoneal or intravenous bolus plus infusion after laparoscopic colon resection in adults (Osborne 2025). Intraperitoneal absorption is first order, parameterised by an absorption half-life, with a fraction of the intraperitoneal dose lost before reaching plasma. Total lidocaine is reconstructed from unbound lidocaine and a postoperative alpha-1-acid glycoprotein (AAG) turnover model whose synthesis is stimulated by a decaying surgical stimulus (fraction unbound = 1 / (1 + AAG / KD)). An exploratory direct inhibitory sigmoid Emax model of unbound lidocaine on the visual analogue pain score at rest is combined with a sigmoid natural-history pain-resolution function of time after the end of surgery. Clearances scale allometrically with total body weight (exponent 3/4), volumes linearly and the absorption half-life with exponent 1/4, standardised to 70 kg; AAG production carries a postmenstrual-age maturation function estimated from paediatric AAG data."
  reference <- "Osborne KW, MacFater WS, Anderson BJ, Svirskis D, Hill AG, Hannam JA. Pharmacokinetics of Intraperitoneal Lidocaine for Sustained Postoperative Analgesia in Adults. Eur J Drug Metab Pharmacokinet. 2025;50(4):295-306. doi:10.1007/s13318-025-00948-1"
  vignette <- "Osborne_2025_lidocaine"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Intraperitoneal doses go to depot_ip, intravenous doses to central.
  dosing <- c("depot_ip", "central")

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Theory-based allometry standardised to 70 kg (Supplementary Methods 1, Equation 11): clearances (CLU, Q, CL_AAG) and the AAG production rate scale with (WT/70)^0.75; volumes (V1, V2, VIP, V_AAG) with WT/70; the absorption half-life with (WT/70)^0.25 (supplementary model code 1, FSZT). Fat-free mass and normal fat mass were tested and not retained (Table S2).",
      source_name = "WT"
    ),
    PAGE = list(
      description = "Postmenstrual age",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = "Weeks (not the register's default months), because the maturation constants PMA50 = 32 weeks and Hill = 8.94 are on that scale (supplementary model code 1, PMAW). Drives the maturation fraction FMAT on the AAG production rate (paper Eq. 5). The AAG model was fitted to paediatric data (Booker 1996, Linnarsson 2023; Table S1); for any adult FMAT is 1 to machine precision, so PAGE = 52.18 * age_years + 40 is adequate.",
      source_name = "PMAW"
    ),
    T_SURG = list(
      description = "Duration of the surgical operation",
      units = "h",
      type = "continuous",
      reference_category = NULL,
      notes = "Time from the lidocaine start (time 0, just before surgical dissection) to the end of surgery. Natural-history pain resolution runs on time after the end of surgery, t - T_SURG, floored at 0 (supplementary model code 2: TNOW = TIME - DURATION). Used only by the pain-score output.",
      source_name = "DURATION"
    )
  )

  covariatesDataExcluded <- list(
    FFM = list(
      description = "Fat-free mass",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested as the size descriptor for CL and V (Supplementary Methods 1, Equation 12 with Ffat = 0) and not retained (Table S2: -1.8 points on CL, +16.2 on V).",
      source_name = "FFM"
    )
  )

  compartmentData <- list(
    depot_ip = list(analyte = "lidocaine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "lidocaine unbound", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "lidocaine unbound", units = "mg", specimen = "plasma", verified = TRUE),
    aag = list(analyte = "alpha-1-acid glycoprotein", units = "umol/L", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 56L,
    n_studies = 5L,
    n_observations = 1268L,
    age_range = "45-89 years (lidocaine RCT)",
    age_median = "74 years (intraperitoneal arm), 68 years (intravenous arm)",
    weight_range = "not reported in the source paper",
    sex_female_pct = NA_real_,
    race_ethnicity = "not reported in the source paper",
    disease_state = "Adults undergoing elective laparoscopic colon resection (lidocaine PK and pain scores); children after major surgery (AAG data); adults after laparoscopic colectomy or cholecystectomy given placebo (natural-history pain scores)",
    dose_range = "Lidocaine 2 mg/kg ideal body weight bolus (maximum 100 mg) before surgical dissection, then 1.5 mg/kg/h for up to 72 h, intraperitoneal (n = 28) or intravenous (n = 28)",
    regions = "New Zealand (lidocaine RCT)",
    notes = "Lidocaine data from the MacFater double-blind RCT (NCT03105193; Osborne 2025 Table 1): 388 unbound and 388 total plasma concentrations and 275 resting VAS pain scores from 56 adults. Pooled with 210 AAG concentrations from 141 children after major surgery (Booker 1996, n = 116; Linnarsson 2023, n = 25; Table S1) and 17 mean placebo-arm pain scores from adults after laparoscopic colectomy (n = 32) or cholecystectomy (n = 10). n_observations counts 766 lidocaine concentrations + 292 pain scores + 210 AAG concentrations. Weight and sex were not reported in this paper."
  )

  ini({
    # Unbound lidocaine PK, standardised to 70 kg. Values from supplementary
    # model code 1 $THETA (final estimates; they match Table 2).
    ltabs <- log(0.543); label("Intraperitoneal absorption half-life (h)") # Code 1 POP_TABS_IP = 0.543; Table 2 t1/2abs = 0.54 h
    lfdepot <- log(0.664); label("Fraction of the intraperitoneal dose reaching plasma (fraction)") # Code 1 POP_FFPM = 0.664; Table 2 F = 0.66
    lcl <- log(121); label("Unbound clearance (L/h/70 kg)") # Table 2 CLU = 121 L/h/70 kg
    lvc <- log(103); label("Central volume of unbound lidocaine (L/70 kg)") # Table 2 V1 = 103 L/70 kg
    lq <- log(195); label("Intercompartmental clearance (L/h/70 kg)") # Table 2 Q = 195 L/h/70 kg
    lvp <- log(555); label("Peripheral volume of unbound lidocaine (L/70 kg)") # Table 2 V2 = 555 L/70 kg
    lvip <- fixed(log(0.1)); label("Volume of the intraperitoneal compartment (L/70 kg)") # Table 2 VIP = 0.1 L/70 kg, held constant (Methods)
    e_wt_cl <- fixed(0.75); label("Allometric exponent on clearances and AAG production (unitless)") # Supplementary Methods 1 Eq. 11; code 1 FSZCL
    e_wt_vc <- fixed(1); label("Allometric exponent on volumes (unitless)") # Supplementary Methods 1; code 1 FSZV
    e_wt_tabs <- fixed(0.25); label("Allometric exponent on the absorption half-life (unitless)") # Code 1 FSZT = FSZV^0.25

    # AAG turnover and lidocaine binding (Table 2 'AAG model' / 'Total binding model')
    lkin_aag <- log(0.0579); label("AAG production rate (umol/h/70 kg)") # Code 1 POP_RATEIN = 0.0579; Table 2 RATEIN_AAG = 0.058
    lcl_aag <- log(0.00627); label("AAG clearance (L/h/70 kg)") # Table 2 CL_AAG = 0.00627 (code 1 initial 0.0063)
    lvc_aag <- fixed(log(3.5)); label("Plasma volume for AAG (L/70 kg)") # Table 2 V_AAG = 3.5 L/70 kg, held constant (code 1: plasma volume of a 70 kg male)
    lthalf_stim <- log(34.64); label("Half-life of the surgical stimulus on AAG synthesis (h)") # Table 2 t1/2 = 34.64 h
    lslope_stim <- log(11.9); label("Slope of AAG synthesis on the surgical stimulus (unitless)") # Code 1 POPSLOPE = 11.9; Table 2 Slope = 12
    lkaag <- log(2.98); label("Lidocaine-AAG dissociation constant KD (umol/L)") # Table 2 KD = 2.98 umol/L
    ltm50_kin_aag <- fixed(log(32)); label("Postmenstrual age at 50% maturation of AAG production (weeks)") # Code 1 POPPMA50 = 32 weeks, held constant
    lhill_kin_aag <- fixed(log(8.94)); label("Hill coefficient of AAG production maturation (unitless)") # Code 1 POPHILL = 8.94, held constant

    # Pain-score PD (Table 3; supplementary model code 2)
    le0 <- fixed(log(10)); label("Baseline pain score without analgesia (VAS cm)") # Table 3 E0 = 10.0 VAS cm, held constant
    lec50 <- fixed(log(0.206)); label("Unbound lidocaine concentration at half-maximal effect (mg/L)") # Code 2 POP_C50 = 0.206; Table 3 C50 = 0.21 mg/L, taken from the literature
    logitemax <- logit(0.593); label("Maximal fractional drug reduction of the pain score, logit scale (unitless)") # Code 2 POP_EMAX = 0.593; Table 3 Emax = 0.59
    lhill <- log(4); label("Hill coefficient of the lidocaine effect (unitless)") # Table 3 HILL = 4
    lemax_pain <- fixed(log(0.983)); label("Maximal fractional natural pain resolution (fraction)") # Code 2 POP_PAINEMX = 0.983; Table 3 Emax,PAIN = 0.98, held constant
    lhill_pain <- log(0.33); label("Hill coefficient of natural pain resolution (unitless)") # Table 3 HILLPAIN = 0.33
    lt50_pain <- log(16.3); label("Time after surgery to half-maximal pain resolution (h)") # Table 3 T50,PAIN = 16.3 h

    # IIV. Table 2 / 3 %PPV = sqrt(omega) * 100. The BLOCK(5) carries FIX on its
    # last element in code 1, which in NONMEM holds the whole block constant.
    # Order: TABS, CLU, VU, Q, V2.
    etaltabs + etalcl + etalvc + etalq + etalvp ~ fixed(c(
      0.385,
      -0.11, 0.058,
      -0.197, 0.00722, 0.195,
      0.117, -0.0239, -0.0662, 0.0827,
      -0.057, 0.0843, -0.0893, 0.0398, 0.233
    ))
    etalkin_aag + etalcl_aag ~ c(0.15, 0.261, 0.513) # Code 1 OMEGA BLOCK(2) PPVRATEIN, PPVCLAAG; Table 2 PPV 38.7% and 71.6%
    etaltm50_kin_aag ~ fixed(0.018) # Code 1 OMEGA PPVPMA50 = 0.018, held constant
    etalec50 ~ 0.255 # Code 2 OMEGA PPVC50 = 0.255; Table 3 PPV 50.5%
    etalogitemax ~ 0.470 # Code 2 OMEGA PPVEMX = 0.470 on the logit scale; Table 3 PPV 68.6%
    etalhill_pain ~ 0.074 # Code 2 OMEGA PPVPAINHL = 0.074; Table 3 PPV 27.2%
    etalt50_pain ~ 9.89 # Code 2 OMEGA PPVPAINT50 = 9.89; Table 3 PPV 314%

    # Residual error. The code's thetas multiply EPS(1) with SIGMA 1 FIX, so each
    # is the SD itself; Table 2 / 3 print sqrt() of the proportional thetas.
    propSd <- 0.271; label("Proportional residual SD, total lidocaine (fraction)") # Code 1 CV_LTOT = 0.271 (Table 2 Errprop,LTOT 52% = sqrt(0.271))
    addSd <- 0.0419; label("Additive residual SD, total lidocaine (mg/L)") # Code 1 SD_LTOT = 0.0419; Table 2 Erradd,LTOT = 0.04
    propSd_Cunbound <- 0.26; label("Proportional residual SD, unbound lidocaine (fraction)") # Code 1 CV_LUN = 0.26 (Table 2 Errprop,LUB 51% = sqrt(0.26))
    addSd_Cunbound <- 0.0304; label("Additive residual SD, unbound lidocaine (mg/L)") # Code 1 SD_LUN = 0.0304; Table 2 Erradd,LUB = 0.03
    propSd_aag <- sqrt(0.0372); label("Proportional residual SD, AAG (fraction)") # Code 1 SIGMA ERRCV = 0.0372 (variance); SD = 0.1929
    propSd_score <- 0.462; label("Proportional residual SD, pain score (fraction)") # Code 2 CV_PDE = 0.462 (Table 3 Errprop,PD 68% = sqrt(0.462))
    addSd_score <- 0.127; label("Additive residual SD, pain score (VAS cm)") # Code 2 SD_PDE = 0.127; Table 3 Erradd,PD = 0.13
  })

  model({
    # Size and maturation (supplementary model code 1)
    fsz_cl <- (WT / 70)^e_wt_cl
    fsz_v <- (WT / 70)^e_wt_vc
    fsz_t <- (WT / 70)^e_wt_tabs
    tm50_kin_aag <- exp(ltm50_kin_aag + etaltm50_kin_aag)
    hill_kin_aag <- exp(lhill_kin_aag)
    fmat <- PAGE^hill_kin_aag / (PAGE^hill_kin_aag + tm50_kin_aag^hill_kin_aag)

    # Unbound lidocaine PK
    tabs <- exp(ltabs + etaltabs) * fsz_t
    ka <- log(2) / tabs
    fdepot <- exp(lfdepot)
    cl <- exp(lcl + etalcl) * fsz_cl
    vc <- exp(lvc + etalvc) * fsz_v
    q <- exp(lq + etalq) * fsz_cl
    vp <- exp(lvp + etalvp) * fsz_v
    vip <- exp(lvip) * fsz_v

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # AAG turnover with a surgical stimulus decaying from 1 at time 0
    kin_aag <- exp(lkin_aag + etalkin_aag) * fsz_cl * fmat
    cl_aag <- exp(lcl_aag + etalcl_aag) * fsz_cl
    vc_aag <- exp(lvc_aag) * fsz_v
    thalf_stim <- exp(lthalf_stim)
    slope_stim <- exp(lslope_stim)
    kaag <- exp(lkaag)
    stim <- exp(-log(2) / thalf_stim * t)
    syn <- 1 + slope_stim * stim

    # The whole intraperitoneal dose leaves the depot at rate ka; only the
    # fraction fdepot reaches the central compartment (code 1 DADT(1), DADT(2)).
    d/dt(depot_ip) <- -ka * depot_ip
    d/dt(central) <- fdepot * ka * depot_ip - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    # AAG concentration (umol/L); code 1 carries the amount A(4) = AAG * V_AAG
    d/dt(aag) <- (kin_aag * syn - cl_aag * aag) / vc_aag
    aag(0) <- kin_aag / cl_aag

    Cunbound <- central / vc
    Cip <- depot_ip / vip
    fu <- 1 / (1 + aag / kaag)
    Cc <- Cunbound / fu

    # Pain score (supplementary model code 2)
    e0 <- exp(le0)
    ec50 <- exp(lec50 + etalec50)
    emax <- expit(logitemax + etalogitemax)
    hill <- exp(lhill)
    emax_pain <- exp(lemax_pain)
    hill_pain <- exp(lhill_pain + etalhill_pain)
    t50_pain <- exp(lt50_pain + etalt50_pain)
    tpost <- t - T_SURG
    if (tpost <= 0) tpost <- 0
    pain_res <- emax_pain * tpost^hill_pain / (tpost^hill_pain + t50_pain^hill_pain)
    drug_eff <- emax * Cunbound^hill / (Cunbound^hill + ec50^hill)
    score <- e0 * (1 - pain_res) * (1 - drug_eff)

    Cc ~ add(addSd) + prop(propSd)
    Cunbound ~ add(addSd_Cunbound) + prop(propSd_Cunbound)
    aag ~ prop(propSd_aag)
    score ~ add(addSd_score) + prop(propSd_score)
  })
}
