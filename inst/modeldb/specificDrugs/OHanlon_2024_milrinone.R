OHanlon_2024_milrinone <- function() {
  description <- "One-compartment population PK model of intravenous milrinone from premature neonates to adolescents (O'Hanlon 2024), pooling 50 children after cardiac surgery with cardiopulmonary bypass and 29 very premature neonates (Paradisis 2007). Clearance scales with fat-free mass (exponent 3/4), a sigmoid postmenstrual-age maturation function and relative renal function (eGFR / normal GFR), and rises over time on milrinone from 63.7% of its asymptote at infusion start (fractional deficit 0.363, half-time 11.9 h). Volume scales linearly with fat-free mass and is 2.06 times its mature value at birth, decaying toward it with a postnatal-age half-time of 0.969 days. Correlated log-normal IIV on CL and V; combined residual error estimated separately for two LC-MS bioanalytical batches and the HPLC-UV neonatal study, with a log-normal random effect on the residual SD."
  reference <- "O'Hanlon CJ, Sumpter A, Anderson BJ, Hannam JA. Time-Varying Clearance in Milrinone Pharmacokinetics from Premature Neonates to Adolescents. Clin Pharmacokinet. 2024;63(5):695-706. doi:10.1007/s40262-024-01372-5"
  vignette <- "OHanlon_2024_milrinone"
  paper_specific_etas <- c("etaRUV")

  units <- list(time = "h", dosing = "ug", concentration = "ug/L")

  compartmentData <- list(
    central = list(analyte = "milrinone", units = "ug", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    FFM = list(
      description = "Fat-free mass",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Size descriptor for both CL (exponent 3/4) and V (exponent 1), standardised to 56.1 kg (the FFM of a 70 kg, 176 cm male). Predicted from total body mass, height, sex and age with the neonate-to-adult FFM model of O'Hanlon et al. (CPT Pharmacometrics Syst Pharmacol 2023), not the adult Janmahasatian formula. The paper's normal fat mass NFM = FFM + Ffat * (WT - FFM) reduces to FFM because Ffat was fixed to 0 for both CL and V in the final model.",
      source_name = "FFMKG"
    ),
    PAGE = list(
      description = "Postmenstrual age",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = "Declared in WEEKS (not the register default of months) because the maturation half-time TM50 = 49.8 weeks is on that scale. Gestational age at birth (40 weeks if described only as term) plus postnatal weeks, per the control stream $INPUT comment 'PMAW ; weeks postmenstrual age GEST+PNAW'.",
      source_name = "PMAW"
    ),
    PNA = list(
      description = "Postnatal age",
      units = "months",
      type = "continuous",
      reference_category = NULL,
      notes = "Canonical units are months; the paper's volume birth effect is in DAYS (control stream PNAD), so model() converts with PNA * 30.4375 days per month. Only the first few days of life matter: the excess volume has a 0.969-day half-time, so it is below 3% of its birth value by PNA 5 days.",
      source_name = "PNAD"
    ),
    RENALFUNC_REL = list(
      description = "Relative renal function, eGFR / normal GFR (1 = normal for age, size and body composition)",
      units = "(dimensionless)",
      type = "continuous",
      reference_category = NULL,
      notes = "O'Hanlon 2024 Eq. 6-8: RF = eGFR / nGFR, with eGFR = creatinine production rate / Jaffe serum creatinine and nGFR from size, PMA and PNA maturation (O'Hanlon et al. 2023 models). In the source data RF uses the post-surgery Day-1 serum creatinine carried forward (control stream RF_D1LOCF). RF was set to 1 for every Paradisis 2007 neonate (all < 2 days old, when serum creatinine reflects maternal values); set RENALFUNC_REL = 1 for any subject whose renal function is assumed normal.",
      source_name = "RF_D1LOCF"
    ),
    ASSAY_BATCH2 = list(
      description = "Second LC-MS bioanalytical batch indicator (1 = sample analysed in batch 2 of the Starship study; 0 = batch 1 or another assay)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Starship LC-MS batch 1)",
      notes = "Residual error only. Constant per patient: the control stream assigns Starship patients with original study ID OID <= 14 to batch 1 and OID > 14 to batch 2. Must be 0 when STUDY_PARADISIS2007 = 1. For simulation of a new population, leave 0 (batch 1).",
      source_name = "OID"
    ),
    STUDY_PARADISIS2007 = list(
      description = "Paradisis 2007 neonatal study indicator (1 = observation from the Paradisis et al. 2007 dose-escalation study in very premature neonates, measured by HPLC-UV; 0 = Starship LC-MS study)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Starship Children's Hospital study)",
      notes = "Residual error only (control stream SID = 2). Selects the Paradisis residual error (CV 0 FIX, SD 31.4 ug/L). For simulation of a new population, leave 0.",
      source_name = "SID"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body mass",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters the paper's normal-fat-mass size descriptor NFM = FFM + Ffat * (WT - FFM) (Eq. 2-3), but Ffat was estimated near 0 (~0.001) for both CL and V and fixed to 0 in the final model (Table 2 'Ffat CL' and 'Ffat V' = 0 FIX; Table S2 model 5 vs 6), so total body mass drops out of the final model. Needed only to compute FFM."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 79L,
    n_studies = 2L,
    n_observations = 369L,
    age_range = "Starship study: PNA 4 days to 16 years (PMA 39.1-769 weeks, 2.5th-97.5th percentile); Paradisis 2007: PNA 0.125-1.25 days, PMA 23.7-27.3 weeks",
    age_median = "Starship study: PNA 0.656 years, PMA 73.8 weeks; Paradisis 2007: PNA 0.75 days, PMA 26.0 weeks",
    weight_range = "Starship study: 3.44-80.9 kg; Paradisis 2007: 0.549-1.17 kg (2.5th-97.5th percentile)",
    weight_median = "Starship study: 8.05 kg; Paradisis 2007: 0.85 kg",
    ffm_median = "Starship study: 5.80 kg (3.02-54.6); Paradisis 2007: 0.808 kg (0.522-1.11)",
    sex_female_pct = 40.5,
    disease_state = "Starship study: children after cardiac surgery supported by cardiopulmonary bypass, milrinone for prophylaxis or treatment of low cardiac output syndrome; Paradisis 2007: very premature neonates (< 28 weeks PMA) given milrinone as prophylaxis of low systemic blood flow after birth.",
    dose_range = "Starship study: IV infusion 0.25-1.25 ug/kg/min started after the end of bypass, titrated to response; Paradisis 2007: 0.25, 0.5 or 0.75 ug/kg/min for 3 h then 0.2 ug/kg/min.",
    renal_function = "Day-1 eGFR median 54.1 mL/min/1.73 m^2 (18.1-93.3) in the Starship study; RF assumed 1 in Paradisis 2007 neonates.",
    regions = "New Zealand (Starship Children's Hospital, Auckland) and Australia (Paradisis 2007).",
    notes = "O'Hanlon 2024 Table 1 and Results. 245 Starship samples (LC-MS, LLOQ 5 ug/L, two analytical batches; 37 BLQ handled with the M3 method) and 124 Paradisis samples (HPLC-UV, LLOQ 20 ug/L, none BLQ). Female percentage pooled from Table 1 (20 of 50 plus 12 of 29). NONMEM 7.5.1, FOCE-I with Laplacian for the M3 likelihood."
  )

  ini({
    # Structural parameters, standardised to NFMstd = 56.1 kg FFM (a 70 kg
    # male) at full maturation, normal renal function, after the time-varying
    # CL deficit has resolved, and with no birth effect on V.
    lcl <- log(17.6); label("Asymptotic clearance for a 56.1 kg FFM (70 kg) mature adult with normal renal function (L/h)") # Table 2 'CL' = 17.6 L/h/70 kg; ctl $THETA POP_CL = 17.6
    lvc <- log(20.8); label("Mature volume of distribution for a 56.1 kg FFM (70 kg) adult (L)") # Table 2 'V' = 20.8 L/70 kg; ctl $THETA POP_V = 20.8

    # Theory-based allometric exponents on normal fat mass (Eq. 1, 13, 14).
    e_ffm_cl <- fixed(0.75); label("Allometric exponent of FFM on CL (unitless)") # Eq. 13 exponent 3/4; ctl FSZ_CL_RF = (NFM/NFM_STD)**0.75
    e_ffm_vc <- fixed(1); label("Allometric exponent of FFM on V (unitless)") # Eq. 14 exponent 1; ctl FSZ_V = NFM/NFM_STD

    # Sigmoid postmenstrual-age maturation of CL (Eq. 4).
    ltmat50 <- log(49.8); label("Postmenstrual age at 50% CL maturation, TM50 (weeks)") # Table 2 'TM50' = 49.8 weeks PMA; ctl TM50_RF = 49.8
    lhill_mat <- log(4.04); label("Hill coefficient of the CL maturation function (unitless)") # Table 2 'HILL' = 4.04; ctl HILL_RF = 4.04

    # Time-varying CL (Eq. 9, 13): cl(t) = cl * (1 - cl_exp_fdef * exp(-cl_exp_kdes * t)).
    cl_exp_fdef <- 0.363; label("Fractional deficit of CL below its asymptote at milrinone start, FCL (fraction)") # Table 2 'FCL' = 0.363; ctl POP_CONCL = 0.363
    lcl_exp_kdes <- log(log(2) / 11.9); label("Rate constant of the rise in CL toward its asymptote (1/h)") # Table 2 'TCL' = 11.9 h (half-time); ctl POP_TCL = 11.9; rate = ln(2) / 11.9

    # Postnatal-age birth effect on V (Eq. 5, 14): V = vc * (1 + vc_pna_famp * exp(-vc_pna_kdes * PNA_days)).
    vc_pna_famp <- 1.06; label("Fractional excess of V over its mature value at birth, FV (fraction)") # Table 2 'FV' = 1.06; ctl POP_CONV_PNA = 1.06
    lvc_pna_kdes <- log(log(2) / 0.969); label("Rate constant of the postnatal decline of excess V (1/day of PNA)") # Table 2 'TVOL' = 0.969 days PNA (half-time); ctl POP_TVOL_PNA = 0.969; rate = ln(2) / 0.969

    # IIV (ctl $OMEGA BLOCK(2)): Table 2 PPV is sqrt(omega) -> 0.359, 0.339;
    # r = 0.114 / sqrt(0.129 * 0.115) = 0.936.
    etalcl + etalvc ~ c(0.129, 0.114, 0.115) # ctl $OMEGA BLOCK(2) 0.129, 0.114, 0.115; Table 2 'PPV CL' 0.359, 'PPV V' 0.339, 'r CL V' 0.936

    # Log-normal random effect on the residual SD (ctl PPV_RUVCP).
    etaRUV ~ 0.0669 # ctl $OMEGA PPV_RUVCP = 0.0669; Table 2 'PPV RUV' = 0.259 = sqrt(0.0669)

    # Combined residual error, estimated separately per sample group (Eq. 12).
    propSd <- 0.382; label("Proportional residual SD, Starship LC-MS batch 1 (fraction)") # Table 2 'CV batch 1' = 0.382; ctl RUV_CVCP1 = 0.382
    addSd <- 34.8; label("Additive residual SD, Starship LC-MS batch 1 (ug/L)") # Table 2 'SD batch 1' = 34.8 ug/L; ctl RUV_SDCP1 = 34.8
    propSd_batch2 <- 0.332; label("Proportional residual SD, Starship LC-MS batch 2 (fraction)") # Table 2 'CV batch 2' = 0.332 (ctl RUV_CVCP2 reads 0.322; Table 2 value used)
    addSd_batch2 <- 7.86; label("Additive residual SD, Starship LC-MS batch 2 (ug/L)") # Table 2 'SD batch 2' = 7.86 ug/L; ctl RUV_SDCP2 = 7.86
    propSd_paradisis <- fixed(0); label("Proportional residual SD, Paradisis 2007 HPLC-UV (fraction)") # Table 2 'CV Paradisis' = 0 FIX; ctl RUV_CVCP3 = 0 FIX
    addSd_paradisis <- 31.4; label("Additive residual SD, Paradisis 2007 HPLC-UV (ug/L)") # Table 2 'SD Paradisis' = 31.4 ug/L; ctl RUV_SDCP3 = 31.4
  })
  model({
    # Size: Ffat = 0, so NFM = FFM and NFMstd = 56.1 kg.
    fsize_cl <- (FFM / 56.1)^e_ffm_cl
    fsize_vc <- (FFM / 56.1)^e_ffm_vc

    # Sigmoid PMA maturation of CL; PAGE in weeks.
    tmat50 <- exp(ltmat50)
    hill_mat <- exp(lhill_mat)
    fmat <- 1 / (1 + (PAGE / tmat50)^(-hill_mat))

    # Birth effect on V over postnatal age in days.
    pna_days <- PNA * 30.4375
    vc_pna_kdes <- exp(lvc_pna_kdes)
    fage_vc <- 1 + vc_pna_famp * exp(-vc_pna_kdes * pna_days)

    cl <- exp(lcl + etalcl) * fsize_cl * fmat * RENALFUNC_REL
    vc <- exp(lvc + etalvc) * fsize_vc * fage_vc

    # Time-varying CL inside the ODE; t is the time since the first milrinone
    # dose (ctl $DES T, 'TIME ; h time since first dose').
    cl_exp_kdes <- exp(lcl_exp_kdes)
    cl_time <- 1 - cl_exp_fdef * exp(-cl_exp_kdes * t)
    kel <- cl * cl_time / vc

    d/dt(central) <- -kel * central

    Cc <- central / vc

    # NONMEM SD = sqrt((CP * CV)^2 + SD^2) * exp(etaRUV), i.e. combined2
    # with both components scaled by exp(etaRUV). The group SDs are selected
    # per observation: batch 1 when both indicators are 0.
    grp_batch1 <- (1 - ASSAY_BATCH2) * (1 - STUDY_PARADISIS2007)
    err_scale <- exp(etaRUV)
    propSd_i <- (propSd * grp_batch1 + propSd_batch2 * ASSAY_BATCH2 + propSd_paradisis * STUDY_PARADISIS2007) * err_scale
    addSd_i <- (addSd * grp_batch1 + addSd_batch2 * ASSAY_BATCH2 + addSd_paradisis * STUDY_PARADISIS2007) * err_scale
    Cc ~ add(addSd_i) + prop(propSd_i)
  })
}
