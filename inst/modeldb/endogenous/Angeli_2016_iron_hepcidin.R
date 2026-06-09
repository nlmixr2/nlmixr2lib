Angeli_2016_iron_hepcidin <- function() {
  description <- "Joint turnover model of serum iron and serum hepcidin during the menstrual cycle in healthy non-menopausal women; both molecules follow first-order turnover with a menses-induced increase in elimination (kloss) shared across iron and hepcidin and a delayed post-menses rebound in synthesis (krelI, krelH) starting on day 2 of the cycle, with serum iron multiplicatively modulating hepcidin synthesis around the iron baseline."
  reference   <- "Angeli A, Laine F, Lavenu A, Ropert M, Lacut K, Gissot V, Sacher-Huvelin S, Jezequel C, Moignet A, Laviolle B, Comets E. Joint Model of Iron and Hepcidin During the Menstrual Cycle in Healthy Women. AAPS J. 2016 May;18(3):490-504. doi:10.1208/s12248-016-9875-4"
  vignette    <- "Angeli_2016_iron_hepcidin"
  units       <- list(
    time          = "day",
    dosing        = "(none -- endogenous joint turnover, no exogenous drug)",
    concentration = "umol/L (serum iron) and nmol/L (serum hepcidin)"
  )

  covariateData <- list(
    CONMED_BIRTHCONTROL = list(
      description        = "Oral hormonal contraceptive use; `1` = currently taking an oral contraceptive, `0` = not on hormonal contraception.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no oral contraceptive). Effect direction: women on contraception eliminate iron about 18% more slowly than women not on contraception (paper Discussion p. 499; Table V).",
      notes              = "Time-fixed per subject in the Angeli 2016 HEPMEN study (61% of 90 women took daily oral contraceptives, Methods p. 491). Encoded with the opposite sign of the paper's reported `beta_NOCONTRA` because the canonical CONMED_BIRTHCONTROL indicator is the inverse of the paper's `NOCONTRA` covariate (1 means on contraception in the canonical encoding; 1 means NOT on contraception in the paper).",
      source_name        = "NOCONTRA (inverted: CONMED_BIRTHCONTROL = 1 - NOCONTRA; canonical effect coefficient = -1 x paper's beta_NOCONTRA)"
    ),
    BMI = list(
      description        = "Body mass index at the inclusion visit.",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = "n/a -- enters as the power-law multiplier `(BMI / tBMI)^e_bmi_krel_iron` on the iron post-menses rebound parameter `krel_iron`, with `tBMI = 22` (a set value close to the population median 22.6 kg/m^2; Angeli 2016 Table I).",
      notes              = "Time-fixed per subject (single inclusion-visit measurement). Range in the Angeli 2016 cohort: 10th-90th percentile 19.3-26.6 kg/m^2 (Table V); mean 22.6 (Table I).",
      source_name        = "BMI"
    ),
    HGB = list(
      description        = "Baseline (end-of-cycle visit) hemoglobin concentration.",
      units              = "g/dL",
      type               = "continuous",
      reference_category = "n/a -- enters as the power-law multiplier `(HGB / tHGB)^e_hgb_krel_iron` on the iron post-menses rebound parameter `krel_iron`, with `tHGB = 13.5` (a set value close to the population median 13.48 g/dL; Angeli 2016 Table I).",
      notes              = "Per-subject 'reference value' taken as the concentration at the last-cycle visit (Angeli 2016 Methods p. 492-493); analysed centrally at the Rennes hospital laboratory for assay consistency. Range: 10th-90th percentile 12.5-14.6 g/dL (Table V); mean 13.48 (Table I).",
      source_name        = "Haemoglobin (baseline; last-visit value)"
    ),
    HIGHAM = list(
      description        = "Higham's pictorial blood-loss assessment chart (PBAC) score; semi-quantitative measure of menstrual blood loss.",
      units              = "(score; unitless)",
      type               = "continuous",
      reference_category = "n/a -- enters as the power-law multiplier `(HIGHAM / tHIGHAM)^e_higham_<param>` on hepcidin synthesis (`ksyn_hep`) and elimination (`kout_hep`) rate constants, with `tHIGHAM = 96.6` (a set value at the population mean / median; Angeli 2016 Table I).",
      notes              = "Per-subject mean of the three most recent menstrual cycles' Higham scores (Angeli 2016 Methods p. 491). Higham score >= 100 is the clinical threshold for menorrhagia (Higham 1990 BJOG). Range: 10th-90th percentile 38.2-160.2 (Angeli 2016 Table V); mean 96.6 with SD 60.5 (Table I).",
      source_name        = "HiS (Higham's score)"
    ),
    FERRITIN_BL = list(
      description        = "Baseline (end-of-cycle visit) serum ferritin concentration.",
      units              = "ug/L",
      type               = "continuous",
      reference_category = "n/a -- enters as the power-law multiplier `(FERRITIN_BL / tFERRITIN)^e_ferritin_bl_<param>` on hepcidin elimination (`kout_hep`) and hepcidin post-menses rebound (`krel_hep`), with `tFERRITIN = 53` (a set value close to the population mean 53.14 ug/L; Angeli 2016 Table I).",
      notes              = "Per-subject 'reference value' taken at the last-cycle visit. Range: 10th-90th percentile 18.6-97.4 ug/L (Angeli 2016 Table V); mean 53.14 with SD 44.3 (Table I).",
      source_name        = "Ferritin (baseline; last-visit value)"
    ),
    HT = list(
      description        = "Subject height at the inclusion visit.",
      units              = "cm",
      type               = "continuous",
      reference_category = "n/a -- enters as the power-law multiplier `(HT / tHT)^e_ht_krel_hep` on the hepcidin post-menses rebound parameter `krel_hep`, with `tHT = 165` (a set value close to the population median 165.4 cm; Angeli 2016 Table I). The very large exponent (32.7) reflects the narrow height range across the cohort (158-173 cm at the 10th-90th percentile).",
      notes              = "Time-fixed per subject (single inclusion-visit measurement). Mean 165.4 cm with SD 6.1 (Table I).",
      source_name        = "Height"
    ),
    DLOSS = list(
      description        = "Per-subject length of the menstrual period in the cycle observed during the study.",
      units              = "days",
      type               = "continuous",
      reference_category = "n/a -- the upper bound of the loss-phase indicator window; the increased elimination kloss adds to both `kout_iron` and `kout_hep` while `t < DLOSS`.",
      notes              = "Fixed to the individual's observed menses length (Angeli 2016 Methods p. 494). Inclusion criterion required menses length between 3 and 5 days (Methods p. 491). For typical-value simulations a value of 4 days is reasonable.",
      source_name        = "dloss"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 90,
    n_studies      = 1,
    age_range      = "19-44 years",
    age_median     = "27.6 years (mean)",
    weight_range   = "mean 61.8 kg (SD 9.2)",
    weight_median  = "61.8 kg (mean reported)",
    sex_female_pct = 100,
    race_ethnicity = "Not reported; French metropolitan cohort recruited at four university hospitals.",
    disease_state  = "Healthy non-menopausal women with regular menstrual cycles and menses length 3-5 days; subjects with low baseline iron concentrations or low haemoglobin were excluded.",
    dose_range     = "(none -- no exogenous drug administered)",
    regions        = "France (Brest, Nantes, Rennes, Tours university hospitals)",
    notes          = "Multicenter observational HEPMEN study (Methods p. 491). 514 paired iron + hepcidin measurements across six fasting morning blood samples per subject spaced over a single menstrual cycle. Sixty-one percent took daily oral contraceptives. Baseline demographics in Angeli 2016 Table I; baseline biological covariates measured at the last cycle visit (end-of-cycle) and used as per-subject reference values in the covariate model."
  )

  ini({
    # ----- Iron turnover parameters (Table III, 'Model with covariates' column) -----
    lksyn_iron  <- log(7.57);       label("Iron zero-order synthesis rate constant (umol/L/d)")  # Angeli 2016 Table III (with covariates); paper symbol ksynI
    lkout_iron  <- log(0.42);       label("Iron first-order elimination rate constant (1/d)")    # Angeli 2016 Table III (with covariates); paper symbol koutI
    lkrel_iron  <- log(2.55);       label("Iron post-menses synthesis increment, active during the rebound window (umol/L/d)")  # Angeli 2016 Table III (with covariates); paper symbol krelI
    lkloss      <- log(0.14);       label("Additional iron AND hepcidin elimination rate during menses (1/d); shared between iron and hepcidin per Methods p. 497")  # Angeli 2016 Table III (with covariates); paper symbol kloss

    # ----- Hepcidin turnover parameters (Table III, 'Model with covariates' column) -----
    lksyn_hep   <- log(2.48);       label("Hepcidin zero-order synthesis rate constant (nmol/L/d)")  # Angeli 2016 Table III (with covariates); paper symbol ksynH
    lkout_hep   <- log(0.82);       label("Hepcidin first-order elimination rate constant (1/d); half-life ~0.85 d")  # Angeli 2016 Table III (with covariates); paper symbol koutH
    lkrel_hep   <- log(0.28);       label("Hepcidin post-menses synthesis increment, active during the rebound window (nmol/L/d)")  # Angeli 2016 Table III (with covariates); paper symbol krelH

    # ----- Iron -> hepcidin coupling (linear multiplicative around iron baseline) -----
    lalpha      <- log(0.03);       label("Linear coupling factor: fractional sensitivity of hepcidin synthesis to (Ir - Ir0)")  # Angeli 2016 Table III (with covariates); paper symbol alpha

    # ----- Rebound window duration -----
    ldrel       <- log(5.74);       label("Duration of the post-menses synthesis-rebound window (days), starting at day 2 of the cycle")  # Angeli 2016 Table III (with covariates); paper symbol drel

    # ----- Covariate effects (Table V; encoded as fixed structural effects per the paper's stepwise selection) -----
    # All continuous-covariate effects are exponents in a power-law multiplier on the rate
    # constant: param = exp(lparam + eta) * (cov / tcov)^e_cov_param.  The reference values
    # tcov for each covariate are documented in covariateData[[<cov>]]$reference_category.
    e_conmed_birthcontrol_kout_iron <- fixed(-0.20);  label("Multiplicative effect of contraception on iron elimination (canonical CONMED_BIRTHCONTROL = 1 means ON contraception; sign flipped from the paper's +0.20 for NOCONTRA)")  # Angeli 2016 Table V; paper reports beta_NOCONTRA = +0.20 for 'No contraception'
    e_bmi_krel_iron                 <- fixed(-3.31);  label("Power exponent for BMI on iron rebound parameter krel_iron; reference tBMI = 22 kg/m^2")            # Angeli 2016 Table V
    e_hgb_krel_iron                 <- fixed( 6.93);  label("Power exponent for baseline hemoglobin on iron rebound parameter krel_iron; reference tHGB = 13.5 g/dL")  # Angeli 2016 Table V
    e_higham_ksyn_hep               <- fixed( 0.66);  label("Power exponent for Higham score on hepcidin synthesis ksyn_hep; reference tHIGHAM = 96.6")           # Angeli 2016 Table V
    e_higham_kout_hep               <- fixed( 0.83);  label("Power exponent for Higham score on hepcidin elimination kout_hep; reference tHIGHAM = 96.6")          # Angeli 2016 Table V
    e_ferritin_bl_kout_hep          <- fixed(-0.60);  label("Power exponent for baseline ferritin on hepcidin elimination kout_hep; reference tFERRITIN = 53 ug/L")  # Angeli 2016 Table V
    e_ht_krel_hep                   <- fixed(32.70);  label("Power exponent for height on hepcidin rebound parameter krel_hep; reference tHT = 165 cm")             # Angeli 2016 Table V
    e_ferritin_bl_krel_hep          <- fixed(-1.95);  label("Power exponent for baseline ferritin on hepcidin rebound parameter krel_hep; reference tFERRITIN = 53 ug/L")  # Angeli 2016 Table V

    # ----- IIV: log-normal, variance on the internal log-scale; omega^2 = log(1 + CV^2) -----
    # CV% from Angeli 2016 Table III 'Variability (RSE)' column for the 'Model with covariates'.
    etalksyn_iron ~ 0.0223   # 15% CV: log(1 + 0.15^2)
    etalkrel_iron ~ 0.2643   # 55% CV: log(1 + 0.55^2)
    etalkloss     ~ 0.6432   # 95% CV: log(1 + 0.95^2)
    etalksyn_hep  ~ 0.0560   # 24% CV: log(1 + 0.24^2)
    etalkout_hep  ~ 0.2814   # 57% CV: log(1 + 0.57^2)
    etalkrel_hep  ~ 0.7232   # 103% CV: log(1 + 1.03^2)
    etalalpha     ~ 0.5636   # 87% CV: log(1 + 0.87^2)
    # Note: koutI IIV was removed from the final model with covariates (0%, Table III);
    # drel had no IIV in any model. Both are intentionally omitted from this block.

    # ----- Residual error (Table III, 'Model with covariates' column) -----
    # Iron uses proportional error only (the additive part was dropped per Methods p. 497;
    # see Table IV model 32, the final selected variant). Hepcidin keeps the combined
    # additive + proportional form (additive part kept for numerical stability per p. 498).
    propSd_Iron <- 0.29;     label("Proportional residual SD for serum iron (fraction)")                      # Angeli 2016 Table III; paper symbol b_Iron
    addSd_Hep   <- 0.17;     label("Additive residual SD for serum hepcidin (nmol/L); kept for numerical stability")  # Angeli 2016 Table III; paper symbol a_Hep
    propSd_Hep  <- 0.47;     label("Proportional residual SD for serum hepcidin (fraction)")                  # Angeli 2016 Table III; paper symbol b_Hep
  })

  model({
    # ----- Reference values for power-law covariate centering -----
    # Set values close to the population median (Angeli 2016 Table I); see
    # covariateData[[<cov>]]$reference_category for the rationale per covariate.
    tBMI      <- 22
    tHGB      <- 13.5
    tHIGHAM   <- 96.6
    tFERRITIN <- 53
    tHT       <- 165

    # ----- Rebound-window start (paper grid search across days 1-5 selected day 2) -----
    tbeg <- 2

    # ----- Typical (baseline) parameter values: combine population means with covariate effects -----
    # The structural baseline rate constants do NOT include the menses-loss or rebound
    # increments; those enter additively below via the schedule indicators.
    ksyn_iron_typ <- exp(lksyn_iron + etalksyn_iron)
    kout_iron_typ <- exp(lkout_iron + e_conmed_birthcontrol_kout_iron * CONMED_BIRTHCONTROL)
    krel_iron_typ <- exp(lkrel_iron + etalkrel_iron) * (BMI / tBMI)^e_bmi_krel_iron * (HGB / tHGB)^e_hgb_krel_iron
    kloss_typ     <- exp(lkloss + etalkloss)

    ksyn_hep_typ <- exp(lksyn_hep + etalksyn_hep) * (HIGHAM / tHIGHAM)^e_higham_ksyn_hep
    kout_hep_typ <- exp(lkout_hep + etalkout_hep) * (HIGHAM / tHIGHAM)^e_higham_kout_hep * (FERRITIN_BL / tFERRITIN)^e_ferritin_bl_kout_hep
    krel_hep_typ <- exp(lkrel_hep + etalkrel_hep) * (HT / tHT)^e_ht_krel_hep * (FERRITIN_BL / tFERRITIN)^e_ferritin_bl_krel_hep

    drel_typ  <- exp(ldrel)
    alpha_typ <- exp(lalpha + etalalpha)

    # ----- Iron baseline (used as the reference point for the iron -> hepcidin coupling) -----
    # Drug-free steady state of the iron equation when neither menses-loss nor rebound is active.
    Ir0 <- ksyn_iron_typ / kout_iron_typ

    # ----- Time-varying schedule (Angeli 2016 Fig. 3 right) -----
    # Menses (loss phase): kloss adds to both kout_iron and kout_hep while 0 <= t < DLOSS.
    # Rebound phase: krel adds to both ksyn while tbeg <= t < tbeg + drel.
    # rxode2 evaluates logical expressions as 1/0 so the indicators multiply cleanly into
    # the additive rate terms.
    loss_ind <- (t < DLOSS)
    rel_ind  <- (t >= tbeg) * (t < tbeg + drel_typ)

    ksyn_iron <- ksyn_iron_typ + krel_iron_typ * rel_ind
    kout_iron <- kout_iron_typ + kloss_typ     * loss_ind
    ksyn_hep  <- ksyn_hep_typ  + krel_hep_typ  * rel_ind
    kout_hep  <- kout_hep_typ  + kloss_typ     * loss_ind

    # ----- ODE system (Angeli 2016 Eq. 4) -----
    # iron: serum iron in umol/L
    # hep:  serum hepcidin in nmol/L
    # Iron drives hepcidin synthesis multiplicatively around the iron baseline Ir0.
    d/dt(iron) <- ksyn_iron - kout_iron * iron
    d/dt(hep)  <- ksyn_hep * (1 + alpha_typ * (iron - Ir0)) - kout_hep * hep

    # ----- Initial conditions (Angeli 2016 Eq. 5): drug-free steady state -----
    iron(0) <- Ir0
    hep(0)  <- ksyn_hep_typ / kout_hep_typ

    # ----- Observations -----
    Iron <- iron
    Hep  <- hep

    Iron ~ prop(propSd_Iron)
    Hep  ~ add(addSd_Hep) + prop(propSd_Hep)
  })
}
