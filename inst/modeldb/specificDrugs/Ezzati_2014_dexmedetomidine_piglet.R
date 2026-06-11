Ezzati_2014_dexmedetomidine_piglet <- function() {
  description <- "Preclinical (newborn piglet). One-compartment IV population PK model of dexmedetomidine in a piglet perinatal-asphyxia model with therapeutic hypothermia (Ezzati 2014). Clearance scales allometrically with body weight (Holford exponent 0.75) standardised to 70 kg, decreases with body temperature centred at 37 C (Ftemp), and is multiplied by a paper-specific factor FAED (= 0.558) in the post-hypoxic-ischemic state; volume scales allometrically with weight (exponent 1)."
  reference   <- "Ezzati M, Broad K, Kawano G, Faulkner S, Hassell J, Fleiss B, Gressens P, Fierens I, Rostami J, Maze M, Sleigh JW, Anderson B, Sanders RD, Robertson NJ. Pharmacokinetics of dexmedetomidine combined with therapeutic hypothermia in a piglet asphyxia model. Acta Anaesthesiol Scand. 2014; 58(6):733-742. doi:10.1111/aas.12318"
  vignette    <- "Ezzati_2014_dexmedetomidine_piglet"
  paper_specific_etas <- c("etalfaed", "etaRUV")

  units <- list(time = "h", dosing = "ug", concentration = "ug/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per piglet in Ezzati 2014. Used for a-priori allometric scaling per the paper Methods: CL scales as (WT/70)^0.75 and V scales as (WT/70)^1; reference weight 70 kg (Holford size standard).",
      source_name        = "WT"
    ),
    BODYTEMP = list(
      description        = "Body (rectal) temperature; time-varying across the cooling / rewarming cycle",
      units              = "degC",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Drives the linear-deviation effect on CL: Effect_TEMP = 1 + Ftemp * (BODYTEMP - 37). Reference 37 degC per Ezzati 2014 PDF page 736 equation. Piglet normothermia 38.5 degC; therapeutic-hypothermia target 33.5 degC. Source column `TEMP` is an existing canonical alias of BODYTEMP per Kloprogge 2013 / 2014 precedent.",
      source_name        = "TEMP"
    ),
    HIE_POST = list(
      description        = "Post hypoxic-ischemic event indicator (1 = subject is in the post-insult state from the documented HI event onward; 0 = no HI event)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no HI event or pre-insult)",
      notes              = "Time-varying within a subject as the indicator flips from 0 to 1 at the HI insult time and remains 1 thereafter. In Ezzati 2014 dosing began 0.5 h or 4 h after the HI insult so HIE_POST = 1 for the entire PK observation window of the 9 HI-exposed piglets; 1 control piglet had HIE_POST = 0 throughout (1.5 ug/kg/h infusion without HI and without hypothermia). The paper's source column `AED` encodes a continuous NTP-integral severity score that the final model dichotomises to a 0/1 state flag.",
      source_name        = "AED"
    )
  )

  population <- list(
    species        = "piglet (newborn, male, age <24 h)",
    n_subjects     = 10L,
    n_studies      = 1L,
    age_range      = "<24 h (mean 22.9 h, SD 1.2 h)",
    age_median     = "mean 22.9 h",
    weight_range   = "1.6-2.0 kg",
    weight_median  = "mean 1.76 kg (SD 0.23)",
    sex_female_pct = 0,
    disease_state  = "Perinatal cerebral hypoxia-ischaemia induced by bilateral common-carotid occlusion + FiO2 reduction to 0.09 for 12.5 min, followed by whole-body cooling to 33.5 degC for 18-24 h. 9 of 10 piglets received the HI insult; 1 control piglet had no HI and no cooling.",
    dose_range     = "1 ug/kg IV loading dose over 20 min + maintenance IV infusion 0.6-10 ug/kg/h for 46-48 h",
    regions        = "United Kingdom (single-centre preclinical study).",
    notes          = "Population is a small preclinical cohort (Ezzati 2014 Table 1 + Methods). All piglets male. 1000-replicate bootstrap and prediction-corrected VPC used for evaluation (Methods). Estimation in NONMEM VII (ADVAN1 TRANS2, FOCE-I; Globomax LLC)."
  )

  ini({
    # Structural PK parameters - reference subject: 70 kg pig, 37 degC, no HI insult.
    lcl   <- log(3.52); label("Standardised clearance CLstd at 70 kg (L/h)")              # Ezzati 2014 Table 2: CLstd = 3.52 L/h/70 kg (95% CI 1.35-9.01); typical 1.76 kg piglet -> 0.126 L/kg/h (Abstract)
    lvc   <- log(236);  label("Standardised central volume Vstd at 70 kg (L)")           # Ezzati 2014 Table 2: Vstd = 236 L/70 kg (95% CI 85.8-497); typical 1.76 kg piglet -> 3.37 L/kg (Abstract). Text page 737 reports "(109.1%)" for V CV; treated as a typo (Table 2 and Abstract agree on 191%).

    # Allometric exponents - fixed a-priori per Holford size standard (Ezzati 2014 Methods, PDF page 736).
    e_wt_cl <- fixed(0.75); label("Body-weight allometric exponent on CL (unitless)")    # Ezzati 2014 Methods: PWR = 0.75 for clearance; fixed.
    e_wt_vc <- fixed(1);    label("Body-weight allometric exponent on V (unitless)")     # Ezzati 2014 Methods: PWR = 1 for distribution volumes; fixed.

    # Temperature effect on CL (linear-deviation, reference 37 degC). Effect_TEMP = 1 + Ftemp * (TEMP - 37) per Ezzati 2014 PDF page 736.
    Ftemp <- 0.0934; label("Linear-deviation temperature coefficient on CL (per degC)")  # Ezzati 2014 Table 2: Ftemp = 0.0934 (95% CI 0.0127-0.244); 32.7% CL reduction at 33.5 degC (= 0.0934 * 3.5).

    # Post hypoxic-ischemic clearance factor (paper-specific name FAED) - multiplicative on CL when HIE_POST = 1.
    # Encoded log-transformed so the eta IIV is well-defined on the positive scale; the paper's NONMEM
    # parameterisation was on linear scale (Table 2 95% CI 0.329-1.21 includes 1).
    lfaed <- log(0.558); label("Log of post-HI clearance multiplicative factor FAED (unitless)") # Ezzati 2014 Table 2: FAED = 0.558 (95% CI 0.329-1.21). Paper page 737 prose says "reducing CL by 55.8%" but equation CL = CLstd * (Wt/70)^0.75 * Effect_TEMP * FAED with FAED = 0.558 mathematically means CL is multiplied by 0.558 (i.e. CL reduced TO 55.8% of pre-insult value = 44.2% reduction). See vignette Errata.

    # Inter-individual variability - log-normal omega^2 = log(CV^2 + 1).
    # CL/V correlated block: r = -0.756 (paper Results page 737), cov = r * sqrt(omegaCL^2 * omegaV^2).
    # CL: %BSV = 46.6 -> omega^2 = log(0.466^2 + 1) = 0.1965
    # V : %BSV = 191  -> omega^2 = log(1.91^2  + 1) = 1.5365
    # cov(CL,V) = -0.756 * sqrt(0.1965 * 1.5365) = -0.4154
    etalcl + etalvc ~ c(0.1965,
                        -0.4154, 1.5365)                                                  # Ezzati 2014 Table 2 (CL %BSV 46.6, V %BSV 191) + Results page 737 (correlation -0.756).

    # FAED IIV: Table 2 %BSV = 2.6% -> omega^2 = log(0.026^2 + 1) = 0.000676 (effectively negligible
    # inter-piglet variability on the post-HI factor; preserved per operator instruction).
    etalfaed ~ 0.000676                                                                    # Ezzati 2014 Table 2 (FAED %BSV = 2.6).

    # Karlsson 1998 additional eta on residual variance (eta-on-epsilon). Encoded as a multiplicative
    # log-normal scaling of both additive and proportional residual SDs in model() so each piglet has
    # its own residual-error magnitude: SD_i = SD_pop * exp(etaRUV).
    etaRUV ~ 0.643                                                                          # Ezzati 2014 Table 2 (eta_RUV variance = 0.643); Karlsson, Jonsson, Wiltse & Wade 1998 J Pharmacokinet Biopharm 26:207-246.

    # Residual error - combined additive + proportional on plasma concentration (ug/L).
    addSd  <- 0.171; label("Additive residual error on Cc (ug/L)")                          # Ezzati 2014 Table 2: ERR_ADD = 0.171 ug/L (95% CI 0.0356-0.657).
    propSd <- 0.26;  label("Proportional residual error on Cc (fraction)")                  # Ezzati 2014 Table 2: ERR_PROP = 26% (95% CI 4.5-48.9).
  })
  model({
    # Individual PK parameters with a-priori allometric weight scaling (reference 70 kg).
    # Temperature: Effect_TEMP = 1 + Ftemp * (TEMP - 37) -- linear deviation centered at 37 degC.
    # Post-HI: FAED applies only when HIE_POST = 1; the (1 + (faed - 1) * HIE_POST) form gives
    # factor = 1 when HIE_POST = 0 and factor = faed when HIE_POST = 1.
    effect_temp <- 1 + Ftemp * (BODYTEMP - 37)
    faed        <- exp(lfaed + etalfaed)
    faed_effect <- 1 + (faed - 1) * HIE_POST

    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * effect_temp * faed_effect
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc

    kel <- cl / vc

    # One-compartment IV disposition (Ezzati 2014 used ADVAN1 TRANS2; no depot).
    d/dt(central) <- -kel * central

    Cc <- central / vc

    # Combined residual error with Karlsson 1998 additional eta-on-epsilon: each piglet's
    # residual-error magnitude is scaled by exp(etaRUV). The add() / prop() error helpers
    # accept only a plain SD name, so the eta-scaled magnitudes are materialised into
    # derived variables `addSd_i` and `propSd_i` first.
    err_scale <- exp(etaRUV)
    addSd_i   <- addSd  * err_scale
    propSd_i  <- propSd * err_scale
    Cc ~ add(addSd_i) + prop(propSd_i)
  })
}
