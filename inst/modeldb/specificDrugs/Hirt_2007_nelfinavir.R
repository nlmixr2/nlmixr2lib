Hirt_2007_nelfinavir <- function() {
  description <- paste(
    "Six-compartment population pharmacokinetic model for nelfinavir and",
    "its M8 metabolite describing placental transfer from maternal plasma",
    "into umbilical (cord) plasma and amniotic fluid (Hirt 2007). Oral",
    "nelfinavir is absorbed first-order with a lag time into the maternal",
    "central compartment. Nelfinavir is then (i) eliminated, (ii) converted",
    "to M8 in a maternal M8 compartment, and (iii) transferred to a cord",
    "nelfinavir compartment. M8 is eliminated from the mother and",
    "transferred to a cord M8 compartment. Both nelfinavir and M8 transfer",
    "from cord to amniotic fluid and are eliminated from amniotic fluid by",
    "first-order rate constants. The distribution volume of M8 in the",
    "mother and the volumes of all cord and amniotic-fluid compartments",
    "were not estimable and are fixed at 1 L per the paper. Covariate",
    "effects: day-of-delivery indicator increases maternal nelfinavir CL",
    "and V each by 92 percent and gates a body-weight effect on CL within",
    "the delivery cohort only (reference 73 kg, exponent 2.81); pregnancy",
    "increases M8 elimination by 67 percent; body weight scales M8",
    "elimination on the full database (reference 63 kg, exponent 1.41);",
    "concomitant NNRTI use increases M8 elimination by 148 percent."
  )
  reference <- paste(
    "Hirt D, Urien S, Jullien V, Firtion G, Chappuy H, Rey E, Pons G,",
    "Mandelbrot L, Treluyer JM. (2007). Pharmacokinetic modelling of the",
    "placental transfer of nelfinavir and its M8 metabolite: a population",
    "study using 75 maternal-cord plasma samples.",
    "Br J Clin Pharmacol 64(5):634-644.",
    "doi:10.1111/j.1365-2125.2007.02885.x."
  )
  vignette <- "Hirt_2007_nelfinavir"
  paper_specific_compartments <- c(
    "mother_m8", "cord_n", "cord_m8", "af_n", "af_m8"
  )

  units <- list(
    time          = "hour",
    dosing        = "mg",
    concentration = "mg/L"
  )

  covariateData <- list(
    WT = list(
      description        = paste(
        "Body weight. Used as a continuous covariate through two power",
        "scalings: on the maternal nelfinavir clearance CL_Nm_No within",
        "the day-of-delivery cohort only ((WT/73)^2.81), and on the M8",
        "maternal elimination rate k_M8m_M8o on the full database",
        "((WT/63)^1.41)."
      ),
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Two reference body weights apply: 73 kg for the delivery-group",
        "CL effect (paper Results: '10 kg increase from mean bodyweight",
        "increased CL_Nm_No 1.44 fold'; Table 1 reports mean weight",
        "73 +/- 14 kg for the 75-woman delivery cohort) and 63 kg for the",
        "M8 elimination effect (paper Results: 'for every 10 kg increase",
        "above the mean weight of 63 kg, k_M8m_M8o was increased by",
        "1.23'). Treat as time-fixed at the value recorded at the",
        "sampling occasion."
      ),
      source_name        = "BW"
    ),
    PREG = list(
      description        = paste(
        "Pregnancy indicator. 1 = pregnant (sampled before the day of",
        "delivery); 0 = non-pregnant. In Hirt 2007 the indicator is",
        "explicitly turned off on the day of delivery (paper Methods:",
        "'On the day of delivery, the coding was zero for pregnancy and",
        "one for delivery'), so PREG and DAY_DELIVERY are mutually",
        "exclusive in the source coding."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-pregnant)",
      notes              = paste(
        "Multiplicative effect (1 + theta_PREG * PREG) on the M8",
        "maternal elimination rate k_M8m_M8o, with theta_PREG = 0.67",
        "(67 percent higher M8 elimination in pregnant non-delivery",
        "women relative to non-pregnant women). The pregnancy effect on",
        "CL_Nm_No was screened but deleted during backward elimination",
        "(OFV penalty of 5 units, below the 7-unit retention threshold)",
        "so the final model carries PREG only on k_M8m_M8o."
      ),
      source_name        = "PREG"
    ),
    DAY_DELIVERY = list(
      description        = paste(
        "Day-of-delivery indicator. 1 = sample collected on the day of",
        "delivery; 0 = otherwise. Per the paper's coding the indicator",
        "is mutually exclusive with PREG = 1. Affects maternal nelfinavir",
        "clearance (+92 percent), maternal nelfinavir distribution",
        "volume (+92 percent), and gates a body-weight power effect on",
        "CL within the delivery cohort only."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not day of delivery)",
      notes              = paste(
        "Multiplicative effect (1 + theta_DEL * DAY_DELIVERY) on both",
        "CL_Nm_No and V with the shared theta_DEL = 1.92 (paper Results:",
        "'the same delivery effect on both V and CL_Nm_No, acting in",
        "fact on bioavailability, was chosen because the use of two",
        "different delivery effects did not improved the OFV'). The body",
        "weight power effect (WT/73)^2.81 on CL_Nm_No applies only when",
        "DAY_DELIVERY = 1 (paper Results: 'an effect of bodyweight was",
        "added to CL_Nm_No in the delivery group only'). Ratified",
        "alongside Hirt_2007_nelfinavir.R as a specific-scope canonical;",
        "future placental-transfer / labour-PK papers can ratify general",
        "scope."
      ),
      source_name        = "DEL"
    ),
    CONMED_NNRTI = list(
      description        = paste(
        "Concomitant non-nucleoside reverse-transcriptase inhibitor",
        "(NNRTI) coadministration indicator. 1 = subject is on any",
        "NNRTI (efavirenz, nevirapine, delavirdine, etravirine, or",
        "rilpivirine) during sample collection; 0 = no NNRTI",
        "coadministration. NNRTIs are CYP3A4 / UGT inducers that",
        "increase clearance of co-administered drugs metabolised by",
        "these pathways."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no NNRTI)",
      notes              = paste(
        "Multiplicative effect (1 + theta_NNRTI * CONMED_NNRTI) on the",
        "M8 maternal elimination rate k_M8m_M8o with theta_NNRTI = 1.48",
        "(148 percent higher M8 elimination in NNRTI-exposed women).",
        "Per Hirt 2007 Table 1, 8 of 75 delivery-cohort women and 13 of",
        "121 non-delivery-cohort women received concomitant NNRTI",
        "therapy. The indicator is a class-level CYP-induction flag",
        "rather than a single-drug indicator."
      ),
      source_name        = "NNRTI"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 196L,
    n_studies      = 2L,
    age_range      = "32-34 years (mean 32.7 +/- 4.4 at delivery; 33.4 +/- 4.5 non-delivery cohorts)",
    weight_range   = "approx 50-90 kg (mean 73 +/- 14 kg at delivery; 65 +/- 15 kg non-delivery cohorts)",
    sex_female_pct = 100,
    disease_state  = paste(
      "Adult women with HIV-1 infection receiving oral nelfinavir as part",
      "of antiretroviral therapy. 75 women were sampled on the day of",
      "delivery (gestational age 31-41 weeks, median 38 weeks); 53",
      "pregnant women were sampled before delivery (mean gestation 32",
      "weeks, range 10-39 weeks); 61 non-pregnant women provided",
      "therapeutic-drug-monitoring samples; 7 additional women were",
      "sampled in both pregnant and non-pregnant states (counted once in",
      "the subject total). Four women co-administered with ritonavir were",
      "excluded."
    ),
    dose_range     = paste(
      "750 mg three times daily or 1250 mg twice daily, oral. 92 percent",
      "of delivery-cohort women and 82 percent of non-delivery-cohort",
      "women were on the twice-daily regimen."
    ),
    regions        = "France (Port Royal Hospital, Paris; Louis Mourier Hospital, Colombes)",
    notes          = paste(
      "Demographics from Hirt 2007 Table 1. The model was fit to 292",
      "nelfinavir / M8 samples in maternal plasma (77 at delivery + 215",
      "from pregnant / non-pregnant cohorts), 77 umbilical-plasma",
      "samples, and 27 amniotic-fluid samples. 8 delivery-cohort women",
      "and 13 non-delivery-cohort women were on concomitant NNRTI",
      "therapy. AGE was screened as a covariate during model building",
      "but was not retained in the final model."
    )
  )

  ini({
    # === Structural PK parameters (Hirt 2007 Table 5, original-dataset mean column) ===
    lka            <- log(0.67)   ; label("Absorption rate ka (1/h)")                                # Table 5 row 'ka'
    ltlag          <- log(0.87)   ; label("Absorption lag time tlag (h)")                            # Table 5 row 'tlag'
    lvc            <- log(557)    ; label("Maternal nelfinavir distribution volume V/F (L)")         # Table 5 row 'V'
    lcl            <- log(39.5)   ; label("Maternal nelfinavir oral clearance CL_Nm_No/F (L/h)")     # Table 5 row 'CL_Nm_No/F'
    lcl_m8m        <- log(0.77)   ; label("M8 formation clearance CL_Nm_M8m/F (L/h)")                # Table 5 row 'CL_Nm_M8m/F'
    lk_m8m_out     <- log(3.41)   ; label("Maternal M8 elimination rate k_M8m_M8o (1/h)")            # Table 5 row 'k_M8m_M8o'
    lcl_nc         <- log(0.058)  ; label("Mother-to-cord nelfinavir clearance CL_Nm_Nc/F (L/h)")    # Table 5 row 'CL_Nm_Nc/F'
    lk_nc_naf      <- log(0.23)   ; label("Cord-to-amniotic-fluid nelfinavir rate k_Nc_Naf (1/h)")   # Table 5 row 'k_Nc_Naf'
    lk_naf_out     <- log(0.36)   ; label("Amniotic-fluid nelfinavir elimination rate k_Naf_No (1/h)") # Table 5 row 'k_Naf_No'
    lk_m8m_m8c     <- log(0.35)   ; label("Mother-to-cord M8 rate k_M8m_M8c (1/h)")                  # Table 5 row 'k_M8m_M8c'
    lk_m8c_m8af    <- log(0.59)   ; label("Cord-to-amniotic-fluid M8 rate k_M8c_M8af (1/h)")         # Table 5 row 'k_M8c_M8af'
    lk_m8af_out    <- log(0.49)   ; label("Amniotic-fluid M8 elimination rate k_M8af_M8o (1/h)")     # Table 5 row 'k_M8af_M8o'

    # === Covariate effects ===
    # Day-of-delivery effect on CL_Nm_No and V (shared theta_DEL = 1.92, multiplicative)
    e_day_delivery_cl_vc <- 1.92   ; label("Day-of-delivery additive multiplier shared by CL_Nm_No and V")  # Table 5 row 'CL_Nm_No/F and V, theta_DEL'
    # Body-weight exponent on CL_Nm_No (delivery cohort only): (WT/73)^2.81
    e_wt_cl              <- 2.81   ; label("Body-weight power exponent on CL_Nm_No within the delivery cohort") # Table 5 row 'CL_Nm_No/F, theta_BW (on DEL)'
    # Pregnancy effect on k_M8m_M8o (multiplicative)
    e_preg_k_m8m_out     <- 0.67   ; label("Pregnancy additive multiplier on k_M8m_M8o")             # Table 5 row 'k_M8m_M8o, theta_PREG'
    # Body-weight exponent on k_M8m_M8o (full database): (WT/63)^1.41
    e_wt_k_m8m_out       <- 1.41   ; label("Body-weight power exponent on k_M8m_M8o")                # Table 5 row 'k_M8m_M8o, theta_BW'
    # Concomitant NNRTI effect on k_M8m_M8o (multiplicative)
    e_conmed_nnrti_k_m8m_out <- 1.48 ; label("Concomitant NNRTI additive multiplier on k_M8m_M8o")   # Table 5 row 'k_M8m_M8o, theta_NNRTI'

    # === IIV (exponential; omega^2 = log(CV^2 + 1) from Hirt 2007 Table 5 % CV) ===
    etalvc          ~ 1.0089   # CV 132% on V/F           -> log(1 + 1.32^2)
    etalcl          ~ 0.2476   # CV 53%  on CL_Nm_No/F    -> log(1 + 0.53^2)
    etalk_m8m_out   ~ 0.4368   # CV 74%  on k_M8m_M8o     -> log(1 + 0.74^2)
    etalcl_nc       ~ 0.1094   # CV 34%  on CL_Nm_Nc/F    -> log(1 + 0.34^2)
    etalk_m8m_m8c   ~ 0.4560   # CV 76%  on k_M8m_M8c     -> log(1 + 0.76^2)

    # === Residual error (paper: additive model in linear concentration space) ===
    # Per Hirt 2007 Table 5 labels: 'FETUS' = amniotic-fluid (matrix surrounding the
    # fetus), 'CORD' = umbilical plasma. The nelfinavir 'FETUS/CORD' notation marks a
    # single SD shared by the cord and amniotic-fluid nelfinavir observations
    # (paper data was too sparse to estimate separately, especially for AF nelfinavir
    # with 56 percent of samples below LOQ per Table 2). M8 carried distinct SDs for
    # the two fetal-side matrices. See vignette Assumptions and deviations for the
    # alternative readings considered.
    addSd            <- 1.07 ; label("Mother nelfinavir additive residual SD (mg/L)")               # Table 5 row 's_NELFI_MOTHER'
    addSd_Cm8        <- 0.24 ; label("Mother M8 additive residual SD (mg/L)")                       # Table 5 row 's_M8_MOTHER'
    addSd_Ccord_n    <- 0.09 ; label("Cord nelfinavir additive residual SD (mg/L)")                 # Table 5 row 's_NELFI_FETUS/CORD' (shared with AF nelfinavir)
    addSd_Caf_n      <- 0.09 ; label("Amniotic-fluid nelfinavir additive residual SD (mg/L)")       # Table 5 row 's_NELFI_FETUS/CORD' (shared with cord nelfinavir)
    addSd_Ccord_m8   <- 0.12 ; label("Cord M8 additive residual SD (mg/L)")                         # Table 5 row 's_M8_CORD'
    addSd_Caf_m8     <- 0.03 ; label("Amniotic-fluid M8 additive residual SD (mg/L)")               # Table 5 row 's_M8_FETUS'
  })

  model({
    # === 1. Individual parameters ===
    ka              <- exp(lka)
    tlag            <- exp(ltlag)
    # Day-of-delivery multiplier on CL_Nm_No and V (shared multiplicative effect)
    del_factor      <- 1 + e_day_delivery_cl_vc * DAY_DELIVERY
    # Body-weight power scaling on CL_Nm_No applies only when DAY_DELIVERY = 1.
    # The piecewise form below evaluates to 1 in the non-delivery cohort and to
    # (WT/73)^2.81 in the delivery cohort, matching the paper's
    # 'in the delivery group only' qualifier.
    bw_cl_factor    <- (1 - DAY_DELIVERY) + DAY_DELIVERY * (WT / 73)^e_wt_cl
    # Pregnancy / NNRTI / body-weight effects on M8 elimination
    preg_factor     <- 1 + e_preg_k_m8m_out * PREG
    nnrti_factor    <- 1 + e_conmed_nnrti_k_m8m_out * CONMED_NNRTI
    bw_m8_factor    <- (WT / 63)^e_wt_k_m8m_out

    cl              <- exp(lcl + etalcl) * del_factor * bw_cl_factor
    vc              <- exp(lvc + etalvc) * del_factor
    cl_m8m          <- exp(lcl_m8m)
    k_m8m_out       <- exp(lk_m8m_out + etalk_m8m_out) * preg_factor * bw_m8_factor * nnrti_factor
    cl_nc           <- exp(lcl_nc + etalcl_nc)
    k_nc_naf        <- exp(lk_nc_naf)
    k_naf_out       <- exp(lk_naf_out)
    k_m8m_m8c       <- exp(lk_m8m_m8c + etalk_m8m_m8c)
    k_m8c_m8af      <- exp(lk_m8c_m8af)
    k_m8af_out      <- exp(lk_m8af_out)

    # === 2. ODE system ===
    # Maternal nelfinavir (central, V = vc): three clearances drain the
    # central compartment (elimination, M8 formation, mother-to-cord).
    d/dt(depot)     <- -ka * depot
    d/dt(central)   <-  ka * depot - (cl + cl_m8m + cl_nc) / vc * central

    # Maternal M8 (V_M8m fixed at 1 L per paper). Formation flux from central
    # nelfinavir = (cl_m8m / vc) * central; elimination + transfer to cord M8
    # are first-order rate constants on the amount.
    d/dt(mother_m8) <- (cl_m8m / vc) * central - (k_m8m_out + k_m8m_m8c) * mother_m8

    # Cord (umbilical) nelfinavir (V fixed at 1 L). Influx from maternal
    # central via mother-to-cord clearance; outflux to amniotic fluid via the
    # k_nc_naf rate constant.
    d/dt(cord_n)    <- (cl_nc / vc) * central - k_nc_naf * cord_n

    # Cord M8 (V fixed at 1 L). Influx from maternal M8 via k_m8m_m8c;
    # outflux to amniotic-fluid M8 via k_m8c_m8af.
    d/dt(cord_m8)   <- k_m8m_m8c * mother_m8 - k_m8c_m8af * cord_m8

    # Amniotic-fluid nelfinavir (V fixed at 1 L).
    d/dt(af_n)      <- k_nc_naf * cord_n - k_naf_out * af_n

    # Amniotic-fluid M8 (V fixed at 1 L).
    d/dt(af_m8)     <- k_m8c_m8af * cord_m8 - k_m8af_out * af_m8

    # === 3. Absorption lag ===
    lag(depot)      <- tlag

    # === 4. Observation variables and additive residual errors ===
    Cc              <- central / vc
    Cm8             <- mother_m8
    Ccord_n         <- cord_n
    Ccord_m8        <- cord_m8
    Caf_n           <- af_n
    Caf_m8          <- af_m8

    Cc              ~ add(addSd)
    Cm8             ~ add(addSd_Cm8)
    Ccord_n         ~ add(addSd_Ccord_n)
    Ccord_m8        ~ add(addSd_Ccord_m8)
    Caf_n           ~ add(addSd_Caf_n)
    Caf_m8          ~ add(addSd_Caf_m8)
  })
}
