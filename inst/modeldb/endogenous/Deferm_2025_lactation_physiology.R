Deferm_2025_lactation_physiology <- function() {
  description <- paste(
    "Lactation system-parameter (population physiology) model for the",
    "postpartum period (Deferm 2025 Front Pharmacol). Meta-analysis of",
    "230 studies (36,689 data points from 20,801 postpartum women)",
    "yielding continuous mathematical functions of postpartum age for",
    "13 maternal and breast-milk parameters that govern drug transfer",
    "into human milk: milk volume, milk pH, milk fat, milk total",
    "protein, milk water content, daily infant milk intake, maternal",
    "haematocrit, alpha-1-acid glycoprotein, human serum albumin,",
    "breast volume, plasma volume, cardiac output and glomerular",
    "filtration rate. This is the time-varying SYSTEM layer intended to",
    "be coupled to a lactation PBPK drug model -- it contains no drug,",
    "no dosing, no compartments and no ODEs. Every output is an",
    "algebraic function of the rxode2 time variable, which the model",
    "interprets as postpartum age in MONTHS (paper symbol PpT), so a",
    "solve over time 0 to 12 traces the full first postpartum year.",
    "Between-subject variability is the constant coefficient of",
    "variation the authors applied per parameter to recover the",
    "observed spread (Simcyp virtual-population style), encoded here as",
    "independent lognormal etas. Three parameters are piecewise: milk",
    "volume (sigmoidal to 6 months, mono-exponential decline",
    "thereafter), daily milk intake (sigmoidal to 1 month,",
    "mono-exponential decline thereafter) and alpha-1-acid glycoprotein",
    "(double-exponential to 1 month, linear thereafter). Valid over 0",
    "to 12 months postpartum only.",
    sep = " "
  )
  reference <- paste(
    "Deferm N, Dinh J, Pansari A, Jamei M, Abduljalil K.",
    "Postpartum changes in maternal physiology and milk composition:",
    "a comprehensive database for developing lactation",
    "physiologically-based pharmacokinetic models.",
    "Front Pharmacol. 2025;16:1517069.",
    "doi:10.3389/fphar.2025.1517069.",
    sep = " "
  )
  vignette <- "Deferm_2025_lactation_physiology"
  units <- list(
    time          = "month (postpartum age, PpT; valid 0 to 12)",
    dosing        = "n/a (no exogenous dosing; maternal / milk physiology model)",
    concentration = paste(
      "n/a (no drug concentration). Outputs carry their own units:",
      "milk_volume L/day, milk_ph unitless, milk_fat g/dL,",
      "milk_protein % w/v, milk_water % w/w, milk_intake L/kg/day,",
      "hct %, agp g/L, hsa g/L, breast_volume L, plasma_volume L,",
      "cardiac_output L/h, mgfr mL/min"
    )
  )

  # No covariate columns: every output is a function of postpartum age
  # alone. Deferm 2025 Results 3.1 explicitly states that parity,
  # ethnicity and delivery type were too inconsistently reported to be
  # carried as covariates, and Discussion notes the same for the timing
  # and quantity of complementary-food introduction (milk volume) and
  # for foremilk-versus-hindmilk sampling (milk fat).
  covariatesDataExcluded <- list(
    PARITY = list(
      description = "Number of previous births",
      units       = "(count)",
      type        = "continuous",
      notes       = paste(
        "Screened but not retained: 'as many subject-specific",
        "characteristics, such as parity, ethnicity, or delivery type,",
        "were not consistently reported, these covariates were not",
        "considered in the analysis' (Deferm 2025 Results 3.1). The",
        "Discussion identifies parity as a likely driver of the wide",
        "early-postpartum cardiac-output spread (Ling et al. 2019)."
      )
    ),
    DELIVERY_CAESAREAN = list(
      description = "Caesarean-section delivery indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened but not retained (Deferm 2025 Results 3.1 and section",
        "3.7): 'Since the type of delivery was not consistently",
        "reported across studies, no distinction was made between",
        "vaginal and caesarean deliveries.' The Discussion attributes",
        "part of the early plasma-volume variability to the higher",
        "blood loss of caesarean delivery."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 20801L,
    n_studies      = 230L,
    age_range      = "20.8 to 40 years (weighted mean 28.59 years)",
    weight_range   = "45 to 100.40 kg (weighted mean 63.93 kg)",
    height_range   = "149 to 173.50 cm (weighted mean 163.78 cm)",
    sex_female_pct = 100,
    race_ethnicity = "Not reported / not analysed (Deferm 2025 Results 3.1: ethnicity inconsistently reported across the 230 pooled studies)",
    disease_state  = paste(
      "Healthy breastfeeding women after an uncomplicated, full-term",
      "pregnancy, taking no medication during or after pregnancy",
      "(Deferm 2025 Methods 2.3 inclusion criteria 1-6). Studies",
      "focused on preterm infants, or in which preterm data could",
      "confound the parameter of interest, were excluded; mixed cohorts",
      "were used only where full-term cases represented at least 90% of",
      "the study."
    ),
    dose_range     = "n/a (no drug administered)",
    postpartum_range = "Immediately after childbirth to 12 months postpartum; 60% of data points fall within the first month and only 251 data points after 7 months",
    n_datapoints   = 36689L,
    regions        = "Not restricted (no language or date restriction applied to the PubMed / Google Scholar search)",
    notes          = paste(
      "Per-parameter contributing evidence (Deferm 2025 sections 3.2",
      "to 3.9, Supplementary Table 1): milk volume 11 studies / 312",
      "points / 763 mothers; milk pH 15 / 790 / 328; milk fat 43 /",
      "5,012 / 2,661; milk total protein 8 / 961 / 481; milk water 8 /",
      "1,152 / 1,005; daily milk intake 30 / 164 / 2,417 infants;",
      "haematocrit 16 / 20,186 / 11,004; alpha-1-acid glycoprotein 6 /",
      "145 / 88; human serum albumin 13 / 1,757 / 1,680; breast volume",
      "6 / 226 / 191; plasma volume 22 / 651 / 560; cardiac output 47 /",
      "2,616 / 1,542; GFR 10 / 165 / 165.",
      "Where observed data did not extend to 12 months (haematocrit,",
      "alpha-1-acid glycoprotein, human serum albumin, plasma volume),",
      "the authors anchored the tail with a single value simulated from",
      "1,000 healthy non-pregnant female volunteers aged 18-45 years in",
      "the Simcyp Simulator V23R2, assumed to represent the 12-month",
      "postpartum level. Functions were fitted by weighted least",
      "squares in Excel Solver with each point weighted by study n."
    )
  )

  ini({
    # ---------------------------------------------------------------
    # All coefficients below are published point estimates of fitted
    # regression / sigmoidal functions. Deferm 2025 reports no standard
    # errors, RSEs or confidence intervals for any coefficient, so every
    # value is encoded as fixed().
    # ---------------------------------------------------------------

    # ---- Milk volume, L/day (Deferm 2025 Eq. 5 and Eq. 6) ----
    # 0 to 6 months:  0.81 * PpT^4.37 / (0.1^4.37 + PpT^4.37)
    # 6 to 12 months: 1.619 * exp(-0.116 * PpT)
    milkvol_max <- fixed(0.81)
    label("Maximum (plateau) daily milk production over 0 to 6 months postpartum (L/day)")  # Eq. 5 numerator coefficient
    milkvol_t50 <- fixed(0.1)
    label("Postpartum age at half-maximal milk production (months)")  # Eq. 5 denominator constant 0.1^4.37
    milkvol_hill <- fixed(4.37)
    label("Hill coefficient on postpartum age in the milk-production rise (unitless)")  # Eq. 5 exponent
    milkvol_a <- fixed(1.619)
    label("Intercept of the mono-exponential milk-volume decline after 6 months (L/day)")  # Eq. 6 coefficient
    milkvol_kdecl <- fixed(0.116)
    label("First-order rate of milk-volume decline after 6 months postpartum (1/month)")  # Eq. 6 exponent

    # ---- Milk pH (Deferm 2025 Eq. 7) ----
    # 0.443 * exp(-13.07 * PpT) + 7.167 * exp(0.0023 * PpT)
    milkph_amp <- fixed(0.443)
    label("Amplitude of the fast colostrum-to-mature-milk pH fall (pH units)")  # Eq. 7 first coefficient
    milkph_kfast <- fixed(13.07)
    label("First-order rate of the early milk-pH fall (1/month)")  # Eq. 7 first exponent
    milkph_base <- fixed(7.167)
    label("Mature-milk baseline pH at time zero of the slow term (pH units)")  # Eq. 7 second coefficient
    milkph_kslow <- fixed(0.0023)
    label("Slow first-order rate of the late milk-pH rise (1/month; positive as printed)")  # Eq. 7 second exponent

    # ---- Milk fat, g/dL (Deferm 2025 Eq. 8) ----
    # 3.69 * (1 + 0.012083 * PpT + 0.000171 * PpT^2)
    milkfat_base <- fixed(3.69)
    label("Milk fat concentration at birth, multiplying the polynomial (g/dL)")  # Eq. 8 leading coefficient
    milkfat_lin <- fixed(0.012083)
    label("Linear fractional coefficient on postpartum age for milk fat (1/month)")  # Eq. 8 linear term
    milkfat_quad <- fixed(0.000171)
    label("Quadratic fractional coefficient on postpartum age for milk fat (1/month^2)")  # Eq. 8 quadratic term

    # ---- Milk total protein, % (Deferm 2025 Eq. 9) ----
    # 1.219 + 1.127 * exp(-5.058 * PpT)
    milkprot_ss <- fixed(1.219)
    label("Steady-state (mature-milk) total protein content (% w/v)")  # Eq. 9 intercept
    milkprot_amp <- fixed(1.127)
    label("Colostrum excess total protein above the mature-milk steady state (% w/v)")  # Eq. 9 amplitude
    milkprot_k <- fixed(5.058)
    label("First-order rate of the milk total-protein decline (1/month)")  # Eq. 9 exponent

    # ---- Milk water content, % (Deferm 2025 section 3.2.5) ----
    milkwater_mean <- fixed(87.5)
    label("Milk water content, constant across the postpartum period (% w/w)")  # Section 3.2.5 weighted mean (no time trend; Figure 6)

    # ---- Daily infant milk intake, L/kg/day (Deferm 2025 Eq. 10 and Eq. 11) ----
    # 0 to 1 month:  0.181 * PpT^2.411 / (0.114^2.411 + PpT^2.411)
    # 1 to 12 months: 0.004 + (0.208 - 0.004) * exp(-0.15 * PpT)
    milkintake_max <- fixed(0.181)
    label("Maximum daily milk intake reached during the first postpartum month (L/kg/day)")  # Eq. 10 numerator coefficient
    milkintake_t50 <- fixed(0.114)
    label("Postpartum age at half-maximal daily milk intake (months)")  # Eq. 10 denominator constant 0.114^2.411
    milkintake_hill <- fixed(2.411)
    label("Hill coefficient on postpartum age in the milk-intake rise (unitless)")  # Eq. 10 exponent
    milkintake_ss <- fixed(0.004)
    label("Asymptotic daily milk intake at late postpartum age (L/kg/day)")  # Eq. 11 additive constant
    milkintake_peak <- fixed(0.208)
    label("Extrapolated daily milk intake at time zero of the declining phase (L/kg/day)")  # Eq. 11 term (0.208 - 0.004)
    milkintake_k <- fixed(0.15)
    label("First-order rate of the daily milk-intake decline after 1 month (1/month)")  # Eq. 11 exponent

    # ---- Maternal haematocrit, % (Deferm 2025 Eq. 12) ----
    # 31.17 + (38.74 - 31.17) * PpT^2.49 / (0.133^2.49 + PpT^2.49)
    hct_birth <- fixed(31.17)
    label("Maternal haematocrit immediately after delivery (%)")  # Eq. 12 baseline
    hct_max <- fixed(38.74)
    label("Maternal haematocrit plateau, at pre-pregnancy level (%)")  # Eq. 12 plateau
    hct_t50 <- fixed(0.133)
    label("Postpartum age at half of the haematocrit recovery (months)")  # Eq. 12 denominator constant 0.133^2.49
    hct_hill <- fixed(2.49)
    label("Hill coefficient on postpartum age in the haematocrit recovery (unitless)")  # Eq. 12 exponent

    # ---- Maternal alpha-1-acid glycoprotein, g/L (Deferm 2025 Eq. 13 and Eq. 14) ----
    # 0 to 1 month:  exp(-1.277 * PpT) - exp(-6.749 * PpT) + 0.6
    # 1 to 12 months: -0.016 * PpT + 0.90
    agp_kslow <- fixed(1.277)
    label("Slow first-order rate of the alpha-1-acid glycoprotein double-exponential (1/month)")  # Eq. 13 first exponent
    agp_kfast <- fixed(6.749)
    label("Fast first-order rate of the alpha-1-acid glycoprotein double-exponential (1/month)")  # Eq. 13 second exponent
    agp_base <- fixed(0.6)
    label("Additive offset of the early alpha-1-acid glycoprotein function (g/L)")  # Eq. 13 additive constant
    agp_slope <- fixed(-0.016)
    label("Linear slope of the alpha-1-acid glycoprotein decline from 1 to 12 months (g/L/month)")  # Eq. 14 slope; SIGN CORRECTED -- see the note in model()
    agp_icept <- fixed(0.90)
    label("Intercept of the linear alpha-1-acid glycoprotein decline from 1 to 12 months (g/L)")  # Eq. 14 intercept

    # ---- Maternal human serum albumin, g/L (Deferm 2025 Eq. 15) ----
    # 32.7 + 12.15 / (1 + exp(-7.16 * (PpT - 0.866)))
    hsa_base <- fixed(32.7)
    label("Maternal human serum albumin immediately after delivery (g/L)")  # Eq. 15 baseline
    hsa_amp <- fixed(12.15)
    label("Rise in maternal human serum albumin from delivery to plateau (g/L)")  # Eq. 15 numerator
    hsa_k <- fixed(7.16)
    label("Steepness of the logistic maternal human serum albumin recovery (1/month)")  # Eq. 15 logistic rate
    hsa_t50 <- fixed(0.866)
    label("Postpartum age at half of the human serum albumin recovery (months)")  # Eq. 15 logistic midpoint

    # ---- Breast volume, L (Deferm 2025 Eq. 16) ----
    # 1.549 - 0.024 * PpT
    breastvol_icept <- fixed(1.549)
    label("Empty (milk-free) breast volume at delivery (L)")  # Eq. 16 intercept
    breastvol_slope <- fixed(0.024)
    label("Linear rate of empty-breast-volume decline (L/month)")  # Eq. 16 slope (subtracted in model())

    # ---- Maternal plasma volume, L (Deferm 2025 Eq. 17) ----
    # 2.67 + 0.106 * 0.133^PpT
    plasmavol_ss <- fixed(2.67)
    label("Pre-pregnancy (asymptotic) maternal plasma volume (L)")  # Eq. 17 additive constant
    plasmavol_amp <- fixed(0.106)
    label("Excess maternal plasma volume at delivery above the pre-pregnancy asymptote (L)")  # Eq. 17 amplitude
    plasmavol_decay <- fixed(0.133)
    label("Per-month geometric decay factor of the excess plasma volume (unitless base, raised to PpT)")  # Eq. 17 base of the 0.133^PpT term

    # ---- Maternal cardiac output, L/h (Deferm 2025 Eq. 18) ----
    # 98.8 * exp(-3.33 * PpT) + 304.4 * exp(-0.00096 * PpT)
    cardiacout_amp <- fixed(98.8)
    label("Amplitude of the fast postpartum cardiac-output decline (L/h)")  # Eq. 18 first coefficient
    cardiacout_kfast <- fixed(3.33)
    label("Fast first-order rate of the cardiac-output decline (1/month)")  # Eq. 18 first exponent
    cardiacout_base <- fixed(304.4)
    label("Pre-pregnancy baseline cardiac output at time zero of the slow term (L/h)")  # Eq. 18 second coefficient
    cardiacout_kslow <- fixed(0.00096)
    label("Slow first-order rate of the late cardiac-output drift (1/month)")  # Eq. 18 second exponent

    # ---- Maternal glomerular filtration rate, mL/min (Deferm 2025 Eq. 19) ----
    # 151.0285 - 57.1898*PpT + 17.1856*PpT^2 - 1.8479*PpT^3 + 0.0661*PpT^4
    mgfr_c0 <- fixed(151.0285)
    label("Constant term of the maternal GFR quartic polynomial (mL/min)")  # Eq. 19 intercept
    mgfr_c1 <- fixed(-57.1898)
    label("First-order coefficient of the maternal GFR quartic polynomial (mL/min/month)")  # Eq. 19 linear term
    mgfr_c2 <- fixed(17.1856)
    label("Second-order coefficient of the maternal GFR quartic polynomial (mL/min/month^2)")  # Eq. 19 quadratic term
    mgfr_c3 <- fixed(-1.8479)
    label("Third-order coefficient of the maternal GFR quartic polynomial (mL/min/month^3)")  # Eq. 19 cubic term
    mgfr_c4 <- fixed(0.0661)
    label("Fourth-order coefficient of the maternal GFR quartic polynomial (mL/min/month^4)")  # Eq. 19 quartic term

    # ---------------------------------------------------------------
    # Between-subject variability. Deferm 2025 applies a constant CV per
    # parameter to recover the spread of the observed data (Methods 2.4,
    # and the closing sentence of each Results subsection); the CV is
    # the only variability the paper reports. Encoded as an independent
    # lognormal eta per parameter with omega^2 = log(1 + CV^2), all
    # fixed because the CVs are stated, not estimated with uncertainty.
    # The paper does not report correlations between parameters, so the
    # etas are independent (documented in the vignette Errata).
    # ---------------------------------------------------------------
    etamilkvol ~ fixed(0.1033685)    # CV 33%   -> log(1 + 0.33^2)  = 0.1033685 (section 3.2.1)
    etamilkph ~ fixed(0.0008996)     # CV 3%    -> log(1 + 0.03^2)  = 0.0008996 (section 3.2.2)
    etamilkfat ~ fixed(0.1283053)    # CV 37%   -> log(1 + 0.37^2)  = 0.1283053 (section 3.2.3)
    etamilkprot ~ fixed(0.0515483)   # CV 23%   -> log(1 + 0.23^2)  = 0.0515483 (section 3.2.4)
    etamilkwater ~ fixed(0.0002250)  # CV 1.5%  -> log(1 + 0.015^2) = 0.0002250 (section 3.2.5)
    etamilkintake ~ fixed(0.0606246) # CV 25%   -> log(1 + 0.25^2)  = 0.0606246 (section 3.2.6)
    etahct ~ fixed(0.0063796)        # CV 8%    -> log(1 + 0.08^2)  = 0.0063796 (section 3.3)
    etaagp ~ fixed(0.0560022)        # CV 24%   -> log(1 + 0.24^2)  = 0.0560022 (section 3.4)
    etahsa ~ fixed(0.0099503)        # CV 10%   -> log(1 + 0.10^2)  = 0.0099503 (section 3.5)
    etabreastvol ~ fixed(0.0099503)  # CV 10%   -> log(1 + 0.10^2)  = 0.0099503 (section 3.6)
    etaplasmavol ~ fixed(0.0167588)  # CV 13%   -> log(1 + 0.13^2)  = 0.0167588 (section 3.7)
    etacardiacout ~ fixed(0.0252778) # CV 16%   -> log(1 + 0.16^2)  = 0.0252778 (section 3.8)
    etamgfr ~ fixed(0.1218636)       # CV 36%   -> log(1 + 0.36^2)  = 0.1218636 (section 3.9)
  })

  model({
    # The rxode2 time variable IS postpartum age in months (paper symbol
    # PpT). Solve over 0 to 12 to trace the first postpartum year; the
    # functions are only supported over that window (Deferm 2025 Methods
    # 2.3 inclusion criterion 5, "data recorded up to 12 months
    # postpartum"). The quartic GFR polynomial in particular diverges
    # rapidly outside it.
    ppt <- time

    # ---- Milk volume, L/day (Eq. 5 below 6 months, Eq. 6 above) ----
    # The two branches are continuous at the 6-month knot: Eq. 5 gives
    # 0.810 L/day and Eq. 6 gives 0.807 L/day there.
    tv_milkvol <-
      (ppt < 6) * (milkvol_max * ppt^milkvol_hill /
                     (milkvol_t50^milkvol_hill + ppt^milkvol_hill)) +
      (ppt >= 6) * (milkvol_a * exp(-milkvol_kdecl * ppt))
    milk_volume <- tv_milkvol * exp(etamilkvol)

    # ---- Milk pH (Eq. 7) ----
    # Printed as a "double exponential decay", but the second exponent
    # is positive (+0.0023 /month), giving the slow late rise from about
    # 7.18 at 2 weeks to about 7.37 at 12 months that Figure 3A shows
    # and that section 3.2.2 describes ("It then slowly increased again
    # as mature milk was produced"). Encoded exactly as printed.
    tv_milkph <- milkph_amp * exp(-milkph_kfast * ppt) +
      milkph_base * exp(milkph_kslow * ppt)
    milk_ph <- tv_milkph * exp(etamilkph)

    # ---- Milk fat, g/dL (Eq. 8) ----
    tv_milkfat <- milkfat_base *
      (1 + milkfat_lin * ppt + milkfat_quad * ppt^2)
    milk_fat <- tv_milkfat * exp(etamilkfat)

    # ---- Milk total protein, % (Eq. 9) ----
    tv_milkprot <- milkprot_ss + milkprot_amp * exp(-milkprot_k * ppt)
    milk_protein <- tv_milkprot * exp(etamilkprot)

    # ---- Milk water content, % (section 3.2.5) ----
    # No postpartum-age trend was detected; the weighted mean is used
    # across the whole period (Figure 6).
    tv_milkwater <- milkwater_mean
    milk_water <- tv_milkwater * exp(etamilkwater)

    # ---- Daily infant milk intake, L/kg/day (Eq. 10 below 1 month, Eq. 11 above) ----
    # Continuous at the 1-month knot: Eq. 10 gives 0.1800 and Eq. 11
    # gives 0.1796 L/kg/day there.
    tv_milkintake <-
      (ppt < 1) * (milkintake_max * ppt^milkintake_hill /
                     (milkintake_t50^milkintake_hill + ppt^milkintake_hill)) +
      (ppt >= 1) * (milkintake_ss +
                      (milkintake_peak - milkintake_ss) *
                      exp(-milkintake_k * ppt))
    milk_intake <- tv_milkintake * exp(etamilkintake)

    # ---- Maternal haematocrit, % (Eq. 12) ----
    tv_hct <- hct_birth + (hct_max - hct_birth) * ppt^hct_hill /
      (hct_t50^hct_hill + ppt^hct_hill)
    hct <- tv_hct * exp(etahct)

    # ---- Maternal alpha-1-acid glycoprotein, g/L (Eq. 13 below 1 month, Eq. 14 above) ----
    # SIGN CORRECTION on the Eq. 14 slope. Equation 14 is printed as
    # "m-AGP (g/L) = 0.016 x PpT + 0.90", i.e. with a POSITIVE slope,
    # which rises from 0.92 g/L at 1 month to 1.09 g/L at 12 months.
    # That contradicts the paper in three independent ways:
    #   (1) section 3.4 prose -- "m-AGP levels were assumed to decline
    #       linearly ... returning to pre-pregnancy levels by 12 months";
    #   (2) Figure 9A -- the simulated mean line falls from about
    #       0.88 g/L at 1 month to about 0.71 g/L at 12 months, exactly
    #       reproduced by a slope of -0.016;
    #   (3) continuity with Eq. 13, which gives 0.878 g/L at the 1-month
    #       knot: the negative slope gives 0.884 there, the printed
    #       positive slope 0.916.
    # The negative slope is therefore used (agp_slope = -0.016) and the
    # deviation is recorded in the vignette Errata.
    tv_agp <-
      (ppt < 1) * (exp(-agp_kslow * ppt) - exp(-agp_kfast * ppt) + agp_base) +
      (ppt >= 1) * (agp_slope * ppt + agp_icept)
    agp <- tv_agp * exp(etaagp)

    # ---- Maternal human serum albumin, g/L (Eq. 15) ----
    tv_hsa <- hsa_base + hsa_amp / (1 + exp(-hsa_k * (ppt - hsa_t50)))
    hsa <- tv_hsa * exp(etahsa)

    # ---- Empty (milk-free) breast volume, L (Eq. 16) ----
    # "Empty" = total breast volume minus the estimated milk it holds
    # (section 3.6), so this is breast tissue only.
    tv_breastvol <- breastvol_icept - breastvol_slope * ppt
    breast_volume <- tv_breastvol * exp(etabreastvol)

    # ---- Maternal plasma volume, L (Eq. 17) ----
    # The 0.133 term is a base raised to PpT (a per-month geometric
    # decay of the residual pregnancy plasma-volume expansion), not a
    # linear coefficient: it takes plasma volume from 2.776 L at
    # delivery to the 2.67 L pre-pregnancy asymptote by month 2, which
    # is what section 3.7 describes and Figure 12 shows.
    tv_plasmavol <- plasmavol_ss + plasmavol_amp * plasmavol_decay^ppt
    plasma_volume <- tv_plasmavol * exp(etaplasmavol)

    # ---- Maternal cardiac output, L/h (Eq. 18) ----
    tv_cardiacout <- cardiacout_amp * exp(-cardiacout_kfast * ppt) +
      cardiacout_base * exp(-cardiacout_kslow * ppt)
    cardiac_output <- tv_cardiacout * exp(etacardiacout)

    # ---- Maternal glomerular filtration rate, mL/min (Eq. 19) ----
    # Quartic in postpartum age; supported over 0 to 12 months only.
    tv_mgfr <- mgfr_c0 + mgfr_c1 * ppt + mgfr_c2 * ppt^2 +
      mgfr_c3 * ppt^3 + mgfr_c4 * ppt^4
    mgfr <- tv_mgfr * exp(etamgfr)
  })
}
