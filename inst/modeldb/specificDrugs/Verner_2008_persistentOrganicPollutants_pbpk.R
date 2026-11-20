Verner_2008_persistentOrganicPollutants_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, lifetime; original implementation in acslX).",
    "Lifetime (0-55 years) toxicokinetics of persistent organic pollutants",
    "(POPs) in women, for breast-cancer exposure assessment (Verner et al.",
    "2008). Flow-limited mass-balance model: a woman is a network of tissue",
    "compartments (liver, adipose, richly perfused, slowly perfused, skin,",
    "mammary, uterine and a pregnancy conceptus) perfused by blood, with the",
    "liver carrying both oral intake and hepatic metabolism and the mammary",
    "tissue carrying lactational excretion. Model time is HOURS SINCE BIRTH,",
    "so age (years) = time/8760 drives every organ volume and blood flow",
    "through the age/body-weight/body-height regressions of Haddad et al.",
    "2006 and Price PS et al. 2003 (Supplemental Table S1); pregnancy adds",
    "uterine, mammary, adipose, placental and fetal growth as functions of a",
    "pregnancy sub-clock, and the conceptus burden leaves the mother at",
    "parturition. Tissue:blood partition coefficients are computed from each",
    "compound's log Kow and the tissue lipid/water composition (Supplemental",
    "Table S3). Three POPs are parameterised by (Kow, half-life): PCB-153",
    "(default), PCB-180 and hexachlorobenzene (HCB). Exposure is a constant",
    "daily oral intake (ng/kg/day) into the liver plus optional contaminated",
    "breast-milk intake in infancy; there are no dose records. Deterministic",
    "typical-value model (no IIV, no residual error): the paper is a forward",
    "simulation tool that was not fitted to individual data, and reports no",
    "random effects. Body burden at birth is set to zero."
  )
  reference <- paste(
    "Verner MA, Charbonneau M, Lopez-Carrillo L, Haddad S. (2008).",
    "Physiologically based pharmacokinetic modeling of persistent organic",
    "pollutants for lifetime exposure assessment: a new tool in breast",
    "cancer epidemiologic studies. Environ Health Perspect 116(7):886-892.",
    "doi:10.1289/ehp.10917. PMCID PMC2453156. Model equations from the main",
    "text (Eqs 1-12); age/body-weight/body-height physiology regressions,",
    "milk-volume/composition equations and tissue composition from the",
    "Supplemental Material (Tables S1-S3). Physiology regressions originate",
    "in Haddad S, Tardif GC, Tardif R. (2006) J Toxicol Environ Health A",
    "69:2095-2136; Price K, Haddad S, Krishnan K. (2003) J Toxicol Environ",
    "Health A 66:417-433; and the pregnancy equations in Gentry PR et al.",
    "(2002) Regul Toxicol Pharmacol 36:51-68; milk equations from Neville MC",
    "et al. (1991) Am J Clin Nutr 54:81-92."
  )
  vignette <- "Verner_2008_persistent_organic_pollutants"

  units <- list(
    time = "h",
    dosing = "ug",
    concentration = "ug/L",
    amount = "ug",
    weight = "kg"
  )

  paper_specific_compartments <- c(
    "liver",
    "fat",
    "rich",
    "slow",
    "skin",
    "mammary",
    "uterine",
    "conceptus",
    "a_intake",
    "a_elim",
    "auc_blood"
  )

  compartmentData <- list(
    liver = list(analyte = "POP", units = "ug", specimen = "tissue", verified = TRUE),
    fat = list(analyte = "POP", units = "ug", specimen = "tissue", verified = TRUE),
    rich = list(analyte = "POP", units = "ug", specimen = "tissue", verified = TRUE),
    slow = list(analyte = "POP", units = "ug", specimen = "tissue", verified = TRUE),
    skin = list(analyte = "POP", units = "ug", specimen = "tissue", verified = TRUE),
    mammary = list(analyte = "POP", units = "ug", specimen = "tissue", verified = TRUE),
    uterine = list(analyte = "POP", units = "ug", specimen = "tissue", verified = TRUE),
    conceptus = list(analyte = "POP", units = "ug", specimen = "tissue", verified = TRUE),
    a_intake = list(analyte = "POP", units = "ug", specimen = "not applicable", verified = TRUE),
    a_elim = list(analyte = "POP", units = "ug", specimen = "not applicable", verified = TRUE),
    auc_blood = list(analyte = "POP", units = "ug*h/L", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 0L,
    n_studies = 0L,
    age_range = "0-55 years (lifetime simulation of a woman)",
    weight_range = paste(
      "Driven by a lifetime body-weight profile, not a fixed value. The",
      "'normal' profile rises along a 50th-percentile female growth curve to",
      "50 kg at 14 years, then linearly to 70 kg at 25 years and holds (Figure",
      "2). Overweight (90 kg) and weight-loss (90 kg then back to 70 kg over",
      "25-35 years) profiles are selectable through WT_PEAK / WT_FINAL."
    ),
    sex_female_pct = 100,
    disease_state = "healthy women; lifetime environmental POP exposure",
    dose_range = paste(
      "No dose records. Constant background oral intake, default 10 ng/kg/day",
      "into the liver (EXPOSURE), optionally preceded by contaminated",
      "breast-milk intake in infancy (CMILK default 2 ug/L for BF_INFANT_MO",
      "months)."
    ),
    regions = "Generic (model framework); illustrative simulations not region-specific.",
    notes = paste(
      "NO SUBJECTS WERE FITTED. The paper is a forward PBPK simulation tool;",
      "'model validation could not be achieved because of the evident lack of",
      "data'. Three compounds were simulated by changing log Kow and half-life:",
      "PCB-153 (log Kow 6.72, t1/2 27.5 y), PCB-180 (7.21, 9.9 y) and HCB",
      "(5.73, 6.0 y)."
    )
  )

  ini({
    # =====================================================================
    # Everything here is a recoded literature constant: nothing was estimated,
    # so every parameter is fixed(). Sources: Verner et al. 2008 main text
    # (Eqs 1-12) and Supplemental Tables S1-S3, unless noted. Partition
    # coefficients are not parameters here -- they are computed in model()
    # from KOW and the tissue composition fractions below (Eq 6).
    # =====================================================================

    # --- Chemical-specific properties (default = PCB-153) ----------------
    # Override KOW and HL_YR through rxSolve(params=) to switch compound:
    # PCB-180  KOW = 10^7.21 = 1.6218e7, HL_YR = 9.9
    # HCB      KOW = 10^5.73 = 5.3703e5, HL_YR = 6.0
    lkow <- fixed(log(5.2481e6)) ; label("Natural log of the n-octanol:water partition coefficient Kow (unitless)") # text: log10 Kow 6.72 for PCB-153 (ATSDR); Kow = 10^6.72 = 5.2481e6
    half_life_yr <- fixed(27.5) ; label("Elimination half-life (years)")                                        # text: 27.5 y for PCB-153 (ATSDR; To-Figueras 2000)

    # --- Tissue composition, lipid/water fractions (Table S3) ------------
    # Used in Eq 6 for tissue:blood partition coefficients. Fat and
    # muscle/skin differ between children (<14 y) and adults (>=18 y), with a
    # linear change over 14-18 y (Table S3); the age interpolation is in
    # model(). Richly perfused composition is used for liver, mammary,
    # uterine and conceptus (no distinct values are tabulated).
    fl_blood  <- fixed(0.006) ; label("Lipid fraction in blood (unitless)")                 # Table S3 blood Fl 0.6%
    fw_blood  <- fixed(0.790) ; label("Water fraction in blood (unitless)")                 # Table S3 blood Fw 79.0%
    fl_fat_c  <- fixed(0.550) ; label("Lipid fraction in adipose tissue, child (unitless)") # Table S3 adipose 0-14 y Fl 55.0%
    fw_fat_c  <- fixed(0.411) ; label("Water fraction in adipose tissue, child (unitless)") # Table S3 adipose 0-14 y Fw 41.1%
    fl_fat_a  <- fixed(0.741) ; label("Lipid fraction in adipose tissue, adult (unitless)") # Table S3 adipose 18-55 y Fl 74.1%
    fw_fat_a  <- fixed(0.212) ; label("Water fraction in adipose tissue, adult (unitless)") # Table S3 adipose 18-55 y Fw 21.2%
    fl_musc_c <- fixed(0.021) ; label("Lipid fraction in muscle/skin, child (unitless)")    # Table S3 muscles and skin 0-14 y Fl 2.1%
    fl_musc_a <- fixed(0.042) ; label("Lipid fraction in muscle/skin, adult (unitless)")    # Table S3 muscles and skin 18-55 y Fl 4.2%
    fw_musc   <- fixed(0.741) ; label("Water fraction in muscle/skin (unitless)")           # Table S3 muscles and skin Fw 74.1%
    fl_rich   <- fixed(0.0368) ; label("Lipid fraction in richly perfused tissue (unitless)") # Table S3 richly perfused Fl 3.68%
    fw_rich   <- fixed(0.781) ; label("Water fraction in richly perfused tissue (unitless)")  # Table S3 richly perfused Fw 78.1%

    # --- Milk composition (Table S2), used for milk:blood partition ------
    fl_milk_0   <- fixed(3.8)    ; label("Milk lipid content intercept (%)")   # Table S2 Fl_milk = 3.8 + 0.0095*Days_P
    fl_milk_slp <- fixed(0.0095) ; label("Milk lipid content slope (%/day postpartum)") # Table S2
    fp_milk_0   <- fixed(0.8)    ; label("Milk protein content intercept (%)") # Table S2 Fp_milk = 0.8 + 0.0004*Days_P
    fp_milk_slp <- fixed(0.0004) ; label("Milk protein content slope (%/day postpartum)") # Table S2
    vmilk_0     <- fixed(1.069)    ; label("Milk volume intercept (L/day)")   # Table S2 Vmilk = (1.069 - 0.001212*Days_P)/24 L/h
    vmilk_slp   <- fixed(0.001212) ; label("Milk volume slope (L/day per day postpartum)") # Table S2

    # --- Exposure and scenario knobs -------------------------------------
    EXPOSURE      <- fixed(10)  ; label("Background oral intake (ng/kg body weight/day)")  # text: background daily exposure 10 ng/kg/day
    CMILK         <- fixed(2)   ; label("Contaminated breast-milk concentration drunk in infancy (ug/L)") # text: milk concentration 2 ug/L
    BF_INFANT_MO  <- fixed(0)   ; label("Breast-feeding duration in infancy (months; 0 = not breast-fed)") # text: 6-month infancy breast-feeding scenario
    AGE_PREG_YR   <- fixed(999) ; label("Age at start of pregnancy (years; >=55 = no pregnancy)")          # text: pregnancies simulated at 20/30/35/40 y
    GEST_WK       <- fixed(40)  ; label("Gestation length (weeks)")                                        # standard term gestation
    LACT_MO       <- fixed(0)   ; label("Lactation duration after birth (months; 0 = no lactation)")       # text: 6- and 12-month lactation scenarios

    # --- Body-weight / body-height lifetime profile ----------------------
    WT_14    <- fixed(50) ; label("Body weight at 14 years (kg)")                      # text/Figure 2: 50 kg at 14 y
    WT_PEAK  <- fixed(70) ; label("Adult peak body weight, reached at 25 years (kg)")  # text/Figure 2: 70 kg (normal), 90 kg (overweight)
    WT_FINAL <- fixed(70) ; label("Final body weight from 35 years (kg)")              # text/Figure 2: weight-loss drops 90->70 over 25-35 y
    HT_ADULT <- fixed(165) ; label("Adult body height (cm)")                           # Figure 2A / Figure S2 (50th-percentile female)
    WT_REF   <- fixed(70) ; label("Reference adult body weight for the half-life/metabolism calibration (kg)") # typical adult at which ATSDR half-lives were measured

    # =====================================================================
    # No IIV and no residual error (deterministic forward-simulation tool).
    # =====================================================================
  })

  model({
    # =====================================================================
    # 1. Master clock and chemical-/compound-level constants.
    # =====================================================================
    age <- time / 8760                 # age in years (8760 h/year)
    kow <- exp(lkow)                   # n-octanol:water partition coefficient

    # Age-interpolated lipid/water fractions for adipose and muscle/skin
    # (Table S3: child <14 y, adult >=18 y, linear over 14-18 y).
    frac_adult <- min(1, max(0, (age - 14) / (18 - 14)))
    fl_fat <- fl_fat_c + (fl_fat_a - fl_fat_c) * frac_adult
    fw_fat <- fw_fat_c + (fw_fat_a - fw_fat_c) * frac_adult
    fl_musc <- fl_musc_c + (fl_musc_a - fl_musc_c) * frac_adult

    # Tissue:blood partition coefficients, Eq 6: P = (Kow*Fl_t + Fw_t)/(Kow*Fl_b + Fw_b).
    denom_b <- kow * fl_blood + fw_blood
    p_fat  <- (kow * fl_fat + fw_fat) / denom_b
    p_slow <- (kow * fl_musc + fw_musc) / denom_b
    p_skin <- p_slow                                   # skin uses muscle/skin composition (Table S3)
    p_rich <- (kow * fl_rich + fw_rich) / denom_b
    p_liver <- p_rich                                  # liver uses richly perfused composition (no distinct Table S3 value)
    p_mam <- p_rich                                    # mammary uses richly perfused PC (stated in Methods)
    p_ut <- p_rich                                     # uterine uses richly perfused PC (stated in Methods)
    p_concep <- p_rich                                 # placenta/fetus use richly perfused PC (stated in Methods)

    # =====================================================================
    # 2. Lifetime body weight BW(age) and body height BH(age).
    #    BH: 50th-percentile female growth to HT_ADULT by ~18 y, then held
    #    (digitised from Figure 2A / Figure S2). BW: 50th-percentile female
    #    growth to WT_14 by 14 y, then the Figure 2 adult profile (linear to
    #    WT_PEAK at 25 y, linear to WT_FINAL at 35 y, held thereafter).
    # =====================================================================
    bh <- 50
    bh <- bh + (75 - 50)   * min(1, max(0, (age - 0) / (1 - 0)))
    bh <- bh + (87 - 75)   * min(1, max(0, (age - 1) / (2 - 1)))
    bh <- bh + (109 - 87)  * min(1, max(0, (age - 2) / (5 - 2)))
    bh <- bh + (138 - 109) * min(1, max(0, (age - 5) / (10 - 5)))
    bh <- bh + (160 - 138) * min(1, max(0, (age - 10) / (14 - 10)))
    bh <- bh + (HT_ADULT - 160) * min(1, max(0, (age - 14) / (18 - 14)))

    bw <- 3.4
    bw <- bw + (9.5 - 3.4)   * min(1, max(0, (age - 0) / (1 - 0)))
    bw <- bw + (12 - 9.5)    * min(1, max(0, (age - 1) / (2 - 1)))
    bw <- bw + (18 - 12)     * min(1, max(0, (age - 2) / (5 - 2)))
    bw <- bw + (26 - 18)     * min(1, max(0, (age - 5) / (8 - 5)))
    bw <- bw + (32 - 26)     * min(1, max(0, (age - 8) / (10 - 8)))
    bw <- bw + (41 - 32)     * min(1, max(0, (age - 10) / (12 - 10)))
    bw <- bw + (WT_14 - 41)  * min(1, max(0, (age - 12) / (14 - 12)))
    bw <- bw + (WT_PEAK - WT_14)   * min(1, max(0, (age - 14) / (25 - 14)))
    bw <- bw + (WT_FINAL - WT_PEAK) * min(1, max(0, (age - 25) / (35 - 25)))

    # =====================================================================
    # 3. Pregnancy sub-clocks (hours since start of pregnancy / since birth).
    # =====================================================================
    t_preg0 <- AGE_PREG_YR * 8760
    gest_h <- GEST_WK * 168
    t_birth <- t_preg0 + gest_h
    preg_on <- (time >= t_preg0) * (time < t_birth)           # pregnant now
    postp_on <- (time >= t_birth) * (time < t_birth + 4380)   # <6 months postpartum
    hours_s <- min(gest_h, max(0, time - t_preg0))            # hours since start of pregnancy (clamped)
    hours_p <- max(0, time - t_birth)                         # hours postpartum
    days_p <- hours_p / 24

    # =====================================================================
    # 4. Basal (non-pregnant) organ volumes (L), Table S1. age-piecewise.
    # =====================================================================
    s_bsa <- bw^0.5150 * bh^0.4220 * 234.9       # body surface area (cm^2)
    v_liver <- 0.0501 * bw^0.780

    rp_young <- (-1.919e-2 * age + 3.193 * (bw^2 / bh)^0.2657 - 1.374) - v_liver
    rp_teen  <- (2.515e-2 * age + 7.619 * (bw^2 / bh)^0.1499 - 6.098) - v_liver
    rp_adult <- (2.331e-3 * age + 0.1253 * bw^0.8477 + bh^0.3821 - 4.725) - v_liver
    v_rich <- rp_young
    v_rich <- v_rich + (rp_teen - rp_young) * (age >= 3)
    v_rich <- v_rich + (rp_adult - rp_teen) * (age >= 18)

    v_tongue <- 1.190e-3 * bw - 4.302e-4
    v_heart <- 1.017e-7 * (bh^0.6862 * bw^0.3561 * 242.7)^1.420
    skm_young <- 9.563e-2 * bw + 1.650e-2 * bh + 9.102e-2 * age - 1.642e-1
    skm_teen  <- 1.629e-1 * bw + 2.603e-2 * bh + 4.661e-1 * age - 3.332
    skm_adult <- 6.780 * (s_bsa / 1e4)^1.629 - 1.492e-3 * age + 3.580
    v_skm <- skm_young
    v_skm <- v_skm + (skm_teen - skm_young) * (age >= 3)
    v_skm <- v_skm + (skm_adult - skm_teen) * (age >= 18)
    v_slow <- v_tongue + v_heart + v_skm

    der_young <- 0.664 * (s_bsa / 1e4)
    der_teen  <- -9.356e-5 - 2.151e-5 * age - 5.058e-1 * (s_bsa / 1e4) +
      1.134e-6 * age^2 + 0.117 * age * (s_bsa / 1e4) - 1.673e-5 * (s_bsa / 1e4)^2
    der_adult <- 1.834 * (s_bsa / 1e4)
    v_dermis <- der_young
    v_dermis <- v_dermis + (der_teen - der_young) * (age >= 10)
    v_dermis <- v_dermis + (der_adult - der_teen) * (age >= 20)
    v_epidermis <- 7.850e-2 * (s_bsa / 1e4)^1.049
    v_skin <- v_dermis + v_epidermis

    v_fat_basal <- 0.91 * bw - (v_liver + v_rich + v_slow + v_skin)
    v_mam_basal <- 0.0062 * bw
    v_ut_basal <- 0.0014 * bw

    # =====================================================================
    # 5. Pregnancy volume changes (Table S1), driven by hours_s / hours_p.
    #    Increase during gestation; linear 6-month (4380 h) return postpartum.
    # =====================================================================
    vf_incr <- bw * 0.09 * exp(-12.90995862 * exp(-0.000797 * hours_s))
    vf_endp <- bw * 0.084026
    vmam_incr <- bw * 0.0065 * exp(-7.44868477 * exp(-0.000678 * hours_s))
    vmam_endp <- bw * 0.007088
    vu_incr <- bw * 0.02 * exp(-4.715669973 * exp(-0.000376 * hours_s))
    vu_endp <- bw * 0.029799

    v_fat <- v_fat_basal + vf_incr * preg_on + vf_endp * (1 - hours_p / 4380) * postp_on
    v_mam <- v_mam_basal + vmam_incr * preg_on + vmam_endp * (1 - hours_p / 4380) * postp_on
    v_ut <- v_ut_basal + vu_incr * preg_on + vu_endp * (1 - hours_p / 4380) * postp_on

    # Conceptus (placenta + fetus) volume (L); held at its birth value after
    # parturition (hours_s clamped at gest_h) so the burden stays trapped.
    v_pla <- 0.85 * exp(-9.434 * exp(-5.23e-4 * hours_s))
    v_fet <- 3.50 * (exp(-16.081 * exp(-5.67e-4 * hours_s)) +
      exp(-140.178 * exp(-7.01e-4 * hours_s)))
    v_concep <- v_pla + v_fet

    # =====================================================================
    # 6. Blood flows (L/h), Table S1.
    # =====================================================================
    qc_basal <- 15.048 * bw^0.7609
    q_liver <- 60.00 * v_liver
    q_slow <- 1.80 * (v_tongue + v_skm) + 57.60 * v_heart
    q_skin <- 9.0 * v_skin
    qf_basal <- 1.80 * v_fat_basal
    q_fat <- qf_basal * v_fat / v_fat_basal
    qmam_basal <- 0.027 * qc_basal
    q_mam <- qmam_basal * v_mam / v_mam_basal
    qu_basal <- 0.0062 * qc_basal
    q_ut <- qu_basal * v_ut / v_ut_basal
    # Placental/conceptus flow is not separately tabulated; set to 3% of basal
    # cardiac output while pregnant (same order as the mammary fraction 0.027).
    # Maternal blood concentration is insensitive to it -- the paper assigns
    # the placental diffusion constant an arbitrary value of 1 for that reason.
    q_concep <- 0.03 * qc_basal * preg_on
    q_rich <- qc_basal - (q_slow + q_skin + qf_basal + q_liver + qu_basal + qmam_basal)
    qc <- qc_basal + (q_fat - qf_basal) + (q_mam - qmam_basal) +
      (q_ut - qu_basal) + q_concep

    # =====================================================================
    # 7. Hepatic metabolism (Eqs 8-12). The age-invariant intrinsic clearance
    #    per litre of liver, CLint_C, is fixed from the adult reference
    #    physiology at which the half-life was measured; the extraction ratio
    #    is then recomputed at the current liver volume and flow.
    # =====================================================================
    # Adult reference volumes (WT_REF kg, HT_ADULT, adult branch of Table S1).
    # WT_REF is the half-life-measurement reference adult and is held fixed so
    # the age-invariant intrinsic clearance does not depend on the simulated
    # body-weight scenario.
    bw_r <- WT_REF
    bh_r <- HT_ADULT
    s_r <- bw_r^0.5150 * bh_r^0.4220 * 234.9
    vl_r <- 0.0501 * bw_r^0.780
    vrich_r <- (2.331e-3 * 30 + 0.1253 * bw_r^0.8477 + bh_r^0.3821 - 4.725) - vl_r
    vskm_r <- 6.780 * (s_r / 1e4)^1.629 - 1.492e-3 * 30 + 3.580
    vslow_r <- (1.190e-3 * bw_r - 4.302e-4) +
      1.017e-7 * (bh_r^0.6862 * bw_r^0.3561 * 242.7)^1.420 + vskm_r
    vskin_r <- 1.834 * (s_r / 1e4) + 7.850e-2 * (s_r / 1e4)^1.049
    vfat_r <- 0.91 * bw_r - (vl_r + vrich_r + vslow_r + vskin_r)
    vmam_r <- 0.0062 * bw_r
    vut_r <- 0.0014 * bw_r
    ql_r <- 60.00 * vl_r
    # Adult partition coefficients (Eq 6, adult composition).
    pfat_r <- (kow * fl_fat_a + fw_fat_a) / denom_b
    pslow_r <- (kow * fl_musc_a + fw_musc) / denom_b
    # Volume of distribution, Eq 11 (the blood term Vb is <0.1% of Vd for
    # log Kow > 5 and is dropped).
    vd_r <- pfat_r * vfat_r + p_rich * (vrich_r + vl_r + vmam_r + vut_r) +
      pslow_r * vslow_r + pslow_r * vskin_r
    cl_r <- (0.6931472 / (half_life_yr * 8760)) * vd_r   # Eq 10, ln(2)/HL(h) * Vd
    eh_r <- cl_r / ql_r                                  # Eq 9
    clint_c <- (eh_r * ql_r / (1 - eh_r)) / vl_r         # Eq 8, age-invariant
    eh <- clint_c * v_liver / (clint_c * v_liver + q_liver)  # Eq 12, current age

    # =====================================================================
    # 8. Concentrations (ug/L) and quasi-steady arterial blood (Eq 1 treats
    #    Ca as the common arterial input). Blood storage is negligible on the
    #    years timescale of POP kinetics, so Ca is the flow-weighted mixed
    #    venous concentration Sum(Qt*Ct/Pt)/Qc.
    # =====================================================================
    cv_liver <- (liver / v_liver) / p_liver
    cv_fat <- (fat / v_fat) / p_fat
    cv_rich <- (rich / v_rich) / p_rich
    cv_slow <- (slow / v_slow) / p_slow
    cv_skin <- (skin / v_skin) / p_skin
    cv_mam <- (mammary / v_mam) / p_mam
    cv_ut <- (uterine / v_ut) / p_ut
    cv_concep <- (conceptus / (v_concep + 1e-8)) / p_concep

    ca <- (q_liver * cv_liver + q_fat * cv_fat + q_rich * cv_rich +
      q_slow * cv_slow + q_skin * cv_skin + q_mam * cv_mam +
      q_ut * cv_ut + q_concep * cv_concep) / qc

    # =====================================================================
    # 9. Exposure fluxes (ug/h).
    # =====================================================================
    # Background oral intake: EXPOSURE ng/kg/day * bw / 24 h, ng -> ug (/1000).
    intake_oral <- EXPOSURE * bw / 24 / 1000
    # Infant breast-milk intake: milk volume (Neville, Table S2) * CMILK,
    # during the first BF_INFANT_MO months of life.
    bf_on <- (age < BF_INFANT_MO / 12)
    vmilk_inf <- (vmilk_0 - vmilk_slp * (time / 24)) / 24      # L/h, Days_P = age in days
    intake_milk <- vmilk_inf * CMILK * bf_on

    # Hepatic metabolism, Eq 2: RAM = Ql * Eh * Ca.
    ram <- q_liver * eh * ca

    # Lactational excretion, Eqs 4-5: RAE = Qmilk * (Cmam/Pmam) * Pmilk.
    lact_on <- (time >= t_birth) * (time < t_birth + LACT_MO * 730.5)
    fl_milk <- (fl_milk_0 + fl_milk_slp * days_p) / 100
    fp_milk <- (fp_milk_0 + fp_milk_slp * days_p) / 100
    fw_milk <- 1 - fl_milk - fp_milk
    p_milk <- (kow * fl_milk + fw_milk) / denom_b
    q_milk <- (vmilk_0 - vmilk_slp * days_p) / 24
    rae <- q_milk * (mammary / v_mam / p_mam) * p_milk * lact_on

    # =====================================================================
    # 10. Mass-balance ODE system (Eqs 1-5). Flow-limited: for each tissue
    #     d/dt(At) = Qt*(Ca - Ct/Pt), with the liver carrying intake and
    #     metabolism and the mammary tissue carrying lactational excretion.
    # =====================================================================
    d/dt(liver) <- q_liver * (ca - cv_liver) - ram + intake_oral + intake_milk
    d/dt(fat) <- q_fat * (ca - cv_fat)
    d/dt(rich) <- q_rich * (ca - cv_rich)
    d/dt(slow) <- q_slow * (ca - cv_slow)
    d/dt(skin) <- q_skin * (ca - cv_skin)
    d/dt(mammary) <- q_mam * (ca - cv_mam) - rae
    d/dt(uterine) <- q_ut * (ca - cv_ut)
    d/dt(conceptus) <- q_concep * (ca - cv_concep)

    d/dt(a_intake) <- intake_oral + intake_milk   # cumulative intake
    d/dt(a_elim) <- ram + rae                      # cumulative metabolism + lactation
    d/dt(auc_blood) <- ca                          # blood-concentration AUC (independent PKNCA cross-check)

    # Body burden at birth is set to 0 (Discussion).
    liver(0) <- 0
    fat(0) <- 0
    rich(0) <- 0
    slow(0) <- 0
    skin(0) <- 0
    mammary(0) <- 0
    uterine(0) <- 0
    conceptus(0) <- 0

    # =====================================================================
    # 11. Outputs. Cc is the blood (venous/arterial) POP concentration, the
    #     quantity the paper plots against age in Figure 3.
    # =====================================================================
    Cc <- ca
    Cfat <- fat / v_fat
    Cliver <- liver / v_liver
    Cmammary <- mammary / v_mam
    BW <- bw
    BH <- bh
    AGEyr <- age
    Eh <- eh
    # Diagnostics for the vignette's exact structural-identity checks.
    Pfat <- p_fat               # Eq 6 partition coefficient at the current age
    PfatAdult <- pfat_r         # Eq 6 partition coefficient at the adult reference
    Prich <- p_rich             # Eq 6 richly-perfused partition coefficient
    Vdist <- vd_r               # Eq 11 adult volume of distribution (L)
    CLint <- clint_c            # Eq 8 age-invariant intrinsic clearance per L liver
    EhRef <- eh_r               # Eq 9 extraction ratio at the HL reference age
    Vliver <- v_liver           # current liver volume (L)
    Qliver <- q_liver           # current hepatic blood flow (L/h)
    IntakeOral <- intake_oral   # background oral intake rate (ug/h)
    RAM <- ram                  # hepatic metabolism rate, Eq 2 (ug/h)
    # Mass balance: cumulative intake minus (body burden + cumulative
    # elimination). The conceptus holds the burden carried off at birth.
    MassBalance <- a_intake -
      (liver + fat + rich + slow + skin + mammary + uterine + conceptus + a_elim)
  })
}
