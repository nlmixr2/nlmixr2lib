Verner_2009_persistentOrganicPollutants_infant <- function() {
  description <- paste(
    "PBPK (whole-body infant submodel; original implementation in acslX).",
    "Infant (birth to 12 months) toxicokinetics of persistent organic",
    "pollutants (POPs) through placental transfer and breast-feeding (Verner",
    "et al. 2009). This is the FIVE-compartment infant submodel (liver,",
    "adipose, richly perfused, poorly perfused, brain) that Verner et al.",
    "2009 added to their previously published maternal lifetime PBPK model",
    "(Verner et al. 2008, shipped separately as",
    "Verner_2008_persistentOrganicPollutants_pbpk). Model time is HOURS SINCE",
    "BIRTH, so age (years) = time/8760 drives every infant organ volume and",
    "blood flow through the sex-specific age/body-weight/body-height",
    "regressions of Table 1. Flow-limited mass-balance: for each tissue",
    "dA/dt = Qt*(Ca - Ct/Pt:b), with the liver carrying breast-milk intake",
    "and hepatic metabolism. Tissue:blood partition coefficients are the",
    "ratio of tissue to blood lipid fractions (Eq 8; valid for log Kow > 4).",
    "Exposure has two routes: (1) placental transfer, encoded as the initial",
    "body burden where every infant tissue's lipid-adjusted concentration at",
    "birth equals the maternal blood lipid-adjusted concentration at delivery",
    "(CMAT_LA); (2) breast-milk ingestion over the exclusive breast-feeding",
    "period (BF_MO months), milk intake volume and milk lipid content as",
    "functions of infant age (Eqs 6-7). The compound is selected by its",
    "liver-volume-adjusted intrinsic clearance CLINT (Table 3); default is",
    "PCB-153. Deterministic typical-value model (no IIV, no residual error):",
    "the paper is a forward simulation tool validated against an Inuit",
    "mother-infant cohort through Spearman correlations, not fitted with",
    "random effects."
  )
  reference <- paste(
    "Verner MA, Ayotte P, Muckle G, Charbonneau M, Haddad S. (2009). A",
    "physiologically based pharmacokinetic model for the assessment of infant",
    "exposure to persistent organic pollutants in epidemiologic studies.",
    "Environ Health Perspect 117(3):481-487. doi:10.1289/ehp.0800047. PMCID",
    "PMC2661921. Infant submodel equations and parameters from the main text",
    "(Eqs 1-8, Tables 1-3). The maternal model and its physiology/milk",
    "regressions were published in Verner MA et al. (2008) Environ Health",
    "Perspect 116(7):886-892, doi:10.1289/ehp.10917 (Supplemental Tables",
    "S1-S3); the infant organ-volume/blood-flow regressions originate in",
    "Haddad S et al. (2006) and Price K et al. (2003)."
  )
  vignette <- "Verner_2009_infant_persistent_organic_pollutants"

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
    "brain",
    "a_intake",
    "a_elim",
    "auc_blood"
  )

  compartmentData <- list(
    liver = list(analyte = "POP", units = "ug", specimen = "tissue", verified = TRUE),
    fat = list(analyte = "POP", units = "ug", specimen = "tissue", verified = TRUE),
    rich = list(analyte = "POP", units = "ug", specimen = "tissue", verified = TRUE),
    slow = list(analyte = "POP", units = "ug", specimen = "tissue", verified = TRUE),
    brain = list(analyte = "POP", units = "ug", specimen = "tissue", verified = TRUE),
    a_intake = list(analyte = "POP", units = "ug", specimen = "not applicable", verified = TRUE),
    a_elim = list(analyte = "POP", units = "ug", specimen = "not applicable", verified = TRUE),
    auc_blood = list(analyte = "POP", units = "ug*h/L", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 0L,
    n_studies = 1L,
    age_range = "0-12 months (infant submodel)",
    weight_range = paste(
      "Driven by a 50th-percentile infant growth curve, not a fixed value.",
      "Female (default): 3.4 kg at birth to 9.0 kg at 12 months; male: 3.5 kg",
      "to 9.6 kg (CDC 50th percentile, Kuczmarski et al. 2000). Selected by",
      "SEXF (1 = female, default)."
    ),
    sex_female_pct = 100,
    disease_state = "healthy infants; pre- and postnatal environmental POP exposure",
    dose_range = paste(
      "No dose records. Two exposure routes: placental transfer at birth",
      "(initial body burden set so every infant tissue's lipid-adjusted",
      "concentration equals the maternal delivery blood lipid-adjusted",
      "concentration CMAT_LA, default 100 ug/kg lipids) and breast-milk",
      "ingestion for the first BF_MO months (default 3)."
    ),
    regions = paste(
      "Model framework is generic; validated against Inuit mothers and infants",
      "from Nunavik (Northern Quebec, Canada; Muckle et al. 2001)."
    ),
    notes = paste(
      "NO SUBJECTS WERE FITTED. Forward PBPK simulation tool; predictions were",
      "compared with measured infant plasma, cord blood and breast-milk levels",
      "through Spearman correlations (r = 0.83-0.96 for PCB-153, PCB-180,",
      "PCB-138, HCB and p,p'-DDE). Seven compounds were simulated by changing",
      "CLINT (Table 3): PCB-153 (default), PCB-180, PCB-138, HCB, beta-HCH,",
      "p,p'-DDE and p,p'-DDT. The maternal lifetime model is shipped",
      "separately as Verner_2008_persistentOrganicPollutants_pbpk and supplies",
      "CMAT_LA and the breast-milk POP concentration."
    )
  )

  ini({
    # =====================================================================
    # Every value here is a recoded literature constant: nothing was
    # estimated, so each parameter is fixed(). Sources: Verner et al. 2009
    # main text (Eqs 1-8, Tables 1-3). Partition coefficients are NOT
    # parameters -- they are computed in model() as the ratio of tissue to
    # blood lipid fractions (Eq 8). Organ volumes and blood flows are
    # likewise computed in model() from the Table 1 regressions.
    # =====================================================================

    # --- Compound-specific intrinsic clearance (Table 3) -----------------
    # Liver-volume-adjusted intrinsic clearance (L/hr/kg liver). Switch
    # compound through rxSolve(params = c(clint = ...)):
    #   PCB-153  0.0082601 (half-life 27.5 y)   <- default
    #   PCB-180  0.0229504 (9.9 y)
    #   PCB-138  0.0144328 (16.3 y)
    #   HCB      0.0395820 (6.0 y)
    #   beta-HCH 0.0309631 (7.6 y)
    #   p,p'-DDE 0.0156840 (15 y)
    #   p,p'-DDT 0.0144328 (5 y)
    clint <- fixed(0.0082601) ; label("Liver-volume-adjusted intrinsic clearance (L/hr/kg liver)") # Table 3, PCB-153

    # --- Infant tissue/blood lipid fractions (Table 2) -------------------
    # Age-dependent fractions are interpolated in model() between the ages
    # tabulated below. Fl is the kg-lipid-per-kg-tissue fraction.
    fl_blood  <- fixed(0.005) ; label("Lipid fraction in infant blood (unitless)")          # Table 2 blood 0-1 y = 0.005
    fl_liver0 <- fixed(0.021) ; label("Lipid fraction in infant liver at birth (unitless)") # Table 2 liver age 0 = 0.021
    fl_liver1 <- fixed(0.041) ; label("Lipid fraction in infant liver at 1 year (unitless)")# Table 2 liver age 1 = 0.041
    fl_fat0   <- fixed(0.347) ; label("Lipid fraction in infant adipose at birth (unitless)")    # Table 2 adipose age 0 = 0.347
    fl_fat05  <- fixed(0.472) ; label("Lipid fraction in infant adipose at 0.5 year (unitless)") # Table 2 adipose age 0.5 = 0.472
    fl_fat9   <- fixed(0.550) ; label("Lipid fraction in infant adipose at 9 years (unitless)")  # Table 2 adipose age 9 = 0.550
    fl_rich   <- fixed(0.018) ; label("Lipid fraction in richly perfused tissue (unitless)")     # Table 2 richly perfused 0-1 y = 0.018
    fl_slow   <- fixed(0.021) ; label("Lipid fraction in poorly perfused tissue (unitless)")     # Table 2 poorly perfused 0-1 y = 0.021
    fl_brain0 <- fixed(0.026) ; label("Lipid fraction in infant brain at birth (unitless)")      # Table 2 brain age 0 = 0.026
    fl_brain15 <- fixed(0.061) ; label("Lipid fraction in infant brain at 1.5 years (unitless)") # Table 2 brain age 1.5 = 0.061

    # --- Breast-milk intake volume (Eq 6) and lipid content (Eq 7) -------
    qmilk_slp <- fixed(-0.114) ; label("Daily milk intake slope (L/kg infant body weight per year of age)") # Eq 6: -0.114*Age + 0.158
    qmilk_int <- fixed(0.158)  ; label("Daily milk intake intercept (L/kg infant body weight per day)")     # Eq 6: 0.158 L/kg at birth
    flmilk_slp <- fixed(0.0034) ; label("Milk lipid fraction slope (per ln-year of age)")  # Eq 7: 0.0034*ln(Age) + 0.0414
    flmilk_int <- fixed(0.0414) ; label("Milk lipid fraction intercept (unitless)")        # Eq 7: reaches 0.0414 at 1 year
    flmilk_birth <- fixed(0.0101) ; label("Milk lipid fraction floor at delivery (unitless)") # text: ~0.0101 at delivery

    # --- Exposure and scenario knobs -------------------------------------
    CMAT_LA <- fixed(100) ; label("Maternal delivery blood lipid-adjusted POP concentration (ug/kg lipids)") # placental transfer + milk driver; illustrative input
    BF_MO   <- fixed(3)   ; label("Exclusive breast-feeding duration (months; 0 = not breast-fed)")          # text/Figure 2 scenario: 3-month breast-feeding
    SEXF    <- fixed(1)   ; label("Infant sex (1 = female, 0 = male)")                                       # Table 1 has sex-specific physiology; default female

    # =====================================================================
    # No IIV and no residual error (deterministic forward-simulation tool).
    # =====================================================================
  })

  model({
    # =====================================================================
    # 1. Master clock and infant 50th-percentile growth curve.
    #    Age in years; body weight BW(age) and body height BH(age) are
    #    piecewise-linear CDC 50th-percentile profiles (Kuczmarski 2000),
    #    sex-selected by SEXF.
    # =====================================================================
    age <- time / 8760

    bw_f <- 3.4 +
      (5.8 - 3.4) * min(1, max(0, (age - 0) / (0.25 - 0))) +
      (7.3 - 5.8) * min(1, max(0, (age - 0.25) / (0.5 - 0.25))) +
      (8.5 - 7.3) * min(1, max(0, (age - 0.5) / (0.75 - 0.5))) +
      (9.0 - 8.5) * min(1, max(0, (age - 0.75) / (1 - 0.75)))
    bw_m <- 3.5 +
      (6.4 - 3.5) * min(1, max(0, (age - 0) / (0.25 - 0))) +
      (7.9 - 6.4) * min(1, max(0, (age - 0.25) / (0.5 - 0.25))) +
      (9.2 - 7.9) * min(1, max(0, (age - 0.5) / (0.75 - 0.5))) +
      (9.6 - 9.2) * min(1, max(0, (age - 0.75) / (1 - 0.75)))
    bw <- SEXF * bw_f + (1 - SEXF) * bw_m

    bh_f <- 49.5 +
      (60.0 - 49.5) * min(1, max(0, (age - 0) / (0.25 - 0))) +
      (65.7 - 60.0) * min(1, max(0, (age - 0.25) / (0.5 - 0.25))) +
      (70.0 - 65.7) * min(1, max(0, (age - 0.5) / (0.75 - 0.5))) +
      (74.0 - 70.0) * min(1, max(0, (age - 0.75) / (1 - 0.75)))
    bh_m <- 49.9 +
      (61.4 - 49.9) * min(1, max(0, (age - 0) / (0.25 - 0))) +
      (67.6 - 61.4) * min(1, max(0, (age - 0.25) / (0.5 - 0.25))) +
      (72.0 - 67.6) * min(1, max(0, (age - 0.5) / (0.75 - 0.5))) +
      (75.7 - 72.0) * min(1, max(0, (age - 0.75) / (1 - 0.75)))
    bh <- SEXF * bh_f + (1 - SEXF) * bh_m

    # =====================================================================
    # 2. Infant organ volumes (L) and blood flows (L/hr), Table 1.
    #    Sex-divergent quantities are computed for each sex and blended by
    #    SEXF; M & F equations that agree are written once.
    # =====================================================================
    s_bsa <- bw^0.515 * bh^0.422 * 234.9              # Table 1: BW^0.515 * BH^0.422 * 234.9
    v_liver <- 0.05012 * bw^0.780                     # Table 1 liver (M & F)
    v_tongue <- 0.00119 * bw - 0.0004302              # Table 1 tongue (M & F)
    v_skin <- 0.664 * (s_bsa / 1e4) + 0.07850 * (s_bsa / 1e4)^1.049 # Table 1 skin (M & F)

    v_heart_f <- 1.017e-7 * (bh^0.6862 * bw^0.3561 * 242.7)^1.420
    v_heart_m <- 1.017e-7 * (bh^0.6640 * bw^0.3851 * 242.7)^1.420
    v_heart <- SEXF * v_heart_f + (1 - SEXF) * v_heart_m

    v_skm_f <- 0.09563 * bw + 0.01650 * bh + 0.09102 * age - 0.1642
    v_skm_m <- 0.09561 * bw + 0.01601 * bh + 0.1097 * age
    v_skm <- SEXF * v_skm_f + (1 - SEXF) * v_skm_m

    v_rich_f <- -0.01919 * age + 3.193 * (bw^2 / bh)^0.2657 - 1.374 - v_liver
    v_rich_m <- -0.01068 * age + 2.038 * (bw^2 / bh)^0.4014 - 0.2046 - v_liver
    v_rich <- SEXF * v_rich_f + (1 - SEXF) * v_rich_m

    v_brain_f <- 10 * (age + 0.226) / (6.521 + 7.514 * age)
    v_brain_m <- 10 * (age + 0.213) / (6.030 + 6.895 * age)
    v_brain <- SEXF * v_brain_f + (1 - SEXF) * v_brain_m

    v_slow <- v_tongue + v_heart + v_skm                 # poorly perfused (Table 1 Vpp)
    v_fat <- 0.91 * bw - (v_liver + v_rich + v_slow + v_skin) # adipose (Table 1 Vat)

    qc_f <- 0.2508 * bw^0.7815 * 60
    qc_m <- 0.2519 * bw^0.7609 * 60
    qc <- SEXF * qc_f + (1 - SEXF) * qc_m

    q_liver_f <- v_liver * 60
    q_liver_m <- 0.84 * v_liver * 60
    q_liver <- SEXF * q_liver_f + (1 - SEXF) * q_liver_m

    q_fat_f <- 0.0300 * v_fat * 60
    q_fat_m <- 0.0209 * v_fat * 60
    q_fat <- SEXF * q_fat_f + (1 - SEXF) * q_fat_m

    q_slow_f <- (0.03 * (v_tongue + v_skm) + 0.96 * v_heart) * 60
    q_slow_m <- (0.03 * (v_tongue + v_skm) + 0.73 * v_heart) * 60
    q_slow <- SEXF * q_slow_f + (1 - SEXF) * q_slow_m

    q_brain <- -0.0024 * age^4 + 0.1305 * age^3 - 2.4822 * age^2 +
      18.025 * age + 15.197                              # Table 1 Qbrain (M & F)
    q_rich <- qc - (q_slow + q_fat + q_liver + q_brain)  # balancing flow (Table 1 Qrp)

    # =====================================================================
    # 3. Age-interpolated tissue lipid fractions (Table 2) and partition
    #    coefficients (Eq 8: P = Fl_tissue / Fl_blood).
    # =====================================================================
    fl_liver <- fl_liver0 + (fl_liver1 - fl_liver0) * min(1, max(0, age / 1))
    fl_fat <- fl_fat0 +
      (fl_fat05 - fl_fat0) * min(1, max(0, age / 0.5)) +
      (fl_fat9 - fl_fat05) * max(0, (age - 0.5)) / (9 - 0.5)
    fl_brain <- fl_brain0 + (fl_brain15 - fl_brain0) * min(1, max(0, age / 1.5))

    p_liver <- fl_liver / fl_blood
    p_fat <- fl_fat / fl_blood
    p_rich <- fl_rich / fl_blood
    p_slow <- fl_slow / fl_blood
    p_brain <- fl_brain / fl_blood

    # =====================================================================
    # 4. Hepatic metabolism (Eqs 4-5). Extraction ratio from the intrinsic
    #    clearance (per L liver) at the current liver volume and flow.
    # =====================================================================
    eh <- (clint * v_liver) / (clint * v_liver + q_liver)   # Eq 5

    # =====================================================================
    # 5. Tissue venous concentrations and quasi-steady arterial blood (Eq 3).
    #    Blood storage is negligible on the years timescale of POP kinetics,
    #    so Ca is the flow-weighted mixed venous concentration.
    # =====================================================================
    cv_liver <- (liver / v_liver) / p_liver
    cv_fat <- (fat / v_fat) / p_fat
    cv_rich <- (rich / v_rich) / p_rich
    cv_slow <- (slow / v_slow) / p_slow
    cv_brain <- (brain / v_brain) / p_brain

    ca <- (q_liver * cv_liver + q_fat * cv_fat + q_rich * cv_rich +
      q_slow * cv_slow + q_brain * cv_brain) / qc

    ram <- eh * q_liver * ca                                 # Eq 4: RAM = Eh*Ql*Ca

    # =====================================================================
    # 6. Breast-milk intake (Eqs 1, 6, 7). Milk POP concentration is the
    #    maternal lipid-adjusted blood level times the age-dependent milk
    #    lipid fraction; the infant drinks a body-weight- and age-dependent
    #    milk volume over the exclusive breast-feeding period.
    # =====================================================================
    age_milk <- max(age, 1e-5)                              # avoid log(0) at birth
    fl_milk <- max(flmilk_birth, flmilk_slp * log(age_milk) + flmilk_int) # Eq 7, floored at the delivery value
    c_milk <- CMAT_LA * fl_milk                             # milk POP concentration (ug/L)
    q_milk <- max(0, qmilk_slp * age + qmilk_int) * bw / 24 # Eq 6, L/kg/day -> L/hr
    bf_on <- (age < BF_MO / 12)
    intake <- c_milk * q_milk * bf_on                       # Eq 1: Intake = Cmilk*Qmilk (ug/hr)

    # =====================================================================
    # 7. Mass-balance ODE system (Eq 2). Flow-limited; the liver carries
    #    breast-milk intake and hepatic metabolism.
    # =====================================================================
    d/dt(liver) <- q_liver * (ca - cv_liver) - ram + intake
    d/dt(fat) <- q_fat * (ca - cv_fat)
    d/dt(rich) <- q_rich * (ca - cv_rich)
    d/dt(slow) <- q_slow * (ca - cv_slow)
    d/dt(brain) <- q_brain * (ca - cv_brain)

    d/dt(a_intake) <- intake      # cumulative breast-milk intake
    d/dt(a_elim) <- ram           # cumulative hepatic metabolism
    d/dt(auc_blood) <- ca         # blood-concentration AUC (PKNCA cross-check)

    # =====================================================================
    # 8. Placental transfer -> initial body burden. Every infant tissue's
    #    lipid-adjusted concentration at birth equals the maternal blood
    #    lipid-adjusted concentration CMAT_LA, so A_t(0) = CMAT_LA*Fl_t*V_t.
    # =====================================================================
    liver(0) <- CMAT_LA * fl_liver0 * v_liver
    fat(0) <- CMAT_LA * fl_fat0 * v_fat
    rich(0) <- CMAT_LA * fl_rich * v_rich
    slow(0) <- CMAT_LA * fl_slow * v_slow
    brain(0) <- CMAT_LA * fl_brain0 * v_brain

    # =====================================================================
    # 9. Outputs. Cc is the infant blood POP concentration (ug/L); the paper
    #    plots the lipid-adjusted blood concentration Cblood_LA = Ca/Fl_blood
    #    (ug/kg lipids) in Figures 2-5.
    # =====================================================================
    Cc <- ca
    Cblood_LA <- ca / fl_blood        # lipid-adjusted blood concentration (ug/kg lipids)
    Cfat <- fat / v_fat
    Cbrain <- brain / v_brain
    Cliver <- liver / v_liver
    BW <- bw
    BH <- bh
    AGEyr <- age
    Eh <- eh
    # Diagnostics for the vignette's exact structural-identity checks.
    Pfat <- p_fat                 # Eq 8 adipose:blood partition coefficient
    Pbrain <- p_brain             # Eq 8 brain:blood partition coefficient
    Qc <- qc                      # cardiac output (L/hr)
    Qrich <- q_rich               # richly-perfused flow (balancing, L/hr)
    Vfat <- v_fat                 # adipose volume (L)
    Vliver <- v_liver             # liver volume (L)
    Vbrain <- v_brain             # brain volume (L)
    Intake <- intake              # breast-milk intake rate (ug/hr)
    Cmilk <- c_milk               # breast-milk POP concentration (ug/L)
    RAM <- ram                    # hepatic metabolism rate, Eq 4 (ug/hr)
    # Mass balance: cumulative intake + initial body burden minus
    # (current body burden + cumulative metabolism).
    MassBalance <- a_intake -
      (liver + fat + rich + slow + brain + a_elim)
  })
}
