Henninger_2026_dndi6148_mouse <- function() {
  description <- paste(
    "Preclinical (mouse, female BALB/c infected with Leishmania major Friedlin:Luc).",
    "Translational target-site PK/PD model for the benzoxaborole antileishmanial",
    "DNDI-6148 in a murine cutaneous-leishmaniasis model. One-compartment plasma PK",
    "with first-order absorption, saturable (Michaelis-Menten) elimination and a",
    "power dose effect on relative oral bioavailability; skin, liver and spleen are",
    "non-mass-transferring effect compartments holding a tissue CONCENTRATION and",
    "equilibrating with plasma at a fixed rate scaled by a tissue-to-plasma",
    "penetration coefficient. Free skin concentration drives a sigmoidal Emax",
    "parasite-killing rate acting on the skin parasite bioluminescence pool, and the",
    "lesion area grows at a rate proportional to the log10 parasite burden while",
    "healing at a constant first-order rate. See",
    "modellib('Henninger_2026_dndi6148_human') for the allometrically scaled human",
    "translation used for efficacious-dose prediction."
  )
  reference <- paste(
    "Henninger RH, Schouten WM, Arana B, Gillon JY, Mowbray CE, Kratz JM,",
    "Van Bocxlaer K, Dorlo TPC. (2026). Translational Pharmacokinetic-Pharmacodynamic",
    "Modeling and Efficacious Human Dose Prediction of DNDI-6148 for the Treatment of",
    "Cutaneous Leishmaniasis. Clin Transl Sci 19(4):e70535.",
    "doi:10.1111/cts.70535."
  )
  vignette <- "Henninger_2026_dndi6148"

  # Lesion AREA in mm^2 is a paper-mechanistic state. It is deliberately NOT the
  # canonical `lesion` compartment, which is registered as a drug-CONCENTRATION
  # state for site-of-action penetration models (Mehta 2023 TB cavitary lesion);
  # here the state is the size of the cutaneous lesion, a PD endpoint carrying no
  # drug at all. `lesion_size` is left paper-specific rather than claiming a new
  # canonical; see the vignette Errata.
  paper_specific_compartments <- c("lesion_size")

  units <- list(time = "h", dosing = "ug", concentration = "ug/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling of Vmax/F (exponent 0.75) and V/F (exponent 1) about the",
        "cohort median mouse weight of 22 g = 0.022 kg (Henninger 2026 Equations 1-2",
        "and Table 1 footnote b). The paper writes WT in grams; the canonical WT",
        "column is in kg, so the reference is 0.022 kg and the RATIO is unchanged."
      ),
      source_name        = "WT"
    ),
    DOSE_DNDI6148_MGKG = list(
      description        = "Administered DNDI-6148 oral dose level",
      units              = "mg/kg (free base)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power effect on relative oral bioavailability normalised to the lowest",
        "studied dose level, F = (DOSE_DNDI6148_MGKG / 6.25)^e_dose_fdepot",
        "(Table 1 footnote c). Time-fixed per animal: each mouse stays on one dose",
        "level for the whole 10-day bid course. Needed as a data column because",
        "rxode2 cannot read the amt of the dose record it is scaling."
      ),
      source_name        = "dose"
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "DNDI-6148", units = "ug",        specimen = "administration site", verified = TRUE),
    central     = list(analyte = "DNDI-6148", units = "ug",        specimen = "plasma",              verified = TRUE),
    skin        = list(analyte = "DNDI-6148", units = "ug/L",      specimen = "tissue",              verified = TRUE),
    liver       = list(analyte = "DNDI-6148", units = "ug/L",      specimen = "tissue",              verified = TRUE),
    spleen      = list(analyte = "DNDI-6148", units = "ug/L",      specimen = "tissue",              verified = TRUE),
    parasites   = list(analyte = "Leishmania major bioluminescence", units = "photons/s", specimen = "tissue", verified = TRUE),
    lesion_size = list(analyte = "cutaneous lesion area",           units = "mm^2",      specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species        = "mouse (female BALB/c, Leishmania major Friedlin:Luc infected)",
    n_subjects     = 30L,
    n_studies      = 1L,
    weight_median  = "22 g (0.022 kg), the median used to normalise the allometric terms",
    sex_female_pct = 100,
    disease_state  = paste(
      "Cutaneous leishmaniasis. Mice were infected in the rump with 4e7 stationary-phase",
      "L. major Friedlin REH promastigotes (0.2 mL subcutaneous) and randomised once they",
      "had developed a mean lesion diameter of 6.12 +/- 0.90 mm and a bioluminescence",
      "signal of 1.35e8 +/- 7.64e7 photons/s (Supplementary methods S1). Treatment started",
      "14 days post-infection."
    ),
    dose_range     = "6.25, 12.5, 25 and 50 mg/kg DNDI-6148 arginine monohydrate (free-base basis) by oral gavage, twice daily for 10 days; dose intervals 8-12 h",
    notes          = paste(
      "36 mice were infected in six groups of six (Supplementary methods S1); the 30",
      "animals in this model are the 24 DNDI-6148-treated mice plus the 6 vehicle",
      "controls. The sixth group (paromomycin 50 mg/kg qd IP) was an active comparator",
      "and is not part of the model. Data comprised 141 plasma, 43 skin (21 infected,",
      "22 non-infected), 23 liver and 23 spleen PK samples, 207 parasite bioluminescence",
      "and 325 lesion-size measurements (Results 3.1). Fitted in NONMEM 7.5 with FOCE-I",
      "using the PPP&D approach; BLQ handled by the M6 method. Disease state",
      "(infected vs non-infected skin) was tested as a categorical covariate on",
      "R_skin-plasma and was NOT significant, so a single skin compartment serves both."
    )
  )

  ini({
    # =================================================================
    # Plasma PK -- Henninger 2026 Table 1, typical values for a 22 g mouse.
    # Only oral data were available, so CL and V are apparent (relative to F).
    # =================================================================
    lka   <- log(2.7)    ; label("Oral absorption rate constant (ka, 1/h)")                                # Table 1: ka = 2.7 [2.0-3.9]
    lvc   <- log(0.025)  ; label("Apparent central volume of distribution for a 22 g mouse (V/F, L)")      # Table 1: V/F = 0.025 [0.022-0.028]
    lvmax <- log(65)     ; label("Apparent maximal elimination rate for a 22 g mouse (Vmax/F, ug/h)")      # Table 1: Vmax/F = 65 [55-77]
    lkm   <- log(8000)   ; label("Apparent Michaelis-Menten constant (Km, ug/L)")                          # Table 1: Km = 8.0e3 [6.5e3-9.5e3]

    # Allometric exponents. Table 1 footnote b gives them without uncertainty
    # under "power exponent of 0.75 for clearance and 1 for volume of
    # distribution"; Equations 1-2 print them as literal constants. The 0.75
    # exponent applies to the clearance-like parameter, which under
    # Michaelis-Menten elimination is Vmax/F (CL at low concentration = Vmax/Km,
    # and Km is not weight-scaled).
    e_wt_vmax <- fixed(0.75) ; label("Allometric exponent on Vmax/F (unitless)")                           # Eq 1 and Table 1 footnote b
    e_wt_vc   <- fixed(1)    ; label("Allometric exponent on V/F (unitless)")                              # Eq 2 and Table 1 footnote b

    # Dose-dependent relative bioavailability, centred on the 6.25 mg/kg group
    # (Table 1 F_rel = 1.0, fixed). Reproduces the Results 3.2.1 values of
    # 79.1 / 62.6 / 49.5 percent at 12.5 / 25 / 50 mg/kg.
    e_dose_fdepot <- -0.34   ; label("Power exponent for dose level on relative oral bioavailability (unitless; reference 6.25 mg/kg)")  # Table 1: theta_dose = -0.34 [-0.39 - -0.28]

    # =================================================================
    # Tissue distribution -- Equation 3. Separate effect compartments with NO
    # mass transfer from plasma; each state holds a tissue CONCENTRATION.
    # =================================================================
    lke0 <- fixed(log(20)) ; label("Plasma-to-tissue equilibration rate constant (k_plasma-tissue, 1/h)")  # Table 1: 20, Fixed -- tissues sampled at a single time point could not identify it (Methods 2.2.2)
    lkp_skin   <- log(0.56); label("Skin-to-plasma penetration coefficient (R_skin-plasma, unitless)")     # Table 1: 0.56 [0.49-0.68]
    lkp_liver  <- log(1.9) ; label("Liver-to-plasma penetration coefficient (R_liver-plasma, unitless)")   # Table 1: 1.9 [1.6-2.2]
    lkp_spleen <- log(0.34); label("Spleen-to-plasma penetration coefficient (R_spleen-plasma, unitless)") # Table 1: 0.34 [0.28-0.41]

    # Not paper-estimated: an in-house measurement the authors state and assume
    # to hold in skin as well as plasma.
    fu <- fixed(0.066)     ; label("Fraction unbound in plasma and in skin (unitless; in-house measurement assumed equal in both matrices)")  # Methods 2.2.3: fu = 6.6% for BALB/c mice, "a similar protein binding for skin tissue was assumed"

    # =================================================================
    # Drug-induced parasite elimination -- Equations 4-5, Table 2.
    # =================================================================
    lemax <- log(0.049)    ; label("Maximal rate of parasite elimination (Emax, 1/h)")                     # Table 2: 0.049 [0.042-0.058]
    lec50 <- log(165)      ; label("Free skin concentration achieving half of Emax (fEC50, ug/L)")         # Table 2: 165 [125-236]
    lhill <- log(1.6)      ; label("Sigmoidicity factor on free skin concentration (gamma, unitless)")     # Table 2: 1.6 [1.1-2.8]

    # Baseline is reported as 8.01 log10 photons/s; the STATE is the linear
    # bioluminescence in photons/s, so the initial condition is 10^8.01.
    # Corroborated by Supplementary methods S1, which reports a pre-treatment
    # signal of 1.35e8 +/- 7.64e7 photons/s (mean 8.13 log10, CV 56.6%) against
    # the model's 8.01 log10 and 60% CV.
    lrbase_parasites <- log(10^8.01) ; label("Baseline skin parasite bioluminescence (photons/s)")         # Table 2: Parasite_base = 8.0 log10 photons/s [7.9-8.1]; Results 3.2.3 gives 8.01

    # =================================================================
    # Lesion healing -- Equation 6, Table 2.
    # =================================================================
    lrbase_lesion <- log(31)     ; label("Baseline lesion area (mm^2)")                                    # Table 2: Lesion_base = 31 [26-36]
    lslope_lesion <- log(0.0029) ; label("Parasite-induced lesion growth slope, per log10 photons/s (1/h)")# Table 2: slope = 0.0029 [0.0010-0.0050]
    lkheal        <- log(0.027)  ; label("Lesion healing rate constant (k_heal, 1/h)")                     # Table 2: 0.027 [0.010-0.044]; Results 3.2.3 gives 0.0272

    # =================================================================
    # Between-subject variability. Table 1 / Table 2 report CV%; the footnote
    # defines CV = sqrt(exp(OMEGA) - 1), so OMEGA = log(CV^2 + 1).
    # =================================================================
    etalkm              ~ 0.004216  # Table 1: BSV on Km = 6.5% CV [3.7-8.6];   log(0.065^2 + 1)
    etalkp_skin         ~ 0.128309  # Table 1: BSV on R_skin-plasma = 37% CV [23-47];   log(0.37^2 + 1)
    etalkp_liver        ~ 0.086178  # Table 1: BSV on R_liver-plasma = 30% CV [21-42];  log(0.30^2 + 1)
    etalkp_spleen       ~ 0.176977  # Table 1: BSV on R_spleen-plasma = 44% CV [34-67]; log(0.44^2 + 1)
    etalrbase_parasites ~ 0.307485  # Table 2: BSV on Parasite_base = 60% CV [45-94];   log(0.60^2 + 1)
    etalrbase_lesion    ~ 0.110620  # Table 2: BSV on Lesion_base = 35% CV [26-43]; Results 3.2.3 gives 34.2%; log(0.342^2 + 1)

    # =================================================================
    # Residual unexplained variability.
    # =================================================================
    addSd  <- fixed(2.5) ; label("Additive residual SD for plasma concentration (ug/L)")                   # Table 1: add_plasma = 2.5, Fixed
    propSd <- 0.36       ; label("Proportional residual SD for plasma concentration (fraction)")           # Table 1: prop_plasma = 36% CV [31-42]

    addSd_log10BLI <- 0.13 ; label("Additive residual SD for parasite bioluminescence (log10 photons/s)")  # Table 2: log10 add_parasite = 0.13; NOTE the printed 95% CI [0.045-0.075] does not contain the point estimate -- see vignette Errata

    addSd_lesionArea  <- fixed(0.39) ; label("Additive residual SD for lesion area (mm^2)")                # Table 2: add_lesion = 0.39, Fixed
    propSd_lesionArea <- 0.37        ; label("Proportional residual SD for lesion area (fraction)")        # Table 2: prop_lesion = 37% CV [30-43]
  })

  model({
    # ---- 1. Derived covariate terms -------------------------------------
    # Equations 1-2. The paper normalises by WT(g)/22; WT is carried in kg, so
    # the reference is 0.022 kg and the ratio is identical.
    wt_ratio <- WT / 0.022

    # Table 1 footnote c: F = (dose / 6.25 mg/kg)^theta_dose.
    fdepot <- (DOSE_DNDI6148_MGKG / 6.25)^e_dose_fdepot

    # ---- 2. Individual parameters ---------------------------------------
    ka   <- exp(lka)
    vc   <- exp(lvc)   * wt_ratio^e_wt_vc
    vmax <- exp(lvmax) * wt_ratio^e_wt_vmax
    km   <- exp(lkm + etalkm)

    ke0        <- exp(lke0)
    kp_skin    <- exp(lkp_skin   + etalkp_skin)
    kp_liver   <- exp(lkp_liver  + etalkp_liver)
    kp_spleen  <- exp(lkp_spleen + etalkp_spleen)

    emax <- exp(lemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    rbase_parasites <- exp(lrbase_parasites + etalrbase_parasites)
    rbase_lesion    <- exp(lrbase_lesion    + etalrbase_lesion)
    slope_lesion    <- exp(lslope_lesion)
    kheal           <- exp(lkheal)

    # ---- 3. Plasma PK ----------------------------------------------------
    Cc <- central / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - vmax * Cc / (km + Cc)

    f(depot) <- fdepot

    # ---- 4. Tissue distribution (Equation 3) -----------------------------
    # dCtissue/dt = k_plasma-tissue * (R_tissue-plasma * Aplasma/Vplasma - Ctissue).
    # No mass is removed from `central`: these are effect compartments, so the
    # tissue states are concentrations (ug/L) rather than amounts.
    d/dt(skin)   <- ke0 * (kp_skin   * Cc - skin)
    d/dt(liver)  <- ke0 * (kp_liver  * Cc - liver)
    d/dt(spleen) <- ke0 * (kp_spleen * Cc - spleen)

    Cskin   <- skin
    Cliver  <- liver
    Cspleen <- spleen

    # ---- 5. Parasite killing (Equations 4-5) -----------------------------
    # Equation 4 is Emax * (Cskin*fu)^gamma / (EC50^gamma + (Cskin*fu)^gamma);
    # EC50 is the FREE-drug value fEC50 because the driver is the free skin
    # concentration (Results 3.2.3, "directly linked to the free skin target
    # site concentration").
    fCskin <- Cskin * fu
    kkill  <- emax * fCskin^hill / (ec50^hill + fCskin^hill)

    d/dt(parasites) <- -kkill * parasites
    parasites(0)    <-  rbase_parasites

    # Observations were made on the log10 bioluminescence scale (Table 2 RUV
    # unit is log10 photons/s).
    log10BLI <- log10(parasites)

    # Fractional reduction from each animal's own baseline. Under Equation 5
    # this is independent of the baseline value.
    parasiteReduction <- 100 * (1 - parasites / rbase_parasites)

    # ---- 6. Lesion dynamics (Equation 6) ---------------------------------
    # Equation 6 prints dAlesion/dt = slope * Aparasite * Alesion - kheal *
    # Alesion. Aparasite there is the LOG10 bioluminescence, not the linear
    # photons/s that Equation 5 propagates: with the linear state (order 1e8)
    # the growth term would exceed kheal by seven orders of magnitude and the
    # lesion would diverge instantly, whereas slope * 8.01 = 0.0232 /h sits
    # alongside kheal = 0.0272 /h, i.e. the lesion is near steady state at the
    # pre-treatment parasite burden as a model with no other growth term
    # requires. Encoding it this way reproduces the four published end-of-
    # treatment lesion reductions (Results 3.2.5); the fully-log10 reading
    # predicts >99% at every dose level and is falsified. See vignette Errata.
    kgrow <- slope_lesion * log10BLI

    d/dt(lesion_size) <- (kgrow - kheal) * lesion_size
    lesion_size(0)    <-  rbase_lesion

    lesionArea <- lesion_size

    # ---- 7. Error model --------------------------------------------------
    Cc         ~ add(addSd) + prop(propSd)
    log10BLI   ~ add(addSd_log10BLI)
    lesionArea ~ add(addSd_lesionArea) + prop(propSd_lesionArea)
  })
}
