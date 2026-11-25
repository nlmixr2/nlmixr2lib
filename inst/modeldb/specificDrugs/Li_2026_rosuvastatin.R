Li_2026_rosuvastatin <- function() {
  description <- paste(
    "One-compartment population PK model with parallel lagged first-order",
    "and zero-order absorption for oral microdose rosuvastatin (OATP1B / BCRP",
    "probe) in Chinese healthy volunteers and type 2 diabetes patients (Li",
    "2026). The dose splits into a first-order fraction (66.8%, lag 0.88 h)",
    "and a zero-order fraction infused into the central compartment over 7.3",
    "h. CL/F falls with the stool relative abundance of the gut taxon",
    "Clostridium cluster XIVb and is 25.9% lower in ABCG2 rs2231142",
    "heterozygotes (BCRP decreased function). Type 2 diabetes is not a",
    "covariate. Each administration needs two dose records with the same",
    "amount: one to depot and one to central with rate = -2."
  )
  reference <- paste(
    "Li Y, Yang J, Zhou E, Geng K, Lei Z, Yao X, Feng F, Jin Y, Kong W, Hu Y,",
    "Liang H, Li H, Yang G-P, Hong T, Liu D.",
    "A microdose cocktail and population pharmacokinetic study suggests",
    "potentially reduced CYP3A and P-gp activities in type 2 diabetes.",
    "Front Pharmacol. 2026;17:1924307. doi:10.3389/fphar.2026.1924307.",
    sep = " "
  )
  vignette <- "Li_2026_microdose_cocktail_t2d"
  units <- list(time = "h", dosing = "ug", concentration = "pg/mL")

  covariateData <- list(
    GUT_CLOSTRIDIUM_XIVB_PCT = list(
      description = "Stool relative abundance of the gut taxon Clostridium cluster XIVb",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "16S rRNA V3-V4 sequencing (primers 338F/806R, MGI-2000) of a",
        "pre-dose stool sample; the RDP genus-level taxon Clostridium_XlVb.",
        "Enters as (Clostridium_XlVb / 0.0628)^-0.0847 on CL/F (Li 2026",
        "Eq. 5), so the effect is 1 at 0.0628%. Figure 3e/f spans about",
        "0-0.42%. The negative exponent makes CL/F infinite at an abundance",
        "of exactly 0, so supply a positive value (e.g. the sequencing",
        "detection floor)."
      ),
      source_name = "Clostridium_XlVb"
    ),
    SNP_ABCG2_RS2231142_HET = list(
      description = "ABCG2 (BCRP) rs2231142 heterozygote indicator (1 = GT, BCRP decreased function; 0 = GG, normal function)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Li 2026 classifies BCRP phenotype from rs2231142 per CPIC: GG",
        "normal function (NF), GT decreased function (DecF), TT poor",
        "function. Only GG and GT were observed (Table S5), so the model's",
        "'BCRP DecF' factor 0.741 (Eq. 5) is the GT heterozygote effect;",
        "the model was not fitted to TT homozygotes."
      ),
      source_name = "BCRP DecF"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "rosuvastatin", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "rosuvastatin", units = "ug", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 27L,
    n_studies = 2L,
    age_range = "18-65 years (inclusion criterion)",
    age_mean = "35.6 years (healthy volunteers); 49.2 years (type 2 diabetes)",
    weight_mean = "58.8 kg (healthy volunteers); 67.1 kg (type 2 diabetes)",
    sex_female_pct = 66.7,
    race_ethnicity = c(Chinese = 100),
    disease_state = paste(
      "14 healthy volunteers and 13 type 2 diabetes patients (HbA1c 7.62 +/-",
      "0.77%, disease duration 1.5-8.5 years, normal hepatic and renal",
      "function; 7 on metformin, 2 acarbose, 2 glyburide, 1 dapagliflozin)."
    ),
    dose_range = paste(
      "Single oral dose of rosuvastatin 50 ug in a fasted microdose cocktail",
      "with midazolam 10 ug, dabigatran etexilate 375 ug, pitavastatin 10 ug",
      "and atorvastatin 100 ug."
    ),
    regions = "China (Peking University Third Hospital, Beijing)",
    renal_function = "Mean eGFR 106.6 (healthy volunteers) and 103.6 (type 2 diabetes) mL/min/1.73 m^2",
    notes = paste(
      "Li 2026 Table 1. Healthy-volunteer data come from the earlier Kong",
      "2025 cocktail study and were sampled at 0.5-24 h; type 2 diabetes",
      "patients were sampled at 1, 4 and 24 h. ABCG2 rs2231142 GT in 7/13",
      "healthy volunteers and 5/13 patients genotyped (Table S6).",
      "Microbiome data were available for 9 healthy volunteers and 8",
      "patients (Figure 3). NCT06273215."
    )
  )

  ini({
    lka <- log(0.0651); label("First-order absorption rate constant of the depot fraction (1/h)") # Table S4 RSV Ka 0.0651 (RSE 16.0%)
    ltlag <- log(0.88); label("Lag time of the first-order absorption, ALAG1 (h)") # Table S4 RSV ALAG1 0.88 (RSE 6.60%)
    lcl <- log(66.9); label("Apparent clearance CL/F at Clostridium_XlVb 0.0628% in BCRP normal function (L/h)") # Table S4 RSV CL/F 66.9 (RSE 9.60%); Eq. 5
    lvc <- log(33.1); label("Apparent volume V/F (L)") # Table S4 RSV VC/F 33.1 (RSE 25.4%)
    lfdepot <- log(0.668); label("Fraction of the dose absorbed first-order through the depot, F1 (fraction)") # Table S4 RSV F1 0.668 (RSE 9.00%)
    ld2 <- log(7.3); label("Duration of the zero-order input into central, D2 (h)") # Table S4 RSV D2 7.3 (RSE 3.10%)

    e_gut_clostridium_xivb_pct_cl <- -0.0847; label("Power exponent of Clostridium_XlVb/0.0628 on CL/F (unitless)") # Table S4 'Clostridium_XlVb on CL/F' -0.0847 (RSE 32.8%); Eq. 5
    e_abcg2_het_cl <- 0.741; label("CL/F ratio of BCRP decreased-function (rs2231142 GT) to normal-function subjects (ratio)") # Table S4 'BCRP DecF on CL/F' 0.741 (RSE 12.3%); Eq. 5

    etalcl ~ 0.134 # Table S4 RSV omega^2 (CL/F) 0.134 (RSE 31.3%)
    etalvc ~ 0.292 # Table S4 RSV omega^2 (VC/F) 0.292 (RSE 53.4%)
    etalfdepot ~ 0.615 # Table S4 RSV omega^2 (F1) 0.615 (RSE 43.6%)

    propSd <- sqrt(0.0695); label("Proportional residual error (fraction)") # Table S4 RSV sigma^2 0.0695 (RSE 17.1%); SD = sqrt(variance)
  })
  model({
    # Eq. 5: CL/F = 66.9 x (Clostridium_XlVb / 0.0628)^-0.0847, x 0.741 in
    # BCRP decreased function
    clos_cl <- (GUT_CLOSTRIDIUM_XIVB_PCT / 0.0628)^e_gut_clostridium_xivb_pct_cl
    bcrp_cl <- e_abcg2_het_cl^SNP_ABCG2_RS2231142_HET

    ka <- exp(lka)
    tlag <- exp(ltlag)
    cl <- exp(lcl + etalcl) * clos_cl * bcrp_cl
    vc <- exp(lvc + etalvc)
    d2 <- exp(ld2)
    # Exponential IIV on F1 (Methods: exponential IIV models) exceeds 1 for
    # etalfdepot > log(1 / 0.668); F1 is capped at 1 so the zero-order
    # fraction 1 - F1 cannot go negative.
    fdepot <- min(1, exp(lfdepot + etalfdepot))

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Figure 2c: DOSE x F1 is absorbed first-order from depot after ALAG1 and
    # DOSE x (1 - F1) enters central zero-order over D2, in parallel. Each
    # administration is two dose records of the full amount (depot, and
    # central with rate = -2); f() splits them.
    f(depot) <- fdepot
    alag(depot) <- tlag
    f(central) <- 1 - fdepot
    dur(central) <- d2

    # ug / L = ng/mL; x 1000 -> pg/mL (the paper's concentration unit)
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
