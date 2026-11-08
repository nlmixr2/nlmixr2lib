Chakravartula_2024_fentanyl_pk_pediatric <- function() {
  description <- paste(
    "Pediatric (ages 2-12 years, plus infants and adults) IV fentanyl",
    "population PK. Extends the adult Mann 2022 three-compartment + biophase",
    "fentanyl model (Mann_2022_fentanyl_iv, Algera 2021 estimates) to",
    "children by allometrically scaling the plasma-compartment clearance on",
    "total body weight (WT/70)^0.75 and the central volume of distribution",
    "on fat-free mass (FFM/56.1) following Okada et al.; the peripheral",
    "volumes and inter-compartmental clearances inherit the adult body-",
    "weight allometry. Because children have smaller clearance and central",
    "volume than adults, the same body-weight-normalized dose yields a",
    "higher plasma Cmax (Chakravartula 2024 Discussion). Outputs plasma Cc",
    "(ng/mL) and effect-site Ce (ng/mL and pM) for downstream consumption",
    "by the Mann 2022 mu-opioid binding layer in a respiratory-depression",
    "overdose chain. No residual error is reported in the source."
  )
  reference <- paste(
    "Chakravartula S, Thrasher B, Mann J, Chaturbedi A, Han X, Dahan A,",
    "Florian J, Strauss D, Li Z. Physiologically based modeling reveals",
    "different risk of respiratory depression after fentanyl overdose",
    "between adults and children. Clin Transl Sci. 2024;17(4):e13780.",
    "doi:10.1111/cts.13780. PMID 38618722.",
    "Pediatric PK scaling: clearance by (WT/70)^0.75 (Methods) and central",
    "volume by fat-free mass FFM/56.1 (Methods, following Okada et al.",
    "Clin Pharmacokinet 2015). Adult structural fentanyl PK inherited from",
    "Mann J et al. Clin Pharmacol Ther. 2022;112(5):1020-1032",
    "(doi:10.1002/cpt.2696), itself importing the Algera MH et al.",
    "Clin Pharmacol Ther. 2021;109(3):637-645 (doi:10.1002/cpt.2027)",
    "intravenous-fentanyl popPK fit."
  )
  vignette <- "Chakravartula_2024_fentanyl_pediatric_respiratory_depression"
  units <- list(
    time = "min",
    dosing = "mg",
    concentration = "ng/mL"
  )

  compartmentData <- list(
    central = list(analyte = "fentanyl", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "fentanyl", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral2 = list(analyte = "fentanyl", units = "mg", specimen = "plasma", verified = FALSE),
    effect = list(analyte = "fentanyl", units = "mg", specimen = "not applicable", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight; drives the allometric clearance scaling (WT/70)^0.75 and the peripheral-volume / inter-compartmental-clearance scaling inherited from the adult Mann 2022 model.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Chakravartula 2024 Methods: clearance P_i = P * (W_i / W)^0.75",
        "with adult reference W = 70 kg. Age-specific typical WT values",
        "(ages 2-12) are listed in the paper's Table S1; e.g. 34.0 kg at",
        "age 10. Infant (6 months) and child (2.7 years) Fig 2 cohorts",
        "use the Singleton-fitted WT values 7.4 kg and 14.9 kg."
      ),
      source_name = "weight"
    ),
    FFM = list(
      description = "Fat-free mass; drives the central-volume allometric scaling (FFM/56.1) per Okada et al., the pediatric refinement over the adult Mann model (which scales central volume on total body weight).",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Chakravartula 2024 Methods: central volume P_i = P * (FFM_i / FFM)",
        "with adult reference FFM = 56.1 kg (Okada et al.). Typical",
        "pediatric FFM values (ages 2-12) are from Table S1 (translating",
        "CDC-growth-chart median weight and BMI to FFM and averaging",
        "sexes); e.g. 25.20 kg at age 10. Infant / child Fig 2 cohorts",
        "use the Singleton-fitted FFM values 6.5 kg and 9.7 kg."
      ),
      source_name = "Fat Free Mass"
    ),
    Q_TOTAL_LPM = list(
      description = "Total cardiac output Qb + Qt feeding the inherited Mann 2022 shock-state Q_Scale feedback that concentrates fentanyl in the central / biophase compartment as overdose-driven hyperperfusion raises cardiac output above baseline.",
      units = "L/min",
      type = "continuous",
      reference_category = "4.87 (baseline, gives Q_Scale ~ 1, no amplification)",
      notes = paste(
        "Inherited verbatim from Mann_2022_fentanyl_iv. Standalone PK use",
        "leaves Q_TOTAL_LPM at the baseline 4.87 L/min so Q_Scale ~ 1 and",
        "the effective central volume is unchanged. When composed with the",
        "pediatric physiology layer, pass that model's Q_total so the FDA",
        "shock-state feedback engages during overdose hyperperfusion.",
        "Source: FDA delaymymod.c Q_0 = 4.87/60 L/s baseline and the",
        "Q_Scale sigmoid (lines 353-368)."
      ),
      source_name = "(none; inherited from Mann 2022 PK layer)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = NA_integer_,
    age_range = "2-12 years (pediatric scaling); Fig 2 validation also covers infants (6 months) and adults",
    weight_range = "7.4 kg (infant) to 70 kg (adult); Table S1 lists age-specific typical WT 12.5-42.0 kg for ages 2-12",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state = paste(
      "Age-dependent allometric extension of an adult healthy-volunteer",
      "intravenous-fentanyl popPK model to pediatric body sizes. The adult",
      "structural fit pools healthy opioid-naive and chronic opioid-user",
      "volunteers (Algera 2021, n = 30, re-tabulated in Mann 2022)."
    ),
    dose_range = paste(
      "Chakravartula 2024 Fig 2 single IV bolus validation: 31.2 ug/kg",
      "(infant), 30.8 ug/kg (2.7-year child), 20.7 ug/kg (adult), after",
      "Singleton et al. Overdose scenarios: 21.5 ug/kg (child) and",
      "23.2 ug/kg (adult) IV push (the adult dose equals 1.625 mg for a",
      "70 kg adult, the reported median fatal community overdose)."
    ),
    regions = NA_character_,
    notes = paste(
      "The adult-to-pediatric mapping changes ONLY the body-size scaling:",
      "clearance stays on total body weight (WT/70)^0.75 but the central",
      "volume moves from the adult Mann weight scaling (WT/70) to Okada",
      "fat-free-mass scaling (FFM/56.1). Peripheral volumes V2, V3 and",
      "inter-compartmental clearances Q2, Q3 retain the adult Mann body-",
      "weight allometry because Chakravartula 2024 specifies only the",
      "plasma-compartment CL and central-volume rules explicitly.",
      "IIV (omega^2 on the log scale) and the biophase equilibration rate",
      "k1 are inherited unchanged from the adult Mann 2022 / Algera 2021",
      "fit; no pediatric-specific variability was estimated.",
      "Effect-site Ce_pM is the downstream input to the Mann 2022 binding",
      "model (fentanyl = ligand 6)."
    ),
    upstream_model = "Mann J et al. Clin Pharmacol Ther. 2022;112(5):1020-1032 (Mann_2022_fentanyl_iv), importing Algera MH et al. Clin Pharmacol Ther. 2021;109(3):637-645."
  )

  ini({
    # Adult structural PK (Mann 2022 Supplement 1 Table S1, IV fentanyl;
    # rate constants converted from per-second to per-minute by *60).
    lcl  <- log(1.26)
    label("Adult clearance (L/min)")                                        # Mann 2022 Table S1: CL = 0.021 L/s = 1.26 L/min
    lvc  <- log(10.5)
    label("Adult central volume V1 (L)")                                    # Mann 2022 Table S1: V1 = 10.5 L
    lvp  <- log(14.4)
    label("Adult first peripheral volume V2 (L)")                           # Mann 2022 Table S1: V2 = 14.4 L
    lvp2 <- log(166)
    label("Adult second peripheral volume V3 (L)")                          # Mann 2022 Table S1: V3 = 166 L
    lq   <- log(2.04)
    label("Adult inter-compartmental clearance Q2 (L/min)")                 # Mann 2022 Table S1: Q2 = 0.034 L/s = 2.04 L/min
    lq2  <- log(2.28)
    label("Adult inter-compartmental clearance Q3 (L/min)")                 # Mann 2022 Table S1: Q3 = 0.038 L/s = 2.28 L/min
    lk1  <- fixed(log(0.18))
    label("Biophase equilibration rate k1 (1/min,, no IIV reported)")       # Mann 2022 Table S1: k1 = 0.003 1/s = 0.18/min; omega^2 = 0

    # Pediatric allometric scaling (Chakravartula 2024 Methods)
    wt_ref  <- fixed(70)
    label("Reference adult total body weight (kg)")                         # Chakravartula 2024 Methods: adult W = 70 kg
    ffm_ref <- fixed(56.1)
    label("Reference adult fat-free mass (kg)")                             # Chakravartula 2024 Methods: adult FFM = 56.1 kg (Okada et al.)
    e_wt_cl <- fixed(0.75)
    label("Body-weight allometric exponent on CL and Q (unitless)")         # Chakravartula 2024 Methods: CL = CL * (W_i/70)^0.75
    e_wt_vp <- fixed(1)
    label("Body-weight allometric exponent on peripheral volumes (unitless)") # Inherited from Mann 2022 for V2, V3 (paper specifies only central V)
    e_ffm_vc <- fixed(1)
    label("Fat-free-mass allometric exponent on central volume (unitless)")  # Chakravartula 2024 Methods: Vc = Vc * (FFM_i/56.1)

    # IIV inherited from Mann 2022 / Algera 2021 (log-scale omega^2).
    etalcl  ~ 0.087
    # Mann 2022 Table S1: omega^2 CL = 0.087 (30% CV)
    etalvc  ~ 0.43
    # Mann 2022 Table S1: omega^2 V1 = 0.43 (73% CV)
    etalvp  ~ 0.52
    # Mann 2022 Table S1: omega^2 V2 = 0.52 (82% CV)
    etalvp2 ~ 0.060
    # Mann 2022 Table S1: omega^2 V3 = 0.060 (25% CV)
    etalq   ~ 0.39
    # Mann 2022 Table S1: omega^2 Q2 = 0.39 (69% CV)
    etalq2  ~ 0.098
    # Mann 2022 Table S1: omega^2 Q3 = 0.098 (32% CV)
  })

  model({
    # Pediatric allometric scaling factors
    wt_ratio  <- WT / wt_ref
    ffm_ratio <- FFM / ffm_ref

    # Individual structural parameters. Clearance and inter-compartmental
    # clearances scale on body weight (0.75); the central volume scales on
    # fat-free mass (Okada); peripheral volumes retain the adult body-weight
    # allometry.
    cl  <- exp(lcl  + etalcl)  * wt_ratio^e_wt_cl
    vc  <- exp(lvc  + etalvc)  * ffm_ratio^e_ffm_vc
    vp  <- exp(lvp  + etalvp)  * wt_ratio^e_wt_vp
    vp2 <- exp(lvp2 + etalvp2) * wt_ratio^e_wt_vp
    q   <- exp(lq   + etalq)   * wt_ratio^e_wt_cl
    q2  <- exp(lq2  + etalq2)  * wt_ratio^e_wt_cl
    k1  <- exp(lk1)

    # Fentanyl free base 336.4 g/mol (matches FDA delaymymod.c Mmass).
    fentanyl_mw_g_mol <- 336.4

    # Inherited Mann 2022 shock-state PK amplification (FDA delaymymod.c
    # lines 358-368). Standalone use leaves Q_TOTAL_LPM at 4.87 so q_scale
    # collapses to ~1.
    q_scale_raw <- 1 + expit(-((1.6 - Q_TOTAL_LPM / 4.87) / 0.05))
    q_scale <- (q_scale_raw < 1) * 1 + (q_scale_raw > 2) * 2 +
               (q_scale_raw >= 1) * (q_scale_raw <= 2) * q_scale_raw
    vc_eff <- vc / q_scale

    # Three-compartment IV disposition with biophase effect compartment.
    d/dt(central)     <- -(cl + q + q2) / vc * central +
                          q  / vp  * peripheral1 +
                          q2 / vp2 * peripheral2
    d/dt(peripheral1) <-  q  / vc * central - q  / vp  * peripheral1
    d/dt(peripheral2) <-  q2 / vc * central - q2 / vp2 * peripheral2
    d/dt(effect)      <-  k1 * (central / vc_eff - effect)

    # Outputs
    Cc    <- (central / vc) * 1000                                          # plasma fentanyl (ng/mL) = (mg/L) * 1000
    Ce    <- effect * 1000                                                  # effect-site fentanyl (ng/mL)
    Ce_pM <- effect * 1e9 / fentanyl_mw_g_mol                               # effect-site fentanyl (pM); feeds Mann_2022_mu_receptor_binding
  })
}
