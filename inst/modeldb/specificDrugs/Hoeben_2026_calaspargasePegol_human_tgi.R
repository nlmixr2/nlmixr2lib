Hoeben_2026_calaspargasePegol_human_tgi <- function() {
  description <- paste(
    "Translational PK-PD tumor-growth-inhibition model of calaspargase pegol in",
    "adults with hepatocellular carcinoma: the population PK model of plasma",
    "asparaginase activity (PAA) fitted in pediatric and adult patients with",
    "acute lymphoblastic leukemia (two compartments with saturable transfer to",
    "the peripheral compartment and parallel linear and saturable elimination;",
    "body surface area, sex and age covariates) drives the Simeoni 2004",
    "tumor-growth-inhibition model whose parameters were estimated in JHH-5",
    "hepatocellular-carcinoma xenograft mice."
  )
  reference <- paste(
    "Hoeben E, Madelain V, Acharya C, Barbier Saint Hilaire P, Tyagi E, Okour M,",
    "Fouliard S. PKPD-Based Translational Modeling of Calaspargase Pegol",
    "Preclinical Activity in Hepatocellular Carcinoma. Eur J Drug Metab",
    "Pharmacokinet. 2026. doi:10.1007/s13318-026-01010-4.",
    "Earlier version of the PK model: FDA Center for Drug Evaluation and",
    "Research, BLA 761102 (Asparlas) Multidiscipline Review, 2018.",
    "Structural TGI formulation from Simeoni M, et al. Cancer Res.",
    "2004;64(3):1094-1101. doi:10.1158/0008-5472.CAN-03-2524.",
    sep = " "
  )
  vignette <- "Hoeben_2026_calaspargasePegol"

  units <- list(time = "h", dosing = "IU", concentration = "IU/L")

  covariateData <- list(
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effects on Vc, CL, VmaxTR and Vmax, normalised to 1.47 m^2, the",
        "median BSA of the 144-subject PK dataset (Supplementary Table S3B).",
        "Neither the paper nor the FDA review states the normalising value; see",
        "the vignette's Assumptions and deviations section. The paper simulated",
        "HCC patients at 2.05 m^2 (men) and 1.81 m^2 (women)."
      ),
      source_name = "BSA"
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "male (SEXF = 0)",
      notes = paste(
        "Table 1 'Female on Vc' = -0.0928, encoded as the proportional shift",
        "Vc * (1 + e_sexf_vc * SEXF). The paper's figures code SEX = 1 male,",
        "SEX = 2 female; SEXF = SEX - 1."
      ),
      source_name = "SEX"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Uncentred exponential effect on Vmax, exp(e_age_vmax * AGE). Neither the",
        "functional form nor a centring value is printed; this form reproduces the",
        "paper's simulated trough PAA (Table 2, Supplementary Table S4) where",
        "centring on the 12.5-year dataset median does not. See the vignette's",
        "Assumptions and deviations section."
      ),
      source_name = "AGE"
    )
  )

  compartmentData <- list(
    central = list(
      analyte = "calaspargase pegol (asparaginase activity)",
      units = "IU",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "calaspargase pegol (asparaginase activity)",
      units = "IU",
      specimen = "tissue",
      verified = TRUE
    ),
    cycling_cells = list(analyte = "proliferating tumor cells", units = "mm^3", specimen = "tumor", verified = TRUE),
    damaged_cells1 = list(
      analyte = "non-replicating dying tumor cells (transit state 1)",
      units = "mm^3",
      specimen = "tumor",
      verified = TRUE
    ),
    damaged_cells2 = list(
      analyte = "non-replicating dying tumor cells (transit state 2)",
      units = "mm^3",
      specimen = "tumor",
      verified = TRUE
    ),
    damaged_cells3 = list(
      analyte = "non-replicating dying tumor cells (transit state 3)",
      units = "mm^3",
      specimen = "tumor",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 144L,
    n_studies = 3L,
    age_range = "1-73 years",
    age_median = "12.5 years",
    weight_range = "8.45-152 kg",
    weight_median = "48.5 kg",
    bsa_range = "0.42-2.80 m^2",
    bsa_median = "1.47 m^2",
    sex_female_pct = 51,
    race_ethnicity = c(White = 81, Asian = 6.2, Black = 4.2, `Pacific Islander` = 1.4, Other = 1.4, Unknown = 5.6),
    disease_state = paste(
      "PK fitted in pediatric and adult patients with acute lymphoblastic",
      "leukemia; applied, by assumption, to adults with hepatocellular carcinoma"
    ),
    dose_range = "PK data: 1000-2500 IU/m^2 IV; HCC simulations: 500-2500 IU/m^2 as a 2-h IV infusion every 3 weeks",
    regions = "North America (AALL07P4, DFCI-11-001) and the CL2-95015-001 adult study",
    notes = paste(
      "PK dataset: 1680 PAA observations in 144 subjects from AALL07P4 (N = 111),",
      "DFCI-11-001 (N = 13) and CL2-95015-001 (N = 20) (Supplementary Tables",
      "S3A-S3C). The tumor-growth-inhibition parameters are those of the JHH-5",
      "xenograft mice (Hoeben_2026_calaspargasePegol_mouse_tgi), carried to humans",
      "under the paper's assumption that tumor growth, potency and the PAA-effect",
      "relationship translate unchanged; no human tumor data were fitted."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Population PK of PAA in pediatric and adult ALL (Table 1 'Human PK
    # Model'). Structure from Methods 'PK Model in Pediatric and Adult
    # Patients with ALL' and Figure 5.
    # ---------------------------------------------------------------------
    lvc <- log(1.78); label("Central volume of distribution Vc at 1.47 m^2, male (L)") # Table 1 Human PK: Vc = 1.78 L
    lcl <- log(0.00245); label("Linear clearance CL at 1.47 m^2 (L/h)") # Table 1 Human PK: CL = 0.00245 L/h
    lvmax_tr <- log(21.4); label("Maximum rate of saturable transfer to the peripheral compartment VmaxTR at 1.47 m^2 (IU/h)") # Table 1 Human PK: VmaxTR = 21.4 U/h
    lkm_tr <- log(58.8); label("PAA at half-maximum peripheral transfer KmTR (IU/L)") # Table 1 Human PK: KmTR = 58.8 U/L
    lk21 <- log(0.0395); label("Peripheral-to-central rate constant k21 (1/h)") # Table 1 Human PK: k21 = 0.0395 1/h
    lvmax <- log(0.583); label("Maximum saturable elimination capacity Vmax at 1.47 m^2, before the age effect (IU/h)") # Table 1 Human PK: Vmax = 0.583 U/h
    lkm <- log(0.00100); label("PAA at half-maximum saturable elimination Km (IU/L)") # Table 1 Human PK: Km = 0.00100 U/L

    e_bsa_vc <- 1.23; label("Power exponent of BSA on Vc (unitless)") # Table 1 Human PK: BSA on Vc = 1.23
    e_sexf_vc <- -0.0928; label("Proportional shift in Vc for women (fraction)") # Table 1 Human PK: Female on Vc = -0.0928
    e_bsa_cl <- 1.19; label("Power exponent of BSA on CL (unitless)") # Table 1 Human PK: BSA on CL = 1.19
    e_bsa_vmax_tr <- 0.726; label("Power exponent of BSA on VmaxTR (unitless)") # Table 1 Human PK: BSA on VmaxTR = 0.726
    e_bsa_vmax <- 1.24; label("Power exponent of BSA on Vmax (unitless)") # Table 1 Human PK: BSA on Vmax = 1.24
    e_age_vmax <- -0.0166; label("Exponential effect of age on Vmax (1/year)") # Table 1 Human PK: AGE on Vmax = -0.0166

    # ---------------------------------------------------------------------
    # Tumor-growth-inhibition parameters carried from the JHH-5 xenograft fit
    # (Table 1 'Mice PKPD TGI Model'), held fixed in the human model.
    # ---------------------------------------------------------------------
    ltumorExpGrowth <- fixed(log(0.00708)); label("Exponential tumor growth rate kg (1/h)") # Table 1 Mice PKPD TGI: kg = 0.00708 1/h
    ltumorLinGrowth <- fixed(log(70.0)); label("Linear tumor growth rate kgl (mm^3/h)") # Table 1 Mice PKPD TGI: kgl = 70.0 mm3/h
    ldrugSlope <- fixed(log(3.22e-5)); label("Slope of the linear PAA kill effect slp (L/IU/h)") # Table 1 Mice PKPD TGI: slp = 3.22 x 10^-5 L/IU/h
    ldamageTransit <- fixed(log(1 / 41.1)); label("Transit rate through the dying-cell stages 1/tau (1/h)") # Table 1 Mice PKPD TGI: tau = 41.1 h
    lrbase_tumor <- fixed(log(152)); label("Tumor volume at time zero W0 (mm^3)") # Table 1 Mice PKPD TGI: W0 = 152 mm3
    psi <- fixed(20); label("Switch sharpness between exponential and linear growth psi (unitless)") # Table 1 Mice PKPD TGI: psi = 20 FIX

    # IIV. Table 1 reports IIV as SD of the log-normal random effect, so the
    # variance is SD^2. CL-k21 covariance = 0.582 * 0.441 * 0.519.
    etalvc ~ 0.0225 # Table 1 Human PK: IIV SD Vc = 0.150
    etalcl + etalk21 ~ c(
      0.194481,
      0.1332076, 0.269361
    ) # Table 1 Human PK: IIV SD CL = 0.441, k21 = 0.519, Corr. CL:k21 = 0.582
    etaltumorExpGrowth ~ fixed(0.050176) # Table 1 Mice PKPD TGI: IIV SD kg = 0.224
    etaldrugSlope ~ fixed(0.128881) # Table 1 Mice PKPD TGI: IIV SD slp = 0.359
    etalrbase_tumor ~ fixed(0.012996) # Table 1 Mice PKPD TGI: IIV SD W0 = 0.114

    # Residual error. The PK model's proportional error (additive on the log
    # scale) is named in Methods but its magnitude is not reported, so it is
    # carried as fixed(0). The tumor-volume error is the xenograft fit's.
    propSd <- fixed(0); label("Proportional residual error on PAA (fraction; ZERO -- not reported in the source)")
    addSd_tumor_vol <- fixed(0.827); label("Additive residual error on tumor volume (mm^3)") # Table 1 Mice PKPD TGI: RUV1 SD additive component = 0.827
    propSd_tumor_vol <- fixed(0.437); label("Proportional residual error on tumor volume (fraction)") # Table 1 Mice PKPD TGI: RUV2 SD proportional component = 0.437
  })

  model({
    # Body-size, sex and age effects. BSA is normalised to 1.47 m^2, the PK
    # dataset median (Supplementary Table S3B); the reference sex is male; the
    # age effect on Vmax is uncentred.
    bsa_ratio <- BSA / 1.47

    vc <- exp(lvc + etalvc) * bsa_ratio^e_bsa_vc * (1 + e_sexf_vc * SEXF)
    cl <- exp(lcl + etalcl) * bsa_ratio^e_bsa_cl
    vmax_tr <- exp(lvmax_tr) * bsa_ratio^e_bsa_vmax_tr
    km_tr <- exp(lkm_tr)
    k21 <- exp(lk21 + etalk21)
    vmax <- exp(lvmax) * bsa_ratio^e_bsa_vmax * exp(e_age_vmax * AGE)
    km <- exp(lkm)

    tumorExpGrowth <- exp(ltumorExpGrowth + etaltumorExpGrowth)
    tumorLinGrowth <- exp(ltumorLinGrowth)
    drugSlope <- exp(ldrugSlope + etaldrugSlope)
    damageTransit <- exp(ldamageTransit)
    rbase_tumor <- exp(lrbase_tumor + etalrbase_tumor)

    # PAA (IU/L); amounts in IU (dose = IU/m^2 x BSA).
    Cc <- central / vc

    # Figure 5: linear (CL) and saturable (Vmax, Km) elimination from the
    # central compartment, saturable transfer to the peripheral compartment
    # (VmaxTR, KmTR) and first-order return (k21).
    transfer <- vmax_tr * Cc / (km_tr + Cc)
    d/dt(central) <- -cl * Cc - vmax * Cc / (km + Cc) - transfer + k21 * peripheral1
    d/dt(peripheral1) <- transfer - k21 * peripheral1

    # Hoeben 2026 Equations 1-2 with the xenograft parameters (Figure 5).
    drugEffectCyclingCells <- drugSlope * Cc
    d/dt(cycling_cells) <- tumorExpGrowth * cycling_cells /
      (1 + (tumorExpGrowth / tumorLinGrowth * cycling_cells)^psi)^(1 / psi) -
      drugEffectCyclingCells * cycling_cells
    d/dt(damaged_cells1) <- drugEffectCyclingCells * cycling_cells - damageTransit * damaged_cells1
    d/dt(damaged_cells2) <- damageTransit * (damaged_cells1 - damaged_cells2)
    d/dt(damaged_cells3) <- damageTransit * (damaged_cells2 - damaged_cells3)
    cycling_cells(0) <- rbase_tumor

    tumor_vol <- cycling_cells + damaged_cells1 + damaged_cells2 + damaged_cells3

    Cc ~ prop(propSd)
    tumor_vol ~ add(addSd_tumor_vol) + prop(propSd_tumor_vol)
  })
}
