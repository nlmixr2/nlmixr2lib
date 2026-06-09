Elkomy_2015_morphine <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model for morphine and its",
    "two glucuronide metabolites M3G and M6G in 20 infants and young",
    "children (3 days - 5.4 years; 3.1 - 18.5 kg) after congenital",
    "heart surgery (Elkomy 2015 AAPS J). Morphine: linear two-",
    "compartment IV disposition with allometric body-weight scaling",
    "(CL and CLD with exponent 0.75 FIXED; VC and VP with exponent",
    "1.0 FIXED) normalised to a reference of 6 kg (study median). Each",
    "metabolite is modeled as a morphine-driven intermediate effect",
    "compartment (rate constant Kint chasing morphine plasma",
    "concentration via dCint/dt = Kint * (Cc - Cint)) feeding an",
    "empirical Emax transduction where metabolite concentration =",
    "Mmax * Cint / (Cint50 + Cint). Estimated glomerular filtration",
    "rate (Schwartz formula) is a covariate: Kint scales linearly with",
    "GFR/70 and Mmax scales as 70/GFR (both exponents FIXED at +/-1",
    "per the paper's covariate analysis). Between-subject random",
    "effects on Kint, Mmax, and Cint50 are SHARED across the M3G and",
    "M6G channels (one eta per group; Table II). Doses are administered",
    "as nmol of morphine equivalents (one nmol of clinical morphine",
    "sulfate yields two nmol of free morphine); concentrations are",
    "output in nmol/L (nM) for morphine, M3G, and M6G."
  )
  reference <- paste(
    "Elkomy MH, Drover DR, Glotzbach KL, Galinkin JL, Frymoyer A,",
    "Su F, Hammer GB. Pharmacokinetics of morphine and its metabolites",
    "in infants and young children after congenital heart surgery.",
    "AAPS J. 2016 Jan;18(1):124-33.",
    "doi:10.1208/s12248-015-9826-5.",
    sep = " "
  )
  vignette <- "Elkomy_2015_morphine"
  units <- list(
    time          = "hour",
    dosing        = "nmol",
    concentration = "nmol/L"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed total body weight. Power-form allometric covariate",
        "on all four morphine PK parameters: CL ~ WT^0.75, CLD (Q) ~",
        "WT^0.75, VC ~ WT^1.0, VP ~ WT^1.0. Exponents FIXED (not",
        "estimated; the paper tested free-exponent estimation and the",
        "OFV reduction was < 3.84 with 1 df, p > 0.05; Results para 2).",
        "Normalised to a reference WT of 6 kg (study median). Age was",
        "NOT identified as a covariate on any morphine parameter once",
        "weight was included (co-linearity r^2 = 0.93 between postnatal",
        "age and WT); the paper attributes the combined size + hepatic-",
        "maturation effect to WT alone (Discussion para 3)."
      ),
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Schwartz-formula estimated glomerular filtration rate",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "eGFR computed via the Schwartz formula (Schwartz 1976) to",
        "account for the effect of age on renal function (paper Methods",
        "and Results, Fig. 3). Covariate on the two metabolite-channel",
        "parameters (SHARED across M3G and M6G): Kint * (GFR/70)",
        "(exponent +1 FIXED) and Mmax * (70/GFR) (exponent -1 FIXED).",
        "The paper fit the GFR exponents first as free parameters,",
        "found them very close to unity, then fixed them with only +2",
        "OFV penalty (Results, final paragraph). Reference value of",
        "70 mL/min/1.73 m^2 represents a typical pediatric normal-",
        "renal-function value. Encoded under the canonical CRCL name in",
        "this register (creatinine-based renal function, BSA-",
        "normalized); the source paper notation is GFR."
      ),
      source_name        = "GFR"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 20,
    n_studies      = 1,
    age_range      = "3 days - 5.4 years",
    age_median     = "1.4 years (SD 1.6)",
    weight_range   = "3.1 - 18.5 kg",
    weight_median  = "7.8 kg (SD 4.3); 6 kg used as the allometric reference",
    sex_female_pct = 50,
    race_ethnicity = NULL,
    disease_state  = paste(
      "Infants and young children admitted to the cardiovascular",
      "intensive care unit (CVICU) after congenital heart surgery.",
      "Cardiac diagnoses (by count): tetralogy of Fallot 12,",
      "atrioventricular septal defect 3, ventricular septal defect 1,",
      "other 4. Exclusion criteria: single-ventricle physiology,",
      "weight < 3 kg, significant renal / hepatic / neurologic",
      "impairment, anticipated mechanical ventilation < 24 h, recent",
      "morphine within 12 h, chronic opioid therapy in the preceding",
      "30 days. Renal function (Schwartz eGFR) range 19 - 158",
      "mL/min/1.73 m^2 (mean 79); hepatic enzymes ALT 25 - 108 U/L,",
      "AST 78 - 659 U/L; gestational age 35 - 40 weeks (mean 38)."
    ),
    dose_range     = paste(
      "IV bolus loading dose 0.15 mg/kg morphine sulfate (range 0.14",
      "- 0.17 mg/kg); nurse-controlled-analgesia (NCA) follow-up IV",
      "bolus doses averaging 0.06 mg/kg (range 0.02 - 0.21 mg/kg).",
      "Mean 18 NCA doses per subject (range 7 - 35) over the 6-h",
      "post-loading sampling period and subsequent random sampling."
    ),
    regions        = "USA (Lucile Packard Children's Hospital, Palo Alto, CA; Children's Hospital Colorado, Aurora, CO)",
    notes          = paste(
      "Demographics from Elkomy 2015 Table I. NONMEM VII with FOCE-",
      "INTER (eta-eps interaction); ADVAN6 TOL5 general ODE solver.",
      "Sequential fitting: morphine 2-compartment model first, then",
      "metabolite (intermediate-effect-compartment + Emax) model",
      "conditional on individual morphine predictions. Concentrations",
      "expressed in molar units (nM) using molecular weights of",
      "758.83, 285.34, and 461.46 g/mol for morphine sulfate,",
      "morphine, and the metabolites, respectively; doses multiplied",
      "by 2 to convert morphine-sulfate moles to free-morphine moles",
      "(one sulfate molecule contains two morphine moieties). 1080",
      "concentration observations across 20 subjects. Bootstrap CIs",
      "from 1000 resampled datasets."
    )
  )

  ini({
    # ================================================================
    # Morphine 2-compartment PK -- allometric body-weight scaling with
    # exponents FIXED at the canonical 0.75 / 1.0 values (paper tested
    # free-exponent estimation and the OFV reduction was < 3.84 with
    # 1 df, p > 0.05; Results para 2). Reference WT = 6 kg (study
    # median).
    # ================================================================
    lcl <- log(8.7);   label("Morphine clearance at WT = 6 kg (L/h)")                         # Table II: theta_CL = 8.7 L/h (%SE 12; bootstrap median 8.6 [7.0, 11.2])
    lvc <- log(18.7);  label("Morphine central volume at WT = 6 kg (L)")                      # Table II: theta_VC = 18.7 L (%SE 21; bootstrap median 18.7 [13.1, 29.0])
    lq  <- log(0.88);  label("Morphine inter-compartmental clearance CLD at WT = 6 kg (L/h)") # Table II: theta_CLD = 0.88 L/h (%SE 15; bootstrap median 0.91 [0.60, 1.4])
    lvp <- log(25.2);  label("Morphine peripheral volume at WT = 6 kg (L)")                   # Table II: theta_VP = 25.2 L (%SE 32; bootstrap median 25.2 [15.0, 61.1])

    e_wt_cl_q  <- fixed(0.75); label("Shared allometric WT exponent on CL and CLD (unitless; FIXED)") # Methods para 4 / Results para 2: exponent 0.75 FIXED
    e_wt_vc_vp <- fixed(1.0);  label("Shared allometric WT exponent on VC and VP (unitless; FIXED)")  # Methods para 4 / Results para 2: exponent 1.0 FIXED

    # ================================================================
    # Metabolite intermediate-effect-compartment + Emax structural
    # model for M3G and M6G. Each metabolite has its own typical-value
    # Kint, Mmax, Cint50 (paper Table II). The GFR covariate is
    # applied as Kint * (GFR / 70) and Mmax * (70 / GFR), with both
    # exponents FIXED to +/-1 (paper: free-exponent estimation was
    # very close to unity; fixed with only +2 OFV penalty; Results,
    # final paragraph).
    # ================================================================
    lkeo_m3g <- log(0.30);  label("Intermediate-compartment rate constant Kint for M3G at GFR = 70 (1/h)") # Table II: theta_Kint,M3G = 0.30 1/h (%SE 14; bootstrap 0.30 [0.21, 0.43])
    lkeo_m6g <- log(0.24);  label("Intermediate-compartment rate constant Kint for M6G at GFR = 70 (1/h)") # Table II: theta_Kint,M6G = 0.24 1/h (%SE 13; bootstrap 0.24 [0.18, 0.32])

    lemax_m3g <- log(1257); label("Maximum M3G concentration Mmax at GFR = 70 (nM)") # Table II: theta_Mmax,M3G = 1257 nM (%SE 14; bootstrap 1254 [950, 1744])
    lemax_m6g <- log(109);  label("Maximum M6G concentration Mmax at GFR = 70 (nM)") # Table II: theta_Mmax,M6G = 109 nM (%SE 15; bootstrap 109 [81.3, 153])

    lec50_m3g <- log(40.4); label("Half-maximal intermediate concentration Cint50 for M3G (nM)") # Table II: Cint,50,M3G = 40.4 nM (%SE 32; bootstrap 40.9 [21.6, 81.4])
    lec50_m6g <- log(21.9); label("Half-maximal intermediate concentration Cint50 for M6G (nM)") # Table II: Cint,50,M6G = 21.9 nM (%SE 27; bootstrap 22.8 [11.9, 43.6])

    # ================================================================
    # IIV -- variances on the log scale, CV = sqrt(value) * 100%.
    # Per Table II, each metabolite-channel parameter (Kint, Mmax,
    # Cint50) has a SINGLE IIV row that is SHARED across the M3G and
    # M6G channels (one eta applied to both metabolites), per the
    # paper's Results-section finding "Between-subject random effects
    # of M3G and M6G were very close. Therefore, inter-subject
    # variability parameters were shared between M3G and M6G."
    # ================================================================
    etalcl ~ 0.20   # Table II IIV on CL = 0.20 (%SE 45) -> CV approx 45%
    etalvc ~ 0.67   # Table II IIV on VC = 0.67 (%SE 82) -> CV approx 82%
    etalq  ~ 0.17   # Table II IIV on CLD = 0.17 (%SE 41) -> CV approx 41%
    etalvp ~ 0.42   # Table II IIV on VP = 0.42 (%SE 65) -> CV approx 65%

    etalkeo  ~ 0.27 # Table II IIV on Kint = 0.27 (%SE 52) -> CV approx 52% (SHARED M3G and M6G)
    etalemax ~ 0.23 # Table II IIV on Mmax = 0.23 (%SE 48) -> CV approx 48% (SHARED M3G and M6G)
    etalec50 ~ 1.26 # Table II IIV on Cint50 = 1.26 (%SE 112) -> CV approx 112% (SHARED M3G and M6G; poorly estimated)

    # ================================================================
    # Residual error -- log-additive in NONMEM (paper Methods:
    # "Concentrations were logarithmically transformed, and an
    # additive residual error model was used") is equivalent to a
    # proportional error model on the linear concentration scale,
    # with proportional-SD ~ sqrt(sigma^2 on log-scale).
    # ================================================================
    propSd     <- sqrt(0.32); label("Proportional residual SD for morphine (fraction)") # Table II: sigma^2(morphine) = 0.32 (%SE 57; bootstrap 0.32 [0.21, 0.50]) -> propSd approx 0.566
    propSd_m3g <- sqrt(0.32); label("Proportional residual SD for M3G (fraction)")      # Table II: sigma^2(M3G) = 0.32 (%SE 57; bootstrap 0.32 [0.19, 0.49]) -> propSd approx 0.566
    propSd_m6g <- sqrt(0.39); label("Proportional residual SD for M6G (fraction)")      # Table II: sigma^2(M6G) = 0.39 (%SE 62; bootstrap 0.38 [0.18, 0.72]) -> propSd approx 0.624
  })

  model({
    # ----------------------------------------------------------------
    # Allometric and renal-function covariate factors. Allometric
    # exponents 0.75 / 1.0 are paper-fixed structural choices (Methods
    # para 4 / Results para 2). GFR linear scaling factors (exponents
    # +/-1) are paper-fixed at the same precision (Results, final
    # paragraph).
    # ----------------------------------------------------------------
    wt_ratio  <- WT / 6
    gfr_ratio <- CRCL / 70

    # ----------------------------------------------------------------
    # Individual morphine PK parameters (Table II structural model).
    # ----------------------------------------------------------------
    cl <- exp(lcl + etalcl) * wt_ratio^e_wt_cl_q
    vc <- exp(lvc + etalvc) * wt_ratio^e_wt_vc_vp
    q  <- exp(lq  + etalq)  * wt_ratio^e_wt_cl_q
    vp <- exp(lvp + etalvp) * wt_ratio^e_wt_vc_vp

    # ----------------------------------------------------------------
    # Individual metabolite-channel parameters. Shared etas
    # (etalkeo / etalemax / etalec50) apply to both M3G and M6G.
    # ----------------------------------------------------------------
    keo_m3g  <- exp(lkeo_m3g  + etalkeo)  * gfr_ratio
    keo_m6g  <- exp(lkeo_m6g  + etalkeo)  * gfr_ratio
    emax_m3g <- exp(lemax_m3g + etalemax) / gfr_ratio
    emax_m6g <- exp(lemax_m6g + etalemax) / gfr_ratio
    ec50_m3g <- exp(lec50_m3g + etalec50)
    ec50_m6g <- exp(lec50_m6g + etalec50)

    # ----------------------------------------------------------------
    # Morphine 2-compartment IV PK. Internal amounts (central,
    # peripheral1) are in nmol of free morphine equivalents; Cc =
    # central / vc gives nmol/L = nM. Doses are administered into
    # `central` as nmol of free morphine.
    # ----------------------------------------------------------------
    d/dt(central)     <- q * peripheral1 / vp - (cl + q) * central / vc
    d/dt(peripheral1) <- q * central / vc - q * peripheral1 / vp

    Cc <- central / vc

    # ----------------------------------------------------------------
    # Metabolite intermediate (effect) compartment chases morphine
    # plasma concentration with rate constant Kint (Eq. 1:
    # dCint/dt = Kint * (Cc - Cint)). Initial condition 0 = no
    # morphine signal before dosing.
    # ----------------------------------------------------------------
    d/dt(effect_m3g) <- keo_m3g * (Cc - effect_m3g)
    d/dt(effect_m6g) <- keo_m6g * (Cc - effect_m6g)

    # ----------------------------------------------------------------
    # Empirical Emax transduction of the intermediate concentration
    # into the observed metabolite plasma concentration (Eq. 2:
    # M = Mmax * Cint / (Cint50 + Cint)). Cc_m3g and Cc_m6g are the
    # observed M3G and M6G plasma concentrations (nM).
    # ----------------------------------------------------------------
    Cc_m3g <- emax_m3g * effect_m3g / (ec50_m3g + effect_m3g)
    Cc_m6g <- emax_m6g * effect_m6g / (ec50_m6g + effect_m6g)

    Cc     ~ prop(propSd)
    Cc_m3g ~ prop(propSd_m3g)
    Cc_m6g ~ prop(propSd_m6g)
  })
}
