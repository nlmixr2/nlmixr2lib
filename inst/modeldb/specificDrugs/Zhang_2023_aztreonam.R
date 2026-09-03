Zhang_2023_aztreonam <- function() {
  description <- paste(
    "Two-compartment intravenous population PK model for aztreonam, as implemented by Zhang 2023 to",
    "simulate steady-state free-drug exposure in a 10,000-subject virtual population spanning creatinine",
    "clearance 10-150 mL/min/1.73 m2, for a pharmacodynamic analysis of mutant selection by the",
    "aztreonam/amoxicillin/clavulanate combination against NDM- and serine-beta-lactamase co-producing",
    "Escherichia coli and Klebsiella pneumoniae. Clearance is a power function of BSA-normalized creatinine",
    "clearance referenced to 100 mL/min/1.73 m2; central volume is a power function of body weight",
    "referenced to 70 kg with a near-quadratic exponent (1.99). Three algebraic observables are exposed:",
    "total plasma concentration Cc, unbound plasma concentration Ccfree = Cc * 0.44 (aztreonam is 56%",
    "protein bound), and epithelial-lining-fluid concentration Celf = Cc * 0.40, which the paper treats as",
    "entirely unbound because no mucin binding was assumed. Every parameter is a FIXED prior transcribed",
    "from the upstream renal-impairment analysis; Zhang 2023 estimated nothing and reported no residual",
    "error, so the residual SD is carried structurally at zero (see vignette Errata).",
    sep = " "
  )
  reference <- paste(
    "Zhang J, Wu M, Diao S, Zhu S, Song C, Yue J, Martins FS, Zhu P, Lv Z, Zhu Y, Yu M, Sy SKB.",
    "Pharmacokinetic/pharmacodynamic evaluation of aztreonam/amoxicillin/clavulanate combination against",
    "New Delhi metallo-beta-lactamase and serine-beta-lactamase co-producing Escherichia coli and",
    "Klebsiella pneumoniae. Pharmaceutics. 2023;15(1):251. doi:10.3390/pharmaceutics15010251",
    "(Section 2.4; Supplementary Materials 'Detailed description of pharmacokinetic models' and the",
    "sample RxODE script).",
    "Structural model and parameter estimates originate from Xu H, Zhou W, Zhou D, Li J, Al-Huniti N.",
    "J Clin Pharmacol. 2017;57(3):336-344, the aztreonam population PK analysis in patients with normal",
    "and impaired renal function; the dosing regimens simulated are those of the REJUVENATE study",
    "(Cornely OA et al. J Antimicrob Chemother. 2020;75(3):618-627).",
    sep = " "
  )
  vignette <- "Zhang_2023_aztreonam_amoxicillin_clavulanate"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central     = list(analyte = "aztreonam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "aztreonam", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance, normalized to 1.73 m2 body surface area",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power effect on clearance, (CRCL/100)^0.43, printed identically in Section 2.4 of the main text",
        "and in the Supplementary Materials equation and RxODE script. Zhang 2023 sampled CRCL from a",
        "uniform distribution over 10-150 mL/min normalized to 1.73 m2 (Section 2.4); the lower bound of",
        "10 excludes hemodialysis patients. The simulation script draws a separate uniform per renal",
        "function category: 51-150, 31-50 and 10-30 mL/min for the three dosing arms of Table 1.",
        "Source variable name CrCL.",
        sep = " "
      ),
      source_name        = "CrCL"
    ),
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power effect on central volume, (WT/70)^1.99. The exponent is near-quadratic rather than the",
        "conventional allometric 1, and is transcribed as printed; both the main text (Section 2.4) and",
        "the Supplementary Materials equation and RxODE script give 1.99. Weight is not a covariate on",
        "clearance, on peripheral volume or on intercompartmental clearance in this model. Zhang 2023 did",
        "not sample weight directly: the Supplementary Materials derive it from a sex-specific height",
        "distribution (see covariatesDataExcluded for HT and SEXF).",
        sep = " "
      ),
      source_name        = "WT"
    )
  )

  # Inputs to the paper's virtual-population weight simulation, not covariates on
  # any model parameter. Documented here so the cohort in the validation vignette
  # can be reconstructed without re-reading the supplement.
  covariatesDataExcluded <- list(
    HT = list(
      description        = "Height",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Not a covariate on any model parameter. Supplementary Materials: heights were drawn as",
        "male N(176.3, 0.17*sqrt(4482)) cm and female N(162.2, 0.16*sqrt(4857)) cm, then converted to",
        "weight by WT = exp(3.28 + 1.92*log(HT/100))*exp(eta), eta ~ N(0, 0.14) for males and",
        "WT = exp(3.49 + 1.45*log(HT/100))*exp(eta), eta ~ N(0, 0.17) for females. The /100 conversion",
        "from cm to m appears only in the RxODE script, not in the typeset equation.",
        sep = " "
      ),
      source_name        = "HT"
    ),
    SEXF = list(
      description        = "Sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Not a covariate on any model parameter. Selects which of the two height and weight equations",
        "above applies. Section 2.4 states the virtual population consisted of males and females in",
        "equal proportion; the script draws SEX as round(runif(n, 0, 1)).",
        sep = " "
      ),
      source_name        = "SEX"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 10000,
    n_studies      = 1,
    sex_female_pct = 50,
    renal_function = "creatinine clearance 10-150 mL/min/1.73 m2, uniformly distributed",
    disease_state  = paste(
      "simulated adults; the source population PK analysis was conducted in patients with normal and",
      "impaired renal function, and the dosing regimens simulated are those studied in complicated",
      "intra-abdominal infection",
      sep = " "
    ),
    dose_range     = paste(
      "3 h intravenous infusions by renal function category (Table 1): 2 g loading dose then 1.5 g q6h",
      "for CRCL > 50-150 mL/min; 2 g loading dose then 750 mg q6h for CRCL > 30-50 mL/min; 2 g loading",
      "dose then 500 mg q8h for CRCL 10-30 mL/min. The high-dose regimens do not exceed 6 g daily.",
      sep = " "
    ),
    notes          = paste(
      "The 10,000 subjects are a VIRTUAL population constructed by Zhang 2023 (Section 2.4 and the",
      "Supplementary Materials script), not an observed cohort: sex 50/50, height sampled from US",
      "anthropometric reference data, weight derived from height, and creatinine clearance sampled",
      "uniformly within each renal function category. No demographic table is reported because no",
      "patients were enrolled for the PK component. The parameter estimates come from the upstream",
      "renal-impairment analysis cited in the reference field.",
      sep = " "
    )
  )

  ini({
    # Every parameter is FIXED: Zhang 2023 re-used the upstream population PK
    # estimates as simulation priors and estimated nothing.

    # Section 2.4 gives 4.73 L/h but the Supplementary Materials give 4.93 L/h in
    # the typeset equation AND in both RxODE code blocks, and 4.93 is the value
    # that reproduces Table 4. See the vignette Errata.
    lcl <- fixed(log(4.93)); label("Clearance at CRCL 100 mL/min/1.73 m2 (L/h)")
    # Suppl. Materials, aztreonam equation and script: CL = 4.93*(CrCL/100)^0.43

    lvc <- fixed(log(7.43)); label("Central volume at 70 kg (L)")
    # Section 2.4 and Suppl. Materials: VC = 7.43*(WT/70)^1.99

    lvp <- fixed(log(6.44)); label("Peripheral volume (L)")
    # Section 2.4: VP is 6.44 L (27.7% CV); Suppl. Materials: VP = 6.44*exp(eta)

    lq  <- fixed(log(9.26)); label("Intercompartmental clearance (L/h)")
    # Section 2.4 and Suppl. Materials: Q = 9.26 L/h, no variability reported

    # Covariate exponents
    e_crcl_cl <- fixed(0.43); label("Exponent of CRCL/100 on clearance (unitless)")
    # Section 2.4 and Suppl. Materials: (CrCL/100)^0.43

    e_wt_vc   <- fixed(1.99); label("Exponent of WT/70 on central volume (unitless)")
    # Section 2.4 and Suppl. Materials: (WT/70)^1.99

    # Protein binding and lung penetration, used to build the free-drug and
    # epithelial-lining-fluid observables the pharmacodynamic analysis runs on.
    fu   <- fixed(0.44); label("Fraction of aztreonam unbound in plasma (unitless)")
    # Section 2.4: plasma protein binding of aztreonam is 56%, so fu = 1 - 0.56.
    # Suppl. Materials script: res$free <- res$value*(1-0.56)

    relf <- fixed(0.40); label("Epithelial lining fluid to TOTAL plasma concentration ratio (unitless)")
    # Section 3.4: a 40% penetration into the ELF was assumed for aztreonam.
    # Section 2.4 says 55%; 40% is the value that reproduces Table 5, and the
    # ratio is applied to TOTAL (not unbound) plasma concentration. See the
    # vignette Errata, which shows the arithmetic for all four candidate readings.

    # IIV. The paper reports these as CV%; the Supplementary Materials script
    # passes the same numbers to rnorm(sd = ...) on the log scale, so they are
    # log-scale SDs and the omega entries below are their squares.
    # Section 2.4: CL 24.1% CV; script sd = 0.241, so 0.241^2 = 0.058081
    etalcl ~ fixed(0.058081)
    # Section 2.4: VC 50.9% CV; script sd = 0.509, so 0.509^2 = 0.259081
    etalvc ~ fixed(0.259081)
    # Section 2.4: VP 27.7% CV; script theta block sd = 0.277, so 0.277^2 = 0.076729
    etalvp ~ fixed(0.076729)
    # No variability is reported on Q anywhere in the paper or supplement, and the
    # script draws no eta for it; carried structurally at zero.
    etalq  ~ fixed(0)

    # Residual error is not reported: Zhang 2023 simulated concentration-time
    # profiles from between-subject variability only and never fitted data.
    propSd <- fixed(0); label("Proportional residual SD (fraction; not reported in the source)")
  })

  model({
    # Individual parameters
    cl <- exp(lcl + etalcl) * (CRCL / 100)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    vp <- exp(lvp + etalvp)
    q  <- exp(lq  + etalq)

    # Two-compartment intravenous system, written as in the Supplementary
    # Materials ODEs (A1 = central, A2 = peripheral1)
    d/dt(central)     <- -central * (cl / vc + q / vc) + peripheral1 * (q / vp)
    d/dt(peripheral1) <-  central * (q / vc)            - peripheral1 * (q / vp)

    # Observables. Cc is total plasma concentration (A1 / VC). Ccfree is the
    # unbound plasma concentration that drives the plasma fT>MIC / fT>MPC
    # analysis. Celf is the epithelial-lining-fluid concentration; it scales
    # TOTAL plasma concentration and is itself treated as entirely unbound,
    # because Section 2.4 assumed no mucin binding in the ELF.
    Cc     <- central / vc
    Ccfree <- Cc * fu
    Celf   <- Cc * relf

    Cc ~ prop(propSd)
  })
}
