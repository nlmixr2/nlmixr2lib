Aoyama_2017_meloxicam_pkpd <- function() {
  description <- "Typical-value PK/PD model for oral meloxicam and serum thromboxane B2 (TXB2) in healthy Korean men, fitted to digitised mean profiles of CYP2C9*1/*1 and *1/*13 carriers after a single 15 mg dose: two-compartment PK with parallel zero-order (into central) and lagged first-order absorption whose parameters are estimated separately per genotype group, and an indirect-response model with Hill-type inhibition of TXB2 production shared by both groups (Aoyama 2017)"
  reference <- paste(
    "Aoyama T, Hasunuma T, Yamaguchi T, Matsumoto Y. Pharmacokinetics and",
    "Pharmacodynamics of Meloxicam in East Asian Populations: The Role of",
    "Ethnicity on Drug Response. CPT Pharmacometrics Syst Pharmacol.",
    "2017;6(12):823-832. doi:10.1002/psp4.12259.",
    sep = " "
  )
  vignette <- "Aoyama_2017_meloxicam"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")
  # Unit note: doses in mg and volumes in L give central / vc in mg/L; the
  # factor 1000 converts to ng/mL (Table 3 IC50 and Figure 3b axis). TXB2 is
  # a percentage of basal serum TXB2 generation (Eq. 7).

  covariateData <- list(
    CYP2C9_S13_COUNT = list(
      description = "Number of CYP2C9*13 alleles; selects the genotype-group PK parameter set",
      units = "(count, 0/1/2 alleles per subject)",
      type = "continuous",
      reference_category = NULL,
      notes = "The source data contain only two groups, CYP2C9*1/*1 (12 subjects) and *1/*13 (9 subjects), and Table 3 reports a separate PK estimate for each. The model therefore uses CYP2C9_S13_COUNT only as a group selector: 0 = the *1/*1 parameter set, 1 = the *1/*13 parameter set. It is not defined for other genotypes (a *1/*3 subject coded 0 here would receive the *1/*1 values) or for 2 *13 alleles.",
      source_name = "CYP2C9 genotype group (*1/*1 vs *1/*13)"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "meloxicam", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "meloxicam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "meloxicam", units = "mg", specimen = "plasma", verified = TRUE),
    txb2 = list(analyte = "thromboxane B2", units = "% of basal value", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 21,
    n_studies = 1,
    age_range = "mean 23.2 +/- 2.6 years (*1/*1) and 24.4 +/- 2.5 years (*1/*13)",
    weight_range = "mean 71.8 +/- 8.1 kg (*1/*1) and 70.7 +/- 7.3 kg (*1/*13)",
    sex_female_pct = 0,
    race_ethnicity = c(Korean = 21),
    disease_state = "Healthy volunteers",
    dose_range = "Single oral 15 mg meloxicam",
    regions = "South Korea",
    cyp2c9_genotype = c(`*1/*1` = 12, `*1/*13` = 9),
    notes = "Source data are the published mean concentration and percent-inhibition-of-TXB2-generation profiles of Bae et al. 2007 (ref. 14), digitised with UN-SCAN-IT and fitted as group means (no between-subject variability). Sampling pre-dose and 1, 2, 3, 4, 5, 6, 8, 10, 12, 24, 48 and 72 h. Percent TXB2 generation = 100% minus the digitised percent inhibition (Methods)."
  )

  ini({
    # Genotype-group PK estimates (Table 3). Each stratum carries its own
    # suffix (_s1s1 = CYP2C9*1/*1, _s1s13 = CYP2C9*1/*13).
    lcl_s1s1 <- log(0.320); label("Apparent clearance CL/F, CYP2C9*1/*1 group (L/h)") # Table 3 'CL (L/h)' *1/*1 0.320 (RSE 1.35%)
    lcl_s1s13 <- log(0.147); label("Apparent clearance CL/F, CYP2C9*1/*13 group (L/h)") # Table 3 'CL (L/h)' *1/*13 0.147 (RSE 2.12%)
    lvc_s1s1 <- log(5.24); label("Apparent central volume Vc/F, CYP2C9*1/*1 group (L)") # Table 3 'Vc (L)' *1/*1 5.24 (RSE 6.35%)
    lvc_s1s13 <- log(3.38); label("Apparent central volume Vc/F, CYP2C9*1/*13 group (L)") # Table 3 'Vc (L)' *1/*13 3.38 (RSE 22.9%)
    lq_s1s1 <- log(3.41); label("Apparent intercompartmental clearance Q/F, CYP2C9*1/*1 group (L/h)") # Table 3 'Q (L/h)' *1/*1 3.41 (RSE 18.5%)
    lq_s1s13 <- log(3.38); label("Apparent intercompartmental clearance Q/F, CYP2C9*1/*13 group (L/h)") # Table 3 'Q (L/h)' *1/*13 3.38 (RSE 22.4%)
    lvp_s1s1 <- log(5.18); label("Apparent peripheral volume Vp/F, CYP2C9*1/*1 group (L)") # Table 3 'Vp (L)' *1/*1 5.18 (RSE 6.89%)
    lvp_s1s13 <- log(4.55); label("Apparent peripheral volume Vp/F, CYP2C9*1/*13 group (L)") # Table 3 'Vp (L)' *1/*13 4.55 (RSE 20.4%)
    lka_s1s1 <- log(5.01); label("First-order absorption rate constant Ka, CYP2C9*1/*1 group (1/h)") # Table 3 'ka (/h)' *1/*1 5.01 (RSE 9.18%)
    lka_s1s13 <- log(5.27); label("First-order absorption rate constant Ka, CYP2C9*1/*13 group (1/h)") # Table 3 'ka (/h)' *1/*13 5.27 (RSE 8.48%)
    ld1_s1s1 <- log(2.99); label("Zero-order input duration DT (= first-order lag), CYP2C9*1/*1 group (h)") # Table 3 'DT (h)' *1/*1 2.99 (RSE 0.140%); footnote 'tlag = DT'
    ld1_s1s13 <- log(2.96); label("Zero-order input duration DT (= first-order lag), CYP2C9*1/*13 group (h)") # Table 3 'DT (h)' *1/*13 2.96 (RSE 0.242%)
    # Table 3 reports F, the ZERO-order fraction; logitffo = logit(1 - F).
    logitffo_s1s1 <- qlogis(1 - 0.653); label("Logit of the first-order dose fraction 1 - F, CYP2C9*1/*1 group (unitless)") # Table 3 'F' *1/*1 0.653 (RSE 2.25%)
    logitffo_s1s13 <- qlogis(1 - 0.713); label("Logit of the first-order dose fraction 1 - F, CYP2C9*1/*13 group (unitless)") # Table 3 'F' *1/*13 0.713 (RSE 8.79%)

    # PD shared by both groups (Eq. 7; Table 3).
    lkout <- log(0.912); label("First-order loss rate constant of percent serum TXB2 generation kout (1/h)") # Table 3 'kout (/h)' 0.912 (RSE 17.2%)
    lic50 <- log(1390); label("Plasma meloxicam concentration halving TXB2 production IC50 (ng/mL)") # Table 3 'IC50 (ng/mL)' 1390 (RSE 10.2%)
    lhill <- log(0.739); label("Sigmoidicity of the inhibition gamma (unitless)") # Table 3 'gamma' 0.739 (RSE 9.62%)

    # Residual error: the fit used weighted residuals (Figure 3b lower panels)
    # but neither the error model nor its magnitude is reported.
    propSd <- fixed(0); label("Proportional residual SD on meloxicam concentration (fraction); not reported") # Table 3: no residual-error row
    addSd_TXB2 <- fixed(0); label("Additive residual SD on percent TXB2 generation (%); not reported") # Table 3: no residual-error row
  })

  model({
    s13 <- CYP2C9_S13_COUNT
    cl <- exp(lcl_s1s1 * (1 - s13) + lcl_s1s13 * s13)
    vc <- exp(lvc_s1s1 * (1 - s13) + lvc_s1s13 * s13)
    q <- exp(lq_s1s1 * (1 - s13) + lq_s1s13 * s13)
    vp <- exp(lvp_s1s1 * (1 - s13) + lvp_s1s13 * s13)
    ka <- exp(lka_s1s1 * (1 - s13) + lka_s1s13 * s13)
    d1 <- exp(ld1_s1s1 * (1 - s13) + ld1_s1s13 * s13)
    ffo <- expit(logitffo_s1s1 * (1 - s13) + logitffo_s1s13 * s13)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Basal percent TXB2 generation is 100%, so kin = 100 * kout (Methods).
    kout <- exp(lkout)
    kin <- 100 * kout
    ic50 <- exp(lic50)
    hill <- exp(lhill)

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Parallel absorption, as in the population PK model (Figure 2): each
    # dose is two records with the same amt, one on cmt = "central" with
    # rate = -2 (zero-order fraction 1 - ffo over d1) and one on
    # cmt = "depot" (first-order fraction ffo, lagged by tlag = DT).
    f(central) <- 1 - ffo
    dur(central) <- d1
    f(depot) <- ffo
    alag(depot) <- d1

    Cc <- 1000 * central / vc

    inh <- Cc^hill / (ic50^hill + Cc^hill)
    d/dt(txb2) <- kin * (1 - inh) - kout * txb2
    txb2(0) <- 100
    TXB2 <- txb2

    Cc ~ prop(propSd)
    TXB2 ~ add(addSd_TXB2)
  })
}
