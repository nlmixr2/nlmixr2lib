vanErp_2016_everolimus <- function() {
  description <- "Semi-physiological population PK model for oral everolimus in adults with advanced thyroid or breast cancer (van Erp 2016): depot plus four absorption transit compartments (ktr = (n + 1) / MAT) feeding a well-stirred liver compartment (hepatic blood flow 90 L/h, liver volume 1 L, liver plasma flow QH * (1 - hematocrit), hepatic extraction from the unbound intrinsic clearance), and two-compartment plasma disposition. Parameters are plasma-based; whole-blood concentration is plasma * (1 + 4.667 * hematocrit) from a fixed 85:15 erythrocyte:plasma concentration ratio, and the unbound plasma concentration (fu = 0.27) drives the simulation-only S6K1-inhibition Imax model of Tanaka 2008 (Imax 97%, IC50 0.05 ng/mL). Correlated IIV on intrinsic clearance and central volume, IIV on MAT, between-occasion variability on bioavailability, and a proportional residual error on whole blood."
  reference <- paste(
    "van Erp NP, van Herpen CM, de Wit D, Willemsen A, Burger DM,",
    "Huitema ADR, Kapiteijn E, ter Heine R. A Semi-Physiological",
    "Population Model to Quantify the Effect of Hematocrit on Everolimus",
    "Pharmacokinetics and Pharmacodynamics in Cancer Patients.",
    "Clin Pharmacokinet. 2016;55(11):1447-1456.",
    "doi:10.1007/s40262-016-0414-3"
  )
  vignette <- "vanErp_2016_everolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    HCT = list(
      description = "Hematocrit -- packed red-blood-cell volume fraction of whole blood",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying (paired with each whole-blood observation in the source",
        "analysis). Supplied in percent and divided by 100 inside model(),",
        "because van Erp 2016 Eqs. 1 and 3 use the volume fraction",
        "(Cwb = Cp * (1 + 4.667 * Ht); QHP = QH * (1 - Ht)). It enters the",
        "liver plasma flow and therefore the hepatic extraction ratio, and",
        "maps plasma to whole-blood concentration. Cohort median 38%",
        "(range 25-49.7%; van Erp 2016 Table 1); the paper's simulations use",
        "20% and 45%."
      ),
      source_name = "Ht"
    ),
    OCC = list(
      description = "Integer pharmacokinetic sampling-occasion index (1 or 2) for the between-occasion variability on bioavailability",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Each patient was sampled on two occasions (THYRRAD days 1 and 15;",
        "INPRES days 14 and 35; van Erp 2016 Methods 'Bioanalysis'), and the",
        "paper's simulation study uses 'two separate occasions' (Table 3",
        "footnote). OCC = 1 for the first occasion and OCC = 2 for the",
        "second; a single-occasion simulation may use OCC = 1 throughout."
      ),
      source_name = "OCC"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "everolimus", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "everolimus", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "everolimus", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "everolimus", units = "mg", specimen = "administration site", verified = TRUE),
    transit4 = list(analyte = "everolimus", units = "mg", specimen = "administration site", verified = TRUE),
    liver = list(analyte = "everolimus", units = "mg", specimen = "tissue", verified = TRUE),
    central = list(analyte = "everolimus", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "everolimus", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 73L,
    n_studies = 2L,
    age_range = "37-80 years",
    age_median = "62 years",
    weight_range = "45-105 kg",
    weight_median = "73.2 kg",
    sex_female_pct = 69.9,
    race_ethnicity = "Not reported (Dutch centres)",
    disease_state = "Adults with progressive unresectable recurrent or metastatic thyroid cancer (THYRRAD phase II, n = 41) or metastatic breast cancer (INPRES phase IV, n = 32)",
    dose_range = "Everolimus 10 mg orally once daily (dose per occasion 10-15 mg)",
    regions = "Netherlands (Leiden University Medical Center; Radboud university medical center)",
    hematocrit = "Median 38% (range 25-49.7%)",
    sampling = "947 whole-blood EDTA concentrations (915 samples stated in the Results text). THYRRAD: days 1 and 15, limited (0-3 h) or extensive (0-8 h) sampling; INPRES: days 14 and 35, extensive (0-8 h) sampling.",
    notes = "Baseline characteristics from van Erp 2016 Table 1 (22 male, 51 female). Strong CYP3A inhibitors / inducers were prohibited. Whole-blood everolimus measured by UPLC-MS/MS (THYRRAD, 2-160 ug/L) or HPLC-MS/MS (INPRES, 1-150 ug/L)."
  )

  ini({
    # Structural plasma-based PK parameters (van Erp 2016 Table 2, 'Final
    # model parameter estimates'; footnote a: 'All pharmacokinetic parameters
    # are shown as the corresponding plasma pharmacokinetic parameters').
    lmat <- log(0.544); label("Mean absorption time through the 4-transit-compartment chain MAT (h)") # Table 2 'MAT through 4 transit compartments (h)' 0.544 (RSE 6.4%)
    lclint <- log(198); label("Apparent intrinsic hepatic clearance CLint/F, plasma-based (L/h)") # Table 2 'CLint/F (L/h)' 198 (RSE 4.3%)
    lvc <- log(207); label("Apparent central volume V3/F, plasma-based (L)") # Table 2 'V3/F (L)' 207 (RSE 5.0%)
    lq <- log(72.1); label("Apparent inter-compartmental clearance Q/F, plasma-based (L/h)") # Table 2 'Q (L/h)' 72.1 (RSE 3.2%)
    lvp <- log(485); label("Apparent peripheral volume V4/F, plasma-based (L)") # Table 2 'V4/F (L)' 485 (RSE 4.2%)

    # Fixed physiological constants of the well-stirred liver model
    # (Methods 'Structural Pharmacokinetic Model', after Gordi 2005).
    lqh <- fixed(log(90)); label("Hepatic blood flow QH (L/h)") # Methods: 'assuming a liver blood flow (QH) of 90 L/h'
    lvliver <- fixed(log(1)); label("Hepatic volume Vliver (L)") # Methods: 'a hepatic volume (Vliver) of 1 L'
    lfu <- fixed(log(0.27)); label("Unbound fraction of everolimus in plasma fu (fraction)") # Table 2 'Unbound fraction in plasma' 0.27 (fixed); Methods Eq. 4
    rbc_plasma_ratio <- fixed(85 / 15); label("Erythrocyte:plasma everolimus concentration ratio (unitless)") # Table 2 'Red blood cell:plasma concentration ratio' 85:15 (fixed); Eq. 1 coefficient 4.667 = 85/15 - 1
    ntr <- fixed(4); label("Number of absorption transit compartments n (count)") # Results: 'a total of four transit compartments'; Eq. 2 ktr = (n + 1)/MAT

    # S6K1-inhibition model used for the paper's simulation study (Eq. 6,
    # taken from Tanaka 2008 J Clin Oncol, van Erp 2016 refs 16 and 22).
    imax_s6k1 <- fixed(97); label("Maximum tumour S6K1 inhibition Imax (%)") # Methods Eq. 6: 'maximum percentage inhibition of S6K1 activity (97 %)'
    ic50_s6k1 <- fixed(0.05); label("Unbound plasma everolimus concentration giving 50% of Imax, IC50 (ng/mL)") # Methods Eq. 6: 'IC50 ... (0.05 ng/mL)'

    # IIV (Table 2; log-normal, CV% converted with omega^2 = log(1 + CV^2)).
    # Correlation CLint-V3 0.347 -> covariance = 0.347 * sqrt(omega2_CL * omega2_V3).
    etalclint + etalvc ~ c(
      0.14091,
      0.045590, 0.12250
    ) # Table 2 IIV intrinsic clearance 38.9%, IIV V3/F 36.1%, correlation 0.347
    etalmat ~ 0.32527 # Table 2 IIV of MAT 62.0%
    etalfu ~ fixed(0.00089960) # Table 2 IIV in unbound fraction 3% CV

    # Between-occasion variability on relative bioavailability F1 (Table 2
    # 15.2%), two occasions sharing one variance (NONMEM SAME).
    etaiov_fdepot_1 ~ 0.022841 # Table 2 inter-occasion variability on F1 15.2%
    etaiov_fdepot_2 ~ fixed(0.022841) # same variance as occasion 1

    # Residual error on the observed whole-blood concentration.
    propSd <- 0.26; label("Proportional residual error on whole-blood concentration (fraction)") # Table 2 'Residual error (%)' 26.0 (RSE 20.8%)
  })

  model({
    # Hematocrit as a volume fraction (HCT is supplied in percent).
    hct_frac <- HCT / 100

    # Individual parameters
    mat <- exp(lmat + etalmat)
    clint <- exp(lclint + etalclint)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)
    qh <- exp(lqh)
    vliver <- exp(lvliver)
    fu <- exp(lfu + etalfu)

    # Between-occasion variability on bioavailability (typical F1 = 1)
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2
    fdepot <- exp(iov_fdepot)

    # Well-stirred liver (Eqs. 3-5)
    qhp <- qh * (1 - hct_frac)
    eh <- clint * fu / (qhp + clint * fu)
    clh <- eh * qhp

    # Rate constants (Methods, list after Fig. 1)
    ktr <- (ntr + 1) / mat
    k20 <- clh / vliver
    k23 <- qhp * (1 - eh) / vliver
    k32 <- qhp / vc
    k34 <- q / vc
    k43 <- q / vp

    # Depot -> 4 transit compartments -> liver -> central <-> peripheral
    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(transit4) <- ktr * transit3 - ktr * transit4
    d/dt(liver) <- ktr * transit4 + k32 * central - (k20 + k23) * liver
    d/dt(central) <- k23 * liver - k32 * central - k34 * central + k43 * peripheral1
    d/dt(peripheral1) <- k34 * central - k43 * peripheral1

    f(depot) <- fdepot

    # Plasma concentration (mg/L * 1000 = ng/mL)
    Cc <- 1000 * central / vc
    # Unbound plasma concentration
    Cu <- fu * Cc
    # Whole-blood concentration, Eq. 1: Cwb = Cp * (1 + 4.667 * Ht)
    Cblood <- Cc * (1 + (rbc_plasma_ratio - 1) * hct_frac)
    # Tumour S6K1 inhibition (%), Eq. 6
    s6k1_inhibition <- imax_s6k1 * Cu / (ic50_s6k1 + Cu)

    Cblood ~ prop(propSd)
  })
}
