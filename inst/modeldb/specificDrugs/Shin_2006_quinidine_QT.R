Shin_2006_quinidine_QT <- function() {
  description <- paste(
    "Population pharmacodynamic Emax model for quinidine-induced QTc",
    "prolongation in 24 healthy Korean (12 M / 12 F) and 13 healthy",
    "Caucasian (7 M / 6 F) adults following a single 20 min IV infusion of",
    "quinidine gluconate 4 mg/kg (base). The Emax form is",
    "QTc(t) = E0 + DeltaEmax * Cc / (EC50 + Cc) with E0 modulated by sex",
    "(additive +34 ms in females; reference category = male) and DeltaEmax",
    "modulated by ethnicity (multiplicative x1.26 in Caucasians;",
    "reference category = Korean) plus an additive +106 ms interaction in",
    "Caucasian females only. EC50 = 3.13 uM (= 1.0155 mg/L using quinidine",
    "MW 324.42 g/mol). Source publication does not fit a popPK model;",
    "the PK driver in this file is a typical-value 1-compartment IV",
    "approximation with CL = 0.3 L/h/kg and Vc = Vss = 2.5 L/kg derived",
    "from the pooled NCA summary statistics in Shin 2006 Table 2",
    "(see vignette Errata).",
    sep = " "
  )

  reference <- paste(
    "Shin JG, Kang WK, Shon JH, Arefayene M, Yoon YR, Kim KA,",
    "Kim DI, Kim DS, Cho KH, Woosley RL, Flockhart DA. (2007).",
    "Possible interethnic differences in quinidine-induced QT",
    "prolongation between healthy Caucasian and Korean subjects.",
    "British Journal of Clinical Pharmacology 63(2):206-215.",
    "doi:10.1111/j.1365-2125.2006.02793.x.",
    "Published OnlineEarly 10 November 2006.",
    sep = " "
  )

  vignette <- "Shin_2006_quinidine_QT"

  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "mg/L"
  )

  covariateData <- list(
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Source paper uses SEX with 1 = male, 0 = female (Shin 2006",
        "Table 3 footnote). Canonical SEXF inverts the value:",
        "SEXF = 1 - SEX. The E0 covariate effect '408 + 34*(1 - SEX)' in",
        "the source paper is then exactly '408 + 34 * SEXF' in canonical",
        "encoding; the +34 ms shift applies to females."
      ),
      source_name        = "SEX (values inverted: SEXF = 1 - SEX)"
    ),
    RACE_WHITE = list(
      description        = "White (Caucasian) race indicator, 1 = White, 0 = non-White (Korean in this cohort).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Korean; fETHN = 1)",
      notes              = paste(
        "Source paper uses ETHN with 1 = Korean, 0 = Caucasian and a",
        "derived multiplier fETHN of 1 for Korean / 1.26 for Caucasian",
        "(Shin 2006 Table 3 footnote). Canonical RACE_WHITE inverts the",
        "encoding: RACE_WHITE = 1 - ETHN. fETHN is computed inside",
        "model() as (1 + 0.26 * RACE_WHITE), so Korean has fETHN = 1 and",
        "Caucasian has fETHN = 1.26. The Cfemale = 106 ms additive",
        "interaction (Caucasian females only) is implemented as",
        "106 * SEXF * RACE_WHITE."
      ),
      source_name        = "ETHN (values inverted: RACE_WHITE = 1 - ETHN)"
    ),
    WT = list(
      description        = "Body weight at baseline (kg).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Used by the typical-value PK driver to scale the per-kg NCA",
        "summary values (CL = 0.3 L/h/kg, Vc = Vss = 2.5 L/kg from",
        "Shin 2006 Table 2 'Total' columns) to subject-level CL and",
        "Vc. Reference weight 70 kg with linear (exponent 1) scaling on",
        "both CL and Vc, equivalent to setting per-kg constants. WT is",
        "time-fixed per subject for this single-occasion dataset."
      ),
      source_name        = "BWT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 37L,
    n_studies      = 1L,
    age_range      = "Korean 21-29 years (M mean 22.1 +/- 1.6; F mean 22.7 +/- 2.4); Caucasian 21-37 years (M mean 26.2 +/- 7.5; F mean 27.7 +/- 3.6)",
    weight_range   = "Korean 49.7-73.8 kg (M mean 66.5 +/- 7.3; F mean 53.4 +/- 3.7); Caucasian 55.2-78.6 kg (M mean 69.8 +/- 8.8; F mean 60.7 +/- 5.5)",
    sex_female_pct = 49,
    race_ethnicity = c(Korean = 65, Caucasian = 35),
    disease_state  = paste(
      "Healthy volunteers screened for the absence of cardiovascular,",
      "hepatic, renal, neurological, and haematological abnormalities;",
      "normal ECG at baseline; not pregnant (urine pregnancy test in",
      "females); not taking oral contraceptives; female subjects",
      "scheduled within 5 days of the cessation of menses to minimise",
      "menstrual-cycle effects; abstained from alcohol, grapefruit",
      "juice, and caffeine for 3 weeks before and during the study."
    ),
    dose_range     = paste(
      "Single 4 mg/kg quinidine gluconate (base) IV infusion over 20 min",
      "in 20 mL normal saline via Harvard infusion pump. Reported",
      "actual doses (Shin 2006 Table 1): Korean M 266.1 +/- 29.4 mg,",
      "Korean F 213.7 +/- 14.6 mg, Caucasian M 279.2 +/- 35.2 mg,",
      "Caucasian F 242.7 +/- 22.2 mg. The cross-over arm used matching",
      "i.v. saline placebo with a 1 month washout."
    ),
    regions        = "Korea (Inje University Busan Paik Hospital) for Korean subjects; United States (Georgetown University Medical Center, Washington DC) for Caucasian subjects.",
    notes          = paste(
      "Randomised, double-blind crossover study (Shin 2006 Methods).",
      "Baseline QTc by sex was: Korean M 402 +/- 9 ms, Korean F 443",
      "+/- 8 ms (sex difference significant within Korean cohort,",
      "P < 0.05), Caucasian M 421 +/- 13 ms, Caucasian F 445 +/- 24 ms",
      "(Shin 2006 Table 2). NCA showed no statistically significant",
      "PK differences between ethnic groups: CLtot 0.28-0.34 L/h/kg,",
      "Vss 2.18-2.85 L/kg (Shin 2006 Table 2 'Total' columns). HPLC",
      "LLOQ 50 ng/mL with 2.5% CV; QT interval measured by Bazett's",
      "formula on 12-lead ECG. The model was developed in NONMEM with",
      "log-normal IIV and a constant-CV residual error."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Typical-value PK driver (Shin 2006 Table 2 'Total' NCA pooled
    # statistics). The paper does NOT fit a popPK model; these values
    # are fixed for simulation purposes only and are flagged as such
    # in the vignette Errata. CL = 0.30 L/h/kg and Vss = 2.5 L/kg are
    # the rounded pooled means; a 1-compartment IV approximation is
    # used because the paper does not report a peripheral volume or
    # inter-compartmental clearance.
    # ------------------------------------------------------------------
    lcl <- fixed(log(0.30 * 70))
    label("Clearance (L/h, reference WT 70 kg; NCA-derived, fixed)")
    # Shin 2006 Table 2 'Total' rows: CLtot = 0.31 (Korean) and 0.29
    # (Caucasian) L/h/kg; pooled ~0.30 L/h/kg. For 70 kg reference,
    # CL = 21 L/h. NOT a popPK fit -- see vignette Errata.

    lvc <- fixed(log(2.5 * 70))
    label("Central volume of distribution (L, reference WT 70 kg; = NCA Vss, fixed)")
    # Shin 2006 Table 2 'Total' rows: Vss = 2.78 (Korean) and 2.40
    # (Caucasian) L/kg; pooled ~2.5 L/kg. 1-cmt approximation uses
    # Vc = Vss; for 70 kg reference, Vc = 175 L. NOT a popPK fit --
    # see vignette Errata.

    # ------------------------------------------------------------------
    # Published PD Emax structural parameters (Shin 2006 Table 3
    # 'Adjusted model' row). Standard errors reported in parentheses
    # below each estimate.
    # ------------------------------------------------------------------
    lE0 <- log(408)
    label("Baseline QTc, male reference (ms)")
    # Shin 2006 Table 3 'Adjusted model' E0 = 408 + 34*(1-SEX); SE 7.9
    # on the male reference value (408). The +34 ms additive female
    # shift is encoded as a separate covariate effect below.

    lEmax <- log(136)
    label("Maximum QTc prolongation DeltaEmax, Korean reference (ms)")
    # Shin 2006 Table 3 'Adjusted model' DeltaEmax = 136*fETHN +
    # Cfemale; SE 18.7 on the Korean reference value (136). The
    # multiplicative ethnicity and additive interaction are encoded as
    # separate covariate-effect terms below.

    lEC50 <- log(3.13)
    label("EC50 (uM)")
    # Shin 2006 Table 3 'Adjusted model' EC50 = 3.13 uM, SE 0.71.
    # Source paper reports EC50 in uM; converted to mg/L inside model()
    # using quinidine MW = 324.42 g/mol so that the simulated Cc (mg/L)
    # can drive the Emax function directly.

    # ------------------------------------------------------------------
    # Covariate-effect coefficients (Shin 2006 Table 3 'Adjusted
    # model' covariate formulae). All three are fixed structural
    # values reported as part of the published equation -- the SEs in
    # Table 3 apply to E0 / DEmax / EC50 typical values, not to the
    # covariate coefficients themselves.
    # ------------------------------------------------------------------
    e_sexf_E0 <- fixed(34)
    label("Additive E0 shift in females (ms)")
    # Shin 2006 Table 3 footnote and Results paragraph after Table 2:
    # E0(ms) = 408 + 34*(1-SEX) with SEX = 1 male / 0 female, so the
    # +34 ms shift applies to female subjects. Implemented here as
    # 34 * SEXF after the source-to-canonical inversion.

    e_white_Emax <- fixed(0.26)
    label("Multiplicative DeltaEmax factor in Caucasians (unitless)")
    # Shin 2006 Table 3 footnote: fETHN = 1 for Koreans, 1.26 for
    # Caucasians. Implemented here as (1 + 0.26 * RACE_WHITE) after the
    # source-to-canonical inversion (RACE_WHITE = 1 - ETHN).

    e_white_sexf_Emax <- fixed(106)
    label("Additive DeltaEmax interaction in Caucasian females (ms)")
    # Shin 2006 Table 3 footnote and Results paragraph after Table 2:
    # Cfemale = 106 for Caucasian female, 0 otherwise. Implemented as
    # 106 * SEXF * RACE_WHITE.

    # ------------------------------------------------------------------
    # Inter-individual variability (log-normal). Variances are read
    # directly from Shin 2006 Table 3 'Adjusted model' omega^2 columns.
    # The paper reports log-normal distributions for the PD typical
    # values; the IIV is applied multiplicatively as exp(eta) around
    # the covariate-adjusted typical value in model().
    # ------------------------------------------------------------------
    etalE0   ~ 0.004   # Shin 2006 Table 3 'Adjusted model' w^2 E0 = 0.004 (SE 0.002).
    etalEmax ~ 0.0002  # Shin 2006 Table 3 'Adjusted model' w^2 Emax = 0.0002 (SE 0.004).
    etalEC50 ~ 0.48    # Shin 2006 Table 3 'Adjusted model' w^2 EC50 = 0.48 (SE 0.20).

    # ------------------------------------------------------------------
    # Residual error. Shin 2006 Methods describe a 'constant
    # coefficient of variation model' with sigma^2 reported in Table 3
    # 'Adjusted model' row as 0.004 (SE 0.003), giving a proportional
    # residual SD of sqrt(0.004) ~= 0.0632 (~6.3% CV on QTc).
    # ------------------------------------------------------------------
    propSd <- sqrt(0.004)
    label("Proportional residual error on QTc (fraction)")
  })

  model({
    # ================================================================
    # 1. Typical-value PK driver: 1-compartment IV with linear
    #    body-weight scaling on CL and Vc. Parameters are FIXED from
    #    NCA pooled summary statistics (Shin 2006 Table 2); this is
    #    NOT a popPK fit. The infusion rate / duration is supplied on
    #    the dose record by the user (see vignette).
    # ================================================================
    cl <- exp(lcl) * (WT / 70)
    vc <- exp(lvc) * (WT / 70)
    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Quinidine plasma concentration in mg/L.
    Cc <- central / vc

    # Convert Cc (mg/L) to uM to match the published EC50 unit.
    # 1 mg/L * (1000 ug/mg) / (324.42 ug/umol) = 3.0824 umol/L = 3.0824 uM.
    Cc_um <- Cc * 1000 / 324.42

    # ================================================================
    # 2. PD covariate-effect derived terms.
    #    fETHN = 1 (Korean) or 1.26 (Caucasian)
    #    Cfem  = 106 (Caucasian female) or 0 (other)
    # ================================================================
    fETHN <- 1 + e_white_Emax * RACE_WHITE
    Cfem  <- e_white_sexf_Emax * SEXF * RACE_WHITE

    # ================================================================
    # 3. Individual PD parameters. Typical values are the
    #    covariate-adjusted baseline / Emax; log-normal IIV is applied
    #    multiplicatively as exp(eta) around the typical value, per the
    #    source paper's lognormal-distribution assumption.
    # ================================================================
    E0_typ   <- exp(lE0) + e_sexf_E0 * SEXF
    E0_ind   <- E0_typ * exp(etalE0)

    Emax_typ <- exp(lEmax) * fETHN + Cfem
    Emax_ind <- Emax_typ * exp(etalEmax)

    EC50_ind <- exp(lEC50 + etalEC50)

    # ================================================================
    # 4. Emax pharmacodynamic model
    #    QTc(t) = E0 + DeltaEmax * Cc / (EC50 + Cc)
    #    Shin 2006 Eq. for E_t in 'Population pharmacodynamic
    #    analysis'. Cc is converted to uM (Cc_um) so the published
    #    EC50 in uM is used verbatim.
    # ================================================================
    QTc <- E0_ind + Emax_ind * Cc_um / (EC50_ind + Cc_um)

    # ================================================================
    # 5. Constant-CV (proportional) residual error on QTc.
    # ================================================================
    QTc ~ prop(propSd)
  })
}
