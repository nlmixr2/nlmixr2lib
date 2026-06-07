You_2013_hCG_GTN <- function() {
  description <- "Population kinetic biomarker decline model for serum human chorionic gonadotrophin (hCG) during 8-day methotrexate (MTX) therapy in low-risk gestational trophoblastic neoplasia (GTN). hCG decays mono-exponentially from an initial amplitude hCG0 toward a non-zero asymptote hCGres (the modelled residual production attributable to MTX-insensitive tumour cells). The same kinetic structure recasts as an indirect-response / turnover model with kout = K and zero-order production kin = K * hCGres, so the canonical lkout / lrbase names are used. No drug PK is modelled (the fixed 8-day MTX/folinic-acid regimen is implicit). Inter-individual variability: Box-Cox transformed eta on the hCG0 amplitude (shape -0.207) and log-normal etas on lkout and lrbase. Proportional residual error. M3 censoring of BLQ titres (< 2 IU/L) used in the original fit but is not reproduced here -- the packaged model is intended for typical-value + IIV simulation."
  reference <- paste(
    "You B, Harvey R, Henin E, Mitchell H, Golfier F, Savage PM, Tod M,",
    "Wilbaux M, Freyer G, Seckl MJ.",
    "Early prediction of treatment resistance in low-risk gestational",
    "trophoblastic neoplasia using population kinetic modelling of hCG",
    "measurements.",
    "Br J Cancer. 2013 May 14;108(9):1810-1816.",
    "doi:10.1038/bjc.2013.123.",
    "Refines the mono-exponential hCG kinetic model of You et al. 2010",
    "(Ann Oncol 21:1643-1650) by adding the hCGres residual-production",
    "term and validates it in an independent 800-patient Charing Cross",
    "Hospital cohort.",
    sep = " "
  )
  vignette <- "You_2013_hCG_GTN"
  paper_specific_etas <- c("eta_hCG0")
  paper_specific_compartments <- c("hCG")

  units <- list(
    time = "day",
    dosing = "IU/L",
    concentration = "IU/L"
  )

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 418L,
    n_studies      = 1L,
    age_range      = "median 31.2 years (IQR 26.1-36.1) in Model set; 30.9 years (IQR 26.1-35.6) in Test set (Table 1)",
    weight_range   = "not reported",
    sex_female_pct = 100,
    race_ethnicity = NULL,
    disease_state  = "Low-risk gestational trophoblastic neoplasia (FIGO 2000 risk score 0-6) including invasive mole, choriocarcinoma, epithelioid trophoblastic and placental site trophoblastic tumours; cases of MTX resistance defined as three static or three rising hCG measurements (within +/-10%) on the 8-day MTX regimen.",
    dose_range     = "8-day methotrexate regimen: 50 mg MTX intramuscular on days 1, 3, 5, 7 with 15 mg oral folinic acid on days 2, 4, 6, 8; cycles repeated every 2 weeks until hCG normalisation plus three consolidation cycles, or MTX resistance (Patients and methods).",
    regions        = "United Kingdom (Charing Cross Hospital Trophoblastic Disease Centre, London, 1991-2011)",
    notes          = "Model parameters in Table 2 estimated from the 418-patient Model data set (195 MTX-resistant, 223 MTX-sensitive); a separate 382-patient Test set validated the predictive 20.44 IU/L hCGres threshold. NONMEM v7 with M3 likelihood for BLQ (< 2 IU/L) hCG measurements over treatment days 0-50; titres beyond day 50 or after second-line switch were censored. Precision of estimates computed from a 100-resampling bootstrap (Table 2 footnote a)."
  )

  ini({
    # Typical-value hCG0 (initial decline amplitude, IU/L). Linear scale because
    # the IIV is implemented via a Box-Cox transformed eta (Table 2 footnote b),
    # not a log-normal eta on a log-transformed theta.
    hCG0_pop <- 14400; label("Typical hCG0 initial decline amplitude (IU/L)")  # Table 2 row 1, THETA(1) = 14400 (RSE 5.7%)

    # K -- first-order decline rate of the hCG above-residual pool (1/day).
    # Recast as the canonical IDR kout: d/dt(hCG) = -K*(hCG - hCGres) is
    # equivalent to d/dt(hCG) = kin - kout*hCG with kout = K and kin = K*hCGres.
    lkout <- log(0.169); label("log decline rate constant K (1/day)")  # Table 2 row 2, THETA(2) = 0.169/day (RSE 1.3%)

    # hCGres -- residual hCG production attributable to MTX-insensitive cells
    # (IU/L); equivalent to the IDR baseline asymptote rbase = kin/kout.
    lrbase <- log(24.7); label("log residual hCG production hCGres (IU/L)")  # Table 2 row 3, THETA(3) = 24.7 IU/L (RSE 11%)

    # Inter-individual variability.
    # eta_hCG0 -- Box-Cox eta on hCG0 (Table 2 footnote b):
    #   BCOX1 = ((exp(ETA(1)))^SHAPE1 - 1) / SHAPE1
    #   hCG0  = THETA(1) * exp(BCOX1)
    # with SHAPE1 estimated at -0.207 (s.e. 7.6%) and ETA(1) ~ N(0, 4.24).
    # The CV of hCG0 reported in Table 2 (152%) is the cohort statistic of the
    # resulting Box-Cox distribution, NOT a log-normal omega; do not back-solve
    # omega^2 = log(CV^2 + 1) here.
    eta_hCG0 ~ 4.24  # Table 2 footnote b, OMEGA(1,1) = 4.24 (RSE 3.5%)

    # log-normal etas on K and hCGres. Table 2 reports them as CV%; omega^2 on
    # the log scale = log(CV^2 + 1).
    #   K     CV 306%: omega^2 = log(3.06^2 + 1) = log(10.3636) = 2.338
    #   hCGres CV 35%: omega^2 = log(0.35^2 + 1) = log(1.1225)  = 0.1156
    etalkout  ~ 2.338  # Table 2 row 2, IIV CV = 306% (RSE 6.7%)
    etalrbase ~ 0.1156  # Table 2 row 3, IIV CV = 35% (RSE 4.9%)

    # Residual error: proportional on the predicted hCG.
    # Paper equation: hCG_ij(t) = (hCG0_i * exp(-K_i * t) + hCGres_i) * (1 + eps1_ij)
    # with eps1_ij ~ N(0, sigma^2). Table 2 reports the residual SD as 34.4%.
    propSd <- 0.344; label("Proportional residual error on hCG (fraction)")  # Table 2 row 1, SIGMA = 34.4% (RSE 1.4%)
  })

  model({
    # Box-Cox shape parameter for the hCG0 individual distribution. Estimated
    # in the original fit and frozen here at the published point estimate.
    # Kept as a local constant rather than an ini() theta because it has no
    # canonical parameter name in nlmixr2lib conventions (Box-Cox shape is
    # paper-specific to this model class).
    shape_hCG0 <- -0.207  # Table 2 footnote b: SHAPE1 = -0.207 (RSE 7.6%)

    # Individual hCG0 via Box-Cox eta (Table 2 footnote b).
    bcox <- ((exp(eta_hCG0))^shape_hCG0 - 1) / shape_hCG0
    hCG0 <- hCG0_pop * exp(bcox)

    # Individual K and hCGres via log-normal etas.
    K      <- exp(lkout  + etalkout)
    hCGres <- exp(lrbase + etalrbase)

    # IDR/turnover-equivalent form of the paper's analytical solution
    #   hCG(t) = hCG0 * exp(-K*t) + hCGres
    # is the single-state ODE
    #   d/dt(hCG) = -K * (hCG - hCGres),  hCG(0) = hCG0 + hCGres
    # which expands to d/dt(hCG) = K*hCGres - K*hCG = kin - kout*hCG with
    # kin = K*hCGres and kout = K (the canonical IDR identification).
    d/dt(hCG) <- -K * (hCG - hCGres)
    hCG(0)    <- hCG0 + hCGres

    hCG ~ prop(propSd)
  })
}
