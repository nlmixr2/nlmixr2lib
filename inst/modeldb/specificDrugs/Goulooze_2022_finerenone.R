Goulooze_2022_finerenone <- function() {
  description <- "Population PKPD turnover model for serum potassium response to finerenone in patients with chronic kidney disease and type 2 diabetes (FIDELIO-DKD Phase III). Indirect-response model with an Emax effect of finerenone steady-state AUC on the potassium dissipation rate Kout, with a linear annual disease-progression slope on serum K (different typical value for active-treatment vs placebo arms). Finerenone PK is upstream (van den Berg 2022) and reduced here to AUCss = DOSE / CL with typical apparent clearance 28.0 L/h."
  reference <- "Goulooze SC, Snelder N, Seelmann A, Horvat-Broecker A, Brinker M, Joseph A, Garmann D, Lippert J, Eissing T. Finerenone Dose-Exposure-Serum Potassium Response Analysis of FIDELIO-DKD Phase III: The Role of Dosing, Titration, and Inclusion Criteria. Clin Pharmacokinet. 2022;61(3):451-462. doi:10.1007/s40262-021-01083-1. PK structural parameters are inherited from the FIDELIO-DKD popPK analysis (van den Berg JHE et al., Clin Pharmacol Drug Dev. 2022) and reduced here to the typical apparent CL = 28.0 L/h reported in Goulooze 2022 Fig 5 caption."
  vignette <- "Goulooze_2022_finerenone"
  units <- list(
    time          = "h",
    dosing        = "(oral finerenone, mg)",
    concentration = "mmol/L (serum potassium; the PD output molecule is potassium, not the dosed finerenone, so the dosing-vs-concentration unit-dimensional check is intentionally a non-applicable PD comparison and the dosing string is parenthesised to skip it)"
  )

  covariateData <- list(
    CRCL = list(
      description        = "Baseline CKD-EPI estimated glomerular filtration rate, BSA-normalised to 1.73 m^2",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed (baseline). Source column EGFR_EPI0. Enters as a power scaling on baseline serum K and on Emax, both referenced to 45 mL/min/1.73 m^2 (close to the FIDELIO-DKD median of 43.0).",
      source_name        = "EGFR_EPI0"
    ),
    UACR = list(
      description        = "Baseline urine albumin-to-creatinine ratio",
      units              = "mg/g",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed (baseline). Source column UACR0. Enters as a centred linear effect on Emax and on the disease-progression slope TSLOPE, both centred at 800 mg/g (close to the FIDELIO-DKD median of 852 mg/g).",
      source_name        = "UACR0"
    ),
    SEXF = list(
      description        = "Female-sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Source paper codes SEX = 1 (male), 2 (female); canonical SEXF maps as SEXF = SEX - 1. Enters as a multiplicative effect on Emax (-14.3% in females) and on baseline K with the coefficient fixed at 0 (paper $THETA(15) = 0 FIX).",
      source_name        = "(SEX - 1)"
    ),
    RACE_JAPANESE = list(
      description        = "Japanese-ancestry race indicator (1 = Japanese, 0 = otherwise)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Japanese)",
      notes              = "Source data column RACA = 3.2 codes Japanese ancestry (other races coded differently); canonical RACE_JAPANESE = (RACA == 3.2). Enters as a percent multiplicative shift on baseline K (-3.62% in Japanese subjects). The Japanese-specific residual-error multiplier (sigma^2 ratio 87.0%) reported in Goulooze 2022 Table 1 is documented in this entry but not reproduced in the typical-value Gaussian-residual approximation used here.",
      source_name        = "(RACA == 3.2)"
    ),
    ON_TREATMENT = list(
      description        = "Active-treatment arm indicator (1 = randomised to finerenone, 0 = randomised to placebo)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (placebo)",
      notes              = "Per-subject time-fixed. Source paper uses TREA = 1 (10 mg starting dose), 2 (20 mg starting dose), 3 (placebo); canonical ON_TREATMENT = (TREA < 3). Switches the disease-progression slope TSLOPE between the placebo typical value (0.00412 / year) and the active-arm typical value (0.00161 / year). Goulooze 2022 supplement applies the active TSLOPE only when TAFD > 0; here ON_TREATMENT is the per-subject randomisation flag and is the canonical generic on-treatment indicator (parallels Lee 2011 Parkinson's disease-progression usage).",
      source_name        = "(TREA < 3)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 10070L,
    n_studies      = 1L,
    age_range      = "adult (FIDELIO-DKD patients with advanced CKD and type 2 diabetes mellitus)",
    weight_range   = NA_character_,
    sex_female_pct = NA_real_,
    race_ethnicity = c(Japanese = NA_real_, NonJapanese = NA_real_),
    disease_state  = "Chronic kidney disease (CKD) with type 2 diabetes mellitus (T2DM); advanced CKD (eGFR 25 to < 75 mL/min/1.73 m^2 at run-in / screening) per FIDELIO-DKD inclusion criteria",
    dose_range     = "Oral finerenone 10 mg or 20 mg once daily, titrated based on serum potassium and eGFR per FIDELIO-DKD protocol; average dose over follow-up 15.1 mg",
    regions        = "Multi-regional Phase III (FIDELIO-DKD, NCT02540993)",
    notes          = "PD analysis informed by 148,384 serum potassium observations (62,401 local-lab and 85,983 central-lab) from 10,070 subjects: 5,674 randomised FIDELIO-DKD participants (2,841 placebo, 2,833 active-treatment of which 2,622 started on 10 mg and 211 on 20 mg) plus 4,396 subjects who participated in run-in and screening but did not meet the serum-potassium inclusion criterion. Median (5th-95th percentile) baseline eGFR 43.0 (26.7-66.9) mL/min/1.73 m^2 and UACR 852 (140-3,366) mg/g (Goulooze 2022 Sect. 3.1). Median follow-up 2.6 years."
  )

  ini({
    # ---------------------------------------------------------------------------
    # Upstream PK (van den Berg 2022 popPK fit on FIDELIO-DKD)
    # Goulooze 2022 does NOT fit PK; it uses individual posthoc PK estimates from
    # the upstream popPK analysis to compute the exposure metric
    #   AUCss = F1 * DOSE / CL
    # with F1 implicit in the apparent clearance (CL/F = 28.0 L/h typical). Only
    # the typical CL is reproduced here, as a fixed structural anchor.
    # ---------------------------------------------------------------------------
    lcl <- fixed(log(28.0)); label("Typical finerenone apparent clearance CL/F (L/h)")  # Goulooze 2022 Fig 5 caption: "a typical finerenone clearance of 28.0 L/h"

    # ---------------------------------------------------------------------------
    # Serum potassium turnover (indirect-response) parameters - typical values
    # All values are the FINAL estimates from Goulooze 2022 Table 1 ("Parameter
    # estimates and uncertainties of the final PKPD serum potassium model"); the
    # initial estimates in the supplement $THETA block are NOT used.
    # ---------------------------------------------------------------------------
    lbaseK     <- log(4.50);     label("Typical baseline serum potassium (mmol/L)")                                  # Goulooze 2022 Table 1: theta_pop,BSL = 4.50 mmol/L (RSE 0.140%)
    lkin       <- log(0.00981);  label("Zero-order serum K production rate Kin (mmol/L/h)")                          # Goulooze 2022 Table 1: theta_pop,Kin = 0.00981 mmol/L/h (RSE 14.2%)
    lemax      <- log(0.0905);   label("Maximum fractional drug effect Emax on Kout (fraction of Kout)")             # Goulooze 2022 Table 1: theta_pop,EMAX = 0.0905 (RSE 16.2%)
    lec50      <- log(0.512);    label("EC50 of finerenone AUCss on Kout (mg*h/L)")                                  # Goulooze 2022 Table 1: theta_pop,EC50 = 0.512 mg*h/L (RSE 33.3%)
    lhill      <- fixed(log(1)); label("Hill coefficient on AUCss (dimensionless)")                                  # Goulooze 2022 supplement $THETA(9): "1 FIX" (HILL fixed; reduces sigmoid Emax to Emax)

    # Disease-progression slopes on serum K, expressed per year
    ltslope_placebo <- log(0.00412); label("Annual fractional progression slope on serum K, placebo arm (1/year)")  # Goulooze 2022 Table 1: theta_pop,TSLOPE,placebo = 0.00412 /year (RSE 14.2%)
    ltslope_active  <- log(0.00161); label("Annual fractional progression slope on serum K, active arm (1/year)")   # Goulooze 2022 Table 1: theta_pop,TSLOPE,active  = 0.00161 /year (RSE 24.9%)

    # ---------------------------------------------------------------------------
    # Covariate effects
    # Sign / form notes:
    #   - EGFR-power: paper centres at EGFR = 45 mL/min/1.73 m^2.
    #   - UACR-linear: paper centres at UACR = 800 mg/g; coefficient units g/mg.
    #   - Japanese: paper expresses as percent shift (theta is the percent value;
    #     model multiplies by JAP/100, so theta = -3.62 -> -3.62% in Japanese).
    #   - SEX on EMAX: female-multiplicative deviation around male reference.
    # ---------------------------------------------------------------------------
    e_egfr_baseK    <- -0.0429;    label("Power exponent of (CRCL / 45) on baseline serum K")                       # Goulooze 2022 Table 1: theta_EGFR,BSL = -0.0429 (RSE 8.02%)
    e_jap_baseK_pct <- -3.62;      label("Japanese percent shift on baseline serum K (%)")                          # Goulooze 2022 Table 1: theta_JAP,BSL  = -3.62%  (RSE 10.1%)
    e_sexf_baseK    <- fixed(0);   label("Female-vs-male linear deviation on baseline serum K")                     # Goulooze 2022 supplement $THETA(15): "0 FIX" (SEX on BSL not retained)
    e_egfr_emax     <- -0.305;     label("Power exponent of (CRCL / 45) on Emax")                                   # Goulooze 2022 Table 1: theta_EGFR,EMAX = -0.305 (RSE 23.6%)
    e_uacr_emax     <- 0.0000931;  label("Linear effect of (UACR - 800) on Emax (g/mg)")                            # Goulooze 2022 Table 1: theta_UACR,EMAX = 9.31e-5 g/mg (RSE 23.7%)
    e_sexf_emax     <- -0.143;     label("Female-vs-male linear deviation on Emax")                                 # Goulooze 2022 Table 1: theta_SEX,EMAX  = -0.143 (RSE 20.1%)
    e_uacr_tslope   <- 0.00114;    label("Linear effect of (UACR - 800) on TSLOPE (g/mg)")                          # Goulooze 2022 Table 1: theta_UACR,TSLOPE = 1.14e-3 g/mg (RSE 19.6%)
    e_baseK_tslope  <- 1.60;       label("Linear effect of (baseK - 4.4) on TSLOPE (L/mmol)")                       # Goulooze 2022 Table 1: theta_BSL,TSLOPE  = 1.60 L/mmol (RSE 10.8%)

    # ---------------------------------------------------------------------------
    # Inter-individual variability
    # Goulooze 2022 supplement applies a Box-Cox-transformed exponential IIV to
    # the typical baseline K:
    #     ETATR = (exp(ETA(1))^bxpar - 1) / bxpar
    #     BSL   = TVBSL * exp(ETATR)
    # with the Box-Cox shape parameter bxpar estimated. Emax carries proportional
    # IIV: EMAX = TVEMAX * (1 + ETA(2)). The two etas share a 2x2 omega block.
    # ---------------------------------------------------------------------------
    bxpar_baseK <- -1.61;          label("Box-Cox shape parameter on the exponential IIV of baseline K")            # Goulooze 2022 Table 1: theta_boxcox,IIV,BSL = -1.61 (RSE 17.4%)
    etalbaseK + etalemax ~ c(0.00717,
                             -0.0385, 1.49)                                                                          # Goulooze 2022 Table 1: omega^2 BSL (exp) = 0.00717 (RSE 2.43%); cov BSL/EMAX = -0.0385 (RSE 10.8%); omega^2 EMAX (prop) = 1.49 (RSE 10.2%)

    # ---------------------------------------------------------------------------
    # Residual error
    # Goulooze 2022 uses a Student-t-distributed proportional residual with df =
    # 6.60 and sigma^2 = 0.00447 (paper Table 1). The Japanese-ancestry sigma^2
    # multiplier (87.0%) is documented in covariateData[[RACE_JAPANESE]]$notes
    # but not reproduced in this Gaussian-residual approximation. The Student-t
    # distribution is approximated by Gaussian here for typical-value simulation;
    # the SD coefficient is sqrt(0.00447) = 0.0668 (~ 6.7% proportional CV).
    # ---------------------------------------------------------------------------
    propSd <- sqrt(0.00447); label("Proportional residual SD (Gaussian approximation of paper's Student-t scalar)")  # Goulooze 2022 Table 1: sigma^2 = 0.00447 (RSE 0.986%); paper df = 6.60 (RSE 1.81%)
  })

  model({
    # Finerenone apparent clearance (typical, upstream-fixed at 28.0 L/h)
    cl <- exp(lcl)

    # -----------------------------------------------------------------------
    # Typical baseline serum potassium with covariate effects
    # Source equation (Goulooze 2022 supplement, $PK block):
    #   TVBSL = THETA(2) * (1 + JAP*THETA(5)/100) * (EGFREPI0/45)^THETA(7) * CV3
    # where CV3 = 1 + (SEX-1)*THETA(15); since THETA(15) = 0 FIX, CV3 == 1.
    # -----------------------------------------------------------------------
    baseK_typ <- exp(lbaseK) *
      (1 + RACE_JAPANESE * e_jap_baseK_pct / 100) *
      (CRCL / 45)^e_egfr_baseK *
      (1 + SEXF * e_sexf_baseK)

    # Box-Cox-transformed individual baseline K.
    # Source (supplement): ETATR = (exp(ETA(1))^BXPAR - 1) / BXPAR;
    #                      BSL    = TVBSL * exp(ETATR).
    etatr <- (exp(etalbaseK)^bxpar_baseK - 1) / bxpar_baseK
    baseK <- baseK_typ * exp(etatr)

    # Turnover constants. Kout is derived so that the unperturbed steady state
    # equals BSL (Kin / Kout = BSL).
    kin  <- exp(lkin)
    kout <- kin / baseK

    # -----------------------------------------------------------------------
    # Emax with covariate-modified typical Emax
    # Source (supplement, $PK block):
    #   EMAX = THETA(4) * (1 + ETA(2)) * CV1 * (EGFREPI0/45)^THETA(13) * CV2
    # with CV1 = 1 + (UACR0 - 800)*THETA(12) and CV2 = 1 + (SEX-1)*THETA(14).
    # -----------------------------------------------------------------------
    emax_typ <- exp(lemax) *
      (1 + (UACR - 800) * e_uacr_emax) *
      (CRCL / 45)^e_egfr_emax *
      (1 + SEXF * e_sexf_emax)
    emax <- emax_typ * (1 + etalemax)

    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # -----------------------------------------------------------------------
    # Disease-progression slope (annual fractional rate)
    # Source (supplement): TSLOPE = THETA(16) * (1+(BSL-4.4)*THETA(17)) * CV4
    #   for placebo, and THETA(19) instead of THETA(16) for active arm when
    #   TAFD > 0. CV4 = 1 + (UACR0 - 800)*THETA(18). The placebo-vs-active
    #   switch is encoded here by the per-subject ON_TREATMENT indicator.
    # -----------------------------------------------------------------------
    tslope_arm <- (1 - ON_TREATMENT) * exp(ltslope_placebo) +
                  ON_TREATMENT       * exp(ltslope_active)
    tslope <- tslope_arm *
      (1 + (baseK - 4.4) * e_baseK_tslope) *
      (1 + (UACR  - 800) * e_uacr_tslope)

    # Cumulative annual progression effect (rxode2 t is in the model time unit
    # declared in units$time = "h"; convert hours to years).
    teff <- tslope * t / (24 * 365.25)

    # -----------------------------------------------------------------------
    # Finerenone exposure metric: AUCss = F * DOSE / CL with F implicit.
    # podo(depot) returns the most recent dose amount entered into the depot
    # compartment, yielding a step-function AUCss that changes at each
    # titration / interruption / re-initiation event.
    # -----------------------------------------------------------------------
    aucSs <- podo(depot) / cl
    eff   <- emax * aucSs^hill / (ec50^hill + aucSs^hill)
    if (aucSs <= 0) eff <- 0

    # -----------------------------------------------------------------------
    # ODE system
    #
    # depot is a virtual dose receiver: doses (cmt = "depot", amt = mg) enter
    # it, and we read the dose amount via podo(depot). The compartment state
    # itself is not used by any other equation, so its derivative is zero;
    # the integrated state grows as the cumulative dose but is irrelevant
    # downstream. For an interruption / discontinuation, the simulation must
    # send an explicit amt = 0 dose event into the depot to clear the metric.
    # -----------------------------------------------------------------------
    d/dt(depot)  <- 0
    d/dt(serumK) <- kin - kout * (1 - eff) * (1 - teff) * serumK
    serumK(0)    <- baseK

    serumK ~ prop(propSd)
  })
}
