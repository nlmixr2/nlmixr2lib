Goulooze_2022_finerenone_uacr <- function() {
  description <- "Population PKPD disease-progression model for the urine albumin-to-creatinine ratio (UACR) response to finerenone in patients with chronic kidney disease and type 2 diabetes (FIDELIO-DKD Phase III). UACR is integrated as a state whose fractional progression rate is corrected by the model-predicted UACR and by UACR over baseline; the finerenone effect is a power function of steady-state daily AUC acting through an effect compartment, and concomitant SGLT2 inhibitor use gives a direct proportional shift. Finerenone PK is upstream (van den Berg 2022) and reduced here to AUCss = DOSE / CL with typical apparent clearance 29.9 L/h."
  reference <- paste(
    "Goulooze SC, Heerspink HJL, van Noort M, Snelder N, Brinker M, Lippert J,",
    "Eissing T. Dose-Exposure-Response Analysis of the Nonsteroidal",
    "Mineralocorticoid Receptor Antagonist Finerenone on UACR and eGFR:",
    "An Analysis from FIDELIO-DKD. Clin Pharmacokinet. 2022;61(7):1023-1037.",
    "doi:10.1007/s40262-022-01124-3.",
    "Parameter values are the final estimates in Table S1 of the Electronic",
    "Supplementary Material; the model equations are the final UACR NONMEM",
    "control stream in the same supplement. The exposure metric AUCtau,md is",
    "computed from the upstream FIDELIO-DKD population PK analysis",
    "(van den Berg P et al., Clin Pharmacokinet. 2022;61(7):1005-1021;",
    "doi:10.1007/s40262-021-01082-2); see modellib('vandenBerg_2021_finerenone').",
    sep = " "
  )
  vignette <- "Goulooze_2022_finerenone_uacr_egfr"
  units <- list(
    time          = "h",
    dosing        = "(oral finerenone, mg)",
    concentration = "mg/g (UACR; the PD output is an endogenous urinary damage biomarker rather than the dosed finerenone, so the dosing-vs-concentration dimensional check is not applicable and the dosing string is parenthesised to skip it)"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    depot  = list(analyte = "finerenone", units = "mg",   specimen = "administration site", verified = TRUE),
    uacr   = list(analyte = "albumin/creatinine ratio", units = "mg/g", specimen = "urine", verified = TRUE),
    effect = list(analyte = "finerenone", units = "mg*h/L", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    UACR = list(
      description        = "Observed baseline urine albumin-to-creatinine ratio",
      units              = "mg/g",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed (baseline). Source column UACR0. Enters as a power scaling (UACR / 850)^0.877 on the typical baseline UACR state. The centring constant 850 mg/g is the control-stream value, close to the FIDELIO-DKD cohort median of 852 mg/g.",
      source_name        = "UACR0"
    ),
    CRCL = list(
      description        = "Observed baseline CKD-EPI estimated glomerular filtration rate, BSA-normalised to 1.73 m^2",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed (baseline). Source column EGFREPI0. Enters twice: as a power scaling (CRCL / 43)^-0.124 on the typical baseline UACR, and as a centred linear shift -0.00257 * (CRCL - 43) on the annual UACR progression rate. The reference 43 mL/min/1.73 m^2 is the FIDELIO-DKD cohort median.",
      source_name        = "EGFREPI0"
    ),
    HEPIMP_MOD = list(
      description        = "Moderate hepatic impairment indicator (1 = likely or certain Child-Pugh B, 0 = likely Child-Pugh A or healthy)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (likely Child-Pugh A or healthy)",
      notes              = "Classification scheme is Child-Pugh Class B, derived in FIDELIO-DKD from laboratory surrogates rather than a full Child-Pugh score: subjects with total bilirubin < 2 mg/dL AND serum albumin > 3.5 g/dL were categorised as likely Child-Pugh A or healthy, all others as likely or certain Child-Pugh B (Goulooze 2022 Table 1 footnote b). Source column CHILDPSC, with values 2 or 3 mapping to HEPIMP_MOD = 1 and values 1 or 5 to 0. Enters as a proportional shift (1 + 0.0943 * HEPIMP_MOD) on the typical baseline UACR.",
      source_name        = "(CHILDPSC %in% c(2, 3))"
    ),
    RACE_ASIAN = list(
      description        = "Asian-race indicator (1 = any Asian race category, 0 = otherwise)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Asian)",
      notes              = "Source column RACEASIA, a decimal-coded race variable; values in [3.0, 3.9] denote Asian race categories and map to RACE_ASIAN = 1. Enters as an additive shift +0.0634 /year on the annual UACR progression rate.",
      source_name        = "(RACEASIA >= 3 & RACEASIA <= 3.9)"
    ),
    RACE_JAPANESE = list(
      description        = "Japanese-heritage race indicator (1 = Japanese, 0 = otherwise)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Japanese)",
      notes              = "Source column RACEASIA = 3.2 codes Japanese heritage (selected in the control stream by 3.15 <= RACEASIA <= 3.25). Enters as a proportional shift (1 - 0.261 * RACE_JAPANESE) on the finerenone drug-effect slope ESLOPE. Japanese subjects are a subset of RACE_ASIAN, so both indicators are 1 for a Japanese subject.",
      source_name        = "(RACEASIA == 3.2)"
    ),
    AGE = list(
      description        = "Subject age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed (baseline). Enters as a power scaling (AGE / 63)^0.864 on the finerenone drug-effect slope ESLOPE; 63 years is the control-stream centring value.",
      source_name        = "AGE"
    ),
    CONMED_SGLT2I = list(
      description        = "Current concomitant SGLT2 inhibitor use indicator (1 = in use, 0 = not in use)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant SGLT2 inhibitor)",
      notes              = "TIME-VARYING within subject. FIDELIO-DKD recorded SGLT2 inhibitor use as a binary variable over time defined as use / no use in the last 5 days, without regard to the specific agent or its dosing schedule (Goulooze 2022 Sect. 2.2.3). Enters as a direct proportional effect exp(-0.212 * (CONMED_SGLT2I - CONMED_SGLT2I_BASE)) on the observed UACR; note it acts on the OBSERVATION, not on the UACR differential equation. 528 of 5674 patients used an SGLT2 inhibitor at some point during the treatment period.",
      source_name        = "FLAGSGLT"
    ),
    CONMED_SGLT2I_BASE = list(
      description        = "Concomitant SGLT2 inhibitor use at treatment start (1 = in use at randomisation, 0 = not in use)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no SGLT2 inhibitor at treatment start)",
      notes              = "Time-fixed per subject. The control stream applies the SGLT2 inhibitor effect to the CHANGE from the treatment-start value, as SGLTEFF * (FLAGSGLT - SGLTSTART), because the estimated baseline UACR of a subject already on an SGLT2 inhibitor at randomisation already reflects that drug's effect. Set to 0 for a subject who is SGLT2i-naive at randomisation.",
      source_name        = "SGLTSTART"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 5674L,
    n_studies      = 1L,
    age_range      = "adults; median 66-67 years (IQR approximately 57-72) across treatment groups",
    weight_range   = "median 84.6-87.6 kg (IQR approximately 72.6-99.4) across treatment groups",
    sex_female_pct = 29.4,
    race_ethnicity = c(White = 63.2, Black = 4.6, Japanese = 7.3, Chinese = 10.3, Other = 14.4),
    disease_state  = "Chronic kidney disease with type 2 diabetes mellitus; baseline eGFR median 42.1-51.3 mL/min/1.73 m^2 and baseline UACR median 661-887 mg/g across treatment groups",
    dose_range     = "Oral finerenone 10 mg or 20 mg once daily (starting dose 10 mg if screening eGFR 25 to < 60 mL/min/1.73 m^2, 20 mg if >= 60), up- and down-titrated per eGFR and serum potassium, versus matching placebo; all in addition to standard of care with a maximally tolerated labelled dose of a renin-angiotensin system inhibitor",
    regions        = "Multi-regional Phase III (FIDELIO-DKD, NCT02540993)",
    notes          = "37,296 UACR observations from 5,674 patients (full analysis set). 528 patients used an SGLT2 inhibitor at some point during the treatment period. Baseline demographics from Goulooze 2022 Table 1; the four columns of that table are placebo without SGLT2i (N = 2553), placebo with SGLT2i (N = 288), finerenone without SGLT2i (N = 2593) and finerenone with SGLT2i (N = 240), and the ranges quoted here span those four groups. Sex percentage is the pooled female fraction across all four groups."
  )

  ini({
    # -------------------------------------------------------------------------
    # Upstream PK (van den Berg 2022 popPK fit on FIDELIO-DKD; reference [10])
    # This paper fits NO PK. Its control stream computes the exposure metric
    # algebraically as AUCT = IF1 * DOSE / ICL, where IF1 and ICL are individual
    # post-hoc bioavailability and clearance supplied as DATA COLUMNS from the
    # upstream fit. At the upstream reference covariates F1 = 1 and CL/F = 29.9
    # L/h (van den Berg 2022 Table 2), so AUCtau,md = DOSE / 29.9 mg*h/L. Only
    # the typical clearance is reproduced here, as a fixed structural anchor.
    # -------------------------------------------------------------------------
    lcl <- fixed(log(29.9)); label("Typical finerenone apparent clearance CL/F, upstream value (L/h)")  # van den Berg 2022 Table 2: CL/F = 29.9 L/h (RSE 3.62%); see modellib('vandenBerg_2021_finerenone')

    # -------------------------------------------------------------------------
    # UACR disease progression - typical values
    # All values are the FINAL estimates from Goulooze 2022 ESM Table S1. The
    # $THETA block of the supplement's control stream holds INITIAL estimates
    # (e.g. TH1 initial 860.473 vs final 866; TH3 initial 0.0990511 vs final
    # 0.137) and is NOT used for values, only for the equation forms.
    # -------------------------------------------------------------------------
    lrbase   <- log(866);      label("Typical baseline UACR (mg/g)")                                       # ESM Table S1: theta_pop,BSLUACR = 866 mg/g (RSE 0.429%)
    prog     <- 0.137;         label("Typical UACR progression rate (1/year)")                             # ESM Table S1: theta_pop,PROG = 0.137 /year (RSE 9.10%)
    cslope1  <- 0.0000720;     label("Correction of progression by model-predicted UACR (g/mg/year)")      # ESM Table S1: theta_pop,CSLOPE1 = 7.20e-5 g/mg/year (RSE 8.39%)
    cslope2  <- 0.182;         label("Correction of progression by log UACR over baseline (1/year)")       # ESM Table S1: theta_pop,CSLOPE2 = 0.182 /year (RSE 21.2%)

    # -------------------------------------------------------------------------
    # Finerenone effect on UACR: a power function of the effect-compartment
    # exposure, attenuated at high UACR through eslope2 (control stream $DES:
    # EFFECT = ESLOPE * AUCEFF**EPOW * EXP((A(1) - 850) * ESLOPE2)).
    # -------------------------------------------------------------------------
    eslope   <- -1.42;         label("Slope of the finerenone effect on log UACR (L/mg/h)")                # ESM Table S1: theta_pop,ESLOPE = -1.42 (RSE 6.61%)
    epow     <- 0.613;         label("Power of the finerenone exposure in the UACR effect (unitless)")     # ESM Table S1: theta_pop,EPOW = 0.613 (RSE 3.99%)
    eslope2  <- -0.000166;     label("Attenuation of the finerenone UACR effect by current UACR (g/mg)")   # ESM Table S1: theta_pop,ESLOPE2 = -1.66e-4 g/mg (RSE 34.5%)
    lke0     <- log(0.000136); label("Equilibration rate of the finerenone UACR effect compartment (1/h)") # ESM Table S1: theta_pop,KE0,UACR = 1.36e-4 /h (RSE 16.1%)

    # -------------------------------------------------------------------------
    # Covariate effects
    #   - power form on baseline UACR and on the drug-effect slope;
    #   - centred linear form on the progression rate;
    #   - proportional (1 + theta * IND) form for the binary indicators.
    # -------------------------------------------------------------------------
    e_uacr_rbase           <- 0.877;    label("Power exponent of (UACR / 850) on baseline UACR")                  # ESM Table S1: theta_UACR0,BSLUACR = 0.877 (RSE 0.546%)
    e_crcl_rbase           <- -0.124;   label("Power exponent of (CRCL / 43) on baseline UACR")                   # ESM Table S1: theta_EGFREPI0,BSLUACR = -0.124 (RSE 11.4%)
    e_hepimp_mod_rbase     <- 0.0943;   label("Proportional effect of Child-Pugh B on baseline UACR")             # ESM Table S1: theta_CHILDPUGH,BSLUACR = 0.0943 (RSE 17.4%)
    e_crcl_prog            <- -0.00257; label("Linear effect of (CRCL - 43) on the progression rate (1/year per mL/min/1.73 m^2)")  # ESM Table S1: theta_EGFREPI0,PROG = -0.00257 (RSE 23.0%)
    e_race_asian_prog      <- 0.0634;   label("Additive effect of Asian race on the progression rate (1/year)")   # ESM Table S1: theta_ASIAN,PROG = 0.0634 /year (RSE 25.8%)
    e_age_eslope           <- 0.864;    label("Power exponent of (AGE / 63) on the drug-effect slope")            # ESM Table S1: theta_AGE,ESLOPE = 0.864 (RSE 19.4%)
    e_race_japanese_eslope <- -0.261;   label("Proportional effect of Japanese race on the drug-effect slope")    # ESM Table S1: theta_JAPAN,ESLOPE = -0.261 (RSE 21.7%)
    e_conmed_sglt2i_uacr   <- -0.212;   label("Effect of current SGLT2 inhibitor use on log UACR")                # ESM Table S1: theta_SGLT2i,UACR = -0.212 (RSE 13.7%)

    # -------------------------------------------------------------------------
    # Inter-individual variability
    # The source applies Box-Cox-transformed exponential IIV to the progression
    # rate and to the drug-effect slope (control stream $PK:
    #   ETATR = (EXP(ETA)**BXPAR - 1) / BXPAR),
    # with plain exponential IIV on baseline UACR. The published OMEGA is a
    # 4x4 block over (BSLUACR, PROG, ERROR, ESLOPE); the ERROR row/column is
    # omitted here because nlmixr2 cannot carry IIV on the residual magnitude
    # (see the residual-error note below). The retained 3x3 sub-block is
    # reproduced exactly, including its off-diagonals.
    #
    # Published values not reproduced here (ESM Table S1): omega^2 ERROR =
    # 0.0402 (RSE 8.83%), BOXCOX ERROR = 5.89 (RSE 8.37%), and the covariances
    # BSL/ERROR = -0.000317 (RSE 406%), PROG/ERROR = -0.00894 (RSE 27.5%) and
    # ERROR/ESLOPE = -0.0254 (RSE 34.8%).
    # -------------------------------------------------------------------------
    boxcox_prog   <- -0.545; label("Box-Cox shape parameter for the IIV on the progression rate")   # ESM Table S1: BOXCOX PROG = -0.545 (RSE 10.2%)
    boxcox_eslope <- -0.955; label("Box-Cox shape parameter for the IIV on the drug-effect slope")  # ESM Table S1: BOXCOX ESLOPE = -0.955 (RSE 14.6%)

    etalrbase + etaprog + etaeslope ~ c(0.0166,
                                        -0.00315, 0.235,
                                        -0.0354,  0.133, 1.03)  # ESM Table S1: omega^2 BSLUACR = 0.0166 (RSE 12.6%); omega^2 PROG = 0.235 (RSE 12.0%); omega^2 ESLOPE = 1.03 (RSE 11.9%); cov BSL/PROG = -0.00315 (RSE 114%); cov BSL/ESLOPE = -0.0354 (RSE 28.7%); cov PROG/ESLOPE = 0.133 (RSE 18.3%)

    # -------------------------------------------------------------------------
    # Residual error
    # The source residual is additive on the log scale (control stream $ERROR:
    #   IPRED = LOG(A(1)) + EFFECT + SGLTEFF*(FLAGSGLT - SGLTSTART)
    #   Y     = IPRED + ERR(1) * VARCOR * EXP(ETA3S)),
    # which is exactly a log-normal residual on the linear-scale prediction, so
    # it is encoded as ~ lnorm(expSd) with expSd = sqrt(sigma^2).
    #
    # Two features of the published residual model are NOT reproducible in
    # nlmixr2, which requires the residual SD to be a bare estimated parameter
    # rather than a model expression:
    #   (1) VARCOR = (IPRED_linear / 850)^theta_IPRED,ERROR with
    #       theta_IPRED,ERROR = -0.227 (RSE 3.06%), a power scaling of the
    #       residual SD by the current prediction; and
    #   (2) EXP(ETA3S), the Box-Cox-transformed IIV on the residual magnitude.
    # Both equal 1 at the reference prediction of 850 mg/g and at the typical
    # value of the eta respectively, so typical-value simulation is unaffected;
    # only the heteroscedasticity and the between-subject spread of the
    # residual are lost. See the vignette's Assumptions and deviations section.
    # -------------------------------------------------------------------------
    expSd <- sqrt(0.141); label("Log-scale residual SD for UACR (log mg/g)")  # ESM Table S1: sigma^2 = 0.141 (RSE 1.53%), residual error additive on log-scale
  })

  model({
    # -----------------------------------------------------------------------
    # Time-unit conversion. The control stream converts annual rates to the
    # model time unit with 365 * 24 hours per year (the eGFR model uses
    # 365.25; both are reproduced as published).
    # -----------------------------------------------------------------------
    hoursPerYear <- 365 * 24

    # Finerenone apparent clearance (typical, upstream-fixed at 29.9 L/h)
    cl <- exp(lcl)

    # -----------------------------------------------------------------------
    # Individual baseline UACR.
    # Control stream $PK:
    #   CV2 = (UACR0/850)**THETA(2); CV6 = (EGFREPI0/43)**THETA(5)
    #   ECPGH = THETA(7)*SCPGH
    #   BSL = THETA(1) * EXP(ETA(1)) * CV2 * CV6 * (1 + ECPGH)
    # -----------------------------------------------------------------------
    rbase <- exp(lrbase + etalrbase) *
      (UACR / 850)^e_uacr_rbase *
      (CRCL / 43)^e_crcl_rbase *
      (1 + e_hepimp_mod_rbase * HEPIMP_MOD)

    # -----------------------------------------------------------------------
    # Individual progression rate, per hour.
    # Control stream $PK:
    #   CV5 = (EGFREPI0-43)*THETA(4); ERACA = THETA(6)*SRACA
    #   TVPROG = THETA(3) + CV5 + ERACA
    #   ETATR2 = (EXP(ETA(2))**BXPAR2 - 1) / BXPAR2
    #   PROG   = (TVPROG + ETATR2) / (365*24)
    # The eta enters ADDITIVELY on the annual rate, which may therefore be
    # negative (a subject whose UACR falls under standard of care).
    # -----------------------------------------------------------------------
    etatrProg <- (exp(etaprog)^boxcox_prog - 1) / boxcox_prog
    progRate  <- (prog +
                    e_crcl_prog * (CRCL - 43) +
                    e_race_asian_prog * RACE_ASIAN +
                    etatrProg) / hoursPerYear

    cslope1Hr <- cslope1 / hoursPerYear
    cslope2Hr <- cslope2 / hoursPerYear

    # -----------------------------------------------------------------------
    # Individual finerenone drug-effect slope.
    # Control stream $PK:
    #   CV14   = (AGE/63)**THETA(18); ERACA2 = THETA(19)*SRACA2
    #   ETATR4 = (EXP(ETA(4))**BXPAR4 - 1) / BXPAR4
    #   ESLOPE = THETA(13) * CV14 * (1 + ERACA2) * EXP(ETATR4)
    # -----------------------------------------------------------------------
    etatrEslope <- (exp(etaeslope)^boxcox_eslope - 1) / boxcox_eslope
    eslopeInd   <- eslope *
      (AGE / 63)^e_age_eslope *
      (1 + e_race_japanese_eslope * RACE_JAPANESE) *
      exp(etatrEslope)

    ke0 <- exp(lke0)

    # -----------------------------------------------------------------------
    # Exposure metric: steady-state daily AUC, AUCtau,md = F * DOSE / CL, with
    # F implicit in the apparent clearance. podo(depot) returns the most recent
    # dose amount entered into the depot compartment, so the metric is a step
    # function that updates at every titration, interruption and re-initiation
    # event. podo() returns NA until the first dose record, so it is floored to
    # 0 -- which is also what the control stream's IF(TAFD.GT.0) guard does and
    # is what makes the placebo arm (no dose records at all) solvable.
    # -----------------------------------------------------------------------
    aucSs <- podo(depot) / cl
    if (is.na(aucSs)) aucSs <- 0

    # -----------------------------------------------------------------------
    # Drug effect on log UACR (control stream $DES / $ERROR).
    # -----------------------------------------------------------------------
    drugEff <- eslopeInd * effect^epow * exp((uacr - 850) * eslope2)
    if (effect <= 0) drugEff <- 0

    # -----------------------------------------------------------------------
    # ODE system.
    # Control stream $DES:
    #   COR  = CSLOPE  * A(1) * EXP(EFFECTD)
    #   COR2 = CSLOPE2 * LOG(A(1)*EXP(EFFECTD)/BSL)
    #   DADT(1) = PROG*A(1) - COR*A(1) - COR2*A(1)
    #   DADT(2) = KE0 * (AUC2 - A(2))
    # depot is a virtual dose receiver whose state is never read: only the most
    # recent dose amount matters, via podo(depot).
    # -----------------------------------------------------------------------
    uacrDrug <- uacr * exp(drugEff)

    d/dt(depot)  <- 0
    d/dt(uacr)   <- progRate * uacr -
      cslope1Hr * uacrDrug * uacr -
      cslope2Hr * log(uacrDrug / rbase) * uacr
    d/dt(effect) <- ke0 * (aucSs - effect)
    uacr(0)      <- rbase

    # -----------------------------------------------------------------------
    # Observation. The SGLT2 inhibitor effect acts on the observation, not on
    # the differential equation, and is applied to the change from the
    # treatment-start SGLT2i status (control stream $ERROR:
    #   IPRED = LOG(A(1)) + EFFECT + SGLTEFF*(FLAGSGLT - SGLTSTART)).
    # -----------------------------------------------------------------------
    uacrObs <- uacrDrug *
      exp(e_conmed_sglt2i_uacr * (CONMED_SGLT2I - CONMED_SGLT2I_BASE))

    uacrObs ~ lnorm(expSd)
  })
}
