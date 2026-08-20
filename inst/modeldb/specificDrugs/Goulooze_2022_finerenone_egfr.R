Goulooze_2022_finerenone_egfr <- function() {
  description <- "Population PKPD disease-progression model for the estimated glomerular filtration rate (eGFR) response to finerenone in patients with chronic kidney disease and type 2 diabetes (FIDELIO-DKD Phase III). eGFR declines at a constant chronic slope that is flattened by an exponential stabilisation function as eGFR approaches 16.1 mL/min/1.73 m^2, and the chronic slope is driven by the model-predicted UACR time course from the companion UACR model rather than by baseline UACR. Finerenone acts twice: an acute fully reversible eGFR decline through a power function of steady-state daily AUC acting via its own effect compartment, and, via its UACR reduction, a sustained flattening of the chronic slope. Finerenone PK is upstream (van den Berg 2022) and reduced here to AUCss = DOSE / CL with typical apparent clearance 29.9 L/h."
  reference <- paste(
    "Goulooze SC, Heerspink HJL, van Noort M, Snelder N, Brinker M, Lippert J,",
    "Eissing T. Dose-Exposure-Response Analysis of the Nonsteroidal",
    "Mineralocorticoid Receptor Antagonist Finerenone on UACR and eGFR:",
    "An Analysis from FIDELIO-DKD. Clin Pharmacokinet. 2022;61(7):1023-1037.",
    "doi:10.1007/s40262-022-01124-3.",
    "Parameter values are the final estimates in Tables S2 (eGFR) and S1",
    "(embedded UACR sub-model) of the Electronic Supplementary Material; the",
    "model equations are the final eGFR NONMEM control stream in the same",
    "supplement. The exposure metric AUCtau,md is computed from the upstream",
    "FIDELIO-DKD population PK analysis (van den Berg P et al.,",
    "Clin Pharmacokinet. 2022;61(7):1005-1021; doi:10.1007/s40262-021-01082-2);",
    "see modellib('vandenBerg_2021_finerenone'). The UACR sub-model is the",
    "companion model in the same paper; see",
    "modellib('Goulooze_2022_finerenone_uacr').",
    sep = " "
  )
  vignette <- "Goulooze_2022_finerenone_uacr_egfr"
  units <- list(
    time          = "h",
    dosing        = "(oral finerenone, mg)",
    concentration = "mL/min/1.73 m^2 (eGFR; the PD output is an endogenous renal-function measure rather than the dosed finerenone, so the dosing-vs-concentration dimensional check is not applicable and the dosing string is parenthesised to skip it)"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    depot   = list(analyte = "finerenone", units = "mg", specimen = "administration site", verified = TRUE),
    egfr    = list(analyte = "glomerular filtration rate", units = "mL/min/1.73 m^2", specimen = "not applicable", verified = TRUE),
    effect1 = list(analyte = "finerenone", units = "mg*h/L", specimen = "not applicable", verified = TRUE),
    uacr    = list(analyte = "albumin/creatinine ratio", units = "mg/g", specimen = "urine", verified = TRUE),
    effect2 = list(analyte = "finerenone", units = "mg*h/L", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    UACR = list(
      description        = "Observed baseline urine albumin-to-creatinine ratio",
      units              = "mg/g",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed (baseline). Source column UACR0. Enters twice: as a power scaling (UACR / 850)^0.877 on the baseline of the embedded UACR sub-model state, and as a power scaling (UACR / 850)^0.306 on the magnitude of the inter-individual variability of the eGFR decline rate. Note the second use is a covariate on a VARIANCE term, not on a typical value.",
      source_name        = "UACR0"
    ),
    CRCL = list(
      description        = "Observed baseline CKD-EPI estimated glomerular filtration rate, BSA-normalised to 1.73 m^2",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed (baseline). Source column EGFREPI0. Enters four times: power scalings (CRCL / 43)^0.882 on the baseline eGFR state, (CRCL / 43)^0.387 on the chronic decline rate, and (CRCL / 43)^-0.560 on the interaction between model-predicted UACR and the decline rate; plus (CRCL / 43)^-0.124 on the baseline of the embedded UACR sub-model. The reference 43 mL/min/1.73 m^2 is the FIDELIO-DKD cohort median.",
      source_name        = "EGFREPI0"
    ),
    POT = list(
      description        = "Observed baseline serum potassium concentration",
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed (baseline). Source column K0. Enters as a power scaling (POT / 4.4)^-1.38 on DSLOPE, the slope of the acute finerenone effect on eGFR; 4.4 mmol/L is the FIDELIO-DKD cohort median. Subjects with lower baseline potassium have a larger acute eGFR decline.",
      source_name        = "K0"
    ),
    HEPIMP_MOD = list(
      description        = "Moderate hepatic impairment indicator (1 = likely or certain Child-Pugh B, 0 = likely Child-Pugh A or healthy)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (likely Child-Pugh A or healthy)",
      notes              = "Classification scheme is Child-Pugh Class B, derived in FIDELIO-DKD from laboratory surrogates rather than a full Child-Pugh score: subjects with total bilirubin < 2 mg/dL AND serum albumin > 3.5 g/dL were categorised as likely Child-Pugh A or healthy, all others as likely or certain Child-Pugh B (Goulooze 2022 Table 1 footnote b). Source column CHILDPSC, with values 2 or 3 mapping to HEPIMP_MOD = 1 and values 1 or 5 to 0. Enters as a power-form multiplier 1.18^HEPIMP_MOD on the chronic eGFR decline rate, and as a proportional shift (1 + 0.0943 * HEPIMP_MOD) on the baseline of the embedded UACR sub-model.",
      source_name        = "(CHILDPSC %in% c(2, 3))"
    ),
    RACE_BLACK = list(
      description        = "Black or African American race indicator (1 = Black or African American, 0 = otherwise)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not Black or African American)",
      notes              = "Source column RACEASIA = 2 codes Black or African American. Enters as a power-form multiplier 1.24^RACE_BLACK on the chronic eGFR decline rate.",
      source_name        = "(RACEASIA == 2)"
    ),
    RACE_ASIAN = list(
      description        = "Asian-race indicator (1 = any Asian race category, 0 = otherwise)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Asian)",
      notes              = "Source column RACEASIA, a decimal-coded race variable; values in [3.0, 3.9] denote Asian race categories. Used only by the embedded UACR sub-model, as an additive shift +0.0634 /year on the UACR progression rate; it has no direct effect on eGFR.",
      source_name        = "(RACEASIA >= 3 & RACEASIA <= 3.9)"
    ),
    RACE_JAPANESE = list(
      description        = "Japanese-heritage race indicator (1 = Japanese, 0 = otherwise)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Japanese)",
      notes              = "Source column RACEASIA = 3.2 codes Japanese heritage. Used only by the embedded UACR sub-model, as a proportional shift (1 - 0.261 * RACE_JAPANESE) on the finerenone UACR drug-effect slope; that reduced UACR effect then propagates to the chronic eGFR slope. Japanese subjects are a subset of RACE_ASIAN, so both indicators are 1 for a Japanese subject.",
      source_name        = "(RACEASIA == 3.2)"
    ),
    AGE = list(
      description        = "Subject age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed (baseline). Used only by the embedded UACR sub-model, as a power scaling (AGE / 63)^0.864 on the finerenone UACR drug-effect slope; that effect then propagates to the chronic eGFR slope.",
      source_name        = "AGE"
    ),
    CONMED_SGLT2I = list(
      description        = "Current concomitant SGLT2 inhibitor use indicator (1 = in use, 0 = not in use)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant SGLT2 inhibitor)",
      notes              = "TIME-VARYING within subject. FIDELIO-DKD recorded SGLT2 inhibitor use as a binary variable over time defined as use / no use in the last 5 days, without regard to the specific agent or its dosing schedule (Goulooze 2022 Sect. 2.2.3). Enters three times: a proportional effect (1 - 0.457 * CONMED_SGLT2I) on the chronic eGFR decline rate, an acute direct effect (1 - 0.0558 * (CONMED_SGLT2I - CONMED_SGLT2I_BASE)) on the observed eGFR, and a proportional effect on the UACR driver of the chronic slope. The dynamics of the onset of the SGLT2 inhibitor effect could not be estimated from these data, so all three are instantaneous.",
      source_name        = "FLAGSGLT"
    ),
    CONMED_SGLT2I_BASE = list(
      description        = "Concomitant SGLT2 inhibitor use at treatment start (1 = in use at randomisation, 0 = not in use)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no SGLT2 inhibitor at treatment start)",
      notes              = "Time-fixed per subject. The control stream applies the acute SGLT2 inhibitor effects to the CHANGE from the treatment-start value, as theta * (FLAGSGLT - SGLTSTART), because the estimated baseline eGFR and baseline UACR of a subject already on an SGLT2 inhibitor at randomisation already reflect that drug's effect. Note the chronic-slope effect (1 - 0.457 * CONMED_SGLT2I) is NOT centred this way in the source and is applied to the raw current-use flag. Set to 0 for a subject who is SGLT2i-naive at randomisation.",
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
    notes          = "78,132 eGFR observations from 5,674 patients (full analysis set). 528 patients used an SGLT2 inhibitor at some point during the treatment period. Baseline demographics from Goulooze 2022 Table 1; the four columns of that table are placebo without SGLT2i (N = 2553), placebo with SGLT2i (N = 288), finerenone without SGLT2i (N = 2593) and finerenone with SGLT2i (N = 240), and the ranges quoted here span those four groups. The final eGFR model was estimated with SAEM followed by an importance-sampling expectation-only step, rather than FOCE, to avoid rounding errors."
  )

  ini({
    # -------------------------------------------------------------------------
    # Upstream PK (van den Berg 2022 popPK fit on FIDELIO-DKD; reference [10])
    # This paper fits NO PK. Its control stream computes the exposure metric
    # algebraically as AUCT = IF1 * DOSE / ICL, where IF1 and ICL are individual
    # post-hoc bioavailability and clearance supplied as DATA COLUMNS from the
    # upstream fit. At the upstream reference covariates F1 = 1 and CL/F = 29.9
    # L/h (van den Berg 2022 Table 2), so AUCtau,md = DOSE / 29.9 mg*h/L.
    # -------------------------------------------------------------------------
    lcl <- fixed(log(29.9)); label("Typical finerenone apparent clearance CL/F, upstream value (L/h)")  # van den Berg 2022 Table 2: CL/F = 29.9 L/h (RSE 3.62%); see modellib('vandenBerg_2021_finerenone')

    # -------------------------------------------------------------------------
    # eGFR disease progression - typical values
    # All values are the FINAL estimates from Goulooze 2022 ESM Table S2. The
    # $THETA block of the supplement's control stream holds INITIAL estimates
    # (e.g. TH1 initial 3.79 vs final 3.77; TH2 initial 3.1 vs final 3.01) and
    # is NOT used for values, only for the equation forms.
    # -------------------------------------------------------------------------
    lrbase_egfr <- log(exp(3.77));       label("Typical baseline eGFR (mL/min/1.73 m^2)")                              # ESM Table S2: theta_pop,BSLEGFR = ln[BSL_eGFR] = 3.77 (RSE 0.0305%), i.e. 43.4 mL/min/1.73 m^2
    pdecline    <- 3.01;                 label("PMAX, chronic linear eGFR decline rate (mL/min/1.73 m^2 per year)")   # ESM Table S2: theta_pop,PMAX = 3.01 (RSE 1.82%)
    lpint       <- fixed(2.78);          label("PINT, log eGFR stabilisation point at which decline stops (log mL/min/1.73 m^2)")  # ESM Table S2: theta_pop,PINT = ln[PINT] = 2.78, fixed to the placebo-model estimate; exp(2.78) = 16.1, matching the 16.2 mL/min/1.73 m^2 quoted in Sect. 3.3
    pslope      <- fixed(-0.300);        label("Slope of the eGFR decline stabilisation function (1.73 m^2 min/mL)")  # ESM Table S2: theta_pop,PSLOPE = -0.300, fixed to the placebo-model estimate

    # -------------------------------------------------------------------------
    # Acute, fully reversible finerenone effect on eGFR (haemodynamic), acting
    # through its own effect compartment. Control stream $DES / $ERROR:
    #   ADEFFECT = A(2)**DPOW * DSLOPE ;  IPRED = A(1)*(1 - ADEFFECT)*...
    # -------------------------------------------------------------------------
    ldslope <- log(0.0475);   label("DSLOPE, slope of the acute finerenone effect on eGFR (L/mg/h)")            # ESM Table S2: theta_pop,DSLOPE = 0.0475 (RSE 5.24%)
    dpow    <- 0.555;         label("DPOW, power of the finerenone exposure in the acute eGFR effect (unitless)")  # ESM Table S2: theta_pop,DPOW = 0.555 (RSE 7.05%)
    lke0    <- log(0.00230);  label("Equilibration rate of the acute finerenone eGFR effect compartment (1/h)")  # ESM Table S2: theta_pop,KE0,eGFR = 0.00230 /h (RSE 14.3%)

    # -------------------------------------------------------------------------
    # UACR-mediated effect on the chronic eGFR slope. This is the paper's
    # central finding: the model-predicted UACR time course, not baseline UACR,
    # drives the chronic decline, so the finerenone (and SGLT2i) effect on UACR
    # propagates to eGFR with no additional UACR-independent drug effect
    # needed (Sect. 3.4).
    # -------------------------------------------------------------------------
    e_uacr_pdecline      <- 0.000859; label("Effect of model-predicted UACR on the eGFR decline rate (mL/min/1.73 m^2 per year per mg/g)")  # ESM Table S2: theta_UACR,PMAX = 8.59e-4 (RSE 2.40%)
    e_crcl_uacr_pdecline <- -0.560;   label("Power exponent of (CRCL / 43) on the UACR-to-eGFR-decline interaction")                        # ESM Table S2: theta_EGFREPI0,INTER = -0.560 (RSE 14.7%)

    # -------------------------------------------------------------------------
    # Covariate effects on the eGFR model
    # -------------------------------------------------------------------------
    e_crcl_rbase_egfr        <- 0.882;   label("Power exponent of (CRCL / 43) on baseline eGFR")                       # ESM Table S2: theta_EGFREPI0,BSLEGFR = 0.882 (RSE 0.433%)
    e_crcl_pdecline          <- 0.387;   label("Power exponent of (CRCL / 43) on the eGFR decline rate")               # ESM Table S2: theta_EGFREPI0,PMAX = 0.387 (RSE 14.5%)
    e_race_black_pdecline    <- 1.24;    label("Power-form multiplier of Black or African American race on the eGFR decline rate")  # ESM Table S2: theta_BAA,PMAX = 1.24 (RSE 3.51%); applied as 1.24^RACE_BLACK
    e_hepimp_mod_pdecline    <- 1.18;    label("Power-form multiplier of Child-Pugh B on the eGFR decline rate")       # ESM Table S2: theta_CHILDPUGH,PMAX = 1.18 (RSE 3.71%); applied as 1.18^HEPIMP_MOD
    e_pot_dslope             <- -1.38;   label("Power exponent of (POT / 4.4) on the acute finerenone eGFR effect slope")  # ESM Table S2: theta_K0,DSLOPE = -1.38 (RSE 22.6%)
    e_conmed_sglt2i_acute    <- -0.0558; label("Acute proportional effect of current SGLT2 inhibitor use on observed eGFR")  # ESM Table S2: theta_SGLT2i,acute,eGFR = -0.0558 (RSE 8.73%)
    e_conmed_sglt2i_pdecline <- -0.457;  label("Proportional effect of current SGLT2 inhibitor use on the eGFR decline rate")  # ESM Table S2: theta_SGLT2i,PMAX = -0.457 (RSE 11.4%)
    e_uacr_iiv_pdecline      <- 0.306;   label("Power exponent of (UACR / 850) on the IIV magnitude of the eGFR decline rate")  # ESM Table S2: theta_UACR0,IIV,PMAX = 0.306 (RSE 6.13%)

    # -------------------------------------------------------------------------
    # Embedded UACR sub-model.
    # In the source, the eGFR run reads the UACR sub-model parameters
    # (BSLUACR, PROGUACR, CSLOPE, CSLOPE2, ESLOPE, EPOW, ESLOPE2, KE0UACR,
    # SGLTUACR) from DATA COLUMNS holding the individual post-hoc estimates of
    # the previously fitted UACR model - they are not re-estimated here. They
    # are therefore reproduced as fixed() parameters from ESM Table S1, with
    # their covariate model, so that this file is self-contained and reproduces
    # the paper's typical-subject simulations (Tables 3 and 4, Figs 3B and 6B).
    # Fitting this model to data would use individual post-hoc UACR parameters
    # instead; see modellib('Goulooze_2022_finerenone_uacr').
    # -------------------------------------------------------------------------
    lrbase_uacr             <- fixed(log(866));      label("Typical baseline UACR of the embedded sub-model (mg/g)")            # ESM Table S1: theta_pop,BSLUACR = 866 mg/g (RSE 0.429%)
    prog_uacr               <- fixed(0.137);         label("Typical UACR progression rate of the embedded sub-model (1/year)")  # ESM Table S1: theta_pop,PROG = 0.137 /year (RSE 9.10%)
    cslope1_uacr            <- fixed(0.0000720);     label("Correction of UACR progression by model-predicted UACR (g/mg/year)")  # ESM Table S1: theta_pop,CSLOPE1 = 7.20e-5 (RSE 8.39%)
    cslope2_uacr            <- fixed(0.182);         label("Correction of UACR progression by log UACR over baseline (1/year)")   # ESM Table S1: theta_pop,CSLOPE2 = 0.182 /year (RSE 21.2%)
    eslope_uacr             <- fixed(-1.42);         label("Slope of the finerenone effect on log UACR (L/mg/h)")                 # ESM Table S1: theta_pop,ESLOPE = -1.42 (RSE 6.61%)
    epow_uacr               <- fixed(0.613);         label("Power of the finerenone exposure in the UACR effect (unitless)")      # ESM Table S1: theta_pop,EPOW = 0.613 (RSE 3.99%)
    eslope2_uacr            <- fixed(-0.000166);     label("Attenuation of the finerenone UACR effect by current UACR (g/mg)")    # ESM Table S1: theta_pop,ESLOPE2 = -1.66e-4 (RSE 34.5%)
    lke0_uacr               <- fixed(log(0.000136)); label("Equilibration rate of the finerenone UACR effect compartment (1/h)")  # ESM Table S1: theta_pop,KE0,UACR = 1.36e-4 /h (RSE 16.1%)
    e_conmed_sglt2i_uacr    <- fixed(-0.212);        label("Effect of current SGLT2 inhibitor use on log UACR")                   # ESM Table S1: theta_SGLT2i,UACR = -0.212 (RSE 13.7%)
    e_uacr_rbase_uacr       <- fixed(0.877);         label("Power exponent of (UACR / 850) on baseline UACR")                     # ESM Table S1: theta_UACR0,BSLUACR = 0.877 (RSE 0.546%)
    e_crcl_rbase_uacr       <- fixed(-0.124);        label("Power exponent of (CRCL / 43) on baseline UACR")                      # ESM Table S1: theta_EGFREPI0,BSLUACR = -0.124 (RSE 11.4%)
    e_hepimp_mod_rbase_uacr <- fixed(0.0943);        label("Proportional effect of Child-Pugh B on baseline UACR")                # ESM Table S1: theta_CHILDPUGH,BSLUACR = 0.0943 (RSE 17.4%)
    e_crcl_prog_uacr        <- fixed(-0.00257);      label("Linear effect of (CRCL - 43) on the UACR progression rate (1/year per mL/min/1.73 m^2)")  # ESM Table S1: theta_EGFREPI0,PROG = -0.00257 (RSE 23.0%)
    e_race_asian_prog_uacr  <- fixed(0.0634);        label("Additive effect of Asian race on the UACR progression rate (1/year)") # ESM Table S1: theta_ASIAN,PROG = 0.0634 /year (RSE 25.8%)
    e_age_eslope_uacr       <- fixed(0.864);         label("Power exponent of (AGE / 63) on the UACR drug-effect slope")          # ESM Table S1: theta_AGE,ESLOPE = 0.864 (RSE 19.4%)
    e_race_japanese_eslope_uacr <- fixed(-0.261);    label("Proportional effect of Japanese race on the UACR drug-effect slope")  # ESM Table S1: theta_JAPAN,ESLOPE = -0.261 (RSE 21.7%)

    # -------------------------------------------------------------------------
    # Inter-individual variability
    # The source applies Box-Cox-transformed exponential IIV to baseline eGFR
    # and to the decline rate (control stream $PK:
    #   ETATR = (EXP(ETA)**BXPAR - 1) / BXPAR),
    # plain exponential IIV to DSLOPE, and plain exponential IIV to PINT with
    # its variance fixed to the placebo-model estimate.
    #
    # The published OMEGA also carries a variance on the residual magnitude
    # (omega^2 ERROR = 0.0898, RSE 3.63%, with BOXCOX ERROR = 1.06, RSE 8.47%),
    # which is omitted because nlmixr2 cannot carry IIV on the residual
    # magnitude - see the residual-error note below.
    #
    # Note the IIV on the decline rate is ADDITIVE on the annual decline rate
    # and is scaled by (UACR / 850)^0.306, so its variance of 8.62 is in units
    # of (mL/min/1.73 m^2 per year)^2 and is large by construction.
    # -------------------------------------------------------------------------
    boxcox_rbase_egfr <- 7.63;  label("Box-Cox shape parameter for the IIV on baseline eGFR")         # ESM Table S2: BOXCOX BSL = 7.63 (RSE 4.12%)
    boxcox_pdecline   <- 0.145; label("Box-Cox shape parameter for the IIV on the eGFR decline rate") # ESM Table S2: BOXCOX PROG = 0.145 (RSE 6.20%)

    etalrbase_egfr + etaldslope ~ c(0.00373,
                                    -0.0453, 0.833)  # ESM Table S2: omega^2 BSL_eGFR = 0.00373 (RSE 3.37%); omega^2 DSLOPE = 0.833 (RSE 8.06%); cov BSL/DSLOPE = -0.0453 (RSE 7.38%)
    etapdecline ~ 8.62                               # ESM Table S2: omega^2 PROG (eGFR decline rate) = 8.62 (RSE 3.57%)
    etalpint    ~ fixed(0.103)                       # ESM Table S2: omega^2 PINT = 0.103, carried over from the placebo model

    # -------------------------------------------------------------------------
    # Residual error
    # Control stream $ERROR:
    #   Y = IPRED*(1 + ERR(1)*EXP(ETATR1)) + ERR(2)*EXP(ETATR1)
    # i.e. a combined proportional-plus-additive residual whose magnitude
    # carries a Box-Cox-transformed IIV. The EXP(ETATR1) multiplier is NOT
    # reproducible in nlmixr2, which requires the residual SD to be a bare
    # estimated parameter rather than a model expression. It equals 1 at the
    # typical value of the eta, so typical-value simulation is unaffected; only
    # the between-subject spread of the residual is lost. See the vignette's
    # Assumptions and deviations section.
    # -------------------------------------------------------------------------
    propSd <- sqrt(0.0101); label("Proportional residual SD for eGFR (fraction)")                # ESM Table S2: sigma^2 PROP = 0.0101 (RSE 1.69%)
    addSd  <- sqrt(2.01);   label("Additive residual SD for eGFR (mL/min/1.73 m^2)")             # ESM Table S2: sigma^2 ADD = 2.01 (RSE 7.71%)
  })

  model({
    # -----------------------------------------------------------------------
    # Time-unit conversions, reproduced exactly as published: the eGFR model
    # converts annual rates with 365.25 days per year, while the embedded UACR
    # sub-model (carried over from the UACR run) uses 365.
    # -----------------------------------------------------------------------
    hoursPerYearEgfr <- 365.25 * 24
    hoursPerYearUacr <- 365 * 24

    # Finerenone apparent clearance (typical, upstream-fixed at 29.9 L/h)
    cl <- exp(lcl)

    # -----------------------------------------------------------------------
    # Individual baseline eGFR.
    # Control stream $PK:
    #   TVBSL  = EXP(THETA(1)) * (EGFREPI0/43)**THETA(6)
    #   ETATR  = (EXP(ETA(2))**BXPAR2 - 1) / BXPAR2
    #   BSL    = TVBSL * EXP(ETATR)
    # -----------------------------------------------------------------------
    etatrRbaseEgfr <- (exp(etalrbase_egfr)^boxcox_rbase_egfr - 1) / boxcox_rbase_egfr
    rbase_egfr     <- exp(lrbase_egfr) *
      (CRCL / 43)^e_crcl_rbase_egfr *
      exp(etatrRbaseEgfr)

    # -----------------------------------------------------------------------
    # Chronic eGFR decline rate, per hour.
    # Control stream $PK:
    #   TVPMAX     = THETA(2) * (EGFREPI0/43)**THETA(5) * THETA(8)**SCPGH * THETA(7)**SRACA1
    #   SGLTonPMAX = 1 + FLAGSGLT*THETA(20)
    #   PMAX       = TVPMAX * SGLTonPMAX / (365.25*24)
    #   ETATR4     = (EXP(ETA(4))**BXPAR4 - 1) / BXPAR4
    #   ETA4PMAX   = (UACR0/850)**THETA(15) * ETATR4 / (365.25*24)
    # The SGLT2i effect on the chronic slope uses the raw current-use flag; it
    # is not centred on the treatment-start value (unlike the acute effects).
    # -----------------------------------------------------------------------
    pdeclineHr <- pdecline *
      (CRCL / 43)^e_crcl_pdecline *
      e_hepimp_mod_pdecline^HEPIMP_MOD *
      e_race_black_pdecline^RACE_BLACK *
      (1 + e_conmed_sglt2i_pdecline * CONMED_SGLT2I) / hoursPerYearEgfr

    etatrPdecline <- (exp(etapdecline)^boxcox_pdecline - 1) / boxcox_pdecline
    etaPdeclineHr <- (UACR / 850)^e_uacr_iiv_pdecline * etatrPdecline / hoursPerYearEgfr

    # -----------------------------------------------------------------------
    # Stabilisation of the decline as eGFR approaches PINT.
    # Control stream $PK:
    #   PINT = EXP(MU_5 + ETA(5)) ; SLOPE = THETA(4)
    #   PA20 = PMAX / EXP((PINT - 20)*SLOPE)   ; floored at 0
    # PA20 is the offset that makes the eGFR derivative exactly 0 at eGFR =
    # PINT for a typical individual; the anchor at 20 mL/min/1.73 m^2 is a
    # numerical convenience in the source and cancels out.
    # -----------------------------------------------------------------------
    pint <- exp(lpint + etalpint)
    pa20 <- pdeclineHr / exp((pint - 20) * pslope)
    if (pa20 < 0) pa20 <- 0

    # -----------------------------------------------------------------------
    # Acute finerenone effect on eGFR.
    # Control stream $PK: DSLOPE = EXP(LOG(THETA(10)) + ETA(3)) * (K0/4.4)**THETA(12)
    # -----------------------------------------------------------------------
    dslope <- exp(ldslope + etaldslope) * (POT / 4.4)^e_pot_dslope
    ke0    <- exp(lke0)

    # -----------------------------------------------------------------------
    # Embedded UACR sub-model parameters (fixed from the UACR fit; see ini()).
    # Identical equation forms to Goulooze_2022_finerenone_uacr.
    # -----------------------------------------------------------------------
    rbase_uacr <- exp(lrbase_uacr) *
      (UACR / 850)^e_uacr_rbase_uacr *
      (CRCL / 43)^e_crcl_rbase_uacr *
      (1 + e_hepimp_mod_rbase_uacr * HEPIMP_MOD)

    progRateUacr <- (prog_uacr +
                       e_crcl_prog_uacr * (CRCL - 43) +
                       e_race_asian_prog_uacr * RACE_ASIAN) / hoursPerYearUacr

    cslope1HrUacr <- cslope1_uacr / hoursPerYearUacr
    cslope2HrUacr <- cslope2_uacr / hoursPerYearUacr

    eslopeIndUacr <- eslope_uacr *
      (AGE / 63)^e_age_eslope_uacr *
      (1 + e_race_japanese_eslope_uacr * RACE_JAPANESE)

    ke0Uacr <- exp(lke0_uacr)

    # -----------------------------------------------------------------------
    # Exposure metric: steady-state daily AUC, AUCtau,md = F * DOSE / CL, with
    # F implicit in the apparent clearance. Both effect compartments are driven
    # by the same exposure with different equilibration rates. podo() returns NA
    # until the first dose record, so it is floored to 0 -- which is also what
    # the control stream's IF(TAFD.GT.0) guard does and is what makes the
    # placebo arm (no dose records at all) solvable.
    # -----------------------------------------------------------------------
    aucSs <- podo(depot) / cl
    if (is.na(aucSs)) aucSs <- 0

    # Acute (haemodynamic) finerenone effect, from effect compartment 1
    acuteEff <- effect1^dpow * dslope
    if (effect1 <= 0) acuteEff <- 0

    # Finerenone effect on log UACR, from effect compartment 2
    uacrEff <- eslopeIndUacr * effect2^epow_uacr * exp((uacr - 850) * eslope2_uacr)
    if (effect2 <= 0) uacrEff <- 0

    # -----------------------------------------------------------------------
    # The UACR value that drives the chronic eGFR slope is the model-predicted
    # (i.e. drug- and SGLT2i-adjusted) UACR, not the raw progression state.
    # Control stream $DES:
    #   A(3)*EXP(EFFECTD + SGLTUACR*(FLAGSGLT - SGLTSTART))
    # -----------------------------------------------------------------------
    uacrDriver <- uacr *
      exp(uacrEff + e_conmed_sglt2i_uacr * (CONMED_SGLT2I - CONMED_SGLT2I_BASE))

    # -----------------------------------------------------------------------
    # ODE system.
    # Control stream $DES:
    #   DADT(1) = -(PMAX*(1 + (UACRDRIVER - 850)*THETA(13)*INT1) + ETA4PMAX
    #               - PA20*EXP((A(1)*(1-ADEFFECT2) - 20)*SLOPE))
    #   DADT(2) = KE0*AUCT - KE0*A(2)
    #   DADT(3) = PROGUACR*A(3) - COR*A(3) - COR2*A(3)
    #   DADT(4) = KE0UACR*(AUCT - A(4))
    # Note the stabilisation term reads the ACUTE-EFFECT-ADJUSTED eGFR, so a
    # patient whose acute decline brings them near PINT decelerates.
    # depot is a virtual dose receiver whose state is never read.
    # -----------------------------------------------------------------------
    int1     <- (CRCL / 43)^e_crcl_uacr_pdecline
    uacrDrug <- uacr * exp(uacrEff)

    d/dt(depot)   <- 0
    d/dt(egfr)    <- -(pdeclineHr * (1 + (uacrDriver - 850) * e_uacr_pdecline * int1) +
                         etaPdeclineHr -
                         pa20 * exp((egfr * (1 - acuteEff) - 20) * pslope))
    d/dt(effect1) <- ke0 * (aucSs - effect1)
    d/dt(uacr)    <- progRateUacr * uacr -
      cslope1HrUacr * uacrDrug * uacr -
      cslope2HrUacr * log(uacrDrug / rbase_uacr) * uacr
    d/dt(effect2) <- ke0Uacr * (aucSs - effect2)

    egfr(0) <- rbase_egfr
    uacr(0) <- rbase_uacr

    # -----------------------------------------------------------------------
    # Observation. The acute finerenone effect and the acute SGLT2i effect are
    # applied to the eGFR state at the observation, not inside the differential
    # equation, which is what makes them fully reversible on discontinuation
    # (control stream $ERROR:
    #   IPRED = A(1)*(1 - ADEFFECT)*(1 + SGLTEGFR*(FLAGSGLT - SGLTSTART))).
    # -----------------------------------------------------------------------
    egfrObs <- egfr * (1 - acuteEff) *
      (1 + e_conmed_sglt2i_acute * (CONMED_SGLT2I - CONMED_SGLT2I_BASE))

    egfrObs ~ prop(propSd) + add(addSd)
  })
}
