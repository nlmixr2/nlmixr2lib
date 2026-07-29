Willmann_2021_rivaroxaban <- function() {
  description <- paste(
    "Pediatric population PK model for rivaroxaban developed on the integrated",
    "EINSTEIN-Jr phase I / I-II / II / III dataset and interim PK from part A",
    "of the UNIVERSE study (524 children, 1988 plasma concentrations, age",
    "birth to <18 years, body weight 2.7-194 kg). Two-compartment disposition",
    "with first-order absorption and first-order elimination from the central",
    "compartment. Body weight enters as estimated allometric scaling on CL, Q,",
    "Vc, and Vp, centred on the 82.48 kg median of the integrated adult popPK",
    "analysis (a shared exponent is used for Vc and Vp). The undiluted",
    "ready-to-use oral suspension has a lower first-order absorption rate",
    "constant ka than the other three formulations (tablet, granules for oral",
    "suspension, and diluted ready-to-use oral suspension), which share a",
    "common ka. Relative oral bioavailability decreases with dose per body",
    "weight following an exponential function carried over from the integrated",
    "adult popPK analysis (anchored to F1 = 1 at 10 mg / 82.48 kg = 0.1213",
    "mg/kg). Inter-individual variability is on CL and F1 only (no IIV on Ka,",
    "Vc, Vp, or Q); residual error is proportional. Age, eGFR (Schwartz and",
    "Rhodin), serum creatinine, comedications (CYP3A4 inhibitors / inducers,",
    "P-gp inhibitors), and Fontan status were tested and not retained."
  )
  reference <- paste(
    "Willmann S, Coboeken K, Zhang Y, Mayer H, Ince I, Mesic E, Thelen K,",
    "Kubitza D, Lensing AWA, Yang H, Zhu P, Mueck W, Drenth HJ, Lippert J.",
    "Population pharmacokinetic analysis of rivaroxaban in children and",
    "comparison to prospective physiologically-based pharmacokinetic",
    "predictions. CPT Pharmacometrics Syst Pharmacol. 2021;10(10):1195-1207.",
    "doi:10.1002/psp4.12688"
  )
  vignette <- "Willmann_2021_rivaroxaban"
  units    <- list(
    time          = "hour",
    dosing        = "mg",
    concentration = "ug/mL"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying (paper Methods 'Covariate analysis' p. 1197 states",
        "bodyweight was treated as time-varying in children). Used for",
        "estimated allometric scaling centred on the 82.48 kg reference (the",
        "median bodyweight of the integrated adult popPK analysis in",
        "reference 19): cl_typ = exp(lcl) * (WT/82.48)^e_wt_cl;",
        "vc_typ = exp(lvc) * (WT/82.48)^e_wt_vc_vp;",
        "vp_typ = exp(lvp) * (WT/82.48)^e_wt_vc_vp;",
        "q_typ  = exp(lq)  * (WT/82.48)^e_wt_q. The volume exponent is",
        "shared across Vc and Vp (Willmann 2021 Table 2 footnote and",
        "control-stream THETA(5))."
      ),
      source_name        = "WGHT"
    ),
    FORM_UNDILUTED_SUSP = list(
      description        = "Undiluted ready-to-use oral suspension formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (tablet, granules for oral suspension, or diluted ready-to-use oral suspension; share the same ka)",
      notes              = paste(
        "1 = subject received the ready-to-use oral suspension administered",
        "undiluted (directly into the mouth); 0 = subject received any of the",
        "three rapidly-absorbed formulations (tablet, granules for oral",
        "suspension, or ready-to-use oral suspension diluted through mixing",
        "with a defined volume of non-sparkling liquid before administration).",
        "Per-dose-occasion indicator (a subject may switch formulations across",
        "study phases). Selects ka via:",
        "ka = exp(lka) for FORM_UNDILUTED_SUSP = 0 and",
        "ka = exp(lka_susp_undiluted) for FORM_UNDILUTED_SUSP = 1.",
        "Source: Willmann 2021 supplement S2 NONMEM control stream",
        "$PK block 'FD = 1 / IF(FORM.EQ.2.AND.DILU.EQ.1) FD = 2 / KA = THETA(1)",
        "/ IF(FD.EQ.2) KA = THETA(6)' (THETA(6) is the undiluted suspension",
        "ka); Willmann 2021 Results p. 1199 'A lower ka was estimated for the",
        "undiluted ready-to-use oral suspension (0.226 1/h) when compared to",
        "the other tested formulations (i.e., tablet, granules for oral",
        "suspension, and diluted ready-to-use oral suspension, ka = 0.799",
        "1/h)'. PROPOSED NEW CANONICAL pending PR review; sibling of the",
        "FORM_* drug-product / preparation-state family."
      ),
      source_name        = "FORM,DILU (composite: FD = 2 if FORM=2 and DILU=1)"
    ),
    DOSE_RIV_MGKG = list(
      description        = "Per-administration rivaroxaban dose per kg body weight",
      units              = "mg/kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Per-dose-record rivaroxaban dose given per administration in mg per",
        "kg of body weight (mg/kg). Used in the dose-dependent relative oral",
        "bioavailability formula carried over from the integrated adult popPK",
        "analysis (reference 19), where DW = DOSE/WGHT and",
        "tf1 = f1min + (f1max - f1min) * exp(-log(2)/d50 * DW),",
        "with f1max = 1.25, f1min = 0.59, and d50 = 14.4/82.48 mg/kg = 0.1746",
        "mg/kg (Willmann 2021 supplement S2 NONMEM $PK block",
        "'F1MAX = 1.25 / F1MIN = 0.59 / DW = DOSE/WGHT / D50 = 14.4/82.48 /",
        "TF1 = F1MIN + (F1MAX - F1MIN) * EXP(-LOG(2)/D50 * DW)'). The",
        "function is anchored to F1 = 1.0 at DW = 10 mg / 82.48 kg = 0.1213",
        "mg/kg (Willmann 2021 Results p. 1199). Distinct from the per-dose-",
        "event amt value because the bioavailability formula is decoupled",
        "from the rxode2 event-table structure -- the user supplies",
        "DOSE_RIV_MGKG as an explicit covariate column on every record so",
        "the bioavailability expression can be evaluated for any record",
        "without back-computing dose from the event table. PROPOSED NEW",
        "CANONICAL pending PR review; sibling of DOSE_RTV_MGKG (Zhang 2012)",
        "and DOSE_PHT_MGKGD (Yukawa 1990)."
      ),
      source_name        = "DW (= DOSE/WGHT in supplement S2 $PK block)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Tested as an additional covariate on CL and on F1 in the forward",
        "inclusion / backward deletion procedure (Willmann 2021 p. 1197).",
        "No effect of age on CL or F1 could be identified by the model",
        "(Willmann 2021 Results p. 1199 'The median age of the pooled",
        "pediatric population was 9.0 years ... No effect of age on CL or F1",
        "could be identified by the model'). Not retained in the final model."
      )
    ),
    EGFR_SCHWARTZ = list(
      description = "Estimated glomerular filtration rate by the Schwartz formula (height and serum creatinine)",
      units       = "mL/min/1.73 m^2",
      type        = "continuous",
      notes       = paste(
        "Tested as a covariate on CL (Willmann 2021 p. 1197 Methods 'Renal",
        "function'). Range in the pooled pediatric population 43.8-456",
        "mL/min/1.73 m^2 (median 150). 'None of the four methods to test",
        "for an influence of renal function led to a significant improvement",
        "of the objective function' (Willmann 2021 Results p. 1199). Not",
        "retained in the final model."
      )
    ),
    EGFR_RHODIN = list(
      description = "Estimated glomerular filtration rate by the Rhodin formula (body size + postmenstrual age, creatinine-independent)",
      units       = "mL/min/1.73 m^2",
      type        = "continuous",
      notes       = paste(
        "Tested as a covariate on CL (Willmann 2021 p. 1197 Methods 'Renal",
        "function'). Not retained in the final model."
      )
    ),
    SCR_ULN_RATIO = list(
      description = "Ratio of individual serum creatinine to the upper limit of normal",
      units       = "(unitless ratio)",
      type        = "continuous",
      notes       = paste(
        "Tested as a covariate on CL (Willmann 2021 p. 1197 Methods 'Renal",
        "function'). Not retained in the final model."
      )
    ),
    SCR_CAT = list(
      description = "Categorical score for serum creatinine (above vs at or below the upper limit of normal)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Tested as a covariate on CL (Willmann 2021 p. 1197 Methods 'Renal",
        "function'). Not retained in the final model."
      )
    ),
    CONMED_CYP3A4_INHIB_WEAK = list(
      description = "Weak CYP3A4 inhibitor co-medication indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Tested exploratively on CL and F1 (only 41 subjects / 7.8% of the",
        "dataset; Willmann 2021 Table 3). No significant effect identified",
        "(Willmann 2021 Results p. 1199). Not retained."
      )
    ),
    CONMED_CYP3A4_INHIB_MOD = list(
      description = "Moderate CYP3A4 inhibitor co-medication indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Tested exploratively on CL and F1 (only 18 subjects / 3.4%;",
        "Willmann 2021 Table 3). No significant effect identified. Not",
        "retained."
      )
    ),
    CONMED_CYP3A4_IND = list(
      description = "CYP3A4 inducer co-medication indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Tested exploratively on CL and F1 (16 subjects / 3.1%; Willmann",
        "2021 Table 3). No significant effect identified. Not retained."
      )
    ),
    FONTAN = list(
      description = "Post-Fontan-procedure patient indicator (UNIVERSE study) vs venous thromboembolism patient (EINSTEIN-Jr)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Tested on CL (Willmann 2021 p. 1199). Showed a small but",
        "statistically significant drop in OFV in univariate forward",
        "inclusion, but not significant in backward elimination under the",
        "more stringent criterion; not retained in the final model."
      )
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 524L,
    n_observations  = 1988L,
    n_studies       = 7L,
    age_range       = "birth to <18 years (median 9.0)",
    age_median      = "9.0 years",
    weight_range    = "2.7-194 kg (median 29.5)",
    weight_median   = "29.5 kg",
    sex_female_pct  = (109 + 47 + 53 + 30 + 8) / 524 * 100,
    race_ethnicity  = paste(
      "Multi-regional pooled cohort spanning EINSTEIN-Jr (phase I, I-II, II,",
      "III) and UNIVERSE phase III part A. Race / ethnicity subgroups",
      "explored graphically (Japanese, Chinese, Asian outside Japan and",
      "China) but not used as a model covariate; no clustering of any",
      "race/ethnicity subgroup at exposure extremes was observed (Willmann",
      "2021 Results p. 1199 and Supplementary Figures S7-S9)."
    ),
    disease_state   = paste(
      "Acute venous thromboembolism (EINSTEIN-Jr) and post-Fontan",
      "thromboprophylaxis (UNIVERSE part A). 512/524 (97.7%) participated",
      "in EINSTEIN-Jr; 12/524 (2.3%) in UNIVERSE."
    ),
    dose_range      = paste(
      "Bodyweight-adjusted, 0.4-20 mg per dose; absolute single or daily",
      "doses 0.4-20 mg. Bodyweight-normalised single doses ~0.1-0.5 mg/kg",
      "(Willmann 2021 Results p. 1199). Once-daily, twice-daily, or",
      "thrice-daily regimens depending on age group / weight band: once-daily",
      "for body weight >=30 kg, twice-daily for 12-<30 kg, thrice-daily for",
      "<12 kg (Willmann 2021 Discussion p. 1201)."
    ),
    regions         = "Multi-regional (EINSTEIN-Jr phase III recruited globally; supplementary subgroup analyses include Japanese, Chinese, and other Asian subjects).",
    age_group_n     = "12-<18 years: 193; 6-<12 years: 135; 2-<6 years: 110; 6 mo-<2 years: 63; birth-<6 mo: 23 (supplement Table S1)",
    n_subjects_under_2_years = 86L,
    formulations    = "Tablet, granules for oral suspension, ready-to-use oral suspension (administered diluted or undiluted). Tablets and granules for oral suspension were the phase III formulations.",
    bloq_handling   = "Not reported in the main paper; standard EINSTEIN-Jr / UNIVERSE PK assay LLOQ procedures apply per the integrated dataset.",
    sampling        = "Sparse PK sampling per study phase; up to 5 samples per subject in older adolescents and 2 samples in <2 yr olds during the single-dose phase I; sparse 1-2 samples per child per study day in the multiple-dose phases (supplement Table S1).",
    notes           = paste(
      "Sex split derived from supplement Table S1 per-age-group counts:",
      "247 female / 277 male (47.1% female). Demographics summarised from",
      "Willmann 2021 Table 1 and supplement Table S1."
    )
  )

  ini({
    # Structural PK parameters - reference subject body weight 82.48 kg
    # (median of the integrated adult popPK analysis in Willmann 2021
    # reference 19). Final estimates from Willmann 2021 Table 2 'Value' column.
    lka <- log(0.799)
    label("First-order absorption rate constant ka for tablets, granules for oral suspension, and diluted ready-to-use oral suspension (1/h)")  # Willmann 2021 Table 2: ka = 0.799 1/h (95% CI 0.655-0.944; %CV 9.21)
    lka_susp_undiluted <- log(0.226)
    label("First-order absorption rate constant ka for the undiluted ready-to-use oral suspension (1/h)")  # Willmann 2021 Table 2: ka = 0.226 1/h (95% CI 0.154-0.297; %CV 16.2)
    lcl <- log(8.02)
    label("Apparent oral clearance CL at 82.48 kg reference (L/h)")  # Willmann 2021 Table 2: CL = 8.02 L/h (95% CI 7.53-8.51; %CV 3.14)
    lvc <- log(53.2)
    label("Central volume of distribution Vc at 82.48 kg reference (L)")  # Willmann 2021 Table 2: Vc = 53.2 L (95% CI 47.2-59.3; %CV 5.77)
    lvp <- log(59.1)
    label("Peripheral volume of distribution Vp at 82.48 kg reference (L)")  # Willmann 2021 Table 2: Vp = 59.1 L (95% CI 29.1-89.1; %CV 25.9)
    lq  <- log(2.50)
    label("Intercompartmental clearance Q at 82.48 kg reference (L/h)")  # Willmann 2021 Table 2: Q = 2.50 L/h (95% CI 1.69-3.31; %CV 16.6)

    # Estimated allometric exponents on body weight (Willmann 2021 Methods
    # 'Covariate analysis' p. 1197 - exponents estimated, not fixed; Table 2
    # 'Value' column; Discussion p. 1201 - 'No efforts were undertaken to
    # fix the fitted allometric exponents to theoretical or published
    # values'). The exponent on Vc and Vp is shared (control stream
    # THETA(5) used for both V2 and V3).
    e_wt_cl    <- 0.481
    label("Estimated allometric WT exponent on CL (unitless)")  # Willmann 2021 Table 2: 0.481 (95% CI 0.434-0.527; %CV 4.96)
    e_wt_vc_vp <- 0.821
    label("Estimated allometric WT exponent shared across Vc and Vp (unitless)")  # Willmann 2021 Table 2: 0.821 (95% CI 0.760-0.881; %CV 3.75); shared exponent for Vc and Vp per supplement S2 control stream THETA(5)
    e_wt_q     <- 0.761
    label("Estimated allometric WT exponent on Q (unitless)")  # Willmann 2021 Table 2: 0.761 (95% CI 0.561-0.961; %CV 13.4)

    # Dose-dependent relative oral bioavailability function carried over
    # from the integrated adult popPK analysis (Willmann 2021 reference 19)
    # after replacing absolute dose by dose per bodyweight and re-anchoring
    # to F1 = 1.0 at 10 mg / 82.48 kg = 0.1213 mg/kg. The three constants
    # f1max, f1min, and d50 are hard-coded a-priori values in the source
    # NONMEM control stream ($PK block, supplement S2): not estimated by
    # the pediatric fit, hence fixed() here. The function reproduces
    # Willmann 2021 Results p. 1199 anchor points: F1(0.12) = 1.00,
    # F1(0.30) = 0.791, F1(0.50) = 0.681. See vignette source-trace table.
    f1max <- fixed(1.25)
    label("Upper asymptote of the dose-dependent relative bioavailability function (unitless)")  # Willmann 2021 supplement S2 $PK: F1MAX = 1.25 (a priori, not estimated)
    f1min <- fixed(0.59)
    label("Lower asymptote of the dose-dependent relative bioavailability function (unitless)")  # Willmann 2021 supplement S2 $PK: F1MIN = 0.59 (a priori, not estimated)
    d50   <- fixed(14.4 / 82.48)
    label("Dose-per-weight half-decline scaling parameter of the bioavailability function (mg/kg)")  # Willmann 2021 supplement S2 $PK: D50 = 14.4/82.48 mg/kg (a priori, not estimated)

    # Log-fixed bioavailability anchor (lfdepot = log(1) = 0). The
    # typical bioavailability function tf1 carries the dose-dependent
    # absolute value; lfdepot is held at the no-op anchor so the IIV
    # parameter etalfdepot pairs with a named fixed-effect by the
    # nlmixr2lib parameter-naming convention.
    lfdepot <- fixed(log(1))
    label("Bioavailability anchor (held at log(1) = 0; the dose-dependent function tf1 carries the actual F1)")  # No-op anchor for the etalfdepot pairing convention.

    # Inter-individual variability - Willmann 2021 Table 2 'Random effects:
    # Interindividual variability'. Exponential IIV (additive on log
    # scale): cl_i = cl_typ * exp(eta_cl), tf1_i = tf1_typ * exp(eta_f1).
    # Diagonal omega matrix in the source (no off-diagonals reported in
    # supplement S2 $OMEGA block); IIV is on CL and F1 only.
    etalcl     ~ 0.0705
    # Willmann 2021 Table 2: omega^2 on CL = 0.0705 (95% CI 0.0453-0.0957);
    # CV calculated as sqrt(exp(omega^2)-1)*100 = 27.0% (paper footnote g)
    etalfdepot ~ 0.0612
    # Willmann 2021 Table 2: omega^2 on F1 = 0.0612 (95% CI 0.0407-0.0818);
    # CV calculated as sqrt(exp(omega^2)-1)*100 = 25.1% (paper footnote g)

    # Residual variability - Willmann 2021 supplement S2 $ERROR block
    # encodes a proportional residual on the linear scale:
    # Y = IPRED + W*EPS(1) with W = F, equivalent to Y = F * (1 + EPS).
    # Table 2 reports sigma^2 = 0.220 (proportional variance);
    # propSd = sqrt(0.220) = 0.4690 (paper footnote h reports
    # sigvar = sqrt(sigma^2)*100 = 46.9%).
    propSd <- 0.4690
    label("Proportional residual SD for rivaroxaban plasma concentration (fraction)")  # Willmann 2021 Table 2: sigma^2 = 0.220; propSd = sqrt(0.220) = 0.4690
  })

  model({
    # Individual PK parameters with allometric scaling.
    # Reference subject: WT = 82.48 kg, FORM_UNDILUTED_SUSP = 0 (tablet,
    # granules for oral suspension, or diluted ready-to-use suspension).
    # ka is selected piecewise from the two structural Ka estimates by
    # the binary FORM_UNDILUTED_SUSP indicator (Willmann 2021 supplement
    # S2 $PK 'IF(FD.EQ.2) KA = THETA(6)').
    ka <- exp(lka) * (1 - FORM_UNDILUTED_SUSP) +
          exp(lka_susp_undiluted) * FORM_UNDILUTED_SUSP
    cl <- exp(lcl + etalcl) * (WT / 82.48)^e_wt_cl
    vc <- exp(lvc)          * (WT / 82.48)^e_wt_vc_vp
    vp <- exp(lvp)          * (WT / 82.48)^e_wt_vc_vp
    q  <- exp(lq)           * (WT / 82.48)^e_wt_q

    # Micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Dose-dependent relative oral bioavailability (Willmann 2021
    # supplement S2 $PK block). DOSE_RIV_MGKG is supplied as a per-record
    # covariate (per-administration mg/kg). Inter-individual variability
    # is applied multiplicatively to the typical function (NONMEM
    # F1 = TF1 * EXP(ETA(2))).
    tf1 <- f1min + (f1max - f1min) * exp(-log(2) / d50 * DOSE_RIV_MGKG)
    f1  <- tf1 * exp(lfdepot + etalfdepot)

    # Two-compartment ODE system: first-order oral absorption from depot
    # into central, first-order distribution central <-> peripheral1, and
    # first-order elimination from central.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Bioavailability factor on the depot dose.
    f(depot) <- f1

    # Rivaroxaban plasma concentration. Dose mg, volume L -> mg/L; multiply
    # by 1000 in the vignette to compare against the paper's reported ug/L
    # exposure metrics (1 mg/L = 1000 ug/L = 1 ug/mL).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
