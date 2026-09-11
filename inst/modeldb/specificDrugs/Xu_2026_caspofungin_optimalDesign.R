Xu_2026_caspofungin_optimalDesign <- function() {
  description <- "Two-compartment population PK model with first-order elimination for intravenous caspofungin in critically ill Chinese children, fitted to the n = 14 intensive-sampling subset of the Xu 2026 PICU study and used to drive the NONMEM $DESIGN / PopED optimisation of the sparse-sampling scheme for the study's second stage. Clearance and intercompartmental clearance scale with body weight through a FIXED exponent of 0.75, and both volumes through a FIXED exponent of 1, standardised to 70 kg. Aspartate aminotransferase carries an estimated power effect of 0.898 on intercompartmental clearance, centred on the subset median of 52.5 U/L -- an effect that did NOT survive into the paper's final 29-patient model. Inter-individual variability on clearance and central volume is strongly correlated (82.4%), with no random effect on Q or V2, and the residual error is purely additive. Companion to Xu_2026_caspofungin.R, which is the paper's final model. Xu 2026 supplemental Table S2 and control stream, n = 14 patients."
  reference <- "Xu N, Shi Y, Ju G, Liu X, Yan G, Zheng Y, Hou S, Xiang X, Lu G, Ouyang D, Zhu X, Wang Y. Population pharmacokinetics of caspofungin in critically ill Chinese children: a prospective observational study. Antimicrob Agents Chemother. 2026;70(2):e01277-25. doi:10.1128/aac.01277-25. PMC12888871. Parameters from supplemental material (aac.01277-25-s0001.docx) Table S2 and the 'Final PPK model for the design time points optimization in NONMEM' control stream. ClinicalTrials.gov NCT04961593."
  vignette <- "Xu_2026_caspofungin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Same assay and matrix as the paper's final model --
  # TOTAL plasma caspofungin by LC-MS/MS, LLOQ 0.05 ug/mL (Methods,
  # 'Caspofungin assay and fungal cultures'). The control stream's $ERROR
  # block confirms the central compartment carries the observed species:
  # IPRED = A(1) / V1.
  compartmentData <- list(
    central     = list(analyte = "caspofungin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "caspofungin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The body-size descriptor of THIS model, in contrast to the paper's",
        "final model (Xu_2026_caspofungin.R), which uses body surface area.",
        "Enters as (WT / 70)^0.75 on CL and Q and (WT / 70)^1 on V1 and V2.",
        "Both exponents were FIXED at the canonical allometric values, not",
        "estimated: Table S2 rows 'The effect of weight on CL / V1 / Q / V2'",
        "read 0.75 FIX, 1 FIX, 0.75 FIX and 1 FIX, with '/' in the RSE and",
        "SIR columns.",
        "The functional form is not inferred -- it is printed twice in the",
        "supplement. Table S2's footnote gives the equations directly",
        "('CL = 0.477 x (weight/70)^0.75; V1 = 11.8 x (weight/70)^1;",
        "Q = 0.512 x (weight//70)^0.75 x (AST/52.5)^0.898; V2 = 19.2 x",
        "(weight/70)^1'; the doubled slash in the Q line is a typographical",
        "slip in the source), and the $PK block of the accompanying control",
        "stream restates it as TVCL = THETA(1) * (WT/70)**0.75 and so on.",
        "The 70 kg standardisation is the conventional adult reference, NOT a",
        "cohort statistic: the n = 14 subset median weight is 15.9 kg (range",
        "4.90-64.0; Table S1), so every patient in the fit sits well below the",
        "reference and the typical values below are extrapolations to a 70 kg",
        "individual rather than descriptions of a study subject. Time-fixed.",
        "Must be strictly positive; it enters a power term."
      ),
      source_name        = "WT"
    ),
    AST = list(
      description        = "Serum aspartate aminotransferase activity",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The only clinical covariate in this model, on intercompartmental",
        "clearance only, as the power term (AST / 52.5)^0.898. Estimated (RSE",
        "16%, SIR median 0.875 [0.655-1.075]), so NOT wrapped in fixed().",
        "The 52.5 U/L divisor is the median AST of the n = 14 subset that",
        "this model was fitted to (Table S1: 52.5 [18.3, 631]); Figure S3's",
        "caption independently confirms it by naming 'AST = 18.3U/L, 52.5U/L",
        "and 631U/L' as the 5th, 50th and 95th percentiles of the included",
        "patients. The control stream hardcodes the same constant:",
        "TVQ = THETA(3) * (WT/70)**0.75 * (AST/52.5)**0.898.",
        "THIS EFFECT DID NOT SURVIVE into the paper's final 29-patient model.",
        "The Discussion of the main text states that 'Other factors,",
        "including AST and ALT levels, did not appear to influence the PK of",
        "caspofungin', and Table 2 of the main text carries no AST term. The",
        "effect is retained here because this model is extracted as the",
        "authors published it -- it is the model that generated the",
        "sparse-sampling design actually used in stage 2, and its sensitivity",
        "to AST was analysed explicitly (Figure S3B). Treat it as a",
        "subset-specific finding, not as a caspofungin class effect.",
        "Reported as U/L in Table S1 and IU/L in the Table S1 row label;",
        "the two are used interchangeably and the values are identical.",
        "Must be strictly positive; it enters a power term. Studied range",
        "18.3-631 U/L in this subset."
      ),
      source_name        = "AST"
    )
  )

  # Covariates screened in the main paper's stepwise covariate modelling but
  # absent from this model. Documentation only -- neither is referenced in
  # model(). BSA is listed because it is the descriptor that REPLACED weight
  # when the fit was extended to all 29 patients, which is the single most
  # important thing to know about this model's relationship to its companion.
  covariatesDataExcluded <- list(
    BSA = list(
      description = "Body surface area computed with the Mosteller formula",
      units       = "m^2",
      type        = "continuous",
      notes       = paste(
        "NOT used in this model, which scales on total body weight. BSA is",
        "the body-size descriptor of the paper's final model",
        "(Xu_2026_caspofungin.R), where it beat weight, lean body weight and",
        "fat-free mass on OFV and AIC (Table S6) and is standardised to 0.79",
        "m^2 with exponents 0.66 and 1. The n = 14 subset median BSA is 0.655",
        "m^2 (range 0.286-1.69; Table S1). A user who wants the BSA-driven",
        "model should load the companion rather than substituting BSA here,",
        "because the exponents and typical values differ."
      )
    ),
    ECMO_STATUS = list(
      description = "Extracorporeal membrane oxygenation support indicator (1 = receiving ECMO)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "NOT a covariate in this model, although 3 of the 14 patients",
        "(21.4%; Table S1) were on ECMO and the control stream's $INPUT",
        "carries no ECMO column at all. The 18.2-fold ECMO effect on the",
        "central volume appears only in the paper's final 29-patient model",
        "(Xu_2026_caspofungin.R). Any ECMO-driven volume expansion in this",
        "subset is therefore absorbed into the 50.3% inter-individual",
        "variability on V1 rather than described structurally."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 14L,
    n_studies      = 1L,
    age_range      = "0.330-16.0 years (Table S1)",
    age_median     = "4.63 years (Table S1)",
    weight_range   = "4.90-64.0 kg (Table S1)",
    weight_median  = "15.9 kg (Table S1). Far below the 70 kg allometric standardisation, which is the conventional adult reference rather than a cohort statistic.",
    height_median  = "97.5 cm (range 54.0-160; Table S1)",
    bsa_median     = "0.655 m^2 (range 0.286-1.69; Table S1)",
    bmi_median     = "18.01 kg/m^2 (range 9.82-25.64; Table S1)",
    sex_female_pct = 64.3,
    race_ethnicity = "Chinese. Single-centre enrolment at the Children's Hospital of Fudan University, Shanghai; the source reports no further race or ethnicity breakdown.",
    disease_state  = paste(
      "Critically ill children in the paediatric intensive care unit treated",
      "with caspofungin, sampled intensively in the first stage of the Xu",
      "2026 study. 3 of 14 (21.4%) were on ECMO (Table S1). The subset is",
      "more hepatically and renally deranged at the median than the full",
      "cohort: total bilirubin median 26.6 umol/L (range 2.10-357) against",
      "8.50 in the full 29, and direct bilirubin 13.1 umol/L (range",
      "1.10-228) against 4.10."
    ),
    hepatic_function = "AST median 52.5 U/L (range 18.3-631) -- the centring value of the Q covariate effect; ALT median 17.0 U/L (range 3.51-223); total bilirubin median 26.6 umol/L (range 2.10-357); direct bilirubin median 13.1 umol/L (range 1.10-228); albumin median 36.0 g/L (range 24.7-49.9); total protein median 61.6 g/L (range 36.8-75.4) (Table S1).",
    renal_function = "Serum creatinine median 27.7 umol/L (range 14.0-223); uric acid median 177 umol/L (range 87.0-736) (Table S1).",
    dose_range     = "Once-daily 1 h intravenous infusion on a BSA-based regimen: loading dose 70 mg/m^2 on day 1 and maintenance 50 mg/m^2 thereafter, each capped at 70 mg.",
    regions        = "China (single centre: Children's Hospital of Fudan University, National Children's Medical Center, Shanghai)",
    sampling       = "Intensive stage-1 sampling, scheduled pre-dose and 1, 2, 4, 8 (if feasible) and 16 h (if feasible) after the sixth dose -- i.e. nominal clock times of 120, 121, 122, 124, 128 and 136 h after the first dose, which is how the $DESIGN evaluation expresses them.",
    notes          = paste(
      "PURPOSE. This is not a base model discarded on the way to the final",
      "fit; it is a separately reported model with its own SIR-based",
      "uncertainty quantification whose job was to drive the optimal-design",
      "calculation for the study's second stage. It was fitted to a different",
      "(smaller) cohort than the final model, uses a different body-size",
      "descriptor, retains a covariate the final model rejects, and has a",
      "different residual-error structure. The supplement titles Table S2",
      "'Parameter estimates of the final caspofungin population",
      "pharmacokinetic model' and heads the control stream 'Final PPK model",
      "for the design time points optimization in NONMEM' -- final for the",
      "design exercise, superseded for the paper's PK conclusions. It is",
      "extracted as a separate file per the replicate-the-author's-structure",
      "policy (independent fits on different cohorts -> N files, 1 vignette).",
      "DESIGN OUTPUT. Driving $DESIGN (GROUPSIZE = 40, FIMDIAG = 2) with this",
      "model gave RSEs of 12% for CL, 6% for V1, 12% for Q and 43% for V2 on",
      "the intensive design (Table S5), closely matched by stochastic",
      "simulation and estimation (n = 1,000). At least four samples per",
      "patient were needed: three inflated the V2 RSE past 50% (Table S3,",
      "Optimisation 3 at 59%). The selected D-optimal times were 119, 121,",
      "126.5 and 144 h after the first dose, relaxed to the practical windows",
      "119-120, 120.5-121.5, 126-127 and 143-144 h, which preserved over 95%",
      "of the information (Figure S4A). Sensitivity analysis over the 5th to",
      "95th percentiles of weight and AST moved the third sampling point by",
      "under 3 h (Figure S3). Transferability to four other published",
      "paediatric caspofungin models kept the median typical-parameter RSE",
      "under 60% and bias within +/-20%, but IIV bias exceeded +/-80% in two",
      "of them (Figure S5).",
      "UNCERTAINTY METHOD. Table S2 reports SIR (sampling importance",
      "resampling) medians and 95% CIs, unlike the main text's Table 2,",
      "which reports a 1,000-sample bootstrap."
    )
  )

  ini({
    # -----------------------------------------------------------------------
    # Structural fixed effects -- Xu 2026 supplemental Table S2, 'Final
    # Estimates' column, each confirmed against the $THETA block of the
    # accompanying control stream. The paired 'SIR median [95% CI]' column is
    # quoted below for confidence in the point estimates but is NOT carried
    # into any omega: it is parameter precision, not between-subject spread.
    #
    # Reference subject: WT = 70 kg, AST = 52.5 U/L. The 70 kg reference is
    # the conventional adult anchor and lies far above every patient in the
    # n = 14 fit (median 15.9 kg), so these typical values are extrapolations
    # to a 70 kg individual. For a 15.9 kg child the model gives
    # CL = 0.477 * (15.9/70)^0.75 = 0.157 L/h and V1 = 11.8 * (15.9/70) =
    # 2.68 L.
    # -----------------------------------------------------------------------

    lcl <- log(0.477)
    label("Clearance at WT = 70 kg (L/h)")
    # Table S2: CL = 0.477 L/h (RSE 18%); SIR median 0.494
    # [0.372-0.625]. Control stream $THETA line 1: '0.477 ; CL'.

    lvc <- log(11.8)
    label("Central volume of distribution V1 at WT = 70 kg (L)")
    # Table S2: V1 = 11.8 L (RSE 15%); SIR median 12.18 [9.66-15.17].
    # Control stream $THETA line 2: '11.8 ; V1'.

    lq <- log(0.512)
    label("Intercompartmental clearance Q at WT = 70 kg and AST = 52.5 U/L (L/h)")
    # Table S2: Q = 0.512 L/h (RSE 23%); SIR median 0.55 [0.32-0.81].
    # Control stream $THETA line 3: '0.512 ; Q'. This is the typical value at
    # BOTH reference points, since the AST term is centred at 52.5 U/L.

    lvp <- log(19.2)
    label("Peripheral volume of distribution V2 at WT = 70 kg (L)")
    # Table S2: V2 = 19.2 L (RSE 41%); SIR median 19.36 [9.98-28.62].
    # Control stream $THETA line 4: '19.2 ; V2'.

    # -----------------------------------------------------------------------
    # Body-weight allometry, both exponents FIXED at the canonical values
    # (Table S2 rows read '0.75 FIX' and '1 FIX' with '/' in the RSE and SIR
    # columns), hence fixed(). ONE exponent is shared between CL and Q and
    # one between V1 and V2 -- not inferred, but printed explicitly in both
    # the Table S2 footnote equations and the control stream's $PK block.
    # -----------------------------------------------------------------------

    e_wt_cl_q <- fixed(0.75)
    label("Allometric exponent on (WT / 70) for CL and Q (unitless)")
    # Table S2 rows 'The effect of weight on CL' and '... on Q' = 0.75 FIX.
    # Control stream: TVCL = THETA(1) * (WT/70)**0.75 and
    # TVQ = THETA(3) * (WT/70)**0.75 * (AST/52.5)**0.898.

    e_wt_vc_vp <- fixed(1)
    label("Allometric exponent on (WT / 70) for V1 and V2 (unitless)")
    # Table S2 rows 'The effect of weight on V 1' and '... on V 2' = 1 FIX.
    # Control stream: TVV1 = THETA(2) * (WT/70)**1 and
    # TVV2 = THETA(4) * (WT/70)**1.

    # -----------------------------------------------------------------------
    # Hepatic-function effect on intercompartmental clearance. Estimated, so
    # no fixed() wrapper. Power form on AST centred at the subset median of
    # 52.5 U/L.
    # -----------------------------------------------------------------------

    e_ast_q <- 0.898
    label("Power exponent on (AST / 52.5) for Q (unitless)")
    # Table S2 row 'The effect of AST on Q' = 0.898 (RSE 16%); SIR median
    # 0.875 [0.655-1.075]. Control stream: (AST/52.5)**0.898. The interval
    # comfortably excludes 0 but includes 1, so the effect is close to
    # directly proportional. Absent from the paper's final 29-patient model
    # -- see covariateData$AST$notes.

    # -----------------------------------------------------------------------
    # Inter-individual variability: exponential on CL and V1, correlated;
    # NONE on Q or V2.
    #
    # This model is the reason the omega scale is unambiguous for the whole
    # paper, because the supplement publishes the %CV table and the raw
    # $OMEGA BLOCK(2) side by side:
    #
    #   Table S2 IIV CL  64.4%   $OMEGA 0.415        sqrt(0.415) = 0.6442
    #   Table S2 IIV V1  50.3%          0.253        sqrt(0.253) = 0.5030
    #   Table S2 Cor     82.4%          0.267 (cov)  0.267 / sqrt(0.415*0.253) = 0.8240
    #
    # All three reproduce the printed percentages exactly, so '%CV' here is
    # the raw omega SD times 100 and NOT the log-normal back-transform
    # omega^2 = log(1 + CV^2), which would print 71.7% for CL. The variances
    # below are therefore taken VERBATIM from $OMEGA rather than
    # back-transformed from the percentages -- no conversion is involved and
    # no rounding is introduced.
    #
    # $OMEGA additionally declares 'IIV_Q' and 'IIV_V2' as '0 FIXED', which
    # is NONMEM's way of writing 'no random effect on this parameter'. They
    # are omitted here rather than encoded as etalq ~ fixed(0) / etalvp ~
    # fixed(0): a zero-variance eta makes the OMEGA matrix singular and
    # rxode2 fails with 'chol(): decomposition failed' when simulating a
    # cohort. Omitting them is the faithful and the safe encoding.
    # -----------------------------------------------------------------------

    etalcl + etalvc ~ c(0.415,
                        0.267, 0.253)
    # Control stream $OMEGA BLOCK(2): '0.415 ; IIV_CL' and '0.267 0.253 ;
    # IIV_V1'. Table S2 prints the same block as IIV CL 64.4% (RSE 15.8%,
    # SHR 0%), IIV V1 50.3% (RSE 25.1%, SHR 2%) and Cor(CL-V1) 82.4%
    # (RSE 34.06%), with SIR medians 69.79% [50.27%-86.46%], 55.41%
    # [36.28%-74.56%] and 79.44% [70.02%-84.44%].

    # -----------------------------------------------------------------------
    # Residual error: ADDITIVE ONLY. Table S2 reports a single 'Add. err,
    # mg/L' row, and the control stream's $SIGMA declares the proportional
    # component '0 FIXED' -- so the combined proportional-plus-additive form
    # of the paper's final model does not apply here. The proportional
    # component is omitted rather than written as propSd <- fixed(0), which
    # would add a term that contributes nothing.
    #
    # Read as an SD in mg/L, following Table S2's own units column and the
    # same convention the main text's Table 2 uses (where the proportional
    # component is printed as a percent, which is only meaningful on the SD
    # scale). Note that the control stream places 1.58 in $SIGMA, where
    # NONMEM expects a VARIANCE; taken literally that would make the SD
    # sqrt(1.58) = 1.257 mg/L. The table's units column is preferred here,
    # and because the stream is a hand-assembled $DESIGN input rather than
    # an estimation run, the likelier explanation is that the authors pasted
    # the reported SD straight into $SIGMA. The discrepancy is 26% on one
    # residual SD and affects no structural parameter; see the vignette
    # Assumptions and deviations.
    # -----------------------------------------------------------------------

    addSd <- 1.58
    label("Additive residual error (mg/L)")
    # Table S2: 'Add. err, mg/L' = 1.58 (RSE 16%, SHR 14%); SIR median 1.62
    # [1.32-2.02]. Control stream $SIGMA: '1.58 ; Add.error', preceded by
    # '0 FIXED ; Prop.error'. Larger than the final model's 0.838 mg/L,
    # which is expected: with no proportional component the additive term
    # has to absorb the concentration-dependent scatter as well.
  })

  model({
    # Body-weight allometry standardised to 70 kg, and the AST effect on Q
    # centred on the n = 14 subset median of 52.5 U/L. Both forms are
    # printed verbatim in the Table S2 footnote and the control stream's
    # $PK block; see covariateData for the quotations.
    wt_cl_factor <- (WT / 70)^e_wt_cl_q
    wt_v_factor  <- (WT / 70)^e_wt_vc_vp
    ast_q_factor <- (AST / 52.5)^e_ast_q

    # Individual disposition parameters. Exponential IIV on CL and V1 only
    # -- $OMEGA fixes the Q and V2 variances to zero, so those parameters
    # are deterministic given the covariates.
    cl <- exp(lcl + etalcl) * wt_cl_factor
    vc <- exp(lvc + etalvc) * wt_v_factor
    q  <- exp(lq)           * wt_cl_factor * ast_q_factor
    vp <- exp(lvp)          * wt_v_factor

    # Linear two-compartment disposition with intravenous input into the
    # central compartment, matching the control stream's ADVAN3 TRANS4
    # (two-compartment, CL / V1 / Q / V2 parameterisation, no absorption
    # compartment). Caspofungin was given as a 1 h infusion; the duration is
    # encoded on the dose record by the user via rate or dur, exactly as the
    # control stream's $INPUT carries a RATE column.
    d/dt(central)     <- -(cl + q) / vc * central + q / vp * peripheral1
    d/dt(peripheral1) <-        q  / vc * central - q / vp * peripheral1

    # Total plasma caspofungin concentration. The control stream's $ERROR
    # block is explicit: IPRED = A(1) / V1. Dose mg, vc L -> Cc mg/L.
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
