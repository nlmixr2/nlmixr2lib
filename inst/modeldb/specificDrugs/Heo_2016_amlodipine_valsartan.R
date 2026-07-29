Heo_2016_amlodipine_valsartan <- function() {
  description <- paste(
    "Joint two-drug population PK/PD model for the antihypertensive",
    "interaction of amlodipine (calcium-channel blocker, unsuffixed",
    "parent) and valsartan (angiotensin-II receptor blocker, sibling-drug",
    "suffix _val) on systolic (SBP) and diastolic (DBP) blood pressure in",
    "healthy adult Korean volunteers receiving a single-dose fixed-dose-",
    "combination tablet of amlodipine 10 mg + valsartan 160 mg. Each drug",
    "has a two-compartment popPK model with zero-order absorption",
    "(duration D1) and theory-based allometric weight scaling on CL, Q",
    "(exponent 0.75) and V1, V2 (exponent 1) at a reference weight of 70",
    "kg. PD uses an effect-compartment Imax model on BP: amlodipine drives",
    "two separate effect compartments (SBP-side and DBP-side, different",
    "Keqs); valsartan drives one shared effect compartment (single Keq).",
    "Imax is fixed at 0.164 (the valsartan monotherapy estimate) and",
    "applies to all four drug/endpoint arms. Combined therapy uses a",
    "proportional interaction term (Heo 2016 Eq 8):",
    "  BP = BSL * (1 - PD_amlo - PD_val - alpha * PD_amlo * PD_val)",
    "with PD_x = Imax * Ce_x / (IC50_x + Ce_x). alpha < 0 = infra-additive,",
    "alpha = 0 = additive, alpha > 0 = synergistic. Estimated alpha =",
    "-0.171 for SBP and -0.0312 for DBP (both infra-additive). Combined-",
    "therapy baselines (BSL_SBP = 117 mmHg, BSL_DBP = 72.8 mmHg) and alpha",
    "are estimated on the 48-subject FDC dataset; monotherapy IC50s and",
    "Keqs are fixed at their monotherapy point estimates (Tables 1 and 2)."
  )
  reference <- paste(
    "Heo YA, Holford N, Kim Y, Son M, Park K.",
    "Quantitative model for the blood pressure-lowering interaction of",
    "valsartan and amlodipine.",
    "Br J Clin Pharmacol. 2017 Jul;83(7):1502-1514.",
    "doi:10.1111/bcp.13082."
  )
  vignette <- "Heo_2016_amlodipine_valsartan"

  paper_specific_compartments <- c("effect_sbp", "effect_dbp")

  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-invariant. Used for theory-based allometric scaling on CL,",
        "Q (exponent 0.75) and V1, V2 (exponent 1) for both drugs, with",
        "reference weight 70 kg (Heo 2016 Methods, Eq 1-2). Cohort mean",
        "68.7 kg, SD 7.63 kg (Table S3)."
      ),
      source_name        = "WT"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age", units = "years", type = "continuous",
      notes = "Screened by stepwise covariate modelling (SCM) but not retained in the final popPK model (Heo 2016 Results: 'Apart from body weight, no other covariate factors were found to be significant (P > 0.001).')."
    ),
    CREAT = list(
      description = "Serum creatinine", units = "mg/dL", type = "continuous",
      notes = "Screened but not retained in the final model."
    ),
    SMOKE = list(
      description = "Smoking status", units = "(binary)", type = "binary",
      notes = "Screened but not retained in the final model."
    ),
    DRINK = list(
      description = "Alcohol drinking status", units = "(binary)", type = "binary",
      notes = "Screened but not retained in the final model."
    ),
    AST = list(
      description = "Aspartate aminotransferase", units = "IU/L", type = "continuous",
      notes = "Screened for amlodipine only (extensive hepatic metabolism) but not retained."
    ),
    ALT = list(
      description = "Alanine aminotransferase", units = "IU/L", type = "continuous",
      notes = "Screened for amlodipine only but not retained."
    ),
    GGT = list(
      description = "Gamma-glutamyl transferase", units = "IU/L", type = "continuous",
      notes = "Screened for amlodipine only but not retained."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 48,
    n_studies      = 1,
    age_range      = "not reported (mean 29 years, SD 5.98; Table S3)",
    age_median     = "29 years (mean)",
    weight_range   = "not reported (mean 68.7 kg, SD 7.63; Table S3)",
    weight_median  = "68.7 kg (mean)",
    sex_female_pct = 0,
    race_ethnicity = c(Asian = 100),
    disease_state  = "healthy volunteers",
    dose_range     = "single oral dose of a fixed-dose-combination tablet of amlodipine 10 mg + valsartan 160 mg (as besylate or orotate; crossover)",
    regions        = "South Korea",
    notes = paste(
      "Combined-administration PK/PD data were taken from a 48-subject",
      "two-way crossover FDC bioequivalence study (Kim 2013 Clin Ther",
      "35:934-940; Heo 2016 ref [7]) at Yonsei University Severance",
      "Hospital. Concentration sampling: 0 (predose) and 0.5, 1, 1.5,",
      "2, 3, 4, 6, 8, 10, 12, 16, 24, 48 h (both drugs) plus 72, 96,",
      "144, 192 h (amlodipine only). BP sampling: 0, 2, 4, 8, 12, 16,",
      "24 h -- only the first 24 h were used for the PD analysis. 816",
      "measurements of amlodipine (0.61% BLQ) and 624 of valsartan",
      "(6.5% BLQ); 336 SBP + 336 DBP with no missing values. Monotherapy",
      "PD baselines and IC50 / Keq estimates come from literature-mean",
      "SBP / DBP time courses in single-dose healthy-volunteer studies",
      "(amlodipine: Kim 2010; valsartan: Czendlik 1997; naive pooled",
      "PD analysis on the mean profiles). PPP&D sequential PKPD approach",
      "(PK parameters fixed at their point estimates during PD",
      "estimation). All fits used NONMEM 7.3 with FOCE + interaction."
    )
  )

  ini({
    # ==================================================================
    # AMLODIPINE (parent, unsuffixed) POPULATION PK
    # Heo 2016 Supplementary Table S4; values for a 70 kg subject.
    # Two-compartment disposition with zero-order absorption directly
    # into central (dose duration = D1) and first-order elimination.
    # ==================================================================
    lcl <- log(39.4); label("Amlodipine apparent clearance CL/F for a 70 kg adult (L/h)")             # Heo 2016 Table S4 THETA CL/F = 39.4 (RSE 3.69%)
    lvc <- log(1620); label("Amlodipine apparent central volume V1/F for a 70 kg adult (L)")          # Heo 2016 Table S4 THETA V1/F = 1620 (RSE 4.17%)
    lq  <- log(45.4); label("Amlodipine apparent intercompartmental clearance Q/F for a 70 kg adult (L/h)") # Heo 2016 Table S4 THETA Q/F = 45.4 (RSE 16.6%)
    lvp <- log(588);  label("Amlodipine apparent peripheral volume V2/F for a 70 kg adult (L)")       # Heo 2016 Table S4 THETA V2/F = 588 (RSE 8.58%)
    ld1 <- log(5.28); label("Amlodipine zero-order absorption duration D1 into central (h)")          # Heo 2016 Table S4 THETA D1 = 5.28 (RSE 3.07%)

    # ==================================================================
    # VALSARTAN (sibling-drug suffix _val) POPULATION PK
    # Heo 2016 Supplementary Table S4; values for a 70 kg subject.
    # Two-compartment disposition with zero-order absorption directly
    # into central (dose duration = D1_val).
    # ==================================================================
    lcl_val <- log(6.18); label("Valsartan apparent clearance CL/F for a 70 kg adult (L/h)")         # Heo 2016 Table S4 THETA CL/F = 6.18 (RSE 5.76%)
    lvc_val <- log(25.9); label("Valsartan apparent central volume V1/F for a 70 kg adult (L)")      # Heo 2016 Table S4 THETA V1/F = 25.9 (RSE 9.73%)
    lq_val  <- log(2.01); label("Valsartan apparent intercompartmental clearance Q/F for a 70 kg adult (L/h)") # Heo 2016 Table S4 THETA Q/F = 2.01 (RSE 14.5%)
    lvp_val <- log(17.4); label("Valsartan apparent peripheral volume V2/F for a 70 kg adult (L)")   # Heo 2016 Table S4 THETA V2/F = 17.4 (RSE 16.0%)
    ld1_val <- log(4.39); label("Valsartan zero-order absorption duration D1 into central (h)")      # Heo 2016 Table S4 THETA D1 = 4.39 (RSE 5.83%)

    # ==================================================================
    # Allometric weight scaling -- fixed canonical exponents applied to
    # both drugs. Reference weight 70 kg (Heo 2016 Methods Eq 1-2).
    # ==================================================================
    e_wt_cl_q  <- fixed(0.75); label("Allometric exponent on CL and Q (unitless)")   # Heo 2016 Eq 1 FsizeCL: (WT/70)^0.75
    e_wt_vc_vp <- fixed(1);    label("Allometric exponent on V1 and V2 (unitless)")  # Heo 2016 Eq 2 FsizeV: (WT/70)^1

    # ==================================================================
    # PD structural parameters
    #
    # Baselines come from the combined-therapy fit (Heo 2016 Table 3;
    # 48-subject FDC dataset). Imax is fixed at 0.164 for both drugs
    # (Heo 2016 Table 2 valsartan monotherapy estimate; the same value
    # was fixed for amlodipine because the monotherapy data were not
    # informative about amlodipine Imax -- Heo 2016 Results, Table 1
    # footer). Monotherapy IC50s and Keqs (Tables 1 and 2) are fixed
    # at their point estimates during the combined-therapy fit; only
    # baselines and alpha are re-estimated on the combined data (Heo
    # 2016 Table 3 footer: "other parameters (Imax, etc.) fixed to
    # values estimated from combined administration PK models and
    # single drug PD models").
    # ==================================================================
    lrbase_sbp <- log(117);  label("Combined-therapy baseline SBP (mmHg)")   # Heo 2016 Table 3 SBP BSL = 117 (RSE 1.10%)
    lrbase_dbp <- log(72.8); label("Combined-therapy baseline DBP (mmHg)")   # Heo 2016 Table 3 DBP BSL = 72.8 (RSE 1.33%)

    limax <- fixed(log(0.164)); label("Maximum fractional inhibition Imax (shared across both drugs and both endpoints)") # Heo 2016 Table 2 valsartan Imax = 0.164 (RSE 102%); fixed at the same value for amlodipine in Table 1 (footer)

    # Amlodipine PD parameters -- separate SBP-side and DBP-side effect
    # compartments (different Keqs and IC50s).
    lec50_sbp <- fixed(log(8.27));  label("Amlodipine SBP IC50 (ng/mL); fixed at monotherapy value in the combined fit")  # Heo 2016 Table 1 SBP IC50 = 8.27 (RSE 32.2%)
    lec50_dbp <- fixed(log(2.97));  label("Amlodipine DBP IC50 (ng/mL); fixed at monotherapy value in the combined fit")  # Heo 2016 Table 1 DBP IC50 = 2.97 (RSE 36.8%)
    lke0_sbp  <- fixed(log(0.211)); label("Amlodipine SBP effect-compartment equilibration rate Keq (1/h); fixed at monotherapy value") # Heo 2016 Table 1 SBP Keq = 0.211 (RSE 38.3%) [teq = ln(2)/Keq = 3.3 h]
    lke0_dbp  <- fixed(log(0.821)); label("Amlodipine DBP effect-compartment equilibration rate Keq (1/h); fixed at monotherapy value") # Heo 2016 Table 1 DBP Keq = 0.821 (RSE 104.7%) [teq = 0.84 h]

    # Valsartan PD parameters -- one effect compartment shared between
    # SBP and DBP (single Keq, single IC50).
    lec50_val <- fixed(log(1200));  label("Valsartan IC50 (ng/mL, shared SBP+DBP); fixed at monotherapy value in the combined fit") # Heo 2016 Table 2 IC50 = 1200 (RSE 53.3%)
    lke0_val  <- fixed(log(0.542)); label("Valsartan effect-compartment equilibration rate Keq (1/h, shared SBP+DBP); fixed at monotherapy value") # Heo 2016 Table 2 Keq = 0.542 (RSE 59.5%) [teq = 1.3 h]

    # Proportional drug-drug interaction terms (Heo 2016 Eq 8).
    # Not log-transformed because alpha can be negative (infra-additive).
    alpha_sbp <- -0.171;  label("Proportional PD interaction term for SBP (dimensionless; < 0 = infra-additive)") # Heo 2016 Table 3 SBP ALPHA = -0.171 (RSE 11.6%; 95% CI -0.218 to -0.143)
    alpha_dbp <- -0.0312; label("Proportional PD interaction term for DBP (dimensionless; < 0 = infra-additive)") # Heo 2016 Table 3 DBP ALPHA = -0.0312 (RSE 57.8%; 95% CI -0.07739 to -0.00283)

    # ==================================================================
    # IIV -- log-normal on the log-transformed parameters. PPV entries
    # in the source tables are sqrt(NONMEM OMEGA); variance = PPV^2
    # (Table 3 caption: "PPV value calculated from sqrt(NONMEM OMEGA
    # estimate)").
    # ==================================================================

    # Amlodipine PK IIV (Heo 2016 Table S4). Text also mentions a CL-V1
    # correlation was estimated but the numeric value is not reported
    # in the paper; encoded as diagonal here and flagged in vignette
    # Assumptions & deviations.
    etalcl ~ 0.239^2 # PPV CL/F  = 0.239 (RSE 12.8%) -> variance 0.057121
    etalvc ~ 0.203^2 # PPV V1/F  = 0.203 (RSE 11.9%) -> variance 0.041209
    etald1 ~ 0.233^2 # PPV D1    = 0.233 (RSE 11.2%) -> variance 0.054289

    # Valsartan PK IIV (Heo 2016 Table S4). D1 was fixed to 0 IIV in
    # the source. Same CL-V1 correlation gap as amlodipine.
    etalcl_val ~ 0.373^2 # PPV CL/F = 0.373 (RSE 9.68%) -> variance 0.139129
    etalvc_val ~ 0.519^2 # PPV V1/F = 0.519 (RSE 12.9%) -> variance 0.269361

    # Baseline BP IIV (Heo 2016 Table 3).
    etalrbase_sbp ~ 0.059^2 # PPV BSL_SBP = 0.059 (RSE 18.8%) -> variance 0.003481
    etalrbase_dbp ~ 0.071^2 # PPV BSL_DBP = 0.071 (RSE 15.9%) -> variance 0.005041

    # PD IIVs on IC50 and Keq (Heo 2016 Table 3). The paper fit SBP
    # and DBP as two separate models with distinct subject-level etas;
    # in this joint encoding each parameter carries its own eta while
    # sharing the population value where the paper did. The reported
    # PPVs on IC50 and Keq are very large with high RSE (up to 100%+)
    # and are essentially unidentified from the 48-subject data;
    # retained for source fidelity.
    etalec50_sbp     ~ 0.468^2 # PPV IC50_1 (aml, SBP model) = 0.468 (RSE 84.3%) -> variance 0.219024
    etalec50_val_sbp ~ 1.732^2 # PPV IC50_2 (val, SBP model) = 1.732 (RSE 43.4%) -> variance 2.999824
    etalke0_sbp      ~ 1.319^2 # PPV AML Keq (SBP model)     = 1.319 (RSE 46.5%) -> variance 1.739761
    etalke0_val_sbp  ~ 0.376^2 # PPV VAL Keq (SBP model)     = 0.376 (RSE 37.9%) -> variance 0.141376

    etalec50_dbp     ~ 1.612^2 # omega IC50_1 (aml, DBP model) = 1.612 (RSE 33.1%) -> variance 2.598544
    etalec50_val_dbp ~ 1.131^2 # omega IC50_2 (val, DBP model) = 1.131 (RSE 32.8%) -> variance 1.279161

    # ==================================================================
    # Residual error -- reported as sqrt(NONMEM SIGMA); used directly
    # as standard-deviation values. Amlodipine and valsartan PK
    # residuals from Heo 2016 Table S4 (additive + proportional
    # combined error). SBP / DBP additive residuals from Heo 2016
    # Table 3 (combined-therapy fit).
    # ==================================================================
    propSd     <- 0.126; label("Amlodipine proportional residual SD (fraction)")   # Heo 2016 Table S4 sigma_prop = 0.126 (RSE 13.8%)
    addSd      <- 0.242; label("Amlodipine additive residual SD (ng/mL)")          # Heo 2016 Table S4 sigma_add = 0.242 (RSE 9.69%)
    propSd_val <- 0.348; label("Valsartan proportional residual SD (fraction)")    # Heo 2016 Table S4 sigma_prop = 0.348 (RSE 5.17%)
    addSd_val  <- 32;    label("Valsartan additive residual SD (ng/mL)")           # Heo 2016 Table S4 sigma_add = 32 (RSE 24.3%)
    addSd_SBP  <- 7.84;  label("Combined-therapy SBP additive residual SD (mmHg)") # Heo 2016 Table 3 SBP sigma_add = 7.84 (RSE 6.36%)
    addSd_DBP  <- 4.46;  label("Combined-therapy DBP additive residual SD (mmHg)") # Heo 2016 Table 3 DBP sigma_add = 4.46 (RSE 5.22%)
  })

  model({
    # ----------------------------------------------------------------
    # 1. Individual PK parameters -- theory-based allometric weight
    #    scaling on CL, Q (exp 0.75) and V1, V2 (exp 1) at reference
    #    70 kg (Heo 2016 Eq 1-2). D1 is not weight-scaled.
    # ----------------------------------------------------------------
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    q  <- exp(lq)           * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp)          * (WT / 70)^e_wt_vc_vp
    d1 <- exp(ld1 + etald1)

    cl_val <- exp(lcl_val + etalcl_val) * (WT / 70)^e_wt_cl_q
    vc_val <- exp(lvc_val + etalvc_val) * (WT / 70)^e_wt_vc_vp
    q_val  <- exp(lq_val)               * (WT / 70)^e_wt_cl_q
    vp_val <- exp(lvp_val)              * (WT / 70)^e_wt_vc_vp
    d1_val <- exp(ld1_val)

    # ----------------------------------------------------------------
    # 2. Individual PD parameters
    # ----------------------------------------------------------------
    bsl_sbp <- exp(lrbase_sbp + etalrbase_sbp)
    bsl_dbp <- exp(lrbase_dbp + etalrbase_dbp)
    imax    <- exp(limax)

    ec50_sbp     <- exp(lec50_sbp + etalec50_sbp)
    ec50_dbp     <- exp(lec50_dbp + etalec50_dbp)
    ec50_val_sbp <- exp(lec50_val + etalec50_val_sbp)
    ec50_val_dbp <- exp(lec50_val + etalec50_val_dbp)

    ke0_sbp <- exp(lke0_sbp + etalke0_sbp)
    ke0_dbp <- exp(lke0_dbp)
    ke0_val <- exp(lke0_val + etalke0_val_sbp)

    # ----------------------------------------------------------------
    # 3. Micro-constants
    # ----------------------------------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    kel_val <- cl_val / vc_val
    k12_val <- q_val  / vc_val
    k21_val <- q_val  / vp_val

    # ----------------------------------------------------------------
    # 4. ODE system -- amlodipine two-compartment PK with zero-order
    #    absorption directly into central (dose duration = d1). Two
    #    effect compartments (SBP-side and DBP-side, different Keqs).
    # ----------------------------------------------------------------
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1
    dur(central)      <- d1

    d/dt(effect_sbp) <- ke0_sbp * (central / vc - effect_sbp)
    d/dt(effect_dbp) <- ke0_dbp * (central / vc - effect_dbp)

    # Valsartan two-compartment PK with zero-order absorption; one
    # effect compartment shared between SBP and DBP effect terms.
    d/dt(central_val)     <- -kel_val * central_val - k12_val * central_val + k21_val * peripheral1_val
    d/dt(peripheral1_val) <-                          k12_val * central_val - k21_val * peripheral1_val
    dur(central_val)      <- d1_val

    d/dt(effect_val) <- ke0_val * (central_val / vc_val - effect_val)

    # ----------------------------------------------------------------
    # 5. Observation variables. central and central_val hold drug
    #    amount (mg); central/vc is mg/L -> multiply by 1000 to get
    #    ng/mL (matches Heo 2016 concentration units in Table S4 and
    #    the ng/mL IC50 units in Tables 1-3). The effect compartments
    #    equilibrate with the plasma concentration (same amount/L
    #    scale as central/vc) so we scale by 1000 in the same way to
    #    get an effect-site plasma-equivalent concentration in ng/mL.
    # ----------------------------------------------------------------
    Cc     <- (central     / vc)     * 1000
    Cc_val <- (central_val / vc_val) * 1000

    ce_amlo_sbp <- effect_sbp * 1000
    ce_amlo_dbp <- effect_dbp * 1000
    ce_val_p    <- effect_val * 1000

    # ----------------------------------------------------------------
    # 6. PD -- effect-compartment Imax model on each endpoint with a
    #    proportional drug-drug interaction (Heo 2016 Eq 3, 6, 8):
    #      PD_x = imax * Ce_x / (IC50_x + Ce_x)
    #      combined = PD_amlo + PD_val + alpha * PD_amlo * PD_val
    #      BP       = BSL * (1 - combined)
    #    Because valsartan's single effect compartment feeds both SBP
    #    and DBP, ce_val_p is reused; the SBP-model and DBP-model
    #    subject-level ec50_val values are propagated via separate
    #    etas (etalec50_val_sbp, etalec50_val_dbp) on the shared
    #    population lec50_val.
    # ----------------------------------------------------------------
    pd_amlo_sbp <- imax * ce_amlo_sbp / (ec50_sbp     + ce_amlo_sbp)
    pd_val_sbp  <- imax * ce_val_p    / (ec50_val_sbp + ce_val_p)
    eff_sbp     <- pd_amlo_sbp + pd_val_sbp + alpha_sbp * pd_amlo_sbp * pd_val_sbp

    pd_amlo_dbp <- imax * ce_amlo_dbp / (ec50_dbp     + ce_amlo_dbp)
    pd_val_dbp  <- imax * ce_val_p    / (ec50_val_dbp + ce_val_p)
    eff_dbp     <- pd_amlo_dbp + pd_val_dbp + alpha_dbp * pd_amlo_dbp * pd_val_dbp

    SBP <- bsl_sbp * (1 - eff_sbp)
    DBP <- bsl_dbp * (1 - eff_dbp)

    # ----------------------------------------------------------------
    # 7. Residual error models
    # ----------------------------------------------------------------
    Cc     ~ add(addSd)     + prop(propSd)
    Cc_val ~ add(addSd_val) + prop(propSd_val)
    SBP    ~ add(addSd_SBP)
    DBP    ~ add(addSd_DBP)
  })
}
