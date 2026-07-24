Kemal_2026_nemtabrutinib <- function() {
  description <- "Two-compartment population PK model for nemtabrutinib (oral BTK inhibitor) in adults with hematologic malignancies including CLL/SLL (Kemal 2026, full covariate model)"
  reference <- "Kemal CC, Zweers TJ, Krekels EHJ, Chatterjee MS. Population Pharmacokinetic Modeling and Exposure-Response Analyses of Nemtabrutinib in Patients With Hematologic Malignancies. CPT Pharmacometrics Syst Pharmacol. 2026;15(5). doi:10.1002/psp4.70257"
  vignette <- "Kemal_2026_nemtabrutinib"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference (median) 73.5 kg; NONMEM supplement code. Applied as power function to CL and Vc with fixed exponents 0.331 and 0.807 respectively (Table 2). Missing values imputed to the reference in the source analysis.",
      source_name        = "BWT"
    ),
    AGE = list(
      description        = "Baseline age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference (median) 68 years; NONMEM supplement code. Log-linear (power) effect on CL only (exponent -0.503). Missing values imputed to the reference.",
      source_name        = "AGE"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference (median) 41.2 g/L; NONMEM supplement code (paper Table 1 reports in g/dL, i.e. 4.12 g/dL). Log-linear (power) effect on CL only (exponent -0.395). Missing values imputed to the reference.",
      source_name        = "BALB"
    ),
    SEXF = list(
      description        = "Biological sex, 1 = female / 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = male",
      notes              = "Fractional additive effect on CL (-0.133) and Vc (-0.0901) for female relative to male reference. NONMEM SEX column: 0 = male, 1 = female.",
      source_name        = "SEX"
    ),
    RACE_ASIAN = list(
      description        = "Race indicator, 1 = Asian / 0 = other",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = White (Other pooled with White per Methods)",
      notes              = "Fractional additive effect on CL (-0.117) and Vc (-0.116). Derived from RACE = 3 in the NONMEM control stream; White (RACE = 1) is the reference and 'Other' was merged with 'White' per Methods section 2.2.2.",
      source_name        = "RACE"
    ),
    RACE_BLACK = list(
      description        = "Race indicator, 1 = Black / African American / 0 = other",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = White (Other pooled with White per Methods)",
      notes              = "Fractional additive effect on CL (0.0533; 95% CI includes zero, RSE 252%) and Vc (0.102). Derived from RACE = 2 in the NONMEM control stream.",
      source_name        = "RACE"
    ),
    DIS_BCELLNHL = list(
      description        = "Disease indication, 1 = B-cell non-Hodgkin lymphoma / 0 = other",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = CLL/SLL (reference disease indication in the model)",
      notes              = "Fractional additive effect on CL (-0.166) and Vc (0.00224). Derived from INDIC = 1 in the NONMEM control stream.",
      source_name        = "INDIC"
    ),
    DIS_WM = list(
      description        = "Disease indication, 1 = Waldenstrom's macroglobulinemia / 0 = other",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = CLL/SLL",
      notes              = "Fractional additive effect on CL (0.0718; RSE 120%) and Vc (0.152). Derived from INDIC = 3 in the NONMEM control stream.",
      source_name        = "INDIC"
    ),
    DIS_OTHER_HEME = list(
      description        = "Disease indication, 1 = other hematologic malignancy (pooled Other + MZL + FL + MCL + Richter's transformation) / 0 = otherwise",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = CLL/SLL",
      notes              = "Fractional additive effect on CL (-0.0244; RSE 171%) and Vc (-0.0200). Derived from INDIC in {4, 5, 6, 7, 8} in the NONMEM control stream, which pools 'Other', MZL (marginal zone lymphoma), FL (follicular lymphoma), MCL (mantle cell lymphoma), and Richter's transformation into a single effect.",
      source_name        = "INDIC"
    ),
    RENALIMP_MILD = list(
      description        = "Renal impairment, 1 = mild (Cockcroft-Gault eGFR 60-89 mL/min) / 0 = otherwise",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = normal renal function (eGFR >= 90 mL/min); mild impairment is a distinct binary indicator from RENALIMP_MOD",
      notes              = "Fractional additive effect on CL (0.0537). Classification per NCI ODWG using Cockcroft-Gault (Table S2). Missing serum creatinine assumed normal in the NONMEM code.",
      source_name        = "RENIMP"
    ),
    RENALIMP_MOD = list(
      description        = "Renal impairment, 1 = moderate (eGFR 30-59 mL/min) / 0 = otherwise",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = normal renal function; mild impairment enters via RENALIMP_MILD",
      notes              = "Fractional additive effect on CL (-0.00186; 95% CI wide, RSE ~3000%). Cockcroft-Gault classification per Table S2. Severe impairment and kidney failure were pooled with the reference in the model (only 11 subjects total).",
      source_name        = "RENIMP"
    ),
    HEPIMP_MILD = list(
      description        = "Hepatic impairment, 1 = mild (NCI ODWG) / 0 = otherwise",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = normal hepatic function; moderate hepatic impairment (n = 8) pooled with reference",
      notes              = "Fractional additive effect on CL (0.00401; 95% CI includes zero, RSE ~1200%). NCI ODWG classification per Table S2.",
      source_name        = "HEPIMP"
    ),
    CONMED_CYP3A4_IND_MOD = list(
      description        = "Concomitant use of a moderate CYP3A4 inducer at the time of the PK observation",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = no concomitant moderate CYP3A4 inducer",
      notes              = "Fractional additive effect on CL (0.0220; RSE 458%). Time-varying regressor. No strong CYP3A4 inducers in the study; weak inducers were not entered into the model.",
      source_name        = "C3A4INDM"
    ),
    CONMED_CYP3A4_INH_STRONG = list(
      description        = "Concomitant use of a strong CYP3A4 inhibitor at the time of the PK observation",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = no concomitant strong CYP3A4 inhibitor",
      notes              = "Fractional additive effect on CL (-0.0119; RSE 774%). Time-varying regressor. Weak and moderate inhibitors were not entered into the model.",
      source_name        = "C3A4INHS"
    ),
    CONMED_PPI = list(
      description        = "Concomitant use of a proton pump inhibitor at the time of the PK observation",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = no concomitant PPI",
      notes              = "Fractional additive effect on F (0.00232; RSE 915%). Time-varying regressor. Examples: esomeprazole, pantoprazole.",
      source_name        = "PPIF"
    ),
    CONMED_H2RA = list(
      description        = "Concomitant use of a histamine H2 receptor antagonist at the time of the PK observation",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = no concomitant H2 antagonist",
      notes              = "Fractional additive effect on F (0.0264; RSE 132%). Time-varying regressor. Examples: famotidine, cimetidine.",
      source_name        = "H2F"
    ),
    CONMED_ANTACID = list(
      description        = "Concomitant use of an antacid at the time of the PK observation",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = no concomitant antacid",
      notes              = "Fractional additive effect on F (-0.0360; RSE 72.8%). Time-varying regressor. Examples: calcium carbonate, magnesium carbonate.",
      source_name        = "ANTACIDF"
    ),
    DOSE = list(
      description        = "Nemtabrutinib dose at the current dosing event, used to identify low-dose (< 30 mg) bioavailability regime",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used only to construct the low-dose indicator: doses < 30 mg receive a -0.151 fractional shift in F relative to doses >= 30 mg. Both the paper and the NONMEM supplement encode this as a dose-level covariate, not a subject-level covariate.",
      source_name        = "DOSE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 578L,
    n_studies      = 2L,
    age_range      = "25-89 years",
    age_median     = "68 years",
    weight_range   = "41.2-147 kg",
    weight_median  = "74.0 kg",
    sex_female_pct = 34.1,
    race_ethnicity = c(White = 85.8, Black = 2.4, Asian = 7.8, Other = 2.4, Missing = 1.6),
    disease_state  = "Hematologic malignancies including CLL/SLL (49.8%), B-cell NHL (9.7%), Waldenstrom's macroglobulinemia (9.7%), and other (30.6%; MZL / FL / MCL / Richter's transformation pooled)",
    dose_range     = "5-80 mg once daily orally (BELLWAVE-001: 5, 10, 15, 20, 30, 45, 65, 75 mg; BELLWAVE-003: 45, 65, 80 mg)",
    regions        = "Not stated at the aggregate level in the source",
    notes          = "Table 1 baseline demographics for the pooled BELLWAVE-001 (Phase 1/2, n = 136) and BELLWAVE-003 (Phase 2, n = 442) analysis population. 5669 non-BLQ observations across 578 patients. Missing continuous covariates imputed to the median; missing / small-N categorical categories pooled with the largest category ('Other' race merged with 'White'; 'Other' disease indication kept as a distinct category)."
  )

  ini({
    # =================================================================
    # Structural parameters -- final estimates from Table 2 of Kemal
    # et al. 2026. Reference subject: 73.5 kg body weight, 68 years,
    # 41.2 g/L albumin, White male with CLL/SLL, normal renal and
    # hepatic function, no CYP3A4 modulators or acid-reducing agents,
    # dose >= 30 mg.
    # NONMEM ADVAN4 TRANS4 (depot + central + one peripheral) with
    # first-order absorption, absorption lag, and first-order
    # elimination.
    # =================================================================
    lcl   <- log(3.33);  label("Apparent clearance CL/F (L/h)")                    # Table 2 point estimate (95% CI 3.03-3.63; RSE 4.59%)
    lvc   <- log(120);   label("Apparent central volume Vc/F (L)")                 # Table 2 point estimate (95% CI 116-124; RSE 1.82%)
    lka   <- log(2.83);  label("First-order absorption rate Ka (1/h)")             # Table 2 point estimate (95% CI 2.43-3.22; RSE 7.16%)
    lq    <- log(0.681); label("Apparent inter-compartmental clearance Q/F (L/h)") # Table 2 point estimate (95% CI 0.126-1.24; RSE 41.6%)
    lvp   <- log(66.9);  label("Apparent peripheral volume Vp/F (L)")              # Table 2 point estimate (95% CI 44.2-89.6; RSE 17.3%)
    lalag <- log(0.494); label("Absorption lag time (h)")                          # Table 2 point estimate (95% CI 0.418-0.569; RSE 7.79%)

    # =================================================================
    # Body weight power exponents on CL and Vc. Fitted first with a
    # simple weight-only model, then FIXED in the full covariate model
    # (paper Section 3.2 and NONMEM code, THETA(8) and THETA(9) FIX).
    # Reference weight 73.5 kg.
    # =================================================================
    e_wt_cl <- fixed(0.331);  label("Body weight power exponent on CL/F (unitless, ref 73.5 kg)")   # Table 2, fixed after initial WT-only fit
    e_wt_vc <- fixed(0.807);  label("Body weight power exponent on Vc/F (unitless, ref 73.5 kg)")   # Table 2, fixed after initial WT-only fit

    # =================================================================
    # Continuous covariate power exponents (log-linear) on CL.
    # Reference values: age 68 years, albumin 41.2 g/L.
    # =================================================================
    e_age_cl <- -0.503; label("Age power exponent on CL/F (unitless, ref 68 y)")           # Table 2 (95% CI -0.774 to -0.231; RSE 27.6%)
    e_alb_cl <- -0.395; label("Baseline albumin power exponent on CL/F (unitless, ref 41.2 g/L)") # Table 2 (95% CI -0.688 to -0.101; RSE 38.0%)

    # =================================================================
    # Categorical covariate fractional effects on CL (multiplicative,
    # applied as CL = TVCL * (1 + effect) for each active indicator;
    # reference categories accumulate the '1'). Full covariate model
    # retains all indicators regardless of statistical significance
    # (paper Section 2.2.2).
    # =================================================================
    e_sexf_cl         <- -0.133;   label("Sex effect on CL/F: female vs male reference (fraction)")             # Table 2 (95% CI -0.198 to -0.0680; RSE 24.9%)
    e_race_black_cl   <-  0.0533;  label("Race effect on CL/F: Black vs White reference (fraction)")            # Table 2 (95% CI -0.210 to 0.317; RSE 252%)
    e_race_asian_cl   <- -0.117;   label("Race effect on CL/F: Asian vs White reference (fraction)")            # Table 2 (95% CI -0.242 to 0.00921; RSE 55.1%)
    e_dis_bcellnhl_cl <- -0.166;   label("Disease effect on CL/F: B-cell NHL vs CLL/SLL reference (fraction)")  # Table 2 (95% CI -0.262 to -0.0701; RSE 29.5%)
    e_dis_wm_cl       <-  0.0718;  label("Disease effect on CL/F: Waldenstrom vs CLL/SLL reference (fraction)") # Table 2 (95% CI -0.0969 to 0.241; RSE 120%)
    e_dis_other_cl    <- -0.0244;  label("Disease effect on CL/F: other heme malignancy vs CLL/SLL (fraction)") # Table 2 (95% CI -0.106 to 0.0576; RSE 171%)
    e_renmild_cl      <-  0.0537;  label("Renal impairment effect on CL/F: mild vs normal (fraction)")          # Table 2 (95% CI -0.0421 to 0.149; RSE 91.0%)
    e_renmod_cl       <- -0.00186; label("Renal impairment effect on CL/F: moderate vs normal (fraction)")      # Table 2 (95% CI -0.111 to 0.107; RSE 2982%)
    e_hepmild_cl      <-  0.00401; label("Hepatic impairment effect on CL/F: mild vs normal (fraction)")        # Table 2 (95% CI -0.0897 to 0.0977; RSE 1192%)
    e_cyp3a4indmod_cl <-  0.0220;  label("CYP3A4 induction effect on CL/F: moderate inducer vs none (fraction)") # Table 2 (95% CI -0.175 to 0.219; RSE 458%)
    e_cyp3a4inhstr_cl <- -0.0119;  label("CYP3A4 inhibition effect on CL/F: strong inhibitor vs none (fraction)") # Table 2 (95% CI -0.193 to 0.169; RSE 774%)

    # Fractional effects on Vc.
    e_sexf_vc         <- -0.0901;  label("Sex effect on Vc/F: female vs male reference (fraction)")             # Table 2 (95% CI -0.124 to -0.0565; RSE 19.1%)
    e_race_black_vc   <-  0.102;   label("Race effect on Vc/F: Black vs White reference (fraction)")            # Table 2 (95% CI -0.0724 to 0.277; RSE 87.2%)
    e_race_asian_vc   <- -0.116;   label("Race effect on Vc/F: Asian vs White reference (fraction)")            # Table 2 (95% CI -0.181 to -0.0516; RSE 28.3%)
    e_dis_bcellnhl_vc <-  0.00224; label("Disease effect on Vc/F: B-cell NHL vs CLL/SLL reference (fraction)")  # Table 2 (95% CI -0.0674 to 0.0719; RSE 1590%)
    e_dis_wm_vc       <-  0.152;   label("Disease effect on Vc/F: Waldenstrom vs CLL/SLL reference (fraction)") # Table 2 (95% CI 0.0717 to 0.232; RSE 26.9%)
    e_dis_other_vc    <- -0.0200;  label("Disease effect on Vc/F: other heme malignancy vs CLL/SLL (fraction)") # Table 2 (95% CI -0.0598 to 0.0198; RSE 102%)

    # Fractional effects on bioavailability F.
    e_dose_lt30_f <- -0.151;  label("Low-dose (<30 mg) effect on F: fractional shift vs >=30 mg reference")     # Table 2 (95% CI -0.225 to -0.0777; RSE 24.8%)
    e_ppi_f       <-  0.00232; label("PPI effect on F (fraction)")                                              # Table 2 (95% CI -0.0392 to 0.0438; RSE 915%)
    e_h2ra_f      <-  0.0264;  label("H2 antagonist effect on F (fraction)")                                    # Table 2 (95% CI -0.0419 to 0.0946; RSE 132%)
    e_antacid_f   <- -0.0360;  label("Antacid effect on F (fraction)")                                          # Table 2 (95% CI -0.0873 to 0.0154; RSE 72.8%)

    # =================================================================
    # Inter-individual variability (log-normal). Table 2 reports CV%;
    # converted to log-scale variance via omega^2 = log(1 + CV^2).
    #   CL/F: 39.8% CV  -> omega^2 = log(1 + 0.398^2) = 0.14710
    #   Vc/F: 17.1% CV  -> omega^2 = log(1 + 0.171^2) = 0.02884
    # NONMEM $OMEGA in the supplement is diagonal (no covariance
    # estimated between ETA(CL) and ETA(Vc)).
    # =================================================================
    etalcl ~ 0.14710   # Table 2, CL/F IIV 39.8% CV (95% CI 36.1-43.3; RSE 4.60%)
    etalvc ~ 0.02884   # Table 2, Vc/F IIV 17.1% CV (95% CI 15.1-18.8; RSE 5.56%)

    # =================================================================
    # Residual error. NONMEM code:
    #   W = SQRT(THETA(32)^2 * IPRED^2 + THETA(33)^2)
    #   Y = IPRED + W * EPS(1), $SIGMA 1 FIX
    # so THETA(32) is the proportional SD (fraction) and THETA(33) is
    # the additive SD (ng/mL). Additive term is very poorly identified
    # (RSE 101%, 95% CI includes zero) but retained in the full model.
    # =================================================================
    propSd <- 0.222; label("Proportional residual error SD (fraction)")   # Table 2 (95% CI 0.208-0.236; RSE 3.19%)
    addSd  <- 3.19;  label("Additive residual error SD (ng/mL)")          # Table 2 (95% CI -3.11 to 9.4; RSE 101% -- retained in the full model)
  })

  model({
    # =================================================================
    # Multiplicative covariate blocks. Each categorical indicator adds
    # its fractional effect to a running '1'; reference categories
    # leave the '1' unchanged. Continuous covariates enter as power
    # functions (ratio-to-reference)^exponent per the paper equation
    # given in the Table 2 footnote.
    # =================================================================
    cov_cl <-
      (WT  / 73.5)^e_wt_cl  *
      (AGE / 68  )^e_age_cl *
      (ALB / 41.2)^e_alb_cl *
      (1 + e_sexf_cl         * SEXF)         *
      (1 + e_race_black_cl   * RACE_BLACK)   *
      (1 + e_race_asian_cl   * RACE_ASIAN)   *
      (1 + e_dis_bcellnhl_cl * DIS_BCELLNHL) *
      (1 + e_dis_wm_cl       * DIS_WM)       *
      (1 + e_dis_other_cl    * DIS_OTHER_HEME) *
      (1 + e_renmild_cl      * RENALIMP_MILD) *
      (1 + e_renmod_cl       * RENALIMP_MOD)  *
      (1 + e_hepmild_cl      * HEPIMP_MILD)   *
      (1 + e_cyp3a4indmod_cl * CONMED_CYP3A4_IND_MOD) *
      (1 + e_cyp3a4inhstr_cl * CONMED_CYP3A4_INH_STRONG)

    cov_vc <-
      (WT / 73.5)^e_wt_vc *
      (1 + e_sexf_vc         * SEXF)           *
      (1 + e_race_black_vc   * RACE_BLACK)     *
      (1 + e_race_asian_vc   * RACE_ASIAN)     *
      (1 + e_dis_bcellnhl_vc * DIS_BCELLNHL)   *
      (1 + e_dis_wm_vc       * DIS_WM)         *
      (1 + e_dis_other_vc    * DIS_OTHER_HEME)

    # Low-dose (< 30 mg) indicator uses the current DOSE regressor; both
    # the paper and NONMEM code define this on the dose value.
    dose_lt30 <- 0
    if (DOSE < 30) dose_lt30 <- 1

    cov_f <-
      (1 + e_dose_lt30_f * dose_lt30) *
      (1 + e_ppi_f       * CONMED_PPI) *
      (1 + e_h2ra_f      * CONMED_H2RA) *
      (1 + e_antacid_f   * CONMED_ANTACID)

    # PK parameters.
    cl     <- exp(lcl + etalcl) * cov_cl
    vc     <- exp(lvc + etalvc) * cov_vc
    ka     <- exp(lka)
    q      <- exp(lq)
    vp     <- exp(lvp)
    alag_d <- exp(lalag)
    fdepot <- cov_f

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot)    <- fdepot
    alag(depot) <- alag_d

    # Concentration units: dose in mg, vc in L -> mg/L = ug/mL;
    # multiply by 1000 to express Cc in ng/mL (matches NONMEM S2 = V2 / 1000).
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
