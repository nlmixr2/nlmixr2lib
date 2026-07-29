Jadhav_2023_bempedoicAcid_ldlc <- function() {
  description <- "Joint population PK + PK/PD model for bempedoic acid and serum low-density lipoprotein cholesterol (LDL-C) in patients with dyslipidemia (Jadhav 2023). The PK layer is the two-compartment / single-transit-absorption popPK model from the companion Jadhav_2023_bempedoicAcid.R file (Table 2). The PD layer (Table 3) is a type 1 indirect-response model in which bempedoic acid inhibits LDL-C production, with covariate effects on the maximum fractional inhibition Imax (sex, body weight, Black race, concomitant statin intensity, concomitant ezetimibe, prior statin therapy) and on baseline LDL-C (concomitant statin intensity, HeFH, type 2 diabetes, prior ezetimibe therapy, prior statin therapy). The authors fit the PD layer sequentially, conditioned on individual post hoc PK parameters from the popPK model."
  reference   <- "Jadhav SB, Amore BM, Bockbrader H, Crass RL, Chapel S, Sasiela WJ, Emery MG. Population pharmacokinetic and pharmacokinetic-pharmacodynamic modeling of bempedoic acid and low-density lipoprotein cholesterol in healthy subjects and patients with dyslipidemia. J Pharmacokinet Pharmacodyn. 2023;50(5):351-364. doi:10.1007/s10928-023-09864-w"
  vignette    <- "Jadhav_2023_bempedoicAcid"
  paper_specific_compartments <- c("LDL")
  units       <- list(time = "h", dosing = "mg", concentration = "ug/mL", ldlc = "mg/dL")

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power covariate on CL/F (0.61) and Vc/F (0.94) in the PK layer (Table 2) and on Imax (0.544) in the PD layer (Table 3). The paper prints the covariate form theta_TV = theta_REF * (x / x_REF)^theta_x but never the reference x_REF; the model file uses the analysis-set medians from Table 1 -- 83.7 kg for the PK layer (popPK analysis set) and 84.5 kg for the PD layer (bempedoic-acid-treated popPK/PD analysis set), following the same per-layer-reference convention as Pu_2021_evinacumab.R. Median normalisation is the standard for a pooled full-covariate popPK model with ESTIMATED (not fixed-allometric) exponents, and it reconciles the model with the paper's own reported mean steady-state exposure of 12.5 ug/mL at 180 mg/day (see the vignette Assumptions and deviations section for the reconstruction and the discarded 70 kg alternative). Only absolute typical values depend on the choice; published covariate ratios are invariant to it.",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power covariate on Vc/F (0.743) in the PK layer (Table 2). Reference age not printed; the model file uses the popPK analysis-set median of 62.0 years (Table 1), i.e. the same PK-layer reference as the companion Jadhav_2023_bempedoicAcid.R, because the PK layer was estimated on that set. Not retained as a covariate in the PD layer.",
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Proportional-shift covariate on CL/F (-0.127) and Vc/F (-0.0895) in the PK layer (Table 2) and on Imax (+0.203) in the PD layer (Table 3). Jadhav 2023 Online Resource 5: females -26.71% (90% CI -27.84, -25.81) vs. males -21.32% (90% CI -21.92, -20.64) LDL-C change from baseline.",
      source_name        = "SEXF"
    ),
    RACE_BLACK = list(
      description        = "Black / African American race indicator, 1 = Black, 0 = other race",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Black; 92.5% of the popPK/PD analysis set was White)",
      notes              = "Proportional-shift covariate on CL/F (-0.143) in the PK layer (Table 2) and on Imax (-0.240) in the PD layer (Table 3). Jadhav 2023 Online Resource 5: Black -20.93% vs. White -23.39% LDL-C change from baseline.",
      source_name        = "RACE_BLACK"
    ),
    DIS_HYPERLIP = list(
      description        = "Hyperlipidemia diagnosis indicator, 1 = participant has hyperlipidemia, 0 = no hyperlipidemia diagnosis",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no hyperlipidemia diagnosis)",
      notes              = "Proportional-shift covariate on CL/F (-0.0945) in the PK layer (Table 2). 97.7% of the bempedoic-acid-treated popPK/PD analysis set carried a hyperlipidemia diagnosis (Table 1). Not retained in the PD layer.",
      source_name        = "Hyperlipidemia"
    ),
    DIS_DIAB = list(
      description        = "Type 2 diabetes mellitus indicator, 1 = participant has T2DM, 0 = no diabetes",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no diabetes)",
      notes              = "Proportional-shift covariate on CL/F (-0.177) in the PK layer (Table 2) and on baseline LDL-C (-0.0661) in the PD layer (Table 3, row 'T2DM'). 21.1% of the bempedoic-acid-treated popPK/PD analysis set had diabetes (Table 1).",
      source_name        = "T2DM"
    ),
    DIS_HEFH = list(
      description        = "Heterozygous familial hypercholesterolemia indicator, 1 = patient has HeFH, 0 = non-HeFH",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-HeFH)",
      notes              = "Proportional-shift covariate on baseline LDL-C (+0.0671) in the PD layer (Table 3). HeFH patients have higher baseline LDL-C than the non-HeFH reference; Online Resource 5 reports steady-state LDL-C 97.38 mg/dL for HeFH vs. 88.53 mg/dL for non-HeFH. Two phase 3 studies enrolled patients with prior ASCVD and/or HeFH on maximally tolerated statin therapy.",
      source_name        = "HeFH"
    ),
    CRCL = list(
      description        = "MDRD-estimated glomerular filtration rate, expressed in absolute mL/min WITHOUT body-surface-area adjustment",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power covariate on CL/F (0.574) in the PK layer (Table 2, row 'eGFR'). Table 1 footnote a states the values are absolute mL/min with no BSA adjustment, so raw mL/min must be supplied. Reference 89.3 mL/min, the popPK analysis-set median (Table 1), i.e. the same PK-layer reference as the companion Jadhav_2023_bempedoicAcid.R because the PK layer was estimated on that set; it sits essentially on the paper's own normal-renal-function threshold of 90 mL/min (Fig. 2). Not retained in the PD layer.",
      source_name        = "eGFR"
    ),
    CONMED_EZE = list(
      description        = "Concomitant ezetimibe indicator, 1 = receiving concomitant ezetimibe, 0 = not receiving ezetimibe",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant ezetimibe)",
      notes              = "Proportional-shift covariate on CL/F (-0.0934) in the PK layer (Table 2) and on Imax (+0.190) in the PD layer (Table 3, row 'Ezetimibe' under 'Covariates of Imax'). Ezetimibe blocks intestinal cholesterol absorption through a pathway independent of ATP-citrate lyase inhibition, so the combination increases the maximum achievable LDL-C reduction: Online Resource 5 reports -29.44% vs. -22.37% LDL-C change from baseline. Distinct from PRIOR_EZE (pre-study ezetimibe therapy), which acts on baseline LDL-C.",
      source_name        = "Ezetimibe"
    ),
    CONMED_SIMVASTATIN = list(
      description        = "Concomitant simvastatin indicator, 1 = receiving concomitant simvastatin, 0 = not receiving simvastatin",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant simvastatin)",
      notes              = "Proportional-shift covariate on Vc/F (-0.154) in the PK layer (Table 2). Not retained in the PD layer, where statin effects enter through the three intensity strata instead.",
      source_name        = "Simvastatin"
    ),
    CONMED_ATORVASTATIN = list(
      description        = "Concomitant atorvastatin indicator, 1 = receiving concomitant atorvastatin, 0 = not receiving atorvastatin",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant atorvastatin; the F1 = 1 anchor)",
      notes              = "Proportional-shift covariate on relative oral bioavailability F1 (+0.142) in the PK layer (Table 2, footnote b). Not retained in the PD layer.",
      source_name        = "Atorvastatin"
    ),
    CONMED_STATIN_LI = list(
      description        = "Concomitant low-intensity statin indicator, 1 = receiving a low-intensity statin regimen, 0 = not on a low-intensity statin",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant statin of this intensity; the model reference is no concomitant statin at all, i.e. all three intensity indicators = 0)",
      notes              = "Proportional-shift covariate on Imax (-0.238) and on baseline LDL-C (-0.159) in the PD layer (Table 3). The three intensity indicators are mutually exclusive: at most one may be 1 for a given patient. Jadhav 2023 Results: low-intensity statin + bempedoic acid gives -23.7% (90% CI -26.3, -21.4) LDL-C change from baseline vs. -30.5% (90% CI -31.6, -29.6) for bempedoic acid without concomitant statin.",
      source_name        = "Low-intensity statin"
    ),
    CONMED_STATIN_MI = list(
      description        = "Concomitant moderate-intensity statin indicator, 1 = receiving a moderate-intensity statin regimen, 0 = not on a moderate-intensity statin",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant statin of this intensity)",
      notes              = "Proportional-shift covariate on Imax (-0.302) and on baseline LDL-C (-0.268) in the PD layer (Table 3). Mutually exclusive with CONMED_STATIN_LI and CONMED_STATIN_HI. Jadhav 2023 Results: -21.8% (90% CI -22.6, -20.9) LDL-C change from baseline.",
      source_name        = "Moderate-intensity statin"
    ),
    CONMED_STATIN_HI = list(
      description        = "Concomitant high-intensity statin indicator, 1 = receiving a high-intensity statin regimen, 0 = not on a high-intensity statin",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant statin of this intensity)",
      notes              = "Proportional-shift covariate on Imax (-0.424) and on baseline LDL-C (-0.293) in the PD layer (Table 3). Mutually exclusive with CONMED_STATIN_LI and CONMED_STATIN_MI. Jadhav 2023 Results: -18.0% (90% CI -18.7, -17.4) LDL-C change from baseline; the larger baseline-LDL-C reduction offsets the smaller fractional change so that absolute steady-state LDL-C is similar to the moderate-intensity group (Online Resource 5: 83.99 vs. 83.18 mg/dL).",
      source_name        = "High-intensity statin"
    ),
    PRIOR_STATIN = list(
      description        = "Prior statin therapy indicator, 1 = patient was on established statin therapy before study entry, 0 = statin-naive at study entry",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no prior statin therapy)",
      notes              = "Proportional-shift covariate on Imax (-0.373) and on baseline LDL-C (-0.296) in the PD layer (Table 3, rows 'Statin prior therapy'). Distinct from the concomitant-intensity indicators: PD covariate screening included 'prior established LMTs' separately from 'concomitant medication (low-, moderate-, or high-intensity statin or ezetimibe)' (Methods, PD model covariate analysis). Online Resource 5: prior-statin patients show -18.83% (90% CI -24.2, -15.34) LDL-C change from baseline vs. -23.32% for statin-naive.",
      source_name        = "Statin prior therapy"
    ),
    PRIOR_EZE = list(
      description        = "Prior ezetimibe therapy indicator, 1 = patient was on established ezetimibe therapy before study entry, 0 = ezetimibe-naive at study entry",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no prior ezetimibe therapy)",
      notes              = "Proportional-shift covariate on baseline LDL-C (-0.0596) in the PD layer (Table 3, row 'Ezetimibe prior therapy'). Distinct from CONMED_EZE (concomitant ezetimibe), which acts on Imax and on CL/F. Online Resource 5 reports both a 'Prior Treatment (Ezetimibe : No Ezetimibe)' and a 'Concomitant Treatment (Ezetimibe : No Ezetimibe)' contrast, confirming they are separate columns in the analysis dataset.",
      source_name        = "Ezetimibe prior therapy"
    ),
    FED = list(
      description        = "Fed-state dose-record indicator, 1 = dose administered with food, 0 = fasted",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted)",
      notes              = "Proportional-shift covariate on Ka (-0.777) in the PK layer (Table 2). Per-dose-record indicator; food slows absorption by 78% without changing the extent of absorption (Discussion).",
      source_name        = "Food"
    ),
    SAMPLE_INTENSIVE = list(
      description        = "Per-observation sampling-intensity indicator, 1 = serial (intensive) PK sampling, 0 = sparse PK sampling",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (sparse PK sampling)",
      notes              = "Record-level indicator that switches the PK proportional residual-error magnitude between 31.9% (serial) and 54.3% (sparse), Jadhav 2023 Table 2. Applies only to the Cc endpoint; the LDL-C endpoint has a single additive-plus-proportional residual (Table 3).",
      source_name        = "SAMPLE_INTENSIVE"
    )
  )

  population <- list(
    species                = "human",
    n_subjects             = 4459L,
    n_subjects_treated     = 2984L,
    n_subjects_placebo     = 1475L,
    n_observations_ldlc    = 27534L,
    n_studies              = 15L,
    phases                 = "Three phase 1, eight phase 2, and four phase 3 studies (a subset of the 22 popPK studies; Online Resource 1)",
    age_range              = "21.0-91.0 years (bempedoic acid arm); 20.0-88.0 years (placebo arm)",
    age_median             = "64.0 years (bempedoic acid arm, mean 63.3, SD 10.2)",
    weight_range           = "42.5-160 kg (bempedoic acid arm)",
    weight_median          = "84.5 kg (bempedoic acid arm, mean 86.0, SD 17.3)",
    sex_female_pct         = 36.4,
    race_ethnicity         = c(White = 92.5, Black = 5.6, Asian = 1.1, NativeAmerican = 0.1, NativeHawaiian = 0.2, Other = 0.5),
    disease_state          = "Patients with hyperlipidemia (97.7%), including cohorts with prior ASCVD and/or heterozygous familial hypercholesterolemia on maximally tolerated statin therapy and cohorts with statin intolerance; 21.1% had type 2 diabetes mellitus. A single phase 1 study (Study 04) contributed 24 healthy subjects.",
    ldlc_baseline          = "Observed median baseline LDL-C 113 mg/dL (mean 122, SD 39.2, range 48.0-422) in the bempedoic acid arm and 110 mg/dL (mean 118, SD 38.2, range 38.0-411) in the placebo arm. Both arms include patients on stable background lipid-modifying therapy at Day 1 and patients with no ongoing LMT.",
    dose_range             = "Oral bempedoic acid 60-240 mg once daily (commercial regimen 180 mg once daily), or matching placebo. Placebo-treated participants entered the popPK/PD model with a bempedoic acid concentration of 0 ug/mL.",
    regions                = "Multi-regional pool of 15 Esperion-sponsored clinical studies.",
    notes                  = "Demographics from Jadhav 2023 Table 1 (popPK/PD columns). The PD layer was fit sequentially, conditioned on individual post hoc PK parameters from the popPK model; for the 989 (33%) bempedoic-acid-treated patients without measurable plasma concentrations, population-predicted concentrations were used. Patients and treatment periods with < 80% compliance were excluded. LDL-C data from the phase 3 studies are minimally informative about the dynamic onset of response because the first post-baseline LDL-C sample was collected on Day 29, near steady state."
  )

  ini({
    # =========================================================================
    # PK layer -- Jadhav 2023 Table 2 (Final PopPK model parameters)
    #
    # Reference participant: male, non-Black, no hyperlipidemia, no T2DM, no
    # concomitant ezetimibe / simvastatin / atorvastatin, fasted dosing,
    # WT = 83.7 kg, AGE = 62 years, eGFR = 89.3 mL/min.
    #
    # Covariate model forms (Jadhav 2023 Methods, "PK model covariate analysis";
    # the PD Methods states the same forms are used for the popPK/PD parameters):
    #   continuous : theta_TV,ij = theta_REF * (x_ij / x_REF)^theta_x
    #   categorical: theta_TV,ij = theta_REF * (1 + theta_x * x_ij)
    # =========================================================================
    lcl  <- log(0.755); label("Apparent systemic clearance CL/F at the reference condition (L/h)")               # Table 2: CL/F = 0.755 L/h (%RSE 2.6)
    lvc  <- log(19.1);  label("Apparent central distribution volume Vc/F at the reference condition (L)")         # Table 2: Vc/F = 19.1 L (%RSE 6.9)
    lka  <- log(1.41);  label("Absorption / transit rate constant Ka under fasted dosing (1/h)")                  # Table 2: Ka = 1.41 1/h (%RSE 5.6)
    lk12 <- log(0.184); label("Central-to-peripheral distribution rate constant K23 (1/h)")                       # Table 2: K23 = 0.184 1/h (%RSE 8.3)
    lk21 <- log(0.156); label("Peripheral-to-central distribution rate constant K32 (1/h)")                       # Table 2: K32 = 0.156 1/h (%RSE 4.3)
    lfdepot <- fixed(log(1)); label("Relative oral bioavailability F1 in the no-atorvastatin reference (fraction; FIXED at 1)")  # Fig. 1 caption: F1 typical value fixed at 1

    e_wt_cl              <- 0.61;    label("Power exponent of (WT / 83.7 kg) on CL/F (unitless)")                                   # Table 2: Body weight = 0.61 (%RSE 9.8)
    e_crcl_cl            <- 0.574;   label("Power exponent of (CRCL / 89.3 mL/min) on CL/F (unitless)")                             # Table 2: eGFR = 0.574 (%RSE 6.3)
    e_sexf_cl            <- -0.127;  label("Proportional shift of female sex on CL/F (unitless)")                                 # Table 2: Female sex = -0.127 (%RSE 16.9)
    e_race_black_cl      <- -0.143;  label("Proportional shift of Black race on CL/F (unitless)")                                 # Table 2: Black race = -0.143 (%RSE 18.0)
    e_dis_hyperlip_cl    <- -0.0945; label("Proportional shift of hyperlipidemia on CL/F (unitless)")                             # Table 2: Hyperlipidemia = -0.0945 (%RSE 28.4)
    e_dis_diab_cl        <- -0.177;  label("Proportional shift of type 2 diabetes mellitus on CL/F (unitless)")                   # Table 2: T2DM = -0.177 (%RSE 16.7)
    e_conmed_eze_cl      <- -0.0934; label("Proportional shift of concomitant ezetimibe on CL/F (unitless)")                      # Table 2: Ezetimibe = -0.0934 (%RSE 28.4)

    e_wt_vc              <- 0.94;    label("Power exponent of (WT / 83.7 kg) on Vc/F (unitless)")                                   # Table 2: Body weight = 0.94 (%RSE 20.8)
    e_age_vc             <- 0.743;   label("Power exponent of (AGE / 62 years) on Vc/F (unitless)")                               # Table 2: Age = 0.743 (%RSE 15.9)
    e_sexf_vc            <- -0.0895; label("Proportional shift of female sex on Vc/F (unitless)")                                 # Table 2: Female sex = -0.0895 (%RSE 33.4)
    e_conmed_simvastatin_vc <- -0.154; label("Proportional shift of concomitant simvastatin on Vc/F (unitless)")                  # Table 2: Simvastatin = -0.154 (%RSE 29.2)

    e_fed_ka             <- -0.777;  label("Proportional shift of fed-state dosing on Ka (unitless)")                             # Table 2: Food on Ka = -0.777 (%RSE 1.4)
    e_conmed_atorvastatin_fdepot <- 0.142; label("Proportional shift of concomitant atorvastatin on relative oral bioavailability F1 (unitless)")  # Table 2: Atorvastatin = 0.142 (%RSE 19.1)

    # PK inter-individual variability (Table 2, "IIV, %CV"); omega^2 = log(CV^2 + 1)
    etalcl ~ 0.084533   # Table 2: IIV CL/F = 29.7 %CV -> log(1 + 0.297^2)
    etalvc ~ 0.693147   # Table 2: IIV Vc/F = 100  %CV -> log(1 + 1.000^2)
    etalka ~ 0.435749   # Table 2: IIV Ka   = 73.9 %CV -> log(1 + 0.739^2)

    # PK residual error (Table 2, "Residual error, %"); log-additive in NONMEM
    # is proportional in nlmixr2's linear space.
    propSdSerial <- 0.319; label("Proportional PK residual error, serial (intensive) sampling (fraction)")   # Table 2: Serial PK sampling = 31.9% (%RSE 1.1)
    propSdSparse <- 0.543; label("Proportional PK residual error, sparse sampling (fraction)")               # Table 2: Sparse PK sampling = 54.3% (%RSE 1.3)

    # =========================================================================
    # PD layer -- Jadhav 2023 Table 3 (Final popPK/PD model parameters)
    #
    # Reference patient: male, non-Black, WT = 84.5 kg, no concomitant statin of
    # any intensity, no concomitant ezetimibe, no prior statin therapy, no prior
    # ezetimibe therapy, non-HeFH, no T2DM.
    #
    # Type 1 indirect-response model (Methods, Eq. for dLDLC/dt):
    #   dLDLC/dt = kin * [1 - Imax * C / (IC50 + C)] - kout * LDLC
    # with kin = kout * baseline LDL-C at drug-free steady state.
    # =========================================================================
    limax  <- log(0.350); label("Log of the maximum fractional inhibition of LDL-C production Imax at the reference condition (unitless)")  # Table 3: Imax = 0.350 (%RSE 5.20; 95% CI 0.314, 0.386)
    lic50  <- log(3.17);  label("Bempedoic acid concentration producing 50% of Imax, IC50 (ug/mL)")                                          # Table 3: IC50 = 3.17 ug/mL (%RSE 17.7; 95% CI 2.07, 4.26)
    lrbase <- log(143);   label("Typical baseline serum LDL-C at the reference condition (mg/dL)")                                           # Table 3: Baseline LDL-C = 143 mg/dL (%RSE 0.800; 95% CI 141, 145)
    lkout  <- log(1 / 85.8); label("Log first-order LDL-C elimination rate constant kout (1/h; = 1 / TURN)")                                 # Table 3: TURN = 85.8 h (%RSE 9.30; 95% CI 70.2, 101). Methods: "LDL-C turnover represented as the inverse of kout (1/kout)". kout = 1/85.8 h = 0.2795 1/day, matching the Discussion's "kout (0.3 day-1)".

    # ---- Covariate effects on Imax (Table 3, "Covariates of Imax") ----
    e_conmed_eze_imax        <-  0.190; label("Proportional shift of concomitant ezetimibe on Imax (unitless)")                    # Table 3: Ezetimibe = 0.190 (%RSE 20.3; 95% CI 0.114, 0.266)
    e_conmed_statin_li_imax  <- -0.238; label("Proportional shift of concomitant low-intensity statin on Imax (unitless)")         # Table 3: Low-intensity statin = -0.238 (%RSE 18.8; 95% CI -0.325, -0.150)
    e_conmed_statin_mi_imax  <- -0.302; label("Proportional shift of concomitant moderate-intensity statin on Imax (unitless)")    # Table 3: Moderate-intensity statin = -0.302 (%RSE 7.50; 95% CI -0.346, -0.257)
    e_conmed_statin_hi_imax  <- -0.424; label("Proportional shift of concomitant high-intensity statin on Imax (unitless)")        # Table 3: High-intensity statin = -0.424 (%RSE 4.50; 95% CI -0.461, -0.387)
    e_sexf_imax              <-  0.203; label("Proportional shift of female sex on Imax (unitless)")                               # Table 3: Female sex = 0.203 (%RSE 16.5; 95% CI 0.137, 0.269)
    e_wt_imax                <-  0.544; label("Power exponent of (WT / 84.5 kg) on Imax (unitless)")                                 # Table 3: Body weight = 0.544 (%RSE 12.6; 95% CI 0.410, 0.679)
    e_race_black_imax        <- -0.240; label("Proportional shift of Black race on Imax (unitless)")                               # Table 3: Black race = -0.240 (%RSE 19.1; 95% CI -0.330, -0.150)
    e_prior_statin_imax      <- -0.373; label("Proportional shift of prior statin therapy on Imax (unitless)")                     # Table 3: Statin prior therapy = -0.373 (%RSE 20.7; 95% CI -0.525, -0.222)

    # ---- Covariate effects on baseline LDL-C (Table 3, "Covariates of baseline LDL-C") ----
    e_conmed_statin_li_rbase <- -0.159;  label("Proportional shift of concomitant low-intensity statin on baseline LDL-C (unitless)")      # Table 3: Low-intensity statin = -0.159 (%RSE 10.0; 95% CI -0.191, -0.128)
    e_conmed_statin_mi_rbase <- -0.268;  label("Proportional shift of concomitant moderate-intensity statin on baseline LDL-C (unitless)") # Table 3: Moderate-intensity statin = -0.268 (%RSE 2.80; 95% CI -0.282, -0.253)
    e_conmed_statin_hi_rbase <- -0.293;  label("Proportional shift of concomitant high-intensity statin on baseline LDL-C (unitless)")     # Table 3: High-intensity statin = -0.293 (%RSE 2.40; 95% CI -0.307, -0.279)
    e_dis_hefh_rbase         <-  0.0671; label("Proportional shift of heterozygous familial hypercholesterolemia on baseline LDL-C (unitless)")  # Table 3: HeFH = 0.0671 (%RSE 27.0; 95% CI 0.0316, 0.103)
    e_dis_diab_rbase         <- -0.0661; label("Proportional shift of type 2 diabetes mellitus on baseline LDL-C (unitless)")              # Table 3: T2DM = -0.0661 (%RSE 13.0; 95% CI -0.0830, -0.0492)
    e_prior_eze_rbase        <- -0.0596; label("Proportional shift of prior ezetimibe therapy on baseline LDL-C (unitless)")               # Table 3: Ezetimibe prior therapy = -0.0596 (%RSE 27.5; 95% CI -0.0917, -0.0274)
    e_prior_statin_rbase     <- -0.296;  label("Proportional shift of prior statin therapy on baseline LDL-C (unitless)")                  # Table 3: Statin prior therapy = -0.296 (%RSE 8.40; 95% CI -0.345, -0.247)

    # ---- PD inter-individual variability (Table 3, "IIV, %CV") ----
    # Only diagonal %CV values are published. The distributional form is not
    # restated in the PD Methods; the log-normal assumption declared for the
    # popPK analysis is carried forward (see vignette Assumptions and deviations).
    etalimax  ~ 0.170385  # Table 3: IIV Imax           = 43.1 %CV -> log(1 + 0.431^2) = 0.170385
    etalrbase ~ 0.055549  # Table 3: IIV Baseline LDL-C = 23.9 %CV -> log(1 + 0.239^2) = 0.055549

    # ---- PD residual error (Table 3, "Residual error") ----
    # Table 3 reports the additive residual as "3.94 g/dL"; the LDL-C scale
    # throughout the paper is mg/dL (baseline 143 mg/dL), so this is a units
    # typographical error in the source and the value is used as 3.94 mg/dL.
    propSd_LDL <- 0.153; label("Proportional residual error on LDL-C (fraction)")   # Table 3: PD proportional = 15.3% (%RSE 0.900; 95% CI 15.0, 15.6)
    addSd_LDL  <- 3.94;  label("Additive residual error on LDL-C (mg/dL)")          # Table 3: PD additive = 3.94 (%RSE 9.10; 95% CI 3.24, 4.64); printed as "g/dL" in Table 3, read as mg/dL
  })

  model({
    # ---- PK layer: individual parameters (Table 2 reference condition) ----
    cl <- exp(lcl + etalcl) *
      (WT / 83.7)^e_wt_cl *
      (CRCL / 89.3)^e_crcl_cl *
      (1 + e_sexf_cl * SEXF) *
      (1 + e_race_black_cl * RACE_BLACK) *
      (1 + e_dis_hyperlip_cl * DIS_HYPERLIP) *
      (1 + e_dis_diab_cl * DIS_DIAB) *
      (1 + e_conmed_eze_cl * CONMED_EZE)

    vc <- exp(lvc + etalvc) *
      (WT / 83.7)^e_wt_vc *
      (AGE / 62)^e_age_vc *
      (1 + e_sexf_vc * SEXF) *
      (1 + e_conmed_simvastatin_vc * CONMED_SIMVASTATIN)

    ka  <- exp(lka + etalka) * (1 + e_fed_ka * FED)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    kel <- cl / vc

    fdepot <- exp(lfdepot) * (1 + e_conmed_atorvastatin_fdepot * CONMED_ATORVASTATIN)

    # ---- PK ODE system (Jadhav 2023 Fig. 1) ----
    d/dt(depot)       <- -ka * depot
    d/dt(transit1)    <-  ka * depot - ka * transit1
    d/dt(central)     <-  ka * transit1 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot) <- fdepot

    # Bempedoic acid plasma concentration in ug/mL (dose mg / volume L).
    Cc <- central / vc

    # ---- PD layer: individual parameters (Table 3 reference condition) ----
    imax <- exp(limax + etalimax) *
      (WT / 84.5)^e_wt_imax *
      (1 + e_sexf_imax * SEXF) *
      (1 + e_race_black_imax * RACE_BLACK) *
      (1 + e_conmed_eze_imax * CONMED_EZE) *
      (1 + e_conmed_statin_li_imax * CONMED_STATIN_LI) *
      (1 + e_conmed_statin_mi_imax * CONMED_STATIN_MI) *
      (1 + e_conmed_statin_hi_imax * CONMED_STATIN_HI) *
      (1 + e_prior_statin_imax * PRIOR_STATIN)

    ic50 <- exp(lic50)

    rbase <- exp(lrbase + etalrbase) *
      (1 + e_conmed_statin_li_rbase * CONMED_STATIN_LI) *
      (1 + e_conmed_statin_mi_rbase * CONMED_STATIN_MI) *
      (1 + e_conmed_statin_hi_rbase * CONMED_STATIN_HI) *
      (1 + e_dis_hefh_rbase * DIS_HEFH) *
      (1 + e_dis_diab_rbase * DIS_DIAB) *
      (1 + e_prior_eze_rbase * PRIOR_EZE) *
      (1 + e_prior_statin_rbase * PRIOR_STATIN)

    # Methods: "At baseline prior to drug exposure, steady-state conditions
    # prevail (dLDLC/dt = 0) and kin is determined as the product of kout and
    # baseline LDL-C."
    kout <- exp(lkout)
    kin  <- kout * rbase

    # ---- Type 1 indirect-response ODE for LDL-C (Jadhav 2023 Fig. 3) ----
    LDL(0)   <- rbase
    d/dt(LDL) <- kin * (1 - imax * Cc / (ic50 + Cc)) - kout * LDL

    # ---- Observation and residual error ----
    # Cc in ug/mL; LDL in mg/dL. The PK residual magnitude is switched per
    # observation by SAMPLE_INTENSIVE (Table 2); the LDL-C endpoint carries a
    # single combined additive + proportional residual (Table 3).
    propSd <- propSdSerial * SAMPLE_INTENSIVE + propSdSparse * (1 - SAMPLE_INTENSIVE)
    Cc  ~ prop(propSd)
    LDL ~ add(addSd_LDL) + prop(propSd_LDL)
  })
}
