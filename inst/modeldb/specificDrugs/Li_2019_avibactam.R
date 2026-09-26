Li_2019_avibactam <- function() {
  description <- paste(
    "Two-compartment IV population PK model for the avibactam component of",
    "ceftazidime-avibactam, fitted to 13,735 plasma concentrations from",
    "2,249 adults pooled across Phase 1, Phase 2 and Phase 3 studies in",
    "complicated intra-abdominal infection (cIAI), complicated urinary tract",
    "infection (cUTI) and nosocomial pneumonia (NP), including",
    "ventilator-associated pneumonia, and spanning renal function from",
    "end-stage renal disease to augmented renal clearance (Li 2019).",
    "Creatinine clearance acts on clearance through a hinged relationship",
    "-- a power function below 80 mL/min and a shallow linear function at or",
    "above it -- with separate absolute or proportional arms for end-stage",
    "renal disease and hemodialysis. Central volume carries body weight,",
    "infection type and mechanical-ventilation effects. Inter-individual",
    "variability is a full 4x4 OMEGA block across CL, Vc, Vp and Q, and",
    "residual variability is switched across three study-phase strata.",
    "Companion to Li_2019_ceftazidime; the two analytes were fitted as",
    "separate models to separate data sets and are combined only in the",
    "joint PK/PD target attainment analysis."
  )
  reference <- paste(
    "Li J, Lovern M, Green ML, Chiu J, Zhou D, Comisar C, Xiong Y, Hing J,",
    "MacPherson M, Wright JG, Riccobene T, Carrothers TJ, Das S.",
    "Ceftazidime-Avibactam Population Pharmacokinetic Modeling and",
    "Pharmacodynamic Target Attainment Across Adult Indications and Patient",
    "Subgroups. Clin Transl Sci. 2019;12(2):151-163. doi:10.1111/cts.12585.",
    "Fixed-effect, random-effect and residual-error estimates are Table 2.",
    "The functional FORMS of every covariate relationship are the final",
    "avibactam NONMEM control stream printed in supplementary Data S1",
    "(CTS-12-151-s001, 'NONMEM control file for avibactam final PopPK",
    "model'); that control stream also settles the end-stage-renal-disease",
    "term, which Table 2 and the Results text mislabel as an absolute",
    "clearance -- see the ini() comment on e_renalimp_esrd_cl.",
    sep = " "
  )
  vignette <- "Li_2019_ceftazidime_avibactam"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    central = list(analyte = "avibactam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "avibactam", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance estimated with the Cockcroft-Gault equation from chronological serum creatinine records (Li 2019 Methods, 'Analysis data and model construction')",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Raw Cockcroft-Gault mL/min, NOT BSA-normalized; same convention as the companion Li_2019_ceftazidime.R and the sibling pair Chen_2025_ceftazidime.R / Chen_2025_avibactam.R. Observed range in the avibactam data set 11-610 mL/min (Li 2019 Methods), which is wider than the ceftazidime data set because the avibactam studies included dedicated renal-impairment and augmented-renal-clearance cohorts. The effect on CL is a hinge at 80 mL/min -- a power function below and a shallow linear function at or above -- and BOTH arms equal 1 at the hinge, so the ini() clearance is the typical value at CRCL = 80 mL/min. This differs from the companion ceftazidime model, whose two-segment linear spline is hinged at 100 mL/min and does not pass through 1 there. The hinge is bypassed entirely for RENALIMP_ESRD subjects, and the sub-80 power arm is additionally bypassed for RENAL_ARC subjects.",
      source_name = "CLCR"
    ),
    WT = list(
      description = "Body weight at baseline (Li 2019 Table 2, theta14 'WT on Vc (WT/70.0)^theta14')",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Centred at 70 kg. The avibactam control stream's $PROBLEM line records the choice explicitly -- 'centre weight at 70 kg rather than 71.2 kg' -- so 70 is a deliberate round reference rather than the cohort median. Li 2019 Results, 'Avibactam' reports the 10th (51 kg) and 90th (95 kg) weight percentiles as giving 29% lower and 39% higher Vc, which reproduce from (WT/70)^1.08. Weight acts on Vc ONLY: the control stream carries a CLWT term but fixes its exponent theta30 to zero ('CL~WT removed THETA(30) FIXED TO ZERO'), so the final model has no allometric clearance term. Subjects with a missing weight were imputed to 70 kg in the source analysis.",
      source_name = "WT"
    ),
    DIS_CIAI = list(
      description = "Complicated intra-abdominal infection infection-type indicator (Li 2019 control stream POP = 3)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = not cIAI; the shared all-zero reference across the DIS_* set is the Phase 1 subject (healthy volunteer or renal-impairment cohort)",
      notes = "786 of 2,249 subjects (34.9%), Li 2019 Results, 'Analysis populations'. The Vc effect is PHASE-DEPENDENT: Phase 2 cIAI subjects take theta9 = 1.92 and Phase 3 cIAI subjects take theta12 = 0.329, and the two are alternatives rather than a stack (control stream: 'IF(POP.EQ.3.AND.PHASE.EQ.2) V1POP=(1 + THETA(9))' / 'IF(POP.EQ.3.AND.PHASE.EQ.3) V1POP=(1 + THETA(12))'). The Phase 2 cohort additionally takes the only cIAI CLEARANCE effect in the model (theta10 = 0.406); Phase 3 cIAI subjects carry no clearance shift. Use STUDY_CIAI_PH2 to select between the two.",
      source_name = "POP = 3"
    ),
    DIS_CUTI = list(
      description = "Complicated urinary tract infection infection-type indicator (Li 2019 control stream POP = 2)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = not cUTI; shared all-zero reference is the Phase 1 subject",
      notes = "705 of 2,249 subjects (31.3%), Li 2019 Results, 'Analysis populations'. Takes the SAME Vc shift (theta11 = 0.434) in Phase 2 and Phase 3 -- the control stream writes the two branches separately but assigns theta11 to both, with the inline comment 'lumped with TH11 now'. Carries no clearance effect. Unlike the companion ceftazidime model, the avibactam model carries NO separate acute-pyelonephritis term, so DIS_AP is not declared here.",
      source_name = "POP = 2"
    ),
    DIS_HABP = list(
      description = "Hospital-acquired bacterial pneumonia infection-type indicator (component of Li 2019 control stream POP = 4, 'NP')",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = not HABP; shared all-zero reference is the Phase 1 subject",
      notes = "Li 2019 models nosocomial pneumonia as a SINGLE level (POP = 4) covering both hospital-acquired and ventilator-associated pneumonia, so DIS_HABP and DIS_VABP share one coefficient and enter model() as their sum, pooled further with Phase 3 cIAI onto theta12. The columns are kept distinct per the DIS_VABP register entry because successor models in this lineage separate them. 413 of 2,249 subjects (18.4%) had NP of either kind; 138 of those had VAP (Li 2019 Table 3). Whether a given subject was ventilated ON THE PK SAMPLING DAY is carried separately by MECH_VENT.",
      source_name = "POP = 4"
    ),
    DIS_VABP = list(
      description = "Ventilator-associated bacterial pneumonia infection-type indicator (component of Li 2019 control stream POP = 4, 'NP')",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = not VABP; shared all-zero reference is the Phase 1 subject",
      notes = "See DIS_HABP: Li 2019 carries one pooled NP coefficient applied to DIS_HABP + DIS_VABP.",
      source_name = "POP = 4"
    ),
    STUDY_CIAI_PH2 = list(
      description = "Phase 2 complicated-intra-abdominal-infection study cohort indicator (Li 2019 control stream POP = 3 AND PHASE = 2)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = any other study, including the Phase 3 cIAI trials",
      notes = "Li 2019 is the founding paper of the lineage this canonical was registered for: its Table 2 rows 'Population effect on Vc (cIAI, phase II)' (theta9) and 'Population effect on CL (cIAI, phase II)' (theta10) are the register entry's quoted Das 2024 source aliases, which Xie 2025 later labels 'Study2002'. IMPORTANT DIFFERENCE FROM Xie_2025_aztreonam_avibactam.R: there the indicator STACKS on top of DIS_CIAI, whereas in Li 2019 the Phase 2 cIAI volume shift REPLACES the Phase 3 one (theta9 instead of theta12, not in addition to it). The clearance shift has no Phase 3 counterpart, so it is a plain addition. By construction this column equals DIS_CIAI * STUDY_CAZAVI_PHASE2.",
      source_name = "POP = 3 AND PHASE = 2"
    ),
    RACE_ASIAN_OTH = list(
      description = "Non-Chinese, non-Japanese Asian race indicator (Li 2019 control stream RCE = 3, abbreviated ASN)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = all other race levels, which the control stream's CLRCE block leaves at a factor of 1",
      notes = "The dominant reference grouping is the pooled non-Asian population. Unlike the companion ceftazidime model, the avibactam model retains only this ONE race level: the control stream fixes the Black, Japanese, Chinese and Other race coefficients (theta21, theta23, theta24, theta25) to zero and carries no race effect on Vc at all. Li 2019 Results, 'Avibactam' reports the retained effect as 'an estimated 8.65% lower CL (translating to a 9.5% increase in AUC) for non-Chinese, non-Japanese Asians'. Do NOT borrow the ceftazidime Chinese or Asian-volume coefficients for avibactam.",
      source_name = "RCE = 3"
    ),
    APACHE_II_SEV = list(
      description = "Elevated-APACHE-II severity stratum indicator (Li 2019 control stream APACHE = 2)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = APACHE II score <= 10, or score missing (the control stream assigns a factor of 1 to both)",
      notes = "The elevated stratum is APACHE II > 10. Li 2019 states the threshold explicitly, which resolves the gap flagged in the APACHE_II_SEV register entry and in the founding model Xie_2025_aztreonam_avibactam.R, where the threshold was unstated in every available source: Methods, 'Selection of covariates' lists 'Acute Physiology and Chronic Health Evaluation version II (APACHE II > 10)' among the markers of systemic disturbance and adds that 'predicted mortality rises steeply for scores > 10 (> 10% mortality), and this represents a reasonable cutoff for defining more severely ill patients'. Missing scores take the reference factor: APACHE II was collected for cIAI and NP patients only, so it is missing for all 648 cUTI patients (Li 2019 Table 3 footnote b). 438 subjects scored > 10 and 677 scored <= 10. Effect is -0.197 on CL, the largest covariate effect on avibactam clearance aside from renal function.",
      source_name = "APACHE"
    ),
    RENALIMP_ESRD = list(
      description = "End-stage renal disease indicator (Li 2019 control stream ESRD = 1)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = any renal function above ESRD",
      notes = "ESRD subjects came from the Phase 1 renal-impairment studies; the avibactam data set spans 'normal renal function to end-stage renal disease' (Li 2019 Methods). The indicator REPLACES the continuous CRCL hinge rather than multiplying it -- the control stream guards both CRCL arms with 'ESRD.EQ.0', so an ESRD subject takes a renal factor of exactly 1 and the ESRD clearance term instead. This is the same 'and not ESRD' discipline the register entry records for Das 2024 and Xie 2025. Mutually exclusive in effect with RRT_HEMODIAL_STATUS, which substitutes a separate absolute clearance.",
      source_name = "ESRD"
    ),
    RRT_HEMODIAL_STATUS = list(
      description = "Intermittent-hemodialysis treatment-status indicator (Li 2019 control stream DIAL = 1)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = not receiving hemodialysis",
      notes = "Selects an ABSOLUTE clearance of 20.8 L/h (Table 2 theta6) in place of the reference-clearance construction, matching the on-dialysis extracorporeal clearance of a small, highly water-soluble, minimally protein-bound molecule. Li 2019 Results, 'Avibactam': 'For patients with ESRD, CL was 0.0678 L/h off dialysis and 20.8 L/h on-dialysis.' The successor Xie_2025_aztreonam_avibactam.R carries the same construction with a dialysis clearance of 17.9 L/h. Note that the PTA simulations for ESRD 'did not account for drug removal through hemodialysis', so the paper's own Table 4 ESRD exposures are an explicit worst case.",
      source_name = "DIAL"
    ),
    RENAL_ARC = list(
      description = "Augmented renal clearance indicator (Li 2019 control stream ARC = 1)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = not augmented; every subject outside study CXL-PK-04 is classified non-ARC",
      notes = "Li 2019 Methods, 'Analysis data and model construction' defines ARC as a MEASURED creatinine clearance >= 140 mL/min by 8-hour urine collection, and states the flag is 'specific to study CXL-PK-04 (Table S1)', with 'subjects in other studies classified as non-ARC'. The definition therefore rests on a measured rather than a Cockcroft-Gault-estimated clearance, and the column is NOT derivable by thresholding CRCL. Structurally the flag does two things (control stream CLCLCR block): it scales the supra-80 CRCL slope by theta13 = 0.992, and it EXCLUDES the subject from the sub-80 power arm, so an ARC subject who nonetheless had CRCL < 80 would take a renal factor of exactly 1. The scaling itself is all but null (a 0.8% change in the slope), and the flag is 0 for every Phase 3 patient and for the whole simulated population of the paper, so it is inert in all reproduced published values.",
      source_name = "ARC"
    ),
    MECH_VENT = list(
      description = "Presence of a ventilator in the hospital room on the day of PK sampling (Li 2019 control stream POP5 = 1, abbreviated NPv)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = no ventilator present on the PK sampling day",
      notes = "Same definition as in the companion Li_2019_ceftazidime.R: 'the presence of a ventilator in the hospital room on the day of PK sampling, which includes patients with VAP or HAP who were ventilated on the day of sampling' (Li 2019 Methods, 'Selection of covariates'). Time-fixed at the PK sampling day, and NOT the same contrast as VABP vs HABP. Raises avibactam Vc by 17.5% (Table 2 theta28) against 29.7% for ceftazidime; it is the least precisely estimated avibactam fixed effect (RSE 53.3%).",
      source_name = "POP5 = 1"
    ),
    STUDY_CAZAVI_PHASE2 = list(
      description = "Phase 2 stratum indicator of the Li 2019 pooled ceftazidime-avibactam analysis (residual-error switch)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 with STUDY_CAZAVI_PHASE3 also 0 selects the Phase 1 stratum",
      notes = "The avibactam residual-error model keeps THREE distinct strata (control stream $ERROR: PH1 / PH2 / PH3 flags, with PH2 set by 'STDY.EQ.2001.OR.STDY.EQ.2002'), unlike the companion ceftazidime model, which pools Phase 2 with Phase 3. Only the Phase 1 stratum carries an additive component.",
      source_name = "PHASE = 2"
    ),
    STUDY_CAZAVI_PHASE3 = list(
      description = "Phase 3 stratum indicator of the Li 2019 pooled ceftazidime-avibactam analysis (residual-error switch)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 with STUDY_CAZAVI_PHASE2 also 0 selects the Phase 1 stratum",
      notes = "Five Phase 3 trials: RECLAIM 1/2, RECLAIM 3, RECAPTURE 1/2, REPRISE and REPROVE (Li 2019 Methods). See STUDY_CAZAVI_PHASE2.",
      source_name = "PHASE = 3"
    )
  )

  # Screened during stepwise covariate selection but NOT retained in the
  # final avibactam model, or carried in the control stream with their
  # coefficient FIXED TO ZERO (Li 2019 Methods, 'Selection of covariates';
  # supplementary Data S1 header comments). Documented here so the paper's
  # covariate screen is preserved without declaring covariates model()
  # never references.
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = male",
      notes = "Screened on both CL and Vc. Both coefficients are carried in the control stream and FIXED TO ZERO (theta26 CLSEX1, theta29 V1SEX1), with the header comment 'Remove effect of SEX from V1 THETA(29)'.",
      source_name = "SEX"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened as a power effect on Vc centred at 47 years; removed at backward elimination ('V1~AGE removed THETA(27) FIXED TO ZERO'). Li 2019 Results: age-related changes in exposure were adequately captured by changes in CrCL.",
      source_name = "AGE"
    ),
    RACE_CHINESE = list(
      description = "Chinese-heritage race indicator",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = non-Chinese",
      notes = "Screened on avibactam CL and fixed to zero (theta24), with the control-stream comment 'clinically insignificant: CL~BLACK, CL~JAPANESE, CL~CHINESE, CL~OTHER, CL~SEX'. The companion ceftazidime model DOES retain a Chinese clearance effect, so the asymmetry is a real feature of the pair.",
      source_name = "RCE = 14"
    ),
    RACE_JAPANESE = list(
      description = "Japanese-heritage race indicator",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = non-Japanese",
      notes = "Screened on avibactam CL and fixed to zero (theta23). See RACE_CHINESE.",
      source_name = "RCE = 13"
    ),
    CONMED_OAT_INHIBITOR = list(
      description = "Concomitant organic anion transporter 1/3 inhibitor (probenecid, cimetidine or diclofenac) coadministration indicator",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = no concomitant OAT1/OAT3 inhibitor",
      notes = "Tested specifically because avibactam is an OAT1/OAT3 substrate in vitro (Li 2019 Methods, 'Selection of covariates'). Not retained: exposures in the 133 patients receiving one differed by <= 25% from the 1,631 who did not (Li 2019 Table 3). Recorded here as documentation only; the name is descriptive and is NOT registered in inst/references/covariate-columns.md, which is permissible because covariatesDataExcluded entries are never referenced by model().",
      source_name = "CMEDOAT"
    ),
    DIS_BACTEREMIA = list(
      description = "Baseline bacteremia indicator",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = no bacteremia at baseline",
      notes = "Screened; not retained. 88 phase III subjects with baseline bacteremia (Li 2019 Table 3).",
      source_name = "BBACTERM"
    ),
    BMI = list(
      description = "Body mass index (obesity status)",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened as obesity status; not retained.",
      source_name = "BMI"
    ),
    WBC = list(
      description = "White blood cell count",
      units = "cells/uL",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened at a 12,000/uL cutoff as a marker of systemic disturbance; not retained.",
      source_name = "WBC"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 2249,
    n_observations = 13735,
    n_studies = 18,
    disease_state = "Adults with complicated intra-abdominal infection (786, 34.9%), complicated urinary tract infection including acute pyelonephritis (705, 31.3%) or nosocomial pneumonia including ventilator-associated pneumonia (413, 18.4%), plus 345 healthy subjects or subjects with renal impairment from Phase 1 studies (15.3%)",
    renal_function = "Estimated Cockcroft-Gault creatinine clearance 11-610 mL/min, spanning end-stage renal disease (with and without hemodialysis) to augmented renal clearance. Augmented renal clearance was defined as a measured creatinine clearance >= 140 mL/min by 8-hour urine collection and was recorded only in study CXL-PK-04.",
    dose_range = "Avibactam 500 mg every 8 hours as a 2-hour intravenous infusion in subjects with CrCL > 50 mL/min, with label-recommended reductions to 250 mg q8h (CrCL 31-50), 187.5 mg q12h (CrCL 16-30), 187.5 mg q24h (CrCL 6-15) and 187.5 mg q48h (ESRD). Given in a fixed 1:4 ratio with ceftazidime.",
    regions = "Global, including dedicated Chinese, Japanese, Korean, Taiwanese and Vietnamese subgroups",
    studies = "11 Phase 1 studies, 2 Phase 2 studies (cIAI and cUTI) and 5 Phase 3 trials: RECLAIM 1/2, RECLAIM 3, RECAPTURE 1/2, REPRISE and REPROVE",
    unbound_fraction = "0.92 for avibactam (Li 2019 Methods, 'Exposure-response analysis': free plasma concentrations were taken to be 92% of total). The model predicts TOTAL plasma concentration.",
    notes = "Estimation used FOCE-INTER in NONMEM 7.2 for model building; the final models were re-estimated with SAEM followed by importance sampling. Correlation between some random-effect parameters was high (-0.36 < r < 0.99), which is reproduced faithfully in the OMEGA block below."
  )

  ini({
    # =====================================================================
    # Structural parameters (Li 2019 Table 2). Values refer to the
    # reference subject: CRCL = 80 mL/min (where both arms of the renal
    # hinge equal 1), WT = 70 kg, Phase 1 subject with no infection,
    # non-Asian, APACHE II <= 10, no ESRD, no dialysis, not ventilated.
    # =====================================================================
    lcl <- log(10.2); label("Avibactam clearance at CRCL = 80 mL/min (L/h)")          # Table 2: theta1 CL = 10.2 L/h (RSE 1.8%)
    lvc <- log(11.1); label("Avibactam central volume of distribution at WT = 70 kg (L)") # Table 2: theta2 Vc = 11.1 L (RSE 9.9%)
    lq <- log(5.44); label("Avibactam inter-compartmental clearance (L/h)")           # Table 2: theta3 Q = 5.44 L/h (RSE 13.9%)
    lvp <- log(6.91); label("Avibactam peripheral volume of distribution (L)")        # Table 2: theta4 Vp = 6.91 L (RSE 6.5%)

    lcl_dial <- log(20.8); label("Avibactam clearance in hemodialysis patients (L/h)")
    # Table 2: theta6 'CL estimate for patients on dialysis' = 20.8 L/h (RSE
    # 9.6%). An ABSOLUTE clearance that replaces the reference-clearance
    # construction for dialysis patients, not a multiplier - hence an lcl_*
    # parameter rather than an e_* covariate effect (control stream: 'IF
    # (DIAL.EQ.0) THEN CLT1 = THETA(1)*CLPP*CLESRD ELSE CLT1 = THETA(6)').

    # =====================================================================
    # Renal-function effects on clearance. Hinged at CRCL = 80 mL/min: a
    # power arm below and a linear arm at or above, both equal to 1 at the
    # hinge.
    # =====================================================================
    e_crcl_cl_lt80 <- 1.05; label("Avibactam CRCL power exponent on CL below 80 mL/min (unitless)")
    # Table 2: theta7 'Power CrCL (< 80) on CL' = 1.05 (RSE 2.4%). Results:
    # 'estimated as a power function of 1.05 indicating an approximately
    # linear relationship'.
    e_crcl_cl_ge80 <- 0.00279; label("Avibactam CRCL linear slope on CL at or above 80 mL/min (per mL/min)")
    # Table 2: theta8 'Linear CrCL (>= 80) on CL' = 0.00279 (RSE 3.7%).
    # Results: 'avibactam CL increased by 27.9% for an increase of 100 mL/min
    # in CrCL over 80 mL/min', which is 100 * 0.00279 exactly.

    e_renalimp_esrd_cl <- 0.0678; label("Multiplicative factor on avibactam CL for end-stage renal disease off dialysis (unitless factor, not a 1+theta shift)")
    # Table 2: theta5, printed as 'CL estimate for patients with ESRD' = 0.0678
    # (RSE 8.3%), and described in Results as 'CL was 0.0678 L/h off dialysis'.
    # BOTH of those readings are wrong: the control stream uses theta5 as a
    # MULTIPLIER ('CLESRD = THETA(5)' feeding 'CLT1 = THETA(1)*CLPP*CLESRD'),
    # giving an effective ESRD clearance of 10.2 * 0.0678 = 0.692 L/h. Three
    # independent checks confirm the multiplicative reading: (i) an absolute
    # 0.0678 L/h implies a terminal half-life of roughly six days, which is
    # not credible for a small renally-cleared molecule; (ii) 0.692 L/h
    # reproduces the paper's own Table 4 ESRD avibactam AUCss,0-24 of
    # 187.5 mg / (0.692 L/h * 2 days) = 136 mg*h/L against the 127-168 printed;
    # (iii) the successor model Xie_2025_aztreonam_avibactam.R estimates the
    # same effect as leaving 7.7% of reference clearance, against 6.78% here.
    # Recorded in the validation vignette's Assumptions and deviations section.

    e_renal_arc_cl <- 0.992; label("Multiplicative scaling of the supra-80 CRCL slope on avibactam CL for augmented renal clearance (unitless factor)")
    # Table 2: theta13 'Scaling factor for CrCL in subjects with ARC,
    # CL = TVCL*(1 + theta8 * theta13 * [CrCL-80])' = 0.992 (RSE 17.4%).
    # A 0.8% softening of the supra-80 slope. The flag is 0 for every Phase 3
    # patient and for the whole simulated population of the paper.

    # =====================================================================
    # Infection-type and population effects, all proportional shifts of
    # the form parameter * (1 + theta) -- the grammar Li 2019 Table 2
    # prints explicitly in its row labels, e.g. 'Population effect on Vc
    # (cUTI), Vc*(1 + theta11)'. Note this differs from the companion
    # ceftazidime model, whose population effects are bare multipliers.
    # =====================================================================
    e_cuti_vc <- 0.434; label("Proportional shift in avibactam Vc for cUTI (fraction)")                      # Table 2: theta11 = 0.434 (RSE 24%)
    e_ciai_habp_vabp_vc <- 0.329; label("Proportional shift in avibactam Vc for Phase 3 cIAI or nosocomial pneumonia (fraction)")
    # Table 2: theta12 'Population effect on Vc (cIAI, phase III, NP)' = 0.329
    # (RSE 28.6%). Table 2 misprints the accompanying formula as
    # 'Vc*(1 + theta11)'; the control stream confirms it is theta12. One
    # coefficient shared by Phase 3 cIAI, HAP and VAP. Results: 'Vc was 32.9%
    # and 43.4% higher for phase III patients with cIAI and patients with NP
    # and patients with cUTI, respectively, compared with healthy subjects.'
    e_ciai_ph2_vc <- 1.92; label("Proportional shift in avibactam Vc for the Phase 2 cIAI study cohort (fraction)")
    # Table 2: theta9 'Population effect on Vc (cIAI, phase II), Vc*(1 + theta9)'
    # = 1.92 (RSE 25.4%). REPLACES the Phase 3 cIAI shift rather than stacking
    # on it - see covariateData$STUDY_CIAI_PH2.
    e_ciai_ph2_cl <- 0.406; label("Proportional shift in avibactam CL for the Phase 2 cIAI study cohort (fraction)")
    # Table 2: theta10 'Population effect on CL (cIAI, phase II), CL*(1 + theta10)'
    # = 0.406 (RSE 23.2%). The ONLY infection-related clearance effect in the
    # avibactam model; Phase 3 cIAI, cUTI and NP subjects carry none.

    e_apache_ii_sev_cl <- -0.197; label("Proportional shift in avibactam CL for the elevated-APACHE-II (score > 10) stratum (fraction)")
    # Table 2: theta15 'APACHE II on CL, CL*(1 + theta15)' = -0.197 (RSE 8.7%).
    # Results: 'The largest covariate effect on CL in phase III patients aside
    # from renal function was a 19.7% decrease for APACHE II score > 10.'
    e_race_asian_oth_cl <- -0.0865; label("Proportional shift in avibactam CL for non-Chinese, non-Japanese Asian race (fraction)")
    # Table 2: theta22 'ASN on CL, CL*(1 + theta22)' = -0.0865 (RSE 20.2%).
    # Results: 'an estimated 8.65% lower CL (translating to a 9.5% increase in
    # AUC) for non-Chinese, non-Japanese Asians'.
    e_mech_vent_vc <- 0.175; label("Proportional shift in avibactam Vc when a ventilator is present on the PK sampling day (fraction)")
    # Table 2: theta28 'NPv on Vc, Vc*(1 + theta28)' = 0.175 (RSE 53.3%).
    # Results: 'Patients with NPv had estimated Vc 17.5% higher than non-NPv
    # patients.'
    e_wt_vc <- 1.08; label("Body-weight power exponent on avibactam Vc, WT/70 (unitless)")
    # Table 2: theta14 'WT on Vc (WT/70.0)^theta14' = 1.08 (RSE 7.8%). Results:
    # weights at the 10th (51 kg) and 90th (95 kg) percentiles give 29% lower
    # and 39% higher Vc, which reproduce from (WT/70)^1.08.

    # =====================================================================
    # Inter-individual variability: a full 4x4 OMEGA block over CL, Vc, Vp
    # and Q, in the control stream's ETA order (ETA1 CL, ETA2 V1, ETA3 V2,
    # ETA4 Q).
    #
    # SCALE. Li 2019 Table 2 footnote b marks these rows 'Reported as
    # variance', so the diagonals are log-scale variances and the
    # off-diagonals covariances, used verbatim. The table's third column
    # carries etashrinkage for these rows, not a second copy of the
    # variance.
    #
    # VERIFICATION. All six correlations printed in Table 2 reproduce from
    # these ten numbers: r(Vc,CL) = 0.125/sqrt(0.349*1.147) = 0.20 (table
    # 0.2); r(Vp,CL) = 0.611/sqrt(0.349*1.494) = 0.85 (0.85);
    # r(Vp,Vc) = -0.426/sqrt(1.147*1.494) = -0.33 (-0.33);
    # r(Q,CL) = 1.231/sqrt(0.349*6.359) = 0.83 (0.83);
    # r(Q,Vc) = -0.978/sqrt(1.147*6.359) = -0.36 (-0.36);
    # r(Q,Vp) = 3.059/sqrt(1.494*6.359) = 0.99 (0.99).
    # =====================================================================
    etalcl + etalvc + etalvp + etalq ~ c(
      0.349,
      0.125, 1.147,
      0.611, -0.426, 1.494,
      1.231, -0.978, 3.059, 6.359
    )

    # =====================================================================
    # Residual unexplained variability, switched across three study-phase
    # strata.
    #
    # SCALE. Unlike the companion ceftazidime model, the avibactam residual
    # parameters are THETAs on the STANDARD-DEVIATION scale with $SIGMA
    # fixed to 1 (control stream $ERROR: 'W1 = SQRT(WA1**2 +
    # (WP1*IPRED)**2)'), and Li 2019 Table 2 does NOT mark these four rows
    # with its 'reported as variance' footnote. They are therefore used
    # verbatim, with no square root taken.
    #
    # UNITS. The control stream sets S1 = V1/1000 with doses in mg, so its
    # additive term is in ng/mL. This model works in mg/L, so the additive
    # standard deviation is divided by 1000. Only the Phase 1 stratum has
    # an additive component.
    # =====================================================================
    propSdPhase1 <- 0.173; label("Proportional residual standard deviation, Phase 1 studies (fraction)")   # Table 2: theta17 'Proportional error, phase I' = 0.173 (RSE 0.1%)
    addSdPhase1 <- 0.0446; label("Additive residual standard deviation, Phase 1 studies (mg/L)")
    # Table 2: theta18 'Additive variability, phase I' = 44.6 ng/mL (RSE 0.5%)
    # = 0.0446 mg/L.
    propSdPhase2 <- 0.492; label("Proportional residual standard deviation, Phase 2 studies (fraction)")   # Table 2: theta19 'Proportional variability, phase II' = 0.492 (RSE 3%)
    propSdPhase3 <- 0.363; label("Proportional residual standard deviation, Phase 3 studies (fraction)")   # Table 2: theta20 'Proportional variability, phase III' = 0.363 (RSE 1.1%)
  })

  model({
    # ------------------------------------------------------------------
    # 1. Renal-function factor on clearance. Three mutually exclusive
    #    regimes, written with 0/1 indicator arithmetic rather than
    #    branches so every arm stays finite for every positive CRCL:
    #
    #      power arm   CRCL <  80, not ESRD, not ARC:  (CRCL/80)^power
    #      linear arm  CRCL >= 80, not ESRD:           1 + slope*arc*(CRCL-80)
    #      otherwise (any ESRD subject, or an ARC subject below 80):  1
    #
    #    The ARC flag both softens the supra-80 slope and excludes the
    #    subject from the sub-80 power arm; both behaviours come straight
    #    from the control stream's CLCLCR block.
    # ------------------------------------------------------------------
    arc_slope <- 1 + (e_renal_arc_cl - 1) * RENAL_ARC

    renal_lo <- (CRCL < 80) * (1 - RENALIMP_ESRD) * (1 - RENAL_ARC)
    renal_hi <- (CRCL >= 80) * (1 - RENALIMP_ESRD)

    renal_cl <- (CRCL / 80)^e_crcl_cl_lt80 * renal_lo +
      (1 + e_crcl_cl_ge80 * arc_slope * (CRCL - 80)) * renal_hi +
      (1 - renal_lo - renal_hi)

    # ------------------------------------------------------------------
    # 2. Clearance. A dialysis patient takes the separate absolute
    #    clearance; everyone else takes the reference clearance scaled by
    #    the Phase 2 cIAI shift and, for ESRD subjects, the ESRD factor.
    #    The renal hinge and the APACHE II / race factors multiply both
    #    arms, exactly as the control stream does
    #    ('CL = CLCLCR*CLT1*EXP(MU_1+ETA(1))' with 'MU_1 = LOG(CLCOV)').
    # ------------------------------------------------------------------
    clt1 <- exp(lcl) *
      (1 + e_ciai_ph2_cl * STUDY_CIAI_PH2) *
      (1 + (e_renalimp_esrd_cl - 1) * RENALIMP_ESRD) *
      (1 - RRT_HEMODIAL_STATUS) +
      exp(lcl_dial) * RRT_HEMODIAL_STATUS

    cl <- renal_cl * clt1 * exp(etalcl) *
      (1 + e_apache_ii_sev_cl * APACHE_II_SEV) *
      (1 + e_race_asian_oth_cl * RACE_ASIAN_OTH)

    # ------------------------------------------------------------------
    # 3. Central volume. The infection-type indicators are mutually
    #    exclusive, and within cIAI the Phase 2 cohort takes theta9
    #    INSTEAD OF theta12, so the two cIAI terms are gated on
    #    STUDY_CIAI_PH2 and its complement.
    # ------------------------------------------------------------------
    ciai_ph3 <- DIS_CIAI - STUDY_CIAI_PH2

    pop_vc <- 1 +
      e_cuti_vc * DIS_CUTI +
      e_ciai_ph2_vc * STUDY_CIAI_PH2 +
      e_ciai_habp_vabp_vc * (ciai_ph3 + DIS_HABP + DIS_VABP)

    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc * pop_vc *
      (1 + e_mech_vent_vc * MECH_VENT)

    vp <- exp(lvp + etalvp)
    q <- exp(lq + etalq)

    # ------------------------------------------------------------------
    # 4. Micro-constants and the two-compartment IV disposition. Dosing is
    #    a zero-order intravenous infusion into the central compartment.
    # ------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ------------------------------------------------------------------
    # 5. Observation. Dose in mg, volumes in L -> concentration in mg/L.
    #    This is the TOTAL plasma concentration; multiply by 0.92 to
    #    obtain the free concentration the paper's 50% fT > 1 mg/L target
    #    is defined on (see population$unbound_fraction).
    #
    #    Residual variability is switched across three strata; both phase
    #    indicators 0 selects the Phase 1 stratum, the only one carrying
    #    an additive component.
    # ------------------------------------------------------------------
    Cc <- central / vc

    phase1 <- 1 - STUDY_CAZAVI_PHASE2 - STUDY_CAZAVI_PHASE3
    propSd <- propSdPhase1 * phase1 +
      propSdPhase2 * STUDY_CAZAVI_PHASE2 +
      propSdPhase3 * STUDY_CAZAVI_PHASE3
    addSd <- addSdPhase1 * phase1

    Cc ~ add(addSd) + prop(propSd)
  })
}
