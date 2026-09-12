Yang_2026_copd_fev1_adipd_mbma <- function() {
  description <- "MBMA. Combined aggregated-data + individual-patient-data (ADIPD) longitudinal model of morning trough forced expiratory volume in 1 second (FEV1) in chronic obstructive pulmonary disease, fit by NONMEM 7.5.1 to 4,137 arm-mean FEV1 observations from 296 published randomized trials (250,543 patients) jointly with individual FEV1 records from 2,241 patients in two fluticasone furoate / vilanterol trials (Yang 2026). FEV1 = baseline - linear disease progression + immediate placebo effect + drug effect + a post-bronchodilator reconciliation term. Twenty-three compounds are carried, each with an Emax or a constant effect in per-arm total daily dose, four class-level effect-onset time courses, a LABA-LAAC power interaction, and background-therapy contributions proportional to the fraction of the arm on each drug class. Random effects are declared at THREE levels -- between-study (eta_study_*), between-arm for aggregated records (eta_arm_base), and between-subject for individual records -- because the source used NONMEM $LEVEL with an interoccasion-like arm random effect; rxode2 draws one level per solve, so see the vignette for how to simulate each level. Aggregated baselines use a normal approximation to the mean of a log-normal to avoid aggregation bias. There is no PK layer: drug effects are driven by per-arm total daily dose supplied as covariate columns."

  reference <- paste(
    "Yang L, Llanos-Paez C, Yang S, Ambery C, Berges A, Kjellsson MC,",
    "Karlsson MO. A Combined Model-Based Meta-Analysis of Aggregated and",
    "Individual FEV1 Data From Randomized COPD Trials.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15(1):e70059.",
    "doi:10.1002/psp4.70059.",
    "Final parameter estimates are in Supporting Information Table S3;",
    "the model equations are in the Supporting Information section",
    "'NONMEM control stream for the combined ADIPD model'.",
    "The structural skeleton and the aggregated-data set are inherited from",
    "Llanos-Paez C, Ambery C, Yang S, Beerahee M, Plan EL, Karlsson MO.",
    "Joint longitudinal model-based meta-analysis of FEV1 and exacerbation",
    "rate in randomized COPD trials.",
    "J Pharmacokinet Pharmacodyn. 2023;50(4):297-314.",
    sep = " "
  )

  vignette <- "Yang_2026_copd_fev1"

  # The observation is absolute FEV1 in litres -- the canonical `FEV1`
  # compartment. For an aggregated-data record it is the ARM-MEAN FEV1; for an
  # individual-patient record it is that subject's FEV1.

  units <- list(
    time          = "week (weeks since randomization; the disease-progression slope and the effect-onset rates are reported per year and per week respectively -- see each label)",
    dosing        = "ug/day (per-arm TOTAL DAILY dose supplied through the CONMED_<drug>_DOSE covariate columns, NOT as rxode2 dose events; this model has no PK layer. Roflumilast, cilomilast and the two oral small molecules are also expressed as ug/day so that one unit serves every column)",
    concentration = "L (FEV1 absolute volume, observation FEV1)"
  )

  covariateData <- list(
    # ---- Record-type and meta-analysis bookkeeping -----------------------
    DTYPE_AGGREGATED = list(
      description        = "Record-type indicator. 1 = an aggregated-data record (one arm-mean FEV1 from a published trial), 0 = an individual-patient-data record (one subject's FEV1).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (individual-patient record).",
      notes              = "Source column DTYPE, coded 1 = aggregated and 2 = individual; re-coded here to a 0/1 binary so it reads as an ordinary indicator. This column switches THREE things at once and is the single most load-bearing covariate in the model: (1) the baseline uses the normal approximation to the mean of a log-normal for aggregated records and the plain log-normal for individual records (paper Section 2.3.2, Equations 1-8); (2) the residual error is additive scaled by 1/sqrt(NARM) for aggregated records and a power model with its own IIV for individual records (paper Equation 10); (3) the disease-progression slope and the vilanterol reference efficacy carry additional between-subject etas on individual records only. Aggregated and individual records share every structural parameter.",
      source_name        = "DTYPE"
    ),
    NARM = list(
      description        = "Number of patients contributing to this study arm.",
      units              = "(count of patients)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column NTRT. Used ONLY on aggregated-data records, in two places: the residual-error weight (additive error scaled by 1/sqrt(NARM), paper Equation 10 'err = eps_AD / sqrt(N_ij)') and the standard deviation of the normal approximation to the arm-mean baseline (paper Equations 6-8, variance divided by N). The arm sizes in the aggregated data range from 18 to 5,724 patients (paper Section 2.3.2), which is what makes the central-limit-theorem approximation valid. Set to 1 on individual-patient records, where it is not referenced.",
      source_name        = "NTRT"
    ),
    MEAS_POSTBD = list(
      description        = "Indicator that this FEV1 record was measured AFTER a short-acting bronchodilator. 1 = post-bronchodilator measurement, 0 = pre-bronchodilator (trough) measurement.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (pre-bronchodilator trough FEV1, the model's primary scale).",
      notes              = "Source column POSTBD. A post-bronchodilator record is predicted from a baseline shifted upward by the median absolute reversibility of 0.18 L and from a long-acting bronchodilator effect reduced by the estimated fraction rel_postbd (paper Section 2.4 note: 'if both FEV1 baseline and FEV1 during treatment were measured post-SABD, a mean absolute reversibility of 0.18 L was added to FEV1 baseline and a fractional reduction in the overall LABD effect was estimated'). Distinct from FEV1_PBD_ANCHOR below, which handles a different set of studies.",
      source_name        = "POSTBD"
    ),
    FEV1_PBD_ANCHOR = list(
      description        = "Reported FEV1 value used by the source's post-bronchodilator baseline reconciliation term, in litres. 0 for every record to which the reconciliation does not apply -- which is every record outside six specific aggregated-data studies, and every record used for simulation.",
      units              = "L",
      type               = "continuous",
      reference_category = "0 (no reconciliation; the correction term vanishes).",
      notes              = "This column exists to reproduce an ESTIMATION-TIME data-reconciliation term and MUST BE SET TO 0 FOR ALL SIMULATION. In the source control stream the term is written 'IF(REF.EQ.93.AND.POSTBD.EQ.0) POSTBDCORR = FEV1 * (1-THETA(53))' and repeated for REF 127, 169, 319, 399 and 438: six aggregated-data studies whose absolute FEV1 was reconstructed as (post-SABD baseline + change from baseline) while the model's own baseline is pre-bronchodilator. The source multiplies the record's OWN OBSERVED value by (1 - postbd_recon) and adds it to the prediction. Using the dependent variable inside the prediction cannot be expressed in rxode2 (the observation is unknown while solving) and is not meaningful for forward simulation, so the observed value is exposed here as an explicit input column instead. Setting it to 0 -- the default and the only sensible simulation choice -- reproduces the published model exactly for every record outside those six studies. The behaviour is recorded in the vignette's Assumptions and deviations section.",
      source_name        = "FEV1 (the DV itself, for REF in {93, 127, 169, 319, 399, 438} with POSTBD = 0)"
    ),
    OCS_NONRESPONDER = list(
      description        = "Indicator that the study enrolled ONLY patients who had not responded to oral corticosteroids. 1 = OCS-non-responder-only study, 0 otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (an unselected study population).",
      notes              = "Source flag OCNR, set by 'IF (REF.EQ.636) OCNR = 1' for a single aggregated-data study, with the control-stream comment 'only non-responders to OCS are included, virtually zero effect'. It zeroes the inhaled-corticosteroid component of the anti-inflammatory effect (rel_cs = 1 - OCS_NONRESPONDER) but leaves the non-steroid anti-inflammatory agents (roflumilast, cilomilast, AZD9668, PH797804) untouched. Carries no estimated parameter.",
      source_name        = "OCNR (REF = 636)"
    ),
    INCL_EXAC_REQUIRED = list(
      description        = "Indicator that the trial's inclusion criteria required a documented history of COPD exacerbations. 1 = exacerbation history required for enrolment, 0 = not required.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no exacerbation-history entry requirement).",
      notes              = "Source column INCL. A STUDY-LEVEL design covariate, not a patient characteristic: it marks trials that enriched for exacerbating patients, who have lower lung function. It carries an estimated effect on baseline (e_incl_exac_base = -0.0213) and additionally enters the fixed-coefficient regressions that impute missing background-therapy fractions. The paper notes this covariate 'only existed as a covariate for AD as all IPD individuals had no exacerbation history' (Section 3.3), so it is 0 on every individual-patient record. Related to but distinct from NEXAC12M, which is a patient-level count of prior exacerbations rather than a study entry criterion.",
      source_name        = "INCL"
    ),

    # ---- Patient / arm characteristics -----------------------------------
    AGE = list(
      description        = "Subject age (individual records) or arm-mean age (aggregated records), at randomization.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Centred at 63.4 years (control stream 'TVB = TVBL*(1 + THETA(33)*(IMPAGE - 63.4))'), the pooled mean across the aggregated and individual data. NOTE this differs from the sibling model Yang_2026_copd_fev1_ipd, which centres at 62 years because it pools only the two individual-patient studies; the two centrings are NOT interchangeable. Set to the sentinel -99 when the arm's mean age was not reported, which activates the source's fixed-coefficient imputation regression (see the model() code and the vignette source-trace table). In the combined model age enters the BASELINE only -- the age effect on vilanterol efficacy found in the individual-patient model was removed as redundant once the baseline-on-efficacy relationship was in place (paper Section 2.4, 'Redundant covariate relationships were removed').",
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Female-sex indicator (individual records) or the arm's female fraction (aggregated records), on a 0-1 scale.",
      units              = "(binary on individual records; fraction 0-1 on aggregated records)",
      type               = "binary",
      reference_category = "0 (male). The source centres the male-coded column at the pooled male fraction 0.671, so the model's typical value corresponds to a mixed-sex arm rather than to either sex.",
      notes              = "Source column SEX is coded 1 = male and is centred at 0.671 as '(1 + THETA(62)*(SEX - 0.671))'. The canonical SEXF is its complement (SEXF = 1 - SEX), so the identical algebra is '(1 + e_sexf_base * (0.329 - SEXF))' with the published coefficient e_sexf_base = +0.276389 carried over unchanged in both value and sign. Check: a female subject gives 1 + 0.276389*(0.329 - 1) = 0.8145 and a male gives 1 + 0.276389*0.329 = 1.0909, a female/male baseline ratio of 0.747 -- matching the 0.752 ratio implied by the sibling individual-patient model's coefficient of -0.248, and the paper's repeated statement that female sex relates to a lower baseline (Abstract, Sections 3.3 and 4). On an aggregated record the column carries the arm's female FRACTION, which is the aggregation of the individual indicator and is exactly what the linear covariate form requires.",
      source_name        = "SEX (1 = male, centred at 0.671)"
    ),
    SMOKE = list(
      description        = "Current-smoker indicator (individual records) or the arm's current-smoker fraction (aggregated records), on a 0-1 scale.",
      units              = "(binary on individual records; fraction 0-1 on aggregated records)",
      type               = "binary",
      reference_category = "Centred at the pooled current-smoker fraction 0.463, so the typical value corresponds to a mixed arm rather than to either smoking status.",
      notes              = "Source column FL_SMOK, coded 1 = current smoker, matching the canonical SMOKE coding exactly (no transformation). Entered as '(1 + THETA(66)*(IMPFL_SMOK - 0.463))' with e_smoke_base = +0.0365392, i.e. a current smoker has a 3.7% higher FEV1 baseline than a non-current smoker -- the paper quotes 3.6% (Section 4) and cautions explicitly that this is NOT causal: patients tend to stop smoking at more severe COPD stages. Set to the sentinel -99 when the arm's smoking fraction was not reported, which activates the source's fixed-coefficient logistic imputation (13.5% of records were missing, per Supporting Information Table S5 discussion).",
      source_name        = "FL_SMOK"
    ),
    DIS_COPD_GOLD = list(
      description        = "GOLD spirometric severity stage as an ordinal 1-4 category (1 = mild, 2 = moderate, 3 = severe, 4 = very severe). Used on INDIVIDUAL-patient records.",
      units              = "(ordinal stage 1-4)",
      type               = "ordinal",
      reference_category = "3 (severe), the cohort median and the model's centring constant.",
      notes              = "Source column FL_COPD; Supporting Information Table S1 footnote 1 defines it as the '% predicted GOLD Stage Category at the screening phase'. Set to the sentinel -99 when missing, which the source replaces with the median stage 3 (0.62% of records; paper Section 2.4). On AGGREGATED records this column is not used -- the arm's severity is the midpoint of DIS_COPD_GOLD_LOW and DIS_COPD_GOLD_HIGH instead. The combined model uses a SINGLE linear slope in the stage, not the hockey-stick of the sibling individual-patient model: that linearization is deliberate and is the paper's remedy for aggregation bias (Section 2.3.1, 'the covariate effects of lowest/highest disease severity on baseline ... was revised to a single covariate (lowest+highest)/2 on baseline to keep consistent with IPD part of the model'). The source keeps the two hockey-stick slopes in the control stream as THETA(64) and THETA(65) but fixes both to zero.",
      source_name        = "FL_COPD"
    ),
    DIS_COPD_GOLD_LOW = list(
      description        = "Lowest GOLD spirometric stage admitted by the trial's inclusion criteria, as an ordinal 1-4 category. Used on AGGREGATED-data records.",
      units              = "(ordinal stage 1-4)",
      type               = "ordinal",
      reference_category = "n/a -- combined with DIS_COPD_GOLD_HIGH into the arm's mean stage, which is centred at 3.",
      notes              = "Source column LOWDS. A published trial reports its inclusion range rather than a per-patient severity distribution, so the aggregated arm's mean severity is taken as the midpoint (DIS_COPD_GOLD_LOW + DIS_COPD_GOLD_HIGH) / 2 (control stream 'IF (DTYPE.EQ.1) IMPFL_COPD=(LOWDSact+HIGHDSact)/2'). Also enters the fixed-coefficient regressions that impute missing background-therapy fractions, where it is centred at 2.",
      source_name        = "LOWDS"
    ),
    DIS_COPD_GOLD_HIGH = list(
      description        = "Highest GOLD spirometric stage admitted by the trial's inclusion criteria, as an ordinal 1-4 category. Used on AGGREGATED-data records.",
      units              = "(ordinal stage 1-4)",
      type               = "ordinal",
      reference_category = "n/a -- combined with DIS_COPD_GOLD_LOW into the arm's mean stage, which is centred at 3.",
      notes              = "Source column HIGHDS. See DIS_COPD_GOLD_LOW. In the background-therapy imputation regressions it is centred at 4.",
      source_name        = "HIGHDS"
    ),

    # ---- Background (non-randomized) therapy in the arm -------------------
    BGTHER_ICS_RUNIN_PCT = list(
      description        = "Percentage of the arm receiving background inhaled corticosteroid during the run-in period, before randomized treatment begins.",
      units              = "% of patients in the arm (0-100)",
      type               = "continuous",
      reference_category = "0 (no background inhaled corticosteroid).",
      notes              = "Source column ICSpRUNIN. Determines the FEV1 already present at time 0 from background therapy, and also the starting point of the corticosteroid effect-onset time course (control stream 'ICSTC = ICSR/100 + (1-ICSR/100)*(1-EXP(-DSCS*TIME))'), so that an arm already on inhaled corticosteroids starts part-way up the onset curve. The background contribution is valued at the fluticasone propionate b.i.d. Emax (control stream 'ICS = DMX9 * ICSR/100 ;based on fluticasone BID'). Sentinel 9999 means completely unknown and is treated as 0; sentinel 7777 means the class was used but the fraction was not reported and activates the source's fixed-coefficient logistic imputation.",
      source_name        = "ICSpRUNIN"
    ),
    BGTHER_LABA_RUNIN_PCT = list(
      description        = "Percentage of the arm receiving a background long-acting beta-2 agonist during the run-in period.",
      units              = "% of patients in the arm (0-100)",
      type               = "continuous",
      reference_category = "0 (no background LABA).",
      notes              = "Source column LABApRUNIN. Valued at the salmeterol effect (control stream 'LABA = DE18 * LABAR/100 ;based on salmeterol') and sets the starting point of the once-daily beta-agonist onset curve. Sentinels 9999 (unknown, treated as 0) and 7777 (imputed) as for BGTHER_ICS_RUNIN_PCT.",
      source_name        = "LABApRUNIN"
    ),
    BGTHER_LAAC_RUNIN_PCT = list(
      description        = "Percentage of the arm receiving a background long-acting anticholinergic during the run-in period.",
      units              = "% of patients in the arm (0-100)",
      type               = "continuous",
      reference_category = "0 (no background LAAC).",
      notes              = "Source column LAACpRUNIN. Valued at the tiotropium HandiHaler effect at its 18 ug/day reference dose (control stream 'LAAC = DE20CD*LAACR/100 ;based on tiotropium HandiHaler', where DE20CD carries no onset time course) and sets the starting point of the once-daily anticholinergic onset curve. Sentinels 9999 and 7777 as above.",
      source_name        = "LAACpRUNIN"
    ),
    BGTHER_ICS_MAINT_PCT = list(
      description        = "Percentage of the arm receiving background inhaled corticosteroid during the randomized maintenance period.",
      units              = "% of patients in the arm (0-100)",
      type               = "continuous",
      reference_category = "0 (no background inhaled corticosteroid).",
      notes              = "Source column ICSpMAINT. Replaces the run-in fraction for every record after time 0 (control stream: the run-in contributions are carried at TIME = 0 and the maintenance contributions at TIME > 0). Sentinels 9999 and 7777 as for the run-in columns; the imputation uses the same regression and the same random effect, so an arm whose run-in and maintenance fractions are both missing receives the same imputed value for both.",
      source_name        = "ICSpMAINT"
    ),
    BGTHER_LABA_MAINT_PCT = list(
      description        = "Percentage of the arm receiving a background long-acting beta-2 agonist during the randomized maintenance period.",
      units              = "% of patients in the arm (0-100)",
      type               = "continuous",
      reference_category = "0 (no background LABA).",
      notes              = "Source column LABApMAINT. See BGTHER_ICS_MAINT_PCT.",
      source_name        = "LABApMAINT"
    ),
    BGTHER_LAAC_MAINT_PCT = list(
      description        = "Percentage of the arm receiving a background long-acting anticholinergic during the randomized maintenance period.",
      units              = "% of patients in the arm (0-100)",
      type               = "continuous",
      reference_category = "0 (no background LAAC).",
      notes              = "Source column LAACpMAINT. See BGTHER_ICS_MAINT_PCT.",
      source_name        = "LAACpMAINT"
    ),

    # ---- Randomized treatment: per-arm total daily dose -------------------
    # Members of the CONMED_<INN>_DOSE family. Every column is 0 outside its
    # own drug's arms; all 23 zero identifies a placebo (background-therapy
    # only) arm. Units are ug/day of TOTAL DAILY dose throughout, matching the
    # source's avdostot data item and the unit of each drug's ED50.
    CONMED_ACLIDINIUM_DOSE = list(
      description        = "Per-arm total daily aclidinium dose.",
      units              = "ug/day",
      type               = "continuous",
      reference_category = "0 (no aclidinium).",
      notes              = "Reference doses: 200 ug/day for the once-daily regimen and 800 ug/day (400 ug b.i.d.) for the twice-daily regimen, which carry SEPARATE reference efficacies and SEPARATE ED50s selected by FORM_ACLIDINIUM_BID. The control stream comments that the q.d. regimen is 'not used clinically, but present in dataset'.",
      source_name        = "avdostot (drgNo = 2)"
    ),
    CONMED_ARFORMOTEROL_DOSE = list(
      description        = "Per-arm total daily arformoterol dose.",
      units              = "ug/day",
      type               = "continuous",
      reference_category = "0 (no arformoterol).",
      notes              = "Reference dose 50 ug/day. Arformoterol carries NO parameters of its own: the source assumes its Emax equals formoterol's and its ED50 is exactly half of formoterol's, because arformoterol is the single active enantiomer of racemic formoterol (control stream 'ED504 = ED5010/2 ;Arformoterol ED50 = 1/2 formoterol ED50' and 'TVDMX4 = TVDMX10 ;Emax drug 4 (Arformoterol) assume equivalent to formoterol').",
      source_name        = "avdostot (drgNo = 4)"
    ),
    CONMED_BECLOMETHASONE_DOSE = list(
      description        = "Per-arm total daily beclomethasone dose.",
      units              = "ug/day",
      type               = "continuous",
      reference_category = "0 (no beclomethasone).",
      notes              = "NO dose-response was estimated: the source fits a single constant effect and the model keys it on this column being non-zero, following the idiom already recorded for CONMED_MTX_DOSE in Mandema_2011_biologicDMARDs_mbma. The column still carries the dose so the arm is self-describing; do NOT read dose-proportionality into it.",
      source_name        = "avdostot (drgNo = 5)"
    ),
    CONMED_BUDESONIDE_DOSE = list(
      description        = "Per-arm total daily budesonide dose.",
      units              = "ug/day",
      type               = "continuous",
      reference_category = "0 (no budesonide).",
      notes              = "Reference dose 320 ug/day, i.e. the 160 ug b.i.d. regimen against which Table S3 quotes the reference efficacy.",
      source_name        = "avdostot (drgNo = 6)"
    ),
    CONMED_CILOMILAST_DOSE = list(
      description        = "Per-arm total daily cilomilast dose.",
      units              = "ug/day",
      type               = "continuous",
      reference_category = "0 (no cilomilast).",
      notes              = "Constant effect, no dose-response; keyed on the column being non-zero. Table S3 labels the estimate 'Drug7.cil.15mcg.Eff.BID', i.e. the 15 mg b.i.d. regimen. A PDE4 inhibitor, so it follows the PDE4 effect-onset time course and is NOT affected by the OCS-non-responder flag.",
      source_name        = "avdostot (drgNo = 7)"
    ),
    CONMED_FLUTICASONEPROPIONATE_DOSE = list(
      description        = "Per-arm total daily fluticasone propionate dose.",
      units              = "ug/day",
      type               = "continuous",
      reference_category = "0 (no fluticasone propionate).",
      notes              = "Constant effect for the twice-daily regimen, no dose-response; keyed on the column being non-zero. This compound is ALSO the reference against which background inhaled-corticosteroid therapy is valued (control stream 'ICS = DMX9 * ICSR/100 ;based on fluticasone BID'), so its Emax is load-bearing far beyond its own arms. Distinct from fluticasone FUROATE (CONMED_FLUTICASONEFUROATE_DOSE), a different molecule with its own dose-response.",
      source_name        = "avdostot (drgNo = 9)"
    ),
    CONMED_FORMOTEROL_DOSE = list(
      description        = "Per-arm total daily formoterol dose.",
      units              = "ug/day",
      type               = "continuous",
      reference_category = "0 (no formoterol).",
      notes              = "Reference dose 18 ug/day, i.e. the 9 ug b.i.d. regimen against which Table S3 quotes the reference efficacy. Formoterol's Emax and ED50 also determine arformoterol's (see CONMED_ARFORMOTEROL_DOSE).",
      source_name        = "avdostot (drgNo = 10)"
    ),
    CONMED_INDACATEROL_DOSE = list(
      description        = "Per-arm total daily indacaterol dose.",
      units              = "ug/day",
      type               = "continuous",
      reference_category = "0 (no indacaterol).",
      notes              = "Reference dose 75 ug/day. Indacaterol is the ONLY beta-agonist that carries the once-daily beta-agonist effect-onset time course in the final model.",
      source_name        = "avdostot (drgNo = 11)"
    ),
    CONMED_MOMETASONE_DOSE = list(
      description        = "Per-arm total daily mometasone dose.",
      units              = "ug/day",
      type               = "continuous",
      reference_category = "0 (no mometasone).",
      notes              = "Constant effect, no dose-response; keyed on the column being non-zero. The twice-daily regimen's effect is the once-daily effect multiplied by rel_mometasone_bid = 0.784, selected by FORM_MOMETASONE_BID.",
      source_name        = "avdostot (drgNo = 13)"
    ),
    CONMED_GLYCOPYRRONIUM_DOSE = list(
      description        = "Per-arm total daily glycopyrronium dose.",
      units              = "ug/day",
      type               = "continuous",
      reference_category = "0 (no glycopyrronium).",
      notes              = "Reference dose 100 ug/day. Carries the once-daily anticholinergic effect-onset time course.",
      source_name        = "avdostot (drgNo = 14)"
    ),
    CONMED_ROFLUMILAST_DOSE = list(
      description        = "Per-arm total daily roflumilast dose.",
      units              = "ug/day",
      type               = "continuous",
      reference_category = "0 (no roflumilast).",
      notes              = "Reference dose 500 ug/day. A PDE4 inhibitor: it follows the PDE4 effect-onset time course and, being a non-steroid anti-inflammatory, is NOT zeroed by the OCS-non-responder flag.",
      source_name        = "avdostot (drgNo = 16)"
    ),
    CONMED_SALMETEROL_DOSE = list(
      description        = "Per-arm total daily salmeterol dose.",
      units              = "ug/day",
      type               = "continuous",
      reference_category = "0 (no salmeterol).",
      notes              = "Constant effect for the twice-daily regimen, no dose-response; keyed on the column being non-zero. Salmeterol is ALSO the reference against which background LABA therapy is valued (control stream 'LABA = DE18 * LABAR/100 ;based on salmeterol'), so its effect is load-bearing beyond its own arms.",
      source_name        = "avdostot (drgNo = 18)"
    ),
    CONMED_TIOTROPIUM_DOSE = list(
      description        = "Per-arm total daily tiotropium dose.",
      units              = "ug/day",
      type               = "continuous",
      reference_category = "0 (no tiotropium).",
      notes              = "Reference dose 18 ug/day for the Spiriva HandiHaler dry-powder inhaler and 5 ug/day for the Respimat soft-mist inhaler, selected by FORM_TIOTROPIUM_SMI. The two devices carry separate reference efficacies; the Respimat ED50 is DERIVED from the HandiHaler Emax and the Respimat reference efficacy rather than estimated (control stream 'ED50120 = REFDSMI*(TVDMX20/THETA(14)-1)'). Open-label administration scales the effect by rel_tiotropium_ol = 0.918, selected by FORM_TIOTROPIUM_OPENLABEL. Tiotropium at its 18 ug/day reference dose is also the reference for background LAAC therapy.",
      source_name        = "avdostot (drgNo = 20)"
    ),
    CONMED_UMECLIDINIUM_DOSE = list(
      description        = "Per-arm total daily umeclidinium dose.",
      units              = "ug/day",
      type               = "continuous",
      reference_category = "0 (no umeclidinium).",
      notes              = "Constant effect, no dose-response; keyed on the column being non-zero. The source assigns the SAME effect to the once-daily and twice-daily regimens (control stream 'TVDMX22BID = TVDMX22'), so no regimen flag is needed.",
      source_name        = "avdostot (drgNo = 22)"
    ),
    CONMED_AZD9668_DOSE = list(
      description        = "Per-arm total daily AZD9668 dose.",
      units              = "ug/day",
      type               = "continuous",
      reference_category = "0 (no AZD9668).",
      notes              = "An investigational neutrophil-elastase inhibitor identified only by its development code in the source. Constant effect, no dose-response, no effect-onset time course; keyed on the column being non-zero. The estimate is very imprecise (RSE 108.5%).",
      source_name        = "avdostot (drgNo = 23)"
    ),
    CONMED_GSK233705_DOSE = list(
      description        = "Per-arm total daily GSK233705 dose.",
      units              = "ug/day",
      type               = "continuous",
      reference_category = "0 (no GSK233705).",
      notes              = "An investigational long-acting anticholinergic identified only by its development code. Constant effect, no dose-response, no onset time course; keyed on the column being non-zero.",
      source_name        = "avdostot (drgNo = 24)"
    ),
    CONMED_VILANTEROL_DOSE = list(
      description        = "Per-arm total daily vilanterol dose.",
      units              = "ug/day",
      type               = "continuous",
      reference_category = "0 (no vilanterol).",
      notes              = "Reference dose 25 ug/day. This is one of the two compounds present in the individual-patient data, and is the only drug effect that carries a between-subject random effect -- which the source applies on individual records only (control stream 'IF(DTYPE.EQ.2) TVDMX25 = (THETA(41)+ETA(26))/...'). That eta's variance was estimated to zero in the combined model.",
      source_name        = "avdostot (drgNo = 25)"
    ),
    CONMED_BEA2180_DOSE = list(
      description        = "Per-arm total daily BEA2180 dose.",
      units              = "ug/day",
      type               = "continuous",
      reference_category = "0 (no BEA2180).",
      notes              = "An investigational long-acting anticholinergic identified only by its development code. Constant effect, no dose-response, no onset time course; keyed on the column being non-zero.",
      source_name        = "avdostot (drgNo = 26)"
    ),
    CONMED_PH797804_DOSE = list(
      description        = "Per-arm total daily PH797804 dose.",
      units              = "ug/day",
      type               = "continuous",
      reference_category = "0 (no PH797804).",
      notes              = "An investigational p38 MAP-kinase inhibitor identified only by its development code. Constant effect, no dose-response, no onset time course; keyed on the column being non-zero. Grouped with the non-steroid anti-inflammatory agents, so it is NOT zeroed by the OCS-non-responder flag.",
      source_name        = "avdostot (drgNo = 27)"
    ),
    CONMED_REVEFENACIN_DOSE = list(
      description        = "Per-arm total daily revefenacin dose.",
      units              = "ug/day",
      type               = "continuous",
      reference_category = "0 (no revefenacin).",
      notes              = "Reference dose 175 ug/day. A long-acting anticholinergic, but the source applies NO effect-onset time course to it.",
      source_name        = "avdostot (drgNo = 29)"
    ),
    CONMED_OLODATEROL_DOSE = list(
      description        = "Per-arm total daily olodaterol dose.",
      units              = "ug/day",
      type               = "continuous",
      reference_category = "0 (no olodaterol).",
      notes              = "Reference dose 5 ug/day for the once-daily regimen, which carries an Emax dose-response. The twice-daily regimen instead carries a CONSTANT effect with no dose-response, selected by FORM_OLODATEROL_BID.",
      source_name        = "avdostot (drgNo = 30)"
    ),
    CONMED_BATEFENTEROL_DOSE = list(
      description        = "Per-arm total daily batefenterol dose.",
      units              = "ug/day",
      type               = "continuous",
      reference_category = "0 (no batefenterol).",
      notes              = "Reference dose 400 ug/day for the once-daily regimen, which carries an Emax dose-response; the twice-daily regimen carries a constant effect, selected by FORM_BATEFENTEROL_BID. Batefenterol is a dual-pharmacology muscarinic-antagonist / beta-2-agonist (MABA) and is the sole member of its own class in the bronchodilator interaction term.",
      source_name        = "avdostot (drgNo = 31)"
    ),
    CONMED_FLUTICASONEFUROATE_DOSE = list(
      description        = "Per-arm total daily fluticasone furoate dose.",
      units              = "ug/day",
      type               = "continuous",
      reference_category = "0 (no fluticasone furoate).",
      notes              = "Reference dose 100 ug/day. The second of the two compounds present in the individual-patient data. Unlike the other inhaled corticosteroids it carries NO effect-onset time course in the source control stream.",
      source_name        = "avdostot (drgNo = 32)"
    ),

    # ---- Regimen and device / blinding selectors -------------------------
    FORM_ACLIDINIUM_BID = list(
      description        = "Aclidinium twice-daily-regimen indicator. 1 = b.i.d., 0 = q.d.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (once-daily aclidinium).",
      notes              = "Source flag derived from dosfreqNoN = 2 for the slot carrying drgNo 2. Selects the b.i.d. reference efficacy (0.0964 L at 800 ug/day) and b.i.d. ED50 in place of the q.d. pair (0.0752 L at 200 ug/day). The control stream notes the two regimens are 'never given in ambiguous combination in dataset'.",
      source_name        = "dosfreqNo (drgNo = 2)"
    ),
    FORM_MOMETASONE_BID = list(
      description        = "Mometasone twice-daily-regimen indicator. 1 = b.i.d., 0 = q.d.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (once-daily mometasone).",
      notes              = "Source flag derived from dosfreqNoN = 2 for the slot carrying drgNo 13. Multiplies the once-daily effect by rel_mometasone_bid = 0.784.",
      source_name        = "dosfreqNo (drgNo = 13)"
    ),
    FORM_OLODATEROL_BID = list(
      description        = "Olodaterol twice-daily-regimen indicator. 1 = b.i.d., 0 = q.d.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (once-daily olodaterol).",
      notes              = "Source flag derived from dosfreqNoN = 2 for the slot carrying drgNo 30. Switches from the once-daily Emax dose-response to a separate constant b.i.d. effect (0.112 L).",
      source_name        = "dosfreqNo (drgNo = 30)"
    ),
    FORM_BATEFENTEROL_BID = list(
      description        = "Batefenterol twice-daily-regimen indicator. 1 = b.i.d., 0 = q.d.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (once-daily batefenterol).",
      notes              = "Source flag derived from dosfreqNoN = 2 for the slot carrying drgNo 31. Switches from the once-daily Emax dose-response to a separate constant b.i.d. effect (0.203 L).",
      source_name        = "dosfreqNo (drgNo = 31)"
    ),
    FORM_TIOTROPIUM_SMI = list(
      description        = "Tiotropium soft-mist-inhaler (Respimat) device indicator. 1 = Respimat, 0 = Spiriva HandiHaler dry-powder inhaler.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (HandiHaler dry-powder inhaler).",
      notes              = "Source flag dosfrmNSM. Selects the 5 ug/day reference dose and the Respimat reference efficacy (0.120 L) with its DERIVED ED50, in place of the 18 ug/day HandiHaler pair (0.122 L). The control stream notes there are 'no occurances of tio being given as HandiHaler with another drug as SMI'.",
      source_name        = "dosfrm1SM / dosfrm2SM / dosfrm3SM"
    ),
    FORM_TIOTROPIUM_OPENLABEL = list(
      description        = "Tiotropium open-label-administration indicator. 1 = open-label, 0 = blinded.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (blinded administration).",
      notes              = "Source flag DRGOLN for the slot carrying drgNo 20. Scales the tiotropium reference efficacy by rel_tiotropium_ol = 0.918 while keeping the blinded ED50 (control stream 'ED5020OL = ED5020'), i.e. open-label tiotropium is estimated to perform about 8% WORSE than the same drug given blinded. The control stream notes that 'blinded tio never given with another drug being OL in dataset'. Strictly a trial-conduct covariate rather than a formulation, but it selects between two variants of one drug's effect in exactly the way the other FORM_ members do.",
      source_name        = "DRGOL1 / DRGOL2 / DRGOL3 (drgNo = 20)"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 252784L,
    n_studies       = 298L,
    age_range       = "Aggregated arm-mean ages across the published trials, pooled mean 63.4 years; the two individual-patient studies span 40-85 years with a mean of 62.",
    age_median      = "63.4 years (the pooled mean, used as the model's centring constant)",
    weight_range    = "not used as a model covariate",
    sex_female_pct  = 32.9,
    disease_state   = "Chronic obstructive pulmonary disease across the GOLD spirometric severity range, in trials of mono-, dual- and triple-therapy with bronchodilators and anti-inflammatories. Endpoint is morning trough FEV1.",
    dose_range      = "Twenty-three compounds at their clinically studied dose ranges, given as monotherapy, dual therapy or triple therapy; see each CONMED_<drug>_DOSE covariate entry for that drug's reference dose. Placebo arms are encoded by all dose columns being zero.",
    regions         = "Multinational; the aggregated data are all published randomized COPD trials meeting the source analysis's criteria up to 24 November 2020.",
    trials_included = "298 studies in total. Aggregated data: 4,137 arm-mean trough FEV1 observations from 298 studies of 250,543 patients, inherited unchanged from the Llanos-Paez 2023 meta-analysis. Individual-patient data: NCT01053988 (n = 1025) and NCT01054885 (n = 1216), both 24-week fluticasone furoate / vilanterol trials. Those two studies were REMOVED from the aggregated data when fitting the combined model so that no observation contributes twice, leaving 296 aggregated studies alongside the 2 individual-patient studies.",
    notes           = "n_subjects is the sum of the 250,543 aggregated patients and the 2,241 individual patients; because the two individual-patient studies are also among the 298 aggregated studies (and were excluded from the aggregated side to avoid duplication), this total counts each patient once. sex_female_pct is 1 - 0.671, using the pooled male fraction 0.671 that the source uses as its sex centring constant; the current-smoker fraction is correspondingly 0.463 and the mean GOLD stage 3. Study durations run from short trials to 19 studies longer than 52 weeks (Supporting Information Figure S4C). The paper reports a between-study shrinkage of 43.4% on the disease-progression slope random effect, so study-level disease-progression draws are weakly informed."
  )

  ini({
    # ==================================================================
    # Structural model, paper Equation 10:
    #
    #   FEV1 = B_ij + PBO_ij - DP_ij + E_ij + PostBD correction + err
    #
    # with B the baseline, PBO an immediate constant placebo effect, DP a
    # linear-in-time disease progression, and E the drug effect. Full
    # equations are in the Supporting Information section 'NONMEM control
    # stream for the combined ADIPD model'.
    #
    # Values are the final estimates of Supporting Information Table S3.
    # Where the control stream's $THETA record carries more significant
    # digits than Table S3's rounded display (it is the converged run, e.g.
    # 1.06908 vs the tabulated 1.07), the full-precision value is used and
    # both are shown in the trailing comment.
    # ==================================================================

    base <- 1.06908 ; label("Typical baseline FEV1 at the reference covariate values (L)")   # Table S3 'Typical baseline FEV1 (L)' = 1.07 (RSE 0.3%); $THETA (1.) = 1.06908
    dps  <- 0.03432 ; label("Typical disease-progression slope, a DECLINE in FEV1 (L/year)") # Table S3 'Disease progression slope (L/year)' = 0.0343 (RSE 7.5%); $THETA (2.) = 0.03432; enters with a minus sign, so a positive value is a decline
    pmx  <- -0.00588548 ; label("Immediate placebo effect, constant for t > 0 (L)")          # Table S3 'Immediate placebo effect (L)' = -0.00589 (RSE 35.0%); $THETA (3.) = -0.00588548

    # The placebo TIME COURSE was linearized away: the published
    # aggregated-data model had a mixture of an Emax(t) onset and an immediate
    # effect, revised to a purely immediate effect to remove a source of
    # aggregation bias (paper Section 2.3.1). The two thetas of the withdrawn
    # mixture survive in the control stream fixed to zero -- $THETA (4.) '0 FIX
    # ; log.Placebo.T50' and $THETA (30.) '0 FIX ; log.P.fast.plac.resp' -- but
    # the equations that consumed them were deleted, so there is nothing for
    # them to multiply and they are not carried here. Neither appears in
    # Table S3.

    # ---- Effect-onset rates ---------------------------------------------
    lonset_ai   <- -0.950958 ; label("log onset rate of the anti-inflammatory effect (log 1/week)")   # Table S3 'Log of onset rate of anti-inflammatory treatments (log (week^-1))' = -0.951 (RSE 12.3%); shared by the corticosteroid and PDE4-inhibitor time courses
    onset_qdba  <-  9.36613  ; label("Onset rate of the once-daily beta-agonist effect (1/week)")     # Table S3 'Onset rate for q.d. LABA bronchodilator (week^-1)' = 9.37 (RSE 8.4%)
    onset_qdac  <- 11.5118   ; label("Onset rate of the once-daily anticholinergic effect (1/week)")  # Table S3 'Onset rate for q.d. LAAC bronchodilators (week^-1)' = 11.5 (RSE 10.0%)

    # ---- Reference efficacies (L at the drug's reference daily dose) -----
    # The source parameterises each dose-response drug by its effect AT A
    # REFERENCE DOSE rather than by Emax, so that
    #   Emax = EffRef / RefDose * (ED50 + RefDose).
    effref_acli_qd <- 0.075184  ; label("Aclidinium q.d. efficacy at 200 ug/day (L)")           # Table S3 'Reference efficacy of aclidinium 200 ug q.d. (L)' = 0.0752 (RSE 16.1%)
    effref_acli_bid <- 0.0964157; label("Aclidinium b.i.d. efficacy at 800 ug/day (L)")         # Table S3 'Reference efficacy of aclidinium 400 ug b.i.d. (L)' = 0.0964 (RSE 6.9%)
    eff_beclo      <- 0.058458  ; label("Beclomethasone constant effect (L)")                   # Table S3 'Efficacy of beclomethasone (L)' = 0.0585 (RSE 18.8%)
    effref_bud     <- 0.0336889 ; label("Budesonide efficacy at 320 ug/day (L)")                # Table S3 'Reference efficacy of budesonide 160 ug b.i.d. (L)' = 0.0337 (RSE 24.3%)
    eff_cil        <- 0.0449212 ; label("Cilomilast constant effect (L)")                       # Table S3 'Efficacy of cilomilast (L)' = 0.0449 (RSE 21.6%)
    eff_flutiprop  <- 0.0420672 ; label("Fluticasone propionate b.i.d. constant effect (L)")    # Table S3 'Efficacy of fluticasone b.i.d. (L)' = 0.0421 (RSE 8.2%); also values background ICS therapy
    effref_form    <- 0.0699607 ; label("Formoterol efficacy at 18 ug/day (L)")                 # Table S3 'Reference efficacy of formoterol 9 ug b.i.d. (L)' = 0.07 (RSE 5.7%)
    effref_indac   <- 0.12618   ; label("Indacaterol efficacy at 75 ug/day (L)")                # Table S3 'Reference efficacy of indacaterol 75 ug q.d. (L)' = 0.126 (RSE 3.7%)
    eff_mome_qd    <- 0.0750943 ; label("Mometasone q.d. constant effect (L)")                  # Table S3 'Efficacy of mometasone q.d. (L)' = 0.0751 (RSE 23.2%)
    effref_glyco   <- 0.12743   ; label("Glycopyrronium efficacy at 100 ug/day (L)")            # Table S3 'Reference efficacy of glycopyrronium 100 ug q.d. (L)' = 0.127 (RSE 4.7%)
    effref_roflu   <- 0.0811084 ; label("Roflumilast efficacy at 500 ug/day (L)")               # Table S3 'Reference efficacy of roflumilast 500 ug q.d. (L)' = 0.0811 (RSE 12.7%)
    eff_salm       <- 0.07838   ; label("Salmeterol b.i.d. constant effect (L)")                # Table S3 'Efficacy of salmeterol b.i.d. (L)' = 0.0784 (RSE 4.6%); also values background LABA therapy
    effref_tio_dpi <- 0.122236  ; label("Tiotropium HandiHaler efficacy at 18 ug/day (L)")      # Table S3 'Reference efficacy of tiotropium (blinded) 18 ug q.d. (Handihaler) (L)' = 0.122 (RSE 2.8%)
    effref_tio_smi <- 0.120457  ; label("Tiotropium Respimat efficacy at 5 ug/day (L)")         # Table S3 'Reference efficacy of tiotropium (blinded) 5 ug q.d. (Respimat) (L)' = 0.12 (RSE 4.3%)
    eff_ume        <- 0.147083  ; label("Umeclidinium constant effect, q.d. and b.i.d. (L)")    # Table S3 'Efficacy of umeclidinium q.d. (L)' = 0.147 (RSE 3.6%)
    eff_azd9668    <- 0.0163939 ; label("AZD9668 constant effect (L)")                          # Table S3 'Efficacy of AZD9668 (L)' = 0.0164 (RSE 108.5%)
    eff_gsk233705  <- 0.190635  ; label("GSK233705 constant effect (L)")                        # Table S3 'Efficacy of GSK23305 (L)' = 0.191 (RSE 21.7%); Table S3 misprints the code as GSK23305, $THETA (40.) reads 'Drug24.GSK.Eff'
    effref_vil     <- 0.109036  ; label("Vilanterol efficacy at 25 ug/day (L)")                 # Table S3 'Reference efficacy of vilanterol 25 ug q.d. (L)' = 0.109 (RSE 3.0%)
    eff_bea2180    <- 0.105779  ; label("BEA2180 constant effect (L)")                          # Table S3 'Efficacy of BEA2180 (L)' = 0.106 (RSE 8.8%)
    eff_ph797804   <- 0.0837006 ; label("PH797804 constant effect (L)")                         # Table S3 'Efficacy of PH797804 (L)' = 0.0837 (RSE 47.6%)
    effref_rev     <- 0.148945  ; label("Revefenacin efficacy at 175 ug/day (L)")               # Table S3 'Reference efficacy of revefenacin 175 ug q.d. (L)' = 0.149 (RSE 12.4%)
    effref_olo_qd  <- 0.0889892 ; label("Olodaterol q.d. efficacy at 5 ug/day (L)")             # Table S3 'Reference efficacy of olodaterol 5 ug q.d. (L)' = 0.089 (RSE 5.5%)
    eff_olo_bid    <- 0.112063  ; label("Olodaterol b.i.d. constant effect (L)")                # Table S3 'Efficacy of olodaterol b.i.d. (L)' = 0.112 (RSE 32.2%)
    effref_bat_qd  <- 0.195297  ; label("Batefenterol q.d. efficacy at 400 ug/day (L)")         # Table S3 'Reference efficacy of batefenterol 400 ug q.d. (L)' = 0.195 (RSE 11.1%)
    eff_bat_bid    <- 0.203238  ; label("Batefenterol b.i.d. constant effect (L)")              # Table S3 'Efficacy of batefenterol b.i.d. (L)' = 0.203 (RSE 13.8%)
    effref_ff      <- 0.0400295 ; label("Fluticasone furoate q.d. efficacy at 100 ug/day (L)")  # $THETA (38.) 'rel.Drug32.FF.Eff.QD' = 0.0400295; Table S3 displays this row as 'Relative reference efficacy of fluticasone 200 ug q.d. compared to b.i.d = 0.04', which mislabels the fluticasone FUROATE reference efficacy as a fluticasone propionate ratio -- the control stream's use 'TVDMX32QD = THETA(38)/REFDFF*(ED5032QD+REFDFF)' with REFDFF = 100 settles it. Recorded in the vignette Errata.

    # ---- Relative efficacies (unitless multipliers) ----------------------
    rel_mome_bid   <- 0.78399  ; label("Mometasone b.i.d. effect relative to q.d. (unitless)")               # Table S3 'Relative efficacy of mometasone b.i.d. compared to q.d.' = 0.784 (RSE 22.7%)
    rel_tio_ol     <- 0.918377 ; label("Open-label tiotropium effect relative to blinded (unitless)")        # Table S3 'Relative efficacy of tiotropium (open-label) 18 ug q.d. (Spiriva) compared to blinded administration' = 0.918 (RSE 4.7%)
    rel_postbd     <- 0.530969 ; label("Fraction of the bronchodilator effect seen in a post-SABD measurement (unitless)") # Table S3 'Fractional bronchodilator effect for postSA bronchodilator measurements' = 0.531 (RSE 12.6%)
    postbd_recon   <- 0.884349 ; label("Post-bronchodilator baseline reconciliation factor (unitless)")      # Table S3 'Post-bronchodilator correction' = 0.884 (RSE 1.3%); applies only via FEV1_PBD_ANCHOR

    # ---- log ED50 values (log ug/day of total daily dose) ----------------
    led50_acli_qd  <- 3.60989  ; label("log ED50 for aclidinium q.d. (log ug/day)")        # Table S3 'Log of ED50 for aclidinium q.d.' = 3.61 (RSE 30.2%)
    led50_acli_bid <- 5.3772   ; label("log ED50 for aclidinium b.i.d. (log ug/day)")      # Table S3 'Log of ED50 for aclidinium b.i.d.' = 5.38 (RSE 11.9%)
    led50_bud      <- 5.68377  ; label("log ED50 for budesonide (log ug/day)")             # Table S3 'Log of ED50 for budesonide' = 5.68 (RSE 19.7%)
    led50_form     <- 1.36718  ; label("log ED50 for formoterol (log ug/day)")             # Table S3 'Log of ED50 for formoterol' = 1.37 (RSE 53.0%); arformoterol's ED50 is half of this
    led50_indac    <- 0.963735 ; label("log ED50 for indacaterol (log ug/day)")            # Table S3 'Log of ED50 for indacaterol' = 0.964 (RSE 189.8%)
    led50_glyco    <- 2.20241  ; label("log ED50 for glycopyrronium (log ug/day)")         # Table S3 'Log of ED50 for glycopyrronium' = 2.2 (RSE 15.4%)
    led50_roflu    <- 5.1029   ; label("log ED50 for roflumilast (log ug/day)")            # Table S3 'Log of ED50 for roflumilast' = 5.1 (RSE 19.2%)
    led50_tio      <- 0.961315 ; label("log ED50 for tiotropium HandiHaler (log ug/day)")  # Table S3 'Log of ED50 for tiotropium (Handihaler)' = 0.961 (RSE 52.0%); the Respimat ED50 is derived from this, not estimated
    led50_vil      <- 1.2732   ; label("log ED50 for vilanterol (log ug/day)")             # Table S3 'Log of ED50 for vilanterol' = 1.27 (RSE 80.3%)
    led50_rev      <- 3.78054  ; label("log ED50 for revefenacin (log ug/day)")            # Table S3 'Log of ED50 for revefenacin' = 3.78 (RSE 14.8%)
    led50_olo      <- 0.465982 ; label("log ED50 for olodaterol q.d. (log ug/day)")        # Table S3 'Log of ED50 for olodaterol' = 0.465 (RSE 105.2%)
    led50_bat      <- 3.19838  ; label("log ED50 for batefenterol q.d. (log ug/day)")      # Table S3 'Log of ED50 for batefenterol' = 3.2 (RSE 19.9%)
    led50_ff       <- 2.61137  ; label("log ED50 for fluticasone furoate q.d. (log ug/day)") # Table S3 'Log of ED50 for fluticasone q.d.' = 2.61 (RSE 139.5%); the control stream names it 'log.ED50.FF.QD' and uses it for drug 32, fluticasone furoate

    # ---- Bronchodilator class interaction --------------------------------
    int_labd <- 1.35496 ; label("LABA / LAAC / MABA interaction exponent (unitless)")  # Table S3 'LABD interaction parameter' = 1.35 (RSE 3.2%); combines the class effects as (sum of class^int)^(1/int), which is INFRA-additive for int > 1

    # ---- Covariate effects on baseline -----------------------------------
    e_age_base       <- -0.0143876 ; label("Fractional change in baseline FEV1 per year of age above 63.4")             # Table S3 'Covariate effect of age on baseline' = -0.0144 (RSE 3.8%)
    e_gold_base      <- -0.402018  ; label("Fractional change in baseline FEV1 per GOLD stage above 3")                 # Table S3 'Covariate effect of disease severity on baseline' = -0.402 (RSE 1.8%); a SINGLE linear slope -- the hockey stick of the individual-patient model was linearized away
    e_incl_exac_base <- -0.0213129 ; label("Fractional change in baseline FEV1 in a trial requiring exacerbation history") # Table S3 'Covariate effect of exacerbation history on baseline' = -0.0213 (RSE 43.1%)
    e_sexf_base      <-  0.276389  ; label("Coefficient on the male-coded, 0.671-centred sex column for baseline FEV1")  # Table S3 'Covariate effect of sex on baseline' = 0.276 (RSE 3.4%); applied in model() as (0.329 - SEXF), see the SEXF covariateData entry
    e_smoke_base     <-  0.0365392 ; label("Coefficient on the 0.463-centred current-smoker column for baseline FEV1")   # Table S3 'Covariate effect of smoke on baseline' = 0.0366 (RSE 24.7%)

    # The individual-patient model's hockey-stick GOLD slopes on baseline and
    # its GOLD slope on the disease-progression slope were both dropped in the
    # combined model; the source keeps them in the control stream fixed to zero.
    e_gold_lo_base <- fixed(0) ; label("Second GOLD slope on baseline below stage 3 (fractional per stage)") # $THETA (64) '0 FIX ; B_FLCOPD1'; the hockey stick was linearized away, Section 2.3.1
    e_gold_hi_base <- fixed(0) ; label("Second GOLD slope on baseline above stage 3 (fractional per stage)") # $THETA (65) '0 FIX ; B_FLCOPD2'; as above
    e_gold_dps     <- fixed(0) ; label("GOLD slope on the disease-progression slope (fractional per stage)") # $THETA (67) '0 FIX ; DPSFL_COPD'; dropped as redundant once baseline drives disease progression

    # ---- Covariate effects of the predicted baseline on drug effects ------
    # Linearized from the published aggregated-data model's step function at a
    # baseline of 1.2 L (paper Section 2.3.1). Table S3 keeps the step-function
    # wording 'Effect of predicted baseline <1.2 L on ...' for these rows; the
    # control stream shows the linear form actually fitted.
    e_base_ai <- 0.633718 ; label("Fractional change in the anti-inflammatory effect per litre of baseline above 1.2 L") # Table S3 'Effect of predicted baseline <1.2 L on anti-inflammatory efficacy' = 0.634 (RSE 25.7%); control stream 'RELAI = 1 + THETA(32)*(B - 1.2)'
    e_base_bd <- 0.412128 ; label("Fractional change in the bronchodilator effect per litre of baseline above 1.2 L")    # Table S3 'Effect of predicted baseline <1.2 L on bronchodilator efficacy' = 0.412 (RSE 15.4%); control stream 'RELBD = 1 + THETA(34)*(B - 1.2)'

    # $THETA (31) '0 FIX ; Incr.LOWDS.MHLA', the medical-history adjustment of
    # the lowest admitted disease-severity class, was estimated to zero and the
    # single line that consumed it is commented out of the final control stream
    # ('LOWDSact = LOWDS'). It is not carried here because there is no
    # surviving equation to attach it to, and it does not appear in Table S3.

    # ---- Residual error ---------------------------------------------------
    # Two different error models, selected by DTYPE_AGGREGATED (paper Eq. 10).
    addSd_FEV1  <- sqrt(0.040742)   ; label("Additive residual SD for an aggregated arm-mean record, before the 1/sqrt(NARM) weight (L)") # Table S3 'Variance of the additive residual error of AD' = 0.0407 (RSE 2.8%); $SIGMA (1.) = 0.040742. NONMEM $SIGMA is a VARIANCE, so the SD is sqrt(0.040742) = 0.2018
    powSd_FEV1  <- sqrt(0.00635275) ; label("Residual-error scale for an individual-patient record (L^(1-powExp_FEV1))") # Table S3 'Variance of the power residual error of IPD' = 0.0063 (RSE 2.3%); $SIGMA (2.) = 0.00635275; SD = sqrt(0.00635275) = 0.0797
    powExp_FEV1 <- 0.629635         ; label("Power of the prediction in the individual-patient residual-error model (unitless)") # Table S3 'power error index for IPD error model' = 0.63 (RSE 4.1%); $THETA (68) = 0.629635

    # ==================================================================
    # Random effects. THREE levels are represented. rxode2 draws one level
    # of random effects per solve, so a simulation that needs more than one
    # level must either solve once per level or supply the etas explicitly
    # -- see the vignette. Level semantics:
    #   eta_study_*  one draw per STUDY   (the source's $LEVEL super-ID)
    #   eta_arm_base one draw per ARM     (aggregated records only)
    #   eta*         one draw per SUBJECT (individual records only)
    #
    # Several of these are N(0,1) carriers whose SCALE is an estimated or
    # fixed theta rather than the variance itself; that is how the source
    # shares one distribution across the aggregated and individual parts
    # (paper Section 2.2). Those etas are fixed(1) and are multiplied by
    # their scale parameter inside model().
    # ==================================================================

    eta_study_base ~ fixed(0.0103)  # $OMEGA (1.) '0.0103 FIX ; ISV.Baseline'; Table S3 'ISV variance for the typical baseline' = 0.0103*, carried over from the aggregated-data-only fit per the paper's chosen estimation strategy
    eta_study_dps  ~ fixed(1)       # $OMEGA (2.) '1 FIX ; ISV.Disease.progression.slope'; Table S3 'ISV for the disease progression slope estimated as a fixed effect' = 1*; the actual scale is exp(lisv_cv_dps)
    eta_study_pmx  ~ fixed(0.00097) # $OMEGA (3.) '0.00097 FIX ; ISV.Placebo.Emax'; Table S3 rounds this to 'ISV variance for the placebo effect' = 0.001*
    eta_study_isv_cv_bd   ~ fixed(1)       # $OMEGA (4.) '1 FIX ; ISV.BD.Eff'; Table S3 'ISV variance for the bronchodilator efficacy estimated as fixed effect' = 1*; the actual scale is isv_cv_bd
    eta_study_isv_cv_ai   ~ fixed(1)       # $OMEGA (5.) '1 FIX ; ISV.AI.Eff'; Table S3 'ISV variance for the anti-inflammatory efficacy estimated as fixed effect' = 1*; the actual scale is isv_cv_ai
    eta_study_age  ~ fixed(4)       # $OMEGA (6.) '4 FIX ; ISV.AGE'; Table S3 'Variance for imputation of age' = 4*; used only by the missing-age imputation
    eta_study_smoke ~ fixed(0.17)   # $OMEGA (27.) '0.17 FIX ; ISV.smoke'; Table S3 'Variance for imputation of smoke' = 0.17*; used only by the missing-smoking-fraction imputation
    eta_study_bgics ~ fixed(0.5)    # $OMEGA (22.) '0.5 FIX ; ISV.bg.ics'; Table S3 'Variance for imputation of ICS background treatment' = 0.5*
    eta_study_bglaba ~ fixed(1)     # $OMEGA (23.) '1 FIX ; ISV.bg.laba'; Table S3 'Variance for imputation of LABA background treatment' = 1*
    eta_study_bglaac ~ fixed(0.1)   # $OMEGA (24.) '0.1 FIX ; ISV.bg.laac'; Table S3 'Variance for imputation of LAAC background treatment' = 0.1*

    # Arm-level baseline random effect for aggregated records. The source
    # needed TEN separate N(0,1) etas, one per arm slot, because NONMEM's
    # $LEVEL cannot nest an arm level inside a study level; only one of the
    # ten is ever active for a given record. Collapsed here to a single
    # arm-level eta, which is exact whenever one rxode2 ID is one study arm.
    eta_arm_base ~ fixed(1)         # $OMEGA (9.)-(18.) 'IAV.Baseline.A1 AD' ... ten BLOCK(1) 1 FIX records; the source's IOV-like method, paper Section 2.2

    # Subject-level effects, active on individual-patient records only.
    etabase        ~ fixed(1)       # $OMEGA (8.) 'BLOCK(1) 1 FIX ; IIV.Baseline.A1 IPD'; the actual scale is cv_base = 0.236954
    etadps         ~ 0.0607811      # $OMEGA (25.) '0.0607811 ; IIV DPS'; Table S3 'IIV variance of disease progression slope of IPD' = 0.0608 (RSE 6.3%)
    etaeffref_vil  ~ fixed(0)       # $OMEGA (26.) '0 FIX ; IIV VIL eff'; estimated to zero in the combined model, unlike the 0.00698 of the individual-patient-only model
    etapowSd_FEV1  ~ 0.128503       # $OMEGA (7.) '0.128503 ; Eta.on.eps for IPD'; Table S3 'IIV variance of residual error of IPD' = 0.129 (RSE 4.7%)

    # Scale parameters for the N(0,1) random-effect carriers above.
    cv_base     <-  0.236954 ; label("Between-subject log-SD of baseline FEV1 (unitless)")                       # Table S3 'IIV CV for baseline' = 0.237 (RSE 1.4%); $THETA (63)
    lisv_cv_dps <- fixed(-0.265) ; label("log between-study SD of the disease-progression slope (log L/year)")   # Table S3 'Log of ISV CV for the disease progression slope' = -0.265*; $THETA (27)
    isv_cv_bd   <- fixed(0.193)  ; label("Between-study SD of the fractional bronchodilator effect (unitless)")  # Table S3 'ISV CV for bronchodilator efficacy' = 0.193*; $THETA (48)
    isv_cv_ai   <- fixed(0.422)  ; label("Between-study SD of the fractional anti-inflammatory effect (unitless)") # Table S3 'ISV CV for anti-inflammatory efficacy' = 0.422*; $THETA (29)
  })

  model({
    # ==================================================================
    # 0. Missing-covariate imputation.
    #
    # The source carries sentinel-coded missing values in the meta-analysis
    # data set and imputes them inside $PRED with FIXED-coefficient
    # regressions inherited from the upstream Llanos-Paez 2023 analysis
    # (paper Section 2.4). The regression coefficients are hardcoded in the
    # control stream and appear in no table; only the random-effect
    # variances are in Table S3. Every branch below is inert when the
    # covariate is observed, so a complete data set never reaches them.
    #
    # The source also imputes three PRIOR-medical-history fractions
    # (ETA(19)-(21), Table S3 rows 'Variance for imputation of ICS / LABA /
    # LAAC medication history'). Those are DEAD CODE in the final model:
    # their only consumer is the lowest-disease-severity adjustment whose
    # coefficient THETA(31) was estimated to zero and whose line is
    # commented out ('LOWDSact = LOWDS'). They are therefore not encoded
    # here; the omission changes no prediction and is noted in the vignette.
    # ==================================================================

    gold_low <- DIS_COPD_GOLD_LOW
    gold_high <- DIS_COPD_GOLD_HIGH

    # Missing arm-mean age: a linear regression on the admitted severity
    # range and the exacerbation-history entry criterion.
    age_imp <- 64 * (1 + 0.017 * (gold_low - 2)) * (1 + 0.017 * (gold_high - 4)) *
      (1 - 0.007 * INCL_EXAC_REQUIRED) + eta_study_age
    age_use <- (AGE > -1) * AGE + (AGE < -1) * age_imp

    # Missing current-smoker fraction: a logistic regression on imputed age.
    lgt_smoke <- -0.165 * (1 + 0.706 * (age_use - 63.7)) + eta_study_smoke
    smoke_imp <- expit(lgt_smoke)
    smoke_use <- (SMOKE > -1) * SMOKE + (SMOKE < -1) * smoke_imp

    # Missing GOLD stage on an individual record is replaced by the median
    # stage 3; an aggregated record uses the midpoint of the admitted range.
    gold_ipd  <- (DIS_COPD_GOLD > 0) * DIS_COPD_GOLD + (DIS_COPD_GOLD < 0) * 3
    gold_use  <- DTYPE_AGGREGATED * (gold_low + gold_high) / 2 +
      (1 - DTYPE_AGGREGATED) * gold_ipd

    # Missing background-therapy fractions. Sentinel 9999 means the class was
    # not used (treated as zero); sentinel 7777 means it was used but the
    # fraction was not reported, and is imputed by a logistic regression on
    # the admitted severity range and the entry criterion. The run-in and
    # maintenance branches share both the regression and the random effect.
    lgt_bgics  <- 0.36 * (1 + 1.6 * (gold_low - 2)) * (1 + 1.6 * (gold_high - 4)) *
      (1 + 0.12 * INCL_EXAC_REQUIRED) + eta_study_bgics
    bgics_imp  <- 100 * expit(lgt_bgics)

    lgt_bglaba <- 0.49 * (1 - 0.84 * (gold_low - 2)) * (1 - 0.84 * (gold_high - 4)) *
      (1 + 0.96 * INCL_EXAC_REQUIRED) + eta_study_bglaba
    bglaba_imp <- 100 * expit(lgt_bglaba)

    # The control stream notes that the entry criterion 'is not supported in
    # regression' for the anticholinergic class, so it carries no INCL term.
    lgt_bglaac <- -5 * (1 - 0.67 * (gold_low - 2)) * (1 - 0.67 * (gold_high - 4)) +
      eta_study_bglaac
    bglaac_imp <- 100 * expit(lgt_bglaac)

    ics_runin  <- (BGTHER_ICS_RUNIN_PCT < 1000) * BGTHER_ICS_RUNIN_PCT +
      (BGTHER_ICS_RUNIN_PCT > 7000) * (BGTHER_ICS_RUNIN_PCT < 8000) * bgics_imp
    laba_runin <- (BGTHER_LABA_RUNIN_PCT < 1000) * BGTHER_LABA_RUNIN_PCT +
      (BGTHER_LABA_RUNIN_PCT > 7000) * (BGTHER_LABA_RUNIN_PCT < 8000) * bglaba_imp
    laac_runin <- (BGTHER_LAAC_RUNIN_PCT < 1000) * BGTHER_LAAC_RUNIN_PCT +
      (BGTHER_LAAC_RUNIN_PCT > 7000) * (BGTHER_LAAC_RUNIN_PCT < 8000) * bglaac_imp

    ics_maint  <- (BGTHER_ICS_MAINT_PCT < 1000) * BGTHER_ICS_MAINT_PCT +
      (BGTHER_ICS_MAINT_PCT > 7000) * (BGTHER_ICS_MAINT_PCT < 8000) * bgics_imp
    laba_maint <- (BGTHER_LABA_MAINT_PCT < 1000) * BGTHER_LABA_MAINT_PCT +
      (BGTHER_LABA_MAINT_PCT > 7000) * (BGTHER_LABA_MAINT_PCT < 8000) * bglaba_imp
    laac_maint <- (BGTHER_LAAC_MAINT_PCT < 1000) * BGTHER_LAAC_MAINT_PCT +
      (BGTHER_LAAC_MAINT_PCT > 7000) * (BGTHER_LAAC_MAINT_PCT < 8000) * bglaac_imp

    # ==================================================================
    # 1. Baseline FEV1
    # ==================================================================
    # Four multiplicative, linear-in-deviation covariate effects. The sex
    # term is written on the source's male-coded, 0.671-centred column,
    # which on the canonical female indicator is (0.329 - SEXF).
    tvb <- base *
      (1 + e_age_base * (age_use - 63.4)) *
      (1 + e_gold_base * (gold_use - 3)) *
      (1 + e_incl_exac_base * INCL_EXAC_REQUIRED) *
      (1 + e_sexf_base * (0.329 - SEXF)) *
      (1 + e_smoke_base * (smoke_use - 0.463)) *
      # Both hockey-stick slopes are fixed to zero in the combined model.
      (1 + (e_gold_lo_base * (gold_use <= 3) + e_gold_hi_base * (gold_use > 3)) * (gold_use - 3))

    # Post-bronchodilator baseline: a shift by the median absolute
    # reversibility of 0.18 L (paper Section 2.4 note).
    tvb_pbd <- tvb + 0.18

    # Study-level baseline random effect, shared by aggregated and
    # individual records.
    study_b <- exp(eta_study_base)

    # An INDIVIDUAL record's baseline is log-normal about the study value.
    # An AGGREGATED record's baseline is the MEAN of that log-normal over
    # NARM patients, which the paper approximates by a normal with the
    # log-normal's mean and its variance divided by NARM (paper Equations
    # 6-8). Using the log-normal itself for the arm mean is what produced
    # the 8-10% aggregation bias the paper set out to remove (Section 3.2).
    ln_mean  <- exp(cv_base * cv_base / 2)
    ln_sd    <- sqrt(exp(cv_base * cv_base) - 1) * ln_mean
    b_ad     <- tvb * study_b * ln_mean + tvb * study_b * ln_sd / sqrt(NARM) * eta_arm_base
    b_ipd    <- tvb * study_b * exp(etabase * cv_base)
    b        <- DTYPE_AGGREGATED * b_ad + (1 - DTYPE_AGGREGATED) * b_ipd

    b_pbd_ad  <- tvb_pbd * study_b * ln_mean +
      tvb_pbd * study_b * ln_sd / sqrt(NARM) * eta_arm_base
    b_pbd_ipd <- tvb_pbd * study_b * exp(etabase * cv_base)
    b_pbd     <- DTYPE_AGGREGATED * b_pbd_ad + (1 - DTYPE_AGGREGATED) * b_pbd_ipd

    # ==================================================================
    # 2. Effect of the predicted baseline on the drug effects
    # ==================================================================
    rel_ai <- 1 + e_base_ai * (b - 1.2)
    rel_bd <- 1 + e_base_bd * (b - 1.2)

    # ==================================================================
    # 3. Disease progression -- proportional to the baseline, normalised
    #    to 1.2 L, with a log-normal study-level effect and an additive
    #    subject-level effect on individual records only.
    # ==================================================================
    dps_i <- dps * b / 1.2 * exp(eta_study_dps * exp(lisv_cv_dps)) *
      (1 + e_gold_dps * (gold_use - 3)) +
      (1 - DTYPE_AGGREGATED) * etadps
    dp <- dps_i * t / 52

    # ==================================================================
    # 4. Placebo -- an immediate constant for t > 0
    # ==================================================================
    plac <- (pmx + eta_study_pmx) * (t > 0)

    # ==================================================================
    # 5. Effect-onset time courses
    #
    # Each is a rise from the fraction of the arm already on that class at
    # run-in toward 1. Only the corticosteroids, the PDE4 inhibitors, and
    # once-daily indacaterol / glycopyrronium / tiotropium carry one.
    # ==================================================================
    onset_ai_rate <- exp(lonset_ai)
    tc_ics  <- ics_runin / 100 + (1 - ics_runin / 100) * (1 - exp(-onset_ai_rate * t))
    # No PDE4 inhibitor was present at run-in anywhere in the data set, so
    # the source drops the run-in offset from this one time course.
    tc_pde4 <- 1 - exp(-onset_ai_rate * t)
    tc_laba <- laba_runin / 100 + (1 - laba_runin / 100) * (1 - exp(-onset_qdba * t))
    tc_laac <- laac_runin / 100 + (1 - laac_runin / 100) * (1 - exp(-onset_qdac * t))

    # ==================================================================
    # 6. Per-drug maximum effects and dose-responses
    #
    # For a dose-response drug, Emax = EffRef / RefDose * (ED50 + RefDose)
    # and the effect is the hyperbolic Emax in the per-arm total daily
    # dose. For a constant-effect drug the effect is a single estimate,
    # keyed on the dose column being non-zero. Every term collapses to
    # zero for an arm that did not receive the drug.
    # ==================================================================
    ed50_acli_qd  <- exp(led50_acli_qd)
    ed50_acli_bid <- exp(led50_acli_bid)
    ed50_bud      <- exp(led50_bud)
    ed50_form     <- exp(led50_form)
    ed50_arform   <- ed50_form / 2
    ed50_indac    <- exp(led50_indac)
    ed50_glyco    <- exp(led50_glyco)
    ed50_roflu    <- exp(led50_roflu)
    ed50_tio      <- exp(led50_tio)
    ed50_vil      <- exp(led50_vil)
    ed50_rev      <- exp(led50_rev)
    ed50_olo      <- exp(led50_olo)
    ed50_bat      <- exp(led50_bat)
    ed50_ff       <- exp(led50_ff)

    emax_acli_qd  <- effref_acli_qd  / 200 * (ed50_acli_qd + 200)
    emax_acli_bid <- effref_acli_bid / 800 * (ed50_acli_bid + 800)
    emax_bud      <- effref_bud      / 320 * (ed50_bud + 320)
    emax_form     <- effref_form     /  18 * (ed50_form + 18)
    # Arformoterol borrows formoterol's Emax outright.
    emax_arform   <- emax_form
    emax_indac    <- effref_indac    /  75 * (ed50_indac + 75)
    emax_glyco    <- effref_glyco    / 100 * (ed50_glyco + 100)
    emax_roflu    <- effref_roflu    / 500 * (ed50_roflu + 500)
    emax_tio_dpi  <- effref_tio_dpi  /  18 * (ed50_tio + 18)
    # The Respimat ED50 is solved from the HandiHaler Emax and the Respimat
    # reference efficacy rather than estimated separately.
    ed50_tio_smi  <- 5 * (emax_tio_dpi / effref_tio_smi - 1)
    emax_tio_smi  <- effref_tio_smi  /   5 * (ed50_tio_smi + 5)
    emax_vil      <- (effref_vil + (1 - DTYPE_AGGREGATED) * etaeffref_vil) / 25 * (ed50_vil + 25)
    emax_rev      <- effref_rev      / 175 * (ed50_rev + 175)
    emax_olo_qd   <- effref_olo_qd   /   5 * (ed50_olo + 5)
    emax_bat_qd   <- effref_bat_qd   / 400 * (ed50_bat + 400)
    emax_ff       <- effref_ff       / 100 * (ed50_ff + 100)

    # Open-label tiotropium keeps the blinded ED50 and scales the efficacy.
    emax_tio_ol   <- rel_tio_ol * effref_tio_dpi / 18 * (ed50_tio + 18)

    # ---- Long-acting beta-agonists ---------------------------------------
    de_arform <- emax_arform * CONMED_ARFORMOTEROL_DOSE / (CONMED_ARFORMOTEROL_DOSE + ed50_arform)
    de_form   <- emax_form   * CONMED_FORMOTEROL_DOSE   / (CONMED_FORMOTEROL_DOSE + ed50_form)
    de_indac  <- emax_indac  * tc_laba * CONMED_INDACATEROL_DOSE / (CONMED_INDACATEROL_DOSE + ed50_indac)
    de_salm   <- eff_salm    * (CONMED_SALMETEROL_DOSE > 0)
    de_vil    <- emax_vil    * CONMED_VILANTEROL_DOSE / (CONMED_VILANTEROL_DOSE + ed50_vil)
    de_olo    <- (1 - FORM_OLODATEROL_BID) * emax_olo_qd * CONMED_OLODATEROL_DOSE /
      (CONMED_OLODATEROL_DOSE + ed50_olo) +
      FORM_OLODATEROL_BID * eff_olo_bid * (CONMED_OLODATEROL_DOSE > 0)

    eff_laba <- de_arform + de_form + de_indac + de_salm + de_vil + de_olo

    # ---- Long-acting anticholinergics ------------------------------------
    de_acli <- (1 - FORM_ACLIDINIUM_BID) * emax_acli_qd * CONMED_ACLIDINIUM_DOSE /
      (CONMED_ACLIDINIUM_DOSE + ed50_acli_qd) +
      FORM_ACLIDINIUM_BID * emax_acli_bid * CONMED_ACLIDINIUM_DOSE /
      (CONMED_ACLIDINIUM_DOSE + ed50_acli_bid)
    de_glyco <- emax_glyco * tc_laac * CONMED_GLYCOPYRRONIUM_DOSE /
      (CONMED_GLYCOPYRRONIUM_DOSE + ed50_glyco)
    # Device and blinding select among three tiotropium variants.
    de_tio_dpi <- emax_tio_dpi * tc_laac * CONMED_TIOTROPIUM_DOSE / (CONMED_TIOTROPIUM_DOSE + ed50_tio)
    de_tio_smi <- emax_tio_smi * tc_laac * CONMED_TIOTROPIUM_DOSE / (CONMED_TIOTROPIUM_DOSE + ed50_tio_smi)
    de_tio_ol  <- emax_tio_ol  * tc_laac * CONMED_TIOTROPIUM_DOSE / (CONMED_TIOTROPIUM_DOSE + ed50_tio)
    de_tio <- FORM_TIOTROPIUM_OPENLABEL * de_tio_ol +
      (1 - FORM_TIOTROPIUM_OPENLABEL) * (FORM_TIOTROPIUM_SMI * de_tio_smi +
                                           (1 - FORM_TIOTROPIUM_SMI) * de_tio_dpi)
    de_ume  <- eff_ume       * (CONMED_UMECLIDINIUM_DOSE > 0)
    de_gsk  <- eff_gsk233705 * (CONMED_GSK233705_DOSE > 0)
    de_bea  <- eff_bea2180   * (CONMED_BEA2180_DOSE > 0)
    de_rev  <- emax_rev * CONMED_REVEFENACIN_DOSE / (CONMED_REVEFENACIN_DOSE + ed50_rev)

    eff_laac <- de_acli + de_glyco + de_tio + de_ume + de_gsk + de_bea + de_rev

    # ---- Dual-pharmacology bronchodilator (its own class) ----------------
    de_bat <- (1 - FORM_BATEFENTEROL_BID) * emax_bat_qd * CONMED_BATEFENTEROL_DOSE /
      (CONMED_BATEFENTEROL_DOSE + ed50_bat) +
      FORM_BATEFENTEROL_BID * eff_bat_bid * (CONMED_BATEFENTEROL_DOSE > 0)
    eff_maba <- de_bat

    # ---- Inhaled corticosteroids -----------------------------------------
    de_beclo     <- eff_beclo     * tc_ics * (CONMED_BECLOMETHASONE_DOSE > 0)
    de_bud       <- emax_bud      * tc_ics * CONMED_BUDESONIDE_DOSE / (CONMED_BUDESONIDE_DOSE + ed50_bud)
    de_flutiprop <- eff_flutiprop * tc_ics * (CONMED_FLUTICASONEPROPIONATE_DOSE > 0)
    de_mome <- ((1 - FORM_MOMETASONE_BID) * eff_mome_qd +
                  FORM_MOMETASONE_BID * rel_mome_bid * eff_mome_qd) *
      tc_ics * (CONMED_MOMETASONE_DOSE > 0)
    # Fluticasone furoate carries no onset time course in the source.
    de_ff <- emax_ff * CONMED_FLUTICASONEFUROATE_DOSE / (CONMED_FLUTICASONEFUROATE_DOSE + ed50_ff)

    eff_ics <- de_beclo + de_bud + de_flutiprop + de_mome + de_ff

    # ---- Non-steroid anti-inflammatories ---------------------------------
    de_cil   <- eff_cil      * tc_pde4 * (CONMED_CILOMILAST_DOSE > 0)
    de_roflu <- emax_roflu   * tc_pde4 * CONMED_ROFLUMILAST_DOSE / (CONMED_ROFLUMILAST_DOSE + ed50_roflu)
    de_azd   <- eff_azd9668  * (CONMED_AZD9668_DOSE > 0)
    de_ph    <- eff_ph797804 * (CONMED_PH797804_DOSE > 0)

    eff_ai_other <- de_cil + de_roflu + de_azd + de_ph

    # A study enrolling only oral-corticosteroid non-responders sees no
    # inhaled-corticosteroid effect, but retains the non-steroid effects.
    rel_cs <- 1 - OCS_NONRESPONDER

    # ==================================================================
    # 7. Background (non-randomized) therapy
    #
    # Valued at fluticasone propionate b.i.d. for corticosteroids,
    # salmeterol for beta-agonists, and tiotropium at its 18 ug/day
    # reference dose for anticholinergics, each scaled by the fraction of
    # the arm on that class. The run-in fractions apply at t = 0 and the
    # maintenance fractions thereafter. Background therapy is present in
    # placebo arms too, so it is NOT gated on active treatment.
    # ==================================================================
    de_tio_ref <- emax_tio_dpi * 18 / (18 + ed50_tio)
    bg_frac_ics  <- (t > 0) * ics_maint  + (t <= 0) * ics_runin
    bg_frac_laba <- (t > 0) * laba_maint + (t <= 0) * laba_runin
    bg_frac_laac <- (t > 0) * laac_maint + (t <= 0) * laac_runin

    bg_ai   <- eff_flutiprop * bg_frac_ics / 100 * rel_cs
    bg_laba <- eff_salm      * bg_frac_laba / 100
    bg_laac <- de_tio_ref    * bg_frac_laac / 100

    # ==================================================================
    # 8. Class totals, interaction, and study-level effect variability
    #
    # Randomized drug effects are expressed only after randomization.
    # ==================================================================
    treat_on <- (t > 0)

    tot_laba <- bg_laba + treat_on * eff_laba
    tot_laac <- bg_laac + treat_on * eff_laac
    tot_maba <- treat_on * eff_maba

    # Infra-additive combination of the three bronchodilator classes. The
    # exponent is estimated at 1.35, so two classes together deliver less
    # than the sum of their separate effects.
    eff_bd_int <- (tot_laba^int_labd + tot_laac^int_labd + tot_maba^int_labd)^(1 / int_labd)
    eff_bd <- eff_bd_int * rel_bd * (1 + eta_study_isv_cv_bd * isv_cv_bd)

    eff_ai_total <- bg_ai + treat_on * (eff_ics * rel_cs + eff_ai_other)
    eff_ai <- eff_ai_total * rel_ai * (1 + eta_study_isv_cv_ai * isv_cv_ai)

    # ==================================================================
    # 9. Prediction
    #
    # A pre-bronchodilator record sees the full bronchodilator effect; a
    # post-bronchodilator record sees a raised baseline and only the
    # fraction rel_postbd of it. The reconciliation term is zero unless
    # FEV1_PBD_ANCHOR is supplied -- see that covariate's entry.
    # ==================================================================
    drug_pre  <- eff_ai + eff_bd
    drug_post <- eff_ai + eff_bd * rel_postbd
    postbd_corr <- (1 - postbd_recon) * FEV1_PBD_ANCHOR

    fev1_pre  <- b     - dp + plac + drug_pre + postbd_corr
    fev1_post <- b_pbd - dp + plac + drug_post

    FEV1 <- (1 - MEAS_POSTBD) * fev1_pre + MEAS_POSTBD * fev1_post

    # ==================================================================
    # 10. Residual error -- two models selected by record type.
    #
    # An aggregated arm mean of NARM patients carries an additive error
    # shrunk by 1/sqrt(NARM); an individual record carries a power error
    # whose magnitude has its own log-normal between-subject variability.
    # The two are combined into one SD so a single endpoint can serve both
    # record types.
    # ==================================================================
    sd_ad  <- addSd_FEV1 / sqrt(NARM)
    sd_ipd <- powSd_FEV1 * FEV1^powExp_FEV1 * exp(etapowSd_FEV1)
    sdFEV1 <- DTYPE_AGGREGATED * sd_ad + (1 - DTYPE_AGGREGATED) * sd_ipd

    FEV1 ~ add(sdFEV1)
  })
}
