# Population PK model for buprenorphine (BPN) after intravenous, sublingual,
# and subcutaneous CAM2038 weekly / monthly depot administration in healthy
# participants and participants with opioid use disorder (Bjornsson 2023,
# Clin Pharmacokinet 62(9):1429-1443; doi:10.1007/s40262-023-01288-6).

Bjornsson_2023_buprenorphine <- function() {
  description <- paste(
    "Three-compartment population pharmacokinetic model for buprenorphine",
    "(BPN) covering four routes of administration in a single jointly-fitted",
    "model: intravenous BPN, sublingual (SL) BPN tablets, and subcutaneous",
    "injection of the CAM2038 FluidCrystal extended-release depot in its",
    "weekly (Q1W) and monthly (Q4W) formulations. Disposition is",
    "three-compartment with first-order elimination from central. Each",
    "extravascular route absorbs through two parallel depot pathways, one",
    "fast and one slow, with the dose split between them by an estimated",
    "logit-transformed fraction: SL BPN uses a sequential zero-order (D_SL1)",
    "then first-order (ka_SL1) pathway with a lag time plus a purely",
    "first-order pathway (ka_SL2); CAM2038 Q1W uses a sequential zero-order",
    "(D_q1w1) then first-order (ka_q1w1) pathway plus a slow first-order",
    "pathway (ka_q1w2, absorption half-life 123 h); CAM2038 Q4W uses two",
    "parallel first-order pathways (ka_q4w1 and ka_q4w2, absorption half-life",
    "418 h). Absorption from the CAM2038 depot is rate-limiting for BPN",
    "elimination (flip-flop kinetics). Bioavailability is 1 for IV and fixed",
    "to 1 for both CAM2038 formulations, and dose-dependent for SL BPN",
    "(F_SL = 0.14 * (Dose / 16 mg)^-0.371, i.e. 18.1 / 14.0 / 12.0 percent at",
    "8 / 16 / 24 mg). Covariate effects retained in the final model are age",
    "and body weight on CL, opioid-use-disorder status on Vc, and",
    "opioid-use-disorder status and female sex on the CAM2038 Q1W fast-pathway",
    "dose fraction; a thigh injection site routes the entire CAM2038 Q1W dose",
    "through the slow pathway. Residual error is additive on log-transformed",
    "concentrations, i.e. proportional on the linear scale.",
    sep = " "
  )
  reference <- paste(
    "Bjornsson M, Acharya C, Strandgarden K, Tiberg F (2023). Population",
    "pharmacokinetic analysis supports initiation treatment and bridging from",
    "sublingual buprenorphine to subcutaneous administration of a",
    "buprenorphine depot (CAM2038) in the treatment of opioid use disorder.",
    "Clinical Pharmacokinetics 62(9):1429-1443.",
    "doi:10.1007/s40262-023-01288-6.",
    sep = " "
  )
  vignette <- "Bjornsson_2023_buprenorphine"

  # The six absorption depots are the source paper's own Figure 1 compartment
  # labels (CMTSL1, CMTSL2, CMTq1w1, CMTq1w2, CMTq4w1, CMTq4w2). They are
  # declared as paper-specific rather than registered as new canonical
  # compartments because four of them name pathways of two proprietary
  # CAM2038 formulations and do not generalise the way the registered
  # `depot_oral` / `depot_im` parallel-route depots do. Keeping the
  # route-and-pathway names (rather than the blessed generic `depot1` ...
  # `depot6`) matters here: a user dosing this four-route model must be able
  # to tell which depot pair belongs to which formulation, and which member
  # of a pair is the fast versus the slow pathway.
  paper_specific_compartments <- c(
    "depot_sl1", "depot_sl2",
    "depot_q1w1", "depot_q1w2",
    "depot_q4w1", "depot_q4w2"
  )

  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Declared explicitly: buildModelDb()'s auto-detection only recognises
  # compartments literally named `depot` / `central`, so it would label this
  # four-route model as central-dosed only. Every extravascular route is dosed
  # into BOTH members of its depot pair with the same amount -- the f()
  # statements in model() split the administered dose between the parallel fast
  # and slow pathways. Dose records for `depot_sl1` and `depot_q1w1` must carry
  # rate = -2 so rxode2 applies the modelled zero-order duration dur().
  dosing <- c(
    "central",
    "depot_sl1", "depot_sl2",
    "depot_q1w1", "depot_q1w2",
    "depot_q4w1", "depot_q4w2"
  )

  # Every ODE state holds an amount of buprenorphine in mg (doses are in mg and
  # the observation Cc converts central/vc from mg/L to the published ng/mL).
  # The six depots are the sublingual and subcutaneous administration sites
  # (Figure 1's CMTSL1 / CMTSL2 / CMTq1w1 / CMTq1w2 / CMTq4w1 / CMTq4w2), so
  # they take the `administration site` specimen term -- none of them was ever
  # sampled. Only the plasma concentration derived from `central` was assayed
  # (Methods 2.1, LC-MS/MS, LLOQ 0.0250 ng/mL).
  compartmentData <- list(
    depot_sl1   = list(analyte = "BPN", units = "mg", specimen = "administration site", verified = FALSE),
    depot_sl2   = list(analyte = "BPN", units = "mg", specimen = "administration site", verified = FALSE),
    depot_q1w1  = list(analyte = "BPN", units = "mg", specimen = "administration site", verified = FALSE),
    depot_q1w2  = list(analyte = "BPN", units = "mg", specimen = "administration site", verified = FALSE),
    depot_q4w1  = list(analyte = "BPN", units = "mg", specimen = "administration site", verified = FALSE),
    depot_q4w2  = list(analyte = "BPN", units = "mg", specimen = "administration site", verified = FALSE),
    central     = list(analyte = "BPN", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "BPN", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral2 = list(analyte = "BPN", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    AGE = list(
      description        = "Age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power-form effect on CL normalised to the pooled cohort median of 35",
        "years (Table 1, All N = 236 column): CL = 52.1 * (AGE / 35)^-0.233.",
        "Bjornsson 2023 Results 3.1.3 verifies the effect size: CL is 60.8 L/h",
        "for an 18-year-old and 45.1 L/h for a 65-year-old (a 26 percent",
        "decrease). Observed range 18-65 years."
      ),
      source_name        = "Age"
    ),
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power-form effect on CL normalised to the pooled cohort median of",
        "72.4 kg: CL = 52.1 * (WT / 72.4)^0.413. The reference value 72.4 kg",
        "is the body weight of the typical individual used for the paper's",
        "simulations (Methods 2.4); Table 1 reports the pooled median as 72.2",
        "kg for the initial model dataset. Bjornsson 2023 Results 3.1.3",
        "verifies the effect size: CL is 44.7 L/h at 50 kg and 59.5 L/h at 100",
        "kg (a 33 percent increase). Observed range 50.4-127 kg. Body weight",
        "was also tested on F_SL and Vc but was removed from the final model",
        "for poor precision (RSE > 50 percent)."
      ),
      source_name        = "WT"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male; the most common category at 61 percent of the pooled cohort, and therefore the paper's reference for the fractional-difference categorical covariate model of Eq. 5)",
      notes              = paste(
        "1 = female, 0 = male. Additive effect of +0.576 on the logit of the",
        "CAM2038 Q1W fast-pathway dose fraction Fq1w1 (Table 3, 'Sex covariate",
        "on Fq1w1', footnote b 'Females versus males'), giving Fq1w1 = 59.8",
        "percent in women versus 45.5 percent in men. Bjornsson 2023 Results",
        "3.1.3 confirms the direction: 'Higher Fq1w1 was shown in women than",
        "in men'. Overall BPN exposure (AUC) is unaffected because Fq1w1 only",
        "redistributes dose between the fast and slow absorption pathways."
      ),
      source_name        = "Sex (Table 1)"
    ),
    DIS_OUD = list(
      description        = "Opioid use disorder patient indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy participant; the most common category at 62 percent of the pooled cohort, and therefore the paper's reference for the fractional-difference categorical covariate model of Eq. 5)",
      notes              = paste(
        "1 = participant with opioid use disorder (the Phase 2 cohorts),",
        "0 = healthy participant (the Phase 1 cohorts). The paper calls this",
        "covariate 'population' and reports both coefficients as 'Patients",
        "with OUD versus healthy volunteers' (Table 3 footnote a). Two effects",
        "are retained in the final model: a fractional difference of +2.69 on",
        "Vc (64.3 L in healthy participants, 64.3 * (1 + 2.69) = 237 L in",
        "participants with OUD, matching the values quoted in Results 3.1.3),",
        "and an additive -0.671 on the logit of the CAM2038 Q1W fast-pathway",
        "dose fraction Fq1w1 (45.5 percent in healthy participants versus 29.9",
        "percent in participants with OUD). Population was also tested on",
        "Fq4w1 but removed for poor precision (RSE > 50 percent)."
      ),
      source_name        = "Population (healthy participants vs participants with OUD)"
    ),
    INJSITE_THIGH = list(
      description        = "Thigh subcutaneous injection-site indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (buttock, abdomen, or upper arm; Bjornsson 2023 Discussion 4.1 reports no differences in absorption parameters between these three sites, and the paper's simulations presume buttock injection)",
      notes              = paste(
        "1 = the CAM2038 Q1W dose was injected into the thigh. Bjornsson 2023",
        "Results 3.1.1 reports that for thigh injections the estimated fast",
        "pathway dose fraction Fq1w1 was close to 0 and was 'subsequently",
        "fixed to 0', so the entire dose is absorbed through the slow pathway",
        "(depot_q1w2, ka_q1w2 = 0.00565 1/h). Encoded here as a multiplicative",
        "gate (1 - INJSITE_THIGH) on Fq1w1 rather than as an estimated",
        "coefficient, because the paper fixed the thigh value at exactly 0.",
        "Applies only to CAM2038 Q1W: the injection-site effect was assessed",
        "univariately on the CAM2038 Q1W absorption parameters (Methods",
        "2.3.1), and Discussion 4.1 states that no injection-site differences",
        "in overall exposure are expected for CAM2038 Q4W. Per-dose-record",
        "covariate."
      ),
      source_name        = "Injection site (buttock / abdomen / thigh / upper arm)"
    ),
    DOSE_BPN_SL_MG = list(
      description        = "Administered sublingual buprenorphine dose",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Sublingual BPN dose in mg, supplied as a data column so the",
        "dose-dependent SL bioavailability can be evaluated inside model():",
        "F_SL = 0.14 * (DOSE_BPN_SL_MG / 16)^-0.371 (Bjornsson 2023 Results",
        "3.1.3). The reference dose is 16 mg, the typical daily maintenance SL",
        "dose. The paper verifies the equation against its own point values:",
        "F_SL is 18.1, 14.0, and 12.0 percent at 8, 16, and 24 mg. Set this",
        "column to the SL dose for records in an SL treatment period; the",
        "value is unused (and may be left at the 16 mg reference) for IV or",
        "CAM2038-only simulations because F_SL scales only the two SL depots.",
        "Not a covariate on any CAM2038 parameter -- BPN PK was",
        "dose-proportional for CAM2038 Q1W (8-32 mg) and Q4W (64-192 mg)."
      ),
      source_name        = "SL BPN dose (mg)"
    )
  )

  covariatesDataExcluded <- list(
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = paste(
        "Not tested. Bjornsson 2023 Results 3.1.2: 'Given the strong",
        "correlation between body weight and body mass index, only body weight",
        "was included in the covariate analysis.' Pooled median 24.4 kg/m^2",
        "(range 16.8-34.9; Table 1)."
      )
    ),
    RACE_BLACK = list(
      description = "Black or African American race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened in the stepwise covariate model but not retained. Race was",
        "collapsed to Black or African American versus all other races",
        "combined because few participants were Asian, Mixed, or Other",
        "(Results 3.1.2); Results 3.1.3 reports 'No additional effects of",
        "population, body weight, sex, age, race, or creatinine clearance on",
        "the model parameters were found'. 22 percent of the pooled cohort",
        "(Table 1)."
      )
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units       = "mL/min",
      type        = "continuous",
      notes       = paste(
        "Screened but not retained (Results 3.1.3). Discussion 4.1: 'CLCR was",
        "not identified as a predictor for CL in individuals with normal or",
        "mild renal function impairment.' Approximately 85 percent of",
        "participants had normal renal function (>= 90 mL/min) and 15 percent",
        "mild impairment (60 to < 90 mL/min); no participant had moderate or",
        "severe impairment. Pooled median 112 mL/min (range 65.1-209; Table",
        "1)."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 252L,
    n_studies      = 4L,
    age_range      = "18-65 years",
    age_median     = "35 years (initial model dataset, N = 236); 39 years (CAM2038 Q4W model-update dataset, N = 96)",
    weight_range   = "50.4-127 kg",
    weight_median  = "72.2 kg (initial model dataset); 75.6 kg (CAM2038 Q4W model-update dataset)",
    sex_female_pct = 39.0,
    race_ethnicity = c(
      White                       = 66,
      `Black or African American` = 22,
      Asian                       = 6,
      Mixed                       = 2,
      Other                       = 5
    ),
    disease_state  = paste(
      "Pooled cohort of healthy participants (147 / 236 = 62 percent of the",
      "initial-model dataset; Phase 1 Trials 1 and 2) and participants with",
      "opioid use disorder (89 / 236 = 38 percent; Phase 2 Trials 3 and 4)."
    ),
    renal_function = paste(
      "Approximately 85 percent normal renal function (creatinine clearance",
      ">= 90 mL/min) and approximately 15 percent mild impairment (60 to < 90",
      "mL/min). No participant had moderate or severe renal impairment",
      "(Methods 2.1)."
    ),
    dose_range     = paste(
      "Intravenous BPN; sublingual BPN tablets (2-32 mg daily); subcutaneous",
      "CAM2038 Q1W 8, 16, 24, or 32 mg weekly; subcutaneous CAM2038 Q4W 64,",
      "96, 128, 160, or 192 mg monthly. Both single-dose and repeat-dose",
      "regimens (Methods 2.1, Supplementary Appendix A)."
    ),
    bmi_range      = "16.8-34.9 kg/m^2 (median 24.4)",
    crcl_range     = "65.1-237 mL/min (median 112-113)",
    regions        = "Not reported by region; four company-sponsored trials (ISRCTN41550730, ISRCTN24987553, NCT02611752, NCT02710526)",
    notes          = paste(
      "The final model analysis included 252 individuals and 10,658 BPN plasma",
      "concentration observations (Results 3.1.1). This is the union of the",
      "initial population PK model dataset (236 participants, 10,260",
      "observations, Trials 1-4) and the 16 additional participants (398",
      "observations) who received 160 mg CAM2038 Q4W and became available",
      "after the initial model was built; the update re-estimated only the",
      "CAM2038 Q4W absorption parameters. Baseline characteristics in Table 1",
      "are reported separately for the initial-model dataset (N = 236) and the",
      "CAM2038 Q4W model-update dataset (N = 96, which overlaps the initial",
      "dataset); no combined N = 252 demographic column is published, so the",
      "ranges recorded here are the union of the two published columns.",
      "Concentrations below the 0.0250 ng/mL LLOQ were excluded from",
      "estimation but retained when producing visual predictive checks.",
      "NONMEM 7.3.0 / 7.5.0, FOCEI, log-transformed both sides. There were no",
      "missing covariate data."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Disposition: three-compartment, first-order elimination from central.
    # All values from Bjornsson 2023 Table 3 (final BPN population PK model,
    # 'Population mean / Estimate' column). Time unit hour, volumes L.
    # NONMEM V1 / V2 / V3 map to canonical vc / vp / vp2 and Q2 / Q3 to
    # canonical q / q2.
    # ---------------------------------------------------------------------
    lcl  <- log(52.1);   label("Plasma clearance CL (L/h)")                                    # Table 3: CL = 52.1 L/h (RSE 1.80%)
    lvc  <- log(64.3);   label("Central volume of distribution Vc (L), healthy reference")     # Table 3: Vc = 64.3 L (RSE 10.3%); healthy participants are the reference category
    lq   <- log(186);    label("Inter-compartmental clearance Q2 to peripheral1 (L/h)")        # Table 3: Q2 = 186 L/h (RSE 2.51%)
    lvp  <- log(130);    label("Volume of second compartment V2 (L)")                          # Table 3: V2 = 130 L (RSE 1.86%)
    lq2  <- log(60.3);   label("Inter-compartmental clearance Q3 to peripheral2 (L/h)")        # Table 3: Q3 = 60.3 L/h (RSE 4.13%)
    lvp2 <- log(1580);   label("Volume of third compartment V3 (L)")                           # Table 3: V3 = 1580 L (RSE 3.88%)

    # ---------------------------------------------------------------------
    # Sublingual BPN absorption: two parallel pathways.
    #   depot_sl1 -- sequential zero-order input over D_SL1 after a lag
    #                t_lag,SL1, then first-order ka_SL1 into central.
    #   depot_sl2 -- first-order ka_SL2 into central.
    # F_SL is the overall SL bioavailability (dose-dependent, see
    # e_dose_fdepot_sl); F_SL1 splits the bioavailable dose between the two
    # pathways. Both are constrained to (0, 1) by a logit transformation
    # (Bjornsson 2023 Eqs. 2-3), so the ini() values are logits.
    # ---------------------------------------------------------------------
    logitfdepot_sl  <- log(0.140 / (1 - 0.140)); label("Logit of sublingual BPN bioavailability F_SL at the 16 mg reference dose (fraction)") # Table 3: F_SL = 14.0% at 16 mg (RSE 3.54%); logit(0.140) = -1.8148
    logitfdepot_sl1 <- log(0.759 / (1 - 0.759)); label("Logit of the fraction of bioavailable SL dose entering depot_sl1 (fraction)")         # Table 3: F_SL1 = 75.9% (RSE 2.65%); logit(0.759) = +1.1472
    ltlag_sl1       <- log(0.171);               label("Lag time to depot_sl1 t_lag,SL1 (h)")                                                  # Table 3: t_lag,SL1 = 0.171 h (RSE 1.85%)
    ld1_sl1         <- log(0.419);               label("Zero-order input duration into depot_sl1 D_SL1 (h)")                                   # Table 3: D_SL1 = 0.419 h (RSE 4.40%)
    lka_sl1         <- log(1.72);                label("First-order absorption rate constant from depot_sl1 ka_SL1 (1/h)")                     # Table 3: ka,SL1 = 1.72 1/h (RSE 4.33%)
    lka_sl2         <- log(0.0875);              label("First-order absorption rate constant from depot_sl2 ka_SL2 (1/h; absorption half-life 7.9 h)") # Table 3: ka,SL2 = 0.0875 1/h (RSE 15.7%); Results 3.1.3 quotes half-life 7.9 h = log(2)/0.0875

    # ---------------------------------------------------------------------
    # CAM2038 Q1W absorption: two parallel pathways.
    #   depot_q1w1 -- sequential zero-order input over D_q1w1, then
    #                 first-order ka_q1w1 (fast pathway).
    #   depot_q1w2 -- first-order ka_q1w2 (slow, depot-matrix biodegradation).
    # Bioavailability was estimated close to the upper boundary of 1 and
    # therefore fixed to 1 (Results 3.1.1), so no F parameter is carried;
    # Fq1w1 splits the whole dose between the two pathways.
    # ---------------------------------------------------------------------
    logitfdepot_q1w1 <- log(0.455 / (1 - 0.455)); label("Logit of the fraction of CAM2038 Q1W dose entering the fast pathway depot_q1w1 (fraction), male healthy reference") # Table 3: Fq1w1 = 45.5% (RSE 8.07%); logit(0.455) = -0.1806; reference = male healthy, non-thigh injection site
    ld1_q1w1         <- log(10.1);                label("Zero-order input duration into depot_q1w1 D_q1w1 (h)")                                                              # Table 3: D_q1w1 = 10.1 h (RSE 4.68%)
    lka_q1w1         <- log(0.0401);              label("First-order absorption rate constant from depot_q1w1 ka_q1w1 (1/h)")                                                # Table 3: ka,q1w1 = 0.0401 1/h (RSE 4.67%)
    lka_q1w2         <- log(0.00565);             label("First-order absorption rate constant from depot_q1w2 ka_q1w2 (1/h; absorption half-life 123 h)")                    # Table 3: ka,q1w2 = 0.00565 1/h (RSE 6.77%); Results 3.1.3 quotes half-life 123 h = log(2)/0.00565

    # ---------------------------------------------------------------------
    # CAM2038 Q4W absorption: two parallel first-order pathways (no
    # zero-order input component). Bioavailability fixed to 1 as for Q1W.
    # These three parameters were re-estimated in the model update that
    # added the 160 mg Q4W dose group (Results 3.1.1).
    # ---------------------------------------------------------------------
    logitfdepot_q4w1 <- log(0.0863 / (1 - 0.0863)); label("Logit of the fraction of CAM2038 Q4W dose entering the fast pathway depot_q4w1 (fraction)") # Table 3: Fq4w1 = 8.63% (RSE 10.6%); logit(0.0863) = -2.3595
    lka_q4w1         <- log(0.0441);                label("First-order absorption rate constant from depot_q4w1 ka_q4w1 (1/h)")                        # Table 3: ka,q4w1 = 0.0441 1/h (RSE 2.84%)
    lka_q4w2         <- log(0.00166);               label("First-order absorption rate constant from depot_q4w2 ka_q4w2 (1/h; absorption half-life 418 h)") # Table 3: ka,q4w2 = 0.00166 1/h (RSE 6.16%); Results 3.1.3 quotes half-life 418 h = log(2)/0.00166

    # ---------------------------------------------------------------------
    # Covariate effects (Table 3, lower block).
    # Continuous covariates on untransformed parameters enter as power
    # models; categorical covariates enter as a fractional difference to the
    # most common category for untransformed parameters, and additively on
    # the logit scale for logit-transformed parameters (Bjornsson 2023
    # Eqs. 4-5). The sublingual dose effect is a power model on F_SL as
    # printed verbatim in Results 3.1.3.
    # ---------------------------------------------------------------------
    e_age_cl              <- -0.233; label("Power exponent for age on CL (unitless; reference 35 years)")                                  # Table 3: Age covariate on CL = -0.233 (RSE 23.4%)
    e_wt_cl               <-  0.413; label("Power exponent for body weight on CL (unitless; reference 72.4 kg)")                           # Table 3: WT covariate on CL = 0.413 (RSE 20.9%)
    e_dis_oud_vc          <-  2.69;  label("Fractional difference in Vc for participants with OUD vs healthy (unitless)")                   # Table 3: Population covariate on Vc = 2.69 (RSE 20.7%); 64.3 * (1 + 2.69) = 237 L per Results 3.1.3
    e_dis_oud_fdepot_q1w1 <- -0.671; label("Additive shift on logit(Fq1w1) for participants with OUD vs healthy (unitless)")                # Table 3: Population covariate on Fq1w1 = -0.671 (RSE 30.3%); footnote a 'Patients with OUD versus healthy volunteers'
    e_sexf_fdepot_q1w1    <-  0.576; label("Additive shift on logit(Fq1w1) for female vs male sex (unitless)")                              # Table 3: Sex covariate on Fq1w1 = 0.576 (RSE 36.1%); footnote b 'Females versus males'
    e_dose_fdepot_sl      <- -0.371; label("Power exponent for sublingual dose on F_SL (unitless; reference 16 mg)")                        # Table 3: Dose dependency on F_SL = -0.371 (RSE 18.4%); Results 3.1.3 F_SL = 0.14 * (Dose/16)^-0.371

    # ---------------------------------------------------------------------
    # Inter-individual variability (Table 3, 'IIV / CV (%)' column).
    # The paper reports one CV (%) column spanning both exponentially
    # distributed parameters (Eq. 1) and logit-transformed parameters
    # (Eqs. 2-3). Because a coefficient of variation is undefined for a
    # logit-normal parameter, the column is read here as 100 * omega, i.e.
    # the standard deviation of eta on the transformed scale, so
    # variance = (CV / 100)^2. See the vignette 'Assumptions and deviations'
    # section: this reading reproduces the geometric CV percentages of
    # Table 4 across all 11 simulated regimens, and the alternative
    # log-normal reading omega^2 = log(1 + CV^2) differs by at most a few
    # percent for the exponential parameters.
    # No IIV is reported for Q2, V2, t_lag,SL1, ka,SL1, D_SL1, ka,q1w1,
    # D_q1w1, or ka,q4w1.
    # ---------------------------------------------------------------------
    etalcl                ~ 0.043681  # Table 3: IIV(CL)     = 20.9% (RSE 5.83%); 0.209^2
    etalvc                ~ 0.617796  # Table 3: IIV(Vc)     = 78.6% (RSE 8.21%); 0.786^2
    etalq2                ~ 0.144400  # Table 3: IIV(Q3)     = 38.0% (RSE 9.14%); 0.380^2
    etalvp2               ~ 0.088804  # Table 3: IIV(V3)     = 29.8% (RSE 6.74%); 0.298^2
    etalogitfdepot_sl     ~ 0.096721  # Table 3: IIV(F_SL)   = 31.1% (RSE 9.70%); 0.311^2
    etalka_sl2            ~ 0.744769  # Table 3: IIV(ka,SL2) = 86.3% (RSE 9.82%); 0.863^2
    etalogitfdepot_sl1    ~ 0.580644  # Table 3: IIV(F_SL1)  = 76.2% (RSE 13.5%); 0.762^2
    etalka_q1w2           ~ 0.324900  # Table 3: IIV(ka,q1w2)= 57.0% (RSE 7.49%); 0.570^2
    etalogitfdepot_q1w1   ~ 0.811801  # Table 3: IIV(Fq1w1)  = 90.1% (RSE 8.12%); 0.901^2
    etalka_q4w2           ~ 0.150544  # Table 3: IIV(ka,q4w2)= 38.8% (RSE 8.45%); 0.388^2
    etalogitfdepot_q4w1   ~ 0.848241  # Table 3: IIV(Fq4w1)  = 92.1% (RSE 8.61%); 0.921^2

    # ---------------------------------------------------------------------
    # Residual unexplained variability. Results 3.1.1: 'The residual
    # unexplained variability (RUV) for BPN was best described by an
    # additive error model on the log-transformed data (approximately
    # proportional on untransformed data)', and Table 3 reports it directly
    # as sigma_prop = 27.7%.
    # ---------------------------------------------------------------------
    propSd <- 0.277; label("Proportional residual error on BPN plasma concentration (fraction)") # Table 3: sigma_prop = 27.7% (RSE 30.3%)
  })

  model({
    # ------------------ Reference values for covariate models -------------
    # Age 35 years and body weight 72.4 kg are the values of the typical
    # individual used for the paper's simulations (Methods 2.4); the
    # sublingual reference dose is 16 mg (Results 3.1.3).
    age_ref     <- 35
    wt_ref      <- 72.4
    dose_sl_ref <- 16

    # ------------------ Individual disposition parameters -----------------
    # CL = 52.1 * (AGE/35)^-0.233 * (WT/72.4)^0.413 (Results 3.1.3).
    cl  <- exp(lcl + etalcl) * (AGE / age_ref)^e_age_cl * (WT / wt_ref)^e_wt_cl

    # Vc carries the OUD fractional difference (Eq. 5): healthy 64.3 L,
    # OUD 64.3 * (1 + 2.69) = 237 L.
    vc  <- exp(lvc + etalvc) * (1 + e_dis_oud_vc * DIS_OUD)

    q   <- exp(lq)
    vp  <- exp(lvp)
    q2  <- exp(lq2  + etalq2)
    vp2 <- exp(lvp2 + etalvp2)

    # ------------------ Sublingual absorption parameters ------------------
    # Typical F_SL is dose-dependent as printed in Results 3.1.3:
    #   F_SL = 0.14 * (Dose / 16 mg)^-0.371
    # (18.1 / 14.0 / 12.0 percent at 8 / 16 / 24 mg). IIV then acts on the
    # logit of that dose-adjusted typical value, per Eqs. 2-3, which keeps
    # every individual F_SL inside (0, 1).
    # (etalogitfdepot_sl cannot be mu-referenced because the dose-dependent
    # typical value is not a bare theta; rxode2 emits a non-mu-reference note
    # for it, which affects estimation efficiency only, not simulation.)
    fdepot_sl_tv <- (1 / (1 + exp(-logitfdepot_sl))) *
      (DOSE_BPN_SL_MG / dose_sl_ref)^e_dose_fdepot_sl
    lgt_sl       <- log(fdepot_sl_tv / (1 - fdepot_sl_tv)) + etalogitfdepot_sl
    fdepot_sl    <- 1 / (1 + exp(-lgt_sl))

    lgt_sl1    <- logitfdepot_sl1 + etalogitfdepot_sl1
    fdepot_sl1 <- 1 / (1 + exp(-lgt_sl1))

    tlag_sl1 <- exp(ltlag_sl1)
    d1_sl1   <- exp(ld1_sl1)
    ka_sl1   <- exp(lka_sl1)
    ka_sl2   <- exp(lka_sl2 + etalka_sl2)

    # ------------------ CAM2038 Q1W absorption parameters -----------------
    # Categorical covariates on the logit-transformed Fq1w1 enter additively
    # on the logit scale: OUD lowers Fq1w1 (45.5 -> 29.9 percent) and female
    # sex raises it (45.5 -> 59.8 percent), matching the directions stated
    # in Results 3.1.3. A thigh injection routes the whole dose through the
    # slow pathway (Fq1w1 fixed to 0; Results 3.1.1).
    lgt_q1w1    <- logitfdepot_q1w1 +
      e_dis_oud_fdepot_q1w1 * DIS_OUD +
      e_sexf_fdepot_q1w1 * SEXF +
      etalogitfdepot_q1w1
    fdepot_q1w1 <- (1 - INJSITE_THIGH) / (1 + exp(-lgt_q1w1))

    d1_q1w1 <- exp(ld1_q1w1)
    ka_q1w1 <- exp(lka_q1w1)
    ka_q1w2 <- exp(lka_q1w2 + etalka_q1w2)

    # ------------------ CAM2038 Q4W absorption parameters -----------------
    lgt_q4w1    <- logitfdepot_q4w1 + etalogitfdepot_q4w1
    fdepot_q4w1 <- 1 / (1 + exp(-lgt_q4w1))
    ka_q4w1     <- exp(lka_q4w1)
    ka_q4w2     <- exp(lka_q4w2 + etalka_q4w2)

    # ------------------ ODE system ----------------------------------------
    # Compartment slots follow the d/dt() declaration order below. Dose IV
    # into 'central'; dose the SL and CAM2038 routes into BOTH depots of the
    # relevant pair with the same amount -- the f() statements below split
    # the administered dose between the parallel pathways, exactly as the
    # source NONMEM model does with per-compartment F.
    d/dt(depot_sl1)  <- -ka_sl1  * depot_sl1
    d/dt(depot_sl2)  <- -ka_sl2  * depot_sl2
    d/dt(depot_q1w1) <- -ka_q1w1 * depot_q1w1
    d/dt(depot_q1w2) <- -ka_q1w2 * depot_q1w2
    d/dt(depot_q4w1) <- -ka_q4w1 * depot_q4w1
    d/dt(depot_q4w2) <- -ka_q4w2 * depot_q4w2
    d/dt(central)    <-  ka_sl1  * depot_sl1  + ka_sl2  * depot_sl2 +
                         ka_q1w1 * depot_q1w1 + ka_q1w2 * depot_q1w2 +
                         ka_q4w1 * depot_q4w1 + ka_q4w2 * depot_q4w2 -
                         cl * central / vc -
                         q  * central / vc + q  * peripheral1 / vp -
                         q2 * central / vc + q2 * peripheral2 / vp2
    d/dt(peripheral1) <- q  * central / vc - q  * peripheral1 / vp
    d/dt(peripheral2) <- q2 * central / vc - q2 * peripheral2 / vp2

    # ------------------ Bioavailability, duration, lag time ---------------
    # SL: the bioavailable fraction F_SL of the dose is split F_SL1 /
    # (1 - F_SL1) between the fast and slow pathways.
    f(depot_sl1)  <- fdepot_sl * fdepot_sl1
    f(depot_sl2)  <- fdepot_sl * (1 - fdepot_sl1)
    # CAM2038 Q1W and Q4W: bioavailability fixed to 1, so the pathway
    # fractions sum to the whole administered dose.
    f(depot_q1w1) <- fdepot_q1w1
    f(depot_q1w2) <- 1 - fdepot_q1w1
    f(depot_q4w1) <- fdepot_q4w1
    f(depot_q4w2) <- 1 - fdepot_q4w1

    # Sequential zero-order then first-order pathways. Dose records for
    # depot_sl1 and depot_q1w1 must carry rate = -2 so rxode2 applies the
    # modelled duration.
    dur(depot_sl1)  <- d1_sl1
    dur(depot_q1w1) <- d1_q1w1
    alag(depot_sl1) <- tlag_sl1

    # ------------------ Observation and residual error --------------------
    # Dose in mg and vc in L give mg/L = ug/mL; multiply by 1000 for the
    # ng/mL units of the published concentrations (Methods 2.1: calibration
    # range 0.0250-10.0 ng/mL).
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
