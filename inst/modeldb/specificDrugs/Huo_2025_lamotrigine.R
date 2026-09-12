Huo_2025_lamotrigine <- function() {
  description <- "One-compartment population PK model with first-order absorption and elimination for lamotrigine (LTG) in 128 Chinese peripregnancy women with epilepsy on lamotrigine monotherapy (Huo 2025 Eqs 1-7, Table 4 'Final Model' column). Ka (1.93 1/h) and apparent volume V/F (68.8 L) were both FIXED from the literature because the therapeutic-drug-monitoring data were almost all steady-state troughs and carried no absorption or distribution information; apparent clearance CL/F was the only structural parameter estimated. CL/F = 2.42 L/h at 59.8 kg and carries an estimated body-weight power exponent of 0.95, an exponential five-level peripregnancy-stage effect using the paper's own Classification C staging (gestational-week nodes at 5, 14 and 28 weeks plus a postpartum level), and an exponential valproate-comedication effect that lowers CL/F by 45% (exp(-0.60)). Residual error is combined proportional plus additive. Fit in Phoenix NLME 8.3 by FOCE-ELS."
  reference   <- "Huo J, Liu Y, Yang J, Chen M, Yang L, Wang L, Zhang D, Liu T, Gao W, Dai H, Mei S, Zhao Z. Dosing Optimization of Lamotrigine in Peripregnancy Epilepsy Through PopPK Modelling and Simulation. Drug Des Devel Ther. 2025;19:10243-10254. doi:10.2147/DDDT.S541597. PMCID PMC12645405. Structural equations from Eqs 1-5 (p 10246); covariate model from Eq 6 (p 10248) and the peripregnancy-stage / inhibitor coefficient block (p 10249); V/F from Eq 7 (p 10249); parameter estimates from Table 4 'Final Model'; peripregnancy staging from Table 1 row C; cohort demographics from Table 2."
  vignette    <- "Huo_2025_lamotrigine"
  units       <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot   = list(analyte = "lamotrigine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "lamotrigine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters CL/F only, as the power term (WT/59.8)^0.95 printed in Huo 2025 Eq 6. The normalisation constant 59.8 kg is NOT the cohort mean body weight, which Table 2 gives as 60.75 +/- 13.47 kg (range 40.00-114.00); the paper never says what 59.8 kg is, and it is most plausibly the cohort MEDIAN weight. It is transcribed exactly as printed in Eq 6 and never adjusted. Note that the paper's own Discussion sentence 'this study set the Vd value as 59.8 L' re-uses this same number as a volume - that sentence is a transcription slip and 59.8 is the weight normalisation constant, not V/F (see the vignette Errata and the lvc source-trace comment). Body weight is time-varying across a pregnancy; the paper does not state whether the per-record weight or a single baseline weight was used, but the strong stage-3/stage-4 CL/F elevation that the Discussion reports (193% / 199% of stage 1) is only reproduced when the pregnancy weight gain is carried in WT alongside the stage effect, so a per-record (time-varying) weight is the reading used here.",
      source_name        = "BW"
    ),
    EGA = list(
      description        = "Maternal estimated gestational age at the observation",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Carries the peripregnancy stage during pregnancy. The paper does not fit gestational age as a continuous covariate; it fits a five-level categorical 'peripregnancy stage' whose levels are gestational-week bands. Huo 2025 Table 1 row C ('Classification C was selected for the model') defines them as stage 1 = GA < 5 weeks, stage 2 = 5 <= GA < 14 weeks, stage 3 = 14 <= GA <= 28 weeks, stage 4 = 28 < GA < delivery, stage 5 = postpartum. Those cutoffs are the paper's own model equation, so the banding is applied inside model() from EGA rather than requiring the user to supply a pre-computed stage column. Per the EGA register entry, gestational-week stratification is carried on EGA rather than by introducing a trimester-indicator canonical. Set EGA to the gestational age in weeks for pregnancy records; for postpartum records EGA is ignored (TPP > 0 selects stage 5). Observation counts per stage: 13 / 63 / 111 / 84 / 22 (Table 2).",
      source_name        = "Peripregnancy stage"
    ),
    TPP = list(
      description        = "Time postpartum (time since delivery)",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used only as the postpartum GATE: TPP = 0 for every pregnancy record and TPP > 0 for every postpartum record, which selects peripregnancy stage 5 (Huo 2025 Table 1 row C). The paper models postpartum as a single step level, NOT as a continuous time-since-delivery decay, so no magnitude of TPP other than 'zero vs positive' affects the prediction. The postpartum observations span 1 to 84 days after delivery (Limitations), which the authors flag as a limitation because LTG CL/F is reported elsewhere to return to prepregnancy levels within 2-4 weeks; a single step level therefore averages across a window in which recovery is still in progress.",
      source_name        = "Postpartum"
    ),
    CONMED_VPA = list(
      description        = "Concomitant valproate (valproic acid) therapy",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = no concomitant enzyme inhibitor (lamotrigine without valproate)",
      notes              = "Huo 2025 calls this covariate 'inhibitor' / 'enzyme inhibitors', but Data Collection defines the enzyme-inhibitor class as containing exactly one drug: 'enzyme inhibitors (valproic acid(VPA))'. Valproate was used by 3.91% (5) of the 128 patients (Table 2), and both the Discussion ('co-administration of VPA could decrease LTG CL/F by 46%') and the Table 6 dosing recommendations name VPA explicitly, so the drug-specific canonical CONMED_VPA is used rather than the pooled CONMED_UGT_INH. Enters CL/F as exp(-0.60 * CONMED_VPA) = 0.549, a 45.1% reduction, which reproduces the paper's stated 46%. Time-varying in principle; the cohort is on chronic maintenance therapy.",
      source_name        = "inhibitor"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened and not retained. Huo 2025 'Other Unexplored Covariates in LTG Metabolism' attributes the absence of an age effect to the narrow age span of the cohort (19-36 years, Table 2)."
    ),
    BSA = list(
      description = "Body surface area",
      units       = "m^2",
      type        = "continuous",
      notes       = "Screened and not retained. Highly correlated with body weight; the Covariate Model section excluded covariate pairs with a correlation coefficient > 0.5 from entering the model simultaneously (Figure S1)."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened and not retained; correlated with body weight, see the BSA note."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened and not retained. Huo 2025 reports that liver and kidney function were within normal limits for the majority of participants, so no hepatic or renal marker reached significance."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened and not retained; see the ALB note."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened and not retained; see the ALB note."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened and not retained; see the ALB note."
    ),
    DBIL = list(
      description = "Direct (conjugated) bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened and not retained; see the ALB note."
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Screened and not retained. Only one patient had a creatinine below the reference range on at least two occasions, so renal function was effectively invariant in the cohort."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened and not retained; see the BUN note."
    ),
    CONMED_CBZ = list(
      description = "Concomitant carbamazepine",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as an enzyme inducer and not retained. Huo 2025 attributes this to sample size: 'only 3 patients received CBZ in this study', versus other studies that did retain a CBZ effect."
    ),
    CONMED_OXC = list(
      description = "Concomitant oxcarbazepine",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as an enzyme inducer and not retained (7 patients). Huo 2025 argues mechanistically that oxcarbazepine's active monohydroxy metabolite is only a weak UGT inducer and so is unlikely to affect a drug cleared by glucuronidation."
    ),
    CONMED_PB = list(
      description = "Concomitant phenobarbital",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as an enzyme inducer and not retained (5 patients)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 128,
    n_studies      = 1,
    n_observations = 293,
    age_range      = "19-36 years; mean 28.24 +/- 3.79 years (Table 2).",
    weight_range   = "40.00-114.00 kg; mean 60.75 +/- 13.47 kg (Table 2).",
    sex_female_pct = 100,
    race_ethnicity = c(Asian = 100),
    disease_state  = "Women with epilepsy (2014 ILAE classification, confirmed by two board-certified neurologists) diagnosed during pregnancy or within 6 months postpartum, receiving lamotrigine monotherapy. Age < 18 years at conception, major psychiatric comorbidity, poor adherence and > 20% missing data were exclusions (Methods, Study Population).",
    dose_range     = "25-800 mg/day oral tablets, typically given as one or two doses per day (Methods, Dosing Regimens and Concentration Measurement).",
    regions        = "China (two centres: Beijing Tiantan Hospital, Capital Medical University; The Second Affiliated Hospital, Zhejiang University School of Medicine), January 2015 - May 2024.",
    co_medication  = "42.19% of patients used at least one concomitant antiseizure medication. Individual comedication rates (Table 2): levetiracetam 17.19% (22), perampanel 4.69% (6), oxcarbazepine 5.47% (7), phenobarbital 3.91% (5), valproic acid 3.91% (5), carbamazepine 2.34% (3), clobazam 2.34% (3), topiramate 0.78% (1), lacosamide 0.78% (1), vigabatrin 0.78% (1). Only valproate (the enzyme-inhibitor class) was retained as a covariate.",
    gestational_stage_distribution = "Observations per peripregnancy stage under Classification C: stage 1 (GA < 5 weeks) 13, stage 2 (5-14 weeks) 63, stage 3 (14-28 weeks) 111, stage 4 (28 weeks to delivery) 84, stage 5 (postpartum) 22 (Table 2). Postpartum samples span 1-84 days after delivery (Limitations).",
    notes          = "Multicentre retrospective therapeutic-drug-monitoring cohort; the majority of the 293 samples are steady-state trough concentrations drawn after at least 7 days of continuous lamotrigine therapy. Concentrations were measured by a validated UPLC-MS/MS assay linear from 1.37 to 20.9 mg/L with an LLOQ of 1.37 mg/L. Baseline demographics are Table 2. No pre-pregnancy therapeutic-drug-monitoring data were available, so the fold-changes reported by the paper are all relative to peripregnancy stage 1 rather than to a true prepregnancy baseline (Limitations). No genotype data were collected."
  )

  ini({
    # Structural parameters. Final-model estimates from Huo 2025 Table 4
    # 'Final Model' column, with the covariate model as printed in Eq 6.

    # Ka FIXED from the literature: "The absence of both the absorption and
    # distribution phases in the sampled data rendered the estimates of the
    # absorption rate constant (Ka) and the apparent distribution volume (V/F)
    # unreliable. Therefore, Ka and V/F were fixed at 1.93 h-1 and 68.8 L,
    # respectively. Only CL/F was estimated during model development." (Base
    # Model, p 10246). The value is a meta-analytic historical Ka taken from
    # Wang et al (reference 19), per the Discussion. Table 4 prints Ka in all
    # three columns with no %RSE and no 95% CI, confirming it was not estimated.
    lka <- fixed(log(1.93)); label("Absorption rate constant (Ka, 1/h)")                      # Base Model p 10246; Table 4 Ka row (1.93, no RSE / CI in any column)

    # V/F FIXED at 68.8 L. Printed three times: the Base Model sentence above,
    # Eq 7 ("Vd(L) = 68.8 fixed", p 10249), and Table 4's Vd row, which reads
    # "68.8 (fixed)" in the Base Model, Final Model AND Bootstrap columns.
    # The Discussion sentence "this study set the Vd value as 59.8 L" conflicts
    # with all three; 59.8 is the body-weight normalisation constant of Eq 6, so
    # that sentence is a transcription slip and 68.8 L is used. See the vignette
    # Errata.
    lvc <- fixed(log(68.8)); label("Apparent central volume of distribution (V/F, L)")        # Eq 7 p 10249; Base Model p 10246; Table 4 Vd row "68.8 (fixed)"

    lcl <- log(2.42);        label("Apparent oral clearance (CL/F, L/h) at 59.8 kg, peripregnancy stage 1, no valproate")  # Table 4 Final Model CL row: 2.42 L/h (%RSE 12.17, 95% CI 1.84-3.00; bootstrap median 2.44)

    # Body weight on CL/F as a power term referenced to 59.8 kg (Eq 6). The
    # exponent was ESTIMATED, not fixed at an allometric 0.75: Table 4 prints
    # both a %RSE and a 95% CI for it, and adding BW to the model dropped the
    # OFV by 19.15 (Table 3, model 3).
    e_wt_cl <- 0.95;         label("Power exponent: WT on CL/F, referenced to 59.8 kg (unitless)")  # Table 4 Final Model 'BW on CL' row: 0.95 (%RSE 22.98, 95% CI 0.52-1.38); exponent and 59.8 kg reference both from Eq 6

    # Peripregnancy-stage effects on CL/F. Eq 6 applies the stage term as
    # exp(peripregnancy stage), and the coefficient block on p 10249 gives the
    # per-level value: stage 1 = 0 (the reference level), stage 2 = 0.28,
    # stage 3 = 0.59, stage 4 = 0.57, stage 5 (postpartum) = -0.33. The level
    # boundaries are Classification C of Table 1. Stage 1 carries no parameter
    # because its coefficient is exactly 0 by construction.
    e_pregstage2_cl  <-  0.28; label("Effect of peripregnancy stage 2 (5-14 weeks GA) on CL/F (log-scale shift)")        # p 10249 coefficient block; Table 4 'Peripregnancy stage2 on CL' 0.28 (%RSE 47.01, 95% CI 0.02-0.53)
    e_pregstage3_cl  <-  0.59; label("Effect of peripregnancy stage 3 (14-28 weeks GA) on CL/F (log-scale shift)")       # p 10249 coefficient block; Table 4 'Peripregnancy stage3 on CL' 0.59 (%RSE 21.83, 95% CI 0.34-0.85)
    e_pregstage4_cl  <-  0.57; label("Effect of peripregnancy stage 4 (28 weeks GA to delivery) on CL/F (log-scale shift)")  # p 10249 coefficient block; Table 4 'Peripregnancy stage4 on CL' 0.57 (%RSE 24.17, 95% CI 0.30-0.84)
    e_postpartum_cl  <- -0.33; label("Effect of peripregnancy stage 5 (postpartum) on CL/F (log-scale shift)")           # p 10249 coefficient block; Table 4 'Peripregnancy stage5 on CL' -0.33 (%RSE -42.98, 95% CI -0.62 to -0.05)

    # Valproate comedication on CL/F, also exponential per Eq 6
    # (exp(inhibitor)). exp(-0.60) = 0.549, i.e. a 45.1% reduction in CL/F,
    # which reproduces the Discussion's "co-administration of VPA could
    # decrease LTG CL/F by 46%".
    e_conmed_vpa_cl <- -0.60; label("Effect of concomitant valproate on CL/F (log-scale shift)")   # p 10249 "inhibitor = -0.60, when comedicated with enzyme inhibitors, otherwise = 0"; Table 4 'Inhibitors on CL' -0.60 (%RSE -30.54, 95% CI -0.95 to -0.24)

    # IIV on CL/F. Eq 4 is the exponential ("index") random-effect model
    # theta_i = theta_TV * exp(eta_i), so eta is log-normal and nlmixr2's ini()
    # takes its VARIANCE. Table 4 reports IIV_CL as a CV% (32.96% in the final
    # model, down from 50.76% in the base model), so
    #   omega^2 = log(CV^2 + 1) = log(0.3296^2 + 1) = 0.10313.
    # No other parameter carries IIV: Ka and V/F were fixed, which the paper
    # explicitly notes "restricted the assessment of interindividual
    # variability in these parameters".
    etalcl ~ 0.10313                                                                          # Table 4 Final Model 'IIV CL (CV%)' = 32.96%; omega^2 = log(0.3296^2 + 1) via Eq 4

    # Residual error. Eq 5 is Cobs = Cpred * (1 + eps1) + eps2, i.e. combined
    # proportional plus additive in linear concentration space, with sigma1 the
    # proportional SD (reported as a CV) and sigma2 the additive SD.
    # Table 4 labels the sigma2 row "additive, mmol/L"; lamotrigine
    # concentrations are in mg/L throughout the paper (LLOQ 1.37 mg/L,
    # calibration range 1.37-20.9 mg/L, therapeutic range 2.5-15 mg/L), so
    # "mmol/L" is a units slip and the additive SD is carried in mg/L. Its
    # magnitude (0.004 mg/L, versus 0.38 in the base model) is negligible
    # against an LLOQ of 1.37 mg/L. See the vignette Errata.
    propSd <- 0.31;  label("Proportional residual error (fraction)")                           # Table 4 Final Model 'sigma1 (multiplicative, CV)' = 0.31 (95% CI 0.27-0.34)
    addSd  <- 0.004; label("Additive residual error (mg/L)")                                   # Table 4 Final Model 'sigma2 (additive)' = 0.004 (95% CI 0.003-0.005); row labelled mmol/L, read as mg/L
  })

  model({
    # Peripregnancy stage under Huo 2025 Classification C (Table 1 row C,
    # "Classification C was selected for the model"), derived from the maternal
    # gestational age and the time since delivery:
    #   stage 1  TPP == 0 and EGA < 5 weeks         reference, coefficient 0
    #   stage 2  TPP == 0 and 5 <= EGA < 14 weeks
    #   stage 3  TPP == 0 and 14 <= EGA <= 28 weeks
    #   stage 4  TPP == 0 and EGA > 28 weeks (to delivery)
    #   stage 5  TPP > 0                            postpartum
    # The postpartum gate is tested first because Classification C makes
    # postpartum a level of the same categorical, not an extension of the
    # gestational-week axis.
    if (TPP > 0) {
      stage_cl <- e_postpartum_cl
    } else if (EGA < 5) {
      stage_cl <- 0
    } else if (EGA < 14) {
      stage_cl <- e_pregstage2_cl
    } else if (EGA <= 28) {
      stage_cl <- e_pregstage3_cl
    } else {
      stage_cl <- e_pregstage4_cl
    }

    # Ka and V/F were fixed, so neither carries an eta (Base Model, p 10246).
    ka <- exp(lka)
    vc <- exp(lvc)

    # Huo 2025 Eq 6:
    #   CL(L/h) = 2.42 * (BW/59.8)^0.95 * exp(peripregnancy stage)
    #                  * exp(inhibitor) * exp(eta_CL)
    cl <- exp(lcl + etalcl) *
      (WT / 59.8)^e_wt_cl *
      exp(stage_cl) *
      exp(e_conmed_vpa_cl * CONMED_VPA)

    kel <- cl / vc

    # One compartment with first-order absorption and first-order elimination,
    # Huo 2025 Eqs 1-3. Eq 2 is written dAc/dt = Ka*Aa - CL*Cc, which is the
    # same as -kel*Ac once Cc = Ac/Vd (Eq 3) is substituted.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Dose in mg and V/F in L give mg/L, the unit of the paper's 2.5-15 mg/L
    # therapeutic reference range (Model-Informed LTG Dosing Regimens).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
