Byrne_2022_imr687 <- function() {
  description <- "One-compartment population PK model with first-order absorption for IMR-687 (a selective PDE9 inhibitor) in healthy subjects and patients with sickle cell disease (SCD), coupled with a repeated time-to-event (RTTE) exposure-response model for vaso-occlusive crisis (VOC) events. The PD hazard uses a saturable (Michaelis-Menten) drug effect on a constant baseline hazard. The model supports forward simulation of typical-value PK and cumulative-VOC hazard at any once-daily dose; the published covariate effects (body weight on CL/F and V/F; capsule formulation, capsule daily dose, and high-fat meal on absorption) are NOT encoded because the source conference poster reports the covariate point estimates without the functional forms or reference values needed to apply them."
  reference <- paste(
    "Byrne R, Ballal R, Ruiz-Garcia A.",
    "A Population Pharmacokinetic Model and Exposure-Response Model of",
    "Repeated Time Event (RTTE) to Justify a Dose Increase in Patients",
    "with Sickle Cell Disease.",
    "American Conference on Pharmacometrics (ACoP) 2022;",
    "Metrum Research Group / Imara Inc.",
    "doi:10.70534/zdmj9414.",
    sep = " "
  )
  vignette <- "Byrne_2022_imr687"
  units <- list(
    time          = "hour",
    dosing        = "mg",
    concentration = "ng/mL"
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "The poster Results-PK-ER text identifies 'body weight on clearance (CL/F) and volume (V/F)' as covariates retained in the final population PK model, but reports no parameter estimate, no functional form (allometric vs linear), and no reference weight. The model file therefore documents WT in covariatesDataExcluded and does NOT scale CL/F or V/F by body weight. A future re-extraction with the underlying NONMEM control stream should promote WT into covariateData and add e_wt_cl / e_wt_vc effects."
    ),
    FORM_CAPSULE = list(
      description = "Capsule formulation indicator: 1 = capsule, 0 = tablet (the per-paper comparator).",
      units       = "(binary)",
      type        = "binary",
      notes       = "The poster Results-PK-ER text reports a 'formulation effect: 0.239 (0.143, 0.401)' on absorption (Ka), with capsule vs tablet as the contrast. The functional form (multiplicative on Ka? log-additive? change in Ka magnitude?) and the reference category (capsule = 0 or tablet = 0) are not stated. The model file therefore documents FORM_CAPSULE in covariatesDataExcluded and does NOT encode it as a Ka effect. Reference category and functional form must be confirmed from the underlying NONMEM control stream before promotion."
    ),
    DOSE = list(
      description = "Capsule daily dose (mg) entered as a covariate effect on absorption (Ka) in the source PPK model. Distinct from the dose-event AMT carried in the dosing records.",
      units       = "mg",
      type        = "continuous",
      notes       = "The poster Results-PK-ER text reports a 'dose effect: -0.602 (-0.908, -0.296)' on absorption (Ka) for the capsule formulation. The functional form, reference dose, and whether the covariate is per-dose or per-subject are not stated. The model file therefore documents DOSE in covariatesDataExcluded and does NOT encode the dose-on-Ka effect; the dose amount still enters the simulation through the standard rxode2 dosing-event AMT. Promotion requires the underlying control stream."
    ),
    FED_HIGHFAT = list(
      description = "High-fat-meal indicator at time of dosing: 1 = high-fat meal, 0 = fasted or low-fat reference.",
      units       = "(binary)",
      type        = "binary",
      notes       = "The poster Results-PK-ER text reports a 'high fat meal effect: 0.176, (0.121, .254)' on absorption (Ka) for the capsule formulation only. The functional form and reference category are not stated. The model file therefore documents FED_HIGHFAT in covariatesDataExcluded and does NOT encode it as a Ka effect. Promotion requires the underlying control stream."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 112L,
    n_studies      = 2L,
    age_range      = NA_character_,
    age_median     = NA_character_,
    weight_range   = NA_character_,
    weight_median  = NA_character_,
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state  = "Healthy adult subjects (Phase 1a, IMR-SCD-101 or related FIH study) and adult patients with sickle cell disease (SCD; Phase 2a study IMR-SCD-102). The Phase 2a SCD cohort enrolled 92 subjects across placebo, 50 mg, 100 mg, 50-100 mg up-titrated, and 100-200 mg up-titrated arms (Figure 2 dose-cohort tabulation). VOC events were the efficacy endpoint.",
    dose_range     = "Phase 2a: oral 50, 100, or 200 mg IMR-687 once daily for 24 weeks; titration arms started at 50 mg or 100 mg and up-titrated to 100 mg or 200 mg respectively. Simulated higher doses considered in the poster: 300, 400, and 600 mg QD.",
    regions        = NA_character_,
    notes          = "The PPK dataset pooled Phase 1a (healthy subjects) and Phase 2a (SCD patients) for a total of 112 subjects. The PD dataset for the RTTE exposure-response analysis comprised 92 Phase 2a subjects (placebo n=30, 50 mg n=15, 100 mg n=12, 50-100 mg n=21, 100-200 mg n=14; per Figure 2). The poster reports baseline demographic distributions only graphically (Figure 2 dot-plot tabulation), so per-field demographic ranges and medians are not reproduced here."
  )

  ini({
    # ----- Structural PK parameters -----
    # Source: poster Results - PK-ER Model paragraph (page 1). Point estimates
    # are typical-value population means under the final PPK model. The poster
    # also lists body-weight effects on CL/F and V/F as part of the final model
    # but does not report parameter values; those effects are not encoded.
    lcl <- log(14.6); label("Apparent oral clearance CL/F (L/h)")  # poster Results-PK-ER text: "CL/F: 14.6 (13.8, 15.4 95% CI) L/h"
    lvc <- log(104) ; label("Apparent central volume of distribution V/F (L)")  # poster Results-PK-ER text: "V/F: 104 (98.8, 109) L"
    lka <- log(4.90); label("First-order absorption rate constant Ka (1/h)")    # poster Results-PK-ER text: "Ka: 4.90 (3.12, 7.70) 1/h"

    # ----- RTTE base hazard and EC50 -----
    # Source: poster Results-PK-ER Model paragraph (page 1, RTTE block).
    # Hazard time units are taken to be 1/day to match the VPC figure's
    # x-axis ("Time (days)"); inside model() the hourly-equivalent rate is
    # computed as lambda / 24 so the constant hazard simulates on the
    # PK time grid (hours).
    llambda <- log(0.0311); label("Baseline VOC hazard (1/day, time unit per published Figure 6 'Time (days)')")  # poster Results-PK-ER text: "base hazard estimate was 0.0311 (0.0140, 0.0481 95%CI)"
    lec50   <- log(574)   ; label("EC50 of IMR-687 average concentration on the VOC hazard reduction (ng/mL)")    # poster Results-PK-ER text: "EC50 for IMR-687 as 574 ng/mL (0, 1266 95%CI)"

    # ----- Inter-individual variability -----
    # Source: poster Results - PK-ER Model paragraph. Reported as %CV on the
    # natural scale; converted to omega^2 = log(CV^2 + 1) for the
    # log-transformed structural parameters.
    etalcl ~ log(0.223^2 + 1)  # 22.3% CV on CL/F (poster Results-PK-ER text)
    etalvc ~ log(0.162^2 + 1)  # 16.2% CV on V/F  (poster Results-PK-ER text)
    etalka ~ log(0.770^2 + 1)  # 77.0% CV on Ka   (poster Results-PK-ER text)

    # IIV on the baseline hazard. The poster's hazard equation defines
    # eta as "exponential individual random effect, N(0, omega^2)", but the
    # value of omega^2 is not reported. No etallambda is declared in ini();
    # downstream users who recover omega^2 from the control stream can add
    # etallambda ~ <var> and re-introduce + etallambda inside model() at
    # the lambda_day line.

    # ----- Residual error -----
    # Proportional residual error on Cc, expressed as %CV in the poster.
    propSd <- 0.409; label("Proportional residual error on observed IMR-687 plasma concentration (fraction; 40.9% CV)")  # poster Results-PK-ER text: "Residual proportional error was 40.9% CV"
  })

  model({
    # ============================================================
    # 1. Individual PK parameters
    # ============================================================
    # NB: the published covariate effects (WT on CL/F and V/F; FORM_CAPSULE,
    # daily DOSE, FED_HIGHFAT on Ka) are not encoded here because the source
    # poster does not provide the functional forms; see covariatesDataExcluded.
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    ka <- exp(lka + etalka)
    kel <- cl / vc

    # ============================================================
    # 2. RTTE structural parameters (with optional IIV on lambda)
    # ============================================================
    # lambda is the baseline VOC hazard with the published time unit of
    # 1/day; convert to 1/h so the constant hazard adds correctly across the
    # PK time grid (hours).
    lambda_day <- exp(llambda)
    lambda_h   <- lambda_day / 24
    ec50       <- exp(lec50)

    # ============================================================
    # 3. PK ODE system (one-compartment with first-order absorption)
    # ============================================================
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ============================================================
    # 4. Observation model (PK)
    # ============================================================
    # central state is in mg, vc in L, so central/vc is mg/L which is
    # numerically equal to ng/mL (1:1).
    Cc <- central / vc
    Cc ~ prop(propSd)

    # ============================================================
    # 5. RTTE hazard and cumulative hazard
    # ============================================================
    # Saturable (Michaelis-Menten) drug effect on hazard, driven by the
    # instantaneous plasma concentration Cc. The poster parameterises this
    # with Cavg(t) (the rolling-window average concentration); for QD dosing
    # at steady state Cc(t) oscillates around Cavg, and the instantaneous Cc
    # closely approximates Cavg once steady state is reached. See the
    # validation vignette's Assumptions and deviations section for the
    # Cc-vs-Cavg discussion.
    effect_voc <- 1 - Cc / (Cc + ec50)
    hazard_voc <- lambda_h * effect_voc
    d/dt(cumhaz) <- hazard_voc
    cumhaz(0) <- 0
    sur_voc <- exp(-cumhaz)
  })
}
