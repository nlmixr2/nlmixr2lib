Han_2013_fluconazole <- function() {
  description <- "One-compartment IV population PK model for fluconazole in adult burn-ICU patients with suspected or confirmed Candida infection, with a piecewise CL covariate model that switches between a fixed CRRT-cohort CL and a Cockcroft-Gault-CrCl / postburn-recency / sepsis-shifted non-CRRT CL plus an additive WT / edema / postburn-recency model on volume (Han 2013)"
  reference <- "Han S, Kim J, Yim H, Hur J, Song W, Lee J, Jeon S, Hong T, Woo H, Yim DS. Population pharmacokinetic analysis of fluconazole to predict therapeutic outcome in burn patients with Candida infection. Antimicrob Agents Chemother. 2013;57(2):1006-1011. doi:10.1128/AAC.01372-12"
  vignette <- "Han_2013_fluconazole"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  paper_specific_etas <- c("etalcl_rrt")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column WT. Time-fixed at the subject's baseline weight (Han 2013 Table 1: mean 65.6 kg, range 40.0-90.0 kg, n = 60). Enters the V model as an additive linear term normalised to the population reference weight of 65 kg: V_typical = (WT/65) * theta5. The reference weight 65 was chosen by the authors to approximate the cohort mean (65.6 kg).",
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column CLCR. Computed by the Cockcroft-Gault equation in raw mL/min (NOT BSA-normalized to mL/min/1.73 m^2). Stored under the canonical CRCL column per inst/references/covariate-columns.md (CRCL accepts raw mL/min when the source paper does not apply BSA normalization, with the per-model description recording the assay form; precedent: Delattre 2010 amikacin, Shekar 2014 meropenem). Han 2013 cohort median is 123.5 mL/min (range 21.6-282.7); the model normalises CLCR to 120 mL/min inside the CL formula so the coefficient theta4 = 0.557 L/h has units of L/h per unit of (CLCR/120). CLCR is conventionally not defined for CRRT-dependent subjects; in Han 2013 the model formula switches off the CLCR-driven term when RRT_CRRT_STATUS = 1 (the CRRT arm uses a fixed CL = theta3 = 1.85 L/h).",
      source_name        = "CLCR"
    ),
    RRT_CRRT_STATUS = list(
      description        = "Continuous renal replacement therapy status indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Source column CRRT. 1 = subject was receiving continuous renal replacement therapy (CRRT) during the modeled PK sampling period; 0 = no CRRT. 15/60 patients were on CRRT in the Han 2013 cohort. Stored under the canonical RRT_CRRT_STATUS column per inst/references/covariate-columns.md (precedent: Shekar 2014 meropenem). The indicator is treated as time-fixed at the subject level; the model uses RRT_CRRT_STATUS as a switch between the non-RRT additive CL formula and the fixed CRRT-cohort CL theta3.",
      source_name        = "CRRT"
    ),
    DIS_SEPSIS = list(
      description        = "Active-sepsis co-condition indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Source column SEPS. 1 = clinically diagnosed sepsis at the start of the modeled PK interval; 0 = no sepsis. 29/60 subjects were septic in the Han 2013 cohort. Treated as time-fixed at the subject level. Enters the non-RRT CL formula as an additive shift theta9 = -0.369 L/h (septic patients have about 0.37 L/h lower CL in the non-RRT arm; consistent with the documented inflammation-driven reduction of glomerular filtration / fluconazole CL in sepsis -- see Pittrow & Penk 1999, cited as Han 2013 reference 18).",
      source_name        = "SEPS"
    ),
    DIS_EDEMA = list(
      description        = "Clinical edema presence indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Source column EDEM. 1 = clinical edema present (puffy face and pitting peripheral edema; Han 2013 Table 1 footnote b); 0 = no clinical edema. 20/60 subjects had edema in the Han 2013 cohort. Treated as time-fixed at the subject level. Enters the V formula as an additive shift theta6 = 13.6 L (edematous patients have about 13.6 L larger V; consistent with third-space fluid expansion increasing the apparent volume of a hydrophilic small molecule).",
      source_name        = "EDEM"
    ),
    DIS_BURN_RECENT = list(
      description        = "Recent-postburn hypermetabolic-phase indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Source column TBI. 1 = within 30 days of burn injury (the acute hypermetabolic phase); 0 = at least 30 days postburn. The 30-day cutoff was chosen by the authors because the postburn hypermetabolic response is maximised between days 7 and 17 after injury and substantially resolved by day 30 (Han 2013 Methods, citing references 12 and 13). The Han 2013 cohort mean time from burn injury was 23 days (range 7-88). Enters both the non-RRT CL formula (theta7 = 0.504 L/h additive shift) and the V formula (theta8 = 9.61 L additive shift) as time-fixed binary indicators -- recent-postburn patients have higher CL and larger V. NOTE: the Han 2013 source column is named TBI but does NOT refer to traumatic brain injury -- it is a binary recoding of time-from-burn-injury at the 30-day threshold.",
      source_name        = "TBI"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 60L,
    n_studies      = 1L,
    age_range      = "20-82 years",
    age_median     = "mean 50.3 years (range 20-82)",
    weight_range   = "40.0-90.0 kg",
    weight_median  = "mean 65.6 kg (range 40.0-90.0)",
    sex_female_pct = 17,
    race_ethnicity = "Not reported (single-centre Korean burn-ICU cohort)",
    disease_state  = "Adult burn-ICU patients (burns 11-95% of total body surface area; mean TBSA 50.3%) with suspected (47/60) or confirmed (13/60) Candida infection. APACHE II mean 7.6 (range 2-26).",
    dose_range     = "Fluconazole 100-400 mg IV infused over 10-40 minutes every 24 h for treatment of suspected or confirmed fungal infection. PK sampling was performed at steady state after 5-17 days (mean 6.5 days) of fluconazole therapy.",
    regions        = "South Korea (single-centre burn ICU at Hangang Sacred Heart Hospital, Hallym University, Seoul)",
    tbsa_range     = "11-95% (mean 50.3% of total body surface area burned)",
    apache_ii      = "mean 7.6 (range 2-26)",
    renal_function = "Cockcroft-Gault creatinine clearance mean 123.5 mL/min (range 21.6-282.7), raw mL/min, not BSA-normalized",
    rrt_status     = "15/60 (25%) on continuous renal replacement therapy (CRRT)",
    sepsis         = "29/60 (48%) septic at PK sampling",
    edema          = "20/60 (33%) with clinical edema (puffy face and pitting peripheral edema)",
    postburn_time  = "Mean 23 days from burn injury to fluconazole administration (range 7-88)",
    notes          = "Baseline demographics per Han 2013 Table 1. 60 burn-ICU patients enrolled between December 2008 and May 2010. Patients younger than 18, pregnant / breastfeeding, or allergic to fluconazole were excluded. 409 fluconazole concentration observations (8 venous samples per subject at 0, 3, 5, 9, 24, 27, 48, and 51 h after the start of infusion, with daily dosing maintained throughout). Plasma fluconazole assayed by LC-MS/MS (LLOQ 0.04 ug/mL; inter- and intra-day precision below 7.42% RSD). NONMEM 7.2 with FOCE-I."
  )

  ini({
    # Fixed effects from Han 2013 Table 3 (Final parameter estimates).
    # Final CL formula:
    #   CL = (1 - CRRT) * (theta1 + (CLCR/120)*theta4 + TBI*theta7 + SEPS*theta9) + CRRT * theta3
    # Final V formula (theta2 not estimable; dropped per Table 3 footnote b):
    #   V  = (WT/65)*theta5 + EDEM*theta6 + TBI*theta8
    # All thetas in linear space; theta1 / theta3 / theta5 log-transformed in ini()
    # for positivity, additive covariate shifts kept on the linear scale because
    # they can be negative (theta9 = -0.369 L/h for sepsis).
    lcl              <- log(0.693); label("CL intercept in nonseptic, >=30 days postburn patients without CRRT (L/h)") # Han 2013 Table 3: theta1 = 0.693 L/h (RSE 28.3%; bootstrap mean 0.622, 95% CI 0.116-1.03)
    lcl_rrt          <- log(1.85);  label("CL in patients on CRRT (L/h)")                                              # Han 2013 Table 3: theta3 = 1.85 L/h (RSE 4.3%; bootstrap mean 1.86, 95% CI 1.70-2.06)
    lvc              <- log(50.3);  label("V proportionality at WT = 65 kg, no edema, >=30 days postburn (L)")         # Han 2013 Table 3: theta5 = 50.3 L per (WT/65) (RSE 6.8%; bootstrap mean 50.6, 95% CI 43.9-59.0)

    # Covariate effects on the non-RRT CL arm (additive, L/h).
    e_crcl_cl            <- 0.557;  label("CL slope per (CLCR/120) on non-RRT arm (L/h)")                              # Han 2013 Table 3: theta4 = 0.557 L/h per unit (CLCR/120) (RSE 28.9%; bootstrap mean 0.593, 95% CI 0.180-0.993)
    e_dis_burn_recent_cl <- 0.504;  label("CL shift if within 30 postburn days (L/h)")                                 # Han 2013 Table 3: theta7 = 0.504 L/h (RSE 36.3%; bootstrap mean 0.580, 95% CI 0.293-0.946)
    e_dis_sepsis_cl      <- -0.369; label("CL shift if sepsis (L/h)")                                                  # Han 2013 Table 3: theta9 = -0.369 L/h (RSE 45.0%; bootstrap mean -0.390, 95% CI -0.0649 to -0.730)

    # Covariate effects on the V model (additive, L).
    e_dis_edema_vc       <- 13.6;   label("V shift if clinical edema (L)")                                             # Han 2013 Table 3: theta6 = 13.6 L (RSE 33.3%; bootstrap mean 13.5, 95% CI 5.35-22.6)
    e_dis_burn_recent_vc <- 9.61;   label("V shift if within 30 postburn days (L)")                                    # Han 2013 Table 3: theta8 = 9.61 L (RSE 43.1%; bootstrap mean 9.44, 95% CI 0.37-17.7)

    # Between-subject variability (Han 2013 Table 3 "Random effects (estimates, %CV)").
    # omega^2 on the internal log-normal scale is log((CV/100)^2 + 1).
    # The paper reports two distinct CL BSV estimates split by CRRT status:
    # omega1^2 applies to the non-RRT CL arm (via etalcl on cl_norrt),
    # omega3^2 applies to the RRT CL arm (via etalcl_rrt on cl_rrt).
    # etalcl_rrt is declared in paper_specific_etas at the top of the model
    # function because the conditional-IIV split-by-CRRT structure does not
    # follow the 1-to-1 eta<->lparam pairing rule.
    etalcl           ~ 0.12000  # log(0.357^2 + 1); Han 2013 Table 3: omega1^2 BSV of CL non-CRRT = 35.7 %CV (RSE 19.8%; bootstrap mean 32.4, 95% CI 24.8-40.4)
    etalvc           ~ 0.02983  # log(0.174^2 + 1); Han 2013 Table 3: omega2^2 BSV of V         = 17.4 %CV (RSE 36.5%; bootstrap mean 15.9, 95% CI 8.91-22.2)
    etalcl_rrt       ~ 0.02559  # log(0.161^2 + 1); Han 2013 Table 3: omega3^2 BSV of CL CRRT    = 16.1 %CV (RSE 30.5%; bootstrap mean 15.0, 95% CI 9.84-20.6)

    # Residual error (proportional only).
    # Han 2013 Methods retained only the proportional component of the combined
    # additive + proportional error model after model development; Han 2013
    # Results confirms "proportional residual error was chosen as the base PK model".
    # Table 3 sigma1^2 = 0.111 is the variance on the proportional NONMEM scale
    # (Y = F * (1 + EPS)); propSd in nlmixr2 is sqrt(sigma1^2) = sqrt(0.111) = 0.3332.
    propSd <- 0.3332; label("Proportional residual error (fraction)") # Han 2013 Table 3: sigma1^2 = 0.111 (RSE 11.0%; bootstrap mean 0.110, 95% CI 0.0883-0.135) -> propSd = sqrt(0.111)
  })
  model({
    # Piecewise CL: additive covariate model on the non-RRT arm, fixed theta3 on
    # the RRT arm. Each arm carries its own multiplicative log-normal IIV (etalcl
    # on non-RRT, etalcl_rrt on RRT) per Han 2013 Table 3.
    cl_norrt <- (exp(lcl) +
                   e_crcl_cl            * (CRCL / 120) +
                   e_dis_burn_recent_cl * DIS_BURN_RECENT +
                   e_dis_sepsis_cl      * DIS_SEPSIS) * exp(etalcl)
    cl_rrt   <- exp(lcl_rrt) * exp(etalcl_rrt)
    cl       <- (1 - RRT_CRRT_STATUS) * cl_norrt + RRT_CRRT_STATUS * cl_rrt

    # Additive covariate model on V (theta2 intercept dropped per Han 2013
    # Table 3 footnote b -- not estimable in the full model). Multiplicative
    # log-normal IIV via etalvc.
    vc <- ((WT / 65) * exp(lvc) +
             e_dis_edema_vc       * DIS_EDEMA +
             e_dis_burn_recent_vc * DIS_BURN_RECENT) * exp(etalvc)

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose in mg, vc in L; central / vc is mg/L = ug/mL (1 mg/L = 1 ug/mL).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
