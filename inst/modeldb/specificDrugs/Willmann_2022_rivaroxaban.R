Willmann_2022_rivaroxaban <- function() {
  description <- paste(
    "Adapted pediatric population PK model for rivaroxaban in 76 post-Fontan",
    "congenital heart disease patients aged 2-8 years (body weight 9.8-25.3 kg)",
    "from the UNIVERSE study. Structural framework is inherited from the",
    "EINSTEIN-Jr popPK model of Willmann 2021 (already extracted as",
    "Willmann_2021_rivaroxaban): 2-compartment disposition with first-order",
    "absorption from a depot, allometric body-weight scaling of CL, Vc, Vp,",
    "and Q centred on the 82.48 kg adult reference weight, with all structural",
    "parameters (ka, Vc, Vp, Q, and the three allometric exponents) fixed to",
    "the Willmann 2021 EINSTEIN-Jr estimates. The Willmann 2021 dose-dependent",
    "relative bioavailability function is replaced by a bimodal age-binned",
    "F1 (post-Fontan patients aged <5 years vs. >=5 years). Apparent CL,",
    "the two F1 values, IIV on CL and F1, and the proportional residual",
    "error are re-estimated on the 76-patient UNIVERSE dataset. The refined",
    "CL is 6.07 L/h at 82.48 kg vs. 8.02 L/h in EINSTEIN-Jr (24% lower in",
    "post-Fontan patients). This model was used to bridge doses for",
    "thromboprophylaxis in post-Fontan patients aged 9 years or older or",
    ">=30 kg (the target extrapolation population), leading to the US label",
    "of 7.5 mg once daily (30-<50 kg) and 10 mg once daily (>=50 kg)."
  )
  reference <- paste(
    "Willmann S, Ince I, Ahsman M, Coboeken K, Zhang Y, Thelen K, Kubitza D,",
    "Zannikos P, Zhou W, Pina LM, Post T, Lippert J. Model-informed bridging",
    "of rivaroxaban doses for thromboprophylaxis in pediatric patients aged",
    "9 years and older with congenital heart disease. CPT Pharmacometrics",
    "Syst Pharmacol. 2022;11(8):1111-1121. doi:10.1002/psp4.12830"
  )
  vignette <- "Willmann_2022_rivaroxaban"
  units    <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ug/L"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    depot       = list(analyte = "rivaroxaban", units = "mg", specimen = "administration site", verified = FALSE),
    central     = list(analyte = "rivaroxaban", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "rivaroxaban", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight (baseline; the UNIVERSE post-Fontan cohort spans 9.8-25.3 kg).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Fixed allometric scaling on CL (exponent 0.481), Vc and Vp",
        "(shared exponent 0.821), and Q (exponent 0.761), centred on the",
        "82.48 kg adult reference weight of the integrated EINSTEIN-Jr",
        "adult popPK analysis. The allometric structural parameters",
        "themselves (ka, Vc, Vp, Q, and the three exponents) are inherited",
        "unchanged from the Willmann 2021 EINSTEIN-Jr popPK model",
        "(doi:10.1002/psp4.12688); the Willmann 2022 paper's Table S2",
        "reports these values without SE (fixed to the EINSTEIN-Jr point",
        "estimates). Only CL is re-estimated for the post-Fontan",
        "population."
      ),
      source_name        = "WGHT"
    ),
    AGE = list(
      description        = "Baseline age (years); used as a binary cut at 5 years to select the two F1 values.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "In the paper's adapted popPK model, F1 is estimated separately for",
        "post-Fontan patients aged >=5 years (F1 = 0.752) versus <5 years",
        "(F1 = 1.20). The cutoff of 5 years was chosen by comparing objective",
        "function values across candidate cut-ages (Willmann 2022 Results",
        "'PopPK model qualification' paragraph 2: 'A cutoff of younger than",
        "5 years of age for CL provided the best improvement in the",
        "objective function value ... two F1 values with a cutoff of younger",
        "than 5 years of age provided the best improvement in OBJF'). AGE is",
        "used only as this binary switch, not as a continuous covariate."
      ),
      source_name        = "AGE"
    )
  )

  covariatesDataExcluded <- list(
    RACE_JAPANESE = list(
      description = "Japanese race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Explored on CL and F1 in the adapted popPK covariate search",
        "(Willmann 2022 Results 'PopPK model qualification' paragraph 3).",
        "'No statistically significant effect of Japanese patients on CL",
        "was detected, but a tendency toward slightly higher F1 in",
        "Japanese patients was seen by the model. The impact of a",
        "separate F1 estimate for Japanese patients aged 2 to younger",
        "than 5 years was tested but showed an inferior OBJF with the",
        "same number of parameters.' Effects in post-Fontan patients aged",
        "5 years and older could not be assessed due to only N=1",
        "Japanese subject in that subgroup. Not retained."
      )
    ),
    EGFR_SCHWARTZ = list(
      description = "Estimated glomerular filtration rate by the Schwartz formula",
      units       = "mL/min/1.73 m^2",
      type        = "continuous",
      notes       = paste(
        "Explored in Willmann 2022 Discussion paragraph 1 ('the estimated",
        "glomerular filtration rate values at baseline in the UNIVERSE",
        "study were in the normal range and very similar for the age",
        "groups 2 to younger than 5 years (geometric mean, 121.8",
        "ml/min/1.73m^2) and 5-8 years (126.9 ml/min/1.73m^2; data on",
        "file)'). Not retained."
      )
    ),
    CONMED_CYP3A4_INHIB = list(
      description = "CYP3A4 inhibitor co-medication indicator (any strength)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "'Concomitant medications that could potentially affect",
        "rivaroxaban exposure (e.g., CYP3A4 inhibitors) were reviewed but",
        "could not explain the findings' (Willmann 2022 Discussion",
        "paragraph 1). Not retained."
      )
    ),
    T_SINCE_FONTAN = list(
      description = "Time between end of Fontan procedure and start of rivaroxaban administration",
      units       = "days",
      type        = "continuous",
      notes       = paste(
        "Explored as a potential covariate for exposure in Willmann 2022",
        "Discussion paragraph 1, motivated by adult VTE-P studies where",
        "major orthopedic surgery reduces rivaroxaban CL by ~20% in the",
        "first 72 h post-surgery. 'However, the start of rivaroxaban",
        "treatment relative to the end of the Fontan procedure could",
        "also not explain the higher plasma exposure in post-Fontan",
        "patients aged between 2 and younger than 5 years.' Not retained."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 76L,
    n_studies      = 1L,
    trial_id       = "UNIVERSE (NCT02846532), post-Fontan cohort",
    age_range      = "2-8 years (76 post-Fontan patients; 52 aged 2 to <5 y, 24 aged 5-8 y)",
    weight_range   = "9.8-25.3 kg",
    sex_female_pct = NA_real_,
    race_ethnicity = paste(
      "8 of the 76 patients were Japanese; 7 were younger than 5 years",
      "of age (Willmann 2022 Data and Methods 'Bridging concept and",
      "model qualification' paragraph 3)."
    ),
    disease_state  = paste(
      "Pediatric congenital heart disease patients who had completed",
      "the Fontan procedure within 4 months. Received rivaroxaban as an",
      "oral suspension formulation for thromboprophylaxis over 12 months",
      "in the UNIVERSE study."
    ),
    dose_range     = paste(
      "Body-weight-adjusted rivaroxaban delivered twice-daily as an",
      "oral suspension in UNIVERSE (post-Fontan cohort 2-8 y, <30 kg);",
      "regimen designed to reproduce adult 10 mg once-daily exposure",
      "for thromboprophylaxis after major orthopedic surgery per the",
      "Fontan-PBPK model of Zhu et al. (Willmann 2022 Introduction",
      "paragraph 4)."
    ),
    regions        = "Multi-regional (UNIVERSE international phase III trial including a Japanese subgroup).",
    n_observations = "Not reported per model in Willmann 2022; UNIVERSE contributed 76 patients to the popPK re-fit; sparse PK sampling per UNIVERSE protocol.",
    upstream_model = paste(
      "The structural PK framework (2-compartment, ka, Vc, Vp, Q, and",
      "the three allometric WT exponents) is inherited unchanged from",
      "the EINSTEIN-Jr popPK model of Willmann 2021 (already extracted",
      "as Willmann_2021_rivaroxaban.R; doi:10.1002/psp4.12688).",
      "Willmann 2022 Table S2 lists the inherited parameters without",
      "standard errors, confirming they were held fixed during the",
      "UNIVERSE re-fit."
    ),
    forward_use    = paste(
      "For the forward extrapolation to post-Fontan patients aged 9-18",
      "years or >=30 kg (the target population of the paper), the",
      "original EINSTEIN-Jr popPK model (Willmann 2021, CL = 8.02 L/h",
      "at 82.48 kg, dose-dependent F1) was used rather than this",
      "adapted model, based on the finding that the adapted CL and F1",
      "estimates were needed only to fit the 2 to <5 year subgroup",
      "(Willmann 2022 Results 'PopPK model qualification' final",
      "paragraph). This adapted model itself describes the UNIVERSE",
      "post-Fontan 2-8 year retrospective fit."
    ),
    notes          = paste(
      "The paper also reports a Fontan-PBPK model with age-dependent",
      "CL adjustment factors (0.53, 0.64, 0.45 for children aged 4",
      "to <5, 3 to <4, and 2 to <3 years respectively; Willmann 2022",
      "Results 'PBPK model qualification' paragraph 1); the PBPK model",
      "is inherited from Zhu et al. 2022 (reference 20; not on disk)",
      "and its structural / physiologic parameters are not reproduced",
      "in Willmann 2022. Only the popPK model is extracted here."
    )
  )

  ini({
    # Structural PK parameters -- reference subject body weight 82.48 kg
    # (median of the integrated EINSTEIN-Jr adult popPK analysis).
    # Values from Willmann 2022 Table S2 'Parameter Name' column. All
    # non-CL structural parameters are inherited unchanged from Willmann
    # 2021 (EINSTEIN-Jr popPK) and are held fixed in the UNIVERSE
    # adaptation -- Table S2 reports these with '-' in the SE / RSE /
    # 95% CI columns, indicating they were not estimated.

    lka <- fixed(log(0.799))
    label("First-order absorption rate constant for tablet, granules, and diluted oral suspension (1/h; Willmann 2021 EINSTEIN-Jr)")
    # Willmann 2022 Table S2: KA (tablet, granules and diluted suspension) = 0.799 h^-1, fixed (SE reported as '-')

    lcl <- log(6.07)
    label("Apparent oral clearance CL for post-Fontan patients at 82.48 kg reference weight (L/h)")
    # Willmann 2022 Table S2: CL subject with WGHT of 82.48 kg = 6.07 L/h (SE 0.356, RSE 5.86%, 95% CI 5.37-6.76); re-estimated in UNIVERSE (vs. 8.02 L/h in Willmann 2021 EINSTEIN-Jr)

    lvc <- fixed(log(53.2))
    label("Central volume of distribution Vc at 82.48 kg reference weight (L; Willmann 2021 EINSTEIN-Jr)")
    # Willmann 2022 Table S2: Vc subject with WGHT of 82.48 kg = 53.2 L, fixed (SE reported as '-')

    lvp <- fixed(log(59.1))
    label("Peripheral volume of distribution Vp at 82.48 kg reference weight (L; Willmann 2021 EINSTEIN-Jr)")
    # Willmann 2022 Table S2: Vp subject with WGHT of 82.48 kg = 59.1 L, fixed (SE reported as '-')

    lq  <- fixed(log(2.50))
    label("Intercompartmental clearance Q at 82.48 kg reference weight (L/h; Willmann 2021 EINSTEIN-Jr)")
    # Willmann 2022 Table S2: Q subject with WGHT of 82.48 kg = 2.50 L/h, fixed (SE reported as '-')

    # Fixed allometric exponents inherited from Willmann 2021. Vc and Vp
    # share a single exponent per the EINSTEIN-Jr control stream
    # (Willmann 2021 supplement S2 THETA(5)); Willmann 2022 Table S2 also
    # reports a single 'Exponent to scale Vc and Vp on WGHT' row.
    e_wt_cl    <- fixed(0.481)
    label("Allometric WT exponent on CL (unitless; Willmann 2021 EINSTEIN-Jr)")
    # Willmann 2022 Table S2: Exponent to scale CL on WGHT = 0.481, fixed
    e_wt_vc_vp <- fixed(0.821)
    label("Allometric WT exponent shared across Vc and Vp (unitless; Willmann 2021 EINSTEIN-Jr)")
    # Willmann 2022 Table S2: Exponent to scale Vc and Vp on WGHT = 0.821, fixed
    e_wt_q     <- fixed(0.761)
    label("Allometric WT exponent on Q (unitless; Willmann 2021 EINSTEIN-Jr)")
    # Willmann 2022 Table S2: Exponent to scale Q on WGHT = 0.761, fixed

    # Age-binned relative bioavailability. In the adapted popPK the
    # Willmann 2021 dose-dependent F1 function is replaced by two
    # constants: F1(>=5 y) = 0.752 (reference) and F1(<5 y) = 1.20.
    # Values relative to the F1 = 1 reference of the EINSTEIN-Jr model
    # (Willmann 2022 Results 'PopPK model qualification' paragraph 2).
    lfdepot <- log(0.752)
    label("Log relative bioavailability F1 for post-Fontan patients aged >=5 years (unitless; reference age bin)")
    # Willmann 2022 Table S2: F1 >=5 years = 0.752 (SE 0.0818, RSE 10.9%, 95% CI 0.592-0.913)

    e_age_lt5y_f1 <- log(1.20 / 0.752)
    label("Log-ratio shift on F1 for post-Fontan patients aged <5 years (unitless)")
    # Willmann 2022 Table S2: F1 <5 years = 1.20 (SE 0.0845, RSE 7.05%, 95% CI 1.03-1.36); shift = log(1.20/0.752) = 0.4675

    # Inter-individual variability -- Willmann 2022 Table S2 'Variability'
    # block; exponential IIV, so the reported omega^2 is the log-normal
    # variance. %CV footnote: sqrt(exp(OM)-1)*100, matching
    # sqrt(exp(0.101)-1)*100 = 32.6% and sqrt(exp(0.161)-1)*100 = 41.8%.
    etalcl     ~ 0.101
    # Willmann 2022 Table S2: omega^2_CL (exponential) = 0.101 (SE 0.0332, RSE 32.8%, %CV 32.6)
    etalfdepot ~ 0.161
    # Willmann 2022 Table S2: omega^2_F1 (exponential) = 0.161 (SE 0.0418, RSE 26.0%, %CV 41.8)

    # Residual error -- Willmann 2022 Table S2 'Residual Error' block
    # reports sigma^2 (proportional) = 0.269 with stDev = sqrt(sigma^2)
    # = 0.519. In nlmixr2 the proportional residual takes an SD, so
    # propSd = sqrt(0.269) = 0.5187.
    propSd <- sqrt(0.269)
    label("Proportional residual SD for rivaroxaban plasma concentration (fraction)")
    # Willmann 2022 Table S2: sigma^2 proportional = 0.269 (SE 0.0219, RSE 8.16%, stDev 0.519)
  })

  model({
    # Individual PK parameters with allometric scaling. Reference
    # subject: WT = 82.48 kg. AGE is used only as a binary switch at
    # 5 years to select the F1 value for the two age bins.

    # Absorption: single ka (fixed) for the tablet / granules / diluted-
    # suspension arm. UNIVERSE post-Fontan patients received the oral
    # suspension formulation.
    ka <- exp(lka)

    cl <- exp(lcl + etalcl)  * (WT / 82.48)^e_wt_cl
    vc <- exp(lvc)           * (WT / 82.48)^e_wt_vc_vp
    vp <- exp(lvp)           * (WT / 82.48)^e_wt_vc_vp
    q  <- exp(lq)            * (WT / 82.48)^e_wt_q

    # Micro-constants for the explicit 2-compartment ODE system.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Age-binned relative bioavailability. AGE < 5 selects the <5 y
    # value (F1 = 1.20); AGE >= 5 uses the reference (F1 = 0.752).
    # IIV is applied multiplicatively (log-additive) to the selected
    # typical value.
    fdepot <- exp(lfdepot + etalfdepot + e_age_lt5y_f1 * (AGE < 5))

    # ODE system.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    # Bioavailability on the depot.
    f(depot) <- fdepot

    # Rivaroxaban plasma concentration. Dose mg, volume L -> mg/L;
    # multiply by 1000 to convert to ug/L (the units the paper's
    # exposure tables use throughout, e.g. Table S1 AUC in ug*h/L
    # and C_max, C_trough in ug/L).
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
