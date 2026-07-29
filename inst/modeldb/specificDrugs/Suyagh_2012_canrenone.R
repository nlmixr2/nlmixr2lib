Suyagh_2012_canrenone <- function() {
  description <- paste(
    "One-compartment population PK model for canrenone, the pharmacologically",
    "active metabolite of intravenous potassium canrenoate (K-canrenoate), in",
    "23 paediatric patients (2 days to 10 years; median weight 4 kg, range",
    "2.16-28.0 kg) receiving K-canrenoate in the NICU / PICU for retained",
    "fluids or congestive heart failure (Suyagh 2012). The K-canrenoate dose",
    "compartment (modelled as 'depot' because only canrenone is measured) is",
    "converted to canrenone by a first-order metabolic transformation rate",
    "kf (paper symbol) = 5.25 1/h. Canrenone disposition is described by",
    "apparent clearance CL/F and apparent central volume V/F (the 'F' factor",
    "absorbs the unknown fraction of K-canrenoate that ultimately reaches",
    "the canrenone compartment; the model assumes the total dose is",
    "converted to canrenone). Bodyweight scales CL/F and V/F by fixed",
    "allometric exponents (0.75 on CL, 1.0 on V) with a reference weight of",
    "70 kg. No other covariate (gestational age, postnatal age, postmenstrual",
    "age, serum creatinine, serum albumin, haematocrit, sex) was retained",
    "in the final model. Residual error is proportional only."
  )
  reference <- paste(
    "Suyagh M, Hawwa AF, Collier PS, Millership JS, Kole P, Millar M,",
    "Shields MD, Halliday HL, McElnay JC (2012). Population pharmacokinetic",
    "model of canrenone after intravenous administration of potassium",
    "canrenoate to paediatric patients.",
    "British Journal of Clinical Pharmacology 74(5):864-872.",
    "doi:10.1111/j.1365-2125.2012.04257.x.",
    sep = " "
  )
  vignette <- "Suyagh_2012_canrenone"
  units <- list(time = "h", dosing = "umol", concentration = "umol/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline body weight used in the allometric covariate model on",
        "CL/F and V/F with fixed exponents 0.75 (CL) and 1.0 (V), scaled to",
        "a reference weight of 70 kg (Suyagh 2012 Methods + Table 3).",
        "Cohort median 4 kg (range 2.16-28.0 kg)."
      ),
      source_name        = "WT"
    )
  )

  covariatesDataExcluded <- list(
    GA = list(
      description = "Gestational age at birth (weeks)",
      units       = "weeks",
      type        = "continuous",
      notes       = "Screened on CL/F and V/F (Suyagh 2012 Step 2-3 forward inclusion / backward elimination); not retained in the final model."
    ),
    PNA = list(
      description = "Postnatal age",
      units       = "days",
      type        = "continuous",
      notes       = "Screened (Suyagh 2012 covariate screening); not retained. Co-linear with weight per Methods discussion."
    ),
    PMA = list(
      description = "Postmenstrual age (= GA + PNA + 2 weeks)",
      units       = "weeks",
      type        = "continuous",
      notes       = "Screened (Suyagh 2012 covariate screening); not retained. Co-linear with weight."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened (Suyagh 2012 covariate screening); not retained in the final model."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = paste(
        "Screened (Suyagh 2012 covariate screening); not retained in the",
        "final model. Canonical units standardized to SI g/L per the",
        "2026-06-19 canonical-register audit. The prior covariateData entry",
        "declared 'mg/dL', which is not a plausible reporting unit for serum",
        "albumin (a typical value of ~4 g/dL = ~40 g/L would be ~4000 mg/dL);",
        "this is read as a transcription error for the standard albumin unit",
        "g/dL. Conversion to SI: g/L = g/dL x 10 (so a typical ~4 g/dL maps to",
        "~40 g/L). No inline conversion is needed because ALB is an excluded",
        "covariate and is not referenced in model()/ini()."
      )
    ),
    HCT = list(
      description = "Haematocrit",
      units       = "%",
      type        = "continuous",
      notes       = "Screened (Suyagh 2012 covariate screening); not retained in the final model."
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened (Suyagh 2012 covariate screening); not retained in the final model. Cohort 15 male : 8 female."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 23L,
    n_studies      = 1L,
    age_range      = "2 days to 10 years (postnatal); PMA 37.3-574.0 weeks",
    age_median     = "PNA median 71 days (range 2-3738 days); 20 of 23 subjects were younger than 1 year (the other three were 2, 6 and 10 years).",
    weight_range   = "2.16-28.0 kg",
    weight_median  = "4 kg",
    sex_female_pct = 8 / 23 * 100,
    race_ethnicity = NA_character_,
    disease_state  = paste(
      "Paediatric patients in the NICU at the Royal Jubilee Maternity",
      "Service, Belfast or the medical and intensive care wards at the",
      "Royal Belfast Hospital for Sick Children. Indications: congestive",
      "heart failure (8/23, 34.8%) or PICU treatment of retained fluids,",
      "e.g. pulmonary oedema due to chronic lung disease (15/23, 65.2%)."
    ),
    dose_range     = paste(
      "Intravenous potassium canrenoate (K-canrenoate). Median dose 10.06",
      "umol (range 3.77-70.43 umol) per dose, given at clinician",
      "discretion. K-canrenoate molar mass 396.6 g/mol; the equivalent mass",
      "dose range is 1.50-27.93 mg per administration."
    ),
    regions        = "Northern Ireland (Belfast), UK",
    notes          = paste(
      "Demographics from Suyagh 2012 Table 1 (n = 23 patients; 101 plasma",
      "canrenone samples; 1-8 samples per patient, median 4). GA at birth",
      "median 37 weeks (range 25-41); GA group: 2 very preterm (< 32 wk),",
      "9 preterm (32-37 wk; all at 36 wk), 12 full-term. Final population",
      "parameter values from Table 3 (final model, fixed allometric",
      "exponents)."
    )
  )

  ini({
    # =====================================================================
    # Final population estimates from Suyagh 2012 Table 3 ("Final
    # pharmacokinetic model" column, p. 868). Estimation by FOCE with
    # interaction in NONMEM VI level 1.0 using ADVAN6 (Methods, p. 866).
    # Reference weight 70 kg; bodyweight enters with fixed allometric
    # exponents of 0.75 on CL/F and 1.0 on V/F (Methods + Results, p. 866
    # and Table 2). The model assumes total conversion of the
    # K-canrenoate dose into canrenone (Methods, p. 866); CL/F and V/F
    # are therefore apparent values relative to this F = 1 assumption.
    # =====================================================================

    # ---- Structural fixed effects ----
    lcl <- log(11.4)
    label("Apparent canrenone clearance CL/F at WT = 70 kg (L/h)")              # Table 3: theta_CL/F = 11.4 L h-1 70 kg-1 (RSE 10.3%)

    lvc <- log(374.2)
    label("Apparent canrenone central volume V/F at WT = 70 kg (L)")           # Table 3: theta_V/F = 374.2 L 70 kg-1 (RSE 20.04%)

    lka <- log(5.25)
    label("Apparent first-order K-canrenoate to canrenone transformation rate kf (1/h)") # Table 3: theta_kf = 5.25 h-1 (RSE 57.84%); paper symbol kf, mapped to canonical lka for the depot -> central transit

    # ---- Allometric exponents (fixed at canonical 0.75 / 1.0) ----
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent for WT on CL/F (unitless, fixed)")               # Methods p. 866 + Table 3: theta_1 = 0.75, fixed (canonical adult allometric value per Anderson / Holford; not estimated in the final model)

    e_wt_vc <- fixed(1.0)
    label("Allometric exponent for WT on V/F (unitless, fixed)")                # Methods p. 866 + Table 3: theta_2 = 1.0, fixed (canonical adult allometric value per Anderson / Holford; not estimated in the final model)

    # ---- IIV (log-scale variances; omega^2 = log(1 + CV^2)) ----
    # The paper reports IIV in CV% on the log-normal P_i = P_pop * exp(eta_i)
    # scale (Methods p. 866). Only diagonal omegas were retained (Results
    # p. 868: full covariance gave negligible MOFV improvement).
    etalcl ~ 0.15581
    label("IIV variance on log CL/F (paper omega_CL/F 41.06% CV per Table 3)")  # log(1 + 0.4106^2) = 0.15581

    etalvc ~ 0.19019
    label("IIV variance on log V/F (paper omega_V/F 45.76% CV per Table 3)")    # log(1 + 0.4576^2) = 0.19019

    # ---- Residual error ----
    # Combined additive + proportional was tested but the additive term
    # collapsed toward zero and was dropped, yielding a proportional-only
    # model (Results p. 868). Table 3 reports residual CV = 34.07%, which
    # is the proportional SD when applied as Cc ~ prop(propSd).
    propSd <- 0.3407
    label("Proportional residual SD on canrenone concentration (fraction)")     # Table 3: sigma_prop CV = 34.07% (RSE 23.08%)
  })

  model({
    # ---- Individual structural parameters ----
    # Fixed-exponent allometric scaling at reference weight 70 kg
    # (Suyagh 2012 Methods p. 866 + Table 3 final model).
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    ka <- exp(lka)                                # = kf (h^-1); no IIV reported in Table 3

    # ---- Micro-constants ----
    kel <- cl / vc                                # canrenone elimination rate (h^-1); 0.0305 h-1 at 70 kg, 0.0621 h-1 at 4 kg

    # ---- ODE system ----
    # K-canrenoate (the IV-dosed parent prodrug) is carried in the depot
    # compartment; only canrenone (central) is measured. Conversion is
    # mole-for-mole at rate ka == kf (paper Methods p. 866 and Figure 1
    # caption). Users dose K-canrenoate amounts in umol via cmt = depot;
    # the depot is structurally an "absorption" compartment but the
    # physical route is intravenous (the depot label reflects its
    # mathematical role, not the administration route).
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ---- Observation and residual error ----
    # Cc is canrenone plasma concentration in umol/L (paper Table 1 units).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
