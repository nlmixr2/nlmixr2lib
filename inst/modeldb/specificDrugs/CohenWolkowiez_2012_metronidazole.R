CohenWolkowiez_2012_metronidazole <- function() {
  description <- "One-compartment IV population PK model for metronidazole in preterm infants (Cohen-Wolkowiez 2012). Clearance scales linearly with body weight (reference 1.5 kg) and as a power function of postmenstrual age (reference 32 weeks); central volume scales linearly with body weight."
  reference <- "Cohen-Wolkowiez M, Ouellet D, Smith PB, et al. Population pharmacokinetics of metronidazole evaluated using scavenged samples from preterm infants. Antimicrob Agents Chemother. 2012;56(4):1828-1837. doi:10.1128/AAC.06071-11"
  vignette <- "CohenWolkowiez_2012_metronidazole"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying; missing weights carried forward up to 7 days. Cohort median (range) 1.495 (0.678-3.850) kg. Reference 1.5 kg used in the linear weight scaling of CL and V (Cohen-Wolkowiez 2012 Table 4).",
      source_name        = "WT"
    ),
    PAGE = list(
      description        = "Postmenstrual age (gestational age in weeks / 4.35 + postnatal age in months)",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Cohort median (range) 32 (24-43) weeks PMA; converted to canonical months as PAGE = PMA_weeks / 4.35. Drives the CL maturation power term; reference is 32 weeks = 32/4.35 months. Reparameterised inside model() as (PAGE * 4.35 / 32)^theta_CL_PMA so the published 32-week reference is preserved.",
      source_name        = "PMA"
    )
  )

  covariatesDataExcluded <- list(
    SCAV = list(
      description = "Indicator (0 = blood draw, 1 = scavenged sample) used in the publication's residual-error / bias structure",
      units       = "(binary)",
      type        = "binary",
      notes       = "Cohen-Wolkowiez 2012 Table 4 reports a multiplicative scavenged-sample bias factor theta_SCAV = 0.713 (95% CI 0.581-0.899) and an elevated proportional residual error for scavenged samples (29.0 vs 13.5 CV% for blood draws). These are sample-quality artifacts of the scavenged sampling design rather than structural drug PK; this model retains only the blood-draw residual error so that simulations represent the underlying drug concentration time-course without the 30% scavenged-sample underestimation. The omitted artifacts are reported verbatim in the vignette 'Assumptions and deviations' section."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Univariable analysis (Cohen-Wolkowiez 2012 Table 3) showed CL = theta_CL * (WT/1.5) * (0.5/SCR)^theta_CL_SCR reduced OFV by 14.3; however, in the multivariable analysis SCR did not improve goodness of fit beyond PMA and was excluded from the final model. The authors attribute the SCR-CL association to a single SCR=4.7 mg/dL outlier and to the strong PMA-SCR correlation in preterm infants."
    ),
    BGA = list(
      description = "Gestational age at birth",
      units       = "weeks",
      type        = "continuous",
      notes       = "Screened during covariate analysis. Effect on CL was captured indirectly through PMA (PMA = BGA + PNA/7), so BGA is not a retained covariate in the final model. Cohort median (range) 27 (22-32) weeks."
    ),
    PNA = list(
      description = "Postnatal age",
      units       = "days",
      type        = "continuous",
      notes       = "Univariable analysis (Cohen-Wolkowiez 2012 Table 3) showed CL = theta_CL * (WT/1.5) * (PNA/57)^theta_CL_PNA reduced OFV by 11.6, but PMA produced a larger OFV drop and was retained in the final model. Cohort median (range) 41 (0-97) days."
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 32L,
    n_studies        = 1L,
    age_range        = "0-97 days postnatal; 22-32 weeks gestational age at birth; 24-43 weeks postmenstrual age",
    age_median       = "PNA 41 days; PMA 32 weeks; BGA 27 weeks",
    weight_range     = "678-3850 g",
    weight_median    = "1495 g",
    sex_female_pct   = 53,
    race_ethnicity   = c(White = 50, Black_or_other = 50),
    disease_state    = "Preterm infants (gestational age at birth <= 32 weeks; <120 days postnatal) receiving intravenous metronidazole as part of routine clinical care in the neonatal intensive care unit. Indications included anaerobic bacteremia, central nervous system infections, complicated intra-abdominal infections (e.g., necrotizing enterocolitis).",
    dose_range       = "Intravenous metronidazole at the dose and frequency prescribed by routine clinical care; cohort median dose 8 mg/kg (range 4-15 mg/kg); median dosing interval 12 h (range 5.9-48 h).",
    regions          = "United States (5 centers)",
    renal_function   = "Serum creatinine median (range) 0.5 (0.1-4.7) mg/dL",
    n_concentrations = 116L,
    n_concentrations_scavenged = 104L,
    notes            = "5-center prospective open-label PK study (Antimicrobial PK in High-Risk Infants trial, Pediatric Pharmacology Research Unit). 116 plasma metronidazole concentrations were used; 104/116 (90%) were scavenged from discarded routine clinical specimens and 12 were timed blood draws. NONMEM 7 / FOCE-I with WINGS for NONMEM 7.03. Baseline demographics per Cohen-Wolkowiez 2012 Table 2; 16 (50%) White, 4 (9%) Hispanic."
  )

  ini({
    # Structural parameters at the reference subject (1.5 kg body weight, 32 weeks postmenstrual age).
    lcl <- log(0.0397); label("Clearance at WT=1.5 kg, PMA=32 weeks (CL, L/h)")  # Cohen-Wolkowiez 2012 Table 4: theta_CL = 0.0397 (RSE 10.9%)
    lvc <- log(1.07);   label("Central volume at WT=1.5 kg (V, L)")              # Cohen-Wolkowiez 2012 Table 4: theta_V = 1.07 (RSE 15.0%)

    # Fixed allometric exponents on body weight (linear scaling, exponent = 1).
    # An estimated body-size exponent was tested by the authors and excluded for
    # lack of improvement in fit and imprecision (Cohen-Wolkowiez 2012 Results,
    # 'Population PK model building').
    e_wt_cl <- fixed(1); label("Linear weight exponent on CL (unitless, fixed at 1)")  # Cohen-Wolkowiez 2012 Table 3 (final model: linear (WT/1.5))
    e_wt_vc <- fixed(1); label("Linear weight exponent on V (unitless, fixed at 1)")   # Cohen-Wolkowiez 2012 Table 3 (final model: linear (WT/1.5))

    # Power exponent on (PMA/32) for CL. Reference PMA 32 weeks = 32/4.35 months.
    e_page_cl <- 2.49; label("Power exponent on (PMA/32) for CL (unitless)")  # Cohen-Wolkowiez 2012 Table 4: theta_CL-PMA = 2.49 (RSE 29.8%)

    # Inter-individual variability (Cohen-Wolkowiez 2012 Table 4).
    # omega^2 = log(CV^2 + 1) for log-normal eta on log-transformed CL.
    # omega^2 on CL: log(1 + 0.425^2) = 0.16608 (42.5 CV%).
    etalcl ~ 0.16608  # 42.5 CV%; Cohen-Wolkowiez 2012 Table 4 (RSE 28.5%)
    # No IIV on V: paper Table 4 reports only IIV on CL; after WT was incorporated
    # the V-IIV nonparametric estimate was close to zero and the parameter was
    # excluded from the final model (Cohen-Wolkowiez 2012 Results).

    # Proportional residual error from blood-draw samples (the clean measurement
    # condition). The published final model also reported an elevated
    # proportional residual error for scavenged samples (29.0 CV%); that and the
    # multiplicative scavenged-sample bias factor theta_SCAV = 0.713 are
    # documented in covariatesDataExcluded$SCAV and in the validation vignette.
    propSd <- 0.135; label("Proportional residual error from blood-draw samples (fraction)")  # Cohen-Wolkowiez 2012 Table 4: blood-draw sigma_1^2 = 13.5 CV% (RSE 24.5%)
  })
  model({
    # PMA reference: 32 weeks = 32/4.35 months (canonical PAGE is in months).
    page_ratio <- PAGE * 4.35 / 32

    # Individual PK parameters.
    cl <- exp(lcl + etalcl) * (WT / 1.5)^e_wt_cl * page_ratio^e_page_cl
    vc <- exp(lvc)          * (WT / 1.5)^e_wt_vc

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose in mg, volume in L -> central/vc has units mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
