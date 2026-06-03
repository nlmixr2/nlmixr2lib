Chen_2017_tacrolimus <- function() {
  description <- "One-compartment population PK model with first-order absorption and absorption lag for low-dose oral tacrolimus (FK506, Prograf 0.5 mg capsules) in Chinese adult and paediatric myasthenia-gravis (MG) patients (Chen 2017). The absorption parameters ka and tlag are fixed at values obtained from a supplementary dataset of healthy volunteers, because the sparse-trough MG dataset is not informative about the absorption phase. Apparent oral clearance CL/F (3.6 L/h typical) is modulated by hematocrit and blood urea nitrogen through a multiplicative power-of-covariate-ratio form referenced to cohort medians (HCT median 38.4 %, exponent 4.31; BUN median 4.2 mmol/L, exponent 1.42). Apparent volume V/F is 1700 L typical with no retained covariate effects (high-dose IV immunoglobulin treatment was tested as a covariate on V/F but did not survive backward elimination). Inter-individual variability is diagonal on CL/F (141.6% CV) and V/F (72.4% CV); no IIV is estimated on ka or tlag. Residual variability is a pure proportional model (35.8% CV) on whole-blood tacrolimus concentrations."
  reference   <- "Chen YS, Liu ZQ, Chen R, Wang L, Huang L, Zhu X, Zhou TY, Lu W, Ma P. Population pharmacokinetic analysis of tacrolimus in Chinese myasthenia gravis patients. Acta Pharmacol Sin. 2017;38(8):1195-1204. doi:10.1038/aps.2016.174"
  vignette    <- "Chen_2017_tacrolimus"
  units       <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    HCT = list(
      description        = "Hematocrit, expressed as a percentage of total blood volume.",
      units              = "%",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying within subject (routine clinical-chemistry panel at each admission). Chen 2017 enters HCT into CL/F as the power-of-ratio (HCT/38.4)^4.31 with the cohort median 38.4% as the reference; the exponent is positive so higher hematocrit increases apparent oral clearance. Cohort median 38.4% (range 29.6 - 55.3%; Table 1; Table 1's column header 'Hemoglobin (g/L)-HCT' is a publication typo -- the printed values are unambiguously percent volume fraction, not hemoglobin g/L, consistent with the canonical-register HCT unit). Tacrolimus partitions strongly into erythrocytes (~95% red-blood-cell-bound in whole blood); the positive exponent is unusual relative to the inverse erythrocyte-binding hypothesis used by transplant-cohort tacrolimus papers and is discussed in the vignette Assumptions and deviations.",
      source_name        = "HCT"
    ),
    BUN = list(
      description        = "Blood urea nitrogen concentration.",
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying within subject (routine clinical-chemistry panel at each admission). Chen 2017 enters BUN into CL/F as the power-of-ratio (BUN/4.2)^1.42 with the cohort median 4.2 mmol/L as the reference; the exponent is positive so higher BUN increases apparent oral clearance. Cohort median 4.2 mmol/L (range 1.7 - 10.4 mmol/L; Table 1). The paper interprets this mechanistically as urea-driven protein carbamylation reducing albumin-binding of the ~99%-bound tacrolimus, freeing more drug for clearance (Chen 2017 Discussion).",
      source_name        = "BUN"
    )
  )

  covariatesDataExcluded <- list(
    PRO = list(
      description = "High-dose intravenous immunoglobulin treatment indicator: 1 = patient received the PRO co-treatment at the dose event; 0 = patient did not.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened in forward inclusion as a covariate on V/F (Table 3 model 5: dOFV = -5.9, p = 0.015) but did NOT survive backward elimination at the 6.64 OFV threshold; not retained in the final model. Discussion notes that with a covariate frequency of only 7.9% the likelihood-ratio test is underpowered and PRO 'has potential to become an influential covariate on V/F'. Documented here for provenance but not referenced in model()."
    ),
    LP = list(
      description = "Lymphocyte percentage.",
      units       = "%",
      type        = "continuous",
      notes       = "Screened as a covariate on CL/F (Table 3 model 4: dOFV = -8.2, p = 0.0042) but was eliminated after a biological-plausibility / clinical-importance review (Chen 2017 Results 'Population pharmacokinetic model'): no published pharmacokinetic mechanism links lymphocyte percentage to tacrolimus clearance, and removing LP further reduced IIV-CL. Documented here for provenance but not referenced in model()."
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 83L,
    n_studies        = 1L,
    age_range        = "2 - 81 years",
    age_median       = "30 years",
    weight_range     = "13.5 - 110 kg",
    weight_median    = "55 kg",
    sex_female_pct   = 60.2,
    sex_distribution = "50 female / 33 male (60.2% female / 39.8% male).",
    race_ethnicity   = "Chinese (single-country: Beijing, China; General Hospital of PLA Rocket Force).",
    disease_state    = "Myasthenia gravis (MG). Disease grading by the Osserman and Genkins modified classification: Class I 14 (16.9%), Class IIA 20 (24.1%), Class IIB 36 (43.4%), Class III 4 (4.8%), Class IV 7 (8.4%), Class V 2 (2.4%). All subjects were clinically diagnosed with MG and were receiving low-dose oral tacrolimus as their primary immunosuppressant; the cohort is operationally distinct from transplant-recipient tacrolimus cohorts because the dose level is markedly lower (institutional therapeutic window 4 - 8 ng/mL trough vs. 5 - 20 ng/mL in transplant indications).",
    dose_range       = "Tacrolimus (Prograf, FK506; Astellas Ireland, Killorglin, Co. Kerry) 0.5 mg capsules administered orally every 12 h. Starting dose 1 mg/day for adults and 0.5 mg/day for paediatric patients; subsequent doses titrated by therapeutic drug monitoring (TDM) in 0.5 mg/day increments. Median observations per patient 2 (range 1 - 14, mean 3).",
    co_medications   = "Half of the dosing events were accompanied by pyridostigmine for symptomatic MG control (DDIA indicator); ~12% of patients were transitioning from hormone-dependent to tacrolimus-dependent immunotherapy and received concomitant corticosteroids (prednisone or methylprednisolone; DDIB indicator). DDIA and DDIB were screened as covariates and were not retained. Antihypertensives and hypoglycemics were present in less than 5% of the cohort and were not included in the DDI covariate model.",
    n_observations   = 253L,
    sampling_design  = "Retrospective TDM cohort. Trough whole-blood tacrolimus concentrations were measured by validated LC-MS/MS (Shim-pack VP-ODS column 5 um, 150 x 4.6 mm; linear range 0.1 - 25 ng/mL; LLOQ 0.1 ng/mL; LOD 0.05 ng/mL; intra-day CV 3.37%, inter-day CV 4.09%). Routine blood and biochemical panels were drawn at each admission. TDM was performed 4 days after the first dose; most patients had 1 - 2 observations per admission and ~4 admissions over the study window. The absorption-phase dataset used to fix ka and tlag came from a separate healthy-volunteer cohort (the paper's Supplementary Figure S1; not included in this 83-subject cohort count).",
    regions          = "China (single-centre, Beijing).",
    notes            = "Retrospective single-centre cohort enrolled January 2011 - May 2015 at the General Hospital of PLA Rocket Force, Beijing. The empirical MG therapeutic window in Chinese hospitals is 4 - 8 ng/mL trough; this is markedly lower than transplant-cohort targets and motivates a dedicated MG popPK model. One individual was excluded for highly unlikely dosing records (suspected human-error data entry) and three observation records were excluded for one individual whose covariate information was entirely missing. BMI and weight were imputed for 35 missing records: nearest-available-admission value for subjects with at least one populated record, otherwise the cohort median."
  )

  ini({
    # Final-model fixed-effect parameter estimates from Chen 2017 Table 2 (column
    # "Final model -- Estimate"). Time in hours, apparent oral clearance CL/F in
    # L/h, apparent oral volume V/F in L, ka in 1/h, tlag in h.
    #
    # Reference subject for the typical-value structural parameters: HCT = 38.4 %
    # and BUN = 4.2 mmol/L (cohort medians, Table 1), which makes both covariate
    # ratios equal to 1.0.
    lcl   <- log(3.6)             ; label("Apparent oral clearance CL/F (L/h) at HCT = 38.4 %, BUN = 4.2 mmol/L")  # Chen 2017 Table 2 final-model CL/F = 3.6 L/h (RSE 39%)
    lvc   <- log(1700)            ; label("Apparent oral volume of distribution V/F (L)")                          # Chen 2017 Table 2 final-model V/F = 1700 L (RSE 9%)
    lka   <- fixed(log(0.502))    ; label("Absorption rate constant ka (1/h); fixed from supplementary dataset")   # Chen 2017 Table 2 final-model Ka = 0.502 /h FIX
    ltlag <- fixed(log(0.346))    ; label("Absorption lag time tlag (h); fixed from supplementary dataset")        # Chen 2017 Table 2 final-model Tlag = 0.346 h FIX

    # Covariate effects on CL/F (Chen 2017 Eq. 3, power-of-ratio form referenced
    # to the cohort medians from Table 1).
    e_bun_cl <- 1.42  ; label("Power exponent of (BUN / 4.2) on CL/F")  # Chen 2017 Table 2 final-model BUN_CL = 1.42 (RSE 16%)
    e_hct_cl <- 4.31  ; label("Power exponent of (HCT / 38.4) on CL/F") # Chen 2017 Table 2 final-model HCT_CL = 4.31 (RSE 19%)

    # Inter-individual variability -- Chen 2017 Methods state IIV was modelled
    # exponentially (log-normal eta). Final-model values reported in Table 2 as
    # CV%; converted to internal log-scale variance via omega^2 = log(1 + CV^2).
    # Chen 2017 Results explicitly note that the CL/F-V/F correlation present in
    # the base model disappeared after covariate inclusion, so the final IIV is
    # diagonal. No IIV is reported on ka or tlag (both fixed structural).
    etalcl ~ 1.100296  # Chen 2017 Table 2 final IIV-CL = 141.6% CV  -> log(1 + 1.416^2) = 1.100296
    etalvc ~ 0.421454  # Chen 2017 Table 2 final IIV-V  =  72.4% CV  -> log(1 + 0.724^2) = 0.421454

    # Residual unexplained variability -- Chen 2017 Methods describe testing
    # CCV (proportional), additive, and combined error models; Table 2 lists
    # only a proportional term in the final model, so the residual is pure
    # proportional with no additive component.
    propSd <- 0.358   ; label("Proportional residual error (fraction)")  # Chen 2017 Table 2 final-model Proportional = 35.8% (RSE 15%)
  })

  model({
    # Individual PK parameters. Chen 2017 Eq. 3 (final-model CL/F):
    #   CL/F_i = theta_CL/F * (HCT_i / HCT_median)^theta_HCT
    #                       * (BUN_i / BUN_median)^theta_BUN * exp(eta_i)
    # V/F retains no covariate effects in the final model (Table 3 backward
    # elimination removed PRO_V).
    cl   <- exp(lcl + etalcl) * (HCT / 38.4) ^ e_hct_cl * (BUN / 4.2) ^ e_bun_cl
    vc   <- exp(lvc + etalvc)
    ka   <- exp(lka)
    tlag <- exp(ltlag)

    kel <- cl / vc

    # One-compartment oral PK with first-order absorption and absorption lag.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    alag(depot) <- tlag

    # Tacrolimus whole-blood concentrations are reported in ng/mL (= ug/L).
    # Doses are in mg and vc is in L, so central/vc has units mg/L; multiply by
    # 1000 to convert to ng/mL.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
