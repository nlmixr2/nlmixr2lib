Yu_2024_ilaprazole <- function() {
  description <- "Two-compartment population PK model with first-order elimination for ilaprazole, a proton-pump inhibitor, after 0.75 h intravenous infusion in Chinese healthy subjects and patients with duodenal ulcer (Yu 2024). Pooled analysis of 1,560 plasma concentrations from 58 subjects across four phase I studies (healthy) and one phase IIa study (duodenal ulcer), fit in Phoenix NLME 8.3 by FOCE-ELS. Female sex lowers clearance (exp(-0.213)) and duodenal-ulcer disease status raises both clearance (exp(0.290)) and peripheral volume (exp(0.356)); peripheral volume additionally scales with body weight by a power of 1.545 around a 60.6 kg median. Typical values are for a healthy male at 60.6 kg. Inter-individual variability on the inter-compartmental clearance (CLp) was fixed in the final model because of 84% eta-shrinkage and no numeric variance was reported, so it is encoded as fixed(0)."
  reference   <- "Yu M, Liu S, Wu X, Wang H. Population pharmacokinetic modeling of ilaprazole in healthy subjects and patients with duodenal ulcer in China. Front Pharmacol. 2024 Jan 10;14:1306222. doi:10.3389/fphar.2023.1306222."
  vignette    <- "Yu_2024_ilaprazole"
  units       <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline (time-fixed) body weight. Enters the final model as a power function on the peripheral volume of distribution only, normalised to the pooled-cohort median of 60.6 kg (Yu 2024 Eq. 9 and Table 2 'Total' column). Pooled-cohort median 60.6 kg with IQR 55.8-65.1 kg (Table 2). Weight was retained for covariate screening in preference to height and BMI because the WT-HT and WT-BMI correlations (R^2 = 0.810 and 0.638) exceeded the 0.5 collinearity threshold (Yu 2024 Section 3.2).",
      source_name        = "WT"
    ),
    SEXF = list(
      description        = "Female sex indicator: 1 = female, 0 = male.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "The source paper's `Sex` covariate is already coded 1 = female (Yu 2024 Table 3 footnote: 'Sex = 1 for female'), so it maps onto the canonical SEXF with no value transformation and the reported coefficient sign is carried through unchanged. Females have lower clearance (exp(-0.213) = 0.808, i.e. 19.2% lower CL), consistent with the paper's simulation result of 24.7% higher AUC0-t in females (Yu 2024 Section 3.3 and Figure 3). 48.3% of the pooled cohort was female (Table 2).",
      source_name        = "Sex"
    ),
    DIS_DUOD_ULCER = list(
      description        = "Duodenal-ulcer disease-state indicator: 1 = patient with duodenal ulcer enrolled in the phase IIa study (CTR20132846), 0 = healthy subject enrolled in one of the four phase I studies.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy subject; the pooled healthy cohort across CTR20132848 / CTR20140147 / CTR20150686 / CTR20150685)",
      notes              = "Time-fixed per subject. Maps 1:1 onto the paper's `Disease status` flag, which Yu 2024 Table 3 footnote codes as 1 = duodenal ulcer with the healthy cohort as reference, so the Table 3 coefficients and typical values are carried through unchanged with no sign flip or re-baselining. Duodenal-ulcer patients have higher clearance (exp(0.290) = 1.336) and a larger peripheral volume (exp(0.356) = 1.428) than healthy subjects, consistent with the paper's simulation result of a 26.92% lower AUC0-t in patients (Yu 2024 Section 3.3 and Figure 3). Enrolment required an ulcer diameter <= 15 mm with no combined ulcer bleeding. Cohort composition 48 healthy / 10 duodenal ulcer (Yu 2024 Section 3.1).",
      source_name        = "Disease status"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 58L,
    n_studies      = 5L,
    n_observations = 1560L,
    age_range      = "Pooled median 25 years (IQR 23-31); healthy-subject study medians 24-25.5 years, duodenal-ulcer study median 46 years (IQR 38.3-50.5)",
    age_median     = "25 years",
    weight_range   = "Pooled median 60.6 kg (IQR 55.8-65.1); per-study medians 58.1-63 kg",
    weight_median  = "60.6 kg",
    height_median  = "168 cm (IQR 160.5-172.3)",
    bmi_median     = "21.6 kg/m^2 (IQR 20.8-22.7)",
    sex_female_pct = 48.3,
    race_ethnicity = "Chinese (all subjects)",
    disease_state  = "48 healthy subjects (four phase I studies: CTR20132848 single-dose 4x4 crossover, CTR20140147 multiple-dose, CTR20150686 high-dose, CTR20150685 loading-dose) and 10 patients with duodenal ulcer (phase IIa CTR20132846; ulcer diameter <= 15 mm, no combined ulcer bleeding, no history of smoking or drinking)",
    renal_function = "Median creatinine clearance 118.1 mL/min (IQR 98.2-124.1) by Cockcroft-Gault; per-study medians 94.8-126 mL/min",
    dose_range     = "5, 10, 20 and 30 mg ilaprazole as intravenous infusion over 0.75 h; single-dose and once-daily multiple-dose regimens including a 20 mg loading dose followed by 10 mg maintenance doses on days 2-3",
    regions        = "China",
    notes          = "Demographics from Yu 2024 Table 2 ('Total' column); study designs and sampling schedules from Table 1. A total of 1,560 valid plasma concentrations were analysed after excluding 5 below-limit-of-quantification and 7 not-detected samples; the assay was linear over 1-1,000 ng/mL with a 1 ng/mL lower limit of quantification. Only the intravenous-infusion data were modelled; oral comparator arms and positive-control drug arms were not included."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters - Yu 2024 Table 3 final-model estimates.
    # The paper's `Disease status` flag (1 = duodenal ulcer) maps 1:1 onto
    # the canonical DIS_DUOD_ULCER, so every value below is the verbatim
    # Table 3 estimate at the paper's own reference state (healthy,
    # male, WT = 60.6 kg) with no sign flip or re-baselining.
    # ------------------------------------------------------------------
    lvc <- log(6.795)                ; label("Central volume of distribution Vc (L)")                                            # Yu 2024 Table 3: V = 6.795 L (RSE 5.230%, 95% CI 6.098-7.492); Eq. 7 has no covariates on V
    lcl <- log(3.394)                ; label("Clearance CL for a healthy male (L/h)")                                            # Yu 2024 Table 3: CL = 3.394 L/h (RSE 3.225%, 95% CI 3.180-3.609)
    lvp <- log(5.544)                ; label("Peripheral volume of distribution Vp for a healthy subject at WT = 60.6 kg (L)")   # Yu 2024 Table 3: Vp = 5.544 L (RSE 7.182%, 95% CI 4.763-6.326)
    lq  <- log(13.086)               ; label("Inter-compartmental clearance CLp (L/h)")                                          # Yu 2024 Table 3: CLp = 13.086 L/h (RSE 16.896%, 95% CI 8.749-17.423); no covariates

    # ------------------------------------------------------------------
    # Covariate effects - Yu 2024 Table 3 "Covariable effect" block,
    # applied exactly as printed in Eqs. 8-9.
    # ------------------------------------------------------------------
    e_sexf_cl           <- -0.213    ; label("Effect of SEXF (female) on log-CL (unitless)")                     # Yu 2024 Table 3 and Eq. 8: Sex effect on CL = -0.213 (RSE -23.398%, 95% CI -0.311 to -0.115), Sex = 1 for female
    e_dis_duod_ulcer_cl <-  0.290    ; label("Effect of DIS_DUOD_ULCER (duodenal ulcer) on log-CL (unitless)")   # Yu 2024 Table 3 and Eq. 8: Disease status effect on CL = 0.290 (RSE 23.751%, 95% CI 0.155-0.425), Disease status = 1 for duodenal ulcer
    e_dis_duod_ulcer_vp <-  0.356    ; label("Effect of DIS_DUOD_ULCER (duodenal ulcer) on log-Vp (unitless)")   # Yu 2024 Table 3 and Eq. 9: Disease status effect on Vp = 0.356 (RSE 14.707%, 95% CI 0.253-0.459), Disease status = 1 for duodenal ulcer
    e_wt_vp             <-  1.545    ; label("Power exponent of WT / 60.6 kg on Vp (unitless)")                  # Yu 2024 Table 3 and Eq. 9: WT effect on Vp = 1.545 (RSE 13.714%, 95% CI 1.129-1.960); median weight 60.6 kg per Eq. 9 and Table 2

    # ------------------------------------------------------------------
    # Inter-individual variability - Yu 2024 Table 3 reports omega^2
    # variances directly for the exponential IIV model of Eq. 1
    # (P_ik = P_popk * exp(eta_ik)), so the values below are variances on
    # the log scale and need no CV% conversion.
    # ------------------------------------------------------------------
    etalvc ~ 0.013                   # Yu 2024 Table 3: omega^2 V = 0.013 (RSE 31.284%, 95% CI 0.005-0.021; eta-shrinkage 34.0%)
    etalvp ~ 0.032                   # Yu 2024 Table 3: omega^2 Vp = 0.032 (RSE 21.763%, 95% CI 0.018-0.046; eta-shrinkage 1.82%)
    etalcl ~ 0.059                   # Yu 2024 Table 3: omega^2 CL = 0.059 (RSE 25.051%, 95% CI 0.030-0.089; eta-shrinkage 24.4%)
    etalq  ~ fixed(0)                # Yu 2024 Table 3: omega^2 CLp reported only as; Section 3.2 states the CLp random effect was because of 84% eta-shrinkage. No numeric variance is reported anywhere in the paper or supplement, so it is encoded as (0) (see vignette Errata).

    # ------------------------------------------------------------------
    # Residual error - multiplicative (proportional) model.
    # ------------------------------------------------------------------
    propSd <- 0.184                  ; label("Proportional residual error (fraction)")   # Yu 2024 Table 3: sigma_mult = 0.184 (RSE 9.131%, 95% CI 0.151-0.217); Section 3.2 states a multiplicative residual error model was used to characterize intra-individual variability
  })

  model({
    # 1. Individual PK parameters, transcribed directly from Yu 2024
    #    Eqs. 7-9 (`Disease status` = DIS_DUOD_ULCER, `Sex` = SEXF):
    #      Eq. 7: V  = V  * exp(eta_V)                                          -- no covariates
    #      Eq. 8: CL = CL * exp(-0.213 * Sex) * exp(0.29 * Disease status) * exp(eta_CL)
    #      Eq. 9: Vp = Vp * (Weight / 60.6)^1.545 * exp(0.356 * Disease status) * exp(eta_Vp)
    vc <- exp(lvc + etalvc)
    cl <- exp(lcl + etalcl + e_sexf_cl * SEXF + e_dis_duod_ulcer_cl * DIS_DUOD_ULCER)
    vp <- exp(lvp + etalvp + e_dis_duod_ulcer_vp * DIS_DUOD_ULCER) * (WT / 60.6)^e_wt_vp
    q  <- exp(lq + etalq)

    # 2. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 3. ODE system - two compartments, first-order elimination from
    #    central. Ilaprazole was given only as a 0.75 h intravenous
    #    infusion in every study contributing to this analysis, so dosing
    #    events target `central` with a rate or duration.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # 4. Observation. Amounts are in mg and volumes in L, so central / vc
    #    is mg/L; multiply by 1000 to report ng/mL as in the source paper
    #    (assay range 1-1,000 ng/mL, Yu 2024 Section 2.2).
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
