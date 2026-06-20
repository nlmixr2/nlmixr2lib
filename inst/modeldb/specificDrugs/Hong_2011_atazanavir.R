Hong_2011_atazanavir <- function() {
  description <- "C0-delinked one-compartment first-order-absorption population PK model with absorption lag-time for orally administered atazanavir (ATV) in HIV-infected adults and pediatric patients (3 months to 21 years), with covariate effects of age (ka), body weight (CL/F, V/F), sex, study-site region (Africa), ritonavir comedication (CL/F and Frel), and capsule-vs-powder formulation (Frel) (Hong 2011)."
  reference <- "Hong Y, Kowalski KG, Zhang J, Zhu L, Horga M, Bertz R, Pfister M, Roy A. Model-based approach for optimization of atazanavir dose recommendations for HIV-infected pediatric patients. Antimicrob Agents Chemother. 2011;55(12):5746-5752. doi:10.1128/aac.00554-11"
  vignette <- "Hong_2011_atazanavir"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at the steady-state PK assessment; used in allometric power scaling on CL/F (exponent 0.600) and V/F (exponent 0.706) with reference 70 kg (Hong 2011 Materials and Methods page 5747).",
      source_name        = "Body wt"
    ),
    AGE = list(
      description        = "Age at PK sampling",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used in power scaling on ka (exponent -0.822) with reference 18 years; the negative exponent encodes faster apparent absorption in younger pediatric subjects relative to adults (Hong 2011 Table 4 and Results page 5749).",
      source_name        = "Age"
    ),
    SEXF = list(
      description        = "Sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Hong 2011 Table 3 reports 104 of 227 subjects female (45.8%). Multiplicative linear-deviation effect on CL/F: cl *= (1 + e_sexf_cl * SEXF) with e_sexf_cl = -0.115 (Table 4), i.e. CL/F is 11.5% lower in females than males. Source paper categorical encoding was Sex (female), which directly maps to canonical SEXF without inversion.",
      source_name        = "Sex (female)"
    ),
    REGION_AFRICA = list(
      description        = "Study-site region indicator: 1 = Africa, 0 = otherwise (pooled North America + Europe in Hong 2011)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (pooled North America + Europe study sites)",
      notes              = "Hong 2011 Table 3 reports 91 of 227 subjects from African study sites; the remaining 136 are from North America (n = 127) and Europe (n = 9), pooled as the reference. Multiplicative linear-deviation effect on CL/F: cl *= (1 + e_region_africa_cl * REGION_AFRICA) with e_region_africa_cl = 0.145 (Table 4), i.e. CL/F is 14.5% higher at African sites than at the pooled North America / Europe reference.",
      source_name        = "Region (Africa)"
    ),
    CONMED_RTV = list(
      description        = "Concomitant low-dose ritonavir (RTV) coadministration indicator (pharmacokinetic booster)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (unboosted atazanavir regimen)",
      notes              = "Hong 2011 Table 3 reports 117 of 227 subjects (52%) on RTV (typically 100 mg QD adult, body-surface-area-scaled in pediatrics). Multiplicative linear-deviation effects: cl *= (1 + e_rtv_cl * CONMED_RTV) with e_rtv_cl = -0.409 (Table 4 -- 40.9% reduction in CL/F when RTV is coadministered, matching paper abstract); fdepot *= (1 + e_rtv_frel * CONMED_RTV) with e_rtv_frel = 1.32 (Table 4 -- 132% higher relative bioavailability with RTV).",
      source_name        = "RTV comedication"
    ),
    FORM_POWDER = list(
      description        = "Atazanavir oral powder formulation indicator (1 = powder, 0 = capsule reference)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (atazanavir capsule reference)",
      notes              = "Hong 2011 Table 3 reports 64 of 227 subjects (28%) on the pediatric powder formulation (all in PACTG1020 groups 1, 2, 5, 6); the remaining 163 are on capsules (Hong 2011 Table 2 stratification). Multiplicative linear-deviation effect on relative bioavailability: fdepot *= (1 + e_form_powder_frel * FORM_POWDER) with e_form_powder_frel = -0.355 (Table 4 -- 35.5% lower bioavailability for the powder relative to the capsule reference).",
      source_name        = "Formulation (powder)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 227L,
    n_studies      = 4L,
    n_observations = 3939L,
    age_range      = "0.33-64 years (adults 22-64, pediatric 0.33-21)",
    age_median     = "9.17 years (pediatric mean) / 36.3 years (adult mean)",
    weight_range   = "2.6-122 kg (adults 49.7-97.1, pediatric 2.6-122)",
    weight_median  = "33.5 kg (pediatric mean) / 70.7 kg (adult mean)",
    sex_female_pct = 45.8,
    race_ethnicity = c(White = 22.9, Black = 66.1, Other = 11.0),
    disease_state  = "HIV infection on antiretroviral therapy; pooled antiretroviral-treatment-naive and -experienced subjects.",
    dose_range     = "ATV 48-1200 mg QD (pediatric 48-600 mg QD on the BSA-titrated PACTG1020 protocol; adults 300, 400, or 600 mg QD) with or without RTV 26-100 mg QD.",
    regions        = "Africa (n = 91), North America (n = 127), Europe (n = 9)",
    notes          = "Pooled adult (3 studies: AI424008 n = 13, AI424089 n = 27, AI424137 n = 11) and pediatric (PACTG1020 n = 176) data. 3,939 ATV plasma concentrations (620 adult, 3,319 pediatric) collected at steady state from 4 clinical studies. The pediatric data span 8 stratification groups by age (infants 3 mo - 2 yr, children 2-13 yr, adolescents 13-21 yr), formulation (capsule vs powder), and RTV comedication (with vs without). See Hong 2011 Tables 1, 2, and 3 for the full demographic and dose-regimen breakdown."
  )

  ini({
    # Structural parameters at the reference covariate set (age = 18 yr, body weight = 70 kg, male,
    # non-Africa study site, atazanavir capsule, no concomitant ritonavir).
    # Reported in Hong 2011 Table 4 ("Final model parameter estimates"), page 5750.
    lka       <- log(2.04);   label("First-order absorption rate constant at reference age 18 yr (ka, 1/h)")             # Table 4: K_a = 2.04 +/- 0.31
    lcl       <- log(34.6);   label("Apparent oral clearance at reference covariates (CL/F, L/h)")                       # Table 4: CL/F = 34.6 (paper Results page 5749 confirms reference is male, 70 kg, non-Africa, no RTV)
    lvc       <- log(266);    label("Apparent central volume of distribution at reference body weight 70 kg (V/F, L)")    # Table 4: V/F = 266 +/- 25
    ltlag     <- log(0.913);  label("Absorption lag time (tlag, h)")                                                      # Table 4: t_lag = 0.913 +/- 0.002
    lfdepot   <- fixed(log(1));  label("Reference relative bioavailability for atazanavir capsule alone (Frel, unitless)") # Table 4: F_rel = 1 (structural anchor; capsule-alone reference)

    # Continuous covariate power exponents (theta_TV * (x/x_REF)^e form)
    e_age_ka            <- -0.822;  label("Age power exponent on ka (reference 18 yr)")                                   # Table 4: Age ~ K_a = -0.822 +/- 0.139
    e_wt_vc             <- 0.706;   label("Body weight power exponent on V/F (reference 70 kg)")                          # Table 4: Body wt ~ V/F = 0.706 +/- 0.082
    e_wt_cl             <- 0.600;   label("Body weight power exponent on CL/F (reference 70 kg)")                         # Table 4: Body wt ~ CL/F = 0.600 +/- 0.083

    # Categorical / binary covariate effects (theta_TV * (1 + e * x) linear-deviation form)
    e_region_africa_cl  <- 0.145;   label("Region (Africa) fractional shift on CL/F (vs pooled North America + Europe reference)")  # Table 4: Region (Africa) ~ CL/F = 0.145 +/- 0.045
    e_sexf_cl           <- -0.115;  label("Female sex fractional shift on CL/F (vs male reference)")                                # Table 4: Sex (female) ~ CL/F = -0.115 +/- 0.035
    e_rtv_cl            <- -0.409;  label("Concomitant ritonavir fractional shift on CL/F (vs no-RTV reference)")                   # Table 4: Comedication (with RTV) ~ CL/F = -0.409 +/- 0.026
    e_form_powder_frel  <- -0.355;  label("Powder formulation fractional shift on Frel (vs capsule reference)")                     # Table 4: Formulation (powder) ~ F_rel = -0.355 +/- 0.100
    e_rtv_frel          <- 1.32;    label("Concomitant ritonavir fractional shift on Frel (vs no-RTV reference)")                    # Table 4: Comedication (with RTV) ~ F_rel = 1.32 +/- 0.24

    # Inter-individual variability (omega^2 = log(CV^2 + 1) for log-normal IIV).
    # Hong 2011 Table 4 reports IIV on ka, ke, and V/F in the source parameterisation;
    # the canonical nlmixr2lib parameterisation uses lcl and lvc instead of lke and lvc.
    # Translating via CL = ke * V (independent IIVs in the source) gives a correlated
    # bivariate block on (etalcl, etalvc):
    #   var(etalvc)         = log(1 + 0.425^2) = 0.1661     (Table 4: IIV V/F = 42.5%)
    #   var(etalcl)         = log(1 + 0.146^2) + log(1 + 0.425^2)
    #                        = 0.0211 + 0.1661 = 0.1872     (combined IIV ke + V/F on CL)
    #   cov(etalcl, etalvc) = log(1 + 0.425^2) = 0.1661     (shared V/F component)
    # This preserves the source paper's marginal CV%s on ke (14.6%), V/F (42.5%) and
    # the derived 45.4% CV on CL/F.
    etalka                ~ 1.3845                          # Table 4: IIV K_a = 173%; log(1 + 1.73^2) = 1.3845
    etalcl + etalvc       ~ c(0.1872, 0.1661, 0.1661)       # Table 4: IIV K_e = 14.6% and IIV V/F = 42.5% translated via CL = ke * V

    # Residual error. Hong 2011 reports two stratum-specific log-transform residual error
    # CV%s in Table 4 (ATV alone: 55.2%; ATV + RTV: 35.3%). nlmixr2 / rxode2 error models
    # require a single residual SD per output; the value below is the ATV + RTV CV%
    # (the regimen used in the recommended pediatric doses targeted by this paper). The
    # ATV-alone residual CV% (55.2%) is documented in the vignette Errata.
    propSd                <- 0.353;  label("Proportional residual error (CV, fraction) at ATV + RTV regimen")     # Table 4: Residual error (%CV), ATV + RTV = 35.3%
  })

  model({
    # Individual PK parameters
    # ka: age-power scaling with reference age 18 yr (Table 4 footnote, Materials and Methods page 5747)
    ka <- exp(lka + etalka) * (AGE / 18)^e_age_ka

    # CL/F: weight-power scaling plus three multiplicative categorical effects (Africa, female, RTV)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl *
          (1 + e_region_africa_cl * REGION_AFRICA) *
          (1 + e_sexf_cl          * SEXF) *
          (1 + e_rtv_cl           * CONMED_RTV)

    # V/F: weight-power scaling (no categorical covariates in the final model)
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc

    # Relative bioavailability: capsule + no-RTV anchor at 1, multiplicative shifts for
    # powder (-35.5%) and concomitant RTV (+132%).
    fdepot <- exp(lfdepot) *
              (1 + e_form_powder_frel * FORM_POWDER) *
              (1 + e_rtv_frel         * CONMED_RTV)

    # Absorption lag time
    tlag <- exp(ltlag)

    # Micro-constant for explicit ODE system
    kel <- cl / vc

    # ODE system: one-compartment with first-order absorption from depot
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Bioavailability and absorption lag on the depot compartment
    f(depot)    <- fdepot
    alag(depot) <- tlag

    # Observation: dose in mg, V/F in L, central amount in mg -> central / vc gives mg/L,
    # which equals ug/mL. Multiply by 1000 to express the readout in ng/mL (the unit Hong
    # 2011 Methods page 5747 reports: "ATV plasma concentration was determined ... HPLC-MS/MS"
    # and Table 4 lists C_0 values in ng/mL).
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
