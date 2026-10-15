Goyal_2022_vancomycin <- function() {
  description <- "Two-compartment IV population PK model for vancomycin in 34 hospitalized pregnant women (Goyal 2022), fitted to routine therapeutic-drug-monitoring (mostly trough) concentrations. Clearance scales linearly with uncapped creatinine clearance (reference 175 mL/min, exponent fixed at 1) and with fat-free mass to the 0.75 power (reference 45 kg); Vc and Vp scale linearly and Q to the 0.75 power with fat-free mass. Vp, Q and the creatinine-clearance exponent were fixed at the non-pregnant literature-based base model values; IIV is estimated on CL only."
  reference <- "Goyal RK, Moffett BS, Gobburu JVS, Al Mohajer M. Population Pharmacokinetics of Vancomycin in Pregnant Women. Front Pharmacol. 2022;13:873439. doi:10.3389/fphar.2022.873439"
  vignette <- "Goyal_2022_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance (raw, NOT BSA-normalized): Cockcroft-Gault for patients >= 19 years, modified Schwartz for patients < 19 years",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Goyal 2022 Methods 'Patients and Data Collection': CRCL 'was calculated by the modified Schwarz equation for patients < 19 years of age and the Cockroft-Gault equation for patients >= 19 years of age'; the result is reported in ml/min (Table 1) with no BSA normalization. NOT capped: the Discussion reports that capping at 120 or 150 ml/min raised the OFV (483 and 489 vs 472) and biased clearance upward, so the final model uses the raw value. Reference 175 mL/min in the power term (Table 2; Supplementary Code S1 '(CrCL/175)'); Table 1 cohort median 176 (range 43-389) ml/min.",
      source_name = "CrCL"
    ),
    FFM = list(
      description = "Fat-free mass: Janmahasatian et al. (2005) formula for patients >= 18 years, Al-Sallami et al. (2015) formula for patients < 18 years",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Goyal 2022 Methods 'Patients and Data Collection' (FFM derivation). Reference 45 kg on all four structural parameters (Table 2; Supplementary Code S1 '(FFM/45)'), equal to the Table 1 cohort median 45 (range 30-60) kg. FFM beat total body weight as the size descriptor by 12 OFV points (Discussion).",
      source_name = "FFM"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = "Goyal 2022 Discussion: TBW tested in place of FFM as the size covariate raised the OFV by 12 points and gave about 8% higher IIV on CL, so FFM was retained. Supplementary Figure S2e: no residual eta-CL trend. Table 1 median 74 (range 43-157) kg."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Goyal 2022 Supplementary Figure S2a: eta-CL vs age showed no significant correlation; not retained. Table 1 median 28 (range 17-38) years."
    ),
    GA = list(
      description = "Gestational age",
      units = "weeks",
      type = "continuous",
      notes = "Goyal 2022 Supplementary Figure S2b: eta-CL vs gestational age showed no significant correlation; not retained. Table 1 median 27 (range 7-40) weeks; 2 first-, 15 second- and 17 third-trimester patients (Results)."
    ),
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      notes = "Goyal 2022 Supplementary Figure S2c: eta-CL vs height showed no significant correlation; not retained (height enters only through the FFM formula). Table 1 median 163 (range 147-173) cm."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 34L,
    n_studies = 1L,
    age_range = "17-38 years",
    age_median = "28 years",
    weight_range = "43-157 kg",
    weight_median = "74 kg",
    height_median = "163 cm (range 147-173)",
    bmi_median = "28 kg/m^2 (range 19-70)",
    ffm_median = "45 kg (range 30-60)",
    sex_female_pct = 100,
    disease_state = "Hospitalized pregnant women receiving intravenous vancomycin (infection type not reported); gestational age 7-40 weeks (median 27): 2 first-, 15 second- and 17 third-trimester patients. Patients on renal replacement therapy were excluded.",
    ga_range = "7-40 weeks (median 27)",
    renal_function = "CRCL median 176 mL/min (range 43-389), uncapped; serum creatinine median 0.56 mg/dL (range 0.27-1.97). 22 patients had serum creatinine 0.4-0.8 mg/dL, 3 below 0.4 and 9 above 0.8 mg/dL.",
    dose_range = "Intravenous vancomycin by routine clinical dosing; median (IQR) total daily dose 3000 mg (2000-4000 mg). Infusion duration not reported.",
    regions = "United States (Texas Children's Hospital, Houston, TX)",
    n_concentrations = 82L,
    notes = "Goyal 2022 Table 1 and Results. Retrospective TDM database study, 1 January 2011 to 31 May 2019. 91 samples from 34 patients; 9 below the 5 mg/L LLOQ were excluded, leaving 82 (most were troughs, at least half within 2 h before a dose). Assay: VITROS VANC competitive enzyme immunoassay (G6P-DH), range 5-50 mg/L, CV < 6%. Fit in Pumas 2.0 with second-order Laplace (LaplaceI); truncated-normal residual error at LLOQ (M2). 1000-sample bootstrap (Table 2)."
  )

  ini({
    # Structural parameters (Goyal 2022 Table 2, 'Final PK Model' column) for
    # the reference patient with CRCL = 175 mL/min and FFM = 45 kg.
    lcl <- log(7.64)
    label("Clearance at CRCL = 175 mL/min and FFM = 45 kg (L/h)") # Goyal 2022 Table 2: CL = 7.64 L/h; bootstrap 95% CI 6.38-9.73
    lvc <- log(67.35)
    label("Central volume at FFM = 45 kg (L)") # Goyal 2022 Table 2: Vc = 67.35 L; bootstrap 95% CI 41.96-112.95
    lq <- fixed(log(9.064))
    label("Intercompartmental clearance at FFM = 45 kg (L/h)") # Goyal 2022 Table 2: Q = 9.06 (Fixed); Supplementary Code S1 constantcoef tvq = 9.064
    lvp <- fixed(log(37.5))
    label("Peripheral volume at FFM = 45 kg (L)") # Goyal 2022 Table 2: Vp = 37.5 (Fixed); Supplementary Code S1 constantcoef tvvp = 37.5

    # Covariate exponents (Goyal 2022 Table 2 'Formula' column and
    # Supplementary Code S1 @pre block).
    e_crcl_cl <- fixed(1)
    label("Power exponent on (CRCL/175) for CL (unitless)") # Goyal 2022 Table 2: theta_CRCL = 1.0 (Fixed); Supplementary Code S1 constantcoef exp_crcl = 1.0
    e_ffm_cl_q <- fixed(0.75)
    label("Allometric exponent on (FFM/45) for CL and Q (unitless)") # Goyal 2022 Table 2 formulas 'CL.(CRCL/175)^theta_CRCL.(FFM/45)^0.75' and 'Q.(FFM/45)^0.75'
    e_ffm_vc_vp <- fixed(1)
    label("Allometric exponent on (FFM/45) for Vc and Vp (unitless)") # Goyal 2022 Table 2 formulas 'Vc.(FFM/45)' and 'Vp.(FFM/45)'

    # IIV: exponential on CL only (Supplementary Code S1: CL = ... * exp(etaCL)).
    # Table 2 prints IIV as 31.9 CV%; omega^2 = log(1 + 0.319^2) = 0.0969.
    etalcl ~ 0.0969 # Goyal 2022 Table 2: IIV on CL 31.9 CV% (shrinkage 0.21)

    # Residual error: proportional, Supplementary Code S1
    # 'truncated(Normal(cp, abs(cp) * sigma_prop), 4.99999, Inf)'.
    propSd <- 0.321
    label("Proportional residual error (fraction)") # Goyal 2022 Table 2: Proportional Error 32.1% (bootstrap 95% CI 18.1-45.8)
  })
  model({
    # Individual parameters (Goyal 2022 Supplementary Code S1 @pre block).
    cl <- exp(lcl + etalcl) * (CRCL / 175)^e_crcl_cl * (FFM / 45)^e_ffm_cl_q
    vc <- exp(lvc) * (FFM / 45)^e_ffm_vc_vp
    vp <- exp(lvp) * (FFM / 45)^e_ffm_vc_vp
    q <- exp(lq) * (FFM / 45)^e_ffm_cl_q

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Pumas Central1Periph1: two-compartment, first-order elimination, IV dosing.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg, volume in L -> mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
