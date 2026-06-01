Zeng_2009_acyclovir <- function() {
  description <- "One-compartment population PK model with first-order absorption for acyclovir in 43 children and young people (age 0.8-19.9 years; weight 7.3-70.2 kg) with malignancy, after intravenous acyclovir (5 mg/kg q8h, 1 h infusion) or oral valacyclovir prodrug (10 mg/kg q12h), developed in NONMEM v5.1.1 (FOCE-I) from 1216 plasma observations. Structural model: first-order absorption (ka) from a depot with bioavailability F (oral valacyclovir delivered as systemic acyclovir), one-compartment disposition with first-order elimination. Allometric body-weight scaling on CL (fixed exponent 0.75) and V (fixed exponent 1) referenced to the cohort median 19.6 kg; CL additionally varies with creatinine clearance via a power function (CRCL/106.7 mL/min/1.73 m^2)^FAC. Inter-individual variability is diagonal on CL, V, ka, and F. Residual error is a combined exponential (proportional after linearization) + additive model. Inter-occasion variability on CL (19.2% CV) and V (30.4% CV) reported by Zeng 2009 Table 3 is NOT encoded structurally here (per the Andrews 2017 / Brooks 2021 tacrolimus precedent) -- the source paper does not define an operational occasion column for the model-library use case."
  reference <- paste(
    "Zeng L, Nath CE, Blair EYL, Shaw PJ, Stephen K, Earl JW,",
    "Coakley JC, McLachlan AJ. (2009).",
    "Population pharmacokinetics of acyclovir in children and young",
    "people with malignancy after administration of intravenous",
    "acyclovir or oral valacyclovir.",
    "Antimicrob Agents Chemother 53(7):2918-2927.",
    "doi:10.1128/AAC.01138-08"
  )
  vignette <- "Zeng_2009_acyclovir"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight (baseline or time-varying)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling on CL (fixed exponent 0.75) and V (fixed exponent 1) with reference weight 19.6 kg (Zeng 2009 covariate-analysis paragraph: 'the population CL and V terms were standardized to 19.6 kg, which represents the median value of weight in this study group'). Cohort range 7.3-70.2 kg (Table 1).",
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Creatinine clearance estimated by the Counahan formula (CrCl = 0.43 * height_cm / serum_creatinine_mg_per_dL), expressed in the canonical mL/min/1.73 m^2 units.",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source paper computes CrCl with the Counahan formula natively in mL/min/1.73 m^2 (Methods, covariate analysis) but reports it in L/h/m^2 in Table 1 and uses 3.7 L/h/m^2 as the standardising reference. Converting to the canonical CRCL register units: 1 L/h/m^2 = (1000/60) * 1.73 = 28.83 mL/min/1.73 m^2, so 3.7 L/h/m^2 corresponds to 106.7 mL/min/1.73 m^2. The covariate effect is a power function `(CRCL/106.7)^e_crcl_cl`. Cohort range 2.0-5.7 L/h/m^2 = 57.7-164.4 mL/min/1.73 m^2 (Table 1).",
      source_name        = "CLCR"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 43L,
    n_studies      = 1L,
    age_range      = "0.8-19.9 years",
    age_median     = "6.3 years",
    weight_range   = "7.3-70.2 kg",
    weight_median  = "19.6 kg",
    sex_female_pct = 41.9,
    disease_state  = "Children and young people with malignancy (acute lymphoblastic leukemia n=16, acute myeloid leukemia n=6, neuroblastoma n=5, Wiskott-Aldrich syndrome n=3, Fanconi's anemia n=2, other diseases n=11) receiving acyclovir prophylactically against herpes simplex virus and varicella-zoster virus reactivation during chemotherapy or hematopoietic stem cell transplantation. Comedication with mycophenolate mofetil in 9/43 patients; tested as a covariate on CL but not retained.",
    dose_range     = "Intravenous acyclovir 5 mg/kg three times daily (1-h infusion); oral valacyclovir 10 mg/kg twice daily. 25 patients received IV only, 7 received oral only, 11 received both at different times.",
    regions        = "Single centre, Children's Hospital at Westmead, Sydney, Australia.",
    sampling_design = "1216 plasma acyclovir concentrations measured by validated HPLC (LOQ 0.1 mg/L; recovery 101%; intra- and inter-day precision <7% over 0.1-60 mg/L). Heavy children (>20 kg): 12-14 samples per dosing interval (intensive). Light children (<20 kg) and ad hoc samples: sparse design. Median 25 samples per patient (range 3-50). IV samples drawn 0-8 h post-infusion end; oral samples 0-12 h post-dose.",
    iov_structure   = "Inter-occasion variability (IOV) was identified on CL (19.2% CV) and V (30.4% CV) in addition to the diagonal IIV (Table 3 final-model). This model file does NOT encode IOV structurally -- the source paper defines an 'occasion' as 7 days for daily-administered patients but does not supply a per-sample OCC column convention for the model-library use case; the nlmixr2lib convention (Andrews 2017 / Brooks 2021 precedent) is to omit IOV when no operational occasion mapping is defined. Downstream users who want to simulate IOV can add an OCC indicator and per-occasion etas in rxode2.",
    notes          = "Demographics from Zeng 2009 Table 1. Acyclovir is eliminated predominantly by renal excretion (glomerular filtration + tubular secretion); approximately 10% of dose is metabolized in the liver. Valacyclovir is the L-valyl ester prodrug of acyclovir, rapidly hydrolyzed after oral administration to release systemic acyclovir; the F estimated here (0.60) is the bioavailability of acyclovir delivered via the oral valacyclovir prodrug, applied only to the depot (oral) route."
  )

  ini({
    # Structural parameters -- Zeng 2009 Table 3 'Final model' column.
    # Reference subject is a 19.6 kg patient with CrCl = 106.7 mL/min/1.73 m^2
    # (= 3.7 L/h/m^2, the cohort median). The four typical-value point
    # estimates were CL = 3.55 L/h, V = 7.36 L, ka = 0.63 1/h, F = 0.60.
    lka     <- log(0.63); label("Absorption rate constant from oral valacyclovir depot (1/h)")             # Zeng 2009 Table 3 Final ka  = 0.63 1/h
    lcl     <- log(3.55); label("Clearance for a 19.6 kg child at CRCL = 106.7 mL/min/1.73 m^2 (L/h)")    # Zeng 2009 Table 3 Final CL  = 3.55 L/h
    lvc     <- log(7.36); label("Central volume of distribution for a 19.6 kg child (L)")                  # Zeng 2009 Table 3 Final V   = 7.36 L
    lfdepot <- log(0.60); label("Bioavailability of acyclovir delivered via oral valacyclovir (unitless)") # Zeng 2009 Table 3 Final F   = 0.60

    # Allometric exponents on body weight, fixed a priori per the Anderson-
    # Holford paediatric size-scaling convention (Zeng 2009 covariate-analysis
    # paragraph: 'CL = theta_1 * (wt/19.6)^0.75 and V = theta_2 * (wt/19.6)';
    # the V exponent is implicitly 1). No RSE/CI reported for the exponents,
    # consistent with values held fixed during estimation.
    e_wt_cl <- fixed(0.75); label("Allometric (WT) exponent on CL (unitless)") # Zeng 2009 covariate-analysis paragraph
    e_wt_vc <- fixed(1.00); label("Allometric (WT) exponent on V (unitless)")  # Zeng 2009 covariate-analysis paragraph

    # Estimated covariate effect: power function of creatinine clearance on
    # CL. Zeng 2009 Table 3 Final 'RF factor' = 0.51 with RSE 27%; bootstrap
    # 95% CI 0.19-0.80. Centred at 106.7 mL/min/1.73 m^2 (= 3.7 L/h/m^2 in
    # the paper's reported units, the cohort median).
    e_crcl_cl <- 0.51; label("Power exponent of CRCL on CL via (CRCL/106.7)^e_crcl_cl (unitless)") # Zeng 2009 Table 3 Final 'RF factor' = 0.51

    # Inter-individual variability. Zeng 2009 reports IIV as CV percentages
    # consistent with the individual-PK formula theta_i = theta_mean *
    # exp(eta_i), eta_i ~ N(0, omega^2) and omega expressed as a percentage
    # on the log scale (i.e. CV% = omega * 100, NOT log(1 + CV^2)) -- the
    # same convention used in the recent Zhi 2018 and Chen 2021 extractions.
    # Back-transformations from Table 3 Final 'IIV (CV) (%)':
    #   CL : 23.6%  -> 0.236^2 = 0.055696
    #   V  : 35.9%  -> 0.359^2 = 0.128881
    #   ka : 58.1%  -> 0.581^2 = 0.337561
    #   F  : 41.8%  -> 0.418^2 = 0.174724
    etalcl     ~ 0.055696  # Zeng 2009 Table 3 Final omega(CL) = 23.6%
    etalvc     ~ 0.128881  # Zeng 2009 Table 3 Final omega(V)  = 35.9%
    etalka     ~ 0.337561  # Zeng 2009 Table 3 Final omega(ka) = 58.1%
    etalfdepot ~ 0.174724  # Zeng 2009 Table 3 Final omega(F)  = 41.8%

    # Residual error -- Zeng 2009 Methods (Base model building):
    #   Y_ij = Yhat_ij * exp(eps1_ij) + eps2_ij,
    #   eps1, eps2 ~ N(0, sigma^2).
    # Table 3 column heading is 'sigma_1' / 'sigma_2' (not sigma^2_1 / sigma^2_2),
    # matching the formula's sigma-as-SD parameterisation; the reported
    # values 0.26 (sigma_1) and 0.10 (sigma_2) are therefore interpreted as
    # standard deviations, NOT variances. For nlmixr2 we linearise the
    # exponential arm (exp(eps) ~ 1 + eps for small eps) to a proportional
    # error and keep the additive arm unchanged, following the Tanaka 2012
    # phenytoin precedent for the same combined-error specification.
    # Concentration units are mg/L; addSd is in the same units.
    propSd <- 0.26; label("Proportional residual SD (fraction)") # Zeng 2009 Table 3 Final sigma_1 = 0.26
    addSd  <- 0.10; label("Additive residual SD (mg/L)")          # Zeng 2009 Table 3 Final sigma_2 = 0.10
  })

  model({
    # Individual PK parameters: allometric scaling on body weight with
    # reference 19.6 kg (Zeng 2009 covariate paragraph and Table 2 model 4),
    # and a power-function CrCl effect on CL centred at 106.7 mL/min/1.73 m^2
    # (= 3.7 L/h/m^2 in the paper's reported units). The reference subject
    # at WT = 19.6 kg, CRCL = 106.7 has typical CL = 3.55 L/h, V = 7.36 L.
    cl     <- exp(lcl     + etalcl)     * (WT / 19.6)^e_wt_cl * (CRCL / 106.7)^e_crcl_cl
    vc     <- exp(lvc     + etalvc)     * (WT / 19.6)^e_wt_vc
    ka     <- exp(lka     + etalka)
    fdepot <- exp(lfdepot + etalfdepot)

    # One-compartment disposition with first-order absorption from depot.
    # IV acyclovir doses enter `central` directly (cmt = central in the user
    # dataset, F = 1 by default); oral valacyclovir doses enter `depot`
    # (cmt = depot) and absorb to central with rate ka and bioavailability
    # fdepot. ka and fdepot are therefore active only on the oral arm.
    kel <- cl / vc
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Bioavailability of acyclovir from oral valacyclovir applied to the
    # depot route only (IV path bypasses depot).
    f(depot) <- fdepot

    # Plasma concentration; dose in mg and V in L -> mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
