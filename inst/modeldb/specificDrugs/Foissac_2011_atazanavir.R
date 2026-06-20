Foissac_2011_atazanavir <- function() {
  description <- paste(
    "One-compartment first-order-absorption population PK model for orally",
    "administered atazanavir in 51 HIV-1-infected children and adolescents",
    "(3-18 years, 13-79 kg) on therapeutic drug monitoring. Body weight is",
    "carried through a fixed-exponent allometric scaling on CL/F (0.75) and",
    "V/F (1.0) referenced to 70 kg. Two binary co-medication indicators enter",
    "linearly on apparent oral clearance: low-dose ritonavir as a PK booster",
    "reduces CL/F (the typical CL/F = 7.1 L/h is the RTV-boosted reference,",
    "and absence of ritonavir multiplies CL/F by 1.80) and concomitant 300 mg",
    "tenofovir disoproxil fumarate increases CL/F by 25%. Between-subject",
    "variability is retained only on CL/F; residual error is proportional",
    "(Foissac 2011)."
  )
  reference <- paste(
    "Foissac F, Blanche S, Dollfus C, Hirt D, Firtion G, Laurent C, Treluyer JM,",
    "Urien S. Population pharmacokinetics of atazanavir/ritonavir in HIV-1-infected",
    "children and adolescents. Br J Clin Pharmacol. 2011;72(6):940-947.",
    "doi:10.1111/j.1365-2125.2011.04035.x."
  )
  vignette <- "Foissac_2011_atazanavir"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling on CL/F (exponent 0.75, FIXED) and V/F (exponent 1, FIXED), both referenced to 70 kg. The model is built on paediatric data (median 52 kg, range 13-79 kg per Table 1) but extrapolates to adults at WT = 70 kg by construction. Paper Table 2 footnote: [Typical value] = [Typical parameter] * (bodyweight/70)^PWR.",
      source_name        = "BW"
    ),
    CONMED_RTV = list(
      description        = "Concomitant low-dose ritonavir (RTV) indicator (pharmacokinetic booster)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (boosted atazanavir/ritonavir; the cohort's most common regimen)",
      notes              = "1 = subject receives low-dose ritonavir (typically 100 mg q.d.) as a pharmacokinetic booster of atazanavir; 0 = unboosted ATV. The Foissac 2011 cohort included 39 ATV/r-treated children and 9 ATV-only-treated children (3 subjects received both regimens successively). The typical CL/F reported in Table 2 (7.1 L/h at 70 kg) is the RTV-positive reference, so the model uses the polarity-flipped indicator (1 - CONMED_RTV) to apply the theta_NO_RTV effect of +0.80 only when ritonavir is absent: CL/F = exp(lcl) * (1 + e_no_rtv_cl * (1 - CONMED_RTV)) * (WT/70)^0.75. Without RTV, CL/F is multiplied by 1.80 (giving 12.8 L/h at 70 kg per Table 2 footnote).",
      source_name        = "RTV"
    ),
    CONMED_TDF = list(
      description        = "Concomitant tenofovir disoproxil fumarate (TDF) indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant TDF)",
      notes              = "1 = subject receives concomitant TDF (typically 300 mg q.d.) as part of the NRTI backbone of the antiretroviral regimen; 0 = no TDF. 21 of the 39 ATV/r-treated children also received TDF (Foissac 2011 Results, Demographic data). TDF increases ATV/r apparent oral clearance by 25%: CL/F = exp(lcl) * (1 + e_tdf_cl * CONMED_TDF) * (WT/70)^0.75. With TDF (and RTV), CL/F is multiplied by 1.25 (giving 8.9 L/h at 70 kg per Table 2 footnote).",
      source_name        = "TDF"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Median 14 years (range 3-18). Screened as a Hill-equation effect on PK parameters per the Methods 'Population pharmacokinetic modelling of ATV/r' paragraph; not retained because the body-weight allometric scaling removed the residual age effect on PK parameters (Discussion paragraph 2)."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units       = "(binary)",
      type        = "binary",
      notes       = "25 girls and 26 boys (Table 1). Screened as a linear effect on the typical value of a given parameter per the Methods paragraph and not retained in the final model. Canonical SEXF (1 = female) matches the source paper's coding.",
      source_name = "SEX"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 51L,
    n_observations = 151L,
    age_range      = "3-18 years",
    age_median     = "14 years",
    weight_range   = "13-79 kg",
    weight_median  = "52 kg",
    sex_female_pct = 49,
    disease_state  = "HIV-1 infected children and adolescents experienced with antiretroviral therapy.",
    dose_range     = paste(
      "Oral atazanavir 150-600 mg q.d. (median 400 mg) given alone or boosted with ritonavir 100 mg q.d.",
      "(ATV/r dose 100-400 mg, median 300 mg). 21 of 51 children additionally received tenofovir disoproxil",
      "fumarate 150-300 mg q.d. (median 300 mg)."
    ),
    regions        = "France (Paris-Necker, Trousseau, Cochin/Saint-Vincent-de-Paul, and Louis Mourier hospitals; APHP)",
    notes          = paste(
      "Routine therapeutic-drug-monitoring data. Median 2 samples per patient (range 1-13) over a median 2.3-month",
      "follow-up (range 0-34 months); 151 ATV plasma concentrations across 51 subjects. Atazanavir assayed by HPLC",
      "with LOQ = 0.10 mg/L (interassay precision and bias < 15% and 5%, respectively, in the 0.05-5 mg/L calibration",
      "range). Four observations below the LOQ (< 3% of the dataset) were handled by setting them to LOQ/2;",
      "the M2 and M3 methods gave equivalent parameter estimates. Population PK fit with NONMEM VI (FOCE-INTERACTION).",
      "Model stability assessed by 1000 bootstrap analyses (Wings for Nonmem); bootstrap confidence intervals",
      "shown in Table 2 are reasonably narrow and exclude zero. 39 children were treated with ATV/r,",
      "9 with ATV alone, and 3 received both regimens successively."
    )
  )

  ini({
    # ===================================================================
    # Structural parameters (Table 2 final-model column).
    # Typical CL/F is reported for the RTV-positive (boosted) reference
    # stratum standardised to 70 kg. The theta_NO_RTV covariate flips on
    # when CONMED_RTV = 0 and scales CL/F by 1 + 0.80 = 1.80.
    # ===================================================================
    lcl <- log(7.1);  label("Apparent oral clearance at WT = 70 kg, CONMED_RTV = 1, CONMED_TDF = 0 (CL/F, L/h)") # Table 2: CL/F = 7.1 L/h (70 kg)^-1, RSE 8%
    lvc <- log(103);  label("Apparent volume of distribution at WT = 70 kg (V/F, L)")                            # Table 2: V/F  = 103 L (70 kg)^-1, RSE 19%
    lka <- log(0.44); label("First-order absorption rate constant (ka, 1/h)")                                    # Table 2: Ka   = 0.44 1/h, RSE 26%

    # ===================================================================
    # Allometric exponents on CL/F and V/F.
    # Paper Methods: 'from allometric scaling theory these are typically
    # 0.75 for clearance and 1 for volume of distribution'. Held FIXED at
    # those theoretical values rather than estimated.
    # ===================================================================
    allo_cl <- fixed(0.75); label("Allometric exponent on CL/F (unitless, FIXED)") # Methods: PWR = 0.75 for CL (allometric theory)
    allo_v  <- fixed(1);    label("Allometric exponent on V/F  (unitless, FIXED)") # Methods: PWR = 1    for V  (allometric theory)

    # ===================================================================
    # Covariate effects on CL/F (linear-deviation form per Table 2 footnote)
    #   CL/F = exp(lcl) * (1 + e_no_rtv_cl * (1 - CONMED_RTV))
    #                   * (1 + e_tdf_cl    * CONMED_TDF)
    #                   * (WT/70)^0.75
    # ===================================================================
    e_no_rtv_cl <- 0.80; label("CL/F fractional change when ritonavir is absent (unitless)") # Table 2: theta_NO_RTV (CL/F) = 0.80, RSE 24%
    e_tdf_cl    <- 0.25; label("CL/F fractional change when TDF is co-administered (unitless)") # Table 2: theta_TDF (CL/F) = 0.25, RSE 37%

    # ===================================================================
    # Between-subject variability. The paper reports BSV expressed as the
    # square root of the omega estimate (Methods paragraph), so BSV = SD of
    # eta and the variance in the nlmixr2 ini() block is BSV^2.
    # Only CL/F retained an IIV in the final model.
    # ===================================================================
    etalcl ~ 0.16^2 # Table 2: BSV(CL/F) = 0.16 (RSE 34%); variance = 0.16^2 = 0.0256

    # ===================================================================
    # Residual variability - proportional only (Results paragraph 2:
    # 'The proportional model for the residual variability ensured a good
    # adequacy between observed and predicted values.').
    # ===================================================================
    propSd <- 0.53; label("Proportional residual error (CV, fraction)") # Table 2: sigma_proportional = 0.53, RSE 15%
  })

  model({
    # Individual structural parameters with allometric body-weight scaling
    # and the two co-medication effects on CL/F (Table 2 footnote).
    cl <- exp(lcl + etalcl) *
          (1 + e_no_rtv_cl * (1 - CONMED_RTV)) *
          (1 + e_tdf_cl    * CONMED_TDF) *
          (WT / 70)^allo_cl
    vc <- exp(lvc) * (WT / 70)^allo_v
    ka <- exp(lka)

    # Micro-constant
    kel <- cl / vc

    # ODE system: one-compartment first-order oral absorption.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Observation (dose in mg, V in L -> central / vc in mg/L) with
    # proportional residual error.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
