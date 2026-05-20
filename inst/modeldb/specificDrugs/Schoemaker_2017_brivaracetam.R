Schoemaker_2017_brivaracetam <- function() {
  description <- "One-compartment population PK model for oral brivaracetam in paediatric patients with epilepsy aged 1 month to 16 years (Schoemaker 2017). First-order absorption, single-compartment distribution, and first-order elimination, with allometric scaling of CL/F (exponent 0.750 fixed) and V/F (exponent 1.00 fixed) on lean body weight normalised to a 50 kg adult typical value. Co-administration of phenobarbital (PB; pooled with primidone), carbamazepine (CBZ), or valproate (VPA) modify apparent oral clearance via linear-additive multiplicative factors."
  reference   <- "Schoemaker R, Wade JR, Stockis A. Brivaracetam population pharmacokinetics in children with epilepsy aged 1 month to 16 years. Eur J Clin Pharmacol. 2017 Jun;73(6):727-733. doi:10.1007/s00228-017-2230-6"
  vignette    <- "Schoemaker_2017_brivaracetam"
  units       <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    LBM = list(
      description        = "Lean body weight (paper notation LBW), calculated from total body weight and body mass index per Janmahasatian et al.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Paper computes LBW from total body weight and BMI (Schoemaker 2017 Methods, page 2 / Data paragraph, citing Janmahasatian 2005). Allometric scaling exponents are fixed at the theoretical values (0.750 on CL/F and 1.00 on V/F) per Anderson & Holford. Reference value 50 kg corresponds to a typical adult lean body weight, chosen so the typical-CL / typical-V estimates match adult-cohort comparisons (Schoemaker 2017 Methods page 2).",
      source_name        = "LBW"
    ),
    CONMED_PB = list(
      description        = "Concomitant phenobarbital (PB) coadministration indicator: 1 = patient is on phenobarbital or primidone, 0 = neither.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant phenobarbital or primidone)",
      notes              = "Source paper pools primidone with phenobarbital because primidone is metabolised to phenobarbital. Linear-additive multiplicative effect on CL/F: cl *= (1 + 0.408 * CONMED_PB); +40.8% in clearance translates to ~29% lower brivaracetam exposure relative to PB-naive patients (Schoemaker 2017 Table 1 and Results paragraph 4).",
      source_name        = "PB"
    ),
    CONMED_CBZ = list(
      description        = "Concomitant carbamazepine (CBZ) coadministration indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant carbamazepine)",
      notes              = "Linear-additive multiplicative effect on CL/F: cl *= (1 + 0.479 * CONMED_CBZ); +47.9% in clearance translates to ~32% lower brivaracetam exposure relative to CBZ-naive patients (Schoemaker 2017 Table 1 and Results paragraph 4).",
      source_name        = "CBZ"
    ),
    CONMED_VPA = list(
      description        = "Concomitant valproate (VPA) coadministration indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant valproate)",
      notes              = "Linear-additive multiplicative effect on CL/F: cl *= (1 - 0.101 * CONMED_VPA); -10.1% in clearance translates to ~11% higher brivaracetam exposure relative to VPA-naive patients (Schoemaker 2017 Table 1). The authors retained the VPA effect in the final model for informational reporting although it did not strictly meet the SCM forward-selection p < 0.01 criterion; they note the apparent VPA-driven exposure rise may be confounded with VPA-associated weight / fat gain (Discussion paragraph 4).",
      source_name        = "VPA"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 96,
    n_studies      = 1,
    age_range      = "1 month to <16 years",
    age_strata     = "29 patients 1 month to <2 years; 26 patients 2 to <6 years; 24 patients 6 to <12 years; 17 patients 12 to <16 years",
    weight_range   = "Paediatric range across the 1 month - 16 years cohort (full demographics in Schoemaker 2017 Supplemental Table 1).",
    weight_median  = "Not reported in the main text; per Supplemental Table 1.",
    sex_female_pct = NA,
    race_ethnicity = "Not detected as significant covariate (race and ethnicity were tested and excluded during SCM).",
    disease_state  = "Paediatric epilepsy (localisation-related, generalised, or undetermined focal / generalised epileptic syndrome per ILAE classification). Patients were receiving 1 to 3 concomitant antiepileptic drugs other than levetiracetam.",
    dose_range     = "Brivaracetam oral solution as weekly increasing doses: 0.4 / 0.8 / 1.6 mg/kg bid for subjects >=8 years, 0.5 / 1.0 / 2.0 mg/kg bid for subjects <8 years. Doses capped at adult 25 / 50 / 100 mg bid for WT >= 50 kg.",
    regions        = "Multicentre (approximately 50 sites planned; specific regions not enumerated in the main text).",
    n_observations = "600 brivaracetam plasma concentrations across 3 sampling occasions (day 7, 14, 21) per subject in early-morning, late-morning, or afternoon brackets plus one optional sample.",
    co_medication  = "Concomitant AEDs: phenobarbital (PB; pooled with primidone), carbamazepine (CBZ), valproate (VPA), phenytoin (PHT; only 1 patient, excluded from SCM testing). Other tested covariates with no detected effect: race, ethnicity, sex, eGFR, non-AED CYP3A inhibitors, non-AED CYP2C19 inhibitors, age, postconceptional age (PCA).",
    notes          = "Trial N01263 (NCT00422422). Hepatic impairment was an exclusion criterion. Demographics from Schoemaker 2017 Results paragraph 1 and Supplemental Table 1. Bootstrap-based parameter uncertainty (n = 1000 replicates, 23 non-converged excluded) is reported alongside the NONMEM estimates in Table 1."
  )

  ini({
    # Structural parameters - final-model NONMEM estimates (Schoemaker 2017
    # Table 1). All clearance and volume terms are apparent (X/F) because
    # the study used oral solution and no IV reference arm; bioavailability
    # is not separately identifiable.
    lka <- log(1.84);  label("First-order absorption rate constant (Ka, 1/h)")               # Schoemaker 2017 Table 1: Ka = 1.84 (95% CI 0.91/2.78); bootstrap median 1.83
    lcl <- log(3.63);  label("Apparent oral clearance at 50 kg LBW reference (CL/F, L/h)")   # Schoemaker 2017 Table 1: CL/F = 3.63 (95% CI 3.42/3.85); bootstrap median 3.62
    lvc <- log(47.8);  label("Apparent central volume at 50 kg LBW reference (V/F, L)")      # Schoemaker 2017 Table 1: V/F = 47.8 (95% CI 43.1/52.5); bootstrap median 47.6

    # Allometric exponents on lean body weight (theoretical values held
    # fixed by the source authors per Anderson & Holford 2008; Schoemaker
    # 2017 Methods page 2).
    e_lbm_cl <- fixed(0.750); label("Allometric LBW exponent on CL/F (unitless, fixed)")     # Schoemaker 2017 Table 1: "Allometric scaling CL/F 0.750 fixed"
    e_lbm_vc <- fixed(1.00);  label("Allometric LBW exponent on V/F (unitless, fixed)")      # Schoemaker 2017 Table 1: "Allometric scaling V/F 1.00 fixed"

    # AED coadministration effects on CL/F. Linear-additive categorical
    # parameterisation matching the PsN SCM "linear" categorical form,
    # i.e. CL = TVCL * (1 + theta * indicator) for each AED.
    e_conmed_pb_cl  <-  0.408; label("Effect of CONMED_PB on CL/F (fraction)")               # Schoemaker 2017 Table 1: CL change with PB = +40.8% (95% CI +19.9%/+65.2%)
    e_conmed_cbz_cl <-  0.479; label("Effect of CONMED_CBZ on CL/F (fraction)")              # Schoemaker 2017 Table 1: CL change with CBZ = +47.9% (95% CI +27.8%/+71.2%)
    e_conmed_vpa_cl <- -0.101; label("Effect of CONMED_VPA on CL/F (fraction)")              # Schoemaker 2017 Table 1: CL change with VPA = -10.1% (95% CI -18.5%/-0.8%)

    # IIV. Source paper states "Exponential models were used to describe
    # the interindividual variability ... IIV was calculated as the
    # square root of the diagonal element in the omega matrix"
    # (Schoemaker 2017 Methods page 2). The Table 1 IIV column therefore
    # reports omega (log-scale SD); ini() takes the variance, so each
    # line below is omega^2.
    etalcl ~ 0.0520        # 0.228^2;  Schoemaker 2017 Table 1: IIV CL = 22.8% (shrinkage 6.1%)
    etalvc ~ 0.0279        # 0.167^2;  Schoemaker 2017 Table 1: IIV V  = 16.7% (shrinkage 45.6%)
    etalka ~ 0.1018        # 0.319^2;  Schoemaker 2017 Table 1: IIV Ka = 31.9% (shrinkage 73.4%)

    # Residual error. Source paper retained the proportional-only model
    # (combined add+prop offered no improvement, Results paragraph 2);
    # values reported as CV% in Schoemaker 2017 Table 1.
    propSd <- 0.234; label("Proportional residual error (fraction)")                         # Schoemaker 2017 Table 1: Residual error CV = 23.4% (95% CI 19.6%/27.1%)
  })

  model({
    # Reference lean body weight (paper LBW) for allometric scaling.
    # Schoemaker 2017 Methods, page 2: "a typical value of 50 kg" used
    # as POPCOV for LBW to allow easy comparison with reported adult
    # values.
    ref_lbm <- 50

    # AED-coadministration multiplier on CL/F. PsN SCM "linear"
    # categorical form: each AED contributes a multiplicative
    # (1 + theta * indicator) factor; the product captures the joint
    # effect for patients on multiple AEDs (Schoemaker 2017 caveats
    # interaction estimates for the single PB+CBZ patient, but the
    # model structure permits the product).
    aed_cl <- (1 + e_conmed_pb_cl  * CONMED_PB) *
              (1 + e_conmed_cbz_cl * CONMED_CBZ) *
              (1 + e_conmed_vpa_cl * CONMED_VPA)

    # Individual parameters.
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (LBM / ref_lbm)^e_lbm_cl * aed_cl
    vc <- exp(lvc + etalvc) * (LBM / ref_lbm)^e_lbm_vc

    # Micro-constant.
    kel <- cl / vc

    # ODE system: one-compartment with first-order absorption.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Observation. Dose in mg, volume in L -> mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
