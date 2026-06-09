Lee_2017_doripenem <- function() {
  description <- "One-compartment IV-infusion population PK model for doripenem in 37 Korean adults with acute infections (pyelonephritis, intra-abdominal infection, neutropenic fever) and CLCR ranging 20-50 or >50 mL/min (Lee 2017). Clearance and central volume scale linearly with body weight (CL/WT = 0.109 L/h/kg, V/WT = 0.280 L/kg at WT=70 kg, CLCR=57 mL/min); CL additionally scales by a power exponent on Cockcroft-Gault creatinine clearance (raw mL/min, reference 57)."
  reference <- "Lee D-H, Kim YK, Jin K, Kang MJ, Joo Y-D, Kim YW, Moon YS, Shin J-G, Kiem S. Population pharmacokinetic analysis of doripenem after intravenous infusion in Korean patients with acute infections. Antimicrob Agents Chemother. 2017;61(5):e02185-16. doi:10.1128/AAC.02185-16"
  vignette <- "Lee_2017_doripenem"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Lee 2017 Table 1: total cohort mean 59.8 kg (SD 12.4); 250-mg group 55.1 kg (SD 13.3); 500-mg group 61.3 kg (SD 11.9). The paper's structural model is per-kg (CL/WT and V/WT), equivalent to fixed linear allometric scaling on WT for both CL and V. Reference WT = 70 kg is used in this implementation to align with the Caucasian comparison in Figure 1 (calculated for a 70-kg patient).",
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Lee 2017 Methods (Population PK analysis): CLCR computed by the Cockcroft-Gault equation in raw mL/min (NOT BSA-normalized). Table 1 cohort mean 66.7 mL/min (SD 34.4); 250-mg group 38.3 mL/min (SD 10.9); 500-mg group 75.9 mL/min (SD 34.5). Reference value 57 mL/min from the structural-model equation CL = 0.109 x WT x (CLCR/57)^0.688 (Table 2 / Results). Stored under canonical CRCL with raw mL/min per inst/references/covariate-columns.md (CRCL accepts raw Cockcroft-Gault when the source paper does not apply BSA normalization).",
      source_name        = "CLCR"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 37L,
    n_studies        = 1L,
    age_range        = "18+ years (adults; eligibility >=18 y)",
    age_median       = "61.7 years (mean, SD 17.9)",
    weight_range     = "approx 30-90 kg (mean 59.8 kg, SD 12.4)",
    weight_median    = "59.8 kg (mean, SD 12.4)",
    sex_female_pct   = 73,
    race_ethnicity   = "Korean (single-country study in the Republic of Korea)",
    disease_state    = "Adult inpatients with acute infections (pyelonephritis n=27, intra-abdominal infection n=9, neutropenic fever n=1) meeting ACCP/SCCM sepsis (n=30) or severe sepsis (n=7) criteria. Patients with septic shock, central-nervous-system infection, pneumonia, severe cardiovascular or hepatic disease, beta-lactam hypersensitivity, pregnancy, or renal replacement therapy were excluded.",
    dose_range       = "Doripenem 250 mg or 500 mg IV infusion over 1 hour, every 8 hours; 250 mg q8h for CLCR <=50 mL/min (n=9), 500 mg q8h for CLCR >50 mL/min (n=28). Blood sampling at predose and 0, 0.5, and 4-6 h after the fourth infusion.",
    regions          = "Republic of Korea (single center, Inje University Haeundae Paik Hospital, Busan)",
    renal_function   = "Cockcroft-Gault CLCR mean 66.7 mL/min (SD 34.4); 250-mg group 38.3 mL/min (SD 10.9); 500-mg group 75.9 mL/min (SD 34.5)",
    apache_ii        = "median 7 (range 0-15)",
    notes            = "Baseline demographics per Lee 2017 Table 1. 148 plasma concentrations (36 from 250-mg group, 112 from 500-mg group). Sparse-sampling design (3 post-dose samples per patient) supported a one-compartment fit even though doripenem PK is more typically described by a two-compartment model in dense-sampling studies (Methods/Discussion). Plasma drug assay validated LC-MS/MS with meropenem internal standard; LLOQ 0.2 ug/mL, linear range 0.2-50 ug/mL. NONMEM 7.3 FOCE-I; final model selected by chi-square OFV reduction (dOFV = -9.89 for CLCR on CL/WT)."
  )

  ini({
    # Structural parameters at the reference subject (WT=70 kg, CLCR=57 mL/min);
    # Lee 2017 Table 2. The paper parameterizes CL/WT and V/WT (per-kg), with
    # CL = 0.109 x WT x (CLCR/57)^0.688 and V = 0.280 x WT. Encoded here at
    # WT = 70 kg to align with the Caucasian comparison in Figure 1.
    lcl <- log(0.109 * 70); label("Clearance at WT=70 kg, CLCR=57 mL/min (L/h)") # Lee 2017 Table 2: CL/WT = 0.109 L/h/kg (RSE 8.57%); typical CL at 70 kg, CLCR=57 = 7.63 L/h
    lvc <- log(0.280 * 70); label("Central volume of distribution at WT=70 kg (L)")  # Lee 2017 Table 2: V/WT = 0.280 L/kg (RSE 9.60%); typical V at 70 kg = 19.6 L

    # Allometric scaling on body weight: the paper's per-kg parameterization
    # (CL/WT, V/WT) is mathematically equivalent to fixed linear allometric
    # exponents (1.0) on WT for both CL and V. Wrapped in fixed() because
    # this is a structural assumption, not an estimated parameter.
    e_wt_cl <- fixed(1); label("WT exponent on CL (fixed, structural)") # Lee 2017 Methods/Results: structural per-kg parameterization CL/WT = constant
    e_wt_vc <- fixed(1); label("WT exponent on V (fixed, structural)")  # Lee 2017 Methods/Results: structural per-kg parameterization V/WT = constant

    # CLCR power-law covariate effect on CL (Lee 2017 Table 2 / Results
    # equation): CL/WT = 0.109 x (CLCR/57)^0.688.
    e_crcl_cl <- 0.688; label("Power exponent on (CRCL/57) for CL") # Lee 2017 Table 2: theta_2 = 0.688 (RSE 22.9%)

    # Inter-individual variability (Lee 2017 Table 2 omega CV%);
    # omega^2 = log(CV^2 + 1) for log-normal etas.
    etalcl ~ 0.26433 # log(0.550^2 + 1); 55.0% CV on CL (RSE 14.4%)
    etalvc ~ 0.20177 # log(0.473^2 + 1); 47.3% CV on V  (RSE 21.6%)

    # Residual error: Lee 2017 Methods defines the "Poisson" error model as
    # Y = Ypred + eps * Ypred with eps ~ N(0, sigma^2), which is the standard
    # proportional residual-error form in nlmixr2 (~ prop()). Table 2 reports
    # sigma_Poisson = 0.633 (RSE 7.50%) directly as the SD.
    propSd <- 0.633; label("Proportional residual error (fraction)") # Lee 2017 Table 2: sigma_Poisson = 0.633 (RSE 7.50%)
  })
  model({
    # Individual PK parameters. CL scales linearly with WT and by power on
    # (CRCL/57); V scales linearly with WT. The per-kg form CL = 0.109 * WT
    # and V = 0.280 * WT is equivalent to the (WT/70)^1 scaling used here
    # with the reference TVCL = 7.63 L/h and TVV = 19.6 L at WT=70 kg.
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * (CRCL / 57)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose in mg, volume in L -> central/vc has units mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
