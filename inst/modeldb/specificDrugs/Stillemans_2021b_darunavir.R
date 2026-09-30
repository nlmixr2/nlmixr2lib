Stillemans_2021b_darunavir <- function() {
  description <- "One-compartment population pharmacokinetic model with first-order absorption and elimination for oral ritonavir- or cobicistat-boosted darunavir in adults with HIV-1 infection, re-estimated on the merged learning and external-validation datasets (Stillemans 2021, Eur J Clin Pharmacol). Reduced covariate model without alpha-1 acid glycoprotein: apparent clearance is lower in women and in CYP3A5 nonexpressers (*3/*3); apparent volume is higher in SLCO3A1 rs8027174 G>T carriers."
  reference <- paste(
    "Stillemans G, Belkhir L, Vandercam B, Vincent A, Haufroid V, Elens L.",
    "Optimal sampling strategies for darunavir and external validation of",
    "the underlying population pharmacokinetic model. Eur J Clin Pharmacol.",
    "2021;77(4):607-616. doi:10.1007/s00228-020-03036-2. PMCID: PMC7935830.",
    sep = " "
  )
  vignette <- "Stillemans_2021b_darunavir"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Enters CL/F as exp(e_sexf_cl * SEXF). Table 2 reports 'theta sex on",
        "CL' = -0.151 (merged set); the exponential form is the one the",
        "underlying model uses (companion paper doi:10.1007/s40262-020-00920-z,",
        "ESM Online Resource 2 'Sex / CL / Exponential'). 33.1% female in the",
        "learning set and 36.0% in the validation set (Table 1)."
      ),
      source_name = "sex"
    ),
    CYP3A5_EXPR = list(
      description = "CYP3A5 expresser status: 1 = at least one CYP3A5*1 allele (*1/*1 or *1/*3), 0 = *3/*3 nonexpresser",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (nonexpresser, register convention); the published typical value is for expressers (CYP3A5_EXPR = 1)",
      notes = paste(
        "The source parameter is 'theta CYP3A5 on CL' = -0.126 (Table 2,",
        "merged set; the ESM Supplementary Material 1 labels the same row",
        "'CYP3A5*3 on CL'), a categorical effect P = theta * (1 + theta_cov)",
        "applied to *3/*3 nonexpressers (companion paper ESM Online",
        "Resource 2 'Categorical'). The indicator is therefore inverted",
        "relative to the canonical column: CL/F is multiplied by",
        "(1 + e_cyp3a5_cl * (1 - CYP3A5_EXPR)). Missing genotypes (2.4% of",
        "the learning set, 11.0% of the validation set) were replaced by the",
        "most frequent genotype for the patient's race (Methods, Validation",
        "dataset)."
      ),
      source_name = "CYP3A5 g.6986A>G (CYP3A5*3)"
    ),
    SNP_SLCO3A1_RS8027174 = list(
      description = "SLCO3A1 rs8027174 (g.91941607G>T) variant indicator: 1 = at least one T allele, 0 = G/G",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (G/G)",
      notes = paste(
        "Enters V/F as (1 + 0.697 * SNP_SLCO3A1_RS8027174) (Table 2 'theta",
        "SLCO3A1 on V' = 0.697, merged set; categorical form per the",
        "companion paper ESM Online Resource 2). No T/T homozygotes were",
        "observed in either dataset (Table 1: G/T 14.2% learning, 11.6%",
        "validation), so the indicator is heterozygote-vs-G/G in this cohort."
      ),
      source_name = "SLCO3A1 g.91941607G>T"
    )
  )

  covariatesDataExcluded <- list(
    AAG = list(
      description = "Plasma alpha-1 acid glycoprotein concentration",
      units = "g/L",
      type = "continuous",
      notes = paste(
        "Retained on CL/F and V/F in the full model of the companion paper",
        "(doi:10.1007/s40262-020-00920-z) but REMOVED from this model because",
        "AAG was not measured in the validation set (Results, External",
        "validation; removal raised the learning-set OFV from 617.307 to",
        "660.683, ESM Supplementary Material 1)."
      )
    )
  )

  compartmentData <- list(
    depot = list(analyte = "darunavir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "darunavir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 291,
    n_studies = 2,
    n_observations = 584,
    age_median = "55 years (IQR 13) learning set; 48 years (IQR 14) validation set (Table 1)",
    sex_female_pct = 34.7,
    race_ethnicity = "Learning set Caucasian 52.8%, African 43.3%; validation set Caucasian 61.6%, African 36.6% (Table 1)",
    disease_state = "HIV-1 infection, adult outpatients on boosted darunavir",
    dose_range = "darunavir 300 mg q12h to 1200 mg q24h; 800 mg q24h in 91.3% (learning) and 67.1% (validation), 600 mg q12h in 7.9% and 29.3% (Table 1)",
    regions = "Belgium (Cliniques universitaires Saint-Luc, Brussels)",
    booster = "cobicistat 85.8% / ritonavir 14.2% in the learning set; ritonavir 100% in the validation set (Table 1)",
    genotypes = "CYP3A5 *3/*3 45.7% (learning) and 47.0% (validation); SLCO3A1 rs8027174 G/T 14.2% and 11.6% (Table 1)",
    notes = paste(
      "Merged learning set (127 patients, 405 samples, including twelve",
      "6-h rich profiles; the companion Clin Pharmacokinet 2021 study) and",
      "external validation set (164 patients, 180 single random-time",
      "samples collected 2012-2016, one excluded for unknown post-intake",
      "time). 51 patients contributed to both studies; their data were",
      "treated as separate occasions (Methods, External validation), so",
      "n_subjects counts 291 patient-study records from about 240",
      "individuals. sex_female_pct is the pooled Table 1 count",
      "(42 + 59) / (127 + 164)."
    )
  )

  ini({
    # Table 2 'Comparison of model parameters with learning, merged, and
    # bootstrapped data', 'Merged set' column: the reduced model (no AAG)
    # re-estimated in NONMEM on the merged learning + validation data.
    lka <- log(0.724)
    label("Absorption rate constant ka (1/h)") # Table 2, 'k a (h -1)' merged set = 0.724
    lcl <- log(12.4)
    label("Apparent clearance CL/F for a male CYP3A5 expresser (L/h)") # Table 2, 'CL/F (l h -1)' merged set = 12.4
    lvc <- log(147)
    label("Apparent volume V/F for an SLCO3A1 G/G subject (L)") # Table 2, 'V/F (l)' merged set = 147

    e_sexf_cl <- -0.151
    label("Exponential female-sex coefficient on CL/F (unitless)") # Table 2, 'theta sex on CL' merged set = -0.151
    e_cyp3a5_cl <- -0.126
    label("Fractional change in CL/F for CYP3A5 nonexpressers (*3/*3) vs expressers (unitless)") # Table 2, 'theta CYP3A5 on CL' merged set = -0.126
    e_snp_slco3a1_rs8027174_vc <- 0.697
    label("Fractional change in V/F for SLCO3A1 rs8027174 T carriers vs G/G (unitless)") # Table 2, 'theta SLCO3A1 on V' merged set = 0.697

    # IIV: exponential random effects; Table 2 reports the SD, entered here
    # as the variance (SD^2). Diagonal, as in the underlying model.
    etalka ~ 0.506944 # Table 2, 'omega ka (sd)' merged set = 0.712 -> 0.712^2
    etalcl ~ 0.061504 # Table 2, 'omega CL (sd)' merged set = 0.248 -> 0.248^2
    etalvc ~ 0.094864 # Table 2, 'omega V (sd)' merged set = 0.308 -> 0.308^2

    # Residual error: exponential plus additive, encoded as proportional
    # plus additive.
    propSd <- 0.334
    label("Exponential (proportional) residual SD (fraction)") # Table 2, 'sigma exponential (sd)' merged set = 0.334
    addSd <- 0.539
    label("Additive residual SD (mg/L)") # Table 2, 'sigma additive (sd)' merged set = 0.539
  })

  model({
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) *
      exp(e_sexf_cl * SEXF) *
      (1 + e_cyp3a5_cl * (1 - CYP3A5_EXPR))
    vc <- exp(lvc + etalvc) *
      (1 + e_snp_slco3a1_rs8027174_vc * SNP_SLCO3A1_RS8027174)

    kel <- cl / vc

    # One-compartment model with first-order absorption and first-order
    # elimination (Methods, Population PK model).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
