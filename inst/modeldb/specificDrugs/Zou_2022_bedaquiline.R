Zou_2022_bedaquiline <- function() {
  description <- "One-compartment population PK model with first-order absorption for oral bedaquiline in 99 Chinese patients with multidrug-resistant pulmonary tuberculosis sampled during the 200 mg three-times-weekly continuation phase (Zou 2022). Apparent clearance decreases with serum gamma-glutamyl transferase (GGT) as a power function centred on the dataset median of 28.9 U/L, and is lowered by an additive 1.4 L/h in patients homozygous GG at AGBL4 rs319952 (reference AA + AG). Typical values CL/F = 4.54 L/h, Vc/F = 227 L, Ka = 0.447 1/h; exponential IIV on CL/F and Vc/F and a proportional residual error."
  reference <- paste(
    "Zou J, Chen S, Rao W, Fu L, Zhang J, Liao Y, Zhang Y, Lv N, Deng G, Yang S, Lin L, Li L, Liu S, Qu J.",
    "Population pharmacokinetic modeling of bedaquiline among multidrug-resistant pulmonary",
    "tuberculosis patients from China.",
    "Antimicrob Agents Chemother. 2022;66(10):e00811-22.",
    "doi:10.1128/aac.00811-22.",
    sep = " "
  )
  vignette <- "Zou_2022_bedaquiline"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "bedaquiline", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "bedaquiline", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    GGT = list(
      description = "Serum gamma-glutamyl transferase activity",
      units = "U/L",
      type = "continuous",
      reference_category = "n/a (power model centred on the dataset median 28.9 U/L)",
      notes = paste(
        "Laboratory value measured within 2 days of the PK sample (Methods, Study subjects).",
        "Cohort mean 36.4 U/L (SD 31.1), range 8-271.8 U/L (Table 1); the centring value",
        "28.9 U/L is the median printed in Results Equations 1 and 2. Missing values were",
        "imputed with the median (Methods). Effect: CL/F scales by (GGT/28.9)^-0.476, a",
        "28% decrease per doubling of GGT. Combined with the additive rs319952 shift, typical",
        "CL/F becomes non-positive for GG subjects once GGT exceeds about 342 U/L, above",
        "the observed maximum; do not simulate GG subjects beyond the observed GGT range."
      ),
      source_name = "GGT"
    ),
    SNP_AGBL4_RS319952_GG = list(
      description = "AGBL4 rs319952 homozygous GG genotype indicator (1 = GG, 0 = AA or AG)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (genotypes AA and AG pooled)",
      notes = paste(
        "Time-fixed per subject (germline genotype). Table 2: AA 31, AG 39, GG 13, not genotyped",
        "16 of 99; the A allele is the reference allele. Missing categorical covariates were",
        "imputed with the most frequent class (Methods), which for rs319952 is AG, so the",
        "16 ungenotyped subjects entered the fit with this indicator = 0. Effect: an",
        "ADDITIVE shift of -1.4 L/h on typical CL/F (Results Equation 2), not a fractional one."
      ),
      source_name = "rs319952 genotype G (vs A&AG)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Mean 38.1 years (SD 14.1), range 11-78 (Table 1). Significant on CL/F in the graphical screen (P = 0.0185, Figure 1); not retained after forward inclusion."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Significant on CL/F in the graphical screen (P = 0.0272, Figure 1); not retained. The sex split of the cohort is not reported."
    ),
    TP = list(
      description = "Serum total protein",
      units = "g/L",
      type = "continuous",
      notes = "Mean 77.5 g/L (SD 5.8), range 65.4-94.6 (Table 1). Significant on Vc/F in forward inclusion (P < 0.01) but removed in backward elimination (P > 0.005)."
    ),
    SNP_AGBL4_RS320003 = list(
      description = "AGBL4 rs320003 genotype",
      units = "(binary)",
      type = "categorical",
      notes = "A 12, G (reference) 69, missing 18 (Table 2). Significant on CL/F in the graphical screen (AA vs GG, P = 0.0268); not retained."
    ),
    SNP_CYP2E1_RS2031920 = list(
      description = "CYP2E1 rs2031920 genotype",
      units = "(binary)",
      type = "categorical",
      notes = "C (reference) 46, T 3, TC 12, missing 38 (Table 2). Significant on CL/F in the graphical screen (CC vs TC+TT, P = 0.028); not retained."
    ),
    SNP_NOS2_RS11080344 = list(
      description = "NOS2 rs11080344 genotype",
      units = "(binary)",
      type = "categorical",
      notes = "C 29, T (reference) 18, TC 36, missing 16 (Table 2). Significant on CL/F in the graphical screen (CC vs TC+TT, P = 0.0205); not retained."
    ),
    SNP_CYP2C9_RS9332096 = list(
      description = "CYP2C9 rs9332096 genotype",
      units = "(binary)",
      type = "categorical",
      notes = "C (reference) 77, CT 6, missing 16 (Table 2). Significant on Vc/F in the graphical screen (CC vs CT, P = 0.00354); not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 99L,
    n_studies = 1L,
    n_observations = "246 plasma bedaquiline concentrations (1-8 per patient; 38 patients contributed a single sample). One sample below the 0.024 ug/mL limit was set to that limit.",
    age_range = "11-78 years",
    age_mean = "38.1 years (SD 14.1)",
    weight_range = "17-90.7 kg",
    weight_mean = "58 kg (SD 12.6)",
    height_mean = "166.4 cm (SD 12.4)",
    sex_female_pct = NA_real_,
    race_ethnicity = c(Asian = 100),
    disease_state = "Multidrug-resistant pulmonary tuberculosis.",
    dose_range = "Approved regimen: 400 mg orally once daily for 2 weeks, then 200 mg three times per week for 22 weeks. Samples were drawn randomly within a dosing interval during weeks 3-24 (200 mg three times weekly; median 20 weeks on therapy). Background drugs: cycloserine 92%, linezolid 97%, pyrazinamide 49%, moxifloxacin 41%, clofazimine 34%, levofloxacin 27%, ethambutol 1%.",
    regions = "China (single centre, Shenzhen Third People's Hospital).",
    notes = "Recruited October 2020 to October 2021. GGT mean 36.4 U/L (range 8-271.8); ALB mean 46.2 g/L (Table 1). AGBL4 rs319952 genotypes AA/AG/GG = 31/39/13, 16 not genotyped (Table 2). NONMEM 7.4, FOCE-I; 1,000-sample bootstrap (98.7% converged) and prediction-corrected VPC. The sex distribution is not reported."
  )

  ini({
    # Structural PK -- Zou 2022 Table 3 'Parameter estimates for the final model',
    # Parameter estimate column (RSE % in parentheses).
    lka <- log(0.447); label("Absorption rate constant Ka (1/h)")                              # Table 3: Ka = 0.447 1/h (RSE 16.6%)
    lcl <- log(4.54);  label("Apparent clearance CL/F for GGT = 28.9 U/L, rs319952 AA/AG (L/h)") # Table 3: CL/F = 4.54 L/h (RSE 5.3%); Results Eq. 1
    lvc <- log(227);   label("Apparent central volume of distribution Vc/F (L)")               # Table 3: Vc/F = 227 L (RSE 16.8%)

    # Covariate effects on CL/F -- Table 3 'Covariate parameter' rows and
    # Results Equations 1 and 2:
    #   AA/AG: CL/F = 4.54 * (GGT/28.9)^-0.476
    #   GG:    CL/F = 4.54 * (GGT/28.9)^-0.476 - 1.4
    e_ggt_cl <- -0.476; label("Power exponent of GGT/28.9 on CL/F (unitless)")                          # Table 3: theta GGT on CL/F = -0.476 (RSE 18.1%)
    e_snp_agbl4_rs319952_gg_cl <- -1.4; label("Additive shift of CL/F for rs319952 GG vs AA/AG (L/h)") # Table 3: theta rs319952 on CL/F = -1.4 (RSE 27.1%); Results: 'CL/F decreased by 1.4 L/h'

    # IIV -- Table 3 'Interindividual variability (%)', exponential model
    # (Methods Eq. 3). The table does not say whether the percentages are CVs
    # or 100 * sqrt(omega^2); read as CVs, omega^2 = log(1 + CV^2).
    etalcl ~ 0.1395 # Table 3: eta CL/F = 38.7% -> log(1 + 0.387^2) = 0.1395
    etalvc ~ 0.5290 # Table 3: eta Vc/F = 83.5% -> log(1 + 0.835^2) = 0.5290

    # Residual error -- Table 3 'Residual variability (%)', proportional model
    # (Results, Model building).
    propSd <- 0.322; label("Proportional residual error (fraction)") # Table 3: eps prop = 32.2% (RSE 6.7%)
  })

  model({
    # Results Eqs. 1-2: power GGT effect, then an additive genotype shift
    # (Table 4's AUCweekly = 600 mg / CL confirms the shift is applied after the
    # power term, not multiplicatively).
    tvcl <- exp(lcl) * (GGT / 28.9)^e_ggt_cl + e_snp_agbl4_rs319952_gg_cl * SNP_AGBL4_RS319952_GG
    cl <- tvcl * exp(etalcl)
    vc <- exp(lvc + etalvc)
    ka <- exp(lka)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
