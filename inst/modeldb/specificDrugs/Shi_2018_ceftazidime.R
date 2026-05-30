Shi_2018_ceftazidime <- function() {
  description <- "One-compartment IV population PK model for ceftazidime in infants 0.1-2.0 years (Shi 2018) with allometric body-weight scaling and a power-form creatinine-clearance effect on clearance."
  reference <- paste(
    "Shi ZR, Chen XK, Tian LY, Wang YK, Zhang GY, Dong L, Jirasomprasert T,",
    "Jacqz-Aigrain E, Zhao W. Population Pharmacokinetics and Dosing",
    "Optimization of Ceftazidime in Infants.",
    "Antimicrob Agents Chemother. 2018;62(4):e02486-17.",
    "doi:10.1128/AAC.02486-17.",
    sep = " "
  )
  vignette <- "Shi_2018_ceftazidime"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Current body weight at the time of pharmacokinetic sampling",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling on CL (fixed exponent 0.75) and V (fixed exponent 1.0) with reference 7.5 kg (cohort median). Shi 2018 Table 2 footnote a.",
      source_name        = "BW"
    ),
    CRCL = list(
      description        = "Creatinine clearance computed from serum creatinine collected within 48 h of pharmacokinetic sampling (Schwartz-style pediatric formula per Shi 2018 Methods reference 24)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column CRCL in raw mL/min (NOT BSA-normalized to mL/min/1.73 m^2). Stored under the canonical CRCL column per inst/references/covariate-columns.md (CRCL accepts raw mL/min when the source paper does not apply BSA normalization, with the per-model description recording the assay form -- same convention as Delattre_2010_amikacin). Reference value 124 mL/min (population median, Shi 2018 Table 1 / Table 2 footnote a). Effect on CL is a power form: (CRCL / 124)^theta_CRCL.",
      source_name        = "CRCL"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 51L,
    n_studies      = 1L,
    age_range      = "0.1-2.0 years",
    age_median     = "0.5 years",
    weight_range   = "3.0-13.0 kg",
    weight_median  = "7.5 kg",
    sex_female_pct = 45,
    race_ethnicity = "Not reported (single-centre Chinese cohort; Children's Hospital of Hebei Province)",
    disease_state  = "Infants 1 month to 2 years with confirmed or suspected bacterial infection (pneumonia, bronchitis, and other miscellaneous pulmonary infections); preterm newborns (gestational age < 37 weeks) excluded",
    dose_range     = "50 mg/kg ceftazidime IV every 12 h administered as a 30-60 min infusion (uniform dose regimen across cohort)",
    regions        = "China (Children's Hospital of Hebei Province, Shijiazhuang)",
    renal_function = "Creatinine clearance median 124 mL/min (range 65-181, raw mL/min); serum creatinine median 24 umol/L (range 16-44)",
    notes          = "Baseline demographics per Shi 2018 Table 1. 51 infants enrolled September 2015 to December 2016. 90 pharmacokinetic samples available (sparse design, 2 per subject at steady state). LLOQ 0.5 ug/mL; 12 BLQ values replaced with 0.25 ug/mL (half-LLOQ) per Shi 2018 Results. Sampling schedule: randomized to either (1) 4-8 h after infusion start and 3-5 min after infusion end, or (2) 8-12 h after infusion start and 1-2 h after infusion end."
  )

  ini({
    # Structural parameters at the reference subject (WT = 7.5 kg, CRCL = 124 mL/min);
    # Shi 2018 Table 2 final-model column.
    lcl <- log(1.30); label("Clearance at 7.5 kg / 124 mL/min CRCL (L/h)") # Shi 2018 Table 2: theta_1 (CL) = 1.30 L/h
    lvc <- log(2.97); label("Volume of distribution at 7.5 kg (L)")        # Shi 2018 Table 2: theta_2 (V) = 2.97 L

    # Allometric exponents (fixed at theoretical values per Shi 2018 Covariate analysis;
    # "this model with fixed allometric coefficients was more fit than that with unfixed
    # coefficients, which caused the significant drop in the objective function value of 23 points").
    allo_cl <- fixed(0.75); label("Allometric exponent on CL (unitless, fixed)") # Shi 2018 Covariate analysis section
    allo_vc <- fixed(1.0);  label("Allometric exponent on V (unitless, fixed)")  # Shi 2018 Covariate analysis section

    # Power-form CRCL effect on CL: F_CRCL = (CRCL / 124)^theta_3.
    e_crcl_cl <- 0.82; label("CRCL exponent on CL (unitless)") # Shi 2018 Table 2: theta_3 = 0.82

    # Inter-individual variability (Shi 2018 Table 2 final-model CV%);
    # omega^2 = log(CV^2 + 1) for log-normal etas.
    etalcl ~ 0.02851 # log(0.170^2 + 1); 17.0% CV on CL (Shi 2018 Table 2)
    etalvc ~ 0.01575 # log(0.126^2 + 1); 12.6% CV on V  (Shi 2018 Table 2)

    # Combined residual error (Shi 2018 Table 2 final-model rows ERR(1) and ERR(2)).
    # ERR(1) = 38.2 is interpreted as the proportional CV (fraction 0.382) per the
    # standard NONMEM $SIGMA output convention for a combined additive + proportional
    # error model (see paper Methods: "A combined additive and proportional model best
    # described residual variability"). ERR(2) = 16.0 is interpreted as the additive
    # SD in concentration units (mg/L equivalent to ug/mL) following the same
    # convention used in Delattre_2010_amikacin and other antibiotic-PK extractions
    # in this library; the table's "Residual variability (%)" header label applies
    # cleanly only to the proportional component. See the validation vignette's
    # Assumptions and deviations section for the documented unit interpretation.
    propSd <- 0.382; label("Proportional residual error (fraction)") # Shi 2018 Table 2 ERR(1) = 38.2%
    addSd  <- 16.0;  label("Additive residual error (mg/L)")         # Shi 2018 Table 2 ERR(2) = 16.0
  })
  model({
    # Individual PK parameters. CL has a body-weight allometric scaling with fixed
    # exponent 0.75 and a power-form CRCL effect; V has body-weight allometric
    # scaling with fixed exponent 1.0.
    cl <- exp(lcl + etalcl) * (WT / 7.5)^allo_cl * (CRCL / 124)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (WT / 7.5)^allo_vc

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose in mg, vc in L -> central/vc has units mg/L (equivalent to ug/mL).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
