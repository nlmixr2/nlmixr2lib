Li_2024_orin1001 <- function() {
  description <- "Two-compartment population PK model with first-order absorption and a lag time for oral ORIN1001 in Chinese patients with advanced solid tumors"
  reference <- paste(
    "Li X, Bo Y, Zeng Q, Diao L, Greene S, Patterson J, Liu L, Yang F (2024).",
    "Population pharmacokinetic model for oral ORIN1001 in Chinese patients",
    "with advanced solid tumors. Front Pharmacol 15:1322557.",
    "doi:10.3389/fphar.2024.1322557.",
    sep = " "
  )
  vignette <- "Li_2024_orin1001"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    TBILI = list(
      description        = "Total serum bilirubin",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline value; power scaling on CL/F normalized to the cohort median",
        "of 10.3 umol/L (Table 1, Eq. 3). Negative exponent: higher bilirubin",
        "lowers apparent clearance, interpreted by the authors as a marker of",
        "reduced hepatic metabolic / excretory capacity (ORIN1001 is mainly",
        "metabolized by the liver via aldehyde reduction). Already in SI units,",
        "so no mg/dL conversion is required.",
        sep = " "
      ),
      source_name        = "TBIL"
    ),
    LBM = list(
      description        = "Lean body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Derived from body weight, height, and sex via the Janmahasatian et al.",
        "(2005) formula, per Methods section 2.4.2. Power scaling on CL/F",
        "(Eq. 3) and on CL2/F (Eq. S1), both normalized to the cohort median of",
        "45.13 kg (Table 1). The paper's source column is LBW; LBM is the",
        "canonical register name for the same quantity.",
        sep = " "
      ),
      source_name        = "LBW"
    ),
    LDH = list(
      description        = "Serum lactate dehydrogenase",
      units              = "IU/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline value; power scaling on V2/F normalized to the cohort median",
        "of 214 IU/L (Table 1, Eq. S3, and the Figure 6 reference-patient",
        "definition). Retained as statistically significant but the authors",
        "found it had no clinically significant effect on steady-state",
        "exposure. The paper reports IU/L, used interchangeably with the",
        "register's canonical U/L.",
        sep = " "
      ),
      source_name        = "LDH"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte  = "ORIN1001",
      units    = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte  = "ORIN1001",
      units    = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte  = "ORIN1001",
      units    = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 25,
    n_studies      = 1,
    n_observations = 471,
    age_range      = "37-72 years",
    age_median     = "57 years",
    weight_range   = "42-80 kg",
    weight_median  = "64 kg",
    sex_female_pct = 48,
    race_ethnicity = c(Asian = 100),
    disease_state  = "advanced solid tumors",
    dose_range     = "100-900 mg orally once daily",
    regions        = "China",
    notes          = paste(
      "Phase I open-label, dose-escalation and dose-expansion basket trial",
      "(NCT05154201) at Peking University Cancer Hospital. Seven dose groups",
      "(100, 200, 300, 400, 500, 650, 900 mg QD). Each subject was followed for",
      "4 days after a single dose (single-dose period) and then through a",
      "21-day once-daily cycle (multiple-dose period). Baseline demographics and",
      "laboratory values are in Table 1; lean body weight median 45.13 kg",
      "(range 30-60.81), total bilirubin median 10.3 umol/L (range 4.80-23.20),",
      "and lactate dehydrogenase median 214 IU/L (range 98-1268).",
      sep = " "
    )
  )

  ini({
    # Structural parameters - final model, Table 4 "Final model / Estimate (%RSE)"
    lka <- log(0.58); label("Absorption rate constant (1/h)")  # Table 4, Ka = 0.58 (17.12% RSE)
    ltlag <- log(0.35); label("Absorption lag time (h)")  # Table 4, Tlag = 0.35 (19.63% RSE)
    lvc <- log(26.21); label("Apparent central volume of distribution, V/F (L)")  # Table 4, V/F = 26.21 (5.39% RSE); Eq. S2 confirms no covariate on V/F
    lvp <- log(26.60); label("Apparent peripheral volume of distribution, V2/F (L)")  # Table 4, V2/F = 26.60 (9.72% RSE)
    lcl <- log(1.07); label("Apparent clearance, CL/F (L/h)")  # Table 4, CL/F = 1.07 (5.31% RSE)
    lq <- log(0.75); label("Apparent intercompartmental clearance, CL2/F (L/h)")  # Table 4, CL2/F = 0.75 (12.57% RSE)

    # Covariate effects - power exponents on the cohort-median-normalized
    # covariate, per the general form in Eq. 2. The four retained effects are
    # CL/F-TBIL, CL/F-LBW, CL2/F-LBW and V2/F-LDH (Table 3 step 9 "Final model";
    # confirmed by the supplementary stepwise log, scenario cstep0323).
    e_tbili_cl <- -0.46; label("Exponent of total bilirubin on CL/F (unitless)")  # Table 4, "TBIL on CL/F" = -0.46 (-22.51% RSE); Eq. 3
    e_lbm_cl <- 1.11; label("Exponent of lean body weight on CL/F (unitless)")  # Table 4, "LBW on CL/F" = 1.11 (18.00% RSE); Eq. 3
    e_lbm_q <- 2.21; label("Exponent of lean body weight on CL2/F (unitless)")  # Table 4, "LBW on CL2/F" = 2.21 (19.30% RSE); Eq. S1
    e_ldh_vp <- 0.99; label("Exponent of lactate dehydrogenase on V2/F (unitless)")  # Table 4, "LDH on V2/F" = 0.99 (19.33% RSE); Eq. S3

    # IIV - exponential model (Eq. 1). Table 4 reports the variances directly
    # (footnote a: "omega V/F, variance of inter-individual variability for
    # V/F"), so these are used on the omega^2 scale without conversion.
    etalka ~ 0.534  # Table 4, omega^2 Ka = 0.534 (26.66% RSE)
    etaltlag ~ 0.438  # Table 4, omega^2 Tlag = 0.438 (19.33% RSE)
    etalvc ~ 0.049  # Table 4, omega^2 V/F = 0.049 (28.32% RSE)
    etalcl ~ 0.067  # Table 4, omega^2 CL/F = 0.067 (11.69% RSE)

    # Residual error - proportional (Results section 3.2: "Exponential model and
    # proportional model were applied to describe the inter- and
    # intra-individual variability, respectively"). Table 4 footnote a defines
    # stdev0 as a standard deviation, so it is used directly as propSd.
    propSd <- 0.197; label("Proportional residual error (fraction)")  # Table 4, stdev0 = 0.197 (7.48% RSE)
  })

  model({
    # 1. Individual parameters. Covariate effects are power functions of the
    #    covariate normalized by its cohort median (Eq. 2); the medians are the
    #    Table 1 values also quoted as the Figure 6 reference patient.
    ka <- exp(lka + etalka)
    tlag <- exp(ltlag + etaltlag)
    cl <- exp(lcl + etalcl) * (TBILI / 10.3)^e_tbili_cl * (LBM / 45.13)^e_lbm_cl  # Eq. 3
    vc <- exp(lvc + etalvc)  # Eq. S2
    q <- exp(lq) * (LBM / 45.13)^e_lbm_q  # Eq. S1
    vp <- exp(lvp) * (LDH / 214)^e_ldh_vp  # Eq. S3

    # 2. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 3. ODE system - first-order absorption from an oral depot with lag time,
    #    two-compartment distribution, first-order elimination (Results 3.2).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 4. Absorption lag time
    alag(depot) <- tlag

    # 5. Observation and error. Doses are in mg and volumes in L, so Cc is in
    #    mg/L = ug/mL; the paper tabulates concentrations in ng/mL (a factor of
    #    1000 larger).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
