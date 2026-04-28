Takeuchi_2023_ozoralizumab <- function() {
  description <- "One-compartment population PK model with first-order absorption for subcutaneous ozoralizumab (anti-TNF VHH NANOBODY) in Japanese patients with rheumatoid arthritis (Takeuchi 2023)"
  reference <- "Takeuchi T, Chino Y, Mano Y, Kawanishi M, Sato Y, Uchida S, Tanaka Y. Population Pharmacokinetics of Ozoralizumab in Patients with Rheumatoid Arthritis. J Clin Pharmacol. 2024;64(4):418-427. doi:10.1002/jcph.2380"
  vignette <- "Takeuchi_2023_ozoralizumab"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL/F (exponent 0.847) and Vd/F (exponent 0.469) with reference (population median) weight 56.65 kg (Takeuchi 2023 final-model equations, p. 422).",
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Estimated glomerular filtration rate (Japanese formula)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL/F (exponent 0.191) with reference (population median) eGFR 85.95 mL/min/1.73 m^2 (Takeuchi 2023 final-model equation, p. 422). Calculated using the Japanese eGFR formula 194 * Scr^-1.094 * age^-0.287 (multiplied by 0.739 for women), per Methods (p. 420). Stored under the canonical CRCL; the source column name is eGFR.",
      source_name        = "eGFR"
    ),
    SEXF = list(
      description        = "Biological sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male) in the canonical column. The paper's own reference category is female (see notes).",
      notes              = "Takeuchi 2023 encodes sex as a male-indicator (Male = 1 for male, 0 for female) with female as the reference category for the published TVCL = 9.20 mL/h and TVVd = 4.91 L. To store under the canonical SEXF (1 = female, 0 = male) while preserving Takeuchi's female-reference TVCL/TVVd, the effect is applied in model() as (1 + e_male_cl * (1 - SEXF)) and (1 + e_male_vc * (1 - SEXF)), so SEXF = 1 yields factor 1 and SEXF = 0 yields the paper's male-vs-female fractional increase. The cohort is 76% female (Table 2), consistent with female being the typical-patient reference.",
      source_name        = "SEX"
    ),
    CONMED_MTX = list(
      description        = "Concomitant methotrexate use indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant methotrexate) in the canonical column. The paper's own reference category is with-MTX (see notes).",
      notes              = "Takeuchi 2023 encodes MTX use as MTX = 1 for yes, 0 for no, with concomitant-MTX as the reference category for the published TVCL = 9.20 mL/h. The OHZORA trial mandated MTX co-administration and the NATSUZORA trial excluded it, so MTX = 1 corresponds to OHZORA enrollment. To store under the canonical CONMED_MTX (1 = on MTX, 0 = not) while preserving Takeuchi's with-MTX-reference TVCL, the effect is applied in model() as (1 + e_nomtx_cl * (1 - CONMED_MTX)), so CONMED_MTX = 1 (on MTX) yields factor 1 and CONMED_MTX = 0 (no MTX) yields a 12.6% increase in CL/F.",
      source_name        = "MTX"
    ),
    ADA_POS = list(
      description        = "Antidrug-antibody status",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ADA-negative)",
      notes              = "Takeuchi 2023 final-model covariate (Table 3). Source column ADA (1 = positive, 0 = negative) renamed to canonical ADA_POS per covariate-columns.md. ADA status is defined per Methods: positive if antibody titer increased by >= 0.95 after dosing or the patient became positive after dosing; negative otherwise.",
      source_name        = "ADA"
    )
  )

  population <- list(
    n_subjects     = 494L,
    n_studies      = 2L,
    age_range      = "21 - 84 years (mean 56, SD 12; Table 2)",
    age_median     = "not reported (mean 56 years reported instead)",
    weight_range   = "35 - 112 kg (mean 59, SD 13; Table 2)",
    weight_median  = "56.65 kg (population median used as the reference weight in the final-model equations, p. 422)",
    sex_female_pct = 76,
    race_ethnicity = c(Japanese = 100),
    disease_state  = "Rheumatoid arthritis",
    dose_range     = "30 or 80 mg subcutaneously every 4 weeks for up to 52 weeks",
    regions        = "Japan",
    studies        = "OHZORA (Phase II/III with concomitant MTX, jRCT2080223971; n = 363) and NATSUZORA (Phase III without MTX, jRCT2080223973; n = 131); 3412 plasma concentrations after exclusion of placebo, criteria violations, BLQ/missing samples, and CWRES > 6 outliers.",
    renal_function = "Baseline eGFR (Japanese formula): mean 88, SD 20, range 35 - 174 mL/min/1.73 m^2; population median 85.95 (Table 2; reference value in final-model equation).",
    co_medication  = "Concomitant MTX in 363 of 494 patients (74%); 0 of 131 NATSUZORA patients on MTX, 363 of 363 OHZORA patients on MTX.",
    ada_status     = "ADA-positive 185 of 494 (37%; 33% in OHZORA, 49% in NATSUZORA per Table 2).",
    notes          = "Baseline demographics and study design per Takeuchi 2023 Tables 1 and 2 (p. 419-420)."
  )

  ini({
    # Structural parameters - reference (typical-patient) values for a 56.65 kg
    # female on concomitant MTX with eGFR 85.95 mL/min/1.73 m^2 and ADA-negative
    # status (Takeuchi 2023 Table 3, p. 421; final-model equations, p. 422).
    # Paper reports CL/F in mL/h; converted to L/h (/ 1000) for ODE consistency
    # with central in mg and Cc in mg/L = ug/mL.
    lka <- log(0.0343);        label("Absorption rate constant (Ka, 1/h)")                                            # Takeuchi 2023 Table 3: Ka = 0.0343 1/h
    lcl <- log(9.20 / 1000);   label("Apparent clearance for the typical patient (CL/F, L/h; paper reports 9.20 mL/h)") # Takeuchi 2023 Table 3: CL/F = 9.20 mL/h
    lvc <- log(4.91);          label("Apparent volume of distribution for the typical patient (Vd/F, L)")             # Takeuchi 2023 Table 3: Vd/F = 4.91 L

    # Covariate effects on CL/F (Takeuchi 2023 final-model equation, p. 422,
    # and Table 3, p. 421). Power form on continuous covariates; (1 + theta)
    # multiplicative form on binary categoricals. e_male_* and e_nomtx_cl are
    # applied in model() as (1 - SEXF) and (1 - CONMED_MTX) so the canonical-
    # column reference-category flips do not change the published coefficients.
    e_wt_cl     <- 0.847;  label("Power exponent of WT on CL/F (unitless)")                                            # Takeuchi 2023 Table 3: WT on CL/F = 0.847
    e_egfr_cl   <- 0.191;  label("Power exponent of CRCL (eGFR) on CL/F (unitless)")                                   # Takeuchi 2023 Table 3: eGFR on CL/F = 0.191
    e_nomtx_cl  <- 0.126;  label("Fractional change in CL/F when not on concomitant MTX (unitless)")                   # Takeuchi 2023 Table 3: MTX on CL/F = 0.126
    e_male_cl   <- 0.117;  label("Fractional change in CL/F for male sex relative to female (unitless)")               # Takeuchi 2023 Table 3: Sex on CL/F = 0.117
    e_ada_cl    <- 0.0967; label("Fractional change in CL/F for ADA-positive status (unitless)")                       # Takeuchi 2023 Table 3: ADA on CL/F = 0.0967

    # Covariate effects on Vd/F (Takeuchi 2023 final-model equation, p. 422,
    # and Table 3, p. 421).
    e_wt_vc     <- 0.469;  label("Power exponent of WT on Vd/F (unitless)")                                            # Takeuchi 2023 Table 3: WT on Vd/F = 0.469
    e_male_vc   <- 0.134;  label("Fractional change in Vd/F for male sex relative to female (unitless)")               # Takeuchi 2023 Table 3: Sex on Vd/F = 0.134

    # Inter-individual variability (log-normal etas, NONMEM exponential).
    # Block correlation between CL/F and Vd/F; eta on Ka was fixed to 0 in the
    # paper because of high shrinkage and is therefore omitted (etalka not in
    # model()). omega^2 values from Table 3.
    etalcl + etalvc ~ c(0.0423,
                        0.00806, 0.0230)                                                                              # Takeuchi 2023 Table 3: omega^2_CL/F = 0.0423, cov_CL:Vd = 0.00806, omega^2_Vd/F = 0.0230

    # Residual error (proportional only). Paper reports sigma^2 = 0.0347
    # (variance on the proportional scale); propSd is sqrt(sigma^2) ~= 0.186.
    propSd <- sqrt(0.0347); label("Proportional residual error (fraction; sqrt of paper's sigma^2 = 0.0347)")          # Takeuchi 2023 Table 3: sigma^2 = 0.0347
  })
  model({
    # Derived "is male" indicator: paper encoded sex as Male = 1 with female as
    # the reference for TVCL/TVVd. (1 - SEXF) reproduces the paper's Male column
    # while keeping SEXF (1 = female) as the canonical storage convention.
    sex_male <- 1 - SEXF
    # Derived "no MTX" indicator: paper encoded MTX = 1 with on-MTX as the
    # reference for TVCL. (1 - CONMED_MTX) reproduces the paper's (1 - MTX) term.
    no_mtx <- 1 - CONMED_MTX

    # Individual PK parameters (Takeuchi 2023 final-model equations, p. 422).
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) *
      (1 + e_ada_cl   * ADA_POS) *
      (CRCL / 85.95)^e_egfr_cl *
      (1 + e_nomtx_cl * no_mtx) *
      (1 + e_male_cl  * sex_male) *
      (WT / 56.65)^e_wt_cl
    vc <- exp(lvc + etalvc) *
      (1 + e_male_vc * sex_male) *
      (WT / 56.65)^e_wt_vc

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Dose in mg and vc in L -> central/vc has units mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
