Ruehs_2021_vericiguat <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption for oral",
    "vericiguat (a soluble guanylate cyclase stimulator) in adults with worsening",
    "chronic heart failure and left ventricular ejection fraction < 45% from the",
    "phase II SOCRATES-REDUCED study (final covariate model). Apparent clearance",
    "is scaled by body weight (fixed exponent 0.75), age, total bilirubin and",
    "creatinine clearance standardized to 70 kg; apparent volume by body weight",
    "(fixed exponent 1) and sex; and the absorption rate constant by body weight",
    "and serum albumin. Relative bioavailability falls stepwise with the",
    "administered dose level (1.08, 1, 0.867, 0.793 at 1.25, 2.5, 5, 10 mg).",
    "Residual error is combined proportional plus additive."
  )

  reference <- paste(
    "Ruehs H, Klein D, Frei M, Grevel J, Austin R, Becker C, Roessig L,",
    "Pieske B, Garmann D, Meyer M.",
    "Population Pharmacokinetics and Pharmacodynamics of Vericiguat in",
    "Patients with Heart Failure and Reduced Ejection Fraction.",
    "Clin Pharmacokinet. 2021;60(11):1407-1421.",
    "doi:10.1007/s40262-021-01024-y"
  )

  vignette <- "Ruehs_2021_vericiguat"

  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying in the source analysis. Allometric scaling of CL/F",
        "(exponent fixed 0.75) and V/F (exponent fixed 1), and an estimated",
        "power effect on ka (1.28), all normalised to 70 kg (Ruehs 2021",
        "Eqs. 2-3 and the Table 2 footnote 'Parameter covariate relations')."
      ),
      source_name = "WGHT"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on CL/F normalised to 68 years, the population median",
        "(Ruehs 2021 Table 2 footnote). The most influential covariate:",
        "adding it reduced the IIV of CL/F by 14.7% (ESM Table 3)."
      ),
      source_name = "AGE"
    ),
    TBILI = list(
      description = "Total serum bilirubin",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying in the source analysis. The paper normalises to",
        "0.6 (the population median; units not printed). The 0.6 reference",
        "is a mg/dL value (0.6 mg/dL = 10.3 umol/L is a typical median total",
        "bilirubin, 0.6 umol/L is not physiological), so the canonical SI",
        "value (umol/L) is converted in model() via tbili_mgdl = TBILI / 17.1."
      ),
      source_name = "BILI"
    ),
    CRCL = list(
      description = paste(
        "Cockcroft-Gault creatinine clearance standardized to a body weight",
        "of 70 kg (NOT BSA-normalized)"
      ),
      units = "mL/min/70 kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying in the source analysis. 'Standardized CLCR' is the",
        "Cockcroft-Gault creatinine clearance [Ruehs 2021 ref 27] standardized",
        "to 70 kg body weight [ref 28, Mould 2002], i.e.",
        "CRCL = CLcr_CG * 70 / WT -- the per-70-kg-body-weight normalisation",
        "(same convention as Suzuki_2024_mycophenolic_acid.R), not the",
        "canonical BSA normalisation. Power effect on CL/F normalised to",
        "100 mL/min (Table 2 footnote). Supplying a raw or BSA-normalized",
        "value would double-count body size against the allometric WT term."
      ),
      source_name = "CRCLST"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying in the source analysis. The paper normalises to 4.0",
        "(the population median; units not printed), a g/dL value; the",
        "canonical SI value (g/L) is converted in model() via",
        "alb_gdl = ALB / 10. Power effect on ka (exponent 2.37)."
      ),
      source_name = "ALB"
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Table 2 footnote: V/F multiplied by SEX = 1 for males and 0.850 for",
        "females, encoded as e_sexf_vc^SEXF."
      ),
      source_name = "SEX"
    ),
    DOSE = list(
      description = "Administered vericiguat dose level at the dose record",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Use case (a) of the DOSE canonical: the per-record administered dose",
        "selects the dose-level relative bioavailability of Ruehs 2021",
        "Table 2 (F = 1.08 for dose <= 1.25 mg, 1 for 2.5 mg [reference],",
        "0.867 for 5 mg, 0.793 for 10 mg). Encoded as a step function of",
        "DOSE with the upper bounds 1.25, 2.5 and 5 mg; doses above 10 mg",
        "(not studied) take the 10 mg value. In SOCRATES-REDUCED the dose was",
        "up-titrated 2.5 -> 5 -> 10 mg at weeks 2 and 4, so DOSE is",
        "time-varying and must be set on every dose record. Place the DOSE",
        "column AFTER the event columns (id, time, evid, amt, cmt) in the",
        "event table passed to rxode2."
      ),
      source_name = "DOSE"
    )
  )

  covariatesDataExcluded <- list(
    RACE = list(
      description = "Race",
      units = "(categorical)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Race on ka was significant at forward inclusion (p < 0.01) but was",
        "eliminated at backward deletion (p < 0.001); race on V/F was also",
        "tested and eliminated (Ruehs 2021 Section 3.3, ESM Table 3 runs",
        "8-9). Not in the final model (Table 2)."
      ),
      source_name = "RACE"
    ),
    NTPROBNP = list(
      description = "Baseline NT-proBNP",
      units = "pg/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Tested on the PK parameters and not significant (Ruehs 2021",
        "Section 3.3), together with atrial fibrillation, NYHA class and",
        "concomitant CYP3A4, P-gp, UGT1A9, BCRP and UGT1A1 inhibitors."
      ),
      source_name = "NT-proBNP"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "vericiguat", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "vericiguat", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 454,
    n_studies = 1,
    age_median = "68 years (covariate reference value)",
    weight_median = "70 kg (covariate reference value)",
    disease_state = paste(
      "Worsening chronic heart failure with left ventricular ejection",
      "fraction < 45%, on guideline-directed standard of care"
    ),
    dose_range = paste(
      "Oral vericiguat once daily for 12 weeks, target doses 1.25, 2.5, 5",
      "and 10 mg (the 5 and 10 mg arms started at 2.5 mg and were",
      "up-titrated at weeks 2 and 4), or placebo"
    ),
    regions = "Multinational (SOCRATES-REDUCED, NCT01951625)",
    notes = paste(
      "Ruehs 2021 Section 3.1 and Fig. 1: 456 randomized patients, 454 with",
      "PK samples; 363 received vericiguat and contributed 3376 eligible",
      "plasma samples. Covariate medians from the Table 2 footnote: age",
      "68 years, weight 70 kg, bilirubin 0.6 mg/dL, standardized CLcr",
      "100 mL/min, albumin 4.0 g/dL. The paper does not print a demographics",
      "table (sex, race and ranges are not reported)."
    )
  )

  ini({
    # Final covariate model, Ruehs 2021 Table 2 (reference patient: male,
    # 70 kg, 68 years, bilirubin 0.6 mg/dL, standardized CLcr 100 mL/min,
    # albumin 4.0 g/dL, 2.5 mg dose so that F = 1).
    lka <- log(1.29); label("Absorption rate constant ka (1/h)") # Table 2: ka = 1.29 1/h (RSE 9.54%)
    lcl <- log(1.24); label("Apparent clearance CL/F (L/h)") # Table 2: CL/F = 1.24 L/h (RSE 3.19%)
    lvc <- log(34.3); label("Apparent volume of distribution V/F (L)") # Table 2: V/F = 34.3 L (RSE 2.12%); the footnote equation's '3.43' is a typo

    e_age_cl <- -0.418; label("Power exponent of age on CL/F (unitless)") # Table 2: theta CL,age = -0.418 (RSE 20.5%)
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL/F (unitless)") # Table 2: theta CL,bodyweight = 0.75, footnote b 'Fixed value'; Eq. 2
    e_crcl_cl <- 0.164; label("Power exponent of standardized creatinine clearance on CL/F (unitless)") # Table 2: theta CL,standardized creatinine clearance = 0.164 (RSE 26.1%)
    e_tbili_cl <- -0.072; label("Power exponent of total bilirubin on CL/F (unitless)") # Table 2: theta CL,bilirubin = -0.072 (RSE 26.1%); footnote equation prints -0.075
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on V/F (unitless)") # Table 2: theta V,bodyweight = 1.00, footnote b 'Fixed value'; Eq. 3
    e_sexf_vc <- 0.850; label("Multiplicative factor on V/F for females (unitless)") # Table 2: theta V,sex = 0.850 (RSE 4.09%); footnote SEX = 1 male, 0.850 female
    e_wt_ka <- 1.28; label("Power exponent of body weight on ka (unitless)") # Table 2: theta ka,bodyweight = 1.28 (RSE 25.2%)
    e_alb_ka <- 2.37; label("Power exponent of serum albumin on ka (unitless)") # Table 2: theta ka,albumin = 2.37 (RSE 26.1%)

    e_dose_1p25mg_fdepot <- 1.08; label("Relative bioavailability at doses <= 1.25 mg vs 2.5 mg (unitless)") # Table 2: theta F,dose<=1.25 mg = 1.08 (RSE 2.62%)
    e_dose_5mg_fdepot <- 0.867; label("Relative bioavailability at 5 mg vs 2.5 mg (unitless)") # Table 2: theta F,dose=5 mg = 0.867 (RSE 1.94%)
    e_dose_10mg_fdepot <- 0.793; label("Relative bioavailability at 10 mg vs 2.5 mg (unitless)") # Table 2: theta F,dose=10 mg = 0.793 (RSE 2.60%)

    # IIV: omega^2 values as printed for the final run (ESM Table 3 run 10:
    # ka 0.867, CL/F 0.061, V/F 0.043), which reproduce the Table 2 CVs via
    # CV = sqrt(exp(omega^2) - 1): 117%, 25.1%, 21.0%. The CL/F-V/F
    # correlation of the base model was removed before the covariate
    # analysis (Section 3.3), so the final IIV is diagonal.
    etalka ~ 0.867 # ESM Table 3 run 10 IIV of ka = 0.867; Table 2 CV 117%
    etalcl ~ 0.061 # ESM Table 3 run 10 IIV of CL/F = 0.061; Table 2 CV 25.1%
    etalvc ~ 0.043 # ESM Table 3 run 10 IIV of V/F = 0.043; Table 2 CV 21.0%

    propSd <- 0.257; label("Proportional residual error (fraction)") # Table 2: proportional error 25.7% = SQRT(SIGMA^2) x 100 (RSE 4.21%)
    addSd <- 7.27; label("Additive residual error (ug/L)") # Table 2: additive error 7.27 ug/L = SQRT(SIGMA^2) (RSE 14.9%)
  })

  model({
    # Unit conversions: canonical SI covariates -> the US-convention units
    # the paper's reference values (0.6 mg/dL, 4.0 g/dL) are expressed in.
    tbili_mgdl <- TBILI / 17.1
    alb_gdl <- ALB / 10

    ka <- exp(lka + etalka) * (alb_gdl / 4.0)^e_alb_ka * (WT / 70)^e_wt_ka
    cl <- exp(lcl + etalcl) * (AGE / 68)^e_age_cl * (tbili_mgdl / 0.6)^e_tbili_cl *
      (CRCL / 100)^e_crcl_cl * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * e_sexf_vc^SEXF * (WT / 70)^e_wt_vc

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Dose-level relative bioavailability (Table 2 footnote), reference
    # 2.5 mg (F = 1). Doses above 10 mg were not studied and take the
    # 10 mg value.
    fdose <- e_dose_10mg_fdepot
    if (DOSE <= 5) fdose <- e_dose_5mg_fdepot
    if (DOSE <= 2.5) fdose <- 1
    if (DOSE <= 1.25) fdose <- e_dose_1p25mg_fdepot
    f(depot) <- fdose

    # Dose in mg, volume in L -> mg/L; x 1000 -> ug/L (the unit of the
    # additive residual error and of ESM Table 4).
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
