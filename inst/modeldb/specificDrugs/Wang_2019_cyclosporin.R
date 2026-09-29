Wang_2019_cyclosporin <- function() {
  description <- "One-compartment first-order-absorption population PK model for oral cyclosporin in Chinese children with pediatric refractory nephrotic syndrome, with allometric body weight on CL/F and V/F and a concomitant-spironolactone effect on CL/F (Wang 2019)"
  reference <- "Wang DD, Chen X, Li ZP. Cyclosporin population pharmacokinetics in pediatric refractory nephrotic syndrome based on real-world studies: Effects of body weight and spirolactone administration. Exp Ther Med. 2019;17(4):3015-3020. doi:10.3892/etm.2019.7325"
  vignette <- "Wang_2019_cyclosporin"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Cyclosporin whole-blood concentrations were measured
  # with the Emit 2000 Cyclosporine Specific assay (Wang 2019 Methods,
  # 'Analytical method'), so the modelled matrix is whole blood.
  compartmentData <- list(
    depot = list(analyte = "cyclosporin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "cyclosporin", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric power scaling on both CL/F and V/F, normalised to a",
        "standard 70 kg adult (Wang 2019 Methods 'Covariate model':",
        "Pi = Pstd * (WTi / WTstd)^PWR with WTstd = 70 kg, PWR = 0.75 for",
        "CL/F and 1 for V/F, following Anderson and Holford, reference 15).",
        "Cohort body weight median 15 kg, range 10-23 kg (Wang 2019",
        "Table I), so every subject is extrapolated well below the 70 kg",
        "reference. The paper does not describe a time-varying weight",
        "record; the retrospective TDM dataset spans June 2014 to June 2018."
      ),
      source_name = "WT"
    ),
    CONMED_SPIRON = list(
      description = "Concomitant spironolactone coadministration indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant spironolactone)",
      notes = paste(
        "11 of 18 patients received concomitant spironolactone (Wang 2019",
        "Table II; the paper spells the drug 'spirolactone' throughout).",
        "Time-fixed per subject in the published analysis. Linear-deviation",
        "multiplicative effect on CL/F:",
        "cl *= (1 + CONMED_SPIRON * e_conmed_spiron_cl) with",
        "e_conmed_spiron_cl = -0.265 (Wang 2019 Results equation and",
        "Table III), i.e. 26.5% lower apparent oral clearance when",
        "spironolactone is coadministered."
      ),
      source_name = "spirolactone"
    )
  )

  # Covariates that Wang 2019 screened in the stepwise covariate search but
  # did not retain in the final model. The Methods ('Covariate model') list
  # every demographic, laboratory and concomitant-medication term that was
  # tested; only body weight and concomitant spironolactone survived the
  # forward-inclusion (dOFV > 3.84) and backward-elimination (dOFV > 6.64)
  # criteria. These entries carry the provenance of that screen without
  # implying a retained effect. Baseline summaries are in Wang 2019 Table I;
  # concomitant-medication counts are in Table II.
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "13 male / 5 female (Wang 2019 Table I); screened, not retained."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Mean 2.79 +/- 0.90, median 2.75 (range 1.18-4.57) years (Wang 2019 Table I); screened, not retained."
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      notes = "Mean 91.28 +/- 8.23, median 94 (range 77-105) cm (Wang 2019 Table I); screened, not retained."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Mean 19.68 +/- 6.82, median 20.1 (range 10.1-33.2) g/L (Wang 2019 Table I); screened, not retained. Markedly hypoalbuminaemic, as expected in nephrotic syndrome."
    ),
    TPRO = list(
      description = "Total serum protein",
      units = "g/L",
      type = "continuous",
      notes = "Mean 45.16 +/- 7.41, median 43.55 (range 33.9-66.0) g/L (Wang 2019 Table I); screened, not retained."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Mean 12.12 +/- 12.71, median 9 (range 1-79) IU/L (Wang 2019 Table I); screened, not retained."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Mean 26.35 +/- 15.84, median 20 (range 11-75) IU/L (Wang 2019 Table I); screened, not retained."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Mean 25.47 +/- 23.72, median 20.5 (range 7-152) umol/L (Wang 2019 Table I); screened, not retained. Patients with diagnosed kidney failure were excluded (Wang 2019 Fig. 1)."
    ),
    BUN = list(
      description = "Blood urea",
      units = "mmol/L",
      type = "continuous",
      notes = "Mean 4.19 +/- 2.05, median 3.65 (range 1.4-11.5); screened, not retained. Wang 2019 Table I prints the unit as umol/l, which is three orders of magnitude away from any physiological blood-urea concentration; the values are consistent with mmol/L and are recorded as such here."
    ),
    TBA = list(
      description = "Total serum bile acids",
      units = "umol/L",
      type = "continuous",
      notes = "Mean 4.79 +/- 3.85, median 3.6 (range 0.3-15.6) umol/L (Wang 2019 Table I); screened, not retained."
    ),
    DBIL = list(
      description = "Direct (conjugated) bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Mean 0.81 +/- 0.83, median 0.6 (range 0.1-4.8) umol/L (Wang 2019 Table I); screened, not retained."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Mean 4.24 +/- 2.54, median 3.5 (range 1.3-12.1) umol/L (Wang 2019 Table I); screened, not retained."
    ),
    HCT = list(
      description = "Hematocrit",
      units = "%",
      type = "continuous",
      notes = "Mean 39.64 +/- 5.17, median 39 (range 30.22-51.9) % (Wang 2019 Table I); screened, not retained. Cyclosporin partitions extensively into erythrocytes, so hematocrit is a mechanistically plausible covariate that this 18-patient dataset was not powered to resolve."
    ),
    HGB = list(
      description = "Hemoglobin",
      units = "g/L",
      type = "continuous",
      notes = "Mean 130.34 +/- 17.53, median 127.5 (range 97-173.1) g/L (Wang 2019 Table I); screened, not retained."
    ),
    CONMED_STEROID = list(
      description = "Concomitant systemic corticosteroid administration indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant systemic corticosteroid)",
      notes = "Prednisolone in 13 of 18 and methylprednisolone in 4 of 18 patients (Wang 2019 Table II); both were screened individually and neither was retained. Pooled here into the canonical systemic-corticosteroid indicator for documentation only."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 18L,
    n_studies = 1L,
    age_range = "1.18-4.57 years",
    age_median = "2.75 years (mean +/- SD 2.79 +/- 0.90)",
    weight_range = "10-23 kg",
    weight_median = "15 kg (mean +/- SD 15.28 +/- 2.95)",
    sex_female_pct = 27.8,
    race_ethnicity = "Chinese",
    disease_state = "Chinese children (<16 years) with pediatric refractory nephrotic syndrome (PRNS); patients with any other serious disease, including kidney failure, were excluded",
    dose_range = "25-80 mg daily orally (liquid solution), subsequently adjusted on clinical response, adverse events and trough TDM concentration",
    regions = "China (Children's Hospital of Fudan University, Shanghai)",
    notes = paste(
      "Retrospective real-world therapeutic-drug-monitoring (TDM) dataset,",
      "June 2014 to June 2018. All modelled cyclosporin whole-blood",
      "concentrations were trough concentrations measured with the Emit 2000",
      "Cyclosporine Specific assay; the paper does not report the number of",
      "concentration records. Because only troughs were available, neither",
      "bioavailability nor an absorption lag time was estimable and the",
      "absorption rate constant was fixed to a literature value, so CL/F and",
      "V/F are apparent (F-scaled) parameters (Wang 2019 Methods, 'PPK",
      "modeling'). Covariates screened but not retained: sex, age, body",
      "height, albumin, globulin, albumin/globulin ratio, alanine",
      "transaminase, aspartate transaminase, creatinine, urea, total protein,",
      "total bile acid, direct bilirubin, total bilirubin, hematocrit,",
      "hemoglobin, mean corpuscular hemoglobin, mean corpuscular hemoglobin",
      "concentration, and concomitant diltiazem (3/18), dipyridamole (2/18),",
      "felodipine (1/18), fosinopril (3/18), methylprednisolone (4/18),",
      "nifedipine (1/18), piperazine ferulate (5/18) and prednisolone (13/18)",
      "(Wang 2019 Methods 'Covariate model'; Tables I and II). Globulin,",
      "albumin/globulin ratio, mean corpuscular hemoglobin and mean",
      "corpuscular hemoglobin concentration have no canonical covariate",
      "column and are recorded here in prose only."
    )
  )

  ini({
    # Structural parameters; allometric reference body weight 70 kg
    lka <- fixed(log(0.68))
    label("Absorption rate constant (Ka, 1/h); from the Ni 2013 pediatric ciclosporin model, Wang 2019 reference 9 (Wang 2019 Methods 'PPK modeling'; Table III)") # Table III 'Ka (h-1) 0.68 (fixed)'
    lcl <- log(80.7)
    label("Apparent oral clearance at WT = 70 kg without concomitant spironolactone (CL/F, L/h); Wang 2019 Table III") # Table III 'CL/F (L/h) 80.7'
    lvc <- log(2030)
    label("Apparent volume of distribution at WT = 70 kg (V/F, L); Wang 2019 Table III") # Table III 'V/F (L) 2030'

    # Allometric exponents fixed by convention (Anderson and Holford; Wang
    # 2019 Methods 'Covariate model', reference 15)
    e_wt_cl <- fixed(0.75)
    label("Allometric (WT/70) exponent on CL/F (unitless); Wang 2019 Methods 'Covariate model'") # Methods 'Covariate model': PWR = 0.75 for CL/F
    e_wt_vc <- fixed(1.0)
    label("Allometric (WT/70) exponent on V/F (unitless); Wang 2019 Methods 'Covariate model'") # Methods 'Covariate model': PWR = 1 for V/F

    # Covariate effect on CL/F
    e_conmed_spiron_cl <- -0.265
    label("Fractional change in CL/F with concomitant spironolactone (unitless); Wang 2019 Table III") # Table III 'theta-spirolactone -0.265'

    # Inter-individual variability. Wang 2019 Table III prints omega on the
    # standard-deviation scale: the Abstract and Discussion quote the same
    # numbers as 'the inter-individual variability in CL/F and V/F was 44.6
    # and 53.1%', i.e. CV% = 100 * omega under the small-omega approximation.
    # Reading them as variances would give CV = sqrt(0.446) = 66.8% and
    # sqrt(0.531) = 72.9%, which contradicts the quoted percentages. The
    # nlmixr2 eta blocks therefore take omega^2.
    etalcl ~ 0.198916
    # omega CL/F = 0.446 (44.6% CV), Wang 2019 Table III; 0.446^2 = 0.198916
    etalvc ~ 0.281961
    # omega V/F = 0.531 (53.1% CV), Wang 2019 Table III; 0.531^2 = 0.281961

    # Residual error: OB = IP * (1 + eps1) + eps2 (Wang 2019 Methods
    # 'Random-effects model'), i.e. combined proportional plus additive.
    # sigma1 and sigma2 are on the same standard-deviation scale as the omega
    # column above; see the vignette Assumptions and deviations section.
    propSd <- 0.117
    label("Proportional residual error (fraction); 11.7% per Wang 2019 Table III") # Table III 'sigma1 0.117'
    addSd <- 8.062
    label("Additive residual error (ng/mL); Wang 2019 Table III") # Table III 'sigma2 8.062'
  })

  model({
    # Individual PK parameters. CL/F carries allometric weight scaling and the
    # linear-deviation spironolactone effect; V/F carries weight scaling only
    # (Wang 2019 Results):
    #   CL/F = theta_CL/F * (WT/70)^0.75 * (1 + spirolactone * theta_spiro)
    #   V/F  = theta_V/F  * (WT/70)
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * (1 + CONMED_SPIRON * e_conmed_spiron_cl)
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Dose in mg and vc in L give central/vc in mg/L = ug/mL; x1000 converts
    # to the ng/mL whole-blood units the assay and Table III report.
    Cc <- central / vc * 1000
    Cc ~ add(addSd) + prop(propSd)
  })
}
