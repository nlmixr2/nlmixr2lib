Ravix_2024_trametinib <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order absorption and linear",
    "elimination for oral trametinib in adults treated for solid tumours",
    "(Ravix 2024, OpTAT study). Because no intravenous data were collected, all",
    "disposition parameters are apparent oral values (CL/F, Vc/F, Q/F, Vp/F) and",
    "no bioavailability term is estimated. The absorption rate constant was not",
    "identifiable from the sparse real-life sampling and was fixed to a literature",
    "value of 0.913 /h, chosen to give a tmax near the labelled 1.5 h.",
    "Apparent clearance carries two linear median-centred covariate effects:",
    "clearance rises with fat-free mass and falls with age. Inter-individual",
    "variability was supported on clearance only; residual error is proportional.",
    "The paper also identified a clearance increase with concomitant dabrafenib,",
    "but a sensitivity analysis restricted to the 11 patients whose dose history",
    "was electronically monitored showed that effect to be confounded by",
    "non-adherence, so it was removed from the final model and is not encoded here.",
    sep = " "
  )
  reference <- paste(
    "Ravix A, Bandiera C, Cardoso E, Lata-Pedreira A, Chtioui H, Decosterd LA,",
    "Wagner AD, Schneider MP, Csajka C, Guidi M.",
    "Population Pharmacokinetics of Trametinib and Impact of Nonadherence on Drug",
    "Exposure in Oncology Patients as Part of the Optimizing Oral Targeted",
    "Anticancer Therapies Study.",
    "Cancers. 2024;16(12):2193.",
    "doi:10.3390/cancers16122193.",
    "Parameter estimates from Table 2; the clearance covariate equation from the",
    "Table 2 footnote; population characteristics from Table 1.",
    sep = " "
  )
  vignette <- "Ravix_2024_trametinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot       = list(analyte = "trametinib", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "trametinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "trametinib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    FFM = list(
      description        = "Fat-free mass",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Median-centred at FFM_median = 46.35 kg (Ravix 2024 Table 1 cohort median",
        "and Table 2 footnote) and applied as a linear deviation,",
        "cl *= (1 + 1.41 * (FFM - 46.35) / 46.35).",
        "IMPORTANT -- this paper does NOT use the Janmahasatian equation that the",
        "register entry names as the usual FFM derivation. Ravix 2024 Methods",
        "Section 2.1 Equations (1) and (2) derive FFM from the Deurenberg",
        "body-fat-percentage equation: FFM (kg) = BW * (1 - BFPbmi / 100), where",
        "BFPbmi (%) = 1.20 * BMI + 0.23 * AGE - 10.8 * SEX - 5.4, with BW in kg,",
        "BMI in kg/m^2, AGE in years and SEX = 1 for men and 0 for women",
        "(Deurenberg P, Weststrate JA, Seidell JC. Br J Nutr. 1991;65(2):105-114).",
        "A downstream user supplying this column must use the Deurenberg form to",
        "stay on the scale the coefficient was estimated on; the two equations do",
        "not give the same number for the same subject. Note that the Deurenberg",
        "equation makes FFM depend on AGE, so the FFM and AGE effects below are not",
        "fully independent covariates. Cohort median (min, max) 46.35 (32, 68) kg;",
        "4% of records had a missing value and were imputed at the population median",
        "(Ravix 2024 Methods Section 2.2.1).",
        sep = " "
      ),
      source_name        = "FFM"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Median-centred at AGE_median = 63 years (Ravix 2024 Table 1 cohort median",
        "and Table 2 footnote) and applied as a linear deviation,",
        "cl *= (1 + (-0.69) * (AGE - 63) / 63), so apparent clearance DECREASES with",
        "age. Cohort median (min, max) 63 (30, 85) years, no missing values.",
        "The linear form is unbounded and would drive clearance to zero at about",
        "154 years, far outside the observed 30-85 year range; do not extrapolate.",
        sep = " "
      ),
      source_name        = "AGE"
    )
  )

  # Covariates that Ravix 2024 screened in the forward-inclusion / backward-
  # elimination covariate analysis (Methods Section 2.2.1) but did NOT retain in
  # the final model. Documentation only -- none is referenced in model(). Recorded
  # here so the provenance of the covariate screen is not lost.
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Significant on CL in univariate analysis (dOFV = -12.51, p < 0.01) but the",
        "effect disappeared once fat-free mass was in the model (Results Section 3.2),",
        "because the Deurenberg FFM derivation already carries a sex term.",
        sep = " "
      )
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened on CL; not significant (dOFV > -3.40, p > 0.05). Cohort median (min, max) 70 (45, 96) kg."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened on CL; not significant. Cohort median (min, max) 25.3 (17.3, 33.8) kg/m^2. Still needed upstream as an input to the Deurenberg FFM derivation."
    ),
    BSA = list(
      description = "Body surface area",
      units       = "m^2",
      type        = "continuous",
      notes       = "Significant on CL in univariate analysis (dOFV = -7.13, p < 0.01) but not retained in the multivariate model. Cohort median (min, max) 1.78 (1.40, 2.23) m^2."
    ),
    CRCL = list(
      description = "Creatinine clearance, Cockcroft-Gault",
      units       = "mL/min/1.73m^2",
      type        = "continuous",
      notes       = "Significant on CL in univariate analysis (dOFV = -11.15, p < 0.01) but not retained. Cohort median (min, max) 84 (42, 164) mL/min/1.73m^2. Consistent with only 19% urinary excretion of trametinib."
    ),
    BUN = list(
      description = "Serum urea",
      units       = "Not reported",
      type        = "continuous",
      notes       = "Listed among the tested covariates (Methods Section 2.2.1) but not significant and not tabulated in Table 1, so neither its units nor its distribution are reported."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened on CL; not significant. Cohort median (min, max) 5 (3, 12) umol/L."
    ),
    ALP = list(
      description = "Alkaline phosphatase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened on CL; not significant. Cohort median (min, max) 91 (49, 873) U/L. Reported by the paper under the French abbreviation PAL."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened on CL; not significant. Cohort median (min, max) 30 (12, 211) U/L. Reported by the paper under the French abbreviation ASAT."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened on CL; not significant. Cohort median (min, max) 26 (8, 152) U/L. Reported by the paper under the French abbreviation ALAT."
    ),
    CONMED_DABRAFENIB = list(
      description = "Concomitant dabrafenib coadministration indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "This is the paper's headline methodological finding, and the reason the",
        "covariate is absent from the final model. Concomitant dabrafenib was",
        "significant on CL in univariate analysis (dOFV = -7.23, p < 0.01) and again",
        "when added to the FFM + age model (dOFV = -9.00, p < 0.01), raising apparent",
        "clearance by about 40% though poorly estimated (RSE 46%). Repeating the",
        "analysis in the 11 patients whose dose history was captured by an electronic",
        "pillbox removed the effect entirely, while the FFM and age effects persisted.",
        "The authors concluded that the apparent dabrafenib effect was a spurious",
        "correlation created by unrecorded missed doses being fitted as high clearance,",
        "and excluded it (Results Section 3.2; Discussion). 22 of 33 patients (67%)",
        "took trametinib with dabrafenib.",
        sep = " "
      )
    ),
    FED = list(
      description = "Fed-state indicator at the time of dosing",
      units       = "(binary)",
      type        = "binary",
      notes       = "Listed among the tested covariates (Methods Section 2.2.1) as 'food intake'; not significant on CL and not retained."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 33L,
    n_studies      = 1L,
    n_observations = 113L,
    age_range      = "30-85 years",
    age_median     = "63 years",
    weight_range   = "45-96 kg",
    weight_median  = "70 kg",
    sex_female_pct = 45,
    disease_state  = paste(
      "Adults treated for a solid tumour (Ravix 2024 Table 1): melanoma 22 (67%),",
      "ovarian cancer 3 (9%), breast cancer 2 (6%), cholangiocarcinoma 2 (6%), and",
      "one patient each (3%) with thyroid carcinoma, gastrointestinal stromal",
      "tumour, hepatocellular carcinoma and ileocecal carcinoma.",
      sep = " "
    ),
    dose_range     = "Trametinib 0.5-2 mg orally once daily; 7 patients (21%) on monotherapy and 22 (67%) with concomitant dabrafenib",
    regions        = "Switzerland (Lausanne University Hospital; OpTAT study, ClinicalTrials.gov NCT04484064)",
    notes          = paste(
      "Real-life therapeutic-drug-monitoring data, sparse and opportunistic: median",
      "(min, max) 3 (1, 11) samples per patient drawn 5 h (0.13, 202 h) after the",
      "reported intake, with a maximum of eight samples per patient. Two patients",
      "consented to a richer profile of eight samples over 24 h. Assay LLOQ 1 ng/mL",
      "by LC-MS/MS. Dose history was electronically captured with a MEMS digital",
      "pillbox in 11 patients ('the full adherence information group'); for the",
      "remaining patients it was reconstructed from consultation notes, assuming",
      "steady state where no information existed. Missing covariate values were",
      "imputed at the population median. Fit in NONMEM 7.4.3 with PsN 4.8.0;",
      "internally validated by a 2000-replicate bootstrap and a prediction-corrected",
      "VPC of 1000 simulations. Relatively few samples were collected beyond 15 h,",
      "which the authors flag as a likely source of bias in the distribution volumes",
      "(Discussion).",
      sep = " "
    )
  )

  ini({
    # ---- Structural fixed effects (Ravix 2024 Table 2, 'Final Model Estimation (RSE %)') ----
    # All disposition parameters are APPARENT oral values: the study collected no
    # intravenous data and did not estimate a bioavailability term, so CL is CL/F,
    # V2 is Vc/F, Q is Q/F and V3 is Vp/F.
    lka <- fixed(log(0.913)); label("First-order absorption rate constant (1/h)")            # Table 2 'ka (h-1) 0.913 fixed'; not estimable from the sparse absorption-phase data. Results 3.2 and the Discussion both attribute the fixed value to reference [38], which the paper's reference list gives as Balakirouchenane et al. (candidate values of 0.4-2 /h were screened, Methods 2.2.1). Yields a median (min, max) tmax of 1.75 h (1.51, 1.84) vs the 1.5 h FDA reference (Results 3.2).
    lcl <- log(3.96);         label("Apparent clearance at the covariate reference (L/h)") # Table 2 'thetaCL (L.h-1) 3.96 (6)'; bootstrap median 3.98 (95% PI 3.52, 4.45). This is the typical value at AGE = 63 y and FFM = 46.35 kg, where both covariate factors equal 1.
    lvc <- log(108);          label("Apparent central volume of distribution (L)")      # Table 2 'V2 (L) 108 (16)'; bootstrap median 101.84 (95% PI 65.42, 143.40)
    lq  <- log(29.4);         label("Apparent inter-compartmental clearance (L/h)")      # Table 2 'Q (L.h-1) 29.4 (30)'; bootstrap median 28.83 (95% PI 12.69, 60.68)
    lvp <- log(286);          label("Apparent peripheral volume of distribution (L)")   # Table 2 'V3 (L) 286 (25)'; bootstrap median 286.48 (95% PI 104.63, 385.40)

    # ---- Covariate effects on apparent clearance (Ravix 2024 Table 2 and its footnote) ----
    # The footnote prints the final clearance equation as a product of two linear
    # median-centred deviation terms:
    #   CL = thetaCL * (1 + thetaAGE * (AGE - AGE_median) / AGE_median)
    #               * (1 + thetaFFM * (FFM - FFM_median) / FFM_median)
    # with AGE_median = 63 years and FFM_median = 46.35 kg. Both coefficients
    # reproduce the paper's own worked numbers exactly, which is what confirms this
    # reading of the equation:
    #   FFM at the observed maximum of 68 kg (Results 3.2, "a 66% increase in CL"):
    #     1.41 * (68 - 46.35) / 46.35 = 0.6586, i.e. +65.9%
    #   AGE at 30 y and 80 y (Results 3.2, "5.39 ... versus 3.22 ... (40% decrease)"):
    #     3.96 * (1 - 0.69 * (30 - 63) / 63) = 5.391 L/h
    #     3.96 * (1 - 0.69 * (80 - 63) / 63) = 3.223 L/h, a 40.2% decrease
    e_age_cl <- -0.69; label("Effect of age on apparent clearance (fractional change per unit relative deviation from 63 years)") # Table 2 'thetaAGE -0.69 (32)'; bootstrap median -0.66 (95% PI -1.25 to -0.20)
    e_ffm_cl <-  1.41; label("Effect of fat-free mass on apparent clearance (fractional change per unit relative deviation from 46.35 kg)") # Table 2 'thetaFFM 1.41 (17)'; bootstrap median 1.41 (95% PI 0.92, 1.90)

    # ---- Inter-individual variability (Ravix 2024 Table 2) ----
    # Supported on clearance only: adding IIV to the other PK parameters gave
    # estimates close to zero and no significant fit improvement (dOFV > -3.4,
    # p > 0.05; Results 3.2). Reported as a CV of 23%, so on the exponential-eta
    # scale omega^2 = log(1 + CV^2) = log(1 + 0.23^2) = 0.0515459.
    etalcl ~ 0.0515459  # Table 2 'IIV_CL (%) 23 (14)'; bootstrap median 22% (95% PI 13, 28). Down from 45% in the covariate-free base model (Results 3.2).

    # ---- Residual error (Ravix 2024 Table 2) ----
    propSd <- 0.20; label("Proportional residual error (fraction)")  # Table 2 'sigma_prop (%) 20 (10)', expressed as CV (%); bootstrap median 20% (95% PI 16, 24). A proportional model was retained over additive and combined alternatives (Methods 2.2.1, Results 3.2).
  })

  model({
    # Amounts are carried in mg and volumes in L, so central / vc is mg/L. The
    # source reports trametinib concentrations in ng/mL (assay LLOQ 1 ng/mL, the
    # therapeutic Cmin target 10.6 ng/mL), so convert: 1 mg/L = 1000 ng/mL.
    mgPerLToNgPerMl <- 1000

    # ---- Covariate reference values (Ravix 2024 Table 2 footnote; both are the
    #      Table 1 cohort medians) ----
    ageMedian <- 63     # years
    ffmMedian <- 46.35  # kg

    # ---- Derived multiplicative covariate factors ----
    # Linear median-centred deviations; each evaluates to 1 for a subject at the
    # cohort median, where exp(lcl) is the typical apparent clearance.
    cl_age_factor <- 1 + e_age_cl * (AGE - ageMedian) / ageMedian
    cl_ffm_factor <- 1 + e_ffm_cl * (FFM - ffmMedian) / ffmMedian

    # ---- Individual PK parameters ----
    # Inter-individual variability is on clearance only.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * cl_age_factor * cl_ffm_factor
    vc <- exp(lvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    # ---- Two-compartment micro-constants ----
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ---- ODE system ----
    # Oral dosing only; no bioavailability term is applied because CL, V2, Q and V3
    # are already apparent (F-scaled) values.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                               k12 * central - k21 * peripheral1

    # ---- Observation and residual error ----
    Cc <- central / vc * mgPerLToNgPerMl
    Cc ~ prop(propSd)
  })
}
