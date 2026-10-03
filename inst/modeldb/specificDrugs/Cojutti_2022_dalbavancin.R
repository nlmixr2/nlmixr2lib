Cojutti_2022_dalbavancin <- function() {
  description <- paste(
    "Two-compartment intravenous population PK model for dalbavancin in adults receiving",
    "long-term, therapeutic-drug-monitoring-guided dalbavancin for subacute or chronic",
    "Gram-positive infections (mostly bone and joint infections, plus endocarditis and",
    "endovascular prosthetic infections). Dalbavancin clearance rises exponentially with",
    "CKD-EPI creatinine clearance; the effect is UNCENTERED, so exp(lcl) is the clearance",
    "intercept extrapolated to CLCR = 0, not a clearance at any physiological renal function.",
    sep = " "
  )
  reference <- paste(
    "Cojutti PG, Tedeschi S, Gatti M, Zamparini E, Meschiari M, Siega PD, Mazzitelli M, Soavi L,",
    "Binazzi R, Erne EM, Rizzi M, Cattelan AM, Tascini C, Mussini C, Viale P, Pea F. Population",
    "Pharmacokinetic and Pharmacodynamic Analysis of Dalbavancin for Long-Term Treatment of",
    "Subacute and/or Chronic Infectious Diseases: The Major Role of Therapeutic Drug Monitoring.",
    "Antibiotics (Basel). 2022;11(8):996. doi:10.3390/antibiotics11080996",
    sep = " "
  )
  vignette <- "Cojutti_2022_dalbavancin"

  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "dalbavancin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "dalbavancin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance estimated by the CKD-EPI equation, BSA-normalized",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Source column CLCR. Methods 4.1: 'Creatinine clearance was estimated by means of the",
        "CKD-EPI equation'; Table 1 and every renal-function class in Results 2.3 are reported in",
        "mL/min/1.73 m^2 (cohort median 93.0, IQR 72.0-104.0, range 3.0-141.0). The effect is",
        "applied UNCENTERED and EXPONENTIALLY, CL = CLpop x exp(beta x CLCR), even though Methods",
        "4.2 calls continuous-covariate effects a 'power function'. The paper's own numbers close",
        "the exponential reading: 0.029 x exp(0.0043 x 93) = 0.0433 L/h against the reported",
        "median individual CL of 0.043 L/h (Results 2.2, repeated in the Discussion), whereas a",
        "power term with exponent 0.0043 would leave CL at about 0.029 L/h for any CLCR and could",
        "not have cut the CL random effect from 31.76 to 26.44 %CV. The same group's 2024",
        "dalbavancin model prints this form verbatim (CL = 0.030 x e^(0.0042 x eGFR); see",
        "modellib('Cojutti_2024_dalbavancin')). Only CL carries a covariate."
      ),
      source_name = "CLCR"
    )
  )

  # Covariates screened (Methods 4.2: 'age, gender, weight, height, serum
  # albumin, serum creatinine and CLCR') but not retained; Results 2.2: 'CLCR was
  # the only covariate significantly associated with dalbavancin clearance'.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened on the PK parameters (Methods 4.2) but not retained. Table 1 median 62 (IQR 51-73) years."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened as 'gender' (Methods 4.2) but not retained. Table 1 male/female 44/25."
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened (Methods 4.2) but not retained; no allometric scaling is applied. Table 1 median 75 (IQR 62-88) kg."
    ),
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      notes = "Screened (Methods 4.2) but not retained. Table 1 median 170 (IQR 165-177) cm."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/dL",
      type = "continuous",
      notes = "Screened (Methods 4.2) but not retained. Table 1 median 3.7 (IQR 3.3-4.0) g/dL."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened (Methods 4.2) but not retained; renal function entered through the derived CKD-EPI CLCR instead. Not summarised in Table 1."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 69L,
    n_studies = 1L,
    n_observations = "289 total dalbavancin plasma concentrations (Results 2.2); median 3 (range 1-19) TDM samples per patient",
    age_range = "19-90 years",
    age_median = "62 years",
    weight_range = "42-143 kg",
    weight_median = "75 kg",
    sex_female_pct = 36.2,
    race_ethnicity = "Not reported; multicentre Italian cohort.",
    disease_state = paste(
      "Adults with documented or suspected Gram-positive infections treated with dalbavancin as",
      "second-line, long-term therapy. Table 1: prosthetic joint infection 26 (37.7%),",
      "osteomyelitis 11 (15.9%), endovascular prosthetic infection 9 (13.0%), endocarditis 7",
      "(10.1%), spondylodiscitis 5 (7.2%), infected pseudoarthrosis non-unions 4 (5.8%), septic",
      "arthritis 1 (1.5%), and 6 patients with two infection sites. 63/69 (91.3%) had a",
      "microbiological isolate; MRSA and methicillin-resistant S. epidermidis made up 55/74 isolates."
    ),
    renal_function = "CKD-EPI CLCR median 93.0 mL/min/1.73 m^2, IQR 72.0-104.0, range 3.0-141.0 (Results 2.1, Table 1).",
    albumin = "Median 3.7 g/dL, range 2.5-4.6 (Results 2.1).",
    dose_range = paste(
      "All patients started with two intravenous doses on days 1 and 8 of 1000 mg or 1500 mg",
      "(Methods 4.1); 32/69 received exactly two doses (27 of them 1500 mg one week apart), 17",
      "received three and 20 received 4-14 doses (Results 2.1). Infusion duration is not stated."
    ),
    regions = "Italy (Bologna, Modena, Udine, Padua, Bergamo, Bolzano)",
    notes = paste(
      "Retrospective study, April 2021 to April 2022 (Ethics Committee 897/2021/Oss/AOUBo).",
      "Total plasma dalbavancin by LC-MS/MS, LLOQ 0.5 mg/L. Estimation by SAEM in Monolix",
      "2021R1; Monte Carlo simulations in Simulx 2020R1. All individual parameters log-normal",
      "with exponential random effects (Methods 4.2)."
    )
  )

  ini({
    # Structural parameters -- Cojutti 2022 Table 2, Final Model column.
    # Monolix estimated them on the log scale (Methods 4.2: 'All individual
    # parameters were log-normally distributed'). The paper's V1 / V2 map onto
    # the canonical vc / vp. CL is the intercept of the uncentered exponential
    # CLCR relationship (see covariateData$CRCL$notes).
    lcl <- log(0.029); label("Clearance intercept at CRCL = 0 (L/h)") # Cojutti 2022 Table 2 Final Model: CL = 0.029 L/h (RSE 11.6%)
    lvc <- log(6.14); label("Central volume of distribution V1 (L)") # Cojutti 2022 Table 2 Final Model: V1 = 6.14 L (RSE 5.26%)
    lq <- log(0.026); label("Intercompartmental clearance Q (L/h)") # Cojutti 2022 Table 2 Final Model: Q = 0.026 L/h (RSE 18.1%)
    lvp <- log(9.52); label("Peripheral volume of distribution V2 (L)") # Cojutti 2022 Table 2 Final Model: V2 = 9.52 L (RSE 19.0%)

    # Covariate effect, applied as cl = exp(lcl) * exp(e_crcl_cl * CRCL).
    e_crcl_cl <- 0.0043; label("Exponential coefficient of CRCL on CL (per mL/min/1.73 m^2; uncentered)") # Cojutti 2022 Table 2 Final Model: beta CLcr-CL = 0.0043 (RSE 28.9%)

    # Inter-individual variability. Table 2 heads this block 'Random Effects
    # (Inter-patient %CV)', so each value is a coefficient of variation of a
    # log-normal parameter and the log-scale variance is omega^2 = log(1 + CV^2).
    # No correlations are reported.
    etalcl ~ 0.067572 # Cojutti 2022 Table 2 Final Model: IIV CL = 26.44 %CV (RSE 13.3%); log(1 + 0.2644^2)
    etalvc ~ 0.025591 # Cojutti 2022 Table 2 Final Model: IIV V1 = 16.10 %CV (RSE 40.1%); log(1 + 0.1610^2)
    etalq ~ 0.230382 # Cojutti 2022 Table 2 Final Model: IIV Q = 50.90 %CV (RSE 32.9%); log(1 + 0.5090^2)
    etalvp ~ 0.129283 # Cojutti 2022 Table 2 Final Model: IIV V2 = 37.15 %CV (RSE 61.7%); log(1 + 0.3715^2)

    # Residual variability -- Table 2 lists a single proportional term b.
    propSd <- 0.3392; label("Proportional residual error (fraction)") # Cojutti 2022 Table 2 Final Model: b (proportional) = 33.92% (RSE 5.97%)
  })

  model({
    # CLCR effect on CL is exponential and uncentered; V1, Q and V2 carry no
    # covariate.
    cl <- exp(lcl + etalcl) * exp(e_crcl_cl * CRCL)
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Intravenous infusion into central; no absorption compartment.
    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Total (not free) plasma dalbavancin, mg/L (dose mg, volume L).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
