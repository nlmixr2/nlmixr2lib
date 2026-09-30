Sukeishi_2022_remdesivir <- function() {
  description <- paste(
    "One-compartment population PK model for GS-441524, the predominant",
    "circulating nucleoside metabolite of intravenous remdesivir, in Japanese",
    "adults hospitalised with COVID-19 including patients with renal",
    "dysfunction and on ECMO (Sukeishi 2022). Remdesivir itself is not",
    "modelled: its half-life is under 1 h, so each remdesivir infusion is",
    "treated as a direct input of GS-441524 into the central compartment,",
    "converted mole-for-mole by the molecular-weight ratio 291.3/602.6",
    "(complete conversion assumed). Clearance scales with absolute",
    "(non-BSA-indexed) eGFR as a power function referenced to the cohort",
    "median of 74.7 mL/min; the volume of distribution is 42.9% lower in",
    "patients aged 75 years or older. Log-normal between-subject",
    "variability on CL and V; proportional residual error."
  )
  reference <- "Sukeishi A, Itohara K, Yonezawa A, Sato Y, Matsumura K, Katada Y, Nakagawa T, Hamada S, Tanabe N, Imoto E, Kai S, Hirai T, Yanagita M, Ohtsuru S, Terada T, Ito I. Population pharmacokinetic modeling of GS-441524, the active metabolite of remdesivir, in Japanese COVID-19 patients with renal dysfunction. CPT Pharmacometrics Syst Pharmacol. 2022;11(1):94-103. doi:10.1002/psp4.12736"
  vignette <- "Sukeishi_2022_remdesivir"

  # Doses are entered as mg of REMDESIVIR (the administered prodrug), as in
  # the source control stream (Supporting Information). The
  # model converts them to mg of GS-441524 through f(central), so the
  # `central` state holds a GS-441524 amount in mg, and central / vc (mg/L)
  # times 1000 is ng/mL, the scale on which Sukeishi 2022 reports serum
  # GS-441524 concentrations (Table 1, Table 3).
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    central = list(analyte = "GS-441524", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Estimated glomerular filtration rate, BSA-normalized (eGFR indexed), from the Japanese-coefficient MDRD equation on standardized serum creatinine.",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Sukeishi 2022 Methods 'Patients and data collection': the clinical laboratory reported eGFR indexed (mL/min/1.73 m^2) from the Modification of Diet in Renal Disease equation for Japanese, and the authors de-indexed it to absolute mL/min as eGFR non-indexed = eGFR indexed / 1.73 * BSA, with BSA from the Du Bois equation. The final model uses the ABSOLUTE value (Table 2: CL = theta_CL * (eGFR non-indexed / 74.7)^theta), which gave a better fit than eGFR indexed, creatinine clearance or serum creatinine (Results 'PopPK model and model evaluation'). This file therefore takes the canonical BSA-normalized CRCL plus BSA and performs the same de-indexing inside model(); supply CRCL on the mL/min/1.73 m^2 scale, not absolute mL/min, or the renal term is rescaled by BSA/1.73. Time-varying in the source: the covariate was the value at each concentration measurement point. Observed eGFR non-indexed range 16.4-147.7 mL/min (median 74.7). In the Supporting Information control stream the absolute value is the data column EGFR2 (CL = THETA(1)*(EGFR2/74.7)**THETA(3)*EXP(ETA(1))).",
      source_name = "eGFRindexed"
    ),
    BSA = list(
      description = "Body surface area (Du Bois equation)",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Used only to de-index eGFR (eGFR non-indexed = CRCL / 1.73 * BSA), per Sukeishi 2022 Methods 'Patients and data collection'. BSA was also screened directly on CL and V but not retained. Cohort median 1.8 m^2 (range 1.24-2.21), Table 1.",
      source_name = "BSA"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters only through the binary indicator AGE >= 75 years on the volume of distribution (Sukeishi 2022 Table 2 footnote: 'AGE is 1 if a patient is 75 years old or more, or 0 if a patient is under 75 years of age'; Methods 'PopPK modeling' tests age as this categorical covariate on Vd). The indicator is derived inside model() from the continuous canonical AGE; the Supporting Information control stream carries it pre-derived in its AGE data column (V = THETA(2)*(1+THETA(4)*AGE)*EXP(ETA(2))). Cohort median 72 years (range 45-97), Table 1.",
      source_name = "AGE"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened on CL and Vd in the stepwise covariate search (Sukeishi 2022 Methods 'PopPK modeling', Table S1) but not retained. Cohort median 66.8 kg (range 36.7-96.3), Table 1."
    ),
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      notes = "Screened on CL and Vd but not retained in the final (eGFR) model; height was retained on CL only in the alternative creatinine-based model that the authors did not select (Results 'PopPK model and model evaluation', Table S2). Cohort median 167.1 cm (range 144-182), Table 1."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "categorical",
      notes = "Screened on CL as a categorical covariate (1 = female) but not retained (Methods 'PopPK modeling'). 10 of 37 patients female, Table 1."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      notes = "Retained on CL in an alternative model (with height and age on CL and age on Vd; OBJ 1474) that the authors rejected in favour of the eGFR non-indexed model (OBJ 1485) 'considering clinical usefulness' (Results 'PopPK model and model evaluation', Table S2). No point estimates for that alternative model are printed in the main text."
    ),
    AST = list(
      description = "Serum aspartate aminotransferase",
      units = "IU/L",
      type = "continuous",
      notes = "Screened on CL but not retained. Cohort median 31 IU/L (range 12-229), Table 1."
    ),
    ALT = list(
      description = "Serum alanine aminotransferase",
      units = "IU/L",
      type = "continuous",
      notes = "Screened on CL but not retained. Cohort median 30 IU/L (range 4-191), Table 1."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/dL",
      type = "continuous",
      notes = "Screened on CL but not retained. Cohort median 2.5 g/dL (range 1.7-4.1), Table 1."
    ),
    ECMO_STATUS = list(
      description = "Extracorporeal membrane oxygenation treatment-status indicator",
      units = "(binary)",
      type = "categorical",
      notes = "Screened on CL and Vd (1 = on ECMO) but not retained; the paper concludes ECMO 'hardly affected' GS-441524 CL and Vd (Abstract; Results). 4 patients, 21 of 190 measurements on ECMO, Table 1."
    ),
    MECH_VENT = list(
      description = "Invasive mechanical ventilation status indicator",
      units = "(binary)",
      type = "categorical",
      notes = "Screened on CL and Vd (1 = on a ventilator) but not retained (Methods 'PopPK modeling'). 82 of 190 measurements on a ventilator, Table 1."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 37L,
    n_studies = 1L,
    n_observations = "190 serum GS-441524 trough concentrations (Table 1).",
    age_range = "45-97 years",
    age_median = "72 years",
    weight_range = "36.7-96.3 kg",
    weight_median = "66.8 kg",
    sex_female_pct = 27.0,
    race_ethnicity = c(Japanese = 100),
    disease_state = "Hospitalised COVID-19 treated with remdesivir; 4 patients on extracorporeal membrane oxygenation (21 measurements) and 82 measurements taken on a ventilator (Table 1).",
    renal_function = "eGFR non-indexed (absolute, BSA-readjusted) median 74.7 mL/min, range 16.4-147.7 across measurement points (Table 1). At the start of remdesivir, 20 patients had eGFR non-indexed >= 60 mL/min, 15 had 30 to < 60 and 2 had < 30 mL/min (Results 'Patients demographics').",
    hepatic_function = "Median AST 31 IU/L (range 12-229), ALT 30 IU/L (4-191), albumin 2.5 g/dL (1.7-4.1) (Table 1).",
    dose_range = "Remdesivir 200 mg intravenously on day 1 then 100 mg once daily from day 2, each infused over 1 h, for eGFR non-indexed >= 30 mL/min; for eGFR non-indexed < 30 mL/min (2 patients), 200 mg on day 1 then 100 mg once every 2 days from day 3 (Methods 'Remdesivir administration').",
    regions = "Single centre, Kyoto University Hospital, Kyoto, Japan; December 2020 to May 2021.",
    notes = "Retrospective study using surplus serum from routine blood tests drawn at trough, so the data carry almost no information on the post-infusion peak. GS-441524 was measured by LC-MS/MS (calibration 5-500 ng/mL). Estimation was FOCE-I in NONMEM 7.5.0 with PsN 5.0.0; the final model was checked with a 500-replicate nonparametric bootstrap and a prediction-corrected VPC."
  )

  ini({
    # Structural parameters at the covariate reference point (eGFR
    # non-indexed = 74.7 mL/min, age < 75 years). The paper calls these CL
    # and Vd without an F or fm divisor; because GS-441524 formation is not
    # observed they are effectively apparent values conditional on the
    # complete mole-for-mole conversion encoded in model().
    lcl <- log(11.8)
    label("GS-441524 clearance CL (L/h) at eGFR non-indexed 74.7 mL/min")         # Sukeishi 2022 Table 2: theta_CL = 11.8 L/h (RSE 4.9%; bootstrap median 11.8, 95% CI 10.6-13.0); Results equation CL = 11.8 * (eGFRnon-indexed / 74.7)^1.09
    lvc <- log(382)
    label("GS-441524 volume of distribution Vd (L) in patients under 75 years")    # Sukeishi 2022 Table 2: theta_V = 382 L (RSE 9.9%; bootstrap median 383, 95% CI 303-461)

    # Covariate effects (Table 2 and the Results final-model equations):
    #   CL (L/h) = 11.8 * (eGFRnon-indexed / 74.7)^1.09
    #   Vd (L)   = 382 * (1 - 0.429 * AGE>=75)
    e_crcl_cl <- 1.09
    label("Power exponent for eGFR non-indexed on CL (unitless)")                  # Sukeishi 2022 Table 2: theta_eGFRnon-indexed,CL = 1.09 (RSE 10.8%; bootstrap median 1.08, 95% CI 0.833-1.37)
    e_age_ge_75_vc <- -0.429
    label("Fractional change in Vd for age >= 75 years (unitless)")                # Sukeishi 2022 Table 2: theta_Age>=75,V = -0.429 (RSE 21.1%; bootstrap median -0.427, 95% CI -0.578 to -0.171)

    # Between-subject variability: exponential model (Results 'PopPK model
    # and model evaluation'). Table 2 reports it as CV%; converted to the
    # log-scale variance with omega^2 = log(CV^2 + 1).
    etalcl ~ 0.065431   # Sukeishi 2022 Table 2: IIV for CL = 26.0 CV% (RSE 9.6%, shrinkage 4.22%) -> log(0.260^2 + 1)
    etalvc ~ 0.112457   # Sukeishi 2022 Table 2: IIV for Vd = 34.5 CV% (RSE 12.8%, shrinkage 24.1%) -> log(0.345^2 + 1)

    # Residual variability: proportional (Results 'PopPK model and model
    # evaluation'), reported as CV% in Table 2.
    propSd <- 0.152
    label("Proportional residual error (fraction)")                                # Sukeishi 2022 Table 2: proportional error = 15.2 CV% (RSE 11.2%, shrinkage 14.5%; bootstrap median 15.0, 95% CI 12.2-18.7)
  })

  model({
    # De-index eGFR to absolute mL/min exactly as the source does (Methods
    # 'Patients and data collection'): eGFRnon-indexed = eGFRindexed /
    # 1.73 m^2 * BSA.
    egfr_abs <- CRCL / 1.73 * BSA

    # Age indicator (Table 2 footnote): 1 if 75 years or older, else 0.
    age_ge_75 <- 0 + (AGE >= 75)

    cl <- exp(lcl + etalcl) * (egfr_abs / 74.7)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (1 + e_age_ge_75_vc * age_ge_75)

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Remdesivir dose (mg) -> GS-441524 amount (mg), mole-for-mole. The
    # main text does not state the dose basis; the control stream in the
    # Supporting Information (ADVAN1 TRANS2) sets
    #   S1 = (V/1000)/(291.3/602.6)
    # with AMT in mg of remdesivir, i.e. every dose is scaled by
    # MW(GS-441524) / MW(remdesivir) = 291.3 / 602.6 before it enters the
    # GS-441524 compartment. Encoding the factor on the dose rather than
    # on the scaling factor gives identical concentrations and keeps
    # `central` a true GS-441524 amount. The paper's own Monte Carlo
    # trough table (Table 3) is reproduced only with this factor (see the
    # vignette). Give the 1-h infusions with `dur = 1`, not `rate = amt`:
    # with a fixed rate, rxode2 applies f() by shortening the infusion to
    # f * amt / rate (about 29 min) instead of scaling the rate.
    mw_ratio <- 291.3 / 602.6
    f(central) <- mw_ratio

    # mg/L -> ng/mL
    Cc <- central / vc * 1000

    Cc ~ prop(propSd)
  })
}
