Wu_2026_bosutinib <- function() {
  description <- paste(
    "Two-compartment population PK model for oral bosutinib (a dual Src and",
    "Bcr-Abl tyrosine kinase inhibitor) in children and adolescents (1-17",
    "years) with Philadelphia chromosome-positive chronic myelogenous",
    "leukemia from the phase I part of the ITCC-054/COG AAML1921 trial.",
    "First-order absorption with a lag time, absolute bioavailability fixed",
    "to the adult value (33.1%), and fixed-exponent allometric scaling of",
    "all clearance (0.75) and volume (1) terms by body weight (reference 70",
    "kg). The adult model of Garrett 2023 was the reference model and the",
    "NONMEM $PRIOR (NWPRI) subroutine supplied informative adult priors on",
    "Q, V3 and ka. Correlated between-subject variability on CL and V2, a",
    "fixed correlated block on Q and V3, variability on ka, three-occasion",
    "inter-occasion variability on CL and a proportional residual error.",
    sep = " "
  )
  reference <- paste(
    "Wu J-H, Bautista F, Brivio E, Willemse M, Kowalski K, Jermain B,",
    "Garrett M, et al. Population",
    "Pharmacokinetics of Bosutinib in Children with Philadelphia",
    "Chromosome-Positive Chronic Myelogenous Leukemia Based on a Phase I",
    "Dose-Finding Trial (ITCC-054/COG AAML1921) and Pediatric Recommended",
    "Phase II Dosage Evaluation. Clin Pharmacokinet. 2026;65(9):1389-1403.",
    "doi:10.1007/s40262-026-01672-y",
    sep = " "
  )
  vignette <- "Wu_2026_bosutinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "bosutinib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "bosutinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "bosutinib", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Fixed-exponent allometric scaling with reference weight 70 kg:",
        "exponent 0.75 on CL and Q, 1 on V2 and V3 (Supplementary Table 4",
        "control stream: TVCL = THETA(1)*CLCOV*(BWT/70)**0.75, TVV2 =",
        "THETA(2)*V2COV*(BWT/70)**1, TVQ = THETA(3)*(BWT/70)**0.75, TVV3 =",
        "THETA(4)*(BWT/70)**1). Baseline weight in the source (BWT);",
        "median 44.8 kg, range 12.0-115 kg (Table 3)."
      ),
      source_name = "BWT"
    ),
    OCC = list(
      description = "Occasion index for inter-occasion variability on clearance",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Supplementary Table 4 control stream: FLAG1 = OCC == 1, FLAG2 =",
        "OCC == 2, FLAG3 = OCC >= 3, IOVCL = FLAG1*ETA(8) + FLAG2*ETA(9) +",
        "FLAG3*ETA(10), with $OMEGA BLOCK(1) 0.221512 followed by two",
        "BLOCK(1) SAME. Decomposed inside model() into oc1, oc2 and oc3",
        "(the last covering every OCC >= 3). Records with OCC < 1 carry no",
        "inter-occasion term. For a single-occasion simulation pass OCC = 1."
      ),
      source_name = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Adult-model power effect of age on V2, (AGE/51.0)^0.178",
        "(Supplementary Table 2). Removed from the pediatric model because",
        "the OFV increase on exclusion was not significant (Results); the",
        "final control stream keeps the term with its exponent fixed to 0",
        "(THETA(7) V2AGE = 0 FIX), i.e. no effect."
      ),
      source_name = "AGE"
    ),
    CRCL_BASE = list(
      description = "Baseline creatinine clearance",
      units = "not stated in the source (adult-model reference value 98.62)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Adult-model power effect of baseline creatinine clearance on CL,",
        "(BCCL/98.62)^0.198 (Supplementary Tables 2 and 3). Removed from the",
        "pediatric model (OFV increase on exclusion not significant, and",
        "bosutinib is minimally renally cleared; Results); the final control",
        "stream keeps the term with its exponent fixed to 0 (THETA(8) CLBCCL",
        "= 0 FIX), i.e. no effect."
      ),
      source_name = "BCCL"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Asian)",
      notes = paste(
        "Adult-model effects of Asian race on CL (1 - 0.184) and V2",
        "(1 - 0.202) (Supplementary Table 2). Removed because no Asian",
        "patients were in the pediatric dataset (Results); absent from the",
        "final control stream."
      ),
      source_name = "ASIA"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 26L,
    n_studies = 1L,
    n_observations = 235L,
    age_range = "1-17 years",
    age_median = "11.5 years",
    weight_range = "12.0-115 kg",
    weight_median = "44.8 kg",
    bsa_range = "0.530-2.38 m^2",
    bsa_median = "1.32 m^2",
    sex_female_pct = 42.3,
    age_groups = c(`infants (1 year)` = 1L, `children (1 to <12 years)` = 12L, `adolescents (12-18 years)` = 13L),
    disease_state = paste(
      "Chronic-phase Philadelphia chromosome-positive chronic myelogenous",
      "leukemia, resistant or intolerant to prior tyrosine kinase inhibitor",
      "therapy (phase I part of ITCC-054/COG AAML1921, NCT04258943)."
    ),
    dose_range = paste(
      "Oral bosutinib once daily at 300 mg/m^2 (max 500 mg; n = 6), 350",
      "mg/m^2 (max 600 mg; n = 11) or 400 mg/m^2 (max 600 mg; n = 9) in",
      "28-day cycles, as tablets and/or capsules (whole, dissolved or",
      "sprinkled on food), fed."
    ),
    regions = "Multicentre ITCC (Innovative Therapies for Children with Cancer) / Children's Oncology Group trial.",
    notes = paste(
      "Demographics from Table 3 (N = 26: 15 male, 11 female). Sampling",
      "(Table 1): cycle 1 day 14 pre-dose and 1, 3, 6, 8 h post dose, cycle",
      "1 day 15 (24 h), and pre-dose on cycle 2, 3 and 4 day 1 and at dose",
      "escalation. Parameters are typical values for a 70-kg individual."
    )
  )

  ini({
    # Structural PK -- Table 4 'Estimates' column (typical values at 70 kg).
    # CL, V2, Q and V3 are systemic (bioavailability-corrected) values because
    # F1 is a separate fixed theta; CL/F = 74.9 / 0.331 = 226 L/h (Results).
    lcl <- log(74.9) ; label("Clearance CL at 70 kg (L/h)") # Table 4 CL = 74.9 L/h (RSE 7.30%)
    lvc <- log(1110) ; label("Central volume of distribution V2 at 70 kg (L)") # Table 4 V2 = 1110 L (RSE 9.00%)
    lq <- log(37.4) ; label("Intercompartmental clearance Q at 70 kg (L/h)") # Table 4 Q = 37.4 L/h (RSE 7.00%); informative adult prior
    lvp <- log(1960) ; label("Peripheral volume of distribution V3 at 70 kg (L)") # Table 4 V3 = 1960 L (RSE 5.60%); informative adult prior
    lka <- log(0.366) ; label("First-order absorption rate constant ka (1/h)") # Table 4 KA = 0.366 1/h (RSE 6.00%); informative adult prior
    ltlag <- log(0.467) ; label("Absorption lag time (h)") # Table 4 LAG TIME = 0.467 h (RSE 11.2%)
    lfdepot <- fixed(log(0.331)) ; label("Absolute oral bioavailability (fraction)") # Table 4 Bioavailability = 33.1% (FIX, adult value); control stream THETA(9) 0.331 FIX

    # Fixed allometric exponents (Methods; Supplementary Table 4 control stream).
    e_wt_cl <- fixed(0.75) ; label("Allometric exponent of body weight on CL (unitless)") # control stream TVCL = ...*(BWT/70)**0.75
    e_wt_vc <- fixed(1) ; label("Allometric exponent of body weight on V2 (unitless)") # control stream TVV2 = ...*(BWT/70)**1
    e_wt_q <- fixed(0.75) ; label("Allometric exponent of body weight on Q (unitless)") # control stream TVQ = ...*(BWT/70)**0.75
    e_wt_vp <- fixed(1) ; label("Allometric exponent of body weight on V3 (unitless)") # control stream TVV3 = ...*(BWT/70)**1

    # Between-subject variability. Table 4 reports IIV as %CV = 100*sqrt(omega^2)
    # (the fixed Q/V3 variances 0.4 are printed as 63.2% = 100*sqrt(0.4)), so
    # variance = (CV/100)^2. All parameters are EXP(MU + ETA).
    etalcl + etalvc ~ c(0.097344, 0.0696, 0.072361) # Table 4 IIV CL 31.2% -> 0.312^2; Covariance (CL, V2) 0.0696; IIV V2 26.9% -> 0.269^2
    etalq + etalvp ~ fixed(c(0.4, 0.1, 0.4)) # Table 4 IIV Q and V3 63.2% (FIX), Covariance (Q, V3) 0.1 (FIX); control stream $OMEGA BLOCK(2) FIX 0.4 / 0.1 0.4
    etalka ~ 0.840889 # Table 4 IIV KA 91.7% -> 0.917^2

    # Inter-occasion variability on CL: one eta per occasion (OCC 1, 2, >= 3)
    # sharing one variance ($OMEGA BLOCK(1) 0.221512 then two BLOCK(1) SAME).
    etaiov_cl_1 ~ 0.221841 # Table 4 IOV CL 47.1% -> 0.471^2
    etaiov_cl_2 ~ fixed(0.221841) # BLOCK(1) SAME: equal to the occasion-1 IOV variance
    etaiov_cl_3 ~ fixed(0.221841) # BLOCK(1) SAME: equal to the occasion-1 IOV variance

    # Residual error: $ERROR Y = F*(1+EPS(1)).
    propSd <- 0.31 ; label("Proportional residual error (fraction)") # Table 4 'Proportional error oral (untransformed variance)' = 0.0961; sqrt(0.0961) = 0.31
  })

  model({
    # Occasion indicators (control stream FLAG1..FLAG3; FLAG3 covers OCC >= 3).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC >= 3)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3

    # Individual parameters. The ETAs on ALAG1 and F1 are 0 FIX in the control
    # stream and are omitted here.
    cl <- exp(lcl + etalcl + iov_cl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q <- exp(lq + etalq) * (WT / 70)^e_wt_q
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vp
    ka <- exp(lka + etalka)
    tlag <- exp(ltlag)
    fdepot <- exp(lfdepot)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fdepot
    alag(depot) <- tlag

    # Dose in mg and V2 in L give mg/L; x1000 gives ng/mL (the control stream
    # uses dose in ug with S2 = V2, i.e. ug/L = ng/mL).
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
