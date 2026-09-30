Nijstad_2021_clofarabine <- function() {
  description <- paste(
    "Two-compartment IV population PK model for clofarabine in children",
    "(0.5-18 years) and one adult receiving clofarabine-fludarabine-busulfan",
    "myeloablative conditioning before allogeneic hematopoietic cell",
    "transplantation (Nijstad 2021; n = 81, 805 plasma concentrations).",
    "Clearance is split into a non-renal arm (24.0 L/h at 70 kg) and a renal",
    "arm (29.8 L/h at 70 kg and normal renal function) scaled by relative",
    "renal function RF = (absolute eGFR in L/h x 70/WT) / 6 L/h, so",
    "CL = (CL_nonrenal + CL_renal x RF) x (WT/70)^0.75. Central and peripheral",
    "volumes scale with WT^1 and Q with WT^0.75 (exponents fixed, 70 kg",
    "reference). IIV on CL, V1 and Q; inter-occasion variability on CL and V2",
    "with each of the four daily doses as its own occasion; proportional",
    "residual error."
  )
  reference <- paste(
    "Nijstad AL, Nierkens S, Lindemans CA, Boelens JJ, Bierings M,",
    "Versluys AB, van der Elst KCM, Huitema ADR. Population pharmacokinetics",
    "of clofarabine for allogeneic hematopoietic cell transplantation in",
    "paediatric patients. Br J Clin Pharmacol. 2021;87(8):3218-3226.",
    "doi:10.1111/bcp.14738"
  )
  vignette <- "Nijstad_2021_clofarabine"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "clofarabine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "clofarabine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Actual body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Source column BW. Allometric size descriptor on all four structural",
        "parameters, referenced to 70 kg, with exponents fixed a priori at",
        "0.75 for CL and Q and 1 for V1 and V2 (Nijstad 2021 Section 3.2 and",
        "Table 2 parameter formulas). WT also enters the relative renal",
        "function RF through the factor 70/BW (Section 2.3 Eq. 2). BSA and",
        "fat-free mass were tested and did not improve the fit over BW",
        "(Section 3.3). Cohort median 36.6 kg (range 6.6-102.9, IQR",
        "20.1-53.5) per Table 1."
      ),
      source_name = "BW"
    ),
    CRCL = list(
      description = "Absolute (not BSA-normalized) estimated glomerular filtration rate, Schwartz (children) or Cockcroft-Gault (adults)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Source term: 'absolute eGFR (in L/h)' (Nijstad 2021 Section 2.3,",
        "Eq. 2). Stored in raw mL/min under the canonical CRCL column",
        "(raw-mL/min variant, Delattre_2010_amikacin.R /",
        "Chen_2023_nemonoxacin.R precedent) and converted to L/h in model()",
        "by x 60/1000. The paper computed eGFR with the Schwartz equation",
        "for girls < 17 years and boys < 14 years and Cockcroft-Gault",
        "otherwise, from the most recent serum creatinine before the",
        "infusion (maximum 10 days), capped at 140 mL/min/1.73 m^2 and",
        "assumed to rise from 35 mL/min/1.73 m^2 at birth to that cap at",
        "1.5 years of age; the capped value was then expressed as an",
        "absolute eGFR. The paper does not print the de-normalization, so",
        "supply the absolute value (for a BSA-normalized estimate, x BSA /",
        "1.73). A 70 kg subject with an absolute eGFR of 100 mL/min (6 L/h)",
        "has RF = 1. Table 1 reports the BSA-normalized renal function:",
        "median 140 mL/min/1.73 m^2 (range 69.3-140, IQR 123.1-140)."
      ),
      source_name = "eGFR"
    ),
    OCC = list(
      description = "Integer-valued occasion index (one occasion per daily clofarabine dose, 1-4)",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Nijstad 2021 Section 2.2: 'interoccasion variability (IOV) was",
        "implemented similarly as IIV, with each dose and subsequent",
        "sampling defined as a separate occasion.' The conditioning",
        "regimen is four once-daily doses (Day -5 to Day -2), so OCC takes",
        "values 1-4, decomposed in model() into binary indicators that",
        "multiplex the per-occasion IOV etas on CL and V2. Observation",
        "records carry the occasion of the dose that preceded them. OCC = 0",
        "or any value outside 1-4 zeros every indicator and gives the",
        "IIV-only parameters."
      ),
      source_name = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      notes = "Tested as the body-size descriptor (allometric exponent 1 on clearances and volumes) and not retained: 'BSA and FFM were evaluated as metrics for body size, but did not improve the model fit over BW' (Nijstad 2021 Section 3.3). The clinical dose was BSA-based (120 mg/m^2 cumulative)."
    ),
    FFM = list(
      description = "Fat-free mass",
      units = "kg",
      type = "continuous",
      notes = "Tested as the body-size descriptor (exponents 0.75 / 1) and not retained (Nijstad 2021 Section 3.3)."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened as a covariate (Section 2.3) and not retained. Renal maturation on CL by the Rhodin sigmoid on postmenstrual age (age in weeks + 40, Hill 3.92, TM50 54.2 weeks, both fixed) was tested and 'did not result in a better fit of the model, so maturation was not included in the final model' (Section 3.3). Age still enters the eGFR derivation (Schwartz vs Cockcroft-Gault switch and the birth-to-1.5-year cap ramp)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 81L,
    n_studies = 1L,
    n_pk_samples = 805L,
    age_range = "0.5-37.8 years (80 paediatric patients 0.5-18 years and one adult of 37.8 years)",
    age_median = "11.1 years (IQR 5.5-14.8)",
    weight_range = "6.6-102.9 kg",
    weight_median = "36.6 kg (IQR 20.1-53.5)",
    sex_female_pct = 37,
    race_ethnicity = "Not reported in source",
    disease_state = paste(
      "Children receiving myeloablative clofarabine-fludarabine-busulfan",
      "conditioning before allogeneic HCT for ALL (49%), AML (35%),",
      "myelodysplastic syndrome (10%), CML (2%) or other (4%); cord blood",
      "(57%) or bone marrow (43%) grafts. Five patients were younger than",
      "12 months."
    ),
    renal_function = "BSA-normalized eGFR (Schwartz / Cockcroft-Gault, capped at 140) median 140 mL/min/1.73 m^2 (range 69.3-140, IQR 123.1-140); no moderate or severe renal impairment",
    dose_range = "Cumulative 120 mg/m^2 clofarabine given as four once-daily 1-h IV infusions (30 mg/m^2/day, Day -5 to Day -2 before HCT), each directly followed by a 1-h fludarabine and a 3-h busulfan infusion",
    concomitant_medications = "Fludarabine (40 mg/m^2 cumulative), busulfan (targeted to a cumulative AUC of 90 mg*h/L), and rabbit ATG in the unrelated-donor setting, with clemastine, paracetamol and prednisolone premedication",
    regions = "Single centre: University Medical Centre Utrecht / Princess Maxima Center for Pediatric Oncology, the Netherlands",
    sampling_window = "Samples drawn during busulfan TDM, mostly 5-8 h after the end of the clofarabine infusion on day 1 or 2 and day 4; a subset to 24 h; from 2016 an additional sample about 1.5 h after the end of the infusion. Median 10 samples per patient (range 3-20). None below the LLOQ of 1 ng/mL.",
    notes = "Retrospective analysis, October 2011 - January 2019 (protocol UMCU 11/063). NONMEM 7.3.0 with FOCE-I; parameter precision by sampling importance resampling."
  )

  ini({
    # Structural parameters: Nijstad 2021 Table 2, for a 70 kg subject. The
    # Table 2 formulas are
    #   CL = (CL_nonrenal + CL_renal x RF) x (BW/70)^0.75
    #   RF = (eGFR (L/h) x 70/BW) / eGFR_STD,  eGFR_STD = 6 L/h
    #   V1 = V1_70kg x (BW/70)^1 ;  V2 = V2_70kg x (BW/70)^1
    #   Q  = Q_70kg  x (BW/70)^0.75
    lcl_nonren <- log(24.0) ; label("Non-renal clearance at 70 kg (L/h)")                                  # Table 2 'CL non-renal,70 kg (L/h)' = 24.0 (95% CI 13.7-34.4)
    lcl_renal  <- log(29.8) ; label("Renal clearance at 70 kg and relative renal function RF = 1 (L/h)")    # Table 2 'CL renal,70 kg (L/h)' = 29.8 (95% CI 23.9-36.1)
    lvc        <- log(268)  ; label("Central volume of distribution V1 at 70 kg (L)")                      # Table 2 'V1 70 kg (L)' = 268 (95% CI 234.8-296.6)
    lvp        <- log(186)  ; label("Peripheral volume of distribution V2 at 70 kg (L)")                   # Table 2 'V2 70 kg (L)' = 186 (95% CI 165.4-210.7)
    lq         <- log(33.2) ; label("Intercompartmental clearance Q at 70 kg (L/h)")                       # Table 2 'Q 70 kg (L/h)' = 33.2 (95% CI 27.5-40.9)

    # Inter-individual variability, exponential (Section 2.2 Eq. 1:
    # Pi = Ppop x exp(eta_i)). Table 2 prints percentages; they are read as
    # CV% and converted with omega^2 = log(CV^2 + 1).
    etalcl ~ 0.031192  # Table 2 'IIV CL (%)' = 17.8 (95% CI 14.6-22.4)
    etalvc ~ 0.015751  # Table 2 'IIV V1 (%)' = 12.6 (95% CI 6.8-18.1)
    etalq  ~ 0.347854  # Table 2 'IIV Q (%)' = 64.5 (95% CI 49.5-83.7)

    # Inter-occasion variability, one occasion per daily dose (Section 2.2).
    # A single IOV variance per parameter is shared by the four occasions:
    # occasion 1 carries the estimate and occasions 2-4 are fixed equal to it
    # (the equivalent of NONMEM $OMEGA BLOCK(1) SAME).
    etaiov_cl_1 ~ 0.009365         # Table 2 'IOV CL (%)' = 9.7 (95% CI 7.8-11.5)
    etaiov_cl_2 ~ fixed(0.009365)  # equal to the occasion-1 IOV variance on CL
    etaiov_cl_3 ~ fixed(0.009365)  # equal to the occasion-1 IOV variance on CL
    etaiov_cl_4 ~ fixed(0.009365)  # equal to the occasion-1 IOV variance on CL
    etaiov_vp_1 ~ 0.142264         # Table 2 'IOV V2 (%)' = 39.1 (95% CI 29.2-53.7)
    etaiov_vp_2 ~ fixed(0.142264)  # equal to the occasion-1 IOV variance on V2
    etaiov_vp_3 ~ fixed(0.142264)  # equal to the occasion-1 IOV variance on V2
    etaiov_vp_4 ~ fixed(0.142264)  # equal to the occasion-1 IOV variance on V2

    # Residual error: proportional (Section 2.2 and Table 2).
    propSd <- 0.083 ; label("Proportional residual error (fraction)")  # Table 2 'Proportional residual error (%)' = 8.3 (95% CI 7.7-8.8)
  })

  model({
    # Occasion indicators for the IOV multiplexing (one occasion per dose).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 + oc4 * etaiov_cl_4
    iov_vp <- oc1 * etaiov_vp_1 + oc2 * etaiov_vp_2 + oc3 * etaiov_vp_3 + oc4 * etaiov_vp_4

    # Relative renal function (Section 2.3 Eq. 2 and Table 2): the absolute
    # eGFR in L/h, standardized to 70 kg, divided by eGFR_STD = 6 L/h
    # (100 mL/min).
    egfr_lh <- CRCL * 60 / 1000
    rf <- egfr_lh * (70 / WT) / 6

    # Individual parameters. Allometric exponents 0.75 (CL, Q) and 1 (V1, V2)
    # were fixed a priori (Section 3.2).
    cl <- (exp(lcl_nonren) + exp(lcl_renal) * rf) * (WT / 70)^0.75 * exp(etalcl + iov_cl)
    vc <- exp(lvc + etalvc) * (WT / 70)
    vp <- exp(lvp + iov_vp) * (WT / 70)
    q <- exp(lq + etalq) * (WT / 70)^0.75

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Clofarabine is given as a 1-h IV infusion into the central compartment.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg over volume in L gives mg/L (x 1000 for the ng/mL of the
    # paper's figures).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
