deJager_2020_factorviii <- function() {
  description <- "One-compartment population PK model for factor VIII (FVIII) activity after perioperative Haemate P / Humate P (plasma-derived VWF/FVIII concentrate) in patients with von Willebrand disease (de Jager 2020). V and CL are allometrically scaled to 70 kg (fixed exponents 1 and 0.75); typical V is 3.28 L and CL 0.037 L/h. CL falls with longer surgery (power exponent -0.416, reference 81 min), with higher time-varying VWF activity (exponent -0.263, reference 1.65 IU/mL) and is 44.5% lower in ASA class III or IV than class II. The measured lowest-ever endogenous FVIII level is added to the prediction. A fitted virtual dose with bioavailability 0.200 (inter-occasion variability 154.5% CV) absorbs pre-dose FVIII levels above that baseline. IIV on V (30.9% CV) and CL (84.1% CV) is correlated (r = 0.47); residual error is proportional (20.3%)."
  reference <- "de Jager NCB, Bukkems LH, Heijdra JM, Hazendonk HCAM, Fijnvandraat K, Meijer K, Eikenboom J, Laros-van Gorkom BAP, Leebeek FWG, Cnossen MH, Mathot RAA; OPTI-CLOT group. One piece of the puzzle: Population pharmacokinetics of FVIII during perioperative Haemate P/Humate P treatment in von Willebrand disease patients. J Thromb Haemost. 2020;18(2):295-305. doi:10.1111/jth.14652. PMID:31557387. PMCID:PMC7028136."
  vignette <- "deJager_2020_factorviii"
  units <- list(time = "h", dosing = "IU", concentration = "IU/mL")

  # Bookkeeping state for the fitted virtual dose (Methods 2.2): a fixed dose
  # given before the pre-dose FVIII sample that represents the patient's
  # endogenous FVIII above the lowest-ever baseline. It decays with the same
  # elimination rate as central so the two states sum to the measured FVIII.
  paper_specific_compartments <- c("virtual")

  compartmentData <- list(
    central = list(analyte = "factor viii", units = "IU", specimen = "plasma", verified = TRUE),
    virtual = list(analyte = "factor viii", units = "IU", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric power scaling at reference 70 kg with exponents fixed a priori at 1 on V and 0.75 on CL (Methods 2.2; Supplement 1 Eq. 4-5; Table 2 footnote). Cohort median 77.0 kg (range 8.8-118.0; Table 1).",
      source_name = "Weight"
    ),
    T_SURG = list(
      description = "Duration of the surgical procedure (operative time)",
      units = "h",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL, exponent -0.416, centred at 81 min (Table 2 footnote: '(Surgery duration_i / 81)^-0.416'). The paper records duration in minutes; the canonical column is in hours, so the model centres on 81 / 60 = 1.35 h, which leaves the effect unchanged. Time-fixed per surgical occasion (a patient with several surgeries carries a different value on each). Cohort median 65 min (range 7-470; Table 1).",
      source_name = "Surgery duration"
    ),
    VWF = list(
      description = "Plasma von Willebrand factor activity (VWF:Act), time-varying",
      units = "IU/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL, exponent -0.263, centred at 1.65 IU/mL (Table 2 footnote: '(VWF:Act_ij / 1.65)^-0.263'); stored in the canonical IU/dL, so the reference is 165 IU/dL. Measured VWF activity after concentrate dosing, carried forward from each measurement to the next (last observation carried forward; Methods 2.3). Interquartile range of all measured levels 0.78-2.21 IU/mL (Results 3.2). Patient baseline (pre-treatment) VWF:Act median 0.15 IU/mL (range 0.0-0.58; Table 1), which is not the reference value.",
      source_name = "VWF:Act"
    ),
    ASA_GE3 = list(
      description = "American Society of Anesthesiologists physical status class III or higher (1 = ASA III or IV, 0 = ASA II)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ASA class II)",
      notes = "Multiplicative factor 0.555 on CL (Table 2 footnote: '0.555^ASAclass3,4'); a 44.5% lower CL in ASA III/IV (Results 3.2: 'a 44% decrease'). The cohort had no ASA class I or V; 126 ASA II, 25 ASA III and 1 ASA IV classifications, with 11.6% missing (Table 1). Time-fixed per surgical occasion.",
      source_name = "ASA class III or IV"
    ),
    FVIII_BL = list(
      description = "Lowest-ever measured endogenous plasma FVIII:C activity of the patient",
      units = "IU/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Added to the predicted exogenous FVIII in the observation (Supplement 1 Eq. 2: 'Y_obs = (baseline FVIII:C_i + Y_pred) * (1 + eps)'). Defined as the lowest FVIII level ever measured in the patient (Methods 2.1), not the pre-dose level, which was often higher (Methods 2.2). Cohort median 0.41 IU/mL (range 0.01-0.97; Table 1). Time-fixed per patient. Set to 0 to simulate drug-attributable FVIII only.",
      source_name = "Baseline FVIII level"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 117L,
    n_studies = 1L,
    n_surgeries = 172L,
    n_observations = 892L,
    age_range = "0.5-82 years (seven children aged 0.5-16 years)",
    age_median = "50 years",
    weight_range = "8.8-118.0 kg",
    weight_median = "77.0 kg",
    sex_female_pct = 67,
    race_ethnicity = "not reported; Netherlands cohort",
    disease_state = "von Willebrand disease (VWD type 1 47%, 2A 27%, 2B 6%, 2M 11%, 2N 5%, 3 5% of occasions) undergoing minor (28%) or major (72%) surgery; ASA class II 83%, III 16%, IV 1%",
    dose_range = "Perioperative Haemate P / Humate P (VWF:FVIII ratio 2.4:1) bolus doses of 5.5-66.1 IU FVIII/kg (median 20.8 IU/kg, median 5 doses per occasion); 4.7% of doses as continuous infusion at 0.19-4.2 IU/h/kg",
    regions = "Netherlands (five academic haemophilia treatment centres, OPTI-CLOT; 2000-2018)",
    notes = "Final model re-estimated on the index data (97 patients, 141 surgeries, 684 FVIII samples) plus the external-validation data (20 patients, 31 surgeries, 208 samples; Erasmus MC Rotterdam). Demographics from Table 1, 'All available data' column. FVIII measured by one-stage clotting assay, LLOQ 0.01 IU/mL; samples taken from 146 h before to 524 h after surgery. Doses and levels after desmopressin were excluded."
  )

  ini({
    # Structural PK parameters at the 70 kg reference (Table 2, 'Final FVIII
    # Estimate' column).
    lvc <- log(3.28); label("Volume of distribution at 70 kg (L)") # Table 2 final: Volume of distribution 3.28 L/70 kg (RSE 4.4%)
    lcl <- log(0.037); label("Clearance at 70 kg and reference covariates (L/h)") # Table 2 final: Clearance 0.037 L/70 kg/h (RSE 9.8%)

    # Bioavailability of the virtual dose (Supplement 1 Eq. 3: F = theta_F * exp(eta_F)).
    lfvirtual <- log(0.200); label("Bioavailability of the virtual pre-dose correction dose (fraction)") # Table 2 final: Bioavailability virtual dose 0.200 (RSE 21.5%)

    # Allometric exponents set a priori (Methods 2.2; Supplement 1 Eq. 4-5).
    e_wt_vc <- fixed(1); label("Allometric exponent of (WT/70) on V (unitless)") # Supplement 1 Eq. 4; Table 2 footnote '(Weight_i/70)^1'
    e_wt_cl <- fixed(0.75); label("Allometric exponent of (WT/70) on CL (unitless)") # Supplement 1 Eq. 5; Table 2 footnote '(Weight_i/70)^0.75'

    # Covariate effects on CL (Table 2 final; equation form in Table 2 footnote).
    e_tsurg_cl <- -0.416; label("Power exponent of surgery duration (relative to 81 min) on CL (unitless)") # Table 2 final: Surgery duration on CL -0.416 (RSE 24.6%)
    e_asa_ge3_cl <- 0.555; label("Multiplicative factor on CL for ASA class III or IV (unitless)") # Table 2 final: ASA class III or IV on CL 0.555 (RSE 15.4%)
    e_vwf_cl <- -0.263; label("Power exponent of VWF activity (relative to 1.65 IU/mL) on CL (unitless)") # Table 2 final: VWF activity on CL -0.263 (RSE 43.2%)

    # IIV, exponential. Table 2 prints %CV = sqrt(exp(omega^2) - 1) * 100, so
    # omega^2 = log(1 + CV^2); covariance = r * sqrt(omega2_CL * omega2_V).
    etalcl + etalvc ~ c(
      0.534902,
      0.103805, 0.0911935
    ) # Table 2 final: IIV CL 84.1 %CV, IIV V 30.9 %CV, correlation V-CL 47 %

    # Inter-occasion variability on the virtual-dose bioavailability, carried
    # as one eta: one draw per simulated subject, i.e. per surgical occasion.
    etaiov_lfvirtual_1 ~ 1.21995 # Table 2 final: IOV Bioavailability virtual dose 154.5 %CV; log(1 + 1.545^2)

    # Proportional residual error on the linear scale (Supplement 1 Eq. 2).
    propSd <- 0.203; label("Proportional residual error (fraction)") # Table 2 final: Proportional residual variability 20.3 % (RSE 7.6%)
  })
  model({
    # Table 2 footnote:
    #   V  = theta_V  * (Weight_i/70)^1 * exp(eta_V)
    #   CL = theta_CL * (Weight_i/70)^0.75 * (Surgery duration_i/81)^-0.416 *
    #        (VWF:Act_ij/1.65)^-0.263 * 0.555^ASAclass3,4 * exp(eta_CL)
    # Surgery duration is centred at 81 min = 1.35 h and VWF:Act at
    # 1.65 IU/mL = 165 IU/dL, matching the canonical column units.
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    cl <- exp(lcl + etalcl) *
      (WT / 70)^e_wt_cl *
      (T_SURG / 1.35)^e_tsurg_cl *
      (VWF / 165)^e_vwf_cl *
      e_asa_ge3_cl^ASA_GE3

    kel <- cl / vc

    # Real Haemate P doses (bolus or infusion) enter central with F = 1. The
    # fitted virtual dose is given to `virtual` before the first pre-dose
    # sample; the paper does not report its amount.
    d/dt(central) <- -kel * central
    d/dt(virtual) <- -kel * virtual
    f(virtual) <- exp(lfvirtual + etaiov_lfvirtual_1)

    # Supplement 1 Eq. 2: observed FVIII = (baseline FVIII + predicted) * (1 + eps).
    # Dose in IU and vc in L give IU/L; divide by 1000 for IU/mL.
    Cc <- (central + virtual) / vc / 1000 + FVIII_BL
    Cc ~ prop(propSd)
  })
}
