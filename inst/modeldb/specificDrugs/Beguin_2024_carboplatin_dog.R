Beguin_2024_carboplatin_dog <- function() {
  description <- paste(
    "Veterinary (client-owned dog).",
    "One-compartment population PK model with linear elimination for",
    "free (ultrafilterable) carboplatin in 16 client-owned dogs with",
    "solid tumours, fitted to 39 concentration-time profiles after a",
    "300 mg/m2 dose given as an approximately 20 min intravenous",
    "infusion. All disposition parameters are body-weight-normalised:",
    "the authors divided each dog's total dose by its own body weight",
    "before fitting, so clearance is L/h/kg, volume is L/kg, and the",
    "dosing amount supplied to this model is per kg of body weight.",
    "Clearance carries two covariates (Beguin 2024 Table 1): an",
    "uncentred power effect of plasma creatinine, CREAT^-0.25, which",
    "is how beta*ln(CREAT) on the log scale re-expresses, and a",
    "multiplicative effect of neutering worth -22% on the log scale.",
    "Inter-individual and inter-occasion variability on volume were",
    "both FIXED at 0.1 by the authors because the data were too sparse",
    "to estimate them. Residual error is the Monolix combined form",
    "sqrt(a^2 + (b*C)^2), which is rxode2's default add() + prop().",
    "The cumulative auc_central state reproduces the AUC output the",
    "authors added by hand to their Monolix model; it is the exposure",
    "driver of the companion model",
    "Beguin_2024_carboplatin_thrombocytopenia_dog."
  )
  reference <- paste(
    "Beguin J, Mahfoudhi S, Uzel M, Rostang A, Ibish C, Ferran AA,",
    "Pelligand L, Hulin A, Kohlhauer M. (2024).",
    "Population pharmacokinetics modelling for clinical dose",
    "adjustment of carboplatin in dogs.",
    "BMC Veterinary Research 20:575.",
    "doi:10.1186/s12917-024-04404-1.",
    sep = " "
  )
  vignette <- "Beguin_2024_carboplatin_dog"

  # Bookkeeping state that integrates Cc so the exposure-response layer can
  # read AUC off the solve, reproducing the AUC output the authors added by
  # hand to their Monolix model. Same idiom as auc_central in
  # Assmus_2025_benznidazole_qpcr.R; not a biological compartment.
  paper_specific_compartments <- c("auc_central")


  # Concentrations are ug/L and amounts are ug per kg of body weight. The
  # paper never states the modelling units directly and its printed
  # "AUC0-inf (mg.h/L)" label is wrong by a factor of 1000; see the
  # vignette Errata for the two independent numerical anchors (the
  # additive residual constant and Dose/CL) that pin the ug/L scale.
  units <- list(
    time          = "h",
    dosing        = "ug",
    concentration = "ug/L"
  )

  compartmentData <- list(
    central     = list(analyte = "carboplatin", units = "ug", specimen = "plasma", verified = TRUE),
    auc_central = list(analyte = "carboplatin", units = "ug*h/L", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    CREAT = list(
      description        = "Plasma creatinine concentration on the day of carboplatin administration",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters UNCENTRED. Beguin 2024 Eq. (3) adds beta * ln(CREAT)",
        "on the log-parameter scale, which is identical to the power",
        "form CREAT^-0.25 used in model(). There is no reference",
        "value: the typical clearance of 6.9 L/h/kg is therefore the",
        "value at CREAT = 1 mg/L, which is far below the observed",
        "range and is an extrapolation, not a typical dog. Cohort",
        "mean 7.32 +/- 1.86 mg/L; dogs with an admission creatinine",
        "above 14 mg/L were excluded. Note the unit is mg/L, not the",
        "mg/dL or umol/L more common in human papers - divide by 10",
        "to obtain mg/dL."
      ),
      source_name        = "plasma creatinine concentration"
    ),
    NEUTERED = list(
      description        = "Neutering status (1 = surgically neutered, i.e. spayed female or castrated male; 0 = sexually intact)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (sexually intact)",
      notes              = paste(
        "Beguin 2024 calls this 'sterilization status' and codes it",
        "0/1 per Eq. (3). Orthogonal to SEXF: the cohort held intact",
        "and neutered animals of both sexes (4 neutered females and 5",
        "neutered males out of 16 dogs), and sex was screened",
        "separately and NOT retained. The authors caution that the",
        "effect is probably confounded with age and creatinine rather",
        "than causal, since neutering is not expected to alter GFR."
      ),
      source_name        = "neutered status"
    ),
    OCC = list(
      description        = "Occasion index; one occasion per carboplatin administration in the same dog",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Beguin 2024 analysed repeat administrations in the same dog",
        "as separate occasions: 4 dogs received 5 administrations, 2",
        "received 3, 3 received 2 and 7 received 1, giving 39",
        "profiles from 16 dogs. Five occasions are encoded here, the",
        "maximum any dog received. Monolix shares a single IOV",
        "variance across occasions; the per-occasion etaiov_* slots",
        "below reproduce that by fixing occasions 2-5 to the occasion",
        "1 value."
      ),
      source_name        = "occasion"
    )
  )

  # Screened by Beguin 2024 and NOT retained in the final model, so they are
  # documented rather than referenced in model(). Body condition score and
  # inclusion centre were screened too but are omitted here because neither
  # has a canonical register entry and neither carries a reported estimate.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Screened as a continuous covariate and NOT retained",
        "(p = 0.32 for V, p = 0.86 for Cl). Notable because the dose",
        "was prescribed on body surface area, BSA(m2) = 0.1 *",
        "WT(kg)^(2/3) per Eq. (1); the authors argue from this null",
        "result that carboplatin should be dosed per kg rather than",
        "per m2. Cohort mean 21.5 +/- 7.8 kg."
      )
    ),
    AGE = list(
      description = "Age at inclusion",
      units       = "years",
      type        = "continuous",
      notes       = "Screened and NOT retained. Cohort mean 11.1 +/- 1.72 years."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened as a categorical covariate and NOT retained; the",
        "separate NEUTERED indicator was retained instead. Cohort as",
        "printed: 14 females and 12 males, which sums to 26 rather",
        "than the 16 dogs stated elsewhere in the same paragraph -",
        "see the vignette Errata."
      )
    )
  )

  population <- list(
    species      = "dog (client-owned, mixed breeds)",
    n_subjects   = 16L,
    n_studies    = 1L,
    n_occasions  = 39L,
    age_range    = "mean 11.1 +/- 1.72 years",
    weight_range = "mean 21.5 +/- 7.8 kg",
    disease_state = paste(
      "Solid tumours treated with adjuvant carboplatin after surgical",
      "resection: carcinoma (4), osteosarcoma (3), melanoma (3),",
      "mammary adenocarcinoma (2), ovarian dysgerminoma (2),",
      "chondrosarcoma (1), fibrosarcoma (1)."
    ),
    dose_range   = "300 mg/m2 (mean 300.4 +/- 7.6 mg/m2, equal to 10.7 +/- 1.0 mg/kg) as an approximately 20 min IV infusion",
    regions      = "France (National Veterinary School of Alfort and Oniris VetAgroBio)",
    renal_function = "Mean plasma creatinine 7.32 +/- 1.86 mg/L; dogs above 14 mg/L at admission were excluded",
    notes        = paste(
      "Prospective two-centre clinical trial, June 2019 to January",
      "2021. 27 dogs enrolled, 11 excluded for protocol deviations,",
      "16 analysed contributing 39 administrations. Exclusions also",
      "covered neutrophils below 1500/uL or thrombocytes below",
      "50 000/uL at admission. Three to four plasma samples per",
      "administration across four sampling windows (0-1 h, 1-2 h,",
      "2-4 h, 4-12 h). Free carboplatin was measured in plasma",
      "ultrafiltrate by LC-MS/MS. Values below the limit of",
      "quantification were treated as interval-censored. See",
      "Beguin 2024 Results, 'Study dogs'."
    )
  )

  ini({
    # --- Structural disposition, body-weight-normalised (Beguin 2024 Table 1).
    lvc <- log(4.05)
    label("Volume of distribution (L/kg)")                                  # Table 1 V = 4.05 L/kg [RSE 6.04%; bootstrap 4.02 (3.65-4.52)]
    lcl <- log(6.9)
    label("Clearance (L/h/kg)")                                             # Table 1 Cl = 6.9 L/h/kg [RSE 24.7%; bootstrap 6.82 (3.94-10.35)]

    # --- Covariate effects on clearance (Beguin 2024 Table 1, Eq. 3).
    e_creat_cl <- -0.25
    label("Exponent of plasma creatinine on clearance (unitless)")          # Table 1 beta_creatinine = -0.25 [RSE 50.3%; bootstrap -0.24 (-0.45, 0.07)]
    e_neutered_cl <- -0.22
    label("Log-scale effect of neutering on clearance (unitless)")          # Table 1 beta_neutered = -0.22 [RSE 34.2%; bootstrap -0.22 (-0.44, -0.07)]

    # --- Inter-individual variability. Table 1 states twice that the random
    # effects are reported as STANDARD DEVIATIONS, so each is squared here.
    etalvc ~ fixed(0.01)                                                    # Table 1 omega_V = 0.1 SD, held constant by the authors per footnote a; 0.1^2
    etalcl ~ 0.005476                                                       # Table 1 omega_Cl = 0.074 SD [RSE 53.1%]; 0.074^2

    # --- Inter-occasion variability over five occasions. Monolix carries one
    # shared variance per parameter across occasions, so occasions 2-5 are
    # fixed to the occasion-1 value rather than estimated separately.
    etaiov_vc_1 ~ fixed(0.01)                                               # Table 1 gamma_V = 0.1 SD, held constant by the authors per footnote a; 0.1^2
    etaiov_vc_2 ~ fixed(0.01)                                               # shared with occasion 1
    etaiov_vc_3 ~ fixed(0.01)                                               # shared with occasion 1
    etaiov_vc_4 ~ fixed(0.01)                                               # shared with occasion 1
    etaiov_vc_5 ~ fixed(0.01)                                               # shared with occasion 1
    etaiov_cl_1 ~ 0.0121                                                    # Table 1 gamma_Cl = 0.11 SD [RSE 28.3%]; 0.11^2
    etaiov_cl_2 ~ fixed(0.0121)                                             # shared with occasion 1
    etaiov_cl_3 ~ fixed(0.0121)                                             # shared with occasion 1
    etaiov_cl_4 ~ fixed(0.0121)                                             # shared with occasion 1
    etaiov_cl_5 ~ fixed(0.0121)                                             # shared with occasion 1

    # --- Residual error. Beguin 2024 Eq. (5) is
    # Obs = C + sqrt(a^2 + (b*C)^2) * eps, the Monolix "combined2" form,
    # which is exactly rxode2's default add() + prop() combination.
    addSd <- 12.28
    label("Additive residual error (ug/L)")                                 # Table 1 a = 12.28 [RSE 17.2%; bootstrap 11.72 (7.3-19.9)]
    propSd <- 0.23
    label("Proportional residual error (fraction)")                         # Table 1 b = 0.23 [RSE 12.2%; bootstrap 0.22 (0.17-0.25)]
  })

  model({
    # 1. Occasion indicators (Beguin 2024 Eq. 2: log(theta_i) carries both
    #    eta_i and eta_occ).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)

    iov_vc <- oc1 * etaiov_vc_1 + oc2 * etaiov_vc_2 + oc3 * etaiov_vc_3 +
      oc4 * etaiov_vc_4 + oc5 * etaiov_vc_5
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 +
      oc4 * etaiov_cl_4 + oc5 * etaiov_cl_5

    # 2. Individual parameters. Eq. (3) is uncentred, so beta * ln(CREAT)
    #    on the log scale becomes the power form CREAT^e_creat_cl here.
    vc <- exp(lvc + etalvc + iov_vc)
    cl <- exp(lcl + etalcl + iov_cl) * CREAT^e_creat_cl *
      exp(e_neutered_cl * NEUTERED)

    kel <- cl / vc

    # 3. One-compartment disposition with linear elimination; dosing is by
    #    IV infusion with no lag time.
    d/dt(central) <- -kel * central

    Cc <- central / vc

    # 4. Cumulative AUC. The authors state that "Area Under the Curve
    #    (AUC0-last) computation was manually added as an additional output
    #    to the Monolix model"; this state reproduces that output and is the
    #    exposure driver of the companion thrombocytopenia model.
    d/dt(auc_central) <- Cc

    Cc ~ add(addSd) + prop(propSd)
  })
}
