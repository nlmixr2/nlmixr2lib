Balakirouchenane_2020_dabrafenib <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model for oral dabrafenib and its",
    "active metabolite hydroxy-dabrafenib in a real-life cohort of adults with",
    "BRAF V600-mutated solid tumours, mostly metastatic melanoma (Balakirouchenane",
    "2020). Dabrafenib is two-compartment with first-order absorption (rate",
    "constant fixed to a literature value) and an absorption lag time, and is",
    "eliminated exclusively by irreversible conversion to hydroxy-dabrafenib,",
    "which is itself two-compartment with first-order elimination. All",
    "disposition parameters are apparent oral values. Dabrafenib apparent",
    "clearance decreases linearly with age (median-centred at 61.2 years) and is",
    "17% lower in women; hydroxy-dabrafenib apparent clearance decreases",
    "linearly with age. Inter-individual variability on dabrafenib CL/F and",
    "V2/F and on hydroxy-dabrafenib CLm/F and V3/F, inter-occasion variability",
    "on dabrafenib CL/F, and proportional residual error on both analytes.",
    sep = " "
  )
  reference <- paste(
    "Balakirouchenane D, Guegan S, Csajka C, Jouinot A, Heidelberger V,",
    "Puszkiel A, Zehou O, Khoudour N, Courlet P, Kramkimel N, Lheure C,",
    "Franck N, Huillard O, Arrondeau J, Vidal M, Goldwasser F, Maubec E,",
    "Dupin N, Aractingi S, Guidi M, Blanchet B.",
    "Population Pharmacokinetics/Pharmacodynamics of Dabrafenib Plus",
    "Trametinib in Patients with BRAF-Mutated Metastatic Melanoma.",
    "Cancers. 2020;12(4):931. doi:10.3390/cancers12040931.",
    "Parameter estimates from Table 2 (final model column); the covariate",
    "equations from the Table 2 footnote; population from Table 1.",
    sep = " "
  )
  vignette <- "Balakirouchenane_2020_dabrafenib_trametinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(
      analyte = "dabrafenib",
      units = "mg (apparent, i.e. amount/F)",
      specimen = "administration site",
      verified = FALSE
    ),
    central = list(
      analyte = "dabrafenib",
      units = "mg (apparent, i.e. amount/F)",
      specimen = "plasma",
      verified = FALSE
    ),
    peripheral1 = list(
      analyte = "dabrafenib",
      units = "mg (apparent, i.e. amount/F)",
      specimen = "plasma",
      verified = FALSE
    ),
    central_ohd = list(
      analyte = "hydroxy-dabrafenib",
      units = "mg dabrafenib-equivalents (apparent, i.e. amount/F)",
      specimen = "plasma",
      verified = FALSE
    ),
    peripheral1_ohd = list(
      analyte = "hydroxy-dabrafenib",
      units = "mg dabrafenib-equivalents (apparent, i.e. amount/F)",
      specimen = "plasma",
      verified = FALSE
    )
  )

  covariateData <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Median-centred at MAGE = 61.2 years (Table 2 footnote; the Table 1",
        "cohort median) and applied as a linear deviation on both dabrafenib",
        "CL/F and hydroxy-dabrafenib CLm/F, e.g.",
        "cl *= (1 + e_age_cl * (AGE - 61.2) / 61.2) with e_age_cl = -0.536.",
        "Table 2 prints the two coefficients without a sign (0.536 and 0.589);",
        "they are encoded NEGATIVE because every statement in the paper has",
        "clearance falling with age: Results 2.2.1 ('CL/F is reduced ... by 55%",
        "when comparing 20-year-old to 90-year-old patients' and 'A similar",
        "decrease (51%) is observed in CLm/F'), the Figure 3 simulated composite",
        "AUC (10,447 ng.h/mL at age 20 vs 19,542 at age 90 in men), and the",
        "Discussion's recommendation of a reduced starting dose in the elderly.",
        "With the negative sign the typical CL/F ratio age 90 / age 20 is",
        "0.748 / 1.361 = 0.55 and the CLm/F ratio is 0.723 / 1.396 = 0.52,",
        "matching the paper's 55% and 51%; the vignette reproduces all four",
        "Figure 3 composite AUC medians to within 1%. Cohort range 20-90 years.",
        sep = " "
      ),
      source_name = "AGE"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "1 = female, 0 = male, matching the paper's own coding ('Sex being 0/1",
        "for man/woman', Table 2 footnote). Enters dabrafenib CL/F as the",
        "multiplicative factor 0.832^SEXF (Table 2 'theta sex/(CL/F) = 0.832');",
        "Results 2.2.1: 'CL/F is reduced by 17% in women vs. men'. 30 of 73",
        "patients (41%) were female (Table 1).",
        sep = " "
      ),
      source_name = "Sex"
    ),
    OCC = list(
      description = "Occasion index for the inter-occasion variability on dabrafenib CL/F",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "The paper estimates an IOV on CL/F of 17.4% (Table 2) but does not",
        "define its occasions; in this routine-care cohort each sampling visit",
        "is the natural occasion. Six occasion slots are provided (the cohort",
        "averaged 424 / 73 = 5.8 dabrafenib samples per patient); an OCC value",
        "outside 1-6 receives no IOV. Set OCC = 1 throughout to simulate a",
        "single occasion.",
        sep = " "
      ),
      source_name = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = "Tested (Methods 4.4.2) but not retained. Median 73.0 (51.7-166.0) kg (Table 1)."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Tested but not retained. Median 25.9 (17.4-44.6) kg/m^2 (Table 1)."
    ),
    FFM = list(
      description = "Fat-free mass (Janmahasatian)",
      units = "kg",
      type = "continuous",
      notes = "Tested but not retained. Median 53.7 (34.2-94.4) kg (Table 1)."
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      notes = "Tested but not retained. Median 1.9 (1.0-3.0) m^2 (Table 1)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Tested but not retained. Median 30.0 (12.0-103.0) U/L (Table 1)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Tested but not retained. Median 28.0 (10.0-116) U/L (Table 1)."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Tested but not retained. Median 6.3 (1.4-47.0) umol/L (Table 1)."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Tested but not retained. Median 43.0 (23.0-49.0) g/L (Table 1)."
    ),
    CRP = list(
      description = "C-reactive protein",
      units = "mg/L",
      type = "continuous",
      notes = "Tested but not retained. Median 19.9 (1.0-248.0) mg/L (Table 1)."
    ),
    CONMED_PPI = list(
      description = "Concomitant proton-pump inhibitor",
      units = "(binary)",
      type = "binary",
      notes = "No effect on dabrafenib or hydroxy-dabrafenib PK (dOFV = 3.5, p > 0.05; Results 2.2.1). 57 of 73 patients (78%) took a PPI (Table 1)."
    ),
    CONMED_TRAMETINIB = list(
      description = "Concomitant trametinib (combination therapy) indicator",
      units = "(binary)",
      type = "binary",
      notes = "Dabrafenib monotherapy vs combination with trametinib did not differ (dOFV = 2.9, p > 0.05; Results 2.2.1). 60 of 73 patients (82%) received the combination."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 73L,
    n_studies = 1L,
    n_observations = "424 dabrafenib/hydroxy-dabrafenib concentration records",
    age_range = "20-90 years (median 61.2)",
    weight_range = "51.7-166.0 kg (median 73.0)",
    sex_female_pct = 41,
    disease_state = paste(
      "Adults with BRAF V600-mutated metastatic solid tumours treated in routine",
      "care: metastatic melanoma 65 (89%), anaplastic thyroid carcinoma 6,",
      "non-small-cell lung carcinoma 2 (Table 1).",
      sep = " "
    ),
    dose_range = paste(
      "Dabrafenib 150 mg orally twice daily (reduced starting doses such as",
      "75 mg twice daily allowed at physician discretion); 60 patients (82%)",
      "also received trametinib 2 mg once daily.",
      sep = " "
    ),
    regions = "France (three Assistance Publique-Hopitaux de Paris hospitals: Cochin, Henri Mondor, Avicenne)",
    notes = paste(
      "Observational multicentre 'real-life' cohort enrolled July 2015 to June",
      "2017; sparse samples drawn during routine visits with exact sampling",
      "times recorded; all subjects sampled at steady state. HPLC-MS/MS assay,",
      "calibration range 10-2000 ng/mL for both analytes; BLQ samples excluded.",
      "Fit in NONMEM 7.4.1 (FOCEI) with PsN 4.8.0 in molar units; residual",
      "errors of the two analytes were correlated (L2 item). 500-replicate",
      "bootstrap and pcVPC for internal validation; no external validation.",
      sep = " "
    )
  )

  ini({
    # ---- Dabrafenib absorption (Table 2, final model) ----
    lka <- fixed(log(1.8)); label("Dabrafenib first-order absorption rate constant (1/h)") # Table 2 'ka (1/h) 1.8 fixed'; Methods 4.4.1: fixed to a previously published value (reference 18) to allow estimation of the other parameters
    ltlag <- log(0.499); label("Dabrafenib absorption lag time (h)") # Table 2 'Tlag (h) 0.499' (RSE 0.2%); bootstrap median 0.499

    # ---- Dabrafenib disposition (apparent values) ----
    # CL/F is the typical value for a MAN aged 61.2 years. Dabrafenib is
    # eliminated only by conversion to hydroxy-dabrafenib (Methods 4.4.1), so
    # CL/F is also the formation clearance of the metabolite.
    lcl <- log(19.3); label("Dabrafenib apparent clearance CL/F, man aged 61.2 years (L/h)") # Table 2 'CL/F (L/h) 19.3' (RSE 7.5%); bootstrap median 19.2
    lvc <- log(39.1); label("Dabrafenib apparent central volume V2/F (L)") # Table 2 'V2/F (L) 39.1' (RSE 13.7%); bootstrap median 38.3
    lq <- log(3.40); label("Dabrafenib apparent intercompartmental clearance Q/F (L/h)") # Table 2 'Q/F (L/h) 3.40' (RSE 21.5%); bootstrap median 3.32
    lvp <- log(18.7); label("Dabrafenib apparent peripheral volume V4/F (L)") # Table 2 'V4/F (L) 18.7' (RSE 20.1%); bootstrap median 18.7

    # ---- Hydroxy-dabrafenib disposition (apparent values) ----
    lcl_ohd <- log(23.2); label("Hydroxy-dabrafenib apparent clearance CLm/F at age 61.2 years (L/h)") # Table 2 'CLm/F (L/h) 23.2' (RSE 5.9%); bootstrap median 22.9
    lvc_ohd <- log(5.11); label("Hydroxy-dabrafenib apparent central volume V3/F (L)") # Table 2 'V3/F (L) 5.11' (RSE 30.5%); bootstrap median 4.99
    lq_ohd <- log(7.21); label("Hydroxy-dabrafenib apparent intercompartmental clearance Qm/F (L/h)") # Table 2 'Qm/F (L/h) 7.21' (RSE 22.3%); bootstrap median 7.02
    lvp_ohd <- log(27.1); label("Hydroxy-dabrafenib apparent peripheral volume V5/F (L)") # Table 2 'V5/F (L) 27.1' (RSE 23.9%); bootstrap median 24.7

    # ---- Covariate effects (Table 2 and its footnote) ----
    # CL_ind/F  = CL/F  * (1 + theta_age,CL  * (AGE - 61.2)/61.2) * theta_sex^Sex * exp(eta)
    # CLm_ind/F = CLm/F * (1 + theta_age,CLm * (AGE - 61.2)/61.2) * exp(eta)
    # Table 2 prints the age coefficients unsigned; they are negative, as the
    # paper's own 20-vs-90-year clearance decreases (55% and 51%) and its
    # Figure 3 AUCs require (see covariateData$AGE$notes).
    e_age_cl <- -0.536; label("Linear age effect on dabrafenib CL/F per unit relative deviation from 61.2 years (unitless)") # Table 2 'theta age/(CL/F) 0.536' (RSE 28.4%), sign from Results 2.2.1
    e_age_cl_ohd <- -0.589; label("Linear age effect on hydroxy-dabrafenib CLm/F per unit relative deviation from 61.2 years (unitless)") # Table 2 'theta age/(CLm/F) 0.589' (RSE 33.6%), sign from Results 2.2.1
    e_sexf_cl <- 0.832; label("Multiplicative factor on dabrafenib CL/F for women vs men (unitless)") # Table 2 'theta sex/(CL/F) 0.832' (RSE 6.4%); bootstrap median 0.829

    # ---- Inter-individual variability (Table 2) ----
    # Reported as CV%; on the exponential-eta scale omega^2 = log(1 + CV^2).
    etalcl ~ 0.025278 # log(1 + 0.160^2); Table 2 'IIV CL/F (%) 16.0'
    etalvc ~ 0.229574 # log(1 + 0.508^2); Table 2 'IIV V2/F (%) 50.8'
    etalcl_ohd ~ 0.056002 # log(1 + 0.240^2); Table 2 'IIV CLm/F (%) 24.0'
    etalvc_ohd ~ 0.203451 # log(1 + 0.475^2); Table 2 'IIV V3/F (%) 47.5'

    # ---- Inter-occasion variability on dabrafenib CL/F (Table 2) ----
    # Occasion-indicator expansion (NONMEM $OMEGA BLOCK(1) SAME style); the
    # first slot carries the estimate and the others repeat it as fixed.
    etaiov_cl_1 ~ 0.029827 # log(1 + 0.174^2); Table 2 'IOV (%) 17.4' on CL/F
    etaiov_cl_2 ~ fixed(0.029827)
    etaiov_cl_3 ~ fixed(0.029827)
    etaiov_cl_4 ~ fixed(0.029827)
    etaiov_cl_5 ~ fixed(0.029827)
    etaiov_cl_6 ~ fixed(0.029827)

    # ---- Residual error (Table 2) ----
    # The two proportional errors were estimated with a correlation of 87.0%
    # (Table 2 'RUV corr (%) 87.0'); correlated residual errors cannot be
    # expressed in nlmixr2, so they are independent here.
    propSd <- 0.487; label("Dabrafenib proportional residual error (fraction)") # Table 2 'RUV of DAB (%) 48.7' (RSE 6.6%)
    propSd_ohd <- 0.531; label("Hydroxy-dabrafenib proportional residual error (fraction)") # Table 2 'RUV of OHD (%) 53.1' (RSE 6.3%)
  })

  model({
    # The analysis was run in molar units (Methods 4.4) with complete
    # conversion of dabrafenib to hydroxy-dabrafenib. Doses are in mg of
    # dabrafenib, so the metabolite states hold dabrafenib-equivalent mg and
    # the hydroxy-dabrafenib concentration is converted to mass units with the
    # molecular-weight ratio. Molecular weights computed from the formulas
    # C23H20F3N5O2S2 (dabrafenib) and C23H20F3N5O3S2 (hydroxy-dabrafenib, one
    # added oxygen); not stated in the paper.
    mwDab <- 519.56 # g/mol
    mwOhd <- 535.56 # g/mol
    mgPerLToNgPerMl <- 1000
    ageMedian <- 61.2 # years; Table 2 footnote 'MAGE = 61.2 years'

    # ---- Inter-occasion indicators ----
    iov_cl <- (OCC == 1) * etaiov_cl_1 + (OCC == 2) * etaiov_cl_2 +
      (OCC == 3) * etaiov_cl_3 + (OCC == 4) * etaiov_cl_4 +
      (OCC == 5) * etaiov_cl_5 + (OCC == 6) * etaiov_cl_6

    # ---- Individual parameters ----
    ka <- exp(lka)
    tlag <- exp(ltlag)
    cl <- exp(lcl + etalcl + iov_cl) *
      (1 + e_age_cl * (AGE - ageMedian) / ageMedian) * e_sexf_cl^SEXF
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)
    cl_ohd <- exp(lcl_ohd + etalcl_ohd) *
      (1 + e_age_cl_ohd * (AGE - ageMedian) / ageMedian)
    vc_ohd <- exp(lvc_ohd + etalvc_ohd)
    q_ohd <- exp(lq_ohd)
    vp_ohd <- exp(lvp_ohd)

    # ---- Micro-constants ----
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    kel_ohd <- cl_ohd / vc_ohd
    k12_ohd <- q_ohd / vc_ohd
    k21_ohd <- q_ohd / vp_ohd

    # ---- ODE system (Figure 1) ----
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(central_ohd) <- kel * central - kel_ohd * central_ohd -
      k12_ohd * central_ohd + k21_ohd * peripheral1_ohd
    d/dt(peripheral1_ohd) <- k12_ohd * central_ohd - k21_ohd * peripheral1_ohd

    lag(depot) <- tlag

    # ---- Observations ----
    Cc <- central / vc * mgPerLToNgPerMl
    Cc_ohd <- central_ohd / vc_ohd * (mwOhd / mwDab) * mgPerLToNgPerMl

    Cc ~ prop(propSd)
    Cc_ohd ~ prop(propSd_ohd)
  })
}
