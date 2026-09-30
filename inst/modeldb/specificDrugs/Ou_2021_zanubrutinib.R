Ou_2021_zanubrutinib <- function() {
  description <- "Two-compartment population PK model with sequential zero-order then first-order absorption for oral zanubrutinib in healthy volunteers and adults with B-cell malignancies"
  reference <- paste(
    "Ou YC, Liu L, Tariq B, Wang K, Jindal A, Tang Z, Gao Y, Sahasranaman S.",
    "Population pharmacokinetic analysis of the BTK inhibitor zanubrutinib in",
    "healthy volunteers and patients with B-cell malignancies.",
    "Clin Transl Sci. 2021;14(2):764-772. doi:10.1111/cts.12948.",
    "PMCID: PMC7993273.",
    sep = " "
  )
  vignette <- "Ou_2021_zanubrutinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    ALT = list(
      description = "Baseline serum alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Baseline (time-fixed). Enters CL/F as a power function normalised to",
        "18 U/L, the typical-patient value of Eq. 5 and the cohort median",
        "(Table S1a: median 18, range 4-197 IU/L overall; 20 IU/L in healthy",
        "volunteers and 18 IU/L in patients). The paper writes the effect as",
        "-0.189 x log(ALT/18) inside the exponent of Eq. 5, i.e.",
        "(ALT/18)^-0.189. Three subjects with missing ALT were imputed to the",
        "population median (Supplementary Material, Handling of Missing",
        "Covariates)."
      ),
      source_name = "ALT"
    ),
    DIS_HEALTHY = list(
      description = "Healthy-volunteer indicator (1 = healthy volunteer, 0 = patient with a B-cell malignancy)",
      units = "(binary)",
      type = "binary",
      reference_category = "patient with a B-cell malignancy (DIS_HEALTHY = 0)",
      notes = paste(
        "Eq. 5 carries two indicator terms, 5.13 x Patient + 4.77 x HV, so",
        "each health-status group has its own typical log CL/F (Table 2:",
        "exp(theta1) = 170 L/h for patients, exp(theta10) = 118 L/h for",
        "healthy volunteers). Encoded here as the patient typical value plus",
        "a log-ratio shift for DIS_HEALTHY = 1. The patient cohort pools",
        "CLL/SLL, MCL, WM and other B-cell malignancies (Table S1b); tumor",
        "type was screened but not retained."
      ),
      source_name = "HV / Patient"
    ),
    OCC = list(
      description = "Occasion index for the between-occasion random effects",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Two occasions (Results, Base model development): OCC = 1 for",
        "records < 7 days after the first dose (single-dose data) and",
        "OCC = 2 for records >= 7 days after the first dose (after repeated",
        "dosing). Selects the per-occasion IOV etas on CL/F, Vc/F and D1.",
        "A single-occasion simulation may use OCC = 1 throughout."
      ),
      source_name = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Baseline age",
      units = "years",
      type = "continuous",
      notes = "Screened (19-90 years) and not a statistically significant covariate on CL/F or Vc/F (Results, Base model development and covariate assessment; Figure S1C)."
    ),
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened (36-144 kg) and not statistically significant on CL/F or Vc/F (Figure S1D)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened and not statistically significant (Results)."
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units = "(binary)",
      type = "binary",
      notes = "Race (Asian, white, other) screened and not statistically significant (Figure S1A)."
    ),
    CRCL_BASE = list(
      description = "Baseline Cockcroft-Gault creatinine clearance",
      units = "mL/min",
      type = "continuous",
      notes = "Mild or moderate renal impairment (CrCL >= 30 mL/min) screened and not statistically significant (Figure S1E)."
    ),
    AST = list(
      description = "Baseline aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened and not statistically significant (Results)."
    ),
    TBILI = list(
      description = "Baseline total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Screened and not statistically significant (Results)."
    ),
    CONMED_PPI = list(
      description = "Concomitant proton-pump inhibitor use",
      units = "(binary)",
      type = "binary",
      notes = "PPI and H2RA use were pooled as one acid-reducing-agent category and were not statistically significant (P > 0.074, ANOVA; Figure S2)."
    ),
    CONMED_H2RA = list(
      description = "Concomitant H2-receptor-antagonist use",
      units = "(binary)",
      type = "binary",
      notes = "Pooled with PPI use as acid-reducing agents; not statistically significant (Figure S2)."
    ),
    TUMTP_MCL = list(
      description = "Mantle cell lymphoma tumor-type indicator",
      units = "(binary)",
      type = "binary",
      notes = "Tumor type (MCL, CLL/SLL, WM, other) passed the P < 0.01 screen but was not retained by the NONMEM stepwise search (Table S2)."
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "zanubrutinib",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "zanubrutinib",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "zanubrutinib",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 632,
    n_studies = 9,
    n_observations = 4925,
    age_median = "64 years (range 19-90); healthy volunteers 43, patients 66",
    weight_median = "75 kg (range 36-144)",
    sex_female_pct = 29.7,
    race_ethnicity = c(White = 67.7, Asian = 23.1, Black = 4.1, Other = 3.0, Missing = 2.1),
    disease_state = paste(
      "90 healthy volunteers (14.2%) and 542 patients with B-cell",
      "malignancies: CLL/SLL 135 (21.4%), Waldenstrom's macroglobulinemia",
      "196 (31.0%), mantle cell lymphoma 70 (11.1%), other 141 (22.3%)"
    ),
    dose_range = "20-320 mg orally once or twice daily (most subjects 160 mg twice daily); the 480 mg single-dose data were excluded",
    hepatic_function = "baseline ALT median 18 IU/L (range 4-197)",
    renal_function = "baseline Cockcroft-Gault CrCL median 86 mL/min (range 13.5-240)",
    co_medication = "proton-pump inhibitors 21.5%, H2-receptor antagonists 9.2%",
    regions = "global phase I-III program (studies AU-003, 1002, 205, 206, 103, 104, 105, 106, 302)",
    notes = paste(
      "Demographics from Tables S1a/S1b; study list from Table 1 of Ou 2021.",
      "LLOQ 1 ng/mL; 6.05% of samples were below it and were omitted.",
      "Nine healthy volunteers dosed at 480 mg (study BGB-3111-106) were",
      "excluded because of less-than-dose-proportional exposure, and two",
      "subjects with extreme PK parameters were excluded."
    )
  )

  ini({
    # Structural parameters (Table 2, final model estimates). The typical
    # values are for a patient with a B-cell malignancy and ALT = 18 U/L
    # (Results, Final population PK model).
    lcl <- log(170); label("Apparent clearance CL/F, patient (L/h)")                          # Table 2, 'exp(theta1) CL/F (L/hour, patient)' = 170; Eq. 5 5.13 (exp = 169)
    lvc <- log(112); label("Apparent central volume Vc/F (L)")                                # Table 2, 'exp(theta2) Vc/F (L)' = 112; Eq. 6 4.72 (exp = 112)
    lq <- log(26.5); label("Apparent intercompartmental clearance Q/F (L/h)")                 # Table 2, 'exp(theta3) Q/F (L/hour)' = 26.5
    lvp <- log(345); label("Apparent peripheral volume Vp/F (L)")                             # Table 2, 'exp(theta4) Vp/F (L)' = 345
    lka <- log(0.526); label("First-order absorption rate constant ka (1/h)")                 # Table 2, 'exp(theta5) Ka (1/hour)' = 0.526
    ld1 <- log(1.128); label("Duration of zero-order input into the depot D1 (h)")            # Table 2, 'exp(theta6) D1 (hour)' = 1.128

    # Covariate effects on CL/F (Eq. 5).
    e_dis_healthy_cl <- log(118 / 170); label("Log-ratio of CL/F in healthy volunteers vs patients (unitless)") # Table 2, 'exp(theta10) CL/F (L/hour, HV)' = 118 vs 170; Eq. 5 4.77 vs 5.13
    e_alt_cl <- -0.189; label("Power exponent on (ALT/18) for CL/F (unitless)")                # Table 2, 'theta11 Influence of ALT on CL/F' = -0.189

    # IIV. Table 2 reports IIV and IOV as percent CV with the variance
    # recovered as (CV/100)^2. The (CL/F, Vc/F) covariance 0.129 settles
    # that convention: with (CV/100)^2 the implied correlation is
    # 0.129 / (0.367 * 0.371) = 0.947, whereas log(1 + CV^2) variances
    # (0.1261, 0.1289) would imply a correlation of 1.01, which is not a
    # valid covariance matrix.
    etalcl + etalvc ~ c(0.1347, 0.129, 0.1376)                                                # Table 2, IIV CL/F 36.7% (0.367^2), IIV Vc/F 37.1% (0.371^2), 'Covariance (CL/F, Vc/F)' = 0.129
    etalq ~ 1.0404                                                                            # Table 2, IIV Q/F 102% (1.02^2)
    etalvp ~ 0.7465                                                                           # Table 2, IIV Vp/F 86.4% (0.864^2)
    etald1 ~ 0.3881                                                                           # Table 2, IIV D1 62.3% (0.623^2)

    # IOV on CL/F, Vc/F and D1 over two occasions (< 7 days and >= 7 days
    # after the first dose). One variance per parameter, shared across
    # occasions (NONMEM $OMEGA BLOCK(1) SAME convention).
    etaiov_cl_1 ~ 0.0818                                                                      # Table 2, IOV CL/F 28.6% (0.286^2)
    etaiov_cl_2 ~ fix(0.0818)                                                                 # same variance as occasion 1 (SAME)
    etaiov_vc_1 ~ 0.4556                                                                      # Table 2, IOV Vc/F 67.5% (0.675^2)
    etaiov_vc_2 ~ fix(0.4556)                                                                 # same variance as occasion 1 (SAME)
    etaiov_d1_1 ~ 0.3881                                                                      # Table 2, IOV D1 62.3% (0.623^2); identical to the IIV D1 row as printed
    etaiov_d1_2 ~ fix(0.3881)                                                                 # same variance as occasion 1 (SAME)

    # Residual error, switched on time after the previous dose (TFDS) at 5 h:
    # additive only for TFDS < 5 h, combined additive + proportional for
    # TFDS >= 5 h (Results, Base model development).
    addSd_early <- 8.69; label("Additive residual error SD, TFDS < 5 h (ng/mL)")              # Table 2, 'theta9 Additive residual error (ng/mL, TFDS < 5 hour)' = 8.69
    addSd_late <- 0.633; label("Additive residual error SD, TFDS >= 5 h (ng/mL)")             # Table 2, 'theta7 Additive residual error (ng/mL, TFDS >= 5 hour)' = 0.633
    propSd_late <- 0.449; label("Proportional residual error SD, TFDS >= 5 h (fraction)")     # Table 2, 'theta8 Proportional residual error (%)' = 44.9
  })

  model({
    # 1. Occasion indicators (OCC = 1: < 7 days, OCC = 2: >= 7 days).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2
    iov_vc <- oc1 * etaiov_vc_1 + oc2 * etaiov_vc_2
    iov_d1 <- oc1 * etaiov_d1_1 + oc2 * etaiov_d1_2

    # 2. Individual parameters. Eq. 5:
    #   CL/F = exp(5.13 x Patient + 4.77 x HV - 0.189 x log(ALT/18) + eta)
    # Checked against the Figure 4 sensitivity analysis: (10/18)^-0.189 =
    # 1.118 (AUCss -10.5%), (37/18)^-0.189 = 0.873 (AUCss +14.6%), and
    # 170/118 = 1.44 (AUCss +43.7% in healthy volunteers).
    cl <- exp(lcl + e_dis_healthy_cl * DIS_HEALTHY + e_alt_cl * log(ALT / 18) + etalcl + iov_cl)
    vc <- exp(lvc + etalvc + iov_vc)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)
    ka <- exp(lka)
    d1 <- exp(ld1 + etald1 + iov_d1)

    # 3. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODEs (Figure 1). The dose enters the depot as a zero-order input of
    #    duration D1 (R1 = Dose/D1), then moves to the central compartment
    #    by first-order absorption. Dose records need rate = -2 so that
    #    dur(depot) is honoured.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    dur(depot) <- d1

    # 5. Observation and error. Dose in mg and volume in L give mg/L; x 1000
    #    converts to ng/mL (assay units).
    Cc <- 1000 * central / vc

    # Time-dependent residual error: additive only before 5 h after the
    # previous dose, combined additive + proportional from 5 h onward.
    # tad() is evaluated once on its own line.
    tfds <- tad()
    late <- tfds >= 5
    addSdTad <- addSd_early * (1 - late) + addSd_late * late
    propSdTad <- propSd_late * late
    Cc ~ add(addSdTad) + prop(propSdTad)
  })
}
