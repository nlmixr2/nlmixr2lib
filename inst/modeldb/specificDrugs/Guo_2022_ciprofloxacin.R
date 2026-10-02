Guo_2022_ciprofloxacin <- function() {
  description <- paste(
    "Two-compartment population PK model for intravenous ciprofloxacin in",
    "critically ill adult ICU patients, developed on individual data pooled",
    "from three Dutch ICU studies (Guo 2022). First-order elimination from the",
    "central compartment; total body weight enters all four structural",
    "parameters allometrically (exponent 0.75 on CL and Q, 1 on V1 and V2,",
    "reference 70 kg) and MDRD eGFR enters CL as a linear effect centered on",
    "58.64 mL/min/1.73 m^2. IIV on CL, V1 and V2, inter-occasion variability",
    "on CL with one occasion per 24 h of therapy, and a study-specific",
    "residual error (combined additive + proportional for study 1,",
    "proportional for studies 2 and 3)."
  )
  reference <- paste(
    "Guo T, Abdulla A, Koch BCP, van Hasselt JGC, Endeman H, Schouten JA,",
    "Elbers PWG, Bruggemann RJM, van Hest RM; Dutch Antibiotic PK/PD",
    "Collaborators. Pooled Population Pharmacokinetic Analysis for Exploring",
    "Ciprofloxacin Pharmacokinetic Variability in Intensive Care Patients.",
    "Clin Pharmacokinet. 2022;61:869-879. doi:10.1007/s40262-022-01114-5"
  )
  vignette <- "Guo_2022_ciprofloxacin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "ciprofloxacin (total)", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ciprofloxacin (total)", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Added a priori to every structural parameter (Methods 2.2.2):",
        "(WT / 70)^0.75 on CL and Q, (WT / 70)^1 on V1 and V2 (Equations 3-6).",
        "Cohort median 80 kg (IQR 69-93), Table 1."
      ),
      source_name = "WGT"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate by the MDRD equation, BSA-normalized",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Linear effect on CL, [1 + 0.008 x (eGFR - 58.64)] (Equation 3), with",
        "58.64 mL/min/1.73 m^2 the median eGFR of the analysis data set",
        "(Results 3.1). Time-varying in the source; the paper tested and",
        "rejected separating the baseline and time-varying effects. Table 1",
        "median 59 (IQR 37-96). The factor falls to zero at eGFR = -66.4, so",
        "it stays positive over the whole physiological range."
      ),
      source_name = "eGFR"
    ),
    OCC = list(
      description = "Integer-valued occasion index for inter-occasion variability on clearance",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Guo 2022 Methods 2.2.1 defines every 24 h of therapy as an occasion",
        "(OCC = 1 for 0-24 h after the first dose, OCC = 2 for 24-48 h, ...).",
        "The paper does not state how many occasions its data spanned, so",
        "this file encodes seven, covering the first week of therapy; the",
        "per-occasion CL etas share the single estimated IOV variance",
        "(occasions 2-7 fixed equal to occasion 1). Records with OCC outside",
        "1..7 carry no IOV. The paper's own exposure simulation covers only",
        "the first 24 h (OCC = 1)."
      ),
      source_name = "OCC"
    ),
    STUDY_EXPAT = list(
      description = "Study 1 indicator (EXPAT study, Erasmus University Medical Center and Maasstad Hospital, Rotterdam; Abdulla 2020)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (study 2, Radboud University Medical Center; Gieling 2020)",
      notes = paste(
        "Selects the study 1 residual error: combined additive (0.151 mg/L)",
        "+ proportional (17.5%) instead of the reference proportional 13.7%",
        "(Table 2, final model). Has no effect on any structural or IIV",
        "parameter; the authors deliberately did not test 'study' as a",
        "covariate (Discussion). Mutually exclusive with STUDY_RDRN."
      ),
      source_name = "Study 1"
    ),
    STUDY_RDRN = list(
      description = "Study 3 indicator (Right Dose Right Now study, Amsterdam UMC location VUmc)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (study 2, Radboud University Medical Center; Gieling 2020)",
      notes = paste(
        "Selects the study 3 proportional residual error of 24.5% instead of",
        "the reference 13.7% (Table 2, final model). Has no effect on any",
        "structural or IIV parameter. Mutually exclusive with STUDY_EXPAT."
      ),
      source_name = "Study 3"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened in the stepwise covariate search (Methods 2.2.2); not retained. Cohort median 67 years (IQR 59-74), Table 1."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "Screened; not retained. 34% female, Table 1."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Screened (time-varying; baseline and time-varying parts also tested separately); not retained. Cohort median 97 umol/L (IQR 69-156), Table 1."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened (time-varying); not retained. Cohort median 25 g/L (IQR 20-28), Table 1."
    ),
    SOFA = list(
      description = "Sequential Organ Failure Assessment score",
      units = "(score)",
      type = "continuous",
      notes = "Screened (time-varying); not retained. Cohort median 10 (IQR 8-13), Table 1.",
      source_name = "SOFA score"
    ),
    RRT_CRRT_STATUS = list(
      description = "Continuous veno-venous hemofiltration",
      units = "(binary)",
      type = "binary",
      notes = "Screened as CVVH support; not retained. 8% of patients on CVVH, Table 1.",
      source_name = "CVVH"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 140L,
    n_studies = 3L,
    n_observations = 1094L,
    age_median = "67 years (IQR 59-74)",
    weight_median = "80 kg (IQR 69-93)",
    sex_female_pct = 34,
    disease_state = paste(
      "Critically ill adult ICU patients receiving intravenous ciprofloxacin:",
      "study 1 mostly respiratory infection (47.6%) and sepsis (19%); study 2",
      "mostly pneumonia (72%); study 3 sepsis / septic shock (55% met the",
      "sepsis-3 criteria for septic shock). Median SOFA 10 (IQR 8-13); 8% on",
      "continuous veno-venous hemofiltration."
    ),
    dose_range = paste(
      "400 mg IV twice or three times daily, infused over 30-60 min (study 1);",
      "400 mg IV twice daily (study 2); 400 mg IV three times daily or",
      "model-based individualized dosing (study 3)."
    ),
    regions = "The Netherlands (Erasmus MC and Maasstad Hospital Rotterdam; Radboud UMC Nijmegen; Amsterdam UMC location VUmc).",
    renal_function = "MDRD eGFR median 59 mL/min/1.73 m^2 (IQR 37-96); serum creatinine median 97 umol/L (IQR 69-156).",
    notes = paste(
      "Pooled individual data from three prospective clinical studies: study 1",
      "(n = 42, 204 samples; Abdulla 2020), study 2 (n = 39, 531 samples;",
      "Gieling 2020) and study 3 (n = 59, 359 samples). Total plasma",
      "concentrations. NONMEM 7.5 FOCE-I. Demographics from Table 1."
    )
  )

  ini({
    # Structural parameters: Table 2 'Final model' column (RSE in brackets),
    # typical values for a 70-kg patient with eGFR = 58.64 mL/min/1.73 m^2
    # (Equations 3-6).
    lcl <- log(14.7); label("Clearance at 70 kg and eGFR 58.64 mL/min/1.73 m^2 (L/h)")    # Table 2, CL = 14.7 L/h (4.1%); Eq. 3
    lvc <- log(61.2); label("Central volume of distribution V1 at 70 kg (L)")               # Table 2, V1 = 61.2 L (7.1%); Eq. 4
    lq <- log(44.9); label("Intercompartmental clearance at 70 kg (L/h)")                   # Table 2, Q = 44.9 L/h (6.7%); Eq. 5
    lvp <- log(71.6); label("Peripheral volume of distribution V2 at 70 kg (L)")            # Table 2, V2 = 71.6 L (6%); Eq. 6

    # Covariate effects. Allometric exponents were fixed a priori
    # (Methods 2.2.2); the eGFR slope is fractional per mL/min/1.73 m^2
    # (Equation 3).
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL (unitless)")     # Methods 2.2.2, power 0.75 on CL; Eq. 3
    e_wt_q <- fixed(0.75); label("Allometric exponent of body weight on Q (unitless)")       # Methods 2.2.2, power 0.75 on intercompartmental CL; Eq. 5
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on V1 (unitless)")        # Methods 2.2.2, power 1 on volumes; Eq. 4
    e_wt_vp <- fixed(1); label("Allometric exponent of body weight on V2 (unitless)")        # Methods 2.2.2, power 1 on volumes; Eq. 6
    e_crcl_cl <- 0.008; label("Linear eGFR effect on CL (fractional change per mL/min/1.73 m^2)") # Table 2, eGFR on CL (linear) = 0.008 (12%); Eq. 3

    # IIV reported as CV%; exponential etas, omega^2 = log(1 + CV^2), the
    # inverse of the paper's own CV formula (Equation 2).
    etalcl ~ 0.201130 # Table 2, IIV CL = 47.2% -> log(1 + 0.472^2)
    etalvc ~ 0.317232 # Table 2, IIV V1 = 61.1% -> log(1 + 0.611^2)
    etalvp ~ 0.198051 # Table 2, IIV V2 = 46.8% -> log(1 + 0.468^2)

    # IOV on CL, one occasion per 24 h (Methods 2.2.1). A single variance was
    # estimated; seven occasions are encoded (see OCC in covariateData).
    etaiov_cl_1 ~ 0.018327 # Table 2, IOV CL = 13.6% -> log(1 + 0.136^2)
    etaiov_cl_2 ~ fixed(0.018327) # same variance as occasion 1
    etaiov_cl_3 ~ fixed(0.018327) # same variance as occasion 1
    etaiov_cl_4 ~ fixed(0.018327) # same variance as occasion 1
    etaiov_cl_5 ~ fixed(0.018327) # same variance as occasion 1
    etaiov_cl_6 ~ fixed(0.018327) # same variance as occasion 1
    etaiov_cl_7 ~ fixed(0.018327) # same variance as occasion 1

    # Study-specific residual error (Table 2, final model). Study 2 is the
    # reference stratum (STUDY_EXPAT = STUDY_RDRN = 0).
    propSd <- 0.137; label("Proportional residual error, study 2 (fraction)")                  # Table 2, Prop Study2 = 13.7% (3.6%)
    propSd_expat <- 0.175; label("Proportional residual error, study 1 (fraction)")            # Table 2, Prop Study1 = 17.5% (9.5%)
    addSd_expat <- 0.151; label("Additive residual error, study 1 (mg/L)")                     # Table 2, Add Study1 = 0.151 mg/L (27.1%)
    propSd_rdrn <- 0.245; label("Proportional residual error, study 3 (fraction)")             # Table 2, Prop Study3 = 24.5% (5.1%)
  })

  model({
    # Occasion indicators (24-h occasions) and the per-occasion CL eta.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    oc6 <- (OCC == 6)
    oc7 <- (OCC == 7)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 +
      oc4 * etaiov_cl_4 + oc5 * etaiov_cl_5 + oc6 * etaiov_cl_6 +
      oc7 * etaiov_cl_7

    # Individual parameters (Equations 3-6).
    cl <- exp(lcl + etalcl + iov_cl) * (WT / 70)^e_wt_cl * (1 + e_crcl_cl * (CRCL - 58.64))
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q <- exp(lq) * (WT / 70)^e_wt_q
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vp

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # IV infusion into central (rate or dur on the dose record).
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Total plasma ciprofloxacin: mg / L = mg/L.
    Cc <- central / vc

    # Study-specific residual error: the additive component exists for
    # study 1 only.
    propSd_study <- propSd * (1 - STUDY_EXPAT - STUDY_RDRN) +
      propSd_expat * STUDY_EXPAT + propSd_rdrn * STUDY_RDRN
    addSd_study <- addSd_expat * STUDY_EXPAT
    Cc ~ add(addSd_study) + prop(propSd_study)
  })
}
