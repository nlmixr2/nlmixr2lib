Neroutsos_2022_busulfan <- function() {
  description <- "Two-compartment population PK model for intravenous busulfan in children undergoing haematopoietic stem cell transplantation (Neroutsos 2022), with a per-patient syringe-pump infusion lag time supplied as data, allometric body-weight scaling (fixed exponents 0.75 on CL and Q, 1 on V1 and V2) relative to a 70 kg reference, correlated inter-individual variability on CL and V1, three-occasion inter-occasion variability on CL, and a proportional residual error."
  reference <- paste(
    "Neroutsos E, Nalda-Molina R, Paisiou A, Zisaki K, Goussetis E,",
    "Spyridonidis A, Kitra V, Grafakos S, Valsami G, Dokoumetzidis A. (2022).",
    "Development of a Population Pharmacokinetic Model of Busulfan in Children",
    "and Evaluation of Different Sampling Schedules for Precision Dosing.",
    "Pharmaceutics 14(3):647. doi:10.3390/pharmaceutics14030647.",
    sep = " "
  )
  vignette <- "Neroutsos_2022_busulfan"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "busulfan", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "busulfan", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling on a 70 kg reference with fixed exponents: 0.75 on CL and Q, 1 on V1 and V2 (Neroutsos 2022 Section 3.2 final-model equations; supplement NONMEM control script: A=(BW/70)**THETA(5), TV1=THETA(2)*(BW/70)**THETA(6), V2=THETA(3)*(BW/70), Q=THETA(4)*(BW/70)**THETA(5), with THETA(5) = 0.75 FIX and THETA(6) = 1 FIX). Cohort mean 30.6 kg (SD 21.6), range 7.38-104 kg (Table 1).",
      source_name = "BW"
    ),
    T_INFUSION_LAG = list(
      description = "Syringe-pump infusion lag time: delay between starting the pump and the drug entering the circulation (h).",
      units = "h",
      type = "continuous",
      reference_category = NULL,
      notes = "Not estimated: supplied per patient in the analysis dataset and applied as NONMEM ALAG1 = TLAG on the infusion into the central compartment (supplement control script; Section 2.2 'The lag-time for each patient came from Table 2 and was included in the dataset'). Table 2 gives the value per body-weight dosing band, in minutes: <9 kg 40 min; 9-16 kg 40 min; 16-23 kg 35 or 25 min; 23-34 kg 20 min; >34 kg 10 or 5 min. The paper does not state which patients of the 16-23 kg and >34 kg bands received the lower value. Supply the value in hours (40 min = 0.6667 h). The lag shifts the whole 2 h infusion; it does not shorten it.",
      source_name = "TLAG"
    ),
    OCC = list(
      description = "Integer-valued occasion index for inter-occasion variability on CL.",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = "Values 1, 2, 3 each select their own IOV eta on log-CL (supplement control script: OCC1 = OC1*ETA(3) + OC2*ETA(4) + OC3*ETA(5), with ETA(4) and ETA(5) declared $OMEGA BLOCK(1) SAME). The control script also defines an OC4 indicator but never uses it, so OCC = 4 (or any other value) carries no IOV. Blood was sampled after the first dose on day 1 and, for most patients, again on day 2 (Section 2.1); the paper does not define the occasions further. Decomposed inside model() into binary indicators oc1..oc3.",
      source_name = "OCC"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 76L,
    n_studies = 1L,
    age_range = "0.5-19 years",
    age_mean = "7.6 years (SD 5.1)",
    weight_range = "7.38-104 kg",
    weight_mean = "30.6 kg (SD 21.6)",
    sex_female_pct = 35.5,
    disease_state = "Paediatric patients undergoing haematopoietic stem cell (bone marrow) transplantation after busulfan-containing conditioning: adrenoleukodystrophy, acute lymphoblastic and myeloid leukaemia, Blackfan-Diamond anaemia, Ewing sarcoma, myelodysplastic syndrome, stage IV neuroblastoma, non-Hodgkin lymphoma, thalassaemia, Wiskott-Aldrich syndrome.",
    dose_range = "Intravenous busulfan (Busilvex) 0.8-1.2 mg/kg by body-weight band (Table 2: <9 kg 1 mg/kg; 9-16 kg 1.2; 16-23 kg 1.1; 23-34 kg 0.95; >34 kg 0.8), 2 h infusion every 6 h for 16 doses, delivered by syringe pump.",
    regions = "Greece ('Agia Sofia' Children's Hospital of Athens, July 2014 - January 2017).",
    renal_function = "CKD-EPI 197 mL/min/1.73 m^2 (SD 41), range 107-346; serum creatinine 0.36 mg/dL (SD 0.15).",
    notes = "596 plasma busulfan concentrations (HPLC-PDA) sampled at nominal 0, 2, 2.5, 4 and 6 h after the start of the first infusion on day 1 and, for most patients, day 2. 49 of 76 patients male. Demographics from Table 1; weight-band dosing and lag times from Table 2. NONMEM 7.3, FOCE with interaction."
  )

  ini({
    # Structural parameters for a 70 kg patient -- Table 3 'Parameter estimates
    # using the final covariate PopPK model', NONMEM Estimate column. Mapping to
    # the supplement control script: THETA(1) = CL, THETA(2) = V1, THETA(3) = V2,
    # THETA(4) = Q.
    lcl <- log(10.7); label("Clearance CL for a 70 kg patient (L/h)") # Table 3 row 'CL (L/h)' = 10.7 (RSE 4.05%)
    lvc <- log(39.5); label("Central volume of distribution V1 for a 70 kg patient (L)") # Table 3 row 'V1 (L)' = 39.5 (RSE 6.84%)
    lq <- log(4.68); label("Intercompartmental clearance Q for a 70 kg patient (L/h)") # Table 3 row 'Q (L/h)' = 4.68 (RSE 15.2%)
    lvp <- log(17.5); label("Peripheral volume of distribution V2 for a 70 kg patient (L)") # Table 3 row 'V2 (L)' = 17.5 (RSE 17.2%)

    # Allometric exponents, fixed (supplement control script $THETA (0 0.75) FIX
    # and (0 1) FIX; Section 3.2 'allometric model of fixed exponents 0.75 and 1').
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent of body weight on CL and Q (unitless)") # control script THETA(5) = 0.75 FIX
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on V1 (unitless)") # control script THETA(6) = 1 FIX

    # IIV. Table 3 reports the random effects as standard deviations on the
    # log scale (CL IIV 0.284 = the abstract's '28%'; the control script's
    # initial $OMEGA 0.0697 has sqrt 0.264, same scale), so omega^2 = SD^2 and
    # cov = r * SD_CL * SD_V1.
    # omega^2(CL) = 0.284^2 = 0.080656; omega^2(V1) = 0.409^2 = 0.167281;
    # cov = 0.679 * 0.284 * 0.409 = 0.078870.
    etalcl + etalvc ~ c(0.080656, 0.078870, 0.167281) # Table 3 rows 'CL IIV' = 0.284, 'V1 IIV' = 0.409, 'Cor. CL-V1' = 0.679; $OMEGA BLOCK(2)

    # IOV on log-CL across three occasions, one shared variance
    # ($OMEGA BLOCK(1) 0.0125 then two BLOCK(1) SAME). omega^2 = 0.105^2.
    etaiov_cl_1 ~ 0.011025 # Table 3 row 'CL IOV' = 0.105 (SD); 0.105^2 = 0.011025
    etaiov_cl_2 ~ fixed(0.011025) # same variance as occasion 1 per $OMEGA BLOCK(1) SAME
    etaiov_cl_3 ~ fixed(0.011025) # same variance as occasion 1 per $OMEGA BLOCK(1) SAME

    # Residual error: proportional only (Y = F + F*EPS(1)); the additive term of
    # the base model was dropped (Section 3.2). Table 3 reports the SD.
    propSd <- 0.126; label("Proportional residual error (fraction)") # Table 3 row 'Prop. RE' = 0.126 (RSE 1.65%)
  })

  model({
    # Occasion indicators for the three IOV etas on log-CL (control script
    # OCC1 = OC1*ETA(3) + OC2*ETA(4) + OC3*ETA(5)).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3

    # Individual parameters, allometric scaling on a 70 kg reference. V2 scales
    # with a hard-coded exponent of 1 in the control script (V2=THETA(3)*(BW/70)).
    cl <- exp(lcl + etalcl + iov_cl) * (WT / 70)^e_wt_cl_q
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q <- exp(lq) * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp) * (WT / 70)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Syringe-pump infusion lag supplied as data (control script ALAG1 = TLAG);
    # doses are 2 h infusions into the central compartment.
    alag(central) <- T_INFUSION_LAG

    # Dose mg, V1 L -> mg/L (control script S1 = V1).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
