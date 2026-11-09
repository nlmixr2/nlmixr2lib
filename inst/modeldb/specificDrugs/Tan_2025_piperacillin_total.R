Tan_2025_piperacillin_total <- function() {
  description <- "Two-compartment IV population PK model for TOTAL plasma piperacillin (given as piperacillin-tazobactam) in critically ill children receiving perioperative prophylaxis (Tan 2025), with allometric body-weight scaling, a bedside-Schwartz creatinine-clearance power effect on CL, and inter-occasion variability on CL for the pre- and post-operative periods (the intraoperative period is the zero-IOV reference). The paper's unbound-piperacillin fit is Tan_2025_piperacillin_unbound."
  reference <- "Tan WR, Irie K, McIntire C, Torres JL, Jones R, Gibson A, Mizuno T, Tang Girdwood S. Model-informed dose optimization for prophylactic piperacillin-tazobactam in perioperative pediatric critically ill patients. Antimicrob Agents Chemother. 2025;69(3):e01227-24. doi:10.1128/aac.01227-24"
  vignette <- "Tan_2025_piperacillin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "piperacillin (total)", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "piperacillin (total)", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling to 70 kg with fixed exponents 0.75 on CL and Q and 1 on V1 and V2 (Methods, 'Population PK modeling'; Table S1 footnote 1). Cohort median 20.0 kg (IQR 12.65-55.33; Table 1).",
      source_name = "BW"
    ),
    CRCL = list(
      description = "Creatinine clearance estimated with the bedside Schwartz equation, BSA-normalized",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Methods: 'creatinine clearance (CrCL) based on bedside Schwartz equation'. Power effect on CL, (CRCL/120)^0.484, where 120 mL/min/1.73 m^2 is the value the Methods call the 'median value ... in the general population' -- it is NOT the cohort median, which is 143 mL/min/1.73 m^2 (IQR 102.81-181.53; Table 1). The paper does not state whether CrCL was handled as a baseline or a time-varying covariate; Table 1 and Figure 1 report baseline CrCL.",
      source_name = "CrCL"
    ),
    OCC = list(
      description = "Surgical-period occasion index for inter-occasion variability on CL",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "OCC = 1 preoperative, OCC = 2 intraoperative, OCC = 3 postoperative.",
        "Methods: 'The IOV is initially set to 0, which corresponds to the",
        "intraoperative period ... Then the model estimated separate IOV terms",
        "for the pre-surgery and post-surgery periods.' The intraoperative",
        "period therefore carries no IOV eta, and the pre- and post-operative",
        "periods carry their own, separately estimated, variances (Table 2).",
        "Decomposed inside model() into the indicators oc1 / oc3. Records with",
        "any other OCC value carry no IOV, which reproduces the intraoperative",
        "setting the paper's dosing simulations target."
      ),
      source_name = "occasion (preoperative / intraoperative / postoperative)"
    )
  )

  covariatesDataExcluded <- list(
    PMA = list(
      description = "Postmenstrual age",
      units = "weeks",
      type = "continuous",
      notes = "Tested on CL for maturation; dOFV 0.00 (Table S1), not retained. Only seven patients (21%) were younger than 2 years (Discussion)."
    ),
    SURG_LIVER_TRANSPLANT = list(
      description = "Liver-transplantation surgery indicator",
      units = "(binary)",
      type = "binary",
      notes = "Tested on CL, dOFV -4.301 in forward inclusion (Table S1); removed in backward elimination (footnote 2). No canonical register column exists; recorded here for provenance only."
    ),
    SURG_TPIAT = list(
      description = "Total pancreatectomy with islet auto-transplantation surgery indicator",
      units = "(binary)",
      type = "binary",
      notes = "Tested on CL, dOFV -1.42 (Table S1); not retained."
    ),
    SURG_ABDOMINAL = list(
      description = "Non-transplant abdominal surgery indicator",
      units = "(binary)",
      type = "binary",
      notes = "Tested on CL, dOFV -0.5 (Table S1); not retained."
    ),
    SURG_ENT = list(
      description = "Ear, nose and throat surgery indicator",
      units = "(binary)",
      type = "binary",
      notes = "Tested on CL, dOFV -1.264 (Table S1); not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 34L,
    n_studies = 1L,
    n_observations = 286L,
    age_median = "7.07 years (IQR 3.05-14.48)",
    weight_median = "20.0 kg (IQR 12.65-55.33)",
    sex_female_pct = 55.9,
    race_ethnicity = "Not reported",
    disease_state = "Critically ill children in the paediatric intensive care unit receiving prophylactic piperacillin-tazobactam around surgery: liver transplantation (15, 44.1%), total pancreatectomy with islet auto-transplantation (2, 5.9%), non-transplant abdominal surgery (9, 26.5%), ear-nose-throat surgery (8, 23.5%).",
    dose_range = "Piperacillin 80-90 mg/kg every 2-3 h for patients below 45 kg; 3000 or 4000 mg for patients above 45 kg; 30-min IV infusions.",
    renal_function = "Baseline bedside-Schwartz CrCL median 143 mL/min/1.73 m^2 (IQR 102.81-181.53)",
    surgery_duration = "median 5.32 h (IQR 2.93-8.89)",
    regions = "USA (Cincinnati Children's Hospital Medical Center, single centre)",
    notes = "Table 1. Opportunistic scavenged sampling of residual clinical blood, 1-19 samples per patient, collected pre-, intra- and post-operatively; samples drawn during the 30-min infusion were excluded. Total and free piperacillin were measured by HPLC in the same samples and fitted as two separate models (286 total, 285 free concentrations)."
  )

  ini({
    # Table 2 (final total-piperacillin model). Typical values refer to a
    # 70 kg patient with CrCL 120 mL/min/1.73 m^2.
    lcl <- log(7.1);  label("Clearance at 70 kg and CrCL 120 mL/min/1.73 m^2 (L/h)") # Table 2 'CL (L/h/70 kg)' = 7.1 (RSE 9%)
    lvc <- log(17.2); label("Central volume of distribution at 70 kg (L)")           # Table 2 'V1 (L/70 kg)' = 17.2 (RSE 16%)
    lq <- log(3.25);  label("Intercompartmental clearance at 70 kg (L/h)")           # Table 2 'Q (L/h/70 kg)' = 3.25 (RSE 74%)
    lvp <- log(4.79); label("Peripheral volume of distribution at 70 kg (L)")        # Table 2 'V2 (L/70 kg)' = 4.79 (RSE 43%)

    e_crcl_cl <- 0.484; label("Power exponent of CrCL/120 on CL (unitless)") # Table 2 'theta for CrCL on CL' = 0.484 (RSE 30%)

    # Allometric exponents fixed per Methods ('a power coefficient of 0.75 for
    # CL and Q, and 1 for V1 and V2').
    e_wt_cl <- fixed(0.75); label("Allometric exponent of WT/70 on CL and Q (unitless)")
    e_wt_vc <- fixed(1);    label("Allometric exponent of WT/70 on V1 and V2 (unitless)")

    # IIV on CL and V1 only; Q and V2 IIV set to 0 (Results). Table 2 reports
    # CV%; omega^2 = log(1 + CV^2).
    etalcl ~ 0.07915 # Table 2 'IIV for CL (CV%)' = 28.7 -> log(1 + 0.287^2)
    etalvc ~ 0.15538 # Table 2 'IIV for V1 (CV%)' = 41 -> log(1 + 0.41^2)

    # IOV on CL, separate variances for the pre- and post-operative periods;
    # the intraoperative period has none (Methods).
    etaiov_cl_1 ~ 0.13289 # Table 2 'IOV for preoperative (CV%)' = 37.7 -> log(1 + 0.377^2)
    etaiov_cl_3 ~ 0.10636 # Table 2 'IOV for postoperative (CV%)' = 33.5 -> log(1 + 0.335^2)

    propSd <- 0.371; label("Proportional residual error (fraction)") # Table 2 'Proportional error (CV%)' = 37.1
  })
  model({
    # Occasion indicators: OCC = 1 preoperative, 3 postoperative; the
    # intraoperative period (OCC = 2) is the zero-IOV reference.
    oc1 <- (OCC == 1)
    oc3 <- (OCC == 3)
    iov_cl <- oc1 * etaiov_cl_1 + oc3 * etaiov_cl_3

    # Final CL equation (Results): CL = 7.1 * (WT/70)^0.75 * (CrCL/120)^0.484
    cl <- exp(lcl + etalcl + iov_cl) * (WT / 70)^e_wt_cl * (CRCL / 120)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q <- exp(lq) * (WT / 70)^e_wt_cl
    vp <- exp(lvp) * (WT / 70)^e_wt_vc

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Total plasma piperacillin: mg / L = mg/L
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
