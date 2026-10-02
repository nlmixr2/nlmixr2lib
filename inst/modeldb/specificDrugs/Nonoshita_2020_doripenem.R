Nonoshita_2020_doripenem <- function() {
  description <- "Two-compartment IV population PK model for doripenem in 21 Japanese adult intensive care unit patients, 9 of them on continuous renal replacement therapy (continuous hemodiafiltration) (Nonoshita 2020). Total clearance is the sum of a body clearance and, when CRRT is running, a CRRT clearance fixed to the filtrate (effluent) flow rate times the 0.919 unbound fraction of doripenem. Body clearance is estimated separately for the CRRT and non-CRRT strata, each with its own power effect of Cockcroft-Gault creatinine clearance centred on its own stratum median (52.75 and 62.25 mL/min). Body weight and serum albumin were screened but not retained."
  reference <- "Nonoshita K, Suzuki Y, Tanaka R, Kaneko T, Ohchi Y, Sato Y, Yasuda N, Goto K, Kitano T, Itoh H. Population pharmacokinetic analysis of doripenem for Japanese patients in intensive care unit. Sci Rep. 2020;10:22148. doi:10.1038/s41598-020-79076-6"
  vignette <- "Nonoshita_2020_doripenem"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = TRUE: checked against Nonoshita 2020 Methods
  # (total plasma doripenem by HPLC-UV; doses in mg) and Figure 1.
  compartmentData <- list(
    central = list(analyte = "doripenem", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "doripenem", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance estimated with the Cockcroft-Gault equation; raw mL/min, NOT BSA-normalized",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Nonoshita 2020 Methods: 'Ccr obtained by the Cockcroft-Gault equation'. Table 1:",
        "68.0 +/- 33.4 mL/min overall, 76.0 +/- 35.8 in the non-CRRT stratum and 57.3 +/- 26.5",
        "in the CRRT stratum; Supplementary Table S1 per-patient range 20.7-155.6 mL/min.",
        "Enters body clearance as a power term centred on the MEDIAN of the subject's own",
        "stratum: (CRCL / 62.25)^0.64 without CRRT and (CRCL / 52.75)^0.42 with CRRT. The",
        "Discussion explains the split: with CRRT running, serum creatinine is cleared by both",
        "the kidney and the circuit, so the Cockcroft-Gault value indexes combined kidney-plus-",
        "circuit function rather than kidney function alone. The 62.25 mL/min non-CRRT median",
        "is the mean of the 62.2 and 62.3 mL/min values in Table S1, consistent with the sixth",
        "and seventh of 12 sorted non-CRRT values."
      ),
      source_name = "Ccr"
    ),
    RRT_CRRT_STATUS = list(
      description = "Continuous renal replacement therapy (continuous hemodiafiltration) status indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no CRRT)",
      notes = paste(
        "Nonoshita 2020 Methods: 9 of 21 patients were receiving CRRT, conducted as continuous",
        "hemodiafiltration through a cellulose triacetate hemofilter. Subject-level status over",
        "the single 8 h sampling interval after the first dose. Switches (a) which stratum-",
        "specific body clearance and creatinine-clearance exponent apply and (b) the CRRT",
        "clearance arm on or off."
      ),
      source_name = "CRRT"
    ),
    RRT_CRRT_EFFLUENT_FLOW = list(
      description = "Filtrate (effluent) flow rate of the continuous hemodiafiltration circuit",
      units = "mL/h",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Nonoshita 2020 Methods: filtrate flow rate QE = 0.6-1.8 L/h, the sum of the dialysate",
        "flow QD (0.3-0.9 L/h) and replacement-fluid flow QS (0.3-0.9 L/h); blood flow",
        "80-100 mL/min. CRRT clearance is CL_CRRT (L/h) = QE (L/h) x 0.919, where 0.919 is the",
        "unbound fraction of doripenem (Hori 2006, the paper's reference 28), so the sieving",
        "coefficient is FIXED to the unbound fraction rather than absorbed into an estimated",
        "parameter. The source records QE in L/h; supply this column in mL/h (QE x 1000). Only",
        "read when RRT_CRRT_STATUS = 1; set to 0 for non-CRRT subjects. Per-patient values are",
        "not published; the Discussion's CRRT-stratum mean CL_CRRT of 1.31 L/h corresponds to a",
        "mean QE of about 1.43 L/h."
      ),
      source_name = "QE"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Nonoshita 2020 Table 2: screened on CLbody (model 5, dOFV 7.06), V1 (model 6, dOFV 3.01) and V2 (model 7). BW on V1 entered the full model (model 8) but was removed in backward elimination (dOFV 0.91, p = 0.34). Table 1: 61.5 +/- 13.6 kg."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Nonoshita 2020 Table 2: screened on V1, V2 and CL_CRRT. Alb on CL_CRRT entered the full model (model 8) but was removed in backward elimination (dOFV 0.92, p = 0.34). Cohort values and units are not reported; the g/L units above are the register's canonical units, not a source-paper statement."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 21L,
    n_studies = 1L,
    n_concentrations = 97L,
    age_range = "15-86 years (Supplementary Table S1)",
    age_mean = "61.8 +/- 18.9 years (Table 1)",
    weight_range = "28.8-92.4 kg (Supplementary Table S1)",
    weight_mean = "61.5 +/- 13.6 kg (Table 1)",
    sex_female_pct = 14.3,
    race_ethnicity = "Japanese (single-centre ICU, Oita University Hospital)",
    disease_state = "Adult intensive care unit inpatients treated with doripenem for severe infection: bacteremia, sepsis / septic shock, infective endocarditis, pneumonia, intraperitoneal infection, urinary-tract infection and postoperative infection. APACHE II 17.6 +/- 6.7 and SOFA 7.5 +/- 2.6 (Table 1). Patients given another carbapenem beforehand were excluded.",
    renal_function = "Cockcroft-Gault creatinine clearance 68.0 +/- 33.4 mL/min (range 20.7-155.6). 9 of 21 patients on continuous hemodiafiltration (cellulose triacetate membrane; QB 80-100 mL/min, QD 0.3-0.9 L/h, QS 0.3-0.9 L/h, QE 0.6-1.8 L/h), 8 of them for a renal indication; Ccr 57.3 +/- 26.5 mL/min in the CRRT stratum and 76.0 +/- 35.8 mL/min in the non-CRRT stratum.",
    dose_range = "Doripenem 500 mg (20 patients) or 250 mg (1 patient) as a 1 h intravenous infusion; PK sampled after the first dose only.",
    regions = "Japan (Oita University Hospital, Oita)",
    notes = "Plasma sampled before and 1, 2, 4, 6 and 8 h after the start of the first infusion; HPLC-UV assay, LOQ 0.5 ug/mL, one sample below LOQ. NONMEM 7.3.0, FOCE-I, ADVAN6. Evaluated by goodness-of-fit plots, 1000-replicate VPC and 1000-sample non-parametric bootstrap (Table 3)."
  )

  ini({
    # Structural parameters: Nonoshita 2020 Table 3 'Population mean' column,
    # also stated in the Abstract and the Results 'final model' paragraph:
    #   CL_total (L/h) = CL_body(non-CRRT) = 3.65 x (Ccr/62.25)^0.64     (no CRRT)
    #                  = CL_body(CRRT) + CL_CRRT
    #                  = 2.49 x (Ccr/52.75)^0.42 + CL_CRRT               (CRRT)
    #   CL_CRRT = QE x 0.919
    #   V1 = 10.04 L; V2 = 8.13 L; Q = 3.53 L/h
    # Body clearance is estimated separately in each CRRT stratum within one
    # joint fit, so it carries symmetric stratum suffixes.
    lcl_offcrrt <- log(3.65); label("Body clearance without CRRT at CRCL = 62.25 mL/min (L/h)") # Table 3: CLbody(non-CRRT) = 3.65 L/h (bootstrap median 3.59, 95% CI 2.63-4.24)
    lcl_oncrrt <- log(2.49); label("Body clearance with CRRT at CRCL = 52.75 mL/min, excluding CRRT clearance (L/h)") # Table 3: CLbody(CRRT) = 2.49 L/h (bootstrap median 2.49, 95% CI 1.58-4.11)
    lvc <- log(10.04); label("Central volume of distribution V1 (L)") # Table 3: V1 = 10.04 L (bootstrap median 7.47, 95% CI 3.18-10.06)
    lvp <- log(8.13); label("Peripheral volume of distribution V2 (L)") # Table 3: V2 = 8.13 L (bootstrap median 8.91, 95% CI 6.36-11.38)
    lq <- log(3.53); label("Intercompartmental clearance Q (L/h)") # Table 3: Q = 3.53 L/h (bootstrap median 4.59, 95% CI 3.60-5.52)

    # Creatinine-clearance power exponents, one per stratum (final-model
    # equation, Abstract and Results). No uncertainty is printed for them.
    e_crcl_cl_offcrrt <- 0.64; label("Power exponent on (CRCL/62.25) for body clearance without CRRT (unitless)") # Results final-model equation: (Ccr/62.25)^0.64
    e_crcl_cl_oncrrt <- 0.42; label("Power exponent on (CRCL/52.75) for body clearance with CRRT (unitless)") # Results final-model equation: (Ccr/52.75)^0.42

    # Unbound fraction of doripenem, used as the fixed sieving coefficient of
    # the CRRT clearance arm (Results: '0.919 represents the non-protein
    # binding rate of DRPM', citing Hori 2006). Taken from the literature, not
    # estimated.
    fu <- fixed(0.919); label("Unbound fraction of doripenem, used as the sieving coefficient for CRRT clearance (fraction)") # Results final-model equation: CL_CRRT = QE x 0.919

    # Inter-individual variability: Table 3 'Population mean / Inter-individual
    # variability (%)' column, exponential error model (Methods). Converted as
    # omega^2 = log(CV^2 + 1).
    # The Results sentence lists the CL values in the opposite order
    # ('CLbody(CRRT), CLbody(non-CRRT) ... 7.3%, 22.2%'); Table 3 and the
    # Discussion ('non-CRRT group: 17.0-7.3%, CRRT group: 29.4-22.2%') agree
    # with each other and are used here.
    etalcl_offcrrt ~ 0.0053149 # Table 3: IIV CLbody(non-CRRT) = 7.3%; log(0.073^2 + 1)
    etalcl_oncrrt ~ 0.048108 # Table 3: IIV CLbody(CRRT) = 22.2%; log(0.222^2 + 1)
    etalvc ~ 0.017274 # Table 3: IIV V1 = 13.2%; log(0.132^2 + 1)
    etalvp ~ 0.056913 # Table 3: IIV V2 = 24.2%; log(0.242^2 + 1)
    etalq ~ 0.016000 # Table 3: IIV Q = 12.7%; log(0.127^2 + 1)

    # Residual error. Methods: 'residual variability was evaluated using the
    # additive error model'. Table 3 prints the final residual row as an
    # estimate of 0.70 (row unit ug/mL) alongside a percentage of 36.5%, which
    # the Results text calls the residual variability. The two cannot be one
    # additive SD, so both are encoded as a combined additive + proportional
    # model; see the vignette's Assumptions and deviations section.
    addSd <- 0.70; label("Additive residual error (mg/L)") # Table 3: final residual variability estimate 0.70 ug/mL
    propSd <- 0.365; label("Proportional residual error (fraction)") # Table 3 and Results: residual variability 36.5%
  })

  model({
    # Stratum-specific body clearance (Results final-model equation), each
    # centred on its own stratum's median Cockcroft-Gault creatinine clearance.
    cl_body_offcrrt <- exp(lcl_offcrrt + etalcl_offcrrt) * (CRCL / 62.25)^e_crcl_cl_offcrrt
    cl_body_oncrrt <- exp(lcl_oncrrt + etalcl_oncrrt) * (CRCL / 52.75)^e_crcl_cl_oncrrt
    cl_body <- (1 - RRT_CRRT_STATUS) * cl_body_offcrrt + RRT_CRRT_STATUS * cl_body_oncrrt

    # CRRT clearance: filtrate flow (mL/h -> L/h) times the unbound fraction.
    cl_crrt <- RRT_CRRT_STATUS * fu * RRT_CRRT_EFFLUENT_FLOW / 1000

    cl <- cl_body + cl_crrt
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    q <- exp(lq + etalq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Figure 1: both CL_body and CL_CRRT leave the central compartment.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Total plasma doripenem: mg / L = mg/L (= ug/mL).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
