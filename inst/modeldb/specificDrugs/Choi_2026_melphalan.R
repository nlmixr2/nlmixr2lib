Choi_2026_melphalan <- function() {
  description <- "Two-compartment IV-infusion population PK model for melphalan in pediatric patients undergoing autologous hematopoietic stem cell transplantation with busulfan-, thiotepa-, etoposide/carboplatin-, BCNU- or fludarabine-containing conditioning regimens (Choi 2026); allometric body-weight scaling on CL, V1 and V2, a serum-creatinine power effect on CL, and a concomitant-busulfan multiplier on CL."
  reference <- "Choi JY, Kim B, Park HJ, Kim BK, Hong KT, Lee S, Lee S, Kang HJ. Population Pharmacokinetics of Melphalan in Pediatric Patients Undergoing Autologous Hematopoietic Stem Cell Transplantation with Various Conditioning Regimens. Eur J Drug Metab Pharmacokinet. 2026. doi:10.1007/s13318-026-01000-6"
  vignette <- "Choi_2026_melphalan"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Choi 2026 Methods 2.3 (melphalan
  # quantified in plasma by LC-MS/MS) and Methods 2.4 (zero-order IV infusion
  # into the central compartment).
  compartmentData <- list(
    central     = list(analyte = "melphalan", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "melphalan", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric power effect on CL, V1 and V2, each normalised to 28 kg (Choi 2026 Section 3.3 final-model equations). The Table 2 footnote states that 28 kg is the rounded-up median weight of the study population (observed median 27.6 kg, range 5.7-84.8 kg). The exponents were ESTIMATED (0.771 on CL, 0.872 on V1, 0.553 on V2), not fixed to the theoretical allometric 0.75 / 1.0 values, and weight was entered on all three parameters before any other covariate was screened (Table S1 footnote a). Q carries no weight effect: the Table 2 row is labelled 'Q (L/h)' with no '/28 kg' qualifier and the Table 2 footnote restricts the allometric exponents to CL, V1 and V2.",
      source_name        = "WT"
    ),
    CREAT = list(
      description        = "Serum creatinine on the day of melphalan administration",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL normalised to 0.51 mg/dL, the cohort median (Choi 2026 Table 1; final-model equation in Section 3.3). Exponent -0.686, i.e. clearance falls as creatinine rises. Observed range 0.34-0.67 mg/dL; all patients had normal renal function, so the covariate is only supported over that narrow interval. Creatinine-derived eGFR (Schwartz, CKiD U25) was also screened but was not retained after backward elimination (Section 4).",
      source_name        = "creatinine"
    ),
    CONMED_BUSULFAN = list(
      description        = "Concomitant busulfan in the conditioning regimen (1 = busulfan-containing regimen, 0 = no busulfan)",
      units              = "(binary)",
      type               = "categorical",
      reference_category = "0 (no concomitant busulfan)",
      notes              = "Multiplicative effect on CL: clearance is multiplied by 0.846 when the regimen contains busulfan (Choi 2026 Section 3.3, equation for CLB; Table 2 row 'CL~RegimenB'). Busulfan was given on days -9 to -6 (BuMel) or -9 to -7 (BuMelThio), i.e. finishing before the melphalan dose, so the flag marks a carried-over regimen effect rather than simultaneous exposure. In the source cohort 11 of 20 patients received a busulfan-containing regimen (7 BuMel, 4 BuMelThio). The paper's own MCMP power analysis (Table S6) supports detection of the effect at this sample size, but the authors describe the model as exploratory and the effect size as imprecise.",
      source_name        = "RegimenB"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "categorical",
      notes       = "Screened on CL in the stepwise analysis (Table S1) but not retained in the final model."
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate (Schwartz and CKiD U25 creatinine-based equations)",
      units       = "mL/min/1.73 m^2",
      type        = "continuous",
      notes       = "Screened on CL (Table S1) but not retained after backward elimination; the authors kept raw serum creatinine instead (Section 4). Cystatin-C-based eGFR was considered but excluded because pre-dose cystatin C was incomplete (Section 2.5)."
    ),
    HCT = list(
      description = "Hematocrit",
      units       = "%",
      type        = "continuous",
      notes       = "Screened on CL (Table S1) but not retained in the final model. Observed median 32.4 % (range 26.7-40.4 %)."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/dL",
      type        = "continuous",
      notes       = "Screened on V1 and V2 (Table S1) but not retained in the final model. Observed median 3.8 g/dL (range 3.3-4.3 g/dL)."
    ),
    PRIOR_RADIATION = list(
      description = "Prior radiation therapy indicator",
      units       = "(binary)",
      type        = "categorical",
      notes       = "Screened on CL (Table S1) but not retained in the final model. 6 of 20 patients (30 %) had prior radiation therapy."
    ),
    AGE = list(
      description = "Age at melphalan infusion",
      units       = "years",
      type        = "continuous",
      notes       = "Not part of the final model. Tested on CL only as a sensitivity analysis for maturational confounding (Section 3.3; Table S3), as a power model (AGE/9.3)^theta with theta between -0.097 and -0.343 depending on the reference model. Adding age to the final model did not remove the busulfan effect (multiplier moved 0.846 -> 0.871), so the authors retained the final model without age."
    ),
    CONMED_ETOPOSIDE = list(
      description = "Concomitant etoposide in the conditioning regimen",
      units       = "(binary)",
      type        = "categorical",
      notes       = "Screened on CL as one of the tested concomitant drugs (Table S1 footnote b) but not retained in the final model."
    ),
    CONMED_CARBOPLATIN = list(
      description = "Concomitant carboplatin in the conditioning regimen",
      units       = "(binary)",
      type        = "categorical",
      notes       = "Screened on CL as one of the tested concomitant drugs (Table S1 footnote b) but not retained in the final model."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 20L,
    n_studies      = 1L,
    n_observations = 140L,
    age_range      = "0.6-17.4 years",
    age_median     = "9.3 years",
    weight_range   = "5.7-84.8 kg",
    weight_median  = "27.6 kg",
    sex_female_pct = 50,
    race_ethnicity = "Not reported; single-centre Korean cohort",
    disease_state  = "Pediatric autologous haematopoietic stem cell transplantation for solid tumours and haematologic disease (osteosarcoma 6, brain tumour 4, Ewing sarcoma 4, non-Hodgkin lymphoma 3, neuroblastoma 1, benign haematologic disease 2)",
    dose_range     = "140 mg/m^2 as a single dose (BuMel, FluMel, BEAM); 50 mg/m^2 once daily on two consecutive days (BuMelThio); 140 mg/m^2 then 70 mg/m^2 on consecutive days (MEC). Every dose given as a 30-minute IV infusion.",
    regions        = "South Korea (Seoul National University Children's Hospital; NCT04937634)",
    renal_function = "All patients had normal renal and hepatic function on the day of conditioning. Serum creatinine median 0.51 mg/dL (range 0.34-0.67); eGFR (Schwartz Cr) median 102.05 mL/min/1.73 m^2 (range 76.10-149.35).",
    co_medication  = "11 of 20 patients received a busulfan-containing regimen (BuMel 7, BuMelThio 4); the remaining 9 received MEC (4), BEAM (3) or FluMel (2).",
    notes          = "Single-centre prospective study, September 2020 to November 2021. Sampling was five points per dosing occasion: pre-dose and 5, 40, 70 and 170 minutes after the end of infusion; 28 occasions x 5 samples = 140 plasma samples. Assay LLOQ 5 ng/mL (0.005 mg/L); all post-dose samples were above LLOQ. The authors describe the model as exploratory owing to the small sample size."
  )

  ini({
    # Structural parameters: Choi 2026 Table 2 ('Final model' estimate column).
    # Final-model covariate equations, Section 3.3:
    #   CL0 (L/h) = 24.9 * (WT/28)^0.771 * (creatinine/0.51)^-0.686
    #   CLB (L/h) = 0.846 * CL0                (concomitant busulfan)
    #   V1  (L)   = 13.4 * (WT/28)^0.872
    #   V2  (L)   = 7.37 * (WT/28)^0.553
    # Q carries no covariate (Table S1: no covariate was tested on Q).
    lcl <- log(24.9);  label("Typical clearance at 28 kg, serum creatinine 0.51 mg/dL and no concomitant busulfan (L/h)")  # Table 2: CL 24.9, RSE 7.7%
    lvc <- log(13.4);  label("Typical central volume of distribution at 28 kg (L)")                                        # Table 2: V1 13.4, RSE 8.7%
    lq  <- log(14.5);  label("Intercompartmental clearance (L/h)")                                                         # Table 2: Q 14.5, RSE 22.9%
    lvp <- log(7.37);  label("Typical peripheral volume of distribution at 28 kg (L)")                                     # Table 2: V2 7.37, RSE 12.2%

    # Covariate-effect parameters, all estimated (Table 2). The weight exponents
    # were estimated rather than fixed at the theoretical allometric values, so
    # they are NOT wrapped in fixed().
    e_wt_cl              <- 0.771;   label("Body-weight power exponent on clearance, reference 28 kg (unitless)")                    # Table 2: CL~WT 0.771, RSE 15.0%
    e_wt_vc              <- 0.872;   label("Body-weight power exponent on central volume, reference 28 kg (unitless)")               # Table 2: V1~WT 0.872, RSE 15.1%
    e_wt_vp              <- 0.553;   label("Body-weight power exponent on peripheral volume, reference 28 kg (unitless)")            # Table 2: V2~WT 0.553, RSE 14.2%
    e_creat_cl           <- -0.686;  label("Serum-creatinine power exponent on clearance, reference 0.51 mg/dL (unitless)")          # Table 2: CL~creatinine -0.686, RSE 14.7%
    e_conmed_busulfan_cl <- 0.846;   label("Multiplicative clearance factor for a busulfan-containing conditioning regimen (unitless)")  # Table 2: CL~RegimenB 0.846, RSE 5.8%

    # Inter-individual variability. Choi 2026 Methods 2.4 gives a log-normal IIV
    # model, theta_i = theta_TV * exp(eta_i). Table 2 reports the IIV rows under
    # a '(%)' qualifier as CL 27.9, V1 39.9, V2 7.0 and the covariance row,
    # which carries no such qualifier, as CL~V1 0.11.
    #
    # Scale: the '(%)' values are 100 * omega (the log-scale SD), NOT a
    # back-transformed CV. Three independent checks agree:
    #  1. Table 3 is a forward simulation from this model. Because AUC is
    #     inversely proportional to CL, the log-scale SD of the simulated
    #     AUC0-24 recovers omega_CL directly: across the six scenarios
    #     log(p95/p05)/(2*1.645) = 0.275, 0.282, 0.275, 0.287, 0.279, 0.285
    #     -- mean 0.280 against the printed 27.9.
    #  2. Reading the '(%)' values as a back-transformed CV
    #     (omega = sqrt(log(1+CV^2))) gives omega_CL 0.2737 and omega_V1 0.3843,
    #     for which the printed covariance 0.11 implies a correlation of 1.046
    #     -- not a valid covariance matrix.
    #  3. For V1 the two readings differ enough to separate: 100*omega = 39.9
    #     (printed exactly) versus CV% = 41.5.
    # So omega^2 = (value/100)^2: CL 0.279^2 = 0.077841, V1 0.399^2 = 0.159201,
    # V2 0.070^2 = 0.0049. The covariance row sits under an unqualified header
    # and is therefore raw NONMEM OMEGA(2,1) = 0.11, giving a CL-V1 correlation
    # of 0.11/(0.279*0.399) = 0.989. That near-unity correlation is corroborated
    # by the size of the improvement the covariance bought: Table S2 model 4
    # gives dOFV = -38.832 on 1 degree of freedom for adding it.
    etalcl + etalvc ~ c(0.077841,
                        0.110000, 0.159201)  # Table 2: 'CL (%)' 27.9, 'CL~V1' 0.11, 'V1 (%)' 39.9
    etalvp ~ 0.004900                        # Table 2: 'V2 (%)' 7.0

    # Residual error: proportional only (Methods 2.4; Table S2 model 3 shows a
    # combined error model gave no OFV improvement over the proportional model).
    # Table 2 prints 'Proportional 0.0829' with no scale qualifier. Read as the
    # raw NONMEM SIGMA variance, so the SD is sqrt(0.0829) = 0.2879. Evidence:
    #  - estimate +/- 1.96 * RSE * estimate = 0.0829 +/- 0.0167 = (0.0662,
    #    0.0996) reproduces the printed bootstrap CI (0.0616, 0.0946) in width,
    #    so the RSE and the estimate share a scale -- and a NONMEM RSE is on the
    #    variance scale.
    #  - Table position: the authors qualified exactly one block of Table 2 with
    #    '(%)' (the IIV rows). The residual row is unqualified and is therefore
    #    raw NONMEM output.
    # Competing reading: if 0.0829 were already an SD, the proportional error
    # would be 8.29%. To adopt that reading, change the value below to 0.0829.
    propSd <- 0.287924;  label("Proportional residual error (fraction)")  # Table 2: 'Proportional' 0.0829 (variance), RSE 10.3%
  })
  model({
    # Typical-value covariate model (Choi 2026 Section 3.3). The busulfan term is
    # written as a power of the binary flag so that CONMED_BUSULFAN = 0 leaves
    # clearance at CL0 and CONMED_BUSULFAN = 1 multiplies it by 0.846.
    cl <- exp(lcl + etalcl) * (WT / 28)^e_wt_cl * (CREAT / 0.51)^e_creat_cl *
      e_conmed_busulfan_cl^CONMED_BUSULFAN
    vc <- exp(lvc + etalvc) * (WT / 28)^e_wt_vc
    q  <- exp(lq)
    vp <- exp(lvp + etalvp) * (WT / 28)^e_wt_vp

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Linear two-compartment disposition with first-order elimination from the
    # central compartment. Melphalan is given as a zero-order IV infusion
    # directly into central (Methods 2.4: "Drug administration was modeled as a
    # zero-order infusion process, with infusion durations specified in the
    # dataset"), so there is no depot and no bioavailability term.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Plasma melphalan concentration: dose mg / volume L -> mg/L. Choi 2026
    # tabulates exposures in ug/L and ug*h/L (Table 3); multiply by 1000.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
