Das_2019_ceftaroline <- function() {
  description <- "Joint two-compartment ceftaroline fosamil (prodrug) plus two-compartment ceftaroline (active moiety) population PK model in 951 adults and adolescents pooled from 21 clinical studies (healthy subjects, renal impairment and end-stage renal disease, and patients with complicated skin and soft tissue infection or community-acquired pneumonia), used to support ceftaroline fosamil 600 mg q8h by 2 h infusion and the EUCAST q8h Staphylococcus aureus breakpoints (Das 2019). Ceftaroline fosamil is converted 100% (molar) to ceftaroline, so the whole prodrug clearance enters the ceftaroline central compartment scaled by the molar-mass ratio. Body weight enters allometrically on all clearances (0.75) and volumes (1) at 70 kg. Ceftaroline clearance carries a power effect of BSA-normalised creatinine clearance below 80 mL/min/1.73 m^2 (replaced by a multiplicative end-stage renal disease effect in subjects on hemodialysis), a power effect of age centred at 36 years and a patient-versus-healthy multiplier; ceftaroline central volume carries a patient-versus-healthy multiplier. A replacement clearance applies during a hemodialysis session. An intramuscular depot (first-order ka fixed) is retained for the one intramuscular study. Combined additive and proportional residual error on each analyte, with one inter-individual random effect scaling the residual magnitude of both."
  reference <- paste(
    "Das S, Li J, Iaconis J, Zhou D, Stone GG, Yan JL, Melnick D.",
    "Ceftaroline fosamil doses and breakpoints for Staphylococcus aureus in complicated",
    "skin and soft tissue infections. J Antimicrob Chemother. 2019;74(2):425-431.",
    "doi:10.1093/jac/dky439.",
    "Final parameter estimates are Supplementary Table S3; the covariate-model forms",
    "(allometric exponents 0.75 / 1 at 70 kg, power and categorical forms) are the",
    "Supplementary Methods 'Covariate analysis'; demographics are Supplementary Table S2;",
    "the probability of target attainment results used for validation are main-text",
    "Tables 2 and 3 and Supplementary Tables S4 to S7.",
    sep = " "
  )
  vignette <- "Das_2019_ceftaroline"
  dosing <- c("depot", "central")
  paper_specific_etas <- c("etaRUV")

  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Baseline body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling of every clearance (exponent 0.75, fixed) and every volume (exponent 1, fixed) of both ceftaroline fosamil and ceftaroline, referenced to 70 kg (Supplementary Methods, Covariate analysis: 'the x values were fixed as 0.75 and 1.0 for the clearance and volume of distribution parameters ... A reference weight value of 70 kg was used'). Range 40-134 kg (Supplementary Table S2).",
      source_name = "WT"
    ),
    AGE = list(
      description = "Age at study entry.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect (AGE/36)^-0.23 on ceftaroline clearance, applied to every subject (Supplementary Table S3 row 'Age effect on CLc ~ (Age/36)'). The same continuous (AGE/36) form is used by the sibling model Riccobene_2016_ceftaroline from the same model lineage. Range 12-93 years.",
      source_name = "Age"
    ),
    CRCL = list(
      description = "Creatinine clearance normalised by body surface area (NCrCL).",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect (NCrCL/80)^0.54 on ceftaroline clearance (Supplementary Table S3 row 'NCrCL effect on CLc ~ (NCrCL/80)'), active only below the 80 mL/min/1.73 m^2 reference (encoded with min(CRCL, 80)) and only in subjects not on hemodialysis. Neither restriction is printed in the article; both are the form of the sibling models Riccobene_2016_ceftaroline and Riccobene_2017_ceftaroline from the same lineage, and the article's own probability-of-target-attainment tables confirm them (see the vignette: the unrestricted form under-predicts every published median fT>MIC for normal renal function, and applying the renal term on top of the ESRD multiplier would cut ESRD clearance about threefold, against the ESRD rows of Supplementary Tables S6 and S7). Range 6.7-467.4 (Supplementary Table S2).",
      source_name = "NCrCL"
    ),
    DIS_HEALTHY = list(
      description = "Healthy-participant indicator (1 = healthy volunteer, 0 = patient).",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (healthy volunteer)",
      notes = "Supplementary Table S3 prints 'Patient status effect on CLc' = 1.16 and 'Patient status effect on Vcc' = 1.59, multiplicative as theta^X (Supplementary Methods). The printed typical values CLc = 8.12 L/h and Vcc = 9.44 L are the healthy-volunteer reference and the multipliers move them to patients, i.e. the effects are applied as base^(1 - DIS_HEALTHY). The orientation follows the row labels, the Supplementary Methods ('theta_x ... the individual being categorised into the group different from the reference group') and the sibling model Riccobene_2017_ceftaroline; the article's patient PTA tables are reproduced within a few points under either orientation (the packaged one somewhat closer) and do not settle it (see vignette). Patients were those with cSSTI, community-acquired pneumonia or suspected infection; 267 of 951 subjects were healthy volunteers.",
      source_name = "patient status"
    ),
    RRT_HEMODIAL_STATUS = list(
      description = "Intermittent-hemodialysis treatment-status indicator (1 = subject with end-stage renal disease on a hemodialysis programme, 0 = not).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not on hemodialysis)",
      notes = "Multiplicative effect 0.34^RRT_HEMODIAL_STATUS on ceftaroline clearance between dialysis sessions (Supplementary Table S3 row 'ESRD effect on CLc'); it replaces the creatinine clearance term, as in the sibling models of the same lineage. The source label is ESRD; the end-stage renal disease subjects in the dataset came from the hemodialysis studies P903-18 and D3720C00012 (Supplementary Table S1), so the ESRD flag and the on-hemodialysis flag coincide.",
      source_name = "ESRD"
    ),
    RRT_HEMODIAL_ACTIVE = list(
      description = "Hemodialysis-active indicator (1 while a dialysis session is running, 0 otherwise).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no dialysis session running)",
      notes = "Time-varying. Supplementary Table S3 prints 'CLc, dialysis' = 11.18 (units not printed; L/h, as the other clearances). Encoded, as in Riccobene_2017_ceftaroline, as a replacement of ceftaroline clearance during the session (no covariate or eta terms); the article does not print the equation (see vignette Assumptions).",
      source_name = "(dialysis period flag, unnamed in the article)"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "ceftaroline fosamil",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "ceftaroline fosamil",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "ceftaroline fosamil",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    central_ceftaroline = list(
      analyte = "ceftaroline",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1_ceftaroline = list(
      analyte = "ceftaroline",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 951L,
    n_studies = 21L,
    age_range = "12-93 years (median 50)",
    weight_range = "40-134 kg (median 73)",
    sex_female_pct = 100 * 450 / 951,
    race_ethnicity = "Caucasian 62.6%, Asian 22.2%, Black 8.2%, Indian 5.3%, other 1.8% (Supplementary Table S2)",
    disease_state = "267 healthy volunteers (including subjects with renal impairment or end-stage renal disease from the phase 1 studies) and 684 patients: 463 with complicated skin and soft tissue infection, 214 with community-acquired pneumonia and 7 with suspected infection.",
    dose_range = "Ceftaroline fosamil single and multiple doses by intravenous infusion, including 600 mg q12h by 1 h and 600 mg q8h by 2 h infusion, and renally adjusted doses; one phase 1 study gave intramuscular doses.",
    regions = "Multinational (14 phase 1, 1 phase 2 and 6 phase 3 studies; Supplementary Table S1)",
    n_observations = "2575 ceftaroline fosamil and 8174 ceftaroline plasma concentrations",
    notes = "Height median 171 cm (114-196), BSA median 1.8 m^2 (1.2-2.6), NCrCL median 95.8 mL/min/1.73 m^2 (6.7-467.4) (Supplementary Table S2). Fitted with FOCE-I in NONMEM 7.2.0. Concentrations were modelled in micromol/L; this implementation doses in mg and reports mg/L (see vignette)."
  )

  ini({
    # ----- Ceftaroline fosamil (prodrug) -- Supplementary Table S3, typical
    # values at the 70 kg allometric reference.
    lcl <- log(268.7); label("Ceftaroline fosamil clearance CLcf at 70 kg (L/h)")                    # Suppl Table S3: CLcf = 268.7 L/h (RSE 4.2%)
    lvc <- log(5.6);   label("Ceftaroline fosamil central volume Vccf at 70 kg (L)")                 # Suppl Table S3: Vccf = 5.6 L (RSE 6.21%)
    lq  <- log(20.39); label("Ceftaroline fosamil intercompartmental clearance Qcf at 70 kg (L/h)")  # Suppl Table S3: Qcf = 20.39 L/h (RSE 4.26%)
    lvp <- log(8);     label("Ceftaroline fosamil peripheral volume Vpcf at 70 kg (L)")              # Suppl Table S3: Vpcf = 8 L (RSE 4.82%)
    lka <- fixed(log(0.53)); label("Ceftaroline fosamil intramuscular absorption rate Kacf (1/h)")   # Suppl Table S3: Kacf = 0.53 h-1 (FIX)

    # ----- Ceftaroline (active moiety) -- Supplementary Table S3. CLc and Vcc
    # are the healthy-volunteer typical values (see covariateData$DIS_HEALTHY).
    lcl_ceftaroline <- log(8.12); label("Ceftaroline clearance CLc at 70 kg in healthy subjects (L/h)")     # Suppl Table S3: CLc = 8.12 L/h (RSE 2.72%)
    lvc_ceftaroline <- log(9.44); label("Ceftaroline central volume Vcc at 70 kg in healthy subjects (L)")  # Suppl Table S3: Vcc = 9.44 L (RSE 2.25%)
    lq_ceftaroline  <- log(8.39); label("Ceftaroline intercompartmental clearance Qc at 70 kg (L/h)")       # Suppl Table S3: Qc = 8.39 L/h (RSE 1.49%)
    lvp_ceftaroline <- log(9.58); label("Ceftaroline peripheral volume Vpc at 70 kg (L)")                   # Suppl Table S3: Vpc = 9.58 L (RSE 0.97%)

    lcl_hemodialysis_ceftaroline <- log(11.18); label("Ceftaroline clearance during a hemodialysis session (L/h)")  # Suppl Table S3: 'CLc, dialysis' = 11.18 (RSE 0.73%; unit not printed)

    # ----- Covariate effects -- Supplementary Methods and Table S3.
    e_wt_cl_q  <- fixed(0.75); label("Allometric exponent shared by all clearances (unitless)")  # Suppl Methods: exponent fixed as 0.75 for clearance, reference 70 kg
    e_wt_vc_vp <- fixed(1);    label("Allometric exponent shared by all volumes (unitless)")     # Suppl Methods: exponent fixed as 1.0 for volume, reference 70 kg

    e_age_cl_ceftaroline  <- -0.23; label("Exponent of (AGE/36) on ceftaroline CL")                                  # Suppl Table S3: 'Age effect on CLc ~ (Age/36)' = -0.23 (RSE 9.74%)
    e_crcl_cl_ceftaroline <- 0.54;  label("Exponent of (CRCL/80) on ceftaroline CL below 80 mL/min/1.73 m^2")       # Suppl Table S3: 'NCrCL effect on CLc ~ (NCrCL/80)' = 0.54 (RSE 5.17%)
    e_rrt_hemodial_status_cl_ceftaroline <- 0.34; label("ESRD (hemodialysis) multiplier on interdialytic ceftaroline CL")  # Suppl Table S3: 'ESRD effect on CLc' = 0.34 (RSE 5.14%)
    e_dis_healthy_cl_ceftaroline <- 1.16; label("Patient multiplier on ceftaroline CL, applied as base^(1 - DIS_HEALTHY)")  # Suppl Table S3: 'Patient status effect on CLc' = 1.16 (RSE 3.07%)
    e_dis_healthy_vc_ceftaroline <- 1.59; label("Patient multiplier on ceftaroline Vc, applied as base^(1 - DIS_HEALTHY)")  # Suppl Table S3: 'Patient status effect on Vcc' = 1.59 (RSE 2.43%)

    # ----- Inter-individual variability -- Supplementary Table S3 prints each
    # as a CV; carried as omega^2 = CV^2 (log-normal etas).
    etalcl ~ 0.8464           # Suppl Table S3: IIV CLcf CV 0.92 -> 0.92^2
    etalvc ~ 0.8464           # Suppl Table S3: IIV Vccf CV 0.92 -> 0.92^2
    etalka ~ fixed(0.1089)    # Suppl Table S3: IIV Kacf CV 0.33 (FIX) -> 0.33^2
    # Block: VAR(CLc) = 0.26^2 = 0.0676; COV = 0.65 * 0.26 * 0.24 = 0.04056
    # (Suppl Table S3 'IIV Vcc:CLc' correlation 0.65); VAR(Vcc) = 0.24^2 = 0.0576
    etalcl_ceftaroline + etalvc_ceftaroline ~ c(
      0.0676,
      0.04056, 0.0576
    )
    etalvp_ceftaroline ~ 0.009025 # Suppl Table S3: IIV Vpc CV 0.095 -> 0.095^2
    etaRUV ~ 0.36                 # Suppl Table S3: 'IIV Eta on Epsilon' CV 0.6 -> 0.6^2; scales the residual SD of both analytes

    # ----- Residual error -- Supplementary Table S3. Additive SDs were printed
    # in micromol/L and are converted to mg/L with the molar masses used in
    # model() (ceftaroline fosamil 684.7, ceftaroline 604.7 g/mol).
    propSd             <- 0.29;      label("Ceftaroline fosamil proportional residual SD (fraction)")  # Suppl Table S3: proportional residual error cf = 0.29
    addSd              <- 0.0821640; label("Ceftaroline fosamil additive residual SD (mg/L)")          # Suppl Table S3: additive residual error cf = 0.12 umol/L x 0.6847 mg/umol
    propSd_ceftaroline <- 0.18;      label("Ceftaroline proportional residual SD (fraction)")          # Suppl Table S3: proportional residual error c = 0.18
    addSd_ceftaroline  <- 0.0380961; label("Ceftaroline additive residual SD (mg/L)")                  # Suppl Table S3: additive residual error c = 0.063 umol/L x 0.6047 mg/umol
  })

  model({
    # Molar masses (g/mol) for the 100% molar conversion of ceftaroline fosamil
    # to ceftaroline; the model was fitted in micromol/L.
    mw_fosamil     <- 684.7
    mw_ceftaroline <- 604.7

    ref_wt   <- 70 # kg, allometric reference (Supplementary Methods)
    ref_age  <- 36 # years (Supplementary Table S3)
    ref_crcl <- 80 # mL/min/1.73 m^2 (Supplementary Table S3)

    allom_cl <- (WT / ref_wt)^e_wt_cl_q
    allom_v  <- (WT / ref_wt)^e_wt_vc_vp

    # Renal term: saturates at the 80 mL/min/1.73 m^2 reference and is
    # replaced by the ESRD multiplier in subjects on hemodialysis.
    crcl_cl <- (min(CRCL, ref_crcl) / ref_crcl)^(e_crcl_cl_ceftaroline * (1 - RRT_HEMODIAL_STATUS))
    age_cl  <- (AGE / ref_age)^e_age_cl_ceftaroline

    # ----- Ceftaroline fosamil
    cl <- exp(lcl + etalcl) * allom_cl
    vc <- exp(lvc + etalvc) * allom_v
    q  <- exp(lq) * allom_cl
    vp <- exp(lvp) * allom_v
    ka <- exp(lka + etalka)

    # ----- Ceftaroline. Clearance is replaced by the dialysis value while a
    # hemodialysis session runs.
    cl_ceftaroline <-
      exp(lcl_ceftaroline + etalcl_ceftaroline) * allom_cl * crcl_cl * age_cl *
        e_rrt_hemodial_status_cl_ceftaroline^RRT_HEMODIAL_STATUS *
        e_dis_healthy_cl_ceftaroline^(1 - DIS_HEALTHY) *
        (1 - RRT_HEMODIAL_ACTIVE) +
      exp(lcl_hemodialysis_ceftaroline) * RRT_HEMODIAL_ACTIVE
    vc_ceftaroline <- exp(lvc_ceftaroline + etalvc_ceftaroline) * allom_v *
      e_dis_healthy_vc_ceftaroline^(1 - DIS_HEALTHY)
    q_ceftaroline  <- exp(lq_ceftaroline) * allom_cl
    vp_ceftaroline <- exp(lvp_ceftaroline + etalvp_ceftaroline) * allom_v

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    kel_ceftaroline <- cl_ceftaroline / vc_ceftaroline
    k12_ceftaroline <- q_ceftaroline / vc_ceftaroline
    k21_ceftaroline <- q_ceftaroline / vp_ceftaroline

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # All ceftaroline fosamil cleared becomes ceftaroline, mol for mol.
    d/dt(central_ceftaroline) <-
      kel * central * mw_ceftaroline / mw_fosamil -
      kel_ceftaroline * central_ceftaroline -
      k12_ceftaroline * central_ceftaroline + k21_ceftaroline * peripheral1_ceftaroline
    d/dt(peripheral1_ceftaroline) <-
      k12_ceftaroline * central_ceftaroline - k21_ceftaroline * peripheral1_ceftaroline

    Cc             <- central / vc
    Cc_ceftaroline <- central_ceftaroline / vc_ceftaroline

    # One eta scales the residual SD of both analytes (Suppl Table S3
    # 'IIV Eta on Epsilon').
    ruv_scale <- exp(etaRUV)
    addSd_i              <- addSd * ruv_scale
    propSd_i             <- propSd * ruv_scale
    addSd_ceftaroline_i  <- addSd_ceftaroline * ruv_scale
    propSd_ceftaroline_i <- propSd_ceftaroline * ruv_scale

    Cc             ~ add(addSd_i) + prop(propSd_i)
    Cc_ceftaroline ~ add(addSd_ceftaroline_i) + prop(propSd_ceftaroline_i)
  })
}
