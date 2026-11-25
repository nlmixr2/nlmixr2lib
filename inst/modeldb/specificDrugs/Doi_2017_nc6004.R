Doi_2017_nc6004 <- function() {
  description <- paste(
    "Population PK/PD model for NC-6004 (cisplatin-incorporating polymeric",
    "micelles) in adults with advanced solid tumours (Doi 2017): a",
    "two-compartment model for micelle-associated (total plasma) platinum",
    "whose clearance releases platinum into a two-compartment free",
    "(ultrafilterable) platinum model, allometric body-weight effects on the",
    "micellar clearance and central volume, and an indirect-response model in",
    "which free platinum inhibits the production of eGFR with a maximum",
    "inhibition that decays exponentially with time since the first dose and",
    "varies between dosing occasions."
  )
  reference <- paste(
    "Doi T, Hamaguchi T, Shitara K, Takahashi N, Kato K, Saitoh H, Yamaguchi",
    "T, Kitano T, Kita Y, Ohtsu A. NC-6004 Phase I study in combination with",
    "gemcitabine for advanced solid tumors and population PK/PD analysis.",
    "Cancer Chemother Pharmacol. 2017;79(3):569-578.",
    "doi:10.1007/s00280-017-3254-4"
  )
  vignette <- "Doi_2017_nc6004"
  units <- list(
    time = "h",
    dosing = paste(
      "mg of platinum. NC-6004 doses are labelled as cisplatin-equivalent",
      "mg/m^2; convert with x 195.08 / 300.05 (= 0.6502, Pt / cisplatin molar",
      "mass) and multiply by BSA before dosing into central_np"
    ),
    concentration = paste(
      "mg/L (= ug/mL) of platinum for Cc_np (total plasma Pt) and Cc (free",
      "Pt); eGFR in mL/min/1.73 m^2"
    )
  )

  compartmentData <- list(
    central_np = list(
      analyte = "platinum in NC-6004 micelles (total plasma Pt)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1_np = list(analyte = "platinum in NC-6004 micelles", units = "mg", specimen = "tissue", verified = TRUE),
    central = list(
      analyte = "free (ultrafilterable) platinum",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(analyte = "free platinum", units = "mg", specimen = "tissue", verified = TRUE),
    egfr = list(
      analyte = "estimated glomerular filtration rate",
      units = "mL/min/1.73 m^2",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effects (WT / 67.8)^0.731 on the micellar clearance CL and",
        "(WT / 67.8)^0.763 on the micellar central volume V1 (Table 3).",
        "67.8 kg is the typical-value weight named in the Table 3 footnote."
      ),
      source_name = "WT"
    ),
    ASSAY_AAS = list(
      description = paste(
        "Total-platinum assay indicator: 1 = atomic absorption spectrometry",
        "(AAS), 0 = inductively coupled plasma (ICP) method"
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ICP); the structural model is on the ICP scale",
      notes = paste(
        "Per observation. Applies to the total-platinum output Cc_np only.",
        "The UK study NC-6004-001 measured total Pt by AAS in 17 patients",
        "(144 points) and by ICP in 5 (42 points); the Japanese study used ICP",
        "throughout. AAS rows read 1.69 x exp(eta) times the ICP prediction",
        "with a 31.6% instead of 17.0% proportional residual error (Table 3",
        "'Difference by method [folds] (AA to ICP)'). Set ASSAY_AAS = 0 to",
        "simulate the ICP readout."
      ),
      source_name = "method (AA / ICP)"
    ),
    OCC = list(
      description = paste(
        "Dosing occasion: k for the interval from the k-th to the (k+1)-th",
        "NC-6004 administration"
      ),
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Carries the inter-occasion variability on Imax (Methods,",
        "'Population PK/PD analysis': occasion = the period through the k-th",
        "to the (k+1)-th administration). Increment OCC on every dose record",
        "and keep it constant until the next dose. Consumed by the ini()",
        "nesting syntax 'etaiov_limax ~ ... | OCC'."
      ),
      source_name = "occasion k"
    )
  )

  covariatesDataExcluded <- list(
    RACE_JAPANESE = list(
      description = "Japanese (1) vs Caucasian (0) patient population",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Caucasian, UK study)",
      notes = paste(
        "Screened on CL and not retained: dOFV 0.014 (P = 0.906) (Results,",
        "'Population PK analysis')."
      ),
      source_name = "ethnic difference"
    ),
    CRCL = list(
      description = "Baseline eGFR",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened as a PK covariate on clearance and volume and not retained",
        "(Methods, 'Population PK analysis'; only body weight was",
        "significant)."
      ),
      source_name = "baseline eGFR"
    ),
    CONMED_IV_FLUIDS = list(
      description = "Prophylactic hydration therapy around NC-6004 dosing",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no hydration)",
      notes = paste(
        "Evaluated preliminarily on the eGFR model: dOFV -2.556 (P < 0.110);",
        "not retained because only five patients received it (Results,",
        "'Population PK/PD analysis')."
      ),
      source_name = "hydration therapy"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 29L,
    n_studies = 2L,
    age_range = "34-72 years (Japanese study, Table 1); 40-80 years (UK study NC-6004-001)",
    weight_range = "Not reported; typical-value reference weight 67.8 kg (Table 3 footnote)",
    sex_female_pct = 31.0,
    race_ethnicity = "12 Japanese (NC-6004-003) and 17 Caucasian (UK NC-6004-001) patients",
    disease_state = paste(
      "Advanced, treatment-refractory solid tumours (Japanese study: 7",
      "carcinoma, 5 neuroendocrine tumour; ECOG 0-1; eGFR >= 60 mL/min/1.73",
      "m^2 at entry)"
    ),
    dose_range = paste(
      "NC-6004 10-120 mg/m^2 (cisplatin equivalent) as a 1-h intravenous",
      "infusion every 21 days; Japanese study 60 or 90 mg/m^2 with",
      "gemcitabine 1000 mg/m^2 on days 1 and 8 from cycle 2"
    ),
    regions = "Japan; United Kingdom",
    notes = paste(
      "PK: 546 concentrations from 29 patients (12 Japanese, 17 UK); free",
      "Pt in 28 patients (252 points). PK/PD (eGFR): 12 Japanese patients",
      "(mean 8.1 points each) and 17 UK patients (mean 13.7 points each).",
      "Sex: 10 M / 2 F (Japanese, Table 1) and 10 M / 7 F (UK, Plummer 2011",
      "Table 2), so 9 / 29 = 31.0% female. NONMEM 7.3.0, FOCE-I."
    )
  )

  ini({
    # Micellar (NC-6004) platinum -- Table 3. IIV CV% on the SD scale
    # (omega = CV / 100): every Table 3 90% CI is symmetric in CV^2
    # (e.g. Qm 54.4 (22.2-73.7): 0.2959 vs CI-midpoint 0.2962).
    lvc_np <- log(2.69)
    label("Central volume of micellar NC-6004 platinum at 67.8 kg, V1 (L)") # Table 3, TV1 = 2.69 L (typical value at 67.8 kg)
    e_wt_vc_np <- 0.763
    label("Power exponent of (WT / 67.8) on V1 (unitless)") # Table 3, Theta(V1-WT) = 0.763
    lvp_np <- log(3.91)
    label("Peripheral volume of micellar NC-6004 platinum, V2 (L)") # Table 3, V2 = 3.91 L (no IIV)
    lcl_np <- log(0.0687)
    label("Clearance of micellar NC-6004 platinum (release into free Pt) at 67.8 kg, CL (L/h)") # Table 3, CL = 0.0687 L/h
    e_wt_cl_np <- 0.731
    label("Power exponent of (WT / 67.8) on CL (unitless)") # Table 3, Theta(CL-WT) = 0.731
    lq_np <- log(0.0117)
    label("Intercompartmental clearance of micellar NC-6004 platinum, Q (L/h)") # Table 3, Q = 0.0117 L/h

    # Free (ultrafilterable) platinum -- Table 3
    lvc <- log(81.6)
    label("Central volume of free platinum, Vm1 (L)") # Table 3, Vm1 = 81.6 L
    lvp <- log(665)
    label("Peripheral volume of free platinum, Vm2 (L)") # Table 3, Vm2 = 665 L
    lcl <- log(6.54)
    label("Clearance of free platinum, CLm (L/h)") # Table 3, CLm = 6.54 L/h
    lq <- log(9.30)
    label("Intercompartmental clearance of free platinum, Qm (L/h)") # Table 3, Qm = 9.30 L/h

    # Total-platinum assay conversion -- Table 3
    e_assay_aas_cc <- 1.69
    label("Multiplicative factor on predicted total Pt for AAS- vs ICP-assayed samples (unitless)") # Table 3, 'Difference by method [folds] (AA to ICP)' = 1.69
    e_assay_aas_propsd <- 1.858824
    label("Multiplicative factor on the total-Pt proportional residual SD for AAS- vs ICP-assayed samples (unitless)") # Table 3, residual CV Total Pt (AA) 31.6% / Total Pt (ICP) 17.0% = 0.316 / 0.170

    # eGFR indirect-response model -- Table 4 (rates printed per day,
    # converted to per hour for the hour time base)
    lrbase <- log(83.4)
    label("Baseline eGFR, GFR0 (mL/min/1.73 m^2)") # Table 4, GFR0 = 83.4
    lkout <- log(0.171 / 24)
    label("First-order loss rate constant of eGFR, Kout (1/h)") # Table 4, Kout = 0.171 /day; / 24
    lic50 <- log(35.3)
    label("Free-Pt concentration for half-maximal inhibition of eGFR production, IC50 (ng/mL)") # Table 4, IC50 = 35.3 ng/mL
    limax <- log(0.851)
    label("Maximum fractional inhibition of eGFR production at the first dose, Imax0 (fraction)") # Table 4, Imax0 = 0.851
    limax_exp_kdes <- log(0.0273 / 24)
    label("First-order rate of decline of Imax with time since the first dose, trend (1/h)") # Table 4, 'Imax decrease trend (/day)' = 0.0273; / 24

    etalvc_np ~ 0.019321 # Table 3, V1 IIV 13.9 CV% (90% CI 7.3-18.2); 0.139^2
    etalcl_np ~ 0.029584 # Table 3, CL IIV 17.2 CV% (8.9-22.7); 0.172^2
    etalq_np ~ 0.040401 # Table 3, Q IIV 20.1 CV% (0-28.5); 0.201^2
    etalvc ~ 0.059536 # Table 3, Vm1 IIV 24.4 CV% (0-36.5); 0.244^2
    etalvp ~ 0.179776 # Table 3, Vm2 IIV 42.4 CV% (24.2-54.9); 0.424^2
    etalcl ~ 0.208849 # Table 3, CLm IIV 45.7 CV% (21.8-60.9); 0.457^2
    etalq ~ 0.295936 # Table 3, Qm IIV 54.4 CV% (22.2-73.7); 0.544^2
    etae_assay_aas_cc ~ 0.044944 # Table 3, AA-to-ICP fold IIV 21.2 CV% (4.8-29.6); 0.212^2
    etalrbase ~ 0.070225 # Table 4, GFR0 IIV 26.5 CV%; 0.265^2
    etalkout ~ 0.237169 # Table 4, Kout IIV 48.7 CV%; 0.487^2
    etalic50 ~ 4.481689 # Table 4, IC50 IIV 211.7 CV%; 2.117^2
    etaiov_limax ~ 0.281961 | OCC # Table 4, IOV (Imax) 53.1 CV%; 0.531^2

    propSd_np <- 0.170
    label("Proportional residual error, total Pt by ICP (fraction)") # Table 3, Total Pt (ICP) 17.0 CV% (15.2-18.7)
    propSd <- 0.316
    label("Proportional residual error, free Pt (fraction)") # Table 3, Free Pt (ICP) 31.6 CV% (24.0-37.7)
    propSd_egfrObs <- 0.088
    label("Proportional residual error, eGFR (fraction)") # Table 4, intra-individual variance 8.8 CV%
  })
  model({
    # Micellar NC-6004 platinum, with allometric body-weight effects
    # (Table 3: V1 = TV1 x (WT/67.8)^Theta, CL = TCL x (WT/67.8)^Theta)
    vc_np <- exp(lvc_np + etalvc_np) * (WT / 67.8)^e_wt_vc_np
    cl_np <- exp(lcl_np + etalcl_np) * (WT / 67.8)^e_wt_cl_np
    vp_np <- exp(lvp_np)
    q_np <- exp(lq_np + etalq_np)

    # Free platinum
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    cl <- exp(lcl + etalcl)
    q <- exp(lq + etalq)

    # eGFR turnover (Methods, Eq. 1; kin = eGFR0 x Kout)
    rbase <- exp(lrbase + etalrbase)
    kout <- exp(lkout + etalkout)
    kin <- rbase * kout
    ic50 <- exp(lic50 + etalic50)
    imax <- exp(limax)
    imax_exp_kdes <- exp(limax_exp_kdes)

    # Methods, Eq. 2: Imax_ijk = Imax0 x exp(-trend x t_j) x exp(eta_ik^IOV),
    # t_j = time of measurement from the first administration
    tfirst <- max(0, tafd())
    imax_t <- imax * exp(-imax_exp_kdes * tfirst) * exp(etaiov_limax)

    k12_np <- q_np / vc_np
    k21_np <- q_np / vp_np
    k12 <- q / vc
    k21 <- q / vp

    # The micellar clearance releases platinum into the free-Pt central
    # compartment (no other source of free Pt is parameterised)
    d/dt(central_np) <- -cl_np / vc_np * central_np - k12_np * central_np + k21_np * peripheral1_np
    d/dt(peripheral1_np) <- k12_np * central_np - k21_np * peripheral1_np
    d/dt(central) <- cl_np / vc_np * central_np - cl / vc * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Free Pt in ng/mL drives the inhibition (IC50 is in ng/mL)
    Cc <- central / vc
    cfree_ngml <- Cc * 1000
    inh <- imax_t * cfree_ngml / (ic50 + cfree_ngml)

    egfr(0) <- rbase
    d/dt(egfr) <- kin * (1 - inh) - kout * egfr

    # Total plasma Pt on the ICP scale; AAS rows read e_assay_aas_cc x exp(eta)
    # times the ICP value with a larger proportional residual error
    cf_aas <- e_assay_aas_cc * exp(etae_assay_aas_cc)
    Cc_np <- central_np / vc_np * cf_aas^ASSAY_AAS
    propSd_np_eff <- propSd_np * e_assay_aas_propsd^ASSAY_AAS

    # Observed eGFR is the turnover state itself
    egfrObs <- egfr

    Cc_np ~ prop(propSd_np_eff)
    Cc ~ prop(propSd)
    egfrObs ~ prop(propSd_egfrObs)
  })
}
