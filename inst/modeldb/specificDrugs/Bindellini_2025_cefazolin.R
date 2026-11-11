Bindellini_2025_cefazolin <- function() {
  description <- "Two-compartment population PK model for intravenous cefazolin in obese (bariatric surgery, BMI >= 35 kg/m^2) and nonobese surgical adults with plasma and subcutaneous adipose-tissue microdialysis data (Bindellini 2025). The model is parameterised on UNBOUND cefazolin: Cu = central/V1 with linear unbound clearance, and the measured total plasma concentration is Cc = Cu + Bmax * Cu / (Kd + Cu) (saturable protein binding). V1 scales linearly and Q to the 0.75 power with lean body weight (fat-free mass, Janmahasatian); V2 follows the LBW/FM approach, with a fraction R scaling with lean body weight and the rest with fat mass. The interstitial-fluid concentration is a tissue factor times the peripheral-compartment concentration, Cisf = TF * C2; microdialysate concentrations are Cisf times the catheter relative recovery RR, which differs between obese and nonobese patients and carries inter-individual and inter-catheter variability, and retrodialysis observes 1 - RR."
  reference <- "Bindellini D, Simon P, Busse D, Michelet R, Petroff D, Aulin LBS, Dorn C, Zeitlinger M, Huisinga W, Wrigge H, Kloft C. Evaluation of the need for dosing adaptations in obese patients for surgical antibiotic prophylaxis: a model-based analysis of cefazolin pharmacokinetics. Br J Anaesth. 2025;134(4):1041-1049. doi:10.1016/j.bja.2024.11.044. PMCID PMC11947583. Parameter estimates from Table 2; model structure from Figure 1 and the Results text; the authors' deposited mrgsolve simulation code (github.com/Kloft-Lab/Bindellini-et-al._CefazolinPK_Obese_Nonobese_mrgsolve_model, popPK_CEZ_run077.cpp, cited in the Supplementary Material) confirms the LBW/FM scaling, the binding equation and Cisf = TF * C2."
  vignette <- "Bindellini_2025_cefazolin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(
      analyte = "cefazolin",
      units = "mg",
      specimen = "plasma",
      verified = TRUE,
      notes = "Unbound cefazolin amount in the central compartment: the model is parameterised in terms of unbound concentrations (Methods), so central / V1 is the UNBOUND plasma concentration Cu, and the total plasma concentration Cc is reconstructed algebraically through the saturable binding equation."
    ),
    peripheral1 = list(
      analyte = "cefazolin",
      units = "mg",
      specimen = "tissue",
      verified = TRUE,
      notes = "Unbound cefazolin amount in the peripheral compartment V2. The paper attributes microdialysate concentrations to this compartment: the subcutaneous adipose interstitial-fluid concentration is Cisf = TF * peripheral1 / V2 (Results; deposited mrgsolve code CISF = (P/V2)*TF)."
    )
  )

  covariateData <- list(
    FFM = list(
      description = "Lean body weight by the Janmahasatian 2005 equation (fat-free mass)",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Reference 64.3 kg, the median lean body weight of all 30 patients (Table 1; also the reference printed in the Table 2 rows for V1, Q and V2). Enters V1 linearly (exponent 1 fixed), Q with exponent 0.75 fixed, and the lean-body-weight part of V2 (fraction f_ffm_vp). Bindellini 2025 cites Janmahasatian et al. (Clin Pharmacokinet 2005) for LBW, the fat-free-mass equation, so it maps to the canonical FFM. Cohort range 34.2-96.1 kg (nonobese median 45.6, obese median 69.3).",
      source_name = "LBW"
    ),
    FM = list(
      description = "Fat mass, total body weight minus Janmahasatian lean body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Reference 41.5 kg, the median fat mass of all 30 patients (Table 1; also the reference printed in the Table 2 V2 row). Enters only the fat-mass part of V2 (fraction 1 - f_ffm_vp = 23.6%) linearly, per the LBW/FM approach (Huisinga et al. 2012) the paper cites. Compute as FM = WT - FFM. Cohort range 15.8-122 kg (nonobese median 25.0, obese median 84.2).",
      source_name = "FM"
    ),
    DIS_OBESE = list(
      description = "Obese study-cohort indicator (1 = obese bariatric-surgery cohort, 0 = nonobese cohort)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Bindellini 2025 enrolled 15 obese patients (BMI >= 35 kg/m^2) undergoing bariatric surgery and 15 age- and sex-matched nonobese patients (18.5 <= BMI <= 30 kg/m^2) undergoing elective abdominal surgery; observed obese BMI 39.5-69.3 kg/m^2. The indicator acts only on the microdialysis relative recovery (RR 23.3% obese vs 41.1% nonobese, Table 2), a microdialysis-technique and tissue property; it does not change any drug parameter, so simulations of plasma and interstitial-fluid exposure are independent of it.",
      source_name = "OBE"
    ),
    MD_CATHETER = list(
      description = "Microdialysis catheter identifier within a patient (1 or 2)",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = "Two microdialysis catheters per patient, one in the subcutaneous adipose tissue of each upper arm (Bindellini 2025 Methods; Figure 1 catheters 1 and 2). Consumed by the ini() nesting syntax 'etaiov_lrr_md ~ ... | MD_CATHETER', which gives each catheter its own relative-recovery random effect (the paper's intercatheter variability, 58.6% CV). Both catheters are sampled over the same collection intervals, so records of the two catheters interleave in time. It has no effect on Cc, Cu or Cisf, so plasma-only records may carry either value, but rxode2 requires an event table to contain at least two distinct MD_CATHETER values (also when the random effects are zeroed).",
      source_name = "catheter"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 30,
    n_studies = 1,
    age_range = "21-65 years (median 45)",
    weight_range = "50.0-200 kg total body weight (nonobese median 78.0, obese median 155)",
    sex_female_pct = 66.7,
    race_ethnicity = "Not reported (single-centre German cohort)",
    disease_state = "Adults undergoing surgery under general anaesthesia: 15 obese patients (BMI >= 35 kg/m^2, median 52.6, range 39.5-69.3) undergoing bariatric surgery and 15 age- and sex-matched nonobese patients (BMI 18.5-30 kg/m^2, median 26.0, range 18.7-29.8) undergoing elective abdominal surgery.",
    renal_function = "eGFR (CKD-EPI) median 119 mL/min (range 82.8-269); nonobese median 103.8, obese median 133.",
    dose_range = "Single 2 g cefazolin dose as a 30-minute intravenous infusion after induction of anaesthesia (surgical prophylaxis).",
    regions = "Germany (Leipzig University; EudraCT 2012-004383-22, DRKS00004776)",
    notes = "Post hoc analysis of the Dorn 2021 study. Plasma total cefazolin sampled at 0.5, 1, 2, 3, 4, 5, 6 and 8 h, unbound (ultrafiltration) at 0.5, 1, 4 and 8 h; microdialysate (CMA 63 catheters, 2 uL/min, both upper arms) collected over 0-0.5, 0.5-1, 1-1.5, 1.5-2, 2-3, 3-4, 4-5, 5-6, 6-7 and 7-8 h; three 15-minute retrodialysis samples per patient after the sampling period. Lean body weight median 64.3 kg (range 34.2-96.1), fat mass median 41.5 kg (range 15.8-122). Fitted in NONMEM 7.4.3; parameter precision by sampling importance resampling (Table 2)."
  )

  ini({
    # Structural disposition parameters, all referenced to UNBOUND
    # cefazolin (Bindellini 2025 Table 2, fixed-effects block). V1, Q and
    # V2 are typical values at LBW (FFM) = 64.3 kg and FM = 41.5 kg.
    lcl <- log(17.9); label("Unbound clearance CL (L/h)") # Table 2 CL = 17.9 L/h (SIR 95% CI 16.6-19.2)
    lvc <- log(22.9); label("Unbound central volume V1 at FFM 64.3 kg (L)") # Table 2 V1 = 22.9 L (SIR 95% CI 19.2-26.5)
    lq <- log(56.8); label("Intercompartmental clearance Q at FFM 64.3 kg (L/h)") # Table 2 Q = 56.8 L/h (SIR 95% CI 49.2-66.9)
    lvp <- log(34.3); label("Unbound peripheral volume V2 at FFM 64.3 kg and FM 41.5 kg (L)") # Table 2 V2 = 34.3 L (SIR 95% CI 32.4-36.9)

    # Body-size scaling (Table 2; fixed exponents marked with the dagger)
    e_ffm_vc <- fixed(1); label("Exponent of (FFM/64.3) on V1 (unitless)") # Table 2 'Exponent V1' = 1, fixed
    e_ffm_q <- fixed(0.75); label("Exponent of (FFM/64.3) on Q (unitless)") # Table 2 'Exponent Q' = 0.75, fixed
    f_ffm_vp <- 0.764; label("Fraction of V2 scaling with lean body weight; 1 - f_ffm_vp scales with fat mass (fraction)") # Table 2 R = 76.4 percent (SIR 95% CI 64.9-86.0)

    # Saturable plasma protein binding (Table 2; Results)
    lbmax_pb <- log(247); label("Maximum binding capacity Bmax (mg/L)") # Table 2 Bmax = 247 mg/L (SIR 95% CI 219-289)
    lkd_pb <- log(65.3); label("Dissociation constant Kd (mg/L)") # Table 2 Kd = 65.3 mg/L (SIR 95% CI 55.8-79.9)

    # Interstitial-fluid tissue factor (Table 2; Results)
    ratio_isf_peripheral1 <- 0.655; label("Tissue factor TF, ISF concentration / peripheral-compartment concentration (fraction)") # Table 2 TF = 65.5 percent (SIR 95% CI 59.4-70.9)

    # Microdialysis relative recovery (Table 2; Results)
    lrr_md <- log(0.411); label("Relative recovery RR of the microdialysis catheter, nonobese patients (fraction)") # Table 2 RRnonobese = 41.1 percent (SIR 95% CI 35.3-48.8)
    lrr_md_obese <- log(0.233); label("Relative recovery RR of the microdialysis catheter, obese patients (fraction)") # Table 2 RRobese = 23.3 percent (SIR 95% CI 19.4-29.1)

    # Inter-individual variability (Table 2, CV%), converted to log-scale
    # variances with omega^2 = log(CV^2 + 1). Kd carries no IIV.
    etalcl ~ 0.046846 # Table 2 IIV CL = 21.9 percent CV
    etalvc ~ 0.14911 # Table 2 IIV V1 = 40.1 percent CV
    etalq ~ 0.2626 # Table 2 IIV Q = 54.8 percent CV
    etalvp ~ 0.015259 # Table 2 IIV V2 = 12.4 percent CV
    etalbmax_pb ~ 0.0084284 # Table 2 IIV Bmax = 9.20 percent CV
    etalrr_md ~ 0.17184 # Table 2 IIV RR = 43.3 percent CV

    # Intercatheter variability on RR (Table 2, 58.6% CV), one random
    # effect per catheter within a patient
    etaiov_lrr_md ~ 0.2952 | MD_CATHETER # Table 2 Intercatheter = 58.6 percent CV

    # Residual error (Table 2, proportional). One plasma error applies to
    # both total and unbound plasma observations; nlmixr2 needs one
    # parameter per output, so propSd and propSd_Cu carry the same value.
    propSd <- 0.357; label("Proportional residual error, total plasma Cc (fraction)") # Table 2 plasma proportional error = 35.7 percent
    propSd_Cu <- 0.357; label("Proportional residual error, unbound plasma Cu (fraction)") # Table 2 plasma proportional error = 35.7 percent (shared with total plasma)
    propSd_Cmd <- 0.521; label("Proportional residual error, microdialysate Cmd (fraction)") # Table 2 microdialysis proportional error = 52.1 percent
    propSd_Crd_ratio <- 0.289; label("Proportional residual error, retrodialysis ratio Crd_ratio (fraction)") # Table 2 retrodialysis proportional error = 28.9 percent
  })

  model({
    # Individual disposition parameters (Table 2; deposited mrgsolve code)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc) * (FFM / 64.3)^e_ffm_vc
    q <- exp(lq + etalq) * (FFM / 64.3)^e_ffm_q
    vp <- exp(lvp + etalvp) * (f_ffm_vp * FFM / 64.3 + (1 - f_ffm_vp) * FM / 41.5)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Plasma: unbound concentration and total via saturable binding
    bmax_pb <- exp(lbmax_pb + etalbmax_pb)
    kd_pb <- exp(lkd_pb)
    Cu <- central / vc
    Cc <- Cu + bmax_pb * Cu / (kd_pb + Cu)

    # Adipose interstitial fluid: Cisf = TF * C2
    Cisf <- ratio_isf_peripheral1 * peripheral1 / vp

    # Microdialysis catheter relative recovery: obese / nonobese typical
    # value, inter-individual and intercatheter random effects
    rr_md <- exp((1 - DIS_OBESE) * lrr_md + DIS_OBESE * lrr_md_obese + etalrr_md + etaiov_lrr_md)

    # Microdialysate concentration. The paper used the integrated
    # dialysate approach (RR times the mean Cisf over the collection
    # interval); here Cisf is evaluated at the observation time, so place
    # each microdialysate record at its collection-interval midpoint.
    Cmd <- rr_md * Cisf

    # Retrodialysis: C_RD / C_RP = 1 - RR
    Crd_ratio <- 1 - rr_md

    Cc ~ prop(propSd)
    Cu ~ prop(propSd_Cu)
    Cmd ~ prop(propSd_Cmd)
    Crd_ratio ~ prop(propSd_Crd_ratio)
  })
}
