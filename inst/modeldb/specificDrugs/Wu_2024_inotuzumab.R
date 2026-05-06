Wu_2024_inotuzumab <- function() {
  description <- "Two-compartment population PK model for inotuzumab ozogamicin in pediatric and adult patients with relapsed/refractory B-cell precursor acute lymphoblastic leukemia (BCP-ALL) and adult patients with B-cell non-Hodgkin's lymphoma (NHL); linear plus time-dependent (target-mediated) clearance with covariate effects on CL_SS, Vc, CL_TIME, and kdes (Wu 2024, ITCC-059 pediatric trial pooled with 11 adult studies)."
  reference <- "Wu JH, Pennesi E, Bautista F, Garrett M, Fukuhara K, Brivio E, et al. Population Pharmacokinetics of Inotuzumab Ozogamicin in Pediatric Relapsed/Refractory B-Cell Precursor Acute Lymphoblastic Leukemia: Results of Study ITCC-059. Clin Pharmacokinet. 2024;63(7):981-997. doi:10.1007/s40262-024-01386-z"
  vignette <- "Wu_2024_inotuzumab"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    LBM = list(
      description        = "Baseline lean body mass",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline. Power effect on CL_SS (exponent 1.05), Vc (0.977), and CL_TIME (0.687) with reference 52.7 kg (population median; Wu 2024 Table 3 footnote). Estimated by Boer's equations for adults and the Peters et al. equation for children (Wu 2024 Table 2 footnote b).",
      source_name        = "LBM"
    ),
    AGE = list(
      description        = "Baseline age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline. Power effect on kdes (exponent -0.296) for BCP-ALL patients only, with reference 60 years (per Wu 2024 Table 3). For NHL patients the AGE effect is gated off via DIS_BCPALL.",
      source_name        = "AGE"
    ),
    DIS_BCPALL = list(
      description        = "B-cell precursor acute lymphoblastic leukemia disease-state indicator (1 = BCP-ALL, 0 = B-cell non-Hodgkin's lymphoma)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (B-cell NHL)",
      notes              = "Wu 2024 calls this the 'ALL effect' and notes it accounts for both disease type (NHL vs ALL) and the corresponding bioanalytical analysis method (ELISA for adult NHL vs HPLC-MS for ALL). Fractional-change (dummy-variable) effects on CL_SS (-0.767) and CL_TIME (-0.362), and gates the BLSTABL and AGE effects on kdes (kdes itself also takes a -0.924 fractional change for BCP-ALL). Source column 'ALL'; renamed to canonical DIS_BCPALL per inst/references/covariate-columns.md.",
      source_name        = "ALL"
    ),
    CONMED_RITUX = list(
      description        = "Concomitant rituximab combination-therapy indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant rituximab)",
      notes              = "Wu 2024 Table 3 footnote b: in the ITCC-059 reanalysis the reference category is 'without rituximab' (RITUX = 0); the original adult Garrett 2019 model used 'with rituximab' as reference. Fractional change on CL_SS: -0.132 (i.e., CL_SS is ~13% lower with concomitant rituximab). Source column 'RITUX'.",
      source_name        = "RITUX"
    ),
    BLSTABL = list(
      description        = "Baseline absolute blast counts in peripheral blood",
      units              = "10^9 counts/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline. Power effect on kdes (exponent -0.0484) for BCP-ALL patients only, with reference 0.352 x 10^9 counts (population median for BCP-ALL; Wu 2024 Table 3). Not applicable to NHL patients (the BLSTABL effect is gated off via DIS_BCPALL); when simulating an NHL patient supply BLSTABL = 0.352 (the reference) so the gated power term evaluates to 1.",
      source_name        = "BLSTABL"
    )
  )

  population <- list(
    n_subjects     = 818L,
    n_studies      = 12L,
    cohorts        = list(
      adult_NHL          = list(n = 531L, observations = 5609L, age_median = "65 years (range 18-92)", LBM_median = "53.6 kg (range 24.4-87.4)"),
      adult_BCP_ALL      = list(n = 234L, observations = 2752L, age_median = "46 years (range 20-79)", LBM_median = "54.4 kg (range 29.4-95.9)"),
      pediatric_BCP_ALL  = list(n = 53L,  observations = 563L,  age_median = "9 years (range 1-17)",   LBM_median = "29.5 kg (range 10.9-69.0)")
    ),
    age_range      = "1-92 years overall (adult NHL 18-92; adult BCP-ALL 20-79; pediatric BCP-ALL 1-17)",
    age_median     = "Adult NHL 65 y; adult BCP-ALL 46 y; pediatric BCP-ALL 9 y (Wu 2024 Table 2)",
    weight_range   = "12.7-154 kg (Wu 2024 Table 2)",
    weight_median  = "73.2 kg (adult NHL); 74.0 kg (adult BCP-ALL); 34.3 kg (pediatric BCP-ALL)",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not reported in detail in Wu 2024 main text.",
    disease_state  = "Relapsed/refractory CD22+ B-cell precursor acute lymphoblastic leukemia (BCP-ALL) and relapsed/refractory B-cell non-Hodgkin's lymphoma (NHL).",
    dose_range     = "Adult approved regimen: 1.8 mg/m^2/cycle (cycle 1) then 1.5 mg/m^2/cycle (subsequent cycles) IV, fractionated weekly on days 1, 8, 15. Pediatric RP2D matches the adult regimen (1.8 mg/m^2/cycle in cycle 1: 0.8 + 0.5 + 0.5 mg/m^2; 1.5 mg/m^2/cycle thereafter: 0.5 + 0.5 + 0.5 mg/m^2). Phase IA pediatric dose level 1 was 1.4 mg/m^2/cycle (0.6 + 0.4 + 0.4) then 1.2 mg/m^2/cycle (0.4 + 0.4 + 0.4).",
    regions        = "Multi-regional (US adult studies via Pfizer; ITCC-059 pediatric study spanned European centers and additional sites internationally).",
    n_observations = 8924L,
    bioanalytic_methods = "Adult NHL: validated ELISA (direct measurement of N-acetyl-gamma-calicheamicin dimethyl hydrazide linked to the InO antibody). Adult BCP-ALL and pediatric BCP-ALL: HPLC-MS/MS for the conjugated calicheamicin released from the ADC, with InO quantitation based on the average drug-to-antibody ratio of the dosing standard. Three separate residual error magnitudes were estimated to capture the assay differences (Wu 2024 Table 3, footnote d).",
    notes          = "Final population PK model from the Wu 2024 (PMID 38907948) reanalysis of the Garrett 2019 adult model, refit on pooled adult + pediatric data and extended with ALL-effect on CL_TIME and AGE-effect on kdes. Used to evaluate exposure at the pediatric RP2D and to support 'no further dose adjustment required' for pediatric BCP-ALL. NONMEM 7.5.0 (SAEM + IMP); Pearl-speaks-NONMEM 5.3.0; SIR for parameter precision."
  )

  ini({
    # Structural parameters (typical population values, Wu 2024 Table 3).
    # Reference covariate values: LBM = 52.7 kg, AGE = 60 years, BLSTABL = 0.352 x 10^9, DIS_BCPALL = 0 (NHL), CONMED_RITUX = 0.
    # Time unit is hour; clearances in L/h; volumes in L; kdes in 1/h.
    # Wu 2024 names CL1 the "linear clearance" and CL2 the "initial value of
    # time-dependent clearance"; total clearance is CL_total = CL1 + CL2 *
    # exp(-kdes * time). Mapping to nlmixr2lib conventions: CL1 -> CL_SS
    # (the steady-state, non-decaying arm) and CL2 -> CL_TIME (the
    # time-varying decay arm).
    lcl <- log(0.130);  label("Linear (steady-state) clearance for an NHL adult (CL_SS, L/h)")               # Wu 2024 Table 3 (CL1)
    lvc    <- log(6.49);   label("Central volume of distribution for an NHL adult (Vc, L)")                     # Wu 2024 Table 3 (V1)
    lcl_time <- log(0.569); label("Initial value of time-dependent clearance for an NHL adult (CL_TIME, L/h)")  # Wu 2024 Table 3 (CL2)
    lkdes  <- log(0.0577); label("Decay coefficient of time-dependent clearance for an NHL adult (kdes, 1/h)")  # Wu 2024 Table 3
    lq     <- log(0.0437); label("Intercompartmental clearance (Q, L/h)")                                       # Wu 2024 Table 3
    lvp    <- log(4.74);   label("Peripheral volume of distribution (Vp, L)")                                   # Wu 2024 Table 3 (V2)

    # Covariate-effect parameters (Wu 2024 Table 3). Continuous covariates enter
    # as power models centered on the reference value; categorical effects enter
    # as additive fractional-change dummy variables of the form (1 + theta * I).
    e_lbm_cl   <-  1.05;    label("Power exponent of LBM on CL_SS (unitless)")                            # Wu 2024 Table 3
    e_lbm_vc      <-  0.977;   label("Power exponent of LBM on Vc (unitless)")                               # Wu 2024 Table 3
    e_lbm_cl_time <-  0.687;   label("Power exponent of LBM on CL_TIME (unitless)")                          # Wu 2024 Table 3
    e_all_cl   <- -0.767;   label("Fractional change in CL_SS for BCP-ALL (vs NHL, unitless)")            # Wu 2024 Table 3
    e_all_cl_time <- -0.362;   label("Fractional change in CL_TIME for BCP-ALL (vs NHL, unitless)")          # Wu 2024 Table 3
    e_all_kdes    <- -0.924;   label("Fractional change in kdes for BCP-ALL (vs NHL, unitless)")             # Wu 2024 Table 3
    e_blstabl_kdes <- -0.0484; label("Power exponent of BLSTABL on kdes for BCP-ALL only (unitless)")        # Wu 2024 Table 3
    e_age_kdes    <- -0.296;   label("Power exponent of AGE on kdes for BCP-ALL only (unitless)")            # Wu 2024 Table 3
    e_ritux_cl <- -0.132;   label("Fractional change in CL_SS for concomitant rituximab (unitless)")      # Wu 2024 Table 3

    # Inter-individual variability. Wu 2024 Table 3 reports CV% for IIV with
    # the convention CV = sqrt(omega^2) (i.e., the OMEGA variance equals the
    # squared CV); this is verified by reproducing the published correlations
    # from the off-diagonal covariances (e.g., 0.136 / sqrt(0.16 * 0.1608)
    # = 0.847 = 84.7% as in Table 3 'CL1 - V1; correlations'). CL_SS, Vc,
    # CL_TIME are reported as a 3x3 correlated block; kdes is independent.
    etalcl + etalvc + etalcl_time ~ c(0.16,
                                         0.136, 0.16080,
                                         0.194, 0.204, 0.54317)  # Wu 2024 Table 3 (CV%: CL_SS 40.0, Vc 40.1, CL_TIME 73.7; covariances 0.136 / 0.194 / 0.204)
    etalkdes ~ 0.35641  # Wu 2024 Table 3 (CV% kdes 59.7)

    # Residual error. Wu 2024 reports the residual SD on log-transformed data,
    # which maps to a log-normal residual error in nlmixr2 (Cc ~ lnorm(SD)) and
    # to proportional error in the small-SD limit. Three population-stratified
    # SDs are estimated (Table 3, footnote d): adult NHL 0.444, adult BCP-ALL
    # 0.612, pediatric BCP-ALL 0.452. The packaged ini() uses the pediatric
    # BCP-ALL value (the focus of Wu 2024); to simulate adult NHL or adult
    # BCP-ALL, override `expSd` via `ini(model, expSd = 0.444)` or
    # `ini(model, expSd = 0.612)` respectively (see vignette).
    expSd <- 0.452; label("Log-scale residual SD for pediatric BCP-ALL population (unitless)")  # Wu 2024 Table 3 (pediatric ALL proportional residual error 0.452; adult NHL 0.444; adult BCP-ALL 0.612)
  })

  model({
    # Individual PK parameters with covariate effects (Wu 2024 Table 3 final
    # model equations, restated; CL1 -> CL_SS, CL2 -> CL_TIME):
    #   CL_SS   = 0.130 * (1 - 0.767 * ALL) * (LBM/52.7)^1.05  * (1 - 0.132 * RITUX)
    #   Vc      = 6.49  * (LBM/52.7)^0.977
    #   CL_TIME = 0.569 * (1 - 0.362 * ALL) * (LBM/52.7)^0.687
    #   kdes_NHL = 0.0577
    #   kdes_ALL = 0.0577 * (1 - 0.924) * (BLSTABL/0.352)^-0.0484 * (AGE/60)^-0.296
    #   Q     = 0.0437
    #   Vp    = 4.74
    cl <- exp(lcl + etalcl) *
      (1 + e_all_cl   * DIS_BCPALL) *
      (LBM / 52.7)^e_lbm_cl *
      (1 + e_ritux_cl * CONMED_RITUX)

    vc <- exp(lvc + etalvc) * (LBM / 52.7)^e_lbm_vc

    cl_time <- exp(lcl_time + etalcl_time) *
      (1 + e_all_cl_time * DIS_BCPALL) *
      (LBM / 52.7)^e_lbm_cl_time

    # kdes covariate effects (BLSTABL, AGE) only apply for BCP-ALL patients in
    # Wu 2024; gating is encoded by multiplying each exponent by DIS_BCPALL so
    # that for NHL the corresponding power term collapses to 1. For NHL
    # patients supply BLSTABL at the BCP-ALL reference (0.352 x 10^9) so the
    # gated power term evaluates to 1 numerically.
    kdes <- exp(lkdes + etalkdes) *
      (1 + e_all_kdes * DIS_BCPALL) *
      (BLSTABL / 0.352)^(e_blstabl_kdes * DIS_BCPALL) *
      (AGE / 60)^(e_age_kdes * DIS_BCPALL)

    q  <- exp(lq)
    vp <- exp(lvp)

    # Time-dependent total clearance (Wu 2024 Methods, Section 2.3):
    #   CL_t    = CL_TIME * exp(-kdes * time)
    #   CL_total = CL_SS + CL_t
    # `time` is the internal integration time in hours, which corresponds to
    # time from the first dose for the event datasets this model expects.
    cl_t   <- cl_time * exp(-kdes * time)
    cl_tot <- cl + cl_t

    # Two-compartment model with IV input (no depot; doses are administered
    # directly into the central compartment as 60-minute IV infusions, see
    # Wu 2024 Methods Section 2.1).
    kel <- cl_tot / vc
    k12 <- q     / vc
    k21 <- q     / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Concentration: dose in mg, vc in L -> central/vc in mg/L = ug/mL.
    # Multiply by 1000 to report Cc in ng/mL (matching the bioanalytical LLOQ
    # of 1.0 ng/mL reported by Wu 2024).
    Cc <- (central / vc) * 1000

    # Log-normal residual error (Wu 2024 reports the residual as 'additive on
    # log-transformed data', i.e., NONMEM Y = LOG(IPRED) + EPS(1) with
    # EPS(1) ~ N(0, sigma^2)). expSd defaults to the pediatric BCP-ALL value;
    # see ini() and the vignette for adult-population alternatives.
    Cc ~ lnorm(expSd)
  })
}
