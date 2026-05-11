Xu_2023_MBG453 <- function() {
  description <- "Two-compartment population PK model for sabatolimab (MBG453, anti-TIM-3 IgG4) with parallel linear and Michaelis-Menten elimination from the central compartment, fit to pooled adult patients with advanced solid tumors and hematologic malignancies (Xu 2023)."
  reference <- "Xu S, Zhang N, Rinne ML, Sun H, Stein AM. Sabatolimab (MBG453) model-informed drug development for dose selection in patients with myelodysplastic syndrome/acute myeloid leukemia and solid tumors. CPT Pharmacometrics Syst Pharmacol. 2023;12(11):1653-1665. doi:10.1002/psp4.12962"
  vignette <- "Xu_2023_MBG453"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form scaling on CL, Vc, and Vp via the paper's log-deviation form exp(beta * log(WT / WT_median)) (Xu 2023, p1657 covariate equations). The paper does not state the numerical median baseline weight; this model uses a working reference of 75 kg. Users with an alternative reference can override e_wt_*.",
      source_name        = "WT0"
    ),
    DIS_AML = list(
      description        = "Acute myeloid leukemia disease indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-AML; the reference complement in Xu 2023 is the advanced-solid-tumor cohort, alongside DIS_MDS = DIS_CMML = 0)",
      notes              = "Exponential effect on CL (Xu 2023 Table 1 row beta_CL,AML = -0.0146, NS p = 0.801). Source column is the categorical DISEASE_abb in the Monolix supplement Appendix S2; the canonical column is the binary as.integer(DISEASE_abb == 'AML').",
      source_name        = "DISEASE_abb"
    ),
    DIS_MDS = list(
      description        = "Myelodysplastic syndrome disease indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-MDS; the reference complement in Xu 2023 is the advanced-solid-tumor cohort, alongside DIS_AML = DIS_CMML = 0)",
      notes              = "Exponential effect on CL (Xu 2023 Table 1 row beta_CL,MDS = -0.149, statistically significant p = 0.0213; ~14% lower CL than the solid-tumor reference). Source column is the categorical DISEASE_abb; the canonical column is the binary as.integer(DISEASE_abb == 'MDS').",
      source_name        = "DISEASE_abb"
    ),
    DIS_CMML = list(
      description        = "Chronic myelomonocytic leukemia disease indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-CMML; the reference complement in Xu 2023 is the advanced-solid-tumor cohort, alongside DIS_AML = DIS_MDS = 0)",
      notes              = "Exponential effect on CL (Xu 2023 Table 1 row beta_CL,CMML = -0.0411, NS p = 0.76). Source column is the categorical DISEASE_abb; the canonical column is the binary as.integer(DISEASE_abb == 'CMML').",
      source_name        = "DISEASE_abb"
    ),
    CONMED_SPART = list(
      description        = "Spartalizumab (PDR001, anti-PD-1) coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no spartalizumab coadministration; the reference complement in Xu 2023 is the union of sabatolimab monotherapy and sabatolimab + hypomethylating-agent combination)",
      notes              = "Exponential effect on CL (Xu 2023 Table 1 row beta_CL,HASPDR = 0.0194, NS p = 0.7). Renamed from the source column HASPDR (the supplement Appendix S2 documents the column as 'this patient HAS received PDR001 [spartalizumab, anti PD-1 mAb]') to the canonical CONMED_SPART per covariate-columns.md.",
      source_name        = "HASPDR"
    )
  )

  population <- list(
    n_subjects     = 444L,
    n_studies      = 2L,
    age_range      = "adults (>= 18 years); median age in the hematologic-malignancy cohort 71-72.5 years across treatment arms (Table S1)",
    age_median     = "not reported pooled across solid tumor + hematologic; ~72 years in the hematologic subgroup",
    weight_range   = "not stated explicitly in the source paper",
    weight_median  = "not stated in the source paper; this model assumes a working reference of 75 kg for the WT-on-CL/V/V2 scaling",
    sex_female_pct = "~40% in the hematologic-malignancy cohort (Table S1: 38.5-42.1% female across arms); not reported pooled with solid tumors",
    race_ethnicity = "predominantly Caucasian in the hematologic-malignancy cohort (Table S1: 68.9-87.9% Caucasian, 0-6.7% Asian, 0-7.7% Black, 1.4-2.2% Other, 7.9-17.8% Unknown); solid-tumor demographics previously reported separately in NCT02608268",
    disease_state  = "Advanced/metastatic solid tumors (n=252) or hematologic malignancies (n=192): acute myeloid leukemia (AML), myelodysplastic syndrome (MDS, including intermediate-, high-, and very high-risk per Revised IPSS), or chronic myelomonocytic leukemia (CMML) ineligible for intensive chemotherapy",
    dose_range     = "Sabatolimab IV 20-1200 mg Q2W or 80-1200 mg Q4W as monotherapy or in combination with spartalizumab (anti-PD-1) and/or hypomethylating agents (decitabine or azacitidine). Specific arms: solid tumor 80-1200 mg Q2W/Q4W single agent, 20-800 mg Q2W or 80-1200 mg Q4W with spartalizumab; heme 400 or 1200 mg Q2W single agent, 240/400/800 mg with HMA, 160/240/400 mg Q2W with spartalizumab +/- decitabine.",
    regions        = "Two pooled studies: NCT02608268 (advanced solid tumors, phase I-Ib/II) and NCT03066648 (hematologic malignancies, phase Ib). Multi-regional but specific regions are not enumerated in the paper.",
    ecog_distribution = "ECOG performance status 0-2 (eligibility cap). In the hematologic-malignancy cohort, 15.4-32.1% PS 0, 58.6-69.2% PS 1, 9.3-15.6% PS 2 across treatment arms (Table S1).",
    notes          = "Pooled PK dataset of 444 patients (252 solid tumors + 192 hematologic malignancies) supports the final model in Xu 2023 Table 1. PK samples were quantified by LC-MS with LLOQ 1 ug/mL; total soluble TIM-3 was quantified by ELISA with LLOQ 1.02 ng/mL (not modeled here -- the sTIM-3 model is described in Xu 2023 as a separately-fit QSS-TMDD model from the prior phase I solid-tumor analysis and is not republished). Treatment-arm composition: 159 sabatolimab monotherapy (133 solid tumor + 26 heme), 130 sabatolimab + spartalizumab (119 solid tumor + 11 heme), 55 sabatolimab + azacitidine (heme), 100 sabatolimab + decitabine +/- spartalizumab (heme)."
  )

  ini({
    # Structural PK parameters - Xu 2023 Table 1 final-model estimates (p1660), reference patient is
    # the typical solid-tumor patient (DIS_AML = DIS_MDS = DIS_CMML = 0, CONMED_SPART = 0) at
    # WT = WT_ref = 75 kg. CL, Q, and Vm were reported in per-hour units in Table 1 and are
    # converted to per-day here (multiplied by 24) because this model keeps time in days.
    lcl  <- log(0.0103 * 24); label("Linear clearance CL (L/day)")                                  # Xu 2023 Table 1: CL = 0.0103 L/h
    lvc  <- log(3.59);        label("Central volume of distribution V (L)")                          # Xu 2023 Table 1: V  = 3.59 L
    lq   <- log(0.0353 * 24); label("Inter-compartmental clearance Q (L/day)")                       # Xu 2023 Table 1: Q  = 0.0353 L/h
    lvp  <- log(2.38);        label("Peripheral volume of distribution V2 (L)")                      # Xu 2023 Table 1: V2 = 2.38 L

    # Parallel Michaelis-Menten elimination from the central compartment (Xu 2023 Methods p1657
    # ODEs and Table 1). Vm reported in (ug/mL)/h; converted to (ug/mL)/day. Km = 0.5 nM was
    # FIXED at 0.074 ug/mL (sabatolimab MW ~150 kDa) per the paper's footnote a, because Km was
    # not estimable from the PK data alone.
    lvmax  <- log(0.0197 * 24);    label("Michaelis-Menten maximum-velocity Vm ((ug/mL)/day)")          # Xu 2023 Table 1: Vm = 0.0197 ug/mL/h
    lkm  <- fixed(log(0.074));   label("Michaelis-Menten constant Km (ug/mL, fixed at 0.5 nM)")       # Xu 2023 Table 1 footnote a (FIXED)

    # Covariate effects (Xu 2023 Table 1 and the covariate equations on p1657).
    # Allometric exponents on baseline weight: applied as (WT/WT_ref)^e_wt_*.
    e_wt_cl  <- 0.743;  label("Power exponent of WT/WT_ref on CL (unitless)")                         # Xu 2023 Table 1: beta_CL,WT0
    e_wt_vc  <- 0.77;   label("Power exponent of WT/WT_ref on Vc (unitless)")                         # Xu 2023 Table 1: beta_V,WT0
    e_wt_vp  <- 0.597;  label("Power exponent of WT/WT_ref on Vp (unitless)")                         # Xu 2023 Table 1: beta_V2,WT0

    # Disease and concomitant-medication effects on CL (exponential form per the paper's
    # equation: CL_i = CLpop * exp(beta_CL,X * X_i)). The reference category is the
    # solid-tumor cohort with no spartalizumab coadministration.
    e_dis_aml_cl       <- -0.0146; label("Exponential coefficient of DIS_AML on CL (unitless; NS p=0.801)") # Xu 2023 Table 1: beta_CL,AML
    e_dis_mds_cl       <- -0.149;  label("Exponential coefficient of DIS_MDS on CL (unitless; p=0.0213)")   # Xu 2023 Table 1: beta_CL,MDS
    e_dis_cmml_cl      <- -0.0411; label("Exponential coefficient of DIS_CMML on CL (unitless; NS p=0.76)") # Xu 2023 Table 1: beta_CL,CMML
    e_coadmin_spart_cl <-  0.0194; label("Exponential coefficient of CONMED_SPART on CL (unitless; NS p=0.7)") # Xu 2023 Table 1: beta_CL,HASPDR

    # Inter-individual variability - Xu 2023 Table 1 reports the SD (omega) of each log-normal
    # eta and the correlation between eta_V and eta_CL. Variance = omega^2; covariance =
    # corr(V,CL) * omega_V * omega_CL = 0.634 * 0.230 * 0.473 = 0.06897. No IIV on Q or Km
    # (paper Methods p1657: random effects on Q and Km were removed due to non-convergence).
    etalcl + etalvc ~ c(0.22373,
                        0.06897, 0.05290)                                                            # Xu 2023 Table 1: omega_CL = 0.473, omega_V = 0.230, corr_V_CL = 0.634
    etalvp ~ 0.11424                                                                                  # Xu 2023 Table 1: omega_V2 = 0.338
    etalvmax ~ 0.41088                                                                                    # Xu 2023 Table 1: omega_Vm = 0.641

    # Residual error - Xu 2023 Methods p1657 combined-error model:
    # LIDV = C + sqrt(a^2 + (b*C)^2) * e, with e ~ N(0, 1). a is the additive SD in ug/mL and
    # b is the proportional SD (fraction). Equivalent to nlmixr2's add(addSd) + prop(propSd).
    addSd  <- 1.41; label("Additive residual error (ug/mL)")                                          # Xu 2023 Table 1: a = 1.41
    propSd <- 0.19; label("Proportional residual error (fraction)")                                   # Xu 2023 Table 1: b = 0.19
  })
  model({
    # Reference baseline weight for the (WT/WT_ref) scaling. The paper does not state the
    # numerical median baseline weight; 75 kg is used here as a working reference and is
    # documented in the model's covariate notes. The exponent estimates are anchored to the
    # actual (unreported) median, so substituting a different reference rescales the typical-
    # value PK parameters but leaves the WT effect-shape intact.
    wt_ref <- 75

    # Individual baseline PK parameters with the Xu 2023 covariate equations (p1657).
    # Reference patient: WT = wt_ref, solid tumor (DIS_AML = DIS_MDS = DIS_CMML = 0),
    # no spartalizumab (CONMED_SPART = 0).
    cl <- exp(lcl + etalcl) *
          (WT / wt_ref)^e_wt_cl *
          exp(e_dis_aml_cl       * DIS_AML +
              e_dis_mds_cl       * DIS_MDS +
              e_dis_cmml_cl      * DIS_CMML +
              e_coadmin_spart_cl * CONMED_SPART)
    vc <- exp(lvc + etalvc) * (WT / wt_ref)^e_wt_vc
    vp <- exp(lvp + etalvp) * (WT / wt_ref)^e_wt_vp
    q    <- exp(lq)
    vmax <- exp(lvmax + etalvmax)
    km   <- exp(lkm)

    # Two-compartment PK with parallel first-order linear and Michaelis-Menten elimination
    # from the central compartment. Concentration in ug/mL (= mg/L); central in mg.
    Cc <- central / vc

    d/dt(central)     <- -(cl / vc) * central -
                          vmax * vc * Cc / (km + Cc) -
                          (q / vc) * central +
                          (q / vp) * peripheral1
    d/dt(peripheral1) <-  (q / vc) * central -
                          (q / vp) * peripheral1

    Cc ~ add(addSd) + prop(propSd)
  })
}
