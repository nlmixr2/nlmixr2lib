Xu_2019_sarilumab <- function() {
  description <- "Two-compartment population PK model for sarilumab in adults with rheumatoid arthritis (Xu 2019), with first-order SC absorption and parallel linear plus Michaelis-Menten (target-mediated) elimination from the central compartment."
  reference <- "Xu C, Su Y, Paccaly A, Kanamaluru V. Population Pharmacokinetics of Sarilumab in Patients with Rheumatoid Arthritis. Clin Pharmacokinet. 2019;58(11):1455-1467. doi:10.1007/s40262-019-00765-1"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight (time-varying)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CLO/F and on Vm, each normalized as WT/71 kg per Xu 2019 Table 3 and the final-model equations for CLO/F and Vm.",
      source_name        = "WT"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male; SEXF term evaluates to 1 and has no effect)",
      notes              = "Xu 2019 codes SEX=1 for female in the final-model equation for CLO/F. This encoding was operator-confirmed (see model extraction task 003 stop-and-ask) based on the paper's narrative that male patients have higher apparent clearance and lower AUC0-14d than female patients.",
      source_name        = "SEX"
    ),
    ADA_POS = list(
      description        = "Anti-drug antibody positivity (time-varying in final model)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ADA-negative; typical patient)",
      notes              = "Time-varying ADA indicator on CLO/F (primary covariate assessment in Xu 2019). The typical patient is ADA-negative.",
      source_name        = "ADA"
    ),
    FORM_DP2 = list(
      description        = "Sarilumab drug product 2 indicator (DP2 formulation)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (DP1 or DP3; DP3 is the commercial product)",
      notes              = "DP2 was used in some phase I studies and the dose-ranging phase II study; it is not the marketed formulation. Affects both Ka and CLO/F. Set to 0 for routine commercial-formulation simulation.",
      source_name        = "DP2"
    ),
    ALBR = list(
      description        = "Serum albumin normalized to the laboratory upper limit of normal",
      units              = "(unitless ratio)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Applied as (ALBR / 0.78)^e_albr_vm on Vm. Reference 0.78 corresponds to a median serum albumin of 38 g/L at a typical ULN of ~48.7 g/L per Xu 2019 final-model narrative.",
      source_name        = "ALBR"
    ),
    CRCL = list(
      description        = "Body-surface-area-normalized creatinine clearance (measured CrCl; CRCL = 1.73 * CrCl / BSA)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Xu 2019 defines the Vm covariate term as (1.73 * CrCl / BSA / 100)^theta13 where CrCl is in mL/min and BSA in m^2; this canonical column carries the precomputed 1.73*CrCl/BSA value with reference 100 mL/min/1.73 m^2. Mapped to the canonical general-scope CRCL covariate (which also accepts MDRD-estimated eGFR in the same units); the measured-CrCl vs estimated-eGFR distinction is documented here in the description.",
      source_name        = "1.73*CrCl/BSA"
    ),
    CRP = list(
      description        = "Baseline (pre-treatment) C-reactive protein; time-fixed per subject",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Applied as (CRP / 14.2)^e_crp_vm on Vm per Xu 2019 Table 3 and the Vm equation. Reference 14.2 mg/L is the median of the Pop-PK dataset. Source column 'BLCRP' (baseline CRP) maps to the canonical general-scope CRP covariate; the baseline-only semantics are documented here in the covariateData entry.",
      source_name        = "BLCRP"
    )
  )

  population <- list(
    n_subjects     = 1770L,
    n_observations = 7676L,
    n_studies      = 12L,
    age_range      = "18-87 years",
    age_median     = "53 years",
    weight_range   = "31.5-176.9 kg",
    weight_median  = "71.0 kg",
    sex_female_pct = 83,
    race_ethnicity = c(White = 88, Black = 3, Asian = 6, Other = 3),
    disease_state  = "Moderate-to-severe rheumatoid arthritis (adults) with inadequate response to methotrexate, TNF inhibitors, or other DMARDs.",
    dose_range     = "50-200 mg SC as single or repeated doses (qw or q2w) across 7 phase I, 1 phase II, and 4 phase III studies. The marketed regimen is 200 mg SC q2w, reducible to 150 mg SC q2w for safety management.",
    regions        = "Multi-regional (North America, EU, Latin America, Japan, and other regions represented across the 12 pooled studies).",
    notes          = "Baseline demographics from Xu 2019 Table 2 (final Pop-PK dataset of 1770 patients and 7676 concentration-time points). Concomitant methotrexate 91%; prior biologics 22% (with 26% unknown); ADA-positive 18%; sarilumab drug product DP1 4%, DP2 20%, DP3 76%. Baseline albumin median 38 g/L, CrCl median 104.8 mL/min, baseline CRP median 14.2 mg/L."
  )

  ini({
    # Structural PK parameters - Xu 2019 Table 3 final-model estimates (reference covariate
    # values: typical 71 kg female, ADA-negative, non-DP2 drug product, ALBR = 0.78, CRCL =
    # 100 mL/min/1.73 m^2, baseline CRP = 14.2 mg/L).
    lka <- log(0.136); label("Absorption rate Ka (1/day)")                                    # Xu 2019 Table 3, Ka row
    lcl <- log(0.260); label("Apparent linear clearance CLO/F (L/day)")                       # Xu 2019 Table 3, CLO/F row
    lvc <- log(2.08);  label("Apparent central volume Vc/F (L)")                              # Xu 2019 Table 3, Vc/F row
    lvp <- log(5.23);  label("Apparent peripheral volume Vp/F (L)")                           # Xu 2019 Table 3, Vp/F row
    lq  <- log(0.156); label("Apparent intercompartmental clearance Q/F (L/day)")             # Xu 2019 Table 3, Q/F row

    # Parallel nonlinear Michaelis-Menten elimination from the central compartment.
    lvm <- log(8.06);  label("Maximum Michaelis-Menten elimination rate Vm (mg/day)")         # Xu 2019 Table 3, Vm row
    lkm <- log(0.939); label("Michaelis-Menten constant Km (mg/L)")                           # Xu 2019 Table 3, Km row

    # Covariate exponents and multiplicative effects - Xu 2019 Table 3 and the final-model
    # equations for Vm, CLO/F, and Ka.
    e_wt_cl    <-  0.885;  label("Power exponent of WT/71 on CLO/F (unitless)")                 # Xu 2019 Table 3: WT effect on CLO/F
    e_wt_vm    <-  0.516;  label("Power exponent of WT/71 on Vm (unitless)")                    # Xu 2019 Table 3: WT effect on Vm
    e_albr_vm  <- -0.844;  label("Power exponent of ALBR/0.78 on Vm (unitless)")                # Xu 2019 Table 3: ALBR effect on Vm
    e_crcl_vm  <-  0.212;  label("Power exponent of CRCL/100 on Vm (unitless)")                 # Xu 2019 Table 3: CrCl effect on Vm
    e_crp_vm   <-  0.0299; label("Power exponent of CRP/14.2 on Vm (unitless)")                 # Xu 2019 Table 3: baseline CRP effect on Vm
    e_dp2_ka   <-  0.663;  label("Multiplier on Ka for drug product DP2 (unitless)")            # Xu 2019 Table 3: DP2 effect on Ka
    e_ada_cl   <-  1.43;   label("Multiplier on CLO/F for ADA-positive (unitless)")             # Xu 2019 Table 3: ADA effect on CLO/F
    e_dp2_cl   <-  1.30;   label("Multiplier on CLO/F for drug product DP2 (unitless)")         # Xu 2019 Table 3: DP2 effect on CLO/F
    e_sexf_cl  <-  0.846;  label("Multiplier on CLO/F for female sex (unitless)")               # Xu 2019 Table 3: SEX effect on CLO/F; SEX=1=female per operator confirmation

    # Inter-individual variability: Xu 2019 Table 3 reports CV% on the linear-parameter scale.
    # Convert to NONMEM-style log-normal variance as omega^2 = log(CV^2 + 1):
    #   Vm   CV 32.4% -> omega^2 = log(0.324^2 + 1) = 0.0998
    #   CLO/F CV 55.3% -> omega^2 = log(0.553^2 + 1) = 0.2669
    #   Vc/F CV 37.3% -> omega^2 = log(0.373^2 + 1) = 0.1302
    #   Ka   CV 32.1% -> omega^2 = log(0.321^2 + 1) = 0.0981
    # Vm-CLO/F correlation -0.566 -> cov = -0.566 * sqrt(0.0998 * 0.2669) = -0.0924
    etalvm + etalcl ~ c(0.0998, -0.0924, 0.2669)  # Xu 2019 Table 3: Vm IIV 32.4% CV, CLO/F IIV 55.3% CV, Vm-CLO/F correlation -0.566
    etalvc ~ 0.1302                                # Xu 2019 Table 3: Vc/F IIV 37.3% CV (shrinkage 64.2%)
    etalka ~ 0.0981                                # Xu 2019 Table 3: Ka IIV 32.1% CV (shrinkage 49.6%)

    # Residual error: Xu 2019 fitted log-transformed concentrations with an additive residual
    # error on the log scale (NONMEM log-EPS), reported variance sigma^2 = 0.395. This maps to
    # a proportional error in linear space with SD = sqrt(0.395) = 0.6285.
    propSd <- 0.6285; label("Proportional residual error (fraction)")   # Xu 2019 Table 3: residual sigma^2 = 0.395 (log-additive)
  })
  model({
    # Individual PK parameters with Xu 2019 covariate models.
    vm <- exp(lvm + etalvm) *
          (WT / 71)^e_wt_vm *
          (ALBR / 0.78)^e_albr_vm *
          (CRCL / 100)^e_crcl_vm *
          (CRP / 14.2)^e_crp_vm
    km <- exp(lkm)
    vc <- exp(lvc + etalvc)
    cl <- exp(lcl + etalcl) *
          (WT / 71)^e_wt_cl *
          e_ada_cl^ADA_POS *
          e_dp2_cl^FORM_DP2 *
          e_sexf_cl^SEXF
    ka <- exp(lka + etalka) * e_dp2_ka^FORM_DP2
    q  <- exp(lq)
    vp <- exp(lvp)

    # Two-compartment PK with parallel linear and Michaelis-Menten elimination from central.
    Cc <- central / vc

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot -
                          (cl / vc) * central -
                          vm * Cc / (km + Cc) -
                          (q / vc) * central +
                          (q / vp) * peripheral1
    d/dt(peripheral1) <-  (q / vc) * central - (q / vp) * peripheral1

    Cc ~ prop(propSd)
  })
}
