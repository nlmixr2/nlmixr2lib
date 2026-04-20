Ma_2020_sarilumab_das28crp <- function() {
  description <- "Indirect-response PK/PD model of sarilumab on the 28-joint disease activity score by C-reactive protein (DAS28-CRP) in adults with rheumatoid arthritis (Ma 2020). Sarilumab inhibits the DAS28-CRP production rate (kin) via a sigmoid Emax function that includes a background DMARD placebo component (PLB). The PK driver is the two-compartment, parallel linear + Michaelis-Menten model of Xu 2019 evaluated at its typical covariate-reference values (adult female, 71 kg, ADA-negative, commercial drug product, ALBR = 0.78, CrCl = 100 mL/min/1.73 m^2, baseline CRP = 14.2 mg/L)."
  reference <- "Ma L, Xu C, Paccaly A, Kanamaluru V. Population Pharmacokinetic-Pharmacodynamic Relationships of Sarilumab Using Disease Activity Score 28-Joint C-Reactive Protein and Absolute Neutrophil Counts in Patients with Rheumatoid Arthritis. Clin Pharmacokinet. 2020;59(11):1451-1466. doi:10.1007/s40262-020-00899-7. PMID: 32451909. PK backbone from Xu C, Su Y, Paccaly A, Kanamaluru V. Population Pharmacokinetics of Sarilumab in Patients with Rheumatoid Arthritis. Clin Pharmacokinet. 2019;58(11):1455-1467. doi:10.1007/s40262-019-00765-1."
  vignette <- "Ma_2020_sarilumab_das28crp"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L", response = "DAS28-CRP score (unitless, 0-10 scale)")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on DAS28-CRP BASE normalized as WT/72.8 (Ma 2020 Table 3 reference weight = median of DAS28-CRP final dataset per paper narrative). The 71 kg reference used for the embedded Xu 2019 PK typical profile is internal to the model and is not exposed through this covariate.",
      source_name        = "WT"
    ),
    CRP = list(
      description        = "Baseline (pre-treatment) C-reactive protein measured by the routine clinical assay; time-fixed per subject",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on BASE and additive log-linear effect on the logit-transformed Emax (Ma 2020 Table 3); reference 15.7 mg/L is the median baseline CRP of the DAS28-CRP dataset per paper narrative. Source column 'CRP' (baseline CRP, standard assay) maps to the canonical general-scope CRP covariate; the baseline-only and standard-assay semantics are documented here in the covariateData entry rather than via a separate CRP canonical.",
      source_name        = "CRP"
    ),
    BLPHYVAS = list(
      description        = "Baseline Physician's Global Assessment of Disease Activity (100-mm VAS)",
      units              = "mm (0-100 VAS)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Power effect on BASE (Ma 2020 Table 3); reference 66 is the median baseline PHYVAS of the DAS28-CRP dataset per paper narrative.",
      source_name        = "BLPHYVAS"
    ),
    BLHAQ = list(
      description        = "Baseline Health Assessment Questionnaire Disability Index (0-3 score)",
      units              = "unitless (0-3 composite)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Power effect on BASE (Ma 2020 Table 3); reference 1.75 is the median baseline HAQ-DI of the DAS28-CRP dataset per paper narrative.",
      source_name        = "BLHAQ"
    ),
    PRICORT = list(
      description        = "Prior corticosteroid treatment indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no prior corticosteroid use)",
      notes              = "Time-fixed per subject. Multiplicative effect on Kout (Ma 2020 Table 3): Kout_i = Kout * 1.26^PRICORT. Paper narrative confirms 0.0333 vs 0.0264 day^-1 for PRICORT = 1 vs 0.",
      source_name        = "PRICORT"
    )
  )

  population <- list(
    n_subjects        = 2082L,
    n_observations    = 17229L,
    n_studies         = 3L,
    age_range         = "adults; mean (SD) 51.6 (12.0) years",
    weight_range      = "mean (SD) 74.6 (18.8) kg; median 72.8 kg (per narrative)",
    sex_female_pct    = 82.1,
    race_ethnicity    = c(Caucasian = 83.6, Other = 16.4),
    disease_state     = "Moderate-to-severely-active rheumatoid arthritis in adults with inadequate response to methotrexate (MTX-IR, 26.1%) or to TNFalpha inhibitors (TNF-IR, 73.9%); all patients received background DMARD therapy (MTX 98.8%).",
    dose_range        = "Sarilumab 100, 150, or 200 mg SC q2w, and 100 or 150 mg SC qw; placebo arm also modelled. Treatment durations 12, 24, and 52 weeks across studies.",
    regions           = "Multi-regional (North America, EU, Latin America, and other regions represented in MOBILITY and TARGET phase II-III programs).",
    baseline_biomarkers = list(
      CRP_mean_sd_mg_L = "24.1 (25.1)",
      CRP_median_mg_L  = 15.7,
      BLIL6_mean_sd_pg_mL = "41.8 (67.2)",
      BLPHYVAS_mean_sd   = "64.6 (16.8)",
      BLPHYVAS_median    = 66,
      BLHAQ_mean_sd      = "1.68 (0.640)",
      BLHAQ_median       = 1.75
    ),
    concomitant_treatment = list(
      methotrexate_pct     = 98.8,
      prior_biologic_pct   = 39.6,
      prior_corticosteroid_pct = 64.6,
      baseline_ACCP_pos_pct = 16.5
    ),
    notes             = "Baseline demographics from Ma 2020 Table 2 (DAS28-CRP final dataset, n=2082 across NCT01061736 Part A/B [MOBILITY phase II/III] and NCT01709578 [TARGET phase III]). Baseline mean DAS28-CRP was not reported in Table 2 for the DAS28-CRP dataset itself (the ANC dataset reported 6.03, consistent with the modelled BASE of 6.06 in Table 3). Sources pooled 17,229 DAS28-CRP observations through week 24."
  )

  ini({
    # --------------------------------------------------------------------------
    # PK backbone (Xu 2019 Table 3) - fixed at typical reference-covariate values.
    # The current model uses these as structural constants so a single file
    # reproduces the full sarilumab PK/PD cascade; individual PK covariates
    # (ADA, drug product, sex, ALBR, CrCl, WT-on-PK, CRP-on-Vm) are omitted
    # here because the Ma 2020 PopPK/PD analysis used sequential individual-PK
    # Bayes estimates as the exposure input. See the vignette's Assumptions and
    # deviations for the discussion of approach (a) vs (b).
    # --------------------------------------------------------------------------
    lka <- log(0.136); label("Absorption rate Ka (1/day)")                                # Xu 2019 Table 3, Ka row
    lcl <- log(0.260); label("Apparent linear clearance CLO/F (L/day)")                   # Xu 2019 Table 3, CLO/F row
    lvc <- log(2.08);  label("Apparent central volume Vc/F (L)")                          # Xu 2019 Table 3, Vc/F row
    lvp <- log(5.23);  label("Apparent peripheral volume Vp/F (L)")                       # Xu 2019 Table 3, Vp/F row
    lq  <- log(0.156); label("Apparent intercompartmental clearance Q/F (L/day)")         # Xu 2019 Table 3, Q/F row
    lvm <- log(8.06);  label("Maximum Michaelis-Menten elimination rate Vm (mg/day)")     # Xu 2019 Table 3, Vm row
    lkm <- log(0.939); label("Michaelis-Menten constant Km (mg/L)")                       # Xu 2019 Table 3, Km row

    # --------------------------------------------------------------------------
    # DAS28-CRP indirect-response PD model (Ma 2020 Table 3, final-model column)
    # Reference covariate values taken from Ma 2020 narrative (median of the
    # DAS28-CRP dataset): CRP 15.7 mg/L, BLPHYVAS 66, BLHAQ 1.75, WT 72.8 kg.
    # --------------------------------------------------------------------------
    lBase  <- log(6.06);  label("Typical DAS28-CRP baseline (unitless score)")            # Ma 2020 Table 3, BASE row (6.06)
    lEmax  <- 0.237;      label("Logit-transformed maximum drug effect on kin (unitless)") # Ma 2020 Table 3, Log(Emax) row; Emax = 1/(1+exp(-lEmax)) = 0.559, matching the paper's stated 55.9% maximum decrease
    lIC50  <- log(2.32);  label("Sarilumab concentration at 50% of Emax (mg/L)")          # Ma 2020 Table 3, IC50 row (2.32 mg/L)
    lKout  <- log(0.0264);label("First-order loss rate constant Kout (1/day)")            # Ma 2020 Table 3, Kout row (0.0264 day^-1)
    lPLB   <- log(0.991); label("Placebo/background DMARD effect, in sarilumab concentration units (mg/L)") # Ma 2020 Table 3, PLB row (0.991 mg/L)
    gamma  <- fixed(1);   label("Hill coefficient for the sigmoidal drug effect (unitless)") # Ma 2020 Table 3, gamma row (1 fixed)

    # Continuous-covariate exponents on BASE (power form; paper narrative describes
    # small effects for each covariate, clinically not meaningful per Ma 2020).
    e_crp_base    <- 0.0564; label("Power exponent of CRP/15.7 on BASE (unitless)")   # Ma 2020 Table 3, CRP on BASE row
    e_blphyvas_base <- 0.105;  label("Power exponent of BLPHYVAS/66 on BASE (unitless)")  # Ma 2020 Table 3, BLPHYVAS on BASE row
    e_blhaq_base    <- 0.0779; label("Power exponent of BLHAQ/1.75 on BASE (unitless)")   # Ma 2020 Table 3, BLHAQ on BASE row
    e_wt_base       <- 0.0522; label("Power exponent of WT/72.8 on BASE (unitless)")      # Ma 2020 Table 3, Weight on BASE row

    # Covariate effect on logit-transformed Emax (additive on log ratio scale).
    e_crp_lemax   <- 0.333;  label("Additive effect of log(CRP/15.7) on lEmax (unitless)") # Ma 2020 Table 3, CRP on Log(Emax) row

    # Binary-covariate multiplier on Kout.
    e_pricort_kout  <- 1.26;   label("Multiplicative effect on Kout for PRICORT = 1 (unitless)") # Ma 2020 Table 3, PRICORT on Kout row; paper narrative 0.0333 vs 0.0264 day^-1

    # --------------------------------------------------------------------------
    # IIV - PK (Xu 2019 Table 3, CV% converted to omega^2 = log(CV^2 + 1))
    #   Vm   CV 32.4% -> omega^2 = log(0.324^2 + 1) = 0.0998
    #   CLO/F CV 55.3% -> omega^2 = log(0.553^2 + 1) = 0.2669
    #   Vc/F CV 37.3% -> omega^2 = log(0.373^2 + 1) = 0.1302
    #   Ka   CV 32.1% -> omega^2 = log(0.321^2 + 1) = 0.0981
    #   Vm-CLO/F correlation -0.566 -> cov = -0.566 * sqrt(0.0998 * 0.2669) = -0.0924
    # --------------------------------------------------------------------------
    etalvm + etalcl ~ c(0.0998, -0.0924, 0.2669)  # Xu 2019 Table 3: Vm IIV 32.4% CV, CLO/F IIV 55.3% CV, Vm-CLO/F correlation -0.566
    etalvc ~ 0.1302                                # Xu 2019 Table 3: Vc/F IIV 37.3% CV
    etalka ~ 0.0981                                # Xu 2019 Table 3: Ka IIV 32.1% CV

    # --------------------------------------------------------------------------
    # IIV - DAS28-CRP PD (Ma 2020 Table 3, CV% converted to omega^2 = log(CV^2 + 1))
    #   BASE     CV  8.05% -> omega^2 = log(0.0805^2 + 1) = 0.00646
    #   lEmax    CV 71.2%  -> omega^2 = log(0.712^2 + 1) = 0.4105
    #   IC50     CV 158%   -> omega^2 = log(1.58^2  + 1) = 1.252
    #   Kout     CV 84.2%  -> omega^2 = log(0.842^2 + 1) = 0.5360
    #   PLB      CV 105%   -> omega^2 = log(1.05^2  + 1) = 0.7431
    # --------------------------------------------------------------------------
    etalBase ~ 0.00646   # Ma 2020 Table 3: BASE IIV 8.05% CV
    etalEmax ~ 0.4105    # Ma 2020 Table 3: Log(Emax) IIV 71.2% CV
    etalIC50 ~ 1.252     # Ma 2020 Table 3: IC50 IIV 158% CV
    etalKout ~ 0.5360    # Ma 2020 Table 3: Kout IIV 84.2% CV
    etalPLB  ~ 0.7431    # Ma 2020 Table 3: PLB IIV 105% CV

    # --------------------------------------------------------------------------
    # Residual errors.
    #   Sarilumab concentration: Xu 2019 log-additive sigma^2 = 0.395 -> propSd = 0.6285
    #   DAS28-CRP: Ma 2020 Table 3 reports an additive term of 0.647 mg/L; units in
    #   the table are ambiguous but the DAS28-CRP score is unitless and the paper's
    #   CWRES/VPC are on DAS28-CRP units. Treated here as additive on DAS28-CRP.
    # --------------------------------------------------------------------------
    CcpropSd    <- 0.6285; label("Proportional residual error on sarilumab concentration (fraction)") # Xu 2019 Table 3: residual sigma^2 = 0.395 (log-additive)
    das28addSd  <- 0.647;  label("Additive residual error on DAS28-CRP (score units)")                # Ma 2020 Table 3, additive residual row
  })
  model({
    # ------------------------------------------------------------------
    # 1. Individual PK parameters (Xu 2019 IIV only; no PK covariates).
    # ------------------------------------------------------------------
    vm <- exp(lvm + etalvm)
    km <- exp(lkm)
    vc <- exp(lvc + etalvc)
    cl <- exp(lcl + etalcl)
    ka <- exp(lka + etalka)
    q  <- exp(lq)
    vp <- exp(lvp)

    # ------------------------------------------------------------------
    # 2. Individual DAS28-CRP PD parameters (Ma 2020 Table 3).
    #    BASE carries four continuous covariates (power form).
    #    Emax on the logit scale carries one continuous covariate
    #      (additive in log-ratio space).
    #    Kout carries one binary covariate (multiplicative).
    # ------------------------------------------------------------------
    Base <- exp(lBase + etalBase) *
            (CRP    / 15.7)^e_crp_base    *
            (BLPHYVAS / 66  )^e_blphyvas_base *
            (BLHAQ    / 1.75)^e_blhaq_base    *
            (WT       / 72.8)^e_wt_base
    lEmax_i <- lEmax + etalEmax + e_crp_lemax * log(CRP / 15.7)
    Emax    <- 1 / (1 + exp(-lEmax_i))
    IC50    <- exp(lIC50 + etalIC50)
    Kout    <- exp(lKout + etalKout) * e_pricort_kout^PRICORT
    PLB     <- exp(lPLB  + etalPLB)
    Kin     <- Kout * Base

    # ------------------------------------------------------------------
    # 3. Two-compartment sarilumab PK with parallel linear + MM clearance
    #    from the central compartment (Xu 2019 Methods / Eq.).
    # ------------------------------------------------------------------
    Cc <- central / vc

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot -
                          (cl / vc) * central -
                          vm * Cc / (km + Cc) -
                          (q  / vc) * central +
                          (q  / vp) * peripheral1
    d/dt(peripheral1) <-  (q  / vc) * central - (q / vp) * peripheral1

    # ------------------------------------------------------------------
    # 4. Indirect-response DAS28-CRP with inhibition of kin.
    #    Effective sarilumab concentration includes the background DMARD
    #    placebo PLB: CeffP = Cc + PLB (Ma 2020 Fig. 1 caption equation).
    # ------------------------------------------------------------------
    CeffP <- Cc + PLB
    Eff   <- Emax * CeffP^gamma / (IC50^gamma + CeffP^gamma)

    das28(0)  <- Base
    d/dt(das28) <- Kin * (1 - Eff) - Kout * das28

    # ------------------------------------------------------------------
    # 5. Observation and error model.
    # ------------------------------------------------------------------
    Cc    ~ prop(CcpropSd)
    das28 ~ add(das28addSd)
  })
}

