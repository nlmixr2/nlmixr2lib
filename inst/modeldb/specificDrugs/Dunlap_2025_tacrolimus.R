Dunlap_2025_tacrolimus <- function() {
  description <- paste0(
    "Two-compartment population pharmacokinetic model for oral immediate-",
    "release tacrolimus in adult allogeneic hematopoietic cell transplant ",
    "(allo-HCT) recipients (Dunlap 2025): first-order absorption with ",
    "bioavailability fixed at 1; allometric (TBW/70 kg) scaling fixed at 0.75 ",
    "on CL/F and Q/F and at 1 on V1/F and V2/F; exponential CYP3A5 ",
    "intermediate / normal metabolizer phenotype effect on CL/F (CYP3A5 IM or ",
    "NM have ~2.14-fold higher CL/F than CYP3A5 PM); exponential reduced-",
    "intensity-conditioning effect on CL/F (RIC recipients have ~37% lower ",
    "CL/F than myeloablative-conditioning recipients); inter-individual ",
    "variability on V1/F, CL/F, and V2/F; and an additive residual error of ",
    "2.51 ng/mL on the linear concentration scale."
  )
  reference <- paste0(
    "Dunlap TC, Zhu J, Weiner DL, Kemper RM, DeVane SC, Ma F, et al. ",
    "A Tacrolimus Population Pharmacokinetic Model for Adult Allogeneic ",
    "Hematopoietic Cell Transplant Recipients Provides Clinical Opportunities ",
    "for Precision Dosing. Clin Pharmacokinet. 2025;64(11):1621-1637. ",
    "doi:10.1007/s40262-025-01529-w."
  )
  vignette <- "Dunlap_2025_tacrolimus"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Time-fixed at baseline in Dunlap 2025 (Table 1 reports a single ",
        "TBW per subject). Allometric power scaling on CL/F and Q/F with ",
        "exponent 0.75 fixed (theoretical), and on V1/F and V2/F with ",
        "exponent 1 fixed; reference 70 kg (Dunlap 2025 Methods 2.4 / Eqs. ",
        "5-8 / Table 2). Cohort median 84 kg [IQR 71-97]."
      ),
      source_name        = "TBW"
    ),
    CYP3A5_EXPR = list(
      description        = paste0(
        "CYP3A5 expresser indicator: 1 if the patient has a CYP3A5 ",
        "intermediate-metabolizer (IM) or normal-metabolizer (NM) phenotype ",
        "(at least one functional CYP3A5*1 allele -- diplotypes *1/*1 or ",
        "*1/*3, *1/*6, *1/*7), 0 if the patient is a poor metabolizer (PM) ",
        "with two non-functional alleles (*3/*3, *3/*6, *3/*7, *6/*6, ",
        "*6/*7, or *7/*7). Time-fixed germline genotype."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP3A5 PM, the reference category in Dunlap 2025 Eq. 6 / Table 2)",
      notes              = paste0(
        "Dunlap 2025 reports a three-level CYP3A5 metabolizer phenotype ",
        "(PM / IM / NM) per CPIC nomenclature (Dunlap 2025 Methods 2.2). ",
        "After fitting both a three-level full covariate model and a ",
        "pooled two-level reduced covariate model, the authors retained the ",
        "single binary `IM or NM` contrast vs PM in the final RCM (Dunlap ",
        "2025 Section 3.2 / Table 2: TVCL/F~CYP3A5 IM or NM = 2.14, no ",
        "separate IM vs NM coefficient). The CYP3A5_EXPR canonical encodes ",
        "the IM-or-NM-vs-PM dichotomy directly: CYP3A5_EXPR = 1 for IM and ",
        "NM (any *1 allele), 0 for PM (no *1 allele). Cohort distribution ",
        "(Dunlap 2025 Table 1): PM 70%, IM 25%, NM 5%; CYP3A5_EXPR = 1 in ",
        "30% of subjects."
      ),
      source_name        = "CYP3A5"
    ),
    COND_RIC = list(
      description        = paste0(
        "Reduced-intensity conditioning regimen indicator: 1 if the patient ",
        "received reduced-intensity conditioning (RIC) chemotherapy prior to ",
        "allo-HCT, 0 if myeloablative conditioning (MAC) was used."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (MAC, myeloablative conditioning)",
      notes              = paste0(
        "Time-fixed per subject (the conditioning regimen is completed ",
        "before transplant and is invariant across the post-transplant ",
        "tacrolimus PK observations). Per Dunlap 2025 Methods 2.4, the ",
        "binary RIC indicator follows the institutional UNCMC classification ",
        "of conditioning intensity. Cohort split (Dunlap 2025 Table 1): MAC ",
        "49%, RIC 51%."
      ),
      source_name        = "RIC"
    )
  )

  population <- list(
    n_subjects        = 290L,
    n_studies         = 2L,
    n_observations    = 906L,
    age_range         = "44-63 years (overall IQR; UNC16-1480 IQR 44-62, UNC19-3328 IQR 44-66)",
    age_median        = "54 years (overall; UNC16-1480 54, UNC19-3328 61)",
    weight_range      = "71-97 kg (overall IQR)",
    weight_median     = "84 kg (overall; UNC16-1480 85, UNC19-3328 79)",
    sex_female_pct    = 43.0,
    race_ethnicity    = c(White = 84.0, Black = 11.0, Other = 5.0),
    disease_state     = paste0(
      "Adult allogeneic hematopoietic cell transplant (allo-HCT) recipients ",
      "treated at the University of North Carolina Medical Center (UNCMC) ",
      "for malignant or non-malignant haematologic disease, receiving oral ",
      "immediate-release tacrolimus capsules as part of acute graft-versus-",
      "host disease (aGVHD) prophylaxis. Primary diagnoses: acute leukemia ",
      "55%, MDS 17%, lymphoma 12%, chronic leukemia 7%, MPN 6%, aplastic ",
      "anaemia 2%, myeloma 1%."
    ),
    dose_range        = paste0(
      "Per UNCMC institutional protocol: oral immediate-release tacrolimus ",
      "0.045 mg/kg twice daily for the first two days starting three days ",
      "prior to transplant (D-3), followed by 0.03 mg/kg twice daily ",
      "thereafter (Methods 2.2). Therapeutic drug monitoring deviations ",
      "from this protocol were retained in the analysis dataset."
    ),
    regions           = "United States (single centre, University of North Carolina Medical Center, Chapel Hill, NC)",
    cyp3a5_distribution = "PM 70% (n=202), IM 25% (n=73), NM 5% (n=15) per CPIC nomenclature; *1/*1, *1/*3, *3/*3 etc. determined by TaqMan assays for CYP3A5*1, *3, *6, *7 (Methods 2.2).",
    conditioning_distribution = "MAC 49% (n=143), RIC 51% (n=147)",
    hla_status        = "Full match 80% (n=231), mismatch 18% (n=54), haploidentical 2% (n=5)",
    donor_status      = "Related 30% (n=88), unrelated 70% (n=202)",
    stem_cell_source  = "Peripheral blood stem cells 94% (n=273), bone marrow 6% (n=16), cord <1% (n=1)",
    n_observations_unc16_1480 = 252L,
    n_observations_unc19_3328 = 652L,
    notes             = paste0(
      "Pooled retrospective (UNC16-1480, n=252) and prospective (UNC19-3328, ",
      "n=38) clinical pharmacology cohorts (clinicaltrials.gov NCT04645667). ",
      "UNC16-1480 contributed a single therapeutic-drug-monitoring trough ",
      "per subject collected on D-1; UNC19-3328 contributed a richer ",
      "sampling schedule (D-2, D-1, D0 trough draws plus D0 post-dose ",
      "samples at 0.5, 1, 2, 4, 6, and 10 h, mean 17.2 samples per subject). ",
      "Bioanalytical assay: tacrolimus by LC-MS/MS in whole blood at the ",
      "CLIA-certified UNCMC McLendon Special Chemistries Lab (analytic ",
      "range 1-40 ng/mL). Patients on concomitant strong CYP3A4 inhibitors ",
      "(e.g. posaconazole) for antifungal prophylaxis, IV tacrolimus, ",
      "tacrolimus oral suspension, or post-transplant cyclophosphamide ",
      "conditioning were excluded. Baseline demographics in Dunlap 2025 ",
      "Table 1."
    )
  )

  ini({
    # ----- Structural PK (Dunlap 2025 Table 2 reduced-covariate-model column) -----
    # Two-compartment first-order oral absorption; F = 1 fixed (Methods 2.3,
    # Section 3.2: "All subsequent model runs were performed with allometric
    # relationships fixed to their theoretical values (i.e., 1 for TVV1/F and
    # TVV2/F, and 0.75 for TVCL/F and TVQ/F)"). Reference subject for the
    # tabulated typical values: 70 kg, CYP3A5 PM phenotype, MAC conditioning
    # (Dunlap 2025 Section 3.2: "The estimated TVKA/F, TVV1/F, TVCL/F, TVV2/F,
    # and TVQ/F for a 70 kg CYP3A5 PM subject receiving MAC chemotherapy were
    # 0.50/h, 150 L, 23 L/h, 1153 L, and 43 L/h").
    lka <- log(0.50);   label("Absorption rate (Ka, 1/h)")  # Dunlap 2025 Table 2 RCM TVKA = 0.50 1/h [95% CI 0.38-0.66]
    lcl <- log(23);     label("Apparent CL/F at reference covariates (TBW=70 kg, CYP3A5 PM, MAC) (L/h)")  # Dunlap 2025 Table 2 RCM TVCL/F = 23 L/h [95% CI 17-30]
    lvc <- log(150);    label("Apparent central volume V1/F at reference TBW=70 kg (L)")  # Dunlap 2025 Table 2 RCM TVV1/F = 150 L [95% CI 104-215]
    lq  <- log(43);     label("Apparent intercompartmental clearance Q/F at reference TBW=70 kg (L/h)")  # Dunlap 2025 Table 2 RCM TVQ/F = 43 L/h [95% CI 37-50]
    lvp <- log(1153);   label("Apparent peripheral volume V2/F at reference TBW=70 kg (L)")  # Dunlap 2025 Table 2 RCM TVV2/F = 1153 L [95% CI 774-1716]

    # ----- Allometric exponents (fixed at theoretical values) -----
    # Dunlap 2025 Table 2 RCM: TVV1/F~TBW = 1 (FIX), TVCL/F~TBW = 0.75 (FIX),
    # TVV2/F~TBW = 1 (FIX), TVQ/F~TBW = 0.75 (FIX); reference 70 kg.
    e_wt_cl_q   <- 0.75;  label("Allometric exponent of (TBW/70) on CL/F and Q/F (unitless; fixed)")  # Dunlap 2025 Table 2 RCM allometry on CL/F and Q/F (FIX)
    e_wt_vc_vp  <- 1.00;  label("Allometric exponent of (TBW/70) on V1/F and V2/F (unitless; fixed)")  # Dunlap 2025 Table 2 RCM allometry on V1/F and V2/F (FIX)

    # ----- Covariate effects on CL/F (Dunlap 2025 Eq. 6 / Table 2 RCM column) -----
    # Eq. 6: CL/F_i = TVCL/F * (TBW/70)^0.75 * (TVCL/F~RIC)^RIC *
    #               (TVCL/F~IM_or_NM)^IMNM * exp(eta_CL).
    # Both covariate factors are reported in the natural domain (i.e., the
    # tabulated value is the multiplicative shift applied when the indicator
    # is 1). In log-domain notation: theta_RIC = log(0.63), theta_IM_or_NM =
    # log(2.14). Reference category is MAC (RIC = 0) and CYP3A5 PM
    # (CYP3A5_EXPR = 0).
    e_cyp3a5_expr_cl <- 2.14;  label("CYP3A5 IM or NM multiplicative factor on CL/F (expressers have ~114% higher CL/F than PM)")  # Dunlap 2025 Table 2 RCM TVCL/F~CYP3A5 IM or NM = 2.14 [95% CI 1.67-2.74]
    e_cond_ric_cl    <- 0.63;  label("Reduced-intensity-conditioning multiplicative factor on CL/F (RIC have ~37% lower CL/F than MAC)")  # Dunlap 2025 Table 2 RCM TVCL/F~RIC = 0.63 [95% CI 0.51-0.77]

    # ----- Inter-individual variability -----
    # Dunlap 2025 Table 2 RCM IIV reported as %CV: V1/F 95%, CL/F 55%, V2/F 66%.
    # KA and Q/F do not carry IIV in the RCM (Eqs. 4 and 8 of Dunlap 2025).
    # Convert to log-scale variance via omega^2 = log(CV^2 + 1):
    #   V1/F  CV 95% -> log(0.95^2 + 1) = 0.65823
    #   CL/F  CV 55% -> log(0.55^2 + 1) = 0.26282
    #   V2/F  CV 66% -> log(0.66^2 + 1) = 0.36185
    # Diagonal block (Dunlap 2025 does not report off-diagonal correlations
    # between IIV terms in the RCM).
    etalvc ~ 0.65823  # Dunlap 2025 Table 2 RCM IIV V1/F = 95% CV (shrinkage 65%)
    etalcl ~ 0.26282  # Dunlap 2025 Table 2 RCM IIV CL/F = 55% CV (shrinkage 27%)
    etalvp ~ 0.36185  # Dunlap 2025 Table 2 RCM IIV V2/F = 66% CV (shrinkage 51%)

    # ----- Residual unexplained variability -----
    # Dunlap 2025 Methods 2.3: "RUV was estimated in the normal domain";
    # Table 2 RCM reports an additive RUV of 2.51 ng/mL (reported as standard
    # deviation per Table 2 footnote), corresponding to nlmixr2 add() on the
    # linear concentration scale.
    addSd <- 2.51;  label("Additive residual error (ng/mL)")  # Dunlap 2025 Table 2 RCM additive RUV = 2.51 ng/mL [95% CI 2.37-2.65]
  })

  model({
    # ----- 1. Derived covariate terms -----
    wt_cl_q  <- (WT / 70) ^ e_wt_cl_q
    wt_vc_vp <- (WT / 70) ^ e_wt_vc_vp

    # ----- 2. Individual PK parameters -----
    # Reference subject: 70 kg, CYP3A5 PM (CYP3A5_EXPR = 0), MAC (COND_RIC = 0).
    # The covariate factors are applied as `theta ^ indicator`, which
    # collapses to 1 at the reference category and to theta when the
    # indicator is 1, reproducing Dunlap 2025 Eq. 6 directly.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * wt_cl_q *
          e_cyp3a5_expr_cl ^ CYP3A5_EXPR *
          e_cond_ric_cl    ^ COND_RIC
    vc <- exp(lvc + etalvc) * wt_vc_vp
    q  <- exp(lq)            * wt_cl_q
    vp <- exp(lvp + etalvp) * wt_vc_vp

    # ----- 3. Micro-constants -----
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ----- 4. ODE system (2-cmt with first-order oral absorption) -----
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # ----- 5. Bioavailability -----
    # F = 1 fixed (Dunlap 2025 Methods 2.3); no explicit f(depot) statement.

    # ----- 6. Observation and error -----
    # Dose in mg, central amount in mg, vc in L -> mg/L; multiply by 1000 to
    # report tacrolimus whole-blood concentration in ng/mL, the units of the
    # bioanalytical assay (Methods 2.2) and Table 2 / Figure 5.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd)
  })
}
