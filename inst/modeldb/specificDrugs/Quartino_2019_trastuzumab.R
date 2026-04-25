Quartino_2019_trastuzumab <- function() {
  description <- "Two-compartment population PK model with parallel linear and Michaelis-Menten nonlinear elimination from the central compartment for intravenous trastuzumab (Herceptin) in patients with metastatic breast cancer, early breast cancer, advanced gastric cancer, or other solid tumors (Quartino 2019)"
  reference <- "Quartino AL, Li H, Kirschbrown WP, Mangat R, Wada DR, Garg A, Jin JY, Lum B. Population pharmacokinetic and covariate analyses of intravenous trastuzumab (Herceptin), a HER2-targeted monoclonal antibody, in patients with a variety of solid tumors. Cancer Chemother Pharmacol. 2019;83(2):329-340. doi:10.1007/s00280-018-3728-z"
  vignette <- "Quartino_2019_trastuzumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. Power effect on linear CL (exponent 0.967); reference 66 kg per Quartino 2019 CL covariate equation (typical patient). Weight-adjusted dosing (mg/kg) already partially compensates for body size; residual CL-vs-WT exponent captures the slight overcompensation described in Quartino 2019 Results 'Assessment of the impact of identified covariates on PK exposure'.",
      source_name        = "Wt"
    ),
    AST = list(
      description        = "Baseline serum aspartate aminotransferase activity",
      units              = "IU/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. Power effect on linear CL (exponent 0.205); reference 24 IU/L per Quartino 2019 CL covariate equation (typical patient; median of the analysis population). Source column 'SGOT' (the legacy clinical-chemistry name for AST) maps to the canonical AST covariate.",
      source_name        = "SGOT"
    ),
    ALB = list(
      description        = "Baseline serum albumin concentration",
      units              = "g/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value (US convention, g/dL). Power effect on linear CL (exponent -0.998, i.e. approximately inversely proportional); reference 4 g/dL per Quartino 2019 CL covariate equation (typical patient). Source column 'ALBU' maps to the canonical ALB covariate.",
      source_name        = "ALBU"
    ),
    LMET = list(
      description        = "Baseline presence of liver metastases",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no liver metastases at baseline)",
      notes              = "Time-fixed baseline indicator. Exponential effect on linear CL (coefficient 0.152, i.e. a ~16.4% CL increase for LMET-positive subjects). Reference LMET = 0 per Quartino 2019 CL covariate equation. Source column 'LMET' maps directly to the canonical LMET covariate.",
      source_name        = "LMET"
    ),
    TUMTP_GC = list(
      description        = "Tumor-type indicator for advanced gastric cancer",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (MBC, EBC, HV, or 'Others' group)",
      notes              = "Quartino 2019 decomposes primary tumor type (TTYPE) into three linear-CL groups (MBC/EBC/HV reference, AGC, Others) and two Vc groups (non-AGC reference, AGC). The AGC indicator enters as a per-group typical-value switch on both linear CL and Vc (different typical values selected by indicator, not a multiplicative exponential / power effect). Derived from the source categorical column TTYPE as TUMTP_GC = as.integer(TTYPE == 'AGC').",
      source_name        = "TTYPE"
    ),
    TUMTP_OTH = list(
      description        = "Tumor-type indicator for 'Others' group (NSCLC and other non-breast, non-gastric solid tumors pooled in Quartino 2019)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (MBC, EBC, HV, or AGC)",
      notes              = "Quartino 2019 pools 107/1582 patients with NSCLC, prostate, ovarian, and other miscellaneous solid tumors into a single 'Others' category for linear CL (distinct typical value theta8 = 0.148 L/day). No Vc effect. Derived from the source categorical column TTYPE as TUMTP_OTH = as.integer(TTYPE == 'Others'). A subject has at most one of TUMTP_GC / TUMTP_OTH set to 1; both zero = MBC/EBC/HV reference.",
      source_name        = "TTYPE"
    )
  )

  population <- list(
    n_subjects       = 1582L,
    n_studies        = 18L,
    n_observations   = 26040L,
    phase_mix        = "Pooled phase I, II, and III trials (all using the innovator trastuzumab Herceptin; no biosimilar data).",
    age_median       = "53 years",
    weight_median    = "66 kg",
    sex_female_pct   = 82.7,
    race_ethnicity   = c(`non-Asian` = 83.4, Asian = 16.6),
    disease_state    = "Metastatic breast cancer (MBC, 810/1582), early breast cancer (EBC, 391/1582), advanced gastric cancer (AGC, 274/1582), non-small cell lung cancer or other solid tumors (107/1582), and healthy volunteers (HV, 6/1582).",
    dose_range       = "Weekly (qw): 4 mg/kg IV loading dose followed by 2 mg/kg maintenance. Every 3 weeks (q3w): 8 mg/kg IV loading dose followed by 6 mg/kg maintenance. 917 patients on q3w, 643 on qw, 28 single-dose. 1188 single-agent; remainder combination with anthracyclines, docetaxel, paclitaxel, cisplatin, or other chemotherapy.",
    regions          = "Global; 18 pooled phase I-III studies.",
    ecog_status      = "94.5% ECOG performance status 0 or 1.",
    reference_subject = "66 kg, AST (SGOT) 24 IU/L, ALB 4 g/dL, no liver metastases, MBC/EBC/HV primary tumor type (the typical patient Quartino 2019 uses for covariate-impact assessment and for the 'typical BC' simulation in Table 2).",
    notes            = "Baseline demographics per Quartino 2019 Results section 'Patient population' and Online Resource 6. The final PK dataset contained 1582 patients and 26,040 serum concentrations across 18 studies; 1588 patients passed PK data handling (316 samples or patients excluded as outliers or BLQ). Among AGC patients, 12.6% had a prior gastrectomy (not a retained model covariate). Shed-antigen ECD-HER2 (SHED) was an exploratory analysis only and is NOT part of the final model."
  )

  ini({
    # Structural parameters (Quartino 2019 Table 1, final FOCEI estimates).
    # Reference subject: 66 kg MBC/EBC/HV patient, AST 24 IU/L, ALB 4 g/dL,
    # no liver metastases. Concentration in the central compartment is
    # Cc = central / vc, with dose in mg and volumes in L -> Cc in mg/L
    # (= ug/mL). Vmax is mg/day and Km is mg/L, so Vmax*Cc/(Km+Cc) is mg/day
    # and matches the linear-elimination flux kel*central (mg/day).
    lcl  <- log(0.127); label("Linear CL for the MBC/EBC/HV reference subject (L/day)")        # Quartino 2019 Table 1, theta1
    lvc  <- log(2.62);  label("Central volume of distribution for non-AGC subject (L)")        # Quartino 2019 Table 1, theta2
    lq   <- log(0.544); label("Intercompartmental clearance Q (L/day)")                        # Quartino 2019 Table 1, theta3
    lvp  <- log(2.97);  label("Peripheral volume of distribution Vp (L)")                      # Quartino 2019 Table 1, theta4
    lvmax <- log(8.81); label("Maximum nonlinear (Michaelis-Menten) elimination rate Vmax (mg/day)") # Quartino 2019 Table 1, theta5
    lkm  <- log(8.92);  label("Michaelis-Menten constant Km (mg/L = ug/mL)")                   # Quartino 2019 Table 1, theta6

    # Tumor-type linear-CL typical values (Quartino 2019 Table 1 and the CL
    # covariate equation in Results). Quartino 2019 models TTYPE as a
    # per-group typical-value switch on CL (three levels: MBC/EBC/HV reference
    # = theta1; AGC = theta9; Others = theta8) rather than a multiplicative
    # covariate effect, so the three linear CLs are log-stored separately.
    lcl_agc <- log(0.176); label("Linear CL for AGC subjects (L/day)")                         # Quartino 2019 Table 1, theta9
    lcl_oth <- log(0.148); label("Linear CL for 'Others' tumor types (L/day)")                 # Quartino 2019 Table 1, theta8

    # Central-volume AGC typical value (Quartino 2019 Table 1 and Vc covariate
    # equation). AGC is the only tumor-type group with a distinct Vc; all
    # non-AGC groups share theta2.
    lvc_agc <- log(3.63);  label("Central volume of distribution for AGC subjects (L)")       # Quartino 2019 Table 1, theta13

    # Covariate effects on linear CL (Quartino 2019 Table 1 and CL covariate
    # equation). Continuous covariates enter as (COV / ref)^exponent with
    # reference WT = 66 kg, AST = 24 IU/L, ALB = 4 g/dL. The LMET indicator
    # enters exponentially: exp(theta12 * LMET).
    e_wt_cl   <-  0.967; label("Power exponent of WT on linear CL (unitless)")                 # Quartino 2019 Table 1, theta7
    e_ast_cl  <-  0.205; label("Power exponent of AST (SGOT) on linear CL (unitless)")         # Quartino 2019 Table 1, theta10
    e_alb_cl  <- -0.998; label("Power exponent of ALB on linear CL (unitless)")                # Quartino 2019 Table 1, theta11
    e_lmet_cl <-  0.152; label("Exponential coefficient of liver-metastases indicator on linear CL (unitless)") # Quartino 2019 Table 1, theta12

    # Inter-individual variability. Quartino 2019 Table 1 reports omega as %CV
    # on log-normal parameters; convert via omega^2 = log(CV^2 + 1). CL and
    # Vc share a correlated block with off-diagonal covariance Omega_CL,Vc =
    # 0.0230 (Table 1 footnote c). IIV on Vp and Km is retained per the paper
    # despite moderate-to-high shrinkage (Vp 22.9%, Km 44.0%).
    etalcl + etalvc ~ c(0.149110,
                        0.0230, 0.058756)   # CL 40.1% CV, Vc 24.6% CV, cov 0.0230 -- Quartino 2019 Table 1
    etalvp ~ 0.219156                       # Vp 49.5% CV  -- Quartino 2019 Table 1
    etalkm ~ 1.075719                       # Km 139%  CV  -- Quartino 2019 Table 1

    # Residual error: combined additive + proportional (Quartino 2019 Table 1).
    # Table 1 reports sigma1 = 19.7% (proportional SD) and sigma2 = 1.38 ug/mL
    # (additive SD). Table 1 footnote d confirms these are SDs on the
    # natural scale: the tabulated RSE column for residual variability is
    # relative to the estimated variance (sigma^2), so the point estimate
    # itself is the SD.
    propSd <- 0.197; label("Proportional residual error (fraction)")                           # Quartino 2019 Table 1, sigma1
    addSd  <- 1.38;  label("Additive residual error (ug/mL)")                                  # Quartino 2019 Table 1, sigma2
  })

  model({
    # Per-subject typical linear CL by tumor-type group (Quartino 2019 CL
    # covariate equation in Results). TUMTP_GC = 1 picks the AGC typical CL;
    # TUMTP_OTH = 1 picks the Others typical CL; both zero = MBC/EBC/HV
    # reference. A subject has at most one indicator set to 1.
    lcl_typ <- lcl +
      (lcl_agc - lcl) * TUMTP_GC +
      (lcl_oth - lcl) * TUMTP_OTH

    # Per-subject typical central volume (Quartino 2019 Vc covariate equation).
    # Only AGC receives a distinct typical Vc; all non-AGC groups share theta2.
    lvc_typ <- lvc + (lvc_agc - lvc) * TUMTP_GC

    # Individual linear CL with weight, AST, and albumin power effects plus a
    # liver-metastases exponential effect (Quartino 2019 CL covariate equation).
    cl <- exp(lcl_typ + etalcl) *
      (WT  / 66)^e_wt_cl *
      (AST / 24)^e_ast_cl *
      (ALB /  4)^e_alb_cl *
      exp(e_lmet_cl * LMET)

    vc   <- exp(lvc_typ + etalvc)
    q    <- exp(lq)
    vp   <- exp(lvp + etalvp)
    vmax <- exp(lvmax)
    km   <- exp(lkm + etalkm)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    Cc <- central / vc

    # Two-compartment IV model with parallel linear and Michaelis-Menten
    # nonlinear elimination from the central compartment (Quartino 2019 Results:
    # "a two-compartment model with parallel linear and nonlinear (Michaelis-
    # Menten) elimination from the central compartment"). Linear flux kel *
    # central and nonlinear flux Vmax * Cc / (Km + Cc) are both in mg/day; dose
    # is mg so the central state carries mg and Cc is mg/L (= ug/mL).
    d/dt(central)     <- -kel * central - vmax * Cc / (km + Cc) - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                          k12 * central - k21 * peripheral1

    Cc ~ add(addSd) + prop(propSd)
  })
}
