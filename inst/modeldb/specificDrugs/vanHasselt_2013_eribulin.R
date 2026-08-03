vanHasselt_2013_eribulin <- function() {
  description <- "Three-compartment IV PK driver coupled with a Friberg-style semi-physiological PD model for eribulin-induced neutropenia in adult patients with late-stage metastatic breast cancer (van Hasselt 2013). The PK layer parameters are FIXED from the Majid 2014 popPK analysis (J Clin Pharmacol 54:1134); van Hasselt 2013 cites the same Eisai popPK analysis as unpublished data on file and does not report PK parameter values in the paper. CL depends on body weight (allometric 0.75), serum albumin, alkaline phosphatase, and total bilirubin; V1, V2, V3 scale linearly with body weight; Q2 and Q3 scale allometrically with body weight. The PD layer (proliferation + three transit compartments + circulating neutrophils + feedback) is estimated on 1579 patients / 23 427 ANC measurements from a pool of 12 phase I / II / III studies (Table 3 final covariate model): ANC0 = 4.03e9 cells/L, MTT = 109 h, Gamma = 0.216, SLOPE = 0.0451 L/ug (linear drug effect). Albumin, bilirubin, alkaline phosphatase, LDH, and G-CSF-receival covariates enter MTT; albumin, AST, and G-CSF-receival covariates enter SLOPE. IIV on ANC0, MTT, Gamma, and SLOPE (diagonal; off-diagonal covariances were not computationally feasible per the paper). Proportional residual error on circulating ANC. Eribulin doses must be supplied in milligrams of eribulin FREE BASE (conversion factor 1.23/1.4 from the mesilate dose)."
  reference <- paste(
    "van Hasselt JGC, Gupta A, Hussein Z, Beijnen JH, Schellens JHM, Huitema ADR. (2013).",
    "Population pharmacokinetic-pharmacodynamic analysis for eribulin mesilate-associated neutropenia.",
    "Br J Clin Pharmacol 76(3):412-424. doi:10.1111/bcp.12143.",
    "PK structure and parameter values FIXED from the Eisai popPK analysis published as",
    "Majid O, Gupta A, Reyderman L, Olivo M, Hussein Z (2014).",
    "Population pharmacometric analyses of eribulin in patients with locally advanced or metastatic breast cancer previously treated with anthracyclines and taxanes.",
    "J Clin Pharmacol 54(10):1134-1143. doi:10.1002/jcph.315;",
    "van Hasselt 2013 refers to this analysis as unpublished data on file (Dr Z. Hussein, Eisai Ltd)",
    "and does not report PK parameter values in the paper itself.",
    "PD structure based on Friberg LE, Henningsson A, Maas H, Nguyen L, Karlsson MO (2002).",
    "Model of chemotherapy-induced myelosuppression with parameter consistency across drugs.",
    "J Clin Oncol 20(24):4713-4721. doi:10.1200/JCO.2002.02.140;",
    "see modellib('Friberg_2002_paclitaxel') and modellib('Kawamura_2018_eribulin')",
    "for the Friberg-family template and the Kawamura 2018 postmarketing-surveillance re-fit of the same PD structure.",
    sep = " "
  )
  vignette <- "vanHasselt_2013_eribulin"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "mg/L",
    anc           = "cells/uL"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling on the Majid 2014 PK parameters with reference 68.7 kg (per the Majid 2014 CL / Q / V equations reproduced in Kawamura 2018 section 2.3): CL, Q2, Q3 scale as (WT/68.7)^0.75; V1, V2, V3 scale linearly as (WT/68.7). van Hasselt 2013 Table 2 reports the pooled-cohort bodyweight median 67.7 kg (IQR 59.0-77.6, footnote scaling value 70 kg for the PK model). WT is used only by the PK layer; the van Hasselt 2013 PD model does not retain bodyweight as a covariate in the final MTT / SLOPE relationships.",
      source_name        = "WT"
    ),
    ALB = list(
      description        = "Serum albumin at baseline (canonical SI unit g/L).",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters BOTH the PK CL term (Majid 2014: (ALB_gdL/4.0)^0.946 with reference 4.0 g/dL) AND the van Hasselt 2013 PD MTT and SLOPE terms (Table 3 final model: (ALB_gdL/4)^0.374 on MTT, (ALB_gdL/4)^0.763 on SLOPE). Canonical register unit is g/L (SI); van Hasselt 2013 Table 2 reports median 3.90 g/dL (IQR 3.6-4.27 g/dL, scaling value 4 g/dL per Table 2 footnote). The model() body applies an inline conversion `alb_gdL <- ALB * 0.1` so the g/L canonical column drives both the PK and PD terms at their published g/dL calibration.",
      source_name        = "ALB"
    ),
    ALP = list(
      description        = "Serum alkaline phosphatase at baseline.",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters BOTH the PK CL term (Majid 2014: (ALP/132)^-0.209 with reference 132 U/L) AND the van Hasselt 2013 PD MTT term (Table 3 final model: (ALP/100)^-0.0337 with scaling value 100 U/L per Table 2 footnote). van Hasselt 2013 Table 2 reports median 118 U/L (IQR 82.0-206 U/L). The two references differ because the PK analysis was fit on a different (earlier) pooled cohort than the PD analysis; both are applied verbatim.",
      source_name        = "ALP"
    ),
    TBILI = list(
      description        = "Total serum bilirubin at baseline (canonical SI unit umol/L).",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters BOTH the PK CL term (Majid 2014: (TBILI_mgdL/0.5)^-0.180 with reference 0.5 mg/dL) AND the van Hasselt 2013 PD MTT term (Table 3 final model: (TBILI_mgdL/2)^-0.046 with scaling value 2 mg/dL per Table 2 footnote). Canonical register unit is umol/L (SI); van Hasselt 2013 Table 2 reports median 0.50 mg/dL (IQR 0.40-0.70 mg/dL). The model() body applies an inline conversion `tbili_mgdL <- TBILI / 17.1` so the umol/L canonical column drives both the PK and PD terms at their published mg/dL calibration. Source column name in van Hasselt 2013 is BILI; the canonical register uses TBILI.",
      source_name        = "BILI"
    ),
    LDH = list(
      description        = "Serum lactate dehydrogenase at baseline.",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "van Hasselt 2013 Table 3 final covariate model: LDH enters MTT as (LDH/238)^-0.0561 with scaling value 238 U/L per Table 2 footnote. Cohort median 328 U/L (IQR 211-486 U/L). Interpreted as a disease-burden / cell-turnover proxy in this metastatic-breast-cancer cohort (see canonical register notes for LDH).",
      source_name        = "LDH"
    ),
    AST = list(
      description        = "Serum aspartate aminotransferase at baseline.",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "van Hasselt 2013 Table 3 final covariate model: AST enters SLOPE as (AST/30)^0.119. The scaling value is not explicitly listed in Table 2 footnote (which enumerates deviations from median for the five covariates WT, ALB, ALP, BILI, LDH); AST scaling therefore defaults to the cohort median 30.0 U/L per Table 2 (IQR 22.0-46.0 U/L). Positive exponent -> elevated AST increases the linear drug-effect slope, i.e. patients with worse liver function are more sensitive to eribulin-induced myelosuppression.",
      source_name        = "AST"
    ),
    CONMED_GCSF = list(
      description        = "Concomitant granulocyte colony-stimulating factor (G-CSF) treatment indicator (1 = received G-CSF; 0 = did not).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant G-CSF treatment)",
      notes              = "van Hasselt 2013 Table 2 counts 382 / 1302 (yes / no) patients receiving G-CSF (22.7% of the cohort). Table 3 final covariate model retains G-CSF as a proportional / dichotomous effect on BOTH MTT and SLOPE per Eq. 10: MTT_ind = MTT_typ * ... * 0.883^CONMED_GCSF (G-CSF shortens MTT by ~12%) and SLOPE_ind = SLOPE_typ * ... * 1.3^CONMED_GCSF (G-CSF increases SLOPE by ~30%). van Hasselt 2013 treats this as a subject-level indicator in the covariate model (the paper does not describe per-cycle time-varying encoding), so users should supply CONMED_GCSF as a per-subject constant (1 for any subject who received G-CSF at any point during the observation period, 0 otherwise). Distinct from CSF1 (colony-stimulating factor 1 / M-CSF, a target-engagement biomarker); see canonical register CONMED_GCSF entry.",
      source_name        = "G-CSF"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 1579L,
    n_studies        = 12L,
    age_range        = "median 58.0 years (IQR 49.0-66.0)",
    weight_range     = "median 67.7 kg (IQR 59.0-77.6 kg)",
    sex_female_pct   = 85.6,
    race_ethnicity   = c(Caucasian = 79.4, `Black/African American` = 5.3, `Asian/Pacific Islander` = 1.5, Japanese = 6.1, `Other/unknown` = 8.4),
    disease_state    = "Late-stage metastatic breast cancer. Pooled data set of 12 phase I, II, and III studies of eribulin mesilate (Table 1 of van Hasselt 2013 lists the individual studies). 1579 patients contributed 23 427 absolute neutrophil count (ANC) measurements to the PD analysis; PK data were available for 428 patients (27%). Prior therapy exposure was high: previous chemotherapy 1527/61 (yes/no), previous Pt-containing chemotherapy 446/1142, previous radiotherapy 1212/375, previous hormonal therapy 923/761. G-CSF was administered in 382/1302 patients (22.7%).",
    dose_range       = "Eribulin mesilate 1.4 mg/m^2 IV per the approved dosing schedule (day 1 and day 8 of a 21-day treatment cycle) evaluated in phase II / III; phase I studies explored 0.6-4 mg/m^2 dose escalation. Doses must be supplied to this model in milligrams of eribulin FREE BASE; conversion from the mesilate dose is D_free_base = D_mesilate * 1.23/1.4 per the Majid 2014 / Kawamura 2018 conversion factor.",
    regions          = "Not tabulated per se; the pooled 12-study data set includes phase I / II / III trials conducted across North America, Europe, and Japan (Japanese subgroup n = 96, per Table 2 ethnicity breakdown).",
    prior_therapy    = "Highly pretreated cohort: 96% received previous chemotherapy, 28% received prior Pt-containing chemotherapy, 76% received prior radiotherapy, 58% received prior hormonal therapy, 10% received blood transfusions. Median 4 prior chemotherapy lines is not explicitly reported in this paper but is characteristic of the Phase III EMBRACE population from which the largest study was drawn.",
    notes            = "Semi-physiological Friberg-style myelosuppression model with linear drug effect (E = SLOPE * C) selected over the E-max alternative on the basis of superior parameter precision and a high correlation between EC50 and Emax fixed effects. First-order estimation method was used (first-order conditional estimation with interaction was not computationally feasible). Off-diagonal IIV covariances and IOV were not computationally feasible. Bootstrap validation n=200 confirmed the point estimates (Table 3 final-covariate-model column). ANC0 covariates identified in the univariate analysis (albumin, sex, blood transfusion, prior Pt-chemotherapy) were NOT retained in the final full model because ANC0 is known per subject at start of therapy and is of less prognostic importance than MTT / gamma / SLOPE. Ethnicity covariate on SLOPE (Asian / Japanese) was tested in the univariate analysis (Table 4) but did not reach the P<0.001 significance threshold and was not retained. The reason G-CSF-on-gamma dropped out of the final model was likely confounding with the retained G-CSF-on-SLOPE effect."
  )

  ini({
    # =========================================================================
    # PK layer -- Majid 2014 three-compartment model with linear elimination,
    # reproduced verbatim from Kawamura 2018 section 2.3. van Hasselt 2013 cites
    # the same Eisai popPK analysis as "unpublished data on file, Dr Z. Hussein
    # (Eisai)" and does not report PK parameter values; the values below are
    # FIXED from the Majid 2014 publication of the same Eisai analysis (see
    # `reference` above and vignette Errata). Reference values: WT 68.7 kg,
    # ALB 4.0 g/dL, ALP 132 U/L, TBILI 0.5 mg/dL.
    # =========================================================================
    lcl  <- fixed(log(3.11)); label("CL: eribulin clearance (L/h)")                                  # Majid 2014 (via Kawamura 2018 section 2.3): CL[L/h] = 3.11 * (WT/68.7)^0.75 * (ALB/4.0)^0.946 * (ALP/132)^-0.209 * (BILI/0.5)^-0.180
    lvc  <- fixed(log(4.06)); label("V1: central volume of distribution (L)")                         # Majid 2014 (via Kawamura 2018 section 2.3): V1[L] = 4.06 * (WT/68.7)
    lq   <- fixed(log(2.64)); label("Q2: intercompartmental clearance central <-> peripheral1 (L/h)") # Majid 2014 (via Kawamura 2018 section 2.3): Q2[L/h] = 2.64 * (WT/68.7)^0.75
    lvp  <- fixed(log(2.42)); label("V2: peripheral1 volume of distribution (L)")                     # Majid 2014 (via Kawamura 2018 section 2.3): V2[L] = 2.42 * (WT/68.7)
    lq2  <- fixed(log(6.60)); label("Q3: intercompartmental clearance central <-> peripheral2 (L/h)") # Majid 2014 (via Kawamura 2018 section 2.3): Q3[L/h] = 6.60 * (WT/68.7)^0.75
    lvp2 <- fixed(log(121));  label("V3: peripheral2 volume of distribution (L)")                     # Majid 2014 (via Kawamura 2018 section 2.3): V3[L] = 121 * (WT/68.7)

    # PK covariate exponents -- all FIXED from Majid 2014.
    e_wt_cl    <- fixed(0.75);   label("Allometric exponent of WT on CL")
    e_alb_cl   <- fixed(0.946);  label("Power exponent of ALB (g/dL scale) on CL")
    e_alp_cl   <- fixed(-0.209); label("Power exponent of ALP on CL")
    e_tbili_cl <- fixed(-0.180); label("Power exponent of TBILI (mg/dL scale) on CL")
    e_wt_q     <- fixed(0.75);   label("Allometric exponent of WT on Q2")
    e_wt_q2    <- fixed(0.75);   label("Allometric exponent of WT on Q3")
    e_wt_vc    <- fixed(1);      label("Linear exponent of WT on V1")
    e_wt_vp    <- fixed(1);      label("Linear exponent of WT on V2")
    e_wt_vp2   <- fixed(1);      label("Linear exponent of WT on V3")

    # =========================================================================
    # PD layer -- van Hasselt 2013 semi-physiological Friberg-style myelosup-
    # pression model. Final covariate model parameter estimates from Table 3
    # ("Full covariate model" block); IIV converted from CV% to log-normal
    # variance via omega^2 = log(1 + CV^2).
    # =========================================================================
    lcirc0 <- log(4030); label("ANC0: baseline circulating neutrophil count (cells/uL)")            # van Hasselt 2013 Table 3 final: ANC0 = 4.03 * 10^9 cells/L = 4030 cells/uL (RSE 1.2%)
    lmtt   <- log(109);  label("MTT: mean transit time through the maturation chain (h)")           # van Hasselt 2013 Table 3 final: MTT = 109 h (RSE 1.8%)
    lgamma <- log(0.216);label("Gamma: feedback exponent on (ANC0/circ) (unitless)")                # van Hasselt 2013 Table 3 final: gamma = 0.216 (RSE 1.4%)
    lslope <- log(0.0451); label("Slope: linear drug-effect coefficient (L/ug, equivalently mL/ng)") # van Hasselt 2013 Table 3 final: SLOPE = 0.0451 L/ug (RSE 3.2%); paper unit "ug l^-1" is the inverse-concentration form L/ug (dimensionally identical to mL/ng)

    # -------------------------------------------------------------------------
    # Covariate effects on MTT -- van Hasselt 2013 Table 3 final covariate
    # model. Functional form per Eq. 10 of the paper:
    #   MTT_ind = MTT_typ * (ALB_gdL/4)^e_alb_mtt * (TBILI_mgdL/2)^e_tbili_mtt
    #           * (ALP/100)^e_alp_mtt * (LDH/238)^e_ldh_mtt
    #           * e_gcsf_mtt^CONMED_GCSF
    # Scaling values (Table 2 footnote): ALB 4 g/dL, TBILI 2 mg/dL, ALP 100
    # U/L, LDH 238 U/L. G-CSF is dichotomous per Eq. 10 second product term.
    # -------------------------------------------------------------------------
    e_alb_mtt   <-  0.374;  label("Power exponent of ALB on MTT (unitless)")                        # van Hasselt 2013 Table 3 final: 0.374 (RSE 23.6%)
    e_tbili_mtt <- -0.046;  label("Power exponent of TBILI on MTT (unitless)")                      # van Hasselt 2013 Table 3 final: -0.046 (RSE 25.7%)
    e_alp_mtt   <- -0.0337; label("Power exponent of ALP on MTT (unitless)")                        # van Hasselt 2013 Table 3 final: -0.0337 (RSE 30.3%)
    e_ldh_mtt   <- -0.0561; label("Power exponent of LDH on MTT (unitless)")                        # van Hasselt 2013 Table 3 final: -0.0561 (RSE 20.5%)
    e_gcsf_mtt  <-  0.883;  label("Multiplicative effect of CONMED_GCSF on MTT (theta^cov form)")   # van Hasselt 2013 Table 3 final: 0.883 (RSE 2.2%)

    # -------------------------------------------------------------------------
    # Covariate effects on SLOPE -- van Hasselt 2013 Table 3 final covariate
    # model. Functional form per Eq. 10 of the paper:
    #   SLOPE_ind = SLOPE_typ * (ALB_gdL/4)^e_alb_slope * (AST/30)^e_ast_slope
    #             * e_gcsf_slope^CONMED_GCSF
    # AST scaling defaults to the cohort median 30 U/L per Table 2 (Table 2
    # footnote enumerates non-median scaling values only for WT/ALB/ALP/BILI/
    # LDH; AST is not asterisked, so scaling = median).
    # -------------------------------------------------------------------------
    e_alb_slope  <- 0.763;  label("Power exponent of ALB on SLOPE (unitless)")                      # van Hasselt 2013 Table 3 final: 0.763 (RSE 18.6%)
    e_ast_slope  <- 0.119;  label("Power exponent of AST on SLOPE (unitless)")                      # van Hasselt 2013 Table 3 final: 0.119 (RSE 24.4%)
    e_gcsf_slope <- 1.3;    label("Multiplicative effect of CONMED_GCSF on SLOPE (theta^cov form)") # van Hasselt 2013 Table 3 final: 1.3 (RSE 8.2%)

    # -------------------------------------------------------------------------
    # Inter-individual variability -- van Hasselt 2013 Table 3 final covariate
    # model. Paper reports CV% for a log-normal IIV structure; converted to the
    # log-eta variance via omega^2 = log(1 + CV^2). Off-diagonal covariances
    # were not computationally feasible per the paper ("Estimation of inter-
    # occasion variability (IOV) and off-diagonal covariances in IIV were not
    # computationally feasible"), so all IIVs are diagonal.
    #   ANC0:  CV% = 37.3 -> omega^2 = log(1 + 0.373^2) = 0.13025
    #   MTT:   CV% = 13.9 -> omega^2 = log(1 + 0.139^2) = 0.01914
    #   Gamma: CV% = 12.2 -> omega^2 = log(1 + 0.122^2) = 0.01477
    #   SLOPE: CV% = 54.0 -> omega^2 = log(1 + 0.540^2) = 0.25610
    # -------------------------------------------------------------------------
    etalcirc0 ~ 0.13025  # van Hasselt 2013 Table 3 final: ANC0 CV% = 37.3
    etalmtt   ~ 0.01914  # van Hasselt 2013 Table 3 final: MTT CV% = 13.9
    etalgamma ~ 0.01477  # van Hasselt 2013 Table 3 final: gamma CV% = 12.2
    etalslope ~ 0.25610  # van Hasselt 2013 Table 3 final: SLOPE CV% = 54.0

    # -------------------------------------------------------------------------
    # Residual unexplained variability -- van Hasselt 2013 Methods "Statistical
    # model development" states residual error on ANC was included as a pro-
    # portional relationship. Table 3 final: proportional CV% = 49.6.
    # -------------------------------------------------------------------------
    propSd <- 0.496; label("Proportional residual SD on circulating ANC (fraction)")                # van Hasselt 2013 Table 3 final: proportional CV% = 49.6 (RSE 2.9%)
  })

  model({
    # ---- Unit conversions for covariate columns supplied in canonical SI units ----
    # The canonical register unit for ALB is g/L (SI) and TBILI is umol/L (SI);
    # the Majid 2014 PK coefficients and the van Hasselt 2013 PD coefficients
    # were calibrated against the US-convention units g/dL and mg/dL. Convert
    # inline so the g/L / umol/L canonical columns drive the model at the
    # published calibration.
    alb_gdL   <- ALB * 0.1     # canonical g/L -> paper g/dL: 1 g/dL = 10 g/L
    tbili_mgdL <- TBILI / 17.1 # canonical umol/L -> paper mg/dL: 1 mg/dL = 17.1 umol/L

    # ---- PK individual parameters (typical-value only; all theta FIXED from Majid 2014) ----
    cl   <- exp(lcl)  * (WT/68.7)^e_wt_cl * (alb_gdL/4.0)^e_alb_cl * (ALP/132)^e_alp_cl * (tbili_mgdL/0.5)^e_tbili_cl
    vc   <- exp(lvc)  * (WT/68.7)^e_wt_vc
    q    <- exp(lq)   * (WT/68.7)^e_wt_q
    vp   <- exp(lvp)  * (WT/68.7)^e_wt_vp
    q2   <- exp(lq2)  * (WT/68.7)^e_wt_q2
    vp2  <- exp(lvp2) * (WT/68.7)^e_wt_vp2

    kel <- cl  / vc
    k12 <- q   / vc
    k21 <- q   / vp
    k13 <- q2  / vc
    k31 <- q2  / vp2

    # ---- PD individual parameters with covariate effects ----
    # van Hasselt 2013 Table 3 final covariate model (Eq. 10 of the paper):
    #   MTT   = MTT_typ * (ALB_gdL/4)^e_alb_mtt * (TBILI_mgdL/2)^e_tbili_mtt
    #         * (ALP/100)^e_alp_mtt * (LDH/238)^e_ldh_mtt
    #         * e_gcsf_mtt^CONMED_GCSF
    #   SLOPE = SLOPE_typ * (ALB_gdL/4)^e_alb_slope * (AST/30)^e_ast_slope
    #         * e_gcsf_slope^CONMED_GCSF
    circ0 <- exp(lcirc0 + etalcirc0)
    mtt   <- exp(lmtt   + etalmtt) *
      (alb_gdL/4)^e_alb_mtt *
      (tbili_mgdL/2)^e_tbili_mtt *
      (ALP/100)^e_alp_mtt *
      (LDH/238)^e_ldh_mtt *
      e_gcsf_mtt^CONMED_GCSF
    gamma <- exp(lgamma + etalgamma)
    slope <- exp(lslope + etalslope) *
      (alb_gdL/4)^e_alb_slope *
      (AST/30)^e_ast_slope *
      e_gcsf_slope^CONMED_GCSF
    ktr   <- 4 / mtt

    # ---- Eribulin plasma concentration ----
    # Dose is supplied in mg of eribulin free base; central holds the amount
    # in mg; vc is in L; so Cc = central/vc is in mg/L. 1 mg/L = 1000 ug/L,
    # so the factor of 1000 in edrug below converts Cc to ug/L (the concen-
    # tration unit implicit in the paper's SLOPE = 0.0451 L/ug).
    Cc <- central / vc

    # ---- Drug effect (linear in concentration) ----
    # van Hasselt 2013 Eq. 6: E = SLOPE * C. Paper's SLOPE has units L/ug
    # (equivalent to mL/ng). Convert Cc from mg/L to ug/L via the factor 1000.
    edrug <- slope * Cc * 1000

    # ---- Three-compartment IV PK (Majid 2014 structure) ----
    d/dt(central)     <- -(kel + k12 + k13) * central + k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    # ---- Friberg myelosuppression chain (van Hasselt 2013 Eqs. 1-6) ----
    # Compartment naming follows the Friberg-family precedent in
    # nlmixr2lib (Friberg_2002_paclitaxel, Kawamura_2018_eribulin):
    #   precursor1 = Prol (proliferating progenitor pool)
    #   precursor2..4 = Transit1..3 (maturation chain)
    #   circ = Circ (circulating neutrophils, observed)
    # Feedback term uses the individual estimated baseline ANC0 (circ0) in
    # the ratio (circ0/circ)^gamma per Eq. 3 of the paper. Under the canonical
    # Friberg parameterisation kprol = kcirc = ktr = 4/MTT, so a single ktr
    # rate constant drives the whole chain including the circulating loss.
    d/dt(precursor1) <- ktr * precursor1 * (1 - edrug) * (circ0 / circ)^gamma - ktr * precursor1
    d/dt(precursor2) <- ktr * precursor1 - ktr * precursor2
    d/dt(precursor3) <- ktr * precursor2 - ktr * precursor3
    d/dt(precursor4) <- ktr * precursor3 - ktr * precursor4
    d/dt(circ)       <- ktr * precursor4 - ktr * circ

    # ---- Initial conditions ----
    # All five Friberg compartments start at the individual estimated baseline
    # ANC0. Because kprol = kcirc = ktr exactly (single ktr rate), the system
    # is at exact steady state at t=0 with no drug: d/dt(circ)(0) = ktr*ANC0 -
    # ktr*ANC0 = 0, d/dt(precursor1)(0) = ktr*ANC0*1*(ANC0/ANC0)^gamma - ktr*
    # ANC0 = 0, etc. No baseline drift.
    precursor1(0) <- circ0
    precursor2(0) <- circ0
    precursor3(0) <- circ0
    precursor4(0) <- circ0
    circ(0)       <- circ0

    # ---- Observations ----
    # ANC is the observed circulating neutrophil count (cells/uL). Proportional
    # residual error per van Hasselt 2013 Table 3 (proportional CV% = 49.6).
    ANC <- circ
    ANC ~ prop(propSd)
  })
}
