Wang_2021_pertuzumab <- function() {
  description <- "Two-compartment population PK model with first-order subcutaneous absorption and bioavailability for pertuzumab (Perjeta) administered either intravenously or as the fixed-dose combination subcutaneous formulation with trastuzumab (PH FDC SC) in patients with HER2-positive early breast cancer in the FeDeriCa study (Wang 2021)"
  reference <- "Wang B, Deng R, Hennig S, Badovinac Crnjevic T, Kaewphluk M, Kagedal M, Quartino AL, Girish S, Li C, Kirschbrown WP. Population pharmacokinetic and exploratory exposure-response analysis of the fixed-dose combination of pertuzumab and trastuzumab for subcutaneous injection in patients with HER2-positive early breast cancer in the FeDeriCa study. Cancer Chemother Pharmacol. 2021;88(3):439-451. doi:10.1007/s00280-021-04296-0"
  vignette <- "Wang_2021_pertuzumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    LBM = list(
      description        = "Baseline lean body weight (LBW)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. Power effect on CL (exponent 1.252), Vc (0.839), and Vp (0.716); reference 45.09 kg per Wang 2021 Results CL/Vc/Vp covariate equations (median LBW of the 489-patient analysis dataset). The paper reports the covariate as 'lean body weight (LBW)'; the canonical nlmixr2lib column is LBM (lean body mass).",
      source_name        = "LBW"
    ),
    ALB = list(
      description        = "Baseline serum albumin concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value (SI convention, g/L). Power effect on CL (exponent -0.629); reference 43.25 g/L per Wang 2021 CL covariate equation (median albumin of the analysis dataset).",
      source_name        = "ALB"
    ),
    RACE_ASIAN = list(
      description        = "Asian region indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Asian region)",
      notes              = "Per-subject indicator for enrollment in an Asian region (1 = Asian region, 0 = non-Asian). Multiplicative effect on CL: typical CL is 1.123-fold higher (12.3%) in Asian patients per Wang 2021 CL covariate equation (theta13 = 0.123). Wang 2021 frames this as a region-of-enrollment covariate ('Asian region') that pools Asian-region cohorts (n=100 in the analysis dataset) against the non-Asian-region reference; the canonical nlmixr2lib column RACE_ASIAN is used here because the paper's region effect is operationally a race / ethnicity proxy (the paper notes the analogous Garg 2014 IV pertuzumab model had a similar effect driven by Japanese subjects).",
      source_name        = "Asian region"
    ),
    ROUTE_IV = list(
      description        = "Indicator for intravenous (IV) administration of pertuzumab",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (subcutaneous; PH FDC SC arm)",
      notes              = "Per-subject dosing-route indicator (1 = P+H IV cohort, 0 = PH FDC SC cohort) used to switch the proportional residual SD between the IV-arm value (`CcpropSdIv` = 0.175) and the SC-arm value (`CcpropSdSc` = 0.155) per Wang 2021 Table 1 theta7 / theta8. Distinct from the rxode2 `cmt` event column (cmt = central for IV doses, cmt = depot for SC doses); ROUTE_IV is the per-subject covariate that picks the residual-error term.",
      source_name        = "IV"
    )
  )

  population <- list(
    n_subjects       = 489L,
    n_observations   = 5180L,
    n_studies        = 1L,
    study            = "FeDeriCa (NCT03493854); randomized, open-label, international, multicenter, phase III non-inferiority trial of PH FDC SC vs P + H IV in the neoadjuvant-adjuvant early breast cancer setting (500 patients randomized; 489 contributed pertuzumab PK to the popPK analysis).",
    age_median       = "Median age not reproduced in the main paper text (Online Resource 1 supplement is not on disk); 43 of 489 patients (8.8%) were >65 years old per Wang 2021 Results.",
    weight_range     = "Lean body weight 5th-95th percentiles 38-53 kg (median 45.09 kg); total body weight not separately tabulated in the main paper text.",
    sex_female_pct   = 100,
    race_ethnicity   = c(Asian = 20.4, `non-Asian` = 79.6),
    disease_state    = "HER2-positive (immunohistochemistry 3+ or in situ hybridization-positive) operable, locally advanced, or inflammatory stage II-IIIC early breast cancer with a primary tumor >2 cm in diameter, or node-positive disease; ECOG performance status 0 or 1; LVEF >=55%.",
    dose_range       = "P + H IV arm (n=246, 50.3%): pertuzumab 840 mg loading dose IV, then 420 mg IV every 3 weeks (q3w) co-administered with trastuzumab 8 mg/kg loading then 6 mg/kg q3w. PH FDC SC arm (n=243, 49.7%): pertuzumab 1200 mg + trastuzumab 600 mg fixed-dose combination SC loading dose in 15 mL, then pertuzumab 600 mg + trastuzumab 600 mg SC maintenance in 10 mL q3w (with 2000 U/mL recombinant human hyaluronidase to aid SC dispersion).",
    regions          = "International, 19 countries / 106 centers; ~20% Asian-region enrollment.",
    sampling         = "Sparse sampling cycles 5 through day 1 of cycle 8. P + H IV arm: pre- and post-dose day 1 of cycles 5, 6, 7, 8, day 15 of cycle 5. PH FDC SC arm: pre-dose day 1 of cycles 5, 6, 7, 8; days 2 and 15 of cycle 5; days 2, 4, 8, and 15 of cycle 7. 5180 evaluable pertuzumab samples (2093 SC + 3087 IV), assayed by validated duplex hybrid immunoaffinity capture LC-MS/MS with LLOQ 100 ng/mL.",
    reference_subject = "Median LBW = 45.09 kg, median albumin = 43.25 g/L, non-Asian region; the typical patient against which Wang 2021 forest plots and exposure ratios are computed.",
    notes            = "Demographics from Wang 2021 Results 'Patients and samples' and Online Resource 1 (not on disk in this worktree); Online Resource 1 covariate-distribution figures (LBW Q1 <42.0, Q4 >48.8 kg) are summarized in Wang 2021 Figs 1-2. The dataset pools P + H IV (intravenous pertuzumab + trastuzumab) and PH FDC SC (subcutaneous fixed-dose combination of pertuzumab + trastuzumab + recombinant human hyaluronidase) cohorts; the same model is fit jointly with route-specific proportional residual error. ER analyses (logistic regression of tpCR / safety endpoints versus model-predicted exposure) are not reproduced in the model file; the popPK structural model and parameter estimates are."
  )

  ini({
    # Structural typical-value parameters at the reference subject:
    # LBW = 45.09 kg, albumin = 43.25 g/L, non-Asian region (Wang 2021 Table 1).
    # Time = day, dose = mg, central volume in L -> Cc = central / vc has units
    # mg/L = ug/mL, which matches the published Ctrough / Cmax units.
    lcl      <- log(0.163);  label("Linear CL at reference covariates (L/day)")                 # Wang 2021 Table 1, theta1
    lvc      <- log(2.77);   label("Central volume of distribution Vc at reference (L)")        # Wang 2021 Table 1, theta2
    lq       <- log(0.616);  label("Inter-compartmental clearance Q (L/day)")                   # Wang 2021 Table 1, theta3
    lvp      <- log(2.49);   label("Peripheral volume of distribution Vp at reference (L)")     # Wang 2021 Table 1, theta4
    lka      <- log(0.348);  label("First-order SC absorption rate ka (1/day)")                 # Wang 2021 Table 1, theta5
    lfdepot  <- log(0.712);  label("Bioavailability of the SC depot (fraction)")                # Wang 2021 Table 1, theta6

    # Covariate effects on linear CL, Vc, Vp.
    # CL: power on LBW (exp 1.252; ref 45.09 kg), power on ALB (exp -0.629; ref
    # 43.25 g/L), and a multiplicative additive Asian-region effect (CL_Asian =
    # CL * (1 + e_asian_cl); the paper writes 'x 1.123 if Asian region' which
    # corresponds to e_asian_cl = 0.123, theta13).
    # Vc and Vp: power on LBW only.
    e_lbw_cl    <-  1.252; label("Power exponent of LBW on linear CL (unitless)")                                # Wang 2021 Table 1, theta10
    e_alb_cl    <- -0.629; label("Power exponent of albumin on linear CL (unitless)")                            # Wang 2021 Table 1, theta9
    e_asian_cl  <-  0.123; label("Multiplicative Asian-region effect on linear CL (unitless; CL_Asian = CL * (1 + e_asian_cl))") # Wang 2021 Table 1, theta13
    e_lbw_vc    <-  0.839; label("Power exponent of LBW on Vc (unitless)")                                       # Wang 2021 Table 1, theta11
    e_lbw_vp    <-  0.716; label("Power exponent of LBW on Vp (unitless)")                                       # Wang 2021 Table 1, theta12

    # Inter-individual variability. Wang 2021 Table 1 reports omega as %CV on
    # log-normal parameters; convert via omega^2 = log(CV^2 + 1). Etas are
    # reported as marginal CVs (no off-diagonal correlation block in Table 1),
    # so a diagonal Omega is used. F is parameterized as
    # F = exp(lfdepot + etalfdepot) per Wang 2021 ('F = 0.712 * exp(eta_F)');
    # individual F can exceed 1 when etalfdepot > log(1/0.712) ~= 0.34 sd, but
    # this is faithful to the published parameterization. Shrinkage on F (50.2%)
    # and Vp (49.7%) is high per Wang 2021 Table 1.
    etalcl     ~ 0.053723   # CL CV% 23.5% -> log(1 + 0.235^2); Wang 2021 Table 1
    etalvc     ~ 0.114349   # Vc CV% 34.8% -> log(1 + 0.348^2); Wang 2021 Table 1
    etalvp     ~ 0.063491   # Vp CV% 25.6% -> log(1 + 0.256^2); Wang 2021 Table 1
    etalfdepot ~ 0.031194   # F  CV% 17.8% -> log(1 + 0.178^2); Wang 2021 Table 1

    # Proportional residual error. Wang 2021 Table 1 reports separate SDs for
    # the SC cohort (theta7 = 0.155) and the IV cohort (theta8 = 0.175);
    # selected per-record by the ROUTE_IV indicator inside model().
    CcpropSdSc <- 0.155;  label("Proportional residual error, SC cohort (fraction)")  # Wang 2021 Table 1, theta7
    CcpropSdIv <- 0.175;  label("Proportional residual error, IV cohort (fraction)")  # Wang 2021 Table 1, theta8
  })

  model({
    # Individual structural parameters with covariate effects (Wang 2021
    # Results: CL / Vc / Vp / Q / ka / F covariate equations).
    cl <- exp(lcl + etalcl) *
            (LBM / 45.09)^e_lbw_cl *
            (ALB / 43.25)^e_alb_cl *
            (1 + e_asian_cl * RACE_ASIAN)
    vc <- exp(lvc + etalvc) * (LBM / 45.09)^e_lbw_vc
    vp <- exp(lvp + etalvp) * (LBM / 45.09)^e_lbw_vp
    q  <- exp(lq)
    ka <- exp(lka)

    # Bioavailability of the SC depot (applies to depot doses only; IV doses
    # bypass the depot via cmt = central in the data).
    fdepot <- exp(lfdepot + etalfdepot)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    Cc <- central / vc

    # Two-compartment PK with first-order SC absorption from depot and
    # IV / SC dosing supported simultaneously: SC doses go to depot (with F),
    # IV doses go directly to central.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1
    f(depot)          <- fdepot

    # Route-specific proportional residual error: pick CcpropSdSc when
    # ROUTE_IV = 0 (SC cohort) and CcpropSdIv when ROUTE_IV = 1 (IV cohort).
    CcpropSd <- CcpropSdSc + (CcpropSdIv - CcpropSdSc) * ROUTE_IV
    Cc ~ prop(CcpropSd)
  })
}
