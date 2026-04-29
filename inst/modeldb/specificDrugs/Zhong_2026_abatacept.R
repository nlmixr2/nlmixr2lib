Zhong_2026_abatacept <- function() {
  description <- "Two-compartment population PK model for abatacept (CTLA4-Ig Fc-fusion) pooled across 9 phase 2/3 studies (Zhong 2026): adults with rheumatoid arthritis, patients aged 2-17 years with polyarticular juvenile idiopathic arthritis, and patients aged 6+ years with hematologic malignancies receiving HLA-matched unrelated-donor HSCT (the ABA2 trial). Final model has zero-order IV infusion, first-order SC absorption, first-order linear elimination, additive plus proportional residual error, allometric weight on CL/VC/VP, hepatic (AST) and renal (cGFR) markers on CL, sex on CL and VC, two HSCT cohort indicators (7-of-8 and 8-of-8 HLA-matched URD) on CL/VC, and a logit-scale SC bioavailability sub-model with weight, age, and pJIA-disease covariates fixed to a previously developed internal JIA PPK model (values match Gandhi 2021)."
  reference <- "Zhong R, Maxwell K, Passarell J, Murthy B, Aras U, Williams D. Model-Informed Abatacept Dose Recommendation in Pediatric Patients With Acute Graft-Versus-Host Disease. J Clin Pharmacol. 2026;66(2):[in-issue]. doi:10.1002/jcph.70156"
  vignette <- "Zhong_2026_abatacept"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 67.9 kg (Zhong 2026 Figure 1 caption: reference patient weighs 67.9 kg). Power effect on CL (exp 0.876), VC (exp 0.712), VP (exp 0.839), and on logit-F (slope -0.506; F1 piece transferred from previous internal JIA PPK model whose Gandhi 2021 published reference is 68 kg, ~0.1 kg discrepancy with Zhong 2026 structural reference flagged in vignette Errata).",
      source_name        = "BWT"
    ),
    AGE = list(
      description        = "Baseline age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 49 years for F1 covariate effect (Zhong 2026 does not state the AGE reference; the F1 sub-model and its parameters are fixed to a previous internal JIA PPK model whose Gandhi 2021 published reference is 49 years, used here for consistency). Not retained on CL, VC, or VP in the final model. Slope effect on logit-F (slope 0.487).",
      source_name        = "AGE"
    ),
    AST = list(
      description        = "Baseline aspartate aminotransferase",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 20 U/L (Zhong 2026 Figure 1 caption: 'baseline AST level of 20 U/L'). Power effect on CL (exp -0.115); not clinically relevant per Zhong 2026 Discussion ('baseline AST and cGFR did not influence abatacept CL').",
      source_name        = "AST"
    ),
    CRCL = list(
      description        = "Calculated glomerular filtration rate (BSA-normalized)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 103 mL/min/1.73 m^2 (Zhong 2026 Figure 1 caption: 'cGFR of 103 mL/min/1.73 m^2'). Zhong 2026 uses 'cGFR' (calculated GFR, BSA-normalized via the Schwartz equation per the supplement); mapped to the canonical CRCL (which accepts either MDRD/CKD-EPI eGFR or BSA-normalized measured CrCl). Power effect on CL (exp 0.279); not clinically relevant per Zhong 2026 Discussion.",
      source_name        = "cGFR"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male; reference subject is male per Zhong 2026 Figure 1 caption 'The reference individual is male')",
      notes              = "Zhong 2026 supplement Methods: 'the reference for the categorical covariates was the mode in the PPK dataset, except for sex where male was used as the reference'. SEXF = 1 for female, 0 for male; the Table 2 'Exponent of sex effect in female' coefficients (-0.0572 on CL, -0.0967 on VC) are applied as `exp(SEXF * coef)` so females have ~5.6% lower CL and ~9.2% lower VC than males. Sign matches Li 2019 (sex on CL: -0.0722; same direction) but reversed from Gandhi 2021 (sex on CL: +0.0674; sign-flipped because the Gandhi 2021 dataset added pJIA data which shifted the female-vs-male CL contrast). The female-on-CL effect is not clinically relevant per Zhong 2026 (within 80%-125% reference range).",
      source_name        = "SEX"
    ),
    DIS_PJIA = list(
      description        = "Polyarticular juvenile idiopathic arthritis disease-state indicator, 1 = pJIA, 0 = adult RA or HM (or other non-pJIA)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-pJIA; the reference complement here is the union of adult RA and HM patients pooled with pJIA in Zhong 2026)",
      notes              = "Zhong 2026 codes a JIA-vs-non-JIA disease indicator (JIA = 1 if pJIA, 0 otherwise). Both JIA studies in the analysis (IM101033, IM101301) are polyarticular JIA per Supplementary Table S1, so JIA = pJIA in this dataset. The covariate enters the bioavailability model only (additive on the logit scale: `logit_F = logit_F_TV + 3.08 * DIS_PJIA + ...`); pJIA patients have substantially higher SC bioavailability than the non-pJIA reference (logit_F shifts from 1.21 at the reference to ~4.29 in pJIA, i.e. F_abs ~ 0.77 -> ~0.987). The F1 sub-model and all F1-covariate effects are fixed in Zhong 2026 to the final estimates from a previously developed internal abatacept JIA PPK model (Methods text). Values are identical to Gandhi 2021 Table 2 for F1, weight-on-F1, age-on-F1, JIA-on-F1, and IIV on F1, suggesting the internal model is the published Gandhi 2021 model or a closely related variant.",
      source_name        = "JIA"
    ),
    HSCT_URD_7OF8 = list(
      description        = "Hematopoietic stem cell transplant from a 7-of-8 HLA-matched unrelated donor (single-allele mismatch); 1 = yes, 0 = otherwise",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not in the 7-of-8-matched HSCT cohort; the reference complement is the union of RA, pJIA, and 8-of-8-matched HSCT patients)",
      notes              = "Zhong 2026 Table 2 / Figure 1: 'Cohort 7/8' from the ABA2 trial (Study IM101311), encoding patients receiving HSCT from a single-allele-mismatched URD. Exponential coefficient -0.326 on CL: HSCT_URD_7OF8 = 1 patients have ~28% lower CL than the reference complement (clinically relevant per Zhong 2026 Discussion). No effect on VC, VP, or KA. The 7-of-8 cohort represents a higher GvHD-risk population because of the single-allele HLA mismatch; the underlying disease biology of the CL effect is not characterized.",
      source_name        = "COHORT7"
    ),
    HSCT_URD_8OF8 = list(
      description        = "Hematopoietic stem cell transplant from an 8-of-8 HLA-matched unrelated donor (full match); 1 = yes, 0 = otherwise",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not in the 8-of-8-matched HSCT cohort; the reference complement is the union of RA, pJIA, and 7-of-8-matched HSCT patients)",
      notes              = "Zhong 2026 Table 2 / Figure 1: 'Cohort 8/8' from the ABA2 trial (Study IM101311), encoding patients receiving HSCT from a fully-HLA-matched URD. Exponential coefficient -0.0934 on CL and +0.257 on VC: HSCT_URD_8OF8 = 1 patients have ~9% lower CL and ~29% higher VC than the reference complement. The CL effect is not clinically relevant; the VC effect is potentially clinically relevant per Zhong 2026 Discussion. No effect on VP or KA.",
      source_name        = "COHORT8"
    )
  )

  population <- list(
    n_subjects     = 904L,
    n_observations = 6355L,
    n_studies      = 9L,
    age_range      = "2-84 years (overall pooled cohort: minimum 2 years, maximum 84 years per Zhong 2026 Table 1)",
    weight_range   = "12-187 kg (overall pooled cohort per Zhong 2026 Table 1)",
    sex_female_pct = 68.6,
    race_ethnicity = "White 83.5%, Black 5.53%, Hispanic/Latin 4.09%, Asian/Pacific Islander 1.88%, Other 4.98% (Zhong 2026 Table 1).",
    disease_state  = "Pooled adult rheumatoid arthritis (n = 386, 42.7%), polyarticular juvenile idiopathic arthritis (pJIA; n = 403, 44.6%, ages 2-17 years), and hematologic malignancies (HM; n = 115, 12.7%, ages 6-76 years receiving HLA-matched URD HSCT in the ABA2 trial).",
    dose_range     = "IV abatacept 0.5-10 mg/kg Q4W (RA + pJIA studies including IM101033) and weight-tiered 500-1000 mg flat-dose IV (RA studies IM101029, IM101031, IM101102); IV 10 mg/kg on Days -1, 5, 14, 28 (HM ABA2 / IM101311); SC abatacept weight-tiered 50/87.5/125 mg QW (pJIA IM101301). Pediatric aGvHD recommended regimen (this paper's deliverable): 15 mg/kg loading dose on Day -1 followed by 12 mg/kg on Days 5, 14, 28 for patients aged 2 to <6 years.",
    regions        = "Multi-regional (9 pooled phase 2/3 studies).",
    reference_subject = "Male, 67.9 kg, baseline AST 20 U/L, calculated GFR 103 mL/min/1.73 m^2, not in HSCT cohort 7/8 or 8/8 (Zhong 2026 Figure 1 caption). The F1 sub-model uses an additional reference of 49-year-old, non-pJIA (DIS_PJIA = 0), inherited from the previously developed internal JIA PPK model whose published equivalent Gandhi 2021 used WT_ref = 68 kg and AGE_ref = 49 years. The 0.1 kg discrepancy between Zhong 2026 structural-reference WT (67.9 kg) and Gandhi 2021 F1-reference WT (68 kg) is documented in the vignette Errata; the impact on F1 predictions is < 0.1% across the model's covariate range.",
    notes          = "Below-LLOQ samples accounted for 2.9% of all collected samples and were excluded from the analysis. The HM cohort is decomposed into HSCT_URD_7OF8 (n = 41) and HSCT_URD_8OF8 (n = 74) per Figure 1 stratification, encoded as two orthogonal binary indicators with the RA + pJIA pool as the reference complement. Zhong 2026 was the first comprehensive PPK characterization of both IV and SC abatacept in pooled RA + pJIA + HSCT-aGvHD data and supported the FDA approval of abatacept for prevention of acute GvHD in pediatric patients aged 2 to < 6 years (15 mg/kg loading dose Day -1 followed by 12 mg/kg on Days 5, 14, 28)."
  )

  ini({
    # Structural PK parameters - Zhong 2026 Table 2 final-model estimates.
    # Paper reports CL, Q in L/h and KA in 1/h; values below convert to per-day
    # (x 24) to match the nlmixr2lib convention of time in days. KA is
    # parameterised with the Gandhi 2021 / Li 2019 flip-flop-prevention scheme
    # (ka = exp(lka + etaKA) + kel) so that exp(lka) is the "additional"
    # absorption rate above kel, matching the BMS author-group convention used
    # for both prior abatacept popPK models (Li 2019, Gandhi 2021); see
    # vignette Errata for the rationale.
    lcl         <- log(0.0230 * 24);    label("Clearance CL (L/day) at reference covariates")                  # Zhong 2026 Table 2: CL_TV = 0.0230 L/h (%RSE 2.32)
    lvc         <- log(3.19);           label("Central volume VC (L) at reference covariates")                 # Zhong 2026 Table 2: VC_TV = 3.19 L (%RSE 2.48)
    lq          <- log(0.0303 * 24);    label("Inter-compartmental clearance Q (L/day)")                       # Zhong 2026 Table 2: Q_TV = 0.0303 L/h (%RSE 5.23); IIV not estimated (NE)
    lvp         <- log(5.07);           label("Peripheral volume VP (L) at reference covariates")              # Zhong 2026 Table 2: VP_TV = 5.07 L (%RSE 3.11)
    lka         <- log(0.00705 * 24);   label("Typical relative absorption rate KA_TV (1/day)")                 # Zhong 2026 Table 2: KA_TV = 0.00705 1/h (%RSE 19.2)
    logitfdepot <- 1.21;                label("Logit-scale typical SC bioavailability F_TV,ref (unitless)")     # Zhong 2026 Table 2: F_TV = 1.21 FIXED (normal-scale F = 0.77); transferred from previous internal JIA PPK model (matches Gandhi 2021)

    # Covariate effects - Zhong 2026 Table 2. Reference subject per Figure 1
    # caption: male, 67.9 kg, AST 20 U/L, cGFR 103 mL/min/1.73 m^2, not in
    # HSCT cohort 7/8 or 8/8. F1 covariates use AGE_ref = 49 yr (inherited
    # from the internal JIA PPK model; not stated explicitly in Zhong 2026
    # but the F1 piece is fixed to that internal model and Gandhi 2021's
    # published reference is 49 years).
    e_wt_cl    <-  0.876;  label("Power exponent of (WT/67.9 kg) on CL (unitless)")                              # Zhong 2026 Table 2: 'Power of weight on CL' = 0.876 (%RSE 3.08)
    e_ast_cl   <- -0.115;  label("Power exponent of (AST/20 U/L) on CL (unitless)")                              # Zhong 2026 Table 2: 'Power of AST effect' = -0.115 (%RSE 24.2)
    e_crcl_cl  <-  0.279;  label("Power exponent of (CRCL/103 mL/min/1.73 m^2) on CL (unitless)")                 # Zhong 2026 Table 2: 'Power of cGFR effect' = 0.279 (%RSE 9.12)
    e_sexf_cl  <- -0.0572; label("Exponential coefficient on CL for female sex (SEXF=1; unitless)")               # Zhong 2026 Table 2: 'Exponent of sex effect in female' on CL = -0.0572 (%RSE 41.2)
    e_co7_cl   <- -0.326;  label("Exponential coefficient on CL for HSCT_URD_7OF8 cohort (unitless)")             # Zhong 2026 Table 2: 'Exponent of GVHD cohort 7/8 effect' on CL = -0.326 (%RSE 16.8)
    e_co8_cl   <- -0.0934; label("Exponential coefficient on CL for HSCT_URD_8OF8 cohort (unitless)")             # Zhong 2026 Table 2: 'Exponent of GVHD cohort 8/8 effect' on CL = -0.0934 (%RSE 36.8)
    e_wt_vc    <-  0.712;  label("Power exponent of (WT/67.9 kg) on VC (unitless)")                              # Zhong 2026 Table 2: 'Power of weight on VC' = 0.712 (%RSE 4.30)
    e_sexf_vc  <- -0.0967; label("Exponential coefficient on VC for female sex (SEXF=1; unitless)")               # Zhong 2026 Table 2: 'Exponent of sex effect in female' on VC = -0.0967 (%RSE 28.6)
    e_co8_vc   <-  0.257;  label("Exponential coefficient on VC for HSCT_URD_8OF8 cohort (unitless)")             # Zhong 2026 Table 2: 'Exponent of GVHD cohort 8/8 effect' on VC = +0.257 (%RSE 14.5; paper labels this 'cohort 8' but Supplementary Table S2 confirms 8/8)
    e_wt_vp    <-  0.839;  label("Power exponent of (WT/67.9 kg) on VP (unitless)")                              # Zhong 2026 Table 2: 'Power of weight on VP' = 0.839 (%RSE 5.24)
    e_wt_f     <- -0.506;  label("Slope of log(WT/67.9 kg) on logit-F (unitless)")                                # Zhong 2026 Table 2: 'Power of weight on F1' = -0.506 FIXED; transferred from internal JIA PPK model (matches Gandhi 2021 'Power of body weight on bioavailability' = -0.506)
    e_age_f    <-  0.487;  label("Slope of log(AGE/49 yr) on logit-F (unitless)")                                 # Zhong 2026 Table 2: 'Power of age on F1' = 0.487 FIXED; transferred from internal JIA PPK model (matches Gandhi 2021 'Power of age on bioavailability' = 0.487)
    e_jia_f    <-  3.08;   label("Additive coefficient on logit-F for pJIA (DIS_PJIA=1; unitless)")               # Zhong 2026 Table 2: 'Exponent of JIA on F1' = 3.08 FIXED; transferred from internal JIA PPK model (matches Gandhi 2021 'Exponent of JIA on bioavailability' = 3.08)

    # Inter-individual variability - Zhong 2026 Supplementary Table S2
    # 'Final Model With Study IM101301' column reports IIV variances
    # directly. Q has no IIV (Table 2 'NE'). All non-F1 IIVs are diagonal
    # (paper reports separate per-parameter variances and shrinkage with
    # no full-block notation). IIV on F1 is fixed to the internal JIA PPK
    # model's value (variance 0.515524 = (0.718 SD FIXED)^2).
    etalcl         ~ 0.0689 ; label("IIV variance on log-CL (Zhong 2026 Suppl. Table S2; 26.7% CV; shrinkage 7.7%)")
    etalvc         ~ 0.0309 ; label("IIV variance on log-VC (Zhong 2026 Suppl. Table S2; 17.7% CV; shrinkage 51.5%)")
    etalvp         ~ 0.164  ; label("IIV variance on log-VP (Zhong 2026 Suppl. Table S2; 42.2% CV; shrinkage 32.8%)")
    etalka         ~ 0.861  ; label("IIV variance on the relative log-KA (Zhong 2026 Suppl. Table S2; 117% CV; shrinkage 59.8%)")
    etalogitfdepot ~ 0.516  ; label("IIV variance on logit-F (Zhong 2026 Suppl. Table S2 'IIV in F1' = 0.516; equivalent to Table 2 SD 0.718 FIXED, variance 0.718^2 = 0.515524; shrinkage 100%)")

    # Residual error - Zhong 2026 Table 2 footnote a reports an
    # additive + proportional model with proportional VARIANCE 0.0724 and
    # additive VARIANCE 5.26E-04 (mg^2/L^2). Convert to SDs for nlmixr2's
    # add() / prop() conventions.
    propSd <- sqrt(0.0724); label("Proportional residual error (fraction)")                # Zhong 2026 Table 2: SIGMA_PROP variance = 0.0724 (%RSE 4.77); SD = 0.269
    addSd  <- sqrt(5.26e-04); label("Additive residual error (mg/L = ug/mL)")               # Zhong 2026 Table 2: SIGMA_ADD variance = 5.26E-04 (%RSE 66.9); SD = 0.0229 mg/L
  })

  model({
    # Individual PK parameters with Zhong 2026 final-model covariate equations
    # (Table 2 / page 6). Reference subject: male, 67.9 kg, AST 20 U/L,
    # cGFR 103 mL/min/1.73 m^2, not in HSCT cohort 7/8 or 8/8. F1 piece uses
    # AGE_ref = 49 yr (inherited from internal JIA PPK model).
    cl <- exp(lcl + etalcl) *
          (WT / 67.9)^e_wt_cl *
          (AST / 20)^e_ast_cl *
          (CRCL / 103)^e_crcl_cl *
          exp(SEXF * e_sexf_cl + HSCT_URD_7OF8 * e_co7_cl + HSCT_URD_8OF8 * e_co8_cl)
    vc <- exp(lvc + etalvc) *
          (WT / 67.9)^e_wt_vc *
          exp(SEXF * e_sexf_vc + HSCT_URD_8OF8 * e_co8_vc)
    q  <- exp(lq)
    vp <- exp(lvp + etalvp) * (WT / 67.9)^e_wt_vp

    # KA is parameterised per the BMS abatacept-author-group convention
    # (Li 2019, Gandhi 2021) to ensure KA > k_el (flip-flop prevention):
    # KA_i = KA_TV * exp(etaKA) + CL_i / VC_i. Zhong 2026 reports KA_TV
    # alongside the prior models' values (Discussion paragraph: "the
    # inclusion of data from patients with GVHD resulted in a 35.3%
    # increase in the KA (0.00521 to 0.00705 1/h)") which strongly
    # implies the same parameterisation.
    kel <- cl / vc
    ka  <- exp(lka + etalka) + kel

    # Logit-scale SC bioavailability with independent IIV. pJIA disease,
    # body weight, and age covariates are added on the logit scale (matches
    # the Gandhi 2021 implementation; see vignette Errata for why the
    # paper's 'multiplicative' equation form is interpreted as additive on
    # the logit). F_abs is bounded in (0, 1) via the inverse logit.
    logit_f <- logitfdepot + etalogitfdepot +
               DIS_PJIA * e_jia_f +
               e_wt_f  * log(WT  / 67.9) +
               e_age_f * log(AGE / 49)
    fdepot  <- 1 / (1 + exp(-logit_f))

    # Two-compartment PK; IV doses go directly to central, SC doses to depot
    # with bioavailability fdepot.
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot) <- fdepot

    # Concentration in mg/L (= ug/mL), matching the addSd unit.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
