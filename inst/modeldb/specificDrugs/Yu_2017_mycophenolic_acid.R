Yu_2017_mycophenolic_acid <- function() {
  description <- "Two-compartment population pharmacokinetic model with first-order oral absorption (no lag) for mycophenolic acid (MPA, the active component of mycophenolate mofetil, MMF) in Chinese adult renal transplant recipients (Yu 2017). Apparent clearance CL/F follows a linear-additive covariate model in body weight and serum creatinine (CL/F = 0.0916 * BW + 0.0417 * Scr + 7.98 L/h); apparent central volume V1/F follows a linear-additive covariate model in the UGT2B7 211G>T (rs7438135) genotype (V1/F = 14.7 + 7.72 * UGT2B7 L) where the paper's ordinal genotype code maps 211GT to 1, 211GG to 2, and 211TT to 3. Residual error is combined proportional plus additive on plasma MPA. Interoccasion variability (13.7% CV on CL/F and V1/F) reported by the paper is documented in the vignette but not encoded in this typical-value model because Karlsson-Sheiner IOV requires per-occasion etas and an OCC data column; the IIV-only encoding remains usable for typical-value and IIV-only simulations."
  reference <- paste(
    "Yu Z-C, Zhou P-J, Wang X-H, Francoise B, Xu D, Zhang W-X, Chen B.",
    "Population pharmacokinetics and Bayesian estimation of mycophenolic acid",
    "concentrations in Chinese adult renal transplant recipients.",
    "Acta Pharmacologica Sinica. 2017;38(11):1566-1579.",
    "doi:10.1038/aps.2017.115.",
    sep = " "
  )
  vignette <- "Yu_2017_mycophenolic_acid"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight at the time of the pharmacokinetic evaluation.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters the CL/F covariate model as a linear additive term: CL/F = 0.0916 * BW + 0.0417 * Scr + 7.98 (Yu 2017 page 1574 / Table 3 footnote). Cohort mean (population group, n = 79) is 58.0 +/- 9.33 kg (range 39-82); cohort mean (full group, n = 118) is 58.3 +/- 9.91 kg (range 36.8-94). No reference-weight normalisation is applied because the paper's model is additive rather than power-form.",
      source_name        = "BW"
    ),
    CREAT = list(
      description        = "Serum creatinine concentration at the time of the pharmacokinetic evaluation.",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters the CL/F covariate model as a linear additive term: CL/F = 0.0916 * BW + 0.0417 * Scr + 7.98 (Yu 2017 page 1574 / Table 3 footnote). Cohort mean (population group, n = 79) is 141.2 +/- 128.5 umol/L (range 61-915); cohort mean (validation group, n = 39) is 129.3 +/- 64.6 umol/L (range 64-328). The paper's CL/F increases with Scr -- a counter-intuitive direction for a renally-influenced apparent clearance that the authors attribute (Discussion page 1576-1577) to higher MPA free fraction under acidosis / uremia / MPAG accumulation, which raises the unbound substrate concentration available for glucuronidation and thereby increases apparent CL/F.",
      source_name        = "Scr"
    ),
    UGT2B7_211GG = list(
      description        = "UGT2B7 211G>T (rs7438135) homozygous G/G (Ala71/Ala71) genotype indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any non-G/G UGT2B7 211 genotype: 211GT or 211TT).",
      notes              = "One of three binary indicators reconstructing Yu 2017's ordinal column 'UGT2B7 genotype' (Table 5; original values 211GT, 211GG, 211TT). The paper encodes the ordinal column as GT = 1, GG = 2, TT = 3 (back-solved from V1/F = 7.72 * UGT2B7 + 14.7 against the Table 5 group means 24.2, 30.8, 36.9 L). The model() block derives the ordinal code as ugt2b7_211_code = UGT2B7_211GT * 1 + UGT2B7_211GG * 2 + UGT2B7_211TT * 3. In the n = 79 population group there were 42 PK evaluations with the 211GG genotype.",
      source_name        = "UGT2B7 genotype"
    ),
    UGT2B7_211GT = list(
      description        = "UGT2B7 211G>T (rs7438135) heterozygous G/T (Ala71/Ser71) genotype indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any non-G/T UGT2B7 211 genotype: 211GG or 211TT).",
      notes              = "Companion to UGT2B7_211GG and UGT2B7_211TT (see UGT2B7_211GG notes for the ordinal reconstruction). The 211GT heterozygous genotype was the most common in Yu 2017's cohort (51 of 101 PK evaluations across the population group).",
      source_name        = "UGT2B7 genotype"
    ),
    UGT2B7_211TT = list(
      description        = "UGT2B7 211G>T (rs7438135) homozygous T/T (Ser71/Ser71) genotype indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any non-T/T UGT2B7 211 genotype: 211GG or 211GT).",
      notes              = "Companion to UGT2B7_211GG and UGT2B7_211GT (see UGT2B7_211GG notes for the ordinal reconstruction). The 211TT homozygous variant was the rarest in Yu 2017's cohort (8 of 101 PK evaluations across the population group).",
      source_name        = "UGT2B7 genotype"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 118L,
    n_studies      = 1L,
    n_population_group = 79L,
    n_validation_group = 39L,
    n_pk_evaluations  = "1 evaluation in all 118 patients; 31 patients had 2 evaluations; 4 patients had 3 evaluations. Total 1172 plasma MPA concentrations (783 in the population group, 389 in the validation group).",
    age_range      = "18-76 years (population group 18-68; validation group 21-76)",
    age_mean       = "41.4 +/- 11.2 years (population group); 44.3 +/- 12.0 years (validation group); 42.5 +/- 11.4 years overall",
    weight_range   = "36.8-94 kg (population group 39-82; validation group 36.8-94)",
    weight_mean    = "58.0 +/- 9.33 kg (population group); 58.9 +/- 11.1 kg (validation group); 58.3 +/- 9.91 kg overall",
    sex_female_pct = 39.8,
    race_ethnicity = "Chinese (single-center cohort at Ruijin Hospital, Shanghai JiaoTong University School of Medicine).",
    disease_state  = "Adult renal transplant recipients receiving triple-immunosuppression with MMF + cyclosporine (or tacrolimus for 12 of 118 patients) + corticosteroids. Pharmacokinetic evaluations were performed 3 to 1460 days after the start of MMF therapy; most occurred within the first month post-transplant.",
    dose_range     = "Oral MMF 1.0 g preoperatively then 2.0 g/day divided BID (target 1.0 g q12h) with clinical adjustment for tolerability. Population-group mean dose 900.1 +/- 177.0 mg per 12 h (range 250-1250); validation-group mean 919.4 +/- 200.7 mg per 12 h (range 500-1500).",
    regions        = "China (Shanghai, single-center).",
    co_medications = "Cyclosporine (CsA, Neoral) in 104 of 118 patients with target C0 200-250 ug/L and C2h 1200-1500 ug/L during month 1, tapered to C0 150-200 ug/L thereafter; tacrolimus (Prograf) in the remaining 12 patients with target trough 10-15 ug/L during week 1, tapered to 5-10 ug/L thereafter; intravenous methylprednisolone 500 mg at surgery, tapered to oral prednisone 5-20 mg/day. CsA was administered 2 h after MMF; tacrolimus was co-administered with MMF.",
    ugt2b7_distribution = "Across the population group (n = 79 patients, 101 PK evaluations) the UGT2B7 211G>T genotype was distributed as 211GT in 51 evaluations (40 patients), 211GG in 42 evaluations (32 patients), and 211TT in 8 evaluations (7 patients).",
    sampling_window= "Full pharmacokinetic profiles used 10 blood samples drawn before (C0) and at 0.5, 1, 1.5, 2, 4, 6, 8, 10, and 12 h after the morning dose. Sparse-sampling occasions used 3 or 4 time-points: 0, 0.5, 2 h or 0, 0.5, 2, and 8 h.",
    notes          = "Single-center retrospective study. Plasma MPA concentrations measured by validated HPLC (LOQ 0.25 mg/L, intraday CV < 6%, interday CV < 8%, accuracy 97.9-101.8%). Patient characteristics from Yu 2017 Table 1."
  )

  ini({
    # Parameter values are the final-model estimates from Yu 2017 Table 3
    # 'With covariates' column (n = 79 population group) plus the
    # additive-covariate formulae stated explicitly on page 1574 just below
    # Table 3:
    #   CL/F = 0.0916 * BW + 0.0417 * Scr + 7.98 = 18.3 L/h (typical)
    #   V1/F = 7.72  * UGT2B7 + 14.7            = 27.9 L   (typical)
    #
    # The paper estimates the structural model in terms of CL/F, V1/F, k12,
    # k21, and ka with exponential IIV on each. This packaged model uses
    # the canonical (CL, Vc, Q, Vp, ka) parameterisation (per the
    # nlmixr2lib parameter-names register) with q = k12 * vc and
    # vp = q / k21 derived in model(). At the typical-value level, this is
    # the same structural model -- the typical values of Q and Vp recover
    # the typical k12 and k21. The IIVs were estimated by the paper on k12
    # and k21 directly; etalq and etalvp here carry the same %CV magnitudes
    # so the simulated marginal variability in q and vp matches what the
    # paper reported for k12 and k21 (the marginal distributions of q and
    # vp are dominated by k12 and k21 when vc is held at its typical
    # value); see vignette Assumptions and deviations for the full
    # discussion of the IIV translation.
    #
    # CV%-to-log-normal-omega translation: omega^2 = log(1 + CV^2). The
    # source CV%s come from Yu 2017 Table 3 column 'IIV [IOV] CV%' with
    # covariates.

    # CL/F = additive linear model in BW and Scr. lcl encodes the LOG of
    # the additive INTERCEPT (7.98 L/h), not the log of typical CL/F. The
    # log() wrapper keeps the intercept positive under any IIV; see
    # vignette Assumptions and deviations for the rationale.
    lcl        <- log(7.98)              ; label("CL/F additive intercept (L/h); typical CL/F is 18.3 L/h at cohort-mean BW and Scr")  # Yu 2017 p1574 formula 'CL/F=0.0916*BW+0.0417*Scr+7.98'
    e_wt_cl    <- 0.0916                  ; label("CL/F additive effect of BW (L/h per kg)")                                              # Yu 2017 Table 3 'a=0.0916 (12.7)'
    e_creat_cl <- 0.0417                  ; label("CL/F additive effect of CREAT (L/h per umol/L of serum creatinine)")                   # Yu 2017 Table 3 'b=0.0417 (34.3)'

    # V1/F = additive linear model in UGT2B7 ordinal code. lvc encodes the
    # LOG of the additive INTERCEPT (14.7 L), not the log of typical V1/F.
    lvc          <- log(14.7)            ; label("V1/F additive intercept (L); typical V1/F is 27.9 L at the cohort-weighted UGT2B7 mix") # Yu 2017 p1574 formula 'V1/F=7.72*UGT2B7+14.7'
    e_ugt2b7_vc  <- 7.72                  ; label("V1/F additive effect per UGT2B7 ordinal-code unit (L); ordinal code GT=1, GG=2, TT=3")  # Yu 2017 Table 3 'a=7.72 (6.85)'

    # Inter-compartmental clearance and peripheral volume. Derived from
    # the paper's k12 = 0.915 1/h and k21 = 0.059 1/h with the typical
    # vc = exp(lvc) + 7.72 * ugt2b7_weighted = 27.9 L:
    #   q  = k12 * vc = 0.915 * 27.9 = 25.5 L/h
    #   vp = q / k21  = 25.5 / 0.059 = 432.7 L
    lq  <- log(25.5)                      ; label("Inter-compartmental clearance Q/F (L/h); derived q = k12 * vc(typical)")              # Yu 2017 Table 3 k12 = 0.915 1/h, with typical vc 27.9 L
    lvp <- log(432.7)                     ; label("Peripheral volume Vp/F (L); derived vp = q / k21")                                    # Yu 2017 Table 3 k21 = 0.059 1/h, with derived q 25.5 L/h

    # Absorption rate constant
    lka <- log(1.89)                      ; label("First-order absorption rate constant ka (1/h)")                                       # Yu 2017 Table 3 ka = 1.89 1/h

    # Inter-individual variability (exponential IIV, log-normal):
    #   omega^2 = log(1 + CV^2)
    # The four primary CL/F, V1/F, k12, k21 IIVs were transferred to
    # etalcl, etalvc, etalq, etalvp respectively (see header comment).
    etalcl ~ 0.110                         # CL/F %CV = 34.2 -> log(1 + 0.342^2) = 0.110
    etalvc ~ 0.0444                        # V1/F %CV = 21.3 -> log(1 + 0.213^2) = 0.0444
    etalq  ~ 0.0929                        # k12  %CV = 31.2 -> log(1 + 0.312^2) = 0.0929 (transferred to etalq; see vignette Assumptions)
    etalvp ~ 1.066                         # k21  %CV = 138  -> log(1 + 1.38^2)  = 1.066  (transferred to etalvp; see vignette Assumptions)
    etalka ~ 0.234                         # ka   %CV = 51.3 -> log(1 + 0.513^2) = 0.234

    # Residual error (Yu 2017 page 1572 / Table 3 'With covariates' column):
    # combined proportional + additive model.
    propSd <- 0.158                       ; label("Proportional residual error on MPA (fraction)")  # Yu 2017 Table 3 sigma_e1 = 15.8% with covariates
    addSd  <- 0.15                        ; label("Additive residual error on MPA (mg/L)")           # Yu 2017 Table 3 sigma_e2 = 0.15 mg/L with covariates
  })

  model({
    # Reconstruct the paper's ordinal UGT2B7 code from the three binary
    # indicators (Yu 2017 Table 5 / page 1574 formula). The ordinal
    # assignment is GT = 1, GG = 2, TT = 3, back-solved from the V1/F
    # group means in Table 5 against the paper's stated linear formula
    # V1/F = 7.72 * UGT2B7 + 14.7.
    ugt2b7_211_code <- UGT2B7_211GT * 1 + UGT2B7_211GG * 2 + UGT2B7_211TT * 3

    # Typical parameter values (no IIV), evaluated at the per-subject
    # covariates. The CL/F and V1/F formulae are LINEAR ADDITIVE per the
    # paper (not power-form). exp(lcl) = 7.98 is the CL/F additive
    # intercept; exp(lvc) = 14.7 is the V1/F additive intercept.
    cl_typical <- exp(lcl) + e_wt_cl * WT + e_creat_cl * CREAT
    vc_typical <- exp(lvc) + e_ugt2b7_vc * ugt2b7_211_code

    # Individual parameters. Inter-individual variability is multiplicative
    # log-normal on the typical-value covariate-conditional CL/F and V1/F,
    # matching Yu 2017's exponential IIV form
    # 'Pj = P_mean * exp(eta_j)' (page 1572).
    cl <- cl_typical * exp(etalcl)
    vc <- vc_typical * exp(etalvc)
    q  <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)
    ka <- exp(lka + etalka)

    # Two-compartment dispositional ODEs with first-order oral absorption.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # MPA plasma concentration in the central compartment.
    Cc <- central / vc

    # Combined proportional + additive residual error on MPA.
    Cc ~ prop(propSd) + add(addSd)
  })
}
