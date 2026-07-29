Aregbe_2012_alvespimycin <- function() {
  description <- "Three-compartment population PK model for the heat shock protein 90 inhibitor 17-DMAG (alvespimycin, NSC 707545) given as a 1 h IV infusion to adult patients with advanced solid tumors (Aregbe 2012), with first-order elimination, log-normal IIV on CL/Q3/V1/V2/V3, and between-occasion variability on Q2 and V1 multiplexed by an OCC indicator across up to five daily dosing occasions."
  reference <- paste(
    "Aregbe AO, Sherer EA, Egorin MJ, Scher HI, Solit DB, Ramanathan RK,",
    "Ramalingam S, Belani CP, Ivy PS, Bies RR. (2012). Population",
    "pharmacokinetic analysis of 17-dimethylaminoethylamino-17-",
    "demethoxygeldanamycin (17-DMAG) in adult patients with solid tumors.",
    "Cancer Chemother Pharmacol 70(2):201-205.",
    "doi:10.1007/s00280-012-1859-1.",
    sep = " "
  )
  vignette <- "Aregbe_2012_alvespimycin"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    OCC = list(
      description        = "Integer-valued occasion / period indicator for between-occasion variability multiplexing.",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Values 1, 2, 3, 4, 5 identify the daily 1 h infusion occasion within subject (the Pittsburgh schedule A used 5 daily doses, schedule B used 3, and the MSKCC cohort used a single dose, so OCC ranges 1-5 across the pooled dataset). Decomposed inside `model()` into binary indicators `oc1` .. `oc5` that multiplex the BOV etas on log-Q (NONMEM Q2) and log-Vc (NONMEM V1).",
      source_name        = "OCC"
    )
  )

  population <- list(
    n_subjects     = 67L,
    n_studies      = 2L,
    age_range      = "28-82 years (median 63)",
    age_median     = "63 years",
    weight_range   = "48.2-136.5 kg (median 80.3)",
    weight_median  = "80.3 kg",
    sex_female_pct = 37,
    disease_state  = "Adult patients with histologically confirmed advanced solid tumors not curable by standard therapies; required adequate hepatic, renal, and cardiac function (ALT/AST <= 1.5 x ULN, normal BUN and creatinine, ECOG <= 2, no QTc prolongation or cardiac comorbidity per the strict 17-AAG-derived exclusion criteria).",
    dose_range     = "1 h IV infusion of 17-DMAG; per-dose range 2.2-413 mg/m^2 (median 36 mg/m^2). Pittsburgh patients followed schedule A (5 daily doses) or schedule B (3 daily doses) under an accelerated dose-titration design; MSKCC patients received a single pre-specified dose.",
    regions        = "United States (University of Pittsburgh Cancer Institute, n = 48; Memorial Sloan-Kettering Cancer Center, n = 19).",
    notes          = "Demographics from Aregbe 2012 Table 1; the cohort was 63% male / 37% female with a median ECOG performance status compatible with phase II eligibility. Baseline laboratory medians: albumin 3.8 g/dL (range 2.6-5.1, missing in 6 subjects), ALT 22 U/L (10-106), AST 25 U/L (12-75), bilirubin 0.5 mg/dL (0.1-3.0), BUN 15 mg/dL (5-70, missing in 2 subjects), creatinine 1.0 mg/dL (0.6-1.8), BSA 1.9 m^2 (1.5-2.6); 1 subject missing demographics other than centre. The paper screened age, albumin, ALT, AST, bilirubin, BUN, BSA, creatinine, weight, and sex by stepwise forward addition / backward elimination but retained no covariate effects in the final model (Aregbe 2012 Results, page 203)."
  )

  ini({
    # Structural PK parameters from Aregbe 2012 Table 2 (final model NONMEM
    # ADVAN11 / TRANS4 estimates). NONMEM names map to nlmixr2 conventions as
    # CL -> lcl, V1 -> lvc, Q2 -> lq, V2 -> lvp, Q3 -> lq2, V3 -> lvp2.
    lcl  <- log(8.4)  ; label("Clearance (L/h)")                               # Aregbe 2012 Table 2: CL = 8.4 (%SE 11.2)
    lvc  <- log(27.4) ; label("Central volume of distribution V1 (L)")         # Aregbe 2012 Table 2: V1 = 27.4 (%SE 11.7)
    lq   <- log(85.1) ; label("Inter-compartmental clearance Q2 (L/h)")        # Aregbe 2012 Table 2: Q2 = 85.1 (%SE 9.6)
    lvp  <- log(66.4) ; label("Peripheral volume of distribution V2 (L)")      # Aregbe 2012 Table 2: V2 = 66.4 (%SE 10.1)
    lq2  <- log(11.6) ; label("Inter-compartmental clearance Q3 (L/h)")        # Aregbe 2012 Table 2: Q3 = 11.6 (%SE 13.1)
    lvp2 <- log(142)  ; label("Peripheral volume of distribution V3 (L)")      # Aregbe 2012 Table 2: V3 = 142 (%SE 13.5)

    # Inter-individual variability (log-normal). The paper reports IIV as a
    # CV%; the internal variance is omega^2 = log(CV^2 + 1). IIV was retained
    # on CL, Q3, V1, V2, V3 in the final model (Aregbe 2012 Table 2, IIV
    # column); no IIV on Q2 (reported as `-`).
    etalcl  ~ 0.2543  # Aregbe 2012 Table 2: IIV CL  53.8% (%SE 22.9), variance = log(0.538^2 + 1)
    etalvc  ~ 0.1034  # Aregbe 2012 Table 2: IIV V1  33.0% (%SE 111.9), variance = log(0.330^2 + 1)
    etalvp  ~ 0.2288  # Aregbe 2012 Table 2: IIV V2  50.7% (%SE 23.7), variance = log(0.507^2 + 1)
    etalq2  ~ 0.4539  # Aregbe 2012 Table 2: IIV Q3  75.8% (%SE 32.0), variance = log(0.758^2 + 1)
    etalvp2 ~ 0.3756  # Aregbe 2012 Table 2: IIV V3  67.5% (%SE 37.3), variance = log(0.675^2 + 1)

    # Between-occasion variability (BOV) on Q2 (NONMEM `$OMEGA BLOCK(1)` reused
    # across occasions per the paper's `BOV on Q2`). Each occasion-1..5 gets an
    # eta drawn from the same single-element variance; nlmixr2 has no `SAME`
    # shortcut (Jonsson_2011_ethambutol pattern), so the first eta is the
    # estimated variance and the remainder are fixed equal to it.
    etaiov_q_1 ~ 0.1010        # Aregbe 2012 Table 2: BOV Q2 32.6% (%SE 37.2), variance = log(0.326^2 + 1)
    etaiov_q_2 ~ fix(0.1010)   # NONMEM $OMEGA BLOCK(1) SAME (occasion 2 shares occasion-1 variance)
    etaiov_q_3 ~ fix(0.1010)   # NONMEM $OMEGA BLOCK(1) SAME (occasion 3 shares occasion-1 variance)
    etaiov_q_4 ~ fix(0.1010)   # NONMEM $OMEGA BLOCK(1) SAME (occasion 4 shares occasion-1 variance)
    etaiov_q_5 ~ fix(0.1010)   # NONMEM $OMEGA BLOCK(1) SAME (occasion 5 shares occasion-1 variance)

    # Between-occasion variability (BOV) on V1 (NONMEM `$OMEGA BLOCK(1)` reused
    # across occasions per the paper's `BOV on V1`). Same `SAME` translation
    # pattern as the Q2 block above. V1 carries both IIV (etalvc above) and
    # BOV (etaiov_vc_1..5) in the final model.
    etaiov_vc_1 ~ 0.3049       # Aregbe 2012 Table 2: BOV V1 59.7% (%SE 31.2), variance = log(0.597^2 + 1)
    etaiov_vc_2 ~ fix(0.3049)  # NONMEM $OMEGA BLOCK(1) SAME (occasion 2 shares occasion-1 variance)
    etaiov_vc_3 ~ fix(0.3049)  # NONMEM $OMEGA BLOCK(1) SAME (occasion 3 shares occasion-1 variance)
    etaiov_vc_4 ~ fix(0.3049)  # NONMEM $OMEGA BLOCK(1) SAME (occasion 4 shares occasion-1 variance)
    etaiov_vc_5 ~ fix(0.3049)  # NONMEM $OMEGA BLOCK(1) SAME (occasion 5 shares occasion-1 variance)

    # Proportional residual error model best described the data per the paper;
    # additive and combined-error structures were tested but not retained.
    propSd <- 0.161 ; label("Proportional residual error (fraction)")          # Aregbe 2012 Table 2: proportional error 16.1% (%SE 2.7)
  })

  model({
    # Decompose the integer-valued occasion column into binary indicators for
    # BOV multiplexing on log-Q and log-Vc. Up to five daily 1 h infusions
    # (Pittsburgh schedule A); patients on schedule B contribute 3 occasions
    # and the MSKCC cohort contributes 1 occasion -- the unused indicator
    # rows simply zero out for those subjects.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)

    iov_q  <- oc1 * etaiov_q_1  + oc2 * etaiov_q_2  + oc3 * etaiov_q_3  + oc4 * etaiov_q_4  + oc5 * etaiov_q_5
    iov_vc <- oc1 * etaiov_vc_1 + oc2 * etaiov_vc_2 + oc3 * etaiov_vc_3 + oc4 * etaiov_vc_4 + oc5 * etaiov_vc_5

    # Individual PK parameters. No covariates are retained in the final
    # model (Aregbe 2012 Results: none of the screened covariates produced
    # an OFV drop greater than the forward-addition / backward-elimination
    # significance threshold). V1 carries both IIV and BOV; Q2 carries BOV
    # only.
    cl  <- exp(lcl  + etalcl)
    vc  <- exp(lvc  + etalvc + iov_vc)
    q   <- exp(lq             + iov_q)
    vp  <- exp(lvp  + etalvp)
    q2  <- exp(lq2  + etalq2)
    vp2 <- exp(lvp2 + etalvp2)

    # Three-compartment IV PK with first-order elimination from the central
    # compartment (NONMEM ADVAN11 / TRANS4 mass-balance ODEs). Dose lands in
    # `central` via the user data set's cmt column; the 1 h infusion duration
    # is encoded on the dose record (rate or dur).
    d/dt(central)     <-  q  / vp  * peripheral1 + q2 / vp2 * peripheral2 -
                          (cl + q + q2) / vc * central
    d/dt(peripheral1) <-  q  / vc  * central     - q  / vp  * peripheral1
    d/dt(peripheral2) <-  q2 / vc  * central     - q2 / vp2 * peripheral2

    # Plasma concentration in the central compartment. Dose units mg, Vc units
    # L -> Cc units mg/L (= ug/mL).
    Cc <- central / vc

    # Proportional residual error on the linear concentration scale.
    Cc ~ prop(propSd)
  })
}
