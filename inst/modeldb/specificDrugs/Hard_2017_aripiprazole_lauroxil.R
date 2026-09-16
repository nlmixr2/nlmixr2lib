Hard_2017_aripiprazole_lauroxil <- function() {
  description <- "Two-compartment population PK model (2MPopPK) for aripiprazole released from the long-acting intramuscular prodrug aripiprazole lauroxil, with lagged zero-order IM input and a first-order oral aripiprazole route, in adults with schizophrenia"
  reference <- paste(
    "Hard ML, Mills RJ, Sadler BM, Wehr AY, Weiden PJ, von Moltke L (2017).",
    "Pharmacokinetic Profile of a 2-Month Dose Regimen of Aripiprazole Lauroxil:",
    "A Phase I Study and a Population Pharmacokinetic Model.",
    "CNS Drugs 31(7):617-624. doi:10.1007/s40263-017-0447-7.",
    "Parameter estimates from Supplemental Table 7 of the Electronic Supplementary Material.",
    "Vp/F, Q/F and the 70 kg weight-centering constant are carried from the earlier",
    "aripiprazole lauroxil PopPK model of Hard ML, Mills RJ, Sadler BM, Turncliff RZ,",
    "Citrome L (2017) J Clin Psychopharmacol 37(3):289-295, doi:10.1097/JCP.0000000000000685",
    "(not currently in nlmixr2lib).",
    sep = " "
  )
  vignette <- "Hard_2017_aripiprazole_lauroxil"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # The typical-value parameters of the 2MPopPK model were estimated on the log
  # scale and the exponentiated values are what Supplemental Table 7 prints
  # (Table 7 footnote '*'); the values below are therefore log()-wrapped.
  # Fixed values are flagged '**' in that table.
  paper_specific_residual_sds <- c("propSdStdyA105", "propSdStdyPrior")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Subject-level. Enters as an allometric power effect on the central apparent volume",
        "Vc/F with the exponent fixed at 1.0 (Hard 2017 CNS Drugs Supplemental Table 7 row",
        "'WT ON VC/F' = 1.00, flagged as fixed). The 2MPopPK paper states only that the",
        "weight effect on Vc/F was 'retained' from the earlier PopPK model and never prints",
        "the centering weight; the earlier model supplies it explicitly -- Hard 2017",
        "J Clin Psychopharmacol Supplementary Data Content, section 'Covariate model',",
        "reads 'The effect of WT on both VC/F and CL/F was evaluated using the power model",
        "with values centered for a weight of 70 kg' and prints VC/F = 268*(WT/70)^1.0,",
        "confirmed by its Supplemental Table 2 footnote 'power effect = VC/F * (WT/70)^1.0'.",
        "Weight does NOT act on CL/F in either final model. Cohort mean weight in the phase I",
        "study was 89.8 kg (Hard 2017 CNS Drugs Sect. 3.1); the earlier PopPK dataset spanned",
        "43-144 kg with a median of 81.7 kg."
      ),
      source_name = "WT"
    ),
    CYP2D6_PM = list(
      description = "CYP2D6 poor-metabolizer phenotype indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (extensive, intermediate, inconclusive or missing CYP2D6 phenotype)",
      notes = paste(
        "Subject-level. Multiplicative effect on apparent clearance: CL/F is 23% lower in",
        "poor metabolizers (Hard 2017 CNS Drugs Supplemental Table 7 row 'CL/F PMs' = 0.767,",
        "95% CI 0.614-0.921, footnote 'In reference to non-PM'). Phenotype was assigned from",
        "genotype by the allele-activity tables in Supplemental Tables 5 and 6: poor",
        "metabolizer = no-activity/no-activity. Of the 700 patients in the 2MPopPK dataset,",
        "25 were poor metabolizers, 416 extensive, 183 intermediate, 3 inconclusive and 73",
        "missing; no ultra-rapid metabolizers were present, so the model offers no guidance",
        "for that phenotype (Hard 2017 CNS Drugs Sect. 4). Intermediate metabolizers are",
        "pooled with the non-PM reference -- the model carries no separate IM term."
      ),
      source_name = "CYP2D6 phenotype (PM vs non-PM)"
    ),
    OCC = list(
      description = "Intramuscular injection occasion number",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Time-varying; 1..7, one occasion per aripiprazole lauroxil injection. Multiplexes the",
        "inter-occasion variability on the zero-order input duration D1 (Hard 2017 CNS Drugs",
        "Supplemental Table 7, 'Inter-occasion Variability' block: D1 IOV variance 0.125,",
        "%RSE 5.23, 35.4% CV). Seven occasions are carried because the 2MPopPK model carried",
        "seven IM dosing depots -- 'Seven IM depots were included to accommodate the maximum",
        "number of AL injections in the clinical studies (seven in study A105)'",
        "(Supplemental Fig. 1 footnote). Set OCC on each IM dose record; oral aripiprazole",
        "records and any IM dose beyond the seventh take OCC = 0, which zeroes every indicator",
        "and leaves that dose with no IOV contribution."
      ),
      source_name = "injection number"
    ),
    STUDY_A105 = list(
      description = "Phase I study A105 (NCT02320032) membership indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the four earlier pooled studies 002, 101, 102 and 003)",
      notes = paste(
        "Subject-level. Selects the proportional residual-error magnitude: the 2MPopPK model",
        "estimated separate residual variances for A105 and for the pooled earlier studies",
        "(Hard 2017 CNS Drugs Supplemental Table 7, 'Residual Variability' block:",
        "sigma^2 prop ARI A105 = 0.0238 -> SD 0.154; sigma^2 prop ARI Non-A105 = 0.0565 ->",
        "SD 0.238). A105 contributed 5077 aripiprazole concentrations from 124 patients and",
        "the earlier studies 576 patients (Hard 2017 CNS Drugs Sect. 2.2). A105 sampled far",
        "more densely (up to 48 samples per patient over 309 days), which is the likely",
        "reason its residual error is smaller."
      ),
      source_name = "Study"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "aripiprazole", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "aripiprazole", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "aripiprazole", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 700,
    n_studies = 5,
    age_range = "18-65 years",
    age_median = "47 years (phase I study A105)",
    weight_range = "43-144 kg (earlier PopPK dataset)",
    weight_median = "81.7 kg (earlier PopPK dataset); mean 89.8 kg in phase I study A105",
    sex_female_pct = 27,
    race_ethnicity = c(Black = 75, White = 25),
    disease_state = "schizophrenia or schizoaffective disorder (chronic stable, plus one acute-exacerbation phase III cohort)",
    dose_range = "aripiprazole lauroxil 441-1064 mg IM q4wk / q6wk / q8wk (aripiprazole-equivalent 300-724 mg), plus oral aripiprazole 10-15 mg once daily",
    regions = "USA",
    cyp2d6 = "25 poor metabolizers, 416 extensive, 183 intermediate, 3 inconclusive, 73 missing; no ultra-rapid metabolizers",
    notes = paste(
      "The 2MPopPK dataset pooled 14,524 aripiprazole concentrations (826 [6%] below the",
      "lower limit of quantification, handled by the M3 method) from 700 patients across five",
      "studies: phase I A105 (NCT02320032, 124 patients contributing 5077 concentrations),",
      "phase I studies 002, 101 and 102, and the pivotal phase III study 003. Study 001 from",
      "the earlier PopPK analysis was excluded because it used a non-commercial formulation.",
      "Demographics quoted above are the phase I (A105) safety population of Hard 2017",
      "CNS Drugs Table 1 (mean age 44.5 years, 73% male, 75% Black/African American,",
      "mean weight 89.8 kg) plus the per-study medians of Supplemental Table 2; the paper",
      "does not tabulate pooled demographics for all 700 patients."
    )
  )

  ini({
    # ---- Structural parameters (Supplemental Table 7) ------------------------
    lka <- log(0.803)
    label("First-order absorption rate constant for oral aripiprazole (1/h)")
    # Supplemental Table 7 row 'Ka (h-1)' = 0.803 (%RSE 29.9; 95% CI 0.333, 1.27)

    lcl <- log(1.898)
    label("Apparent clearance of aripiprazole in CYP2D6 non-poor-metabolizers (CL/F, L/h)")
    # Supplemental Table 7 row 'CL/F (L/hr)' = 1.898 (%RSE 2.57; 95% CI 1.80, 1.99)

    lvc <- log(317)
    label("Apparent central volume of distribution at 70 kg (VC/F, L)")
    # Supplemental Table 7 row 'VC/F (L)' = 317 (%RSE 2.25; 95% CI 303, 331)

    lvp <- fixed(log(2122))
    label("Apparent peripheral volume of distribution, from Hard 2017 J Clin Psychopharmacol (VP/F, L)")
    # Supplemental Table 7 row 'VP/F (L)' = 2122, flagged '**' fixed and footnote 'a'
    # 'Fixed at estimate from previous final model'; matches that model's 2122 L exactly.

    lq <- fixed(log(0.423))
    label("Apparent inter-compartmental clearance, from Hard 2017 J Clin Psychopharmacol (Q/F, L/h)")
    # Supplemental Table 7 row 'Q/F (L/hr)' = 0.423, flagged '**' fixed and footnote 'a';
    # matches the earlier model's 0.423 L/h exactly.

    ld1 <- log(1043)
    label("Duration of the zero-order aripiprazole input after an intramuscular aripiprazole lauroxil injection (D1, h)")
    # Supplemental Table 7 row 'D1 (hr)' = 1043 (%RSE 2.09; 95% CI 1000, 1086);
    # 1043 h = 43.5 days, the 'duration of absorption of 43 days' of the Abstract.

    ltlag <- log(77.5)
    label("Lag time before aripiprazole appears in the central compartment after an intramuscular injection (ALAG, h)")
    # Supplemental Table 7 row 'ALAG (hr)' = 77.5 (%RSE 4.07; 95% CI 71.3, 83.7);
    # 77.5 h = 3.2 days, the lag time of the Abstract.

    lc0 <- fixed(log(0.915))
    label("Baseline (pre-first-dose) aripiprazole plasma concentration (ng/mL)")
    # Supplemental Table 7 row 'ARI(0) (ng/mL)' = 0.915, flagged '**' fixed; ESM Sect. 4
    # 'ARI(0) was fixed to a value estimated during base model development (0.915 ng/mL)'.
    # Carries quantifiable pre-dose concentrations from the pre-randomization oral
    # aripiprazole tolerability test doses rather than excluding or baseline-correcting them.

    lfdepot_im <- log(0.571)
    label("Bioavailability of intramuscular aripiprazole lauroxil relative to oral aripiprazole (FIM, unitless)")
    # Supplemental Table 7 row 'FIM' = 0.571 (%RSE 2.55; 95% CI 0.542, 0.599)

    lfdepot <- fixed(log(1))
    label("Bioavailability of oral aripiprazole, the reference route (FPO, unitless)")
    # Supplemental Table 7 row 'FPO' = 1.00, flagged '**' fixed, footnote 'b' 'Fixed as 1.0
    # (as reference for FPO ...)'. Its IIV was also fixed to zero (ESM Sect. 4), so no eta.

    # ---- Covariate effects ---------------------------------------------------
    e_wt_vc <- fixed(1)
    label("Allometric exponent on body weight for VC/F, centred at 70 kg (unitless)")
    # Supplemental Table 7 row 'WT ON VC/F' = 1.00, flagged '**' and footnote 'b'
    # 'at allometric value for WT on VC/F'. Centering weight 70 kg from the earlier model.

    e_cyp2d6_pm_cl <- 0.767
    label("Multiplicative effect of CYP2D6 poor-metabolizer status on CL/F (unitless)")
    # Supplemental Table 7 row 'CL/F PMs' = 0.767 (%RSE 10.2; 95% CI 0.614, 0.921),
    # footnote 'c' 'In reference to non-PM'; a 23% reduction in CL/F.

    # ---- Inter-individual variability ---------------------------------------
    # Supplemental Table 7 'Inter-individual Variability' Value column holds the
    # omega^2 VARIANCE on the log scale. Confirmed against the table's own CV%
    # column and its footnote 'd': CV = sqrt(exp(omega^2) - 1) * 100 for the rows
    # flagged 'd' -- e.g. Ka sqrt(exp(2.67) - 1) = 3.67 -> 367%, VC/F
    # sqrt(exp(0.182) - 1) = 0.447 -> 44.7%. FIM's 32.2% carries no 'd' and is
    # sqrt(0.104) = 0.322, the small-variance approximation the ESM Sect. 1 says
    # was used when the variance is below 0.15.
    etalka ~ 2.67 # Supplemental Table 7 Ka IIV variance 2.67 (95% CI 1.34, 4.00)
    etalcl ~ 0.366 # Supplemental Table 7 CL/F IIV variance 0.366 (95% CI 0.302, 0.430)
    etalvc ~ 0.182 # Supplemental Table 7 VC/F IIV variance 0.182 (95% CI 0.149, 0.215)
    etalvp ~ 2.99 # Supplemental Table 7 VP/F IIV variance 2.99 (95% CI 2.13, 3.85); estimated even though VP/F itself is fixed
    etalq ~ 1.13 # Supplemental Table 7 Q/F IIV variance 1.13 (95% CI 0.895, 1.37); estimated even though Q/F itself is fixed
    etald1 ~ 0.318 # Supplemental Table 7 D1 IIV variance 0.318 (95% CI 0.250, 0.386)
    etaltlag ~ 0.344 # Supplemental Table 7 ALAG IIV variance 0.344 (95% CI 0.276, 0.412)
    etalc0 ~ 5.89 # Supplemental Table 7 ARI(0) IIV variance 5.89 (95% CI 5.00, 6.78); estimated even though ARI(0) itself is fixed
    etalfdepot_im ~ 0.104 # Supplemental Table 7 FIM IIV variance 0.104 (95% CI 0.0803, 0.128)

    # ---- Inter-occasion variability on D1 ------------------------------------
    # One eta per IM injection occasion, sharing a single estimated variance --
    # the NONMEM 'OMEGA BLOCK(1) SAME' idiom. Occasions 2-7 are fixed() to the
    # occasion-1 estimate so only one variance is free, as the paper reports.
    etaiov_d1_1 ~ 0.125 # Supplemental Table 7 'Inter-occasion Variability' D1 variance 0.125 (%RSE 5.23; 95% CI 0.112, 0.345; 35.4% CV)
    etaiov_d1_2 ~ fixed(0.125) # shared variance
    etaiov_d1_3 ~ fixed(0.125) # shared variance
    etaiov_d1_4 ~ fixed(0.125) # shared variance
    etaiov_d1_5 ~ fixed(0.125) # shared variance
    etaiov_d1_6 ~ fixed(0.125) # shared variance
    etaiov_d1_7 ~ fixed(0.125) # shared variance

    # ---- Residual error ------------------------------------------------------
    propSdStdyA105 <- 0.15427
    label("Proportional residual SD for phase I study A105 (fraction)")
    # Supplemental Table 7 'sigma2prop ARI A105' = 0.0238 (%RSE 2.30; 95% CI 0.0227,
    # 0.0249); sqrt(0.0238) = 0.15427, matching the table's own 15.4% CV.

    propSdStdyPrior <- 0.23770
    label("Proportional residual SD for the four earlier pooled studies (fraction)")
    # Supplemental Table 7 'sigma2prop ARI Non-A105' = 0.0565 (%RSE 1.98; 95% CI 0.0543,
    # 0.0587); sqrt(0.0565) = 0.23770, matching the table's own 23.8% CV.
  })

  model({
    # ---- 1. Derived covariate terms -----------------------------------------
    # Inter-occasion variability on the zero-order input duration: exactly one
    # indicator is 1 on an IM dose record carrying OCC in 1..7; OCC = 0 (oral
    # records, or an eighth-and-later injection) zeroes them all.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    oc6 <- (OCC == 6)
    oc7 <- (OCC == 7)
    iov_d1 <-
      oc1 * etaiov_d1_1 + oc2 * etaiov_d1_2 + oc3 * etaiov_d1_3 +
      oc4 * etaiov_d1_4 + oc5 * etaiov_d1_5 + oc6 * etaiov_d1_6 +
      oc7 * etaiov_d1_7

    # ---- 2. Individual parameters -------------------------------------------
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * e_cyp2d6_pm_cl^CYP2D6_PM
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    vp <- exp(lvp + etalvp)
    q <- exp(lq + etalq)
    d1 <- exp(ld1 + etald1 + iov_d1)
    tlag <- exp(ltlag + etaltlag)
    c0 <- exp(lc0 + etalc0)
    fdepot <- exp(lfdepot)
    fdepot_im <- exp(lfdepot_im + etalfdepot_im)

    # ---- 3. Micro-constants --------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- 4. ODE system -------------------------------------------------------
    # Supplemental Fig. 1: seven IM depots and one oral depot all feed a single
    # aripiprazole central compartment, which exchanges with one peripheral
    # compartment and is cleared by CL/F. The seven IM depots are a NONMEM
    # book-keeping device -- the figure's footnote says a new depot was added per
    # injection purely to accommodate up to seven overlapping inputs, and the
    # model carries no rate constant out of an IM depot. Each depot's
    # contribution is a zero-order input of duration D1 starting ALAG after the
    # injection, so the seven depots superpose exactly onto a single lagged
    # modelled-duration input to `central`, which is how it is encoded here.
    # IM dose records must therefore use cmt = "central" with rate = -2.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Baseline aripiprazole carried into the model as an initial condition.
    # c0 is a concentration in ng/mL and vc a volume in L, so the /1000 converts
    # ng/mL * L = ug to the mg in which the compartment amounts are held.
    central(0) <- c0 * vc / 1000

    # ---- 5. Route-specific input ---------------------------------------------
    f(depot) <- fdepot # oral aripiprazole, the FPO = 1 reference route
    f(central) <- fdepot_im # intramuscular aripiprazole lauroxil, relative to oral
    dur(central) <- d1
    alag(central) <- tlag

    # ---- 6. Observation and error --------------------------------------------
    # central is in mg and vc in L, so central/vc is ug/mL; *1000 gives ng/mL.
    Cc <- 1000 * central / vc
    propSdStudy <- propSdStdyA105 * STUDY_A105 + propSdStdyPrior * (1 - STUDY_A105)
    Cc ~ prop(propSdStudy)
  })
}
