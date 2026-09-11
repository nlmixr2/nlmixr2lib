Prohn_2021_letermovir_healthy <- function() {
  description <- "Four-compartment population PK model for letermovir given orally or intravenously to healthy phase I participants over a 30-960 mg dose range. Absorption is a Savic analytical transit chain (3.58 transit compartments) whose mean transit time increases with dose. Both elimination and the intercompartmental clearance to the fast peripheral compartment are concentration-dependent, CL = EAI * CLmax / (1 + Cc/KMcl) and Q1 = Q1max / (1 + Cc/KMq), which produces the greater-than-proportional rise in exposure observed over the dose range. Clearance is additionally auto-induced through a turnover enzyme pool driven by plasma concentration, which reproduces the observed fall in trough concentrations on repeated dosing. Body weight scales CLmax and all four volumes, and Asian participants have a 28.1 percent lower volume of distribution."
  reference <- paste(
    "Prohn M, Viberg A, Zhang D, Dykstra K, Davis C, Macha S, Sabato P,",
    "de Alwis D, Iwamoto M, Fancourt C, Cho CR (2021). Population",
    "pharmacokinetics of letermovir following oral and intravenous",
    "administration in healthy participants and allogeneic hematopoietic",
    "cell transplantation recipients. CPT Pharmacometrics Syst Pharmacol",
    "10(3):255-267. doi:10.1002/psp4.12593.",
    "This file encodes the healthy participant (phase I) model of Table 2 and",
    "Figure 1a; the HSCT recipient (phase III) model of Table 3 and Figure 1b",
    "is a separate, independently fitted model and is encoded in",
    "Prohn_2021_letermovir_hsct.R.",
    sep = " "
  )
  vignette <- "Prohn_2021_letermovir"
  # Figure 1b of the source defines Cp = A_central / Vc * 1000, and both
  # Michaelis-Menten constants of this model are tabulated in ng/mL, so the
  # observation is scaled from the mg/L that mg and L give directly.
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot       = list(analyte = "letermovir", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "letermovir", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "letermovir", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "letermovir", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral3 = list(analyte = "letermovir", units = "mg", specimen = "plasma", verified = TRUE),
    enzyme      = list(analyte = "relative amount of the clearing enzyme system", units = "(fraction of baseline)", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power-form scaling normalised at 65.8 kg, which is the reference",
        "weight printed in the Figure 1a equations rather than a rounded",
        "standard; the cohort median was 66 kg (Table 1). Two separately",
        "ESTIMATED exponents rather than fixed allometric 0.75 / 1 values:",
        "0.566 on CLmax and 0.667 on the volume of distribution. The volume",
        "exponent applies to Vd as a whole, defined by the source as",
        "V1 + V2 + V3 + V4, and Figure 1a writes the covariate model as",
        "V(1..4) = V(1..4)base * (1 + Vdjpn) * (WT/65.8)^Vdwt, so the same",
        "factor multiplies all four volumes. Body weight was NOT retained in",
        "the companion HSCT recipient (phase III) model.",
        sep = " "
      ),
      source_name        = "WT"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Asian, predominantly White)",
      notes              = paste(
        "Enters the volume of distribution as a LINEAR fractional shift, not",
        "a log-additive or power one: Figure 1a prints",
        "V(1..4) = V(1..4)base * (1 + Vdjpn) * (WT/65.8)^Vdwt with",
        "Vdjpn = -0.281, so the factor is 1 - 0.281 = 0.719 in Asian",
        "participants and exactly 1 otherwise. The Discussion confirms the",
        "magnitude and the direction: 'Asian participants had a 28.1% lower",
        "Vd compared with White participants'. The Results quote a 27.6%",
        "reduction from the same coefficient; the 28.1% figure matches the",
        "tabulated -0.281 exactly and is used here. The covariate was",
        "originally fitted as Japanese ethnicity (n = 30) and later replaced",
        "by Asian ethnicity (n = 33), 'which slightly improved the model' --",
        "hence the jpn subscript surviving in the figure's notation for what",
        "is an Asian-race indicator.",
        sep = " "
      ),
      source_name        = "Vd-jpn / ASIAN"
    ),
    DOSE = list(
      description        = "Administered letermovir dose per administration",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Use case (a) of the DOSE canonical: the per-subject assigned dose",
        "level entering a covariate effect on the mean transit time,",
        "normalised at 240 mg. Mean transit time INCREASES with dose, i.e.",
        "absorption slows at higher doses. The functional form was settled",
        "against the source's own derived value rather than read off the",
        "figure, because the Figure 1a rendering lost the grouping of its",
        "final term. Results: 'MTT in the TCAM model for a 240 mg oral dose",
        "estimated to be 1.4 h', against a tabulated TVMTT of 1.04 h and",
        "MTTdose of 0.344. Of the three readings the figure could bear, only",
        "MTT = TVMTT * (1 + DOSE/240 * MTTdose) reproduces it:",
        "1.04 * (1 + 1 * 0.344) = 1.398 h. A power form",
        "TVMTT * (DOSE/240)^MTTdose returns TVMTT itself, 1.04 h, and the",
        "literal rendering TVMTT * (1 + DOSE/240) * MTTdose returns 0.716 h;",
        "both are excluded. Observed dose levels were 30-960 mg.",
        sep = " "
      ),
      source_name        = "DOSE"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Tested on clearance and volume of distribution in the stepwise covariate search and not retained (Supplementary Information, 'Healthy participant (phase I model) covariate analysis'). Median 30 years, range 18-59 (Table 1).",
      source_name = "AGE"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested and not retained: 'although most of the healthy participants were female, gender was not found to have a significant covariate effect on PK, and individual parameter estimates stratified by gender were not significantly different.' 254 of 280 participants (91%) were female, because most phase I studies recruited only women following preclinical testicular toxicity in rats.",
      source_name = "SEX"
    ),
    SNP_SLCO1B1_RS4149056 = list(
      description = "OATP1B1 (SLCO1B1) rs4149056 genotype",
      units       = "(genotype)",
      type        = "categorical",
      notes       = "Tested and not retained. 224 participants had single nucleotide polymorphism data for OATP1B1 rs4149056 and rs2306283 and UGT1A1 rs4148323; 'these functional variants had no statistically significant effect on letermovir exposure and they were not included in the final model'. No coefficient is reported, so the effect cannot be encoded even as a fixed zero with provenance.",
      source_name = "OATP1B1 rs4149056"
    ),
    SNP_SLCO1B1_RS2306283 = list(
      description = "OATP1B1 (SLCO1B1) rs2306283 genotype",
      units       = "(genotype)",
      type        = "categorical",
      notes       = "Tested and not retained; see SNP_SLCO1B1_RS4149056. No coefficient reported.",
      source_name = "OATP1B1 rs2306283"
    ),
    SNP_UGT1A1_RS4148323 = list(
      description = "UGT1A1 rs4148323 genotype",
      units       = "(genotype)",
      type        = "categorical",
      notes       = "Tested and not retained; see SNP_SLCO1B1_RS4149056. No coefficient reported.",
      source_name = "UGT1A1 rs4148323"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 280,
    n_studies      = 12,
    age_range      = "median 30 years, range 18-59 (Table 1)",
    weight_range   = "median 66 kg, range 45-99 (Table 1)",
    sex_female_pct = 91,
    race_ethnicity = "Reported only as the Asian / non-Asian contrast retained in the model. The covariate was first fitted as Japanese ethnicity (n = 30) and replaced in the final model by Asian ethnicity (n = 33); the remainder are predominantly White. Counts by race are not otherwise tabulated.",
    disease_state  = "Healthy volunteers.",
    dose_range     = "30-960 mg, single and multiple doses, orally or intravenously",
    regions        = "Not reported; one of the pooled phase I studies was conducted in an Asian population.",
    notes          = paste(
      "Pooled across 12 phase I studies. 9020 concentration observations,",
      "6391 (71%) after oral and 2629 (29%) after intravenous dosing, and",
      "4680 (52%) after single and 4340 (48%) after multiple dosing. 174",
      "observations (1.9%) below the 1 ng/mL lower limit of quantification",
      "were excluded. Most phase I studies recruited only female participants",
      "because of preclinical testicular toxicity in rats, so 91% of the",
      "pooled cohort is female. Fitted in NONMEM 7.3; parameter precision was",
      "assessed by a 500-replicate trial-stratified non-parametric bootstrap.",
      "The source notes that the model does not describe the full",
      "concentration-time profile at the extremes of the dose range (30 mg and",
      "960 mg), although goodness-of-fit plots showed no structural bias, and",
      "that it fits well over 60-720 mg.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters -- Prohn 2021 Table 2, "Estimate" column.
    # Table 2 reports arithmetic values; they are carried as logs here so the
    # log-normal IIV terms add on the natural-log scale.
    # ------------------------------------------------------------------
    lclmax <- log(12.3)
    label("Maximal clearance at a reference 65.8 kg body weight (L/h)")  # Table 2 'Clearance Vmax, L/h' CLmax 12.3 (RSE 2.60%, bootstrap 95% CI 11.3-13.3)
    lvc <- log(7.46)
    label("Central volume of distribution at reference (L)")  # Table 2 'Central volume of distribution, L' V1 7.46 (RSE 6.30%, CI 6.94-7.89)
    lvp <- log(61.6)
    label("Fast peripheral volume of distribution at reference (L)")  # Table 2 'Peripheral volume, L' V2 61.6 (RSE 4.80%, CI 55.3-68.2)
    lvp2 <- log(12.1)
    label("Second peripheral volume of distribution at reference (L)")  # Table 2 'Peripheral volume, L' V3 12.1 (RSE 4.20%, CI 11.4-12.9)
    lvp3 <- log(19.0)
    label("Third peripheral volume of distribution at reference (L)")  # Table 2 'Peripheral volume, L' V4 19.0 (RSE 4.00%, CI 17.4-21.2)
    lq <- log(4.39)
    label("Maximal intercompartmental clearance to peripheral1 (L/h)")  # Table 2 'Intercompartment clearance Vmax, L/h' Q1max 4.39 (RSE 4.90%, CI 3.93-4.95)
    lq2 <- log(31.3)
    label("Intercompartmental clearance to peripheral2 (L/h)")  # Table 2 'Intercompartment clearance, L/h' Q2 31.3 (RSE 12.3%, CI 27.2-36.6)
    lq3 <- log(4.91)
    label("Intercompartmental clearance to peripheral3 (L/h)")  # Table 2 'Intercompartment clearance, L/h' Q3 4.91 (RSE 4.90%, CI 4.49-5.39)
    lfdepot <- log(0.938)
    label("Oral bioavailability (fraction)")  # Table 2 'Bioavailability' F1 0.938 (RSE 2.10%, CI 0.906-0.974); Results quote "93.8%"

    # Saturable-elimination and saturable-distribution constants. Both are
    # tabulated in ng/mL, which is why the observation below is scaled to
    # ng/mL rather than left in mg/L.
    lkm_cl <- log(2680)
    label("Michaelis-Menten constant for clearance (ng/mL)")  # Table 2 'Michaelis-Menten constant, ng/ml' KMCL 2.68e3 (RSE 4.50%, CI 2.20e3-3.41e3)
    lkm_q <- log(5630)
    label("Michaelis-Menten constant for intercompartmental clearance (ng/mL)")  # Table 2 'Michaelis-Menten constant, ng/ml' KMQ1 5.63e3 (RSE 10.3%, CI 3.81e3-7.32e3)

    # Transit-compartment absorption (Savic parameterisation). NTR is not an
    # integer: it is an estimated continuous chain length, which the analytical
    # Savic form accommodates through the gamma function.
    lntr <- log(3.58)
    label("Number of transit compartments (unitless)")  # Table 2 'Number of transit compartments' NTR 3.58 (RSE 1.70%, CI 3.12-4.09)
    lmtt <- log(1.04)
    label("Mean transit time intercept (h)")  # Table 2 'Mean transit time, h' MTT 1.04 (RSE 4.20%, CI 0.949-1.14)

    # Auto-induction of clearance. Figure 1a:
    #   d(EAI)/dt = kin * (1 + IMAG * Cc / 1000) - kout * EAI
    # Table 2 reports kout and IMAG but NOT kin. kin is not missing data: the
    # enzyme pool is a unitless relative amount whose baseline is kin/kout, and
    # clearance is CL = EAI * CLmax / (1 + Cc/KMcl), so CLmax is identifiable
    # as "the maximal clearance" -- the quantity Table 2 tabulates and to which
    # the body-weight covariate is referenced -- only if the pool sits at
    # exactly 1 in the absence of drug. That forces kin = kout, and the same
    # equality is the established idiom for this construct in the registry
    # (Wicha_2018_rifampicin.R writes d/dt(enzyme) <- kenz * (1 + eff) -
    # kenz * enzyme with enzyme(0) <- 1). kout is therefore the only estimated
    # rate constant and kin is derived from it in model().
    lkout <- log(0.00783)
    label("Turnover rate constant of the enzyme induction pool (1/h)")  # Table 2 'Turnover rate induction,/h' Kout 0.00783 (RSE not applicable, bootstrap CI 0.00580-0.00966)
    imag <- 0.0829
    label("Slope of the concentration-driven induction effect (per ug/mL)")  # Table 2 'Slope of induction effect' IMAG 0.0829 (RSE not calculated, bootstrap CI 0.0672-0.102)

    # ------------------------------------------------------------------
    # Covariate effects
    # ------------------------------------------------------------------
    e_wt_clmax <- 0.566
    label("Power exponent of body weight on CLmax, normalised at 65.8 kg (unitless)")  # Table 2 'Weight effect on clearance Vmax' CLmax-wt 0.566 (RSE 26.4%, CI 0.376-0.741)
    e_wt_vd <- 0.667
    label("Power exponent of body weight on all four volumes, normalised at 65.8 kg (unitless)")  # Table 2 'Weight effect on Vd' Vd-wt 0.667 (RSE 10.9%, CI 0.452-0.854)
    e_asian_vd <- -0.281
    label("Fractional effect of Asian race on all four volumes (unitless)")  # Table 2 'Asian effect on Vd' Vd-jpn -0.281 (RSE 8.80%, CI -0.334 to -0.227); enters as (1 + e_asian_vd * RACE_ASIAN), i.e. 28.1% lower Vd
    e_dose_mtt <- 0.344
    label("Fractional effect of dose on mean transit time, normalised at 240 mg (unitless)")  # Table 2 'Dose effect on MTT' MTTdose 0.344 (RSE 11.5%, CI 0.258-0.433)

    # ------------------------------------------------------------------
    # Inter-individual variability. As in the companion HSCT model, Table 2's
    # IIV column holds the NONMEM $OMEGA elements, i.e. log-scale VARIANCES,
    # and is carried verbatim. The scale is settled for the phase III model by
    # that model's reported exposure prediction intervals (see
    # Prohn_2021_letermovir_hsct.R); Table 2 uses the same unlabelled,
    # unitless column style and the same "IIV, <parameter>" row labels, so the
    # same reading applies. Note the contrast with the residual-error row
    # immediately below it, which IS labelled with a unit ("%").
    #
    # No IIV was estimated on V3, V4, Q2, Q3 or NTR.
    # ------------------------------------------------------------------
    etalclmax ~ 0.0647
    label("IIV on maximal clearance")  # Table 2 'IIV, clearance Vmax' 0.0647 (RSE 9.90%, CI 0.0533-0.0766); CV 25.7%
    etalvc ~ 0.0642
    label("IIV on central volume")  # Table 2 'IIV, central volume' 0.0642 (RSE 30.7%, CI 0.0387-0.0977); CV 25.6%
    etalmtt ~ 0.0841
    label("IIV on mean transit time")  # Table 2 'IIV, mean transit time' 0.0841 (RSE 14.2%, CI 0.0700-0.0976); CV 29.4%
    etalq ~ 0.155
    label("IIV on maximal intercompartmental clearance")  # Table 2 'IIV, intercompartmental clearance' 0.155 (RSE 16.0%, CI 0.126-0.203); CV 40.9%
    etalvp ~ 0.224
    label("IIV on fast peripheral volume")  # Table 2 'IIV, peripheral volume' 0.224 (RSE 11.4%, CI 0.184-0.273); CV 50.1%
    etalfdepot ~ 0.0323
    label("IIV on bioavailability")  # Table 2 'IIV, bioavailability' 0.0323 (RSE 13.3%, CI 0.0232-0.0425); CV 18.1%

    # Residual error. Unlike the IIV rows this one carries an explicit unit in
    # its row label -- Table 2 heads it "Proportional residual variability, %"
    # -- so 0.283 is a 28.3% proportional error on the SD scale, not a
    # variance that would imply a 53% error.
    propSd <- 0.283
    label("Proportional residual error (fraction)")  # Table 2 'Proportional residual variability, %' 0.283 (RSE 0.300%, CI 0.271-0.294)
  })

  model({
    # --- Covariate model (Prohn 2021 Figure 1a) --------------------------
    # V(1..4) = V(1..4)base * (1 + Vdjpn) * (WT/65.8)^Vdwt. One shared factor
    # multiplies all four volumes; the Asian term is a linear fractional shift
    # and the weight term is a power.
    vdcov <- (1 + e_asian_vd * RACE_ASIAN) * (WT / 65.8)^e_wt_vd

    vc  <- exp(lvc  + etalvc) * vdcov
    vp  <- exp(lvp  + etalvp) * vdcov
    vp2 <- exp(lvp2)          * vdcov
    vp3 <- exp(lvp3)          * vdcov

    # CLmax = CLmbase * (WT/65.8)^CLMwt
    clmax <- exp(lclmax + etalclmax) * (WT / 65.8)^e_wt_clmax
    q1max <- exp(lq + etalq)
    q2    <- exp(lq2)
    q3    <- exp(lq3)
    km_cl <- exp(lkm_cl)
    km_q  <- exp(lkm_q)

    fdepot <- exp(lfdepot + etalfdepot)

    # MTT = TVMTT * (1 + DOSE/240 * MTTdose); see covariateData$DOSE for how
    # the form was identified from the source's own 1.4 h value at 240 mg.
    mtt <- exp(lmtt + etalmtt) * (1 + DOSE / 240 * e_dose_mtt)
    ntr <- exp(lntr)
    ktr <- (ntr + 1) / mtt

    kout <- exp(lkout)
    # Forced equal to kout so the enzyme pool has a baseline of exactly 1 and
    # CLmax retains the meaning Table 2 assigns it; see the ini() comment.
    kin <- kout

    # --- Plasma concentration --------------------------------------------
    # Amounts are in mg and volumes in L, so central/vc is mg/L; the factor of
    # 1000 converts to the ng/mL in which both Michaelis-Menten constants are
    # tabulated. Cc is computed before the ODE right-hand sides so the
    # concentration-dependent terms below read the current value.
    Cc <- central / vc * 1000

    # --- Concentration-dependent clearance and distribution (Figure 1a) ---
    # CL = EAI * CLmax / (1 + Cc/KMcl) and Q1 = Q1max / (1 + Cc/KMq). Both
    # fall as concentration rises, which is what produces the greater than
    # proportional increase in exposure with dose.
    cl <- enzyme * clmax / (1 + Cc / km_cl)
    q  <- q1max / (1 + Cc / km_q)

    # --- Savic analytical transit-chain absorption ------------------------
    # Figure 1a: Dose -> A_tr1 -> ... -> A_trn -> Central, all at rate Ktr,
    # with Ktr = (NTR + 1) / MTT. The chain empties directly into the central
    # compartment; there is no separate first-order absorption step. The
    # analytical Savic 2007 input rate,
    #   rate = F1 * Dose * Ktr * (Ktr * t)^NTR * exp(-Ktr * t) / gamma(NTR + 1),
    # reproduces the chain without integrating its members, and accommodates
    # the non-integer NTR = 3.58 through the gamma function. It is written out
    # here rather than delegated to rxode2's transit() macro because the
    # chain's terminus is central rather than a first-order depot.
    #
    # depot exists only so that podo(depot) and tad(depot) have a dose record
    # to read; f(depot) <- 0 keeps the dose from also entering as a bolus, and
    # podo(depot) returns the raw dose amount regardless of f(depot). Route is
    # therefore selected entirely by the event table: oral doses go to depot
    # and are delivered through the chain with bioavailability F1, intravenous
    # doses go straight to central and bypass both.
    tdose  <- tad(depot)
    absorb <- fdepot * podo(depot) * ktr *
      exp(ntr * log(ktr * tdose) - ktr * tdose - lgamma(ntr + 1))

    # --- ODE system --------------------------------------------------------
    d/dt(depot)       <- 0
    d/dt(central)     <- absorb - cl * central / vc +
      q  * (peripheral1 / vp  - central / vc) +
      q2 * (peripheral2 / vp2 - central / vc) +
      q3 * (peripheral3 / vp3 - central / vc)
    d/dt(peripheral1) <- q  * (central / vc - peripheral1 / vp)
    d/dt(peripheral2) <- q2 * (central / vc - peripheral2 / vp2)
    d/dt(peripheral3) <- q3 * (central / vc - peripheral3 / vp3)

    # Enzyme induction pool, driven by plasma concentration in ug/mL (hence
    # the division of the ng/mL Cc by 1000). Baseline 1 = uninduced.
    d/dt(enzyme)      <- kin * (1 + imag * Cc / 1000) - kout * enzyme
    enzyme(0)         <- 1

    f(depot) <- 0

    Cc ~ prop(propSd)
  })
}
