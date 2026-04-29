Hu_2014_bapineuzumab <- function() {
  description <- "Two-compartment population PK model for bapineuzumab in adults with mild-to-moderate Alzheimer's disease following IV administration (Hu 2014, reduced model)"
  reference <- "Hu C, Yu G, Tomaszewski EN, et al. Confirmatory population pharmacokinetic analysis for bapineuzumab phase 3 studies in patients with mild to moderate Alzheimer's disease. J Clin Pharmacol. 2015;55(2):221-229. doi:10.1002/jcph.393"
  vignette <- "Hu_2014_bapineuzumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline; power scaling on CL and Vc with reference 70 kg standardized weight (Hu 2014 abstract and Reduced Covariate Model paragraph). Mean +/- SD body weight in the analysis population was 72.6 +/- 15.2 kg (Table 1).",
      source_name        = "WT"
    ),
    RACE_WHITE = list(
      description        = "White (Caucasian) race indicator (1 = Caucasian, 0 = non-Caucasian)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (Caucasian; reference for the typical-value parameters CL = 0.17 L/day and Vc = 3.13 L)",
      notes              = "Multiplicative effect on CL: non-Caucasian subjects have 15% higher CL than the Caucasian reference (Hu 2014 Table 2 'Race on CL' = 1.15). Source paper used a binary Caucasian-vs-non-Caucasian dichotomy without further decomposition. Renamed from the source descriptor 'race' (Caucasian vs non-Caucasian) to canonical RACE_WHITE per inst/references/covariate-columns.md; the typical-value reference is the Caucasian (RACE_WHITE = 1) subgroup, which is the inverse of the Lin 2024 use of the same canonical column.",
      source_name        = "RACE"
    )
  )

  population <- list(
    n_subjects     = 1458,
    n_studies      = 2,
    age_range      = "Mean +/- SD 73.0 +/- 8.8 years (Hu 2014 Table 1, overall)",
    weight_range   = "Mean +/- SD 72.6 +/- 15.2 kg (Hu 2014 Table 1, overall)",
    sex_female_pct = 30.7,
    race_ethnicity = "Predominantly Caucasian (Study 301: 93.6% Caucasian per Hu 2014 Methods; Study 302 was similarly composed). Non-Caucasian subjects pooled into the model's non-Caucasian indicator.",
    disease_state  = "Mild-to-moderate Alzheimer's disease (baseline MMSE 16-26)",
    dose_range     = "0.5 or 1.0 mg/kg IV (1-hour infusion) every 13 weeks for up to six infusions; an additional 141 patients in Study 301 received 2.0 mg/kg before the protocol was amended and they were transitioned to 1.0 mg/kg",
    regions        = "Study 301 (APOE4 non-carriers) and Study 302 (APOE4 carriers); both Phase 3 multinational trials (ELN115727-301 and ELN115727-302).",
    apoe4_status   = "Study 301 enrolled APOE*E4 non-carriers; Study 302 enrolled APOE*E4 carriers. APOE*E4 carrier status was tested as a covariate and had no meaningful effect on PK.",
    n_observations = 8040,
    n_blq_excluded = 100,
    bapineuzumab_subjects_only = TRUE,
    notes          = "n_subjects = 1458 reflects the bapineuzumab-treated subjects whose serum samples (8040 measurements) were analyzed (Hu 2014 abstract and Results). Hu 2014 Table 1 reports n = 1937 for the wider covariate-evaluation dataset that includes placebo subjects. Sex breakdown: Study 301 enrolled 69.3% male per Hu 2014 Methods; Study 302 sex split is not reported in the trimmed text but the paper's Table 1 covers continuous variables only. The 30.7% female value here is computed from the Study 301 male fraction and is documented as an approximation in the vignette."
  )

  ini({
    # Structural parameters for bapineuzumab. Reference subject for the typical
    # values: Caucasian (RACE_WHITE = 1), 70 kg standardized body weight
    # (Hu 2014 abstract; Reduced Covariate Model paragraph: "the typical
    # values for CL and Vc in a Caucasian subject with a standardized body
    # weight of 70 kg as 0.17 L/day ... and 3.13 L").
    lcl     <- log(0.17);   label("Bapineuzumab clearance for Caucasian, 70 kg reference (CL, L/day)")    # Hu 2014 Table 2: CL reduced model estimate 0.17
    lvc     <- log(3.13);   label("Bapineuzumab central volume of distribution for 70 kg reference (Vc, L)")  # Hu 2014 Table 2: Vc reduced model estimate 3.13
    lq      <- log(0.871);  label("Bapineuzumab inter-compartmental clearance (Q, L/day)")                # Hu 2014 Table 2: Q reduced model estimate 0.871
    lvp     <- log(3.61);   label("Bapineuzumab peripheral volume of distribution (Vp, L)")               # Hu 2014 Table 2: Vp reduced model estimate 3.61

    # Power exponents on body weight (reference 70 kg). Both estimated in the
    # reduced covariate model. The Sensitivity Analysis paragraph notes that
    # fixing these to allometric values 0.75 (CL) and 1.0 (Vc) gave similar
    # results, but the reduced model uses the estimated values.
    e_wt_cl   <- 0.64;  label("Power exponent: body weight on CL, (WT/70)^e_wt_cl (unitless)")  # Hu 2014 Table 2: Weight on CL reduced model estimate 0.64
    e_wt_vc   <- 0.78;  label("Power exponent: body weight on Vc, (WT/70)^e_wt_vc (unitless)")  # Hu 2014 Table 2: Weight on Vc reduced model estimate 0.78

    # Race effect on CL. Source paper (Eq. 2) parameterizes discrete covariate
    # effects as u_i = b^X_i * u with X_i = 0 for the reference patient and
    # X_i = 1 for the comparator. Caucasian (RACE_WHITE = 1) is the reference
    # so the indicator used in the model is (1 - RACE_WHITE). The reduced
    # model estimate b = 1.15 corresponds to a 15% higher CL in non-Caucasian
    # subjects; e_nonwhite_cl = 0.15 is the fractional-change form used here.
    e_nonwhite_cl <- 0.15;  label("Fractional increase in CL for non-Caucasian vs Caucasian reference (unitless)")  # Hu 2014 Table 2: Race on CL reduced model estimate 1.15

    # Inter-individual variability for the log-normal random effects
    # P_i = theta * exp(eta_i). Hu 2014 Table 2 footnote explicitly defines
    # "BSV, between-subject variability, calculated as (variance)^(1/2)*100%"
    # so the reported BSV percentages are the square root of the NONMEM
    # internal omega^2 (not approximate CV%). Therefore:
    #   omega^2(CL) = (28.1/100)^2 = 0.078961
    #   omega^2(Vc) = (33.0/100)^2 = 0.108900
    #   cov(eta_CL, eta_Vc) = r * sqrt(omega^2(CL) * omega^2(Vc))
    #                       = 0.494 * sqrt(0.078961 * 0.108900) = 0.045809.
    etalcl + etalvc ~ c(0.078961, 0.045809, 0.108900)  # Hu 2014 Table 2 reduced model: BSV(CL) 28.1, BSV(Vc) 33.0 (defined as sqrt(variance)*100%), correlation 0.494

    # Residual variability. Hu 2014 Methods state "log-transform-both-sides
    # approach for all model runs"; Methods Eq. for residual error is
    # Y_ij = C_ij * exp(e_ij) with e_ij ~ N(0, sigma^2). The reduced model
    # value 0.413 (Table 2 footnote: "additive error in the log domain")
    # maps to a proportional error in linear space with propSd = 0.413.
    propSd <- 0.413; label("Bapineuzumab proportional residual error (fraction)")  # Hu 2014 Table 2: residual error reduced model estimate 0.413
  })
  model({
    # Multiplicative race effect on CL relative to Caucasian reference
    # (Hu 2014 Eq. 2 with X_i = 1 for non-Caucasian).
    race_cl <- 1 + e_nonwhite_cl * (1 - RACE_WHITE)

    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * race_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q  <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Concentration: dose in mg, volume in L -> mg/L = ug/mL (matches the
    # Hu 2014 ELISA assay range 1.98-147.2 ng/mL with 1:4 minimum dilution
    # and the y-axis scale on Figure 2).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
