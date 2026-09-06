Duke_2024_cefazolin <- function() {
  description <- "Two-compartment population PK model for intravenous cefazolin in infected Indigenous Australian adults with end-stage kidney disease on three-times-weekly intermittent high-flux haemodialysis, receiving a 2 g three-times-weekly post-dialysis regimen. PK is parameterised on unbound drug: the central state carries unbound cefazolin and an explicit second-order albumin-binding exchange (k1 on / k2 off) against a capacity bmax derived from serum albumin carries the bound drug, so total and unbound plasma concentrations are both model outputs. Clearance is replaced (not augmented) by a 41-fold higher dialytic clearance while a session is running, gated by the time-varying RRT_HEMODIAL_ACTIVE covariate; interdialytic clearance falls with the number of months the patient has been established on haemodialysis (T_HEMODIAL_INIT) through an inverse-power relationship, a surrogate for the progressive loss of residual renal function. Estimated with the Pmetrics non-parametric adaptive grid (NPAG). Duke 2024, n = 16 subjects, 130 paired total-and-unbound plasma samples."
  reference <- "Duke C, Parker SL, Zam BB, Chiong F, Sajiv C, Pawar B, Ashok A, Cooper BP, Tong SYC, Janson S, Wallis SC, Roberts JA, Tsai D. Population pharmacokinetics of unbound cefazolin in infected hospitalized patients requiring intermittent high-flux haemodialysis: can a three-times-weekly post-dialysis dosing regimen provide optimal treatment? J Antimicrob Chemother. 2024;79(11):2980-2989. doi:10.1093/jac/dkae318"
  vignette <- "Duke_2024_cefazolin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Mapped from the Table S3 Pmetrics model file:
  # X(1) = unbound drug in the central compartment, X(2) = albumin-bound
  # drug, X(3) = peripheral drug (Figure S1 compartment diagram).
  compartmentData <- list(
    central     = list(analyte = "unbound cefazolin", units = "mg", specimen = "plasma", verified = FALSE),
    complex     = list(analyte = "albumin-bound cefazolin", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "cefazolin", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    T_HEMODIAL_INIT = list(
      description        = "Time the patient has been established on intermittent haemodialysis therapy for end-stage kidney disease, measured from the initiation of that therapy",
      units              = "month",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject in the source analysis. The only covariate retained in the final model (Duke 2024 Results: TOH 'was the only covariate retained in the final pharmacokinetic model'). Enters as an inverse-power effect on the interdialytic clearance arm, CL = CLnHD * (59 / TOH)^0.28 (Duke 2024 Results equation 'When dialysis is off', reproduced verbatim as 'CL=CLnHD*(59/TOH)**0.28' in the Table S3 Pmetrics model file). The 59-month reference is the cohort median TOH (Table 1: 59 months, IQR 24.3-120), so the effect is centred rather than arbitrary -- unlike the bilirubin reference in the sibling model Tsai_2023_ceftriaxone.R. Clearance and TOH followed an inverse-power relationship with r^2 = 0.433. The Discussion reads TOH as 'a surrogate for the incremental reduction in the residual renal function from the initiation of haemodialysis therapy', which is why clearance FALLS as TOH rises; the authors state this covariate had not previously been included in a popPK model for patients requiring intermittent haemodialysis. Must be strictly positive: the covariate enters as a denominator, so TOH = 0 is undefined. The paper's own dosing simulations (Table 4) span TOH = 6, 12, 24, 36 and 60 months.",
      source_name        = "TOH"
    ),
    ALB = list(
      description        = "Serum albumin concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject in the source analysis. Not a covariate on any structural PK parameter; instead it sets the albumin-binding capacity of the central compartment through the Table S3 secondary variable Bmax1 = Alb * Vc * 4.1 (mg). The 4.1 mg cefazolin per g albumin constant encodes the paper's Bmax equation Bmax = Alb * N * (MCFZ / MAlb) * 1000 with N = 0.6 binding sites per albumin molecule, MCFZ = 455 g/mol and MAlb = 66500 g/mol: 0.6 * 455 / 66500 * 1000 = 4.105, rounded to 4.1 in the model file. Cohort median 38.5 g/L (IQR 35.5-40); the Discussion notes the absence of hypoalbuminaemia (< 24 g/L) in this cohort and attributes the unusually high unbound fraction to competitive displacement by uraemia (median pre-dialysis urea 19.6 mmol/L) and by heparin-induced free fatty acids instead.",
      source_name        = "Alb"
    ),
    RRT_HEMODIAL_ACTIVE = list(
      description        = "Haemodialysis-active indicator (1 while an intermittent high-flux haemodialysis session is running, 0 in the interdialytic interval)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no dialysis session running)",
      notes              = "Time-varying within subject. Implemented in the source as the Pmetrics conditional '&IF (HDx.EQ.1) CL=CLHD' (Table S3 secondary variables), i.e. the dialytic clearance REPLACES the interdialytic clearance arm for the duration of the session rather than being added to it. This is the opposite composition rule from the additive dialysis-arm precedents (Veinstein_2013_gentamicin.R, Eyler_2014_ertapenem.R, Jacobs_2016_colistin.R, Dohmann_2025_piperacillin.R) and follows the replacement precedent already set by the same group's Tsai_2023_ceftriaxone.R; it is encoded here as the paper wrote it. Because CLHD replaces CL entirely, the T_HEMODIAL_INIT covariate does not act during a dialysis session. Unbound cefazolin clearance was 41-fold higher during dialysis (16.36 vs 0.40 L/h), which the authors attribute to the high-flux membranes used (FX80 / FX100 / FX120, Fresenius). Doses in this study were given post-dialysis (slow push over 5 min at the completion of the session), so RRT_HEMODIAL_ACTIVE = 0 at the dosing times of the observed data. Median dialysis session duration was 4.0 h (Table 2).",
      source_name        = "HDx"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 16L,
    n_studies        = 1L,
    n_samples        = 130L,
    age_median       = "51 years (IQR 38.8-62.3); full range not reported",
    weight_median    = "69.5 kg (IQR 58.5-76.3); full range not reported",
    sex_female_pct   = 87.5,
    race_ethnicity   = "100% Indigenous Australian (an explicit inclusion criterion). The Discussion notes that this population commences haemodialysis considerably younger than their non-Indigenous counterparts, which is why vein preservation -- and therefore a post-dialysis regimen requiring no separate cannulation -- carries particular weight.",
    disease_state    = "Adults with end-stage kidney disease established on three-times-weekly intermittent high-flux haemodialysis, treated with cefazolin for an active infection or for surgical prophylaxis. Indications as printed in Table 1: line-associated cellulitis (4), wound infection (3), abscess (3), bacteraemia (3), diabetic foot infection (2), periorbital cellulitis (2), surgical prophylaxis (1). Baseline laboratory values (median, IQR): albumin 38.5 g/L (35.5-40), pre-dialysis urea 19.6 mmol/L (16.2-22.8), total bilirubin 7.5 umol/L (6-12), ALP 233 U/L (179-307), GGT 142 U/L (71-192), ALT 9 U/L (5.75-18.75). No adverse drug reactions were reported.",
    renal_function   = "End-stage kidney disease requiring three-times-weekly intermittent haemodialysis. Residual renal function could not be quantified because serum creatinine in maintenance-dialysis patients is dominated by time since the last session (stated study limitation); the authors used months since initiation of haemodialysis (TOH, median 59, IQR 24.3-120) as its surrogate instead. Dialysers were high-flux throughout: FX80 in 5 subjects, FX100 in 10, FX120 in 1 (ultrafiltration coefficients 59, 73 and 87 mL/h/mmHg; surface areas 1.8, 2.2 and 2.5 m^2). Dialysis parameters (Table S4, mean +/- SD): blood flow rate 348 +/- 47 mL/min, ultrafiltration volume 2939 +/- 819 mL, Kt/V 1.69 +/- 0.37, recirculation 11.4 +/- 1.9%. No subject received haemodiafiltration.",
    dose_range       = "2 g cefazolin (Cefazolin-AFT) reconstituted in 10 mL water-for-injection and injected through the arteriovenous fistula or central line as a slow push over 5 min at the completion of each dialysis session, three times weekly.",
    regions          = "Australia (renal dialysis unit of a remote Northern Territory hospital, Alice Springs)",
    protein_binding  = "Measured directly rather than assumed: median unbound fraction 0.38 (IQR 0.32-0.46), roughly double the 0.21 reported for healthy volunteers. Median pre-dialysis unbound trough was 35.7 mg/L (IQR 27.5-45.7) over a 2-day interval and 17.7 mg/L (IQR 13.5-31.4) over a 3-day interval; the lowest pre-dialysis unbound concentration observed in the whole study was 9.1 mg/L. The unbound fraction was higher immediately before dialysis than immediately after (mean 36.5% +/- 6.8% versus 24.5% +/- 14.3%).",
    notes            = "Prospective single-centre population PK study. 260 concentrations (130 total, 130 unbound) from 16 patients. Plasma sampled over two dosing or dialysis intervals: directly before dialysis, immediately after dialysis, then 5, 15, 60 and 1440 min after the dose, then at 48 h or immediately before the next dialysis session (whichever came first), and again before the next session when the interval was 72 h. Total and unbound cefazolin assayed 1-500 mg/L by validated UHPLC-MS/MS (Table S1); the unbound fraction was isolated by ultrafiltration at 37 C with Centrifree devices. Exclusion criteria: pregnancy, cephalosporin allergy, or a requirement for more frequent dialysis."
  )

  ini({
    # Structural parameters: Duke 2024 Table 3, 'Mean' column of the
    # Pmetrics NPAG non-parametric population distribution. Table 3 also
    # reports a 'Median' column; the mean is used as the typical value
    # here (it is the column the Abstract, Results and Discussion quote)
    # and the median is noted per line. Only the seven primary variables
    # of the Table S3 model file are estimated; V and the two half-lives
    # in Table 3 carry footnote 'a' ("Data not available as the entry is
    # manually calculated") and are derived quantities, not model
    # parameters.
    lcl <- log(0.40)
    label("Interdialytic (dialysis-off) clearance CLnHD (L/h)")
    # Duke 2024 Table 3: CLnHD mean 0.40, SD 0.19, CV 46.00%, median 0.39 L/h

    lcl_hemodialysis <- log(16.36)
    label("Intradialytic (dialysis-on) clearance CLHD (L/h)")
    # Duke 2024 Table 3: CLHD mean 16.36, SD 4.26, CV 26.04%, median 16.31 L/h

    lvc <- log(6.51)
    label("Central volume of distribution Vc (L)")
    # Duke 2024 Table 3: Vc mean 6.51, SD 1.30, CV 20.02%, median 6.74 L

    lk1 <- log(2.17)
    label("Second-order cefazolin-albumin association rate constant Kon (L/mg/h)")
    # Duke 2024 Table 3: Kon mean 2.17, SD 0.56, CV 25.86%, median 2.39 L/mg/h

    lk2 <- log(92.13)
    label("First-order cefazolin-albumin dissociation rate constant Koff (1/h)")
    # Duke 2024 Table 3: Koff mean 92.13, SD 15.15, CV 16.44%, median 89.83 1/h
    # Implied KD = Koff / Kon = 92.13 / 2.17 = 42.46 mg/L (Methods: KD = 1/KA = Koff/Kon)

    lk12 <- log(4.01)
    label("Central-to-peripheral rate constant Kcp (1/h)")
    # Duke 2024 Table 3: Kcp mean 4.01, SD 2.54, CV 63.23%, median 3.51 1/h

    lk21 <- log(1.72)
    label("Peripheral-to-central rate constant Kpc (1/h)")
    # Duke 2024 Table 3: Kpc mean 1.72, SD 1.48, CV 86.22%, median 1.17 1/h

    # Time-on-haemodialysis effect on the interdialytic clearance arm. The
    # exponent is hard-coded in the Table S3 model file rather than
    # reported as an estimated parameter in Table 3, so it is encoded as
    # fixed().
    e_t_hemodial_init_cl <- fixed(0.28)
    label("Inverse-power exponent of months on haemodialysis on interdialytic CL (unitless)")
    # Duke 2024 Table S3 secondary variables: CL=CLnHD*(59/TOH)**0.28;
    # same form printed in Results ('When dialysis is off'). Equivalent to
    # (TOH / 59)^-0.28. Supported by the reported inverse-power fit,
    # r^2 = 0.433.

    # Interindividual variability. Pmetrics NPAG estimates a discrete
    # non-parametric distribution rather than a parametric omega matrix;
    # Table 3 summarises that distribution by its mean, SD and CV%. The
    # CV% is carried here into a log-normal random effect using the
    # standard omega^2 = log(CV^2 + 1) identity. This is a parametric
    # APPROXIMATION of a non-parametric distribution (see vignette
    # 'Assumptions and deviations'); it is required to reproduce the
    # paper's own Monte Carlo PTA simulations, which sample the
    # population distribution.
    #   CLnHD : 46.00% CV -> omega^2 = log(0.4600^2 + 1) = 0.191942
    #   CLHD  : 26.04% CV -> omega^2 = log(0.2604^2 + 1) = 0.065608
    #   Vc    : 20.02% CV -> omega^2 = log(0.2002^2 + 1) = 0.039298
    #   Kon   : 25.86% CV -> omega^2 = log(0.2586^2 + 1) = 0.064733
    #   Koff  : 16.44% CV -> omega^2 = log(0.1644^2 + 1) = 0.026669
    #   Kcp   : 63.23% CV -> omega^2 = log(0.6323^2 + 1) = 0.336332
    #   Kpc   : 86.22% CV -> omega^2 = log(0.8622^2 + 1) = 0.555831
    etalcl              ~ 0.191942  # Duke 2024 Table 3 (CLnHD, CV 46.00%)
    etalcl_hemodialysis ~ 0.065608  # Duke 2024 Table 3 (CLHD,  CV 26.04%)
    etalvc              ~ 0.039298  # Duke 2024 Table 3 (Vc,    CV 20.02%)
    etalk1              ~ 0.064733  # Duke 2024 Table 3 (Kon,   CV 25.86%)
    etalk2              ~ 0.026669  # Duke 2024 Table 3 (Koff,  CV 16.44%)
    etalk12             ~ 0.336332  # Duke 2024 Table 3 (Kcp,   CV 63.23%)
    etalk21             ~ 0.555831  # Duke 2024 Table 3 (Kpc,   CV 86.22%)

    # Residual error. Table S3 '#Error' block gives one assay-error
    # polynomial per output equation, identical for both:
    #   0.3, 0.1, 0, 0   ->  SD = 0.3 + 0.1 * conc  (C2 = C3 = 0)
    # so each output carries a 0.3 mg/L additive plus 10% proportional
    # term. The C1 = 0.1 slope is consistent with the Table S1 assay
    # validation (total-cefazolin precision 4.1-5.3%, unbound 3.6-6.3%).
    # NOTE: Pmetrics multiplies this assay polynomial by an estimated
    # noise-inflation factor gamma; the Table S3 file sets the gamma
    # STARTING value 'G=2', and the paper does not report the final
    # estimated gamma anywhere. The assay polynomial is therefore carried
    # here unscaled (equivalent to gamma = 1), which is the minimum-
    # assumption reading of the on-disk file, matching the sibling
    # extraction Tsai_2023_ceftriaxone.R. See vignette 'Assumptions and
    # deviations'.
    addSd <- 0.3
    label("Additive residual error on total Cc (mg/L)")
    # Duke 2024 Table S3 #Error, output 2 (total): C0 = 0.3
    propSd <- 0.1
    label("Proportional residual error on total Cc (fraction)")
    # Duke 2024 Table S3 #Error, output 2 (total): C1 = 0.1
    addSd_Cunbound <- 0.3
    label("Additive residual error on unbound Cunbound (mg/L)")
    # Duke 2024 Table S3 #Error, output 1 (unbound): C0 = 0.3
    propSd_Cunbound <- 0.1
    label("Proportional residual error on unbound Cunbound (fraction)")
    # Duke 2024 Table S3 #Error, output 1 (unbound): C1 = 0.1
  })

  model({
    # Stoichiometric constant for the albumin-binding capacity, carried
    # exactly as hard-coded in the Table S3 secondary variable
    # Bmax1 = Alb * Vc * 4.1. Units: mg cefazolin bound per g albumin.
    # It encodes the Methods equation Bmax = Alb * N * (MCFZ/MAlb) * 1000
    # with N = 0.6 binding sites per albumin molecule --
    #   0.6 * 455 (cefazolin g/mol) / 66500 (albumin g/mol) * 1000 = 4.105
    # -- which the model file rounds to 4.1.
    bmax_per_g_alb <- 4.1

    # Reference time on haemodialysis for the inverse-power clearance
    # covariate (months). Hard-coded in the Table S3 model file; it is the
    # Table 1 cohort median TOH of 59 months.
    toh_ref <- 59

    # Individual parameters.
    cl              <- exp(lcl + etalcl) * (toh_ref / T_HEMODIAL_INIT)^e_t_hemodial_init_cl
    cl_hemodialysis <- exp(lcl_hemodialysis + etalcl_hemodialysis)
    vc              <- exp(lvc + etalvc)
    k1              <- exp(lk1 + etalk1)
    k2              <- exp(lk2 + etalk2)
    k12             <- exp(lk12 + etalk12)
    k21             <- exp(lk21 + etalk21)

    # Dialysis REPLACES the interdialytic clearance arm rather than adding
    # to it (Table S3: '&IF (HDx.EQ.1) CL=CLHD'). Note this differs from
    # the additive dialysis-arm convention used by Veinstein 2013 /
    # Eyler 2014 / Jacobs 2016 / Dohmann 2025; it follows the same group's
    # Tsai 2023 ceftriaxone model and is encoded as Duke 2024 wrote it.
    cl_total <- (1 - RRT_HEMODIAL_ACTIVE) * cl + RRT_HEMODIAL_ACTIVE * cl_hemodialysis
    kel      <- cl_total / vc

    # Albumin-binding capacity of the central compartment, as a MASS (mg)
    # rather than a concentration -- so it is directly comparable with the
    # bound-drug amount held in the 'complex' state (Table S3 Bmax1).
    bmax <- ALB * vc * bmax_per_g_alb

    # ODE system, transcribed from the Table S3 '#Differential equations'
    # block. X(1) -> central (unbound drug), X(2) -> complex (albumin-bound
    # drug), X(3) -> peripheral1. Elimination and inter-compartmental
    # distribution act on unbound drug only; the bound state exchanges
    # solely with central.
    #   XP(1) = RATEIV(1) - (Ke + Kcp)*X(1) - (Kon/Vc)*(Bmax1-X(2))*X(1)
    #                     + Koff*X(2) + Kpc*X(3)
    #   XP(2) =             (Kon/Vc)*(Bmax1-X(2))*X(1) - Koff*X(2)
    #   XP(3) =  Kcp*X(1) - Kpc*X(3)
    # The dose enters 'central' (Pmetrics RATEIV(1)) via the event table.
    # At binding equilibrium this system reproduces the paper's Methods
    # relation Ctotal = Cunbound + Bmax * Cunbound / (KD + Cunbound) with
    # KD = Koff / Kon and Bmax = Alb * 4.1 mg/L.
    d/dt(central) <- -(kel + k12) * central -
      (k1 / vc) * (bmax - complex) * central + k2 * complex + k21 * peripheral1
    d/dt(complex) <-
      (k1 / vc) * (bmax - complex) * central - k2 * complex
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Output equations (Table S3 '#Output equations').
    #   Y(1) = X(1)/Vc          -> unbound plasma concentration
    #   Y(2) = (X(2)+X(1))/Vc   -> total plasma concentration
    Cunbound <- central / vc
    Cc       <- (complex + central) / vc

    Cc       ~ add(addSd) + prop(propSd)
    Cunbound ~ add(addSd_Cunbound) + prop(propSd_Cunbound)
  })
}
