vanSteeg_2007_atenolol_noniso_rat <- function() {
  description <- paste(
    "Preclinical (rat).",
    "Three-compartment IV-infusion PK of S(-)-atenolol in conscious male Wistar-Kyoto (WKY) rats",
    "linked to a DIRECT-effect Emax model for the reduction of RESTING heart rate (beats/min).",
    "This is the NON-ISOPRENALINE (vehicle) group -- the negative control for the companion",
    "isoprenaline-tachycardia model modellib('vanSteeg_2007_atenolol_iso_rat').",
    "It is deliberately the paper's UNDER-DETERMINED arm and is packaged as such: the drug effect",
    "on resting heart rate (Emax = -43 beats/min) is barely larger than the ~30 beats/min",
    "spontaneous variation in rat heart rate, so EC50 is estimated with 97% RSE, and neither a",
    "Hill coefficient nor an effect-compartment ke0 could be estimated at all -- n is FIXED at 1",
    "and there is NO effect compartment, unlike the isoprenaline arm.",
    "The PK was fitted jointly across both groups (no PK difference was found), so the six PK",
    "parameters are identical to those of the isoprenaline model; only the PD differs.",
    "Note the two non-standard variability scales, both taken verbatim from the paper:",
    "PK inter-individual variability is EXPONENTIAL (Equation 1) whereas PD inter-individual",
    "variability is ADDITIVE in beats/min (Equation 5), and the PD residual error is additive in",
    "beats/min despite the paper's Equation 6 printing a proportional form (see the ini() notes).",
    sep = " "
  )

  reference <- paste(
    "van Steeg TJ, Freijer J, Danhof M, de Lange ECM.",
    "Pharmacokinetic-pharmacodynamic modelling of S(-)-atenolol in rats:",
    "reduction of isoprenaline-induced tachycardia as a continuous pharmacodynamic endpoint.",
    "British Journal of Pharmacology. 2007;151(3):356-366.",
    "doi:10.1038/sj.bjp.0707234. PMCID: PMC2013984.",
    "S(-)-atenolol PK parameters are Table 3; PD parameters are Table 4,",
    "'Value (CV) non-ISO' column.",
    "Companion model files from the same paper:",
    "modellib('vanSteeg_2007_atenolol_iso_rat'), modellib('vanSteeg_2007_isoprenaline_rat').",
    sep = " "
  )

  vignette <- "vanSteeg_2007_atenolol_isoprenaline_rat"

  units <- list(time = "min", dosing = "ng", concentration = "ng/mL")

  # No covariates. The paper fits absolute (not weight-normalised) volumes and clearances,
  # so body weight enters only when converting the reported mg/kg dose into the absolute
  # infused amount. See population$notes.
  covariateData <- list()

  compartmentData <- list(
    central     = list(analyte = "S(-)-atenolol", units = "ng", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "S(-)-atenolol", units = "ng", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "S(-)-atenolol", units = "ng", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "rat (Wistar-Kyoto, WKY; male)",
    n_subjects     = 8L,
    n_studies      = 1L,
    age_range      = "not reported; acclimatised at least 5 days, instrumented 7 days before the experiment",
    weight_range   = "294 g +/- 49 (mean +/- SD across the 39 WKY rats of this study)",
    sex_female_pct = 0,
    disease_state  = paste(
      "Normotensive (WKY) conscious rats at resting heart rate. The continuous IV infusion",
      "carried vehicle only -- 0.1% (w/w) sodium metabisulphite in saline -- rather than",
      "isoprenaline. No disease model.",
      sep = " "
    ),
    dose_range     = "S(-)-atenolol 5 mg/kg as a 15-minute IV infusion (single dose), on a background of continuous IV vehicle.",
    regions        = "The Netherlands (Leiden University)",
    biomarkers     = paste(
      "Heart rate captured continuously from a femoral-artery pressure signal (P10EZ-1 transducer,",
      "Spike 2 acquisition) throughout the 8-hour experiment and used directly as the PD endpoint.",
      "Plasma S(-)-atenolol by reversed-phase HPLC with fluorescence detection (ex 235 nm,",
      "em 300 nm) after liquid-liquid extraction; calibration range 40-20000 ng/mL, LOQ 40 ng/mL;",
      "intra- and inter-assay variability 4% and 11%. Serial arterial samples at 0, 5, 10, 15,",
      "17.5, 20, 22.5, 27.5, 32.5, 40, 55, 70, 90, 120, 180, 240, 360 and 480 min.",
      sep = " "
    ),
    notes = paste(
      "SUBJECT COUNT. n_subjects = 8 is the non-isoprenaline group. The PD of the two groups was",
      "fitted SEPARATELY ('In contrast to the PK data, the PD data were analysed separately, since",
      "the baseline heart rate differed between the isoprenaline and the non-isoprenaline group'),",
      "which is why this paper contributes two S(-)-atenolol model files rather than one stratified",
      "file. The PK, by contrast, was fitted on the pooled 17 rats of both groups, so the Table 3",
      "PK block is identical in the two files.",
      "WHY THIS ARM IS UNDER-DETERMINED, and what that means for users. This model is packaged",
      "because it is the paper's stated NEGATIVE RESULT, not because it is a well-identified PD",
      "model: 'it is shown that the measurement of the heart-rate effects of S(-)-atenolol without",
      "isoprenaline-induced tachycardia resulted in imprecise estimation of the potency of the drug",
      "(EC50)'. EC50 = 27.9 ng/mL carries a 97.1% RSE; 'it was not possible to estimate a Hill",
      "coefficient (n) and keo in the non-isoprenaline group, as incorporation of these parameters",
      "in the PD model resulted in increased imprecision in the EC50'. The mechanism is that the",
      "drug effect and the noise are the same size: Emax = -43 beats/min against a ~30 beats/min",
      "spontaneous variation in rat heart rate from circadian rhythm, movement and stress, and an",
      "inter-animal SD in baseline of 37 beats/min that 'approximates the maximal drug effect'.",
      "Simulations from this model are correspondingly wide; use the isoprenaline arm for any",
      "quantitative potency question.",
      "DOSE UNITS: the paper reports a 5 mg/kg dose but fits absolute volumes (mL) and clearances",
      "(mL/min), so a body weight is needed to build an event table. The validation vignette doses",
      "the paper's own typical 294 g rat (Discussion), for which 5 mg/kg is 1.47 mg = 1.47e6 ng.",
      "The body weight is NOT a model covariate.",
      sep = " "
    )
  )

  ini({
    # ==================================================================
    # PHARMACOKINETICS -- van Steeg 2007 Table 3 ("Population parameter
    # estimates including variabilities and BS replicates for CL, V1, V2,
    # V3 and Q2, Q3"), 'Value (CV)' column. Standard three-compartment
    # model (Methods/Pharmacokinetics). Fitted on the pooled 17 rats of
    # both PD groups, so this block is byte-identical in the companion
    # isoprenaline model file.
    #
    # NONMEM V1/V2/V3 and Q2/Q3 map onto the nlmixr2lib canonicals as
    # V1 -> vc, V2 -> vp, V3 -> vp2, Q2 -> q, Q3 -> q2.
    #
    # The '(CV)' column is the RELATIVE STANDARD ERROR of the estimate,
    # not a between-animal CV: Results states "All PK parameters were
    # estimated precisely with coefficients of variation ranging between
    # 3 and 36%", which brackets exactly the Table 3 range (3.4% to
    # 36.7%).
    # ==================================================================
    lcl  <- log(11.7) ; label("Clearance CL (mL/min)")                                       # Table 3 CL = 11.7 mL/min (RSE 3.4%; bootstrap 11.7)
    lvc  <- log(115)  ; label("Central volume of distribution V1 (mL)")                      # Table 3 V1 = 115 mL (RSE 7.4%; bootstrap 114)
    lq   <- log(15.0) ; label("Inter-compartmental clearance to peripheral1, Q2 (mL/min)")   # Table 3 Q2 = 15.0 mL/min (RSE 7.2%; bootstrap 15.2)
    lvp  <- log(173)  ; label("Volume of the second compartment V2 (mL)")                    # Table 3 V2 = 173 mL (RSE 7.7%; bootstrap 172)
    lq2  <- log(8.50) ; label("Inter-compartmental clearance to peripheral2, Q3 (mL/min)")   # Table 3 Q3 = 8.50 mL/min (RSE 5.3%; bootstrap 8.50)
    lvp2 <- log(849)  ; label("Volume of the third compartment V3 (mL)")                     # Table 3 V3 = 849 mL (RSE 4.2%; bootstrap 848)

    # ------------------------------------------------------------------
    # PK inter-individual variability. Exponential, per Equation 1:
    #   P_i = theta * exp(eta_i),  eta_i ~ N(0, omega^2)
    # "Interindividual variability was identified for clearance (CL), the
    # volume of the second compartment (V2) and the volume of the third
    # compartment (V3) ... A significant correlation was obtained for CL
    # and V3 and this correlation was estimated in the final model."
    #
    # Table 3 reports VARIANCES, which the covariance row settles: an
    # off-diagonal of 0.026 alongside diagonals of 0.033 and 0.023
    # implies a correlation of 0.026 / sqrt(0.033 * 0.023) = 0.944, which
    # is admissible; reading the diagonals as SDs instead would imply
    # 0.026 / (0.033 * 0.023) = 34.3, which is impossible. The block is
    # positive definite (determinant 8.3e-5), so no numerical nudge is
    # needed.
    # ------------------------------------------------------------------
    # Block entries in source order (a comment INSIDE the c(...) breaks the
    # rxode2 ini() parser, so they are listed here instead):
    #   var(etalcl)             = Table 3 omega^2_CL = 0.033 (RSE 9.8%; bootstrap 0.032)  -> CV 18.2%
    #   cov(etalcl, etalvp2)    = Table 3 omega^2_CL,V3 = 0.026 (RSE 31.9%; bootstrap 0.024)
    #   var(etalvp2)            = Table 3 omega^2_V3 = 0.023 (RSE 36.7%; bootstrap 0.022)  -> CV 15.2%
    etalcl + etalvp2 ~ c(0.033, 0.026, 0.023)
    etalvp ~ 0.170    # Table 3 omega^2_V2 = 0.170 (RSE 28.2%; bootstrap 0.169) -> CV 43.2%

    # ------------------------------------------------------------------
    # PK residual error. Proportional, per Equation 2:
    #   C_obs = C_pred * (1 + eps),  eps ~ N(0, sigma^2)
    # Table 3 reports the VARIANCE sigma^2 = 0.027, so the SD that
    # nlmixr2's prop() wants is its square root (0.164, i.e. 16.4%).
    # The row is labelled "sigma^2_PD" even though Table 3 is the PK
    # table; that subscript is a copy-paste slip carried across all four
    # tables of the paper (see the vignette Errata). Its placement under
    # the "Residual error" heading is unambiguous.
    # ------------------------------------------------------------------
    propSd <- sqrt(0.027) ; label("Proportional residual error on plasma S(-)-atenolol (fraction)")  # Table 3 sigma^2 = 0.027 (RSE 9.5%; bootstrap 0.027)

    # ==================================================================
    # PHARMACODYNAMICS -- van Steeg 2007 Table 4, 'Value (CV) non-ISO'
    # column ("Population estimates of PD parameters and variabilities
    # for S(-)-atenolol with (ISO) and without (non-ISO)
    # isoprenaline-induced tachycardia").
    #
    # DIRECT-effect Emax on the PLASMA concentration, i.e. Equation 3
    # driven by Cp rather than by the Equation 4 effect compartment:
    #   E = E0 + Emax * Cp^n / (EC50^n + Cp^n)
    # There is deliberately NO effect compartment here. Table 4 records
    # ke0 as "NE" (not estimated) in this column, and the Results explain
    # why: "it was not possible to estimate a Hill coefficient (n) and
    # keo in the non-isoprenaline group, as incorporation of these
    # parameters in the PD model resulted in increased imprecision in the
    # EC50."
    #
    # SIGN OF Emax. Table 4 prints Emax with a leading minus sign that
    # the PDF's symbol font drops in every text extraction (the same
    # dropout that turns "+/-" into "7" and "S(-)-atenolol" into
    # "S()-atenolol" throughout this paper); the column alignment in the
    # PDF shows the extra leading character on the Emax row and not on
    # the E0 row. The prose is explicit: "the E max was a reduction of
    # 43 +/- 18 b.p.m." and "Baseline heart rate in the non-isoprenaline
    # group was 362 +/- 20 b.p.m. and the E max was a reduction of
    # 43 +/- 17.7 b.p.m."
    #
    # E0 and Emax are NOT log-transformed. The paper models the PD
    # parameters ADDITIVELY (Equation 5, P_i = theta + eta_i), unlike the
    # exponential PK Equation 1, so the eta must sit on the natural scale
    # of beats/min -- and Emax is negative, which a log transform cannot
    # represent at all. The Table 4 magnitudes confirm the scale and the
    # Discussion corroborates them: omega^2_E0 = 1380 is an SD of 37.1
    # beats/min, and "the inter animal variability in baseline heart rate
    # approximates the maximal drug effect" (43 beats/min).
    # ==================================================================
    e0    <- 362        ; label("Baseline resting heart rate E0 (beats/min)")                                            # Table 4 non-ISO E0 = 362 beats/min (RSE 5.5%); Results text 362 +/- 20
    emax  <- -43.0      ; label("Maximum S(-)-atenolol effect on heart rate, Emax (beats/min; negative = reduction)")    # Table 4 non-ISO Emax = -43.0 beats/min (RSE 41.2%); Results text reports the Emax as a reduction of 43 +/- 17.7 b.p.m.
    lec50 <- log(27.9)  ; label("Plasma concentration at half-maximal effect, EC50 (ng/mL)")                             # Table 4 non-ISO EC50 = 27.9 ng/mL (RSE 97.1%); Results text 28 +/- 27
    lhill <- fixed(log(1)) ; label("Hill coefficient n of the S(-)-atenolol sigmoid (unitless)")                         # Table 4 non-ISO n = NE (fixed at n = 1) -- not estimable in this arm

    # ------------------------------------------------------------------
    # PD inter-individual variability, ADDITIVE in beats/min per
    # Equation 5. "Inter individual variability was identified for E0 and
    # E max for this group", so EC50 carries no eta. The two etas are
    # reported as separate diagonal entries with no covariance row, so
    # they are encoded uncorrelated.
    # ------------------------------------------------------------------
    etae0   ~ 1380      # Table 4 non-ISO omega^2_E0 = 1380 (beats/min)^2 (RSE 55.4%) -> SD 37.1 beats/min
    etaemax ~ 913       # Table 4 non-ISO omega^2_Emax = 913 (beats/min)^2 (RSE 55.9%) -> SD 30.2 beats/min

    # ------------------------------------------------------------------
    # PD residual error, ADDITIVE in beats/min. Table 4 reports the
    # VARIANCE sigma^2 = 250, so the SD is 15.8 beats/min.
    #
    # The paper's Equation 6 prints a PROPORTIONAL form for the PD
    # residual, but that equation is a verbatim copy-paste of the PK
    # Equation 2 -- its own explanatory sentence still reads "in which
    # C_obs,ij is the jth observed CONCENTRATION", inside the
    # Pharmacodynamics section. The reported magnitudes settle it:
    # sigma^2 = 250 read proportionally would be a residual CV of 1581%,
    # whereas read additively it is 15.8 beats/min on a 319-362
    # beats/min signal (4-5%), matching a continuously recorded
    # heart-rate trace. All three of this paper's PD residual variances
    # (250 here, 409 in Table 2 and 747 in the ISO column) behave the
    # same way. Encoded additively; recorded in the vignette Errata.
    # ------------------------------------------------------------------
    addSd_hr <- sqrt(250) ; label("Additive residual error on heart rate (beats/min)")   # Table 4 non-ISO sigma^2 = 250 (beats/min)^2 (RSE 23.2%)
  })

  model({
    # ---- 1. Individual parameters ----
    # PK: exponential etas (Equation 1). PD: additive etas (Equation 5).
    cl   <- exp(lcl  + etalcl)
    vc   <- exp(lvc)
    q    <- exp(lq)
    vp   <- exp(lvp  + etalvp)
    q2   <- exp(lq2)
    vp2  <- exp(lvp2 + etalvp2)

    ec50 <- exp(lec50)
    hill <- exp(lhill)

    e0_i   <- e0   + etae0
    emax_i <- emax + etaemax

    # ---- 2. Micro-constants ----
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # ---- 3. Three-compartment PK ----
    # Amount form, so rxode2's own zero-order infusion machinery supplies
    # the infusion rate; dose `central` with dur = 15 min.
    d/dt(central)     <- -(kel + k12 + k13) * central + k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1) <-   k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-   k13 * central - k31 * peripheral2

    Cc <- central / vc

    # ---- 4. Direct-effect Emax on the PLASMA concentration ----
    # No effect compartment: ke0 is "NE" in this arm of Table 4. The
    # `^hill` power is retained with n fixed at 1 so the structure stays
    # the paper's Equation 3 rather than being silently collapsed.
    hr <- e0_i + emax_i * Cc^hill / (ec50^hill + Cc^hill)

    # ---- 5. Observations ----
    Cc ~ prop(propSd)
    hr ~ add(addSd_hr)
  })
}
