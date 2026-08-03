Jensen_2023_lngIus52mg <- function() {
  description <- paste(
    "Eight-year population PK / SHBG-turnover model for levonorgestrel (LNG)",
    "released from the LNG-releasing intrauterine system (LNG-IUS) 52 mg",
    "(Mirena). In-vivo release from the reservoir is described as a mixture",
    "of a first-order term (rate coefficient c12), a constant zero-order term",
    "(c13), and a time-dependent first-order term with time-decay constant t1,",
    "giving DADT(depot) = -c12 * depot - c13 * (1 + depot / (t1 + t)). Two-",
    "compartment LNG disposition treats only unbound drug as distributable /",
    "eliminable via K20 * fuLNG * central and K23 * fuLNG * central, where",
    "fuLNG is the closed-form free-fraction solution to reversible LNG binding",
    "to albumin (constant, KDA = 18209 nmol/L, ALB = 700000 nmol/L) and to",
    "SHBG (KDS = 1.82 nmol/L). SHBG serum concentration is modelled with an",
    "indirect-response turnover (zero-order kin, first-order kout) whose",
    "synthesis is linearly inhibited by a delay-compartment-smoothed unbound",
    "LNG signal (delay time-constant tau). Body weight enters as an",
    "allometric-style power on apparent (free) LNG clearance and on the SHBG",
    "baseline. Bioavailability of the loaded LNG reservoir is fixed via a",
    "logit anchor (F1 ~ 0.971 of the 52 mg reservoir). All rate constants are",
    "expressed per hour."
  )
  reference <- paste(
    "Jensen JT, Reinecke I, Post TM, Lukkari-Lax E, Hofmann BM.",
    "Extended use of levonorgestrel-releasing intrauterine system (LNG-IUS)",
    "52 mg: A population pharmacokinetic approach to estimate in vivo",
    "levonorgestrel release rates and systemic exposure including comparison",
    "with two other LNG-IUSs.",
    "Contraception. 2023 May;121:109954.",
    "doi:10.1016/j.contraception.2023.109954.",
    "Upstream 5-year integrated LNG contraceptive popPK meta-analysis:",
    "Reinecke I, Hofmann B, Mesic E, Drenth HJ, Garmann D.",
    "J Clin Pharmacol. 2018 Dec;58(12):1639-1654. doi:10.1002/jcph.1288."
  )
  vignette <- "Jensen_2023_lngIus52mg"
  units <- list(
    time          = "h",
    dosing        = "mg (levonorgestrel loaded in the intrauterine reservoir)",
    concentration = "ng/L (total LNG in plasma; also derived unbound LNG in ng/L) and nmol/L (SHBG in serum)"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight (kg). Continuous covariate on apparent clearance of LNG and on SHBG baseline via power-law scaling.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Reference weight is 68 kg (medWT), the pooled median across the",
        "8-year popPK dataset (median weight column in Supplementary Table 1:",
        "70 kg MET, 67 kg Phase 2 Study 308901, 62 kg Phase 3 Study 89532;",
        "NONMEM control stream fixes medWT = 68). Missing / zero WT is imputed",
        "to 68 in the source NONMEM code (`IF (WGHTC.EQ.-999.OR.WGHTC.EQ.0)",
        "WT=medWT`) so the packaged model expects WT populated; downstream",
        "simulation callers should substitute the reference 68 kg when WT is",
        "unknown for a given subject. Body weight is the only covariate the",
        "8-year model retained (Supplementary Section 1.5)."
      ),
      source_name        = "WGHTC"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 920L,
    n_studies      = 3L,
    studies        = c(
      "NCT02985541 (Mirena Extension Trial, MET, Phase 3, US 54 centers, 2016-2021)",
      "NCT00185380 (Phase 2 Study 308901, LNG-IUS 52 mg, up to ~3 y)",
      "Leiras Study Report 02-89532-07 (Phase 3 Study 89532, LNG-IUS 52 mg, up to ~5 y)"
    ),
    age_range      = "18-40 years (medians 29-33 across the three studies)",
    weight_range   = "39-164 kg (medians 62 / 67 / 70 kg across studies; pooled median 68 kg)",
    sex_female_pct = 100,
    race_ethnicity = c(White = 75.4, Black = 14.1, Asian = 2.5, AmericanIndianAlaskaNative = 0.6, Multiple = 3.9, NotReported = 3.6),
    disease_state  = "Premenopausal fertile women using LNG-IUS 52 mg (Mirena) for contraception, with or without heavy menstrual bleeding.",
    dose_range     = "Single insertion of a 52 mg LNG intrauterine reservoir (F1 ~ 0.971 of the loaded amount enters the release pool); no additional dosing over 8 years of continuous use.",
    regions        = "United States (MET); global for the pooled Phase 2 / 3 upstream studies.",
    notes          = paste(
      "Race percentages are baseline demographics for the 362 MET subjects",
      "starting extended treatment in year 6 (Supplementary Section 2.1);",
      "pooled race distribution across all 3 studies not reported. Total",
      "observation counts included in the 8-year popPK fit: 1806 plasma LNG",
      "concentrations, 1567 serum SHBG concentrations, and 262 residual-",
      "content measurements. LNG was assayed by validated LC-MS/MS (LLOQ",
      "10 ng/L); SHBG by dissociation-enhanced lanthanide fluoroimmunoassay",
      "(LLOQ 10 nmol/L); residual LNG content by LC after dichloromethane",
      "extraction (Supplementary Section 1.3)."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Release parameters. The paper's release-rate formulation combines a
    # first-order term (c12 * A_depot), a zero-order term (c13), and a
    # time-dependent first-order term (c13 * A_depot / (t1 + t)). Both c12
    # and c13 carry an explicit 10^-6 scale in the NONMEM control stream
    # (Supplementary Appendix, $PK: TC12 = THETA(1)*10**-6). The scale is
    # applied inside model() so that the ini() typical values match the
    # NONMEM $THETA numbers verbatim.
    # ------------------------------------------------------------------
    lc12   <- log(7.98)          ; label("First release rate coefficient c12 (x 1e-6 per hour); THETA(1) in Jensen 2023 Supp Appendix")   # Supp Appendix $THETA TH1 = 7.98
    lc13   <- fixed(log(303))    ; label("Second release rate coefficient c13 (x 1e-6 mg/h); THETA(11) FIX in Jensen 2023 Supp Appendix")  # Supp Appendix $THETA TH11 = 303 FIX
    lt1    <- log(58.8)          ; label("Release time-dependency parameter t1 (hour); THETA(12) in Jensen 2023 Supp Appendix")               # Supp Appendix $THETA TH12 = 58.8

    # ------------------------------------------------------------------
    # LNG disposition. V2, V3, and Q3 are all FIXED in the NONMEM control
    # stream (upstream 5-year integrated popPK values from Reinecke 2018);
    # only apparent free clearance is re-estimated. Note that "clearance"
    # here is the apparent free-drug clearance CL/F -- the elimination in
    # the ODE is CL * fuLNG * (A / V) so total-conc-based clearance is
    # CL * fuLNG ~ 242 * 0.016 = 3.87 L/h.
    # ------------------------------------------------------------------
    lvc    <- fixed(log(20.7))   ; label("Apparent central volume V2 (L); THETA(2) FIX in Jensen 2023 Supp Appendix")     # Supp Appendix $THETA TH2 = 20.7 FIX
    lcl    <- log(242)           ; label("Apparent (free) clearance CL/F of LNG (L/h); THETA(3) in Jensen 2023 Supp Appendix")     # Supp Appendix $THETA TH3 = 242
    lvp    <- fixed(log(4690))   ; label("Apparent peripheral volume V3 (L); THETA(4) FIX in Jensen 2023 Supp Appendix")  # Supp Appendix $THETA TH4 = 4690 FIX
    lq     <- fixed(log(600))    ; label("Apparent intercompartmental clearance Q3 (L/h); THETA(5) FIX in Jensen 2023 Supp Appendix")  # Supp Appendix $THETA TH5 = 600 FIX

    # ------------------------------------------------------------------
    # SHBG delayed-inhibition parameters. tau is the effect-compartment
    # equilibration time-constant (delay for the LNG signal driving SHBG
    # inhibition); ri is the linear inhibition slope.
    # ------------------------------------------------------------------
    ltau   <- fixed(log(13.7))   ; label("Delay time-constant tau for LNG effect on SHBG (hour); THETA(6) FIX in Jensen 2023 Supp Appendix")  # Supp Appendix $THETA TH6 = 13.7 FIX
    lri    <- fixed(log(0.232))  ; label("Linear inhibition slope ri of SHBG synthesis by delayed LNG (L/nmol); THETA(7) FIX in Jensen 2023 Supp Appendix")  # Supp Appendix $THETA TH7 = 0.232 FIX

    # ------------------------------------------------------------------
    # SHBG turnover. SBL is the drug-free baseline SHBG concentration (kin
    # is derived as SBL * kout so that the SHBG steady state without LNG
    # equals SBL). kout is FIXED to the upstream Reinecke 2018 estimate.
    # ------------------------------------------------------------------
    lrbase_shbg <- log(51.8)             ; label("SHBG baseline concentration SBL (nmol/L); THETA(8) in Jensen 2023 Supp Appendix")   # Supp Appendix $THETA TH8 = 51.8
    lkout_shbg  <- fixed(log(0.00313))   ; label("SHBG first-order elimination rate constant kout (1/h); THETA(9) FIX in Jensen 2023 Supp Appendix")  # Supp Appendix $THETA TH9 = 0.00313 FIX

    # ------------------------------------------------------------------
    # Bioavailability of the loaded reservoir. F1 = plogis(logitfdepot) so
    # logitfdepot = 3.51 gives F1 = exp(3.51)/(1+exp(3.51)) = 0.9710. NONMEM
    # code: F1 = EXP(TF1)/(1+EXP(TF1)); THETA(10) FIX at 3.51.
    # ------------------------------------------------------------------
    logitfdepot <- fixed(3.51)   ; label("Logit-transformed bioavailability of LNG-IUS 52 mg reservoir into depot; THETA(10) FIX in Jensen 2023 Supp Appendix")  # Supp Appendix $THETA TH10 = 3.51 FIX

    # ------------------------------------------------------------------
    # Body-weight effect. Power-law scaling with reference weight medWT =
    # 68 kg for both apparent clearance and SHBG baseline. NONMEM:
    # CO1 = (WT/68)^THETA(13); CO2 = (WT/68)^THETA(14).
    # ------------------------------------------------------------------
    e_wt_cl          <- 0.728    ; label("Body-weight power exponent on apparent clearance of LNG (WT / 68 kg)^e_wt_cl; THETA(13) in Jensen 2023 Supp Appendix")           # Supp Appendix $THETA TH13 = 0.728
    e_wt_rbase_shbg  <- -0.987   ; label("Body-weight power exponent on SHBG baseline (WT / 68 kg)^e_wt_rbase_shbg; THETA(14) in Jensen 2023 Supp Appendix")               # Supp Appendix $THETA TH14 = -0.987

    # ------------------------------------------------------------------
    # Inter-individual variability. NONMEM $OMEGA BLOCK(2) reports
    # variance-covariance entries verbatim (log-normal exponential-eta
    # model on both CL and SBL). CL IIV variance 0.0345 corresponds to
    # sqrt(exp(0.0345) - 1) ~ 18.7% CV; SHBG-baseline IIV variance 0.173
    # corresponds to sqrt(exp(0.173) - 1) ~ 43.5% CV. Off-diagonal
    # covariance -0.0505 implies a correlation of
    # -0.0505 / sqrt(0.0345 * 0.173) = -0.654 (strong negative correlation
    # of individual CL and SBL, matching the paper's Discussion: 'the
    # positive correlation with LNG clearance and the negative correlation
    # with SHBG baseline concentration').
    # ------------------------------------------------------------------
    etalcl + etalrbase_shbg ~ c(0.0345,
                                -0.0505, 0.173)

    # ------------------------------------------------------------------
    # Residual error. NONMEM $SIGMA has proportional error on LNG and on
    # SHBG (Y = IPRED + W*EPS where W = IPRED for CMT 3/6), and pure
    # additive error on the residual-content observation (Y = IPRED +
    # EPS for CMT 7). Values are variances in NONMEM; the standard
    # deviations used by nlmixr2 are sqrt(variance).
    # ------------------------------------------------------------------
    propSd            <- sqrt(0.0329)  ; label("Proportional residual SD on plasma LNG concentration (fraction of Cc); $SIGMA(1) in Jensen 2023 Supp Appendix, sqrt(0.0329) ~= 0.181")   # Supp Appendix $SIGMA(1) = 0.0329
    propSd_shbg       <- sqrt(0.0344)  ; label("Proportional residual SD on serum SHBG concentration (fraction of shbg); $SIGMA(2) in Jensen 2023 Supp Appendix, sqrt(0.0344) ~= 0.186")  # Supp Appendix $SIGMA(2) = 0.0344
    addSd_iusResidual <- sqrt(1.52)    ; label("Additive residual SD on residual LNG content in explanted device (mg); $SIGMA(3) in Jensen 2023 Supp Appendix, sqrt(1.52) ~= 1.233")     # Supp Appendix $SIGMA(3) = 1.52
  })

  model({
    # ------------------------------------------------------------------
    # Physicochemical constants for the closed-form free-fraction of LNG.
    # Values are the NONMEM $PK block constants; MW LNG = 312.5 g/mol,
    # albumin dissociation constant KDA = 18209 nmol/L, albumin
    # concentration ALB = 700000 nmol/L, SHBG dissociation constant
    # KDS = 1.82 nmol/L (Supplementary Table, "Free fraction of LNG").
    # ------------------------------------------------------------------
    MWLNG <- 312.5
    KDS   <- 1.82
    KDA   <- 18209
    ALB   <- 700000
    KALB  <- 1 + ALB / KDA

    # Reference weight (median across the pooled 8-year popPK dataset).
    wtRef <- 68

    # ------------------------------------------------------------------
    # Structural release parameters (with the 10^-6 scale factor applied
    # here so the ini() values match the NONMEM $THETA numbers verbatim).
    # ------------------------------------------------------------------
    c12   <- exp(lc12) * 1e-6
    c13   <- exp(lc13) * 1e-6
    t1    <- exp(lt1)

    # LNG disposition individual parameters.
    vc <- exp(lvc)
    cl <- exp(lcl + e_wt_cl * log(WT / wtRef) + etalcl)
    vp <- exp(lvp)
    q  <- exp(lq)
    k20 <- cl / vc
    k23 <- q  / vc
    k32 <- q  / vp

    # SHBG parameters.
    tau       <- exp(ltau)
    ri        <- exp(lri)
    rbase_shbg <- exp(lrbase_shbg + e_wt_rbase_shbg * log(WT / wtRef) + etalrbase_shbg)
    kout_shbg <- exp(lkout_shbg)
    kin_shbg  <- rbase_shbg * kout_shbg

    # Bioavailability of the LNG reservoir into the depot compartment.
    fdepot <- exp(logitfdepot) / (1 + exp(logitfdepot))

    # ------------------------------------------------------------------
    # Free fraction of LNG (closed-form solution to reversible binding).
    # The molar concentration in the central compartment is
    # A3nM = (central [mg] / vc [L]) * (10^6 [ng/mg]) / MWLNG [g/mol] =
    # (central / vc) * (10^6 / MWLNG) [nmol/L]. The +1e-6 guard mirrors
    # the NONMEM code and keeps A3nM strictly positive early in the run
    # when the depot has released essentially nothing.
    # ------------------------------------------------------------------
    A3nM <- (central / vc) * (1e6 / MWLNG) + 1e-6
    temp1 <- KALB * KDS + shbg - A3nM
    temp2 <- temp1 / (KALB * 2)
    temp3 <- temp2^2 + (KDS / KALB) * A3nM
    fuLNG <- (1 / A3nM) * (-temp2 + sqrt(temp3))

    # Inhibition of SHBG synthesis by delayed LNG (capped at 1 as in the
    # source NONMEM code).
    inh_raw <- ri * effect
    inh     <- (inh_raw <= 1) * inh_raw + (inh_raw > 1) * 1

    # Release inputs from the IUS device to central circulation.
    input1 <- c12 * depot
    input2 <- c13 * (1 + depot / (t1 + t + 1e-6))

    # ------------------------------------------------------------------
    # ODE system (Supplementary Table, "The 8-year popPK model"):
    #   depot        : LNG remaining in the IUS device (mg)
    #   central      : LNG in the plasma central compartment (mg)
    #   peripheral1  : LNG in the peripheral tissue compartment (mg)
    #   effect       : delayed molar LNG signal driving SHBG inhibition
    #                  (nmol/L; equivalent to Sheiner effect compartment)
    #   shbg         : SHBG in serum (nmol/L)
    # ------------------------------------------------------------------
    d/dt(depot)       <- -input1 - input2
    d/dt(central)     <-  input1 + input2 - (k20 + k23) * fuLNG * central + k32 * peripheral1
    d/dt(peripheral1) <-  k23 * fuLNG * central - k32 * peripheral1
    d/dt(effect)      <- (1 / tau) * (A3nM - effect)
    d/dt(shbg)        <-  kin_shbg * (1 - inh) - kout_shbg * shbg

    f(depot) <- fdepot
    shbg(0)  <- rbase_shbg

    # ------------------------------------------------------------------
    # Observations. Total LNG plasma concentration Cc is central amount
    # (mg) rescaled to ng/L via 10^6 factor. Unbound LNG concentration
    # is derived for reporting only (no separate residual error attached).
    # Residual LNG content in the device is simply the depot amount (mg).
    # ------------------------------------------------------------------
    Cc          <- (central / vc) * 1e6
    CcUnbound   <- fuLNG * Cc
    iusResidual <- depot

    Cc          ~ prop(propSd)
    shbg        ~ prop(propSd_shbg)
    iusResidual ~ add(addSd_iusResidual)
  })
}
