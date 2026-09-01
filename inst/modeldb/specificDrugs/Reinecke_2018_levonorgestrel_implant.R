Reinecke_2018_levonorgestrel_implant <- function() {
  description <- paste(
    "Final subdermal-implant model of the Reinecke 2018 integrated",
    "levonorgestrel (LNG) population PK analysis, for Jadelle (2 rods of",
    "75 mg LNG each). In-vivo release from the implant uses the same two",
    "processes as LNG-IUS 12 and LNG-IUS 8: a zero-order term c12 plus a",
    "time-dependent first-order term c13 * depot / (t1 + t). Released LNG",
    "enters the comprehensive-model two-compartment disposition in which",
    "only unbound drug is eliminated or distributed, the free fraction fuLNG",
    "being the closed-form solution of reversible binding to SHBG",
    "(KDS = 1.82 nmol/L) and albumin (KDA = 18209 nmol/L, constant",
    "ALB = 700000 nmol/L). SHBG is retained as a latent turnover state",
    "because it drives the free fraction, but no SHBG data were available",
    "for the implant studies, so its baseline is fixed to the comprehensive",
    "model's intrauterine value and it carries neither inter-individual",
    "variability nor a residual error. Only release parameters,",
    "bioavailability, inter-individual variability on clearance and the two",
    "residual errors were estimated; everything else is fixed to the final",
    "comprehensive model."
  )
  reference <- paste(
    "Reinecke I, Hofmann B, Mesic E, Drenth HJ, Garmann D.",
    "An integrated population pharmacokinetic analysis to characterize",
    "levonorgestrel pharmacokinetics after different administration routes.",
    "J Clin Pharmacol. 2018 Dec;58(12):1639-1654. doi:10.1002/jcph.1288.",
    "Parameter values from Supplemental Table S2, 'final implant model'",
    "block; model equations from Supplemental Table S3a. Parameters not",
    "listed in that block are fixed to the final comprehensive model,",
    "available as modellib('Reinecke_2018_levonorgestrel_lngIus20') and its",
    "LNG-IUS 12 / LNG-IUS 8 companions."
  )
  vignette <- "Reinecke_2018_levonorgestrel_contraceptives"
  units <- list(
    time          = "h",
    dosing        = "mg (levonorgestrel loaded in the subdermal implant)",
    concentration = "ng/L (total and unbound LNG in serum); SHBG in nmol/L; residual implant content in mg"
  )

  compartmentData <- list(
    depot       = list(analyte = "levonorgestrel", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "levonorgestrel", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "levonorgestrel", units = "mg", specimen = "serum", verified = TRUE),
    effect      = list(analyte = "levonorgestrel", units = "nmol/L", specimen = "not applicable", verified = TRUE),
    shbg        = list(analyte = "SHBG", units = "nmol/L", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight (kg). Power-function covariate on clearance of unbound levonorgestrel (increasing with weight) and on the SHBG baseline concentration (decreasing with weight). Both exponents are fixed to the final comprehensive model.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline body weight, centred on 66 kg. The Implant Model section",
        "states this explicitly: 'Because only partial information on body",
        "weight was available, the covariate effect implemented in the",
        "comprehensive model was applied to the implant population (centered",
        "around the median used in the comprehensive model, 66 kg). Missing",
        "body weights were replaced by the median in the comprehensive",
        "model.' The implant cohort's own median weight is 58 kg (Table 1),",
        "which is NOT the centering value. No covariate analysis was",
        "performed in the implant step because no age data were available",
        "and weight data were incomplete."
      ),
      source_name        = "WGHT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 471L,
    n_studies      = 3L,
    studies        = "Phase 3 studies 39, 41 and 42 (subdermal implant Jadelle, comprehensive records of levonorgestrel concentrations and residual content over 5 years of intended treatment)",
    weight_range   = "39-87 kg (median 58 kg, among the 197 subjects with concentration data)",
    sex_female_pct = 100,
    disease_state  = "Healthy premenopausal women using a levonorgestrel-releasing subdermal implant for contraception",
    dose_range     = "Single subdermal placement of 2 rods containing 75 mg levonorgestrel each (150 mg total), indicated for up to 5 years",
    notes          = paste(
      "Table 1 splits the implant cohort into 274 subjects contributing 275",
      "residual-content measurements and a separate 197 subjects",
      "contributing 1634 levonorgestrel concentrations; the two sets are",
      "disjoint, which the Discussion flags as a limitation on the",
      "bioavailability estimate. No age data and only partial weight data",
      "were available, so no covariate analysis was performed in this step.",
      "No SHBG concentrations were measured. Nine outliers (|weighted",
      "residual| > 6) were excluded during model development. Epsilon",
      "shrinkage was below 5% for both outputs but eta shrinkage on",
      "clearance was high (35.8%), so individual clearance estimates from",
      "this model should be treated with caution."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # In-vivo release from the implant, using the same branch of Table S3a
    # as LNG-IUS 12 and LNG-IUS 8:
    #   input1 = C12
    #   input2 = C13 * A(2) / (T1 + T)
    # Table S2 reports IC12 in ng/h; the 1e-6 conversion to the mg/h scale
    # of the device compartment is applied inside model(). IC13 is
    # dimensionless.
    # ------------------------------------------------------------------
    lc12 <- log(1320)    ; label("Zero-order release rate IC12 (ng/h)")                          # Table S2, final implant model, IC12 = 1320 (1260-1380) ng/h
    lc13 <- log(0.0246)  ; label("Second release process coefficient IC13 (dimensionless)")      # Table S2, final implant model, IC13 = 0.0246 (0.0145-0.0347)
    lt1  <- log(3640)    ; label("Time parameter IT1 for the drop in release process 2 (h)")     # Table S2, final implant model, IT1 = 3640 (1990-5290)

    # Absolute bioavailability on the logit scale: TF = 0.664 gives
    # F = 0.660, matching the value Table S2 prints in parentheses. The
    # Discussion notes this is lower than expected and may be biased by the
    # residual-content and concentration data coming from different
    # subjects.
    logitfdepot <- 0.664 ; label("Logit-transformed absolute bioavailability of the implant reservoir (logit units)")  # Table S2, final implant model, F_Implant = 0.664 (0.521-0.807), i.e. 0.660 (0.627-0.691)

    # ------------------------------------------------------------------
    # Everything below is fixed to the final comprehensive model
    # (Table S2 footnotes g and i).
    # ------------------------------------------------------------------
    lvc <- fixed(log(20.7))  ; label("Central volume of distribution Vc (L)")                          # Table S2, final comprehensive model, Vc = 20.7; fixed here (footnote i)
    lcl <- fixed(log(243))   ; label("Clearance of unbound levonorgestrel CL (L/h)")                   # Table S2, final comprehensive model, CL = 243; fixed here (footnote i)
    lvp <- fixed(log(4690))  ; label("Peripheral volume of distribution Vp (L)")                       # Table S2, final comprehensive model, Vp = 4690; fixed here (footnote i)
    lq  <- fixed(log(600))   ; label("Intercompartmental clearance of unbound levonorgestrel Q (L/h)") # Table S2, final comprehensive model, Q = 600; fixed here (footnote i)

    lrbase_shbg <- fixed(log(51.5))    ; label("Baseline SHBG serum concentration SBL at 66 kg (nmol/L)")            # Table S2, final implant model, SBL = 51.5, fixed to the comprehensive SBL_IUS (footnote g)
    lkout_shbg  <- fixed(log(0.00313)) ; label("First-order elimination rate constant of SHBG KOUT (1/h)")           # Table S2, final comprehensive model, KOUT = 0.00313; fixed here (footnote i)
    ltau        <- fixed(log(13.7))    ; label("Delay time constant tau of the levonorgestrel effect on SHBG (h)")   # Table S2, final comprehensive model, tau = 13.7; fixed here (footnote i)
    lri         <- fixed(log(0.232))   ; label("Linear inhibition factor RI of SHBG synthesis by delayed LNG (L/nmol)") # Table S2, final comprehensive model, R_I = 0.232; fixed here (footnote i)

    e_wt_cl         <- fixed( 0.823)   ; label("Body-weight power exponent on clearance, (WT/66)^e_wt_cl")                 # Table S2, final implant model, C_WGHT,CL = 0.823, fixed to the comprehensive value (footnote g)
    e_wt_rbase_shbg <- fixed(-0.950)   ; label("Body-weight power exponent on the SHBG baseline, (WT/66)^e_wt_rbase_shbg") # Table S2, final implant model, C_WGHT,SBL = -0.950, fixed to the comprehensive value (footnote g)

    # ------------------------------------------------------------------
    # Inter-individual variability. Only clearance carries IIV here: the
    # SHBG-baseline variance and the covariance were both fixed to zero
    # because no SHBG data were available (Table S2 footnote j), so the
    # comprehensive model's two-by-two block collapses to a single eta.
    # log(1 + 0.351^2) = 0.116183.
    # ------------------------------------------------------------------
    etalcl ~ 0.116183   # Table S2, final implant model, IIV_CL 35.1% (30.2-39.5); IIV_SBL and COV_IIV fixed to 0 (footnote j)

    # ------------------------------------------------------------------
    # Residual error. No proportional error was estimated for SHBG (fixed
    # to zero, footnote j), so SHBG is a latent state here and carries no
    # endpoint. The additive residual-content entry is on the variance
    # scale, matching the other Table S2 additive rows; sqrt(19.1) = 4.37
    # mg on a 150 mg implant.
    # ------------------------------------------------------------------
    propSd                <- 0.205         ; label("Proportional residual error SD on total levonorgestrel concentration (fraction)")  # Table S2, final implant model, SIGMA_LNG (CV) = 20.5% (19.4-21.6)
    addSd_implantResidual <- sqrt(19.1)    ; label("Additive residual error SD on residual implant content (mg)")                      # Table S2, final implant model, SIGMA Implant res. content = 19.1 (15.7-22.5) as a variance; sqrt = 4.37 mg
  })

  model({
    # Physicochemical constants for the closed-form free fraction of LNG
    # (Table S3a, "Free fraction of LNG").
    MWLNG <- 312.5      # Table S3a: molecular weight LNG, g/mol
    KDS   <- 1.82       # Table S3a: dissociation constant SHBG, nmol/L
    KDA   <- 18209      # Table S3a: dissociation constant albumin, nmol/L
    ALB   <- 700000     # Table S3a: albumin concentration, nmol/L
    KALB  <- 1 + ALB / KDA

    wtRef <- 66

    # Release parameters. c12 converts ng/h to the mg/h scale of the
    # implant compartment; c13 is dimensionless.
    c12 <- exp(lc12) * 1e-6
    c13 <- exp(lc13)
    t1  <- exp(lt1)

    vc <- exp(lvc)
    cl <- exp(lcl + e_wt_cl * log(WT / wtRef) + etalcl)
    vp <- exp(lvp)
    q  <- exp(lq)

    k20 <- cl / vc
    k23 <- q / vc
    k32 <- q / vp

    tau        <- exp(ltau)
    ri         <- exp(lri)
    rbase_shbg <- exp(lrbase_shbg + e_wt_rbase_shbg * log(WT / wtRef))
    kout_shbg  <- exp(lkout_shbg)
    kin_shbg   <- rbase_shbg * kout_shbg

    fdepot <- exp(logitfdepot) / (1 + exp(logitfdepot))

    A3nM  <- (central / vc) * (1e6 / MWLNG) + 1e-6
    temp1 <- KALB * KDS + shbg - A3nM
    temp2 <- temp1 / (KALB * 2)
    temp3 <- temp2^2 + (KDS / KALB) * A3nM
    fuLNG <- (1 / A3nM) * (-temp2 + sqrt(temp3))

    inh <- ri * effect

    # ------------------------------------------------------------------
    # Implant release branch of Table S3a. Note where the bioavailability
    # enters: for the implant it scales the flux ARRIVING in the central
    # compartment, while the implant compartment itself depletes at the
    # full release rate. This is the difference Table S2 footnote e refers
    # to ("implementation for IUS differs from implementation for the
    # implant for numerical reasons"), and it is pinned by the paper's own
    # numbers -- the implant is loaded with 150 mg (2 rods of 75 mg),
    # Supplemental Figure S1d shows the residual-content fit starting at
    # 150 mg rather than at 0.660 * 150 = 99 mg, and only the full 150 mg
    # reproduces the Table 4 release rate of 53.0 ug/day at 24 days
    # (1320e-6 + 0.0246 * 150 / (3640 + 576) = 2.19e-3 mg/h = 52.7
    # ug/day). Applying F to the dose instead would give 45.5 ug/day.
    # ------------------------------------------------------------------
    input1 <- c12
    input2 <- c13 * depot / (t1 + t)

    d/dt(depot)       <- -input1 - input2
    d/dt(central)     <-  fdepot * (input1 + input2) - (k20 + k23) * fuLNG * central + k32 * peripheral1
    d/dt(peripheral1) <-  k23 * fuLNG * central - k32 * peripheral1
    d/dt(effect)      <- (1 / tau) * A3nM - (1 / tau) * effect
    d/dt(shbg)        <-  kin_shbg * (1 - inh) - kout_shbg * shbg

    shbg(0) <- rbase_shbg

    Cc              <- (central / vc) * 1e6
    CcUnbound       <- fuLNG * Cc
    implantResidual <- depot
    releaseRate     <- (input1 + input2) * 1000 * 24

    Cc              ~ prop(propSd)
    implantResidual ~ add(addSd_implantResidual)
  })
}
