Reinecke_2018_levonorgestrel_progestinOnlyPill <- function() {
  description <- paste(
    "Oral arm of the Reinecke 2018 final COMPREHENSIVE levonorgestrel (LNG)",
    "population PK model, as applied to the progestin-only pill",
    "Microlut/Norgeston (30 ug LNG once daily). Structurally identical to the",
    "final iv-oral model -- two-compartment LNG disposition in which only",
    "unbound drug is eliminated or distributed, the free fraction fuLNG being",
    "the closed-form solution of reversible binding to SHBG (KDS = 1.82",
    "nmol/L) and albumin (KDA = 18209 nmol/L, constant ALB = 700000 nmol/L),",
    "coupled to an SHBG indirect-response turnover state inhibited by a",
    "delay-compartment-smoothed molar LNG signal -- but carrying the",
    "comprehensive fit's parameter values, including the separately estimated",
    "iv/oral SHBG baseline of 76.4 nmol/L and the body-weight effect on",
    "clearance that only the intrauterine data could identify. This is the",
    "parameterisation the paper used to derive the progestin-only-pill",
    "exposure statistics in Table 3. Use",
    "modellib('Reinecke_2018_levonorgestrel_ivOral') for the final iv-oral",
    "model itself."
  )
  reference <- paste(
    "Reinecke I, Hofmann B, Mesic E, Drenth HJ, Garmann D.",
    "An integrated population pharmacokinetic analysis to characterize",
    "levonorgestrel pharmacokinetics after different administration routes.",
    "J Clin Pharmacol. 2018 Dec;58(12):1639-1654. doi:10.1002/jcph.1288.",
    "Parameter values from Supplemental Table S2, 'final comprehensive",
    "model' block; model equations from Supplemental Table S3a. The",
    "intrauterine arms of the same fit are",
    "modellib('Reinecke_2018_levonorgestrel_lngIus20') and its LNG-IUS 12 /",
    "LNG-IUS 8 companions."
  )
  vignette <- "Reinecke_2018_levonorgestrel_contraceptives"
  units <- list(
    time          = "h",
    dosing        = "mg (levonorgestrel)",
    concentration = "ng/L (total and unbound LNG in serum/plasma); SHBG in nmol/L"
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
      description        = "Body weight (kg). Power-function ('allometric') covariate on the absolute oral bioavailability of levonorgestrel (decreasing with weight) and on the clearance of unbound levonorgestrel (increasing with weight).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline body weight, centred on 66 kg -- the pooled",
        "intrauterine-system median that the comprehensive fit uses",
        "throughout. The comprehensive F_oral row of Supplemental Table S2 is",
        "LABELLED 'median[WGHT]=65kg', but that label is inconsistent with",
        "the derived bioavailabilities printed in the same row block: 0.751 *",
        "(62/65)^-1.09 = 0.791, whereas the table reports 0.804 at 62 kg,",
        "0.798 at 62.4 kg and 0.818 at 61 kg. Those three values are all",
        "reproduced to three decimal places by 0.751 * (WT/66)^-1.09, so 66",
        "kg is the reference weight and '65kg' is a typographical error. 66",
        "kg is also the value the F_LNG-IUS rows carry ('median[WGHT]=66kg,",
        "IUS data') and the value the Implant Model section names as 'the",
        "median used in the comprehensive model'. The body-weight effect on",
        "the SHBG baseline was identified on, and applies only to, the",
        "intrauterine baseline SBL_IUS, so it is not carried here."
      ),
      source_name        = "WGHT"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline (years). Screened as a power-function covariate on levonorgestrel clearance and on the SHBG baseline concentration.",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened in the comprehensive covariate step but not retained; only",
        "body weight reached the P <= .001 forward-inclusion criterion",
        "(Results, 'Comprehensive Model'). No point estimate is reported."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 21L,
    n_studies      = 2L,
    studies        = c(
      "Phase 1 study 15687 (progestin-only pill Microlut/Norgeston, 30 ug LNG once daily, dense LNG and SHBG profiles over 1 treatment cycle, 21 women)",
      "Phase 1 study 92085 (single-dose intravenous 90 ug LNG; included in every model-development step so that absolute bioavailability is estimable)"
    ),
    age_range      = "31-45 years (median 42 years in study 15687)",
    weight_range   = "52.5-82.9 kg (median 62.4 kg in study 15687)",
    sex_female_pct = 100,
    disease_state  = "Healthy premenopausal women",
    dose_range     = "30 ug levonorgestrel orally once daily",
    notes          = paste(
      "Demographics for the progestin-only pill from Table 1. The",
      "comprehensive fit these parameters come from was estimated jointly on",
      "intravenous, oral (progestin-only pill plus supportive single-dose",
      "data) and intrauterine data -- 2884 women across the three",
      "intrauterine devices plus the 51 women contributing oral and",
      "intravenous data -- so the disposition parameters here are identical",
      "to those of the LNG-IUS arms. Study 15687 contributed 599 LNG and 621",
      "SHBG concentrations; LNG was assayed there by LC-MS/MS. Eta and",
      "epsilon shrinkage for the comprehensive model were below 18%."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Levonorgestrel disposition, shared verbatim with the intrauterine
    # arms of the comprehensive fit. cl and q are clearances of UNBOUND
    # drug: the ODEs multiply them by the free fraction fuLNG.
    # ------------------------------------------------------------------
    lka  <- log(2.06)   ; label("Oral first-order absorption rate constant Ka (1/h)")             # Table S2, final comprehensive model, Ka = 2.06 (1.72-2.40)
    lvc  <- log(20.7)   ; label("Central volume of distribution Vc (L)")                          # Table S2, final comprehensive model, Vc = 20.7 (18.7-22.7)
    lcl  <- log(243)    ; label("Clearance of unbound levonorgestrel CL (L/h)")                   # Table S2, final comprehensive model, CL = 243 (240-246)
    lvp  <- log(4690)   ; label("Peripheral volume of distribution Vp (L)")                       # Table S2, final comprehensive model, Vp = 4690 (3910-5470)
    lq   <- log(600)    ; label("Intercompartmental clearance of unbound levonorgestrel Q (L/h)") # Table S2, final comprehensive model, Q = 600 (534-666)

    # ------------------------------------------------------------------
    # SHBG turnover. The comprehensive fit estimated the SHBG baseline
    # separately for iv/oral (76.4 nmol/L) and intrauterine (51.5 nmol/L)
    # treatment; the iv/oral value applies here. tau, RI and KOUT could not
    # be estimated once the IUS data were added and were fixed to the final
    # iv-oral estimates (Table S2 footnote f).
    # ------------------------------------------------------------------
    lrbase_shbg <- log(76.4)            ; label("Baseline SHBG serum concentration SBL_iv,oral (nmol/L)")           # Table S2, final comprehensive model, SBL_iv,oral = 76.4 (66.0-86.8)
    lkout_shbg  <- fixed(log(0.00313))  ; label("First-order elimination rate constant of SHBG KOUT (1/h)")         # Table S2, final comprehensive model, KOUT = 0.00313, fixed to the final iv-oral value (footnote f)
    ltau        <- fixed(log(13.7))     ; label("Delay time constant tau of the levonorgestrel effect on SHBG (h)") # Table S2, final comprehensive model, tau = 13.7, fixed to the final iv-oral value (footnote f)
    lri         <- fixed(log(0.232))    ; label("Linear inhibition factor RI of SHBG synthesis by delayed LNG (L/nmol)") # Table S2, final comprehensive model, R_I = 0.232, fixed to the final iv-oral value (footnote f)

    # ------------------------------------------------------------------
    # Absolute oral bioavailability at the reference weight of 66 kg, and
    # the two body-weight power exponents that apply to an oral subject.
    # The oral arm carries the body-weight effect on clearance because CL
    # is one shared parameter across the whole comprehensive fit; the
    # body-weight effect on the SHBG baseline was identified on, and
    # applies only to, SBL_IUS, so it is absent here.
    # ------------------------------------------------------------------
    lfdepot     <- log(0.751)  ; label("Absolute oral bioavailability F_oral at 66 kg (fraction)")                  # Table S2, final comprehensive model, F_oral = 0.751 (0.699-0.803)
    e_wt_fdepot <- -1.09       ; label("Body-weight power exponent on oral bioavailability, (WT/66)^e_wt_fdepot")   # Table S2, final comprehensive model, C_WGHT,Foral = -1.09 (-1.47 to -0.714)
    e_wt_cl     <-  0.823      ; label("Body-weight power exponent on clearance, (WT/66)^e_wt_cl")                  # Table S2, final comprehensive model, C_WGHT,CL = 0.823 (0.770-0.876)

    # ------------------------------------------------------------------
    # Inter-individual variability; diagonals inverted from the reported
    # CV(%) via log(1 + CV^2) per Table S2 footnote d, off-diagonal
    # reported untransformed (correlation -0.788).
    # ------------------------------------------------------------------
    etalcl + etalrbase_shbg ~ c(0.039221,
                                -0.0655, 0.176237)   # Table S2, final comprehensive model: IIV_CL 20.0% (17.9-22.0), IIV_SBL 43.9% (42.0-45.7), COV_IIV -0.0655 (-0.0741 to -0.0569)

    propSd      <- 0.235   ; label("Proportional residual error SD on total levonorgestrel concentration (fraction)")  # Table S2, final comprehensive model, SIGMA_LNG (CV) = 23.5% (21.6-25.3)
    propSd_shbg <- 0.229   ; label("Proportional residual error SD on SHBG concentration (fraction)")                  # Table S2, final comprehensive model, SIGMA_SHBG (CV) = 22.9% (21.5-24.3)
  })

  model({
    # ------------------------------------------------------------------
    # Physicochemical constants for the closed-form free fraction of LNG
    # (Supplemental Table S3a, "Free fraction of LNG"). KDS and KDA are
    # annotated in the source as unpublished data; ALB is fixed to "an
    # appropriate physiological value" (Discussion: physiological range
    # 530000-830000 nmol/L, i.e. 3.5-5.5 g/dL).
    # ------------------------------------------------------------------
    MWLNG <- 312.5      # Table S3a: molecular weight LNG, g/mol
    KDS   <- 1.82       # Table S3a: dissociation constant SHBG, nmol/L
    KDA   <- 18209      # Table S3a: dissociation constant albumin, nmol/L
    ALB   <- 700000     # Table S3a: albumin concentration, nmol/L
    KALB  <- 1 + ALB / KDA

    # Reference body weight for the oral-bioavailability covariate.
    wtRef <- 66

    # Individual disposition parameters.
    ka <- exp(lka)
    vc <- exp(lvc)
    cl <- exp(lcl + e_wt_cl * log(WT / wtRef) + etalcl)
    vp <- exp(lvp)
    q  <- exp(lq)

    k20 <- cl / vc
    k23 <- q / vc
    k32 <- q / vp

    # SHBG turnover parameters.
    tau        <- exp(ltau)
    ri         <- exp(lri)
    rbase_shbg <- exp(lrbase_shbg + etalrbase_shbg)
    kout_shbg  <- exp(lkout_shbg)
    kin_shbg   <- rbase_shbg * kout_shbg

    # Absolute oral bioavailability with the body-weight power function.
    fdepot <- exp(lfdepot + e_wt_fdepot * log(WT / wtRef))

    # ------------------------------------------------------------------
    # Free fraction of levonorgestrel. A3nM is the molar total LNG
    # concentration in the central compartment: central is in mg and vc in
    # L, so (central / vc) is mg/L and the 1e6 / MWLNG factor converts it
    # to nmol/L exactly as written in Table S3a. The 1e-6 guard keeps
    # A3nM strictly positive so the 1 / A3nM factor is finite before the
    # first dose reaches the central compartment.
    # ------------------------------------------------------------------
    A3nM  <- (central / vc) * (1e6 / MWLNG) + 1e-6
    temp1 <- KALB * KDS + shbg - A3nM
    temp2 <- temp1 / (KALB * 2)
    temp3 <- temp2^2 + (KDS / KALB) * A3nM
    fuLNG <- (1 / A3nM) * (-temp2 + sqrt(temp3))

    # Delayed linear inhibition of SHBG synthesis by levonorgestrel.
    inh <- ri * effect

    # ------------------------------------------------------------------
    # ODE system (Table S3a). Only unbound LNG leaves the central
    # compartment, hence the fuLNG factor on both k20 and k23.
    # ------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - (k20 + k23) * fuLNG * central + k32 * peripheral1
    d/dt(peripheral1) <-  k23 * fuLNG * central - k32 * peripheral1
    d/dt(effect)      <- (1 / tau) * A3nM - (1 / tau) * effect
    d/dt(shbg)        <-  kin_shbg * (1 - inh) - kout_shbg * shbg

    f(depot) <- fdepot
    shbg(0)  <- rbase_shbg

    # ------------------------------------------------------------------
    # Observations. Cc is the total LNG serum/plasma concentration in
    # ng/L: (central [mg] / vc [L]) * 1e6 ng/mg. CcUnbound is reported by
    # the paper (Tables S4, S6, S8) and is derived here without a separate
    # residual error, because the paper attaches residual error to the
    # measured total LNG and SHBG concentrations only.
    # ------------------------------------------------------------------
    Cc        <- (central / vc) * 1e6
    CcUnbound <- fuLNG * Cc

    Cc   ~ prop(propSd)
    shbg ~ prop(propSd_shbg)
  })
}
