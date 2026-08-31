Reinecke_2018_levonorgestrel_iusRelease <- function() {
  description <- paste(
    "Release-only ('in vivo release rate') model of the Reinecke 2018",
    "integrated levonorgestrel (LNG) analysis, covering all three",
    "levonorgestrel-releasing intrauterine systems -- LNG-IUS 20 (Mirena),",
    "LNG-IUS 12 (Kyleena) and LNG-IUS 8 (Jaydess/Skyla) -- in one model. The",
    "paper fitted this model to the levonorgestrel residual-content data",
    "ALONE, independently of levonorgestrel disposition, so that the release",
    "parameters are not influenced by absorption of released drug into serum;",
    "it is the model that generated the published in-vivo release rates. It",
    "carries no pharmacokinetic disposition and no SHBG: the only state is the",
    "device reservoir, and the outputs are the residual device content and the",
    "instantaneous release rate. Two release processes act in parallel. For",
    "LNG-IUS 20 the first is first-order in the remaining content (c12 *",
    "depot) and the second is c13 * (1 + depot / (t1 + t)); for LNG-IUS 12 and",
    "LNG-IUS 8 the first is zero-order (c12) and the second is c13 * depot /",
    "(t1 + t). Unlike the comprehensive model, which estimates the",
    "time-dependency parameter t1 directly per device, this fit reparameterises",
    "it as t1 = fc1 * (c13 - fc2) with fc1 shared across all three devices and",
    "fc2 shared between LNG-IUS 12 and LNG-IUS 8 (and zero for LNG-IUS 20).",
    "The device is selected by the mutually exclusive indicators",
    "FORM_LNG_IUS20, FORM_LNG_IUS12 and FORM_LNG_IUS8, exactly one of which",
    "must be 1. For the levonorgestrel disposition that the released drug",
    "enters, and for the subdermal implant (for which the paper developed no",
    "separate release model), see the companion models named in the reference."
  )
  reference <- paste(
    "Reinecke I, Hofmann B, Mesic E, Drenth HJ, Garmann D.",
    "An integrated population pharmacokinetic analysis to characterize",
    "levonorgestrel pharmacokinetics after different administration routes.",
    "J Clin Pharmacol. 2018 Dec;58(12):1639-1654. doi:10.1002/jcph.1288.",
    "Parameter values from Supplemental Table S1 ('final intrauterine system",
    "release model'); release equations from Supplemental Table S3a. The",
    "corresponding arms of the final comprehensive model, which embed the same",
    "release structure inside the full levonorgestrel/SHBG model, are",
    "modellib('Reinecke_2018_levonorgestrel_lngIus20'),",
    "modellib('Reinecke_2018_levonorgestrel_lngIus12') and",
    "modellib('Reinecke_2018_levonorgestrel_lngIus8'); the subdermal implant is",
    "modellib('Reinecke_2018_levonorgestrel_implant')."
  )
  vignette <- "Reinecke_2018_levonorgestrel_contraceptives"
  units <- list(
    time          = "h",
    dosing        = "mg (levonorgestrel loaded into the device reservoir)",
    concentration = "not applicable (no disposition); residual content in mg, release rate in ug/day"
  )

  compartmentData <- list(
    # The reservoir of the inserted device: levonorgestrel that has not yet been
    # released into the uterine cavity, hence "administration site" rather than
    # a biological matrix.
    depot = list(analyte = "levonorgestrel", units = "mg", specimen = "administration site", verified = TRUE)
  )

  covariateData <- list(
    FORM_LNG_IUS20 = list(
      description        = "Indicator for the levonorgestrel-releasing intrauterine system 20 (LNG-IUS 20, Mirena): 1 = this device, 0 = otherwise.",
      units              = "(binary)",
      type               = "categorical",
      reference_category = "0 (a device other than LNG-IUS 20)",
      notes              = paste(
        "Exactly one of FORM_LNG_IUS20, FORM_LNG_IUS12 and FORM_LNG_IUS8 must",
        "be 1; the three are mutually exclusive. LNG-IUS 20 is the only device",
        "whose first release process is first-order in the remaining reservoir",
        "content, and the only one whose fc2 is zero, so this indicator",
        "switches structure as well as selecting parameter values",
        "(Supplemental Table S3a; Supplemental Table S1 footnote d)."
      ),
      source_name        = "not a data column in the paper; the paper coded the device as a NONMEM branch on the treatment identifier"
    ),
    FORM_LNG_IUS12 = list(
      description        = "Indicator for the levonorgestrel-releasing intrauterine system 12 (LNG-IUS 12, Kyleena): 1 = this device, 0 = otherwise.",
      units              = "(binary)",
      type               = "categorical",
      reference_category = "0 (a device other than LNG-IUS 12)",
      notes              = paste(
        "Mutually exclusive with FORM_LNG_IUS20 and FORM_LNG_IUS8. LNG-IUS 12",
        "shares the zero-order first release process and the fc2 value with",
        "LNG-IUS 8."
      ),
      source_name        = "not a data column in the paper; the paper coded the device as a NONMEM branch on the treatment identifier"
    ),
    FORM_LNG_IUS8 = list(
      description        = "Indicator for the levonorgestrel-releasing intrauterine system 8 (LNG-IUS 8, Jaydess/Skyla): 1 = this device, 0 = otherwise.",
      units              = "(binary)",
      type               = "categorical",
      reference_category = "0 (a device other than LNG-IUS 8)",
      notes              = paste(
        "Mutually exclusive with FORM_LNG_IUS20 and FORM_LNG_IUS12. LNG-IUS 8",
        "shares the zero-order first release process and the fc2 value with",
        "LNG-IUS 12."
      ),
      source_name        = "not a data column in the paper; the paper coded the device as a NONMEM branch on the treatment identifier"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 2645L,
    n_studies      = 3L,
    studies        = c(
      "Phase 2 study 308901 and Phase 3 study 89532 (LNG-IUS 20; 64 + 87 residual-content measurements)",
      "Phase 3 study 310442, LNG-IUS 12 arm (854 residual-content measurements over 5 years)",
      "Phase 3 study 310442, LNG-IUS 8 arm (763 residual-content measurements over 3 years)"
    ),
    age_range      = "18-40 years (median 28 years across the pooled intrauterine data, N = 2790)",
    weight_range   = "39-160 kg (median 66 kg across the pooled intrauterine data, N = 2790)",
    sex_female_pct = 100,
    disease_state  = "Healthy premenopausal women using a levonorgestrel-releasing intrauterine system for contraception",
    dose_range     = paste(
      "A single insertion of one device. The reservoir loading is not stated in",
      "the paper; it is 52 mg (LNG-IUS 20), 19.5 mg (LNG-IUS 12) and 13.5 mg",
      "(LNG-IUS 8) -- see the vignette Errata for how these are recovered from",
      "the paper's own Table 4."
    ),
    notes          = paste(
      "Subject counts are the sums of the intrauterine rows of Table 1",
      "(239 + 94 for LNG-IUS 20, 1306 for LNG-IUS 12, 1245 for LNG-IUS 8",
      "gives 2884 women), restricted here to the 2043 residual-content",
      "measurements that this fit used; the paper does not report how many",
      "distinct women contributed a residual-content measurement, so the",
      "2645 figure is the count of measurements' parent rows and should be",
      "read as approximate. No covariates and no inter-individual variability",
      "were estimated in this fit -- only the release parameters and one",
      "additive residual error per device (Supplemental Table S1). All",
      "parameters were identified with a CV below 17%."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Release process 1. For LNG-IUS 20 this is a FIRST-ORDER rate
    # coefficient (units 1e-6 per hour); for LNG-IUS 12 and LNG-IUS 8 it is
    # a ZERO-ORDER rate (ng/h). Both are scaled inside model() onto the mg/h
    # scale of the reservoir so that ini() carries the tabulated numbers
    # verbatim.
    # ------------------------------------------------------------------
    lc12_ius20 <- log(16.1)   ; label("First-order release rate coefficient c12, LNG-IUS 20 (x 1e-6 per hour)")  # Table S1, C12 LNG-IUS 20 = 16.1 (15.7-16.5) x 1e-6/h
    lc12_ius12 <- log(310)    ; label("Zero-order release rate c12, LNG-IUS 12 (ng/h)")                          # Table S1, C12 LNG-IUS 12 = 310 (303-317) ng/h
    lc12_ius8  <- log(225)    ; label("Zero-order release rate c12, LNG-IUS 8 (ng/h)")                           # Table S1, C12 LNG-IUS 8 = 225 (223-227) ng/h

    # ------------------------------------------------------------------
    # Release process 2. For LNG-IUS 20 c13 is tabulated in units of 1e-6;
    # for LNG-IUS 12 and LNG-IUS 8 it is dimensionless.
    # ------------------------------------------------------------------
    lc13_ius20 <- log(68.0)   ; label("Second release process coefficient c13, LNG-IUS 20 (x 1e-6)")             # Table S1, C13 LNG-IUS 20 = 68.0 (59.7-76.3) x 1e-6
    lc13_ius12 <- log(0.0877) ; label("Second release process coefficient c13, LNG-IUS 12 (dimensionless)")      # Table S1, C13 LNG-IUS 12 = 0.0877 (0.0654-0.110)
    lc13_ius8  <- log(0.0186) ; label("Second release process coefficient c13, LNG-IUS 8 (dimensionless)")       # Table S1, C13 LNG-IUS 8 = 0.0186 (0.0173-0.0199)

    # ------------------------------------------------------------------
    # Time-dependency of release process 2. This fit does NOT estimate t1
    # per device; it estimates two shared correlation factors and derives
    # t1 = fc1 * (c13 - fc2), with fc2 = 0 for LNG-IUS 20 (Table S1
    # footnote d). The derived values Table S1 reports -- 4.37 h, 4580 h and
    # 141 h -- are reproduced by this relation to three significant digits,
    # which is what pins the parameterisation.
    # ------------------------------------------------------------------
    lfc1 <- log(64300)        ; label("Shared correlation factor fc1 in t1 = fc1 * (c13 - fc2) (h)")             # Table S1, FC1 = 64300 (59700-68900)
    lfc2 <- log(0.0164)       ; label("Correlation factor fc2 for LNG-IUS 12 and LNG-IUS 8 in t1 = fc1 * (c13 - fc2)")  # Table S1, FC2 LNG-IUS 8 & 12 = 0.0164 (0.0155-0.0173)

    # ------------------------------------------------------------------
    # Additive residual error on the residual device content, estimated
    # separately per device. Table S1 reports these on the STANDARD
    # DEVIATION scale, unlike the corresponding rows of Table S2, which are
    # variances -- see the vignette Errata for the cross-check that pins the
    # two scales against each other.
    # ------------------------------------------------------------------
    addSd_ius20 <- 1.43       ; label("Additive residual error SD on LNG-IUS 20 residual content (mg)")          # Table S1, SIGMA LNG-IUS 20 res. content = 1.43 (1.19-1.65)
    addSd_ius12 <- 0.435      ; label("Additive residual error SD on LNG-IUS 12 residual content (mg)")          # Table S1, SIGMA LNG-IUS 12 res. content = 0.435 (0.392-0.473)
    addSd_ius8  <- 0.234      ; label("Additive residual error SD on LNG-IUS 8 residual content (mg)")           # Table S1, SIGMA LNG-IUS 8 res. content = 0.234 (0.191-0.269)
  })

  model({
    # ------------------------------------------------------------------
    # Device selection. Exactly one indicator is 1. The 1e-6 factors put
    # every release coefficient on the mg/h scale of the reservoir: c12 is
    # tabulated as 1e-6/h for LNG-IUS 20 and as ng/h (= 1e-6 mg/h) for the
    # other two, and c13 is tabulated as 1e-6 for LNG-IUS 20 and as a
    # dimensionless factor for the other two.
    # ------------------------------------------------------------------
    c12 <- 1e-6 * (FORM_LNG_IUS20 * exp(lc12_ius20) +
                   FORM_LNG_IUS12 * exp(lc12_ius12) +
                   FORM_LNG_IUS8  * exp(lc12_ius8))

    c13 <- FORM_LNG_IUS20 * exp(lc13_ius20) * 1e-6 +
           FORM_LNG_IUS12 * exp(lc13_ius12) +
           FORM_LNG_IUS8  * exp(lc13_ius8)

    # fc2 is zero for LNG-IUS 20 and shared by the other two devices.
    fc1 <- exp(lfc1)
    fc2 <- (1 - FORM_LNG_IUS20) * exp(lfc2)
    t1  <- fc1 * (c13 - fc2)

    # ------------------------------------------------------------------
    # The two release processes (Supplemental Table S3a). LNG-IUS 20 takes
    # the "If LNG-IUS 20" branch; LNG-IUS 12 and LNG-IUS 8 take the
    # "If LNG-IUS 12, LNG-IUS 8, implant" branch. As printed, the LNG-IUS 20
    # form adds a dimensionless 1 to an amount-over-time ratio; the paper
    # describes this release approach as "purely descriptive" (Discussion),
    # and it is transcribed here exactly as written.
    # ------------------------------------------------------------------
    input1 <- FORM_LNG_IUS20 * c12 * depot + (1 - FORM_LNG_IUS20) * c12
    input2 <- c13 * (FORM_LNG_IUS20 + depot / (t1 + t))

    d/dt(depot) <- -input1 - input2

    # ------------------------------------------------------------------
    # Outputs. iusResidual is the levonorgestrel left in the device (mg),
    # the quantity actually measured in the explanted devices; releaseRate
    # converts the instantaneous release (mg/h) to the ug/day of Table 4.
    # ------------------------------------------------------------------
    iusResidual <- depot
    releaseRate <- (input1 + input2) * 1000 * 24

    # Device-specific additive residual error on the measured content.
    addSd_sel <- FORM_LNG_IUS20 * addSd_ius20 +
                 FORM_LNG_IUS12 * addSd_ius12 +
                 FORM_LNG_IUS8  * addSd_ius8

    iusResidual ~ add(addSd_sel)
  })
}
