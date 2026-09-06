Ravva_2015_varenicline <- function() {
  description <- paste(
    "Population PK-PD model of the effect of a single 2 mg oral dose of varenicline on",
    "abstinence-induced cigarette craving in adult smokers (Ravva 2015). Craving is the mean",
    "of the five Smoking Urges Scale items on a 0-100 scale. Two-compartment varenicline PK",
    "with first-order absorption and a lag time, carried forward unestimated from the Ravva",
    "2009 pooled population analysis, supplies the drug concentration; a hypothetical placebo",
    "'kinetic' system, into which a dummy unit dose is given at every dosing time and which",
    "transfers at kon and empties at koff, supplies the placebo driver. The craving response",
    "is the sum of a baseline, a linear placebo term and a linear varenicline term, all on the",
    "logit scale, back-transformed onto 0-100.",
    sep = " "
  )
  reference <- paste(
    "Ravva P, Gastonguay MR, Faessel HM, Lee TC, Niaura R. Pharmacokinetic-pharmacodynamic",
    "modeling of the effect of varenicline on nicotine craving in adult smokers.",
    "Nicotine Tob Res. 2015;17(1):106-113. doi:10.1093/ntr/ntu154",
    sep = " "
  )
  vignette <- "Ravva_2015_varenicline"
  # Stated explicitly because the automatic detector only finds depot and
  # central. Every dosing time carries TWO records: the varenicline dose into
  # depot, and the dimensionless dummy dose of 1 into depot_placebo that drives
  # the hypothetical placebo kinetic system (Ravva 2015 Methods). A simulation
  # that omits the depot_placebo record silently drops the entire placebo term.
  dosing <- c("depot", "depot_placebo")
  units <- list(
    time          = "h",
    dosing        = paste(
      "mg (varenicline, oral); the placebo 'kinetic' system takes a dimensionless dummy dose",
      "of 1 into depot_placebo at every dosing time",
      sep = " "
    ),
    concentration = paste(
      "ng/mL (varenicline plasma concentration Cc); the observation smokingurges is a",
      "Smoking Urges Scale mean (1-5) score on a 0-100 scale",
      sep = " "
    )
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = TRUE for the placebo states because their
  # role is stated explicitly in Ravva 2015 Figure 1 and the Methods.
  compartmentData <- list(
    depot           = list(analyte = "varenicline", units = "mg", specimen = "administration site", verified = FALSE),
    central         = list(analyte = "varenicline", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1     = list(analyte = "varenicline", units = "mg", specimen = "plasma", verified = FALSE),
    depot_placebo   = list(analyte = "placebo response dummy dose", units = "dimensionless", specimen = "not applicable", verified = TRUE),
    central_placebo = list(analyte = "placebo response dummy dose", units = "dimensionless", specimen = "not applicable", verified = TRUE)
  )

  # The PK layer's covariates are those of the Ravva 2009 pooled model, which
  # Ravva 2015 Results explicitly re-uses: "The previously established
  # varenicline population PK model, which included the covariate effects of
  # renal function on clearance (CL/F) and body weight on volume of
  # distribution (V/F) ... was used to predict the individual
  # concentration-time profiles based on the subjects' characteristics in this
  # study". Notes are carried from Ravva_2009_varenicline.R, the on-disk
  # sibling extraction of that model.
  covariateData <- list(
    CRCL = list(
      description        = "Estimated creatinine clearance by the Cockcroft-Gault formula (raw, NOT BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL/F with reference 100 mL/min (Ravva 2009 Table 4 final-model q_CRCL). Carried into Ravva 2015 unchanged: the PK layer was not re-estimated here, only used as MAP-Bayesian priors (Ravva 2015 Methods, 'Pharmacokinetic and Pharmacodynamic Analyses'). The Ravva 2015 cohort is 40 healthy adult smokers, so CRCL sits at the upper end of the Ravva 2009 range.",
      source_name        = "CLcr"
    ),
    WT = list(
      description        = "Total body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on V2/F (estimated exponent), and on V3/F and Q/F (fixed allometric exponents 1 and 0.75) with reference 70 kg (Ravva 2009 Table 4). Ravva 2015 Results reports mean body weight 77 kg (range 59-95) for males and 72 kg (range 58-90) for females.",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Subject age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on V2/F with reference 45 years (Ravva 2009 Table 4 final-model q_AGE). Ravva 2015 Results reports a mean age of 36 years (range 18-63).",
      source_name        = "AGE"
    ),
    RACE_BLACK = list(
      description        = "Black / African American race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = White (reference race for the Ravva 2009 typical individual)",
      notes              = "Power-of-categorical-indicator form on CL/F and V2/F (Ravva 2009 Table 4 q_Black). Ravva 2015 Results reports that all but one of the 40 subjects (97.5%) were White, so this indicator is 0 for essentially the whole Ravva 2015 cohort; it is retained because the PK layer is the Ravva 2009 model in full.",
      source_name        = "Race (Black)"
    ),
    RACE_OTHER = list(
      description        = "Composite 'Other' race indicator pooling Hispanic, Asian, and Other (Ravva 2009 grouping)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = White (reference race for the Ravva 2009 typical individual)",
      notes              = "Power-of-categorical-indicator form on CL/F and V2/F (Ravva 2009 Table 4 q_Other). As with RACE_BLACK, essentially 0 across the Ravva 2015 cohort (97.5% White).",
      source_name        = "Race (Hispanic + Asian + Other)"
    )
  )

  # Screened but not retained. Ravva 2015 Results, 'Relationship Between
  # Varenicline Concentration and Mean (1-5) Craving Score': "Several covariate
  # factors ... were available for inclusion in the population PK-PD model. A
  # graphical inspection of potential covariate-parameter relationships was not
  # indicative of any obvious trends. Further attempts to explain random
  # variability by univariate inclusion of covariate factors directly into the
  # nonlinear mixed-effects model revealed no improvements in goodness-of-fit
  # ... No further covariate model building was conducted." No point estimate
  # is reported for any of them, so none can be encoded.
  covariatesDataExcluded <- list(
    NICOTINE = list(
      description = "Study (not treatment period) baseline plasma nicotine concentration",
      units       = "ng/mL",
      type        = "continuous",
      notes       = "Screened in the Ravva 2015 covariate search and not retained; mean 12.3 ng/mL, range 0-29."
    ),
    COTININE = list(
      description = "Study (not treatment period) baseline plasma cotinine concentration",
      units       = "ng/mL",
      type        = "continuous",
      notes       = "Screened in the Ravva 2015 covariate search and not retained; mean 235 ng/mL, range 74-430."
    ),
    CO_EXHALED = list(
      description = "Study (not treatment period) baseline exhaled carbon monoxide",
      units       = "ppm",
      type        = "continuous",
      notes       = "Screened in the Ravva 2015 covariate search and not retained; mean 20 ppm, range 10-38."
    ),
    SMOKING_YEARS = list(
      description = "Years of smoking",
      units       = "years",
      type        = "continuous",
      notes       = "Screened in the Ravva 2015 covariate search and not retained; subjects had smoked since age 17 (range 12-31) at a mean age of 36 years."
    ),
    CIGS_PER_DAY = list(
      description = "Average number of cigarettes smoked per day",
      units       = "cigarettes/day",
      type        = "continuous",
      notes       = "Screened in the Ravva 2015 covariate search and not retained; mean approximately 21 per day, range 16-40."
    ),
    FTND_Q1 = list(
      description = "Fagerstrom Test for Nicotine Dependence question 1 score (time to first cigarette after waking)",
      units       = "(ordered category)",
      type        = "categorical",
      notes       = "Screened in the Ravva 2015 covariate search and not retained; cohort distribution <5 min: 9, 6-30 min: 23, 31-60 min: 6, >60 min: 2."
    ),
    SEQUENCE = list(
      description = "Crossover treatment sequence (varenicline-then-placebo vs placebo-then-varenicline)",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Screened in the Ravva 2015 covariate search and not retained."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened in the Ravva 2015 covariate search as a PD covariate and not retained; the randomized cohort was 21 male and 19 female. Sex still enters the PK layer indirectly through the Cockcroft-Gault CRCL and through body weight."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 40L,
    n_studies      = 1L,
    age_range      = "18-63 years",
    age_mean       = "36 years",
    weight_range   = "59-95 kg (males); 58-90 kg (females)",
    weight_mean    = "77 kg (males); 72 kg (females)",
    sex_female_pct = 47.5,
    race_ethnicity = c(White = 97.5, Other = 2.5),
    disease_state  = paste(
      "Otherwise healthy adult smokers aged 18-65 not currently intending to quit, smoking at",
      "least 20 cigarettes/day, or 11-19 cigarettes/day with a first cigarette within 30 min",
      "of waking. Overnight smoking and food abstinence, confirmed biochemically with exhaled",
      "CO < 15 ppm.",
      sep = " "
    ),
    dose_range     = "Single 2 mg oral dose (2 x 1 mg tablets) of varenicline versus placebo",
    regions        = "Single center, United States (Center for Behavioral Medicine, The Miriam Hospital)",
    samples        = paste(
      "305 varenicline plasma concentrations, 387 placebo-period craving responses and 390",
      "varenicline-period craving responses. Only the abstinence window from time 0 (just",
      "before dosing) to 4 h postdose, before the cue reactivity session, was analysed.",
      sep = " "
    ),
    notes          = paste(
      "Randomized, double-blind, placebo-controlled, two-period crossover with a 7-day washout",
      "(Ravva 2015 Methods, 'Study Design' and 'Study Procedures'; demographics from Results,",
      "'Subject Disposition and Smoking History'). Smoking history: about 21 cigarettes/day",
      "(range 16-40), smoking since age 17 (range 12-31), a mean Fagerstrom Test for Nicotine",
      "Dependence total score of 5.6 (SE 0.31) out of 10, and 7 of 40 subjects reporting prior",
      "bupropion (Zyban) quit attempts. Two subjects did not complete both crossover periods;",
      "their partial data were retained in the mixed-effects analysis.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Varenicline PK layer.
    #
    # Ravva 2015 did NOT estimate any PK parameter. Methods, 'Pharmacokinetic
    # and Pharmacodynamic Analyses': "The final PK model and population
    # parameters previously obtained in a varenicline pooled analysis were
    # used as prior information in a maximum a posteriori (MAP) Bayesian
    # analysis (POSTHOC option) of the new PK data from this study" -- that
    # pooled analysis is reference 30, Ravva 2009, which is packaged here as
    # inst/modeldb/specificDrugs/Ravva_2009_varenicline.R. Every fixed effect
    # and the whole IIV block below are therefore transcribed verbatim from
    # that on-disk sibling model (its own source trace points at Ravva 2009
    # Table 4) and wrapped in fixed(), because in the Ravva 2015 analysis they
    # are prior constants rather than estimated quantities.
    # ------------------------------------------------------------------
    lka     <- fixed(log(1.69)); label("First-order absorption rate constant (Ka, 1/h)")            # Ravva 2009 Table 4 final-model Ka, via Ravva_2009_varenicline.R
    lcl     <- fixed(log(10.4)); label("Apparent clearance (CL/F, L/h) at reference covariates")    # Ravva 2009 Table 4 final-model q_CL
    lvc     <- fixed(log(337));  label("Apparent central volume of distribution (V2/F, L)")         # Ravva 2009 Table 4 final-model q_V2
    lvp     <- fixed(log(78.1)); label("Apparent peripheral volume of distribution (V3/F, L)")      # Ravva 2009 Table 4 final-model q_V3
    lq      <- fixed(log(2.08)); label("Apparent intercompartmental clearance (Q/F, L/h)")          # Ravva 2009 Table 4 final-model q_Q
    ltlag   <- fixed(log(0.43)); label("Absorption lag time (Alag, h)")                             # Ravva 2009 Table 4 final-model q_Alag

    e_crcl_cl   <- fixed(0.54); label("Power exponent of CRCL/100 on CL/F (unitless)")                     # Ravva 2009 Table 4 q_CRCL on CL/F
    e_black_cl  <- fixed(1.16); label("Multiplicative factor on CL/F for Black vs White (unitless)")       # Ravva 2009 Table 4 q_Black on CL/F
    e_other_cl  <- fixed(1.11); label("Multiplicative factor on CL/F for Other vs White (unitless)")       # Ravva 2009 Table 4 q_Other on CL/F

    e_wt_vc     <- fixed(0.77); label("Power exponent of WT/70 on V2/F (unitless)")                        # Ravva 2009 Table 4 q_WT on V2/F
    e_age_vc    <- fixed(0.13); label("Power exponent of AGE/45 on V2/F (unitless)")                       # Ravva 2009 Table 4 q_AGE on V2/F
    e_black_vc  <- fixed(0.92); label("Multiplicative factor on V2/F for Black vs White (unitless)")       # Ravva 2009 Table 4 q_Black on V2/F
    e_other_vc  <- fixed(0.71); label("Multiplicative factor on V2/F for Other vs White (unitless)")       # Ravva 2009 Table 4 q_Other on V2/F

    e_wt_vp     <- fixed(1);    label("Allometric exponent of WT/70 on V3/F (unitless)")                   # Ravva 2009 Table 4 q_WT on V3/F = 1 (Fixed)
    e_wt_q      <- fixed(0.75); label("Allometric exponent of WT/70 on Q/F (unitless)")                    # Ravva 2009 Table 4 q_WT on Q/F = 0.75 (Fixed)

    # Ravva 2009 Table 4 interindividual-variance BLOCK on Ka, CL and V2, in
    # that row/column order; entries are the lower-triangular variances and
    # covariances. Fixed here for the same reason as the fixed effects: in
    # Ravva 2015 this block is the MAP-Bayesian prior, not an estimate.
    etalka + etalcl + etalvc ~ fixed(c(
       0.49,
      -0.009, 0.061,
       0.24,  0.006, 0.25
    ))

    # ------------------------------------------------------------------
    # Placebo 'kinetic' system.
    #
    # Ravva 2015 Methods: "The time-varying behavior of the placebo response
    # was modeled as a separate hypothetical kinetic system with a 'dummy'
    # unit dose introduced at times of placebo or varenicline dosing."
    # Figure 1 draws it as a depot box feeding a C_PBO box at kon, with koff
    # leaving C_PBO. Supplementary Table S1 names kon "rate constant for onset
    # of placebo response" and koff "rate constant for offset of placebo
    # response".
    #
    # Both are wrapped in fixed() on the authority of Table S1 footnote a:
    # "Precision of the k_on and k_off parameter estimates was obtained when
    # fitting the placebo scores alone. These parameters were later fixed when
    # modeling the drug treatment effect." Their IIV is fixed for the same
    # reason -- the individual-specific values carried into the combined fit
    # come from the placebo-only run, so neither the typical value nor the
    # spread was re-estimated at the step this file encodes.
    #
    # Log-transformed because Table S1 reports the IIV of both as a percent
    # CV, i.e. log-normal, in contrast to the three effect parameters below
    # whose IIV is reported as an SD on the parameter's own scale.
    # ------------------------------------------------------------------
    lkon_placebo  <- fixed(log(0.0112)); label("Onset rate constant of the placebo 'kinetic' system (kon, 1/h)")  # Ravva 2015 Supplementary Table S1: k_on = 0.0112 /h (%RSE 142)
    lkoff_placebo <- fixed(log(0.130));  label("Offset rate constant of the placebo 'kinetic' system (koff, 1/h)") # Ravva 2015 Supplementary Table S1: k_off = 0.130 /h (%RSE 65.4)

    etalkon_placebo  ~ fixed(0.97987)  # Table S1 k_on IIV %CV = 129 (%RSE 154);  omega^2 = log(1 + 1.29^2) = 0.97987
    etalkoff_placebo ~ fixed(0.88216)  # Table S1 k_off IIV %CV = 119 (%RSE 154); omega^2 = log(1 + 1.19^2) = 0.88216

    # ------------------------------------------------------------------
    # Linear logit-scale effect model.
    #
    # Ravva 2015 Methods: "The concentrations of varenicline and placebo in
    # the central compartment of each kinetic system were linked to effect by
    # a linear PD model with slope parameters, DSLP and PSLP, respectively.
    # Response was described as the net effect of the baseline (E0), placebo
    # (PSLP*C_PBO), and drug treatment (DSLP*C_DRUG)." Supplementary Table S1
    # footnote b puts all three on the logit scale: "Parameters are presented
    # in the logit scale as data were modeled using logistic transformation to
    # constrain between 0 and 100."
    #
    # All three carry ADDITIVE interindividual variability, per Methods
    # ("Inter-individual and residual random effects were characterized by an
    # additive error model") and per Table S1, which reports their IIV as an
    # SD on the parameter's own scale rather than as a %CV. An additive eta is
    # also the only encoding compatible with the two negative slopes: it lets
    # a subject's individual slope take either sign, which is exactly what the
    # per-subject IPRED panels of Figure 3 show.
    # ------------------------------------------------------------------
    logite0       <-  0.991; label("Baseline craving response E0 on the logit scale (unitless)")                       # Ravva 2015 Supplementary Table S1: E0 = 0.991 (%RSE 29.6); 100*expit(0.991) = 72.9 on the Smoking Urges Scale
    slope_placebo <- -3.32;  label("Placebo slope PSLP (logit per unit of the placebo kinetic-system amount)")         # Ravva 2015 Supplementary Table S1: PSLP = -3.32 (%RSE 50.6)
    slope_drug    <- -0.192; label("Varenicline slope DSLP (logit per ng/mL of plasma varenicline)")                   # Ravva 2015 Supplementary Table S1: DSLP = -0.192 (%RSE 33.5)

    etalogite0       ~ 3.2041     # Table S1 E0 IIV SD = 1.79 (%RSE 32.0);   variance = 1.79^2  = 3.2041
    etaslope_placebo ~ 16.1604    # Table S1 PSLP IIV SD = 4.02 (%RSE 85.8); variance = 4.02^2  = 16.1604
    etaslope_drug    ~ 0.157609   # Table S1 DSLP IIV SD = 0.397 (%RSE 40.8); variance = 0.397^2 = 0.157609

    # Residual error on the craving observation, additive on the logit scale.
    # Ravva 2015 Supplementary Table S1 residual-error row: sigma_(3,3) add,PD
    # SD = 1.06 (%RSE 17.0). Encoded as logitNorm on 0-100 so the additive
    # error acts on the same logit scale the effect parameters live on, which
    # is what the footnote-b logistic transformation implies.
    addSd <- 1.06; label("Additive residual SD of the craving response on the logit scale")  # Ravva 2015 Supplementary Table S1: sigma_(3,3) add,PD SD = 1.06 (%RSE 17.0)
  })
  model({
    # ---------------- Varenicline PK (Ravva 2009 final model) ----------------
    ka  <- exp(lka + etalka)
    cl  <- exp(lcl + etalcl) *
           (CRCL / 100)^e_crcl_cl *
           e_black_cl^RACE_BLACK *
           e_other_cl^RACE_OTHER
    vc  <- exp(lvc + etalvc) *
           (WT / 70)^e_wt_vc *
           (AGE / 45)^e_age_vc *
           e_black_vc^RACE_BLACK *
           e_other_vc^RACE_OTHER
    vp  <- exp(lvp) * (WT / 70)^e_wt_vp
    q   <- exp(lq)  * (WT / 70)^e_wt_q

    tlag <- exp(ltlag)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ---------------- Placebo 'kinetic' system ----------------
    kon_placebo  <- exp(lkon_placebo  + etalkon_placebo)
    koff_placebo <- exp(lkoff_placebo + etalkoff_placebo)

    # ---------------- ODE system ----------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # A dimensionless dummy dose of 1 is administered into depot_placebo at
    # every varenicline or placebo dosing time (Ravva 2015 Methods). The
    # contents of central_placebo are the paper's C_PBO driver.
    d/dt(depot_placebo)   <- -kon_placebo * depot_placebo
    d/dt(central_placebo) <-  kon_placebo * depot_placebo - koff_placebo * central_placebo

    alag(depot) <- tlag

    # Plasma varenicline (dose in mg, volumes in L -> mg/L; reported in ng/mL).
    Cc <- 1000 * central / vc

    # Individual effect parameters. Each eta enters ADDITIVELY (Ravva 2015
    # Methods; Table S1 reports these three IIVs as SDs on the parameter's own
    # scale). Written one per line, following the Zhang_2023_brazikumab_crp.R
    # `iplac_i <- iplac + etaiplac` idiom, so nlmixr2 can mu-reference them.
    logite0_i       <- logite0 + etalogite0
    slope_placebo_i <- slope_placebo + etaslope_placebo
    slope_drug_i    <- slope_drug + etaslope_drug

    # Net craving response: baseline + placebo + drug, summed on the logit
    # scale and back-transformed onto the 0-100 Smoking Urges Scale.
    smokingurges <- expit(
      logite0_i +
        slope_placebo_i * central_placebo +
        slope_drug_i * Cc,
      0, 100
    )
    smokingurges ~ logitNorm(addSd, 0, 100)
  })
}
