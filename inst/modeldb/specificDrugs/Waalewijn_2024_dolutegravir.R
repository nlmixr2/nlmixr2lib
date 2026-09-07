Waalewijn_2024_dolutegravir <- function() {
  description <- paste(
    "Two-compartment population PK model for once-daily oral dolutegravir in",
    "42 African children weighing >14 kg taking dolutegravir with a low-fat",
    "breakfast in the CHAPAS-4 trial (Waalewijn 2024). Absorption is a Savic",
    "analytical Erlang transit chain (mean transit time 1.08 h, 10.1 transit",
    "compartments) feeding a depot that empties at first-order rate ka into",
    "the central compartment. Disposition is two-compartment with apparent",
    "clearance and volumes scaled allometrically on total body weight (fixed",
    "exponents 0.75 on CL and Q, 1.0 on Vc and Vp) referenced to 30 kg.",
    "One covariate effect is retained: an emtricitabine/tenofovir alafenamide",
    "(FTC/TAF) nucleos(t)ide backbone reduces relative bioavailability by",
    "19.6% relative to the standard-of-care backbone (lamivudine/abacavir or",
    "lamivudine/zidovudine). Between-subject variability is a single eta on",
    "clearance (21.1% CV); between-occasion variability is carried on relative",
    "bioavailability, ka and mean transit time across three dosing occasions,",
    "with the bioavailability BOV inflated 3.35-fold on the occasion whose",
    "dose was taken at home and not observed by study staff. Residual error",
    "is combined proportional (9.1%) and additive (0.07 mg/L).",
    sep = " "
  )
  reference <- paste(
    "Waalewijn H, Wasmann RE, Bamford A, Gibb DM, McIlleron HM, Colbers A,",
    "Burger DM, Denti P, and the CHAPAS-4 trial team.",
    "Population Pharmacokinetics of Dolutegravir in African Children:",
    "Results From the CHAPAS-4 Trial.",
    "J Pediatric Infect Dis Soc. 2024;13(10):533-536.",
    "doi:10.1093/jpids/piae076.",
    sep = " "
  )
  vignette <- "Waalewijn_2024_dolutegravir"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric size descriptor on all four disposition parameters, with",
        "exponents fixed at 0.75 for the flow parameters (CL, Q) and 1.0 for",
        "the volume parameters (Vc, Vp), referenced to a 30 kg child",
        "(Waalewijn 2024 Supplementary Table S3 footnote b; NONMEM $PK",
        "TVWT = 30, ALLMCL_WT = (WT/TVWT)**0.75, ALLMV_WT = (WT/TVWT)).",
        "Total body weight and fat-free mass were both tested as body size",
        "descriptors and total body weight fitted best (Results paragraph 2);",
        "fat-free mass is therefore documented in covariatesDataExcluded.",
        sep = " "
      ),
      source_name        = "WT"
    ),
    CONMED_TAF = list(
      description        = "Emtricitabine/tenofovir alafenamide (FTC/TAF) nucleos(t)ide backbone indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (standard-of-care backbone: lamivudine/abacavir or lamivudine/zidovudine)",
      notes              = paste(
        "Time-fixed per subject (CHAPAS-4 randomised each child to one",
        "backbone). Applied as a linear additive effect on relative",
        "bioavailability: F *= (1 + e_taf_fdepot * CONMED_TAF) with",
        "e_taf_fdepot = -0.196, i.e. F = 0.804 on FTC/TAF versus 1 (fixed) on",
        "the SOC reference (Supplementary Table S3 'Effect of FTC/TAF on",
        "bioavailability' -19.6%; NONMEM $PK BB_BIO = THETA(14) = 0.804 when",
        "BB == 1, THETA(4) = 1 FIX otherwise).",
        "IMPORTANT -- the indicator flags the WHOLE FTC/TAF fixed-dose",
        "backbone, not isolated tenofovir alafenamide exposure: both arms",
        "receive a two-drug nucleos(t)ide backbone, so the estimate is the",
        "contrast FTC/TAF versus 3TC/ABC or 3TC/ZDV and cannot be attributed",
        "to TAF alone. The paper explicitly leaves the mechanism unresolved",
        "(Discussion paragraph 2) and notes that adult FTC/TAF + dolutegravir",
        "studies did not reproduce the reduction.",
        "The two SOC backbones were pooled only after being tested",
        "separately: relative to 3TC/ABC, the 3TC/ZDV effect on",
        "bioavailability was -8.6% [-30.9% to +19.2%] (95% CI spanning zero)",
        "while FTC/TAF was -23.3% [-37.7% to -6.7%] (Supplementary Table S4).",
        "Waalewijn 2024 also tested TAF and tenofovir AUC as continuous",
        "replacements for this categorical indicator; neither improved the fit",
        "(Results paragraph 3), so they are documented in",
        "covariatesDataExcluded rather than carried here.",
        sep = " "
      ),
      source_name        = "BB"
    ),
    OCC = list(
      description        = "Dosing-occasion index for between-occasion variability",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Three occasions, matching the NONMEM $PK IF(OCC==1/2/3) blocks that",
        "select the per-occasion BOV etas on relative bioavailability, ka and",
        "mean transit time. Occasion 1 is the occasion whose dose was taken at",
        "home on the day before the pharmacokinetic assessment and was NOT",
        "observed by study staff -- the occasion the pre-dose sample belongs",
        "to. Occasions 2 and 3 are doses observed by study staff on the",
        "assessment day. Waalewijn 2024 estimates a single shared BOV",
        "magnitude across occasions ($OMEGA BLOCK(1) SAME) plus a scaling",
        "factor of 3.35 (95% CI 2.25-4.86) that multiplies the bioavailability",
        "eta on occasion 1 only, giving 22.1% CV for observed doses and about",
        "75% CV for the unobserved pre-dose occasion (Results paragraph 5;",
        "Supplementary Table S3). Set OCC constant across a dosing interval;",
        "the value on a record selects that record's occasion etas.",
        sep = " "
      ),
      source_name        = "OCC"
    )
  )

  # Covariates that Waalewijn 2024 screened but did not retain in the final
  # model. Documented for provenance only; none is referenced in model().
  covariatesDataExcluded <- list(
    FFM = list(
      description = "Fat-free mass",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Tested against total body weight as the allometric body size",
        "descriptor; total body weight fitted better and was retained",
        "(Results paragraph 2). The control stream retains the FFM machinery",
        "but the final $PK block scales on WT only. The sex-specific constants",
        "it carries (WHSMAX 42.92 male / 37.99 female, WHS50 30.93 male /",
        "35.98 female, with height in metres squared) are Janmahasatian's,",
        "not Al-Sallami's; the reference value is TVFFM = 23.9 kg. Cohort",
        "median 22.6 kg (Supplementary Table S2).",
        sep = " "
      )
    ),
    CRCL = list(
      description = "Creatinine clearance estimated with the Schwartz formula",
      units       = "mL/min",
      type        = "continuous",
      notes       = paste(
        "Tested on clearance as a renal-function descriptor; did not",
        "significantly improve the fit (Results paragraph 4). Cohort median",
        "113 mL/min, range 71.8-247 (Supplementary Table S2). The control",
        "stream imputes missing values to the cohort median 113.",
        sep = " "
      )
    ),
    ALT = list(
      description = "Serum alanine transaminase",
      units       = "U/L",
      type        = "continuous",
      notes       = paste(
        "Liver-function biomarker tested on clearance; not retained",
        "(Results paragraph 4). Control stream imputes missing values to",
        "28.0 U/L.",
        sep = " "
      )
    ),
    AST = list(
      description = "Serum aspartate transaminase",
      units       = "U/L",
      type        = "continuous",
      notes       = paste(
        "Liver-function biomarker tested on clearance; not retained",
        "(Results paragraph 4). Control stream imputes missing values to",
        "37.6 U/L.",
        sep = " "
      )
    ),
    FORM_DTG_DT = list(
      description = "Dolutegravir 25 mg dispersible-tablet formulation indicator (reference: 50 mg film-coated tablet)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Tested on the absorption parameters including bioavailability and did",
        "not significantly improve the fit (Results paragraph 4). This is a",
        "notable negative result: Chandasana 2024 and the ODYSSEY trial both",
        "report substantially higher dispersible-tablet bioavailability, and",
        "Waalewijn 2024 attributes the absence of a formulation effect here to",
        "every child having dosed with a low-fat breakfast (Discussion",
        "paragraph 3). Because formulation is confounded with weight band in",
        "this trial (25 mg DT for 14 to <20 kg, 50 mg FCT above 20 kg), the",
        "null result is specific to the fed condition.",
        sep = " "
      )
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "dolutegravir", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "dolutegravir", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "dolutegravir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 42L,
    n_studies      = 1L,
    n_observations = 358L,
    age_range      = "5.46-15.5 years",
    age_median     = "10.9 years",
    weight_range   = "15.9-53.0 kg",
    weight_median  = "27.8 kg",
    sex_female_pct = 47.6,
    disease_state  = paste(
      "Children living with HIV taking second-line antiretroviral therapy in",
      "the CHAPAS-4 trial (ISRCTN22964075), randomised to an FTC/TAF",
      "nucleos(t)ide backbone (n = 21) or a standard-of-care backbone",
      "(3TC/ZDV n = 12, 3TC/ABC n = 9).",
      sep = " "
    ),
    dose_range     = paste(
      "Once-daily dolutegravir per WHO weight-band recommendations: 25 mg as",
      "five 5 mg dispersible tablets for 14 to <20 kg (n = 10), one 50 mg",
      "film-coated tablet for >20 kg (n = 32). All doses on the",
      "pharmacokinetic assessment day were taken with a standardised low-fat",
      "breakfast (250 kcal, 5% fat) and observed by study staff.",
      sep = " "
    ),
    regions        = "Uganda, Zambia, and Zimbabwe.",
    sampling       = paste(
      "Steady-state intensive sampling pre-dose and at 1, 2, 4, 6, 8, 12, and",
      "24 h post-dose, with an additional 0.5 h sample for children taking",
      "TAF. Dolutegravir assayed by validated LC-MS/MS.",
      sep = " "
    ),
    notes          = paste(
      "Enrolment January 2019 to March 2021. 42 children contributed 358",
      "dolutegravir concentrations, 2 of them below the limit of",
      "quantification and handled by Beal's M6 method. Fit in NONMEM 7.5 with",
      "PsN 5.3.1; parameter uncertainty from sampling importance resampling.",
      "Baseline demographics from Supplementary Table S2 (Total column).",
      sep = " "
    )
  )

  ini({
    # ---- Structural fixed effects ------------------------------------------
    # All disposition values are the typical values for the 30 kg reference
    # child in Waalewijn 2024 Supplementary Table S3 (final estimates with
    # sampling-importance-resampling 95% CIs). Note that the $THETA block of
    # the supplement's control stream holds INITIAL estimates, not these final
    # ones -- e.g. $THETA(1) CL = 1.21 L/h and $THETA(2) V = 13.2 L both fall
    # far outside the Table S3 SIR confidence intervals below -- so Table S3 is
    # the authority for every value in this block.
    lcl  <- log(0.722); label("Apparent oral clearance CL/F at 30 kg (L/h)")                         # Table S3 'Clearance (L/h)' 0.722 (0.645, 0.797)
    lvc  <- log(6.66);  label("Apparent central volume of distribution Vc/F at 30 kg (L)")           # Table S3 'Central distribution volume (L)' 6.66 (5.84, 7.46)
    lq   <- log(0.278); label("Apparent inter-compartmental clearance Q/F at 30 kg (L/h)")           # Table S3 'Intercompartmental clearance (L/h)' 0.278 (0.157, 0.457)
    lvp  <- log(2.00);  label("Apparent peripheral volume of distribution Vp/F at 30 kg (L)")        # Table S3 'Peripheral distribution volume (L)' 2.00 (1.45, 2.71)
    lka  <- log(0.950); label("First-order absorption rate constant ka, depot to central (1/h)")     # Table S3 'Absorption rate constant (h-1)' 0.950 (0.704, 1.31)
    lmtt <- log(1.08);  label("Mean transit time through the Savic transit chain (h)")               # Table S3 'Mean transit time (h)' 1.08 (0.837, 1.33)
    lnn  <- log(10.1);  label("Number of transit compartments N (Savic chain, non-integer allowed) (unitless)") # Table S3 'Number of transit compartments' 10.1 (7.66, 14.5)

    lfdepot <- fixed(log(1)); label("Relative bioavailability F on the SOC-backbone reference (unitless)") # Table S3 'Relative bioavailability [F]' = 1 fixed; NONMEM $THETA(4) (1) FIX

    # ---- Allometric exponents (fixed a priori) -----------------------------
    # Supplementary Methods 'Population pharmacokinetic analysis' paragraph 2:
    # "We incorporated fixed allometric scaling exponents of 0.75 for clearance
    # and 1.0 for volume in disposition parameters". Reference weight 30 kg
    # (Table S3 footnote b; NONMEM $PK TVWT = 30).
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL and Q (unitless)")
    e_wt_vc <- fixed(1.0);  label("Allometric exponent on Vc and Vp (unitless)")

    # ---- Covariate effect on relative bioavailability ----------------------
    # Linear additive form F *= (1 + e_taf_fdepot * CONMED_TAF), so the FTC/TAF
    # arm carries F = 1 - 0.196 = 0.804, matching NONMEM $THETA(14) = 0.804
    # exactly (the one structural THETA in the control stream that equals its
    # Table S3 final value).
    e_taf_fdepot <- -0.196; label("Factor change in F on an FTC/TAF backbone (fraction; -19.6% with CONMED_TAF = 1)") # Table S3 'Effect of FTC/TAF on bioavailability' -19.6% (-30.8%, -8.13%)

    # ---- Between-subject variability ---------------------------------------
    # Table S3 footnote c defines the reported CV% on the log-scale SD:
    # %CV = sqrt(omega^2) * 100%, so omega^2 = (CV/100)^2 directly -- NOT the
    # log(1 + CV^2) transform used by papers that report the exponentiated CV.
    # Clearance is the only parameter carrying BSV; every other $OMEGA in the
    # control stream's BSV block ($OMEGA 2-9: V, KA, BIO, V3, Q, V4, Q2, MTT)
    # is 0 FIX.
    etalcl ~ 0.044521  # Table S3 'Clearance' BSV 21.1% (16.0, 27.0) -> 0.211^2 = 0.044521

    # ---- Between-occasion variability --------------------------------------
    # Three occasions (NONMEM $PK IF(OCC==1/2/3)), with $OMEGA BLOCK(1) SAME
    # sharing one variance across occasions for each of F, ka and MTT.
    #
    # Bioavailability: the control stream multiplies the OCCASION 1 eta by
    # THETA(13) = 3.35 (BOVBIO = ETA(13) * PDBOVBIO), which scales that
    # occasion's standard deviation rather than adding a fixed-effect shift.
    # Occasion 1 is the unobserved home dose taken the day before the PK visit.
    # nlmixr2lib has no canonical name for a per-occasion eta-SD multiplier, so
    # the scaling is folded exactly into occasion 1's variance instead of being
    # carried as a separate estimated parameter:
    #   observed occasions : SD = 0.221          -> omega^2 = 0.048841
    #   occasion 1         : SD = 0.221 * 3.35   -> omega^2 = 0.740350^2 = 0.548118
    # 0.221 * 3.35 = 0.740 reproduces the paper's "about 75%" for the pre-dose
    # samples (Results paragraph 5). This is a re-parameterisation, not an
    # approximation -- the two forms are numerically identical -- but it means
    # the 3.35 folds cannot be re-estimated as a single parameter. See the
    # vignette 'Assumptions and deviations' section.
    etaiov_fdepot_1 ~ fixed(0.548118)  # Table S3 BOV F 22.1% x scaling factor 3.35 (2.25, 4.86); (0.221*3.35)^2
    etaiov_fdepot_2 ~ 0.048841         # Table S3 'Bioavailability of observed doses' BOV 22.1% (17.2, 28.6) -> 0.221^2
    etaiov_fdepot_3 ~ fixed(0.048841)  # NONMEM $OMEGA BLOCK(1) SAME

    etaiov_ka_1 ~ 0.966289         # Table S3 'Absorption rate constant' BOV 98.3% (77.8, 121) -> 0.983^2
    etaiov_ka_2 ~ fixed(0.966289)  # NONMEM $OMEGA BLOCK(1) SAME
    etaiov_ka_3 ~ fixed(0.966289)  # NONMEM $OMEGA BLOCK(1) SAME

    etaiov_mtt_1 ~ 0.567009         # Table S3 'Mean transit time' BOV 75.3% (59.7, 95.0) -> 0.753^2
    etaiov_mtt_2 ~ fixed(0.567009)  # NONMEM $OMEGA BLOCK(1) SAME
    etaiov_mtt_3 ~ fixed(0.567009)  # NONMEM $OMEGA BLOCK(1) SAME

    # ---- Residual error ----------------------------------------------------
    # Combined proportional + additive, matching the control stream's $ERROR
    # W = SQRT(ADD**2 + PROP**2). Table S3 reports the assembled additive term
    # (the control stream builds it as THETA(6) + 0.2 * censoring threshold).
    propSd <- 0.091; label("Proportional residual error (fraction)")  # Table S3 'Proportional error (%)' 9.1 (7.9-10.7)
    addSd  <- 0.07;  label("Additive residual error (mg/L)")          # Table S3 'Additive error (mg/L)' 0.07 (0.04-0.10)
  })

  model({
    # --- 1. Occasion indicators and between-occasion random effects ---------
    # Mutually exclusive; a record with OCC outside 1-3 carries no BOV.
    occ1 <- (OCC == 1)
    occ2 <- (OCC == 2)
    occ3 <- (OCC == 3)

    iov_fdepot <- occ1 * etaiov_fdepot_1 + occ2 * etaiov_fdepot_2 + occ3 * etaiov_fdepot_3
    iov_ka     <- occ1 * etaiov_ka_1     + occ2 * etaiov_ka_2     + occ3 * etaiov_ka_3
    iov_mtt    <- occ1 * etaiov_mtt_1    + occ2 * etaiov_mtt_2    + occ3 * etaiov_mtt_3

    # --- 2. Individual PK parameters ----------------------------------------
    # Allometric scaling on total body weight referenced to 30 kg. Only CL
    # carries between-subject variability.
    cl <- exp(lcl + etalcl) * (WT / 30)^e_wt_cl
    vc <- exp(lvc)          * (WT / 30)^e_wt_vc
    q  <- exp(lq)           * (WT / 30)^e_wt_cl
    vp <- exp(lvp)          * (WT / 30)^e_wt_vc

    ka  <- exp(lka  + iov_ka)
    mtt <- exp(lmtt + iov_mtt)
    nn  <- exp(lnn)

    # Relative bioavailability: SOC reference F = 1 (fixed), FTC/TAF arm 0.804.
    fdepot <- exp(lfdepot + iov_fdepot) * (1 + e_taf_fdepot * CONMED_TAF)

    # --- 3. Micro-constants --------------------------------------------------
    kel <- cl / vc
    k23 <- q  / vc
    k32 <- q  / vp

    # --- 4. ODE system -------------------------------------------------------
    # Dose is administered to depot (cmt = "depot"). Absorption is the Savic
    # 2007 analytical Erlang transit chain, written out by hand rather than
    # through the rxode2 transit() macro. The hand-expanded form is a literal
    # transcription of the control stream and is numerically identical to
    # rxode2's transit() built-in; it is used because transit() silently
    # evaluates to zero inside an nlmixr2 model function under rxode2 5.1.7
    # (see the vignette 'Assumptions and deviations' section).
    #
    # Control stream $PK / $DES:
    #   KTR     = (NN+1)/MTT
    #   PIZZA   = LOG(BIO*PD*KTR + 1E-12) - GAMLN(NN+1)   [PD = dose amount]
    #   TEMPO   = T - TDOS                                 [time after dose]
    #   TRANSIT = EXP(PIZZA + NN*LOG(KTR*TEMPO) - KTR*TEMPO)
    #
    # The 1e-12 offset is the control stream's own guard (its comment: "without
    # +0.00001, it won't work with ETAs in bioavailability"); here it also keeps
    # log() finite on records that precede the first dose.
    ktr        <- (nn + 1) / mtt
    tempo      <- tad(depot)
    transit_in <- exp(log(fdepot * podo(depot) * ktr + 1e-12) - lgamma(nn + 1) +
                        nn * log(ktr * tempo) - ktr * tempo)

    # The depot then empties at first-order rate ka into central; disposition
    # is two-compartment (the control stream's V4 / Q2 second peripheral
    # compartment is 0 FIX and is commented out of $MODEL).
    d/dt(depot)       <- transit_in - ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k23 * central + k32 * peripheral1
    d/dt(peripheral1) <-                              k23 * central - k32 * peripheral1

    # --- 5. Bioavailability --------------------------------------------------
    # Suppress the dose bolus on depot so the analytical transit() chain is the
    # only input pathway, exactly as the control stream does with F1 = 0.
    f(depot) <- 0

    # --- 6. Observation and residual error -----------------------------------
    # Dose in mg / volume in L -> mg/L, the unit Table S3 reports the additive
    # error in. NONMEM $ERROR IPRED = A(2)/V.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
