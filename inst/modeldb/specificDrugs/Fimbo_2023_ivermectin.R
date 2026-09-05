# Two-compartment population PK of oral ivermectin with a Savic analytical
# transit-absorption chain, allometric scaling on CL/F and Vc/F, a relative
# bioavailability increase in the 3 mg dose pool, and a two-class mixture on
# mean transit time, fitted to 468 individuals dosed during a lymphatic
# filariasis mass-drug-administration campaign in Tanga, Tanzania (Fimbo 2023,
# CPT Pharmacometrics Syst Pharmacol 12(12):1884-1896;
# doi:10.1002/psp4.13038).

Fimbo_2023_ivermectin <- function() {
  description <- paste(
    "Two-compartment population PK model for oral ivermectin (IVM) in 468",
    "individuals aged 5-78 years and weighing more than 15 kg who received a",
    "single height-based IVM dose (3, 6, 9, or 12 mg) together with",
    "albendazole 400 mg during a lymphatic filariasis mass drug administration",
    "(MDA) campaign in the Mkinga district, Tanga region, Tanzania (1404",
    "plasma samples drawn predose and at 2, 4, and 6 h). Absorption is the",
    "Savic analytical transit chain: a gamma-distributed input with mean",
    "transit time MTT = 1.523 h through NN = 6 transit compartments (fixed)",
    "feeds an absorption compartment that empties into the central",
    "compartment at ka = 0.708 /h. Body weight enters as allometric scaling",
    "with exponents fixed at 0.75 on CL/F and 1 on Vc/F, referenced to 70 kg;",
    "intercompartmental clearance Q/F and peripheral volume Vp/F are NOT",
    "allometrically scaled, which is what makes the model's terminal",
    "half-life fall with increasing body weight (66 h in the 3 mg pool vs 38",
    "h in the 12 mg pool). Dose group is the only retained covariate:",
    "relative bioavailability is 48.2% higher in the 3 mg dose pool than in",
    "the 6, 9, and 12 mg pools. A two-class NONMEM mixture on MTT identifies",
    "a 16.1% subpopulation whose mean transit time is 97.0% longer than the",
    "rest of the population; it is supplied here as the per-subject binary",
    "MIX_LONG_MTT. Sex, age, renal and hepatic laboratory markers and ABCB1 /",
    "CYP3A4 / CYP3A5 / CYP2C9 / CYP2C19 / CYP2J2 genotypes were screened and",
    "none were retained. Inter-individual variability is estimated on CL/F,",
    "Vc/F, Q/F, Vp/F, MTT and F (none on ka, which NONMEM fixed to zero);",
    "residual variability is purely proportional (26.2%), the additive term",
    "having been fixed to zero.",
    sep = " "
  )
  reference <- paste(
    "Fimbo AM, Mlugu EM, Kitabi EN, Kulwa GS, Iwodyah MA, Mnkugwe RH,",
    "Kunambi PP, Malishee A, Kamuhabwa AAR, Minzi OM, Aklillu E (2023).",
    "Population pharmacokinetics of ivermectin after mass drug administration",
    "in lymphatic filariasis endemic communities of Tanzania.",
    "CPT Pharmacometrics Syst Pharmacol 12(12):1884-1896.",
    "doi:10.1002/psp4.13038.",
    "Structural model adapted, via the NONMEM PRIOR subroutine, from",
    "Duthaler U et al. (2019) Br J Clin Pharmacol 85(3):626-633; all",
    "parameter values encoded here are Fimbo 2023's own final estimates",
    "(Table 2), not Duthaler's priors.",
    sep = " "
  )
  vignette <- "Fimbo_2023_ivermectin"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ng/mL"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed (single-dose study). Allometric scaling is applied to",
        "CL/F (exponent 0.75) and Vc/F (exponent 1) ONLY, both referenced to",
        "70 kg; Q/F and Vp/F carry no weight effect. Both exponents were",
        "FIXED, not estimated (Fimbo 2023 Results, 'Population",
        "pharmacokinetic model': 'Allometric scaling was applied to clearance",
        "(CL/F) and volume (Vc/F) with fixed allometric exponents, i.e.,",
        "CLi = CLpop x (WT/70)^0.75, and Vci = Vcpop x (WT/70)'; Data S3",
        "control stream 'TVCL = THETA(1) * (WT/70)**0.75' and 'TVV2 =",
        "THETA(2) * (WT/70)'). Weight spans roughly 15-90 kg; the study",
        "excluded individuals below 15 kg, who do not take part in the MDA",
        "programme. Mean (SD) weight per dose pool is 19.92 (3.28) kg at 3",
        "mg, 25.32 (3.97) kg at 6 mg, 52.26 (15.17) kg at 9 mg and 61.73",
        "(10.79) kg at 12 mg (Fimbo 2023 Table 1), so weight and dose pool",
        "are strongly confounded by construction of the height-based dosing",
        "poles."
      ),
      source_name        = "WT"
    ),
    DOSE_3MG = list(
      description        = paste(
        "1 = the subject is in the 3 mg ivermectin dose pool (height 90-119",
        "cm); 0 = the subject is in the 6, 9, or 12 mg pool. Subject-level",
        "(single-dose study, one dose pool per subject)."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the pooled 6, 9, and 12 mg dose pools, 456 of the 468 subjects)",
      notes              = paste(
        "The only covariate retained in the final model. Relative",
        "bioavailability is 48.22% higher in the 3 mg pool:",
        "F = (1 + 0.4822)^DOSE_3MG (Fimbo 2023 Table 2 row 'Proportional",
        "increase in bioavailability for the 3 mg dose' = 0.4822, RSE 40%,",
        "bootstrap median 0.48, 95% CI 0.06-1.28; Data S3 control stream",
        "'TVBIO = 1 * (1 + THETA(10))**DOSE3MG' with DOSE3MG = 1 when",
        "DOSEGRP == 3). Fimbo 2023 Results: dose was treated as a CATEGORICAL",
        "covariate rather than a continuous one because the ETA_F1 vs DOSE",
        "plot showed a step function -- the medians for the 6, 9, and 12 mg",
        "pools lie on a horizontal line and only the 3 mg pool sits above it.",
        "Only 12 of 468 subjects are in the 3 mg pool, which is why the RSE",
        "is 40% and the bootstrap CI is wide. For simulation of the routine",
        "MDA population set DOSE_3MG = 0 unless the 3 mg (smallest-children)",
        "pool is specifically being reproduced."
      ),
      source_name        = "DOSE3MG (derived in the control stream from DOSEGRP == 3)"
    ),
    MIX_LONG_MTT = list(
      description        = paste(
        "1 = the subject is classified to the long-mean-transit-time",
        "subpopulation (MTT 97.0% longer than the rest of the population);",
        "0 = the subject is classified to the reference MTT subpopulation."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the reference-MTT majority class, 83.9% of the source cohort)",
      notes              = paste(
        "Not a measured clinical covariate -- this is the per-subject latent",
        "class index from Fimbo 2023's NONMEM $MIX block (Data S3 control",
        "stream: 'NSPOP=2 / P(1) = THETA(11) / P(2) = 1 - P(1)', then",
        "'MIXPOP1 = 0 / IF(MIXNUM.EQ.1) MIXPOP1=1 / MTTMIX = (1 +",
        "THETA(12))**MIXPOP1' and 'TVMTT = THETA(5) * MTTMIX'). The",
        "population probability of MIX_LONG_MTT = 1 is THETA(11) = 0.1606",
        "(Fimbo 2023 Table 2 row 'Proportion of a subpopulation with a",
        "different typical MTT', RSE 30%, bootstrap median 0.15, 95% CI",
        "0.07-0.30), and the MTT multiplier is (1 + 0.9695) = 1.9695 (Table 2",
        "row 'Proportional difference in MTT between the subpopulation and",
        "the rest of the population'). For typical-value simulation set",
        "MIX_LONG_MTT = 0 (the majority phenotype). For population",
        "simulation draw MIX_LONG_MTT ~ Bernoulli(0.1606) per subject.",
        "Adding the mixture improved the fit by dOFV = -17 (Fimbo 2023",
        "Results and Table S3 run15 vs run11). Fimbo 2023 offers no",
        "mechanistic explanation for the two classes; the label is latent."
      ),
      source_name        = "MIXNUM / MIXEST (NONMEM $MIX class index; MIX_LONG_MTT = as.integer(MIXNUM == 1))"
    )
  )

  # Covariates that Fimbo 2023 screened but did NOT retain in the final model.
  # Documented here for provenance; they are deliberately absent from
  # covariateData because they are never referenced in model().
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened via covariate-vs-eta plots (Fimbo 2023 Table 1 lists sex",
        "among the tested covariates) and not retained. Fimbo 2023",
        "Discussion: 'Baseline characteristics did not show significant",
        "influence on IVM PKs'."
      )
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened (Fimbo 2023 Table 1) and not retained."
    ),
    EGFR = list(
      description = "Estimated glomerular filtration rate",
      units       = "mL/min/1.73m^2",
      type        = "continuous",
      notes       = paste(
        "Screened (Fimbo 2023 Table 1 'eGFR category'; control-stream $INPUT",
        "columns EGFRMD / EGFRKD / EGFR) and not retained. Ivermectin is",
        "hepatically cleared, so a renal effect was not expected."
      )
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "units/L",
      type        = "continuous",
      notes       = "Screened (Fimbo 2023 Table 1 'ALT category') and not retained."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "units/L",
      type        = "continuous",
      notes       = "Screened (Fimbo 2023 Table 1 'AST category') and not retained."
    ),
    SNP_ABCB1_RS1045642 = list(
      description = "ABCB1 c.3435C>T genotype",
      units       = "(count of variant alleles)",
      type        = "categorical",
      notes       = paste(
        "Formally tested on F1 in Fimbo 2023 Table S3 run17: coefficient",
        "-0.1741 with RSE 89%, dOFV = -1.03 vs run15, not significant and",
        "not retained."
      )
    ),
    SNP_ABCB1_RS3842 = list(
      description = "ABCB1 rs3842 genotype",
      units       = "(count of variant alleles)",
      type        = "categorical",
      notes       = paste(
        "Formally tested on F1 in Fimbo 2023 Table S3 run19: coefficient",
        "0.05182 with RSE 156%, dOFV = -0.44 vs run15, not significant and",
        "not retained."
      )
    ),
    SNP_CYP2C9 = list(
      description = "CYP2C9 *2 / *3 genotype",
      units       = "(categorical)",
      type        = "categorical",
      notes       = paste(
        "Formally tested on F1 in Fimbo 2023 Table S3 run16: coefficient",
        "0.1828 with RSE 121%, dOFV = -0.81 vs run15, not significant and",
        "not retained."
      )
    ),
    SNP_CYP2J2 = list(
      description = "CYP2J2 *7 genotype",
      units       = "(categorical)",
      type        = "categorical",
      notes       = paste(
        "Formally tested on MTT in Fimbo 2023 Table S3 run18: coefficient",
        "-0.2425 with RSE 53%, dOFV = -3.27 vs run15 -- short of the 3.84",
        "retention threshold, so not retained."
      )
    )
  )

  compartmentData <- list(
    depot = list(
      analyte  = "ivermectin",
      units    = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte  = "ivermectin",
      units    = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte  = "ivermectin",
      units    = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 468,
    n_studies      = 1,
    n_observations = 1404,
    age_range      = "5-78 years",
    weight_range   = "roughly 15-90 kg (individuals below 15 kg are excluded from the MDA programme)",
    sex_female_pct = 40.1,
    race_ethnicity = c(Black = 100),
    disease_state  = paste(
      "MDA-eligible residents of lymphatic filariasis endemic communities;",
      "participants were prescreened for circulating filarial antigen and",
      "microfilaraemia but were treated regardless of infection status, so",
      "the cohort is the general at-risk population rather than a",
      "confirmed-infected one. Pregnant women and children under 5 years were",
      "excluded."
    ),
    dose_range     = paste(
      "single oral ivermectin 3, 6, 9, or 12 mg assigned by height pole",
      "(90-119, 119-139, 139-159, and >159 cm respectively; roughly 150-200",
      "ug/kg), co-administered with albendazole 400 mg"
    ),
    regions        = "Mkinga district, Tanga region, north-eastern Tanzania",
    notes          = paste(
      "Baseline demographic, clinical and biochemical characteristics are in",
      "Fimbo 2023 Table 1, stratified by dose pool (3 mg n = 12, 6 mg n = 49,",
      "9 mg n = 159, 12 mg n = 248); genotype characteristics are in Table",
      "S2. Sex was recorded for 466 of 468 subjects (279 male, 187 female);",
      "the 40.1% female figure is 187/466. Weight and dose pool are strongly",
      "confounded because dosing was by height pole. NOTE: Fimbo 2023",
      "Methods states participants were 'aged between 5- and 78 years old',",
      "but the per-weight-band age ranges in Table S1 extend to 91 years;",
      "the Methods range is recorded here and the discrepancy is flagged in",
      "the validation vignette."
    )
  )

  ini({
    # =========================================================================
    # Structural disposition parameters. All are Fimbo 2023's own final
    # estimates for a 70 kg individual, from Table 2 (column 'Estimates
    # (RSE)'), reproduced identically in Data S2 Table S3 row 'run15'.
    # =========================================================================
    lcl <- log(7.698)
    label("Apparent oral clearance CL/F at 70 kg (L/h)")
    # Fimbo 2023 Table 2 row 'CL (L/h)' = 7.698 (RSE 7%); bootstrap median 7.7 (95% CI 7.59-8.09)

    lvc <- log(146.1)
    label("Apparent central volume of distribution Vc/F at 70 kg (L)")
    # Fimbo 2023 Table 2 row 'Vc (L)' = 146.1 (RSE 11%); bootstrap median 144.38 (95% CI 10.79-168.02)

    lq <- log(20.42)
    label("Apparent intercompartmental clearance Q/F (L/h; not weight-scaled)")
    # Fimbo 2023 Table 2 row 'Q (L/h)' = 20.42 (RSE 12%); bootstrap median 20.13 (95% CI 15.59-23.63)

    lvp <- log(207.1)
    label("Apparent peripheral volume of distribution Vp/F (L; not weight-scaled)")
    # Fimbo 2023 Table 2 row 'Vp (L)' = 207.1 (RSE 11%); bootstrap median 208.72 (95% CI 199.35-249)

    # =========================================================================
    # Absorption. Savic analytical transit chain feeding an absorption
    # (depot) compartment that empties into central at ka.
    # =========================================================================
    lmtt <- log(1.523)
    label("Mean transit time MTT through the NN-compartment Savic chain (h), reference subpopulation")
    # Fimbo 2023 Table 2 row 'MTT (h)' = 1.523 (RSE 5%); bootstrap median 1.53 (95% CI 1.33-1.98)

    lka <- log(0.7082)
    label("First-order absorption rate constant ka from the depot to the central compartment (1/h)")
    # Fimbo 2023 Table 2 row 'K a (/h)' = 0.7082 (RSE 17%); bootstrap median 0.7 (95% CI 0.15-0.91)

    nn_fix <- fixed(6)
    label("Number of Savic transit compartments NN (integer, unitless)")
    # Fimbo 2023 Table 2 row 'NN' = '6 (Fixed to this value)'; Data S3 control stream '$THETA (6) FIX ; N' and 'LNFAC = 6.579251 ; log NN factorial, NN=number of transit compartments [log(6!)]'

    lfdepot <- fixed(log(1))
    label("Reference relative oral bioavailability F (unitless; 1 because no intravenous data were available)")
    # Fimbo 2023 Results, 'Population pharmacokinetic model': 'we estimated the IIV for bioavailability (IIV-F1), whereas F1 was fixed to 1'; Data S3 control stream 'TVBIO = 1 * (1 + THETA(10))**DOSE3MG'

    # =========================================================================
    # Allometric exponents. Both FIXED by the authors, not estimated, and
    # applied to CL/F and Vc/F only.
    # =========================================================================
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent on CL/F, reference 70 kg (unitless)")
    # Fimbo 2023 Results: 'Allometric scaling was applied to clearance (CL/F) and volume (Vc/F) with fixed allometric exponents, i.e., CLi = CLpop x (WT/70)^0.75'; Data S3 'TVCL = THETA(1) * (WT/70)**0.75'

    e_wt_vc <- fixed(1)
    label("Allometric exponent on Vc/F, reference 70 kg (unitless)")
    # Fimbo 2023 Results: '... and Vci = Vcpop x (WT/70)' -- the printed equation carries no exponent on the Vc term, i.e. exponent 1; Data S3 'TVV2 = THETA(2) * (WT/70)'

    # =========================================================================
    # Covariate effects.
    # =========================================================================
    e_dose_3mg_fdepot <- 0.4822
    label("Proportional increase in relative bioavailability F in the 3 mg dose pool (fraction)")
    # Fimbo 2023 Table 2 row 'Proportional increase in bioavailability for the 3 mg dose' = 0.4822 (RSE 40%); bootstrap median 0.48 (95% CI 0.06-1.28). Enters as F = (1 + 0.4822)^DOSE_3MG per Data S3 'TVBIO = 1 * (1 + THETA(10))**DOSE3MG'

    e_mix_long_mtt_mtt <- 0.9695
    label("Proportional increase in MTT in the long-transit-time mixture subpopulation (fraction)")
    # Fimbo 2023 Table 2 row 'Proportional difference in MTT between the subpopulation and the rest of the population' = 0.9695 (RSE 12%); bootstrap median 1.02 (95% CI 0.8-1.28). Enters as MTT = MTT_pop * (1 + 0.9695)^MIX_LONG_MTT per Data S3 'MTTMIX = (1+THETA(12))**MIXPOP1'

    # =========================================================================
    # Inter-individual variability.
    #
    # Fimbo 2023 Table 2 reports these as 'Interindividual variability for
    # <par> (%CV)' with the footnote 'Inter-individual variabilities are on a
    # standard deviation scale (%CV)'. The Methods separately define
    # CV% = 100 * sqrt(exp(omega^2) - 1), which would make the printed number
    # a CV rather than an omega. Data S3 settles it: the control stream's
    # $OMEGA initial estimates are exactly the SQUARES of the corresponding
    # Table S3 run10 entries (0.201 -> 0.4479, 0.0758 -> 0.2753, 0.247 ->
    # 0.4971, 0.172 -> 0.4147, 0.16 -> 0.4002, 0.0922 -> 0.3036). The printed
    # values are therefore omega on the log (eta) scale, NOT the CV, and the
    # nlmixr2 variances below are the printed values squared.
    # =========================================================================
    etalcl ~ 0.198738
    # Fimbo 2023 Table 2 row 'Interindividual variability for CL (%CV)' = 0.4458 (RSE 15%); variance = 0.4458^2

    etalvc ~ 0.111489
    # Fimbo 2023 Table 2 row 'Interindividual variability for Vc (%CV)' = 0.3339 (RSE 21%); variance = 0.3339^2

    etalq ~ 0.235807
    # Fimbo 2023 Table 2 row 'Interindividual variability for Q (%CV)' = 0.4856 (RSE 11%); variance = 0.4856^2

    etalvp ~ 0.170156
    # Fimbo 2023 Table 2 row 'Interindividual variability for Vp (%CV)' = 0.4125 (RSE 16%); variance = 0.4125^2

    etalmtt ~ 0.067029
    # Fimbo 2023 Table 2 row 'Interindividual variability for MTT (%CV)' = 0.2589 (RSE 11%); variance = 0.2589^2

    etalfdepot ~ 0.077284
    # Fimbo 2023 Table 2 row 'Interindividual variability for F (%CV)' = 0.278 (RSE 15%); variance = 0.278^2; shrinkage 22.4% per Fimbo 2023 Results

    # No IIV on ka: Data S3 control stream '$OMEGA BLOCK(1) 0 FIX ; EKA'.

    # =========================================================================
    # Residual unexplained variability.
    # =========================================================================
    propSd <- 0.262
    label("Proportional residual error (fraction)")
    # Fimbo 2023 Table 2 row 'Proportional residual error' = 0.262 (RSE 3%); bootstrap median 0.26 (95% CI 0.24-0.3). The additive term is absent by construction: Data S3 declares '$THETA (0) FIX ; ADD' and forms W = SQRT(ADD**2 + PROP**2*IPRED**2), so the error model is purely proportional.
  })

  model({
    # --- 1. Individual PK parameters ---------------------------------------
    # Allometry is applied to CL/F and Vc/F only; Q/F and Vp/F carry no
    # weight effect. This asymmetry is deliberate in the source model and is
    # what makes the terminal half-life fall as body weight rises.
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q  <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    # Mean transit time carries the two-class mixture multiplier.
    mtt <- exp(lmtt + etalmtt) * (1 + e_mix_long_mtt_mtt)^MIX_LONG_MTT

    # No IIV on ka (NONMEM fixed the corresponding OMEGA to zero).
    ka <- exp(lka)

    nn <- nn_fix

    # Relative bioavailability is raised in the 3 mg dose pool.
    fdepot <- exp(lfdepot + etalfdepot) * (1 + e_dose_3mg_fdepot)^DOSE_3MG

    # --- 2. Micro-constants -------------------------------------------------
    kel <- cl / vc
    k23 <- q  / vc
    k32 <- q  / vp

    # --- 3. ODE system ------------------------------------------------------
    # Savic analytical transit absorption. The rxode2 builtin
    # transit(nn, mtt, fdepot) emits the gamma-density input rate to the
    # depot from the most recent dose, with ktr = (nn + 1) / mtt computed
    # internally; this is the exact analogue of the Data S3 $DES line
    #   DADT(1) = EXP(LOG(BIO*PODO+X)+LOG(KTR+X)+NN*LOG(KTR*T+X)-KTR*T-LNFAC)
    #             - KA*A(1)
    # with LNFAC = log(NN!) and KTR = (NN+1)/MTT defined in $PK. The depot
    # then empties into central at first-order rate ka.
    d/dt(depot)       <- transit(nn, mtt, fdepot) - ka * depot
    d/dt(central)     <- ka * depot - kel * central - k23 * central + k32 * peripheral1
    d/dt(peripheral1) <-                             k23 * central - k32 * peripheral1

    # --- 4. Bioavailability -------------------------------------------------
    # Suppress the standard dose bolus on depot so the analytical transit()
    # chain is the only input pathway. This is the direct analogue of the
    # Data S3 $PK line 'F1 = 0' whose comment reads: ';;Prevent NONMEM from
    # administering dose into absorption compartment as this will be done
    # through transit compartments'.
    f(depot) <- 0

    # --- 5. Observation and residual error ----------------------------------
    # Doses are in mg and volumes in L, so central/vc is mg/L; the factor of
    # 1000 converts to the reported ng/mL. This is the Data S3 scaling
    # 'S2 = V2/1000 ; dose in mg, DV in ng/mL' with IPRED = A(2)/S2.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
