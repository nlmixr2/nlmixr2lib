Wang_2022_aripiprazole <- function() {
  description <- "Population PK model for aripiprazole after oral tablet dosing and after intramuscular injection of aripiprazole once monthly (AOM, aripiprazole monohydrate extended-release injectable suspension) in healthy adults and in adults with schizophrenia or schizoaffective disorder (Wang 2022). Three-compartment linear disposition fed by two parallel absorption routes: oral doses enter `depot` as a zero-order input of fixed rate R1 = 9.33 mg/h and are then absorbed first-order with Ka = 0.540 /h (the paper's 'sigmoid absorption', which produces a short effective lag followed by a steeper rise in absorbed amount), while AOM doses enter `depot2` and are absorbed first-order with a much slower rate constant so that the AOM disposition is absorption-rate limited (average absorption half-life about 28 days versus an average terminal elimination half-life of 7.5 days). AOM bioavailability is 48% higher than oral. The AOM absorption rate constant decreases as a power function of body mass index and is 34.6% higher in men, so the AOM absorption half-life at BMI 28 kg/m^2 is about 32 days in women and 24 days in men. Apparent oral clearance is halved in CYP2D6 poor metabolizers, reduced 51.1% by a concomitant CYP2D6 inhibitor and 23.7% by a concomitant CYP3A4 inhibitor. The proportional residual error magnitude switches between the rich-sampling phase 1 studies and the sparse-sampling phase 3 study. The companion time-to-relapse exposure-response model from the same paper is modellib('Wang_2022_aripiprazole_relapse')."
  reference <- paste(
    "Wang X, Raoufinia A, Bihorel S, Passarell J, Mallikaarjun S, Phillips L.",
    "Population Pharmacokinetic Modeling and Exposure-Response Analysis for",
    "Aripiprazole Once Monthly in Subjects With Schizophrenia.",
    "Clin Pharmacol Drug Dev. 2022;11(2):150-164. doi:10.1002/cpdd.1022.",
    sep = " "
  )
  vignette <- "Wang_2022_aripiprazole"

  # Study-stratified proportional residual magnitudes; the residual actually
  # applied to Cc is the canonical `propSd`, selected inside model() from these
  # two paper-specific SDs. Same shape as Chan_2023_nirmatrelvir.R and
  # Friberg_2012_voriconazole.R.
  paper_specific_residual_sds <- c("propSdPhase1", "propSdPhase3")

  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Oral doses enter `depot` and AOM (intramuscular) doses enter `depot2`.
  # buildModelDb()'s dosing heuristic only looks for states literally named
  # `depot` and `central`, so without this field the registry would report
  # `depot,central` -- true for `depot`, wrong for `central` (never dosed) and
  # silently missing the AOM route entirely.
  dosing <- c("depot", "depot2")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = TRUE means checked against the source paper
  # (Wang 2022 Table 1 and Results, "Base Model Development").
  compartmentData <- list(
    depot       = list(analyte = "aripiprazole", units = "mg", specimen = "administration site",   verified = TRUE),
    depot2      = list(analyte = "aripiprazole", units = "mg", specimen = "administration site",   verified = TRUE),
    central     = list(analyte = "aripiprazole", units = "mg", specimen = "plasma",                               verified = TRUE),
    peripheral1 = list(analyte = "aripiprazole", units = "mg", specimen = "not applicable",                       verified = TRUE),
    peripheral2 = list(analyte = "aripiprazole", units = "mg", specimen = "not applicable",                       verified = TRUE)
  )

  covariateData <- list(
    BMI = list(
      description        = "Body mass index.",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required input for AOM (intramuscular) dosing; enters the AOM absorption rate constant as a power function normalised to 28 kg/m^2 (Wang 2022 Table 1 footnote, 'Related equations': AOM Ka = 0.000904 * (BMI/28)^-0.975 * (1 + 0.346 * Male)). The exponent is negative, so the AOM absorption rate constant falls and the AOM absorption half-life lengthens as BMI increases. Observed BMI range in the popPK analysis population was 15 to 61 kg/m^2 (Wang 2022 Discussion). BMI was not retained on any disposition parameter, and the paper reports that the model-predicted steady-state exposures showed no trend with BMI across that range. Time-fixed per subject in this encoding.",
      source_name        = "BMI"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (female) -- see notes; the source parameterises the shift on men",
      notes              = "The source paper's covariate equation is written on a male indicator: AOM Ka = 0.000904 * (BMI/28)^-0.975 * (1 + 0.346 * Male) (Wang 2022 Table 1 footnote). The canonical column is SEXF (1 = female), so the effect is applied here as (1 + e_sexf_ka_im * (1 - SEXF)), i.e. women carry the reference AOM absorption rate constant and men are 34.6% faster. This reproduces the paper's own reported AOM absorption half-lives at BMI 28 kg/m^2 exactly: 0.693/0.000904 = 766.6 h = 31.9 days for women (paper: about 32 days) and 0.693/(0.000904 * 1.346) = 569.7 h = 23.7 days for men (paper: about 24 days). Same (1 - SEXF) inversion idiom as Sathe_2024_sacituzumab.R, Kuchimanchi_2024_dostarlimab.R and Zufferey_2018_fondaparinux.R. Sex was screened on oral Ka, CL/F and Vc/F and was not retained on any of them.",
      source_name        = "Male"
    ),
    CYP2D6_PM = list(
      description        = "CYP2D6 poor-metabolizer phenotype indicator, 1 = poor metabolizer, 0 = extensive metabolizer.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP2D6 extensive metabolizer)",
      notes              = "Selects between the two separately estimated apparent oral clearances, CL/F = 3.71 L/h for extensive metabolizers and 1.88 L/h for poor metabolizers (Wang 2022 Table 1; Table 1 footnote equation CL/F = (3.71 * EM + 1.88 * PM) * (1 - 0.511 * CYP2D6) * (1 - 0.237 * CYP3A4)). Only two phenotype levels appear in the source model, so the companion CYP2D6_EM / CYP2D6_IM / CYP2D6_UM indicators are not needed: EM is simply CYP2D6_PM = 0. In the source analysis, phenotype was measured for part of the cohort and imputed for the remainder with a NONMEM mixture model that assigned the most probable status assuming 90% extensive metabolizers in the general population (Wang 2022 Methods, 'Population Pharmacokinetic Model Development'); the mixture step is an estimation device for subjects with unknown status and is not part of the final structural model, so this extraction carries the phenotype as an ordinary observed covariate. Of the 663 subjects in the model development data set, 621 were classified extensive and 42 poor metabolizers (Wang 2022 Results, 'Model Simulations'). Time-fixed per subject (germline genotype-derived phenotype).",
      source_name        = "PM / EM"
    ),
    CONMED_CYP2D6_INH = list(
      description        = "Concomitant CYP2D6 inhibitor coadministration indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant CYP2D6 inhibitor)",
      notes              = "Multiplies apparent oral clearance by (1 - 0.511), i.e. a 51.1% reduction, matching the paper's statement that CL/F in the presence of strong CYP2D6 inhibitors is approximately half that of CYP2D6 extensive metabolizers (Wang 2022 Abstract and Table 1). The = 1 category is populated exclusively by the STRONG CYP2D6 inhibitors administered in the dedicated phase 1 drug-interaction studies (Wang 2022 Discussion: 'phase 1 trials designed to evaluate the effect of CYP2D6 and CYP3A4 inhibitors'); the paper does not name the individual agents, so the strength tier of the pooled category is 'strong only' and no weaker inhibitors contribute. The class-level canonical is used rather than a hypothetical CONMED_CYP2D6_INH_STRONG because the source's own model term and Table 1 row are unstratified ('proportional change in CL/F for CYP2D6 inhibitor'). Distinct from CYP2D6_PM, which is intrinsic phenotype rather than a drug-drug interaction; the two enter multiplicatively, so a poor metabolizer on a CYP2D6 inhibitor gets both reductions. Per-record time-varying in a fixed-sequence interaction study. In the model development cohort, 13 CYP2D6 extensive metabolizers took a CYP2D6 inhibitor (Wang 2022 Results, 'Model Simulations').",
      source_name        = "CYP2D6"
    ),
    CONMED_CYP3A4_INH = list(
      description        = "Concomitant CYP3A4 inhibitor coadministration indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant CYP3A4 inhibitor)",
      notes              = "Multiplies apparent oral clearance by (1 - 0.237), i.e. a 23.7% reduction, matching the paper's statement that CL/F in the presence of strong CYP3A4 inhibitors is about 24% lower than in CYP2D6 extensive metabolizers (Wang 2022 Abstract, Results and Table 1). As for CONMED_CYP2D6_INH, the = 1 category comprises only the STRONG CYP3A4 inhibitors used in the dedicated phase 1 drug-interaction studies; the paper does not name them. The class-level canonical is used rather than CONMED_CYP3A4_INH_STRONG so that the two inhibitor covariates in this model are encoded symmetrically and match the unstratified Table 1 row labels. The effect is consistent with the independent Koue 2007 oral-aripiprazole popPK analysis cited in the Wang 2022 Discussion, which found a 14% CL/F reduction with itraconazole coadministration. In the model development cohort, 25 CYP2D6 extensive metabolizers took a CYP3A4 inhibitor (Wang 2022 Results, 'Model Simulations').",
      source_name        = "CYP3A4"
    ),
    STUDY_ARI_PHASE3 = list(
      description        = "Study-phase indicator for the observation record, 1 = the sparse-sampling phase 3 study 31-07-246, 0 = one of the four rich-sampling phase 1 studies (31-98-206, 31-98-207, CN138020, 31-05-244).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (phase 1 observation record)",
      notes              = "Record-level study-design property rather than a subject-level covariate: it selects the proportional residual-error magnitude only (24.23% for phase 1, 28.11% for phase 3; Wang 2022 Table 1, 'Phase 1 RV (%CV)' and 'Phase 3 RV (%CV)'). The phase 3 study contributed only predose samples plus single samples on days 7, 14 and 28 after dosing (Wang 2022 Results, 'Base Model Development'), which is the usual reason a sparse arm carries the larger residual. Follows the STUDY_<drug>_<phase> convention established by STUDY_NMV_PHASE23, STUDY_NIPOCALIMAB_PHASE1, STUDY_FARLETUZUMAB_PHASE2 and STUDY_POSA_PHASE3. Not collinear with route: the phase 1 studies contributed both oral and AOM records, and the phase 3 study contributed both oral lead-in and AOM records.",
      source_name        = "study phase"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 663L,
    n_studies      = 5L,
    n_observations = "6153 aripiprazole plasma concentration records in the model development data set (Wang 2022 Methods, 'Population PK Analysis Source Data').",
    age_range      = "Adults; the per-study age distribution is in Wang 2022 Table S2 (supplement not on disk -- see the vignette Errata).",
    weight_range   = "Not reported in the main text; per-study body-weight summaries are in Wang 2022 Table S2 (supplement not on disk).",
    bmi_range      = "15 to 61 kg/m^2 in the popPK analysis population (Wang 2022 Discussion, 'Population PK Analysis').",
    sex_female_pct = NA_real_,
    race_ethnicity = "Self-reported race category (interpreted per US FDA guidance) was screened as a covariate on CL/F and was not retained. The distribution is in Wang 2022 Table S3 (supplement not on disk).",
    disease_state  = "52 healthy subjects (oral aripiprazole only) and 611 subjects with schizophrenia or schizoaffective disorder (Wang 2022 Methods).",
    dose_range     = "Oral aripiprazole tablets and AOM intramuscular injections. The model development data set contained no AOM dose below 300 mg; the phase 3 arms used 400 mg AOM with one permitted reduction to 300 mg. Doses of 50/25 mg AOM (study 31-07-247) were shown NOT to be adequately described by this model and were excluded from the source paper's own downstream analyses.",
    regions        = "Not reported in the main text.",
    cyp2d6_status  = "621 CYP2D6 extensive metabolizers and 42 poor metabolizers; 13 extensive metabolizers were taking a CYP2D6 inhibitor and 25 a CYP3A4 inhibitor (Wang 2022 Results, 'Model Simulations').",
    notes          = "Pooled from five studies: phase 1 studies 31-98-206, 31-98-207, CN138020 and 31-05-244 (serial PK sampling) and phase 3 study 31-07-246 (sparse PK sampling). All AOM doses were given in the gluteus maximus except in CN138020 (13 subjects) where they were given in the nondominant arm or the midlateral thigh; injection site was screened on AOM Ka and was not retained. A sixth study, phase 3 study 31-07-247, was used only for external validation. The final model met the external-validation criteria for the 400/300 mg AOM arm (median %PPE -6.8%, median |%PPE| 29.2%) but not for the 50/25 mg arm (median |%PPE| 41.3%, 75th percentile 74.6%), so this model must NOT be extrapolated to AOM doses below 300 mg. Model development used NONMEM VI level 2.0 with FOCE-I; all analyses were run in 2011. Minimum objective function value 48892.907."
  )

  ini({
    # ==================================================================
    # All values from Wang 2022 Table 1 ("Final Model - Phase 1 and
    # Phase 3 Oral + AOM Data") and its footnote "Related equations".
    #
    # The final model estimated only relative bioavailability for AOM,
    # AOM Ka, CL/F, Vc/F, the related IIV terms and the residual-error
    # terms; every other parameter was carried over fixed from the base
    # model built before the sparse phase 3 data were added (Wang 2022
    # Results, "Final PopPK Model"). Those carried-over parameters are
    # wrapped in fixed() here and are shown as "Fixed" in the Table 1
    # %RSE column.
    # ==================================================================

    # ---- Oral absorption: zero-order input into depot, then first-order ----
    # The paper calls this "sigmoid absorption": the oral dose is
    # delivered into `depot` at a constant rate R1 and is simultaneously
    # absorbed first-order into `central`, which delays the appearance of
    # measurable concentrations and gives an S-shaped absorbed-amount
    # curve. Applied via rate(depot) <- r1 in model(); dose records must
    # carry rate = -1 so rxode2 uses the modelled rate. At a 10 mg oral
    # dose the zero-order input therefore lasts 10 / 9.33 = 1.07 h.
    # Naming follows Brekkan_2018_pegfilgrastim.R, the existing model
    # that uses the same zero-order-into-depot construction.
    lr1        <- fixed(log(9.33))     ; label("Zero-order input rate R1 of the oral dose into the depot compartment (mg/h)")   # Wang 2022 Table 1: R1 = 9.33 mg/h (Fixed)
    lka_oral   <- fixed(log(0.540))    ; label("First-order absorption rate constant out of the oral depot (1/h)")              # Wang 2022 Table 1: Ka oral first-order absorption rate = 0.540 /h (Fixed)

    # ---- AOM (intramuscular) absorption ----
    lka_im     <- log(0.000904)        ; label("First-order absorption rate constant out of the AOM injection-site depot at BMI 28 kg/m^2 in women (1/h)")  # Wang 2022 Table 1: IM Ka AOM first-order absorption rate = 0.000904 /h (%RSE 5.3)
    e_bmi_ka_im <- -0.975              ; label("Power exponent on (BMI/28) for the AOM absorption rate constant (unitless)")     # Wang 2022 Table 1: IM Ka power for (BMI/28) = -0.975 (%RSE 11.5)
    e_sexf_ka_im <- 0.346              ; label("Proportional increase in the AOM absorption rate constant in men, applied on (1 - SEXF) (unitless)")  # Wang 2022 Table 1: IM Ka proportional shift for men = 0.346 (%RSE 28.9)
    lfdepot_im <- log(1.48)            ; label("Bioavailability of the AOM injection relative to oral dosing (unitless)")        # Wang 2022 Table 1: F2 relative bioavailability for AOM = 1.48 (%RSE 4.9)

    # ---- Disposition ----
    # CL/F is estimated separately in each CYP2D6 phenotype stratum, so
    # both estimates carry a stratum suffix and neither keeps the bare
    # `lcl` name (parameter-names.md, "Stratum-suffixed parameters";
    # same shape as Shen_2024_vancomycin.R). Oral bioavailability is the
    # implicit anchor F = 1, so these are apparent (CL/F) values.
    lcl_em     <- log(3.71)            ; label("Apparent oral clearance in CYP2D6 extensive metabolizers (L/h)")                # Wang 2022 Table 1: CL/F clearance for EM = 3.71 L/h (%RSE 4.0)
    lcl_pm     <- log(1.88)            ; label("Apparent oral clearance in CYP2D6 poor metabolizers (L/h)")                     # Wang 2022 Table 1: CL/F clearance for PM = 1.88 L/h (%RSE 6.9)
    e_cyp2d6_inh_cl <- fixed(-0.511)   ; label("Proportional change in apparent oral clearance with a concomitant CYP2D6 inhibitor (unitless)")  # Wang 2022 Table 1: CL/F proportional change for CYP2D6 inhibitor = -0.511 (Fixed)
    e_cyp3a4_inh_cl <- fixed(-0.237)   ; label("Proportional change in apparent oral clearance with a concomitant CYP3A4 inhibitor (unitless)")  # Wang 2022 Table 1: CL/F proportional change for CYP3A4 inhibitor = -0.237 (Fixed)
    lvc        <- log(93.4)            ; label("Apparent central volume of distribution (L)")                                   # Wang 2022 Table 1: Vc/F central volume = 93.4 L (%RSE 8.8)
    lq         <- fixed(log(0.591))    ; label("Apparent intercompartmental clearance to the first peripheral compartment (L/h)")  # Wang 2022 Table 1: Q1/F intercompartmental CL/F = 0.591 L/h (Fixed)
    lvp        <- fixed(log(118))      ; label("Apparent first peripheral volume of distribution (L)")                          # Wang 2022 Table 1: Vp1/F peripheral volume = 118 L (Fixed)
    lq2        <- fixed(log(28.8))     ; label("Apparent intercompartmental clearance to the second peripheral compartment (L/h)")  # Wang 2022 Table 1: Q2/F second intercompartmental CL/F = 28.8 L/h (Fixed)
    lvp2       <- fixed(log(134))      ; label("Apparent second peripheral volume of distribution (L)")                         # Wang 2022 Table 1: Vp2/F second peripheral volume = 134 L (Fixed)

    # ---- Interindividual variability ----
    # Table 1 reports IIV as a percent coefficient of variation, so the
    # internal variances are omega^2 = log(CV^2 + 1):
    #   CL      38.34 %CV -> log(0.3834^2 + 1) = 0.137146
    #   Vc     124.50 %CV -> log(1.2450^2 + 1) = 0.936103
    #   Ka oral 65.88 %CV -> log(0.6588^2 + 1) = 0.360480
    #   Ka AOM  55.59 %CV -> log(0.5559^2 + 1) = 0.269282
    # A single CL IIV term is shared by both CYP2D6 phenotype strata:
    # Table 1 reports the IIV on the "CL/F: clearance for EM" row only
    # and leaves the PM row's IIV cell blank.
    #
    # No off-diagonal (correlation) terms: after the covariates were
    # added "it was determined that off-diagonal elements of the
    # covariance matrix (describing the covariance of IIV terms) could
    # not be successfully estimated" (Wang 2022 Results, "Covariate
    # Analysis"). The r = 0.901 reported in the Table 1 footnote is the
    # correlation between the ESTIMATES of CL for extensive metabolizers
    # and of AOM relative bioavailability -- an estimation-precision
    # correlation, not an IIV covariance -- so it is deliberately not
    # encoded as a block.
    #
    # Source traces sit ABOVE each line rather than trailing it: a trailing
    # comment on an omega line is read as that eta's label, and one of these
    # traces has to carry the word Fixed from the Table 1 %RSE column.
    #
    # Wang 2022 Table 1: IIV on CL/F = 38.34 %CV (%RSE 6.9)
    etalcl      ~ 0.137146
    # Wang 2022 Table 1: IIV on Vc/F = 124.50 %CV (%RSE 15.2)
    etalvc      ~ 0.936103
    # Wang 2022 Table 1: IIV on oral Ka = 65.88 %CV, %RSE column reads Fixed
    etalka_oral ~ fixed(0.360480)
    # Wang 2022 Table 1: IIV on IM Ka = 55.59 %CV (%RSE 8.2)
    etalka_im   ~ 0.269282

    # ---- Residual variability ----
    # Two proportional residual magnitudes, selected per observation
    # record by study phase; see the STUDY_ARI_PHASE3 covariate entry.
    propSdPhase1 <- 0.2423             ; label("Proportional residual SD for phase 1 observation records (fraction)")   # Wang 2022 Table 1: Phase 1 RV = 24.23 %CV (%RSE 8.4)
    propSdPhase3 <- 0.2811             ; label("Proportional residual SD for phase 3 observation records (fraction)")   # Wang 2022 Table 1: Phase 3 RV = 28.11 %CV (%RSE 4.7)
  })

  model({
    # ---- Oral absorption ----
    r1        <- exp(lr1)
    ka_oral   <- exp(lka_oral + etalka_oral)

    # ---- AOM (intramuscular) absorption ----
    # Wang 2022 Table 1 footnote, "Related equations":
    #   AOM Ka = 0.000904 * (BMI/28)^-0.975 * (1 + 0.346 * Male)
    # The canonical SEXF column is 1 = female, so the paper's `Male`
    # indicator is (1 - SEXF). Women therefore carry the reference rate
    # constant and men are 34.6% faster.
    ka_im     <- exp(lka_im + etalka_im) *
      (BMI / 28)^e_bmi_ka_im *
      (1 + e_sexf_ka_im * (1 - SEXF))
    fdepot_im <- exp(lfdepot_im)

    # ---- Disposition ----
    # Wang 2022 Table 1 footnote, "Related equations":
    #   CL/F = (3.71 * EM + 1.88 * PM) * (1 - 0.511 * CYP2D6) * (1 - 0.237 * CYP3A4)
    # EM and PM are mutually exclusive indicators, so EM = 1 - CYP2D6_PM.
    cl  <- (exp(lcl_em) * (1 - CYP2D6_PM) + exp(lcl_pm) * CYP2D6_PM) *
      (1 + e_cyp2d6_inh_cl * CONMED_CYP2D6_INH) *
      (1 + e_cyp3a4_inh_cl * CONMED_CYP3A4_INH) *
      exp(etalcl)
    vc  <- exp(lvc + etalvc)
    q   <- exp(lq)
    vp  <- exp(lvp)
    q2  <- exp(lq2)
    vp2 <- exp(lvp2)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    d/dt(depot)       <- -ka_oral * depot
    d/dt(depot2)      <- -ka_im * depot2
    d/dt(central)     <- ka_oral * depot + ka_im * depot2 -
      (kel + k12 + k13) * central + k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # Oral doses are delivered into `depot` at the fixed zero-order rate
    # R1; rxode2 reads rate() at solve time for dose records flagged
    # rate = -1 and converts the bolus into a constant-rate input of
    # duration amt / r1. Oral bioavailability is the implicit F = 1
    # anchor, so only the AOM depot carries an f().
    rate(depot) <- r1
    f(depot2)   <- fdepot_im

    # Amounts are in mg and volumes in L, so central/vc is mg/L; the
    # factor 1000 converts to the ng/mL used throughout the paper.
    Cc <- 1000 * central / vc

    # Study-phase-dependent proportional residual magnitude
    # (Wang 2022 Table 1, "Phase 1 RV" and "Phase 3 RV" rows).
    propSd <- propSdPhase1 * (1 - STUDY_ARI_PHASE3) + propSdPhase3 * STUDY_ARI_PHASE3
    Cc ~ prop(propSd)
  })
}
