Wang_2019_oral1cmt_base <- function() {
  description <- paste(
    "Methodology reference. Base (no-covariate) data-generating",
    "one-compartment population PK model for a HYPOTHETICAL drug, from",
    "Wang and Liu 2019, the simulation study that introduced the modified",
    "missing dose method (MDM2) and the compartment initialization method",
    "(CIM) for handling missing dosing records in population PK analysis.",
    "There is no real molecule and no real patients: the authors DEFINED",
    "these population parameters and simulated 100 virtual subjects per",
    "replicate, then deleted dosing records from the simulated data sets to",
    "compare six analysis methods (OM, PDM, MDM, MDM2, CIM and the complete",
    "data ideal method IDM). Every value here is therefore an author-chosen",
    "simulation constant, not an estimate, and all are encoded with fixed().",
    "Structure is first-order absorption into a depot, one-compartment",
    "distribution and linear elimination, matching the authors' NONMEM",
    "ADVAN2 TRANS2 subroutine. Fraction absorbed is 1 (complete absorption),",
    "so no bioavailability term is applied. Residual error is combined",
    "proportional plus additive, which is nlmixr2's default combined2",
    "root-sum-square form and reproduces the authors' NONMEM error line",
    "Y = F * (1 + ERR(1)) + ERR(2) exactly, because the two SIGMA elements",
    "are independent. Doses may be given into depot (oral) or directly into",
    "central (i.v. bolus); the paper simulated both routes from this one",
    "structure. The companion model Wang_2019_oral1cmt_covariate adds the",
    "body weight and sex covariate effects the authors simulated separately.",
    sep = " "
  )
  reference <- paste(
    "Wang Y, Liu X. Handling Missing Dosing History in Population",
    "Pharmacokinetic Modeling: An Extension to MDM Method.",
    "CPT Pharmacometrics Syst Pharmacol. 2019;8(1):39-49.",
    "doi:10.1002/psp4.12374. PMCID PMC6363138.",
    "Structural parameter values transcribed from Table 1 (parameter set A)",
    "and cross-checked against the authors' deposited data-generating NONMEM",
    "control stream in Code Examples S1 (supplementary file PSP4-8-39-s003",
    ".txt, $PROB RUN# 011_1, the $SIMULATION problem), which uses exactly",
    "these values.",
    sep = " "
  )
  vignette <- "Wang_2019_missing_dosing_history"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. The drug is hypothetical, so it has no analyte name.
  # The paper states only that "blood samples were collected during visits and
  # later analyzed for drug concentrations (ng/mL)" (Simulated trial design);
  # it never says plasma, serum or whole blood, so specimen is unverified.
  compartmentData <- list(
    depot = list(
      analyte = "hypothetical drug",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "hypothetical drug",
      units = "mg",
      specimen = "plasma",
      verified = FALSE
    )
  )

  covariateData <- list()

  # Body weight and sex WERE simulated as subject descriptors for every
  # virtual cohort (Table S1), but the base model deliberately carries no
  # covariate terms -- the paper simulated and analysed "base models and
  # covariate models" separately (PK data simulation). They are documented
  # here so the cohort provenance is preserved without raising
  # unused-covariate warnings. The covariate effects live in
  # Wang_2019_oral1cmt_covariate.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Simulated as Normal with mean (SD) 70 (10) kg for men and 60 (10) kg",
        "for women (Table S1), but not entered into the base model.",
        sep = " "
      ),
      source_name = "WT"
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "categorical",
      notes = paste(
        "Balanced 50 men / 50 women in every simulated data set (PK data",
        "simulation; Table S1), but not entered into the base model. The",
        "source column SEX is a MALE indicator (the control stream comments",
        "'SEX=1 IF MALE, 0 IF FEMALE'), so SEXF = 1 - SEX.",
        sep = " "
      ),
      source_name = "SEX"
    )
  )

  population <- list(
    species = "None (methodology paper; a hypothetical drug simulated in virtual subjects, not a fit of any real molecule).",
    n_subjects = 100L,
    n_studies = 1L,
    weight_range = "Normally distributed, mean (SD) 70 (10) kg for men and 60 (10) kg for women (Table S1).",
    sex_female_pct = 50,
    disease_state = "N/A (Monte Carlo simulation study).",
    dose_range = paste(
      "Three consecutive doses at 0, 12 and 24 h (Table 2). The nominal",
      "prescribed dose is 8 mg every 12 h, but each administered dose was",
      "drawn independently from a Uniform(5, 10) mg distribution so that the",
      "true dose amount could not be reconstructed from the prescription --",
      "this is what makes the prescribed dose method (PDM) biased. Given as",
      "an oral or an i.v. bolus dose.",
      sep = " "
    ),
    regions = "N/A",
    scope_note = paste(
      "Filed under inst/modeldb/pharmacokinetics/ rather than specificDrugs/",
      "because there is no drug. This follows the precedent set by",
      "Beal_2001_iv1cmt_bql and the Schoning_2026_oral1cmt_* family, the other",
      "methodology-reference toy models in the library. The file stem uses the",
      "structural descriptor oral1cmt in the slot where a drug name would",
      "normally go, as Beal_2001 does with iv1cmt.",
      sep = " "
    ),
    notes = paste(
      "Sampling times 11.83, 12.25, 12.75, 13, 17, 23.83, 24.25, 24.75, 25,",
      "27, 31 and 35 h, giving 1,200 observations per 100-subject data set",
      "(Table 2). Table 1 lists five parameter sets A-E that differ ONLY in",
      "clearance (CL = 15, 20, 25, 30 and 35 L/h); every other value is shared",
      "across the five sets. Set A (CL = 15 L/h) is encoded here because it is",
      "the set the authors' own deposited data-generating control stream",
      "actually runs ($THETA 15 ;[CL, L/h] in Code Examples S1). Change lcl to",
      "log(20), log(25), log(30) or log(35) for sets B-E. The paper ran",
      "5 parameter sets x 2 dosing routes x 4 missing-record cases x 100",
      "replicates, i.e. 8,000 model runs per analysis method.",
      sep = " "
    )
  )

  ini({
    # Every value is an author-chosen simulation constant (Table 1; deposited
    # $SIMULATION control stream in Code Examples S1), never an estimate, so
    # all are fixed(). The fitted arms of the paper (OM/PDM/MDM/MDM2/CIM/IDM)
    # publish only relative estimation error and standardized RMSE summaries
    # -- no point estimates -- so they are not encodable as models.
    lka <- fixed(log(2)); label("Absorption rate constant (1/h)") # Table 1 ka = 2 (1/h), shared by all five sets; control stream $THETA 2 ;[KA, 1/h]
    lcl <- fixed(log(15)); label("Clearance (L/h)") # Table 1 set A CL = 15 L/h; control stream $THETA 15 ;[CL, L/h]
    lvc <- fixed(log(100)); label("Central volume of distribution (L)") # Table 1 V = 100 L, shared by all five sets; control stream $THETA 100 ;[V, L]

    # Table 1 gives the BSV variances directly as 0.04 (annotated '20% CV'),
    # and the control stream confirms them as $OMEGA diagonal elements on an
    # exponential eta: CL=TVCL*EXP(ETA(1)), V=TVV*EXP(ETA(2)). No eta on ka.
    # NOTE: Table 1's two omega ROW LABELS are transposed relative to their
    # descriptions ('omega2-V' is described as 'Variance of BSV on CL' and
    # vice versa). Both values are 0.04, so the transposition has no numerical
    # consequence; the control stream settles the mapping.
    etalcl ~ fixed(0.04) # variance = 0.2^2; Table 1 'Variance of BSV on CL' = 0.04 (20% CV); control stream $OMEGA 0.04 omega(1,1)
    etalvc ~ fixed(0.04) # variance = 0.2^2; Table 1 'Variance of BSV on V' = 0.04 (20% CV); control stream $OMEGA 0.04 omega(2,2)

    propSd <- fixed(0.1)
    label("Proportional residual error (fraction)") # sqrt of Table 1 sigma1^2 = 0.01 (annotated '10% CV'); control stream $SIGMA 0.01 ;[P] sigma(1,1)
    addSd <- fixed(1)
    label("Additive residual error (ng/mL)") # sqrt of Table 1 sigma2^2 = 1 (ng/mL); control stream $SIGMA 1 ;[A] sigma(2,2)
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    kel <- cl / vc

    # NONMEM ADVAN2 TRANS2: depot is compartment 1, central is compartment 2.
    # F1 is never set in the control stream, so fraction absorbed is 1
    # ("Complete absorption was assumed and fraction absorbed F was fixed to
    # 1", PK data simulation; Table 1 row F = '1 fixed'). No f(depot) term.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Amounts are mg and vc is L, so central/vc is mg/L; the factor 1000
    # converts to the ng/mL reporting scale of the paper. This reproduces the
    # control stream's scaling statement S2 = V/1000 exactly, since NONMEM
    # forms the prediction as F = A(2)/S2 = 1000 * A(2)/V.
    Cc <- central / vc * 1000

    # NONMEM Y = F * (1 + ERR(1)) + ERR(2) with a DIAGONAL $SIGMA: the two
    # error terms are independent, so Var(Y|F) = F^2 * sigma1^2 + sigma2^2 and
    # sd(Y|F) = sqrt((propSd * F)^2 + addSd^2). That root-sum-square form is
    # nlmixr2's DEFAULT combined2, so no combined1() modifier is used here.
    Cc ~ prop(propSd) + add(addSd)
  })
}
