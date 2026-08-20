Alshehri_2023_ivermectin <- function() {
  description <- paste0(
    "Two-compartment population PK model for ivermectin (IVM) in 56 adults ",
    "(32 Wuchereria bancrofti-infected, 24 uninfected) in Cote d'Ivoire given a ",
    "single 200 ug/kg oral dose as part of ivermectin + diethylcarbamazine + ",
    "albendazole (IDA) triple-drug therapy for lymphatic filariasis, taken after ",
    "a high-fat breakfast (724 plasma samples). Absorption is a zero-order input ",
    "into the depot compartment (TK0 = 3.74 h) preceded by a lag time (Tlag = ",
    "0.757 h), followed by first-order transfer into the central compartment (Ka ",
    "= 0.718 1/h) and linear elimination. Body weight enters as allometric ",
    "scaling with exponents fixed at 1 on Vc/F and Vp/F and 0.75 on CL/F and ",
    "Q/F, centered on the study population's body weight. Sex is the only ",
    "covariate retained by the stepwise selection: Vp/F is 53% lower in men than ",
    "in women (424.3 L in women, 201.8 L in men). Lymphatic filariasis infection ",
    "status was screened and had no effect on any PK parameter. Residual ",
    "variability is combined additive (0.461 ng/mL) plus proportional (22.8%)."
  )
  reference <- paste(
    "Alshehri A, Chhonker YS, Bala V, Edi C, Bjerum CM, Koudou BG, John LN,",
    "Mitja O, Marks M, King CL, Murry DJ (2023). Population pharmacokinetic",
    "model of ivermectin in mass drug administration against lymphatic",
    "filariasis. PLoS Negl Trop Dis 17(6):e0011319.",
    "doi:10.1371/journal.pntd.0011319",
    sep = " "
  )
  vignette <- "Alshehri_2023_ivermectin"
  units <- list(
    time          = "h",
    dosing        = "ug",
    concentration = "ng/mL"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    depot       = list(analyte = "ivermectin", units = "ug", specimen = "administration site", verified = FALSE),
    central     = list(analyte = "ivermectin", units = "ug", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "ivermectin", units = "ug", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Time-fixed (single-dose study). Allometric scaling is applied to all ",
        "four disposition parameters with exponents FIXED by the authors: 1 on ",
        "Vc/F and Vp/F, 0.75 on CL/F and Q/F (Alshehri 2023 Methods 'Covariate ",
        "analysis': 'Fixed exponent values of 1 and 0.75 were applied to the ",
        "apparent volume of distribution and clearance of central and ",
        "peripheral compartments, respectively'; S1 Text PML dVdWeight = 1, ",
        "dV2dWeight = 1, dCldWeight = 0.75, dQdWeight = 0.75). The S1 Text PML ",
        "centers on `mean(Weight)`, the arithmetic mean body weight of the ",
        "analysis dataset, which is NOT reported anywhere in the paper or its ",
        "supplements. This model file centers on 61.6 kg, the median body ",
        "weight reported in Alshehri 2023 Table 1 (range 51-135 kg); see the ",
        "vignette 'Assumptions and deviations' section."
      ),
      source_name        = "Weight"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste0(
        "REFERENCE-CATEGORY INVERSION relative to the source. Alshehri 2023 ",
        "codes its `Sex` covariate as 1 = male / 0 = female (Methods ",
        "'Covariate analysis': 'categorical covariate value of 0 (female sex ",
        "or -ve infection status) or 1 (male sex or +ve infection status)'), ",
        "i.e. the source column is an SEXM-style indicator. The canonical ",
        "column is SEXF (1 = female), so the model file re-expresses the ",
        "effect against the male indicator SEXM = 1 - SEXF and keeps the ",
        "paper's coefficient sign unchanged: Vp/F = tvVp * exp(e_sexf_vp * (1 ",
        "- SEXF)) with e_sexf_vp = -0.74325 (S1 Text PML dV2dSex1). The ",
        "structural typical value `lvp` is therefore the WOMEN's value. ",
        "Arithmetic check against the paper's own reported endpoints: women ",
        "424.33 L (exponent term = 0) and men 424.334 * exp(-0.743251) = ",
        "201.8 L, vs the 424.3 L / 200.4 L pair quoted in Alshehri 2023 ",
        "Results; the 0.7% difference on the male value comes from the paper ",
        "quoting the rounded -0.74 coefficient. Sex is the only covariate ",
        "retained in the final model, and it acts only on Vp/F."
      ),
      source_name        = "Sex (1 = male, 0 = female)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = paste0(
        "Screened in the ETA plots and in the stepwise forward-addition / ",
        "backward-elimination covariate selection, but not retained. Alshehri ",
        "2023 Results: 'Other tested covariates did not show any trend in the ",
        "ETA value distribution of PK parameters' and 'sex was the only ",
        "covariate from the tested covariates to decrease the OFV'. Median 40 ",
        "years (range 18-66), Table 1."
      ),
      source_name = "Age"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "mg/dL",
      type        = "continuous",
      notes       = paste0(
        "Screened but not retained (see AGE note). Median 1.1 mg/dL (range ",
        "0.6-1.6), Table 1. Creatinine clearance was deliberately NOT tested: ",
        "Alshehri 2023 Methods 'Covariate analysis' states 'Creatinine ",
        "clearance (CrCl) was not included in the covariate analysis since IVM ",
        "is mainly eliminated into feces with negligible renal excretion'."
      ),
      source_name = "Creatine level"
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = paste0(
        "Screened but not retained (see AGE note). Median 25 U/L (range ",
        "14-67), Table 1."
      ),
      source_name = "ALT"
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = paste0(
        "Screened but not retained (see AGE note). Median 29 U/L (range ",
        "15-53), Table 1."
      ),
      source_name = "AST"
    ),
    WBANCROFTI_POS = list(
      description = paste0(
        "1 = Wuchereria bancrofti-infected (lymphatic filariasis positive, ",
        "microfilaria plasma level >= 50 Mf/mL at screening), 0 = uninfected ",
        "community control."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (uninfected)",
      notes       = paste0(
        "Screened but not retained -- this is the paper's headline negative ",
        "finding. Alshehri 2023 Abstract: 'The final model identifies that the ",
        "PK parameters of IVM are not affected by LF infection'; S2 Fig shows ",
        "the ETA box plots by infection status. 32 of 56 subjects (57%) were ",
        "infected (Table 1). DOCUMENTATION-ONLY NAME: because the covariate is ",
        "not referenced in model(), `WBANCROFTI_POS` is not ratified in ",
        "inst/references/covariate-columns.md; it is spelled to match the ",
        "existing <PATHOGEN>_POS canonical family (HIV_POS, TB_POS, HCV_POS) ",
        "so that a future model retaining a filarial-infection effect can ",
        "adopt it directly."
      ),
      source_name = "Infection status (1 = +ve, 0 = -ve)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 56L,
    n_studies      = 1L,
    age_range      = "18-66 years (eligibility 18-70 years)",
    age_median     = "40 years",
    weight_range   = "51-135 kg",
    weight_median  = "61.6 kg",
    sex_female_pct = 43,
    race_ethnicity = "not reported; participants resident in the Agboville district, Cote d'Ivoire",
    disease_state  = paste0(
      "32 (57%) treatment-naive Wuchereria bancrofti-infected adults ",
      "(microfilaria plasma level >= 50 Mf/mL) and 24 (43%) uninfected adults"
    ),
    dose_range     = paste0(
      "single oral dose of ivermectin 200 ug/kg co-administered with ",
      "diethylcarbamazine 6 mg/kg and albendazole 400 mg (IDA triple-drug ",
      "therapy), taken after a high-fat breakfast"
    ),
    regions        = "Agboville district, Cote d'Ivoire",
    n_observations = paste0(
      "724 plasma ivermectin concentrations (pre-dose and 1, 2, 3, 4, 6, 8, ",
      "12, 24, 36, 48, 72 h and 7 days post-dose); all above the 0.1 ng/mL ",
      "LLOQ; the 168 h sample was missing for 4 participants"
    ),
    exclusions     = paste0(
      "renal or hepatic disease; ALT, AST or creatinine > 1.5x ULN; ",
      "hemoglobin < 7 g/dL; positive pregnancy test; interacting concomitant ",
      "medication within one week; urinary tract infection; albendazole or ",
      "ivermectin within the past two years"
    ),
    external_validation = paste0(
      "286 additional ivermectin concentrations from 25 participants (12 male ",
      "/ 13 female, median age 27.5 years [18-59], median weight 62 kg ",
      "[46-93]) receiving the same 200 ug/kg dose in an azithromycin-IDA ",
      "drug-interaction trial (S2 Table); used for a VPC-based external ",
      "validation with structural and error parameters fixed"
    ),
    notes          = paste0(
      "Data are from the open-label IDA triple-drug therapy cohort study ",
      "registered as NCT02845713 and NCT03664063. Model built in Phoenix NLME ",
      "8.3 (Certara) using FOCE-I; the full PML control code for the final ",
      "model is reproduced in S1 Text. Baseline demographics are Alshehri 2023 ",
      "Table 1; final parameter estimates with bootstrap CIs are Table 2."
    )
  )

  ini({
    # ========================================================================
    # Structural PK. Alshehri 2023 Table 2 ('Final model estimations of
    # population PK parameters of IVM') and S1 Text (Phoenix PML of the final
    # model). Values below use the full precision reported in the S1 Text
    # `fixef()` statements; Table 2 quotes the same values rounded.
    #
    # UNIT CONVENTION: the S1 Text PML carries volumes in mL and clearances in
    # mL/h with the dose in ug, so that C = A1/V lands in ug/mL. This model
    # file instead carries volumes in L and clearances in L/h with the dose in
    # ug, so that Cc = central/vc lands in ug/L == ng/mL -- the unit the paper
    # reports throughout (LLOQ 0.1 ng/mL; Cmax ~59 ng/mL). Only the
    # volume/clearance values are rescaled by 1e-3; no structural change.
    # ========================================================================
    lka <- log(0.717722387163157)
    label("First-order absorption rate constant Ka (1/h)")
    # S1 Text PML: fixef(tvKa = 0.717722387163157).
    # Table 2: Ka = 0.71 1/h (CV 13.6%); bootstrap median 0.72 (0.53-0.99).

    ld1 <- log(3.73958407653678)
    label("Zero-order duration of the dose input into the depot compartment TK0 (h)")
    # S1 Text PML: fixef(tvTK0 = 3.73958407653678); dosepoint(Aa, duration = TK0).
    # Table 2: Tk0 = 3.73 h (CV 5.3%); bootstrap median 3.73 (3.33-4.11).

    ltlag <- log(0.756753148154855)
    label("Absorption lag time Tlag (h)")
    # S1 Text PML: fixef(tvTlag = 0.756753148154855); dosepoint(Aa, tlag = Tlag).
    # Table 2: Tlag = 0.75 h (CV 9.5%); bootstrap median 0.76 (0.63-0.93).

    lcl <- log(7.02663375480997)
    label("Apparent clearance CL/F (L/h) at the reference body weight")
    # S1 Text PML: fixef(tvCl = 7026.63375480997) mL/h = 7.0266 L/h.
    # Table 2: CL/F = 7.02 L/h (CV 6.7%); bootstrap median 7.01 (6.17-7.95).

    lvc <- log(138.07174209207)
    label("Apparent central volume of distribution Vc/F (L) at the reference body weight")
    # S1 Text PML: fixef(tvV = 138071.74209207) mL = 138.07 L.
    # Table 2: Vc/F = 138 L (CV 8.4%); bootstrap median 137.17 (117-159.50).

    lq <- log(9.11317023736889)
    label("Apparent intercompartmental clearance Q/F (L/h) at the reference body weight")
    # S1 Text PML: fixef(tvQ = 9113.17023736889) mL/h = 9.1132 L/h.
    # Table 2: Q/F = 9.11 L/h (CV 10.4%); bootstrap median 9.04 (7.43-10.88).

    lvp <- log(424.334354384267)
    label("Apparent peripheral volume of distribution Vp/F (L) in WOMEN at the reference body weight")
    # S1 Text PML: fixef(tvV2 = 424334.354384267) mL = 424.33 L.
    # Table 2: Vp/F = 424.33 L (CV 8.4%); bootstrap median 437.94 (368.45-530.87).
    # This is the FEMALE typical value: in the PML the sex term is
    # exp(dV2dSex1 * (Sex == 1)), which is exp(0) = 1 for women (Sex == 0).
    # Alshehri 2023 Results: "The model estimation of Vp/F was 424.3 liters in
    # women and 200.4 liters in men."

    # ========================================================================
    # Allometric scaling on body weight. All four exponents are FIXED by the
    # authors, not estimated. Alshehri 2023 Methods 'Covariate analysis':
    # "Fixed exponent values of 1 and 0.75 were applied to the apparent volume
    # of distribution and clearance of central and peripheral compartments,
    # respectively." Confirmed by S1 Text PML fixef() statements, none of which
    # carry bounds or standard errors for these four terms; Table 2 reports no
    # CV% for any of them.
    # ========================================================================
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent of body weight on CL/F (unitless)")
    # S1 Text PML: fixef(dCldWeight = 0.75).

    e_wt_q <- fixed(0.75)
    label("Allometric exponent of body weight on Q/F (unitless)")
    # S1 Text PML: fixef(dQdWeight = 0.75).

    e_wt_vc <- fixed(1)
    label("Allometric exponent of body weight on Vc/F (unitless)")
    # S1 Text PML: fixef(dVdWeight = 1).

    e_wt_vp <- fixed(1)
    label("Allometric exponent of body weight on Vp/F (unitless)")
    # S1 Text PML: fixef(dV2dWeight = 1).

    # ========================================================================
    # Sex effect on Vp/F -- the only covariate retained by the stepwise
    # forward-addition / backward-elimination procedure.
    # ========================================================================
    e_sexf_vp <- -0.74325100598213
    label("Log-ratio of Vp/F for men vs women (unitless); multiplies the male indicator (1 - SEXF)")
    # S1 Text PML: fixef(dV2dSex1 = -0.74325100598213), entering as
    # exp(dV2dSex1 * (Sex == 1)) where the paper's Sex == 1 denotes MALE.
    # Table 2: "Sex (male) on Vp/F" = -0.74 (CV 12.4%); bootstrap median -0.77.
    # Re-expressed onto the canonical SEXF orientation (1 = female) as
    # exp(e_sexf_vp * (1 - SEXF)); the coefficient sign is unchanged because
    # (1 - SEXF) is exactly the paper's male indicator. Fold change in men
    # = exp(-0.743251) = 0.4756, i.e. Vp/F is 52.4% lower in men, matching the
    # Abstract's "53% lower in men than in women".

    # ========================================================================
    # Inter-individual variability. Alshehri 2023 Table 2 'Random-effect
    # parameters' block reports omega^2 directly (the CV% column is derived as
    # sqrt(exp(omega^2) - 1) * 100, per the Table 2 footnote), so the values
    # below are used as variances without a CV -> variance conversion. The same
    # seven variances appear in the S1 Text PML ranef() statement.
    #
    # DEVIATION: the published final model used a FULL BLOCK (non-diagonal)
    # omega -- Alshehri 2023 Results: "A model with non-diagonal covariance
    # improved the model fitting by minimizing OFV (dOFV = -144)". The
    # off-diagonal covariance estimates are not reported in Table 2, in S1
    # Table, or in the S1 Text PML (whose ranef(block(...)) line carries only
    # the seven diagonal entries). Per the standing "never invent variances"
    # policy the off-diagonals are left at zero here, giving a diagonal omega.
    # Marginal (per-parameter) IIV magnitudes are therefore exactly as
    # published; only the between-parameter correlation is absent. See the
    # vignette 'Assumptions and deviations' section.
    # ========================================================================
    etalvc ~ 0.23
    # Table 2: Vc/F omega^2 = 0.23 (50.8% CV), shrinkage 8.1%; bootstrap 0.22 (49.6%).
    # S1 Text PML ranef(block(nV, ...)) first entry = 0.23.

    etalcl ~ 0.25
    # Table 2: CL/F omega^2 = 0.25 (53.2% CV), shrinkage 2.9%; bootstrap 0.25 (53.2%).

    etalka ~ 0.55
    # Table 2: Ka omega^2 = 0.55 (85.6% CV), shrinkage 19.6%; bootstrap 0.51 (81.5%).

    etalvp ~ 0.17
    # Table 2: Vp/F omega^2 = 0.17 (43% CV), shrinkage 13.9%; bootstrap 0.18 (44.4%).
    # Alshehri 2023 Results: including sex on Vp/F "decreased IIV in Vp/F by 43%
    # from 0.30 to 0.17 compared to the base model".

    etalq ~ 0.19
    # Table 2: Q/F omega^2 = 0.19 (45.7% CV), shrinkage 11.1%; bootstrap 0.19 (45.7%).

    etaltlag ~ 0.35
    # Table 2: Tlag omega^2 = 0.35 (64.7% CV), shrinkage 7.7%; bootstrap 0.35 (64.7%).

    etald1 ~ 0.11
    # Table 2: Tk0 omega^2 = 0.11 (34% CV), shrinkage 23.7%; bootstrap 0.1 (32.4%).

    # ========================================================================
    # Residual variability -- combined additive plus proportional. Alshehri
    # 2023 Results: "The combined additive-proportional error model best
    # explained the residual variability." The S1 Text PML writes this in
    # Phoenix's scaled form
    #   observe(CObs = C + CEps * sqrt(1 + C^2 * (CMultStdev/sigma())^2))
    # with error(CEps = 0.46110892757607), i.e. sigma() = 0.46110892757607.
    # The residual variance is therefore
    #   Var = CEps^2 + C^2 * CMultStdev^2,
    # which is exactly nlmixr2's `add(addSd) + prop(propSd)` with
    # addSd = CEps and propSd = CMultStdev. Epsilon shrinkage was 18%.
    # ========================================================================
    addSd <- 0.46110892757607
    label("Additive residual error (ng/mL)")
    # S1 Text PML: error(CEps = 0.46110892757607).
    # Table 2: additive residual error = 0.46 ng/mL (CV 8.9%); bootstrap 0.45 (0.37-0.64).

    propSd <- 0.228212446759966
    label("Proportional residual error (fraction)")
    # S1 Text PML: fixef(tvCMultStdev = 0.228212446759966).
    # Table 2: proportional residual error = 0.22 (CV 5.4%); bootstrap 0.22 (0.19-0.24).
  })

  model({
    # ----------------------------------------------------------------------
    # 1. Reference body weight for the allometric terms.
    #
    #    The S1 Text PML centers every allometric term on `mean(Weight)`, the
    #    arithmetic mean body weight of the analysis dataset. That mean is not
    #    reported in the paper or in any supplement; Alshehri 2023 Table 1
    #    reports only the MEDIAN, 61.6 kg (range 51-135 kg), which is used
    #    here. At WT = 61.6 kg every allometric factor is exactly 1, so the
    #    ini() values above reproduce the Table 2 typical values verbatim.
    # ----------------------------------------------------------------------
    wtref <- 61.6  # kg; Alshehri 2023 Table 1 median body weight

    # ----------------------------------------------------------------------
    # 2. Individual parameters. IIV is exponential (log-normal) on every
    #    structural parameter, per Alshehri 2023 Methods 'Base model':
    #    "Interindividual variability (IIV) in estimated PK parameters was
    #    assumed to be log-normally distributed".
    # ----------------------------------------------------------------------
    ka   <- exp(lka + etalka)
    d1   <- exp(ld1 + etald1)
    tlag <- exp(ltlag + etaltlag)

    cl <- exp(lcl + etalcl) * (WT / wtref)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / wtref)^e_wt_vc
    q  <- exp(lq  + etalq)  * (WT / wtref)^e_wt_q

    # Sex enters only Vp/F. (1 - SEXF) is the paper's male indicator, so the
    # PML coefficient dV2dSex1 carries over unchanged; lvp is the women's
    # typical value.
    vp <- exp(lvp + etalvp) * (WT / wtref)^e_wt_vp * exp(e_sexf_vp * (1 - SEXF))

    # ----------------------------------------------------------------------
    # 3. Micro-constants. The S1 Text PML uses
    #    cfMicro(A1, Cl/V, Cl2/V, Q/V2, first = (Aa = Ka)), i.e. Phoenix's
    #    cfMicro(A1, Ke, K12, K21, first = ...) macro. `Cl2` is never declared
    #    as an stparm() in the PML and `Q` is; the third argument is therefore
    #    the standard K12 = Q/V, matching K21 = Q/V2 in the fourth argument.
    # ----------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ----------------------------------------------------------------------
    # 4. ODE system. Two-compartment disposition with a separate absorption
    #    (depot) compartment. Alshehri 2023 Fig 2 and Results: "A two-
    #    compartment model with zero-order dose input into the absorption
    #    compartment with a lag time function followed by first-order
    #    absorption and linear elimination best described the IVM disposition."
    #    S1 Text PML: dosepoint(Aa, tlag = Tlag, duration = TK0) with
    #    first = (Aa = Ka); `Aa` is the depot compartment here.
    # ----------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central +
                          k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # ----------------------------------------------------------------------
    # 5. Dosing modifiers on the depot compartment: the dose is delivered as a
    #    zero-order input over d1 hours, beginning tlag hours after the dose
    #    record. Requires rate = -2 (modeled duration) on the dose event.
    # ----------------------------------------------------------------------
    dur(depot)  <- d1
    alag(depot) <- tlag

    # ----------------------------------------------------------------------
    # 6. Observation. Dose in ug and vc in L give ug/L == ng/mL.
    # ----------------------------------------------------------------------
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
