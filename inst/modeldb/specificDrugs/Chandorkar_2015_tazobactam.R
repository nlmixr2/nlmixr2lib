Chandorkar_2015_tazobactam <- function() {
  description <- paste(
    "Two-compartment population PK model for tazobactam given as a 1-hour",
    "intravenous infusion in a fixed 2:1 ceftolozane:tazobactam combination,",
    "fitted to 4,249 plasma concentrations from 243 adults pooled across the",
    "subset of studies in which tazobactam was co-administered: Phase 1",
    "studies in healthy volunteers, studies spanning mild to severe renal",
    "impairment, and one Phase 2 study in patients with complicated",
    "intra-abdominal infection (cIAI) (Chandorkar 2015). Elimination is first",
    "order. Baseline Cockcroft-Gault creatinine clearance acts on clearance",
    "through a power function, and cIAI raises the central volume.",
    "Between-subject variability is diagonal on CL and Vc only; the",
    "peripheral parameters carry none, and the residual error is purely",
    "proportional. Companion model to Chandorkar_2015_ceftolozane, which the",
    "same paper fitted separately because the two drugs do not interact."
  )
  reference <- paste(
    "Chandorkar G, Xiao A, Mouksassi MS, Hershberger E, Krishna G.",
    "Population pharmacokinetics of ceftolozane/tazobactam in healthy",
    "volunteers, subjects with varying degrees of renal function and patients",
    "with bacterial infections. J Clin Pharmacol. 2015;55(2):230-239.",
    "doi:10.1002/jcph.395.",
    "All fixed-effect, random-effect and residual-error estimates are Table 3",
    "panel [B]. No supplement was deposited with the article (EuropePMC",
    "hasSuppl 'N'); Table 3B, the Results narrative and Figure 1B together",
    "report the complete final model.",
    sep = " "
  )
  vignette <- "Chandorkar_2015_ceftolozane_tazobactam"
  units    <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Doses are in mg and the paper reports plasma
  # concentrations in ug/mL, numerically identical to mg/L.
  compartmentData <- list(
    central     = list(analyte = "tazobactam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "tazobactam", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = paste(
        "Baseline creatinine clearance estimated by the Cockcroft-Gault",
        "formula; raw mL/min, NOT body-surface-area normalized"
      ),
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Chandorkar 2015 Table 3B: 'CL (L/h), 18.0 (3.39)*(CrCL/115)^0.67",
        "(11.1)'. Reference value 115 mL/min; power exponent 0.67 (RSE",
        "11.1%). The Results restate it -- 'the effect of baseline CrCL on CL",
        "showing a power function of 0.67 (ie, [CrCL/115]^0.67)' -- as does",
        "the Discussion: 'CL of tazobactam = 18.0 L/h * (CrCL/115)^0.67'.",
        "NOTE THE DIFFERENT REFERENCE from the companion ceftolozane model,",
        "which centres on 109 mL/min. The two analytes were fitted to",
        "different subject sets (376 vs 243) with different median renal",
        "function, so 115 and 109 are not interchangeable; using the",
        "ceftolozane centring here would misstate tazobactam CL by about 4%.",
        "ASSAY FORM. Methods: Cockcroft-Gault on ACTUAL body weight with the",
        "0.85 female multiplier, in raw mL/min, NOT normalized to 1.73 m^2.",
        "BASELINE (time-fixed) per subject; serum creatinine was stable over",
        "the short treatment duration (Results, 'Datasets').",
        "RANGE FITTED. Table 2 gives tazobactam-cohort means 100.4 mL/min",
        "(range 19-238) without infection and 105 (range 41-309) with cIAI;",
        "Figure 1B labels the observed range as 19.1-308.5 mL/min. The",
        "renal-impairment categories quoted in the paper are descriptive",
        "strata only -- the model uses the continuous column."
      ),
      source_name        = "CrCL"
    ),
    DIS_CIAI = list(
      description        = "Complicated intra-abdominal infection cohort indicator (1 = cIAI)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy volunteer or renal-impairment subject with no infection)",
      notes              = paste(
        "Chandorkar 2015 Table 3B row 'Vc (L), With cIAI, x1.47 (21.9)'.",
        "Printed as an exp(beta) factor under the paper's stated",
        "parameterisation for categorical covariates, so the ini()",
        "coefficient is log(1.47). The Discussion restates it: 'Vc of",
        "tazobactam = 14.2 L with a multiplicative factor of 1.47 for",
        "patients with cIAI ... the Vc of tazobactam is about 47% larger in",
        "cIAI patients than in healthy subjects'.",
        "VOLUME ONLY. Unlike the companion ceftolozane model, tazobactam",
        "carries NO infection effect on clearance -- Table 3B's CL row has a",
        "single unstratified entry, and the Results name only 'the effect of",
        "cIAI infection on Vc and of CrCL on CL' as retained.",
        "NO cUTI TERM EXISTS FOR THIS ANALYTE, and this is a data fact rather",
        "than an omission: Results, 'note there were no tazobactam data from",
        "cUTI patients'. Table 1 confirms the cUTI study (Umeh 2010,",
        "NCT00921024) dosed ceftolozane 1000 mg q8h alone. Do not borrow the",
        "companion model's cUTI coefficients.",
        "77 of the 243 tazobactam subjects had cIAI (Table 2), receiving",
        "ceftolozane/tazobactam 1000/500 mg q8h."
      ),
      source_name        = "cIAI"
    )
  )

  # Covariates the paper screened but did not retain in the final tazobactam
  # model. Documented so the provenance of the covariate screen survives,
  # without declaring covariateData entries that model() never references.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Actual total body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "SCREENED, REACHED THE STEPWISE PROCEDURE, BUT NOT RETAINED IN THE",
        "FINAL MODEL. The Results report the forward-addition trace 'The",
        "dMOF2 was -103.02 (P = .001) when the effect of cIAI infection was",
        "included in the model and -109.73 (P = .01) when the effect of",
        "weight on Vc was also included', so a weight-on-Vc term was tested",
        "and improved the objective function. It nonetheless does not appear",
        "in the final model: Table 3B's Vc rows are a bare '14.2' and",
        "'x1.47' with no '*(weight/74)' term (contrast Table 3A, where the",
        "ceftolozane Vc rows do carry it), and the Discussion's final",
        "equation is 'Vc of tazobactam = 14.2 L with a multiplicative factor",
        "of 1.47 for patients with cIAI' -- no weight. The Results also list",
        "only 'the effect of cIAI infection on Vc (note there were no",
        "tazobactam data from cUTI patients) and of CrCL on CL' as the",
        "retained effects, and state 'No trends were noted between other",
        "covariates tested and tazobactam PK.' The backward-elimination step",
        "used a stricter threshold (P < .001) than the forward step",
        "(P < .01), which is consistent with a term admitted at P = .01",
        "forward being dropped again. NO COEFFICIENT IS PUBLISHED for it in",
        "any case, so it could not be carried even if it had been retained.",
        "Weight still enters this model INDIRECTLY through the",
        "Cockcroft-Gault CRCL column, which is computed on actual body",
        "weight. See the vignette Errata."
      )
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened as an intrinsic covariate (Methods) and rejected. Results:",
        "'No trends were noted between other covariates tested and tazobactam",
        "PK.' Cohort range 18-86 years."
      )
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened and rejected. Sex nonetheless enters the model INDIRECTLY",
        "through the Cockcroft-Gault CRCL column, which multiplies by 0.85",
        "for female subjects."
      )
    ),
    RACE_WHITE = list(
      description = "White race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened and rejected. The tazobactam cohort was 81.9% white without",
        "infection and 98.7% with cIAI (Table 2), so the analysis had little",
        "power to resolve a race effect."
      )
    ),
    DOSE_TAZOBACTAM_MG = list(
      description = "Administered tazobactam dose",
      units       = "mg",
      type        = "continuous",
      notes       = paste(
        "Screened as an extrinsic covariate to test for dose-dependent",
        "(nonlinear) PK and rejected, supporting the linear structural model.",
        "Doses spanned 250-1500 mg tazobactam (Table 1). Discussion: 'the PK",
        "of ceftolozane/tazobactam is dose-proportional and linear'."
      )
    ),
    CONMED_CEFTOLOZANE = list(
      description = "Co-administration of ceftolozane indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened as the drug-drug-interaction covariate and rejected.",
        "Discussion: 'no drug-drug interaction was observed between",
        "ceftolozane and tazobactam'. This is the finding that licenses",
        "fitting -- and extracting -- ceftolozane and tazobactam as two",
        "independent models. In practice every tazobactam subject in this",
        "analysis also received ceftolozane, since the drug is only given in",
        "the fixed 2:1 combination."
      )
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 243L,
    n_studies        = 10L,
    n_concentrations = 4249L,
    age_range        = "18-86 years",
    age_median       = "means 43.7 years (no infection) and 47.0 years (cIAI); medians not reported",
    weight_range     = "50-145 kg",
    weight_median    = "means 74.1 kg (no infection) and 78.0 kg (cIAI); medians not reported",
    sex_female_pct   = 42.8,
    race_ethnicity   = c(White = 87.2, Other = 12.8),
    disease_state    = paste(
      "Pooled healthy adults, adults with mild to severe renal impairment,",
      "and hospitalized patients with complicated intra-abdominal infection.",
      "Table 2: 166 of 243 subjects had no infection and 77 had cIAI. NO cUTI",
      "PATIENTS contributed tazobactam data, because the Phase 2 cUTI study",
      "dosed ceftolozane alone."
    ),
    dose_range       = paste(
      "All doses given as 1-hour intravenous infusions of the fixed 2:1",
      "ceftolozane:tazobactam combination. Healthy volunteers received single",
      "or multiple (q8h or q12h) doses of ceftolozane/tazobactam 500/250,",
      "1000/500, 1500/750, 2000/1000 or 3000/1500 mg. Renal-impairment",
      "subjects received a single 1000/500 mg (mild/moderate) or 500/250 mg",
      "(severe) dose. cIAI patients received 1000/500 mg q8h (Table 1). The",
      "tazobactam component is therefore 250-1500 mg."
    ),
    renal_function   = paste(
      "Spans severe renal impairment to augmented clearance. Estimated CrCL",
      "means 100.4 mL/min (range 19-238) without infection and 105 mL/min",
      "(range 41-309) with cIAI; Figure 1B gives the observed range as",
      "19.1-308.5 mL/min. Category counts (Table 2, no infection / cIAI):",
      "normal 137/48, mild 17/26, moderate 6/3, severe 6/0. No",
      "end-stage-renal-disease or dialysis subjects were enrolled."
    ),
    bmi_range        = "18-51 kg/m^2 (means 26.0 without infection, 26.6 with cIAI)",
    pkpd_target      = paste(
      "Not fitted in this paper. Discussion: 'the PD driver for tazobactam is",
      "thought to be the percentage of time above a threshold concentration",
      "(%T>threshold)'. The paper positions this PK model as the input to a",
      "later probability-of-target-attainment analysis rather than performing",
      "one."
    ),
    notes            = paste(
      "sex_female_pct is the pooled tazobactam data set: (70 + 34) / 243 =",
      "42.8% (Table 2). race_ethnicity is likewise pooled white (136 + 76) /",
      "243 = 87.2%; the paper reports only the white percentage, so the",
      "remainder is recorded as Other rather than broken out. n_studies is",
      "the 10 studies of the pooled analysis (Table 1); the tazobactam data",
      "come from the subset in which tazobactam was co-administered.",
      "Estimation was first-order conditional estimation-extended least",
      "squares (FOCE-ELS) in Phoenix NLME 1.2 (Certara), NOT NONMEM.",
      "Model evaluation used goodness-of-fit diagnostics, a visual predictive",
      "check with 1,000 replicates, and 1,000-sample nonparametric bootstrap",
      "resampling, with bootstrap-versus-final differences under 4%.",
      "Seventeen samples from 13 subjects had CWRES > 4 and were RETAINED;",
      "the paper notes that excluding them would have reduced the BSV of CL",
      "and Vc by 42.5% and 32.5% respectively, so the IIV magnitudes below",
      "are the outlier-inclusive ones."
    )
  )

  ini({
    # =====================================================================
    # STRUCTURAL PARAMETERS -- Chandorkar 2015 Table 3B, "Population
    # estimates (RSE %)" column. Values refer to the reference subject:
    # CrCL = 115 mL/min, no infection.
    # =====================================================================
    lcl <- log(18.0); label("Tazobactam clearance at CrCL = 115 mL/min (L/h)")                     # Table 3B: CL = 18.0 L/h (RSE 3.39%)
    lvc <- log(14.2); label("Tazobactam central volume of distribution, no infection (L)")         # Table 3B: Vc, No infection = 14.2 L (RSE 4.45%)
    lq  <- log(3.13); label("Tazobactam inter-compartmental clearance (L/h)")                      # Table 3B: CL2 = 3.13 L/h (RSE 4.59%)
    lvp <- log(4.29); label("Tazobactam peripheral volume of distribution (L)")                    # Table 3B: Vp = 4.29 L (RSE 2.61%)
    # Vp is ESTIMATED here (RSE 2.61%), unlike the companion ceftolozane
    # model where Table 3A prints "2.88 (fixed)".

    # =====================================================================
    # COVARIATE EFFECTS. Methods: continuous covariates enter as a power
    # model standardized by the median; categorical covariates as "a linear
    # model with an exponentiated factor relative to the reference", which
    # is what Table 3B prints as "x1.47". The categorical ini() coefficient
    # is therefore the log of the printed factor.
    # =====================================================================
    e_crcl_cl <- 0.67; label("Power exponent of baseline CrCL on tazobactam CL, CRCL/115 (unitless)")
    # Table 3B: 18.0*(CrCL/115)^0.67, RSE 11.1% on the exponent. Cross-checked
    # against the Figure 1B tornado panel, whose six printed relative-CL
    # values (19.1 -> 0.29, 15 -> 0.25, 30 -> 0.40, 50 -> 0.57, 90 -> 0.84,
    # 308.5 -> 1.93) are reproduced by (CrCL/115)^0.67 to within the panel's
    # two-decimal rounding.

    e_ciai_vc <- log(1.47); label("Effect of cIAI on tazobactam Vc (log multiplicative factor)")   # Table 3B: Vc, With cIAI = x1.47 (RSE 21.9%)

    # =====================================================================
    # BETWEEN-SUBJECT VARIABILITY -- Table 3B, "BSV % (RSE %)" column.
    #
    # Results: "a two-compartmental structural model with a diagonal
    # variance (BSV) for CL and Vc". Diagonal, and the CL2 / Vp rows read
    # "Fixed at 0", so only CL and Vc carry IIV.
    #
    # SCALE. Methods: the variance component "assumed a log-normal
    # distribution of PK parameters", so the %CV column converts as
    # omega^2 = log(1 + CV^2). The reported RSEs corroborate an SD-like
    # quantity rather than a variance: for a variance estimated on 243
    # subjects the asymptotic RSE is about sqrt(2/243) = 9.1%, against about
    # 4.5% for the corresponding SD, and the printed RSEs are 4.98% and
    # 6.14%.
    # =====================================================================
    etalcl ~ 0.224745  # BSV on CL 50.2% CV (RSE 4.98%, shrinkage 4.68%): log(1 + 0.502^2)
    etalvc ~ 0.243436  # BSV on Vc 52.5% CV (RSE 6.14%, shrinkage 11.5%): log(1 + 0.525^2)

    # =====================================================================
    # RESIDUAL VARIABILITY -- Table 3B. Purely proportional: Results, "a
    # proportional model for unexplained residual variability", and Table 3B
    # has no additive row, unlike Table 3A.
    # =====================================================================
    propSd <- 0.260; label("Proportional residual error (fraction)")  # Table 3B: Proportional error = 26.0% (RSE 1.64%)
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual PK parameters.
    #
    # Discussion: "The final equations for CL and Vc of tazobactam were:
    # CL of tazobactam = 18.0 L/h * (CrCL/115)^0.67 and Vc of tazobactam =
    # 14.2 L with a multiplicative factor of 1.47 for patients with cIAI."
    # Clearance carries no infection term and volume carries no renal term.
    # ------------------------------------------------------------------
    cl <- exp(lcl + etalcl) * (CRCL / 115)^e_crcl_cl
    vc <- exp(lvc + etalvc + e_ciai_vc * DIS_CIAI)

    q  <- exp(lq)
    vp <- exp(lvp)

    # ------------------------------------------------------------------
    # 2. Micro-constants.
    # ------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ------------------------------------------------------------------
    # 3. Two-compartment disposition. Tazobactam is given as a zero-order
    #    1-hour intravenous infusion straight into the central compartment;
    #    the infusion is expressed in the event table, not here.
    # ------------------------------------------------------------------
    d/dt(central)     <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <-   k12 * central        - k21 * peripheral1

    # ------------------------------------------------------------------
    # 4. Observation. Dose in mg, volumes in L -> mg/L, numerically
    #    identical to the ug/mL the paper reports. Total plasma
    #    concentrations; no unbound fraction is reported in this paper.
    # ------------------------------------------------------------------
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
