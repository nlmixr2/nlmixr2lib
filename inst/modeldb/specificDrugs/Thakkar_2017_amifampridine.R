Thakkar_2017_amifampridine <- function() {
  description <- paste(
    "Joint parent-metabolite population PK + fractional-Emax PD model",
    "for 3,4-diaminopyridine (3,4-DAP, amifampridine) free base and its",
    "N-acetyl metabolite 3-Ac DAP in 49 adults with Lambert-Eaton",
    "myasthenia (Thakkar 2017). Two-compartment parent + one-compartment",
    "metabolite with Fm fixed to 1 (all parent clearance forms",
    "metabolite). Body weight is allometrically scaled on CL/F and",
    "CLm/F3ACDAP (exponent 0.75 fixed) and linearly on Vp/F (exponent 1",
    "fixed), all with reference weight 82 kg. Serum creatinine acts on",
    "CLm/F3ACDAP through (0.8/SCR)^0.7 with median SCR 0.8 mg/dL. The PD",
    "submodel describes the Triple Timed Up and Go (3TUG) score in",
    "seconds via a fractional-inhibitory Emax equation Effect = E0 *",
    "(1 - Emax * Cp / (EC50 + Cp)) where Cp is the parent 3,4-DAP",
    "plasma concentration in ng/mL."
  )
  reference <- paste(
    "Thakkar N, Guptill JT, Ales K, Jacobus D, Jacobus L, Peloquin C,",
    "Cohen-Wolkowiez M, Gonzalez D, for the DAPPER Study Group.",
    "Population Pharmacokinetics/Pharmacodynamics of 3,4-Diaminopyridine",
    "Free Base in Patients With Lambert-Eaton Myasthenia.",
    "CPT Pharmacometrics Syst Pharmacol. 2017;6(9):625-634.",
    "doi:10.1002/psp4.12218. ClinicalTrials.gov NCT01511978.",
    sep = " "
  )
  vignette <- "Thakkar_2017_amifampridine"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed total body weight. Allometric scaling on CL/F and",
        "CLm/F3ACDAP (exponent 0.75 fixed) and linear scaling on Vp/F",
        "(exponent 1 fixed); reference weight 82 kg (population median",
        "from Thakkar 2017 Table 1)."
      ),
      source_name        = "TBW"
    ),
    CREAT = list(
      description        = "Serum creatinine",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed serum creatinine. Power-form covariate on",
        "metabolite clearance: CLm = CLm_pop * (0.8/CREAT)^0.7 (median",
        "SCR 0.8 mg/dL from Thakkar 2017 Table 1; the implementation",
        "applies the ratio so that higher SCR decreases CLm). Source",
        "column name SCR is the canonical CREAT alias per",
        "inst/references/covariate-columns.md."
      ),
      source_name        = "SCR"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 49,
    n_studies      = 1,
    age_range      = "23-83 years",
    age_median     = "60 years",
    weight_range   = "45.8-131.5 kg",
    weight_median  = "82.6 kg",
    sex_female_pct = 53,
    race_ethnicity = c(White = 94, Other = 6),
    disease_state  = "Adults with Lambert-Eaton myasthenia on chronic 3,4-DAP free base treatment",
    dose_range     = "10-30 mg single oral doses (median 20 mg); total daily dose 30-100 mg (median 80 mg)",
    regions        = "United States (DAPPER multicenter trial)",
    notes          = paste(
      "Demographics from Thakkar 2017 Table 1. Single multicenter",
      "double-blind placebo-controlled withdrawal phase II study",
      "(ClinicalTrials.gov NCT01511978). PK analysis n = 49 (1270",
      "samples). PD analysis n = 32 randomized patients (1091 3TUG",
      "data points). 1 patient (2%) reported Hispanic ethnicity.",
      "NAT2 acetylator status was available for only 11/49 subjects",
      "(5 slow, 5 intermediate, 1 rapid) and did not retain",
      "statistical significance in the final model. The paper reports",
      "no significant covariate effects in the PD submodel."
    )
  )

  ini({
    # =====================================================================
    # PK structural parameters - Thakkar 2017 Table 2 Final model column
    # Reference weight 82 kg (population median); reference SCR 0.8 mg/dL.
    # CL/F and CLm/F3ACDAP are apparent clearances (oral); Vc/F, Vp/F and
    # Vm/F3ACDAP are apparent volumes. Fm is fixed to 1 so all parent CL
    # forms metabolite.
    # =====================================================================
    lka      <- log(0.9);   label("3,4-DAP absorption rate constant Ka (1/h)")                   # Table 2 Final: KA = 0.9 /h (RSE 8.9%)
    lcl      <- log(90);    label("3,4-DAP apparent clearance CL/F at 82 kg (L/h)")              # Table 2 Final: CL/F = 90 L/h (RSE 7.4%)
    lvc      <- log(24);    label("3,4-DAP apparent central volume Vc/F (L)")                    # Table 2 Final: Vc/F = 24 L (RSE 33.3%)
    lq       <- log(111);   label("3,4-DAP apparent intercompartmental clearance Q/F (L/h)")     # Table 2 Final: Q/F = 111 L/h (RSE 12.1%)
    lvp      <- log(669);   label("3,4-DAP apparent peripheral volume Vp/F at 82 kg (L)")        # Table 2 Final: Vp/F = 669 L (RSE 17.6%)
    lcl_acdap <- log(20.5); label("3-Ac DAP apparent clearance CLm/F3ACDAP at 82 kg, SCR 0.8 mg/dL (L/h)")  # Table 2 Final: CLm/F3ACDAP = 20.5 L/h (RSE 6%)
    lvc_acdap <- log(36);   label("3-Ac DAP apparent volume Vm/F3ACDAP (L)")                     # Table 2 Final: Vm/F3ACDAP = 36 L (RSE 10.7%)

    # Allometric / SCR exponents
    e_wt_cl       <- fixed(0.75); label("Allometric exponent of WT on CL/F (unitless)")          # Methods Eq 3: fixed at 0.75 per allometric convention
    e_wt_vp       <- fixed(1);    label("Linear exponent of WT on Vp/F (unitless)")              # Methods Eq 4: fixed at 1 (volume scales linearly with WT)
    e_wt_cl_acdap <- fixed(0.75); label("Allometric exponent of WT on CLm/F3ACDAP (unitless)")   # Methods Eq 3: fixed at 0.75 per allometric convention
    e_creat_cl_acdap <- 0.7;      label("Power exponent of SCR ratio on CLm/F3ACDAP (unitless)") # Table 2 Final: SCR exponent = 0.7 (RSE 31.8%); applied as (0.8/SCR)^0.7

    # =====================================================================
    # PK IIV - omega^2 = log(CV^2 + 1) for log-normal eta on log-transformed
    # parameters; CV values from Thakkar 2017 Table 2 "Between subject
    # variability (%CV)" block. CL and Vc are correlated (r = 0.53;
    # Table 2 footnote d) - covariance derived as 0.53 * sqrt(var_CL * var_Vc).
    # =====================================================================
    etalcl + etalvc ~ c(0.2082,
                        0.2926, 1.4649)
    # Table 2 Final: CV(CL/F) = 48.1% -> var = log(1 + 0.481^2) = 0.2082;
    # CV(Vc/F) = 182.4% -> var = log(1 + 1.824^2) = 1.4649;
    # correlation 0.53 (footnote d) -> cov = 0.53 * sqrt(0.2082 * 1.4649) = 0.2926.
    etalka       ~ 0.1239
    # Table 2 Final: CV(Ka) = 36.2% -> var = log(1 + 0.362^2) = 0.1239
    etalq        ~ 0.2122
    # Table 2 Final: CV(Q/F) = 48.7% -> var = log(1 + 0.487^2) = 0.2122
    etalvp       ~ 0.6989
    # Table 2 Final: CV(Vp/F) = 89.8% -> var = log(1 + 0.898^2) = 0.6989
    etalcl_acdap ~ 0.1389
    # Table 2 Final: CV(CLm/F3ACDAP) = 38.6% -> var = log(1 + 0.386^2) = 0.1389

    # =====================================================================
    # PK residual error - proportional, separate magnitudes for parent and
    # metabolite (Thakkar 2017 Table 2 Residual error block).
    # =====================================================================
    propSd       <- 0.348; label("3,4-DAP proportional residual error (fraction)")         # Table 2 Final: 34.8% (RSE 16.7%)
    propSd_acdap <- 0.201; label("3-Ac DAP proportional residual error (fraction)")        # Table 2 Final: 20.1% (RSE 12.7%)

    # =====================================================================
    # PD structural parameters - Thakkar 2017 Table 3 Final model column.
    # Effect (sec) = E0 * (1 - Emax * Cp / (EC50 + Cp)) with EC50 in ng/mL
    # of parent 3,4-DAP plasma concentration; the paper's Results section
    # quotes the explicit form Effect = 18.2 * (1 - 0.816 * Cp / (29.8 + Cp)).
    # Fractional Emax is constrained to (0, 1) via a logit transformation.
    # =====================================================================
    logitemax <- log(0.816 / (1 - 0.816)); label("Logit of fractional Emax (unitless)")     # Table 3 Final: Fractional Emax = 0.816 (RSE 16%); logit(0.816) = 1.4895
    lec50     <- log(29.8);                label("3,4-DAP EC50 for 3TUG response (ng/mL)")  # Table 3 Final: EC50 = 29.8 ng/mL (RSE 36%; 273 nM equivalent)
    le0       <- log(18.2);                label("Baseline 3TUG in absence of drug (sec)")  # Table 3 Final: E0 = 18.2 sec (RSE 11.5%)

    # =====================================================================
    # PD IIV - logit-scale variance for Emax (paper reports the variance
    # directly as 2.93 on logit scale); log-normal etas for EC50 and E0
    # with CV% converted via omega^2 = log(1 + CV^2).
    # =====================================================================
    etalogitemax ~ 2.93
    # Table 3 Final: logit-scale variance for fractional Emax = 2.93 (RSE 71%; bootstrap 0.35 - 12.85)
    etalec50     ~ 0.5760
    # Table 3 Final: CV(EC50) = 88.3% -> var = log(1 + 0.883^2) = 0.5760
    etale0       ~ 0.4109
    # Table 3 Final: CV(E0) = 71.3% -> var = log(1 + 0.713^2) = 0.4109

    # =====================================================================
    # PD residual error - proportional on observed 3TUG (sec).
    # =====================================================================
    propSd_tug3 <- 0.214; label("3TUG proportional residual error (fraction)")              # Table 3 Final: 21.4% (RSE 31%)
  })
  model({
    # ------------------------------------------------------------
    # Individual PK parameters - reference 82 kg, SCR 0.8 mg/dL.
    # Body weight enters allometrically on CL/F and CLm/F3ACDAP
    # (exponent 0.75 fixed) and linearly on Vp/F (exponent 1 fixed).
    # Serum creatinine enters CLm/F3ACDAP through (0.8/CREAT)^0.7.
    # ------------------------------------------------------------
    ka       <- exp(lka + etalka)
    cl       <- exp(lcl + etalcl) * (WT / 82)^e_wt_cl
    vc       <- exp(lvc + etalvc)
    q        <- exp(lq  + etalq)
    vp       <- exp(lvp + etalvp) * (WT / 82)^e_wt_vp
    cl_acdap <- exp(lcl_acdap + etalcl_acdap) *
                (WT / 82)^e_wt_cl_acdap *
                (0.8 / CREAT)^e_creat_cl_acdap
    vc_acdap <- exp(lvc_acdap)

    # ------------------------------------------------------------
    # Micro-constants and ODE system. Fm is fixed to 1 in the source,
    # so the parent CL flux feeds the metabolite compartment in full.
    # Concentrations are computed in ng/mL via the dose-mg / volume-L
    # to mg/L scale and then multiplied by 1000 (1 mg/L = 1000 ng/mL).
    # ------------------------------------------------------------
    kel  <- cl / vc
    k12  <- q  / vc
    k21  <- q  / vp
    kelm <- cl_acdap / vc_acdap

    d/dt(depot)         <- -ka * depot
    d/dt(central)       <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1)   <-  k12 * central - k21 * peripheral1
    d/dt(central_acdap) <-  kel * central - kelm * central_acdap

    Cc       <- 1000 * central       / vc        # ng/mL parent
    Cc_acdap <- 1000 * central_acdap / vc_acdap  # ng/mL metabolite (apparent, relative to F3ACDAP)

    # ------------------------------------------------------------
    # PD - fractional-inhibitory Emax on 3TUG response (seconds).
    # Emax is bounded to (0, 1) via a logit transformation of the
    # population logit plus a normal-distributed individual deviation.
    # ------------------------------------------------------------
    logit_emax_i <- logitemax + etalogitemax
    emax <- exp(logit_emax_i) / (1 + exp(logit_emax_i))
    ec50 <- exp(lec50 + etalec50)
    e0   <- exp(le0   + etale0)
    tug3 <- e0 * (1 - emax * Cc / (ec50 + Cc))

    # ------------------------------------------------------------
    # Observations. Proportional residual error on each output.
    # ------------------------------------------------------------
    Cc       ~ prop(propSd)
    Cc_acdap ~ prop(propSd_acdap)
    tug3     ~ prop(propSd_tug3)
  })
}
