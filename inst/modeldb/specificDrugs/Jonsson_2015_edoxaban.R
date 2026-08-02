Jonsson_2015_edoxaban <- function() {
  description <- "Joint parent plus metabolite population PK model for edoxaban and its main metabolite M4 in adults with normal kidney function through severe renal impairment, from a dedicated renal impairment study. Edoxaban disposition is two-compartment with Savic 2007 analytical transit-compartment absorption (non-integer NN, with the absorption rate constant constrained equal to the transit rate constant ktr) and absolute oral bioavailability estimated on the logit scale; edoxaban clearance is split into a renal arm that is linear in creatinine clearance and drives a urinary-excretion compartment, and a non-renal arm that is assumed to form all of the M4 metabolite. M4 is described by a one-compartment model with first-order formation. Fixed allometric body-weight scaling (0.75 on all clearance terms that are not a function of kidney function, 1 on all volumes) is applied around a 70 kg reference subject."
  reference <- paste(
    "Jonsson S., Simonsson U. S. H., Miller R., Karlsson M. O. (2015).",
    "Population pharmacokinetics of edoxaban and its main metabolite in a",
    "dedicated renal impairment study.",
    "The Journal of Clinical Pharmacology 55(11):1268-1279.",
    "doi:10.1002/jcph.541.",
    sep = " "
  )
  vignette <- "Jonsson_2015_edoxaban"

  # The original NONMEM analysis was run in molar units throughout: doses in
  # nmol and plasma / urine concentrations in nmol/L (nM). Jonsson 2015 Methods
  # "Data": "By using concentrations in nM and dose amounts in nmol, the
  # difference in molecular weights for edoxaban and M4 was taken into
  # account." Working in molar units is what makes the 1:1 metabolite-formation
  # stoichiometry in d/dt(central_m4) below correct without a molecular-weight
  # ratio. To dose in mass units, convert with the edoxaban free-base
  # molecular weight of 548.0 g/mol (15 mg = 27372 nmol).
  units <- list(time = "hour", dosing = "nmol", concentration = "nmol/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight, used for fixed allometric scaling around a 70 kg reference subject",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed body weight. Jonsson 2015 Results 'Final Population PK Model':",
        "'Introduction of allometrically scaled body weight on disposition parameters",
        "(CLNR, central and peripheral volume of distribution for edoxaban [Vc, Vp],",
        "intercompartmental clearance for edoxaban [Q], CLM4, and VM4) with exponents",
        "fixed to 1 and 3/4 for volume and clearance terms, respectively, resulted in an",
        "improved fit'. Applied as (WT/70)^0.75 on cl_nonren, q and cl_m4 and as",
        "(WT/70)^1 on vc, vp and vc_m4. NOTE: the renal clearance arm cl_renal is",
        "deliberately NOT allometrically scaled, because allometry was applied only to",
        "'all clearance terms not being a function of kidney function' (Jonsson 2015",
        "Methods 'Modeling'); cl_renal is instead a pure linear function of CRCL.",
        "Observed body weight in the four renal-function groups was 74.4 (58.9-89.4),",
        "76.5 (60.3-91.0), 78.6 (58.0-90.0) and 71.7 (56.0-95.0) kg (Jonsson 2015",
        "Results 'Data'); study eligibility required 55-110 kg."
      ),
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Creatinine clearance estimated by the Cockcroft-Gault formula, NOT body-surface-area normalised",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Raw (non-BSA-normalised) Cockcroft-Gault creatinine clearance in mL/min, the",
        "single kidney-function covariate of the model and the variable used to allocate",
        "subjects to renal-function groups. Jonsson 2015 Methods 'Data': 'Subjects",
        "classified based on CLcr, estimated by the Cockcroft-Gault formula'. Enters the",
        "model as a linear, non-centred, through-the-origin term on the renal clearance",
        "arm: cl_renal = slopeCLcr * CLcr (Jonsson 2015 Table 1 footnote d: 'typical",
        "CL = CLNR,typ * (WT/70)^(3/4) + CLR, where CLR is the slopeClcr * CLcr').",
        "Because the relationship passes through the origin the paper has no reference",
        "or centring value; a subject with CLcr = 0 has cl_renal = 0 and eliminates",
        "edoxaban by the non-renal route only. This implementation anchors the same",
        "straight line at the paper's typical CLcr = 100 mL/min so the arm can carry the",
        "canonical lcl_renal / etalcl_renal names (an exact reparameterisation, since",
        "slope * CLcr == (slope * 100) * (CLcr/100)); recover the published slope as",
        "exp(lcl_renal)/100 = 0.109. Renal-function groups (Jonsson 2015",
        "Methods 'Data'): normal CLcr > 80, mild 50 <= CLcr <= 80, moderate",
        "30 <= CLcr < 50, severe CLcr < 30 mL/min (not on dialysis). Observed group",
        "means (ranges) were 94.6 (83.0-123.0), 64.7 (54.0-77.0), 42.0 (33.0-49.0) and",
        "21.8 (14.0-27.0) mL/min (Jonsson 2015 Results 'Data'). The units are plain",
        "mL/min rather than the register's default mL/min/1.73 m^2; the same raw",
        "Cockcroft-Gault usage is already registered for CRCL via the CLCR source alias",
        "(precedent: Delattre_2010_amikacin.R)."
      ),
      source_name        = "CRCL"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 32L,
    n_studies      = 1L,
    age_range      = "30.0-67.0 years across the four renal-function groups (study eligibility 18-75 years)",
    age_median     = "group means 50.1, 56.8, 50.8 and 53.1 years for normal, mild, moderate and severe renal impairment",
    weight_range   = "56.0-95.0 kg across the four renal-function groups (study eligibility 55-110 kg)",
    weight_median  = "group means 74.4, 76.5, 78.6 and 71.7 kg for normal, mild, moderate and severe renal impairment",
    sex_female_pct = 43.8,
    race_ethnicity = c(White = 100),
    disease_state  = paste(
      "Open-label, parallel-group, single-dose dedicated renal impairment study in",
      "subjects with varying degrees of kidney function, allocated by Cockcroft-Gault",
      "creatinine clearance to normal kidney function (CLcr > 80 mL/min, n = 8), mild",
      "renal impairment (50 <= CLcr <= 80 mL/min, n = 8), moderate renal impairment",
      "(30 <= CLcr < 50 mL/min, n = 8) or severe renal impairment (CLcr < 30 mL/min,",
      "not on dialysis, n = 8). Subjects were otherwise free of medical conditions",
      "affecting study-drug PK."
    ),
    dose_range     = "Single oral 15 mg edoxaban tablet (27372 nmol) taken after a light breakfast with 240 mL water",
    regions        = "Not stated in the publication",
    renal_function = paste(
      "Full spectrum from normal kidney function to severe renal impairment. A fifth",
      "study group of subjects with end-stage renal disease undergoing peritoneal",
      "dialysis was enrolled but EXCLUDED from this analysis (Jonsson 2015 Methods",
      "'Data': 'group 5 not included in current analysis'); the author's control stream",
      "excludes them with IGNORE=(CRCL.EQ.0). The model therefore should not be",
      "extrapolated to dialysis patients."
    ),
    notes          = paste(
      "Baseline demographics from Jonsson 2015 Results 'Data'. Sex distribution",
      "(male/female) was 6/2, 4/4, 5/3 and 3/5 in the normal, mild, moderate and severe",
      "groups respectively, i.e. 18 male / 14 female overall (43.8% female). All",
      "subjects were White (RACE = 1 in the author's analysis dataset). Observations:",
      "360 edoxaban plasma, 159 edoxaban urine and 294 M4 plasma concentrations; 56",
      "(13.5%) edoxaban and 83 (22.0%) M4 plasma samples were below the lower limit of",
      "quantification and were EXCLUDED from the fit rather than handled by the M3",
      "method, which was numerically unstable in this analysis. Blood sampling predose",
      "and 0.5, 1, 2, 3, 4, 6, 8, 10, 12, 24, 36, 48 and 72 h postdose; urine",
      "collections -12 to 0, 0 to 4, 4 to 8, 8 to 12, 12 to 24 and 24 to 48 h postdose.",
      "Estimation was FOCE with interaction in NONMEM 7.2; precision from MATRIX=S and",
      "from a 200-sample bootstrap stratified by renal-function group."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # All values are the FINAL population estimates from Jonsson 2015
    # Table 1, column "Estimate (%RSE)". The author-supplied NONMEM control
    # stream (run164-share.mod) is the authoritative source for the model
    # STRUCTURE encoded in model() below, but its $THETA / $OMEGA / $SIGMA
    # entries are INITIAL estimates, not final ones -- e.g. it carries
    # slopeCLcr = 0.0862453 whereas the final estimate in Table 1 is 0.109.
    # Where the control stream's initial value happens to agree with the
    # paper to full precision (NN, IIV NN, and the three residual-error
    # terms, all of which were FIXed or converged), the extra digits are
    # noted in the trailing comment.
    # ---------------------------------------------------------------------

    # -- Absorption: Savic 2007 analytical transit-compartment model. The
    #    number of transit compartments NN is non-integer and was FIXED;
    #    the absorption rate constant was not separately identifiable and
    #    was constrained equal to ktr = (NN + 1) / MTT.
    lmtt        <- log(0.731)      ; label("Mean transit time MTT (h, log scale)")                                    # Jonsson 2015 Table 1 'Mean transit time (MTT), h = 0.731' (RSE 16.3%); control stream initial 0.735371
    lnn         <- fixed(log(8.08)); label("Number of transit compartments NN (unitless, log scale)")           # Jonsson 2015 Table 1 'Number of transit compartments (NN) = 8.08 (ne)'; footnote c 'ne, not estimated. Parameter estimate was fixed'; control stream $THETA 8.07871 FIX
    logitfdepot <- log(0.723 / (1 - 0.723)); label("Absolute oral bioavailability F (logit scale; F_pop = 0.723)")     # Jonsson 2015 Table 1 'Absolute oral bioavailability (F) = 0.723' (RSE 8.0%); logit form per control stream PSI = LOG(TVF1/(1-TVF1)); control stream initial 0.739588

    # -- Edoxaban disposition. Clearance is split into a non-renal arm
    #    (which is assumed to form ALL of the M4 metabolite) and a renal arm
    #    (which is a linear function of CRCL and drives urinary excretion).
    #    Jonsson 2015 Table 1 footnote d: 'typical CL = CLNR,typ *
    #    (WT/70)^(3/4) + CLR, where CLR is the slopeClcr * CLcr'.
    lcl_nonren  <- log(10.1)       ; label("Non-renal (M4-forming) edoxaban clearance CLNR at 70 kg (L/h, log scale)") # Jonsson 2015 Table 1 'Nonrenal clearance (CLNR), L/h = 10.1' (RSE 19.4%); control stream initial 10.0954
    lvc         <- log(95.4)       ; label("Edoxaban central volume of distribution Vc at 70 kg (L, log scale)")       # Jonsson 2015 Table 1 'Central volume of distribution (Vc), L = 95.4' (RSE 11.6%); control stream initial 97.0195
    lvp         <- log(54.3)       ; label("Edoxaban peripheral volume of distribution Vp at 70 kg (L, log scale)")    # Jonsson 2015 Table 1 'Peripheral volume of distribution (Vp), L = 54.3' (RSE 16.3%); control stream initial 55.336
    lq          <- log(5.19)       ; label("Edoxaban intercompartmental clearance Q at 70 kg (L/h, log scale)")        # Jonsson 2015 Table 1 'Intercompartmental clearance (Q), L/h = 5.19' (RSE 13.1%); control stream initial 5.24162

    # -- M4 metabolite disposition (one compartment, first-order formation).
    #    Both are APPARENT values: because the fraction of edoxaban
    #    converted to M4 was not identifiable, the model assumes all
    #    non-renally-cleared edoxaban becomes M4 (fm,M4 = 1 - fe =
    #    CLNR/CL), so CLM4 and VM4 are overestimates of the true values
    #    by a factor 1/fm,M4 (Jonsson 2015 Methods 'Modeling').
    lcl_m4      <- log(113)        ; label("Apparent M4 clearance CLM4 at 70 kg (L/h, log scale)")                     # Jonsson 2015 Table 1 'Apparent clearance M4 (CLM4), L/h = 113' (RSE 16.2%); control stream initial 113.391
    lvc_m4      <- log(78.1)       ; label("Apparent M4 volume of distribution VM4 at 70 kg (L, log scale)")           # Jonsson 2015 Table 1 'Apparent volume of distribution M4 (VM4), L = 78.1' (RSE 19.9%); control stream initial 79.2843

    # -- Renal clearance arm. The paper parameterises this as a linear,
    #    non-centred, through-the-origin slope on creatinine clearance:
    #    CLR = slopeCLcr * CLcr with slopeCLcr = 0.109 L/h per mL/min
    #    (Jonsson 2015 Table 1 footnote d). Because the relationship passes
    #    through the origin, anchoring it at any CLcr is an exact
    #    reparameterisation: CLR = (slopeCLcr * CLcr_ref) * (CLcr/CLcr_ref).
    #    It is anchored here at the paper's own typical subject, CLcr = 100
    #    mL/min (Jonsson 2015 Abstract: 'For a typical subject (70 kg;
    #    CLcr, 100 mL/min) ... clearance ... was estimated to ... 21.0
    #    L/h'), so that (a) the parameter carries the canonical
    #    `lcl_renal` name and pairs with the canonical `etalcl_renal` IIV
    #    that Table 1 itself labels 'Renal clearance (CLR), %CV', and (b)
    #    the published slope stays literally visible in the expression
    #    below. At the anchor, cl_renal = 10.9 L/h, so total CL = 10.1 +
    #    10.9 = 21.0 L/h, reproducing the Abstract exactly. Recover the
    #    paper's slope as exp(lcl_renal)/100. NOTE this arm carries NO
    #    allometric weight term (see e_wt_cl_q below).
    lcl_renal   <- log(0.109 * 100); label("Renal clearance CLR at CLcr = 100 mL/min (L/h, log scale)")                # Jonsson 2015 Table 1 'CLcr on CL (slopeCLcr) = 0.109' L/h per mL/min (RSE 11.1%), anchored at the paper's typical CLcr = 100 mL/min; control stream initial slope 0.0862453

    # -- Covariate effects.
    #    Fixed allometric exponents, applied to the clearance terms that are
    #    NOT a function of kidney function (cl_nonren, q, cl_m4) and to all
    #    volume terms (vc, vp, vc_m4).
    e_wt_cl_q       <- fixed(0.75) ; label("Allometric exponent on CLNR, Q and CLM4 (unitless)")                # Jonsson 2015 Results 'Final Population PK Model': exponents 'fixed to 1 and 3/4 for volume and clearance terms'; control stream $THETA 0.75 FIX
    e_wt_vc_vp      <- fixed(1)    ; label("Allometric exponent on Vc, Vp and VM4 (unitless)")                  # Jonsson 2015 Results 'Final Population PK Model': exponents 'fixed to 1 and 3/4 for volume and clearance terms'; control stream $THETA 1 FIX

    # ---------------------------------------------------------------------
    # Inter-individual variability. Jonsson 2015 Results: 'Interindividual
    # variability was included on CLNR, CLR, Vc, Vp, F, MTT, NN, CLM4, and
    # VM4. For CLM4 and VM4 a correlation was estimated. Interindividual
    # variability for NN was fixed in the model.' There is deliberately NO
    # IIV on Q (absent from that list and from the control stream $PK,
    # where Q = TVQ3 with no ETA).
    #
    # BACK-TRANSFORM. Table 1 reports IIV as '%CV' (footnote f: 'CV,
    # coefficient of variation'), but that %CV is the APPROXIMATE
    # standard-deviation-scale CV, i.e. %CV = 100 * sqrt(omega^2), NOT the
    # exact log-normal 100 * sqrt(exp(omega^2) - 1). This is pinned down
    # unambiguously by NN, whose IIV was FIXED and therefore cannot have
    # moved between the control stream and the final fit: the control
    # stream has $OMEGA 1.50033 FIX and Table 1 reports 122% CV, and
    # sqrt(1.50033) = 1.2249 -> 122% whereas sqrt(exp(1.50033) - 1) =
    # 1.8663 -> 187%. Footnote a corroborates that IIV precision is
    # 'reported on the approximate standard deviation scale'. Every
    # variance below is therefore (%CV / 100)^2.
    # ---------------------------------------------------------------------
    etalcl_nonren     ~ 0.031329   # var = 0.177^2 ; Jonsson 2015 Table 1 'CLNR, %CV = 17.7' (RSE 45.3%); control stream initial 0.032937
    etalcl_renal      ~ 0.073441   # var = 0.271^2 ; Jonsson 2015 Table 1 'Renal clearance (CLR), %CV = 27.1' (RSE 35.0%); control stream initial 0.0588206
    etalvc            ~ 0.027225   # var = 0.165^2 ; Jonsson 2015 Table 1 'Vc, %CV = 16.5' (RSE 32.7%); control stream initial 0.0281244
    etalvp            ~ 0.123201   # var = 0.351^2 ; Jonsson 2015 Table 1 'Vp, %CV = 35.1' (RSE 23.2%); control stream initial 0.123084
    etalmtt           ~ 0.283024   # var = 0.532^2 ; Jonsson 2015 Table 1 'MTT, %CV = 53.2' (RSE 27.9%); control stream initial 0.284762
    etalnn            ~ fixed(1.4884) # var = 1.22^2 ; Jonsson 2015 Table 1 'NN, %CV = 122 (ne)' -- fixed, footnote c; control stream $OMEGA 1.50033 FIX

    # IIV on F is on the LOGIT scale (control stream: VF1 = EXP(PSI +
    # ETA(5)); BIO = VF1/(1 + VF1)). Table 1 reports 0.129 for F, which per
    # footnote g is the delta-method standard deviation on the NATURAL F
    # scale: 'Because of logit transformation, variability is standard
    # deviation (SD) = F * (1-F) * vF'. Inverting for the logit-scale SD:
    # sqrt(omega^2) = 0.129 / (0.723 * (1 - 0.723)) = 0.644127, so
    # omega^2 = 0.414900. Cross-check with the control stream's initial
    # $OMEGA 0.446911: F * (1-F) * sqrt(0.446911) = 0.1339, consistent with
    # the reported 0.129 after final estimation.
    etalogitfdepot    ~ 0.414900   # Jonsson 2015 Table 1 'F = 0.129' (RSE 43.1%) with footnote g; control stream initial 0.446911

    # Correlated IIV on the two M4 disposition parameters. Table 1 reports
    # 'Correlation CLM4-VM4 = 0.848', so the covariance is
    # 0.848 * sqrt(0.248004 * 0.180625) = 0.179479.
    etalcl_m4 + etalvc_m4 ~ c(0.248004, 0.179479, 0.180625)  # var_clm4 = 0.498^2, cov from corr 0.848, var_vcm4 = 0.425^2 ; Jonsson 2015 Table 1 'CLM4, %CV = 49.8', 'Correlation CLM4-VM4 = 0.848', 'VM4, %CV = 42.5'; control stream $OMEGA BLOCK(2) initials 0.241464 / 0.175534 / 0.17453

    # ---------------------------------------------------------------------
    # Residual variability: proportional on all three measured entities
    # (Jonsson 2015 Results: 'The residual variability for each entity
    # (edoxaban plasma and urine concentrations, M4 plasma concentrations)
    # was modeled as proportional residual error models'). The reported
    # '%CV' equals 100 * sqrt(sigma^2), which matches the control stream
    # $SIGMA values to three digits.
    #
    # NOT ENCODED: Table 1 also reports a correlation of 0.232 between the
    # edoxaban-plasma and M4-plasma residual errors (a $SIGMA BLOCK(2) in
    # the control stream, estimated with the NONMEM L2 data item). nlmixr2
    # has no idiomatic encoding for cross-endpoint residual correlation, so
    # the two proportional errors are independent here. Precedent:
    # Svensson_2014_bedaquiline.R drops the analogous 55% BDQ/M2 residual
    # correlation. See the vignette's Assumptions and deviations section.
    # ---------------------------------------------------------------------
    propSd        <- 0.154 ; label("Edoxaban plasma proportional residual error (SD, fraction)")   # Jonsson 2015 Table 1 'Proportional residual error, edoxaban plasma concentrations, %CV = 15.4' (RSE 5.9%); control stream $SIGMA 0.023842, sqrt = 0.1544
    propSd_m4     <- 0.184 ; label("M4 plasma proportional residual error (SD, fraction)")         # Jonsson 2015 Table 1 'Proportional residual error, M4 plasma concentrations, %CV = 18.4' (RSE 8.3%); control stream $SIGMA 0.0337782, sqrt = 0.1838
    propSd_Aurine <- 0.254 ; label("Edoxaban urine proportional residual error (SD, fraction)")    # Jonsson 2015 Table 1 'Proportional residual error, edoxaban urine concentrations, %CV = 25.4' (RSE 8.1%); control stream $SIGMA 0.0649608, sqrt = 0.2549
  })

  model({
    # 1. Fixed allometric size scaling around a 70 kg reference subject.
    #    Jonsson 2015 Results: 'Typical Vc = Vc,typ * (WT/70)^1' and
    #    'Typical Q = Q_typ * (WT/70)^(3/4)'.
    allcl <- (WT / 70)^e_wt_cl_q
    allv  <- (WT / 70)^e_wt_vc_vp

    # 2. Individual absorption parameters. ktr is the transit rate constant
    #    of the Savic 2007 chain; the depot empties into central at this
    #    same rate because the absorption rate constant was not separately
    #    identifiable (Jonsson 2015 Results: 'The transit chain was
    #    therefore allowed to deliver the amount of drug into the central
    #    compartment without a separate absorption rate constant, that is,
    #    the absorption rate constant was set to the same value as the rate
    #    constant governing the transit through the compartments (ktr)').
    mtt    <- exp(lmtt + etalmtt)
    nn     <- exp(lnn  + etalnn)
    ktr    <- (nn + 1) / mtt

    # 3. Absolute oral bioavailability, logit-normal so every individual F
    #    stays inside (0, 1).
    fdepot <- expit(logitfdepot + etalogitfdepot)

    # 4. Individual disposition parameters. The renal arm carries no
    #    allometric term (allometry was applied only to clearance terms not
    #    driven by kidney function); it is a pure linear function of CRCL
    #    through the origin. Total edoxaban clearance is the sum of the two
    #    arms, matching Table 1 footnote d.
    cl_nonren <- exp(lcl_nonren + etalcl_nonren) * allcl
    cl_renal  <- exp(lcl_renal  + etalcl_renal)  * (CRCL / 100)
    cl        <- cl_nonren + cl_renal
    vc        <- exp(lvc  + etalvc) * allv
    vp        <- exp(lvp  + etalvp) * allv
    q         <- exp(lq)            * allcl
    cl_m4     <- exp(lcl_m4 + etalcl_m4) * allcl
    vc_m4     <- exp(lvc_m4 + etalvc_m4) * allv

    # 5. ODE system (author's control stream $MODEL / $DES, ADVAN13:
    #    DEPOT / EDOCP / EDOPER / MET / EDOCU).
    #    - transit() is the rxode2 built-in for the Savic 2007 analytical
    #      gamma-density transit input; it reproduces the control stream's
    #      hand-coded DADT(1) input term, whose LNFAC is the Stirling
    #      approximation of lgamma(NN + 1).
    #    - The non-renal clearance flux leaves edoxaban central and enters
    #      M4 central 1:1 in MOLAR units: all edoxaban not renally excreted
    #      is assumed to form M4 (fm,M4 = 1 - fe = CLNR/CL).
    #    - The renal clearance flux leaves edoxaban central and accumulates
    #      in the urine compartment.
    d/dt(depot)       <- transit(nn, mtt, fdepot) - ktr * depot
    d/dt(central)     <- ktr * depot -
                           cl_nonren * central / vc -
                           cl_renal  * central / vc -
                           q * central / vc + q * peripheral1 / vp
    d/dt(peripheral1) <- q * central / vc - q * peripheral1 / vp
    d/dt(central_m4)  <- cl_nonren * central / vc - cl_m4 * central_m4 / vc_m4
    d/dt(urine)       <- cl_renal  * central / vc

    # 6. Bioavailability. The dose amount is read by transit() via podo(),
    #    which delivers F * dose smoothly through the transit chain, so the
    #    bolus contribution into depot must be suppressed. This mirrors the
    #    control stream's 'F1 = 0' with BIO carried inside the DADT(1)
    #    input expression.
    f(depot) <- 0

    # 7. Observations. All three outputs are declared BEFORE any residual
    #    error line so that the ODE states they reference are not parsed
    #    into the error block.
    #    Cc     : edoxaban plasma concentration (nmol/L)
    #    Cc_m4  : M4 plasma concentration (nmol/L)
    #    Aurine : CUMULATIVE amount of edoxaban excreted in urine (nmol).
    #             The paper fits interval urine CONCENTRATION (amount
    #             excreted during a collection interval divided by the
    #             recorded interval urine volume UVOL, with the urine
    #             compartment reset at each interval). A proportional
    #             residual error is invariant to that known-volume
    #             division, so propSd_Aurine transfers unchanged; the
    #             interval amount is recovered by differencing Aurine
    #             across collection boundaries. See the vignette's
    #             Assumptions and deviations section.
    Cc     <- central    / vc
    Cc_m4  <- central_m4 / vc_m4
    Aurine <- urine

    Cc     ~ prop(propSd)
    Cc_m4  ~ prop(propSd_m4)
    Aurine ~ prop(propSd_Aurine)
  })
}
