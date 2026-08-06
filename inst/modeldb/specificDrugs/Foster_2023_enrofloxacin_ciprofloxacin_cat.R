Foster_2023_enrofloxacin_ciprofloxacin_cat <- function() {
  description <- paste(
    "Veterinary (domestic cat). Joint parent + metabolite population PK model for",
    "enrofloxacin and its active metabolite ciprofloxacin after a single 5 mg/kg",
    "intravenous enrofloxacin dose infused over 30 minutes in 34 client-owned cats",
    "hospitalised for clinical illness and spanning normal to markedly reduced kidney",
    "function (Foster 2023 Figure S2 and Table 3). Each analyte gets one plasma",
    "compartment: enrofloxacin is eliminated by a renal clearance into urine plus a",
    "metabolic formation clearance that generates ciprofloxacin, and ciprofloxacin is",
    "eliminated by its own clearance into urine. Retained covariates are body weight",
    "on the enrofloxacin volume and blood urea nitrogen on the metabolic formation",
    "clearance -- the paper's finding that more azotemic cats form ciprofloxacin",
    "faster, which the authors themselves argue is unlikely to be causal because",
    "serum creatinine and SDMA were not significant. This is a SECOND, independent",
    "fit reported alongside the paper's primary total-fluoroquinolone model",
    "(Foster_2023_enrofloxacin_totalFluoroquinolone_cat), not a sub-model of it.",
    "Two caveats carried from the source: the volumes and clearances are ABSOLUTE",
    "(mL, mL/h) despite Table 3's 'mL/kg' labels, and the published clearance pair",
    "reproduces a ciprofloxacin fraction of total fluoroquinolone of only about 2%",
    "against the ~18% the same paper reports -- see the lvc_cipro / lcl_met comments",
    "and the vignette Errata. Between-subject variances and the residual-error",
    "magnitude are not reported, so every eta and residual SD is fixed(0) and the",
    "model simulates typical-value profiles unless the user supplies variances.",
    sep = " "
  )
  reference <- paste(
    "Foster JD, Abouraya M, Papich MG, Muma NA. (2023).",
    "Population pharmacokinetic analysis of enrofloxacin and its active metabolite",
    "ciprofloxacin after intravenous injection to cats with reduced kidney function.",
    "Journal of Veterinary Internal Medicine 37(6):2230-2240.",
    "doi:10.1111/jvim.16866.",
    sep = " "
  )
  vignette <- "Foster_2023_enrofloxacin_cat"
  units <- list(
    time = "h",
    dosing = "ug",
    concentration = "ug/mL"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Both states are plasma compartments in Foster 2023
  # Figure S2 ("Volume Enro / Concentration Enro" and "Volume Cipro /
  # Concentration Cipro"); the urine sinks of that figure are not carried
  # as states because no urine data were collected.
  compartmentData <- list(
    central = list(
      analyte = "enrofloxacin", units = "ug", specimen = "plasma", verified = TRUE
    ),
    central_cipro = list(
      analyte = "ciprofloxacin", units = "ug", specimen = "plasma", verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed baseline body weight. Retained on the enrofloxacin volume of",
        "distribution only (Table 3 dVd Enrofloxacin-Weight; Figure 5 row 'nVparent,",
        "Weight' shows the positive trend). Foster 2023 states the centring value",
        "only for the companion total-fluoroquinolone model ('body weight was",
        "centered around the mean population body weight, 3.8 kg', Results 3.2); the",
        "same 3.8 kg centring is applied here because it is the same cohort, the same",
        "software and the same covariate. Entered as a power of the ratio",
        "(WT / 3.8)^exponent, Phoenix NLME's standard continuous-covariate form.",
        "Cohort median 3.8 kg, range 1.8-8.8 kg.",
        sep = " "
      ),
      source_name        = "weight (Foster 2023 Table 3 covariate name dVd Enrofloxacin-Weight)"
    ),
    BUN = list(
      description        = "Blood urea nitrogen concentration",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed baseline BUN. Retained on the metabolic formation clearance of",
        "ciprofloxacin from enrofloxacin (Table 3 'dCl Enrofloxacin metabolism to",
        "ciprofloxacin BUN'; Figure 5 row 'nCl_EnroCipro, BUN' shows the positive",
        "trend). Foster 2023 does NOT state a centring value for BUN; the cohort mean",
        "of 137.6 mg/dL (Results section 3) is used because the paper centres its only",
        "other continuous covariate on the population mean. Entered as a power of the",
        "ratio (BUN / 137.6)^exponent. Cohort mean 137.6 +/- 103.0 mg/dL; the Figure 5",
        "x-axes span roughly 0-400 mg/dL. The authors caution in the Discussion that",
        "BUN is a poor GFR marker and that the association is unlikely to be",
        "mechanistic, since serum creatinine and SDMA were not significant on the same",
        "parameter.",
        sep = " "
      ),
      source_name        = "BUN (Foster 2023 Table 3)"
    )
  )

  # Screened in the stepwise covariate search for this model (Foster 2023
  # Results 3.4, Figure 5: age, BUN, creatinine, SDMA and weight against the
  # empirical Bayes etas of all five structural parameters) but not retained,
  # and with no published point estimate. BUN and WT are absent from this
  # list because they ARE retained (see covariateData above).
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age", units = "years", type = "continuous",
      notes = "Screened against all five etas (Figure 5 column 1); not retained. Cohort mean 11.4 +/- 4.3 years."
    ),
    CREAT = list(
      description = "Serum creatinine concentration", units = "mg/dL", type = "continuous",
      notes = "Screened against all five etas (Figure 5 column 3); not retained. Cohort mean 7.5 +/- 6.0 mg/dL. The Discussion leans on this non-significance to argue the BUN effect on metabolism is not GFR-mediated."
    ),
    SDMA = list(
      description = "Symmetric dimethylarginine concentration", units = "ug/dL", type = "continuous",
      notes = "Screened against all five etas (Figure 5 column 4); not retained. Cohort mean 39.8 +/- 27.0 ug/dL. SDMA has no canonical entry in inst/references/covariate-columns.md; none is proposed here because the covariate is documentation only and carries no published point estimate."
    )
  )

  population <- list(
    species        = "cat (domestic, client-owned)",
    n_subjects     = 34L,
    n_studies      = 1L,
    age_range      = "mean 11.4 +/- 4.3 years",
    weight_range   = "1.8-8.8 kg (median 3.8 kg)",
    weight_median  = "3.8 kg",
    sex_female_pct = 52.9,
    disease_state  = paste(
      "Hospitalised for clinical illness with suspected bacterial infection and",
      "variable kidney function: 9 cats with normal kidney function (median serum",
      "creatinine 1.6 mg/dL), 15 with moderate kidney dysfunction (median 5.6) and",
      "10 with severe kidney dysfunction (median 13.3).",
      sep = " "
    ),
    renal_function = "Serum creatinine mean 7.5 +/- 6.0 mg/dL; BUN mean 137.6 +/- 103.0 mg/dL; SDMA mean 39.8 +/- 27.0 ug/dL",
    dose_range     = "Single 5 mg/kg enrofloxacin (Baytril 2.27%, Elanco) diluted 1:1 with sterile saline and infused intravenously over 30 minutes; mean total dose 20.5 +/- 7.4 mg",
    sampling       = "Sparse: 3 plasma samples per cat in the 24 h after dosing (98 samples from 34 cats); enrofloxacin and ciprofloxacin quantified separately by validated HPLC, which is what makes this parent + metabolite fit possible.",
    regions        = "United States (single referral hospital, enrolment 2019-01-01 to 2021-09-01)",
    notes          = paste(
      "Foster 2023 Results 3.4. Table 2 is the base model of this structure and",
      "Table 3 the final covariate model; the two differ only in the two added",
      "covariates and in rounding. This model was fitted to the separately measured",
      "enrofloxacin and ciprofloxacin concentrations, whereas the companion",
      "Foster_2023_enrofloxacin_totalFluoroquinolone_cat model was fitted to their",
      "sum. No urine was collected, so both 'clearance into urine' parameters are",
      "plasma-side estimates of routes the authors label renal.",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # UNITS -- same erratum as the companion total-fluoroquinolone model
    # =====================================================================
    # Tables 2 and 3 print "(mL/kg)" and "(mL/kg/h)" but the parameters are
    # ABSOLUTE. With the mean total dose of 20.5 mg = 20500 ug the absolute
    # reading gives C_enro(0) = 20500 / 5797.5 = 3.5 ug/mL, consistent with
    # Figure 1's total-fluoroquinolone concentrations of ~3-8 ug/mL just
    # after the infusion; the weight-normalised reading gives 5000 / 5797.5
    # = 0.86 ug/mL, four-fold too low. The elimination half-life implied by
    # the absolute reading, log(2) * 5797.5 / (499.1 + 0.13) = 8.05 h, also
    # matches the ~9 h terminal slope of Figure 1. Doses are therefore
    # entered in ug and volumes are in mL. See vignette Errata.

    # =====================================================================
    # Enrofloxacin (parent) -- Foster 2023 Table 3
    # =====================================================================
    lvc <- log(5797.5); label("Enrofloxacin volume of distribution for a 3.8 kg cat (mL)")  # Table 3: theta Vd Enrofloxacin = 5797.5 (SE 313.5, CV% 5.4); Table 2 base model gives the same 5797.5 (SE 323.5)
    lcl_nonmet <- log(499.1); label("Enrofloxacin clearance into urine (mL/h)")             # Table 3: Clearance of enrofloxacin into urine = 499.1 (SE 42.6, CV% 8.5); identical in Table 2

    # Metabolic formation clearance of ciprofloxacin from enrofloxacin.
    # Table 3 prints 0.1 and Table 2 prints 0.13 for the same parameter with
    # the same SE (0.03) and the same CV% (22.9). 0.03 / 0.13 = 23.1%, which
    # rounds to the printed 22.9%, whereas 0.03 / 0.1 = 30.0% does not, so
    # the underlying estimate is 0.13 and Table 3 simply prints it to one
    # decimal. The Table 2 value is used here.
    lcl_met <- log(0.13); label("Metabolic formation clearance of ciprofloxacin from enrofloxacin, for the mean BUN of 137.6 mg/dL (mL/h)")  # Table 2: Metabolism of enrofloxacin to ciprofloxacin = 0.13 (SE 0.03, CV% 22.9); Table 3 prints the same parameter rounded to 0.1 with the same SE and CV%

    # =====================================================================
    # Ciprofloxacin (metabolite) -- Foster 2023 Table 3
    # =====================================================================
    # Two properties of this pair are worth flagging before use.
    # (a) 22.6 mL is not a physically possible distribution volume for a
    #     ~4 kg cat (plasma volume alone is roughly 180 mL). It is a
    #     fitted apparent volume that also absorbs the fraction of the
    #     enrofloxacin dose actually converted, which the sparse design
    #     cannot identify separately from the formation clearance.
    # (b) At quasi-steady state this model gives
    #       C_cipro / C_enro = cl_met / cl_cipro = 0.13 / 4.6 = 0.028,
    #     i.e. ciprofloxacin is ~2.8% of total fluoroquinolone, whereas
    #     Table 1 of the same paper reports a ciprofloxacin:quinolone ratio
    #     of 0.18 (SE 0.02). That ratio depends only on cl_met / cl_cipro,
    #     so no choice of volume reconciles the two; matching 0.18 would
    #     need cl_met / cl_cipro near 0.22. The published values are
    #     encoded unchanged -- tuning them to hit the reported ratio is
    #     not permitted -- and the vignette quantifies the gap.
    lvc_cipro <- log(22.6); label("Ciprofloxacin apparent volume of distribution (mL)")  # Table 3: Vd ciprofloxacin = 22.6 (SE 4.2, CV% 18.7); Table 2 base model gives 22.7
    lcl_cipro <- log(4.6); label("Ciprofloxacin clearance into urine (mL/h)")            # Table 3: Clearance of ciprofloxacin into urine = 4.6 (SE 0.8, CV% 16.9); identical in Table 2

    # =====================================================================
    # Covariate effects -- Foster 2023 Abstract and Table 3
    # =====================================================================
    # Both are entered as powers of the ratio to a centring value, Phoenix
    # NLME's standard continuous-covariate form. The centred-exponential
    # alternative is impossible for either: exp(0.7 * (8.8 - 3.8)) = 33-fold
    # on the volume of the largest cat, and exp(0.51 * (400 - 137.6)) = 10^58
    # on the formation clearance at the top of Figure 5's BUN axis.
    e_wt_vc <- 0.7; label("Body-weight power exponent on the enrofloxacin volume (unitless)")  # Table 3: dVd Enrofloxacin-Weight = 0.7 (SE 0.15, CV% 20.0)
    # The Abstract prints this effect to two decimals ("effect 0.51, SE 0.08")
    # where Table 3 rounds it to 0.5; 0.08 / 0.51 = 15.7% reproduces the
    # printed CV% of 15.5 better than 0.08 / 0.5 = 16.0%, so 0.51 is used.
    e_bun_cl_met <- 0.51; label("Blood-urea-nitrogen power exponent on the metabolic formation clearance (unitless)")  # Abstract: "Blood urea nitrogen concentration influenced the metabolic generation of ciprofloxacin from enrofloxacin (effect 0.51, SE 0.08, P < .01)"; Table 3 dCl ... BUN = 0.5 (SE 0.08)

    # =====================================================================
    # Between-subject variability
    # =====================================================================
    # Foster 2023 Methods paragraph 7 specifies exponential IIV, and Figure 5
    # plots the empirical Bayes etas of all five structural parameters
    # (nVparent, nCl_Enro_urine, nCl_EnroCipro, nVcipro, nCl_cipro_urine), so
    # each parameter carried an eta. No omega, omega^2, SD or CV of any eta
    # is published, and the bootstrap supplement (Table S2) covers only the
    # companion total-fluoroquinolone model. Per the standing policy for
    # unreported IIV with structural values present, each eta is fixed(0)
    # rather than invented.
    etalvc ~ fixed(0)
    etalcl_nonmet ~ fixed(0)
    etalcl_met ~ fixed(0)
    etalvc_cipro ~ fixed(0)
    etalcl_cipro ~ fixed(0)

    # =====================================================================
    # Residual unexplained variability
    # =====================================================================
    # Methods paragraph 7: multiplicative (proportional in nlmixr2 terms),
    # Cobs_ij = Cpred_ij * (1 + eps_ij). No sigma value is published for
    # either analyte, so both are fixed(0) and simulations are noise-free
    # typical-value profiles.
    propSd <- fixed(0); label("Proportional residual error on enrofloxacin (fraction; not published)")           # Methods paragraph 7 states the multiplicative form; no sigma value is given
    propSd_cipro <- fixed(0); label("Proportional residual error on ciprofloxacin (fraction; not published)")    # Methods paragraph 7 states the multiplicative form; no sigma value is given
  })

  model({
    # 1. Individual parameters. Body weight scales the enrofloxacin volume
    #    and BUN scales the metabolic formation clearance; the two urinary
    #    clearances and the ciprofloxacin volume carry no covariate
    #    (Foster 2023 Table 3).
    vc <- exp(lvc + etalvc) * (WT / 3.8)^e_wt_vc
    cl_nonmet <- exp(lcl_nonmet + etalcl_nonmet)
    cl_met <- exp(lcl_met + etalcl_met) * (BUN / 137.6)^e_bun_cl_met
    vc_cipro <- exp(lvc_cipro + etalvc_cipro)
    cl_cipro <- exp(lcl_cipro + etalcl_cipro)

    # 2. ODE system (Foster 2023 Figure S2). Enrofloxacin leaves its plasma
    #    compartment by two parallel routes -- clearance into urine and
    #    metabolic conversion -- and the conversion flux is the input to the
    #    ciprofloxacin compartment, which is drained by its own clearance
    #    into urine. The intravenous dose (30-minute infusion) goes into
    #    `central`.
    d/dt(central) <- -(cl_nonmet + cl_met) / vc * central
    d/dt(central_cipro) <- cl_met / vc * central - cl_cipro / vc_cipro * central_cipro

    # 3. Observations. Amounts are in ug and volumes in mL, so both
    #    concentrations come out in ug/mL. Their sum is the "total
    #    fluoroquinolone" concentration that the companion model describes
    #    directly.
    Cc <- central / vc
    Cc_cipro <- central_cipro / vc_cipro

    Cc ~ prop(propSd)
    Cc_cipro ~ prop(propSd_cipro)
  })
}
