Melander_2025_cetirizine <- function() {
  description <- "One-compartment population PK model with first-order absorption describing cetirizine concentrations in the breast milk of lactating women, built to predict milk exposure and the relative infant dose (RID) in a breastfed infant. Only breast-milk samples were collected (no maternal plasma), so every disposition parameter is an APPARENT milk parameter that lumps oral bioavailability, systemic disposition and plasma-to-milk transfer into a single empirical milk compartment; the model describes the milk concentration-time course and must not be read as systemic PK. Inter-individual variability on apparent milk clearance only. The duration of the breast-milk pumping session used to collect the sample is the single retained covariate and increases the apparent milk volume of distribution, consistent with fat-content enrichment of milk over a longer expression."
  reference <- paste(
    "Melander E, Nielsen EI, Lindqvist A, Hovd M, Gandia P, Panchaud A,",
    "Guidi M, Annaert P, Baranczewski P, Spigset O, Nordeng H.",
    "Population pharmacokinetic modelling of cetirizine concentrations in",
    "human breast milk - A contribution from the ConcePTION project.",
    "Basic Clin Pharmacol Toxicol. 2025;136(1):e14100.",
    "doi:10.1111/bcpt.14100.",
    "All parameter values are the final estimates of Table 1; the structural",
    "ODE and the covariate-parameter relationship are the display equations",
    "of Section 3 (Results).",
    sep = " "
  )
  vignette <- "Melander_2025_cetirizine"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  covariateData <- list(
    T_PUMP = list(
      description        = "Duration of the breast-milk pumping (expression) session used to collect the milk sample",
      units              = "h",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The only covariate retained by the stepwise covariate-modelling",
        "search (Melander 2025 Results section 3 and Table 2); maternal age,",
        "maternal body weight, BMI, breastfeeding exclusivity and infant age",
        "were all screened and rejected. Enters the apparent milk volume of",
        "distribution as a linear effect centred at T_PUMP = 0.22 h",
        "(13.2 min), the value at which the covariate factor is 1 and Vm",
        "equals the tabulated typical value of 19.9 L. The paper reports the",
        "covariate in HOURS in the equation text ('TPUMP is the recorded",
        "pumping duration in the individual in hours') while the cohort",
        "summary is quoted in minutes (mean 14.9 min, range 2-50 min, i.e.",
        "0.033-0.833 h). Source column TPUMP.",
        sep = " "
      ),
      source_name        = "TPUMP"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Maternal age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened by the stepwise covariate model and not retained (Melander 2025 Methods section 2.3 and Table 2, 'No significant effects were found from maternal characteristics'). Cohort mean 30 years (range 22-40)."
    ),
    WT = list(
      description = "Maternal body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened by the stepwise covariate model and not retained (Melander 2025 Results section 3, 'No significant effects were found from maternal characteristics such as body weight or BMI'). Cohort mean 78 kg (range 53-110). Retained here as documentation because maternal weight is nevertheless needed to compute the relative infant dose, which normalises the maternal dose per kg of maternal body weight."
    ),
    BMI = list(
      description = "Maternal body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened by the stepwise covariate model and not retained (Melander 2025 Results section 3 and Table 2). No cohort summary reported."
    ),
    AGE_INFANT = list(
      description = "Postnatal age of the breastfed infant",
      units       = "months",
      type        = "continuous",
      notes       = "Screened as a proxy for milk maturity and not retained (Melander 2025 Methods section 2.3, 'Age of the child gave an indication of the maturity of the milk, which potentially could impact milk composition', and Results section 3, 'nor was the age of the infant found to be a significant covariate'). Cohort mean 8.2 months (range 8.3 weeks to 21 months)."
    ),
    WT_INFANT = list(
      description = "Body weight of the breastfed infant",
      units       = "kg",
      type        = "continuous",
      notes       = "Not a covariate on any PK parameter. Reported for the cohort (mean 8.3 kg, range 3.7-11.9 kg) and used outside the PK model in the relative-infant-dose calculation, where the assumed daily milk intake is 150 mL per kg of infant body weight (Melander 2025 Methods section 2.4); because both the absolute infant dose and the maternal dose in the RID ratio are expressed per kg of the respective body weight, infant weight cancels out of the RID itself."
    ),
    BF_EXCLUSIVE = list(
      description = "Exclusive-breastfeeding indicator (yes / no)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened by the stepwise covariate model as a breastfeeding-pattern covariate and not retained (Melander 2025 Methods section 2.3, 'breastfeeding exclusivity (yes/no, categorical variable)'). No point estimate is reported, so no effect can be encoded. Not registered in inst/references/covariate-columns.md because the model does not use it."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "cetirizine", units = "mg", specimen = "administration site", verified = TRUE),
    milk  = list(analyte = "cetirizine", units = "mg", specimen = "milk", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 35,
    n_studies      = 1,
    n_observations = 205,
    age_range      = "22-40 years",
    age_mean       = "30 years",
    weight_range   = "53-110 kg",
    weight_mean    = "78 kg",
    sex_female_pct = 100,
    disease_state  = "Healthy lactating (breastfeeding) women taking cetirizine or levocetirizine for allergic conditions; no disease state was required for enrolment.",
    dose_range     = "Cetirizine 10 mg once daily in 32 of 35 women; 20 mg once daily in 2 women; levocetirizine 5 mg once daily in 1 woman. All women had dosed for at least two days before the sampling day so that milk sampling was at steady state.",
    regions        = "Norway (mothers resident in Norway; samples analysed at Uppsala University, Sweden). Regional Committee for Medical and Health Research Ethics in South-East Norway, REK no. 232351.",
    infant_partner = "Each mother had one breastfed infant: mean age 8.2 months (range 8.3 weeks to 21 months), mean body weight 8.3 kg (range 3.7-11.9 kg). Infant characteristics enter the relative-infant-dose calculation, not the PK model.",
    feeding_pattern = "Relative infant dose assumes a daily milk intake of 150 mL per kg of infant body weight, per the FDA clinical-lactation-study guidance (Melander 2025 Methods section 2.4).",
    notes          = "Baseline demographics from Melander 2025 Methods section 2.1; the parent lactation study is Nordeng et al. Each woman provided 4-6 milk samples, nominally pre-dose (time 0) and 2, 4, 8, 12 and 24 h after dose intake, self-collected at home with a supplied electric breast pump. The women were instructed to pump until the breast felt empty and to retain 20 mL of the full pumped volume, so each sample is a whole-breast average rather than fore- or hind-milk. Mean pumping duration 14.9 min (range 2-50 min). Assay: LC-MS/MS with LLOQ 0.39 ug/L and limit of detection 0.04 ug/L; 17 detectable samples from three women fell between the two limits and were carried into the fit with their own additive residual-error term rather than being censored."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters. All values are the FINAL estimates of
    # Melander 2025 Table 1 ("Final parameter estimates of cetirizine
    # breast milk popPK model, n = 205 samples"). Every disposition
    # parameter is apparent: the study collected milk only, so
    # bioavailability, systemic disposition and plasma-to-milk transfer
    # are all folded into these numbers (Discussion, "all pharmacokinetic
    # parameters are apparent parameters").
    #
    # The volume below is the value at the covariate reference
    # T_PUMP = 0.22 h; see e_t_pump_vc and the model() block.
    # ------------------------------------------------------------------
    lka <- log(0.095); label("Apparent first-order absorption rate constant into breast milk, ka (1/h)")                                       # Table 1 'ka (h-1) 0.095, RSE 4.5%, 95% CI 0.088-0.110'
    lcl <- log(25.2);  label("Apparent clearance of cetirizine from breast milk, CL/F_milk (L/h)")                                             # Table 1 'Apparent Clearance milk (L/h) 25.2, RSE 6.7%, 95% CI 22.8-29.7'
    lvc <- log(19.9);  label("Apparent volume of distribution of the milk compartment at the reference pumping duration of 0.22 h, Vm (L)")    # Table 1 'Apparent VM (L) 19.9, RSE 11.6%, 95% CI 15.6-24.8'

    # ------------------------------------------------------------------
    # Covariate effect of pumping duration on the apparent milk volume.
    #
    # SOURCE EQUATION AS PRINTED (Melander 2025, Results section 3):
    #     VTpump = (1 + VTpump1) x (TPUMP - 0.22)
    #     Vm     = VTpump x TVV
    # As printed this is arithmetically impossible: it makes Vm exactly
    # zero at TPUMP = 0.22 h and NEGATIVE for every shorter pumping
    # duration, and 13 of the study's pumping durations are shorter than
    # 13.2 min (cohort range 2-50 min). It also destroys the reported
    # typical value: at the cohort mean pumping duration of 14.9 min it
    # returns Vm = 2.5 L rather than the tabulated 19.9 L.
    #
    # ENCODED HERE (misplaced-parenthesis correction, the standard
    # centred linear NONMEM covariate form):
    #     VTpump = 1 + VTpump1 x (TPUMP - 0.22)
    #     Vm     = VTpump x TVV
    # Three independent checks confirm this reading; the arithmetic is
    # reproduced in the vignette's source-trace section:
    #   (1) at the reference TPUMP = 0.22 h the factor is exactly 1, so
    #       Vm = 19.9 L, the value Table 1 tabulates;
    #   (2) with ka = 0.095/h, CL = 25.2 L/h and Vm = 19.9 L the
    #       analytical single-dose peak is Cmax = 30.56 ug/L, matching
    #       Table 1's calculated Cmax of 30.6 ug/L to three significant
    #       figures, and Tmax = 2.21 h against a reported 2.53 +/- 1.09 h;
    #   (3) sweeping T_PUMP over the observed 2-50 min range gives
    #       Tmax from 1.04 h to 4.68 h, which brackets the reported Tmax
    #       range of 1.23-5.83 h once the 12.2% CV on clearance is added.
    #       The alternative power reading, VTpump = (1 + VTpump1)^(TPUMP
    #       - 0.22), spans only 1.82-4.07 h and cannot reach the reported
    #       lower bound of 1.23 h at any pumping duration in range.
    # Table 1 labels this coefficient's units "(L)"; under the encoded
    # form it is a fractional change per hour of pumping, i.e. 1/h. The
    # "(L)" label appears to be carried down from the Vm row above it.
    # ------------------------------------------------------------------
    e_t_pump_vc <- 3.47; label("Linear effect of pumping duration on the apparent milk volume of distribution, per hour of pumping centred at 0.22 h (1/h)")  # Table 1 'Covariate effect of pumping duration on VM 3.47, RSE 24.1%, 95% CI 1.86-5.20'

    # ------------------------------------------------------------------
    # Inter-individual variability. IIV was estimated on apparent milk
    # clearance only ("IIV was included in CL", Results section 3); no IIV
    # was reported on ka or on Vm. Table 1 reports it already converted to
    # a coefficient of variation, so the internal variance is
    # omega^2 = log(CV^2 + 1) = log(0.122^2 + 1) = 0.01477.
    # ------------------------------------------------------------------
    etalcl ~ 0.01477   # Table 1 'Inter-individual variability in Clearance (%CV) 12.2, RSE 21.1%, shrinkage 5.3%'

    # ------------------------------------------------------------------
    # Residual error. The source fitted a combined additive-plus-
    # proportional model "with the additive part only included for BLQ
    # data" (Results section 3): the additive term exists to downweight
    # the 17 detectable-but-below-LLOQ samples so they could be used
    # rather than censored (dOFV -194).
    #
    # nlmixr2 has no per-record switch on a residual-error component, so
    # both components are applied to every record here. The approximation
    # is close in the direction that matters: the additive term is
    # negligible where concentrations are high (at the 30.6 ug/L peak it
    # raises the residual SD from 10.7 to 11.5 ug/L, +8%) and dominates
    # only at the low concentrations where the BLQ samples actually live,
    # which is the behaviour the source model encodes explicitly. The
    # deviation is confined to non-BLQ records near the LLOQ.
    #
    # Both values are read as standard deviations, not variances: Table 1
    # gives the additive term in ug/L rather than (ug/L)^2, and reports
    # the IIV already back-transformed to %CV, so the table presents
    # interpretable scales throughout rather than raw NONMEM $SIGMA /
    # $OMEGA variances.
    # ------------------------------------------------------------------
    propSd <- 0.35; label("Proportional residual error (fraction)")                       # Table 1 'Proportional error 0.35, RSE 9.1%, shrinkage 4.3%, 95% CI 0.30-0.44'
    addSd  <- 4.28; label("Additive residual error, fitted to below-LLOQ samples (ug/L)") # Table 1 'Additive error for BLQ samples (ug/L) 4.28, RSE 115.5%, shrinkage 4.3%'
  })

  model({
    # ---- 1. Individual parameters -------------------------------------
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)

    # Apparent milk volume with the centred linear pumping-duration
    # effect. 0.22 h (13.2 min) is the centring value carried in the
    # published covariate equation; at that pumping duration the factor
    # is 1 and vc equals the tabulated 19.9 L.
    vc <- exp(lvc) * (1 + e_t_pump_vc * (T_PUMP - 0.22))

    # ---- 2. Micro-constants (1/h) -------------------------------------
    kel <- cl / vc

    # ---- 3. ODE system ------------------------------------------------
    # Melander 2025 Results section 3 prints the milk compartment as
    #   dA(cet)/dt = ka x Dose - ke x A(cet)
    # which is the usual shorthand for a first-order input whose driving
    # term is the undepleted dose. The depot state below supplies that
    # input in its depleting form, ka x depot, so the pair reproduces the
    # standard one-compartment first-order-absorption solution. That is
    # the solution the paper itself used for its "calculated parameters
    # from model estimates": it returns Cmax = 30.56 ug/L against Table
    # 1's 30.6 ug/L. Amounts are in mg.
    d/dt(depot) <- -ka * depot
    d/dt(milk)  <-  ka * depot - kel * milk

    # ---- 4. Observation and residual error ----------------------------
    # milk is in mg and vc in L, so milk/vc is mg/L; the factor 1000
    # converts to the ug/L reporting scale of Table 1 and of the assay.
    Cmilk <- 1000 * milk / vc

    Cmilk ~ add(addSd) + prop(propSd)
  })
}
