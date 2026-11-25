Joshi_2017_vancomycin_rat <- function() {
  description <- paste(
    "Preclinical (rat).",
    "Two-compartment population PK model of vancomycin in male Sprague-Dawley",
    "rats after intraperitoneal (i.p.) injection, with first-order absorption",
    "from the peritoneal cavity (depot), an estimated fraction absorbed, and",
    "intercompartmental transfer parameterised directly as the rate constants",
    "Kcp and Kpc. Fitted non-parametrically with the NPAG algorithm in Pmetrics",
    "to richly sampled plasma concentrations from 10 rats given 150 mg/kg once",
    "daily for 1 or 3 days. The paper prints the means and SDs of the NPAG",
    "marginal distributions without units; volume is in litres for the studied",
    "250-350 g rats and the dose is the absolute amount in mg (see the vignette",
    "for how the units were settled). Each parameter carries an independent",
    "inter-individual variability derived from its marginal SD; residual",
    "variability is fixed(0) because the Pmetrics error model was not reported.",
    sep = " "
  )
  reference <- paste(
    "Joshi MD, O'Donnell JN, Venkatesan N, Chang J, Nguyen H, Rhodes NJ,",
    "Pais G, Chapman RL, Griffin B, Scheetz MH. High-Performance Liquid",
    "Chromatography Method for Rich Pharmacokinetic Sampling Schemes in",
    "Translational Rat Toxicity Models With Vancomycin.",
    "Clin Transl Sci. 2017;10(6):496. doi:10.1111/cts.12484.",
    "PMID 28675684; PMCID PMC5698807.",
    "All parameter values are Table 3 ('Pharmacokinetic parameters from the",
    "two-compartment model based on plasma sampling'); no supplement was",
    "deposited. The unit convention follows the same group's predecessor rat",
    "vancomycin model, Rhodes NJ et al. Antimicrob Agents Chemother.",
    "2016;60(10):5742-5751 (doi:10.1128/AAC.00591-16), Table 2, which gives",
    "volume in 'liter/0.3 kg'.",
    sep = " "
  )
  vignette <- "Joshi_2017_vancomycin_rat"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # No covariates: Table 3 carries no covariate term and the Methods describe
  # none. Body weight enters only when converting the mg/kg dose into the
  # absolute mg amount the model is dosed with.
  covariateData <- list()

  compartmentData <- list(
    depot = list(analyte = "vancomycin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "rat (Sprague-Dawley)",
    n_subjects = 10L,
    n_studies = 1L,
    age_range = "not reported",
    weight_range = "250-350 g",
    sex_female_pct = 0,
    disease_state = "healthy (translational vancomycin nephrotoxicity model)",
    dose_range = paste(
      "150 mg/kg/day vancomycin (clinical grade, Hospira) as a single daily",
      "intraperitoneal injection for 1 day (n = 5) or 3 days (n = 5)"
    ),
    regions = "United States (Midwestern University, Downers Grove, IL)",
    notes = paste(
      "Methods, 'Pharmacokinetic study': ten male Sprague-Dawley rats were",
      "selected from a larger study (Rhodes 2016) to illustrate the HPLC assay.",
      "Blood (0.25 mL) was drawn from a jugular-vein catheter at 0 (predose),",
      "10, 30, 60, 120, 240 and 1,440 min after the FINAL dose. Plasma was",
      "assayed by HPLC-UV (LLOQ 3 mg/L, linear to 75 mg/L). Estimation was",
      "NONPARAMETRIC (NPAG in Pmetrics); Table 3 summarises the marginal",
      "support-point distributions as mean and SD, so the joint density and",
      "every parameter correlation are not recoverable."
    )
  )

  ini({
    # =====================================================================
    # STRUCTURAL PARAMETERS -- Joshi 2017 Table 3, 'Mean' column.
    #
    # Each value is the MEAN of the NPAG marginal support-point
    # distribution, and the 'SD' column is the spread of that distribution
    # across rats (inter-individual variability, not uncertainty). Following
    # the library's Pmetrics precedents (Sime_2019_tazobactam.R,
    # Hughes_2024_vancomycin_nonparametric.R), each mean is carried as the
    # MEDIAN of a log-normal marginal.
    #
    # UNITS are not printed in Table 3. Time is hours (Results report the
    # half-life as 1.7 h, and the typical-value terminal half-life computed
    # from these rate constants is 1.58 h). Volume is litres for the
    # studied 250-350 g rats, and the model is dosed with the absolute
    # amount in mg (150 mg/kg x body weight): the same group's Rhodes 2016
    # rat model reports V as 'liter/0.3 kg', and a simulated cohort under
    # this reading lands within ~30% of the Results' median Cmax and AUC
    # while the alternative (L/kg with a mg/kg dose) overshoots by a factor
    # of 2-2.5. See the vignette.
    # =====================================================================
    lkel <- log(0.911843)
    label("Elimination rate constant from the central compartment, Ke (1/h)")
    # Table 3, Ke: mean 0.911843, SD 0.785514
    lvc <- log(0.185885)
    label("Central volume of distribution, V, for a 250-350 g rat (L)")
    # Table 3, V: mean 0.185885, SD 0.14991
    lka <- log(3.946906)
    label("First-order absorption rate constant from the peritoneal cavity, Ka (1/h)")
    # Table 3, Ka: mean 3.946906, SD 2.009931
    lk12 <- log(45.33389)
    label("Transfer rate constant central -> peripheral1, KCP (1/h)")
    # Table 3, KCP: mean 45.33389, SD 26.78453
    lk21 <- log(42.25276)
    label("Transfer rate constant peripheral1 -> central, KPC (1/h)")
    # Table 3, KPC: mean 42.25276, SD 31.89025
    logitfdepot <- log(0.459455 / (1 - 0.459455))
    label("Logit of the fraction of the i.p. dose absorbed, FA1 (fraction)")
    # Table 3, FA1: mean 0.459455, SD 0.156658. Held on the logit scale so
    # simulated fractions stay inside (0, 1); logit(0.459455) = -0.162537.

    # =====================================================================
    # INTER-INDIVIDUAL VARIABILITY -- Joshi 2017 Table 3, 'SD' column.
    #
    # CV = SD / mean, converted to a log-normal variance exactly as
    # omega^2 = log(CV^2 + 1). For FA1 (logit scale) the SD is carried
    # through the delta method: var(logit F) = (SD / (F * (1 - F)))^2.
    # The marginals are encoded as INDEPENDENT because the NPAG joint
    # distribution is not published.
    # =====================================================================
    etalkel ~ 0.555096 # Table 3 Ke: CV 0.785514/0.911843 = 86.1%; log(CV^2 + 1)
    etalvc ~ 0.501010 # Table 3 V: CV 0.14991/0.185885 = 80.6%; log(CV^2 + 1)
    etalka ~ 0.230578 # Table 3 Ka: CV 2.009931/3.946906 = 50.9%; log(CV^2 + 1)
    etalk12 ~ 0.299421 # Table 3 KCP: CV 26.78453/45.33389 = 59.1%; log(CV^2 + 1)
    etalk21 ~ 0.450851 # Table 3 KPC: CV 31.89025/42.25276 = 75.5%; log(CV^2 + 1)
    etalogitfdepot ~ 0.397883 # Table 3 FA1: (0.156658 / (0.459455 * 0.540545))^2

    # =====================================================================
    # RESIDUAL UNEXPLAINED VARIABILITY is not reported. The Methods name
    # NPAG in Pmetrics but give neither the error model (gamma or lambda)
    # nor an assay-error polynomial. The Results quote only the Bayesian
    # posterior fit statistics (R^2 0.98, bias -0.25 mg/L, imprecision
    # 1.05 mg^2/L^2), which describe individual predictions, not a
    # population residual. Carried as fixed(0) rather than invented.
    # =====================================================================
    propSd <- fixed(0)
    label("Proportional residual SD (fraction; not reported in the source)")
    addSd <- fixed(0)
    label("Additive residual SD (mg/L; not reported in the source)")
  })

  model({
    # Individual primary parameters (Table 3 parameterisation).
    kel <- exp(lkel + etalkel)
    vc <- exp(lvc + etalvc)
    ka <- exp(lka + etalka)
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)
    fdepot <- expit(logitfdepot + etalogitfdepot)

    # Clearance, intercompartmental clearance and peripheral volume derived
    # from the rate constants, so that rxode2's linear-compartment
    # recognition (keyed on cl / vc / q / vp / ka) and the explicit ODEs
    # below describe the same system.
    cl <- kel * vc
    q <- k12 * vc
    vp <- q / k21

    # Methods, 'Pharmacokinetic model': 'An absorption term (i.e., Ka), a
    # bioavailability term (Fa), and a peritoneal compartment were included
    # to account for intraperitoneal absorption.' The peritoneal cavity is
    # the depot; FA1 is the fraction of the i.p. dose that reaches plasma.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    f(depot) <- fdepot

    # Dose in mg, volume in L -> mg/L (the paper's 'mg/L' axes, Figures 3-4).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
