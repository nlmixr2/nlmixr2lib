Kang_2026_ticagrelor <- function() {
  description <- paste(
    "Joint parent-metabolite one-compartment population PK model for ticagrelor and its",
    "active metabolite AR-C124910XX (TAM) in adults with acute coronary syndrome supported by",
    "veno-arterial extracorporeal membrane oxygenation (VA-ECMO) (Kang 2026). First-order",
    "absorption from a depot into the ticagrelor central compartment; ticagrelor leaves the",
    "central compartment by two parallel first-order routes sharing the same apparent clearance",
    "CL/F - a non-metabolic route (fraction 1 - fm) and a metabolic route (fraction fm) that",
    "forms TAM in its own one-compartment space with its own apparent clearance CLM/F and volume",
    "VM/F. Ka and fm were fixed for identifiability, and both additive residual error terms were",
    "fixed to preliminary estimates. Final-model covariates: a binary ECMO treatment-status",
    "indicator reducing CL/F 0.35-fold and expanding Vd/F 2.74-fold, plus an ECMO circuit",
    "blood-flow power term 0.518^(Q_ECMO/3) on Vd/F that applies only while on ECMO.",
    sep = " "
  )
  reference <- paste(
    "Kang S, Min KL, Yang S, Hahn J, Kim D, Jin BH, Chae SU, Bae SK, Wi J, Chang MJ.",
    "Population Pharmacokinetics of Ticagrelor during Veno-Arterial ECMO in Acute Coronary",
    "Syndrome: Model-Informed Dosing Simulations.",
    "Clin Pharmacol Ther. 2026;120(1):193-202. doi:10.1002/cpt.70282",
    sep = " "
  )
  vignette <- "Kang_2026_ticagrelor"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Amount units are mg because Appendix S2 sets S2 = Vd/F
  # and S3 = VM/F with no unit-conversion factor, so the state / volume ratio
  # is the modelled plasma concentration directly (see the units block).
  compartmentData <- list(
    depot       = list(analyte = "ticagrelor",             units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "ticagrelor",             units = "mg", specimen = "plasma",              verified = TRUE),
    central_tam = list(analyte = "AR-C124910XX (TAM)",     units = "mg", specimen = "plasma",              verified = TRUE)
  )

  covariateData <- list(
    ECMO_STATUS = list(
      description        = "Binary indicator: 1 = patient currently supported by VA-ECMO, 0 = patient off VA-ECMO (after weaning) (time-varying within subject)",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "Source column ECMO, coded ON-ECMO = 1 / OFF-ECMO = 0 (Kang 2026 Methods,",
        "'Population PK modeling', categorical-covariate list). ECMO modality: veno-arterial only;",
        "the cohort used a Capiox SP-101 centrifugal pump with a Capiox EBS X-coated conduit",
        "(Terumo) and a Sechrist air-oxygen mixer, femoral-femoral peripheral cannulation with",
        "17-Fr arterial and 21-Fr venous cannulas (Methods, 'ECMO system').",
        "Time-varying within subject: paired PK sampling was performed during ECMO (>= 24 h after",
        "initiation) and again after weaning, and Methods states 'all available data were included",
        "in the population PK analysis with ECMO status coded per occasion', so the per-record",
        "indicator switches at decannulation. 127 ON-ECMO observations from 19 patients and 98",
        "OFF-ECMO observations from 13 patients.",
        "Effect on CL/F is a power-form multiplier: CL/F = 12.2 * 0.35^ECMO_STATUS (a 65% clearance",
        "reduction on ECMO; 4.27 vs 12.2 L/h).",
        "Effect on Vd/F is a power-form multiplier 2.74^ECMO_STATUS, and it additionally gates the",
        "Q_ECMO flow term - see that entry.",
        "Mechanistic discussion (Kang 2026 Discussion): the CL/F reduction is attributed to",
        "end-organ hypoperfusion in cardiogenic shock reducing intrinsic hepatic metabolic capacity",
        "(ticagrelor is a low-extraction CYP3A4 substrate), plus possible circuit sequestration of a",
        "highly lipophilic, >99% protein-bound drug; the Vd/F expansion is attributed to",
        "hemodilution from resuscitation fluids, transfusions, albumin and circuit priming",
        "(all 20 patients were transfused and 18 received albumin) and to circuit sequestration.",
        sep = " "
      ),
      source_name        = "ECMO"
    ),
    Q_ECMO = list(
      description        = "Blood flow rate delivered through the VA-ECMO circuit",
      units              = "L/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column 'ECMO flow rate' (Kang 2026 Methods, 'Data collection and sample analysis'",
        "and the continuous-covariate list under 'Population PK modeling'; Table 2 row",
        "'theta ECMO flow rate on Vd/F').",
        "Enters Vd/F as a power term normalised to a 3 L/min reference:",
        "0.518^(Q_ECMO / 3). The 3 L/min denominator is printed in the final-model equation itself",
        "(Kang 2026 Results, 'Population PK modeling and evaluation') and is not a cohort median",
        "reported in Table 1 - Table 1 does not tabulate flow rate at all.",
        "Because 0.518 < 1, HIGHER flow gives a SMALLER Vd/F: at 3 L/min Vd/F is 222.8 L on ECMO,",
        "and reducing flow to 1 L/min raises it to 345.5 L.",
        "The flow term is nested inside the ECMO indicator in the NONMEM control stream",
        "(Appendix S2: V = THETA(3)*EXP(ETA(2))*THETA(10)**ECMO*(THETA(11)**(LPM/3))**ECMO), so it",
        "applies only while the patient is cannulated; off ECMO the circuit flow is zero and the",
        "term is 1 either way. See the model file comment at the vc assignment and the vignette",
        "Errata for the (numerically equivalent) main-text form.",
        "Time-varying: flow was deliberately titrated during support and stepped down during",
        "weaning - the paper's own Monte Carlo scenarios taper 3.4 L/min to 1 L/min in 0.2-0.3",
        "L/min steps every 6 h. Kang 2026 Discussion cautions that flow rate may be confounded by",
        "disease severity and that circuit extraction may be time- and flow-dependent, so the",
        "flow effect 'should be interpreted cautiously and considered exploratory'.",
        sep = " "
      ),
      source_name        = "ECMO flow rate"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 20L,
    n_studies      = 1L,
    age_range      = "36-88 years",
    age_median     = "59 years",
    weight_range   = "58-110 kg",
    weight_median  = "70.2 kg (Results); 70.8 kg ON-ECMO and 70.95 kg OFF-ECMO (Table 1)",
    sex_female_pct = 10,
    race_ethnicity = "Not reported (single-centre South Korean cohort at Severance Hospital, Seoul; presumed predominantly Korean)",
    disease_state  = paste(
      "Adults (> 19 years) with acute coronary syndrome receiving ticagrelor during and after",
      "VA-ECMO support in a coronary intensive care unit. 19 STEMI / 1 NSTEMI; 19 underwent PCI;",
      "15 of 20 had cardiac arrest. Median body mass index 25.1 kg/m^2 (20.7-35.9).",
      "Exclusions: pregnancy and strong CYP inducers (carbamazepine, phenytoin, rifampin,",
      "St. John's wort) or inhibitors (clarithromycin, itraconazole, ketoconazole,",
      "lopinavir/ritonavir).",
      sep = " "
    ),
    dose_range     = "Ticagrelor 180 mg oral loading dose followed by 90 mg twice daily maintenance (standard ACS regimen; the observed data come from this regimen only - the reduced 45-135 mg regimens in the paper are simulated, not observed).",
    regions        = "South Korea (single-centre prospective observational cohort, coronary intensive care unit, Severance Hospital, Yonsei University College of Medicine, Seoul; October 2015 - April 2018).",
    renal_function = "Median serum creatinine 1.30 mg/dL ON-ECMO (range 0.79-4.89) and 1.72 mg/dL OFF-ECMO (1.16-4.38); median BUN 20.8 vs 43.8 mg/dL. 5 patients received CRRT while ON-ECMO and 6 OFF-ECMO; CRRT was screened as a covariate and not retained. BUN, creatinine and uric acid were not evaluated for patients on CRRT.",
    hepatic_function = "Median total bilirubin 2.1 mg/dL ON-ECMO (0.9-6.3) and 1.4 mg/dL OFF-ECMO (0.4-24.2); median albumin 2.9 vs 2.95 g/dL; median total protein 4.9 vs 5.55 g/dL. Median AST 155 vs 63.5 IU/L.",
    ecmo_duration  = "Median 5.87 days (range 1.82-16.26)",
    notes          = paste(
      "Prospective observational cohort. 225 ticagrelor and 225 AR-C124910XX plasma",
      "concentrations (127 ON-ECMO, 98 OFF-ECMO) from 20 patients; paired sampling at pre-dose and",
      "1, 2, 3, 6, 8 and 12 h post-dose during ECMO (>= 24 h after initiation) and after weaning,",
      "drawn from an indwelling arterial line. Not all patients contributed both occasions.",
      "Validated LC-MS/MS quantification. NONMEM 7.4 with ADVAN6 (TOL = 3) and FOCE-I, Pirana",
      "2.9.7, PsN and Xpose4 in R 4.4.0. Model evaluation by goodness-of-fit plots, nonparametric",
      "bootstrap (n = 5,000) and prediction-corrected VPC (n = 1,000) stratified by ECMO status.",
      "Covariates screened but NOT retained: height, weight, time to start ECMO, duration of ECMO,",
      "ECMO pump speed, BUN, serum creatinine, uric acid, total protein, albumin, total bilirubin,",
      "sex, smoking and CRRT.",
      sep = " "
    )
  )

  ini({
    # ---- Structural fixed effects (Kang 2026 Table 2, 'Population estimations' column) ----
    # CL/F and Vd/F are the OFF-ECMO reference values; the ECMO multipliers below
    # move them to their ON-ECMO values.
    lcl     <- log(12.2); label("Apparent ticagrelor clearance CL/F, OFF-ECMO (L/h)")                        # Kang 2026 Table 2: theta CL/F = 12.2 (RSE 16.7%; bootstrap median 12.01, 95% CI 8.99-16.24)
    lvc     <- log(157);  label("Apparent ticagrelor volume of distribution Vd/F, OFF-ECMO (L)")             # Kang 2026 Table 2: theta Vd/F = 157 (RSE 17.4%; bootstrap median 152.93, 95% CI 101.04-291.14)
    lcl_tam <- log(8.8);  label("Apparent AR-C124910XX clearance CLM/F (L/h)")                               # Kang 2026 Table 2: theta CLM/F = 8.8 (RSE 10.7%; bootstrap median 8.83, 95% CI 7.27-10.76)
    lvc_tam <- log(29.5); label("Apparent AR-C124910XX volume of distribution VM/F (L)")                     # Kang 2026 Table 2: theta VM/F = 29.5 (RSE 42.7%; bootstrap median 28.86, 95% CI 12.63-102.22)

    # Ka and fm were FIXED, not estimated: Kang 2026 Methods ('Population PK
    # modeling') states both were fixed "due to identifiability limitations in
    # this VA-ECMO dataset, including limited absorption-phase sampling", and
    # Table 2 prints them as "0.533 FIX" and "0.224 FIX" with no RSE, no
    # bootstrap median and no bootstrap CI. Ka came from a preliminary base-model
    # estimation on these data; fm was taken from published non-ECMO ticagrelor
    # popPK models (Kang 2026 refs 23-24, Li 2016 and Roshammar 2017).
    # Appendix S3 reports the fixed-parameter sensitivity analysis: the parent
    # estimates are stable across plausible Ka (12.2-12.8 L/h for CL/F over
    # Ka 0.533-2) and completely invariant to fm (CL/F 12.2 and Vd/F 157 at every
    # fm from 0.1 to 0.9, identical OFV -1142.42), while CLM/F and VM/F scale
    # in proportion to fm with an essentially unchanged CLM/VM ratio.
    lka <- fixed(log(0.533)); label("First-order absorption rate constant Ka (1/h, : held at a preliminary base-model estimate for identifiability)")     # Kang 2026 Table 2: theta Ka = 0.533 FIX
    fm  <- fixed(0.224);      label("Fraction of ticagrelor converted to AR-C124910XX (unitless, : held at a published non-ECMO value for identifiability)")  # Kang 2026 Table 2: theta FMET = 0.224 FIX

    # ---- Covariate effects (Kang 2026 Table 2 and the final-model equations in
    # Results, 'Population PK modeling and evaluation') ----
    # Stored on the paper's own printed scale so that the model() expressions
    # reproduce the published equations verbatim:
    #   CL/F = 12.2 * 0.35^ECMO
    #   Vd/F = 157  * 0.518^(flow/3) * 2.74^ECMO
    # All three are power-form multipliers raised to an indicator or a scaled
    # continuous covariate, so a value below 1 is a decrement and above 1 an
    # increment.
    e_ecmo_cl   <- 0.35;  label("Power-form ECMO_STATUS multiplier on CL/F (unitless)")                      # Kang 2026 Table 2: theta ECMO on CL/F = 0.35 (RSE 6.5%; bootstrap median 0.356, 95% CI 0.302-0.433)
    e_ecmo_vc   <- 2.74;  label("Power-form ECMO_STATUS multiplier on Vd/F (unitless)")                      # Kang 2026 Table 2: theta ECMO on Vd/F = 2.74 (RSE 34.5%; bootstrap median 2.82, 95% CI 1.28-5.30)
    e_qecmo_vc <- 0.518; label("Power-form multiplier on Vd/F for (Q_ECMO / 3 L/min) (unitless)")           # Kang 2026 Table 2: theta ECMO flow rate on Vd/F = 0.518 (RSE 32%; bootstrap median 0.526, 95% CI 0.248-0.839)

    # ---- Inter-individual variability ----
    # Exponential (log-normal) IIV on CL/F, Vd/F and CLM/F per Methods
    # ('Interindividual variability (IIV) was modeled exponentially and included
    # on CL/F, Vd/F, and CLM/F') and Appendix S2 ('IIV was modeled using
    # exponential random effects'). Table 2 prints omega^2 directly, so the
    # values below are used as-is with no CV-to-variance conversion. The scale is
    # confirmed by Table 2 footnote b, %CV = sqrt(exp(omega^2) - 1) * 100, which
    # reproduces every printed %CV from the printed omega^2:
    #   sqrt(exp(0.441) - 1) = 74.4%, sqrt(exp(0.661) - 1) = 96.8%,
    #   sqrt(exp(0.156) - 1) = 41.1%.
    # No off-diagonal covariances were reported, so the etas are independent.
    etalcl     ~ 0.441  # Kang 2026 Table 2, row 'omega^2 CL/F'   = 0.441 (RSE 34.7%) [74.4% CV]; bootstrap median 0.420, 95% CI 0.188-0.684
    etalvc     ~ 0.661  # Kang 2026 Table 2, row 'omega^2 Vd/F'   = 0.661 (RSE 34%)   [96.8% CV]; bootstrap median 0.624, 95% CI 0.310-1.06
    etalcl_tam ~ 0.156  # Kang 2026 Table 2, row 'omega^2 CLM/F'  = 0.156 (RSE 47.74%) [41.1% CV]; bootstrap median 0.144, 95% CI 0.051-0.315

    # ---- Residual unexplained variability ----
    # Additive on the linear concentration scale for BOTH analytes. Kang 2026
    # Methods: additive, proportional and combined error models were compared,
    # "the proportional component was negligible and did not meaningfully improve
    # model diagnostics; therefore, an additive error model was selected for both
    # parent and metabolite", and "residual error parameters were fixed to the
    # preliminary estimates, with residual variance fixed at 0.1" - hence fixed()
    # on both lines and the 'FIX' flag in Table 2.
    #
    # SCALE. Appendix S2 writes the error as Y = IPRED + W*EPS(1) with
    # W = THETA(7) for CMT 2 (parent) and THETA(8) for CMT 3 (metabolite), and
    # Var(EPS(1)) fixed at 0.1. Table 2 footnote c states the same in words:
    # "Residual variance was fixed at 0.1; therefore, residual SD = sqrt(0.1) x W".
    # So the additive residual SD is W * sqrt(0.1), i.e. 1.65 * sqrt(0.1) = 0.5217
    # and 0.2 * sqrt(0.1) = 0.06325.
    #
    # UNITS. These SDs are in the model's own concentration units, mg/L, not the
    # ng/mL used for reporting in the paper's figures and target range. Appendix
    # S2 sets S2 = Vd/F and S3 = VM/F with no 1000-fold conversion factor, so with
    # doses in mg and volumes in L the modelled concentration is in mg/L. Three
    # independent checks agree: (i) 0.5217 and 0.06325 are ~30% and ~33% of the
    # respective typical steady-state concentrations on the 90 mg q12h observed
    # regimen in mg/L (1.76 and 0.191 mg/L), a plausible residual magnitude,
    # whereas on a ng/mL scale they would be ~0.03% and would be
    # arithmetically incompatible with the reported OFV of -1142.42;
    # (ii) their ratio 8.25 tracks the parent/metabolite concentration ratio 9.2;
    # (iii) 90 mg / 157 L = 0.573 mg/L reproduces the reported ticagrelor
    # exposure range. The vignette converts to ng/mL for comparison with the
    # paper's figures.
    addSd     <- fixed(1.65 * sqrt(0.1)); label("Additive residual error for ticagrelor (mg/L, : from the preliminary model)")          # Kang 2026 Table 2: sigma ADD  = 1.65 FIX; footnote c, residual SD = sqrt(0.1) x W
    addSd_tam <- fixed(0.2  * sqrt(0.1)); label("Additive residual error for AR-C124910XX (mg/L, : from the preliminary model)")        # Kang 2026 Table 2: sigma ADDM = 0.2 FIX;  footnote c, residual SD = sqrt(0.1) x W
  })

  model({
    # ---- 1. Individual PK parameters ----
    ka <- exp(lka)

    # CL/F = 12.2 * 0.35^ECMO (Kang 2026 Results final-model equation;
    # Appendix S2 $PK: CL = THETA(2)*EXP(ETA(1))*THETA(9)**ECMO).
    cl <- exp(lcl + etalcl) * e_ecmo_cl^ECMO_STATUS

    # Vd/F. Appendix S2 $PK writes
    #   V = THETA(3)*EXP(ETA(2)) * THETA(10)**ECMO * (THETA(11)**(LPM/3))**ECMO
    # i.e. the ECMO-flow power term is itself raised to the ECMO indicator, so it
    # applies only while the patient is cannulated. That control-stream form is
    # encoded here because it is the code that was actually run and because it is
    # robust to however Q_ECMO is carried on OFF-ECMO records. The main-text
    # equation prints the flow term ungated,
    #   Vd/F = 157 * 0.518^(flow/3) * 2.74^ECMO,
    # which is algebraically identical whenever the delivered circuit flow is 0
    # off ECMO - the physically correct coding, and the one this model expects.
    # Both forms reproduce every ON-ECMO number the paper reports:
    #   157 * 2.74 * 0.518^(3.0/3) = 222.8 L at 3 L/min on ECMO,
    #   157 * 2.74 * 0.518^(1.0/3) = 345.5 L at 1 L/min on ECMO.
    # The paper's "81.3 L OFF-ECMO" is the same 222.8 L divided by the 2.74-fold
    # ECMO effect (222.8 / 2.74 = 81.3), i.e. a matched-flow counterfactual that
    # isolates the volume expansion; see the vignette Errata.
    vc <- exp(lvc + etalvc) * e_ecmo_vc^ECMO_STATUS *
          (e_qecmo_vc^(Q_ECMO / 3))^ECMO_STATUS

    cl_tam <- exp(lcl_tam + etalcl_tam)
    vc_tam <- exp(lvc_tam)

    # ---- 2. Micro-constants (Kang 2026 Figure 1 caption and Appendix S1) ----
    # Ticagrelor leaves `central` by two parallel first-order routes that share
    # the same apparent clearance CL/F: k20 is the non-metabolic route and k23
    # the metabolic route forming TAM. Their sum is CL/F / (Vd/F), so total
    # apparent parent clearance is CL/F and AUCinf(parent) = Dose / (CL/F).
    # k30 is the TAM elimination rate constant.
    k20 <- (1 - fm) * cl / vc
    k23 <- fm       * cl / vc
    k30 <- cl_tam / vc_tam

    # ---- 3. ODE system (Appendix S1 / Appendix S2 $DES) ----
    # DADT(1) = -KA*A(1)
    # DADT(2) =  KA*A(1) - K20*A(2) - K23*A(2)
    # DADT(3) =  K23*A(2) - K30*A(3)
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - k20 * central - k23 * central
    d/dt(central_tam) <-  k23 * central - k30 * central_tam

    # ---- 4. Observations ----
    # Appendix S1: IPRED_parent = A(2)/(Vd/F) = A(2)/S2 and
    #              IPRED_metab  = A(3)/(VM/F) = A(3)/S3.
    # Bioavailability F is not estimated, so no f(depot) term: the whole dose
    # enters `depot` and F is absorbed into the apparent CL/F and Vd/F.
    # Doses in mg and volumes in L give concentrations in mg/L.
    Cc     <- central     / vc
    Cc_tam <- central_tam / vc_tam

    Cc     ~ add(addSd)
    Cc_tam ~ add(addSd_tam)
  })
}
