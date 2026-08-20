Centanni_2025_sunitinib_thrombocytopenia <- function() {
  description <- "Semi-physiological Friberg-Karlsson myelosuppression model for the absolute thrombocyte (platelet) count during sunitinib treatment in adults with imatinib-resistant gastrointestinal stromal tumours (GIST). A self-renewing proliferating progenitor pool feeds five transit compartments reflecting megakaryocyte maturation and a circulating-platelet pool, with an Emax drug-effect function driven directly by the sunitinib daily AUC inhibiting proliferation and a (Circ0/circ)^gamma feedback term. This is the one model newly developed in Centanni 2025; it was added to the previously published Hansson 2013 sunitinib/GIST pharmacometric framework, whose other components are packaged separately (see reference). The model has no PK ODE: it consumes the current daily dose (DOSE) and the individual posthoc upstream popPK clearance (CLI) as data covariates and forms the exposure driver AUC = DOSE / CLI. NOTE: three parameters in the source Table S2 are encoded at 1/1000 of the printed value; see the vignette Errata and the inline comments in ini() for the full justification."
  reference <- paste(
    "Centanni M, Nijhuis J, Karlsson MO, Friberg LE.",
    "Comparative Analysis of Traditional and Pharmacometric-Based",
    "Pharmacoeconomic Modeling in the Cost-Utility Evaluation of",
    "Sunitinib Therapy.",
    "PharmacoEconomics. 2025;43(1):31-43.",
    "doi:10.1007/s40273-024-01438-z.",
    "The thrombocyte model parameters are in Online Resource 1, Table S2;",
    "its structure is given there only by reference to the ANC model.",
    "Parent pharmacometric framework (not re-estimated by Centanni 2025):",
    "Hansson EK, Ma G, Amantea MA, French J, Milligan PA, Friberg LE,",
    "Karlsson MO. PKPD modeling of predictors for adverse effects and",
    "overall survival in sunitinib-treated patients with GIST.",
    "CPT Pharmacometrics Syst Pharmacol. 2013;2(11):e85.",
    "doi:10.1038/psp.2013.62; see",
    "modellib('Hansson_2013_sunitinib_myelosuppression') [ANC],",
    "modellib('Hansson_2013_sunitinib_dbp'),",
    "modellib('Hansson_2013_sunitinib_hfs'),",
    "modellib('Hansson_2013_sunitinib_os').",
    "Friberg myelosuppression backbone:",
    "Friberg LE et al. J Clin Oncol 2002;20(24):4713-4721,",
    "doi:10.1200/JCO.2002.02.140.",
    sep = " "
  )
  vignette <- "Centanni_2025_sunitinib_thrombocytopenia"
  units <- list(time = "h", dosing = "mg", concentration = "10^9 cells/L (platelets)")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Matches the sibling
  # Hansson_2013_sunitinib_myelosuppression compartmentData shape, with the
  # transit chain extended from three to five stages per Online Resource 1.
  # verified = FALSE means NOT checked against the source paper (the source
  # gives the chain only by reference to the ANC model).
  compartmentData <- list(
    prol     = list(analyte = "megakaryocyte progenitor cells", units = "10^9 cells/L", specimen = "blood cell", verified = FALSE),
    transit1 = list(analyte = "maturing megakaryocyte precursors (stage 1)", units = "10^9 cells/L", specimen = "blood cell", verified = FALSE),
    transit2 = list(analyte = "maturing megakaryocyte precursors (stage 2)", units = "10^9 cells/L", specimen = "blood cell", verified = FALSE),
    transit3 = list(analyte = "maturing megakaryocyte precursors (stage 3)", units = "10^9 cells/L", specimen = "blood cell", verified = FALSE),
    transit4 = list(analyte = "maturing megakaryocyte precursors (stage 4)", units = "10^9 cells/L", specimen = "blood cell", verified = FALSE),
    transit5 = list(analyte = "maturing megakaryocyte precursors (stage 5)", units = "10^9 cells/L", specimen = "blood cell", verified = FALSE),
    circ     = list(analyte = "absolute thrombocyte count", units = "10^9 cells/L", specimen = "whole blood", verified = FALSE)
  )

  covariateData <- list(
    DOSE = list(
      description        = "Current administered sunitinib daily dose (mg) carried as a time-varying data column. Set to 0 during a dose interruption or for untreated subjects so the derived AUC = DOSE / CLI becomes 0 and the drug effect vanishes.",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Centanni 2025 simulated a continuous 37.5 mg once-daily sunitinib regimen (Methods 2.1), with protocol dose reductions to 25, 12.5 or 0 mg on unacceptable adverse events (Methods 2.1, citing the prescribing information). For typical-cohort vignette simulations the value is held at 37.5 mg, matching the paper's base-case comparator arm.",
      source_name        = "DOSE"
    ),
    CLI = list(
      description        = "Individual posthoc total plasma clearance (L/h) of sunitinib from the upstream popPK fit. Per-subject, time-fixed. Used only to form the exposure driver AUC = DOSE / CLI.",
      units              = "L/h",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required input. Centanni 2025 did not re-estimate the PK; it reuses the Hansson 2013 framework, whose upstream popPK is Houk et al. 2009 Clin Cancer Res 15:2497-2506 (not packaged in nlmixr2lib at extraction time). The sibling Hansson_2013_sunitinib_myelosuppression, Hansson_2013a_sunitinib and Hansson_2013c_sunitinib model files use a typical-value reference of 32.819 L/h; the same value is used in this model's vignette so the exposure driver is consistent across the framework.",
      source_name        = "CL"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1000L,
    n_studies      = 1L,
    age_range      = "adults with GIST; Centanni 2025 generated virtual patients from the covariate distributions in Methods section 2.1 rather than reporting an observed baseline-demographics table",
    weight_range   = "not reported in the source",
    sex_female_pct = NA_real_,
    race_ethnicity = "not reported in the source",
    disease_state  = "imatinib-resistant gastrointestinal stromal tumours (GIST)",
    dose_range     = "sunitinib 37.5 mg orally once daily, continuous dosing, with protocol dose reductions to 25, 12.5 or 0 mg on unacceptable adverse events; comparator arm received no treatment",
    regions        = "not reported in the source",
    biomarkers     = "absolute thrombocyte (platelet) count, evaluated every 6 weeks in the simulation per established clinical protocols (Methods 2.1). Unacceptable thrombocytopenia was defined as a platelet count < 50 x 10^9/L, aggregated per 6-week cycle.",
    notes          = "IMPORTANT PROVENANCE NOTE. Centanni 2025 is a health-economics methods-comparison study; the N=1,000 above is the size of the simulated two-arm virtual trial that the pharmacoeconomic comparison ran on, NOT an observed cohort for this PD model. The thrombocyte model itself is described in Online Resource 1 as 'developed' and is the one component of the framework that Centanni 2025 newly estimated (main text Methods 2.2: 'The pharmacometric-based models were not re-estimated, except for the thrombocytopenia model, which was newly developed'). The source does not state which patient dataset the thrombocyte model was estimated on; by construction of the framework it is the pooled four-study imatinib-resistant GIST sunitinib dataset (n = 303) underlying Hansson 2013. Treat the estimation-population description as inherited rather than reported."
  )

  ini({
    # ------------------------------------------------------------------
    # All structural parameters are from Centanni 2025 Online Resource 1,
    # Table S2 ("Absolute Thrombocyte Count Model"). The source gives the
    # model structure only by reference: "The structural model is similar
    # to the published ANC model, with the exception that five transit
    # compartment were used instead of three. Additionally, sunitinib area
    # under the curve (AUC) was the best predictor of thrombocyte changes.
    # An Emax function described the relationship between AUC and
    # thrombocytes." No equation is printed for this model.
    #
    # SUSPECTED UNIT / SCALE ERRATUM IN THE SOURCE TABLE S2 -- READ THIS.
    # Three of the five structural parameters are encoded here at 1/1000 of
    # the value printed in Table S2. The printed values cannot be solved and
    # cannot reproduce the paper's own reported thrombocytopenia incidence;
    # the corrected values can. The three corrections are the SAME factor of
    # 1000 in the SAME direction within one table, which is consistent with a
    # single systematic scale/unit error in the typesetting of Table S2
    # rather than with three independent estimation results. Each corrected
    # line is annotated individually below, and the full falsification
    # evidence is tabulated in the vignette's "Assumptions and deviations"
    # section. Encoded per operator ruling on task sidecar
    # oare_PMC11724784 request-001 (2026-08-14).
    # ------------------------------------------------------------------

    # Baseline circulating platelet count. Printed value used as-is.
    lcirc0 <- log(323);    label("Baseline circulating thrombocyte count Circ0 (10^9 cells/L)")   # Table S2 'Baseline thrombocyte count (x 10^9/liter)' = 323 (RSE 2.3%)

    # Mean transit time through the proliferation -> 5-transit chain.
    # Printed value used as-is.
    lmtt   <- log(175);    label("Mean transit time MTT through the progenitor -> transit maturation chain (h)")  # Table S2 'MTT (hours)' = 175 (RSE 4.2%)

    # ERRATUM 1 of 3. Table S2 prints Emax = 497. In the Friberg
    # proliferation term ktr*prol*(1 - Emax*f(AUC))*(Circ0/circ)^gamma the
    # drug effect Emax*f(AUC) must stay below 1 or the proliferation rate
    # becomes negative and the system is unsolvable. At the paper's own
    # 37.5 mg/day with the framework clearance of 32.819 L/h the Emax
    # function evaluates to a drug effect of 492.6 with the printed Emax --
    # far outside the solvable region. Encoded at 0.497, which additionally
    # sits essentially on top of the sibling ANC model's Emax of 0.520
    # (Hansson 2013 Table 2) in the same framework by the same authors.
    lemax  <- log(0.497);  label("Maximum fractional inhibition Emax of megakaryocyte-progenitor proliferation by sunitinib AUC (unitless)")  # Table S2 'Emax' printed 497 (RSE 33.2%); encoded 0.497 -- suspected 1000x scale erratum, see block comment above

    # ERRATUM 2 of 3. Table S2 prints EC50 = 10300 with the unit
    # "ng * hour / L", i.e. 0.0103 mg*h/L. The framework's exposure driver
    # is the sunitinib daily AUC = DOSE / CLI, which at 37.5 mg/day and
    # 32.819 L/h is 1.14 mg*h/L -- 111-fold above an EC50 of 0.0103 mg*h/L,
    # so the Emax function would sit 99.1% saturated. A saturated EC50 also
    # makes the separately reported EC50 IIV of 66.9% CV inoperative, which
    # is internal evidence against the printed unit. Reading the same
    # printed number 10300 as ng*hour/mL gives 10.3 mg*h/L, which places the
    # 37.5 mg/day exposure at 10% of the Emax curve -- a position where a
    # 67% CV on EC50 genuinely matters. Encoded as 10.3 mg*h/L.
    lec50  <- log(10.3);   label("Sunitinib daily AUC producing half-maximal inhibition of platelet proliferation EC50 (mg*h/L)")  # Table S2 'EC50 (ng * hour / L)' printed 10300 (RSE 34.6%); encoded 10.3 mg*h/L (= 10300 ng*h/mL) -- suspected 1000x unit erratum, see block comment above

    # ERRATUM 3 of 3. Table S2 prints gamma = 92.3. gamma is the feedback
    # exponent in (Circ0/circ)^gamma; the feedback sensitivity is
    # proportional to gamma, so gamma = 92.3 is an extremely strong feedback
    # gain that drives the 7-stage chain far outside its stable region (the
    # solution diverges within the paper's own 104-week horizon). It is also
    # self-defeating at steady state: the analytic steady state of this
    # chain is circ_ss/Circ0 = (1 - edrug)^(1/gamma), so gamma = 92.3 moves
    # platelets by only 0.06% under any attainable drug effect, contradicting
    # the thrombocytopenia events the paper reports. Encoded at 0.0923, which
    # is stable and is the same order as the sibling ANC model's gamma of
    # 0.362 (Hansson 2013 Table 2).
    lgamma <- log(0.0923); label("Feedback exponent gamma on (Circ0 / circ) driving compensatory progenitor proliferation (unitless)")  # Table S2 'gamma' printed 92.3 (RSE 8.6%); encoded 0.0923 -- suspected 1000x scale erratum, see block comment above

    # ------------------------------------------------------------------
    # Inter-individual variability. Table S2 reports IIV as CV%, with the
    # footnote "calculated as the sqrt(exp(omega^2) - 1)", so the
    # back-transform to the nlmixr2 variance scale is
    # omega^2 = log((CV/100)^2 + 1).
    #   Circ0 CV = 33.3% -> omega2 = log(0.333^2 + 1) = 0.105161
    #   MTT   CV = 17.7% -> omega2 = log(0.177^2 + 1) = 0.030848
    #   EC50  CV = 66.9% -> omega2 = log(0.669^2 + 1) = 0.369880
    # Emax and gamma are reported with "-" in the IIV column, i.e. no IIV
    # was estimated for them; they carry no eta here.
    # Table S2 reports no correlation between the random effects, so the
    # omega matrix is diagonal (contrast the sibling ANC model, where
    # Hansson 2013 reported a 90% ANC0-Emax correlation).
    # ------------------------------------------------------------------
    etalcirc0 ~ 0.105161  # Table S2 Circ0 IIV CV = 33.3% (RSE 5%)
    etalmtt   ~ 0.030848  # Table S2 MTT IIV CV = 17.7% (RSE 9.4%)
    etalec50  ~ 0.369880  # Table S2 EC50 IIV CV = 66.9% (RSE 11.3%)

    # ------------------------------------------------------------------
    # Residual error. Table S2 row 'Residual error (x 10^9/liter)' = 0.268
    # (RSE 6%). The printed value and its printed unit are encoded here
    # literally, as an additive residual SD of 0.268 x 10^9/L. Unlike the
    # three parameters above, this reading is solvable and is therefore NOT
    # treated as an erratum -- the standing rule is to deviate from a printed
    # value only where the printed value cannot reproduce the source. Note
    # for downstream users that 0.268 x 10^9/L against a 323 x 10^9/L
    # baseline is a residual CV of only 0.08%, which is very tight for a
    # population PD residual; an alternative reading (a proportional /
    # transformed-scale residual of 26.8%) is discussed in the vignette
    # Assumptions and deviations section, where both readings are scored
    # against the paper's reported thrombocytopenia incidence.
    # ------------------------------------------------------------------
    addSd <- 0.268; label("Additive residual error SD on the absolute thrombocyte count (10^9 cells/L)")  # Table S2 'Residual error (x 10^9/liter)' = 0.268 (RSE 6%)
  })

  model({
    # ---- 1. Exposure driver: sunitinib daily AUC (mg*h/L) ----
    # mg / (L/h) = mg*h/L. Online Resource 1: "sunitinib area under the
    # curve (AUC) was the best predictor of thrombocyte changes". This is
    # the same DOSE / CLI daily-AUC convention the sibling Hansson 2013
    # sunitinib framework files use.
    auc <- DOSE / CLI

    # ---- 2. Individual structural parameters ----
    circ0 <- exp(lcirc0 + etalcirc0)
    mtt   <- exp(lmtt + etalmtt)
    emax  <- exp(lemax)
    ec50  <- exp(lec50 + etalec50)
    gamma <- exp(lgamma)

    # ---- 3. Emax drug effect on proliferation, and the feedback term ----
    # Online Resource 1: "An Emax function described the relationship
    # between AUC and thrombocytes."
    edrug <- emax * auc / (ec50 + auc)
    feed  <- (circ0 / circ)^gamma

    # ---- 4. Friberg myelosuppression chain with FIVE transit compartments ----
    # ktr = (n_transit + 1) / MTT with n_transit = 5 per Online Resource 1
    # ("five transit compartment were used instead of three"), following the
    # Friberg 2002 / Wahlby 2004 convention also used by the sibling
    # Hansson_2013_sunitinib_myelosuppression file (which uses 4 / mtt for
    # its three-transit chain).
    ktr <- 6 / mtt

    d/dt(prol)     <- ktr * prol * (1 - edrug) * feed - ktr * prol
    d/dt(transit1) <- ktr * prol     - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(transit4) <- ktr * transit3 - ktr * transit4
    d/dt(transit5) <- ktr * transit4 - ktr * transit5
    d/dt(circ)     <- ktr * transit5 - ktr * circ

    prol(0)     <- circ0
    transit1(0) <- circ0
    transit2(0) <- circ0
    transit3(0) <- circ0
    transit4(0) <- circ0
    transit5(0) <- circ0
    circ(0)     <- circ0

    # ---- 5. Observation: absolute thrombocyte count (10^9 cells/L) ----
    # The circulating-platelet state IS the observed endpoint; `circ` is the
    # registered canonical compartment name for the circulating-cell pool of
    # a Friberg chain, so it is used directly as the observation variable.
    circ ~ add(addSd)
  })
}
