vanErp_2010_sunitinib <- function() {
  description <- paste0(
    "One-compartment population PK model for oral sunitinib in cancer ",
    "patients with a mechanism-specific grapefruit-juice (GJ) drug-",
    "interaction module. Sunitinib is absorbed first-order (ka, tlag) into ",
    "a single central compartment with linear elimination (CL/F, Vd/F). A ",
    "paper-specific intestinal CYP3A4-activity state (baseline 1, recovery ",
    "first-order with t1/2 = 23 h fixed from Greenblatt 2003) is fully ",
    "depleted to 0 by each GJ ingestion event. The relative bioavailability ",
    "is F = 1 + deltaF * (1 - cyp3a4), so simultaneous GJ + sunitinib ",
    "intake gives F = 1.11 (deltaF = 0.11) and the GJ-induced increase in ",
    "sunitinib exposure decays back to baseline with the CYP3A4 recovery ",
    "half-life (8.9% at 7 h, 5.3% at 24 h, 1.3% at 72 h, 0.07% at 1 week ",
    "after the last GJ dose). No covariates were retained in the final ",
    "model. Eight metastatic-cancer patients (1 female / 7 male, age 41-78 ",
    "years) on chronic sunitinib 25-50 mg once daily contributed 268 ",
    "plasma concentrations."
  )
  reference <- paste(
    "van Erp NP, Baker SD, Zandvliet AS, Ploeger BA, den Hollander M,",
    "Chen Z, den Hartigh J, Konig-Quartel JMC, Guchelaar H-J, Gelderblom H.",
    "Marginal increase of sunitinib exposure by grapefruit juice.",
    "Cancer Chemother Pharmacol. 2011;67(3):695-703.",
    "doi:10.1007/s00280-010-1367-0.",
    "CYP3A4 recovery half-life (23 h) fixed from",
    "Greenblatt DJ et al., Clin Pharmacol Ther. 2003;74(2):121-129.",
    sep = " "
  )
  vignette <- "vanErp_2010_sunitinib"
  paper_specific_compartments <- c("cyp3a4")

  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list()

  covariatesDataExcluded <- list()

  population <- list(
    species        = "human",
    n_subjects     = 8L,
    n_observations = 268L,
    n_studies      = 1L,
    age_range      = "41-78 years (median 54)",
    sex_female_pct = 12.5,
    race_ethnicity = NULL,
    disease_state  = paste0(
      "Adult patients with advanced solid tumours for which sunitinib is ",
      "registered first-line (metastatic renal cell carcinoma) or second-",
      "line (gastrointestinal stromal tumour), or who had exhausted other ",
      "treatment options. WHO performance status <= 2; adequate bone-",
      "marrow, renal and hepatic function (creatinine clearance >= 60 ",
      "mL/min, bilirubin <= 1.75x ULN)."
    ),
    dose_range     = paste0(
      "Oral sunitinib 25, 37.5 or 50 mg once daily in a 4-weeks-on / ",
      "2-weeks-off cycle. PK sampling occurred at steady state (day 14-20 ",
      "without GJ; day 28 with GJ). On days 25-27, patients took 200 mL ",
      "grapefruit juice three times a day; sunitinib was co-administered ",
      "with the morning GJ on day 28."
    ),
    regions        = "Netherlands (Leiden University Medical Center).",
    notes          = paste0(
      "Cohort and dosing details from van Erp 2011 Table 1 and Methods. ",
      "Baseline serum creatinine median 77 uM (range 56-122); total ",
      "bilirubin median 9 uM (range 6-15); ALT median 39 U/L (range ",
      "18-68); Hb median 8.7 mM (range 7.0-9.4); WBC median 5.5x10^9/L ",
      "(range 3.5-38.2); thrombocytes median 196x10^9/L (range 149-318). ",
      "No covariates were tested or retained in the final population PK ",
      "model. Sample size of 8 was sized for a paired two-sided test with ",
      "80% power to detect a 25% change in sunitinib exposure."
    )
  )

  ini({
    # -----------------------------------------------------------------------
    # Sunitinib disposition -- van Erp 2011 Table 2 final-model estimates.
    # The reported parameters are apparent (relative to oral bioavailability
    # under no-GJ conditions, which is implicitly anchored to F = 1).
    # -----------------------------------------------------------------------
    lka    <- log(0.468); label("Sunitinib absorption rate ka (1/h)")                                      # Table 2 (ka 0.468 1/h, RSE 27.6%)
    lcl    <- log(50.5);  label("Sunitinib apparent clearance CL/F (L/h)")                                 # Table 2 (Cl/F 50.5 L/h, RSE 20.6%)
    lvc    <- log(3210);  label("Sunitinib apparent central volume of distribution Vd/F (L)")              # Table 2 (Vd/F 3210 L, RSE 7.8%)
    ltlag  <- log(0.487); label("Sunitinib absorption lag time (h)")                                      # Table 2 (Absorption lag time 0.487 h, RSE 6.8%)

    # -----------------------------------------------------------------------
    # Grapefruit-juice CYP3A4 interaction parameters.
    #   ldeltaF -- log of the relative-bioavailability increment under
    #              full intestinal CYP3A4 inhibition. The published
    #              Relative F of 1.11 (Table 2; RSE 70%; profile-likelihood
    #              95% CI 1.042-1.182) corresponds to deltaF = 0.11.
    #   lkdeg   -- log of the first-order CYP3A4-activity recovery rate
    #              constant, FIXED at log(2)/23 h^-1 from Greenblatt 2003
    #              [ref 27 of the paper], as stated in the Methods.
    # -----------------------------------------------------------------------
    ldeltaF <- log(0.11);             label("Relative-F increment under full intestinal CYP3A4 inhibition (unitless)")  # Table 2 (Relative F 1.11; deltaF = 1.11 - 1 = 0.11)
    lkdeg   <- fixed(log(log(2)/23)); label("Intestinal CYP3A4-activity recovery rate (1/h)")              # Methods (recovery half-life set to 23 h from Greenblatt 2003)

    # -----------------------------------------------------------------------
    # Inter-individual variability -- van Erp 2011 Table 2 (exponential / log-
    # normal IIV; CV% reported). Convert to omega^2 via log(1 + CV^2):
    #   ka:  CV = 63.9% -> omega^2 = log(1 + 0.639^2) = 0.342  (RSE 42.9%)
    #   Cl:  CV = 67.9% -> omega^2 = log(1 + 0.679^2) = 0.379  (RSE 42.7%)
    # No IIV was estimated on Vd/F, Relative F, lag time, or residual error
    # ("nd" in Table 2).
    # -----------------------------------------------------------------------
    etalka ~ 0.342
    etalcl ~ 0.379

    # -----------------------------------------------------------------------
    # Residual error -- van Erp 2011 Table 2: proportional residual error
    # 16.3% (RSE 22.9%). Encoded as the proportional SD on the linear scale
    # (propSd = 0.163), with NONMEM INTER between IIV and residual error.
    # -----------------------------------------------------------------------
    propSd <- 0.163; label("Proportional residual SD on sunitinib (fraction)")  # Table 2 (proportional residual error 16.3%, RSE 22.9%)
  })

  model({
    # -------------------------------------------------------------------
    # Individual PK parameters.
    # -------------------------------------------------------------------
    ka   <- exp(lka + etalka)
    cl   <- exp(lcl + etalcl)
    vc   <- exp(lvc)
    tlag <- exp(ltlag)

    # -------------------------------------------------------------------
    # Grapefruit-juice / CYP3A4 interaction module.
    #
    # The state `cyp3a4` is the fractional intestinal CYP3A4 activity:
    # 1 = full baseline activity, 0 = fully depleted. Recovery follows
    # first-order regeneration toward 1 at the rate kdeg (t1/2 = 23 h,
    # fixed). Each grapefruit-juice ingestion event fully depletes
    # CYP3A4 activity to 0 -- in the vignette this is implemented by
    # bolus replacement events (evid = 5) on the `cyp3a4` compartment.
    #
    # The relative bioavailability of sunitinib is then
    #
    #   F = 1 + deltaF * (1 - cyp3a4)
    #
    # so simultaneous GJ + sunitinib (cyp3a4 = 0) gives F = 1 + deltaF =
    # 1.11 (Table 2), and the GJ-induced uplift decays back toward 1 as
    # the CYP3A4 activity regenerates. The model assumes no GJ effect on
    # CL, Vd, ka, or tlag (van Erp 2011 Discussion: only intestinal
    # CYP3A4 is inhibited, so the GJ effect is restricted to absorption).
    # -------------------------------------------------------------------
    kdeg   <- exp(lkdeg)
    deltaF <- exp(ldeltaF)

    d/dt(cyp3a4) <- kdeg * (1 - cyp3a4)
    cyp3a4(0)    <- 1

    Frel <- 1 + deltaF * (1 - cyp3a4)

    # -------------------------------------------------------------------
    # PK ODE system: 1-compartment with first-order absorption from a
    # depot. Bioavailability and absorption lag time are applied to the
    # depot compartment. CL/Vc gives the elimination rate constant.
    # -------------------------------------------------------------------
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - cl / vc * central

    f(depot)    <- Frel
    alag(depot) <- tlag

    # -------------------------------------------------------------------
    # Observation. Volumes are in L and amounts in mg, giving plasma in
    # mg/L; multiply by 1000 to obtain ng/mL (the units in which the
    # paper reports concentrations and the Table 2 derived parameters).
    # -------------------------------------------------------------------
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
