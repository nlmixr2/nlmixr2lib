AbdRahman_2020_chloroquine_wholeblood <- function() {
  description <- paste(
    "Semi-mechanistic population PK/PD model of chloroquine and its active",
    "metabolite desethylchloroquine, driven by WHOLE BLOOD drug",
    "concentrations, in 24 healthy adults inoculated with blood-stage",
    "Plasmodium vivax in a volunteer infection study (Abd-Rahman 2020). PK is",
    "a joint parent-metabolite model: chloroquine follows a two-compartment",
    "disposition with first-order absorption and elimination, and a fixed",
    "fraction FM = 0.18 of chloroquine elimination forms desethylchloroquine,",
    "which itself follows a two-compartment disposition with first-order",
    "elimination. Concentrations and doses are in molar units (chloroquine MW",
    "319.8, desethylchloroquine MW 291.8 g/mol). Whole-blood apparent volumes",
    "and clearances are about sixfold lower than the plasma values because",
    "chloroquine partitions into erythrocytes. PD is a delayed-effect",
    "(effect-compartment) parasite-clearance model: chloroquine in the effect",
    "compartment equilibrates with whole blood at rate ke0 and drives a",
    "sigmoid Emax parasite-killing rate kkill; the log10 P. vivax parasitaemia",
    "grows at a net rate kgrow and declines by kkill. Whole-blood EC50 was",
    "fixed to a literature value and the Hill coefficient fixed at 2.5 because",
    "no subject recrudesced. See the companion plasma model",
    "modellib('AbdRahman_2020_chloroquine_plasma').",
    sep = " "
  )
  reference <- paste(
    "Abd-Rahman AN, Marquart L, Gobeau N, Kummel A, Simpson JA, Chalon S,",
    "Mohrle JJ, McCarthy JS.",
    "Population Pharmacokinetics and Pharmacodynamics of Chloroquine in a",
    "Plasmodium vivax Volunteer Infection Study.",
    "Clin Pharmacol Ther. 2020;108(5):1055-1066.",
    "doi:10.1002/cpt.1893.",
    sep = " "
  )
  vignette <- "AbdRahman_2020_chloroquine"
  units <- list(time = "h", dosing = "umol", concentration = "umol/L")

  # The log10-parasitaemia state is paper-mechanistic (Abd-Rahman 2020
  # Supplementary Material S1, "Population pharmacokinetic-pharmacodynamic
  # modelling": dPL/dt = kgrow - kkill with PL the log-transformed parasite
  # count) and has no canonical density-scale counterpart, so it is declared
  # here rather than as the canonical `parasites` burden pool.
  paper_specific_compartments <- c(
    "parasitemia_log10"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Chloroquine and desethylchloroquine amounts are in
  # micromoles because the paper converted all concentrations and doses to
  # molar units (Methods, "Pharmacokinetic modelling").
  compartmentData <- list(
    depot = list(analyte = "chloroquine", units = "umol", specimen = "administration site", verified = TRUE),
    central = list(analyte = "chloroquine", units = "umol", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "chloroquine", units = "umol", specimen = "whole blood", verified = TRUE),
    central_dcq = list(analyte = "desethylchloroquine", units = "umol", specimen = "whole blood", verified = TRUE),
    peripheral1_dcq = list(analyte = "desethylchloroquine", units = "umol", specimen = "whole blood", verified = TRUE),
    effect = list(
      analyte = "chloroquine",
      units = "umol/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    parasitemia_log10 = list(
      analyte = "Plasmodium vivax parasitaemia, log10-transformed",
      units = "log10(parasites/mL)",
      specimen = "whole blood",
      verified = TRUE
    )
  )

  # Age, sex and body weight were screened during covariate model building
  # but none was retained in the final model (Results, "Population
  # pharmacokinetic modelling": "No covariate was included in the final PK
  # model"; Supplementary Material S1 reports the screened effects fell within
  # the +/-20% clinical-relevance region). They are recorded as documented-
  # but-unused so the covariate screen is preserved without convention
  # warnings for unreferenced covariates.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Screened on chloroquine and desethylchloroquine apparent clearance",
        "and central volume as a power function; not retained (narrow",
        "body-weight range; 95% CI of the effect within the +/-20%",
        "clinical-relevance region). Supplementary Material S1."
      )
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened during covariate model building; not retained."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Screened as a categorical effect on relative bioavailability, CL and",
        "central volume; not retained (95% CI within the +/-20% region).",
        "Supplementary Material S1."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 24L,
    n_studies = 1L,
    age_mean_sd = "25.6 (6.6) years",
    age_range = "19-44 years",
    weight_mean_sd = "73.3 (10.8) kg",
    weight_range = "57.2-99.5 kg",
    sex_female_pct = 45.8,
    race_ethnicity = c(White = 87.5, `Indigenous aboriginal` = 4.2, Latino = 4.2, Asian = 4.2),
    disease_state = paste(
      "Healthy malaria-naive adults inoculated intravenously with ~564 viable",
      "blood-stage Plasmodium vivax-infected erythrocytes on day 0 (induced",
      "blood-stage malaria volunteer infection study)."
    ),
    dose_range = paste(
      "Chloroquine phosphate tablets (Avloclor); total chloroquine base dose",
      "1.55 g for >= 60 kg adults (or 25 mg base/kg for < 60 kg) over 3 days.",
      "Per-subject base doses were 620 mg at 0 h, then 310 mg at 6, 24 and",
      "48 h (chloroquine phosphate converted to base by multiplying by 0.62).",
      "Treatment commenced on day 8 (cohort 1) or day 10 (cohorts 2 and 3)."
    ),
    regions = "Australia (QIMR Berghofer, Brisbane)",
    trial_registration = "ACTRN12616000174482",
    notes = paste(
      "Phase Ib, single-centre, open-label trial in three cohorts of eight",
      "subjects. Chloroquine and desethylchloroquine were measured by",
      "LC-MS/MS in plasma and whole blood (whole-blood LLOQ 1.0 ug/L for",
      "chloroquine and 1.75 ug/L for desethylchloroquine); P. vivax",
      "parasitaemia by 18S rRNA qPCR (limit of detection 10 parasites/mL).",
      "Whole blood was used as the surrogate matrix for parasite killing",
      "because whole-blood and erythrocyte concentrations were highly",
      "correlated. PK was fit separately to plasma and to whole blood data;",
      "this file is the whole-blood model. Modelling used Monolix 4.3.3 with",
      "left-censoring of below-limit data. Structural PK estimates are",
      "Abd-Rahman 2020 Table 2 (whole-blood columns); PD estimates are",
      "Table 3."
    )
  )

  ini({
    # =========================================================================
    # Chloroquine PK -- Abd-Rahman 2020 Table 2, whole-blood-sample columns.
    # IIV values in Table 2 are Monolix omega estimates (standard deviations
    # of the log-normal random effects); the nlmixr2 variances below are their
    # squares. Concentrations and doses are molar (chloroquine MW 319.8 g/mol).
    # =========================================================================
    lka     <- log(0.574); label("Chloroquine first-order absorption rate constant ka (1/h)")                 # Table 2 whole-blood row 'k a (1/h) = 0.574 (95% CI 0.378-0.770)'
    lcl     <- log(8.96);  label("Apparent chloroquine clearance CL_CQ/F (L/h)")                              # Table 2 whole-blood row 'CL CQ/F (L/h) = 8.96 (95% CI 8.10-9.82)'
    lvc     <- log(560);   label("Apparent chloroquine central volume Vc_CQ/F (L)")                           # Table 2 whole-blood row 'V c CQ/F (L) = 560 (95% CI 456-664)'
    lq      <- log(38.5);  label("Apparent chloroquine intercompartmental clearance Q1_CQ/F (L/h)")           # Table 2 whole-blood row 'Q 1 CQ/F (L/h) = 38.5 (95% CI 30.1-46.9)'
    lvp     <- log(1230);  label("Apparent chloroquine peripheral volume Vp1_CQ/F (L)")                       # Table 2 whole-blood row 'V p1 CQ/F (L) = 1230 (95% CI 1079-1381)'
    lfdepot <- fixed(log(1)); label("Chloroquine relative bioavailability Frel (fraction)")                   # Table 2 whole-blood row 'F rel (%) = 100 (fixed)'; typical F fixed to 100%, only its IIV estimated

    # =========================================================================
    # Desethylchloroquine PK -- Abd-Rahman 2020 Table 2, whole-blood columns.
    # Formed from chloroquine with fixed fraction FM = 0.18; molar units
    # (desethylchloroquine MW 291.8 g/mol). Apparent volumes and clearances
    # were estimated assuming FM = 0.18 (Results, "Population pharmacokinetic
    # modelling").
    # =========================================================================
    lcl_dcq <- log(4.42);  label("Apparent desethylchloroquine clearance CL_DCQ/F (L/h)")                     # Table 2 whole-blood row 'CL DCQ/F (L/h) = 4.42 (95% CI 3.87-4.97)'
    lvc_dcq <- log(16.1);  label("Apparent desethylchloroquine central volume Vc_DCQ/F (L)")                  # Table 2 whole-blood row 'V c DCQ/F (L) = 16.1 (95% CI 11.6-20.6)'
    lq_dcq  <- log(4.46);  label("Apparent desethylchloroquine intercompartmental clearance Q1_DCQ/F (L/h)")  # Table 2 whole-blood row 'Q 1 DCQ/F (L/h) = 4.46 (95% CI 3.40-5.52)'
    lvp_dcq <- log(259);   label("Apparent desethylchloroquine peripheral volume Vp1_DCQ/F (L)")              # Table 2 whole-blood row 'V p1 DCQ/F (L) = 259 (95% CI 196-322)'

    fm <- fixed(0.18); label("Fraction of chloroquine elimination forming desethylchloroquine (fraction)")    # Results, "Population pharmacokinetic modelling": 'it was assumed that 18% of chloroquine is converted to desethylchloroquine' (from urinary-recovery literature, references 9 and 12); fixed for identifiability

    # =========================================================================
    # PD -- delayed-effect parasite-clearance model. Abd-Rahman 2020 Table 3.
    # Emax, kgrow and the baseline log10 parasitaemia are shared with the
    # plasma model; EC50 and ke0 are whole-blood-specific.
    # =========================================================================
    lemax  <- log(0.213);  label("Maximum chloroquine parasite-killing rate EmaxCQ (1/h)")                    # Table 3 row 'E maxCQ (1/h) = 0.213 (95% CI 0.196-0.230)'; shared across matrices
    lec50  <- fixed(log(0.28)); label("Chloroquine effect-site EC50 in whole blood, EC50_CQ_wholeblood (umol/L)") # Table 3 row 'EC 50CQ whole blood (umol/L) = 0.28 (fixed)'; fixed to a previously reported relapse-based value
    lhill  <- fixed(log(2.5));   label("Hill coefficient of the chloroquine killing sigmoid in whole blood (unitless)") # Table 3 row 'Hill CQ whole blood = 2.5 (fixed)'
    lke0   <- log(0.0288); label("Chloroquine effect-compartment equilibration rate ke0 in whole blood (1/h)")     # Table 3 row 'k e0CQ whole blood (1/h) = 0.0288 (95% CI 0.0212-0.0364)'
    lkgrow <- log(0.059);  label("Net P. vivax growth rate constant kgrow (1/h)")                             # Table 3 row 'k grow (1/h) = 0.059 (95% CI 0.056-0.062)'; shared across matrices
    plbase <- -3.36;       label("Baseline log10 P. vivax parasitaemia PLbase (log10 parasites/mL)")          # Table 3 row 'PL base (log 10 parasite/mL) = -3.36 (95% CI -3.97 to -2.75)'; normally distributed, shared across matrices

    # =========================================================================
    # Inter-individual variability. PK: Table 2 whole-blood omega estimates
    # (standard deviations), squared for nlmixr2 variances. Whole blood adds
    # IIV on Vc_DCQ (absent for plasma). PD: Table 3 omega estimates squared;
    # the baseline PLbase random effect is NORMAL (added directly, not
    # log-normal) per Methods. EC50 and ke0 IIV were fixed.
    # =========================================================================
    etalka     ~ 0.218089   # Table 2 whole-blood 'omega k a = 0.467'; variance = 0.467^2
    etalfdepot ~ 0.046225   # Table 2 whole-blood 'omega F rel = 0.215'; variance = 0.215^2
    etalcl     ~ 0.013689   # Table 2 whole-blood 'omega CL CQ/F = 0.117'; variance = 0.117^2
    etalvc     ~ 0.051984   # Table 2 whole-blood 'omega Vc CQ/F = 0.228'; variance = 0.228^2
    etalcl_dcq ~ 0.047961   # Table 2 whole-blood 'omega CL DCQ/F = 0.219'; variance = 0.219^2
    etalvc_dcq ~ 0.147456   # Table 2 whole-blood 'omega Vc DCQ/F = 0.384'; variance = 0.384^2

    etaplbase  ~ 2.1609            # Table 3 'omega PL base = 1.47'; variance = 1.47^2; NORMAL baseline random effect
    etalkgrow  ~ 0.018769          # Table 3 'omega k grow = 0.137'; variance = 0.137^2
    etalemax   ~ 0.037249          # Table 3 'omega E maxCQ = 0.193'; variance = 0.193^2
    etalec50   ~ fixed(0.09)       # Table 3 'omega EC 50CQ whole blood = 0.3', not estimated; variance = 0.3^2
    etalke0    ~ fixed(0.01)       # Table 3 'omega k e CQ whole blood = 0.1', not estimated; variance = 0.1^2

    # =========================================================================
    # Residual error. Table 2 reports proportional residual SDs for the drug
    # concentrations (Monolix 'epsilon prop'); Table 3 an additive residual SD
    # on the log10 parasitaemia measurements.
    # =========================================================================
    propSd     <- 0.221; label("Proportional residual error on whole-blood chloroquine concentration (fraction)")            # Table 2 whole-blood row 'epsilon prop CQ = 0.221 (95% CI 0.199-0.243)'
    propSd_dcq <- 0.253; label("Proportional residual error on whole-blood desethylchloroquine concentration (fraction)")    # Table 2 whole-blood row 'epsilon prop DCQ = 0.253 (95% CI 0.228-0.278)'
    addSd_parasitemia_log10 <- 0.379; label("Additive residual error on log10 parasitaemia (log10 parasites/mL)")            # Table 3 row 'epsilon add whole blood = 0.379 (95% CI 0.357-0.401)'
  })

  model({
    # -----------------------------------------------------------------------
    # Individual PK parameters (log-normal IIV; Frel fixed to 1 with IIV).
    # -----------------------------------------------------------------------
    ka     <- exp(lka + etalka)
    cl     <- exp(lcl + etalcl)
    vc     <- exp(lvc + etalvc)
    q      <- exp(lq)
    vp     <- exp(lvp)
    fdepot <- exp(lfdepot + etalfdepot)

    cl_dcq <- exp(lcl_dcq + etalcl_dcq)
    vc_dcq <- exp(lvc_dcq + etalvc_dcq)
    q_dcq  <- exp(lq_dcq)
    vp_dcq <- exp(lvp_dcq)

    # -----------------------------------------------------------------------
    # Individual PD parameters. PLbase carries a NORMAL random effect (added
    # directly); the others are log-normal.
    # -----------------------------------------------------------------------
    emax     <- exp(lemax + etalemax)
    ec50     <- exp(lec50 + etalec50)
    hill     <- exp(lhill)
    ke0      <- exp(lke0 + etalke0)
    kgrow    <- exp(lkgrow + etalkgrow)
    plbase_i <- plbase + etaplbase

    # -----------------------------------------------------------------------
    # Micro-constants. Molar 1:1 conversion of chloroquine to
    # desethylchloroquine, so the metabolite formation flux is
    # fm * kel * central (umol chloroquine eliminated by metabolism = umol
    # desethylchloroquine formed).
    # -----------------------------------------------------------------------
    kel     <- cl / vc
    kel_dcq <- cl_dcq / vc_dcq

    # -----------------------------------------------------------------------
    # PK ODEs: chloroquine two-compartment with first-order absorption, and
    # desethylchloroquine two-compartment formed from chloroquine central.
    # -----------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <- ka * depot - kel * central -
      q / vc * central + q / vp * peripheral1
    d/dt(peripheral1) <- q / vc * central - q / vp * peripheral1

    d/dt(central_dcq) <- fm * kel * central - kel_dcq * central_dcq -
      q_dcq / vc_dcq * central_dcq + q_dcq / vp_dcq * peripheral1_dcq
    d/dt(peripheral1_dcq) <- q_dcq / vc_dcq * central_dcq -
      q_dcq / vp_dcq * peripheral1_dcq

    f(depot) <- fdepot

    # Molar concentrations (umol/L).
    Cc     <- central / vc
    Cc_dcq <- central_dcq / vc_dcq

    # -----------------------------------------------------------------------
    # PD: effect-compartment delay and sigmoid Emax parasite killing.
    #   dCe/dt = ke0 * (Cc - Ce)                        (biophase, umol/L)
    #   kkill  = Emax * Ce^Hill / (EC50^Hill + Ce^Hill) (1/h)
    #   dPL/dt = kgrow - kkill        with PL = log10 parasitaemia
    # The supplement writes dPL/dt = kgrow - kkill with PL the log-transformed
    # parasite count. Because the paper's own secondary parameters use natural
    # e (PMR48 = exp(kgrow*48) = 17.4 and PCt1/2 = ln(2)/(Emax-kgrow) = 4.5 h),
    # kgrow and Emax are natural-log growth/kill rates. Storing the state in
    # log10 (so PLbase and the additive residual apply directly) makes the
    # ODE dPL10/dt = (kgrow - kkill)/ln(10), which is algebraically identical
    # to the natural-log implementation and reproduces every reported
    # secondary parameter. See vignette source trace.
    # -----------------------------------------------------------------------
    ce    <- effect
    kkill <- emax * ce^hill / (ec50^hill + ce^hill)

    d/dt(effect) <- ke0 * (Cc - effect)
    d/dt(parasitemia_log10) <- (kgrow - kkill) / log(10)
    parasitemia_log10(0) <- plbase_i

    Cc     ~ prop(propSd)
    Cc_dcq ~ prop(propSd_dcq)
    parasitemia_log10 ~ add(addSd_parasitemia_log10)
  })
}
