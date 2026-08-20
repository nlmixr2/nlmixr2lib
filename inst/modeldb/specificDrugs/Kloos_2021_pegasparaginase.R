Kloos_2021_pegasparaginase <- function() {
  description <- "One-compartment population PK model with time-dependent (split-point) clearance for intravenous PEGasparaginase in pediatric acute lymphoblastic leukemia patients treated on the Dutch Childhood Oncology Group ALL-11 protocol (Kloos 2021). Clearance and volume of distribution are normalized to body surface area; clearance is constant for the first 12.7 days after a dose and then increases linearly with time after dose as the polyethylene glycol moiety is hydrolyzed. Clearance is 38 percent higher during an active infection and 11-19 percent lower outside induction, encoded as multiplicative treatment-phase factors with protocol 1A as the reference. Inter-individual variability on clearance is shared with the volume of distribution through a scaled eta; inter-occasion variability on clearance uses one occasion per administered dose. Combined proportional and additive residual error."
  reference   <- "Kloos RQH, Mathot R, Pieters R, van der Sluis IM. Individualized dosing guidelines for PEGasparaginase and factors influencing the clearance: a population pharmacokinetic model. Haematologica. 2021;106(5):1254-1261. doi:10.3324/haematol.2019.242289. Erratum (corrected Table 5 dose-adjustment nomogram): Haematologica. 2023;108(9):2558. doi:10.3324/haematol.2023.283685"
  vignette    <- "Kloos_2021_pegasparaginase"
  units       <- list(time = "day", dosing = "IU", concentration = "IU/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    central = list(analyte = "pegasparaginase", units = "IU", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    BSA = list(
      description        = "Body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Computed by the Mosteller formula (Kloos 2021 Methods, Population pharmacokinetic analysis, citing Mosteller 1987). Index-dataset median 0.76 m^2 (IQR 0.65-1.05; Supplemental Table 2); measured at the start of PEGasparaginase therapy. BSA normalizes BOTH clearance and volume of distribution (no exponent; a linear proportionality, not an allometric power term) -- Kloos 2021 Table 2 final-model equations and Online Supplementary Appendix, Supplemental Results, PK analysis: 'Cl = Theta1 * e^(eta + eta_IOV) * BSA' and 'Vd = Theta3 * e^(Theta4 * eta) * BSA'. Normalization by BSA reduced unexplained inter-individual variability in clearance from 29.6 percent to 24.1 percent (Kloos 2021 Results, Structural model). The Monte Carlo simulations that generated the dosing guidelines spanned BSA 0.52-2.3 m^2 (Supplemental Results, Simulations).",
      source_name        = "BSA"
    ),
    DIS_INFECT_ACTIVE = list(
      description        = "Active clinical infection episode indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no active infection)",
      notes              = "Time-varying per record. Kloos 2021 defines an infection as 'fever (>38 degrees Celsius) and hospital admission or prescription of antibiotics' (Online Supplementary Appendix, Supplemental Table 2 footnote; also Methods, Covariate analysis). 19 of 92 index-dataset patients had at least one infection (Supplemental Table 2). Multiplies clearance by 1.38 (Table 2 final model); the authors attribute the increase to activation of the mononuclear phagocyte system, which clears PEGylated asparaginase (Discussion). ICU admission was screened alongside infection in the multivariate analysis but did not improve the fit on top of infection and treatment phase (delta OFV -2.6; Supplemental Results, PK analysis) and is therefore not a covariate of this model.",
      source_name        = "INFECTION"
    ),
    TRTPH_1B = list(
      description        = "DCOG ALL-11 protocol 1B treatment-phase indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (protocol 1A induction, the model reference phase)",
      notes              = "Time-varying per record. Protocol 1B follows induction protocol 1A and comprises cyclophosphamide, cytarabine and 6-mercaptopurine plus a 1,500 IU/m^2 PEGasparaginase dose at day 40 (Supplemental Table 1). Its clearance effect was FIXED to 1 (i.e., pooled with the 1A reference) because only two patients were treated as high risk and the association between the high-risk blocks and PEGasparaginase clearance could not be estimated reliably (Kloos 2021 Results, Covariate analysis; Table 2 lists the 1B row as '1 (fix)'). Retained as an explicit column so the model's phase mapping is complete and the fixed status is visible in ini().",
      source_name        = "Treatment phase 1B"
    ),
    TRTPH_M = list(
      description        = "DCOG ALL-11 protocol M treatment-phase indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (protocol 1A induction, the model reference phase)",
      notes              = "Time-varying per record. Protocol M comprises 6-mercaptopurine plus high-dose methotrexate (Supplemental Table 1); 69 of 816 index-dataset levels fall in this phase (Table 1). Multiplies clearance by 0.87 (Table 2 final model).",
      source_name        = "Treatment phase M"
    ),
    TRTPH_MR_INTENS = list(
      description        = "DCOG ALL-11 medium-risk-group intensification treatment-phase indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (protocol 1A induction, the model reference phase)",
      notes              = "Time-varying per record. Medium-risk intensification comprises dexamethasone, vincristine, 6-mercaptopurine and (in TEL/AML1-negative patients) doxorubicin, or methotrexate in place of doxorubicin for TEL/AML1-positive patients (Figure 1, Supplemental Table 1); 168 of 816 index-dataset levels fall in this phase (Table 1). Multiplies clearance by 0.89 (Table 2 final model).",
      source_name        = "Treatment phase MR intensification"
    ),
    TRTPH_MR_MAINT = list(
      description        = "DCOG ALL-11 medium-risk-group maintenance treatment-phase indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (protocol 1A induction, the model reference phase)",
      notes              = "Time-varying per record. Medium-risk maintenance comprises dexamethasone, vincristine, methotrexate and 6-mercaptopurine (Figure 1, Supplemental Table 1); 250 of 816 index-dataset levels fall in this phase, the largest single phase (Table 1). Multiplies clearance by 0.81 (Table 2 final model).",
      source_name        = "Treatment phase MR maintenance"
    ),
    TRTPH_SR_IV = list(
      description        = "DCOG ALL-11 standard-risk-group protocol IV treatment-phase indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (protocol 1A induction, the model reference phase)",
      notes              = "Time-varying per record. Standard-risk protocol IV comprises dexamethasone and vincristine plus a single individualized PEGasparaginase dose at day 1 (Figure 1, Supplemental Table 1); 38 of 816 index-dataset levels fall in this phase (Table 1). Multiplies clearance by 0.81 (Table 2 final model).",
      source_name        = "Treatment phase SR protocol IV"
    ),
    OCC = list(
      description        = "Integer occasion index for inter-occasion variability on clearance",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Kloos 2021 Online Supplementary Appendix, Supplemental Methods, Population PK analysis defines the occasion explicitly: 'inter-individual variability and inter-occasion variability, with an occasion defined as administration of a new dose'. This model instantiates three occasions, matching the three fixed 1,500 IU/m^2 PEGasparaginase doses given biweekly during induction on protocols 1A (days 12 and 26) and 1B (day 40) per Supplemental Table 1 and Kloos 2021 Methods, Patients and treatment protocol. Data assemblers set OCC = 1, 2, 3, ... incrementing at each new administration; records with OCC outside 1-3 receive no IOV contribution.",
      source_name        = "occasion (administration of a new dose)"
    )
  )

  # Covariates screened by Kloos 2021 but NOT retained in the final model.
  # Documented for provenance only; these names are deliberately absent from
  # model() (see Kloos 2021 Table 3 univariate analysis and Online
  # Supplementary Appendix, Supplemental Results, PK analysis).
  covariatesDataExcluded <- list(
    WBC = list(
      description        = "Leukocyte count",
      units              = "10^9/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened in the univariate analysis (Kloos 2021 Table 3: effect -0.09, 95% CI -0.13 to -0.05, delta OFV -22.9) as a power covariate centred on the cohort median 2.4 x 10^9/L (Supplemental Methods). Not retained: the supplement records that leukocytes 'had large relative standard errors and the 95% confidence interval included 0' (Supplemental Results, PK analysis). Index-dataset median 2.4 x 10^9/L (IQR 1.5-4.0), 8 percent of measurements missing (Supplemental Table 2)."
    ),
    CREAT = list(
      description        = "Serum creatinine",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened in the univariate analysis (Kloos 2021 Table 3: effect -0.21, 95% CI -0.36 to -0.07, delta OFV -17.0). Not retained for the same reason as WBC (large RSE / CI including 0; Supplemental Results, PK analysis). Index-dataset median 27 umol/L (IQR 22-33), 62 percent of measurements missing (Supplemental Table 2)."
    ),
    ADA_POS = list(
      description        = "Anti-asparaginase antibody status (anti-native E. coli asparaginase and anti-PEGasparaginase antibodies, reported as extinction / optical-density values)",
      units              = "(binary; source reported continuous optical density)",
      type               = "binary",
      reference_category = "0 (antibody negative)",
      notes              = "Screened in the univariate analysis as continuous extinction values (Kloos 2021 Table 3: anti-native E. coli asparaginase antibodies effect 0.05, 95% CI -0.01 to 0.11, delta OFV -12.9; anti-PEGasparaginase antibodies effect 0.04, 95% CI -0.01 to 0.10, delta OFV -7.3). Neither was retained (large RSE, 95% CI including 0; Supplemental Results, PK analysis). Only 2 of 92 index patients developed an allergy and 2 developed silent inactivation (Supplemental Table 2), so the cohort carried little antibody signal. Registered here under the canonical binary ADA_POS name even though the source screened a continuous optical-density readout, because the concept excluded is anti-drug-antibody status."
    )
  )

  population <- list(
    species             = "human",
    n_subjects          = 92L,
    n_studies           = 1L,
    n_observations      = 816L,
    age_range           = "1-18 years (eligibility); index-dataset median 4.8 years, IQR 3.3-8.2",
    weight_range        = "Index-dataset median 19.2 kg, IQR 14.9-29.3",
    bsa_range           = "Index-dataset median 0.76 m^2, IQR 0.65-1.05",
    sex_female_pct      = 44.6,
    race_ethnicity      = "Not reported (single-center Dutch cohort with additional trough samples from other Dutch pediatric oncology centers).",
    disease_state       = "Newly diagnosed pediatric acute lymphoblastic leukemia treated per the Dutch Childhood Oncology Group (DCOG) ALL-11 protocol between November 2014 and May 2017. Risk stratification after induction: standard risk 13 (14 percent), medium risk 76 (83 percent), high risk 2 (2 percent), not stratified 1 (1 percent) (Supplemental Table 2).",
    dose_range          = "Three fixed doses of 1,500 IU/m^2 biweekly during induction (protocols 1A and 1B), then individualized doses guided by therapeutic drug monitoring of asparaginase activity; standard-risk patients received one individualized dose in protocol IV, medium-risk patients up to 14 individualized doses during intensification and maintenance. Each dose was administered intravenously over 1 hour.",
    regions             = "Netherlands (Sophia Children's Hospital - Erasmus MC, Rotterdam, plus trough samples from other Dutch pediatric oncology centers via the national TDM program).",
    bioanalytic_methods = "Asparaginase activity in serum by the L-aspartic beta-hydroxamate (AHA) assay with photometric detection of indooxine at 690 nm; lower limit of quantification 10 IU/L. Values below the LLQ were excluded from the analysis (M2 / M3 methods gave unstable runs).",
    validation_cohort   = "An independent validation dataset of 28 patients / 405 samples, obtained by randomly selecting 25 percent of the total 120-patient population, was used for external validation (Table 1, Supplemental Table 2).",
    notes               = "Baseline demographics from Kloos 2021 Supplemental Table 2; sample distribution from Table 1. NONMEM 7.2 with FOCE+I on log-transformed asparaginase activity levels. One-compartment model; addition of a second compartment did not improve the fit, and first-order, zero-order and Michaelis-Menten elimination all failed to describe the data, so a split-point time-dependent clearance model after Wurthwein 2017 was adopted. The final model was assessed by bootstrap (1000 replicates), prediction-corrected VPC, and external validation. Trial registration NL50250.078.14 (CCMO register)."
  )

  ini({
    # ---- Structural PK ----------------------------------------------------
    # Kloos 2021 Table 2, 'Final model' column. Clearance and volume are
    # reported per m^2 of body surface area; BSA multiplies both inside
    # model(), so these are the values for a 1 m^2 patient.
    lcl <- log(0.084); label("Clearance per m^2 BSA during the first 12.7 days after a dose (CL, L/day/m^2)")  # Kloos 2021 Table 2 final model (CL 0.084 L/day/m^2, RSE 4.4%; bootstrap 0.084, 95% CI 0.078-0.090)
    lvc <- log(0.94);  label("Volume of distribution per m^2 BSA (Vd, L/m^2)")                                 # Kloos 2021 Table 2 final model (Vd 0.94 L/m^2, RSE 4.5%; bootstrap 0.94, 95% CI 0.87-1.01)

    # ---- Time-dependent clearance ----------------------------------------
    # PEGasparaginase clearance is constant for the first `tsplit` days after
    # a dose and then rises linearly with time after dose, as the PEG moiety
    # is hydrolyzed from the molecule and the residual E. coli asparaginase
    # (with a linker attached) is cleared far faster (Kloos 2021 Discussion).
    #
    # IMPORTANT -- ADDITIVE vs MULTIPLICATIVE SLOPE. The printed equations in
    # Table 2 and in the Online Supplementary Appendix give the MULTIPLICATIVE
    # form Cl_ind = 1 + Theta2 * (TAD - split point). This model implements the
    # ADDITIVE form CL = Theta1 + Theta2 * (TAD - split point) instead, i.e. it
    # treats the printed 'Cl_ind = 1 + Theta2 * (TAD - split)' as missing a
    # / Theta1. The additive reading is the one consistent with the paper's own
    # numeric results: the Table 2 row label gives this parameter's units as
    # L/day/m^2/day (absolute clearance per day, meaningless in the
    # multiplicative form, where Theta2 would be 1/day); the abstract, Results
    # (Structural model) and Results (Covariate analysis) all phrase it as CL
    # 'increasing by 0.082 L/day/m^2/day'; and the additive form exactly
    # reproduces the three half-lives the authors publish for the structural
    # model (8.50 / 4.14 / 2.74 days vs the stated 8.5 / 4.1 / 2.7) where the
    # multiplicative form misses them by 1.9-2.7x. See the vignette's
    # 'Assumptions and deviations' section for the full evidence table. This
    # reading was ratified by the operator (sidecar request-001 / response-001,
    # question q1, option A).
    lcl_time <- log(0.082); label("Slope of the linear clearance increase after the split point (L/day/m^2/day)")  # Kloos 2021 Table 2 final model ('Slope CL_ind' 0.082 L/day/m^2/day, RSE 20.5%; bootstrap 0.080, 95% CI 0.052-0.115)
    tsplit   <- 12.7;       label("Time after dose at which clearance starts to increase (split point, days)")     # Kloos 2021 Table 2 final model (split point 12.7 days, RSE 0.2%; bootstrap 12.7, 95% CI 11.8-13.1)

    # ---- Covariate effects on clearance ----------------------------------
    # All covariate effects are multiplicative on CL and are written in the
    # paper as <estimate>^<indicator> (Kloos 2021 Table 2 final-model
    # equations). Because every indicator is binary, <estimate>^1 = <estimate>
    # and <estimate>^0 = 1, so the reference level (protocol 1A, no infection)
    # collapses each factor to 1.
    e_infect_cl <- 1.38; label("Multiplicative effect of an active infection on CL (unitless)")  # Kloos 2021 Table 2 final model (Infection 1.38, RSE 10.5%; bootstrap 1.38, 95% CI 1.15-1.67)

    # Treatment phase, protocol 1A (induction) = reference. Protocol 1B is
    # FIXED to 1 (pooled with the 1A reference) because only two patients were
    # treated as high risk and the effect was not estimable -- Kloos 2021
    # Table 2 prints '1 (fix)' for this row and Results (Covariate analysis)
    # states 'The CL in phase 1B was equal to 1A ... could not be estimated
    # reliably and was fixed to one.'
    e_phase1b_cl      <- fixed(1); label("Multiplicative effect of DCOG ALL-11 protocol 1B on CL relative to 1A (unitless)")                  # Kloos 2021 Table 2 final model (1B '1 (fix)')
    e_phasem_cl       <- 0.87;     label("Multiplicative effect of DCOG ALL-11 protocol M on CL relative to 1A (unitless)")                   # Kloos 2021 Table 2 final model (M 0.87, RSE 5.2%; bootstrap 0.87, 95% CI 0.80-0.95)
    e_phasemrint_cl   <- 0.89;     label("Multiplicative effect of medium-risk intensification on CL relative to 1A (unitless)")              # Kloos 2021 Table 2 final model (MRG intens. 0.89, RSE 5.2%; bootstrap 0.88, 95% CI 0.82-0.98)
    e_phasemrmaint_cl <- 0.81;     label("Multiplicative effect of medium-risk maintenance on CL relative to 1A (unitless)")                  # Kloos 2021 Table 2 final model (MRG maint. 0.81, RSE 3.9%; bootstrap 0.81, 95% CI 0.75-0.86)
    e_phasesriv_cl    <- 0.81;     label("Multiplicative effect of standard-risk protocol IV on CL relative to 1A (unitless)")                # Kloos 2021 Table 2 final model (SRG protocol IV 0.81, RSE 6.5%; bootstrap 0.81, 95% CI 0.73-0.90)

    # ---- Shared eta between CL and Vd ------------------------------------
    # IIV on Vd could not be estimated independently because it correlated
    # completely with the CL-Vd correlation (Kloos 2021 Online Supplementary
    # Appendix, Supplemental Results, PK analysis). Instead the supplement's
    # structural equation 'Vd = Theta3 * e^(Theta4 * eta) * BSA' scales the CL
    # eta onto Vd by Theta4, printed in Table 2 as 'Fractional increase
    # IIV-Vd'. Same scaled-shared-eta encoding as Hirt_2009_efavirenz.R
    # (vc_eta_scale), Prytula_2016_tacrolimus.R (q_eta_scale),
    # Sharma_2016_hydroxyprogesteroneCaproate.R and
    # Leshinsky_2017_caspofungin_cat.R (scale_etalvc).
    scale_etalvc <- 1.26; label("Scaling of the CL inter-individual eta onto Vd (unitless)")  # Kloos 2021 Table 2 final model ('Fractional increase IIV-Vd' 1.26, RSE 10.2%; bootstrap 1.28, 95% CI 0.98-1.51)

    # ---- Inter-individual variability ------------------------------------
    # Table 2 reports IIV on CL as 19.7% CV. Exponential IIV, so
    # omega^2 = log(1 + CV^2) = log(1 + 0.197^2) = 0.038077.
    etalcl ~ 0.038077  # Kloos 2021 Table 2 final model (IIV CL 19.7% CV, RSE 13.2%; bootstrap 19.5%, 95% CI 14.6-24.0) converted as log(1 + 0.197^2)

    # ---- Inter-occasion variability --------------------------------------
    # Table 2 reports IOV on CL as 23.6% CV; pi^2 = log(1 + 0.236^2) = 0.054197.
    # The supplement defines an occasion as the administration of a new dose.
    # nlmixr2 has no NONMEM `$OMEGA BLOCK(1) SAME` shortcut, so occasions 2 and
    # 3 get their own etas with the variance fixed equal to the occasion-1
    # estimate (Jonsson_2011_ethambutol.R / Chen_2023_nemonoxacin.R precedent).
    etaiov_cl_1 ~ 0.054197       # Kloos 2021 Table 2 final model (IOV CL 23.6% CV, RSE 10.3%; bootstrap 23.1%, 95% CI 18.5-27.5) converted as log(1 + 0.236^2)
    etaiov_cl_2 ~ fix(0.054197)  # SAME-equivalent: equal to the occasion-1 IOV variance
    etaiov_cl_3 ~ fix(0.054197)  # SAME-equivalent: equal to the occasion-1 IOV variance

    # ---- Residual error ---------------------------------------------------
    # Kloos 2021 Online Supplementary Appendix, Supplemental Results, PK
    # analysis: 'a combined model of proportional and additive error was
    # superior'. Table 2 final model reports the additive term as 20.2 IU/L
    # and the proportional term as 17.0%.
    propSd <- 0.170; label("Proportional residual error (fraction)")  # Kloos 2021 Table 2 final model (Proportional 17.0%, RSE 9.2%; bootstrap 16.8%, 95% CI 14.2-19.5)
    addSd  <- 20.2;  label("Additive residual error (IU/L)")          # Kloos 2021 Table 2 final model (Additional 20.2 IU/L, RSE 18.2%; bootstrap 20.2, 95% CI 12.7-298)
  })

  model({
    # ---- Inter-occasion variability on clearance -------------------------
    # One occasion per administered dose (Supplemental Methods). Records with
    # OCC outside 1-3 contribute no IOV.
    oc1    <- (OCC == 1)
    oc2    <- (OCC == 2)
    oc3    <- (OCC == 3)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3

    # ---- Time-dependent clearance ----------------------------------------
    # Before the split point the induced arm is switched off (cl_ind = 1);
    # after it, CL rises linearly at exp(lcl_time) per day. ADDITIVE reading
    # of the slope -- see the ini() note above.
    #
    # Written in ratio form so the population value stays mu-referenced
    # (nlmixr2 requires the eta to sit inside exp() alongside its theta). This
    # is algebraically identical to the additive form:
    #   CL_per_m2 = Theta1 + Theta2 * (TAD - split)
    #             = Theta1 * (1 + (Theta2 / Theta1) * (TAD - split))
    # so cl_ind = 1 + (0.082 / 0.084) * (TAD - 12.7) = 1 + 0.9762 * (TAD - 12.7),
    # which is exactly the supplement's printed 'Cl_ind = 1 + Theta2 * (TAD -
    # split point)' with the missing / Theta1 restored.
    #
    # `split_point <- tsplit` and `slope_ratio <- exp(lcl_time - lcl)` are
    # deliberately given their own simple lines. rxode2's mu-reference parser
    # fails with "mu-ref err: subscript out of bounds" when an ini() population
    # parameter is used in a non-mu-referenced position inside a larger
    # expression while any eta is present; assigning the parameter to a plain
    # local first is the documented work-around ("try putting the
    # mu-referenced expression on a simple line").
    split_point <- tsplit
    slope_ratio <- exp(lcl_time - lcl)

    # ---- Covariate factors on clearance ----------------------------------
    # Written exactly as Kloos 2021 Table 2: <estimate>^<binary indicator>.
    # Protocol 1A (all treatment-phase indicators 0, no infection) is the
    # reference and collapses every factor to 1.
    infect_factor <- e_infect_cl^DIS_INFECT_ACTIVE
    phase_factor  <- e_phase1b_cl^TRTPH_1B *
      e_phasem_cl^TRTPH_M *
      e_phasemrint_cl^TRTPH_MR_INTENS *
      e_phasemrmaint_cl^TRTPH_MR_MAINT *
      e_phasesriv_cl^TRTPH_SR_IV

    # ---- Individual parameters -------------------------------------------
    # Pre-split clearance (L/day) = Theta1 * BSA * covariate factors *
    # exp(IIV + IOV). The etas multiply the whole clearance, including the
    # induced arm (supplement: 'Cl (after 13 days) = Theta1 * e^(eta +
    # eta_IOV) * BSA * Cl_ind').
    cl0 <- exp(lcl + etalcl + iov_cl) * BSA * infect_factor * phase_factor

    # Vd (L) = Theta3 * exp(Theta4 * eta_CL) * BSA. The CL eta is shared and
    # scaled rather than an independent eta (see the ini() note).
    vc <- exp(lvc + scale_etalvc * etalcl) * BSA

    # One-compartment model with intravenous input; each dose is infused into
    # the central compartment over 1 hour (Kloos 2021 Methods, Patients and
    # treatment protocol). The infusion duration is supplied by the event
    # table (rate / dur), not by the model.
    #
    # The time-after-dose term is written INLINE in the d/dt right-hand side,
    # as `t - tlast()`, rather than being pre-computed into a local via
    # `tad()`. This is load-bearing, not stylistic. rxode2 evaluates the
    # standalone algebraic assignments of model() ONCE PER OUTPUT RECORD and
    # holds them constant across the integration interval, so a clearance
    # built from `tad()` in a local variable is a step function of the
    # observation grid: with a single trough observation 14 days after an
    # 800 IU/m^2 dose the solved trough is 50 IU/L, versus 227 IU/L on a dense
    # grid -- a 4.5-fold error that shrinks only as the grid is refined. Terms
    # appearing directly inside `d/dt(...)` are instead evaluated at every
    # solver step, so `t` advances continuously while `tlast()` (the time of
    # the most recent dose) correctly stays constant within a dosing interval.
    # Written this way the solution is exact and independent of the output
    # grid. Do not "simplify" this back into a `tad()` local.
    d/dt(central) <- -(cl0 / vc) * central *
      (1 + slope_ratio * (t - tlast() - split_point) *
         ((t - tlast()) > split_point))

    # Record-level reporting copies of the induced clearance, the total
    # clearance and the elimination rate constant. These are evaluated at each
    # output record (where `tad()` is exact) and are returned as solve columns
    # for diagnostics -- e.g. the instantaneous half-life ln(2) * vc / cl. They
    # do NOT feed the ODE above.
    tad_split <- tad() - split_point
    cl_ind    <- 1 + slope_ratio * tad_split * (tad_split > 0)
    cl        <- cl0 * cl_ind
    kel       <- cl / vc

    # Asparaginase activity level: dose in IU, vc in L -> IU/L, matching the
    # AHA assay output.
    Cc <- central / vc

    # Combined proportional + additive residual error.
    Cc ~ prop(propSd) + add(addSd)
  })
}
