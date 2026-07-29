Han_2015_decitabine <- function() {
  description <- paste(
    "Two-compartment IV population pharmacokinetic model coupled with two",
    "parallel Friberg-style myelosuppression PD chains (absolute neutrophil",
    "count, ANC, and platelet count, PC) for decitabine post-transplant",
    "maintenance in adult Korean patients with higher-risk myelodysplastic",
    "syndrome or secondary acute myeloid leukemia (Han 2015). The platelet",
    "feedback baseline rises asymptotically over cycles per the paper's",
    "IMP / IMK extension (BASE_P_t = BASE_P + IMP * (1 - exp(-IMK * t)));",
    "the neutrophil chain uses a time-invariant baseline. PK parameters",
    "are body-surface-area-normalized (per m^2): doses must be supplied",
    "in mg/m^2 and central-compartment concentrations are returned in",
    "mg/L (= ug/mL = 1000 ng/mL). PD outputs ANC and PLT are in",
    "10^9 cells/L."
  )
  reference <- paste(
    "Han S, Kim Y-J, Lee J, Jeon S, Hong T, Park G-J, Yoon J-H, Yahng S-A,",
    "Shin S-H, Lee S-E, Eom K-S, Kim H-J, Min C-K, Lee S, Yim D-S. (2015).",
    "Model-based adaptive phase I trial design of post-transplant",
    "decitabine maintenance in myelodysplastic syndrome.",
    "Journal of Hematology & Oncology 8:118.",
    "doi:10.1186/s13045-015-0208-3 (PMID 26482429).",
    "ClinicalTrials.gov NCT01277484."
  )
  vignette <- "Han_2015_decitabine"
  units <- list(
    time          = "hour",
    dosing        = "mg/m^2",
    concentration = "mg/L",
    anc           = "10^9 cells/L",
    plt           = "10^9 cells/L"
  )

  # No model-referenced covariates. PK parameters are reported on a per-m^2
  # body-surface-area-normalized scale (Han 2015 Table 3), so the dataset
  # is expected to carry doses already converted to mg/m^2 (the protocol
  # records doses as 4-11 mg/m^2/day, 60-minute IV infusion); BSA is not
  # used inside model(). Covariate analysis on sex, age, baseline body
  # weight, body-surface area, and clinical-lab variables found no
  # meaningful covariate effect (Methods + Results: "No meaningful
  # covariate was found in either the patient demographic or clinical
  # variables").
  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 15L,
    n_studies      = 1L,
    age_range      = "19-64 years",
    weight_range   = NA_character_,
    sex_female_pct = 40,
    race_ethnicity = c(Asian = 100),
    disease_state  = paste0(
      "Adult Korean patients with higher-risk myelodysplastic syndrome ",
      "(MDS; intermediate-2 or high IPSS risk) or secondary acute myeloid ",
      "leukemia (sAML) evolving from MDS, receiving decitabine maintenance ",
      "after allogeneic hematopoietic stem cell transplantation (allo-HSCT). ",
      "Eligibility: age <= 65; PC > 30,000/mm^3 and ANC > 1000/mm^3 ",
      "maintained > 7 days without transfusion or growth factors; no ",
      "grade III/IV acute GVHD; ECOG 0-2; no renal or hepatic impairment. ",
      "Decitabine started on days 42-90 post-transplant (median 86 days)."
    ),
    dose_range     = paste0(
      "Decitabine 4-12 mg/m^2/day (Cycle 2-4 range 1.5-12 mg/m^2/day from ",
      "Table 2) as a 60-minute IV infusion for 5 consecutive days, cycles ",
      "repeated every 4 weeks up to Cycle 12. Initial Cycle-1 dose ",
      "5 mg/m^2/day for cohort 1; subsequent cohort initial doses ",
      "estimated by cohort dose estimation (CDE) at 4, 5, 5.5, and ",
      "5 mg/m^2/day for cohorts 2-5. Doses for Cycles 2-4 individualized ",
      "via PK-PD adaptive titration (IDT); the Cycle-4 dose was maintained ",
      "for all subsequent cycles. The actual dosing interval was 34.5 ",
      "+/- 8.7 days (mean +/- SD)."
    ),
    regions        = "Republic of Korea (Seoul St. Mary's Hospital, The Catholic University of Korea).",
    co_medication  = paste0(
      "GVHD prophylaxis: calcineurin inhibitor (cyclosporine for related ",
      "donors, tacrolimus for unrelated donors) plus short-course ",
      "methotrexate. Antithymocyte globulin given to all patients prior ",
      "to transplant. At decitabine initiation: acute GVHD grade 0-2 in ",
      "all patients, mild chronic GVHD in 1 patient."
    ),
    notes          = paste0(
      "Baseline demographics from Han 2015 Tables 1 and 2. 16 patients ",
      "were enrolled (9 male / 7 female); the popPK-PD analysis excluded ",
      "1 female patient whose dose-limiting toxicity factor was platelet ",
      "count due to immune thrombocytopenia after transplantation (managed ",
      "with steroids), leaving 15 patients in the mixed-effect analysis ",
      "(8 male / 6 female assumed; sex_female_pct = 40 reflects 6/15). ",
      "Donor type for the 16 enrolled: matched sibling donor (n=7), ",
      "matched unrelated donor (n=8), partially matched unrelated donor ",
      "(n=1). WHO diagnoses: RAEB-1 (n=1), RAEB-2 (n=10), AML evolving ",
      "from MDS (n=5). Final analysis dataset: 95 PK observations and ",
      "622 PD observations (311 ANC + 311 PC) across the 15 subjects ",
      "(Han 2015 Results, 'Overall mixed-effect PK-PD analysis'). PK ",
      "proportional residual error variance sigma^2_PK = 0.441 (= 66% CV); ",
      "PD additive residual error variances sigma^2_PD,P = 25000 /mm^3 ",
      "squared and sigma^2_PD,N = 754 /mm^3 squared (Table 3). Successful ",
      "bootstrap convergence proportion was 78.8% (PK) and 78.0% (PD)."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # PK structural parameters - Han 2015 Table 3, "Pharmacokinetic
    # parameters" block. All values are body-surface-area-normalized
    # (per m^2): CL in L/h per m^2, volumes in L per m^2. With doses
    # supplied in mg/m^2 (BSA-normalized amount), the central-compartment
    # state holds mg/m^2 and `central / vc` therefore yields mg/L. NOTE:
    # the paper's Discussion claims a predicted Cmax of 66.0 ng/mL after
    # a 1-h infusion of 5 mg/m^2, but the published Table 3 parameters
    # imply ~51 ng/mL (asymptotic infusion concentration k0/CL = 5/87.8
    # = 57 ng/mL ceiling). The validation vignette reports both values
    # and flags the paper-internal discrepancy.
    # ---------------------------------------------------------------------
    lcl <- log(87.8); label("CL: decitabine clearance (L/h per m^2)")                          # Table 3
    lvc <- log(18.5); label("Vc: central compartment volume (L per m^2)")                       # Table 3
    lvp <- log(22.9); label("Vp: peripheral compartment volume (L per m^2)")                    # Table 3
    lq  <- log(13.1); label("Q: intercompartmental clearance (L/h per m^2)")                    # Table 3

    # ---------------------------------------------------------------------
    # PD - neutrophil chain - Han 2015 Table 3, "Pharmacodynamic
    # parameters for neutrophil" block. Friberg-style transit-compartment
    # model with feedback (Wallin 2009 reference; standard 1 proliferating
    # + 3 transit + 1 circulating chain with ktr applied uniformly across
    # proliferation, transits, and circulating-pool clearance). The
    # neutrophil baseline is time-invariant (no asymptotic recovery
    # structure for ANC: the paper notes "baseline cell count increase
    # was not meaningful" for neutrophils).
    # ---------------------------------------------------------------------
    lktr_anc   <- log(0.0132);   label("ktr_N: rate constant of inter-compartmental neutrophil movement (1/h)")  # Table 3
    lslope_anc <- log(0.263);    label("SLOPE_N: linear drug-effect slope on ANC proliferation (L/mg)")          # Table 3
    lbase_anc  <- log(3.24);     label("BASE_N: baseline neutrophil count (10^9 cells/L)")                       # Table 3 (3240/mm^3 = 3.24 * 10^9/L)
    lgamma_anc <- log(0.193);    label("GAMMA_N: feedback exponent on (BASE_N / circ_anc) (unitless)")           # Table 3

    # ---------------------------------------------------------------------
    # PD - platelet chain - Han 2015 Table 3, "Pharmacodynamic parameters
    # for platelet" block. Same Friberg structure as the ANC chain but
    # with a time-dependent baseline that increases asymptotically over
    # cycles: BASE_P_t = BASE_P + IMP * (1 - exp(-IMK * TIME)). The
    # time-dependent BASE_P_t substitutes for BASE_P everywhere it
    # appears in the Friberg feedback term (Han 2015 Results, "Overall
    # mixed-effect PK-PD analysis").
    # ---------------------------------------------------------------------
    lktr_plt   <- log(0.0244);    label("ktr_P: rate constant of inter-compartmental platelet movement (1/h)")    # Table 3
    lslope_plt <- log(0.0656);    label("SLOPE_P: linear drug-effect slope on PC proliferation (L/mg)")           # Table 3
    lbase_plt  <- log(49.2);      label("BASE_P: baseline platelet count at t=0 (10^9 cells/L)")                  # Table 3 (49200/mm^3 = 49.2 * 10^9/L)
    lgamma_plt <- log(0.304);     label("GAMMA_P: feedback exponent on (BASE_P_t / circ_plt) (unitless)")         # Table 3
    limp       <- log(55.0);      label("IMP: maximum platelet count recovery expected over cycles (10^9 cells/L)") # Table 3 (55000/mm^3)
    limk       <- log(0.000530);  label("IMK: rate constant for asymptotic platelet count recovery (1/h)")         # Table 3

    # ---------------------------------------------------------------------
    # IIV - Han 2015 Table 3, "Between-subject variability" columns
    # (estimates reported as CV%). Converted to log-eta variance via
    # omega^2 = log(CV^2 + 1) for log-normal individual parameters. NE
    # (not estimated) entries in Table 3 for Vc, Vp, Q, ktr_N, ktr_P,
    # GAMMA_P, and IMK are encoded by simply omitting the corresponding
    # etas (rather than fixing to 0), so the typical-value simulation
    # reproduces these as fixed structural parameters with no individual
    # spread.
    # ---------------------------------------------------------------------
    etalcl        ~ 0.04478;  label("IIV variance on log-CL (Table 3: 21.4% CV)")            # log(0.214^2 + 1) = 0.04478
    etalslope_anc ~ 0.28565;  label("IIV variance on log-SLOPE_N (Table 3: 57.5% CV)")       # log(0.575^2 + 1) = 0.28565
    etalbase_anc  ~ 0.17326;  label("IIV variance on log-BASE_N (Table 3: 43.5% CV)")        # log(0.435^2 + 1) = 0.17326
    etalgamma_anc ~ 0.14430;  label("IIV variance on log-GAMMA_N (Table 3: 39.4% CV)")       # log(0.394^2 + 1) = 0.14430
    etalslope_plt ~ 0.06396;  label("IIV variance on log-SLOPE_P (Table 3: 25.7% CV)")       # log(0.257^2 + 1) = 0.06396
    etalbase_plt  ~ 0.74320;  label("IIV variance on log-BASE_P (Table 3: 105% CV)")         # log(1.05^2 + 1)  = 0.74320
    etalimp       ~ 0.47798;  label("IIV variance on log-IMP (Table 3: 78.3% CV)")           # log(0.783^2 + 1) = 0.47798

    # ---------------------------------------------------------------------
    # Residual error - Han 2015 Table 3, "Residual error" block. Reported
    # as variances on the native concentration / cell-count scales.
    #   sigma^2_PK   = 0.441   -> propSd    = sqrt(0.441)   = 0.664 (proportional, fraction of Cc)
    #   sigma^2_PD,P = 25000 (/mm^3)^2 -> addSd_PLT = sqrt(25000) /mm^3 = 158.1 /mm^3
    #                                  = 0.158 * 10^9/L (canonical unit)
    #   sigma^2_PD,N = 754   (/mm^3)^2 -> addSd_ANC = sqrt(754) /mm^3   = 27.46 /mm^3
    #                                  = 0.0275 * 10^9/L (canonical unit)
    # ---------------------------------------------------------------------
    propSd    <- 0.664;  label("Proportional residual SD on plasma decitabine concentration (fraction)")  # Table 3, sqrt(0.441)
    addSd_ANC <- 0.0275; label("Additive residual SD on circulating neutrophil count (10^9 cells/L)")     # Table 3, sqrt(754)  /mm^3 / 1000
    addSd_PLT <- 0.158;  label("Additive residual SD on circulating platelet count (10^9 cells/L)")        # Table 3, sqrt(25000) /mm^3 / 1000
  })

  model({
    # --------------------------------------------------------------------
    # Individual PK parameters - log-normal IIV applied to CL only (the
    # only PK BSV that Han 2015 could estimate; Vc, Vp, Q have NE in
    # Table 3 BSV column). All parameters remain on the per-m^2 scale.
    # --------------------------------------------------------------------
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    vp <- exp(lvp)
    q  <- exp(lq)

    # PK micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # --------------------------------------------------------------------
    # Individual PD parameters - ANC chain
    # --------------------------------------------------------------------
    ktr_anc   <- exp(lktr_anc)
    slope_anc <- exp(lslope_anc + etalslope_anc)
    base_anc  <- exp(lbase_anc  + etalbase_anc)
    gamma_anc <- exp(lgamma_anc + etalgamma_anc)

    # --------------------------------------------------------------------
    # Individual PD parameters - PLT chain (with time-varying baseline)
    # --------------------------------------------------------------------
    ktr_plt   <- exp(lktr_plt)
    slope_plt <- exp(lslope_plt + etalslope_plt)
    base_plt  <- exp(lbase_plt  + etalbase_plt)
    gamma_plt <- exp(lgamma_plt)
    imp       <- exp(limp + etalimp)
    imk       <- exp(limk)

    # Time-dependent platelet baseline (Han 2015 second formula in
    # "Overall mixed-effect PK-PD analysis": BASE_P_t = BASE_P +
    # IMP * (1 - exp(-IMK * TIME))). Uses the rxode2 reserved variable
    # `t` for current simulation time. At t = 0, baseline_plt_t reduces
    # to base_plt, which is consistent with the per-subject initial
    # condition below.
    baseline_plt_t <- base_plt + imp * (1 - exp(-imk * t))

    # --------------------------------------------------------------------
    # Concentration in central compartment. Doses are mg/m^2 (BSA-
    # normalized); vc is L per m^2; the m^2 cancels and Cc is returned
    # in mg/L (= ug/mL).
    # --------------------------------------------------------------------
    Cc <- central / vc

    # Drug effect and feedback (Friberg form: linear drug effect on
    # proliferation, power-law feedback from the circulating pool toward
    # the baseline).
    edrug_anc <- 1 - slope_anc * Cc
    feed_anc  <- (base_anc / circ_anc)^gamma_anc
    edrug_plt <- 1 - slope_plt * Cc
    feed_plt  <- (baseline_plt_t / circ_plt)^gamma_plt

    # --------------------------------------------------------------------
    # ODE system
    # --------------------------------------------------------------------
    # PK - 2-compartment IV
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Friberg ANC chain: proliferation + 3 transits + circulating
    d/dt(precursor1_anc) <- ktr_anc * precursor1_anc * edrug_anc * feed_anc - ktr_anc * precursor1_anc
    d/dt(precursor2_anc) <- ktr_anc * precursor1_anc - ktr_anc * precursor2_anc
    d/dt(precursor3_anc) <- ktr_anc * precursor2_anc - ktr_anc * precursor3_anc
    d/dt(precursor4_anc) <- ktr_anc * precursor3_anc - ktr_anc * precursor4_anc
    d/dt(circ_anc)       <- ktr_anc * precursor4_anc - ktr_anc * circ_anc

    # Friberg PLT chain: proliferation + 3 transits + circulating, time-
    # varying baseline (baseline_plt_t enters via feed_plt only; the
    # transit kinetics themselves use the structural ktr_plt).
    d/dt(precursor1_plt) <- ktr_plt * precursor1_plt * edrug_plt * feed_plt - ktr_plt * precursor1_plt
    d/dt(precursor2_plt) <- ktr_plt * precursor1_plt - ktr_plt * precursor2_plt
    d/dt(precursor3_plt) <- ktr_plt * precursor2_plt - ktr_plt * precursor3_plt
    d/dt(precursor4_plt) <- ktr_plt * precursor3_plt - ktr_plt * precursor4_plt
    d/dt(circ_plt)       <- ktr_plt * precursor4_plt - ktr_plt * circ_plt

    # --------------------------------------------------------------------
    # Initial conditions - both Friberg chains at steady state at t = 0
    # (no drug before treatment initiation). The PLT chain uses base_plt
    # because baseline_plt_t(0) = base_plt + IMP * (1 - exp(0)) =
    # base_plt.
    # --------------------------------------------------------------------
    precursor1_anc(0) <- base_anc
    precursor2_anc(0) <- base_anc
    precursor3_anc(0) <- base_anc
    precursor4_anc(0) <- base_anc
    circ_anc(0)       <- base_anc

    precursor1_plt(0) <- base_plt
    precursor2_plt(0) <- base_plt
    precursor3_plt(0) <- base_plt
    precursor4_plt(0) <- base_plt
    circ_plt(0)       <- base_plt

    # --------------------------------------------------------------------
    # Observations
    # --------------------------------------------------------------------
    ANC <- circ_anc
    PLT <- circ_plt

    Cc  ~ prop(propSd)
    ANC ~ add(addSd_ANC)
    PLT ~ add(addSd_PLT)
  })
}
