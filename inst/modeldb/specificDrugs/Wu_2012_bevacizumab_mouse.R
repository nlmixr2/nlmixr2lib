Wu_2012_bevacizumab_mouse <- function() {
  description <- paste(
    "Preclinical (mouse).",
    "Mechanistic two-compartment population PK model for bevacizumab",
    "labelled with the near-infrared dye IRDye 800CW, fitted",
    "simultaneously to plasma profiles after a single 0.45 mg/kg",
    "intravenous (penile vein) and subcutaneous (front footpad) dose and",
    "to the draining axillary lymph-node concentrations that follow the",
    "subcutaneous dose, in male SKH-1 mice. The subcutaneous input is",
    "split at the injection site between a direct blood-capillary route",
    "and a lymphatic route through the draining node: a fraction flnode",
    "(the paper's Frc) of the absorbed dose passes through the lymph-node",
    "compartment and returns to plasma with rate constant",
    "k_lnode_central (the paper's ka2), while the complementary fraction",
    "1 - flnode enters plasma directly. Disposition is parameterised in",
    "micro-constant form exactly as the authors report it (Vc, k10, k12,",
    "k21) rather than as clearances. Every dose, volume and rate in the",
    "source is normalised to body weight, so the model is coded on a",
    "per-kilogram basis: state amounts are ug/kg and volumes are mL/kg,",
    "which puts central/vc and lnode/v_lnode directly in ug/mL.",
    "The model was fitted by naive pooling of Bailer-method sacrificial",
    "sampling means, so the source reports neither between-subject",
    "variability nor a residual-error model; both residual SDs are",
    "therefore fixed at zero and the model is deterministic.",
    sep = " ")
  reference <- paste(
    "Wu F, Tamhane M, Morris ME.",
    "Pharmacokinetics, Lymph Node Uptake, and Mechanistic PK Model of",
    "Near-Infrared Dye-Labeled Bevacizumab After IV and SC Administration",
    "in Mice.",
    "AAPS J. 2012;14(2):252-261.",
    "doi:10.1208/s12248-012-9342-9.",
    sep = " ")
  vignette <- "Wu_2012_bevacizumab_mouse"
  units <- list(
    time          = "h",
    dosing        = "ug",
    concentration = "ug/mL"
  )

  # Per-kg normalisation: Table II reports Vc and VLN in mL/kg and the dose
  # is 0.45 mg/kg, so state amounts are ug/kg and amount / volume is
  # ug/mL -- the units of Table I (Cmax in ug/mL) and Figs. 6 and 7.
  compartmentData <- list(
    depot       = list(analyte = "bevacizumab", units = "ug", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "bevacizumab", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "bevacizumab", units = "ug", specimen = "tissue", verified = TRUE),
    # Axillary lymph node draining the front-footpad SC injection site.
    # Assayed as node homogenate by ELISA, hence "tissue" rather than
    # "lymph" (which is the fluid): Methods, "Plasma and Lymph Node
    # Concentrations Determined by ELISA".
    lnode       = list(analyte = "bevacizumab", units = "ug", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species        = "mouse (SKH-1, male)",
    n_subjects     = 63L,
    n_studies      = 1L,
    weight_range   = "25-30 g",
    sex_female_pct = 0,
    disease_state  = "Healthy (no disease model; bevacizumab does not bind murine VEGF-A)",
    dose_range     = paste(
      "Single 0.45 mg/kg dose of bevacizumab conjugated to IRDye 800CW",
      "(dye:protein ratio 3.5:1; 10.5 nmol/kg of dye), given either",
      "intravenously by penile-vein injection or subcutaneously into the",
      "front footpad. The dose was capped by the 15 uL maximum footpad",
      "injection volume."
    ),
    regions        = "USA (University at Buffalo)",
    notes          = paste(
      "Destructive sacrificial sampling: three animals were killed at each",
      "of 5 min (IV only), 15 min, 30 min, 1, 2, 4, 8, 24, 72, 168 and",
      "288 h, giving 11 IV and 10 SC timepoints (n = 63 animals in total).",
      "Bevacizumab was quantified in plasma and in axillary lymph-node",
      "homogenate by a validated human-IgG ELISA (assay range 0.82-200",
      "ng/mL); plasma was additionally quantified to 24 h by near-infrared",
      "fluorescence imaging (LLOQ 0.5 ug/mL). All concentrations used for",
      "the model fit were the ELISA ones. Node concentrations were",
      "quantifiable only through 8 h, and at 8 h in only one of three",
      "animals. Lymph-node uptake after the IV dose was negligible and is",
      "not modelled. See Wu 2012 Methods, 'Animal Studies' and",
      "'Compartmental Analysis and Modeling'."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural disposition -- Wu 2012 Table II ("PK Parameter Estimates
    # for Bevacizumab 800CW After IV and SC Administration (0.45 mg/kg) to
    # Mice"), Mean column. The CV% column is the relative standard error
    # of each estimate from the naive-pooled fit, NOT a variance
    # component, so no IIV is encoded (see the IIV block below).
    #
    # The paper parameterises disposition on micro-constants, so k10 is
    # carried as the canonical kel and Vc as lvc; no clearance or
    # inter-compartmental clearance is reported. Derived: CL = kel * vc =
    # 0.0162 * 319 = 5.17 mL/h/kg.
    # ------------------------------------------------------------------
    lvc  <- log(319)    ; label("Central volume of distribution Vc (log mL/kg)")                  # Table II: Vc = 319 mL/kg (CV 9.59%)
    lkel <- log(0.0162) ; label("Elimination rate constant from central k10 (log 1/h)")           # Table II: k10 = 0.0162 /h (CV 24.5%)
    lk12 <- log(0.0314) ; label("Rate constant central -> peripheral1 k12 (log 1/h)")             # Table II: k12 = 0.0314 /h (CV 40.1%)
    lk21 <- log(0.0131) ; label("Rate constant peripheral1 -> central k21 (log 1/h)")             # Table II: k21 = 0.0131 /h (CV 57.0%)

    # ------------------------------------------------------------------
    # Subcutaneous absorption and the lymphatic limb -- Wu 2012 Table II
    # and Eqs. (3), (5), (6).
    #
    # ka1 is the single first-order rate constant emptying the SC site;
    # Fig. 1's legend states it "includes two absorption processes,
    # systemic absorption and lymphatic absorption", which the model then
    # splits by Frc. It is therefore the canonical lka.
    #
    # NAMING. The lymph-node limb had no canonical coverage before this
    # extraction. The four names below were ratified by the operator on
    # 2026-09-02 (sidecar oare_PMC3326166 request-001 / response-001,
    # questions q1-q4, all option A):
    #   Frc -> lflnode  (fraction routed to the `lnode` compartment, named
    #                    by destination exactly as fdepot names the
    #                    fraction routed to `depot`; the register's
    #                    f_lymph is a DIFFERENT quantity, lymph flow as a
    #                    fraction of plasma flow)
    #   VLN -> lv_lnode (member of the lv_<space> family founded by lv_elf)
    #   ka2 -> lk_lnode_central (member of the k_<from>_<to> directional-
    #                    transfer family; the paper calls it an absorption
    #                    rate constant but structurally it is a transfer
    #                    between two modelled states)
    #   A_LN,sc -> the already-registered `lnode` compartment, whose role
    #                    text is broadened by this extraction.
    # ------------------------------------------------------------------
    lka             <- log(4.82)    ; label("First-order absorption rate constant from the SC site ka1 (log 1/h)")            # Table II: ka1 = 4.82 /h (CV 40.9%)
    lflnode         <- log(0.00964) ; label("Fraction of the absorbed SC dose routed through the lymph node Frc (log)")       # Table II: Frc = 0.00964 (CV 19.6%)
    lk_lnode_central <- log(0.723)  ; label("Rate constant lymph node -> central ka2 (log 1/h)")                              # Table II: ka2 = 0.723 /h (CV 9.64%)

    # ------------------------------------------------------------------
    # Fixed parameters -- Wu 2012 Methods, "Compartmental Analysis and
    # Modeling": "BIO was fixed as 1 based on the NCA results, and volume
    # of distribution for the lymph node compartment, VLN, was fixed as
    # 0.33 mL/kg based on the actual weight of axillary lymph nodes
    # (assuming 1 g/ml specific density) and the average weight of SKH-1
    # mice." Table II marks both "(Fixed)".
    #
    # BIO is the paper's bioavailability and is carried as the canonical
    # lfdepot. Note that Eqs. (3), (5) and (6) place BIO on the depot RATE
    # rather than on the dose: it multiplies the depot efflux term and
    # both inflow terms identically, so the system is mass-conserving for
    # any BIO and BIO scales the absorption RATE rather than the absorbed
    # FRACTION. The equations are encoded exactly as printed. This is
    # immaterial here because BIO is fixed at 1; see the vignette
    # "Assumptions and deviations" section.
    # ------------------------------------------------------------------
    lv_lnode <- fixed(log(0.33)) ; label("Lymph-node volume of distribution VLN (log mL/kg)")  # Table II: VLN = 0.33 mL/kg (Fixed); Methods: axillary-node weight at 1 g/mL over mean SKH-1 body weight
    lfdepot  <- fixed(log(1))    ; label("Bioavailability of the SC dose BIO (log fraction)")  # Table II: BIO = 1 (Fixed); Methods: fixed from the NCA bioavailability estimate

    # ------------------------------------------------------------------
    # IIV -- none. The plasma and lymph-node data were obtained by
    # destructive sacrificial sampling (three animals per timepoint) and
    # fitted by naive pooling of the Bailer-method means, so no
    # between-subject variance component is estimable and none is
    # reported. The Table II CV% column is the relative standard error of
    # each point estimate. No eta terms are encoded rather than inventing
    # variances.
    #
    # Residual error -- not reported. The source states no residual-error
    # model and Table II carries no sigma. Both residual SDs are fixed at
    # zero so the endpoint structure is declared without inventing an
    # error magnitude; the model is deterministic.
    # ------------------------------------------------------------------
    propSd         <- fixed(0) ; label("Plasma proportional residual SD (fraction; ZERO - not reported in source)")      # Wu 2012 reports no residual-error model (naive-pooled Bailer-method fit)
    propSd_Clnode  <- fixed(0) ; label("Lymph-node proportional residual SD (fraction; ZERO - not reported in source)")  # Wu 2012 reports no residual-error model (naive-pooled Bailer-method fit)
  })

  model({
    # 1. Individual parameters. No IIV (see ini()).
    vc              <- exp(lvc)
    kel             <- exp(lkel)
    k12             <- exp(lk12)
    k21             <- exp(lk21)
    ka              <- exp(lka)
    flnode          <- exp(lflnode)
    k_lnode_central <- exp(lk_lnode_central)
    v_lnode         <- exp(lv_lnode)
    fdepot          <- exp(lfdepot)

    # 2. ODE system -- Wu 2012 Eqs. (1)-(6), transcribed verbatim.
    #
    # The paper writes two parallel sets of equations, one for the IV arm
    # (Eqs. 1-2) and one for the SC arm (Eqs. 3-6), with the same
    # disposition constants in both. They are the same system: the IV
    # equations are the SC equations with an empty depot. The single
    # system below therefore reproduces both arms -- dose `central` for
    # the IV arm, dose `depot` for the SC arm.
    #
    # Eq. (6)  dAsc/dt    = -BIO * ka1 * Asc
    # Eq. (3)  dAs,sc/dt  = -(k12 + k10) * As,sc + k21 * At,sc
    #                       + ka2 * ALN,sc + BIO * (1 - Frc) * ka1 * Asc
    # Eq. (4)  dAt,sc/dt  =  k12 * As,sc - k21 * At,sc
    # Eq. (5)  dALN,sc/dt =  BIO * Frc * ka1 * Asc - ka2 * ALN,sc
    #
    # (The PDF's display equations are set in a symbol font that
    # substitutes glyphs for the operators; they recover exactly as above
    # from `pdftotext -layout`. See the vignette source-trace table.)
    d/dt(depot)       <- -fdepot * ka * depot
    d/dt(central)     <- -(k12 + kel) * central + k21 * peripheral1 +
      k_lnode_central * lnode + fdepot * (1 - flnode) * ka * depot
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(lnode)       <- fdepot * flnode * ka * depot - k_lnode_central * lnode

    # 3. Observations -- Wu 2012 Eq. (7): As = Cs * Vc and
    #    ALN,sc = CLN,sc * VLN. Amount in ug/kg over volume in mL/kg gives
    #    ug/mL, the units of Table I and Figs. 6-7.
    Cc     <- central / vc
    Clnode <- lnode / v_lnode

    Cc     ~ prop(propSd)
    Clnode ~ prop(propSd_Clnode)
  })
}
