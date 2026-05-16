FehlingKaschek_2019_trastuzumab_skbr3 <- function() {
  description <- "In vitro (SKBR3 cell line). Mechanistic ODE model of trastuzumab-induced HER2 receptor internalization with two cell-membrane phenotypes (ruffled vs flat); Model B of Fehling-Kaschek 2019, no recycling or degradation."
  reference <- "Fehling-Kaschek M, Peckys DB, Kaschek D, Timmer J, de Jonge N. Mathematical modeling of drug-induced receptor internalization in the HER2-positive SKBR3 breast cancer cell-line. Sci Rep. 2019 Sep 5;9(1):12709. doi:10.1038/s41598-019-49019-x. PMID: 31488874; PMCID: PMC6728336."
  vignette <- "FehlingKaschek_2019_trastuzumab_skbr3"
  units <- list(
    time = "min",
    dosing = "dimensionless",
    concentration = "fraction of initial total membrane HER2"
  )

  covariateData <- list()

  population <- list(
    species = "in vitro (SKBR3 cell line)",
    n_subjects = NA_integer_,
    n_studies = 1,
    cell_line = "SKBR3 (human HER2-overexpressing breast carcinoma; ATCC)",
    n_cells_per_condition = "approximately 200 (range 66-302); Tables 1 and 4",
    drug_concentration = "10 ug/mL trastuzumab in serum-free DMEM",
    affibody_concentration = "0.2 uM anti-HER2 (ZHER2:477)2 Affibody with biotin-streptavidin Qdot-655 label",
    temperature = "37 deg C continuous; cells never cooled below room temperature",
    disease_state = "HER2-overexpressing breast cancer (cell-line model system)",
    regions = "in vitro",
    notes = "Two-population model partitions HER2 between flat membrane regions (NF) and ruffled regions (NR) across the SKBR3 population. The estimated fraction in flat regions (NF(0) = 0.62) exceeds the approximately 6 percent share of completely flat cells because bulk ruffled cells also contain flat membrane sub-areas. The receptor counts are normalized so that NR(0) + NF(0) = 1; absolute counts in the experiments are approximately 2e6 HER2 per cell (Methods)."
  )

  ini({
    # Mechanistic rate constants. All times are minutes. Values are Model B
    # final estimates from Table 3 (right column).
    kact_R_T0 <- 0.4;          label("Effective trastuzumab activation rate on ruffled membrane (kact_R * T0, 1/min)")  # Table 3 Model B
    kact_F_T0 <- fixed(0.41);  label("Effective trastuzumab activation rate on flat membrane (kact_F * T0, 1/min); paper holds equal to kact_R per footnote *") # Table 3 Model B footnote *
    kdiss     <- 0.05;         label("Trastuzumab dissociation rate (kdiss, 1/min)")              # Table 3 Model B (5.0e-2 /min)
    kon_A0    <- fixed(1.0);   label("Effective Affibody association rate (kon * A0, 1/min); non-identifiable in Model B, fixed at Model A value per footnote **") # Table 3 Model B footnote **
    koff      <- 6.8e-4;       label("Affibody dissociation rate (koff, 1/min)")                  # Table 3 Model B (6.8e-4 /min)
    kint_R_T  <- 49;           label("Internalization rate of trastuzumab-bound HER2 on ruffled membrane (kint_RT, 1/min); only lower bound identifiable") # Table 3 Model B
    kint_F_T  <- 0.0126;       label("Internalization rate of trastuzumab-bound HER2 on flat membrane (kint_FT, 1/min)") # Table 3 Model B
    nf0       <- 0.62;         label("Initial fraction of receptors in flat membrane regions (NF(0), unitless)") # Table 3 Model B
    sc        <- 535;          label("Observation scaling factor sc (arbitrary QD-fluorescence units per fractional receptor)") # Table 3 Model B

    # Experimental-design switches: when the drug and Affibody bath are
    # available. Encoded as fixed values that the user overrides per
    # simulation via rxSolve(params = ...). Default is the standard 60 min
    # trastuzumab incubation with a 10 min Affibody pre-incubation
    # (Methods: Trastuzumab incubation and labeling of HER2).
    drug_start_min  <- fixed(0);    label("Time trastuzumab bath turns on (min); user-overridable")
    drug_stop_min   <- fixed(60);   label("Time trastuzumab bath turns off (min); user-overridable")
    affib_start_min <- fixed(-10);  label("Time Affibody bath turns on (min); user-overridable. Negative value implements the 10 min pre-incubation described in Methods.")
    affib_stop_min  <- fixed(1e6);  label("Time Affibody bath turns off (min); user-overridable. Default value keeps Affibody present indefinitely (no washout).")
  })

  model({
    # Region populations sum to 1 by normalization (Results: "the freedom of
    # scale is used to fix NR(t=0) + NF(t=0) to one").
    nr0 <- 1 - nf0

    # Time-gated availability of trastuzumab and Affibody bath. The paper
    # implements these as boolean switches (Results, recycling model
    # paragraph: "T = NT * s_on with [switch definition]"); reproduced here
    # via the time variable. The product of two comparisons evaluates to
    # 1 when the bath is on and 0 otherwise.
    trast_on <- (time >= drug_start_min) * (time < drug_stop_min)
    affib_on <- (time >= affib_start_min) * (time < affib_stop_min)

    flux_actR <- kact_R_T0 * trast_on
    flux_actF <- kact_F_T0 * trast_on
    flux_kon  <- kon_A0    * affib_on

    # ------------------------------------------------------------------
    # Membrane states - ruffled regions (Table 2 + Model B reductions:
    # kprod = 0, kdeg = 0, krec = 0, kint(unbound) = 0).
    # ------------------------------------------------------------------
    d/dt(nm_r)   <- -flux_actR * nm_r   + kdiss * nmt_r  - flux_kon * nm_r   + koff * nma_r
    d/dt(nmt_r)  <-  flux_actR * nm_r   - kdiss * nmt_r  - kint_R_T * nmt_r  - flux_kon * nmt_r  + koff * nmta_r
    d/dt(nma_r)  <- -flux_actR * nma_r  + kdiss * nmta_r + flux_kon * nm_r   - koff * nma_r
    d/dt(nmta_r) <-  flux_actR * nma_r  - kdiss * nmta_r - kint_R_T * nmta_r + flux_kon * nmt_r  - koff * nmta_r

    # ------------------------------------------------------------------
    # Membrane states - flat regions
    # ------------------------------------------------------------------
    d/dt(nm_f)   <- -flux_actF * nm_f   + kdiss * nmt_f  - flux_kon * nm_f   + koff * nma_f
    d/dt(nmt_f)  <-  flux_actF * nm_f   - kdiss * nmt_f  - kint_F_T * nmt_f  - flux_kon * nmt_f  + koff * nmta_f
    d/dt(nma_f)  <- -flux_actF * nma_f  + kdiss * nmta_f + flux_kon * nm_f   - koff * nma_f
    d/dt(nmta_f) <-  flux_actF * nma_f  - kdiss * nmta_f - kint_F_T * nmta_f + flux_kon * nmt_f  - koff * nmta_f

    # ------------------------------------------------------------------
    # Internal states - ruffled regions. With Model B reductions, the
    # only source of internal trastuzumab-bound receptors is the
    # internalization flux kint_R_T * nmt_r (or nmta_r). Recycling /
    # degradation / unbound internalization are all zero.
    # ------------------------------------------------------------------
    d/dt(ni_r)   <-  kdiss * nit_r                       + koff * nia_r
    d/dt(nit_r)  <- -kdiss * nit_r   + kint_R_T * nmt_r  + koff * nita_r
    d/dt(nia_r)  <-  kdiss * nita_r                      - koff * nia_r
    d/dt(nita_r) <- -kdiss * nita_r  + kint_R_T * nmta_r - koff * nita_r

    # ------------------------------------------------------------------
    # Internal states - flat regions
    # ------------------------------------------------------------------
    d/dt(ni_f)   <-  kdiss * nit_f                       + koff * nia_f
    d/dt(nit_f)  <- -kdiss * nit_f   + kint_F_T * nmt_f  + koff * nita_f
    d/dt(nia_f)  <-  kdiss * nita_f                      - koff * nia_f
    d/dt(nita_f) <- -kdiss * nita_f  + kint_F_T * nmta_f - koff * nita_f

    # Initial conditions: at t = 0 the system is in equilibrium without
    # drug or Affibody (Methods: pre-experiment serum starvation), so all
    # receptors are unbound on the membrane.
    nm_r(0)   <- nr0
    nm_f(0)   <- nf0
    nmt_r(0)  <- 0
    nmt_f(0)  <- 0
    nma_r(0)  <- 0
    nma_f(0)  <- 0
    nmta_r(0) <- 0
    nmta_f(0) <- 0
    ni_r(0)   <- 0
    ni_f(0)   <- 0
    nit_r(0)  <- 0
    nit_f(0)  <- 0
    nia_r(0)  <- 0
    nia_f(0)  <- 0
    nita_r(0) <- 0
    nita_f(0) <- 0

    # Observation: Affibody-bound receptors on the membrane (only species
    # that contribute to QD fluorescence). Reproduces N_obs definition in
    # the Extended-model section ("the observed receptors on the membrane
    # thereby consist of contributions of both populations...").
    nobs <- sc * (nma_r + nmta_r + nma_f + nmta_f)

    # Diagnostic outputs (not used by the paper's fit but useful for
    # mass-balance and steady-state checks in the validation vignette).
    nm_total  <- nm_r + nmt_r + nma_r + nmta_r + nm_f + nmt_f + nma_f + nmta_f
    ni_total  <- ni_r + nit_r + nia_r + nita_r + ni_f + nit_f + nia_f + nita_f
    n_total   <- nm_total + ni_total
  })
}
