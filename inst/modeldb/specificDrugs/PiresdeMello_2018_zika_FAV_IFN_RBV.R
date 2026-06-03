PiresdeMello_2018_zika_FAV_IFN_RBV <- function() {
  description <- paste(
    "In vitro (Vero cells). Translational mechanism-based pharmacodynamic",
    "(MBM) model of Zika virus replication and inhibition by favipiravir",
    "(FAV), interferon alpha (IFN), and ribavirin (RBV) as monotherapy and",
    "in two-drug combinations. Eight-state model: uninfected (U) and",
    "infected (I) host cells, five intracellular virus transit compartments",
    "(vi1..vi5) capturing maturation delay, and extracellular virus (vextra)",
    "as the observation output (log10 PFU/mL). IFN inhibits cellular",
    "infection via a sigmoidal Hill function; FAV and RBV both inhibit the",
    "vi4 -> vi5 maturation transit; RBV additionally causes first-order",
    "cytotoxicity to both uninfected and infected host cells. FAV+RBV",
    "antagonism is encoded via a competitive-interaction factor PSI",
    "(= 1 monotherapy, = 1.37 combination). Drug concentrations are static",
    "covariates -- the in vitro experiment fixes nominal concentrations",
    "for the 4-day window. All parameters fixed at the Table 1 point",
    "estimates; the between-curve CVs reported in Table 1 are not encoded",
    "as etas (typical-value mechanism)."
  )
  reference <- paste(
    "Pires de Mello CP, Tao X, Kim TH, Bulitta JB, Rodriquez JL, Pomeroy JJ,",
    "Brown AN. (2018). Zika virus replication is substantially inhibited by",
    "novel favipiravir and interferon alpha combination regimens.",
    "Antimicrob Agents Chemother 62:e01983-17.",
    "doi:10.1128/AAC.01983-17."
  )
  vignette <- "PiresdeMello_2018_zika_FAV_IFN_RBV"

  paper_specific_compartments <- c(
    "uninfected", "infected",
    "vi1", "vi2", "vi3", "vi4", "vi5",
    "vextra"
  )

  units <- list(
    time          = "hour",
    dosing        = "static covariates (uM FAV, ug/mL RBV, IU/mL IFN) -- not administered events",
    concentration = "log10(PFU/mL) for the model observation Cc"
  )

  covariateData <- list(
    CONC_FAV_UM = list(
      description        = "Static favipiravir extracellular concentration in the in vitro time-kill assay (uM)",
      units              = "uM",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-invariant per Pires de Mello 2018 Materials and methods",
        "(Antiviral evaluations): combination assays used a 6-by-6 checkerboard",
        "format. FAV tested at 0, 31.25, 62.5, 125, 250, and 500 uM (Figure 1A).",
        "Set to 0 in regimens without FAV. Drives the Hill-Imax inhibition of",
        "the intracellular virus vi4 -> vi5 transit (Eq 10) and the combination",
        "competitive-interaction term with RBV (Eq 13). Paper-specific covariate",
        "not in inst/references/covariate-columns.md because the canonical",
        "concentration concept in nlmixr2lib is a state-derived plasma",
        "concentration (Cc), not a static exogenous-drug-concentration",
        "covariate used to drive an in vitro PD model."
      ),
      source_name        = "C_FAV (Pires de Mello 2018 Eq 10)"
    ),
    CONC_RBV_UGML = list(
      description        = "Static ribavirin extracellular concentration in the in vitro time-kill assay (ug/mL)",
      units              = "ug/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-invariant. RBV tested at 0, 1, 10, 100, and 1,000 ug/mL",
        "(Figure 1B). Set to 0 in regimens without RBV. Drives three effects:",
        "the Hill-Imax inhibition of vi4 -> vi5 transit (Eq 11), the FAV+RBV",
        "antagonism factor (Eq 13), and the first-order cytotoxicity rate",
        "k_cytotox on uninfected and infected host cells and on the",
        "intracellular virus chain (Eq 12). Paper-specific covariate."
      ),
      source_name        = "C_RBV (Pires de Mello 2018 Eq 11)"
    ),
    CONC_IFN_IUML = list(
      description        = "Static interferon alpha extracellular concentration in the in vitro time-kill assay (IU/mL)",
      units              = "IU/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-invariant. IFN tested at 0, 10, 100, 1,000, and 10,000 IU/mL",
        "(Figure 1C). Set to 0 in regimens without IFN. Drives the Hill-Imax",
        "inhibition of cellular infection (Eq 9), entering as a multiplicative",
        "factor on the second-order infection rate kinfect * INH_IFN * vextra",
        "* U. Paper-specific covariate."
      ),
      source_name        = "C_IFN (Pires de Mello 2018 Eq 9)"
    )
  )

  population <- list(
    species        = "in vitro (Vero cells)",
    n_subjects     = NA_integer_,
    n_studies      = 1L,
    cell_line      = "Vero (ATCC CCL-81, American Type Culture Collection)",
    virus_strain   = "Zika virus, 2015 human Puerto Rican strain PRVABC59 (BEI Resources)",
    inoculum       = paste(
      "Multiplicity of infection 0.01 PFU/cell on confluent Vero cell",
      "monolayers. Initial total uninfected cells = 10^6.30 ~ 1.995e6",
      "cells/mL; initial infected cells = 10^3.38 ~ 2399 cells/mL",
      "(Table 1 Log_U fixed; Log_I estimated)."
    ),
    disease_state  = paste(
      "Zika virus in vitro time-course infection experiments in Vero CCL-81",
      "cells cultured in Eagle's MEM with 5% fetal bovine serum and 1%",
      "penicillin-streptomycin at 37 C in 5% CO2. Plaque assay limit of",
      "detection 100 PFU/mL (log10 = 2). Single-agent assays sampled daily",
      "for 4 days; combination assays sampled on day 3 posttreatment",
      "(peak viral burden). Three independent samples per regimen.",
      "Beal M3 method used in S-ADAPT to handle the BLQ samples at time zero."
    ),
    dose_range     = paste(
      "Static drug concentrations (no PK dosing) as monotherapy and in all",
      "two-drug combinations of FAV, RBV, and IFN. FAV 0-500 uM; RBV",
      "0-1,000 ug/mL; IFN 0-10,000 IU/mL. Combination experiments used a",
      "6-by-6 checkerboard."
    ),
    notes          = paste(
      "Single experiment per drug regimen with three independent replicate",
      "samples. Between-curve CVs reported in Table 1 (e.g. 0.0841 on",
      "log10 k_infect, 0.365 on Log_I, 0.491 on IC50_RBV) reflect day-to-day",
      "biological variability of separate in vitro curves; they are NOT",
      "encoded as etas in this model file -- the registry uses the typical-",
      "value mechanism. See the validation vignette Errata for the rationale.",
      "Model fit and parameter estimation performed in S-ADAPT (version 1.57)",
      "via Monte Carlo parametric expectation maximization; simulations in",
      "the original paper used Berkeley Madonna (version 8.23.3.0)."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Host-cell dynamics (Pires de Mello 2018 Table 1, Eqs 1-2).
    # log10kinfect is kept on the log10 scale (matching the paper's reported
    # form) and back-transformed in model() via kinfect <- 10^log10kinfect.
    # ---------------------------------------------------------------------
    log10kinfect  <- fixed(-4.10)
    label("log10 of 2nd-order virus-host infection rate constant (mL/PFU/h)")
    # Table 1 row 1: log10(k_infect) = -4.10 (RSE 2.39 percent).

    lksyn         <- fixed(log(9.35))
    label("ksyn -- virus synthesis rate constant in infected cells (1/h)")
    # Table 1 row 2: k_syn = 9.35 1/h (RSE 7.09 percent).

    lktr          <- fixed(log(0.125))
    label("ktr -- intracellular-virus transit rate constant (1/h)")
    # Table 1 row 3: T_Delay = 5/ktr = 40.0 h (RSE 2.56 percent) -> ktr = 5/40 = 0.125 1/h.

    lkdeath       <- fixed(log(1 / 70.5))
    label("kdeath -- first-order death rate constant of infected host cells (1/h)")
    # Table 1 row 4: MST_Infected = 1/k_death = 70.5 h (RSE 8.92 percent) -> k_death = 1/70.5 1/h.

    lklossvirus   <- fixed(log(1 / 14.3))
    label("klossvirus -- first-order loss rate constant for extracellular virus (1/h)")
    # Table 1 row 5: MST_Virus = 1/k_loss,virus = 14.3 h (RSE 10.4 percent) -> k_loss = 1/14.3 1/h.

    log10U0       <- fixed(6.30)
    label("log10 initial number of uninfected host cells (cells/mL)")
    # Table 1 row 6: Log_U = 6.30 (fixed); U(0) = 10^6.30 ~ 1.995e6 cells/mL.

    log10I0       <- fixed(3.38)
    label("log10 initial number of infected host cells (cells/mL)")
    # Table 1 row 7: Log_I = 3.38 (RSE 2.66 percent); I(0) = 10^3.38 ~ 2399 cells/mL.

    # ---------------------------------------------------------------------
    # Drug effects: FAV (Pires de Mello 2018 Table 1, Eq 10).
    # ---------------------------------------------------------------------
    imax_fav      <- fixed(0.9999)
    label("Imax_FAV -- maximum extent of inhibition by FAV on vi4->vi5 transit (fraction)")
    # Table 1: Imax_FAV = 0.9999 (95 percent CI ~0.9992-1.00) on normal scale.

    lic50_fav     <- fixed(log(41.7))
    label("IC50_FAV -- FAV concentration causing 50 percent of Imax_FAV (uM)")
    # Table 1: IC50_FAV = 41.7 uM (RSE 2.55 percent).

    lhill_fav     <- fixed(log(2.79))
    label("Hill_FAV -- Hill coefficient for FAV inhibition (unitless)")
    # Table 1: Hill_FAV = 2.79 (RSE 4.53 percent).

    # ---------------------------------------------------------------------
    # Drug effects: RBV (Pires de Mello 2018 Table 1, Eqs 11-12).
    # ---------------------------------------------------------------------
    imax_rbv      <- fixed(0.954)
    label("Imax_RBV -- maximum extent of inhibition by RBV on vi4->vi5 transit (fraction)")
    # Table 1: Imax_RBV = 0.954 (95 percent CI ~0.924-0.973).

    lic50_rbv     <- fixed(log(7.86))
    label("IC50_RBV -- RBV concentration causing 50 percent of Imax_RBV (ug/mL)")
    # Table 1: IC50_RBV = 7.86 ug/mL (RSE 9.99 percent).

    lhill_rbv     <- fixed(log(2.90))
    label("Hill_RBV -- Hill coefficient for RBV inhibition (unitless)")
    # Table 1: Hill_RBV = 2.90 (RSE 16.5 percent).

    lsmax_rbv     <- fixed(log(1 / 11.9))
    label("Smax_RBV -- maximum RBV cytotoxicity rate constant (1/h)")
    # Table 1: MST_TOX = 1/Smax_RBV = 11.9 h (RSE 7.21 percent) -> Smax_RBV = 1/11.9 1/h.

    lsc50_rbv     <- fixed(log(150))
    label("SC50_RBV -- RBV concentration causing 50 percent of Smax_RBV (ug/mL)")
    # Table 1: SC50_RBV = 150 ug/mL (RSE 11.3 percent).

    lhill_rbvtox  <- fixed(log(4.16))
    label("Hill_RBVTOX -- Hill coefficient for RBV cytotoxicity (unitless)")
    # Table 1: Hill_RBVTOX = 4.16 (RSE 12.5 percent).

    # ---------------------------------------------------------------------
    # Drug effects: IFN (Pires de Mello 2018 Table 1, Eq 9).
    # ---------------------------------------------------------------------
    imax_ifn      <- fixed(0.99997)
    label("Imax_IFN -- maximum extent of inhibition by IFN on cellular infection (fraction)")
    # Table 1: Imax_IFN = 0.99997 (95 percent CI ~0.99990-1.00).

    lic50_ifn     <- fixed(log(4.12))
    label("IC50_IFN -- IFN concentration causing 50 percent of Imax_IFN (IU/mL)")
    # Table 1: IC50_IFN = 4.12 IU/mL (RSE 15.9 percent).

    lhill_ifn     <- fixed(log(2.00))
    label("Hill_IFN -- Hill coefficient for IFN inhibition (unitless, fixed)")
    # Table 1: Hill_IFN = 2.00 (fixed; informed by IFN monotherapy modelling).

    # ---------------------------------------------------------------------
    # FAV+RBV combination antagonism factor (Pires de Mello 2018 Eq 13).
    # PSI = 1 in monotherapy; PSI = 1.37 in any FAV+RBV combination.
    # Eq 13 already reduces to the mono Hill form (Eq 10 / Eq 11) when the
    # other drug concentration is 0, so the same single expression handles
    # mono and combination input without a branch.
    # ---------------------------------------------------------------------
    psi           <- fixed(1.37)
    label("PSI -- FAV+RBV competitive-interaction antagonism factor (unitless)")
    # Table 1: PSI = 1 (monotherapy) or 1.37 (combination), RSE 8.21 percent.

    # ---------------------------------------------------------------------
    # Residual error (Pires de Mello 2018 Table 1 footer; System outputs
    # and residual-error model section: additive error on log10 PFU/mL).
    # ---------------------------------------------------------------------
    addSd         <- fixed(0.333)
    label("addSd -- additive residual SD on log10 PFU/mL scale (Beal M3 method for BLQ)")
    # Table 1 residual-error row: SDin = 0.333 (RSE 4.67 percent).
  })

  model({
    # ================================================================
    # 1. Back-transform structural parameters from log / log10 scale.
    # ================================================================
    kinfect     <- 10^log10kinfect
    ksyn        <- exp(lksyn)
    ktr         <- exp(lktr)
    kdeath      <- exp(lkdeath)
    klossvirus  <- exp(lklossvirus)
    ic50_fav    <- exp(lic50_fav)
    ic50_rbv    <- exp(lic50_rbv)
    ic50_ifn    <- exp(lic50_ifn)
    hill_fav    <- exp(lhill_fav)
    hill_rbv    <- exp(lhill_rbv)
    hill_ifn    <- exp(lhill_ifn)
    smax_rbv    <- exp(lsmax_rbv)
    sc50_rbv    <- exp(lsc50_rbv)
    hill_rbvtox <- exp(lhill_rbvtox)

    # ================================================================
    # 2. Initial conditions (Pires de Mello 2018 Table 1).
    #    U(0) = 10^Log_U = 10^6.30; I(0) = 10^Log_I = 10^3.38.
    #    Intracellular virus chain (vi1..vi5) and extracellular virus
    #    (vextra) default to 0 (Materials and methods, Viral replication).
    # ================================================================
    uninfected(0) <- 10^log10U0
    infected(0)   <- 10^log10I0

    # ================================================================
    # 3. Drug-effect terms.
    #    INH_IFN (Eq 9), INH for FAV+RBV (Eq 13), k_cytotox (Eq 12).
    #    The +1e-30 floors guard CONC^Hill / 0^0 limits when a drug
    #    concentration is exactly 0; Hill > 0 always so the floor only
    #    affects the double-zero limit.
    # ================================================================
    cfav <- CONC_FAV_UM   + 1e-30
    crbv <- CONC_RBV_UGML + 1e-30
    cifn <- CONC_IFN_IUML + 1e-30

    inh_ifn <- 1 - imax_ifn * cifn^hill_ifn /
      (cifn^hill_ifn + ic50_ifn^hill_ifn)

    # Combination antagonism (Eq 13). Per Table 1, PSI = 1 in monotherapy
    # (either drug absent) and PSI = 1.37 only when both FAV and RBV are
    # present. The indicator form below evaluates psi_eff = psi when both
    # CONC_FAV_UM and CONC_RBV_UGML are strictly positive, and psi_eff = 1
    # otherwise -- so a FAV-only or RBV-only regimen recovers the mono
    # Hill IC50 without antagonism scaling.
    psi_eff <- 1 + (psi - 1) * (CONC_FAV_UM > 0) * (CONC_RBV_UGML > 0)
    xfav <- (cfav / (psi_eff * ic50_fav))^hill_fav
    xrbv <- (crbv / (psi_eff * ic50_rbv))^hill_rbv
    inh  <- 1 - (imax_fav * xfav + imax_rbv * xrbv) /
      (xfav + xrbv + 1)

    k_cytotox <- smax_rbv * crbv^hill_rbvtox /
      (crbv^hill_rbvtox + sc50_rbv^hill_rbvtox)

    # ================================================================
    # 4. ODE system (Pires de Mello 2018 Eqs 1-8).
    #    Eq 1: dU/dt    = -kinfect * INH_IFN * vextra * U - k_cytotox * U
    #    Eq 2: dI/dt    =  kinfect * INH_IFN * vextra * U - (kdeath + k_cytotox) * I
    #    Eq 3: dvi1/dt  =  ksyn * I - ktr * vi1   - (kdeath + k_cytotox) * vi1
    #    Eq 4: dvi2/dt  =  ktr * (vi1 - vi2)      - (kdeath + k_cytotox) * vi2
    #    Eq 5: dvi3/dt  =  ktr * (vi2 - vi3)      - (kdeath + k_cytotox) * vi3
    #    Eq 6: dvi4/dt  =  ktr * (vi3 - vi4*INH)  - (kdeath + k_cytotox) * vi4
    #    Eq 7: dvi5/dt  =  ktr * (vi4*INH - vi5)  - (kdeath + k_cytotox) * vi5
    #    Eq 8: dvextra/dt = ktr * vi5 - klossvirus * vextra
    # ================================================================
    d/dt(uninfected) <- -kinfect * inh_ifn * vextra * uninfected -
      k_cytotox * uninfected
    d/dt(infected)   <-  kinfect * inh_ifn * vextra * uninfected -
      (kdeath + k_cytotox) * infected
    d/dt(vi1)        <-  ksyn * infected - ktr * vi1 -
      (kdeath + k_cytotox) * vi1
    d/dt(vi2)        <-  ktr * (vi1 - vi2) -
      (kdeath + k_cytotox) * vi2
    d/dt(vi3)        <-  ktr * (vi2 - vi3) -
      (kdeath + k_cytotox) * vi3
    d/dt(vi4)        <-  ktr * (vi3 - vi4 * inh) -
      (kdeath + k_cytotox) * vi4
    d/dt(vi5)        <-  ktr * (vi4 * inh - vi5) -
      (kdeath + k_cytotox) * vi5
    d/dt(vextra)     <-  ktr * vi5 - klossvirus * vextra

    # ================================================================
    # 5. Observation: log10 PFU/mL of extracellular virus.
    #    Plaque assay LoD is 100 PFU/mL (log10 = 2); the Beal M3 method
    #    is used in NONMEM/S-ADAPT to handle the BLQ samples at time zero.
    #    A small floor (1e-30) avoids -Inf when vextra is exactly 0 at
    #    t = 0; observed data are sampled from t = 24 h onward, by which
    #    point vextra is well above the floor.
    # ================================================================
    Cc <- log10(vextra + 1e-30)
    Cc ~ add(addSd)
  })
}
