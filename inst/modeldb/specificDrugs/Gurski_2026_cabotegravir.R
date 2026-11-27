Gurski_2026_cabotegravir <- function() {
  description <- paste(
    "QSP (within-host HIV-1 viral dynamics with Darwinian trait evolution).",
    "Long-acting cabotegravir (CAB-LA) pre-exposure prophylaxis during early",
    "HIV-1 infection, from Gurski et al. 2026. Plasma cabotegravir declines",
    "mono-exponentially after each injection (superposition across",
    "injections) and blocks the integration step (early-to-productive",
    "T-cell transition) through an Emax model. Healthy CD4+ T cells are",
    "infected either by free virions or by cell-to-cell transfer, giving two",
    "early-infected classes that progress to a common productively infected",
    "class releasing free virus. Two continuous mean resistance traits (free",
    "virion and cell-to-cell) lower infectivity as exp(-u^2/2) and escape the",
    "drug as exp(-u^2/(2 sigma^2)), and evolve along the gradient of the",
    "Malthusian fitness r (Lande's equation), where r is the dominant",
    "eigenvalue of the Jacobian at the disease-free equilibrium, computed in",
    "the model. Deterministic typical-value model with no IIV and no",
    "residual error; the evolutionary speeds were not reported and are",
    "back-solved defaults."
  )
  reference <- paste(
    "Gurski K, Kang Y, Ma Y. A Darwinian model for the evolution of drug",
    "resistance to long-acting PrEP during an early HIV infection. PLoS One.",
    "2026;21(8):e0355833. doi:10.1371/journal.pone.0355833.",
    "PMID 42623409, PMCID PMC13492818. PK/PD from Eq 3, structure from",
    "Eq 5 and Eq 31, parameter values from Table 2, initial conditions from",
    "Table 1, trait functions from Table 4, fitness definition from",
    "Section 8 (Eq 19, Eq 32-34)."
  )
  vignette <- "Gurski_2026_cabotegravir"

  # Time in days. One dose record of amt = 1 is one CAB-LA injection;
  # central holds the remaining injection-equivalents and Cc converts them
  # to the plasma concentration of Eq 3 (5.04 ug/mL per injection at the
  # time of the injection). Cell and virion states are concentrations per
  # microlitre of blood, as in Table 1.
  units <- list(time = "day", dosing = "injection", concentration = "ug/mL")

  dosing <- "central"

  paper_specific_compartments <- c(
    "tcell_uninfected",
    "tcell_early_free",
    "tcell_early_c2c",
    "tcell_infected",
    "trait_free",
    "trait_c2c"
  )

  covariateData <- list()

  compartmentData <- list(
    central = list(
      analyte = "cabotegravir",
      units = "injection",
      specimen = "plasma",
      verified = TRUE
    ),
    tcell_uninfected = list(
      analyte = "healthy CD4+ T cells",
      units = "cells/uL",
      specimen = "blood cell",
      verified = TRUE
    ),
    tcell_early_free = list(
      analyte = "early infected CD4+ T cells, infected by free virions (pre-integration)",
      units = "cells/uL",
      specimen = "blood cell",
      verified = TRUE
    ),
    tcell_early_c2c = list(
      analyte = "early infected CD4+ T cells, infected by cell-to-cell transfer (pre-integration)",
      units = "cells/uL",
      specimen = "blood cell",
      verified = TRUE
    ),
    tcell_infected = list(
      analyte = "productively infected CD4+ T cells",
      units = "cells/uL",
      specimen = "blood cell",
      verified = TRUE
    ),
    virus = list(
      analyte = "free HIV-1 virions",
      units = "virions/uL",
      specimen = "plasma",
      verified = TRUE
    ),
    trait_free = list(
      analyte = "mean free-virion resistance trait u_k",
      units = "(unitless)",
      specimen = "not applicable",
      verified = TRUE
    ),
    trait_c2c = list(
      analyte = "mean cell-to-cell resistance trait u_w",
      units = "(unitless)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = NA_integer_,
    disease_state = paste(
      "Theoretical HIV-negative adult on bimonthly long-acting cabotegravir",
      "PrEP who is exposed to HIV-1 (acute infection phase only, no latent",
      "reservoir)"
    ),
    dose_range = paste(
      "CAB-LA injections at day 0 and day 14, then every 56 days (Fig 3,",
      "WHO 2022 guideline as described by the paper); exposure assumed at",
      "the third injection"
    ),
    regions = "Not applicable (theoretical model)",
    notes = paste(
      "No cohort was fitted. T-cell and viral parameters are literature",
      "values (Table 2) except the cell-to-cell infection rate omega_m, fitted",
      "by nudging and Nelder-Mead to the smoothed mean viral load of 46",
      "acute HIV-1 infections from the RV217 cohort. The PK/PD curve (Eq 3)",
      "was fitted in the authors' earlier work to the single-injection human",
      "plasma data of Shaik et al. The resistance-breadth parameters",
      "sigma_k and sigma_w are scenario values (0.5, 1, 2 in the paper)."
    )
  )

  ini({
    # ---- CAB-LA PK/PD (Eq 3; fitted in the authors' earlier work) ----------
    lcinj <- fixed(log(5.04)); label("Plasma cabotegravir concentration immediately after one injection (ug/mL)") # Eq 3, C(t) = 5.04 exp(-0.0363 t)
    lkel <- fixed(log(0.0363)); label("Cabotegravir elimination rate constant (1/day)") # Eq 3, decay constant 0.0363 day^-1
    lic50 <- fixed(log(0.018)); label("Cabotegravir concentration for 50% inhibition of integration IC50 (ug/mL)") # Eq 3, IC50 = 0.018 (Hill coefficient n = 1, Eq 1-2)

    # ---- Within-host HIV-1 dynamics (Table 2) ------------------------------
    lambda_t <- fixed(60); label("Recruitment rate of healthy T cells beta (cells/uL/day)") # Table 2, beta = 60 cells per day per uL
    delta_t <- fixed(0.05); label("Natural death rate of T cells mu (1/day)") # Table 2, mu = 0.05 day^-1
    delta_c2c <- fixed(0.25); label("Death rate of early cells infected by cell-to-cell transfer mu_c (1/day)") # Table 2, mu_c = 0.25 day^-1
    delta_burst <- fixed(1); label("Bursting (virus-induced death) rate of productively infected T cells delta (1/day)") # Table 2, delta = 1 day^-1
    n_burst <- fixed(100); label("Virions released per burst of an infected cell N (virions/cell)") # Table 2, N = 100 per infected cell
    cl_virus <- fixed(25); label("Free virion clearance rate c (1/day)") # Table 2, c = 25 day^-1
    k_t <- fixed(1.6); label("Progression rate from early to productive infection phi (1/day)") # Table 2, phi = 24/15 day^-1 (calculated)
    beta_t <- fixed(2.5e-4); label("Maximum free-virion infection rate k_m (uL/virion/day)") # Table 2, k_m = 2.5e-4 per virions per uL
    omega_c2c <- 2.93e-3; label("Maximum cell-to-cell infection rate omega_m (uL/cell/day)") # Table 2, omega_m = 2.93e-3 per infected cell per uL (fitted)

    # ---- Darwinian trait dynamics (Section 8, Tables 3-4, Eq 19/31) --------
    sigma_free <- 2; label("Breadth of CAB-LA effectiveness against the free-virion trait sigma_k (unitless)") # Section 8, scenario values 0.5, 1, 2 (Figs 12-18)
    sigma_c2c <- 2; label("Breadth of CAB-LA effectiveness against the cell-to-cell trait sigma_w (unitless)") # Section 8, scenario values 0.5, 1, 2 (Figs 10-18)
    s_free <- 0.018; label("Speed of evolution of the free-virion trait s_k (1/day per unit fitness gradient)") # NOT REPORTED; back-solved by the maintainers from Fig 16
    s_c2c <- 0.024; label("Speed of evolution of the cell-to-cell trait s_w (1/day per unit fitness gradient)") # NOT REPORTED; back-solved by the maintainers from Fig 16
    trait_free0 <- 0.1; label("Mean free-virion trait u_k at HIV exposure (unitless)") # Figs 15-16, starting traits uk = uw = 0.1
    trait_c2c0 <- 0.1; label("Mean cell-to-cell trait u_w at HIV exposure (unitless)") # Figs 15-16, starting traits uk = uw = 0.1

    # ---- Initial conditions at HIV exposure (Table 1) ----------------------
    tcell0 <- fixed(1000); label("Healthy T cells at exposure H(0) (cells/uL)") # Table 1, H = 10^3 cells/uL
    tcell_infected0 <- fixed(1e-3); label("Productively infected T cells at exposure I(0) (cells/uL)") # Table 1, I = 10^-3 cells/uL
    virus0 <- fixed(100); label("Free virions at exposure V(0) (virions/uL)") # Table 1, V = 100 virions/uL
  })

  model({
    # ---- 1. PK/PD (Eq 1, 3, 4) -----------------------------------------------
    cinj <- exp(lcinj)
    kel <- exp(lkel)
    ic50 <- exp(lic50)
    d/dt(central) <- -kel * central
    Cc <- cinj * central
    # Eq 1 with n = 1: Feff = 1 - IC50 / (IC50 + C)
    feff <- Cc / (ic50 + Cc)

    # ---- 2. Trait-dependent rates (Table 4) ---------------------------------
    uk <- trait_free
    uw <- trait_c2c
    kinf <- beta_t * exp(-uk^2 / 2)
    winf <- omega_c2c * exp(-uw^2 / 2)
    esc_k <- exp(-uk^2 / (2 * sigma_free^2))
    esc_w <- exp(-uw^2 / (2 * sigma_c2c^2))
    d_free <- k_t * (1 - feff * esc_k)
    d_c2c <- k_t * (1 - feff * esc_w)

    # ---- 3. Malthusian fitness r (Section 8) ---------------------------------
    # r is the dominant (Perron) eigenvalue of the Jacobian of
    # (Ef, Ec, I, V) at the disease-free state H0 = beta / mu. Its
    # characteristic polynomial (equivalent to Eq 32) is
    #   P(l) = (l + dd)(l + aa)(l + ee)(l + bb) - q1 (l + bb)
    #          - q2 (l + aa)(l + ee)
    # The Jacobian is Metzler, so every root has real part <= r and P is
    # increasing and convex on (r, Inf); Newton's method started from an
    # upper bound therefore decreases monotonically to r. The start is the
    # Collatz-Wielandt bound max_i (J x)_i / x_i with x = (1, 1, 1, N delta / c).
    h0 <- lambda_t / delta_t
    aa <- d_free + delta_t
    bb <- d_c2c + delta_c2c
    dd <- delta_burst + delta_t
    ee <- cl_virus
    gg <- h0 * n_burst * delta_burst
    q1 <- d_free * kinf * gg
    q2 <- d_c2c * winf * h0
    lam <- max(max(kinf * gg / ee - aa, winf * h0 - bb), max(d_free + d_c2c - dd, 0)) + 0.01
    # 20 unrolled Newton steps: over 0.1 <= sigma <= 5, 0 <= u <= 6 and
    # 0 <= Feff <= 1 the error against eigen() is < 1e-11 after 15 steps.
    pl <- (lam + dd) * (lam + aa) * (lam + ee) * (lam + bb) - q1 * (lam + bb) - q2 * (lam + aa) * (lam + ee)
    dpl <- (lam + aa) * (lam + ee) * (lam + bb) + (lam + dd) * (lam + ee) * (lam + bb) + (lam + dd) * (lam + aa) * (lam + bb) + (lam + dd) * (lam + aa) * (lam + ee) - q1 - q2 * (2 * lam + aa + ee)
    lam <- lam - pl / dpl
    pl <- (lam + dd) * (lam + aa) * (lam + ee) * (lam + bb) - q1 * (lam + bb) - q2 * (lam + aa) * (lam + ee)
    dpl <- (lam + aa) * (lam + ee) * (lam + bb) + (lam + dd) * (lam + ee) * (lam + bb) + (lam + dd) * (lam + aa) * (lam + bb) + (lam + dd) * (lam + aa) * (lam + ee) - q1 - q2 * (2 * lam + aa + ee)
    lam <- lam - pl / dpl
    pl <- (lam + dd) * (lam + aa) * (lam + ee) * (lam + bb) - q1 * (lam + bb) - q2 * (lam + aa) * (lam + ee)
    dpl <- (lam + aa) * (lam + ee) * (lam + bb) + (lam + dd) * (lam + ee) * (lam + bb) + (lam + dd) * (lam + aa) * (lam + bb) + (lam + dd) * (lam + aa) * (lam + ee) - q1 - q2 * (2 * lam + aa + ee)
    lam <- lam - pl / dpl
    pl <- (lam + dd) * (lam + aa) * (lam + ee) * (lam + bb) - q1 * (lam + bb) - q2 * (lam + aa) * (lam + ee)
    dpl <- (lam + aa) * (lam + ee) * (lam + bb) + (lam + dd) * (lam + ee) * (lam + bb) + (lam + dd) * (lam + aa) * (lam + bb) + (lam + dd) * (lam + aa) * (lam + ee) - q1 - q2 * (2 * lam + aa + ee)
    lam <- lam - pl / dpl
    pl <- (lam + dd) * (lam + aa) * (lam + ee) * (lam + bb) - q1 * (lam + bb) - q2 * (lam + aa) * (lam + ee)
    dpl <- (lam + aa) * (lam + ee) * (lam + bb) + (lam + dd) * (lam + ee) * (lam + bb) + (lam + dd) * (lam + aa) * (lam + bb) + (lam + dd) * (lam + aa) * (lam + ee) - q1 - q2 * (2 * lam + aa + ee)
    lam <- lam - pl / dpl
    pl <- (lam + dd) * (lam + aa) * (lam + ee) * (lam + bb) - q1 * (lam + bb) - q2 * (lam + aa) * (lam + ee)
    dpl <- (lam + aa) * (lam + ee) * (lam + bb) + (lam + dd) * (lam + ee) * (lam + bb) + (lam + dd) * (lam + aa) * (lam + bb) + (lam + dd) * (lam + aa) * (lam + ee) - q1 - q2 * (2 * lam + aa + ee)
    lam <- lam - pl / dpl
    pl <- (lam + dd) * (lam + aa) * (lam + ee) * (lam + bb) - q1 * (lam + bb) - q2 * (lam + aa) * (lam + ee)
    dpl <- (lam + aa) * (lam + ee) * (lam + bb) + (lam + dd) * (lam + ee) * (lam + bb) + (lam + dd) * (lam + aa) * (lam + bb) + (lam + dd) * (lam + aa) * (lam + ee) - q1 - q2 * (2 * lam + aa + ee)
    lam <- lam - pl / dpl
    pl <- (lam + dd) * (lam + aa) * (lam + ee) * (lam + bb) - q1 * (lam + bb) - q2 * (lam + aa) * (lam + ee)
    dpl <- (lam + aa) * (lam + ee) * (lam + bb) + (lam + dd) * (lam + ee) * (lam + bb) + (lam + dd) * (lam + aa) * (lam + bb) + (lam + dd) * (lam + aa) * (lam + ee) - q1 - q2 * (2 * lam + aa + ee)
    lam <- lam - pl / dpl
    pl <- (lam + dd) * (lam + aa) * (lam + ee) * (lam + bb) - q1 * (lam + bb) - q2 * (lam + aa) * (lam + ee)
    dpl <- (lam + aa) * (lam + ee) * (lam + bb) + (lam + dd) * (lam + ee) * (lam + bb) + (lam + dd) * (lam + aa) * (lam + bb) + (lam + dd) * (lam + aa) * (lam + ee) - q1 - q2 * (2 * lam + aa + ee)
    lam <- lam - pl / dpl
    pl <- (lam + dd) * (lam + aa) * (lam + ee) * (lam + bb) - q1 * (lam + bb) - q2 * (lam + aa) * (lam + ee)
    dpl <- (lam + aa) * (lam + ee) * (lam + bb) + (lam + dd) * (lam + ee) * (lam + bb) + (lam + dd) * (lam + aa) * (lam + bb) + (lam + dd) * (lam + aa) * (lam + ee) - q1 - q2 * (2 * lam + aa + ee)
    lam <- lam - pl / dpl
    pl <- (lam + dd) * (lam + aa) * (lam + ee) * (lam + bb) - q1 * (lam + bb) - q2 * (lam + aa) * (lam + ee)
    dpl <- (lam + aa) * (lam + ee) * (lam + bb) + (lam + dd) * (lam + ee) * (lam + bb) + (lam + dd) * (lam + aa) * (lam + bb) + (lam + dd) * (lam + aa) * (lam + ee) - q1 - q2 * (2 * lam + aa + ee)
    lam <- lam - pl / dpl
    pl <- (lam + dd) * (lam + aa) * (lam + ee) * (lam + bb) - q1 * (lam + bb) - q2 * (lam + aa) * (lam + ee)
    dpl <- (lam + aa) * (lam + ee) * (lam + bb) + (lam + dd) * (lam + ee) * (lam + bb) + (lam + dd) * (lam + aa) * (lam + bb) + (lam + dd) * (lam + aa) * (lam + ee) - q1 - q2 * (2 * lam + aa + ee)
    lam <- lam - pl / dpl
    pl <- (lam + dd) * (lam + aa) * (lam + ee) * (lam + bb) - q1 * (lam + bb) - q2 * (lam + aa) * (lam + ee)
    dpl <- (lam + aa) * (lam + ee) * (lam + bb) + (lam + dd) * (lam + ee) * (lam + bb) + (lam + dd) * (lam + aa) * (lam + bb) + (lam + dd) * (lam + aa) * (lam + ee) - q1 - q2 * (2 * lam + aa + ee)
    lam <- lam - pl / dpl
    pl <- (lam + dd) * (lam + aa) * (lam + ee) * (lam + bb) - q1 * (lam + bb) - q2 * (lam + aa) * (lam + ee)
    dpl <- (lam + aa) * (lam + ee) * (lam + bb) + (lam + dd) * (lam + ee) * (lam + bb) + (lam + dd) * (lam + aa) * (lam + bb) + (lam + dd) * (lam + aa) * (lam + ee) - q1 - q2 * (2 * lam + aa + ee)
    lam <- lam - pl / dpl
    pl <- (lam + dd) * (lam + aa) * (lam + ee) * (lam + bb) - q1 * (lam + bb) - q2 * (lam + aa) * (lam + ee)
    dpl <- (lam + aa) * (lam + ee) * (lam + bb) + (lam + dd) * (lam + ee) * (lam + bb) + (lam + dd) * (lam + aa) * (lam + bb) + (lam + dd) * (lam + aa) * (lam + ee) - q1 - q2 * (2 * lam + aa + ee)
    lam <- lam - pl / dpl
    pl <- (lam + dd) * (lam + aa) * (lam + ee) * (lam + bb) - q1 * (lam + bb) - q2 * (lam + aa) * (lam + ee)
    dpl <- (lam + aa) * (lam + ee) * (lam + bb) + (lam + dd) * (lam + ee) * (lam + bb) + (lam + dd) * (lam + aa) * (lam + bb) + (lam + dd) * (lam + aa) * (lam + ee) - q1 - q2 * (2 * lam + aa + ee)
    lam <- lam - pl / dpl
    pl <- (lam + dd) * (lam + aa) * (lam + ee) * (lam + bb) - q1 * (lam + bb) - q2 * (lam + aa) * (lam + ee)
    dpl <- (lam + aa) * (lam + ee) * (lam + bb) + (lam + dd) * (lam + ee) * (lam + bb) + (lam + dd) * (lam + aa) * (lam + bb) + (lam + dd) * (lam + aa) * (lam + ee) - q1 - q2 * (2 * lam + aa + ee)
    lam <- lam - pl / dpl
    pl <- (lam + dd) * (lam + aa) * (lam + ee) * (lam + bb) - q1 * (lam + bb) - q2 * (lam + aa) * (lam + ee)
    dpl <- (lam + aa) * (lam + ee) * (lam + bb) + (lam + dd) * (lam + ee) * (lam + bb) + (lam + dd) * (lam + aa) * (lam + bb) + (lam + dd) * (lam + aa) * (lam + ee) - q1 - q2 * (2 * lam + aa + ee)
    lam <- lam - pl / dpl
    pl <- (lam + dd) * (lam + aa) * (lam + ee) * (lam + bb) - q1 * (lam + bb) - q2 * (lam + aa) * (lam + ee)
    dpl <- (lam + aa) * (lam + ee) * (lam + bb) + (lam + dd) * (lam + ee) * (lam + bb) + (lam + dd) * (lam + aa) * (lam + bb) + (lam + dd) * (lam + aa) * (lam + ee) - q1 - q2 * (2 * lam + aa + ee)
    lam <- lam - pl / dpl
    pl <- (lam + dd) * (lam + aa) * (lam + ee) * (lam + bb) - q1 * (lam + bb) - q2 * (lam + aa) * (lam + ee)
    dpl <- (lam + aa) * (lam + ee) * (lam + bb) + (lam + dd) * (lam + ee) * (lam + bb) + (lam + dd) * (lam + aa) * (lam + bb) + (lam + dd) * (lam + aa) * (lam + ee) - q1 - q2 * (2 * lam + aa + ee)
    lam <- lam - pl / dpl
    fitness <- lam

    # dr/du by implicit differentiation of P(r(u)) = 0 (equivalent to
    # Eq 33-34): dr/du = -(dP/du) / P'(r).
    rr <- fitness
    dpdr <- (rr + aa) * (rr + ee) * (rr + bb) + (rr + dd) * (rr + ee) * (rr + bb) + (rr + dd) * (rr + aa) * (rr + bb) + (rr + dd) * (rr + aa) * (rr + ee) - q1 - q2 * (2 * rr + aa + ee)
    dp_dfree <- (rr + dd) * (rr + ee) * (rr + bb) - kinf * gg * (rr + bb) - q2 * (rr + ee)
    dp_dc2c <- (rr + dd) * (rr + aa) * (rr + ee) - q1 - winf * h0 * (rr + aa) * (rr + ee)
    dp_dkinf <- -d_free * gg * (rr + bb)
    dp_dwinf <- -d_c2c * h0 * (rr + aa) * (rr + ee)
    ddfree_duk <- k_t * feff * esc_k * uk / sigma_free^2
    ddc2c_duw <- k_t * feff * esc_w * uw / sigma_c2c^2
    dkinf_duk <- -uk * kinf
    dwinf_duw <- -uw * winf
    drdu_free <- -(dp_dfree * ddfree_duk + dp_dkinf * dkinf_duk) / dpdr
    drdu_c2c <- -(dp_dc2c * ddc2c_duw + dp_dwinf * dwinf_duw) / dpdr

    # ---- 4. ODEs (Eq 31; Eq 5 is the special case u_k = u_w = 0) ------------
    tcell_uninfected(0) <- tcell0
    tcell_infected(0) <- tcell_infected0
    virus(0) <- virus0
    trait_free(0) <- trait_free0
    trait_c2c(0) <- trait_c2c0

    d/dt(tcell_uninfected) <- lambda_t - (kinf * virus + winf * tcell_infected + delta_t) * tcell_uninfected
    d/dt(tcell_early_free) <- kinf * virus * tcell_uninfected - (d_free + delta_t) * tcell_early_free
    d/dt(tcell_early_c2c) <- winf * tcell_infected * tcell_uninfected - (d_c2c + delta_c2c) * tcell_early_c2c
    d/dt(tcell_infected) <- d_free * tcell_early_free + d_c2c * tcell_early_c2c - (delta_burst + delta_t) * tcell_infected
    d/dt(virus) <- n_burst * delta_burst * tcell_infected - cl_virus * virus
    d/dt(trait_free) <- s_free * drdu_free
    d/dt(trait_c2c) <- s_c2c * drdu_c2c

    # ---- 5. Outputs ----------------------------------------------------------
    # Plasma viral load in copies/mL (Fig 7 and Fig 15-18 plot V x 1000:
    # 100 virions/uL at exposure is drawn at 1e5 vRNA/mL).
    vl <- 1000 * virus
  })
}
