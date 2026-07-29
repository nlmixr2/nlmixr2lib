Brillac_2025_isatuximab <- function() {
  description <- "Two-compartment population PK model with linear elimination for isatuximab in pediatric and adult patients with relapsed/refractory acute leukemias (Brillac 2025)"
  reference <- "Brillac C, Semiond D, Oprea C, Baruchel A, Zwaan CM, Nguyen L. Selection of isatuximab dosing regimen in pediatric patients with leukemia using population pharmacokinetics. Cancer Chemother Pharmacol. 2025;95:116. doi:10.1007/s00280-025-04832-2"
  vignette <- "Brillac_2025_isatuximab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Body weight is the only retained covariate; allometric power scaling applied to CL, V1, Q, and V2 with reference 38 kg (population median in the pooled adult + pediatric analysis dataset; Brillac 2025 Table 1).",
      source_name        = "WT"
    )
  )

  population <- list(
    n_subjects     = 79L,
    n_studies      = 2L,
    age_range      = "1.4-74 years (pediatric 1.4-17.0; adults 16-74)",
    age_median     = "8.0 years (pediatric only; adults: not reported individually)",
    weight_range   = "8.8-108 kg (pediatric 8.8-108; adults 46-93)",
    weight_median  = "38 kg (pooled population median; pediatric median 32.5 kg)",
    sex_female_pct = 38.5,
    race_ethnicity = "not reported",
    disease_state  = "Relapsed/refractory acute leukemia: pediatric AML (n = 26), pediatric B-ALL (n = 27), pediatric T-ALL (n = 12), and adult R/R T-ALL or T-cell lymphoblastic lymphoma (n = 14)",
    dose_range     = "20 mg/kg IV infusion; QW for 4 doses then Q2W (ALL design: D1, 8, 15, 22, 29, 43, 57; AML design: D1, 8, 15 of C1 and optional C2)",
    regions        = "Multicenter (ISLAY NCT02999633 and ISAKIDS NCT03860844)",
    notes          = "Final pop-PK dataset of 674 plasma concentrations (555 pediatric + 119 adult; ~8 observations per patient). Pediatric subset baseline characteristics and weight/age distributions reported in Results / Patients section and Figure 1 of Brillac 2025; the only patient < 24 months in the analysis was a single 1.4-year-old infant. Sex split based on the 61.5% male statement for the 65 pediatric patients evaluable for PK; treated as the population value because adult-cohort sex is not separately reported. PopPK fit with Monolix 2021R2 SAEM (proportional error)."
  )

  ini({
    # Structural parameters - typical values at the population reference
    # weight of 38 kg (the pooled-dataset median; Brillac 2025 Table 1).
    # CL and Q are reported in L/h in Table 1; converted to L/day (x 24)
    # because this model keeps time in days. Brillac 2025 Results section
    # also reports the converted values directly: "typical CL was 0.00556
    # L/h (0.133 L/day) and typical Q was 0.0358 L/h (0.859 L/day)".
    lcl <- log(0.00556 * 24); label("Linear clearance from central compartment at 38 kg reference (CL, L/day)") # Brillac 2025 Table 1: CL = 0.00556 L/h
    lvc <- log(1.98);         label("Central volume of distribution at 38 kg reference (V1, L)")                # Brillac 2025 Table 1: V1 = 1.98 L
    lq  <- log(0.0358 * 24);  label("Intercompartmental clearance at 38 kg reference (Q, L/day)")               # Brillac 2025 Table 1: Q  = 0.0358 L/h
    lvp <- log(2.20);         label("Peripheral volume of distribution at 38 kg reference (V2, L)")             # Brillac 2025 Table 1: V2 = 2.20 L

    # Body-weight allometric exponents on each PK parameter (Brillac 2025
    # Table 1 and Results: "CL = 0.00556 * (WT/38)^0.833; V1 = 1.98 *
    # (WT/38)^0.821; Q = 0.0358 * (WT/38)^0.85 (fixed); V2 = 2.20 *
    # (WT/38)^0.72"). The exponent on Q was held fixed because it was
    # poorly estimated (RSE 136% on the free fit).
    e_wt_cl <- 0.833;        label("Allometric exponent of WT on CL (unitless)")              # Brillac 2025 Table 1: beta_CL_log(WT/MedWT) = 0.833
    e_wt_vc <- 0.821;        label("Allometric exponent of WT on V1 (unitless)")              # Brillac 2025 Table 1: beta_V1_log(WT/MedWT) = 0.821
    e_wt_q  <- fixed(0.85);  label("Allometric exponent of WT on Q (fixed; unitless)")        # Brillac 2025 Table 1: beta_Q_log(WT/MedWT)  = 0.85 (fixed)
    e_wt_vp <- 0.72;         label("Allometric exponent of WT on V2 (unitless)")              # Brillac 2025 Table 1: beta_V2_log(WT/MedWT) = 0.72

    # Inter-individual variability (Monolix exponential model:
    # parameter_i = TV * exp(eta_i), eta_i ~ N(0, omega^2)).
    # Brillac 2025 Methods state: "Use of an exponential model implies a
    # log-normal distribution for the parameters and omega is thus an
    # approximate coefficient of variation." The percentages reported in
    # Table 1 ("Interindividual variability (%)") are therefore the SD of
    # eta on the log scale; the variance entered here is omega^2 = (%/100)^2.
    # IIV is diagonal (no covariances) per Methods.
    etalcl ~ 0.388  # Brillac 2025 Table 1: omega(CL) = 62.3% -> 0.623^2
    etalvc ~ 0.163  # Brillac 2025 Table 1: omega(V1) = 40.4% -> 0.404^2
    etalq  ~ 0.257  # Brillac 2025 Table 1: omega(Q)  = 50.7% -> 0.507^2
    etalvp ~ 0.243  # Brillac 2025 Table 1: omega(V2) = 49.3% -> 0.493^2

    # Residual error (proportional model;
    # Cp_ij = Cp_pred,ij * (1 + epsilon_ij)).
    propSd <- 0.257; label("Proportional residual error (fraction)") # Brillac 2025 Table 1: residual error proportional = 25.7%
  })
  model({
    # Individual PK parameters with body-weight allometric scaling
    # to the population median weight of 38 kg.
    cl <- exp(lcl + etalcl) * (WT / 38)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 38)^e_wt_vc
    q  <- exp(lq  + etalq)  * (WT / 38)^e_wt_q
    vp <- exp(lvp + etalvp) * (WT / 38)^e_wt_vp

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment model with IV (central-compartment) dosing and
    # linear elimination from the central compartment. Brillac 2025
    # Methods give the published ODE system:
    #   dA1/dt = K * Dose + k21*A2 - (CL/V1 + k12) * A1
    #   dA2/dt = k12 * A1 - k21 * A2
    # with K an infusion rate constant; here the infusion is supplied via
    # the event table's RATE column on the central-compartment dose row.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Concentration: dose in mg, vc in L -> central / vc has units mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
