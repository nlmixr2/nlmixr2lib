Pan_2019_rituximab <- function() {
  description <- "Kinetic-pharmacodynamic (K-PD) model of rituximab on CD19+ B lymphocytes in paediatric autoimmune disease: a first-order decaying virtual dose drives a saturable (Emax) increase in the first-order loss of a CD19+ turnover pool, with methotrexate and cyclophosphamide comedication increasing Emax (no rituximab concentrations were measured)"
  reference <- paste(
    "Pan S, Yu H, Surti A, Cheng I, Marks SD, Brogan PA, Eleftheriou D,",
    "Standing JF (2019). Pharmacodynamics of rituximab on B lymphocytes in",
    "paediatric patients with autoimmune diseases.",
    "Br J Clin Pharmacol 85(8):1790-1797. doi:10.1111/bcp.13970.",
    "Parameter values are the final estimates printed in the NONMEM model file",
    "of the paper's supplementary material (Section 2), which agree with Table 2",
    "to its printed precision.",
    sep = " "
  )
  vignette <- "Pan_2019_rituximab"
  dosing <- "depot_kpd"
  units <- list(time = "day", dosing = "mg", concentration = "10^6 cells/L")

  # Model structure (Pan 2019 Eqs 1-2 and the supplement $DES block):
  #   DADT(1) = -KE * A(1)
  #   DG      = (EMAX * A(1)) / (ED50 + A(1))
  #   DADT(2) = KIN - KOUT * (1 + DG) * A(2)
  # A(1) is the residual rituximab dose (mg; every dose adds its amount to
  # A(1)), A(2) the CD19+ lymphocyte count (10^6/L). KIN = BSL * KOUT and
  # A_0(2) = BSL, so the count starts at the estimated baseline.
  compartmentData <- list(
    depot_kpd = list(
      analyte = "rituximab (virtual K-PD dose)",
      units = "mg",
      specimen = "not applicable",
      verified = TRUE
    ),
    bcell = list(analyte = "CD19+ B lymphocytes", units = "10^6 cells/L", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    CONMED_MTX = list(
      description = "Concomitant methotrexate (1 = yes, 0 = no)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant methotrexate)",
      notes = paste(
        "Multiplies Emax by (1 + 0.657) when 1 (supplement $PK: CCOV1 = THETA(8)",
        "if MTX == 1; EMAX = THETA(4) * (1 + CCOV) * (1 + CCOV1)). 25.6% of the",
        "cohort received methotrexate (Table 1). The paper does not say whether",
        "the flag was time-varying or per-patient.",
        sep = " "
      ),
      source_name = "MTX"
    ),
    CONMED_CYCLOPHOSPHAMIDE = list(
      description = "Concomitant cyclophosphamide (1 = yes, 0 = no)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant cyclophosphamide)",
      notes = paste(
        "Multiplies Emax by (1 + 0.381) when 1 (supplement $PK: CCOV = THETA(7)",
        "if CTX == 1). 43.6% of the cohort received cyclophosphamide (Table 1).",
        "The methotrexate and cyclophosphamide factors multiply when both are 1.",
        sep = " "
      ),
      source_name = "CTX"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Tested as an asymptotic decrease in kin with age; not significant (Pan 2019 Methods and Results)."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "Tested as a higher kin in females; not significant (Pan 2019 Methods and Results)."
    ),
    CONMED_AZA = list(
      description = "Concomitant azathioprine (1 = yes, 0 = no)",
      units = "(binary)",
      type = "binary",
      notes = "Tested on Emax and ED50; not significant (Pan 2019 Results). 33.3% of the cohort received it (Table 1)."
    ),
    CONMED_MYCOPHENOLATE = list(
      description = "Concomitant mycophenolate mofetil (1 = yes, 0 = no)",
      units = "(binary)",
      type = "binary",
      notes = "Tested on Emax and ED50; not significant (Pan 2019 Results). 48.7% of the cohort received it (Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 39,
    n_studies = 1,
    n_observations = 258,
    age_range = "3.9-18 years",
    age_median = "13.2 years",
    sex_female_pct = 87.2,
    disease_state = paste(
      "Paediatric autoimmune disease treated with rituximab: systemic lupus",
      "erythematosus (22), granulomatosis with polyangiitis (6), juvenile",
      "dermatomyositis (4), anti-SRP positive myositis (3), and one each of",
      "juvenile idiopathic arthritis, lupus overlap syndrome, anti-neutrophil",
      "cytoplasmic antibody-positive vasculitis (MPA) and polyarteritis nodosa",
      sep = " "
    ),
    dose_range = paste(
      "2 x 750 mg/m^2 IV on days 1 and 15, or 4 x 375 mg/m^2 IV on days 1, 8,",
      "15 and 22; per-infusion dose median 1000 mg (range 413-1000 mg)",
      sep = " "
    ),
    baseline_cd19 = "170 x 10^6/L median (range 10-940) (Table 1)",
    co_medication = paste(
      "All on prednisolone and methylprednisolone; methotrexate 25.6%,",
      "cyclophosphamide 43.6%, azathioprine 33.3%, mycophenolate 48.7%,",
      "tacrolimus 5.1% (Table 1)",
      sep = " "
    ),
    regions = "United Kingdom (Great Ormond Street Hospital, London; 2010-2016)",
    notes = paste(
      "Retrospective electronic-health-record data; no rituximab concentrations",
      "were measured. 108 of 258 CD19+ counts were reported as zero and were set",
      "to 5 x 10^6/L (half the lowest detectable count; M6 method). Demographics",
      "from Pan 2019 Table 1.",
      sep = " "
    )
  )

  ini({
    # Structural parameters: final estimates from the supplement NONMEM file
    # ($THETA, 'values set to final estimates'); Table 2 prints them rounded.
    lkel <- log(0.0359064); label("First-order decay rate of the rituximab effect (virtual dose) KE (1/day)") # Supplement $THETA(1) = 0.0359064; Table 2 ke 0.036 (RSE 22.7%)
    lrbase <- log(266.449); label("Baseline CD19+ lymphocyte count BSL (10^6/L)") # Supplement $THETA(2) = 266.449; Table 2 Baseline 266.4 (RSE 37.8%)
    lkout <- log(0.0160064); label("First-order loss rate of CD19+ lymphocytes KOUT (1/day)") # Supplement $THETA(3) = 0.0160064; Table 2 kout 0.02 (RSE 41.0%)
    lemax <- log(35.1538); label("Maximal fractional increase in CD19+ loss rate EMAX (unitless)") # Supplement $THETA(4) = 35.1538; Table 2 Emax 35.2 (RSE 11.8%)
    led50 <- log(0.811751); label("Residual rituximab dose giving half-maximal effect ED50 (mg)") # Supplement $THETA(5) = 0.811751; Table 2 ED50 0.81 mg (RSE 190.3%)

    e_conmed_cyclophosphamide_emax <- 0.381171; label("Fractional increase in EMAX with cyclophosphamide (unitless)") # Supplement $THETA(7) CCOV = 0.381171; Table 2 covCYC 0.38 (RSE 74.4%)
    e_conmed_mtx_emax <- 0.657277; label("Fractional increase in EMAX with methotrexate (unitless)") # Supplement $THETA(8) CCOV1 = 0.657277; Table 2 covMTX 0.66 (RSE 37.4%)

    # Exponential IIV; supplement $OMEGA variances. Table 2 BSV% = 100*sqrt(omega).
    etalkel ~ 0.54598 # Supplement $OMEGA BSV_KE = 0.54598; Table 2 BSV 73.88%
    etalrbase ~ 0.96308 # Supplement $OMEGA BSV_BSL = 0.96308; Table 2 BSV 98.13%
    etalkout ~ 0.803053 # Supplement $OMEGA BSV_KOUT = 0.803053; Table 2 BSV 89.61%

    # $ERROR: SD = SQRT(0.0001 + CV * IPRED**2) with $SIGMA 1 FIX, so THETA(6)
    # CV = 0.351255 is the proportional VARIANCE: propSd = sqrt(0.351255).
    # The 0.0001 additive variance is a fixed constant: addSd = sqrt(0.0001).
    propSd <- 0.592668; label("Proportional residual error on CD19+ count (fraction)") # Supplement $THETA(6) CV = 0.351255 -> sqrt = 0.5927; Table 2 'Prop residual error (cv%)' 59.3 (RSE 14.4%)
    addSd <- fixed(0.01); label("Additive residual error on CD19+ count (10^6/L)") # Supplement $ERROR constant 0.0001 (variance) -> sqrt = 0.01
  })

  model({
    kel <- exp(lkel + etalkel)
    rbase <- exp(lrbase + etalrbase)
    kout <- exp(lkout + etalkout)
    emax <- exp(lemax) *
      (1 + e_conmed_cyclophosphamide_emax * CONMED_CYCLOPHOSPHAMIDE) *
      (1 + e_conmed_mtx_emax * CONMED_MTX)
    ed50 <- exp(led50)
    kin <- rbase * kout

    # Virtual (K-PD) rituximab dose in mg; doses are given into depot_kpd.
    d/dt(depot_kpd) <- -kel * depot_kpd
    drugEffect <- emax * depot_kpd / (ed50 + depot_kpd)

    d/dt(bcell) <- kin - kout * (1 + drugEffect) * bcell
    bcell(0) <- rbase

    # Observation: the CD19+ lymphocyte count itself (10^6/L).
    bcell ~ add(addSd) + prop(propSd) + combined2()
  })
}
