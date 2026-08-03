Kang_2023_pyronaridine_hamster_pbpk <- function() {
  description <- paste(
    "Preclinical (golden hamster). PBPK (minimal, five-compartment;",
    "WinNonlin 8.3). Oral pyronaridine disposition with the antiviral",
    "target tissues -- lung and trachea -- carried as explicit",
    "perfusion-rate-limited organ compartments and all nontarget tissues",
    "lumped into a rest-of-body compartment (Kang et al. 2023,",
    "Pharmaceutics 15:838). Whole blood is the sampling matrix because",
    "pyronaridine partitions strongly into erythrocytes (reported",
    "blood-to-plasma ratio 4.9-17.8), so no blood-to-plasma conversion is",
    "applied. Tissue uptake follows the well-stirred perfusion-limited",
    "form Q_tissue * (C_blood - C_tissue / K_tissue), and a first-order",
    "bidirectional exchange (k_tl, k_lt) links lung and trachea directly.",
    "Physiological volumes and blood flows (Table 1) are fixed; the seven",
    "biochemical parameters (Table 2) were fitted to 343 pooled",
    "pyronaridine measurements by naive pooled-data analysis, so the model",
    "carries no between-subject variability and is intended for",
    "typical-value simulation of blood, lung and trachea profiles after",
    "single or daily oral dosing."
  )
  reference <- paste(
    "Kang DW, Kim KM, Kim JH, Cho H-Y. Application of Minimal",
    "Physiologically-Based Pharmacokinetic Model to Simulate Lung and",
    "Trachea Exposure of Pyronaridine and Artesunate in Hamsters.",
    "Pharmaceutics. 2023;15(3):838. doi:10.3390/pharmaceutics15030838.",
    sep = " "
  )
  vignette <- "Kang_2023_pyronaridine_artesunate_hamster_pbpk"
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL")

  # No covariates: the Kang 2023 minimal PBPK model fixes every
  # physiological volume and blood flow at the study-cohort mean
  # (102.01 g hamster) and fits the biochemical parameters to pooled
  # data. Body weight is therefore baked into the Table 1 constants
  # rather than entering as a scaling covariate.
  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    depot   = list(analyte = "pyronaridine", units = "ug", specimen = "administration site", verified = FALSE),
    blood   = list(analyte = "pyronaridine", units = "ug", specimen = "whole blood", verified = FALSE),
    lung    = list(analyte = "pyronaridine", units = "ug", specimen = "tissue", verified = FALSE),
    trachea = list(analyte = "pyronaridine", units = "ug", specimen = "tissue", verified = FALSE),
    other   = list(analyte = "pyronaridine", units = "ug", specimen = "tissue", verified = FALSE)
  )

  covariateData <- list()

  population <- list(
    species = "hamster (male golden / Syrian, Janvier Labs)",
    n_subjects = 108L,
    n_studies = 1L,
    age_range = "adult",
    weight_range = "102.01 +/- 5.72 g (mean total body volume 102.01 mL)",
    sex_female_pct = 0,
    disease_state = paste(
      "Healthy (uninfected) male golden hamsters, 12 h dark-light cycle,",
      "23 +/- 3 degC, 55 +/- 15% relative humidity. Study conducted under",
      "IACUC protocol SP2021-14 (Shin Poong Pharm. Co., Ltd., Seoul,",
      "Republic of Korea). Sampling was destructive and sparse: 4 animals",
      "per blood / tissue time point, each animal bled 1-3 times within",
      "72 h, with lung and trachea harvested after the final bleed."
    ),
    dose_range = paste(
      "Oral pyronaridine-artesunate fixed-dose combination once daily for",
      "3 days: low-dose group (n = 60) 180/60 mg/kg and high-dose group",
      "(n = 48) 360/120 mg/kg pyronaridine/artesunate. This model file",
      "covers the pyronaridine arm."
    ),
    regions = "Republic of Korea",
    notes = paste(
      "343 of the 530 total measurements were pyronaridine (blood, lung,",
      "trachea). Model fitted in WinNonlin 8.3 by a naive pooled-data",
      "approach, so no between-subject variability and no residual-error",
      "estimates are reported; see the vignette Assumptions and deviations",
      "section. The minimal PBPK structure was adapted from the Jermain",
      "et al. ivermectin COVID-19 repurposing model (reference 11 of the",
      "paper)."
    )
  )

  ini({
    # All seven parameters below are the fitted biochemical parameters of
    # Kang 2023 Table 2, pyronaridine block. The paper reports point
    # estimates only -- no standard errors, RSEs or confidence intervals
    # are given for any biochemical parameter, and the normalized
    # sensitivity coefficients (Figure 4a) are the paper's only
    # parameter-precision diagnostic.

    lka <- log(0.03)
    label("Absorption rate constant ka (1/h)")                                  # Kang 2023 Table 2, pyronaridine (ka = 0.03 1/h)

    lcl <- log(0.21)
    label("Apparent total clearance CL/F from blood (L/h)")                     # Kang 2023 Table 2, pyronaridine (CL/F = 0.21 L/h)

    lkp_lung <- log(26.06)
    label("Lung-to-blood partition coefficient Klung (unitless)")               # Kang 2023 Table 2, pyronaridine (Klung = 26.06)

    lkp_trachea <- log(8.67)
    label("Trachea-to-blood partition coefficient Ktrachea (unitless)")         # Kang 2023 Table 2, pyronaridine (Ktrachea = 8.67)

    lkp_other <- log(5.25e-7)
    label("Rest-of-body-to-blood partition coefficient Krest (unitless)")       # Kang 2023 Table 2, pyronaridine (Krest = 5.25e-7)

    lk_trachea_lung <- log(1.01)
    label("First-order trachea-to-lung transfer rate constant ktl (1/h)")       # Kang 2023 Table 2, pyronaridine (ktl = 1.01 1/h)

    lk_lung_trachea <- log(0.92)
    label("First-order lung-to-trachea transfer rate constant klt (1/h)")       # Kang 2023 Table 2, pyronaridine (klt = 0.92 1/h)

    # Kang 2023 used a naive pooled-data least-squares fit in WinNonlin
    # 8.3 and reports no residual-error model and no variance parameters
    # of any kind. nlmixr2 model definitions require a residual-error
    # term, so propSd is a fixed placeholder for syntactic completeness
    # only and must NOT be read as an estimate. Same convention as
    # An_2012_mitoxantrone_mouse_pbpk. See the vignette Assumptions and
    # deviations section.
    propSd <- fixed(0.10)
    label("Proportional residual error placeholder (fraction)")                 # not reported in Kang 2023; placeholder only
  })

  model({
    # ===== Hamster physiology (Kang 2023 Table 1), fixed =====
    # Table 1 reports volumes in mL and flows in mL/h. They are written
    # here in L and L/h so that they are dimensionally consistent with
    # CL/F (L/h, Table 2) with no conversion factor, and so that an
    # amount in ug divided by a volume in L gives ug/L = ng/mL, the unit
    # the paper reports pyronaridine concentrations in.
    v_blood <- 7.20 / 1000        # Kang 2023 Table 1 (Vblood = 7.20 mL, ref 22)
    v_lung <- 0.48 / 1000         # Kang 2023 Table 1 (Vlung = 0.48 mL, experimental)
    v_trachea <- 0.06 / 1000      # Kang 2023 Table 1 (Vtrachea = 0.06 mL, experimental)
    v_other <- 94.27 / 1000       # Kang 2023 Table 1 (Vrest = 94.27 mL = Vtotal 102.01 - Vblood - Vlung - Vtrachea)
    q_co <- 1181.28 / 1000        # Kang 2023 Table 1 (Qco = 1181.28 mL/h, ref 23)
    q_trachea <- 24.81 / 1000     # Kang 2023 Table 1 (Qtrachea = 24.81 mL/h = 2.1% of Qco, refs 24-26)
    q_other <- 1156.47 / 1000     # Kang 2023 Table 1 (Qrest = 1156.47 mL/h = Qco - Qtrachea)

    # ===== Individual (typical) biochemical parameters =====
    ka <- exp(lka)
    cl <- exp(lcl)
    kp_lung <- exp(lkp_lung)
    kp_trachea <- exp(lkp_trachea)
    kp_other <- exp(lkp_other)
    ktl <- exp(lk_trachea_lung)
    klt <- exp(lk_lung_trachea)

    # ===== Concentrations (ng/mL = ug/L) =====
    # States hold amount in ug; volumes are in L.
    c_blood <- blood / v_blood
    c_lung <- lung / v_lung
    c_trachea <- trachea / v_trachea
    c_other <- other / v_other

    # ===== ODE system (Kang 2023 Section 2.3, pyronaridine block) =====
    # Transcribed term-by-term from the printed equations on page 4:
    #   d(Aa)/dt       = -(ka * Aa)
    #   d(Ablood)/dt   = (ka * Aa)
    #                    + Clung*Qco/Klung + Ctrachea*Qtrachea/Ktrachea
    #                    + Crest*Qrest/Krest
    #                    - Cblood*Qco - Cblood*Qtrachea - Cblood*Qrest
    #                    - Cblood*CL/F
    #   d(Alung)/dt    = Cblood*Qco - Clung*Qco/Klung
    #                    + Atrachea*ktl - Alung*klt
    #   d(Atrachea)/dt = Cblood*Qtrachea - Ctrachea*Qtrachea/Ktrachea
    #                    + Alung*klt - Atrachea*ktl
    #   d(Arest)/dt    = Cblood*Qrest - Crest*Qrest/Krest
    # Note that blood perfuses lung at the full cardiac output AND
    # trachea and rest-of-body in parallel, so the total flow leaving
    # blood is Qco + Qtrachea + Qrest. That is the topology the authors
    # printed and drew (Figure 1a) and it is reproduced verbatim here.
    d/dt(depot) <- -ka * depot
    d/dt(blood) <-
      ka * depot +
      c_lung * q_co / kp_lung +
      c_trachea * q_trachea / kp_trachea +
      c_other * q_other / kp_other -
      c_blood * q_co -
      c_blood * q_trachea -
      c_blood * q_other -
      c_blood * cl
    d/dt(lung) <-
      c_blood * q_co -
      c_lung * q_co / kp_lung +
      trachea * ktl -
      lung * klt
    d/dt(trachea) <-
      c_blood * q_trachea -
      c_trachea * q_trachea / kp_trachea +
      lung * klt -
      trachea * ktl
    # The paper does not print the definition Crest = Arest / Vrest, but
    # Vrest is tabulated (Table 1) and every other compartment uses the
    # same amount/volume relation, so it is the only reading consistent
    # with the rest of the equation set. See the vignette Errata.
    d/dt(other) <-
      c_blood * q_other -
      c_other * q_other / kp_other

    # ===== Observations =====
    # Cc is the whole-blood concentration (the paper's sampling matrix
    # for pyronaridine). Clung and Ctrachea are the target-tissue
    # concentrations reproduced in Figures 2 and 6.
    Cc <- c_blood
    Clung <- c_lung
    Ctrachea <- c_trachea

    Cc ~ prop(propSd)
  })
}
