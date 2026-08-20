Kang_2023_artesunate_hamster_pbpk <- function() {
  description <- paste(
    "Preclinical (golden hamster). PBPK (minimal, joint parent +",
    "metabolite; WinNonlin 8.3). Oral artesunate and its active",
    "metabolite dihydroartemisinin, each described by a five-compartment",
    "minimal PBPK model in which lung and trachea -- the antiviral target",
    "tissues -- are explicit perfusion-rate-limited organs and all",
    "nontarget tissues are lumped into a rest-of-body compartment (Kang",
    "et al. 2023, Pharmaceutics 15:838). Artesunate is treated as a",
    "prodrug that is completely hydrolysed to dihydroartemisinin in the",
    "blood compartment, so the parent apparent clearance CL/F is",
    "simultaneously the metabolite formation rate; the two subsystems are",
    "therefore one jointly fitted model and are kept in a single file.",
    "Plasma was the assay matrix and is recovered from blood via the",
    "fixed blood-to-plasma ratio 0.75 for both compounds. States carry",
    "molar amounts (nmol) so that the 1:1 artesunate to",
    "dihydroartemisinin conversion is mass-balanced. Physiological",
    "volumes and blood flows (Table 1) are fixed; the biochemical",
    "parameters (Table 2) were fitted to 60 artesunate and 127",
    "dihydroartemisinin pooled measurements by naive pooled-data",
    "analysis, so the model carries no between-subject variability and is",
    "intended for typical-value simulation. Two Table 2 / equation",
    "anomalies are reproduced verbatim -- see the vignette Errata."
  )
  reference <- paste(
    "Kang DW, Kim KM, Kim JH, Cho H-Y. Application of Minimal",
    "Physiologically-Based Pharmacokinetic Model to Simulate Lung and",
    "Trachea Exposure of Pyronaridine and Artesunate in Hamsters.",
    "Pharmaceutics. 2023;15(3):838. doi:10.3390/pharmaceutics15030838.",
    sep = " "
  )
  vignette <- "Kang_2023_pyronaridine_artesunate_hamster_pbpk"
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

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
    depot              = list(analyte = "artesunate", units = "nmol", specimen = "administration site", verified = FALSE),
    blood              = list(analyte = "artesunate", units = "nmol", specimen = "blood cell", verified = FALSE),
    lung               = list(analyte = "artesunate", units = "nmol", specimen = "tissue", verified = FALSE),
    trachea            = list(analyte = "artesunate", units = "nmol", specimen = "tissue", verified = FALSE),
    other              = list(analyte = "artesunate", units = "nmol", specimen = "tissue", verified = FALSE),
    blood_dihydroart   = list(analyte = "dihydroartemisinin", units = "nmol", specimen = "blood cell", verified = FALSE),
    lung_dihydroart    = list(analyte = "dihydroartemisinin", units = "nmol", specimen = "tissue", verified = FALSE),
    trachea_dihydroart = list(analyte = "dihydroartemisinin", units = "nmol", specimen = "tissue", verified = FALSE),
    other_dihydroart   = list(analyte = "dihydroartemisinin", units = "nmol", specimen = "tissue", verified = FALSE)
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
      "covers the artesunate / dihydroartemisinin arm (60 or 120 mg/kg",
      "artesunate)."
    ),
    regions = "Republic of Korea",
    notes = paste(
      "60 of the 530 total measurements were artesunate and 127 were",
      "dihydroartemisinin. Model fitted in WinNonlin 8.3 by a naive",
      "pooled-data approach, so no between-subject variability and no",
      "residual-error estimates are reported; see the vignette",
      "Assumptions and deviations section. Both compounds were",
      "essentially eliminated within about 2 h and did not accumulate on",
      "a 24 h dosing interval, so no steady state is reached and the",
      "paper reports tissue-to-blood AUC ratios rather than",
      "steady-state Cavg ratios."
    )
  )

  ini({
    # Fitted biochemical parameters of Kang 2023 Table 2, "Artesunate and
    # dihydroartemisinin" block. The paper reports point estimates only;
    # no standard errors, RSEs or confidence intervals are given, and the
    # normalized sensitivity coefficients (Figure 4b, 4c) are the only
    # parameter-precision diagnostic.

    # --- Artesunate (parent) ---
    lka <- log(1.74)
    label("Absorption rate constant ka (1/h)")                                        # Kang 2023 Table 2, artesunate (ka = 1.74 1/h)

    lcl <- log(2517.70)
    label("Apparent artesunate clearance CL/F from blood (L/h)")                      # Kang 2023 Table 2, artesunate (CL/F = 2517.70 L/h); this is simultaneously the dihydroartemisinin formation clearance

    lkp_lung <- log(10.33)
    label("Artesunate lung-to-blood partition coefficient Klung (unitless)")          # Kang 2023 Table 2, artesunate (Klung = 10.33). NOTE: this value does not reproduce the paper's own reported artesunate lung-to-blood AUC ratio of 3.34 and is numerically identical to CLm/F in the same table; see vignette Errata E2

    lkp_trachea <- log(1.48)
    label("Artesunate trachea-to-blood partition coefficient Ktrachea (unitless)")    # Kang 2023 Table 2, artesunate (Ktrachea = 1.48)

    lkp_other <- log(1.32)
    label("Artesunate rest-of-body-to-blood partition coefficient Krest (unitless)")  # Kang 2023 Table 2, artesunate (Krest = 1.32)

    lk_trachea_lung <- log(1.50)
    label("Artesunate trachea-to-lung transfer rate constant ktl (1/h)")              # Kang 2023 Table 2, artesunate (ktl = 1.50 1/h)

    lk_lung_trachea <- log(0.34)
    label("Artesunate lung-to-trachea transfer rate constant klt (1/h)")              # Kang 2023 Table 2, artesunate (klt = 0.34 1/h)

    # --- Dihydroartemisinin (active metabolite) ---
    lcl_dihydroart <- log(10.33)
    label("Apparent dihydroartemisinin clearance CLm/F (L/h)")                        # Kang 2023 Table 2, dihydroartemisinin (CLm/F = 10.33 L/h)

    lkp_lung_dihydroart <- log(0.34)
    label("Dihydroartemisinin lung-to-blood partition coefficient Klung,m (unitless)")           # Kang 2023 Table 2, dihydroartemisinin (Klung,m = 0.34)

    lkp_trachea_dihydroart <- log(1.08)
    label("Dihydroartemisinin trachea-to-blood partition coefficient Ktrachea,m (unitless)")     # Kang 2023 Table 2, dihydroartemisinin (Ktrachea,m = 1.08); used only by the blood equation -- the printed trachea equation uses Klung,m instead, see vignette Errata E1

    lkp_other_dihydroart <- log(1.21)
    label("Dihydroartemisinin rest-of-body-to-blood partition coefficient Krest,m (unitless)")   # Kang 2023 Table 2, dihydroartemisinin (Krest,m = 1.21)

    lk_trachea_lung_dihydroart <- log(6.98)
    label("Dihydroartemisinin trachea-to-lung transfer rate constant ktl,m (1/h)")    # Kang 2023 Table 2, dihydroartemisinin (ktl,m = 6.98 1/h)

    lk_lung_trachea_dihydroart <- log(0.35)
    label("Dihydroartemisinin lung-to-trachea transfer rate constant klt,m (1/h)")    # Kang 2023 Table 2, dihydroartemisinin (klt,m = 0.35 1/h)

    # Note on Ktrachea,m: it appears in the printed dihydroartemisinin
    # BLOOD equation (as Ctrachea,m * Qtrachea / Ktrachea,m, mirroring the
    # parent) but NOT in the printed dihydroartemisinin TRACHEA equation,
    # whose efflux term is written with Clung,m and Klung,m. Both are
    # reproduced exactly as printed, which means the metabolite trachea
    # efflux out of the trachea and the matching return into blood are
    # written with different coefficients. See vignette Errata E1.

    # Kang 2023 used a naive pooled-data least-squares fit in WinNonlin
    # 8.3 and reports no residual-error model and no variance parameters
    # of any kind. nlmixr2 model definitions require a residual-error
    # term, so the two propSd values below are fixed placeholders for
    # syntactic completeness only and must NOT be read as estimates.
    # Same convention as An_2012_mitoxantrone_mouse_pbpk. See the
    # vignette Assumptions and deviations section.
    propSd <- fixed(0.10)
    label("Artesunate proportional residual error placeholder (fraction)")            # not reported in Kang 2023; placeholder only
    propSd_dihydroart <- fixed(0.10)
    label("Dihydroartemisinin proportional residual error placeholder (fraction)")    # not reported in Kang 2023; placeholder only
  })

  model({
    # ===== Hamster physiology (Kang 2023 Table 1), fixed =====
    # Table 1 reports volumes in mL and flows in mL/h. They are written
    # here in L and L/h so that they are dimensionally consistent with
    # CL/F and CLm/F (L/h, Table 2) with no conversion factor, and so
    # that an amount in nmol divided by a volume in L gives nmol/L, the
    # unit the paper reports artesunate and dihydroartemisinin
    # concentrations in (Table 3 footnote 1).
    v_blood <- 7.20 / 1000        # Kang 2023 Table 1 (Vblood = 7.20 mL, ref 22)
    v_lung <- 0.48 / 1000         # Kang 2023 Table 1 (Vlung = 0.48 mL, experimental)
    v_trachea <- 0.06 / 1000      # Kang 2023 Table 1 (Vtrachea = 0.06 mL, experimental)
    v_other <- 94.27 / 1000       # Kang 2023 Table 1 (Vrest = 94.27 mL = Vtotal 102.01 - Vblood - Vlung - Vtrachea)
    q_co <- 1181.28 / 1000        # Kang 2023 Table 1 (Qco = 1181.28 mL/h, ref 23)
    q_trachea <- 24.81 / 1000     # Kang 2023 Table 1 (Qtrachea = 24.81 mL/h = 2.1% of Qco, refs 24-26)
    q_other <- 1156.47 / 1000     # Kang 2023 Table 1 (Qrest = 1156.47 mL/h = Qco - Qtrachea)

    # Blood-to-plasma partition coefficients, fixed from the literature
    # (Kang 2023 Table 1, ref 21). Both compounds share the value 0.75,
    # which is why plasma rather than whole blood was the assay matrix.
    kbp <- 0.75                   # Kang 2023 Table 1 (Kb:p artesunate = 0.75)
    kbp_dihydroart <- 0.75        # Kang 2023 Table 1 (Kb:p,m dihydroartemisinin = 0.75)

    # ===== Individual (typical) biochemical parameters =====
    ka <- exp(lka)
    cl <- exp(lcl)
    kp_lung <- exp(lkp_lung)
    kp_trachea <- exp(lkp_trachea)
    kp_other <- exp(lkp_other)
    ktl <- exp(lk_trachea_lung)
    klt <- exp(lk_lung_trachea)

    cl_dihydroart <- exp(lcl_dihydroart)
    kp_lung_dihydroart <- exp(lkp_lung_dihydroart)
    kp_trachea_dihydroart <- exp(lkp_trachea_dihydroart)
    kp_other_dihydroart <- exp(lkp_other_dihydroart)
    ktl_dihydroart <- exp(lk_trachea_lung_dihydroart)
    klt_dihydroart <- exp(lk_lung_trachea_dihydroart)

    # ===== Concentrations (nmol/L) =====
    # States hold amount in nmol; volumes are in L.
    c_blood <- blood / v_blood
    c_lung <- lung / v_lung
    c_trachea <- trachea / v_trachea
    c_other <- other / v_other

    c_blood_dihydroart <- blood_dihydroart / v_blood
    c_lung_dihydroart <- lung_dihydroart / v_lung
    c_trachea_dihydroart <- trachea_dihydroart / v_trachea
    c_other_dihydroart <- other_dihydroart / v_other

    # ===== Artesunate (parent) ODEs =====
    # Kang 2023 Section 2.3, "Parent-metabolite PBPK model for
    # artesunate" (page 5), transcribed term-by-term. Structurally
    # identical to the pyronaridine model.
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
    d/dt(other) <-
      c_blood * q_other -
      c_other * q_other / kp_other

    # ===== Dihydroartemisinin (metabolite) ODEs =====
    # Kang 2023 Section 2.3, "Parent-metabolite PBPK model for
    # dihydroartemisinin" (page 5). Formation is Cblood * CL/F, i.e. the
    # whole parent elimination flux appears as metabolite, matching the
    # paper's statement that artesunate is "thoroughly converted to
    # dihydroartemisinin in the blood compartment by the CL/F".
    d/dt(blood_dihydroart) <-
      c_blood * cl +
      c_lung_dihydroart * q_co / kp_lung_dihydroart +
      c_trachea_dihydroart * q_trachea / kp_trachea_dihydroart +
      c_other_dihydroart * q_other / kp_other_dihydroart -
      c_blood_dihydroart * q_co -
      c_blood_dihydroart * q_trachea -
      c_blood_dihydroart * q_other -
      c_blood_dihydroart * cl_dihydroart
    d/dt(lung_dihydroart) <-
      c_blood_dihydroart * q_co -
      c_lung_dihydroart * q_co / kp_lung_dihydroart +
      trachea_dihydroart * ktl_dihydroart -
      lung_dihydroart * klt_dihydroart
    # As printed, the efflux term of the dihydroartemisinin TRACHEA
    # equation is (Clung,m * Qtrachea) / Klung,m -- the lung
    # concentration and the LUNG partition coefficient, not the trachea
    # ones -- even though the matching return term in the metabolite
    # BLOOD equation above is written (Ctrachea,m * Qtrachea) /
    # Ktrachea,m, and even though the parallel pyronaridine and
    # artesunate trachea equations both use (Ctrachea * Qtrachea) /
    # Ktrachea. This reads as a typesetting slip, but it is the printed
    # form that reproduces the paper's reported dihydroartemisinin
    # trachea-to-blood AUC ratio of 0.15 (as-printed gives 0.136; the
    # "corrected" form gives 1.06), so it is reproduced verbatim.
    # See vignette Errata E1.
    d/dt(trachea_dihydroart) <-
      c_blood_dihydroart * q_trachea -
      c_lung_dihydroart * q_trachea / kp_lung_dihydroart +
      lung_dihydroart * klt_dihydroart -
      trachea_dihydroart * ktl_dihydroart
    d/dt(other_dihydroart) <-
      c_blood_dihydroart * q_other -
      c_other_dihydroart * q_other / kp_other_dihydroart

    # ===== Observations =====
    # Blood is the model's native matrix; plasma is recovered by
    # dividing by Kb:p, matching the paper's Cplasma = Cblood / Kb:p.
    Cc <- c_blood
    Cplasma <- c_blood / kbp
    Clung <- c_lung
    Ctrachea <- c_trachea

    Cc_dihydroart <- c_blood_dihydroart
    Cplasma_dihydroart <- c_blood_dihydroart / kbp_dihydroart
    Clung_dihydroart <- c_lung_dihydroart
    Ctrachea_dihydroart <- c_trachea_dihydroart

    Cc ~ prop(propSd)
    Cc_dihydroart ~ prop(propSd_dihydroart)
  })
}
