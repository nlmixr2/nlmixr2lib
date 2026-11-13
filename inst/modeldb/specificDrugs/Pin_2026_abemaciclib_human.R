Pin_2026_abemaciclib_human <- function() {
  description <- paste(
    "QSP. Human small-intestinal epithelium. Pin 2026 translation of the",
    "abemaciclib toxicity quantified in human small-intestinal organoids to",
    "clinical epithelial injury and diarrhea risk. The four-compartment",
    "crypt-villus epithelium of Gall 2023 (stem, proliferative and",
    "non-proliferative crypt cells, villus cells) with villus-to-crypt",
    "feedback and arrested/dead crypt and villus compartments is driven by",
    "the organoid kill function of abemaciclib (CDK4/6 inhibitor),",
    "evaluated at the free drug concentration in the crypt cells, supplied",
    "as the time-varying covariate CEFFECT (uM). The authors used the",
    "Simcyp-predicted free enterocyte concentration (ileum I, M-ADAM) as the",
    "surrogate, with the free plasma concentration as the alternative; the",
    "PBPK layer is not included. Villus injury drives a nonlinear",
    "multinomial-logistic model of the probability of diarrhea of any grade",
    "and of grade 3 or higher. Deterministic typical-value model (no IIV, no",
    "residual error). Time is in hours."
  )
  reference <- paste(
    "Pin C, Maheshvare M D, Gall L, Hall A, Coen M, Bianco A, Liu S,",
    "Pryde K. (2026). Intestinal Organoid-Based Mathematical Modeling",
    "Predicts Clinical Gastrointestinal Toxicity of Oral Oncology Drugs.",
    "CPT Pharmacometrics Syst Pharmacol 15(5):e70256.",
    "doi:10.1002/psp4.70256.",
    "Model code deposited as Supporting Information Data S2 and at BioModels",
    "(MODEL2212120003)."
  )
  vignette <- "Pin_2026_intestinal_organoid_gi_toxicity"

  paper_specific_compartments <- c(
    "stem",
    "prolif",
    "nonprolif",
    "villus",
    "carr",
    "cdead",
    "varr",
    "vdead"
  )

  units <- list(
    time = "h",
    dosing = "n/a (the crypt free drug exposure is supplied externally via the time-varying covariate CEFFECT; the model has no dosing compartment)",
    concentration = "uM for CEFFECT; cells (one-dimensional crypt-villus strip) for the epithelial states"
  )

  covariateData <- list(
    CEFFECT = list(
      description = "Time-varying free abemaciclib concentration in the proliferative crypt cells, the site of action of the kill (and arrest) functions.",
      units = "uM",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Not measurable directly. Pin 2026 bounded it by two surrogates",
        "predicted with Simcyp V24 for the clinical regimen: the free",
        "enterocyte concentration in the ileum I segment (M-ADAM), the",
        "preferred surrogate, and the free plasma concentration. Supply the",
        "profile at each record and interpolate linearly. The paper's",
        "symbol is Drug."
      ),
      source_name = "Drug"
    )
  )

  compartmentData <- list(
    stem = list(analyte = "crypt stem cells", units = "cells", specimen = "tissue", verified = TRUE),
    prolif = list(
      analyte = "proliferative crypt progenitor cells",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    nonprolif = list(
      analyte = "non-proliferative (differentiated) crypt cells",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    villus = list(analyte = "healthy villus epithelial cells", units = "cells", specimen = "tissue", verified = TRUE),
    carr = list(analyte = "drug-arrested crypt cells", units = "cells", specimen = "tissue", verified = TRUE),
    cdead = list(analyte = "dead crypt cells", units = "cells", specimen = "tissue", verified = TRUE),
    varr = list(analyte = "drug-arrested villus cells", units = "cells", specimen = "tissue", verified = TRUE),
    vdead = list(analyte = "dead villus cells", units = "cells", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human (small-intestinal epithelium)",
    n_subjects = "Typical-value translational model; not fitted to individual patients. The diarrhea sub-model was fitted to 13 published clinical-study arms of the six drugs (Pin 2026 Table 2).",
    n_studies = 11L,
    disease_state = "Patients with HR+/HER2- breast cancer (abemaciclib monotherapy)",
    dose_range = "Simulated clinical regimen: 150 mg and 200 mg orally twice daily",
    regions = "Translational (human small-intestine parameters from the literature; clinical diarrhea incidences from published trials)",
    notes = paste(
      "Clinical comparator for abemaciclib: 327 (150 mg BID) and 132 (200 mg BID) patients [Table 2, ref 34].",
      "The epithelial cell-dynamics parameters are human small-intestine",
      "literature values shared by all six Pin 2026 human models; the kill",
      "function is the abemaciclib hSIO estimate (Pin_2026_abemaciclib_organoid)",
      "carried over unchanged. The diarrhea sub-model was fitted once (SAS",
      "NLMIXED, binomial) to the diarrhea counts of all 13 arms."
    )
  )

  ini({
    # Human epithelial cell dynamics. Values from the authors' deposited
    # MATLAB code (Data S2, GI_human_Integr_Abema.m), which reproduces the
    # Table 2 villus injury; Table S3 prints a different set (see vignette).
    lambda_s <- fixed(log(2) / (21 / 12 * 29.7))
    label("Stem-cell proliferation rate lambda_S (1/h)") # Table S3 and deposited code: ln(2)/(21/12 x 29.7)
    lambda_p <- fixed(log(2) / 29.7)
    label("Proliferative-cell proliferation rate lambda_P (1/h)") # Table S3 and deposited code: ln(2)/29.7
    tau_dead <- fixed(log(2) / 12)
    label("Detachment rate of dead cells tau_Dead (1/h)") # deposited code tau_A = log(2)/12 (Table S3 prints 1/12)
    gam <- fixed(3.41394)
    label("Villus-to-crypt proliferation feedback strength gamma (unitless)") # deposited code gamma = 3.41394 (Table S3 prints 5.31)
    mm <- fixed(1.468566)
    label("Crypt-to-villus transfer shape coefficient m (unitless)") # deposited code m = 1.468566 (Table S3 prints 1.46)
    sss <- fixed(2)
    label("Steady-state stem-cell number s_ss (cells)") # deposited code S_ss = 2 (Table S3 prints 4)
    pss <- fixed(19.446)
    label("Steady-state proliferative-progenitor number p_ss (cells)") # deposited code P_ss = 19.446 (Table S3 prints 17.446)
    dss <- fixed(11.154)
    label("Steady-state non-proliferative crypt-cell number d_ss (cells)") # Table S3 and deposited code: 11.154
    vss <- fixed(4 * 32.6)
    label("Steady-state villus-cell number v_ss (cells)") # Table S3 and deposited code: 4 x 32.6

    # Abemaciclib toxicity function -- hSIO estimates (Table 1, Abemaciclib
    # column) carried unchanged into the human model.
    karrest <- fixed(0)
    label("Maximal drug-induced cell-cycle arrest k_arrest (unitless)") # Table 1, Abemaciclib column: 0
    c50arrest <- fixed(1)
    label("Free drug concentration for half-maximal arrest C_50arrest (uM)") # inactive while k_arrest = 0; value 1 is the deposited-code placeholder (GI_human_Integr_Abema.m), not reported in Table 1
    hillarrest <- fixed(1)
    label("Hill coefficient of the arrest function (unitless)") # inactive while k_arrest = 0; deposited-code placeholder 1
    kkill <- fixed(1.90)
    label("Maximal drug-induced killing of proliferative cells k_kill (unitless)") # Table 1, Abemaciclib column: 1.90
    c50kill <- fixed(3.51)
    label("Free drug concentration for half-maximal killing C_50kill (uM)") # Table 1, Abemaciclib column: 3.51
    hillkill <- fixed(1)
    label("Hill coefficient of the kill function (unitless)") # Table 1: 1, assumed

    # Probability of diarrhea (Eqs E16-E17). Estimates from the SAS NLMIXED
    # output in the Supporting Information, which match the deposited code
    # and reproduce Table 2; Table S3 prints different values.
    g0 <- 1.7535
    label("Logit intercept of grade-0 (no) diarrhea g_0 (unitless)") # Supporting Information SAS 'Parameter Estimates': g0 = 1.7535 (Table S3 prints 2.94)
    g2 <- 6.2109
    label("Logit intercept of grade <= 2 diarrhea g_2 (unitless)") # SAS 'Parameter Estimates': g12 = 6.2109 (Table S3 prints 7.047)
    r_inj <- -5.1074
    label("Slope of the diarrhea logits on the injury Hill term r (unitless)") # SAS 'Parameter Estimates': r_inj = -5.1074 (Table S3 prints -5.56)
    inj50_0 <- 0.4818
    label("Villus injury at half-maximal grade-0 Hill term I_0 (fraction)") # SAS 'Parameter Estimates': Dam0_50 = 0.4818 (Table S3 prints 0.2627)
    inj50_2 <- 0.3998
    label("Villus injury at half-maximal grade <= 2 Hill term I_2 (fraction)") # SAS 'Parameter Estimates': Dam12_50 = 0.3998 (Table S3 prints 0.22)
    hill_inj <- fixed(5)
    label("Hill coefficient of the injury terms h (unitless)") # Table S3 and SAS code: hill = 5, fixed
  })
  model({
    # Steady-state-derived transfer and shedding rates (Table S3).
    delta_sp <- lambda_s
    delta_pd <- lambda_p * (1 + gam / 2) + delta_sp * sss / pss
    delta_dv <- delta_pd * pss / dss
    tau_v <- delta_dv * dss / vss

    # Drug-induced arrest and killing of proliferative cells (Eqs 5-6).
    farrest <- karrest * CEFFECT^hillarrest / (CEFFECT^hillarrest + c50arrest^hillarrest)
    fkill <- kkill * CEFFECT^hillkill / (CEFFECT^hillkill + c50kill^hillkill)

    # Villus-to-crypt feedback and crypt-to-villus crowding. The crowding term
    # follows the deposited code; Eq E15 prints the arrested and dead cells in
    # the numerator instead (see vignette).
    feedback <- gam * vss / (vss + villus)
    crowd <- ((stem + prolif + nonprolif) / (sss + pss + dss + carr + cdead))^mm

    # Epithelial cell dynamics (Eqs E7-E14).
    d/dt(stem) <- lambda_s * stem * (1 - farrest - fkill) - delta_sp * stem * stem / sss
    d/dt(prolif) <- delta_sp * stem * stem / sss +
      lambda_p * prolif * (1 + feedback - farrest - fkill) - delta_pd * prolif
    d/dt(nonprolif) <- delta_pd * prolif - delta_dv * nonprolif * crowd
    d/dt(villus) <- delta_dv * nonprolif * crowd - tau_v * villus
    d/dt(carr) <- farrest * (lambda_s * stem + lambda_p * prolif) - delta_dv * carr * crowd
    d/dt(cdead) <- fkill * (lambda_s * stem + lambda_p * prolif) -
      delta_dv * cdead * crowd - tau_dead * cdead
    d/dt(varr) <- delta_dv * carr * crowd - tau_v * varr
    d/dt(vdead) <- delta_dv * cdead * crowd - tau_v * vdead - tau_dead * vdead

    stem(0) <- sss
    prolif(0) <- pss
    nonprolif(0) <- dss
    villus(0) <- vss

    # Villus injury, the fraction of missing villus cells (Eq E16).
    injury <- max(0, min(1, 1 - (villus + varr + vdead) / vss))

    # Cumulative-logit diarrhea model (Eq E17): P(grade <= k) for k = 0, 2.
    hill0 <- injury^hill_inj / (injury^hill_inj + inj50_0^hill_inj)
    hill2 <- injury^hill_inj / (injury^hill_inj + inj50_2^hill_inj)
    prob_diarrhea <- 1 - expit(g0 + r_inj * hill0)
    prob_diarrhea_grade3 <- 1 - expit(g2 + r_inj * hill2)
  })
}
