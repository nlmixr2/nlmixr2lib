Bilous_2019_nsclc_bm_dormancy_patient1 <- function() {
  description <- paste(
    "Mechanistic disease model of primary tumor growth and brain-metastasis",
    "dissemination with dormancy in non-small cell lung cancer, individual",
    "fit to patient 1 (EGFR-mutated adenocarcinoma) of Bilous 2019. Time 0 is",
    "the first cancer cell. The primary tumor (in cells) grows by Gompertz",
    "kinetics until diagnosis and then follows a tumor growth inhibition",
    "model (exponential growth p, log-kill kkill that decays with half-life",
    "thalf_res). Brain metastases are born at rate mu * tumor^gamma (Iwata",
    "framework), stay dormant at one cell for tau_dorm days and then grow by",
    "the same Gompertz law as the untreated primary tumor; metastases are",
    "unaffected by treatment. The Gompertz parameters are computed from the",
    "primary tumor diameter at diagnosis (TUMSZ), the adenocarcinoma doubling",
    "time and a 1e12-cell carrying capacity. Outputs are the cumulative",
    "number of metastases born and the number visible on MRI (diameter >= 3",
    "mm). Single-patient deterministic fit: no IIV and no drug input.",
    sep = " "
  )
  reference <- paste(
    "Bilous M, Serdjebi C, Boyer A, Tomasini P, Pouypoudat C, Barbolosi D,",
    "Barlesi F, Chomy F, Benzekry S. Quantitative mathematical modeling of",
    "clinical brain metastasis dynamics in non-small cell lung cancer.",
    "Sci Rep. 2019;9:13018. doi:10.1038/s41598-019-49407-3.",
    "Equations are from Materials and Methods (equations 1-5 and the dormancy",
    "model); parameter values are from Table 2 (patient 1), Table 1 and",
    "Supplementary Table S1.",
    sep = " "
  )
  vignette <- "Bilous_2019_nsclc_brain_metastasis"

  units <- list(
    time = "day",
    dosing = "none (no drug input; treatment acts from the diagnosis time td computed in model())",
    concentration = "cells (tumor size) and count (number of brain metastases); there is no drug concentration"
  )

  paper_specific_compartments <- c("metastases", "metastases_visible")

  compartmentData <- list(
    tumor = list(
      analyte = "primary lung tumor size",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    metastases = list(
      analyte = "cumulative number of brain metastases born, N(t)",
      units = "count",
      specimen = "not applicable",
      verified = TRUE
    ),
    metastases_visible = list(
      analyte = "number of brain metastases above the imaging visibility threshold",
      units = "count",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    TUMSZ = list(
      description = "Largest diameter of the primary lung tumor at diagnosis",
      units = "mm",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Single-lesion largest diameter on CT, held constant per subject.",
        "Converted to cells assuming a sphere and 1e6 cells/mm^3; it fixes",
        "the Gompertz parameters (with the histology doubling time and the",
        "carrying capacity) and the age of the tumor at diagnosis td.",
        "Patient 1 value: 36.0 mm (Bilous 2019 Table 2, S_d)."
      ),
      source_name = "S_d"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1L,
    n_studies = 1L,
    disease_state = "unresectable EGFR-mutated lung adenocarcinoma (non-small cell lung cancer) with brain metastases",
    dose_range = "erlotinib, then docetaxel, erlotinib re-challenge and gemcitabine; therapy is not resolved by agent in the model",
    regions = "France (Institut Bergonie, Bordeaux)",
    notes = paste(
      "Patient 1 of Bilous 2019, the model-development patient: 10 primary",
      "tumor size measurements and 47 brain-metastasis size measurements over",
      "6 MRI time points, follow-up to 48 months after diagnosis (Materials",
      "and Methods, Patient data). The paper reports no demographics for this",
      "patient."
    )
  )

  ini({
    # Gompertz calibration of the untreated primary tumor (Materials and
    # Methods, Pre-clinical primary tumor growth): carrying capacity
    # exp(alpha / beta) = 1e12 cells and the histology-dependent volume
    # doubling time at the diagnosis size.
    ltsmax <- fixed(log(1e12))
    label("Natural log of the Gompertz carrying capacity (log cells)") # Methods: carrying capacity exp(alpha0p/betap) = 10^12 cells
    tdouble <- fixed(201)
    label("Primary tumor volume doubling time at diagnosis, adenocarcinoma (day)") # Supplementary Table S1 'Average' adenocarcinoma = 201 days; Methods text 201 days
    cell_density <- fixed(1e6)
    label("Tumor cell density (cells/mm^3)") # Methods: '1 mm3 = 10^6 cells'
    diam_vis <- fixed(3)
    label("Imaging visibility threshold for a brain metastasis (mm diameter)") # Table 2 'Visibility threshold' = 3 mm (Fixed); Figure 5A 'Imaging detection limit (3 mm)'

    # Tumor growth inhibition model of the primary tumor after diagnosis
    # (equation 1), fitted to the patient 1 primary tumor sizes.
    lp <- log(5.72e-4)
    label("Primary tumor exponential growth rate during treatment alpha1 (1/day)") # Table 2 alpha1 patient 1 = 5.72e-4 day-1 (RSE 23.6%)
    lkkill <- log(4.46e-3)
    label("Primary tumor log-kill rate at treatment start kappa (1/day)") # Table 2 kappa patient 1 = 4.46e-3 day-1 (RSE 20.6%)
    lthalf_res <- log(149)
    label("Half-life of the treatment effect (resistance) t_res (day)") # Table 2 t_res patient 1 = 149 day (RSE 1.05e-3%)

    # Dormancy model of brain metastasis, fitted to the patient 1 brain
    # metastasis sizes (equations 2-5 and the dormancy model).
    lmu <- log(2.00e-12)
    label("Cellular metastatic potential mu (1/(cell^gamma * day))") # Table 2 mu patient 1 = 2.00e-12 cell-1 day-1 (RSE 2.83%)
    gamma <- fixed(1)
    label("Exponent of the primary tumor size in the dissemination rate gamma (unitless)") # Table 2 gamma patient 1 = 1 (best of the grid 0, 0.4, 0.5, 2/3, 1; Methods)
    ltau_dorm <- log(133)
    label("Dormancy duration of a newborn brain metastasis tau (day)") # Table 2 tau patient 1 = 133 day (RSE 4.23%)

    # Additive residual error on the number of metastases larger than a
    # given size (Methods, y = f + sigma * eps). Not printed; computed from
    # the Table 1 dormancy objective SSE = 4.93 over N = 47 brain-metastasis
    # sizes and P = 3 fitted parameters as sqrt(SSE / (N - P)).
    addSd_metastases_visible <- 0.3347
    label("Additive residual SD on the number of visible brain metastases (count)") # Table 1 dormancy patient 1 SSE 4.93; sqrt(4.93 / (47 - 3))
  })
  model({
    # Spherical tumors with 1e6 cells/mm^3: cells = density * pi / 6 * d^3.
    s_diag <- cell_density * pi / 6 * TUMSZ^3
    s_vis <- cell_density * pi / 6 * diam_vis^3

    # Gompertz parameters of the untreated primary tumor. With
    # alpha0p = betap * ltsmax, the Methods doubling-time equation
    # DT(Sd) = -1/betap * ln((betap ln(2 Sd) - alpha0p) / (betap ln(Sd) - alpha0p))
    # has the closed-form solution below. Table 2 'Computed' values for
    # patient 1 (TUMSZ = 36 mm): alpha0p = 0.0284, betap = 1.03e-3 /day.
    beta_gomp <- -log((log(2 * s_diag) - ltsmax) / (log(s_diag) - ltsmax)) / tdouble
    alpha_gomp <- beta_gomp * ltsmax

    # Age of the primary tumor at diagnosis, inverting
    # Sp(t) = exp(alpha0p / betap * (1 - exp(-betap * t))) at Sp = Sd
    # (Table 2 'T_d' = 5.3 years for patient 1).
    td <- -log(1 - log(s_diag) / ltsmax) / beta_gomp

    p <- exp(lp)
    kkill <- exp(lkkill)
    thalf_res <- exp(lthalf_res)
    res <- log(2) / thalf_res
    mu <- exp(lmu)
    tau_dorm <- exp(ltau_dorm)

    # Primary tumor: Gompertz before diagnosis, equation 1 (TGI) after.
    if (t < td) {
      rate_tumor <- alpha_gomp - beta_gomp * log(tumor)
    } else {
      rate_tumor <- p - kkill * exp(-res * (t - td))
    }
    tumor(0) <- 1
    d/dt(tumor) <- rate_tumor * tumor

    # Equations 2 and 4: dissemination rate d = mu * Sp^gamma and the
    # cumulative number of metastases N(t), each born at one cell.
    d/dt(metastases) <- mu * tumor^gamma

    # Dormancy model: f_dorm(t, s) = N(t - t(s) - tau), with t(s) the
    # Gompertz time for a metastasis (same alpha0, beta as the primary
    # tumor; Table 2) to grow from one cell to s. The visible count is
    # f_dorm(t, s_vis). N at the lagged time is accumulated from the
    # closed-form primary-tumor size at t - lag_vis (Gompertz, then the
    # integrated equation 1), which equals the tumor state solved above.
    t_vis <- -log(1 - beta_gomp / alpha_gomp * log(s_vis)) / beta_gomp
    lag_vis <- t_vis + tau_dorm
    t_lag <- t - lag_vis
    if (t_lag <= 0) {
      dissem_lag <- 0
    } else if (t_lag < td) {
      dissem_lag <- mu * exp(alpha_gomp / beta_gomp * (1 - exp(-beta_gomp * t_lag)))^gamma
    } else {
      dissem_lag <- mu * (s_diag * exp(p * (t_lag - td) - kkill / res * (1 - exp(-res * (t_lag - td)))))^gamma
    }
    d/dt(metastases_visible) <- dissem_lag

    # Derived outputs: primary tumor largest diameter (mm) and the time
    # after diagnosis (day).
    tumor_diam <- (6 * tumor / (pi * cell_density))^(1 / 3)
    time_after_diag <- t - td

    metastases_visible ~ add(addSd_metastases_visible)
  })
}
