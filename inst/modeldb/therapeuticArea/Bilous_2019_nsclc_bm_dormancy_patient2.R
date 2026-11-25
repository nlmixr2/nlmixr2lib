Bilous_2019_nsclc_bm_dormancy_patient2 <- function() {
  description <- paste(
    "Mechanistic disease model of primary tumor growth and brain-metastasis",
    "dissemination with dormancy in non-small cell lung cancer, individual",
    "fit to patient 2 (EGFR wild-type squamous cell carcinoma) of Bilous 2019. Time 0 is",
    "the first cancer cell. The primary tumor (in cells) grows by Gompertz",
    "kinetics until diagnosis and then follows a tumor growth inhibition",
    "model (exponential growth p, log-kill kkill that decays with half-life",
    "thalf_res). Brain metastases are born at rate mu * tumor^gamma (Iwata",
    "framework), stay dormant at one cell for tau_dorm days and then grow by",
    "the same Gompertz law as the untreated primary tumor; metastases are",
    "unaffected by treatment. The Gompertz parameters are computed from the",
    "primary tumor diameter at diagnosis (TUMSZ), the squamous-cell doubling",
    "time and a 1e12-cell carrying capacity. Outputs are the cumulative",
    "number of metastases born and the number visible on MRI (diameter >= 5",
    "mm). Single-patient deterministic fit: no IIV and no drug input.",
    sep = " "
  )
  reference <- paste(
    "Bilous M, Serdjebi C, Boyer A, Tomasini P, Pouypoudat C, Barbolosi D,",
    "Barlesi F, Chomy F, Benzekry S. Quantitative mathematical modeling of",
    "clinical brain metastasis dynamics in non-small cell lung cancer.",
    "Sci Rep. 2019;9:13018. doi:10.1038/s41598-019-49407-3.",
    "Equations are from Materials and Methods (equations 1-5 and the dormancy",
    "model); parameter values are from Table 2 (patient 2), Table 1 and",
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
        "Patient 2 value: 53.7 mm (Bilous 2019 Table 2, S_d)."
      ),
      source_name = "S_d"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1L,
    n_studies = 1L,
    disease_state = "unresectable EGFR wild-type squamous cell lung carcinoma (non-small cell lung cancer) with brain metastases",
    dose_range = "routine systemic therapy (agents not reported); therapy is not resolved by agent in the model",
    regions = "France (University Hospital of Marseille)",
    notes = paste(
      "Patient 2 of Bilous 2019, the validation patient: 11 primary tumor",
      "size measurements and 16 brain-metastasis size measurements over 4 MRI",
      "time points (Materials and Methods, Patient data). The paper reports no",
      "demographics for this patient."
    )
  )

  ini({
    # Gompertz calibration of the untreated primary tumor (Materials and
    # Methods, Pre-clinical primary tumor growth): carrying capacity
    # exp(alpha / beta) = 1e12 cells and the histology-dependent volume
    # doubling time at the diagnosis size.
    ltsmax <- fixed(log(1e12))
    label("Natural log of the Gompertz carrying capacity (log cells)") # Methods: carrying capacity exp(alpha0p/betap) = 10^12 cells
    tdouble <- fixed(104)
    label("Primary tumor volume doubling time at diagnosis, squamous cell carcinoma (day)") # Supplementary Table S1 'Average' squamous cell carcinoma = 104 days (Methods text attributes 104 days to undifferentiated carcinoma; Table S1 gives 91 for that column)
    cell_density <- fixed(1e6)
    label("Tumor cell density (cells/mm^3)") # Methods: '1 mm3 = 10^6 cells'
    # Table 2 prints 3 mm for both patients, but the Methods define the
    # visible size as 5 mm and the patient 2 fit (Figure S7A) bins visible
    # metastases from 5 mm. With 5 mm the model reproduces the Figure S7B
    # visible-count steps (first at about 11.5 months after diagnosis);
    # with 3 mm every step comes about 1.5 months early.
    diam_vis <- fixed(5)
    label("Imaging visibility threshold for a brain metastasis (mm diameter)") # Methods 'S_vis ... taken here to be 5 mm in diameter'; Figure S7A size bins start at 5 mm; Table 2 prints 3 mm

    # Tumor growth inhibition model of the primary tumor after diagnosis
    # (equation 1), fitted to the patient 2 primary tumor sizes.
    lp <- log(0.0101)
    label("Primary tumor exponential growth rate during treatment alpha1 (1/day)") # Table 2 alpha1 patient 2 = 0.0101 day-1 (RSE 5.28e-3%)
    lkkill <- log(0.0439)
    label("Primary tumor log-kill rate at treatment start kappa (1/day)") # Table 2 kappa patient 2 = 0.0439 day-1 (RSE 1803%, poorly identified per Results)
    lthalf_res <- log(35.4)
    label("Half-life of the treatment effect (resistance) t_res (day)") # Table 2 t_res patient 2 = 35.4 day (RSE 0.998%)

    # Dormancy model of brain metastasis, fitted to the patient 2 brain
    # metastasis sizes (equations 2-5 and the dormancy model).
    lmu <- log(1.02e-12)
    label("Cellular metastatic potential mu (1/(cell^gamma * day))") # Table 2 mu patient 2 = 1.02e-12 cell-1 day-1 (RSE 71.9%)
    gamma <- fixed(1)
    label("Exponent of the primary tumor size in the dissemination rate gamma (unitless)") # Table 2 gamma patient 2 = 1 (best of the grid 0, 0.4, 0.5, 2/3, 1; Methods)
    ltau_dorm <- log(171)
    label("Dormancy duration of a newborn brain metastasis tau (day)") # Table 2 tau patient 2 = 171 day (RSE 33.3%)

    # Additive residual error on the number of metastases larger than a
    # given size (Methods, y = f + sigma * eps). Not printed; computed from
    # the Table 1 dormancy objective SSE = 1.71 over N = 16 brain-metastasis
    # sizes and P = 3 fitted parameters as sqrt(SSE / (N - P)).
    addSd_metastases_visible <- 0.3627
    label("Additive residual SD on the number of visible brain metastases (count)") # Table 1 dormancy patient 2 SSE 1.71; sqrt(1.71 / (16 - 3))
  })
  model({
    # Spherical tumors with 1e6 cells/mm^3: cells = density * pi / 6 * d^3.
    s_diag <- cell_density * pi / 6 * TUMSZ^3
    s_vis <- cell_density * pi / 6 * diam_vis^3

    # Gompertz parameters of the untreated primary tumor. With
    # alpha0p = betap * ltsmax, the Methods doubling-time equation
    # DT(Sd) = -1/betap * ln((betap ln(2 Sd) - alpha0p) / (betap ln(Sd) - alpha0p))
    # has the closed-form solution below. Table 2 'Computed' values for
    # patient 2 (TUMSZ = 53.7 mm): alpha0p = 0.0858, betap = 3.10e-3 /day.
    beta_gomp <- -log((log(2 * s_diag) - ltsmax) / (log(s_diag) - ltsmax)) / tdouble
    alpha_gomp <- beta_gomp * ltsmax

    # Age of the primary tumor at diagnosis, inverting
    # Sp(t) = exp(alpha0p / betap * (1 - exp(-betap * t))) at Sp = Sd
    # (Table 2 'T_d' = 2.1 years for patient 2).
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
