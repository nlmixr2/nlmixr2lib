McKenna_2017_doxorubicin_mdamb468_pd <- function() {
  description <- "In vitro (MDA-MB-468 human triple-negative breast cancer cell line). Doxorubicin treatment-response (pharmacodynamic) model for cell number over ~30 days after a timed exposure (McKenna 2017 Eqs 4, 5A, 5B and 10). Cells grow logistically at proliferation rate kp; after treatment at t = 0 the growth rate is reduced either by a constant death rate kd,A (Eq 5A, apoptosis-like, state cells_a) or by a transient death rate kd,B * r * t * exp(1 - r * t) that allows regrowth (Eq 5B, mitotic-catastrophe-like, state cells_b), each with its own carrying capacity. The reported cell number is the AIC-weighted average cells = w_A * cells_a + w_B * cells_b (Eq 10). The paper fitted kd,A, kd,B, r, theta_A, theta_B and the weights separately for each concentration x exposure condition and showed that they collapse onto one curve per cell line when plotted against the maximum bound intracellular concentration CB,max. Here each of them is a piecewise-linear function of log2(CBMAX_DOX_NM) through 11 knots on the two-fold dilution lattice, digitised from Fig 5, Fig S7 and Fig S8 and averaged across the 6, 12 and 24 h exposure series: the maintainers' stand-in for the paper's local linear regression and binomial logistic weight model, whose coefficients are not printed. CBMAX_DOX_NM is supplied by the user, computed with the companion PK model McKenna_2017_doxorubicin_mdamb468_pk. CBMAX_DOX_NM = 0 gives the untreated control (no drug effect, control carrying capacity). Initial cell number is set by dosing the same count into both cells_a and cells_b at t = 0. Deterministic; no residual error. Sibling cell lines: McKenna_2017_doxorubicin_sum149pt_pd, McKenna_2017_doxorubicin_mdamb231_pd, McKenna_2017_doxorubicin_mdamb453_pd."
  reference <- "McKenna MT, Weis JA, Barnes SL, Tyson DR, Miga MI, Quaranta V, Yankeelov TE. A Predictive Mathematical Modeling Approach for the Study of Doxorubicin Treatment in Triple Negative Breast Cancer. Sci Rep. 2017;7:5725. doi:10.1038/s41598-017-05902-z"
  vignette <- "McKenna_2017_doxorubicin"

  units <- list(
    time = "h",
    dosing = "cells (initial cell number, dosed into both cells_a and cells_b at t = 0)",
    concentration = "cells per well (cell number; CBMAX_DOX_NM is in nM)"
  )

  # cells_a / cells_b are the two candidate solutions N_TC,A and N_TC,B that
  # Eq 10 averages; the averaged output uses the registered `cells` name.
  paper_specific_compartments <- c("cells_a", "cells_b")
  dosing <- c("cells_a", "cells_b")

  compartmentData <- list(
    cells_a = list(
      analyte = "MDA-MB-468 cells (Eq 5A solution)",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    cells_b = list(
      analyte = "MDA-MB-468 cells (Eq 5B solution)",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    CBMAX_DOX_NM = list(
      description = "Maximum DNA-bound intracellular doxorubicin concentration CB,max of the treatment",
      units = "nM",
      type = "continuous",
      reference_category = NULL,
      notes = "Summary statistic of one treatment condition, constant within a simulation. Computed by running McKenna_2017_doxorubicin_mdamb468_pk forward with the treatment's medium concentration time course and taking the bound concentration after the free pool has washed out; simulate at least 48 h past the end of exposure (for a single pulse, CB,max = 0.01368 nM per nM*h of medium AUC). The knots span 0.802 to 820.8 nM; values outside that range are clamped to the end knots. 0 = untreated control.",
      source_name = "CB,max"
    )
  )

  population <- list(
    species = "in vitro (MDA-MB-468 human triple-negative breast cancer cell line, H2B-mRFP labelled)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    disease_state = "triple-negative breast cancer (MDA-MB-468, Lehmann basal-like 1 subtype). About 2,500 H2B-mRFP-labelled cells per well of a 96-well plate, grown for at least 3 days before treatment, then imaged daily by fluorescence microscopy for at least 30 days; nuclei counted in ImageJ. Media refreshed every 3 days.",
    dose_range = "9.77-2500 nM (two-fold dilution series, nine concentrations) applied for 6, 12 or 24 h then removed by media replacement (medium AUC 59-60,000 nM*h); six replicate wells per condition plus untreated controls.",
    cell_line = "MDA-MB-468 (ATCC), transduced with histone H2B-mRFP (Addgene plasmid 18982).",
    notes = "kp and the control carrying capacity were fitted to pre-treatment and untreated-control counts (Fig S3) and printed in the Results; kp was then fixed while the drug-effect parameters were fitted per condition with a smoothness penalty across conditions (Eq 7; Supplement). Figs 5, S7 and S8 show those per-condition fits against the cell line's simulated CB,max."
  )

  ini({
    lkgrow   <- log(0.0118) ; label("Proliferation rate constant kp (1/h)")  # Results, Doxorubicin treatment response: kp = 1.18 x 10^-2 hr^-1 for MDA-MB-468
    lcellmax <- log(16400) ; label("Carrying capacity of untreated control wells theta (cells)")  # Results, Doxorubicin treatment response: theta = 1.64 x 10^4 for MDA-MB-468 (Fig S3)

    # Knot 0 is the 9.77 nM x 6 h condition (2500 nM / 2^8, medium AUC 58.6 nM*h); knot n is the condition with 2^n times its AUC. Its CB,max is the companion PK model's value, 0.01368 nM per nM*h x 58.59 nM*h, which matches the Fig 5 abscissae.
    cbmax_n00 <- fixed(0.8016) ; label("CB,max at knot 0 (nM)")  # derived: McKenna_2017_doxorubicin_mdamb468_pk CB,max for 9.77 nM x 6 h (back-solved from the Figs 3-4 printed CB,max values)

    # Step death rate kd,A (Eq 5A) -- Fig 5 kd,A panel, MDA-MB-468; knot n sits at CB,max = cbmax_n00 * 2^n
    kda_n00    <- 0.00243   ; label("Step death rate kd,A (Eq 5A) at knot 0 (1/h)")  # Fig 5 kd,A panel, 6 h point (digitised)
    kda_n01    <- 0.00279   ; label("Step death rate kd,A (Eq 5A) at knot 1 (1/h)")  # Fig 5 kd,A panel, mean of the 6 and 12 h points (digitised)
    kda_n02    <- 0.00455   ; label("Step death rate kd,A (Eq 5A) at knot 2 (1/h)")  # Fig 5 kd,A panel, mean of the 6, 12 and 24 h points (digitised)
    kda_n03    <- 0.00762   ; label("Step death rate kd,A (Eq 5A) at knot 3 (1/h)")  # Fig 5 kd,A panel, mean of the 6, 12 and 24 h points (digitised)
    kda_n04    <- 0.0108    ; label("Step death rate kd,A (Eq 5A) at knot 4 (1/h)")  # Fig 5 kd,A panel, mean of the 6, 12 and 24 h points (digitised)
    kda_n05    <- 0.0133    ; label("Step death rate kd,A (Eq 5A) at knot 5 (1/h)")  # Fig 5 kd,A panel, mean of the 6, 12 and 24 h points (digitised)
    kda_n06    <- 0.0131    ; label("Step death rate kd,A (Eq 5A) at knot 6 (1/h)")  # Fig 5 kd,A panel, mean of the 6, 12 and 24 h points (digitised)
    kda_n07    <- 0.0129    ; label("Step death rate kd,A (Eq 5A) at knot 7 (1/h)")  # Fig 5 kd,A panel, mean of the 6, 12 and 24 h points (digitised)
    kda_n08    <- 0.0125    ; label("Step death rate kd,A (Eq 5A) at knot 8 (1/h)")  # Fig 5 kd,A panel, mean of the 6, 12 and 24 h points (digitised)
    kda_n09    <- 0.0123    ; label("Step death rate kd,A (Eq 5A) at knot 9 (1/h)")  # Fig 5 kd,A panel, mean of the 12 and 24 h points (digitised)
    kda_n10    <- 0.0116    ; label("Step death rate kd,A (Eq 5A) at knot 10 (1/h)")  # Fig 5 kd,A panel, 24 h point (digitised)

    # Transient death-rate amplitude kd,B (Eq 5B) -- Fig 5 kd,B panel, MDA-MB-468; knot n sits at CB,max = cbmax_n00 * 2^n
    kdb_n00    <- 0.00258   ; label("Transient death-rate amplitude kd,B (Eq 5B) at knot 0 (1/h)")  # Fig 5 kd,B panel, 6 h point (digitised)
    kdb_n01    <- 0.00297   ; label("Transient death-rate amplitude kd,B (Eq 5B) at knot 1 (1/h)")  # Fig 5 kd,B panel, mean of the 6 and 12 h points (digitised)
    kdb_n02    <- 0.00514   ; label("Transient death-rate amplitude kd,B (Eq 5B) at knot 2 (1/h)")  # Fig 5 kd,B panel, mean of the 6, 12 and 24 h points (digitised)
    kdb_n03    <- 0.00979   ; label("Transient death-rate amplitude kd,B (Eq 5B) at knot 3 (1/h)")  # Fig 5 kd,B panel, mean of the 6, 12 and 24 h points (digitised)
    kdb_n04    <- 0.0129    ; label("Transient death-rate amplitude kd,B (Eq 5B) at knot 4 (1/h)")  # Fig 5 kd,B panel, mean of the 6, 12 and 24 h points (digitised)
    kdb_n05    <- 0.0165    ; label("Transient death-rate amplitude kd,B (Eq 5B) at knot 5 (1/h)")  # Fig 5 kd,B panel, mean of the 6, 12 and 24 h points (digitised)
    kdb_n06    <- 0.0166    ; label("Transient death-rate amplitude kd,B (Eq 5B) at knot 6 (1/h)")  # Fig 5 kd,B panel, mean of the 6, 12 and 24 h points (digitised)
    kdb_n07    <- 0.0163    ; label("Transient death-rate amplitude kd,B (Eq 5B) at knot 7 (1/h)")  # Fig 5 kd,B panel, mean of the 6, 12 and 24 h points (digitised)
    kdb_n08    <- 0.0155    ; label("Transient death-rate amplitude kd,B (Eq 5B) at knot 8 (1/h)")  # Fig 5 kd,B panel, mean of the 6, 12 and 24 h points (digitised)
    kdb_n09    <- 0.0154    ; label("Transient death-rate amplitude kd,B (Eq 5B) at knot 9 (1/h)")  # Fig 5 kd,B panel, mean of the 12 and 24 h points (digitised)
    kdb_n10    <- 0.0148    ; label("Transient death-rate amplitude kd,B (Eq 5B) at knot 10 (1/h)")  # Fig 5 kd,B panel, 24 h point (digitised)

    # Transient-effect rate constant r (Eq 5B) -- Fig 5 r panel, MDA-MB-468; knot n sits at CB,max = cbmax_n00 * 2^n
    ron_n00    <- 0.00993   ; label("Transient-effect rate constant r (Eq 5B) at knot 0 (1/h)")  # Fig 5 r panel, 6 h point (digitised)
    ron_n01    <- 0.00988   ; label("Transient-effect rate constant r (Eq 5B) at knot 1 (1/h)")  # Fig 5 r panel, mean of the 6 and 12 h points (digitised)
    ron_n02    <- 0.0101    ; label("Transient-effect rate constant r (Eq 5B) at knot 2 (1/h)")  # Fig 5 r panel, mean of the 6, 12 and 24 h points (digitised)
    ron_n03    <- 0.00681   ; label("Transient-effect rate constant r (Eq 5B) at knot 3 (1/h)")  # Fig 5 r panel, mean of the 6, 12 and 24 h points (digitised)
    ron_n04    <- 0.00468   ; label("Transient-effect rate constant r (Eq 5B) at knot 4 (1/h)")  # Fig 5 r panel, mean of the 6, 12 and 24 h points (digitised)
    ron_n05    <- 0.00296   ; label("Transient-effect rate constant r (Eq 5B) at knot 5 (1/h)")  # Fig 5 r panel, mean of the 6, 12 and 24 h points (digitised)
    ron_n06    <- 0.00303   ; label("Transient-effect rate constant r (Eq 5B) at knot 6 (1/h)")  # Fig 5 r panel, mean of the 6, 12 and 24 h points (digitised)
    ron_n07    <- 0.00277   ; label("Transient-effect rate constant r (Eq 5B) at knot 7 (1/h)")  # Fig 5 r panel, mean of the 6, 12 and 24 h points (digitised)
    ron_n08    <- 0.00252   ; label("Transient-effect rate constant r (Eq 5B) at knot 8 (1/h)")  # Fig 5 r panel, mean of the 6, 12 and 24 h points (digitised)
    ron_n09    <- 0.00241   ; label("Transient-effect rate constant r (Eq 5B) at knot 9 (1/h)")  # Fig 5 r panel, mean of the 12 and 24 h points (digitised)
    ron_n10    <- 0.0021    ; label("Transient-effect rate constant r (Eq 5B) at knot 10 (1/h)")  # Fig 5 r panel, 24 h point (digitised)

    # Carrying capacity under Eq 5A theta_A -- Fig S7 theta (5A) panel, MDA-MB-468; knot n sits at CB,max = cbmax_n00 * 2^n
    thetaa_n00 <- 27000     ; label("Carrying capacity under Eq 5A theta_A at knot 0 (cells)")  # Fig S7 theta (5A) panel, 6 h point (digitised)
    thetaa_n01 <- 27200     ; label("Carrying capacity under Eq 5A theta_A at knot 1 (cells)")  # Fig S7 theta (5A) panel, mean of the 6 and 12 h points (digitised)
    thetaa_n02 <- 28600     ; label("Carrying capacity under Eq 5A theta_A at knot 2 (cells)")  # Fig S7 theta (5A) panel, mean of the 6, 12 and 24 h points (digitised)
    thetaa_n03 <- 30300     ; label("Carrying capacity under Eq 5A theta_A at knot 3 (cells)")  # Fig S7 theta (5A) panel, mean of the 6, 12 and 24 h points (digitised)
    thetaa_n04 <- 28400     ; label("Carrying capacity under Eq 5A theta_A at knot 4 (cells)")  # Fig S7 theta (5A) panel, mean of the 6, 12 and 24 h points (digitised)
    thetaa_n05 <- 29700     ; label("Carrying capacity under Eq 5A theta_A at knot 5 (cells)")  # Fig S7 theta (5A) panel, mean of the 6, 12 and 24 h points (digitised)
    thetaa_n06 <- 28700     ; label("Carrying capacity under Eq 5A theta_A at knot 6 (cells)")  # Fig S7 theta (5A) panel, mean of the 6, 12 and 24 h points (digitised)
    thetaa_n07 <- 29500     ; label("Carrying capacity under Eq 5A theta_A at knot 7 (cells)")  # Fig S7 theta (5A) panel, mean of the 6, 12 and 24 h points (digitised)
    thetaa_n08 <- 29100     ; label("Carrying capacity under Eq 5A theta_A at knot 8 (cells)")  # Fig S7 theta (5A) panel, mean of the 6, 12 and 24 h points (digitised)
    thetaa_n09 <- 28100     ; label("Carrying capacity under Eq 5A theta_A at knot 9 (cells)")  # Fig S7 theta (5A) panel, mean of the 12 and 24 h points (digitised)
    thetaa_n10 <- 28500     ; label("Carrying capacity under Eq 5A theta_A at knot 10 (cells)")  # Fig S7 theta (5A) panel, 24 h point (digitised)

    # Carrying capacity under Eq 5B theta_B -- Fig S7 theta (5B) panel, MDA-MB-468; knot n sits at CB,max = cbmax_n00 * 2^n
    thetab_n00 <- 23900     ; label("Carrying capacity under Eq 5B theta_B at knot 0 (cells)")  # Fig S7 theta (5B) panel, 6 h point (digitised)
    thetab_n01 <- 23700     ; label("Carrying capacity under Eq 5B theta_B at knot 1 (cells)")  # Fig S7 theta (5B) panel, mean of the 6 and 12 h points (digitised)
    thetab_n02 <- 23300     ; label("Carrying capacity under Eq 5B theta_B at knot 2 (cells)")  # Fig S7 theta (5B) panel, mean of the 6, 12 and 24 h points (digitised)
    thetab_n03 <- 24100     ; label("Carrying capacity under Eq 5B theta_B at knot 3 (cells)")  # Fig S7 theta (5B) panel, mean of the 6, 12 and 24 h points (digitised)
    thetab_n04 <- 11800     ; label("Carrying capacity under Eq 5B theta_B at knot 4 (cells)")  # Fig S7 theta (5B) panel, mean of the 6, 12 and 24 h points (digitised)
    thetab_n05 <- 11500     ; label("Carrying capacity under Eq 5B theta_B at knot 5 (cells)")  # Fig S7 theta (5B) panel, mean of the 6, 12 and 24 h points (digitised)
    thetab_n06 <- 11400     ; label("Carrying capacity under Eq 5B theta_B at knot 6 (cells)")  # Fig S7 theta (5B) panel, mean of the 6, 12 and 24 h points (digitised)
    thetab_n07 <- 11200     ; label("Carrying capacity under Eq 5B theta_B at knot 7 (cells)")  # Fig S7 theta (5B) panel, mean of the 6, 12 and 24 h points (digitised)
    thetab_n08 <- 11400     ; label("Carrying capacity under Eq 5B theta_B at knot 8 (cells)")  # Fig S7 theta (5B) panel, mean of the 6, 12 and 24 h points (digitised)
    thetab_n09 <- 11600     ; label("Carrying capacity under Eq 5B theta_B at knot 9 (cells)")  # Fig S7 theta (5B) panel, mean of the 12 and 24 h points (digitised)
    thetab_n10 <- 11500     ; label("Carrying capacity under Eq 5B theta_B at knot 10 (cells)")  # Fig S7 theta (5B) panel, 24 h point (digitised)

    # AIC weight of Eq 5B w_B -- Fig S8 panel, MDA-MB-468; knot n sits at CB,max = cbmax_n00 * 2^n
    wb_n00     <- 1         ; label("AIC weight of Eq 5B w_B at knot 0 (unitless)")  # Fig S8 panel, 6 h point (digitised)
    wb_n01     <- 0.865     ; label("AIC weight of Eq 5B w_B at knot 1 (unitless)")  # Fig S8 panel, mean of the 6 and 12 h points (digitised)
    wb_n02     <- 1         ; label("AIC weight of Eq 5B w_B at knot 2 (unitless)")  # Fig S8 panel, mean of the 6, 12 and 24 h points (digitised)
    wb_n03     <- 1         ; label("AIC weight of Eq 5B w_B at knot 3 (unitless)")  # Fig S8 panel, mean of the 6, 12 and 24 h points (digitised)
    wb_n04     <- 1         ; label("AIC weight of Eq 5B w_B at knot 4 (unitless)")  # Fig S8 panel, mean of the 6, 12 and 24 h points (digitised)
    wb_n05     <- 0.183     ; label("AIC weight of Eq 5B w_B at knot 5 (unitless)")  # Fig S8 panel, mean of the 6, 12 and 24 h points (digitised)
    wb_n06     <- 0         ; label("AIC weight of Eq 5B w_B at knot 6 (unitless)")  # Fig S8 panel, mean of the 6, 12 and 24 h points (digitised)
    wb_n07     <- 0.00667   ; label("AIC weight of Eq 5B w_B at knot 7 (unitless)")  # Fig S8 panel, mean of the 6, 12 and 24 h points (digitised)
    wb_n08     <- 0.367     ; label("AIC weight of Eq 5B w_B at knot 8 (unitless)")  # Fig S8 panel, mean of the 6, 12 and 24 h points (digitised)
    wb_n09     <- 1         ; label("AIC weight of Eq 5B w_B at knot 9 (unitless)")  # Fig S8 panel, mean of the 12 and 24 h points (digitised)
    wb_n10     <- 1         ; label("AIC weight of Eq 5B w_B at knot 10 (unitless)")  # Fig S8 panel, 24 h point (digitised)
  })

  model({
    kgrow   <- exp(lkgrow)
    cellmax <- exp(lcellmax)

    # Position of the treatment on the knot lattice: knot n sits at
    # cbmax_n00 * 2^n, so s_cb = log2(CBMAX_DOX_NM / cbmax_n00), clamped to the
    # digitised range [0, 10]. h01..h10 are the hinge fractions of each
    # knot interval, giving a piecewise-linear interpolation in log2(CBMAX_DOX_NM).
    s_cb <- log(max(CBMAX_DOX_NM, 1e-12) / cbmax_n00) / log(2)
    s_cb <- min(max(s_cb, 0), 10)
    h01 <- min(max(s_cb - 0, 0), 1)
    h02 <- min(max(s_cb - 1, 0), 1)
    h03 <- min(max(s_cb - 2, 0), 1)
    h04 <- min(max(s_cb - 3, 0), 1)
    h05 <- min(max(s_cb - 4, 0), 1)
    h06 <- min(max(s_cb - 5, 0), 1)
    h07 <- min(max(s_cb - 6, 0), 1)
    h08 <- min(max(s_cb - 7, 0), 1)
    h09 <- min(max(s_cb - 8, 0), 1)
    h10 <- min(max(s_cb - 9, 0), 1)

    kd_a    <- kda_n00 +
      (kda_n01 - kda_n00) * h01 +
      (kda_n02 - kda_n01) * h02 +
      (kda_n03 - kda_n02) * h03 +
      (kda_n04 - kda_n03) * h04 +
      (kda_n05 - kda_n04) * h05 +
      (kda_n06 - kda_n05) * h06 +
      (kda_n07 - kda_n06) * h07 +
      (kda_n08 - kda_n07) * h08 +
      (kda_n09 - kda_n08) * h09 +
      (kda_n10 - kda_n09) * h10
    kd_b    <- kdb_n00 +
      (kdb_n01 - kdb_n00) * h01 +
      (kdb_n02 - kdb_n01) * h02 +
      (kdb_n03 - kdb_n02) * h03 +
      (kdb_n04 - kdb_n03) * h04 +
      (kdb_n05 - kdb_n04) * h05 +
      (kdb_n06 - kdb_n05) * h06 +
      (kdb_n07 - kdb_n06) * h07 +
      (kdb_n08 - kdb_n07) * h08 +
      (kdb_n09 - kdb_n08) * h09 +
      (kdb_n10 - kdb_n09) * h10
    ron     <- ron_n00 +
      (ron_n01 - ron_n00) * h01 +
      (ron_n02 - ron_n01) * h02 +
      (ron_n03 - ron_n02) * h03 +
      (ron_n04 - ron_n03) * h04 +
      (ron_n05 - ron_n04) * h05 +
      (ron_n06 - ron_n05) * h06 +
      (ron_n07 - ron_n06) * h07 +
      (ron_n08 - ron_n07) * h08 +
      (ron_n09 - ron_n08) * h09 +
      (ron_n10 - ron_n09) * h10
    theta_a <- thetaa_n00 +
      (thetaa_n01 - thetaa_n00) * h01 +
      (thetaa_n02 - thetaa_n01) * h02 +
      (thetaa_n03 - thetaa_n02) * h03 +
      (thetaa_n04 - thetaa_n03) * h04 +
      (thetaa_n05 - thetaa_n04) * h05 +
      (thetaa_n06 - thetaa_n05) * h06 +
      (thetaa_n07 - thetaa_n06) * h07 +
      (thetaa_n08 - thetaa_n07) * h08 +
      (thetaa_n09 - thetaa_n08) * h09 +
      (thetaa_n10 - thetaa_n09) * h10
    theta_b <- thetab_n00 +
      (thetab_n01 - thetab_n00) * h01 +
      (thetab_n02 - thetab_n01) * h02 +
      (thetab_n03 - thetab_n02) * h03 +
      (thetab_n04 - thetab_n03) * h04 +
      (thetab_n05 - thetab_n04) * h05 +
      (thetab_n06 - thetab_n05) * h06 +
      (thetab_n07 - thetab_n06) * h07 +
      (thetab_n08 - thetab_n07) * h08 +
      (thetab_n09 - thetab_n08) * h09 +
      (thetab_n10 - thetab_n09) * h10
    w_b     <- wb_n00 +
      (wb_n01 - wb_n00) * h01 +
      (wb_n02 - wb_n01) * h02 +
      (wb_n03 - wb_n02) * h03 +
      (wb_n04 - wb_n03) * h04 +
      (wb_n05 - wb_n04) * h05 +
      (wb_n06 - wb_n05) * h06 +
      (wb_n07 - wb_n06) * h07 +
      (wb_n08 - wb_n07) * h08 +
      (wb_n09 - wb_n08) * h09 +
      (wb_n10 - wb_n09) * h10

    # Untreated control: no drug effect and the control carrying capacity
    if (CBMAX_DOX_NM <= 0) {
      kd_a    <- 0
      kd_b    <- 0
      theta_a <- cellmax
      theta_b <- cellmax
      w_b     <- 1
    }
    w_a <- 1 - w_b

    # Eqs 5A and 5B for t >= 0 (t = 0 is the start of drug application)
    kd_t_a <- kd_a
    kd_t_b <- kd_b * ron * t * exp(1 - ron * t)

    # Eq 4 with each death-rate function
    d/dt(cells_a) <- (kgrow - kd_t_a) * cells_a * (1 - cells_a / theta_a)
    d/dt(cells_b) <- (kgrow - kd_t_b) * cells_b * (1 - cells_b / theta_b)

    # Eq 10: AIC-weighted average of the two solutions
    cells <- w_a * cells_a + w_b * cells_b
  })
}
