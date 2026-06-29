Frohlich_2018_mRNA_translation <- function() {
  description <- "QSP (mechanistic, single-cell, in vitro). Two-state gene-expression model with ribosomal rate-limited translation (model ii of Frohlich et al. 2018) for eGFP / d2eGFP fluorescent reporter expression after mRNA transfection in HuH7 hepatoma cells, fit by multi-experiment nonlinear mixed-effects modelling (MEMOIR + AMICI)."
  reference   <- paste(
    "Frohlich F, Reiser A, Fink L, Woschee D, Ligon T,",
    "Theis FJ, Radler JO, Hasenauer J.",
    "Multi-experiment nonlinear mixed effect modeling of",
    "single-cell translation kinetics after transfection.",
    "npj Syst Biol Appl. 2018;4:42.",
    "doi:10.1038/s41540-018-0079-7.",
    "Final population parameter values for the chosen model",
    "(model ii, ribosomal translation) were obtained from",
    "the authors' Zenodo supplementary-code deposit",
    "(doi:10.5281/zenodo.1228899, results_ribo.mat;",
    "see vignette for the recovery details).",
    sep = " "
  )
  vignette    <- "Frohlich_2018_mRNA_translation"
  units       <- list(time = "hour", dosing = "normalized mRNA mass (m0 = 1)", concentration = "log fluorescence intensity (a.u.)")

  covariateData <- list(
    STUDY_d2eGFP = list(
      description        = "1 = cell transfected with destabilized eGFP (d2eGFP, ~6.6 h protein half-life via C-terminal PEST sequence); 0 = cell transfected with eGFP (~22.8 h protein half-life). Cohort indicator selecting between the two reporter constructs in the Frohlich 2018 multi-experiment NLME analysis. All structural and ribosomal-binding parameters are shared between cohorts; only the protein degradation rate (gamma_eGFP for STUDY_d2eGFP = 0 vs gamma_d2eGFP for STUDY_d2eGFP = 1) and its IIV variance differ.",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Cell-level (time-fixed) indicator. The two cohorts were measured in parallel using the same transfection / imaging / quantification pipeline; the only experimental difference is the mRNA construct (eGFP vs d2eGFP). Frohlich 2018 Methods (Plasmid vectors and mRNA production).",
      source_name        = "EXPERIMENT (1 = eGFP, 2 = d2eGFP in the source dataset)"
    )
  )

  population <- list(
    species       = "in vitro (human hepatoma HuH7 cell line)",
    n_subjects    = 630,
    n_studies     = 1,
    age_range     = "n/a (immortalised cell line; 4-h adhesion + 30-h imaging per channel)",
    disease_state = "Healthy proliferating cells; no perturbation other than mRNA lipoplex transfection.",
    dose_range    = "Single mRNA lipoplex addition (~0.5 ug/uL mRNA in OptiMEM with Lipofectamine(TM) 2000 at 2.5 uL per 1 ug mRNA), 1-h incubation followed by washout. The structural model is encoded on the normalised mRNA mass m0; one 'dose' = the per-cell mRNA mass entering the cytoplasm.",
    regions       = "Germany (LMU Munich / Helmholtz Zentrum Muenchen).",
    notes         = "Single-cell time-lapse fluorescence microscopy of micropatterned protein arrays (30 um x 30 um fibronectin squares on PLL-g-PEG passivated coverslips). >= 200 cells per experimental condition per replicate; the primary fit (results_ribo.mat, 14-Nov-2016) uses N_eGFP = 236 and N_d2eGFP = 394 single cells (sum 630). Imaging interval 10 min over 30 h; observations start ~2 h after lipoplex addition (eGFP grid: 2.17-32 h; d2eGFP grid: 2.0-32 h). The output is log fluorescence intensity y = log(GFP + offset), so the residual error is additive on the log scale (paper Methods, Data acquisition and quantitative image analysis). Parameter values are the population fixed-effect means (beta) and diagonal random-effect variances (D) from results_ribo.mat in the authors' Zenodo deposit (doi:10.5281/zenodo.1228899); MEMOIR uses log10(parameter) = beta + b internally, so all variances are converted to the natural-log scale used by nlmixr2 via Var[ln(p)] = ln(10)^2 * Var[log10(p)] = 5.302 * D_paper (see in-file comments). Offset has no inter-cell IIV in the source model (only a common fixed effect per the per-experiment phi mapping in experiments_transfection_ribo.m); the C_offset value in the deposited parameter vector is unconstrained and not used here. Residual error is reported in the source experiment definition as sigma_noise = 0.3 on the log-fluorescence observable (initialisation in experiments_transfection_ribo.m); per-cell sigmas are estimated as nuisance parameters in the inner MEMOIR likelihood (estim_sigma = true in optimize_transfection.m) and not retained at the population level, so 0.3 is the reported population-level value carried forward here -- see vignette Assumptions."
  )

  ini({
    # Population means -- Frohlich 2018 Zenodo deposit (doi:10.5281/zenodo.1228899),
    # results_ribo.mat, parameters_MEM.MS.par[, 1] (best of 200 multistarts, logPost = 95747.24).
    # Values stored as log10(parameter); converted here to natural log via ln() of the
    # 10^(stored) linear value. Source-trace cross-checks against paper Table S2:
    #   ln(2) / 0.8096 = 0.857 h  ~ paper's mRNA half-life 0.8 h.
    #   ln(2) / 0.0303 = 22.87 h  ~ paper's eGFP half-life 22.8 h.
    #   ln(2) / 0.1055 =  6.57 h  ~ paper's d2eGFP half-life 6.6 h.
    ldelta_mrna    <- log(0.80958);    label("mRNA degradation rate (1/h)")                                 # results_ribo.mat M_delta1: log10 = -0.0917
    lkdeg_egfp     <- log(0.03031);    label("eGFP protein degradation rate (1/h)")                          # results_ribo.mat M_pbeta:  log10 = -1.5184
    lkdeg_d2egfp   <- log(0.10546);    label("d2eGFP protein degradation rate (1/h)")                        # results_ribo.mat M_pbetad2: log10 = -0.9769
    lk2_m0_scale   <- log(6.198e8);    label("Translation rate x m0 x fluorescence scale (a.u./(h * complex_normunit))")  # results_ribo.mat M_k2_m0_scale: log10 = 8.7923
    lt0            <- log(0.87294);    label("Effective transfection-onset time (h after start of imaging)") # results_ribo.mat M_t0:     log10 = -0.0590
    loffset        <- log(8.1489);     label("Fluorescence background offset (a.u.)")                        # results_ribo.mat M_offset: log10 =  0.9111
    lk1_m0         <- log(2010.07);    label("Ribosome-mRNA binding rate constant x m0 (1/(h * ribosome_normunit))")  # results_ribo.mat M_k1_m0:  log10 = 3.3032
    lfracr0_m0     <- log(6.235e-7);   label("Initial free-ribosome concentration / m0 ratio (ribosome_normunit)")     # results_ribo.mat M_frac_R0_m0: log10 = -6.2052
    lk2            <- log(0.58577);    label("Ribosome-mRNA complex catalytic rate (translation + dissociation, 1/h)") # results_ribo.mat M_k2:     log10 = -0.2323

    # Inter-individual (cell-to-cell) variability.
    # MEMOIR parameterises diag(D) on the log10 scale via D = diag(exp(C)) (xi2D.m
    # `diag-matrix-logarithm` case); convert to natural-log variance using
    # Var[ln(p)] = ln(10)^2 * Var[log10(p)] = 5.30190 * exp(C_*).
    # Source: results_ribo.mat parameters_MEM.MS.par[, 1] entries 10-18.
    etaldelta_mrna   ~ 5.30190 * 0.15494   # exp(-1.8647) = 0.15494 -> 0.8214
    etalkdeg_egfp    ~ 5.30190 * 0.21267   # exp(-1.5480) = 0.21267 -> 1.1276
    etalkdeg_d2egfp  ~ 5.30190 * 0.06133   # exp(-2.7914) = 0.06133 -> 0.3252
    etalk2_m0_scale  ~ 5.30190 * 0.14599   # exp(-1.9242) = 0.14599 -> 0.7741
    etalt0           ~ 5.30190 * 0.02848   # exp(-3.5587) = 0.02848 -> 0.1510
    etalk1_m0        ~ 5.30190 * 10.5161   # exp( 2.3529) = 10.5161 -> 55.7556
    etalfracr0_m0    ~ 5.30190 * 0.01972   # exp(-3.9263) = 0.01972 -> 0.1045
    etalk2           ~ 5.30190 * 0.29135   # exp(-1.2332) = 0.29135 -> 1.5446
    # No IIV on offset: in the source per-experiment phi mapping
    # (experiments_transfection_ribo.m), offset uses beta(6) alone with no
    # random-effect b -- it is a population-common parameter only.

    # Residual error on the log-fluorescence observable.
    # Frohlich 2018 source code (project/models/experiments_transfection_ribo.m
    # lines 21-22, 45-46): Model.exp{s}.sym.sigma_noise = sym(0.3). Per-cell
    # sigmas are estimated as inner-loop nuisance parameters in MEMOIR
    # (estim_sigma = true in optimize_transfection.m); the deposited code
    # carries 0.3 as the population-level value. Wrapped in fixed() because
    # population-level sigma is not estimated.
    addSd_logfluor <- fixed(0.3); label("Additive SD on log(GFP + offset) (a.u. on log scale)")  # experiments_transfection_ribo.m sigma_noise = sym(0.3)
  })

  model({
    # Individual single-cell parameters (linear scale).
    delta_mrna  <- exp(ldelta_mrna + etaldelta_mrna)
    # Protein degradation rate: cohort-specific value AND cohort-specific IIV
    # (only the chosen cohort's eta is observed per cell; the other contributes
    # zero via the STUDY_d2eGFP gating but is still drawn from N(0, var)).
    lkdeg_gfp_i <- lkdeg_egfp    * (1 - STUDY_d2eGFP) + lkdeg_d2egfp     * STUDY_d2eGFP +
                   etalkdeg_egfp * (1 - STUDY_d2eGFP) + etalkdeg_d2egfp  * STUDY_d2eGFP
    kdeg_gfp    <- exp(lkdeg_gfp_i)
    k2_m0_scale <- exp(lk2_m0_scale + etalk2_m0_scale)
    t0          <- exp(lt0          + etalt0)
    offset      <- exp(loffset)
    k1_m0       <- exp(lk1_m0       + etalk1_m0)
    fracr0_m0   <- exp(lfracr0_m0   + etalfracr0_m0)
    k2          <- exp(lk2          + etalk2)

    # ODE system -- Frohlich 2018 model (ii), as encoded in
    # project/models/transfection_ribo_syms.m. State 'mrna' is the free
    # cytoplasmic mRNA (normalised by m0 = 1 at the bolus event); 'gfp' is the
    # protein concentration; 'ribo' is the free ribosome concentration. The
    # bound mRNA-ribosome complex is the conservation residual
    # (fracr0_m0 - ribo) -- ribosome total = free + bound = fracr0_m0.
    d/dt(mrna) <- -delta_mrna * mrna - k1_m0 * mrna * ribo + k2 * (fracr0_m0 - ribo)
    d/dt(gfp)  <-  k2_m0_scale * (fracr0_m0 - ribo) - kdeg_gfp * gfp
    d/dt(ribo) <- -k1_m0 * mrna * ribo + k2 * (fracr0_m0 - ribo)

    # Initial conditions: pre-transfection, only free ribosomes exist.
    mrna(0) <- 0
    gfp(0)  <- 0
    ribo(0) <- fracr0_m0

    # Bolus injection (Frohlich 2018 syms file: +dirac(t - t0) on the mrna ODE,
    # i.e. mrna jumps by 1 at t = t0 in the normalised model). Encoded as a
    # dosing event of amount 1 to the 'mrna' compartment in the source data;
    # the per-cell onset time t0 is realised via alag(mrna) = t0.
    alag(mrna) <- t0

    # Observable: log fluorescence intensity (paper: y = log(scale * GFP + offset);
    # in the identifiable transformed model the scale factor is absorbed into
    # k2_m0_scale, so the observed quantity is log(gfp + offset)).
    logfluor <- log(gfp + offset)
    logfluor ~ add(addSd_logfluor)
  })
}
