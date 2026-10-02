# Population quantitative systems pharmacology (QSP) model of CD19-specific
# CAR-T cell immunotherapy (axicabtagene ciloleucel) in relapsed/refractory
# large B-cell non-Hodgkin lymphoma, published by Mueller-Schoell et al.
# 2021 (Cancers 13:2782). Five ODE species: naive (TN), central-memory
# (TCM), effector-memory (TEM) and terminally differentiated effector
# (TEff) CAR-T cell phenotypes plus CD19+ metabolic tumour volume. The four
# phenotypes follow a progressive-differentiation chain (TN -> TCM -> TEM ->
# TEff); TN/TCM/TEM expand upon tumour contact via a Michaelis-Menten term
# (Vmax1, KM1) with the RESPECTIVE T-cell concentration in the denominator,
# undergo homeostatic proliferation (kp1/kp2/kp3, fixed to literature) and
# apoptosis (ke1..ke4). CD19+ tumour volume grows logistically (k5, K0,
# both fixed / unidentifiable) and is killed by each phenotype (Vmax5,1..4,
# KM5). A two-subpopulation mixture on the baseline maximum expansion rate
# (reference vs low-expansion) plus two covariates on that rate -- a prior
# autologous stem cell transplant (ASCT, encoded HCT_PRIOR) and the day-7
# CD4+/CD8+ CAR-T cell ratio (CD4CD8_RATIO, on the reference population
# only) -- explain 2/3 of the interindividual variability in expansion.
#
# All structural equations are written out in the paper (Equations 5-11)
# and all final parameter estimates are in Table 1. Per-patient covariates
# (baseline metabolic tumour volume, ASCT, day-7 CD4/CD8 ratio, mixture
# class) are in Supplementary Table S1. Observed / predicted CAR-T cell
# kinetic parameters (Cmax, Tmax, AUC0-28d) are in Supplementary Table S2.
#
# The model has no explicit dosing event: the authors discounted the
# initial distribution phase and imputed a low initial concentration of
# 0.1 cells/uL per phenotype (Section 3.1.4), encoded here as the state
# initial condition init_cart. CD19+ tumour volume starts at the measured
# per-patient baseline metabolic tumour volume (TUM_VOL_TOTAL).

MuellerSchoell_2021_axicabtageneCiloleucel_qsp <- function() {
  description <- "QSP (cellular kinetics / tumour dynamics). Population quantitative systems pharmacology model of CD19-specific CAR-T cell immunotherapy (axicabtagene ciloleucel) in 19 relapsed/refractory large B-cell non-Hodgkin lymphoma patients. Five ODE species: naive (TN), central-memory (TCM), effector-memory (TEM) and effector (TEff) CAR-T cell phenotypes on a progressive-differentiation chain, plus CD19+ metabolic tumour volume. TN/TCM/TEM expand upon tumour contact via Michaelis-Menten terms; the tumour grows logistically and is killed by each phenotype. A two-class mixture (reference vs low-expansion subpopulation) on the baseline maximum expansion rate carries a prior-ASCT effect (both classes) and a day-7 CD4+/CD8+ CAR-T cell ratio effect (reference class only). Residual error is log-transform-both-sides (log-normal) per species."
  reference <- paste(
    "Mueller-Schoell A, Puebla-Osorio N, Michelet R, Green MR, Kuenkele A,",
    "Huisinga W, Strati P, Chasen B, Neelapu SS, Yee C, Kloft C (2021).",
    "Early Survival Prediction Framework in CD19-Specific CAR-T Cell",
    "Immunotherapy Using a Quantitative Systems Pharmacology Model.",
    "Cancers 13(11):2782. doi:10.3390/cancers13112782.",
    sep = " "
  )
  vignette <- "MuellerSchoell_2021_axicabtageneCiloleucel_qsp"

  # The four CAR-T cell phenotype states are paper-mechanistic QSP
  # compartments with no canonical analogue in
  # inst/references/compartment-names.md (see that file's "Paper-specific
  # compartments" section). `tumor` is canonical.
  paper_specific_compartments <- c("tn", "tcm", "tem", "teff")

  units <- list(
    time = "day",
    dosing = "cells/uL (imputed initial CAR-T cell concentration per phenotype; set as state initial condition, no explicit dose event)",
    concentration = "cells/uL (CAR-T cell density in blood); mL (CD19+ metabolic tumour volume)"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. The four CAR-T states are peripheral-blood cell
  # densities (cells/uL); the tumour state is a metabolic tumour VOLUME in
  # mL. verified = TRUE: analyte and specimen confirmed against the source.
  compartmentData <- list(
    tn = list(
      analyte = "naive CD19-specific CAR-T cells",
      units = "cells/uL",
      specimen = "blood cell",
      verified = TRUE
    ),
    tcm = list(
      analyte = "central-memory CD19-specific CAR-T cells",
      units = "cells/uL",
      specimen = "blood cell",
      verified = TRUE
    ),
    tem = list(
      analyte = "effector-memory CD19-specific CAR-T cells",
      units = "cells/uL",
      specimen = "blood cell",
      verified = TRUE
    ),
    teff = list(
      analyte = "effector CD19-specific CAR-T cells",
      units = "cells/uL",
      specimen = "blood cell",
      verified = TRUE
    ),
    tumor = list(analyte = "CD19+ metabolic tumour volume", units = "mL", specimen = "tumor", verified = TRUE)
  )

  covariateData <- list(
    TUM_VOL_TOTAL = list(
      description = "Baseline CD19+ metabolic tumour volume",
      units = "mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Per-subject baseline measurement by [18F]FDG PET-CT, used as the initial condition of the CD19+ tumour state (tumor(0) <- TUM_VOL_TOTAL). Values from Supplementary Table S1 (range 2.54-3555 mL; cohort median 85.7 mL). Held constant per individual.",
      source_name = "Baseline metabolic tumour volume"
    ),
    HCT_PRIOR = list(
      description = "Prior autologous stem cell transplant (ASCT) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no previous ASCT)",
      notes = "1 = the patient had a previous autologous stem cell transplant before CAR-T infusion, 0 = no previous ASCT (Supplementary Table S1 'Previous ASCT' column; 7 of 19 patients = 1). Enters the baseline maximum expansion rate of BOTH the reference and the low-expansion subpopulation via the fractional-change term (1 + ASCTVmax1 * HCT_PRIOR); the single ASCTVmax1 effect is shared across the two classes (Equations 10 and 11).",
      source_name = "Previous ASCT [yes/no]"
    ),
    CD4CD8_RATIO = list(
      description = "Ratio of CD4+ to CD8+ CAR-T cells measured on day 7 post-infusion",
      units = "(dimensionless)",
      type = "continuous",
      reference_category = NULL,
      notes = "Day-7 peripheral-blood ratio of CD4+ to CD8+ CAR-T cells (Supplementary Table S1 'ratio of CAR-T cells at day seven'; cohort range 0.0512-8.94). Enters ONLY the reference-population baseline maximum expansion rate as a power term CD4CD8_RATIO^CD4CD8exp, referenced to a ratio value of 1 (Equation 10, Table 1). The paper measured the ratio at several time points but implemented day 7 because it was available for 18 of 19 patients. Not influential on the low-expansion subpopulation (Equation 11).",
      source_name = "CAR+ CD4/CD8 day7"
    ),
    MIX_LOW_EXPANSION = list(
      description = "Latent mixture-model class indicator: low-expansion vs reference-expansion subpopulation",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (reference-expansion population)",
      notes = "Per-subject latent-class index of the paper's two-class $MIXTURE on the baseline maximum expansion rate Vmax1,base. 1 = low-expansion subpopulation (reduced Vmax1,base, no CD4/CD8 effect, no retained IIV; Equation 11), 0 = reference-expansion population (Equation 10). The estimated population proportion of the reference class is MIXP = 0.803 (Table 1, RSE 11%), so the low-expansion class probability is 0.197 (4 of 19 patients). For typical-value simulation set MIX_LOW_EXPANSION = 0 (dominant reference phenotype); for population simulation draw MIX_LOW_EXPANSION ~ Bernoulli(0.197) per subject.",
      source_name = "Low expansion-subpopulation [yes/no]"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 19L,
    n_studies = 1L,
    age_range = "21-80 years (10-year age bands, Supplementary Table S1)",
    sex_female_pct = 31.6,
    disease_state = "relapsed/refractory large B-cell non-Hodgkin lymphoma (DLBCL, PMBCL, transformed follicular lymphoma)",
    dose_range = "single IV infusion of axicabtagene ciloleucel, target 2e6 CAR-T cells/kg (preceded by fludarabine + cyclophosphamide lymphodepletion on days -5/-4/-3)",
    regions = "USA (MD Anderson Cancer Center)",
    notes = "24 patients treated; 5 excluded from the model analysis (1 with no active disease, 2 missing baseline metabolic tumour volume, 2 with failed flow-cytometry staining), leaving n = 19. Baseline demographics and per-patient covariates in Supplementary Table S1."
  )

  ini({
    # -----------------------------------------------------------------
    # Baseline maximum expansion rate per mL tumour volume, per (sub)pop
    # (cells/uL) * day^-1 * mL^-1 -- Table 1, log-transformed
    # -----------------------------------------------------------------
    lvmax1_base_ref <- log(0.00846); label("Baseline max expansion rate per mL tumour volume, reference population (cells/uL/day/mL)")  # Table 1 'Vmax1,base,ref' = 0.00846 (RSE 36%)
    lvmax1_base_low <- log(0.000700); label("Baseline max expansion rate per mL tumour volume, low-expansion subpopulation (cells/uL/day/mL)")  # Table 1 'Vmax1,base,low' = 0.000700 (RSE 17%)

    # Covariate effects on the baseline maximum expansion rate
    e_hct_prior_vmax1 <- 2.53; label("Fractional change in Vmax1,base due to a previous ASCT (unitless)")  # Table 1 'ASCTVmax1' = 2.53 (RSE 31%); enters (1 + e_hct_prior_vmax1 * HCT_PRIOR)
    e_cd4cd8_vmax1 <- -0.385; label("Power exponent for day-7 CD4+/CD8+ CAR-T cell ratio on Vmax1,base,ref (unitless)")  # Table 1 'CD4/CD8 exp' = -0.385 (RSE 45%); referenced to a ratio of 1

    # Michaelis-Menten half-max expansion concentration (TN/TCM/TEM)
    lkm1 <- log(1.13); label("TN/TCM/TEM concentration at half-maximum expansion (cells/uL)")  # Table 1 'KM1' = 1.13 (RSE 22%)

    # Differentiation rate constants (1/day)
    lk12 <- log(0.140); label("Differentiation rate constant TN -> TCM (1/day)")  # Table 1 'k12' = 0.140 (RSE 9%)
    lk23 <- log(0.191); label("Differentiation rate constant TCM -> TEM (1/day)")  # Table 1 'k23' = 0.191 (RSE 11%)
    lk34 <- log(0.355); label("Differentiation rate constant TEM -> TEff (1/day)")  # Table 1 'k34' = 0.355 (RSE 13%)

    # Death rate constants (1/day)
    lke4 <- log(0.518); label("Death rate constant for TEff (1/day)")  # Table 1 'ke4' = 0.518 (RSE 13%)
    # ke1/ke2/ke3 fixed = 0.0104 = ke4 * 0.02, the long-lived:short-lived
    # death-rate ratio from Stein et al. 2019 [ref 21]; unidentifiable from
    # the rapid-expansion-phase data (Section 3.1.4). Table 1 rows ke1/ke2/ke3.
    ke1 <- fixed(0.0104); label("Death rate constant for TN (1/day)")  # Table 1 'ke1' = 0.0104 (derived from ke4)
    ke2 <- fixed(0.0104); label("Death rate constant for TCM (1/day)")  # Table 1 'ke2' = 0.0104 (derived from ke4)
    ke3 <- fixed(0.0104); label("Death rate constant for TEM (1/day)")  # Table 1 'ke3' = 0.0104 (derived from ke4)

    # Homeostatic proliferation rate constants (1/day), fixed to literature
    # values [ref 47]; unidentifiable during the rapid expansion phase.
    kp1 <- fixed(0.0005); label("Homeostatic proliferation rate constant for TN (1/day)")  # Table 1 'kp1' = 0.0005 (literature [47])
    kp2 <- fixed(0.007); label("Homeostatic proliferation rate constant for TCM (1/day)")  # Table 1 'kp2' = 0.007 (literature [47])
    kp3 <- fixed(0.007); label("Homeostatic proliferation rate constant for TEM (1/day)")  # Table 1 'kp3' = 0.007 (literature [47])

    # Maximum tumour-killing rates by each phenotype
    # [mL * day^-1 * (cells/uL)^-1]. Only Vmax5,2 (TCM) was estimated; the
    # others are fixed at absolute values derived from the Vmax5,2 estimate
    # and digitised in-vitro killing-capacity ratios (Schmueck-Henneresse
    # et al. [ref 49]); Table 1.
    lvmax5_2 <- log(4.04); label("Maximum tumour-killing rate by TCM (mL/day per cells/uL)")  # Table 1 'Vmax5,2' = 4.04 (RSE 39%)
    vmax5_1 <- fixed(2.57); label("Maximum tumour-killing rate by TN (mL/day per cells/uL)")  # Table 1 'Vmax5,1' = 2.57 (derived from Vmax5,2)
    vmax5_3 <- fixed(3.78); label("Maximum tumour-killing rate by TEM (mL/day per cells/uL)")  # Table 1 'Vmax5,3' = 3.78 (derived from Vmax5,2)
    vmax5_4 <- fixed(4.24); label("Maximum tumour-killing rate by TEff (mL/day per cells/uL)")  # Table 1 'Vmax5,4' = 4.24 (derived from Vmax5,2)

    lkm5 <- log(276); label("Metabolic tumour volume at half-maximum killing rate (mL)")  # Table 1 'KM5' = 276 (RSE 33%)

    # Tumour logistic-growth parameters, fixed (not identifiable: data
    # contained tumour volumes only in the presence of CAR-T cells).
    k5 <- fixed(0.0023); label("Proliferation rate constant of metabolic tumour volume (1/day)")  # Table 1 'k5' = 0.0023 (fixed)
    k0 <- fixed(5000); label("Maximum tumour volume observable / carrying capacity (mL)")  # Table 1 'K0' = 5000 (fixed)

    # Imputed initial CAR-T cell concentration per phenotype (Section 3.1.4;
    # a 10-fold change had minor impact only on Tmax, Supplementary Fig S5).
    init_cart <- fixed(0.1); label("Imputed initial CAR-T cell concentration per phenotype (cells/uL)")  # Section 3.1.4 imputed dose = 0.1 cells/uL

    # -----------------------------------------------------------------
    # Interindividual variability (exponential; reported as %CV in Table 1)
    # Variance on the log scale = log(1 + CV^2).
    # -----------------------------------------------------------------
    etalvmax1_base_ref ~ 1.17865  # Table 1 'IIV Vmax1,base,ref' = 150% CV (RSE 19%); log-scale variance = log(1 + 1.50^2); reference class only
    etalvmax5_2 ~ 2.3442  # Table 1 'IIV Vmax5,2' = 307% CV (RSE 19%); log-scale variance = log(1 + 3.07^2)

    # -----------------------------------------------------------------
    # Residual unexplained variability -- log-transform-both-sides
    # (log-normal). Reported as %CV in Table 1; log-scale SD = sqrt(log(1 + CV^2)).
    # -----------------------------------------------------------------
    expSd_pred_tn <- 0.547332; label("Residual SD for TN, log scale (LTBS)")  # Table 1 'RUV TN' = 59.1% CV (RSE 11%)
    expSd_pred_tcm <- 0.743415; label("Residual SD for TCM, log scale (LTBS)")  # Table 1 'RUV TCM' = 85.9% CV (RSE 9%)
    expSd_pred_tem <- 0.944456; label("Residual SD for TEM, log scale (LTBS)")  # Table 1 'RUV TEM' = 120% CV (RSE 9%)
    expSd_pred_teff <- 0.635942; label("Residual SD for TEff, log scale (LTBS)")  # Table 1 'RUV TEff' = 70.6% CV (RSE 10%)
    expSd_pred_tumor <- 0.917957; label("Residual SD for CD19+ tumour volume, log scale (LTBS)")  # Table 1 'RUV CD19+ tumour' = 115% CV (RSE 12%)
  })

  model({
    # -----------------------------------------------------------------
    # 1. Mixture-class selector (1 = low-expansion subpopulation)
    # -----------------------------------------------------------------
    mix_low <- MIX_LOW_EXPANSION

    # -----------------------------------------------------------------
    # 2. Baseline maximum expansion rate per mL tumour volume (Vmax1)
    #    Reference population (Eq 10): prior-ASCT fractional change and a
    #    day-7 CD4/CD8 ratio power term; IIV on the reference baseline only.
    #    Low-expansion subpopulation (Eq 11): prior-ASCT change only.
    # -----------------------------------------------------------------
    vmax1_ref <- exp(lvmax1_base_ref + etalvmax1_base_ref) *
      (1 + e_hct_prior_vmax1 * HCT_PRIOR) *
      CD4CD8_RATIO^e_cd4cd8_vmax1
    vmax1_low <- exp(lvmax1_base_low) *
      (1 + e_hct_prior_vmax1 * HCT_PRIOR)
    vmax1 <- mix_low * vmax1_low + (1 - mix_low) * vmax1_ref

    # -----------------------------------------------------------------
    # 3. Remaining individual parameters
    # -----------------------------------------------------------------
    km1 <- exp(lkm1)
    km5 <- exp(lkm5)
    k12 <- exp(lk12)
    k23 <- exp(lk23)
    k34 <- exp(lk34)
    ke4 <- exp(lke4)
    vmax5_2i <- exp(lvmax5_2 + etalvmax5_2)

    # -----------------------------------------------------------------
    # 4. ODE system (Equations 5-9). Expansion of TN/TCM/TEM is limited by
    #    the RESPECTIVE T-cell concentration in the denominator.
    # -----------------------------------------------------------------
    d/dt(tn) <- vmax1 * tumor * tn / (km1 + tn) + kp1 * tn - k12 * tn - ke1 * tn
    d/dt(tcm) <- vmax1 * tumor * tcm / (km1 + tcm) + kp2 * tcm + k12 * tn - k23 * tcm - ke2 * tcm
    d/dt(tem) <- vmax1 * tumor * tem / (km1 + tem) + kp3 * tem + k23 * tcm - k34 * tem - ke3 * tem
    d/dt(teff) <- k34 * tem - ke4 * teff
    d/dt(tumor) <- k5 * (1 - tumor / k0) * tumor -
      vmax5_1 * tn * tumor / (km5 + tumor) -
      vmax5_2i * tcm * tumor / (km5 + tumor) -
      vmax5_3 * tem * tumor / (km5 + tumor) -
      vmax5_4 * teff * tumor / (km5 + tumor)

    # -----------------------------------------------------------------
    # 5. Initial conditions: imputed 0.1 cells/uL per phenotype; CD19+
    #    tumour starts at the measured baseline metabolic tumour volume.
    # -----------------------------------------------------------------
    tn(0) <- init_cart
    tcm(0) <- init_cart
    tem(0) <- init_cart
    teff(0) <- init_cart
    tumor(0) <- TUM_VOL_TOTAL

    # -----------------------------------------------------------------
    # 6. Observations (each species observed directly). Residual error is
    #    log-transform-both-sides (log-normal), per Equation 3.
    # -----------------------------------------------------------------
    pred_tn <- tn
    pred_tcm <- tcm
    pred_tem <- tem
    pred_teff <- teff
    pred_tumor <- tumor
    pred_tn ~ lnorm(expSd_pred_tn)
    pred_tcm ~ lnorm(expSd_pred_tcm)
    pred_tem ~ lnorm(expSd_pred_tem)
    pred_teff ~ lnorm(expSd_pred_teff)
    pred_tumor ~ lnorm(expSd_pred_tumor)
  })
}
