Wen_2026_salbutamol_pbpk <- function() {
  description <- paste0(
    "PBPK (minimal pulmonary, population). Population minimal physiologically ",
    "based pharmacokinetic model for orally inhaled salbutamol in healthy ",
    "adult volunteers, fitted jointly to plasma and to three lung sampling ",
    "techniques (bronchosorption, bronchial biopsy and bronchoalveolar ",
    "lavage). The inhaled dose splits into a lung-deposited fraction (logit ",
    "inter-individual variability) apportioned across bronchial, bronchiolar ",
    "and alveolar deposition depots, and a swallowed fraction absorbed ",
    "first-order from the gut with fixed oral bioavailability. Each lung ",
    "region has an epithelial lining fluid (ELF), epithelium and interstitium ",
    "layer linked by passive diffusion governed by one effective permeability ",
    "and the unbound lung:plasma partition coefficient; each interstitium ",
    "exchanges with plasma at its regional blood flow. Systemic disposition ",
    "is one-compartment with log-normal inter-individual variability on ",
    "clearance. Outputs: plasma (Cc), bronchial ELF (Cbr, bronchosorption), ",
    "mean bronchial tissue (Cbr_tissue, biopsy) and mean ELF over all three ",
    "regions (Celf, BAL)."
  )
  reference <- paste(
    "Wen H, Sadiq MW, Friden M, Hohlfeld JM, Friberg LE, Svensson EM.",
    "Informing Sampling Design for Lung Distribution Studies Using a",
    "Pulmonary Population Minimal PBPK Model.",
    "Clin Pharmacokinet. 2026;65(8):1257-1272.",
    "doi:10.1007/s40262-026-01664-y",
    sep = " "
  )
  vignette <- "Wen_2026_salbutamol_pbpk"
  units <- list(time = "h", dosing = "ug", concentration = "ug/L")

  # The authors' NONMEM data set doses every inhalation into all four
  # deposition depots (the supplement code splits it with F1-F4), so a 200 ug
  # inhalation is four dose records of 200 ug, one per depot.
  dosing <- c(
    "depot_lung_bronchial",
    "depot_lung_bronchiolar",
    "depot_lung_alveolar",
    "depot_oral"
  )

  # The regional lung depots and the three tissue layers of each lung region
  # are paper-mechanistic states with no canonical role; the same elf_ / epi_ /
  # sub_ layer naming is used by the authors' rat model Wen_2025_salbutamol_rat_pbpk.
  paper_specific_compartment_pattern <- c("^depot_lung_", "^elf_", "^epi_", "^sub_")

  compartmentData <- list(
    depot_lung_bronchial = list(
      analyte = "salbutamol",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ), # ESM code ABSORB1
    depot_lung_bronchiolar = list(
      analyte = "salbutamol",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ), # ESM code ABSORB2
    depot_lung_alveolar = list(analyte = "salbutamol", units = "ug", specimen = "administration site", verified = TRUE), # ESM code ABSORB3
    depot_oral = list(analyte = "salbutamol", units = "ug", specimen = "administration site", verified = TRUE), # ESM code ABSORB4 (swallowed fraction)
    elf_bronchial = list(analyte = "salbutamol", units = "ug", specimen = "epithelial lining fluid", verified = TRUE), # ESM code BB_ELF
    epi_bronchial = list(analyte = "salbutamol", units = "ug", specimen = "tissue", verified = TRUE), # ESM code BB_EPI
    sub_bronchial = list(analyte = "salbutamol", units = "ug", specimen = "tissue", verified = TRUE), # ESM code BB_sub (interstitium)
    elf_bronchiolar = list(analyte = "salbutamol", units = "ug", specimen = "epithelial lining fluid", verified = TRUE), # ESM code Bar_ELF
    epi_bronchiolar = list(analyte = "salbutamol", units = "ug", specimen = "tissue", verified = TRUE), # ESM code Bar_EPI
    sub_bronchiolar = list(analyte = "salbutamol", units = "ug", specimen = "tissue", verified = TRUE), # ESM code Bar_sub
    elf_alveolar = list(analyte = "salbutamol", units = "ug", specimen = "epithelial lining fluid", verified = TRUE), # ESM code AV_ELF
    epi_alveolar = list(analyte = "salbutamol", units = "ug", specimen = "tissue", verified = TRUE), # ESM code AV_EPI
    sub_alveolar = list(analyte = "salbutamol", units = "ug", specimen = "tissue", verified = TRUE), # ESM code AV_sub
    central = list(analyte = "salbutamol", units = "ug", specimen = "plasma", verified = TRUE) # ESM code CENTRAL
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = 13L,
    n_studies = 1L,
    age_range = "18-50 years (inclusion criterion)",
    disease_state = "Healthy, non-smoking volunteers with normal lung function",
    dose_range = "Single 200 ug dose of commercially available salbutamol inhaled via a spacer device",
    regions = "Germany (Hannover Medical School / Fraunhofer ITEM)",
    notes = paste0(
      "Inhalation cohort of ClinicalTrials.gov NCT03524066 (Sadiq et al., the ",
      "paper's reference 4). Each participant had two bronchoscopy sessions ",
      "(1 and 6 h, 2 and 9 h, 3 and 12 h, or 4 and 18 h post-dose), each giving ",
      "two bronchosorption samples, one BAL and two bronchial biopsies; plasma ",
      "was sampled at 0.25, 0.5, 1, 2, 4, 8, 12 and 24 h. 101 plasma and 100 ",
      "lung samples (32 bronchosorption, 32 biopsy, 36 BAL). Baseline ",
      "demographics (weight, sex split) are not restated in Wen 2026 ",
      "(Methods, Study Design and Data Description)."
    )
  )

  ini({
    # Systemic PK (Table 3, 'Systemic PK parameters'). Estimated from plasma
    # data alone in the first modelling step, then fixed for the lung step.
    lcl <- log(32.8); label("Plasma clearance CL (L/h)") # Table 3: Clearance 32.8 L/h (RSE 6.87%)
    lvc <- log(211); label("Volume of distribution V (L)") # Table 3: Volume of distribution 211 L (RSE 4.55%)
    lka_oral <- log(2.01); label("First-order absorption rate constant of the swallowed fraction ka (1/h)") # Table 3: ka 2.01 /h (RSE 31.9%)
    lkdiss <- fixed(log(1000)); label("Dissolution rate constant from lung deposition depots into ELF kdiss (1/h)") # Table 3: kdiss 1000 /h FIX; Methods: 'fixed as 1000/h, assuming instantaneous dissolution'
    logitfdepot_lung <- fixed(0); label("Logit of the total lung deposition fraction F1 (unitless)") # Table 3: F1 (lung) 0.500 FIX; logit(0.5) = 0
    lfdepot_oral <- fixed(log(0.418)); label("Oral bioavailability of the swallowed fraction F2 (fraction)") # Table 3: F2 (oral) 0.418 FIX

    # Regional deposition, as fractions of the lung-deposited dose (typical
    # path model, Results); the alveolar fraction is the remainder.
    fdep_bronchial <- fixed(0.0992); label("Fraction of the lung dose deposited in the bronchial region (fraction)") # ESM NONMEM code: FBB = 0.0992 (Results text: 9.90%)
    fdep_bronchiolar <- fixed(0.226); label("Fraction of the lung dose deposited in the bronchiolar region (fraction)") # ESM NONMEM code: FBR = 0.226 (Results text: 22.6%)

    # Compound-specific parameters (Table 3, 'Compound-specific parameters').
    lbp <- fixed(log(0.96)); label("Blood-to-plasma concentration ratio R (unitless)") # Table 3: B:P ratio 0.96 FIX
    lfuelf <- fixed(log(0.74)); label("Fraction unbound in epithelial lining fluid (fraction)") # Table 3: fu_ELF 0.74 FIX
    lfu <- fixed(log(0.65)); label("Fraction unbound in plasma (fraction)") # Table 3: fu 0.65 FIX

    # Pulmonary PK parameters (Table 3, 'Pulmonary PK parameters').
    lkpulung <- log(11.0); label("Unbound lung-to-plasma partition coefficient Kp,u,lung (unitless)") # Table 3: Kp,u,lung 11.0 (RSE 8.58%)
    lpeff <- log(0.543); label("Effective permeability Peff (L/h per m2 of surface; reported as m/h)") # Table 3: Peff 0.543 m/h (RSE 11.6%)

    # IIV. Table 3 reports 'CV%' for both random effects; a CV is undefined for
    # the logit-scale F1 eta, so the column is 100 * sqrt(omega^2), and the
    # same reading is applied to CL.
    etalcl ~ 0.193^2 # Table 3: 'CV% Clearance' 19.3 (RSE 37.0%) [shrinkage 4.15%]
    etalogitfdepot_lung ~ 1.285^2 # Table 3: 'CV% LogitF1' 128.5 (RSE 15.9%) [shrinkage 7.04%]

    # Residual error: additive on log-transformed concentrations (Methods,
    # Data Description and Data Handling), one per sampling technique.
    expSd <- 0.188; label("Residual error SD on the log scale, plasma (unitless)") # Table 3: Prop error - plasma 18.8% (RSE 8.18%)
    expSd_Cbr <- 1.09; label("Residual error SD on the log scale, bronchosorption (unitless)") # Table 3: Prop error - bronchosorption 109% (RSE 22.3%)
    expSd_Cbr_tissue <- 0.848; label("Residual error SD on the log scale, biopsy (unitless)") # Table 3: Prop error - biopsy 84.8% (RSE 20.0%)
    expSd_Celf <- 1.46; label("Residual error SD on the log scale, BAL-derived ELF (unitless)") # Table 3: Prop error - BAL 146% (RSE 20.1%)
  })

  model({
    # Individual parameters
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    ka_oral <- exp(lka_oral)
    kdiss <- exp(lkdiss)
    fdepot_lung <- expit(logitfdepot_lung + etalogitfdepot_lung)
    fdepot_oral <- exp(lfdepot_oral)
    fdep_alveolar <- 1 - fdep_bronchial - fdep_bronchiolar

    bp <- exp(lbp)
    fuelf <- exp(lfuelf)
    fu <- exp(lfu)
    kpulung <- exp(lkpulung)
    peff <- exp(lpeff)
    kplung <- kpulung * fu # Eq 5: Kp,lung = Kp,u,lung x fu,p

    # Adult human lung physiology (ESM NONMEM code; Table S1). Volumes in L
    # (mL x 0.001), surface areas in m2 (cm2 x 0.0001), blood flows in L/h.
    # In the code a flux is peff * area * concentration with no further unit
    # factor, so peff carries L/h per m2.
    v_elf_bronchial <- 0.33 * 0.001
    v_epi_bronchial <- 1.74 * 0.001
    v_sub_bronchial <- 207.04 * 0.001
    v_elf_bronchiolar <- 1.65 * 0.001
    v_epi_bronchiolar <- 4.25 * 0.001
    v_sub_bronchiolar <- 66.25 * 0.001
    v_elf_alveolar <- 0.69 * 0.001
    v_epi_alveolar <- 3.58 * 0.001
    v_sub_alveolar <- 244.45 * 0.001
    a_elf_bronchial <- 305.51 * 0.0001
    a_epi_bronchial <- 317.61 * 0.0001 # ESM code SBB_EPI; Table S1 prints the epithelium volume (1.74) in this cell
    a_elf_bronchiolar <- 2772.37 * 0.0001
    a_epi_bronchiolar <- 2896.77 * 0.0001 # ESM code SBar_EPI; Table S1 prints 4.25
    a_elf_alveolar <- 700000 * 0.0001
    a_epi_alveolar <- 700000 * 0.0001 # ESM code SAV_EPI; Table S1 prints 3.58
    q_bronchial <- 1.05
    q_bronchiolar <- 8.07
    q_alveolar <- 5.2 * 60 # cardiac output, 5.2 L/min

    Cc <- central / vc

    # Passive diffusion (Eqs 2-4): ELF -> epithelium, epithelium ->
    # interstitium, interstitium -> plasma.
    j_elf_bronchial <- peff * a_elf_bronchial * (fuelf * elf_bronchial / v_elf_bronchial - epi_bronchial / v_epi_bronchial / kpulung)
    j_epi_bronchial <- peff * a_epi_bronchial * (epi_bronchial / v_epi_bronchial / kpulung - sub_bronchial / v_sub_bronchial / kpulung)
    j_bld_bronchial <- q_bronchial * (bp * sub_bronchial / v_sub_bronchial / kplung - Cc)

    j_elf_bronchiolar <- peff * a_elf_bronchiolar * (fuelf * elf_bronchiolar / v_elf_bronchiolar - epi_bronchiolar / v_epi_bronchiolar / kpulung)
    j_epi_bronchiolar <- peff * a_epi_bronchiolar * (epi_bronchiolar / v_epi_bronchiolar / kpulung - sub_bronchiolar / v_sub_bronchiolar / kpulung)
    j_bld_bronchiolar <- q_bronchiolar * (bp * sub_bronchiolar / v_sub_bronchiolar / kplung - Cc)

    j_elf_alveolar <- peff * a_elf_alveolar * (fuelf * elf_alveolar / v_elf_alveolar - epi_alveolar / v_epi_alveolar / kpulung)
    j_epi_alveolar <- peff * a_epi_alveolar * (epi_alveolar / v_epi_alveolar / kpulung - sub_alveolar / v_sub_alveolar / kpulung)
    j_bld_alveolar <- q_alveolar * (bp * sub_alveolar / v_sub_alveolar / kplung - Cc)

    d/dt(depot_lung_bronchial) <- -kdiss * depot_lung_bronchial
    d/dt(depot_lung_bronchiolar) <- -kdiss * depot_lung_bronchiolar
    d/dt(depot_lung_alveolar) <- -kdiss * depot_lung_alveolar
    d/dt(depot_oral) <- -ka_oral * depot_oral

    d/dt(elf_bronchial) <- kdiss * depot_lung_bronchial - j_elf_bronchial
    d/dt(epi_bronchial) <- j_elf_bronchial - j_epi_bronchial
    d/dt(sub_bronchial) <- j_epi_bronchial - j_bld_bronchial
    d/dt(elf_bronchiolar) <- kdiss * depot_lung_bronchiolar - j_elf_bronchiolar
    d/dt(epi_bronchiolar) <- j_elf_bronchiolar - j_epi_bronchiolar
    d/dt(sub_bronchiolar) <- j_epi_bronchiolar - j_bld_bronchiolar
    d/dt(elf_alveolar) <- kdiss * depot_lung_alveolar - j_elf_alveolar
    d/dt(epi_alveolar) <- j_elf_alveolar - j_epi_alveolar
    d/dt(sub_alveolar) <- j_epi_alveolar - j_bld_alveolar

    d/dt(central) <- ka_oral * depot_oral + j_bld_bronchial + j_bld_bronchiolar + j_bld_alveolar - cl / vc * central

    # Dose split (ESM code F1-F4): each depot receives the full inhaled dose
    # and keeps its own fraction.
    f(depot_lung_bronchial) <- fdep_bronchial * fdepot_lung
    f(depot_lung_bronchiolar) <- fdep_bronchiolar * fdepot_lung
    f(depot_lung_alveolar) <- fdep_alveolar * fdepot_lung
    f(depot_oral) <- fdepot_oral * (1 - fdepot_lung)

    # Sampled matrices (Methods, Population-Based PBPK Modelling; ESM $ERROR):
    # bronchosorption = bronchial ELF; biopsy = mean over the bronchial ELF,
    # epithelium and interstitium; BAL = mean ELF over all three regions.
    Cbr <- elf_bronchial / v_elf_bronchial
    Cbr_tissue <- (elf_bronchial + epi_bronchial + sub_bronchial) /
      (v_elf_bronchial + v_epi_bronchial + v_sub_bronchial)
    Celf <- (elf_bronchial + elf_bronchiolar + elf_alveolar) /
      (v_elf_bronchial + v_elf_bronchiolar + v_elf_alveolar)

    Cc ~ lnorm(expSd)
    Cbr ~ lnorm(expSd_Cbr)
    Cbr_tissue ~ lnorm(expSd_Cbr_tissue)
    Celf ~ lnorm(expSd_Celf)
  })
}
