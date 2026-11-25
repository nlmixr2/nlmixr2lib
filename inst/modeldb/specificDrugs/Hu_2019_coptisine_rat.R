Hu_2019_coptisine_rat <- function() {
  description <- paste(
    "Preclinical (rat). Mechanism-based PK-PD model of intravenous coptisine",
    "in male Sprague-Dawley rats challenged with LPS (100 ug/kg IV at time 0).",
    "Two-compartment coptisine PK with a lung leaf compartment; plasma",
    "TNF-alpha produced by a two-stage LPS-driven zero-order input (0-0.33 h",
    "and 0.33-1 h after LPS) that coptisine inhibits linearly, with",
    "Michaelis-Menten elimination; then a lung TNF-alpha -> lung iNOS mRNA ->",
    "lung iNOS -> plasma iNOS -> NO precursor (iNOS^delta) -> plasma NO",
    "stimulation cascade. Typical values only (between-animal and residual",
    "variances are not reported). The coptisine PK and plasma TNF-alpha",
    "layers reproduce the paper; the lung TNF-alpha / iNOS / NO cascade as",
    "printed does NOT reproduce the paper's own data (lung TNF-alpha about",
    "1e3-fold low, iNOS and NO peaks at about 3.5 h and 4.5 h instead of",
    "7-8 h) -- see the validation vignette Errata before using that half."
  )
  reference <- paste(
    "Hu Y, Wang L, Xiang L, Wu J, Huang W, Xu C, Meng X, Wang P.",
    "Pharmacokinetic-Pharmacodynamic Modeling for Coptisine Challenge of",
    "Inflammation in LPS-Stimulated Rats. Sci Rep. 2019;9:1450.",
    "doi:10.1038/s41598-018-38164-4"
  )
  vignette <- "Hu_2019_coptisine_rat"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL"
  )

  # Lung TNF-alpha, the iNOS mRNA pool, plasma iNOS and the NO precursor
  # pool are states of this paper's inflammatory cascade with no
  # register entry; plasma TNF-alpha (tnf), lung iNOS (inos) and plasma
  # NO (no) use the registered inflammatory-mediator names.
  paper_specific_compartments <- c("tnf_lung", "inos_mrna", "inos_plasma", "no_pre")

  compartmentData <- list(
    central = list(analyte = "coptisine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "coptisine", units = "mg", specimen = "tissue", verified = TRUE),
    lung = list(analyte = "coptisine", units = "ng/g", specimen = "tissue", verified = TRUE),
    tnf = list(analyte = "TNF-alpha", units = "pg/mL", specimen = "plasma", verified = TRUE),
    tnf_lung = list(analyte = "TNF-alpha", units = "pg/g", specimen = "tissue", verified = TRUE),
    inos_mrna = list(analyte = "iNOS mRNA", units = NA_character_, specimen = "tissue", verified = TRUE),
    inos = list(analyte = "iNOS", units = NA_character_, specimen = "tissue", verified = TRUE),
    inos_plasma = list(analyte = "iNOS", units = NA_character_, specimen = "plasma", verified = TRUE),
    no_pre = list(analyte = "nitric oxide precursor", units = "umol/L", specimen = "not applicable", verified = TRUE),
    no = list(analyte = "nitric oxide (nitrite + nitrate)", units = "umol/L", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "rat (male Sprague-Dawley)",
    n_subjects = 152L,
    n_studies = 4L,
    age_range = "not reported",
    weight_range = "200 +/- 20 g",
    sex_female_pct = 0,
    disease_state = paste(
      "Acute systemic inflammation induced by Escherichia coli O55:B5",
      "lipopolysaccharide, 100 ug/kg IV via the tail vein at time 0."
    ),
    dose_range = paste(
      "Coptisine 1.94, 3.87 or 7.74 mg/kg IV bolus 0.5 h after LPS, or",
      "vehicle (LPS-only arm); untreated controls sampled at time 0."
    ),
    regions = "China (Chengdu University of TCM)",
    notes = paste(
      "Four studies fitted simultaneously in Monolix 2016R1: 152 rats and",
      "1118 observations of plasma and lung coptisine, plasma and lung",
      "TNF-alpha, lung iNOS expression (Western blot, normalised to the",
      "time-0 control) and plasma NO (nitrite + nitrate) (Hu 2019 Table 2",
      "and Methods). Studies 3 and 4 were terminal (urethane anaesthesia,",
      "abdominal-artery blood and lungs, n = 2-3 per group per time).",
      "No covariates were tested."
    )
  )

  ini({
    # Hu 2019 Table 1 is the only parameter table. Its 'CV%' column is the
    # relative standard error of each estimate, not between-animal
    # variability. The between-animal variances (log-normal, Methods 'Data
    # Analysis') and the additive residual variances are not reported, so the
    # model carries no etas and fixed-zero residual SDs. All values are the
    # printed estimates, unchanged.

    # Coptisine PK (Eq 1-3)
    lvc <- log(2630)
    label("Central volume of distribution V1 per animal (mL)") # Table 1 'V1' = 2.63x10^3 (RSE 10.29%); printed 'L/kg', read as mL per rat -- see vignette Errata
    lkel <- log(0.458)
    label("Elimination rate constant k10 (1/h)") # Table 1 'k10' = 4.58x10^-1 (RSE 9.99%)
    lk12 <- log(2.05)
    label("Central-to-peripheral rate constant k12 (1/h)") # Table 1 'k12' = 2.05 (RSE 13.39%)
    lk21 <- log(0.668)
    label("Peripheral-to-central rate constant k21 (1/h)") # Table 1 'K21' = 6.68x10^-1 (RSE 9.51%)
    lk13 <- log(0.0121)
    label("Plasma-to-lung coptisine transfer rate constant k13 (1/h)") # Table 1 'k13' = 1.21x10^-2 (RSE 16.92%)
    lk31 <- log(3.77)
    label("Lung coptisine loss rate constant k31 (1/h)") # Table 1 'k31' = 3.77 (RSE 17.16%)

    # Plasma TNF-alpha (Eq 4, 11, 12)
    lkin_tnf1 <- log(645.73)
    label("TNF-alpha zero-order formation rate 0-0.33 h after LPS (pg/mL/h)") # Table 1 'k0,0-0.33h' = 645.73 (RSE 29.18%)
    lkin_tnf2 <- log(5881.32)
    label("TNF-alpha zero-order formation rate 0.33-1 h after LPS (pg/mL/h)") # Table 1 'k0,0.33-1h' = 5881.32 (RSE 31.11%)
    lslope_tnf <- log(8.83e-4)
    label("Linear coptisine inhibition constant on TNF-alpha formation InhibitCon (mL/ng)") # Table 1 'InhibitCon' = 8.83x10^-4 (RSE 47.33%); Results text '8.83 x 10^-4 mL/ng'
    lvmax_tnf <- log(5990)
    label("Maximal TNF-alpha Michaelis-Menten elimination rate Vm (pg/mL/h)") # Table 1 'Vm' = 5.99x10^3 (RSE 55.85%); printed unit 'h-1', Eq 4 requires pg/mL/h
    lkm_tnf <- log(5180)
    label("TNF-alpha concentration at half-maximal elimination km (pg/mL)") # Table 1 'km' = 5.18x10^3 (RSE 61.92%)

    # Lung TNF-alpha -> iNOS -> NO cascade (Eq 5-10)
    lk_tnf_lung <- log(3.54e-7)
    label("Plasma-to-lung TNF-alpha transfer rate constant kTNFapl (1/h)") # Table 1 'kTNFapl' = 3.54x10^-7 (RSE 177.12%)
    lkout_tnf_lung <- log(0.234)
    label("Lung TNF-alpha elimination rate constant koutTNFa (1/h)") # Table 1 'koutTNFa' = 2.34x10^-1 (RSE 7.09%)
    lkin_inos_mrna <- log(3.38e-3)
    label("Lung-TNF-alpha-driven iNOS mRNA production rate constant kiNOSmRNA (1/h)") # Table 1 'kiNOSmRNA' = 3.38x10^-3 (RSE 10.85%)
    lkout_inos_mrna <- log(2.36)
    label("iNOS mRNA elimination rate constant koutiNOSmRNA (1/h)") # Table 1 'koutiNOSmRNA' = 2.36 (RSE 7.75%)
    lkin_inos <- log(0.299)
    label("iNOS production rate constant from iNOS mRNA kiNOS (1/h)") # Table 1 'kiNOS' = 2.99x10^-1 (RSE 10.33%)
    lkout_inos <- log(3.72)
    label("Lung iNOS elimination rate constant koutiNOS (1/h)") # Table 1 'koutiNOS' = 3.72 (RSE 10.06%)
    lk_inos_lp <- log(23.41)
    label("Lung-to-plasma iNOS transfer rate constant kiNOSlp (1/h)") # Table 1 'kiNOSlp' = 23.41 (RSE 19.96%)
    lkout_inosp <- log(1.92)
    label("Plasma iNOS elimination rate constant koutiNOSp (1/h)") # Table 1 'koutiNOSp' = 1.92 (RSE 6.25%)
    lkin_nopre <- log(772.17)
    label("NO-precursor production rate constant from plasma iNOS kiNOSNO (1/h)") # Table 1 'kiNOSNO' = 772.17 (RSE 21.58%)
    ldelta <- log(1.23)
    label("Amplification exponent on plasma iNOS for NO-precursor production (unitless)") # Table 1 'Delta' = 1.23 (RSE 29.71%); exponent per the Eq 9 MathML superscript
    lkin_no <- log(354.68)
    label("NO production rate constant from the NO precursor kinNO (1/h)") # Table 1 'kinNO' = 354.68 (RSE 63.97%)
    lkout_no <- log(3.46)
    label("Plasma NO elimination rate constant koutNO (1/h)") # Table 1 'koutNO' = 3.46 (RSE 11.13%)

    # Residual error: additive on every output (Methods 'Data Analysis'),
    # variances not reported -> fixed at zero.
    addSd <- fixed(0)
    label("Additive residual SD, plasma coptisine (ng/mL; not reported)") # Methods: additive residual model; value not reported
    addSd_Clung <- fixed(0)
    label("Additive residual SD, lung coptisine (ng/g; not reported)") # Methods: additive residual model; value not reported
    addSd_tnf <- fixed(0)
    label("Additive residual SD, plasma TNF-alpha (pg/mL; not reported)") # Methods: additive residual model; value not reported
    addSd_tnf_lung <- fixed(0)
    label("Additive residual SD, lung TNF-alpha (pg/g; not reported)") # Methods: additive residual model; value not reported
    addSd_inos <- fixed(0)
    label("Additive residual SD, lung iNOS ratio (not reported)") # Methods: additive residual model; value not reported
    addSd_no <- fixed(0)
    label("Additive residual SD, plasma NO (umol/L; not reported)") # Methods: additive residual model; value not reported
  })

  model({
    vc <- exp(lvc)
    kel <- exp(lkel)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    k13 <- exp(lk13)
    k31 <- exp(lk31)

    kin_tnf1 <- exp(lkin_tnf1)
    kin_tnf2 <- exp(lkin_tnf2)
    slope_tnf <- exp(lslope_tnf)
    vmax_tnf <- exp(lvmax_tnf)
    km_tnf <- exp(lkm_tnf)

    k_tnf_lung <- exp(lk_tnf_lung)
    kout_tnf_lung <- exp(lkout_tnf_lung)
    kin_inos_mrna <- exp(lkin_inos_mrna)
    kout_inos_mrna <- exp(lkout_inos_mrna)
    kin_inos <- exp(lkin_inos)
    kout_inos <- exp(lkout_inos)
    k_inos_lp <- exp(lk_inos_lp)
    kout_inosp <- exp(lkout_inosp)
    kin_nopre <- exp(lkin_nopre)
    delta <- exp(ldelta)
    kin_no <- exp(lkin_no)
    kout_no <- exp(lkout_no)

    # Coptisine: dose in mg per rat, V1 in mL per rat, Cc in ng/mL.
    Cc <- 1e6 * central / vc
    # Amount of coptisine in V1 (ng), the driver of lung uptake (Eq 3 read
    # on the amount basis that Fig 2 implies by placing cop in V1).
    amt_cop_ng <- Cc * vc

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    # Lung coptisine concentration (ng/g). Eq 3 has no return term in
    # Eq 1: the lung is a leaf compartment.
    d/dt(lung) <- k13 * amt_cop_ng - k31 * lung
    Clung <- lung

    # LPS (100 ug/kg IV) is given at time 0, so `t` is time after LPS.
    # Two-stage zero-order TNF-alpha formation, zero after 1 h (Table 1,
    # supplement Fig S1).
    kin_tnf <- 0
    if (t < 0.33) {
      kin_tnf <- kin_tnf1
    }
    if (t >= 0.33 & t < 1) {
      kin_tnf <- kin_tnf2
    }
    # Eq 11-12: linear inhibition of the formation rate by plasma coptisine.
    # Valid for Cc < 1 / slope_tnf (about 1133 ng/mL).
    inh_cop <- slope_tnf * Cc
    d/dt(tnf) <- kin_tnf * (1 - inh_cop) - vmax_tnf * tnf / (km_tnf + tnf)

    # Eq 5 read on the amount basis (TNF-alpha in V1, pg): lung TNF-alpha
    # concentration in pg/g.
    amt_tnf_pg <- tnf * vc
    d/dt(tnf_lung) <- k_tnf_lung * amt_tnf_pg - kout_tnf_lung * tnf_lung

    # Eq 6-10. The paper writes no baselines; every state starts at zero.
    d/dt(inos_mrna) <- kin_inos_mrna * tnf_lung - kout_inos_mrna * inos_mrna
    d/dt(inos) <- kin_inos * inos_mrna - kout_inos * inos
    d/dt(inos_plasma) <- k_inos_lp * inos - kout_inosp * inos_plasma
    # max() keeps the power term defined if the solver steps plasma iNOS
    # fractionally below zero near t = 0.
    d/dt(no_pre) <- kin_nopre * max(inos_plasma, 0)^delta - kin_no * no_pre
    d/dt(no) <- kin_no * no_pre - kout_no * no

    Cc ~ add(addSd)
    Clung ~ add(addSd_Clung)
    tnf ~ add(addSd_tnf)
    tnf_lung ~ add(addSd_tnf_lung)
    inos ~ add(addSd_inos)
    no ~ add(addSd_no)
  })
}
