Krippendorff_2012_zalutumumab_cyno <- function() {
  description <- paste(
    "Preclinical (cynomolgus monkey). QSP cell-level PK/PD model of the",
    "anti-EGFR IgG1 antibody zalutumumab: the plasma-interstitial PK model of",
    "Lammerts van Bueren et al. (2006) with its Michaelis-Menten target term",
    "replaced by a cell-level EGF receptor trafficking model (free surface",
    "receptor, internalised free receptor, EGF-receptor complex and",
    "drug-receptor complex per cell) scaled to the whole body by the number of",
    "target-expressing cells. Normal cells plus a tumour-cell population (1% of",
    "the normal cell number) with 10-fold elevated EGFR from either increased",
    "receptor synthesis or decreased internalisation. Deterministic: every",
    "parameter is taken from in vitro data, the inherited PK model or an",
    "explicit calculation, so no etas and no error model are encoded.",
    "Outputs plasma concentration, residual EGFR activation and receptor",
    "saturation for normal and tumour cells.",
    sep = " "
  )
  reference <- paste(
    "Krippendorff BF, Oyarzun DA, Huisinga W. Predicting the F(ab)-mediated",
    "effect of monoclonal antibodies in vivo by combining cell-level kinetic",
    "and pharmacokinetic modelling. J Pharmacokinet Pharmacodyn.",
    "2012;39(2):125-139. doi:10.1007/s10928-012-9243-7. PK backbone",
    "(Table 1) inherited from Lammerts van Bueren JJ et al. Cancer Res.",
    "2006;66:7630-7638. doi:10.1158/0008-5472.CAN-05-4010.",
    sep = " "
  )
  vignette <- "Krippendorff_2012_zalutumumab_cyno"

  units <- list(time = "h", dosing = "mg/kg", concentration = "ug/mL")

  covariateData <- list()

  # Per-cell receptor-trafficking states of the cell-level model (Eqs. 33-40).
  # Each is a number of molecules per cell, not a drug amount, so none fits
  # the amount-based `target` / `complex` canonicals.
  paper_specific_compartments <- c(
    "target_normal",
    "target_int_normal",
    "complex_egf_normal",
    "complex_normal",
    "target_tumor",
    "target_int_tumor",
    "complex_egf_tumor",
    "complex_tumor"
  )

  compartmentData <- list(
    central = list(analyte = "zalutumumab", units = "mg/kg", specimen = "plasma", verified = TRUE),
    isf = list(analyte = "zalutumumab", units = "mg/kg", specimen = "tissue", verified = TRUE),
    target_normal = list(
      analyte = "free surface EGFR, normal cell",
      units = "molecules/cell",
      specimen = "tissue",
      verified = TRUE
    ),
    target_int_normal = list(
      analyte = "internalised free EGFR, normal cell",
      units = "molecules/cell",
      specimen = "tissue",
      verified = TRUE
    ),
    complex_egf_normal = list(
      analyte = "EGF-EGFR complex (activated receptor), normal cell",
      units = "molecules/cell",
      specimen = "tissue",
      verified = TRUE
    ),
    complex_normal = list(
      analyte = "zalutumumab-EGFR complex, normal cell",
      units = "molecules/cell",
      specimen = "tissue",
      verified = TRUE
    ),
    target_tumor = list(
      analyte = "free surface EGFR, tumour cell",
      units = "molecules/cell",
      specimen = "tumor",
      verified = TRUE
    ),
    target_int_tumor = list(
      analyte = "internalised free EGFR, tumour cell",
      units = "molecules/cell",
      specimen = "tumor",
      verified = TRUE
    ),
    complex_egf_tumor = list(
      analyte = "EGF-EGFR complex (activated receptor), tumour cell",
      units = "molecules/cell",
      specimen = "tumor",
      verified = TRUE
    ),
    complex_tumor = list(
      analyte = "zalutumumab-EGFR complex, tumour cell",
      units = "molecules/cell",
      specimen = "tumor",
      verified = TRUE
    )
  )

  population <- list(
    species = "cynomolgus monkey (Macaca fascicularis)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    disease_state = "healthy (no tumour); tumour cells are a hypothetical in silico population",
    dose_range = "Single IV bolus 2, 20 and 40 mg/kg zalutumumab",
    notes = paste(
      "No parameter was fitted. The plasma-interstitial PK structure and",
      "values (Table 1) come from the cynomolgus monkey compartment model of",
      "Lammerts van Bueren et al. (2006); the receptor-trafficking rate",
      "constants (Table 2) were measured in human fibroblasts (Starbuck and",
      "Lauffenburger 1992; Knauer et al. 1984), with recycling x4 and free",
      "receptor internalisation x1/4 for monkeys; zalutumumab affinity is",
      "from Peipp et al. (2008, Table 4). The prediction was checked against",
      "plasma concentrations from the monkey study (Fig. 3a, data courtesy of",
      "Genmab, number of animals not reported) and its receptor-saturation",
      "findings (Table 3). All amounts and volumes are per kg body weight.",
      sep = " "
    )
  )

  ini({
    # ---- PK backbone (Table 1; Lammerts van Bueren 2006 monkey model) ----
    lvc <- fixed(log(70)); label("Plasma volume Vpla (mL/kg)") # Table 1: V pla = 70 ml/kg
    lvisf <- fixed(log(35)); label("Interstitial volume Vint (mL/kg)") # Table 1: V int = 35 ml/kg
    lk12 <- fixed(log(0.043)); label("Plasma-to-interstitial transfer rate constant kpi (1/h)") # Table 1: k pi = 0.043 1/h
    lk21 <- fixed(log(0.043)); label("Interstitial-to-plasma transfer rate constant kip (1/h)") # Table 1: k ip = 0.043 1/h
    lkel <- fixed(log(0.0055)); label("Rate constant kel giving the target-independent clearance CLlin = Vpla * kel (1/h)") # Table 1: kel = 0.0055 1/h; CLlin = V pla * k el
    bmax_pk <- fixed(2); label("Whole-body binding capacity Bmax,PK of the Michaelis-Menten model (mg/kg)") # Table 1: B max,PK = 2 (unit printed as mg/h/kg; see vignette)

    # ---- EGF receptor system (Table 2, rate constants in 1/min as printed) ----
    kon_egf <- fixed(7.2e-2); label("EGF-receptor association rate constant konL (1/(nM*min))") # Table 2: konL = 7.2e-2 1/(nM min)
    koff_egf <- fixed(0.34); label("EGF-receptor dissociation rate constant koffL (1/min)") # Table 2: koffL = 0.34 1/min
    kint_target <- fixed(0.03); label("Free surface receptor internalisation rate constant kdegR, human (1/min)") # Table 2: kdegR = 0.03 1/min
    kint_egf <- fixed(0.03); label("EGF-receptor complex internalisation rate constant kdegRL (1/min)") # Table 2: kdegRL = 0.03 1/min
    ksyn_target <- fixed(130); label("Receptor expression rate in normal cells ksynR,N (receptors/min per cell)") # Table 2: kR,N = 130 receptors/min per cell
    krecy_target <- fixed(5.8e-2); label("Internalised free receptor recycling rate constant krecyRi, human (1/min)") # Table 2: krecyRi = 5.8e-2 1/min
    kdeg_target_int <- fixed(2.2e-3); label("Internalised free receptor degradation rate constant kdegRi (1/min)") # Table 2: kdegRi = 2.2e-3 1/min
    kint <- fixed(0.005); label("Drug-receptor complex internalisation rate constant kdegRC (1/h)") # Table 2: kdegRC = 0.005 1/h
    mw <- fixed(148000); label("Antibody molecular weight MWmAbs (g/mol)") # Table 2: MWmAbs = 148000 Dalton (g/mol)
    f_recy_cyno <- fixed(4); label("Monkey multiplier on krecyRi (unitless)") # Table 2 footnote: krecyRi multiplied by 4 for monkeys
    f_kint_cyno <- fixed(0.25); label("Monkey multiplier on kdegR (unitless)") # Table 2 footnote: kdegR multiplied by 1/4 for monkeys
    egf <- fixed(2.36e-3); label("Interstitial EGF concentration L (nM)") # Theoretical section, after Eq. 6: L = 2.36e-3 nM
    navo <- fixed(6.02e23); label("Avogadro constant Navo (1/mol)") # After Eq. 11: Navo = 6.02e23 1/mol

    # ---- Drug-receptor binding (zalutumumab) ----
    kd_invitro <- fixed(7); label("In vitro zalutumumab-EGFR dissociation constant KD (nM)") # Table 4: zalutumumab 7e-9 M
    sfiviv <- fixed(2.1); label("In vitro-in vivo scaling factor SFiviv as reported, equal to KM,cell/KM,PK (unitless)") # Results, Evaluation against in vivo data: SFiviv = 2.1; equals KM,cell/KM,PK from Tables 1-2 (Eq. 18 orientation is its reciprocal)
    f_kint_drug <- fixed(1); label("Multiplier on kdegRC for a hypothetical antibody; Ncell keeps its zalutumumab value (unitless)") # Results, Downregulation: 5- and 10-fold kdegRC in Fig. 4b-d; 1 = zalutumumab

    # ---- Tumour-cell population (Eqs. 30-40) ----
    frac_tumor <- fixed(0.01); label("Number of tumour cells as a fraction of normal cells NT/NN (unitless)") # Theoretical section, normal and tumour cells: tumour cells set to 1% of normal cells
    fold_ksyn_tumor <- fixed(10); label("Tumour/normal ratio of receptor synthesis ksynR,T/ksynR,N (unitless)") # Fig. 5 caption: 10-fold increased receptor expression
    fold_kint_tumor <- fixed(1); label("Tumour/normal ratio of internalisation kdegR,T/kdegR,N = kdegRL,T/kdegRL,N (unitless)") # Fig. 5 caption: 0.1 for the decreased-internalisation tumour; 1 = unaltered
  })

  model({
    # Unit bookkeeping: time in h, receptor rate constants converted from 1/min.
    sf_unit <- 1e9 / navo # nmol per molecule (SFunit = 1e9/Navo)
    vc <- exp(lvc)
    visf <- exp(lvisf)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    kel <- exp(lkel)
    q_pi <- vc * k12 # q_pi = Vpla * kpi (Table 1)
    q_ip <- visf * k21 # q_ip = Vint * kip (Table 1)
    cl_lin <- vc * kel # CLlin = Vpla * kel (Table 1)

    kon_egf_h <- kon_egf * 60
    koff_egf_h <- koff_egf * 60
    kdeg_int_h <- kdeg_target_int * 60
    krecy_h <- krecy_target * 60 * f_recy_cyno
    kint_r_n <- kint_target * 60 * f_kint_cyno
    kint_rl_n <- kint_egf * 60
    ksyn_n <- ksyn_target * 60
    kint_r_t <- kint_r_n * fold_kint_tumor
    kint_rl_t <- kint_rl_n * fold_kint_tumor
    ksyn_t <- ksyn_n * fold_ksyn_tumor

    # Drug binding: konC = konL (Table 2); affinity changes act through koffC
    # with konC unchanged (Results, Affinity and dose). The KD scaling follows
    # Eq. 18 (KM,PK = SFiviv' KM,cell, so SFiviv' = 1/2.1): it brings the cell
    # model's KM down to the KM,PK of the inherited PK model, which the
    # reported 2.1 = KM,cell/KM,PK would move the other way.
    kon <- kon_egf_h
    # Ncell (Eq. 17) is a property of the receptor system, derived once from
    # zalutumumab's kdegRC; drug-property scans change kdegRC only here.
    kint_drug <- kint * f_kint_drug
    kd <- kd_invitro / sfiviv
    koff <- kd * kon

    # Number of target-expressing cells per kg from Eq. 17,
    # Ncell = kdegRC * Bmax,PK / (ksynR * SFunit), Bmax,PK converted to nmol/kg.
    bmax_pk_nmol <- bmax_pk / mw * 1e6
    n_normal <- kint * bmax_pk_nmol / (ksyn_n * sf_unit)
    n_tumor <- frac_tumor * n_normal

    # Pre-treatment steady state of each cell type (Methods)
    egf_term_n <- egf * kon_egf_h * kint_rl_n / (koff_egf_h + kint_rl_n)
    egf_term_t <- egf * kon_egf_h * kint_rl_t / (koff_egf_h + kint_rl_t)
    r_ss_n <- ksyn_n / (kint_r_n * kdeg_int_h / (krecy_h + kdeg_int_h) + egf_term_n)
    r_ss_t <- ksyn_t / (kint_r_t * kdeg_int_h / (krecy_h + kdeg_int_h) + egf_term_t)
    rl_ss_n <- kon_egf_h * egf * r_ss_n / (koff_egf_h + kint_rl_n)
    rl_ss_t <- kon_egf_h * egf * r_ss_t / (koff_egf_h + kint_rl_t)

    target_normal(0) <- r_ss_n
    target_int_normal(0) <- kint_r_n * r_ss_n / (krecy_h + kdeg_int_h)
    complex_egf_normal(0) <- rl_ss_n
    target_tumor(0) <- r_ss_t
    target_int_tumor(0) <- kint_r_t * r_ss_t / (krecy_h + kdeg_int_h)
    complex_egf_tumor(0) <- rl_ss_t

    # Interstitial free drug concentration in nM (states are mg/kg, volumes mL/kg)
    c_isf_nm <- isf / visf / mw * 1e9

    # Drug-receptor binding fluxes per cell (molecules/h)
    bind_n <- kon * target_normal * c_isf_nm - koff * complex_normal
    bind_t <- kon * target_tumor * c_isf_nm - koff * complex_tumor

    # Eqs. 30-32, written for amounts (mg/kg); the cell terms convert
    # molecules -> nmol (SFunit) -> mg (x MW / 1e6).
    d/dt(central) <- -q_pi * central / vc + q_ip * isf / visf - cl_lin * central / vc
    d/dt(isf) <- q_pi * central / vc - q_ip * isf / visf -
      (n_normal * bind_n + n_tumor * bind_t) * sf_unit * mw / 1e6

    # Eqs. 33-36 (normal cells)
    d/dt(target_normal) <- ksyn_n - kon_egf_h * target_normal * egf - kon * target_normal * c_isf_nm +
      koff_egf_h * complex_egf_normal + koff * complex_normal - kint_r_n * target_normal +
      krecy_h * target_int_normal
    d/dt(target_int_normal) <- kint_r_n * target_normal - krecy_h * target_int_normal -
      kdeg_int_h * target_int_normal
    d/dt(complex_egf_normal) <- kon_egf_h * egf * target_normal - koff_egf_h * complex_egf_normal -
      kint_rl_n * complex_egf_normal
    d/dt(complex_normal) <- bind_n - kint_drug * complex_normal

    # Eqs. 37-40 (tumour cells)
    d/dt(target_tumor) <- ksyn_t - kon_egf_h * target_tumor * egf - kon * target_tumor * c_isf_nm +
      koff_egf_h * complex_egf_tumor + koff * complex_tumor - kint_r_t * target_tumor +
      krecy_h * target_int_tumor
    d/dt(target_int_tumor) <- kint_r_t * target_tumor - krecy_h * target_int_tumor -
      kdeg_int_h * target_int_tumor
    d/dt(complex_egf_tumor) <- kon_egf_h * egf * target_tumor - koff_egf_h * complex_egf_tumor -
      kint_rl_t * complex_egf_tumor
    d/dt(complex_tumor) <- bind_t - kint_drug * complex_tumor

    # Residual receptor activity (Eq. 26, percent) and receptor saturation (Eq. 25)
    rl_pct_normal <- 100 * complex_egf_normal / rl_ss_n
    rl_pct_tumor <- 100 * complex_egf_tumor / rl_ss_t
    sat_normal <- complex_normal / (target_normal + complex_egf_normal + complex_normal)
    sat_tumor <- complex_tumor / (target_tumor + complex_egf_tumor + complex_tumor)

    # Plasma zalutumumab concentration (mg/mL -> ug/mL)
    Cc <- 1000 * central / vc
  })
}
