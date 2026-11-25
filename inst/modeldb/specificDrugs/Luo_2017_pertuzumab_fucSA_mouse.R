Luo_2017_pertuzumab_fucSA_mouse <- function() {
  description <- "Preclinical (mouse, ICR, male). Two-compartment IV-bolus PK model for native (fucosylated, sialylated; Fuc+SA+) pertuzumab, a humanized anti-HER2 IgG1 monoclonal antibody, after a single 1 mg/kg intravenous bolus in ICR mice (Luo 2017). Deterministic typical-value model (no IIV, no residual error): the paper reports per-group biexponential (two-compartment) point estimates only, with between-animal spread given as mean +/- SD across n = 5 mice rather than a fitted variance model. Companion glycoform / asialofetuin-protection models: Luo_2017_pertuzumab_nofucSA_mouse, Luo_2017_pertuzumab_fucNoSA_mouse, Luo_2017_pertuzumab_nofucNoSA_mouse, Luo_2017_pertuzumab_fucNoSA_asf_mouse, Luo_2017_pertuzumab_nofucNoSA_asf_mouse."
  reference <- paste(
    "Luo C, Chen S, Xu N, Wang C, Sai WB, Zhao W, Li YC, Hu XJ, Tian H, Gao XD, Yao WB.",
    "Glycoengineering of pertuzumab and its impact on the pharmacokinetic/pharmacodynamic properties.",
    "Sci Rep. 2017;7:46347. doi:10.1038/srep46347. PMCID PMC5387714.",
    sep = " "
  )
  vignette <- "Luo_2017_pertuzumab_glycoengineering"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  # linCmt() auto-solves the two-compartment disposition; the only named state
  # is the central compartment, which holds pertuzumab measured in mouse serum
  # by a HER2-binding ELISA (Methods, 'In vitro pharmacokinetics of pertuzumab').
  compartmentData <- list(
    central = list(analyte = "pertuzumab", units = "mg", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "mouse (ICR, male)",
    n_subjects = 5L,
    age_range = "6 weeks",
    sex_female_pct = 0,
    disease_state = "Healthy (non-tumour-bearing) ICR mice.",
    dose_range = "1 mg/kg pertuzumab as a single intravenous (tail-vein) bolus.",
    glycoform = "Fuc+SA+ (native pertuzumab expressed in wild-type CHO-K1 cells: fucosylated and sialylated).",
    notes = paste(
      "Luo 2017 Table 2, row 'PertuzumabFuc+SA+'. Free pertuzumab in serum was",
      "measured by a HER2-specific ELISA at 0.5, 1, 2, 4, 8 and 24 h and 2, 3, 4,",
      "5, 7, 9, 14, 20 and 26 days after dosing, and the concentration-time data",
      "were fit to a two-compartment model C(t) = A*exp(-alpha*t) + B*exp(-beta*t)",
      "in WinNonlin 5.2. All PK parameters are expressed per kilogram body weight.",
      sep = " "
    )
  )

  ini({
    # Two-compartment IV-bolus disposition, native pertuzumab (Fuc+SA+).
    # Luo 2017 Table 2 reports (per kg): CL = 18.1 mL/day/kg, V1 = 57.3 mL/kg,
    # Vss = 372.1 mL/kg, distribution half-life alpha_HL = 0.08 day and terminal
    # half-life beta_HL = 14.7 day. Values are converted from mL to L (/1000)
    # so that, with a 1 mg (= 1 mg/kg) bolus, Cc = amount/vc is in mg/L = ug/mL.
    #   lcl : CL            = 18.1 mL/day/kg -> 0.01810 L/day/kg (Table 2)
    #   lvc : V1 (central)  = 57.3 mL/kg     -> 0.05730 L/kg     (Table 2)
    #   lvp : Vp = Vss - V1 = 314.8 mL/kg    -> 0.31480 L/kg     (Table 2, Vss - V1)
    #   lq  : Q derived from CL, V1, alpha_HL, beta_HL via the standard
    #         two-compartment relations (k10 = CL/V1; alpha+beta = k10+k12+k21;
    #         alpha*beta = k10*k21) = 407.0 mL/day/kg -> 0.40700 L/day/kg.
    #         The Q + Vp pair reproduces the reported Vss to within <1%.
    lcl <- log(0.01810); label("Clearance (L/day/kg)") # Luo 2017 Table 2, Fuc+SA+ CL
    lvc <- log(0.05730); label("Central volume V1 (L/kg)") # Luo 2017 Table 2, Fuc+SA+ V1
    lq <- log(0.40700); label("Intercompartmental clearance Q (L/day/kg)") # Luo 2017 Table 2, derived from CL/V1/alpha_HL/beta_HL
    lvp <- log(0.31480); label("Peripheral volume Vp (L/kg)") # Luo 2017 Table 2, Fuc+SA+ Vss - V1
  })

  model({
    cl <- exp(lcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)
    # Two-compartment IV-bolus disposition (dose delivered into the central
    # compartment). Deterministic typical-value prediction: no residual-error
    # endpoint because Luo 2017 reports no fitted measurement-error model.
    Cc <- linCmt()
  })
}
