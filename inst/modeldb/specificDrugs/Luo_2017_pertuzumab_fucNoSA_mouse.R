Luo_2017_pertuzumab_fucNoSA_mouse <- function() {
  description <- "Preclinical (mouse, ICR, male). Two-compartment IV-bolus PK model for fucosylated, desialylated (Fuc+SA-) pertuzumab, a humanized anti-HER2 IgG1 monoclonal antibody (CHO-K1 material enzymatically desialylated with sialidase A), after a single 1 mg/kg intravenous bolus in ICR mice (Luo 2017). Deterministic typical-value model (no IIV, no residual error): the paper reports per-group biexponential (two-compartment) point estimates only, with between-animal spread given as mean +/- SD across n = 5 mice rather than a fitted variance model. Companion glycoform / asialofetuin-protection models: Luo_2017_pertuzumab_fucSA_mouse, Luo_2017_pertuzumab_nofucSA_mouse, Luo_2017_pertuzumab_nofucNoSA_mouse, Luo_2017_pertuzumab_fucNoSA_asf_mouse, Luo_2017_pertuzumab_nofucNoSA_asf_mouse."
  reference <- paste(
    "Luo C, Chen S, Xu N, Wang C, Sai WB, Zhao W, Li YC, Hu XJ, Tian H, Gao XD, Yao WB.",
    "Glycoengineering of pertuzumab and its impact on the pharmacokinetic/pharmacodynamic properties.",
    "Sci Rep. 2017;7:46347. doi:10.1038/srep46347. PMCID PMC5387714.",
    sep = " "
  )
  vignette <- "Luo_2017_pertuzumab_glycoengineering"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

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
    glycoform = "Fuc+SA- (fucosylated CHO-K1 pertuzumab desialylated in vitro with Arthrobacter ureafaciens sialidase A).",
    notes = paste(
      "Luo 2017 Table 2, row 'PertuzumabFuc+SA-'. Free pertuzumab in serum was",
      "measured by a HER2-specific ELISA and the concentration-time data were fit",
      "to a two-compartment model C(t) = A*exp(-alpha*t) + B*exp(-beta*t) in",
      "WinNonlin 5.2. All PK parameters are expressed per kilogram body weight.",
      sep = " "
    )
  )

  ini({
    # Two-compartment IV-bolus disposition, fucosylated desialylated pertuzumab
    # (Fuc+SA-). Luo 2017 Table 2 (per kg): CL = 19.9 mL/day/kg, V1 = 55.1 mL/kg,
    # Vss = 367.4 mL/kg, alpha_HL = 0.05 day, beta_HL = 13.1 day. Converted mL->L
    # (/1000); with a 1 mg (= 1 mg/kg) bolus, Cc = amount/vc is in mg/L = ug/mL.
    # lq: Q derived from CL, V1, alpha_HL, beta_HL via the standard two-compartment
    # relations; the Q + Vp pair reproduces the reported Vss to within <1%.
    lcl <- log(0.01990); label("Clearance (L/day/kg)") # Luo 2017 Table 2, Fuc+SA- CL
    lvc <- log(0.05510); label("Central volume V1 (L/kg)") # Luo 2017 Table 2, Fuc+SA- V1
    lq <- log(0.63500); label("Intercompartmental clearance Q (L/day/kg)") # Luo 2017 Table 2, derived from CL/V1/alpha_HL/beta_HL
    lvp <- log(0.31230); label("Peripheral volume Vp (L/kg)") # Luo 2017 Table 2, Fuc+SA- Vss - V1
  })

  model({
    cl <- exp(lcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)
    Cc <- linCmt()
  })
}
