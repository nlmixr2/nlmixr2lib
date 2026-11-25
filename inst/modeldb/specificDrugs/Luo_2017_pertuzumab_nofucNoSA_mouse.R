Luo_2017_pertuzumab_nofucNoSA_mouse <- function() {
  description <- "Preclinical (mouse, ICR, male). Two-compartment IV-bolus PK model for nonfucosylated, desialylated (Fuc-SA-) pertuzumab, a humanized anti-HER2 IgG1 monoclonal antibody (CHOFUT8-/- material enzymatically desialylated with sialidase A; the fastest-cleared glycoform), after a single 1 mg/kg intravenous bolus in ICR mice (Luo 2017). Deterministic typical-value model (no IIV, no residual error): the paper reports per-group biexponential (two-compartment) point estimates only, with between-animal spread given as mean +/- SD across n = 5 mice rather than a fitted variance model. Companion glycoform / asialofetuin-protection models: Luo_2017_pertuzumab_fucSA_mouse, Luo_2017_pertuzumab_nofucSA_mouse, Luo_2017_pertuzumab_fucNoSA_mouse, Luo_2017_pertuzumab_fucNoSA_asf_mouse, Luo_2017_pertuzumab_nofucNoSA_asf_mouse."
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
    glycoform = "Fuc-SA- (nonfucosylated CHOFUT8-/- pertuzumab desialylated in vitro with Arthrobacter ureafaciens sialidase A).",
    notes = paste(
      "Luo 2017 Table 2, row 'PertuzumabFuc-SA-'. Free pertuzumab in serum was",
      "measured by a HER2-specific ELISA and the concentration-time data were fit",
      "to a two-compartment model C(t) = A*exp(-alpha*t) + B*exp(-beta*t) in",
      "WinNonlin 5.2. This desialylated glycoform showed the highest clearance and",
      "shortest half-life. All PK parameters are expressed per kilogram body weight.",
      sep = " "
    )
  )

  ini({
    # Two-compartment IV-bolus disposition, nonfucosylated desialylated
    # pertuzumab (Fuc-SA-). Luo 2017 Table 2 (per kg): CL = 20.7 mL/day/kg,
    # V1 = 20.8 mL/kg, Vss = 313.5 mL/kg, alpha_HL = 0.01 day, beta_HL = 10.6 day.
    # Converted mL->L (/1000); with a 1 mg (= 1 mg/kg) bolus, Cc = amount/vc is in
    # mg/L = ug/mL. The reported alpha_HL is given to a single significant figure
    # (0.01 day); the derived Q is correspondingly imprecise but the Q + Vp pair
    # still reproduces the reported Vss to within <1%.
    # lq: Q derived from CL, V1, alpha_HL, beta_HL via the standard two-compartment relations.
    lcl <- log(0.02070); label("Clearance (L/day/kg)") # Luo 2017 Table 2, Fuc-SA- CL
    lvc <- log(0.02080); label("Central volume V1 (L/kg)") # Luo 2017 Table 2, Fuc-SA- V1
    lq <- log(1.32770); label("Intercompartmental clearance Q (L/day/kg)") # Luo 2017 Table 2, derived from CL/V1/alpha_HL/beta_HL
    lvp <- log(0.29270); label("Peripheral volume Vp (L/kg)") # Luo 2017 Table 2, Fuc-SA- Vss - V1
  })

  model({
    cl <- exp(lcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)
    Cc <- linCmt()
  })
}
