Dias_2022_tobramycin_rat <- function() {
  description <- paste(
    "Preclinical (rat, Wistar male). Three-compartment population PK model",
    "for intravenous tobramycin in healthy rats and in rats with acute or",
    "chronic biofilm-forming Pseudomonas aeruginosa lung infection, fit",
    "jointly to total plasma concentrations and to unbound lung and",
    "epithelial-lining-fluid (ELF) concentrations measured by",
    "microdialysis. The ODE states carry UNBOUND drug amounts: central",
    "(plasma), peripheral1 and a third lung compartment (the lung",
    "interstitial space sampled by the lung microdialysis probe). Unbound",
    "ELF concentration is the unbound lung concentration times a",
    "distribution factor (Dfactor), and total plasma is unbound plasma",
    "divided by the unbound fraction 0.89. Chronic infection (alginate",
    "beads carrying P. aeruginosa ATCC 27853) has its own clearance and",
    "central volume; every intratracheally inoculated group (acute",
    "PA14 infection, chronic infection and sterile blank alginate beads)",
    "shares a larger lung volume than healthy rats.",
    sep = " "
  )
  reference <- paste(
    "Dias BB, Carreno F, Helfer VE, Garzella PMB, de Lima DMF, Barreto F,",
    "de Araujo BV, Dalla Costa T. Probability of Target Attainment of",
    "Tobramycin Treatment in Acute and Chronic Pseudomonas aeruginosa Lung",
    "Infection Based on Preclinical Population Pharmacokinetic Modeling.",
    "Pharmaceutics. 2022;14(6):1237. doi:10.3390/pharmaceutics14061237.",
    "PMCID: PMC9228144. Parameter estimates: Table 1. Structural ODEs:",
    "Equations 1-3. Lung, ELF and unbound-fraction observation equations:",
    "Results text after Equation 3. Group sizes: Supplementary Table S2.",
    "Human-to-rat dose scaling: Supplementary Equations S1-S2 and Table S1.",
    sep = " "
  )
  vignette <- "Dias_2022_tobramycin_rat"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "tobramycin (unbound)", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "tobramycin (unbound)", units = "mg", specimen = "tissue", verified = TRUE),
    lung = list(analyte = "tobramycin (unbound)", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    DIS_PSEUDOMONAS_LUNG_ACUTE = list(
      description = "Acute Pseudomonas aeruginosa lung infection (preclinical): 1 = rat inoculated intratracheally with planktonic P. aeruginosa PA14 and studied 7 days later; 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy, not inoculated)",
      source_name = "acute infection (acutely infected group)",
      notes = paste(
        "Inoculum 100 uL of a 10^9 CFU/mL PA14 suspension; PK experiments 7 days after inoculation (Section 2.3).",
        "Enters the model only through the shared inoculated-lung volume (V3infected); clearance and central volume are those of healthy rats.",
        "Mutually exclusive with DIS_PSEUDOMONAS_LUNG_CHRONIC and ALGINATE_BEAD_BLANK; all three are 0 for a healthy rat."
      )
    ),
    DIS_PSEUDOMONAS_LUNG_CHRONIC = list(
      description = "Chronic Pseudomonas aeruginosa lung infection (preclinical): 1 = rat inoculated intratracheally with alginate beads carrying P. aeruginosa ATCC 27853 and studied 14 days later; 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not chronically infected)",
      source_name = "chronic infection (chronically infected group)",
      notes = paste(
        "Alginate-bead model mimicking the chronic mucoid infection of cystic-fibrosis patients; 50 uL of beads, PK experiments 14 days after inoculation (Section 2.3).",
        "Selects the chronic-infection clearance and central volume (CLchronic, V1chronic) and the inoculated-lung volume (V3infected).",
        "Mutually exclusive with DIS_PSEUDOMONAS_LUNG_ACUTE and ALGINATE_BEAD_BLANK."
      )
    ),
    ALGINATE_BEAD_BLANK = list(
      description = "Sterile alginate-bead control (preclinical): 1 = rat inoculated intratracheally with blank alginate beads that carry no bacteria; 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no blank-bead inoculation)",
      source_name = "blank-bead group",
      notes = paste(
        "50 uL of beads prepared by the same procedure as the chronic-infection beads but without bacteria (Section 2.3); the control for the effect of alginate itself.",
        "Enters the model only through the shared inoculated-lung volume (V3infected); plasma PK is that of healthy rats.",
        "Mutually exclusive with DIS_PSEUDOMONAS_LUNG_ACUTE and DIS_PSEUDOMONAS_LUNG_CHRONIC."
      )
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened in the exploratory analysis of the base model (Supplementary Figure S4) and not retained because it showed no correlation with the model parameters (Section 2.5). Rats weighed 200-250 g, and the parameters are absolute (L, L/h) for a rat of that size."
    )
  )

  population <- list(
    species = "rat (Wistar, male)",
    n_subjects = 71L,
    n_studies = 1L,
    age_range = "Adult (age not reported)",
    weight_range = "200-250 g",
    sex_female_pct = 0,
    disease_state = paste(
      "Four groups: healthy; acute P. aeruginosa PA14 lung infection (7 days",
      "after intratracheal inoculation); chronic P. aeruginosa ATCC 27853",
      "lung infection established with bacteria-laden alginate beads (14",
      "days after inoculation); and a blank-bead control inoculated with",
      "sterile alginate beads."
    ),
    dose_range = "Single 10 mg/kg intravenous bolus through the femoral vein.",
    regions = "Brazil (Federal University of Rio Grande do Sul, Porto Alegre)",
    notes = paste(
      "1231 observations from 71 rats. Supplementary Table S2",
      "(animals/observations): plasma healthy 6/68, acute 5/59, chronic",
      "6/72, blank bead 4/40; lung microdialysis healthy 6/142, acute 5/120,",
      "chronic 7/106, blank bead 5/111; ELF microdialysis healthy 12/231,",
      "acute 6/91, chronic 3/61, blank bead 5/111. Plasma sampled at 0.08",
      "to 12 h; microdialysate collected every 30 min for 12 h and corrected",
      "for the in vivo probe relative recovery (27.7% lung, 26.6% ELF).",
      "NONMEM 7.4, FOCE-I."
    )
  )

  ini({
    # Structural parameters: Table 1 (estimate, %RSE). Values are absolute
    # for a 200-250 g rat; body weight was not retained as a covariate.
    # CL and V1 were "estimated separately for the chronic infected group"
    # (Results), so each carries a stratum suffix.
    lcl_nonchronic <- log(0.047); label("Clearance, healthy / acute / blank-bead rats (L/h)") # Table 1: CL = 0.047 L/h (RSE 12%)
    lcl_chronic <- log(0.085); label("Clearance, chronically infected rats (L/h)") # Table 1: CLchronic = 0.085 L/h (RSE 23%)
    lvc_nonchronic <- log(0.055); label("Central volume, healthy / acute / blank-bead rats (L)") # Table 1: V1 = 0.055 L (RSE 16%)
    lvc_chronic <- log(0.323); label("Central volume, chronically infected rats (L)") # Table 1: V1chronic = 0.323 L (RSE 17%)
    lq <- log(0.030); label("Intercompartmental clearance, central to peripheral (L/h)") # Table 1: Q1 = 0.030 L/h (RSE 8%)
    lvp <- log(0.154); label("Peripheral volume (L)") # Table 1: V2 = 0.154 L (RSE 16%)
    lq_lung <- log(0.370); label("Intercompartmental clearance, central to lung (L/h)") # Table 1: Q2 = 0.370 L/h (RSE 5%)
    lv_lung_healthy <- log(0.083); label("Lung compartment volume, healthy rats (L)") # Table 1: V3 = 0.083 L (RSE 29%)
    lv_lung_inoc <- log(0.130); label("Lung compartment volume, acute / chronic / blank-bead rats (L)") # Table 1: V3infected = 0.130 L (RSE 22%)
    lr_elf_lung <- log(0.36); label("Distribution factor, unbound ELF to unbound lung concentration ratio (unitless)") # Table 1: Dfactor = 0.36 (RSE 10%)
    fu <- fixed(0.89); label("Unbound fraction of tobramycin in rat plasma (fraction)") # Results after Eq 3: estimates divided by 0.89, the unbound fraction; Methods 2.5: protein binding 11%

    # IIV: Table 1 reports %CV for exponential IIV; omega^2 = log(CV^2 + 1).
    etalcl ~ 0.5339 # Table 1: omega CL = 84 %CV -> log(0.84^2 + 1)
    etalvc ~ 0.3075 # Table 1: omega V1 = 60 %CV -> log(0.60^2 + 1)
    etalv_lung ~ 1.0280 # Table 1: omega V3 = 134 %CV -> log(1.34^2 + 1)
    etalr_elf_lung ~ 0.1697 # Table 1: omega Dfactor = 43 %CV -> log(0.43^2 + 1)

    # Residual error: log-additive, "described separately for plasma and
    # microdialysate data" (Methods 2.5). The paper estimated ONE
    # microdialysis SD for lung and ELF together; nlmixr2 needs a separate
    # parameter per endpoint, so the single estimate is carried twice.
    expSd <- 0.152; label("Log-additive residual SD, total plasma (log scale)") # Table 1: plasma log-additive error = 0.152 (RSE 5%)
    expSd_Clung <- 0.313; label("Log-additive residual SD, lung microdialysate (log scale)") # Table 1: microdialysis log-additive error = 0.313 (RSE 3%), shared with ELF
    expSd_Celf <- 0.313; label("Log-additive residual SD, ELF microdialysate (log scale)") # Table 1: microdialysis log-additive error = 0.313 (RSE 3%), shared with lung
  })

  model({
    # Group indicators are mutually exclusive; every intratracheally
    # inoculated group shares V3infected (Results: "Acute, chronic, and
    # blank-bead were included as covariates in V3").
    inoc <- DIS_PSEUDOMONAS_LUNG_ACUTE + DIS_PSEUDOMONAS_LUNG_CHRONIC + ALGINATE_BEAD_BLANK

    cl <- exp(lcl_nonchronic * (1 - DIS_PSEUDOMONAS_LUNG_CHRONIC) + lcl_chronic * DIS_PSEUDOMONAS_LUNG_CHRONIC + etalcl)
    vc <- exp(lvc_nonchronic * (1 - DIS_PSEUDOMONAS_LUNG_CHRONIC) + lvc_chronic * DIS_PSEUDOMONAS_LUNG_CHRONIC + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)
    q_lung <- exp(lq_lung)
    v_lung <- exp(lv_lung_healthy * (1 - inoc) + lv_lung_inoc * inoc + etalv_lung)
    r_elf_lung <- exp(lr_elf_lung + etalr_elf_lung)

    # Equations 1-3 (unbound amounts).
    d/dt(central) <- -(q / vc + q_lung / vc + cl / vc) * central + q / vp * peripheral1 + q_lung / v_lung * lung
    d/dt(peripheral1) <- -q / vp * peripheral1 + q / vc * central
    d/dt(lung) <- -q_lung / v_lung * lung + q_lung / vc * central

    # Unbound plasma, unbound lung (A3/V3) and unbound ELF ((A3/V3) * Dfactor).
    Cu <- central / vc
    Clung <- lung / v_lung
    Celf <- Clung * r_elf_lung
    # Observed plasma concentrations are total: "the estimates were divided
    # by 0.89, the unbound fraction of the drug" (Results after Eq 3).
    Cc <- Cu / fu

    Cc ~ lnorm(expSd)
    Clung ~ lnorm(expSd_Clung)
    Celf ~ lnorm(expSd_Celf)
  })
}
