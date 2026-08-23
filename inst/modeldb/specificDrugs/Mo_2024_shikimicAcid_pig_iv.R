Mo_2024_shikimicAcid_pig_iv <- function() {
  description <- paste(
    "Preclinical (pig). Two-compartment intravenous pharmacokinetic model for",
    "shikimic acid (SA) in growing Landrace x Large White pigs, coupled to five",
    "direct-effect sigmoid-Emax pharmacodynamic models for the immune-enhancing",
    "effect of SA on plasma complement and immunoglobulin (Mo 2024). SA was",
    "given as a single 2 mg/kg injection into the marginal ear vein. The paper",
    "reports the intravenous disposition as the biexponential C = A*exp(-alpha*t)",
    "+ B*exp(-beta*t) (Table 10); the central volume, clearance,",
    "intercompartmental clearance and peripheral volume carried here are derived",
    "from the mean A, alpha, B and beta and reproduce that equation. The",
    "pharmacodynamic layer links the SA plasma concentration to the absolute",
    "change from the predose (0 h) level of complement C3 and C4 and of",
    "immunoglobulin G, A and M through E = Emax * C^gamma / (EC50^gamma +",
    "C^gamma) (eq 2; parameters in Table 11), with no effect compartment because",
    "Mo 2024 found no hysteresis between concentration and effect. All",
    "disposition and effect states are expressed per kg body weight (volumes in",
    "mL/kg, clearances in mL/h/kg, amounts in ng/kg). Mo 2024 fitted each pig",
    "individually in Phoenix WinNonlin and reported only the mean and SD of the",
    "individual estimates, so no between-subject variability or residual-error",
    "model is available; every parameter is fixed at the published mean and the",
    "residual SDs are fixed at zero. The companion intragastric model is",
    "modellib('Mo_2024_shikimicAcid_pig_oral')."
  )
  reference <- paste(
    "Mo K, Shen Y, Su D, Lv L, Du J, Ding H, Huang X (2024).",
    "Pharmacokinetic-Pharmacodynamic Modeling of the Immune-Enhancing Effect of",
    "Shikimic Acid in Growing Pigs. J Agric Food Chem 72:26224-26235.",
    "doi:10.1021/acs.jafc.4c09250"
  )
  vignette <- "Mo_2024_shikimicAcid"
  units <- list(time = "h", dosing = "ng/kg", concentration = "ng/mL")

  compartmentData <- list(
    central     = list(analyte = "shikimic acid", units = "ng", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "shikimic acid", units = "ng", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species       = "pig (Landrace x Large White)",
    n_subjects    = 6,
    n_studies     = 1,
    age_median    = "approximately 90 days",
    sex           = "male and female",
    disease_state = "healthy growing pigs",
    dose_range    = "shikimic acid 2 mg/kg as a single marginal-ear-vein injection",
    regions       = "China",
    notes         = paste(
      "Mo 2024 Section 2.2 and Table 1: six growing pigs of comparable health",
      "status and genetic background, approximately 90 days old, allocated to a",
      "two-period two-sequence crossover (sequence A n = 3, sequence B n = 3)",
      "with a 7-day washout; every pig therefore contributes both an",
      "intragastric and an intravenous profile. The same six pigs support the",
      "companion model Mo_2024_shikimicAcid_pig_oral. Body weights are given in",
      "supplementary Table S1, which is not part of the open-access record; the",
      "paper states only that weight did not differ between periods (P = 0.320).",
      "Plasma was sampled predose and at 0.083, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4,",
      "6, 8, 12 and 24 h (Section 2.3); complement and immunoglobulin were",
      "measured by ELISA at 0, 0.5, 1, 1.5, 3, 6 and 12 h (Section 2.5). Both",
      "sexes were studied and the authors note sex may affect SA disposition,",
      "but no sex covariate was fitted. No race/ethnicity data apply."
    )
  )

  ini({
    # Two-compartment intravenous disposition. Mo 2024 Table 10 reports the
    # biexponential coefficients and exponents of C = A*exp(-alpha*t) +
    # B*exp(-beta*t) (mean of the six individual fits): A = 7227.76 ng/mL,
    # alpha = 1.56 /h, B = 941.91 ng/mL, beta = 0.23 /h, for the 2 mg/kg
    # (2e6 ng/kg) intravenous dose. The four parameters below are the standard
    # transformation of that mean biexponential:
    #   vc  = Dose / (A + B)
    #   k21 = (A*beta + B*alpha) / (A + B); k10 = alpha*beta / k21
    #   k12 = alpha + beta - k21 - k10
    #   cl  = k10*vc; q = k12*vc; vp = q/k21
    # Solving these values back to a biexponential returns A = 7227.65 ng/mL,
    # alpha = 1.5600 /h, B = 941.95 ng/mL, beta = 0.2300 /h, i.e. Table 10 to
    # rounding. They also reproduce the secondary parameters Mo 2024 tabulates
    # independently: cl 229.14 vs 233.20 mL/h/kg reported, vc + vp 545.4 vs Vss
    # 574.10 mL/kg reported, and vp 300.58 vs the 323.11 mL/kg reported in the
    # "V_Z" column of Table 10 (that column is the peripheral volume, not the
    # terminal-phase volume: Vss - V_Z equals Dose/(A + B) for all six pigs to
    # within 0.01 mL/kg, whereas cl/beta runs 660-1913 mL/kg, i.e. 2.3-4.9 times
    # the V_Z of the same pig).
    lvc <- fixed(log(244.81)); label("Central volume of distribution (mL/kg)")  # derived from Table 10 mean A, B (Dose/(A+B))
    lcl <- fixed(log(229.14)); label("Elimination clearance (mL/h/kg)")  # derived from Table 10 mean A, alpha, B, beta
    lq <- fixed(log(115.23)); label("Intercompartmental clearance (mL/h/kg)")  # derived from Table 10 mean A, alpha, B, beta
    lvp <- fixed(log(300.58)); label("Peripheral volume of distribution (mL/kg)")  # derived from Table 10 mean A, alpha, B, beta

    # Direct-effect sigmoid-Emax models for the absolute change from the
    # predose (0 h) plasma level, E = Emax * C^gamma / (EC50^gamma + C^gamma)
    # (Mo 2024 eq 2). Parameters are the mean of the six individual intravenous
    # fits in Table 11. Table 11 expresses every Emax in ug/mL even where the
    # assay is reported in other units in Figure 5 (C4 in ng/mL, IgG in g/L);
    # the ug/mL values are carried here as tabulated, so dC4 is returned in
    # ug/mL (0.064 ug/mL = 64 ng/mL) and dIgG in ug/mL (9801.86 ug/mL =
    # 9.80 g/L). C in the equation is the SA plasma concentration in ng/mL,
    # matching the EC50 units.
    lemax_c3 <- fixed(log(352.95)); label("Maximum effect Emax on the change in complement C3 (ug/mL)")  # Table 11, C3
    lec50_c3 <- fixed(log(700.45)); label("Half-maximal SA concentration EC50 for complement C3 (ng/mL)")  # Table 11, C3
    lhill_c3 <- fixed(log(3.77)); label("Sigmoidicity gamma for complement C3 (unitless)")  # Table 11, C3

    lemax_c4 <- fixed(log(0.064)); label("Maximum effect Emax on the change in complement C4 (ug/mL)")  # Table 11, C4
    lec50_c4 <- fixed(log(836.61)); label("Half-maximal SA concentration EC50 for complement C4 (ng/mL)")  # Table 11, C4
    lhill_c4 <- fixed(log(3.89)); label("Sigmoidicity gamma for complement C4 (unitless)")  # Table 11, C4

    lemax_igg <- fixed(log(9801.86)); label("Maximum effect Emax on the change in immunoglobulin G (ug/mL)")  # Table 11, IgG
    lec50_igg <- fixed(log(514.22)); label("Half-maximal SA concentration EC50 for immunoglobulin G (ng/mL)")  # Table 11, IgG
    lhill_igg <- fixed(log(5.15)); label("Sigmoidicity gamma for immunoglobulin G (unitless)")  # Table 11, IgG

    lemax_iga <- fixed(log(469.30)); label("Maximum effect Emax on the change in immunoglobulin A (ug/mL)")  # Table 11, IgA
    lec50_iga <- fixed(log(561.88)); label("Half-maximal SA concentration EC50 for immunoglobulin A (ng/mL)")  # Table 11, IgA
    lhill_iga <- fixed(log(6.97)); label("Sigmoidicity gamma for immunoglobulin A (unitless)")  # Table 11, IgA

    lemax_igm <- fixed(log(1252.05)); label("Maximum effect Emax on the change in immunoglobulin M (ug/mL)")  # Table 11, IgM
    lec50_igm <- fixed(log(686.24)); label("Half-maximal SA concentration EC50 for immunoglobulin M (ng/mL)")  # Table 11, IgM
    lhill_igm <- fixed(log(3.75)); label("Sigmoidicity gamma for immunoglobulin M (unitless)")  # Table 11, IgM

    addSd <- fixed(0); label("Additive residual SD on shikimic acid plasma concentration (ng/mL; not reported)")  # Mo 2024 reports no residual-error model
    addSd_dC3 <- fixed(0); label("Additive residual SD on the change in complement C3 (ug/mL; not reported)")  # Mo 2024 reports no residual-error model
    addSd_dC4 <- fixed(0); label("Additive residual SD on the change in complement C4 (ug/mL; not reported)")  # Mo 2024 reports no residual-error model
    addSd_dIgG <- fixed(0); label("Additive residual SD on the change in immunoglobulin G (ug/mL; not reported)")  # Mo 2024 reports no residual-error model
    addSd_dIgA <- fixed(0); label("Additive residual SD on the change in immunoglobulin A (ug/mL; not reported)")  # Mo 2024 reports no residual-error model
    addSd_dIgM <- fixed(0); label("Additive residual SD on the change in immunoglobulin M (ug/mL; not reported)")  # Mo 2024 reports no residual-error model
  })

  model({
    vc <- exp(lvc)
    cl <- exp(lcl)
    q <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central)     <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc

    # Absolute change from the predose level of each immune marker; a direct
    # effect of the plasma concentration, with no effect compartment (Mo 2024
    # Section 3.5 and Discussion: no hysteresis, Emax and Cmax synchronous).
    emax_c3 <- exp(lemax_c3)
    ec50_c3 <- exp(lec50_c3)
    hill_c3 <- exp(lhill_c3)
    dC3 <- emax_c3 * Cc^hill_c3 / (ec50_c3^hill_c3 + Cc^hill_c3)

    emax_c4 <- exp(lemax_c4)
    ec50_c4 <- exp(lec50_c4)
    hill_c4 <- exp(lhill_c4)
    dC4 <- emax_c4 * Cc^hill_c4 / (ec50_c4^hill_c4 + Cc^hill_c4)

    emax_igg <- exp(lemax_igg)
    ec50_igg <- exp(lec50_igg)
    hill_igg <- exp(lhill_igg)
    dIgG <- emax_igg * Cc^hill_igg / (ec50_igg^hill_igg + Cc^hill_igg)

    emax_iga <- exp(lemax_iga)
    ec50_iga <- exp(lec50_iga)
    hill_iga <- exp(lhill_iga)
    dIgA <- emax_iga * Cc^hill_iga / (ec50_iga^hill_iga + Cc^hill_iga)

    emax_igm <- exp(lemax_igm)
    ec50_igm <- exp(lec50_igm)
    hill_igm <- exp(lhill_igm)
    dIgM <- emax_igm * Cc^hill_igm / (ec50_igm^hill_igm + Cc^hill_igm)

    Cc ~ add(addSd)
    dC3 ~ add(addSd_dC3)
    dC4 ~ add(addSd_dC4)
    dIgG ~ add(addSd_dIgG)
    dIgA ~ add(addSd_dIgA)
    dIgM ~ add(addSd_dIgM)
  })
}
