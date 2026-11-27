Quinney_2013_midazolam_ketoconazole_ddi_pbpk <- function() {
  description <- paste(
    "PBPK (semi-physiological, CYP3A drug-drug interaction). Midazolam",
    "(victim) and ketoconazole (perpetrator) solved in one system: the",
    "semi-physiological model of Chien et al. 2006 as implemented in the",
    "R code published with Quinney 2013, with the in vitro CYP3A",
    "inhibition constants of Quinney 2013 (mixed inhibition with",
    "independent substrate, inhibitor and metabolite binding, M-BIM:",
    "Kic 19.2 nmol/L, Kin 39.8 nmol/L). Midazolam has a gut lumen, gut",
    "wall, portal vein, liver, central and peripheral compartment; hepatic",
    "and gut-wall CYP3A metabolism is Michaelis-Menten and is inhibited by",
    "unbound ketoconazole in a hepatic enzyme-site compartment and in the",
    "gut lumen. Ketoconazole is absorbed straight into its liver",
    "compartment, distributes to a central and peripheral compartment and",
    "exchanges saturably with the enzyme site. The competitive and",
    "noncompetitive variants the paper compares are special cases of the",
    "shipped mixed form (kin very large; kic equal to kin), so the",
    "vignette reproduces all three by overriding kic and kin. Deterministic:",
    "the paper simulated a single typical subject, so there is no",
    "between-subject variability and no residual-error model."
  )
  reference <- paste(
    "Quinney SK, Knopp S, Chang C, Hall SD, Li L. Integration of in vitro",
    "binding mechanism into the semiphysiologically based pharmacokinetic",
    "interaction model between ketoconazole and midazolam.",
    "CPT Pharmacometrics Syst Pharmacol. 2013;2(9):e75.",
    "doi:10.1038/psp.2013.50.",
    "The ODEs and every structural constant are transcribed from the",
    "authors' R code in the Supplementary Data (PSP-2013-0039-s06.txt, the",
    "three inhibition-model ODE functions; PSP-2013-0039-s07.txt, the",
    "parameter vector and dosing script). Kic and Kin are Table 3 (BIM",
    "columns). The structural model and the parameter values labelled",
    "'Chien, 2006' in that code are from Chien JY, Lucksiri A, Ernest CS,",
    "Gorski JC, Wrighton SA, Hall SD. Drug Metab Dispos. 2006;34:1208-1219.",
    "doi:10.1124/dmd.105.008730.",
    sep = " "
  )
  vignette <- "Quinney_2013_midazolam_ketoconazole"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")
  # Oral midazolam into `depot`, intravenous midazolam into `central`, oral
  # ketoconazole into `depot_ketoconazole`.
  dosing <- c("depot", "central", "depot_ketoconazole")

  covariateData <- list()

  # The midazolam states use the canonical bare names where one exists
  # (depot, liver, central, peripheral1). The gut wall and portal vein have no
  # bare canonical, and the ketoconazole states carry a `_ketoconazole` suffix
  # because the two drugs are solved in one system (the Nishiyama 2019
  # metformin-cimetidine precedent). `enzyme_site_ketoconazole` is the
  # ketoconazole amount in the hepatic CYP3A enzyme compartment ('Ae' in the
  # authors' code), whose unbound concentration inhibits hepatic midazolam
  # metabolism.
  paper_specific_compartments <- c(
    "gut_wall",
    "portal_vein",
    "depot_ketoconazole",
    "liver_ketoconazole",
    "enzyme_site_ketoconazole",
    "central_ketoconazole",
    "peripheral1_ketoconazole"
  )

  compartmentData <- list(
    depot = list(analyte = "midazolam", units = "mg", specimen = "administration site", verified = TRUE),
    gut_wall = list(analyte = "midazolam", units = "mg", specimen = "tissue", verified = TRUE),
    portal_vein = list(analyte = "midazolam", units = "mg", specimen = "whole blood", verified = TRUE),
    liver = list(analyte = "midazolam", units = "mg", specimen = "tissue", verified = TRUE),
    central = list(analyte = "midazolam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "midazolam", units = "mg", specimen = "tissue", verified = TRUE),
    depot_ketoconazole = list(
      analyte = "ketoconazole",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    liver_ketoconazole = list(analyte = "ketoconazole", units = "mg", specimen = "tissue", verified = TRUE),
    enzyme_site_ketoconazole = list(analyte = "ketoconazole", units = "mg", specimen = "tissue", verified = TRUE),
    central_ketoconazole = list(analyte = "ketoconazole", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_ketoconazole = list(analyte = "ketoconazole", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 1L,
    disease_state = "healthy adults (typical-subject simulation)",
    dose_range = paste(
      "Ketoconazole 200 mg orally once daily for 6 days; midazolam 1 mg",
      "intravenous bolus alone and with the sixth ketoconazole dose"
    ),
    notes = paste(
      "Quinney 2013 simulated a single typical subject with the mean",
      "parameters ('A single subject was simulated for each inhibition",
      "type using the mean parameters'). The observed comparator is the",
      "healthy-volunteer intravenous midazolam / ketoconazole study of",
      "Lucksiri 2005 (PhD thesis, Purdue University; Quinney 2013 ref. 16),",
      "with concentrations normalised to a 1 mg midazolam dose. The in vitro",
      "inhibition constants come from a re-analysis of six published human",
      "liver microsome studies (Quinney 2013 Table 2, studies A-F)."
    )
  )

  ini({
    # ---------------------------------------------------------------
    # Midazolam disposition (s07 parameter vector, para[1]-para[8])
    # ---------------------------------------------------------------
    lka <- fixed(log(1.17))
    label("Midazolam oral absorption rate constant, gut lumen to gut wall (1/h)") # s07 para[1] ka = 1.17 h-1
    lkpv <- fixed(log(1400))
    label("Midazolam transfer rate constant, gut wall to portal vein (1/h)") # s07 para[2] kpv = 1400 h-1
    lvc <- fixed(log(43))
    label("Midazolam central volume V4 (L)") # s07 para[3] V4 = 43 L
    lvp <- fixed(log(88.4))
    label("Midazolam peripheral volume V5 (L)") # s07 para[4] V5 = 88.4 L, Chien 2006
    lq <- fixed(log(14.73))
    label("Midazolam intercompartmental clearance Q4 (L/h)") # s07 para[5] Q4 = 14.73 L/h
    vmax_gut <- fixed(2549)
    label("Midazolam gut-wall CYP3A maximum metabolic rate Vmaxg (mg/h)") # s07 para[6] Vmaxg = 2549 mg/h, Chien 2006
    vmax_liver <- fixed(2549 + 206)
    label("Midazolam hepatic maximum metabolic rate Vmaxh (mg/h)") # s07 para[7] Vmaxh = 2549+206 mg/h, Chien 2006
    km <- fixed(5.8 * 0.362)
    label("Midazolam Michaelis constant Km (mg/L)") # s07 para[8] km = 5.8*0.362 mg/L, Chien 2006
    fm <- fixed(0.93)
    label("Fraction of hepatic midazolam metabolism by CYP3A, fm (fraction)") # s07 para[18] fm = 0.93
    fu <- fixed(0.04)
    label("Midazolam fraction unbound in liver, fu,MDZ (fraction)") # s07 para[22] fum = 0.04, Chien 2006; Methods 'assumed to be 0.04'
    lcl_renal <- fixed(log(0.06))
    label("Midazolam renal clearance (L/h)") # s07 para[33] CLren = 0.06 L/h, Chien 2006

    # ---------------------------------------------------------------
    # Ketoconazole disposition (s07 para[9]-para[16], para[24]-para[25],
    # para[30], para[34]-para[35])
    # ---------------------------------------------------------------
    lka_ketoconazole <- fixed(log(1.10))
    label("Ketoconazole oral absorption rate constant, gut lumen to liver (1/h)") # s07 para[9] ka1 = 1.10, Chien 2006
    lvc_ketoconazole <- fixed(log(32.1))
    label("Ketoconazole central volume V51 (L)") # s07 para[10] V51 = 32.1 L, Chien 2006
    lvp_ketoconazole <- fixed(log(19.6))
    label("Ketoconazole peripheral volume V61 (L)") # s07 para[11] V61 = 19.6 L, Chien 2006
    lq_ketoconazole <- fixed(log(2.39))
    label("Ketoconazole intercompartmental clearance Q5 (L/h)") # s07 para[14] Q5 = 2.39 L/h, Chien 2006
    vmax_ketoconazole <- fixed(41.1)
    label("Ketoconazole hepatic maximum metabolic rate Vmaxh1 (mg/h)") # s07 para[15] Vmaxh1 = 41.1 mg/h, Chien 2006
    km_ketoconazole <- fixed(0.810)
    label("Ketoconazole hepatic Michaelis constant km1 (mg/L)") # s07 para[16] km1 = 0.810 mg/L, Chien 2006
    fu_cl_ketoconazole <- fixed(1)
    label("Ketoconazole fraction unbound applied to its hepatic clearance, fukCL (fraction)") # s07 para[24] fukCL = 1, Chien 2006
    fu_eff_ketoconazole <- fixed(1)
    label("Ketoconazole fraction unbound applied to liver-to-enzyme-site transfer, fukEff (fraction)") # s07 para[25] fukEff = 1, Chien 2006
    k_enzyme_ketoconazole <- fixed(1)
    label("Ketoconazole first-order transfer rate constant, liver to enzyme site k2e (1/h)") # s07 para[13] k2e = 1 h-1, Chien 2006
    tmax_enzyme_ketoconazole <- fixed(1)
    label("Ketoconazole maximum saturable return rate, enzyme site to liver Tmax (mg/h)") # s07 para[34] Tmax = 1 mg/h, Chien 2006
    keff_enzyme_ketoconazole <- fixed(5)
    label("Ketoconazole half-saturation concentration of the enzyme-site return keff (mg/L)") # s07 para[12] keff = 5 mg/L, Chien 2006
    fu_inh_ketoconazole <- fixed(0.01)
    label("Ketoconazole fraction unbound at the hepatic enzyme site, fu,KTZ (fraction)") # s07 para[23] fukddi = 0.01, Chien 2006; Methods 'assumed to be ... 0.01'

    # ---------------------------------------------------------------
    # CYP3A inhibition constants -- Table 3, M-BIM column (unbound,
    # nmol/L); converted to mg/L with the ketoconazole molecular weight
    # the code uses (s07 para[37], para[38]: x*531.4/1000/1000).
    # ---------------------------------------------------------------
    kic <- fixed(19.2)
    label("Ketoconazole competitive inhibition constant Kic (nmol/L)") # Table 3 BIM, M: Kic 19.2 nmol/L; s07 para[38]
    kin <- fixed(39.8)
    label("Ketoconazole noncompetitive inhibition constant Kin (nmol/L)") # Table 3 BIM, M: Kin 39.8 nmol/L; s07 para[37]
    mw_ketoconazole <- fixed(531.4)
    label("Ketoconazole molecular weight used for the nmol/L to mg/L conversion (g/mol)") # s07 para[37]-[38]: 531.4

    # ---------------------------------------------------------------
    # Physiology (s07 para[26]-para[32])
    # ---------------------------------------------------------------
    v_portal <- fixed(0.07)
    label("Portal vein volume Vpv (L)") # s07 para[26] Vpv = 0.07 L, Ito 2003
    q_portal <- fixed(67.5)
    label("Portal vein blood flow Qpv (L/h)") # s07 para[27] Qpv = 67.5 L/h
    qh <- fixed(90)
    label("Hepatic blood flow Qh (L/h)") # s07 para[28] Qh = 90 L/h
    v_liver <- fixed(1.5)
    label("Liver volume for midazolam Vh (L)") # s07 para[29] Vh = 1.5 L, Ito 2003
    v_liver_ketoconazole <- fixed(0.825)
    label("Liver volume for ketoconazole Vhk (L)") # s07 para[30] Vhk = 0.825 L, Chien 2006
    v_enzyme_ketoconazole <- fixed(0.825)
    label("Hepatic enzyme-site volume Ve (L)") # s07 para[35] Ve = 0.825 L, Chien 2006
    v_gut_wall <- fixed(1)
    label("Gut-wall volume Vgw (L)") # s07 para[32] Vgw = 1 L
  })
  model({
    ka <- exp(lka)
    kpv <- exp(lkpv)
    vc <- exp(lvc)
    vp <- exp(lvp)
    q <- exp(lq)
    ka_ketoconazole <- exp(lka_ketoconazole)
    vc_ketoconazole <- exp(lvc_ketoconazole)
    vp_ketoconazole <- exp(lvp_ketoconazole)
    q_ketoconazole <- exp(lq_ketoconazole)
    cl_renal <- exp(lcl_renal)

    # Inhibition constants in mg/L (s07: Ki * 531.4 / 1000 / 1000).
    kic_mgl <- kic * mw_ketoconazole / 1e6
    kin_mgl <- kin * mw_ketoconazole / 1e6

    # Inhibitor concentrations. Gut wall: ketoconazole amount still in the
    # gut lumen divided by the gut-wall volume, with no binding correction
    # (s06: y[7] / Vgw). Liver: unbound ketoconazole at the enzyme site
    # (s06: fukddi * y[10] / Ve).
    c_inh_gut <- depot_ketoconazole / v_gut_wall
    c_inh_liver <- fu_inh_ketoconazole * enzyme_site_ketoconazole / v_enzyme_ketoconazole

    # Mixed-inhibition Michaelis-Menten metabolism (s06
    # mixed.inhibition.model). The CYP3A fraction fm is inhibited; the
    # remaining 1 - fm of hepatic metabolism is not. In the code the
    # unbound fraction multiplies only the substrate term of the
    # denominator, i.e. Km is an unbound-concentration constant.
    met_gut <- vmax_gut * gut_wall /
      (km * v_gut_wall * (1 + c_inh_gut / kic_mgl) + gut_wall * (1 + c_inh_gut / kin_mgl))
    met_liver_cyp3a <- fm * vmax_liver * fu * liver /
      (km * v_liver * (1 + c_inh_liver / kic_mgl) + fu * liver * (1 + c_inh_liver / kin_mgl))
    met_liver_other <- (1 - fm) * vmax_liver * fu * liver / (km * v_liver + fu * liver)

    # Midazolam (s06 y[1]-y[6]). The code multiplies these derivatives by
    # the route flags modose / midose; dosing the depot (oral) or central
    # (intravenous) compartment here is equivalent.
    d/dt(depot) <- -ka * depot
    d/dt(gut_wall) <- ka * depot - kpv * gut_wall - met_gut
    d/dt(portal_vein) <- kpv * gut_wall - q_portal * (portal_vein / v_portal - central / vc)
    d/dt(liver) <- q_portal * (portal_vein / v_portal - central / vc) +
      qh * (central / vc - liver / v_liver) - met_liver_cyp3a - met_liver_other
    d/dt(central) <- -qh * (central / vc - liver / v_liver) -
      q * (central / vc - peripheral1 / vp) - cl_renal * central / vc
    d/dt(peripheral1) <- q * (central / vc - peripheral1 / vp)

    # Ketoconazole (s06 y[7], y[9]-y[12]; y[8] has a zero derivative in the
    # code and is omitted). The oral dose is absorbed directly into the
    # liver compartment.
    d/dt(depot_ketoconazole) <- -ka_ketoconazole * depot_ketoconazole
    d/dt(liver_ketoconazole) <- ka_ketoconazole * depot_ketoconazole +
      qh * (central_ketoconazole / vc_ketoconazole - liver_ketoconazole / v_liver_ketoconazole) -
      vmax_ketoconazole * fu_cl_ketoconazole * liver_ketoconazole /
        (v_liver_ketoconazole * km_ketoconazole + liver_ketoconazole * fu_cl_ketoconazole) -
      fu_eff_ketoconazole * k_enzyme_ketoconazole * liver_ketoconazole +
      tmax_enzyme_ketoconazole * enzyme_site_ketoconazole /
        (v_enzyme_ketoconazole * keff_enzyme_ketoconazole + enzyme_site_ketoconazole)
    d/dt(enzyme_site_ketoconazole) <- fu_eff_ketoconazole * k_enzyme_ketoconazole * liver_ketoconazole -
      tmax_enzyme_ketoconazole * enzyme_site_ketoconazole /
        (v_enzyme_ketoconazole * keff_enzyme_ketoconazole + enzyme_site_ketoconazole)
    d/dt(central_ketoconazole) <- -qh * (central_ketoconazole / vc_ketoconazole - liver_ketoconazole / v_liver_ketoconazole) -
      q_ketoconazole * (central_ketoconazole / vc_ketoconazole - peripheral1_ketoconazole / vp_ketoconazole)
    d/dt(peripheral1_ketoconazole) <- q_ketoconazole *
      (central_ketoconazole / vc_ketoconazole - peripheral1_ketoconazole / vp_ketoconazole)

    # Outputs in ug/L (mg/L x 1000), the units of Quinney 2013 Figure 3.
    Cc <- 1000 * central / vc
    Cc_ketoconazole <- 1000 * central_ketoconazole / vc_ketoconazole
  })
}
