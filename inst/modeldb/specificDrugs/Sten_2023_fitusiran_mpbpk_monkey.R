Sten_2023_fitusiran_mpbpk_monkey <- function() {
  description <- paste(
    "Preclinical (cynomolgus monkey). Minimal PBPK-PD model of GalNAc-conjugated siRNA disposition and gene silencing",
    "(Ayyar 2021 mPBPK-PD, developed on fitusiran data) in the corrected and unambiguously",
    "defined form published by Sten 2023 (ESM MATLAB implementation), monkey parameter set.",
    "18 ODEs: SC depot, plasma, lumped remainder tissue, kidney vascular / tissue / bound pool,",
    "liver vascular / interstitial / bound pool, ASGPR receptor binding with internalization,",
    "endosomal cleavage and recycling, endosomal escape to the cytoplasm, RISC loading, and",
    "RISC-stimulated degradation of the target (antithrombin) mRNA driving plasma antithrombin",
    "protein through an indirect-response (power gamma) translation term. Deterministic",
    "typical-value simulator: no between-subject variability or residual error.",
    sep = " "
  )
  reference <- paste(
    "Sten S, Cardilin T, Antonsson M, Gennemark P. Plasma Pharmacokinetics of",
    "N-Acetylgalactosamine-Conjugated Small-Interfering Ribonucleic Acids",
    "(GalNAc-Conjugated siRNAs). Clin Pharmacokinet. 2023;62(12):1661-1672.",
    "doi:10.1007/s40262-023-01314-7. Original model: Ayyar VS, Song D, Zheng S, Carpenter T,",
    "Heald DL. Minimal Physiologically Based Pharmacokinetic-Pharmacodynamic (mPBPK-PD) Model",
    "of N-Acetylgalactosamine-Conjugated Small Interfering RNA Disposition and Gene Silencing",
    "in Preclinical Species and Humans. J Pharmacol Exp Ther. 2021;379(2):134-146.",
    "doi:10.1124/jpet.121.000805",
    sep = " "
  )
  vignette <- "Sten_2023_galnacSirna"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Sten 2023 writes every state except the depot and the two bound pools as a
  # concentration, so the states keep the paper's units rather than amounts.
  paper_specific_compartments <- c(
    "remainder",
    "is_liver",
    "complex_endo",
    "target_endo",
    "cyto",
    "risc",
    "protein"
  )

  compartmentData <- list(
    depot = list(analyte = "GalNAc-conjugated siRNA", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "GalNAc-conjugated siRNA", units = "ng/mL", specimen = "plasma", verified = TRUE),
    remainder = list(analyte = "GalNAc-conjugated siRNA", units = "ng/mL", specimen = "tissue", verified = TRUE),
    kidney_vas = list(analyte = "GalNAc-conjugated siRNA", units = "ng/mL", specimen = "plasma", verified = TRUE),
    kidney = list(analyte = "GalNAc-conjugated siRNA", units = "ng/g", specimen = "tissue", verified = TRUE),
    kidney_deep = list(analyte = "GalNAc-conjugated siRNA", units = "ng", specimen = "tissue", verified = TRUE),
    liver_vas = list(analyte = "GalNAc-conjugated siRNA", units = "ng/mL", specimen = "plasma", verified = TRUE),
    is_liver = list(analyte = "GalNAc-conjugated siRNA", units = "nmol/L", specimen = "tissue", verified = TRUE),
    liver_deep = list(analyte = "GalNAc-conjugated siRNA", units = "nmol", specimen = "tissue", verified = TRUE),
    target = list(analyte = "free ASGPR (cell surface)", units = "nmol/L", specimen = "tissue", verified = TRUE),
    complex = list(
      analyte = "GalNAc-siRNA - ASGPR complex (cell surface)",
      units = "nmol/L",
      specimen = "tissue",
      verified = TRUE
    ),
    complex_endo = list(
      analyte = "GalNAc-siRNA - ASGPR complex (endosomal)",
      units = "nmol/L",
      specimen = "tissue",
      verified = TRUE
    ),
    target_endo = list(analyte = "free ASGPR (endosomal)", units = "nmol/L", specimen = "tissue", verified = TRUE),
    liver_endo = list(analyte = "free siRNA (endosomal)", units = "nmol/L", specimen = "tissue", verified = TRUE),
    cyto = list(analyte = "free siRNA (cytoplasmic)", units = "nmol/L", specimen = "tissue", verified = TRUE),
    risc = list(analyte = "RISC-loaded siRNA", units = "nmol/L", specimen = "tissue", verified = TRUE),
    mrna_liver = list(
      analyte = "antithrombin mRNA",
      units = "percent of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    protein = list(
      analyte = "antithrombin protein",
      units = "percent of baseline",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list()

  population <- list(
    species = "cynomolgus monkey",
    n_subjects = NA_integer_,
    n_studies = NA_integer_,
    weight_range = "6.2 kg (single typical body weight used by the ESM code for allometric scaling and mg/kg dosing)",
    disease_state = "Not applicable (mechanistic parameterization; Ayyar 2021 calibrated it on published fitusiran data)",
    dose_range = "Validation simulations: Ayyar 2021 Fig. 4 as regenerated in the Sten 2023 ESM (1, 3, 10 and 30 mg/kg SC in a 6.2 kg monkey)",
    notes = paste(
      "Parameter set of the cynomolgus-monkey branch of the ESM MATLAB function. The baseline plasma antithrombin is 110% (protein0 = 110), assigned by Sten 2023 from Ayyar 2021 Fig. 4.",
      "Sten 2023 made six corrections to the published Ayyar 2021 equations so that they",
      "regenerate Ayyar's figures: mRNA0 = protein0 = 100%; V_k,BC,vas = V_l,BC,vas = 0 in the",
      "tissue-concentration outputs; k_esc = k_deg,D * f_esc; k_syn = R_tot * k_deg;",
      "a corrected liver-interstitial equation with V_l,isf in mL; and k_deg,D = 0.012 1/h in all",
      "species. No parameter was re-estimated by Sten 2023; every value is carried from the",
      "corrected model, so all are fixed.",
      sep = " "
    )
  )

  ini({
    # ASGPR receptor kinetics - Sten 2023 ESM, MATLAB function Ayyar_mPBPKPD, branch species == 3 (cyno)
    rtot <- fixed(339); label("Baseline ASGPR density in liver Rtot (nM)") # ESM code: Rtot = 339
    kdeg <- fixed(0.04); label("ASGPR degradation rate constant kdeg (1/h)") # ESM code: kdeg = 0.04 (correction 4: ksyn = Rtot*kdeg)
    koff <- fixed(1.32); label("GalNAc-ASGPR dissociation rate constant koff (1/h)") # ESM code: koff = 1.32
    kon <- fixed(0.53); label("GalNAc-ASGPR association rate constant kon (1/(nM*h))") # ESM code: kon = 0.53
    kint_ref <- fixed(2.4); label("Internalization rate constant of bound ASGPR in the 28 g reference mouse kint (1/h)") # ESM code: kint = 2.4 (mouse), scaled as 2.4/(28/BW_g)^-0.25 in the other species
    krec <- fixed(13.8); label("Endosomal ASGPR recycling rate constant krec (1/h)") # ESM code: krec = 13.8
    kdeg_rendo <- fixed(1.53); label("Endosomal ASGPR degradation rate constant kdegR (1/h)") # ESM code: kdegR = 1.53
    kcle <- fixed(1.32); label("GalNAc cleavage rate constant from the siRNA in the endosome kcle (1/h)") # ESM code: kcle = 1.32

    # Absorption and plasma binding
    lka <- fixed(log(0.35)); label("Subcutaneous absorption rate constant ka (1/h)") # ESM code: ka = 0.35
    lfdepot <- fixed(log(33 / 100)); label("Apparent subcutaneous bioavailability F (fraction)") # ESM code: F = 33 (percent)
    fu <- fixed(0.15); label("Plasma fraction unbound fu (fraction)") # ESM code: fu = 0.15

    # Kidney
    ps_ref <- fixed(49.6); label("Kidney permeability-surface flow in the 28 g reference mouse PS (mL/h)") # ESM code: PS = 49.6 (mouse), scaled as 49.6/(28/BW_g)^0.75
    fuk <- fixed(0.07); label("Kidney tissue fraction unbound fuk (fraction)") # ESM code: fuk = 0.07
    kassk <- fixed(9.9); label("Kidney binding association rate constant kassk (1/h)") # ESM code: kassk = 9.9
    kdisk <- fixed(0.026); label("Kidney binding dissociation rate constant kdisk (1/h)") # ESM code: kdisk = 0.026

    # Liver uptake, binding and intracellular processing
    clup_in_ref <- fixed(207.2); label("Liver uptake clearance into the interstitium in the 28 g reference mouse CLuplivin (mL/h)") # ESM code: CLuplivin = 207.2 (mouse), scaled as 207.2/(28/BW_g)^0.75
    clup_eff_ref <- fixed(0.0025); label("Liver release clearance from the interstitium in the 28 g reference mouse CLupliveff (mL/h)") # ESM code: CLupliveff = 0.0025 (mouse), scaled as 0.0025/(28/BW_g)^0.75
    kassl <- fixed(0.0032); label("Liver binding association rate constant kassl (1/h)") # ESM code: kassl = 0.0032
    kdisl <- fixed(0.0034); label("Liver binding dissociation rate constant kdisl (1/h)") # ESM code: kdisl = 0.0034
    kpremainder <- fixed(0.1); label("Partition coefficient of the remainder tissue Kprest (unitless)") # ESM code: Kprest = 0.1
    fesc <- fixed(0.01); label("Fraction of endosomal siRNA escaping into the cytosol fesc (fraction)") # ESM code: fesc = 0.01 (correction 3: kesc = kdegD*fesc)
    kdeg_endo <- fixed(0.012); label("Endosomal siRNA degradation rate constant kdegD (1/h)") # ESM code: kdegD = 0.012 (correction 6: same value in all species)
    kon_risc <- fixed(0.00023); label("siRNA-RISC association rate constant konRISC (1/(nM*h))") # ESM code: konRISC = 0.00023
    koff_risc <- fixed(1e-7); label("siRNA-RISC dissociation rate constant koffRISC (1/h)") # ESM code: koffRISC = 1E-7
    risctot <- fixed(30); label("Total RISC concentration in hepatocytes RISCtot (nM)") # ESM code: RISCtot = 30
    kdeg_risc <- fixed(0.005); label("Degradation rate constant of loaded RISC kDR (1/h)") # ESM code: kDR = 0.005
    kdegc <- fixed(0.1); label("Degradation rate constant of free cytosolic siRNA kdegc (1/h)") # ESM code: kdegc = 0.1

    # Target (antithrombin) mRNA and protein
    kdeg_mrna <- fixed(0.06); label("Degradation rate constant of target mRNA kdegmRNA (1/h)") # ESM code: kdegmRNA = 0.06
    kdeg_protein <- fixed(0.01); label("Degradation rate constant of target protein kdegprotein (1/h)") # ESM code: kdegprotein = 0.01
    smax <- fixed(14.2); label("Maximal stimulation of target mRNA degradation Smax (unitless)") # ESM code: Smax = 14.2
    sc50 <- fixed(5.3); label("RISC-bound siRNA concentration at half-maximal stimulation SC50 (nM)") # ESM code: SC50 = 5.3; nM
    gamma_prot <- fixed(1.5); label("Power of relative mRNA on protein synthesis gamma (unitless)") # ESM code: gamma = 1.5
    mrna0 <- fixed(100); label("Baseline target mRNA (percent)") # ESM code: mRNA0 = 100 (correction 1)
    protein0 <- fixed(110); label("Baseline target protein (percent)") # ESM code: protein0 = 110

    # Physiology (monkey)
    qremainder <- fixed(10358); label("Plasma flow to the remainder tissue Qrest (mL/h)") # ESM code: Qrest = 10358
    vremainder <- fixed(1087); label("Extracellular volume of the remainder tissue Vrest (mL)") # ESM code: Vrest = 1087
    vplasma <- fixed(187); label("Plasma volume in circulation Vp (mL)") # ESM code: Vp = 187
    qkid <- fixed(3237); label("Plasma flow to the kidney Qkidney (mL/h)") # ESM code: Qkidney = 3237
    vk <- fixed(24.4); label("Interstitial and cellular volume of the kidney Vk (mL)") # ESM code: Vk = 24.4
    vkvas <- fixed(1.5); label("Plasma volume in the kidney vasculature Vkvas (mL)") # ESM code: Vkvas = 1.5
    qliver <- fixed(8838); label("Plasma flow to the liver Qliver (mL/h)") # ESM code: Qliver = 8838
    vlc <- fixed(120); label("Cellular volume of the liver Vlc (mL)") # ESM code: Vlc = 120
    vlisf <- fixed(37.4); label("Interstitial volume of the liver Vlisf (mL)") # ESM code: Vlisf = 37.4
    vlvas <- fixed(15.9); label("Plasma volume in the liver vasculature Vlvas (mL)") # ESM code: Vlvas = 15.9
    gfr <- fixed(1138); label("Glomerular filtration rate GFR (mL/h)") # ESM code: GFR = 1138
    bw <- fixed(6.2); label("Body weight of the simulated animal or human (kg)") # ESM figure scripts: BW = 6.2 kg
    bw_ref <- fixed(0.028); label("Body weight of the reference mouse for allometric scaling (kg)") # ESM code: the 28 (g) in the scaling factor (28/BW_g)
    e_bw_cl <- fixed(0.75); label("Allometric exponent on PS, CLuplivin and CLupliveff (unitless)") # ESM code: ^0.75
    e_bw_kint <- fixed(-0.25); label("Allometric exponent on kint (unitless)") # ESM code: ^-0.25
    mw <- fixed(15000); label("Molecular weight of the siRNA duplex MW (g/mol)") # ESM code: MW = 15000
  })

  model({
    # Allometric scaling from the 28 g reference mouse (ESM code, all non-mouse branches)
    kint <- kint_ref * (bw / bw_ref)^e_bw_kint
    ps <- ps_ref * (bw / bw_ref)^e_bw_cl
    clup_in <- clup_in_ref * (bw / bw_ref)^e_bw_cl
    clup_eff <- clup_eff_ref * (bw / bw_ref)^e_bw_cl

    ka <- exp(lka)
    fdepot <- exp(lfdepot)
    ksyn <- rtot * kdeg # correction 4
    kesc <- kdeg_endo * fesc # correction 3
    ksyn_mrna <- mrna0 * kdeg_mrna
    ksyn_protein <- protein0 * kdeg_protein
    fvol <- vlc * 0.8 / vlisf # ESM code: Fvol = Vlc*0.8/Vlisf

    # The depot holds mg; the ESM code doses in ng (Aa = dose * 1e6), so the input
    # to plasma converts mg to ng before dividing by Vp in mL.
    d/dt(depot) <- -ka * depot
    # Plasma (ng/mL)
    d/dt(central) <- fdepot * ka * depot * 1e6 / vplasma -
      (qremainder + qliver + qkid) * central / vplasma +
      qremainder * remainder / (kpremainder * vplasma) +
      qliver * liver_vas / vplasma + qkid * kidney_vas / vplasma
    # Remainder tissue (ng/mL)
    d/dt(remainder) <- (qremainder * central - qremainder * remainder / kpremainder) / vremainder
    # Kidney vascular (ng/mL)
    d/dt(kidney_vas) <- (qkid * (central - kidney_vas) - ps * (fu * kidney_vas - fuk * kidney) -
      gfr * fu * kidney_vas) / vkvas
    # Kidney tissue (ng/g)
    d/dt(kidney) <- ps * (fu * kidney_vas - fuk * kidney) / vk - kassk * fuk * kidney +
      kdisk * kidney_deep / vk
    # Kidney bound pool (ng)
    d/dt(kidney_deep) <- kassk * fuk * kidney * vk - kdisk * kidney_deep
    # Liver vascular (ng/mL)
    d/dt(liver_vas) <- (qliver * (central - liver_vas) - clup_in * fu * liver_vas +
      clup_eff * is_liver * mw / 1000) / vlvas
    # Liver interstitial fluid (nM); correction 5 (Vlisf in mL)
    d/dt(is_liver) <- (clup_in * fu * liver_vas / (mw / 1000) - clup_eff * is_liver) / vlisf -
      kon * is_liver * target + koff * complex - kassl * is_liver +
      kdisl * liver_deep / (vlisf / 1000)
    # Liver bound pool (nmol); correction 5
    d/dt(liver_deep) <- kassl * is_liver * (vlisf / 1000) - kdisl * liver_deep
    # Free cell-surface ASGPR (nM)
    d/dt(target) <- ksyn - kdeg * target - kon * is_liver * target + koff * complex +
      krec * target_endo * fvol
    # siRNA-bound cell-surface ASGPR (nM)
    d/dt(complex) <- kon * is_liver * target - koff * complex - kint * complex
    # siRNA-bound endosomal ASGPR (nM)
    d/dt(complex_endo) <- kint * complex / fvol - kcle * complex_endo
    # Free endosomal ASGPR (nM)
    d/dt(target_endo) <- kcle * complex_endo - krec * target_endo - kdeg_rendo * target_endo
    # Free endosomal siRNA (nM)
    d/dt(liver_endo) <- kcle * complex_endo - kesc * liver_endo - kdeg_endo * liver_endo
    # Free cytoplasmic siRNA (nM)
    d/dt(cyto) <- fesc * kdeg_endo * liver_endo - kdegc * cyto -
      kon_risc * cyto * (risctot - risc) + koff_risc * risc
    # RISC-loaded siRNA (nM)
    d/dt(risc) <- kon_risc * cyto * (risctot - risc) - koff_risc * risc - kdeg_risc * risc
    # Target mRNA in liver (percent of baseline)
    d/dt(mrna_liver) <- ksyn_mrna - kdeg_mrna * (1 + smax * risc / (sc50 + risc)) * mrna_liver
    # Target protein in plasma (percent of baseline)
    d/dt(protein) <- ksyn_protein * (mrna_liver / mrna0)^gamma_prot - kdeg_protein * protein

    target(0) <- rtot
    mrna_liver(0) <- mrna0
    protein(0) <- protein0

    # Outputs of the ESM figure scripts. Correction 2 sets the blood-cell vascular
    # volumes V_k,BC,vas and V_l,BC,vas to 0 in the tissue-concentration denominators.
    Cc <- central
    Ckidney <- (kidney * vk + kidney_deep) / (vkvas + vk)
    Cliver <- ((is_liver + complex) * vlisf / 1000 +
      (complex_endo + liver_endo + cyto) * vlc * 0.8 / 1000 + liver_deep) * mw /
      (vlisf + vlc + vlvas)
    # ESM: 6.65 kDa for the RISC-bound antisense strand converts nM to ng/g
    Crisc <- risc * 6.65
  })
}
