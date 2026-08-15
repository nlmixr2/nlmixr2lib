Wen_2025_salbutamol_rat_pbpk <- function() {
  description <- paste0(
    "PBPK (mechanistic pulmonary, 24-generation lung; rat). Translational ",
    "pulmonary physiologically based pharmacokinetic model for salbutamol in ",
    "male Wistar Han rats after a single intratracheal instillation. The lung ",
    "is resolved into 24 airway generations (1-16 tracheobronchial, 17-24 ",
    "alveolar), each split into epithelial lining fluid, epithelium and ",
    "subepithelium, giving 72 lung states coupled to a three-compartment ",
    "systemic disposition model. Per-generation airway geometry (number, ",
    "length, diameter, ELF and epithelial thickness) comes from the Yeh rat ",
    "lung cast; passive diffusion between the three tissue layers is driven by ",
    "a single effective permeability and the unbound lung:plasma partition ",
    "coefficient, and each subepithelium exchanges with plasma at a local blood ",
    "flow apportioned from bronchial flow (tracheobronchial) or cardiac output ",
    "(alveolar). Effective permeability and the unbound partition coefficient ",
    "were the two parameters estimated from rat plasma and bronchoalveolar-lavage ",
    "ELF profiles; every other value is a fixed physiological or literature input. ",
    "Instilled dose is apportioned across generations by the geometric deposition ",
    "rule of Appendix S1 Eq S13. Dissolution and mucociliary clearance are inactive ",
    "for this arm because the instillate is a solution of a freely water-soluble ",
    "drug (14.1 g/L) with no particle size reported (Wen 2025)."
  )
  reference <- paste(
    "Wen H, Sadiq MW, Friberg LE, Svensson EM.",
    "Translational physiologically based pharmacokinetic modeling to predict",
    "human pulmonary kinetics after lung delivery.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(4):796-806.",
    "doi:10.1002/psp4.13316",
    sep = " "
  )
  vignette <- "Wen_2025_pulmonary_pbpk"
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

  # The 72 airway-tissue states are genuinely paper-mechanistic and do not map
  # onto any canonical compartment role.
  paper_specific_compartment_pattern <- c("^elf_", "^epi_", "^sub_")

  # Issue #482: every ODE state integrates a drug AMOUNT in nmol. The three
  # sub-compartments of each airway generation are the epithelial lining fluid
  # (elf_NN), the epithelium (epi_NN) and the subepithelium (sub_NN) of Wen 2025
  # Figure 1b. Generations 01-16 are tracheobronchial, 17-24 alveolar.
  compartmentData <- list(
    elf_01 = list(analyte = "salbutamol", units = "nmol", specimen = "epithelial lining fluid", verified = TRUE), # generation 1 (tracheobronchial)
    epi_01 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    sub_01 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    elf_02 = list(analyte = "salbutamol", units = "nmol", specimen = "epithelial lining fluid", verified = TRUE), # generation 2 (tracheobronchial)
    epi_02 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    sub_02 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    elf_03 = list(analyte = "salbutamol", units = "nmol", specimen = "epithelial lining fluid", verified = TRUE), # generation 3 (tracheobronchial)
    epi_03 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    sub_03 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    elf_04 = list(analyte = "salbutamol", units = "nmol", specimen = "epithelial lining fluid", verified = TRUE), # generation 4 (tracheobronchial)
    epi_04 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    sub_04 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    elf_05 = list(analyte = "salbutamol", units = "nmol", specimen = "epithelial lining fluid", verified = TRUE), # generation 5 (tracheobronchial)
    epi_05 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    sub_05 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    elf_06 = list(analyte = "salbutamol", units = "nmol", specimen = "epithelial lining fluid", verified = TRUE), # generation 6 (tracheobronchial)
    epi_06 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    sub_06 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    elf_07 = list(analyte = "salbutamol", units = "nmol", specimen = "epithelial lining fluid", verified = TRUE), # generation 7 (tracheobronchial)
    epi_07 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    sub_07 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    elf_08 = list(analyte = "salbutamol", units = "nmol", specimen = "epithelial lining fluid", verified = TRUE), # generation 8 (tracheobronchial)
    epi_08 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    sub_08 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    elf_09 = list(analyte = "salbutamol", units = "nmol", specimen = "epithelial lining fluid", verified = TRUE), # generation 9 (tracheobronchial)
    epi_09 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    sub_09 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    elf_10 = list(analyte = "salbutamol", units = "nmol", specimen = "epithelial lining fluid", verified = TRUE), # generation 10 (tracheobronchial)
    epi_10 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    sub_10 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    elf_11 = list(analyte = "salbutamol", units = "nmol", specimen = "epithelial lining fluid", verified = TRUE), # generation 11 (tracheobronchial)
    epi_11 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    sub_11 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    elf_12 = list(analyte = "salbutamol", units = "nmol", specimen = "epithelial lining fluid", verified = TRUE), # generation 12 (tracheobronchial)
    epi_12 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    sub_12 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    elf_13 = list(analyte = "salbutamol", units = "nmol", specimen = "epithelial lining fluid", verified = TRUE), # generation 13 (tracheobronchial)
    epi_13 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    sub_13 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    elf_14 = list(analyte = "salbutamol", units = "nmol", specimen = "epithelial lining fluid", verified = TRUE), # generation 14 (tracheobronchial)
    epi_14 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    sub_14 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    elf_15 = list(analyte = "salbutamol", units = "nmol", specimen = "epithelial lining fluid", verified = TRUE), # generation 15 (tracheobronchial)
    epi_15 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    sub_15 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    elf_16 = list(analyte = "salbutamol", units = "nmol", specimen = "epithelial lining fluid", verified = TRUE), # generation 16 (tracheobronchial)
    epi_16 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    sub_16 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    elf_17 = list(analyte = "salbutamol", units = "nmol", specimen = "epithelial lining fluid", verified = TRUE), # generation 17 (alveolar)
    epi_17 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    sub_17 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    elf_18 = list(analyte = "salbutamol", units = "nmol", specimen = "epithelial lining fluid", verified = TRUE), # generation 18 (alveolar)
    epi_18 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    sub_18 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    elf_19 = list(analyte = "salbutamol", units = "nmol", specimen = "epithelial lining fluid", verified = TRUE), # generation 19 (alveolar)
    epi_19 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    sub_19 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    elf_20 = list(analyte = "salbutamol", units = "nmol", specimen = "epithelial lining fluid", verified = TRUE), # generation 20 (alveolar)
    epi_20 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    sub_20 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    elf_21 = list(analyte = "salbutamol", units = "nmol", specimen = "epithelial lining fluid", verified = TRUE), # generation 21 (alveolar)
    epi_21 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    sub_21 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    elf_22 = list(analyte = "salbutamol", units = "nmol", specimen = "epithelial lining fluid", verified = TRUE), # generation 22 (alveolar)
    epi_22 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    sub_22 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    elf_23 = list(analyte = "salbutamol", units = "nmol", specimen = "epithelial lining fluid", verified = TRUE), # generation 23 (alveolar)
    epi_23 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    sub_23 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    elf_24 = list(analyte = "salbutamol", units = "nmol", specimen = "epithelial lining fluid", verified = TRUE), # generation 24 (alveolar)
    epi_24 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    sub_24 = list(analyte = "salbutamol", units = "nmol", specimen = "tissue", verified = TRUE),
    central     = list(analyte = "salbutamol", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "salbutamol", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "salbutamol", units = "nmol", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Sets total lung volume (V_lung = 0.004126 * WT L; Appendix S1 section ",
        "3.1) and cardiac output (Q_CO = 20.77 * WT L/h; Appendix S1 section ",
        "3.1), which together fix the subepithelial volumes and every local ",
        "blood flow. Airway number, length and diameter are NOT scaled with ",
        "body weight - they are the fixed Yeh lung-cast morphometry of a 330 g ",
        "rat (Appendix S1 Table S2). The salbutamol arm used rats of 350 g ",
        "(Wen 2025 section 2.4), i.e. WT = 0.35."
      ),
      source_name        = "BW"
    )
  )

  population <- list(
    species        = "rat (male Wistar Han)",
    n_subjects     = NA_integer_,
    n_studies      = 1L,
    weight_median  = "350 g",
    disease_state  = "Healthy",
    dose_range     = "100 nmol/kg single dose by intratracheal instillation",
    notes          = paste0(
      "Plasma and bronchoalveolar-lavage-obtained ELF concentration-time data ",
      "were digitised from Boger and Friden (2019) J Aerosol Med Pulm Drug Deliv ",
      "32:1-12 (Wen 2025 Appendix S1 Table S1); the number of animals is not ",
      "restated in Wen 2025. The model is deterministic - no inter-individual ",
      "random effects and no residual-error model were reported, because the two ",
      "estimated parameters were fitted to mean profiles by Monte Carlo search ",
      "followed by Levenberg-Marquardt in MoBi (Wen 2025 section 2.4)."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # The two parameters estimated from the rat data (Wen 2025 section 3.2).
    # Both are FIXED here at their published point estimates.
    # ------------------------------------------------------------------
    lpeff    <- fixed(log(1.18e-5)); label("Effective permeability, all generations (cm/s)")          # Wen 2025 section 3.2 (95% CI 9.93e-6 to 1.36e-5)
    lkpulung <- fixed(log(8.83));    label("Unbound lung:plasma partition coefficient (unitless)")    # Wen 2025 section 3.2 (95% CI 6.96 to 10.68)

    # ------------------------------------------------------------------
    # Fixed drug parameters (Wen 2025 Table 1).
    # ------------------------------------------------------------------
    lkplung <- fixed(log(4.66)); label("Total lung:plasma partition coefficient (unitless)")          # Wen 2025 Table 1, from Boger and Friden
    lfuelf  <- fixed(log(0.80)); label("Fraction unbound in epithelial lining fluid (unitless)")      # Wen 2025 Table 1; Appendix S1 Eq S21 from fu,p = 0.71
    lbp     <- fixed(log(0.96)); label("Blood-to-plasma ratio (unitless)")                            # Wen 2025 Table 1, from Boger and Friden

    # ------------------------------------------------------------------
    # Fixed systemic disposition (Wen 2025 Table 1, three-compartment model
    # for salbutamol in the rat; peripheral1 is the slowly and peripheral2
    # the rapidly equilibrating compartment).
    # ------------------------------------------------------------------
    lcl  <- fixed(log(2.31)); label("Systemic clearance (L/h)")                                       # Wen 2025 Table 1
    lvc  <- fixed(log(0.26)); label("Central volume of distribution (L)")                             # Wen 2025 Table 1 (Vc)
    lvp  <- fixed(log(1.03)); label("Slowly equilibrating peripheral volume (L)")                     # Wen 2025 Table 1 (V2)
    lq   <- fixed(log(0.69)); label("Slow distribution clearance (L/h)")                              # Wen 2025 Table 1 (Q2)
    lvp2 <- fixed(log(0.51)); label("Rapidly equilibrating peripheral volume (L)")                    # Wen 2025 Table 1 (V3)
    lq2  <- fixed(log(2.24)); label("Fast distribution clearance (L/h)")                              # Wen 2025 Table 1 (Q3)

    # ------------------------------------------------------------------
    # Fixed rat physiology (Appendix S1 section 3.1).
    # ------------------------------------------------------------------
    lqcopkg    <- fixed(log(20.77));    label("Cardiac output per kg body weight (L/h/kg)")           # Appendix S1 section 3.1
    fqbr       <- fixed(0.02);          label("Bronchial fraction of cardiac output (unitless)")      # Appendix S1 section 3.1
    lvlungpkg  <- fixed(log(0.004126)); label("Total lung volume per kg body weight (L/kg)")          # Appendix S1 section 3.1

    # ------------------------------------------------------------------
    # Residual error. Wen 2025 reports no residual-error model - the two
    # estimated parameters were fitted to mean observed profiles - so these
    # are FIXED placeholders that exist only so the model has a declared
    # endpoint for each observed matrix.
    # ------------------------------------------------------------------
    propSd     <- fixed(0.10); label("Proportional residual error placeholder, plasma (fraction)")    # not reported in Wen 2025
    propSd_Celf <- fixed(0.10); label("Proportional residual error placeholder, mean ELF (fraction)")  # not reported in Wen 2025
  })

  model({
    # ==================================================================
    # 1. Individual parameters
    # ==================================================================
    peff    <- exp(lpeff)      # cm/s
    kpulung <- exp(lkpulung)
    kplung  <- exp(lkplung)
    fuelf   <- exp(lfuelf)
    bp      <- exp(lbp)

    cl  <- exp(lcl)
    vc  <- exp(lvc)
    vp  <- exp(lvp)
    q   <- exp(lq)
    vp2 <- exp(lvp2)
    q2  <- exp(lq2)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # Permeability in solver units: cm/s * 3600 s/h / 1000 cm3/L gives the
    # L/h delivered per cm2 of surface per unit concentration (nmol/L).
    pk <- peff * 3.6

    # Rat physiology scaled to body weight (Appendix S1 section 3.1).
    qco   <- exp(lqcopkg) * WT      # L/h, cardiac output, perfuses the alveolar region
    qbr   <- fqbr * qco             # L/h, bronchial flow, perfuses the tracheobronchial region
    vlung <- exp(lvlungpkg) * WT    # L, total lung volume

    # Subepithelial volume is whatever total lung volume is left after the
    # ELF and epithelium are accounted for (Appendix S1 Eq S8), apportioned
    # across generations in proportion to epithelial volume (Eq S9). The two
    # subtracted constants are sum(velf) and sum(vepi) over generations 1-24.
    vrem <- vlung - 3.697271e-05 - 1.226157e-04

    # Plasma concentration drives every subepithelium-to-blood exchange, so it
    # is formed before the fluxes.
    Cc <- central / vc

    # ==================================================================
    # 2. Airway geometry
    # ==================================================================
    # -- Airway geometry (Appendix S1 Table S2 + Eqs S1-S9). Every literal
    #    below is reproduced from the Table S2 row for that generation by the
    #    equation named in its comment; the validation vignette recomputes all
    #    of them from Table S2 and asserts equality with these values.
    #    velf<i>  = N*L*pi*h_elf*(D - h_elf)          (Eqs S5, S6)  [L]
    #    vepi<i>  = N*L*pi*h_epi*(D + h_epi)          (Eqs S5, S7)  [L]
    #    rsub<i>  = vepi<i> / sum(vepi)               (Eq S9)       [-]
    #    aelf<i>  = pi*D*L*N   + fi*(3870 - sum)      (Eqs S1,S3,S4)[cm2]
    #    aepi<i>  = pi*L*N*(D+2*h_epi) + fi*(3870-s)  (Eqs S2,S3,S4)[cm2]
    #    wq<i>    = blood-flow weight                 (Eqs S10-S12) [-]
    velf01 <- 1.429205e-06; vepi01 <- 6.918782e-06; rsub01 <- 5.642656e-02  # gen 01: N=1, L=2.68 cm, D=0.34 cm, h_elf=5 um, h_epi=24 um
    aelf01 <- 2.862619e+00; aepi01 <- 2.903033e+00; wq01 <- 6.372900e-02  # TB region
    velf02 <- 6.502861e-07; vepi02 <- 1.701256e-06; rsub02 <- 1.387470e-02  # gen 02: N=2, L=0.715 cm, D=0.29 cm, h_elf=5 um, h_epi=13 um
    aelf02 <- 1.302818e+00; aepi02 <- 1.314499e+00; wq02 <- 1.895834e-02  # TB region
    velf03 <- 4.948008e-07; vepi03 <- 1.295304e-06; rsub03 <- 1.056393e-02  # gen 03: N=3, L=0.4 cm, D=0.263 cm, h_elf=5 um, h_epi=13 um
    aelf03 <- 9.914866e-01; aepi03 <- 1.001288e+00; wq03 <- 1.607652e-02  # TB region
    velf04 <- 2.799159e-07; vepi04 <- 7.342505e-07; rsub04 <- 5.988226e-03  # gen 04: N=5, L=0.176 cm, D=0.203 cm, h_elf=5 um, h_epi=13 um
    aelf04 <- 5.612141e-01; aepi04 <- 5.684021e-01; wq04 <- 1.170327e-02  # TB region
    velf05 <- 4.247433e-07; vepi05 <- 1.116565e-06; rsub05 <- 9.106218e-03  # gen 05: N=8, L=0.208 cm, D=0.163 cm, h_elf=5 um, h_epi=13 um
    aelf05 <- 8.521005e-01; aepi05 <- 8.656922e-01; wq05 <- 2.117833e-02  # TB region
    velf06 <- 3.434907e-07; vepi06 <- 9.051174e-07; rsub06 <- 7.381742e-03  # gen 06: N=14, L=0.117 cm, D=0.134 cm, h_elf=5 um, h_epi=13 um
    aelf06 <- 6.895545e-01; aepi06 <- 7.029339e-01; wq06 <- 1.953528e-02  # TB region
    velf07 <- 5.045319e-07; vepi07 <- 1.331058e-06; rsub07 <- 1.085553e-02  # gen 07: N=23, L=0.114 cm, D=0.123 cm, h_elf=5 um, h_epi=13 um
    aelf07 <- 1.013182e+00; aepi07 <- 1.034599e+00; wq07 <- 3.018975e-02  # TB region
    velf08 <- 6.375234e-07; vepi08 <- 1.684320e-06; rsub08 <- 1.373657e-02  # gen 08: N=28, L=0.13 cm, D=0.112 cm, h_elf=5 um, h_epi=13 um
    aelf08 <- 1.280764e+00; aepi08 <- 1.310497e+00; wq08 <- 4.015795e-02  # TB region
    velf09 <- 9.552130e-07; vepi09 <- 2.530860e-06; rsub09 <- 2.064058e-02  # gen 09: N=65, L=0.099 cm, D=0.095 cm, h_elf=5 um, h_epi=13 um
    aelf09 <- 1.920534e+00; aepi09 <- 1.973096e+00; wq09 <- 6.522031e-02  # TB region
    velf10 <- 1.347733e-06; vepi10 <- 3.577024e-06; rsub10 <- 2.917264e-02  # gen 10: N=109, L=0.091 cm, D=0.087 cm, h_elf=5 um, h_epi=13 um
    aelf10 <- 2.711047e+00; aepi10 <- 2.792067e+00; wq10 <- 9.563780e-02  # TB region
    velf11 <- 2.150357e-06; vepi11 <- 5.720783e-06; rsub11 <- 4.665620e-02  # gen 11: N=184, L=0.096 cm, D=0.078 cm, h_elf=5 um, h_epi=13 um
    aelf11 <- 4.328461e+00; aepi11 <- 4.472743e+00; wq11 <- 1.594520e-01  # TB region
    velf12 <- 2.462555e-06; vepi12 <- 6.568468e-06; rsub12 <- 5.356955e-02  # gen 12: N=309, L=0.073 cm, D=0.07 cm, h_elf=5 um, h_epi=13 um
    aelf12 <- 4.960543e+00; aepi12 <- 5.144792e+00; wq12 <- 1.900037e-01  # TB region
    velf13 <- 3.529285e-06; vepi13 <- 9.463394e-06; rsub13 <- 7.717930e-02  # gen 13: N=521, L=0.075 cm, D=0.058 cm, h_elf=5 um, h_epi=13 um
    aelf13 <- 7.119949e+00; aepi13 <- 7.439119e+00; wq13 <- 2.894933e-01  # TB region
    velf14 <- 4.008782e-06; vepi14 <- 1.080966e-05; rsub14 <- 8.815886e-02  # gen 14: N=877, L=0.06 cm, D=0.049 cm, h_elf=5 um, h_epi=13 um
    aelf14 <- 8.100220e+00; aepi14 <- 8.530027e+00; wq14 <- 3.449103e-01  # TB region
    velf15 <- 4.529929e-06; vepi15 <- 1.237500e-05; rsub15 <- 1.009251e-01  # gen 15: N=1477, L=0.055 cm, D=0.036 cm, h_elf=5 um, h_epi=13 um
    aelf15 <- 9.187462e+00; aepi15 <- 9.851001e+00; wq15 <- 4.197537e-01  # TB region
    velf16 <- 2.666234e-06; vepi16 <- 7.572106e-06; rsub16 <- 6.175478e-02  # gen 16: N=2487, L=0.035 cm, D=0.02 cm, h_elf=5 um, h_epi=13 um
    aelf16 <- 5.469199e+00; aepi16 <- 6.180194e+00; wq16 <- 2.770329e-01  # TB region
    velf17 <- 3.738588e-06; vepi17 <- 1.078073e-05; rsub17 <- 8.792290e-02  # gen 17: N=4974, L=0.029 cm, D=0.017 cm, h_elf=5 um, h_epi=13 um
    aelf17 <- 7.703757e+00; aepi17 <- 8.881979e+00; wq17 <- 2.231492e-01  # alveolar region
    velf18 <- 8.746889e-08; vepi18 <- 4.811915e-07; rsub18 <- 3.924387e-03  # gen 18: N=9948, L=0.025 cm, D=0.016 cm, h_elf=0.07 um, h_epi=0.384 um
    aelf18 <- 1.827621e+01; aepi18 <- 1.832319e+01; wq18 <- 9.960133e-03  # alveolar region; alveolarisation f=0.002
    velf19 <- 1.443195e-07; vepi19 <- 7.940927e-07; rsub19 <- 6.476272e-03  # gen 19: N=19896, L=0.022 cm, D=0.015 cm, h_elf=0.07 um, h_epi=0.384 um
    aelf19 <- 4.083983e+01; aepi19 <- 4.089988e+01; wq19 <- 1.643684e-02  # alveolar region; alveolarisation f=0.007
    velf20 <- 2.448976e-07; vepi20 <- 1.347797e-06; rsub20 <- 1.099204e-02  # gen 20: N=39792, L=0.02 cm, D=0.014 cm, h_elf=0.07 um, h_epi=0.384 um
    aelf20 <- 9.275470e+01; aepi20 <- 9.281653e+01; wq20 <- 2.789791e-02  # alveolar region; alveolarisation f=0.02
    velf21 <- 4.653054e-07; vepi21 <- 2.560814e-06; rsub21 <- 2.088488e-02  # gen 21: N=79584, L=0.019 cm, D=0.014 cm, h_elf=0.07 um, h_epi=0.384 um
    aelf21 <- 2.686368e+02; aepi21 <- 2.685460e+02; wq21 <- 5.300603e-02  # alveolar region; alveolarisation f=0.07
    velf22 <- 8.816313e-07; vepi22 <- 4.852069e-06; rsub22 <- 3.957135e-02  # gen 22: N=159168, L=0.018 cm, D=0.014 cm, h_elf=0.07 um, h_epi=0.384 um
    aelf22 <- 5.273855e+02; aepi22 <- 5.271720e+02; wq22 <- 1.004325e-01  # alveolar region; alveolarisation f=0.139
    velf23 <- 1.665304e-06; vepi23 <- 9.165019e-06; rsub23 <- 7.474589e-02  # gen 23: N=318336, L=0.017 cm, D=0.014 cm, h_elf=0.07 um, h_epi=0.384 um
    aelf23 <- 1.052320e+03; aepi23 <- 1.051790e+03; wq23 <- 1.897058e-01  # alveolar region; alveolarisation f=0.282
    velf24 <- 3.330607e-06; vepi24 <- 1.833004e-05; rsub24 <- 1.494918e-01  # gen 24: N=636672, L=0.017 cm, D=0.014 cm, h_elf=0.07 um, h_epi=0.384 um
    aelf24 <- 1.862083e+03; aepi24 <- 1.861570e+03; wq24 <- 3.794116e-01  # alveolar region; alveolarisation f=0.48

    # ==================================================================
    # 3. Mass transfer (Wen 2025 Eqs 1-3; Appendix S1 Eqs S10-S12)
    #    jelf: ELF -> epithelium, jepi: epithelium -> subepithelium,
    #    jbld: subepithelium -> plasma.
    # ==================================================================
    jelf01 <- pk * aelf01 * (fuelf * elf_01 / velf01 - epi_01 / vepi01 / kpulung)
    jepi01 <- pk * aepi01 * (epi_01 / vepi01 / kpulung - sub_01 / (rsub01 * vrem) / kpulung)
    jbld01 <- qbr * wq01 * (bp * sub_01 / (rsub01 * vrem) / kplung - Cc)

    jelf02 <- pk * aelf02 * (fuelf * elf_02 / velf02 - epi_02 / vepi02 / kpulung)
    jepi02 <- pk * aepi02 * (epi_02 / vepi02 / kpulung - sub_02 / (rsub02 * vrem) / kpulung)
    jbld02 <- qbr * wq02 * (bp * sub_02 / (rsub02 * vrem) / kplung - Cc)

    jelf03 <- pk * aelf03 * (fuelf * elf_03 / velf03 - epi_03 / vepi03 / kpulung)
    jepi03 <- pk * aepi03 * (epi_03 / vepi03 / kpulung - sub_03 / (rsub03 * vrem) / kpulung)
    jbld03 <- qbr * wq03 * (bp * sub_03 / (rsub03 * vrem) / kplung - Cc)

    jelf04 <- pk * aelf04 * (fuelf * elf_04 / velf04 - epi_04 / vepi04 / kpulung)
    jepi04 <- pk * aepi04 * (epi_04 / vepi04 / kpulung - sub_04 / (rsub04 * vrem) / kpulung)
    jbld04 <- qbr * wq04 * (bp * sub_04 / (rsub04 * vrem) / kplung - Cc)

    jelf05 <- pk * aelf05 * (fuelf * elf_05 / velf05 - epi_05 / vepi05 / kpulung)
    jepi05 <- pk * aepi05 * (epi_05 / vepi05 / kpulung - sub_05 / (rsub05 * vrem) / kpulung)
    jbld05 <- qbr * wq05 * (bp * sub_05 / (rsub05 * vrem) / kplung - Cc)

    jelf06 <- pk * aelf06 * (fuelf * elf_06 / velf06 - epi_06 / vepi06 / kpulung)
    jepi06 <- pk * aepi06 * (epi_06 / vepi06 / kpulung - sub_06 / (rsub06 * vrem) / kpulung)
    jbld06 <- qbr * wq06 * (bp * sub_06 / (rsub06 * vrem) / kplung - Cc)

    jelf07 <- pk * aelf07 * (fuelf * elf_07 / velf07 - epi_07 / vepi07 / kpulung)
    jepi07 <- pk * aepi07 * (epi_07 / vepi07 / kpulung - sub_07 / (rsub07 * vrem) / kpulung)
    jbld07 <- qbr * wq07 * (bp * sub_07 / (rsub07 * vrem) / kplung - Cc)

    jelf08 <- pk * aelf08 * (fuelf * elf_08 / velf08 - epi_08 / vepi08 / kpulung)
    jepi08 <- pk * aepi08 * (epi_08 / vepi08 / kpulung - sub_08 / (rsub08 * vrem) / kpulung)
    jbld08 <- qbr * wq08 * (bp * sub_08 / (rsub08 * vrem) / kplung - Cc)

    jelf09 <- pk * aelf09 * (fuelf * elf_09 / velf09 - epi_09 / vepi09 / kpulung)
    jepi09 <- pk * aepi09 * (epi_09 / vepi09 / kpulung - sub_09 / (rsub09 * vrem) / kpulung)
    jbld09 <- qbr * wq09 * (bp * sub_09 / (rsub09 * vrem) / kplung - Cc)

    jelf10 <- pk * aelf10 * (fuelf * elf_10 / velf10 - epi_10 / vepi10 / kpulung)
    jepi10 <- pk * aepi10 * (epi_10 / vepi10 / kpulung - sub_10 / (rsub10 * vrem) / kpulung)
    jbld10 <- qbr * wq10 * (bp * sub_10 / (rsub10 * vrem) / kplung - Cc)

    jelf11 <- pk * aelf11 * (fuelf * elf_11 / velf11 - epi_11 / vepi11 / kpulung)
    jepi11 <- pk * aepi11 * (epi_11 / vepi11 / kpulung - sub_11 / (rsub11 * vrem) / kpulung)
    jbld11 <- qbr * wq11 * (bp * sub_11 / (rsub11 * vrem) / kplung - Cc)

    jelf12 <- pk * aelf12 * (fuelf * elf_12 / velf12 - epi_12 / vepi12 / kpulung)
    jepi12 <- pk * aepi12 * (epi_12 / vepi12 / kpulung - sub_12 / (rsub12 * vrem) / kpulung)
    jbld12 <- qbr * wq12 * (bp * sub_12 / (rsub12 * vrem) / kplung - Cc)

    jelf13 <- pk * aelf13 * (fuelf * elf_13 / velf13 - epi_13 / vepi13 / kpulung)
    jepi13 <- pk * aepi13 * (epi_13 / vepi13 / kpulung - sub_13 / (rsub13 * vrem) / kpulung)
    jbld13 <- qbr * wq13 * (bp * sub_13 / (rsub13 * vrem) / kplung - Cc)

    jelf14 <- pk * aelf14 * (fuelf * elf_14 / velf14 - epi_14 / vepi14 / kpulung)
    jepi14 <- pk * aepi14 * (epi_14 / vepi14 / kpulung - sub_14 / (rsub14 * vrem) / kpulung)
    jbld14 <- qbr * wq14 * (bp * sub_14 / (rsub14 * vrem) / kplung - Cc)

    jelf15 <- pk * aelf15 * (fuelf * elf_15 / velf15 - epi_15 / vepi15 / kpulung)
    jepi15 <- pk * aepi15 * (epi_15 / vepi15 / kpulung - sub_15 / (rsub15 * vrem) / kpulung)
    jbld15 <- qbr * wq15 * (bp * sub_15 / (rsub15 * vrem) / kplung - Cc)

    jelf16 <- pk * aelf16 * (fuelf * elf_16 / velf16 - epi_16 / vepi16 / kpulung)
    jepi16 <- pk * aepi16 * (epi_16 / vepi16 / kpulung - sub_16 / (rsub16 * vrem) / kpulung)
    jbld16 <- qbr * wq16 * (bp * sub_16 / (rsub16 * vrem) / kplung - Cc)

    jelf17 <- pk * aelf17 * (fuelf * elf_17 / velf17 - epi_17 / vepi17 / kpulung)
    jepi17 <- pk * aepi17 * (epi_17 / vepi17 / kpulung - sub_17 / (rsub17 * vrem) / kpulung)
    jbld17 <- qco * wq17 * (bp * sub_17 / (rsub17 * vrem) / kplung - Cc)

    jelf18 <- pk * aelf18 * (fuelf * elf_18 / velf18 - epi_18 / vepi18 / kpulung)
    jepi18 <- pk * aepi18 * (epi_18 / vepi18 / kpulung - sub_18 / (rsub18 * vrem) / kpulung)
    jbld18 <- qco * wq18 * (bp * sub_18 / (rsub18 * vrem) / kplung - Cc)

    jelf19 <- pk * aelf19 * (fuelf * elf_19 / velf19 - epi_19 / vepi19 / kpulung)
    jepi19 <- pk * aepi19 * (epi_19 / vepi19 / kpulung - sub_19 / (rsub19 * vrem) / kpulung)
    jbld19 <- qco * wq19 * (bp * sub_19 / (rsub19 * vrem) / kplung - Cc)

    jelf20 <- pk * aelf20 * (fuelf * elf_20 / velf20 - epi_20 / vepi20 / kpulung)
    jepi20 <- pk * aepi20 * (epi_20 / vepi20 / kpulung - sub_20 / (rsub20 * vrem) / kpulung)
    jbld20 <- qco * wq20 * (bp * sub_20 / (rsub20 * vrem) / kplung - Cc)

    jelf21 <- pk * aelf21 * (fuelf * elf_21 / velf21 - epi_21 / vepi21 / kpulung)
    jepi21 <- pk * aepi21 * (epi_21 / vepi21 / kpulung - sub_21 / (rsub21 * vrem) / kpulung)
    jbld21 <- qco * wq21 * (bp * sub_21 / (rsub21 * vrem) / kplung - Cc)

    jelf22 <- pk * aelf22 * (fuelf * elf_22 / velf22 - epi_22 / vepi22 / kpulung)
    jepi22 <- pk * aepi22 * (epi_22 / vepi22 / kpulung - sub_22 / (rsub22 * vrem) / kpulung)
    jbld22 <- qco * wq22 * (bp * sub_22 / (rsub22 * vrem) / kplung - Cc)

    jelf23 <- pk * aelf23 * (fuelf * elf_23 / velf23 - epi_23 / vepi23 / kpulung)
    jepi23 <- pk * aepi23 * (epi_23 / vepi23 / kpulung - sub_23 / (rsub23 * vrem) / kpulung)
    jbld23 <- qco * wq23 * (bp * sub_23 / (rsub23 * vrem) / kplung - Cc)

    jelf24 <- pk * aelf24 * (fuelf * elf_24 / velf24 - epi_24 / vepi24 / kpulung)
    jepi24 <- pk * aepi24 * (epi_24 / vepi24 / kpulung - sub_24 / (rsub24 * vrem) / kpulung)
    jbld24 <- qco * wq24 * (bp * sub_24 / (rsub24 * vrem) / kplung - Cc)

    # ==================================================================
    # 4. ODE system
    # ==================================================================
    d/dt(elf_01) <- -jelf01
    d/dt(epi_01) <-  jelf01 - jepi01
    d/dt(sub_01) <-  jepi01 - jbld01
    d/dt(elf_02) <- -jelf02
    d/dt(epi_02) <-  jelf02 - jepi02
    d/dt(sub_02) <-  jepi02 - jbld02
    d/dt(elf_03) <- -jelf03
    d/dt(epi_03) <-  jelf03 - jepi03
    d/dt(sub_03) <-  jepi03 - jbld03
    d/dt(elf_04) <- -jelf04
    d/dt(epi_04) <-  jelf04 - jepi04
    d/dt(sub_04) <-  jepi04 - jbld04
    d/dt(elf_05) <- -jelf05
    d/dt(epi_05) <-  jelf05 - jepi05
    d/dt(sub_05) <-  jepi05 - jbld05
    d/dt(elf_06) <- -jelf06
    d/dt(epi_06) <-  jelf06 - jepi06
    d/dt(sub_06) <-  jepi06 - jbld06
    d/dt(elf_07) <- -jelf07
    d/dt(epi_07) <-  jelf07 - jepi07
    d/dt(sub_07) <-  jepi07 - jbld07
    d/dt(elf_08) <- -jelf08
    d/dt(epi_08) <-  jelf08 - jepi08
    d/dt(sub_08) <-  jepi08 - jbld08
    d/dt(elf_09) <- -jelf09
    d/dt(epi_09) <-  jelf09 - jepi09
    d/dt(sub_09) <-  jepi09 - jbld09
    d/dt(elf_10) <- -jelf10
    d/dt(epi_10) <-  jelf10 - jepi10
    d/dt(sub_10) <-  jepi10 - jbld10
    d/dt(elf_11) <- -jelf11
    d/dt(epi_11) <-  jelf11 - jepi11
    d/dt(sub_11) <-  jepi11 - jbld11
    d/dt(elf_12) <- -jelf12
    d/dt(epi_12) <-  jelf12 - jepi12
    d/dt(sub_12) <-  jepi12 - jbld12
    d/dt(elf_13) <- -jelf13
    d/dt(epi_13) <-  jelf13 - jepi13
    d/dt(sub_13) <-  jepi13 - jbld13
    d/dt(elf_14) <- -jelf14
    d/dt(epi_14) <-  jelf14 - jepi14
    d/dt(sub_14) <-  jepi14 - jbld14
    d/dt(elf_15) <- -jelf15
    d/dt(epi_15) <-  jelf15 - jepi15
    d/dt(sub_15) <-  jepi15 - jbld15
    d/dt(elf_16) <- -jelf16
    d/dt(epi_16) <-  jelf16 - jepi16
    d/dt(sub_16) <-  jepi16 - jbld16
    d/dt(elf_17) <- -jelf17
    d/dt(epi_17) <-  jelf17 - jepi17
    d/dt(sub_17) <-  jepi17 - jbld17
    d/dt(elf_18) <- -jelf18
    d/dt(epi_18) <-  jelf18 - jepi18
    d/dt(sub_18) <-  jepi18 - jbld18
    d/dt(elf_19) <- -jelf19
    d/dt(epi_19) <-  jelf19 - jepi19
    d/dt(sub_19) <-  jepi19 - jbld19
    d/dt(elf_20) <- -jelf20
    d/dt(epi_20) <-  jelf20 - jepi20
    d/dt(sub_20) <-  jepi20 - jbld20
    d/dt(elf_21) <- -jelf21
    d/dt(epi_21) <-  jelf21 - jepi21
    d/dt(sub_21) <-  jepi21 - jbld21
    d/dt(elf_22) <- -jelf22
    d/dt(epi_22) <-  jelf22 - jepi22
    d/dt(sub_22) <-  jepi22 - jbld22
    d/dt(elf_23) <- -jelf23
    d/dt(epi_23) <-  jelf23 - jepi23
    d/dt(sub_23) <-  jepi23 - jbld23
    d/dt(elf_24) <- -jelf24
    d/dt(epi_24) <-  jelf24 - jepi24
    d/dt(sub_24) <-  jepi24 - jbld24

    # Systemic disposition. Every subepithelium drains into the central
    # compartment (Wen 2025 Figure 1); the three-compartment structure and its
    # rate constants are the MoBi passive transports k_e, k_12/k_21 and
    # k_13/k_31 of the published project file.
    jbldtot <- 
      jbld01 + jbld02 + jbld03 + jbld04 + jbld05 + jbld06 + jbld07 + jbld08 + jbld09 +
      jbld10 + jbld11 + jbld12 + jbld13 + jbld14 + jbld15 + jbld16 + jbld17 + jbld18 +
      jbld19 + jbld20 + jbld21 + jbld22 + jbld23 + jbld24

    d/dt(central)     <- jbldtot - kel * central -
                         k12 * central + k21 * peripheral1 -
                         k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # ==================================================================
    # 5. Deposition of the instilled dose
    # ==================================================================
    # -- Deposition after intratracheal instillation (Eq S13):
    #    df_i = df_1 * 0.985^(i-1), normalised over generations 1-24.
    #    Gives 70.6% tracheobronchial / 29.4% alveolar, reproducing the
    #    29.3% alveolar figure quoted in Appendix S1 section 3.2.
    f(elf_01) <- 4.930580e-02
    f(elf_02) <- 4.856621e-02
    f(elf_03) <- 4.783772e-02
    f(elf_04) <- 4.712015e-02
    f(elf_05) <- 4.641335e-02
    f(elf_06) <- 4.571715e-02
    f(elf_07) <- 4.503139e-02
    f(elf_08) <- 4.435592e-02
    f(elf_09) <- 4.369058e-02
    f(elf_10) <- 4.303522e-02
    f(elf_11) <- 4.238970e-02
    f(elf_12) <- 4.175385e-02
    f(elf_13) <- 4.112754e-02
    f(elf_14) <- 4.051063e-02
    f(elf_15) <- 3.990297e-02
    f(elf_16) <- 3.930442e-02
    f(elf_17) <- 3.871486e-02
    f(elf_18) <- 3.813414e-02
    f(elf_19) <- 3.756212e-02
    f(elf_20) <- 3.699869e-02
    f(elf_21) <- 3.644371e-02
    f(elf_22) <- 3.589706e-02
    f(elf_23) <- 3.535860e-02
    f(elf_24) <- 3.482822e-02

    # ==================================================================
    # 6. Observations
    #    Cc    - plasma concentration (nmol/L)
    #    Celf  - mean ELF concentration over all 24 generations, the quantity
    #            compared against bronchoalveolar-lavage data (Wen 2025 Fig 2)
    #    Clung - total lung concentration (all three layers / lung volume)
    #    Cbr   - mean ELF concentration over generations 4-9 ("bronchial ELF")
    # ==================================================================
    aelftot <- 
      elf_01 + elf_02 + elf_03 + elf_04 + elf_05 + elf_06 + elf_07 + elf_08 + elf_09 +
      elf_10 + elf_11 + elf_12 + elf_13 + elf_14 + elf_15 + elf_16 + elf_17 + elf_18 +
      elf_19 + elf_20 + elf_21 + elf_22 + elf_23 + elf_24
    aepitot <- 
      epi_01 + epi_02 + epi_03 + epi_04 + epi_05 + epi_06 + epi_07 + epi_08 + epi_09 +
      epi_10 + epi_11 + epi_12 + epi_13 + epi_14 + epi_15 + epi_16 + epi_17 + epi_18 +
      epi_19 + epi_20 + epi_21 + epi_22 + epi_23 + epi_24
    asubtot <- 
      sub_01 + sub_02 + sub_03 + sub_04 + sub_05 + sub_06 + sub_07 + sub_08 + sub_09 +
      sub_10 + sub_11 + sub_12 + sub_13 + sub_14 + sub_15 + sub_16 + sub_17 + sub_18 +
      sub_19 + sub_20 + sub_21 + sub_22 + sub_23 + sub_24

    Celf  <- aelftot / 3.697271e-05
    Clung <- (aelftot + aepitot + asubtot) / vlung
    Cbr   <- (elf_04 + elf_05 + elf_06 + elf_07 + elf_08 + elf_09) /
             (velf04 + velf05 + velf06 + velf07 + velf08 + velf09)

    Cc   ~ prop(propSd)
    Celf ~ prop(propSd_Celf)
  })
}
