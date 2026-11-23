Madrasi_2014_tenofovir_qsp <- function() {
  description <- paste(
    "QSP. Mechanistic model of tenofovir uptake into and phosphorylation",
    "inside peripheral blood mononuclear cells (PBMCs) after oral tenofovir",
    "disoproxil fumarate in adults: a two-compartment plasma tenofovir PK",
    "model with first-order absorption, saturable glutathione-driven",
    "antiporter uptake into the PBMC pool, saturable transporter efflux and",
    "a first-order MRP4 leak back to plasma, optional passive diffusion",
    "(permeability-surface product fixed at 0 in the final model), and",
    "sequential Michaelis-Menten phosphorylation by adenylate kinase 2 to",
    "tenofovir monophosphate and by nucleoside diphosphate kinase A to",
    "tenofovir diphosphate, which leaves the PBMC by first-order loss. All",
    "parameters are literature values or values the authors chose to",
    "reproduce saturable uptake; nothing was fitted. Between-subject",
    "variability is on plasma CL, Vc and Vp only."
  )
  reference <- paste(
    "Madrasi K, Burns RN, Hendrix CW, Fossler MJ, Chaturvedula A. Linking",
    "the population pharmacokinetics of tenofovir and its metabolites with",
    "its cellular uptake and metabolism. CPT Pharmacometrics Syst Pharmacol.",
    "2014;3(11):e147. doi:10.1038/psp.2014.46"
  )
  vignette <- "Madrasi_2014_tenofovir_qsp"

  # Doses are entered as milligrams of tenofovir disoproxil fumarate (TDF) and,
  # exactly as in the authors' MATLAB and NONMEM code (Supplement S8 / S9, where
  # a 300 mg tablet is the amount 0.3 g), the TDF mass passes into the plasma
  # tenofovir pool scaled only by the bioavailability Fbio = 0.32 -- no
  # TDF-to-tenofovir molecular-weight conversion is applied. All amount states
  # are carried in mg (the code carried g); every rate law below is converted
  # accordingly so that the authors' numerical outputs are reproduced exactly.
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(
      analyte = "tenofovir disoproxil fumarate",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(analyte = "tenofovir", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "tenofovir", units = "mg", specimen = "tissue", verified = TRUE),
    pbmc = list(analyte = "tenofovir", units = "mg", specimen = "blood cell", verified = TRUE),
    pbmc_tfvmp = list(analyte = "tenofovir monophosphate", units = "mg", specimen = "blood cell", verified = TRUE),
    pbmc_tfvdp = list(analyte = "tenofovir diphosphate", units = "mg", specimen = "blood cell", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    disease_state = "Healthy adults (HIV-negative), the population of the MTN-001 oral tenofovir disoproxil fumarate arm used as the external comparison.",
    dose_range = "Tenofovir disoproxil fumarate 75, 150, 300 and 600 mg orally, single dose and once daily for 10 days (deterministic simulations); 300 mg once daily for 14 days (population simulation).",
    notes = paste(
      "Not fitted to data. The model was assembled from literature values",
      "(Supplementary Table S1) and simulated deterministically in MATLAB and",
      "as a 100-subject population in NONMEM 7.2 (Methods 'Model validation').",
      "The population simulation was compared against steady-state plasma",
      "tenofovir and PBMC tenofovir-diphosphate observations from the oral arm",
      "of MTN-001 (Hendrix 2013, PLoS One 8:e55013), assuming 100 percent",
      "adherence. Plasma PK typical values are from Duwal 2012 (PLoS One",
      "7:e40382); between-subject variances are those printed in the authors'",
      "NONMEM code, attributed in the text to Baheti 2011."
    )
  )

  ini({
    # Plasma tenofovir PK. All values are literature values held fixed by the
    # authors (Supplementary Table S1, references Duwal 2012); none was estimated.
    lka <- fixed(log(1)); label("Absorption rate constant Ka (1/h)") # Table S1 'Ka' = 1 /h; NONMEM code KA=1
    lcl <- fixed(log(29.28)); label("Tenofovir clearance Cl (L/h)") # Table S1 'Cl' = 29.28 L/h; NONMEM $THETA(1) = 29.28
    lvc <- fixed(log(244)); label("Central volume Vc (L)") # Table S1 'Vc' = 244; NONMEM $THETA(2) = 244
    lvp <- fixed(log(464.54)); label("Peripheral volume Vp (L)") # Table S1 'Vp' = 464.54; NONMEM $THETA(3) = 464.54
    lq <- fixed(log(71.41)); label("Intercompartmental clearance Q (L/h)") # Table S1 'Q' = 71.41 L/h; NONMEM code Q=71.41
    lfdepot <- fixed(log(0.32)); label("Bioavailability Fbio (fraction)") # Table S1 'Fbio' = 0.32; NONMEM code FB=0.32
    fu <- fixed(0.93); label("Fraction unbound in plasma fu (fraction)") # Table S1 'fu' = 0.93; NONMEM code FU=0.93

    # PBMC physiology
    n_pbmc <- fixed(0.9e6); label("Number of PBMCs per mL of blood N (cells/mL)") # Table S1 'N' = 0.9 x 10^6 cells/ml
    v_blood <- fixed(5); label("Blood volume Vbl (L)") # Table S1 'Vbl' = 5 L; MATLAB code Vbl=5000 mL
    v_pbmc_cell <- fixed(0.28); label("Volume of one PBMC Vcell (pL)") # Table S1 'Vcell' = 0.28 pL; NONMEM code VC1 = 2.8e-13 L

    # Transporter uptake and efflux across the PBMC membrane
    vmax_pbmc <- fixed(1.77); label("Maximum rate of tenofovir uptake into PBMCs Vmax (g/(L*h))") # Table S1 'Vmax' = 1.77 g/(L*h); NONMEM code VMAX=1.77
    km_pbmc <- fixed(0.154); label("Michaelis constant of tenofovir PBMC uptake Km (uM)") # Table S1 'Km' assumed for model = 0.154 uM (770 uM lowered 5000-fold); NONMEM code KM=0.000044229 g/L = 0.154 uM x 287.2 g/mol
    keff_pbmc <- fixed(3.4e10); label("First-order MRP4 leak of tenofovir from PBMCs back to plasma kleak (1/h)") # Table S1 'kleak' assumed for model = 3.4 x 10^10 /h; NONMEM code KL=34000000000
    gsh_int <- fixed(3100); label("Intracellular glutathione concentration GSHi (uM)") # Table S1 'GSHi' = 3.1 mM; NONMEM code GSHI=3100
    gsh_ext <- fixed(20); label("Extracellular glutathione concentration GSHe (uM)") # Table S1 'GSHe' = 20 uM; NONMEM code GSHO=20
    ps_pbmc <- fixed(0); label("Total permeability-surface area product for passive diffusion PS_PBMC (L/h)") # Supplement text: diffusion eliminated by choosing PS_PBMC as 0; Table S1 'PScell' assumed for model = 0

    # Intracellular phosphorylation
    kcat_ak2 <- fixed(2.4); label("Turnover number of adenylate kinase 2 for tenofovir kcat1 (1/s)") # Table S1 'kcat1' = 2.4 s-1; NONMEM code KCAT1=2.4*3600
    km_ak2 <- fixed(3); label("Michaelis constant of adenylate kinase 2 for tenofovir Km1 (mM)") # Table S1 'Km1' = 3 mM; NONMEM code KM1=0.003 M
    enz_ak2 <- fixed(56); label("Adenylate kinase 2 concentration in PBMCs E1 (nM)") # Table S1 'E1' = 56 nM; NONMEM code E1=0.000000056 M
    kcat_ndpk <- fixed(0.12); label("Turnover number of nucleoside diphosphate kinase A for tenofovir monophosphate kcat2 (1/s)") # Table S1 'kcat2' = 0.12 s-1; NONMEM code KCAT2=0.12*3600
    km_ndpk <- fixed(0.29); label("Michaelis constant of nucleoside diphosphate kinase A for tenofovir monophosphate Km2 (mM)") # Table S1 'Km2' = 0.29 mM; NONMEM code KM2=0.00029 M
    enz_ndpk <- fixed(287.7); label("Nucleoside diphosphate kinase A concentration in PBMCs E2 (nM)") # Table S1 'E2' = 287.7 nM; NONMEM code E2=0.0000002877 M
    lkel_tfvdp <- fixed(log(0.0144375)); label("First-order loss of tenofovir diphosphate from PBMCs kout (1/h)") # NONMEM and MATLAB code KO / kout = 0.0144375 /h; Table S1 'kout' = 0.0144 /h (rounded)

    # Between-subject variability on the plasma PK, as printed in the authors'
    # NONMEM code ($OMEGA, diagonal); the text attributes it to Baheti 2011.
    etalcl ~ fixed(0.36) # NONMEM $OMEGA 0.36 on CL
    etalvc ~ fixed(0.79) # NONMEM $OMEGA 0.79 on V1
    etalvp ~ fixed(0.42) # NONMEM $OMEGA 0.42 on V2

    # Residual error, as printed in the authors' NONMEM code ($SIGMA, variances).
    propSd <- fixed(0.427785); label("Proportional residual error on plasma tenofovir (fraction)") # NONMEM $SIGMA 0.183 -> sqrt(0.183) = 0.427785
    propSd_Cpbmc_tfvdp <- fixed(0.752994); label("Proportional residual error on PBMC tenofovir diphosphate (fraction)") # NONMEM $SIGMA 0.567 -> sqrt(0.567) = 0.752994
  })
  model({
    # Molecular weights (g/mol), Table S1 and NONMEM code M1 / M2 / M3.
    mw_tfv <- 287.2
    mw_tfvmp <- 367.2
    mw_tfvdp <- 447.18

    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    q <- exp(lq)

    kel <- cl / vc
    kel_tfvdp <- exp(lkel_tfvdp)
    k12 <- q / vc
    k21 <- q / vp

    # Total PBMC volume (L): cells per mL x blood volume in mL x one-cell volume
    # in L (NONMEM code N = 4500000000, VC2 = VC1 * N).
    v_pbmc <- n_pbmc * v_blood * 1000 * v_pbmc_cell * 1e-12

    # Glutathione gradient driving the antiporter, (GSHi - GSHe) / GSHe
    # (Supplement Eq. S1h; NONMEM code GR).
    gsh_ratio <- (gsh_int - gsh_ext) / gsh_ext

    # Unit conversions for amounts in mg: Vmax g/(L*h) -> mg/(L*h); uptake Km
    # uM -> mg/L (uM x 1e-3 = mM, and mM x g/mol = mg/L); kcat 1/s -> 1/h;
    # enzyme nM -> mM. A PBMC amount in mg divided by (g/mol x L) is mM, the
    # unit of Km1 / Km2.
    vmax_mg <- vmax_pbmc * 1000
    km_mgl <- km_pbmc * 1e-3 * mw_tfv
    kcat1 <- kcat_ak2 * 3600
    kcat2 <- kcat_ndpk * 3600
    enz1 <- enz_ak2 * 1e-6
    enz2 <- enz_ndpk * 1e-6

    # Saturable antiporter influx from plasma and saturable efflux from PBMCs
    # (Supplement Eq. S1h, S1i as coded in Supplement S9: the unionised-fraction
    # factor (1 - F) = 1 - 4.2e-5 is omitted in the code).
    uptake <- vmax_mg * fu * central / (km_mgl + fu * central / vc) * gsh_ratio
    efflux <- vmax_mg * fu * pbmc / (km_mgl + fu * pbmc / v_pbmc)
    # Passive diffusion (Supplement Eq. S2a / S2c); zero in the final model.
    diffusion <- fu * ps_pbmc / vc * (central - pbmc)

    # Michaelis-Menten phosphorylation (Supplement Eq. S1f, S1g; NONMEM code).
    phos1 <- kcat1 * enz1 * pbmc / (km_ak2 + pbmc / (mw_tfv * v_pbmc))
    phos2 <- kcat2 * enz2 * pbmc_tfvmp / (km_ndpk + pbmc_tfvmp / (mw_tfvmp * v_pbmc))

    d/dt(depot) <- -ka * depot
    f(depot) <- exp(lfdepot)
    # The MRP4 leak returns tenofovir to plasma (Supplement Eq. S2a and the
    # MATLAB code); the printed NONMEM $DES omits the +KL*A(4) term.
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1 -
      uptake + efflux + keff_pbmc * pbmc - diffusion
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(pbmc) <- uptake - efflux - keff_pbmc * pbmc + diffusion - phos1
    d/dt(pbmc_tfvmp) <- phos1 - phos2
    d/dt(pbmc_tfvdp) <- phos2 - kel_tfvdp * pbmc_tfvdp

    # Plasma tenofovir (ng/mL): MATLAB CTFV = Y(1) * 1e6 / V1 with Y in g.
    Cc <- 1000 * central / vc
    # PBMC tenofovir monophosphate (ng/mL of PBMC volume): MATLAB TFVMP.
    Cpbmc_tfvmp <- 1000 * pbmc_tfvmp / v_pbmc
    # PBMC tenofovir diphosphate as the authors report it ('fmol/10^6 cells'):
    # MATLAB TFVDP = Y(5) * 1e15 / (M * Nu) with Nu = 1e6, NONMEM
    # C5 = A(6) * 1e9 / M3, i.e. the total amount in fmol divided by 1e6.
    Cpbmc_tfvdp <- 1e6 * pbmc_tfvdp / mw_tfvdp

    Cc ~ prop(propSd)
    Cpbmc_tfvdp ~ prop(propSd_Cpbmc_tfvdp)
  })
}
