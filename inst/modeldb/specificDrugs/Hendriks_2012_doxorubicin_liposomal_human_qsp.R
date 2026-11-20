Hendriks_2012_doxorubicin_liposomal_human_qsp <- function() {
  description <- paste(
    "QSP. Hendriks 2012 multiscale model of conventional and pegylated",
    "liposomal doxorubicin (PLD) delivery to tumor-cell DNA, scaled from the",
    "mouse model to human physiology (70 kg, 1.8 m2, 20 g tumor) with PK",
    "and tumor-transport parameters re-estimated from clinical data. Free",
    "doxorubicin follows a two-compartment PK model and PLD a one-compartment",
    "model (an optional liposome peripheral compartment is switched off by",
    "k12_lipo = 0) sharing the blood volume as central compartment;",
    "doxorubicin is released from liposomes in blood,",
    "tissue interstitium and cells. A perfused, well-mixed tumor (capillary,",
    "interstitium, cells) and heart (same structure) receive drug by blood",
    "flow and transvascular flux; liposomes enter cells by fluid-phase",
    "uptake. Inside cells doxorubicin partitions into the outer membrane",
    "leaflet, is flipped to the inner leaflet, enters the cytosol and binds",
    "DNA reversibly (modified Eytan-Kuchel cell model). Deterministic",
    "typical-value model (no IIV, no residual error). Amounts in nmol of",
    "doxorubicin; liposome and empty-lipid states are in doxorubicin-",
    "equivalent nmol (20,000 doxorubicin molecules per liposome). Time in",
    "minutes."
  )
  reference <- paste(
    "Hendriks BS, Reynolds JG, Klinz SG, Geretti E, Lee H, Leonard SC,",
    "Gaddy DF, Espelin CW, Nielsen UB, Wickham TJ. (2012). Multiscale",
    "kinetic modeling of liposomal doxorubicin delivery quantifies the role",
    "of tumor and drug-specific parameters in local delivery to tumors.",
    "CPT Pharmacometrics Syst Pharmacol 1(11):e15. doi:10.1038/psp.2012.16.",
    "Structure and fixed physiological constants follow the authors'",
    "SimBiology model (Supplementary Data, SBML PSP-2012-0011-s02) and its",
    "human parameter and reaction workbook (supplementary Excel file);",
    "fitted values follow Hendriks 2012 Table 1, except the doxorubicin",
    "tumor efflux flux, taken from the model code (Table 1 prints it",
    "ten-fold too high)."
  )
  vignette <- "Hendriks_2012_doxorubicin_liposomal"

  # Free doxorubicin enters in `central`; pegylated liposomal doxorubicin in
  # `central_lipo` (as doxorubicin-equivalent nmol).
  dosing <- c("central", "central_lipo")

  # Mechanistic tumor / heart sub-compartments of the Hendriks 2012 model
  # with no canonical analogue in inst/references/compartment-names.md.
  # cap = capillary blood, int = interstitium, cell_* = tumor or heart cells
  # (ol / il = outer / inner membrane leaflet, cyt = cytosol, dna = free
  # DNA binding sites, doxdna = DNA-bound doxorubicin, lipo / lipid =
  # internalised liposomes / empty lipid).
  paper_specific_compartments <- c(
    "central_lipo",
    "central_lipid",
    "peripheral1_lipo",
    "cap_tumor_dox",
    "cap_tumor_lipo",
    "cap_tumor_lipid",
    "int_tumor_dox",
    "int_tumor_lipo",
    "int_tumor_lipid",
    "cell_tumor_ol",
    "cell_tumor_il",
    "cell_tumor_cyt",
    "cell_tumor_dna",
    "cell_tumor_doxdna",
    "cell_tumor_lipo",
    "cell_tumor_lipid",
    "cap_heart_dox",
    "cap_heart_lipo",
    "cap_heart_lipid",
    "int_heart_dox",
    "int_heart_lipo",
    "int_heart_lipid",
    "cell_heart_ol",
    "cell_heart_il",
    "cell_heart_cyt",
    "cell_heart_dna",
    "cell_heart_doxdna",
    "cell_heart_lipo",
    "cell_heart_lipid"
  )

  units <- list(
    time = "min",
    dosing = "nmol",
    concentration = "nM"
  )

  compartmentData <- list(
    central = list(analyte = "doxorubicin", units = "nmol", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "doxorubicin", units = "nmol", specimen = "tissue", verified = TRUE),
    central_lipo = list(
      analyte = "liposome-encapsulated doxorubicin",
      units = "nmol",
      specimen = "whole blood",
      verified = TRUE
    ),
    central_lipid = list(
      analyte = "empty liposomal lipid (doxorubicin-equivalent)",
      units = "nmol",
      specimen = "whole blood",
      verified = TRUE
    ),
    peripheral1_lipo = list(
      analyte = "liposome-encapsulated doxorubicin",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    cap_tumor_dox = list(analyte = "doxorubicin", units = "nmol", specimen = "tumor", verified = TRUE),
    cap_tumor_lipo = list(
      analyte = "liposome-encapsulated doxorubicin",
      units = "nmol",
      specimen = "tumor",
      verified = TRUE
    ),
    cap_tumor_lipid = list(
      analyte = "empty liposomal lipid (doxorubicin-equivalent)",
      units = "nmol",
      specimen = "tumor",
      verified = TRUE
    ),
    int_tumor_dox = list(analyte = "doxorubicin", units = "nmol", specimen = "tumor", verified = TRUE),
    int_tumor_lipo = list(
      analyte = "liposome-encapsulated doxorubicin",
      units = "nmol",
      specimen = "tumor",
      verified = TRUE
    ),
    int_tumor_lipid = list(
      analyte = "empty liposomal lipid (doxorubicin-equivalent)",
      units = "nmol",
      specimen = "tumor",
      verified = TRUE
    ),
    cell_tumor_ol = list(
      analyte = "doxorubicin (outer membrane leaflet)",
      units = "nmol",
      specimen = "tumor",
      verified = TRUE
    ),
    cell_tumor_il = list(
      analyte = "doxorubicin (inner membrane leaflet)",
      units = "nmol",
      specimen = "tumor",
      verified = TRUE
    ),
    cell_tumor_cyt = list(analyte = "doxorubicin (cytosol)", units = "nmol", specimen = "tumor", verified = TRUE),
    cell_tumor_dna = list(analyte = "free DNA binding sites", units = "nmol", specimen = "tumor", verified = TRUE),
    cell_tumor_doxdna = list(analyte = "DNA-bound doxorubicin", units = "nmol", specimen = "tumor", verified = TRUE),
    cell_tumor_lipo = list(
      analyte = "internalised liposome-encapsulated doxorubicin",
      units = "nmol",
      specimen = "tumor",
      verified = TRUE
    ),
    cell_tumor_lipid = list(
      analyte = "internalised empty liposomal lipid (doxorubicin-equivalent)",
      units = "nmol",
      specimen = "tumor",
      verified = TRUE
    ),
    cap_heart_dox = list(analyte = "doxorubicin", units = "nmol", specimen = "tissue", verified = TRUE),
    cap_heart_lipo = list(
      analyte = "liposome-encapsulated doxorubicin",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    cap_heart_lipid = list(
      analyte = "empty liposomal lipid (doxorubicin-equivalent)",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    int_heart_dox = list(analyte = "doxorubicin", units = "nmol", specimen = "tissue", verified = TRUE),
    int_heart_lipo = list(
      analyte = "liposome-encapsulated doxorubicin",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    int_heart_lipid = list(
      analyte = "empty liposomal lipid (doxorubicin-equivalent)",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    cell_heart_ol = list(
      analyte = "doxorubicin (outer membrane leaflet)",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    cell_heart_il = list(
      analyte = "doxorubicin (inner membrane leaflet)",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    cell_heart_cyt = list(analyte = "doxorubicin (cytosol)", units = "nmol", specimen = "tissue", verified = TRUE),
    cell_heart_dna = list(analyte = "free DNA binding sites", units = "nmol", specimen = "tissue", verified = TRUE),
    cell_heart_doxdna = list(analyte = "DNA-bound doxorubicin", units = "nmol", specimen = "tissue", verified = TRUE),
    cell_heart_lipo = list(
      analyte = "internalised liposome-encapsulated doxorubicin",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    cell_heart_lipid = list(
      analyte = "internalised empty liposomal lipid (doxorubicin-equivalent)",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    )
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = "Literature data; typical-value model scaled from the mouse model (no individual patients modelled jointly).",
    n_studies = paste(
      "Free doxorubicin plasma PK compiled from the published literature",
      "(Hendriks 2012 Figure 3g); PLD plasma PK and tumor deposition from",
      "Harrington 2001 (111In-labelled PLD, various solid tumors) and Northfelt",
      "1996 (Kaposi sarcoma); tumor Ktrans values from three DCE-MRI studies"
    ),
    weight_range = "Typical 70 kg adult with 1.8 m2 body surface area (model code)",
    disease_state = "Adults with solid tumors (breast cancer, Kaposi sarcoma and other tumor types)",
    dose_range = "Doxorubicin 10-60 mg/m2 IV; pegylated liposomal doxorubicin 40-50 mg/m2 IV (simulation scenarios)",
    regions = "Merrimack Pharmaceuticals, Cambridge MA, USA (model development)",
    notes = paste(
      "Human scale-up of the mouse model (Hendriks 2012 Methods, 'Human",
      "scale-up'): blood volume, body weight, tumor blood flow and dose were",
      "set to human values; free-doxorubicin PK, liposome elimination and",
      "the doxorubicin tumor efflux flux were estimated from clinical data;",
      "all cellular, release and liposome transport parameters were carried",
      "from the mouse model. tvf_in_dox was fixed from DCE-MRI Ktrans values",
      "and tvf_in_lipo was estimated per patient (Figure 3l)."
    )
  )

  ini({
    # ==================================================================
    # Free doxorubicin plasma PK (two-compartment; central = blood volume).
    # Hendriks 2012 Table 1, 'Human parameters' (Fit).
    # ==================================================================
    lkel <- log(0.0820022)
    label("Doxorubicin elimination rate constant from blood kel_dox (1/min)") # Table 1 human kel_dox = 0.0820 1/min (Fit); model code 0.0820022
    lk12 <- log(0.0475101)
    label("Doxorubicin blood-to-peripheral rate constant k12_dox (1/min)") # Table 1 human k12_dox = 0.0475 1/min (Fit); model code 0.0475101
    lk21 <- log(0.00124661)
    label("Doxorubicin peripheral-to-blood rate constant k21_dox (1/min)") # Table 1 human k21_dox = 0.00125 1/min (Fit); model code 0.00124661

    # ==================================================================
    # Liposome PK and doxorubicin release.
    # ==================================================================
    kel_lipo <- 1.67e-4
    label("Liposome elimination rate constant from blood kel_lipo (1/min)") # Table 1 human kel_lipo = 1.67E-04 1/min (ref 8, Harrington); Supplementary human workbook 1.68862e-4
    k12_lipo <- fixed(0)
    label("Liposome blood-to-peripheral rate constant (1/min)") # Supplementary human workbook k12_lipo = 0 (one-compartment liposome PK; Methods 'Human scale-up')
    k21_lipo <- fixed(0.00341)
    label("Liposome peripheral-to-blood rate constant (1/min)") # Supplementary human workbook k21_lipo = 0.00341 1/min (inactive while k12_lipo = 0)
    krel_plasma <- fixed(5.465e-5)
    label("Doxorubicin release rate from liposomes in plasma and capillaries krel_plasma (1/min)") # Methods 5.47e-5 1/min; Table 1 5.47E-05; model code 5.465e-5 (ref 31)
    krel_tissue <- fixed(1.77e-4)
    label("Doxorubicin release rate from liposomes in tissue interstitium krel_tumor (1/min)") # Methods / Table 1 krel_tumor = 1.77E-04 1/min (ref 31)
    krel_cell <- fixed(9.6e-4)
    label("Doxorubicin release rate from internalised liposomes krel_cell (1/min)") # Methods / Table 1 krel_cell = 9.6E-04 1/min (12-h half-life)
    ke_fluidphase <- fixed(4e-14)
    label("Fluid-phase liposome uptake rate per cell ke_fluidphase (L/min/cell)") # Table 1 ke_fluidphase = 4E-14 L/min/cell (ref 32)
    kdeg_lipo <- fixed(0.0042)
    label("Exocytosis rate of internalised liposomes and lipid kdeg_lipo (1/min)") # Supplementary human workbook / model code kdeg_lipo = 0.0042 1/min
    doxperliposome <- fixed(20000)
    label("Doxorubicin molecules per liposome (count)") # Table 1 doxperliposome = 20,000; Methods stoichiometry 20,000:1

    # ==================================================================
    # Tumor transvascular flux per unit capillary surface area. Table 1
    # reports cm/min; the model works in dm (1 cm/min = 0.1 dm/min).
    # ==================================================================
    tvf_in_tumor_dox <- fixed(3.63e-4)
    label("Doxorubicin transvascular flux capillary to tumor interstitium tvf_in_dox (dm/min)") # Table 1 human tvf_in_dox = 3.63E-03 cm/min (DCE-MRI Ktrans, refs 22-24)
    tvf_out_tumor_dox <- 8.45e-5
    label("Doxorubicin transvascular flux tumor interstitium to capillary tvf_out_dox (dm/min)") # Supplementary human workbook / model code 8.45e-5 dm/min (Northfelt fit); Table 1 prints 8.45E-03 cm/min, see vignette
    tvf_in_tumor_lipo <- 2.64e-7
    label("Liposome transvascular flux capillary to tumor interstitium tvf_in_lipo (dm/min)") # Table 1 tvf_in_lipo = 2.64E-06 cm/min; model code 2.64e-7 dm/min (per-patient estimates in Figure 3l)
    tvf_out_tumor_lipo <- fixed(7.14e-7)
    label("Liposome transvascular flux tumor interstitium to capillary tvf_out_lipo (dm/min)") # Table 1 tvf_out_lipo = 7.14E-06 cm/min (assumed equal in human and mouse)

    # ==================================================================
    # Cellular doxorubicin transport (modified Eytan-Kuchel model).
    # ==================================================================
    kf_ol <- fixed(3.3e-6)
    label("Doxorubicin mass-transfer coefficient medium/cytosol to membrane kf_ol (dm/min)") # Table 1 kf_ol = 3.30E-05 cm/min (median cell-line fit)
    kflippase_in <- fixed(1)
    label("Inward flippase rate constant kflippase_in (1/min)") # Table 1 kflippase_in = 1 1/min (ref 20)
    kflippase_out <- fixed(0.37)
    label("Outward flippase rate constant kflippase_out (1/min)") # Table 1 kflippase_out = 0.37 1/min (median cell-line fit)
    pc_dox <- fixed(40000)
    label("Doxorubicin membrane/aqueous partition coefficient (unitless)") # Table 1 doxpartitioncoefficient = 40,000
    frxnvol_ol <- fixed(5e-5)
    label("Outer-leaflet volume as a fraction of cell volume (unitless)") # Table 1 frxnvol_outerleaflet = 0.00005
    frxnvol_il <- fixed(5e-5)
    label("Inner-leaflet volume as a fraction of cell volume (unitless)") # Table 1 frxnvol_innerleaflet = 0.00005
    kd_dna <- fixed(1000)
    label("Doxorubicin-DNA dissociation constant Kd_DNA (nM)") # Table 1 Kd_DNA = 1.00E03 nmol/l
    kr_dna <- fixed(1200)
    label("Doxorubicin-DNA complex dissociation rate constant kr_DNA (1/min)") # Table 1 kr_DNA = 1200 1/min
    nucleotides_per_cell <- fixed(4.58e7)
    label("Effective nucleotides available for doxorubicin binding per cell (count/cell)") # Table 1 nucleotides_per_cell_tumor = 4.58E+07 (also used for heart in model code)
    sites_per_nucleotide <- fixed(0.18)
    label("Doxorubicin binding sites per nucleotide (unitless)") # Table 1 bindingsites_per_nucleotide = 0.18

    # ==================================================================
    # Human and tumor physiology. Lengths in dm, volumes in L, flows per
    # kg tissue. Tumor geometry is carried from the mouse model.
    # ==================================================================
    bloodvolume <- fixed(4.85)
    label("Human blood volume (L)") # Supplementary human workbook / model code bloodvolume = 4.85 L
    tissuedensity <- fixed(1)
    label("Tissue density (kg/L)") # Table 1 tissuedensity = 1 kg/l
    mw_dox <- fixed(580)
    label("Doxorubicin molecular weight used for unit conversion (g/mol)") # model code MW_dox = 580 g/mol
    tumorweight <- fixed(0.02)
    label("Tumor weight (kg)") # Supplementary human workbook / model code tumorweight = 0.02 kg
    q_tumor <- fixed(0.0282)
    label("Blood flow into tumor per kg tumor Qtumor (L/min/kg)") # Table 1 human Qtumor = 2.82E-02 l/min/kg
    l_tumor <- fixed(0)
    label("Lymph flow out of tumor per kg tumor Ltumor (L/min/kg)") # Table 1 Ltumor = 0 l/min/kg
    mvd_tumor <- fixed(250)
    label("Microvessel density in tumor MVD_tumor (vessels/mm^2)") # Table 1 MVD_tumor = 250 #/mm2
    capradius3_tumor <- fixed(3e-5)
    label("Tumor capillary radius class 1 (dm)") # Supplementary human workbook capillaryradius3_tumor = 3e-5 dm (3 um)
    capradius5_tumor <- fixed(5e-5)
    label("Tumor capillary radius class 2 (dm)") # Supplementary human workbook capillaryradius5_tumor = 5e-5 dm (5 um); Table 1 capillaryradius_tumor = 5 um
    capradius7_tumor <- fixed(7e-5)
    label("Tumor capillary radius class 3 (dm)") # Supplementary human workbook capillaryradius7_tumor = 7e-5 dm (7 um)
    capradius23_tumor <- fixed(2.3e-4)
    label("Tumor arteriole radius, non-permeable (dm)") # Supplementary human workbook capillaryradius23_tumor = 2.3e-4 dm (23 um)
    radfrxn3_tumor <- fixed(0.1)
    label("Fraction of tumor vessels in radius class 1 (unitless)") # Supplementary human workbook radiusfrxn3_tumor = 0.1
    radfrxn5_tumor <- fixed(0.6)
    label("Fraction of tumor vessels in radius class 2 (unitless)") # Supplementary human workbook radiusfrxn5_tumor = 0.6
    radfrxn7_tumor <- fixed(0.2)
    label("Fraction of tumor vessels in radius class 3 (unitless)") # Supplementary human workbook radiusfrxn7_tumor = 0.2
    radfrxn23_tumor <- fixed(0.1)
    label("Fraction of tumor vessels that are arterioles (unitless)") # Supplementary human workbook radiusfrxn23_tumor = 0.1
    cellfrxn_tumor <- fixed(0.412)
    label("Cellular fraction of the non-vascular tumor volume (unitless)") # Table 1 tumor_cell2interstitium_ratio = 0.412 (model code tumor_cell2nonvascular_ratio)
    cellradius_tumor <- fixed(8.5e-5)
    label("Tumor cell radius (dm)") # Table 1 cellradius_tumor = 8.5 um

    # ==================================================================
    # Heart physiology and doxorubicin / liposome heart transport (heart
    # sub-model of the authors' code; Supplementary human workbook).
    # ==================================================================
    heartweight <- fixed(0.33)
    label("Heart weight (kg)") # Supplementary human workbook heartweight = 0.33 kg
    q_heart <- fixed(0.72727273)
    label("Blood flow into heart per kg heart Qheart (L/min/kg)") # Supplementary human workbook Qheart = 0.727273 L/min/kg
    l_heart <- fixed(7.5e-5)
    label("Lymph flow out of heart per kg heart Lheart (L/min/kg)") # Supplementary human workbook Lheart = 7.5e-5 L/min/kg
    mvd_heart <- fixed(1650)
    label("Microvessel density in heart MVD_heart (vessels/mm^2)") # Supplementary human workbook MVD_heart = 1650 #/mm2
    capradius_heart <- fixed(2.2e-5)
    label("Heart capillary radius (dm)") # Supplementary human workbook capillaryradius_heart = 2.2e-5 dm
    cellfrxn_heart <- fixed(0.849)
    label("Cellular fraction of the non-vascular heart volume (unitless)") # Supplementary human workbook heart_cell2nonvascular_ratio = 0.849
    cellradius_heart <- fixed(1.56e-4)
    label("Cardiomyocyte radius (dm)") # Supplementary human workbook cellradius_heart = 0.000156 dm
    tvf_in_heart_dox <- fixed(3.63e-4)
    label("Doxorubicin transvascular flux capillary to heart interstitium (dm/min)") # Supplementary human workbook tvf_IN_heart_dox = 3.63e-4 dm/min (same as tumor)
    tvf_out_heart_dox <- fixed(2.0862069e-4)
    label("Doxorubicin transvascular flux heart interstitium to capillary (dm/min)") # Supplementary human workbook tvf_OUT_heart_dox = 2.086e-4 dm/min (mouse in/out ratio 1.74)
    tvf_in_heart_lipo <- fixed(2.65e-8)
    label("Liposome transvascular flux capillary to heart interstitium (dm/min)") # Supplementary human workbook tvf_IN_heart_lipo = 2.65e-8 dm/min
    tvf_out_heart_lipo <- fixed(0)
    label("Liposome transvascular flux heart interstitium to capillary (dm/min)") # Supplementary human workbook tvf_OUT_heart_lipo = 0 dm/min
    kdeg_dox <- fixed(0.014)
    label("Doxorubicin degradation rate in heart interstitium kdeg_dox (1/min)") # Supplementary human workbook kdeg_dox = 0.014 1/min
    residual_blood_frxn <- fixed(0.5)
    label("Fraction of heart capillary blood counted in heart-tissue readouts (unitless)") # Supplementary human workbook ResidualBloodFrxn = 0.5
  })
  model({
    kel <- exp(lkel)
    k12 <- exp(lk12)
    k21 <- exp(lk21)

    # --- Unit constants: Avogadro's number per nmol; mm per dm. ----------
    nav <- 6.022e14
    mm_per_dm <- 100

    kf_dna <- kr_dna / kd_dna
    kf_il <- kf_ol

    # --- Tumor geometry (assignment rules of the authors' code). ----------
    # The effective capillary radius excludes the non-permeable arterioles
    # from the numerator but not from the fraction sum; the vascular volume
    # includes them.
    radfrxn_sum <- radfrxn3_tumor + radfrxn5_tumor + radfrxn7_tumor + radfrxn23_tumor
    capradius_eff_tumor <- (capradius3_tumor * radfrxn3_tumor + capradius5_tumor * radfrxn5_tumor +
      capradius7_tumor * radfrxn7_tumor) / radfrxn_sum
    capsa_tumor <- mvd_tumor * mm_per_dm^2 * 2 * pi * capradius_eff_tumor
    fvasc_tumor <- mvd_tumor * mm_per_dm^2 * pi * (capradius3_tumor^2 * radfrxn3_tumor +
      capradius5_tumor^2 * radfrxn5_tumor + capradius7_tumor^2 * radfrxn7_tumor +
      capradius23_tumor^2 * radfrxn23_tumor) / radfrxn_sum
    vcap_tumor <- tumorweight / tissuedensity * fvasc_tumor
    vint_tumor <- tumorweight / tissuedensity * (1 - fvasc_tumor) * (1 - cellfrxn_tumor)
    vcell_tumor <- tumorweight / tissuedensity * (1 - fvasc_tumor) * cellfrxn_tumor
    cellvol_tumor <- 4 / 3 * pi * cellradius_tumor^3
    cellsa_tumor <- 4 * pi * cellradius_tumor^2
    ncell_tumor <- vcell_tumor / cellvol_tumor

    # --- Heart geometry. ------------------------------------------------
    capsa_heart <- mvd_heart * mm_per_dm^2 * 2 * capradius_heart * pi
    fvasc_heart <- mvd_heart * mm_per_dm^2 * pi * capradius_heart^2
    vcap_heart <- heartweight / tissuedensity * fvasc_heart
    vint_heart <- heartweight / tissuedensity * (1 - fvasc_heart) * (1 - cellfrxn_heart)
    vcell_heart <- heartweight / tissuedensity * (1 - fvasc_heart) * cellfrxn_heart
    cellvol_heart <- 4 / 3 * pi * cellradius_heart^3
    cellsa_heart <- 4 * pi * cellradius_heart^2
    ncell_heart <- vcell_heart / cellvol_heart

    # --- Concentrations (nM) in each space. -----------------------------
    c_dox <- central / bloodvolume
    c_lipo <- central_lipo / bloodvolume
    c_lipid <- central_lipid / bloodvolume
    c_cap_tumor_dox <- cap_tumor_dox / vcap_tumor
    c_cap_tumor_lipo <- cap_tumor_lipo / vcap_tumor
    c_cap_tumor_lipid <- cap_tumor_lipid / vcap_tumor
    c_int_tumor_dox <- int_tumor_dox / vint_tumor
    c_int_tumor_lipo <- int_tumor_lipo / vint_tumor
    c_int_tumor_lipid <- int_tumor_lipid / vint_tumor
    c_cell_tumor_ol <- cell_tumor_ol / vcell_tumor
    c_cell_tumor_il <- cell_tumor_il / vcell_tumor
    c_cell_tumor_cyt <- cell_tumor_cyt / vcell_tumor
    c_cell_tumor_dna <- cell_tumor_dna / vcell_tumor
    c_cap_heart_dox <- cap_heart_dox / vcap_heart
    c_cap_heart_lipo <- cap_heart_lipo / vcap_heart
    c_cap_heart_lipid <- cap_heart_lipid / vcap_heart
    c_int_heart_dox <- int_heart_dox / vint_heart
    c_int_heart_lipo <- int_heart_lipo / vint_heart
    c_int_heart_lipid <- int_heart_lipid / vint_heart
    c_cell_heart_ol <- cell_heart_ol / vcell_heart
    c_cell_heart_il <- cell_heart_il / vcell_heart
    c_cell_heart_cyt <- cell_heart_cyt / vcell_heart
    c_cell_heart_dna <- cell_heart_dna / vcell_heart

    # --- Tumor fluxes (nmol/min). --------------------------------------
    flow_in_tumor <- q_tumor * tumorweight
    flow_out_tumor <- (q_tumor - l_tumor) * tumorweight
    lymph_tumor <- l_tumor * tumorweight
    tv_tumor_dox <- (tvf_in_tumor_dox * c_cap_tumor_dox - tvf_out_tumor_dox * c_int_tumor_dox) *
      capsa_tumor * tumorweight
    tv_tumor_lipo <- (tvf_in_tumor_lipo * c_cap_tumor_lipo - tvf_out_tumor_lipo * c_int_tumor_lipo) *
      capsa_tumor * tumorweight
    tv_tumor_lipid <- (tvf_in_tumor_lipo * c_cap_tumor_lipid - tvf_out_tumor_lipo * c_int_tumor_lipid) *
      capsa_tumor * tumorweight
    part_tumor_ol <- kf_ol * cellsa_tumor * ncell_tumor *
      (c_int_tumor_dox - c_cell_tumor_ol / pc_dox / frxnvol_ol)
    flip_tumor <- kflippase_in / frxnvol_ol * cell_tumor_ol - kflippase_out / frxnvol_il * cell_tumor_il
    part_tumor_il <- kf_il * cellsa_tumor * ncell_tumor *
      (c_cell_tumor_cyt - c_cell_tumor_il / pc_dox / frxnvol_il)
    bind_tumor <- kf_dna * c_cell_tumor_cyt * c_cell_tumor_dna * vcell_tumor - kr_dna * cell_tumor_doxdna
    uptake_tumor_lipo <- ke_fluidphase * ncell_tumor * c_int_tumor_lipo

    # --- Heart fluxes (nmol/min). --------------------------------------
    flow_in_heart <- q_heart * heartweight
    flow_out_heart <- (q_heart - l_heart) * heartweight
    lymph_heart <- l_heart * heartweight
    tv_heart_dox <- (tvf_in_heart_dox * c_cap_heart_dox - tvf_out_heart_dox * c_int_heart_dox) *
      capsa_heart * heartweight
    tv_heart_lipo <- (tvf_in_heart_lipo * c_cap_heart_lipo - tvf_out_heart_lipo * c_int_heart_lipo) *
      capsa_heart * heartweight
    tv_heart_lipid <- (tvf_in_heart_lipo * c_cap_heart_lipid - tvf_out_heart_lipo * c_int_heart_lipid) *
      capsa_heart * heartweight
    part_heart_ol <- kf_ol * cellsa_heart * ncell_heart *
      (c_int_heart_dox - c_cell_heart_ol / pc_dox / frxnvol_ol)
    flip_heart <- kflippase_in / frxnvol_ol * cell_heart_ol - kflippase_out / frxnvol_il * cell_heart_il
    part_heart_il <- kf_il * cellsa_heart * ncell_heart *
      (c_cell_heart_cyt - c_cell_heart_il / pc_dox / frxnvol_il)
    bind_heart <- kf_dna * c_cell_heart_cyt * c_cell_heart_dna * vcell_heart - kr_dna * cell_heart_doxdna
    uptake_heart_lipo <- ke_fluidphase * ncell_heart * c_int_heart_lipo

    # --- Blood: free doxorubicin (two-compartment), liposomes, lipid. --
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 + krel_plasma * central_lipo -
      flow_in_tumor * c_dox + flow_out_tumor * c_cap_tumor_dox + lymph_tumor * c_int_tumor_dox -
      flow_in_heart * c_dox + flow_out_heart * c_cap_heart_dox + lymph_heart * c_int_heart_dox
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(central_lipo) <- -kel_lipo * central_lipo - krel_plasma * central_lipo -
      k12_lipo * central_lipo + k21_lipo * peripheral1_lipo -
      flow_in_tumor * c_lipo + flow_out_tumor * c_cap_tumor_lipo + lymph_tumor * c_int_tumor_lipo -
      flow_in_heart * c_lipo + flow_out_heart * c_cap_heart_lipo + lymph_heart * c_int_heart_lipo
    d/dt(peripheral1_lipo) <- k12_lipo * central_lipo - k21_lipo * peripheral1_lipo
    d/dt(central_lipid) <- krel_plasma * central_lipo - kel_lipo * central_lipid -
      flow_in_tumor * c_lipid + flow_out_tumor * c_cap_tumor_lipid + lymph_tumor * c_int_tumor_lipid -
      flow_in_heart * c_lipid + flow_out_heart * c_cap_heart_lipid + lymph_heart * c_int_heart_lipid

    # --- Tumor capillary. ----------------------------------------------
    d/dt(cap_tumor_dox) <- flow_in_tumor * c_dox - flow_out_tumor * c_cap_tumor_dox - tv_tumor_dox +
      krel_plasma * cap_tumor_lipo
    d/dt(cap_tumor_lipo) <- flow_in_tumor * c_lipo - flow_out_tumor * c_cap_tumor_lipo - tv_tumor_lipo -
      krel_plasma * cap_tumor_lipo
    d/dt(cap_tumor_lipid) <- flow_in_tumor * c_lipid - flow_out_tumor * c_cap_tumor_lipid - tv_tumor_lipid +
      krel_plasma * cap_tumor_lipo

    # --- Tumor interstitium. -------------------------------------------
    d/dt(int_tumor_dox) <- tv_tumor_dox - lymph_tumor * c_int_tumor_dox + krel_tissue * int_tumor_lipo -
      part_tumor_ol
    d/dt(int_tumor_lipo) <- tv_tumor_lipo - lymph_tumor * c_int_tumor_lipo - krel_tissue * int_tumor_lipo -
      uptake_tumor_lipo + kdeg_lipo * cell_tumor_lipo
    d/dt(int_tumor_lipid) <- tv_tumor_lipid - lymph_tumor * c_int_tumor_lipid + krel_tissue * int_tumor_lipo +
      kdeg_lipo * cell_tumor_lipid

    # --- Tumor cells. ---------------------------------------------------
    d/dt(cell_tumor_ol) <- part_tumor_ol - flip_tumor
    d/dt(cell_tumor_il) <- flip_tumor + part_tumor_il
    d/dt(cell_tumor_cyt) <- -part_tumor_il - bind_tumor + krel_cell * cell_tumor_lipo
    d/dt(cell_tumor_dna) <- -bind_tumor
    d/dt(cell_tumor_doxdna) <- bind_tumor
    d/dt(cell_tumor_lipo) <- uptake_tumor_lipo - krel_cell * cell_tumor_lipo - kdeg_lipo * cell_tumor_lipo
    d/dt(cell_tumor_lipid) <- krel_cell * cell_tumor_lipo - kdeg_lipo * cell_tumor_lipid
    cell_tumor_dna(0) <- nucleotides_per_cell * sites_per_nucleotide / nav * ncell_tumor

    # --- Heart capillary. ----------------------------------------------
    d/dt(cap_heart_dox) <- flow_in_heart * c_dox - flow_out_heart * c_cap_heart_dox - tv_heart_dox +
      krel_plasma * cap_heart_lipo
    d/dt(cap_heart_lipo) <- flow_in_heart * c_lipo - flow_out_heart * c_cap_heart_lipo - tv_heart_lipo -
      krel_plasma * cap_heart_lipo
    d/dt(cap_heart_lipid) <- flow_in_heart * c_lipid - flow_out_heart * c_cap_heart_lipid - tv_heart_lipid +
      krel_plasma * cap_heart_lipo

    # --- Heart interstitium. -------------------------------------------
    d/dt(int_heart_dox) <- tv_heart_dox - lymph_heart * c_int_heart_dox + krel_tissue * int_heart_lipo -
      part_heart_ol - kdeg_dox * int_heart_dox
    d/dt(int_heart_lipo) <- tv_heart_lipo - lymph_heart * c_int_heart_lipo - krel_tissue * int_heart_lipo -
      uptake_heart_lipo + kdeg_lipo * cell_heart_lipo
    d/dt(int_heart_lipid) <- tv_heart_lipid - lymph_heart * c_int_heart_lipid + krel_tissue * int_heart_lipo +
      kdeg_lipo * cell_heart_lipid

    # --- Heart cells. ---------------------------------------------------
    d/dt(cell_heart_ol) <- part_heart_ol - flip_heart
    d/dt(cell_heart_il) <- flip_heart + part_heart_il
    d/dt(cell_heart_cyt) <- -part_heart_il - bind_heart + krel_cell * cell_heart_lipo
    d/dt(cell_heart_dna) <- -bind_heart
    d/dt(cell_heart_doxdna) <- bind_heart
    d/dt(cell_heart_lipo) <- uptake_heart_lipo - krel_cell * cell_heart_lipo - kdeg_lipo * cell_heart_lipo
    d/dt(cell_heart_lipid) <- krel_cell * cell_heart_lipo - kdeg_lipo * cell_heart_lipid
    cell_heart_dna(0) <- nucleotides_per_cell * sites_per_nucleotide / nav * ncell_heart

    # --- Readouts (repeated assignments of the authors' code). ---------
    # nmol x mw_dox / 1e3 = ug; tissue weight kg x 1e3 = g.
    ug_per_g_tumor <- mw_dox / 1e3 / (tumorweight * 1e3)
    ug_per_g_heart <- mw_dox / 1e3 / (heartweight * 1e3)
    cell_tumor_dox <- cell_tumor_ol + cell_tumor_il + cell_tumor_cyt + cell_tumor_doxdna
    cell_heart_dox <- cell_heart_ol + cell_heart_il + cell_heart_cyt + cell_heart_doxdna

    # Total plasma doxorubicin, free plus liposome-encapsulated (PlasmaDox_nM).
    Cc <- c_dox + c_lipo
    # Free and liposome-encapsulated plasma doxorubicin, and total plasma
    # lipid label in doxorubicin-equivalents (nM).
    Cc_free <- c_dox
    Cc_lipo <- c_lipo
    Cc_lipid_total <- c_lipo + c_lipid
    # Tumor readouts (ug doxorubicin per g tumor), capillary blood excluded
    # (TumorTotalDox_ug_per_g, TumorTotalLipid_ugDox_per_g,
    # TumorDNABoundDox_ug_per_g, TumorTotalCellDox_ug_per_g).
    tumor_total_dox <- (int_tumor_dox + int_tumor_lipo + cell_tumor_dox + cell_tumor_lipo) * ug_per_g_tumor
    tumor_total_dox_cap <- tumor_total_dox + (cap_tumor_dox + cap_tumor_lipo) * ug_per_g_tumor
    tumor_total_lipid <- (int_tumor_lipo + int_tumor_lipid + cell_tumor_lipo) * ug_per_g_tumor
    tumor_dna_dox <- cell_tumor_doxdna * ug_per_g_tumor
    tumor_cell_dox <- cell_tumor_dox * ug_per_g_tumor
    # Tumor concentrations (nM doxorubicin) in each space (Supplementary
    # Figure S2B); cell values are per unit cell volume.
    c_cell_tumor_free <- (cell_tumor_ol + cell_tumor_il + cell_tumor_cyt) / vcell_tumor
    c_cell_tumor_bound <- cell_tumor_doxdna / vcell_tumor
    # Liposome-equivalents delivered per tumor cell (TumorCellLiposomes_npc).
    tumor_cell_liposomes <- (cell_tumor_dox + cell_tumor_lipo + cell_tumor_lipid) / doxperliposome * nav / ncell_tumor
    # Heart readouts (ug doxorubicin per g heart), with the residual
    # capillary blood fraction (HeartTotalDox_ug_per_g,
    # HeartDNABoundDox_ug_per_g).
    heart_total_dox <- (int_heart_dox + int_heart_lipo + cell_heart_dox + cell_heart_lipo) * ug_per_g_heart +
      residual_blood_frxn * (cap_heart_dox + cap_heart_lipo) * ug_per_g_heart
    heart_dna_dox <- cell_heart_doxdna * ug_per_g_heart
  })
}
