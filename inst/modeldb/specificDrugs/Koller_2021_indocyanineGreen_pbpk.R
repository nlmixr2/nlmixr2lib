Koller_2021_indocyanineGreen_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, 12 ODEs, SBML/libroadrunner). Indocyanine green",
    "(ICG) liver-function test in healthy adults, liver cirrhosis and",
    "partial hepatectomy, built by Koller, Grzegorzewski, Tautenhahn and",
    "Koenig (2021) from 29 curated clinical studies. An intravenous dose",
    "enters venous plasma and circulates through lung, arterial plasma,",
    "gastrointestinal tract, portal vein, liver and a lumped rest-of-body",
    "plasma pool by irreversible blood-flow transport. ICG is eliminated",
    "only by the liver: saturable (Michaelis-Menten) hepatic uptake with",
    "competitive bilirubin inhibition, saturable biliary export into a",
    "bile pool, and first-order transport from bile to feces, all scaled",
    "by the functional liver tissue volume. Cirrhosis severity (paper",
    "parameter f_cirrhosis, supplied here as HEPFUNC_REL = 1 -",
    "f_cirrhosis) removes the same fraction of parenchymal tissue and",
    "shunts the same fraction of arterial and portal blood past the liver;",
    "partial hepatectomy (resection_rate) removes a fraction of liver",
    "volume. As in the paper's simulations, the liver-plasma bilirubin",
    "that inhibits uptake is held at a fixed amount, so its concentration",
    "rises with (75 / WT) / (1 - resection_rate). The model is",
    "deterministic - the authors fit five hepatic",
    "transport parameters for a typical individual and report no",
    "between-subject variability or residual error."
  )
  reference <- paste(
    "Koller A, Grzegorzewski J, Tautenhahn HM, Koenig M. Prediction of",
    "Survival After Partial Hepatectomy Using a Physiologically Based",
    "Pharmacokinetic Model of Indocyanine Green Liver Function Tests.",
    "Front Physiol. 2021;12:730418. doi:10.3389/fphys.2021.730418.",
    "Fixed physiological, fitted and scan parameters from Table 2; model",
    "structure from Section 3.1 and Figure 1. The paper does not print",
    "rate laws; they were taken from the model archive the paper cites as",
    "the version used (Section 3.1): Koenig M, Koller A. Indocyanine green",
    "physiological based pharmacokinetics model (PBPK) 1.0.0. Zenodo.",
    "doi:10.5281/zenodo.5552405 (file models/icg_body_flat.xml). See the",
    "vignette Errata for two fractional-volume values where Table 2 and",
    "that archive disagree.",
    sep = " "
  )
  vignette <- "Koller_2021_indocyanineGreen_pbpk"
  units <- list(time = "min", dosing = "mg", concentration = "mg/L")

  # The dose enters `depot_iv`, which empties into venous plasma with a
  # half-life equal to the injection time (source SBML reaction iv_icg).
  dosing <- c("depot_iv")

  paper_specific_compartments <- c(
    "gut_plasma",
    "liver_plasma",
    "lung_plasma",
    "rest_plasma",
    "hepatic_vein",
    "bile"
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Scales every absolute organ volume (fractional volumes in L/kg",
        "times body weight) and, through cardiac output per body weight",
        "(COBW), every blood flow. Table 2 reference value BW = 75 kg;",
        "Section 3.2 states that study-specific body weights were used",
        "where reported and 75 kg otherwise. Body weight also scales the",
        "liver-plasma bilirubin concentration by 75 / WT, reproducing the",
        "fixed bilirubin amount of the paper's simulations (see the",
        "vignette Errata)."
      ),
      source_name = "BW"
    ),
    HEPFUNC_REL = list(
      description = "Relative liver function as a fraction of normal",
      units = "(dimensionless)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "1 = healthy liver. The source parameterises liver disease with the",
        "OPPOSITE orientation, as the cirrhosis degree f_cirrhosis on",
        "[0, 1), so this column is its complement: f_cirrhosis = 1 -",
        "HEPFUNC_REL. Section 3.1.4 couples the fraction of lost",
        "parenchymal tissue (f_tissue_loss) and the fraction of arterial",
        "and portal blood shunted past the liver (f_shunts) in lockstep to",
        "f_cirrhosis. Section 3.3 / Figure 5E map Child-Turcotte-Pugh",
        "classes onto f_cirrhosis as control = 0, mild (CTP-A) = 0.41,",
        "moderate (CTP-B) = 0.70 and severe (CTP-C) = 0.82, i.e.",
        "HEPFUNC_REL = 1, 0.59, 0.30 and 0.18. Section 2.6 Eq. 2 estimates",
        "an individual f_cirrhosis from a measured preoperative ICG-R15 as",
        "f_cirrhosis = 0.312 * ln(1.693 * R15) + 0.861. Must be > 0:",
        "HEPFUNC_REL = 0 removes all functional liver tissue and the",
        "hepatic transport terms divide by zero. Same parameterisation as",
        "Nemitz_2026_dapagliflozin_pbpk, from the same modelling group."
      ),
      source_name = "f_cirrhosis"
    )
  )

  compartmentData <- list(
    depot_iv = list(analyte = "indocyanine green", units = "mg", specimen = "administration site", verified = TRUE),
    venous = list(analyte = "indocyanine green", units = "mg", specimen = "plasma", verified = TRUE),
    lung_plasma = list(analyte = "indocyanine green", units = "mg", specimen = "plasma", verified = TRUE),
    arterial = list(analyte = "indocyanine green", units = "mg", specimen = "plasma", verified = TRUE),
    rest_plasma = list(analyte = "indocyanine green", units = "mg", specimen = "plasma", verified = TRUE),
    gut_plasma = list(analyte = "indocyanine green", units = "mg", specimen = "plasma", verified = TRUE),
    portal = list(analyte = "indocyanine green", units = "mg", specimen = "plasma", verified = TRUE),
    liver_plasma = list(analyte = "indocyanine green", units = "mg", specimen = "plasma", verified = TRUE),
    hepatic_vein = list(analyte = "indocyanine green", units = "mg", specimen = "plasma", verified = TRUE),
    liver = list(analyte = "indocyanine green", units = "mg", specimen = "tissue", verified = TRUE),
    bile = list(analyte = "indocyanine green", units = "mg", specimen = "bile", verified = TRUE),
    a_feces = list(analyte = "indocyanine green", units = "mg", specimen = "faeces", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 29L,
    age_range = "adults; not tabulated per study in the source",
    weight_range = paste(
      "study-specific mean body weight where reported (Table 1, 58.2-83.4",
      "kg); 75 kg reference individual otherwise (Section 3.2)"
    ),
    sex_female_pct = NA_real_,
    disease_state = paste(
      "Healthy volunteers (model calibration), patients with liver",
      "cirrhosis of Child-Turcotte-Pugh class A-C, and patients undergoing",
      "partial hepatectomy (validation). The survival classification used",
      "141 Japanese hepatectomy patients (109 survivors, 32",
      "non-survivors; Wakabayashi 2004, Seyama and Kokudo 2009)."
    ),
    dose_range = paste(
      "Intravenous bolus 0.25-5.0 mg/kg or 10-20 mg, and constant",
      "infusions 0.08-2.0 mg/min with or without a priming dose (Table 1)."
    ),
    regions = "international; the hepatectomy classification data are Japanese",
    notes = paste(
      "Deterministic typical-individual model. Section 4.2 states that",
      "all simulations used one identical parameter set, not accounting",
      "for inter-individual or inter-study differences, so no IIV and no",
      "residual-error model are reported and none are invented here.",
      "Data were curated into PK-DB (identifiers in Table 1). Five hepatic",
      "transport parameters were fitted by weighted least squares",
      "(250 local optimisation runs) to plasma, bile and extraction-ratio",
      "time courses from 13 studies in healthy subjects (Table 1 'Fit'",
      "column); all cirrhosis and hepatectomy data were validation only.",
      "The model targets the first 15-30 min of ICG elimination and does",
      "not capture the slower second elimination phase (Section 4.2)."
    )
  )

  ini({
    covbw            <- fixed(0.833333333333333); label("Cardiac output per unit body weight (mL/s/kg)")      # Table 2 COBW = 0.83 (SBML COBW = 0.833333333333333)
    f_cardiac_output <- fixed(1);                 label("Cardiac output scaling factor (unitless)")          # Table 2 scan parameter f_cardiac_output = 1
    f_bloodflow      <- fixed(1);                 label("Hepatic blood flow scaling factor (unitless)")      # Table 2 scan parameter f_bloodflow = 1
    hct              <- fixed(0.51);              label("Hematocrit (fraction)")                             # Table 2 HCT = 0.51
    fblood           <- fixed(0.02);              label("Blood-vessel fraction of organ volume (unitless)")  # Table 2 Fblood = 0.02

    fvgi <- fixed(0.0297);  label("Gastrointestinal tract fractional tissue volume (L/kg)")  # SBML FVgi = 0.0297; Table 2 prints 0.0171 - see vignette Errata
    fvbi <- fixed(0.00071); label("Bile fractional volume (L/kg)")                           # Table 2 FVbi = 0.00071
    fvli <- fixed(0.021);   label("Liver fractional tissue volume (L/kg)")                   # Table 2 FVli = 0.021
    fvlu <- fixed(0.0076);  label("Lung fractional tissue volume (L/kg)")                    # SBML FVlu = 0.0076; Table 2 prints 0.0297 - see vignette Errata
    fvve <- fixed(0.0587);  label("Venous blood fractional volume (L/kg)")                   # Table 2 FVve = 0.0587
    fvar <- fixed(0.0184);  label("Arterial blood fractional volume (L/kg)")                 # Table 2 FVar = 0.0184
    fvpo <- fixed(0.001);   label("Portal vein fractional volume (L/kg)")                    # Table 2 FVpo = 0.001
    fvhv <- fixed(0.001);   label("Hepatic vein fractional volume (L/kg)")                   # Table 2 FVhv = 0.001

    fqgi <- fixed(0.19);  label("Gastrointestinal fractional blood flow (unitless)")       # Table 2 FQgi = 0.19
    fqh  <- fixed(0.255); label("Hepatic venous-side fractional blood flow (unitless)")    # Table 2 FQh = 0.255
    fqlu <- fixed(1);     label("Lung fractional blood flow (unitless)")                   # Section 2.4 'fractional blood flow through the lung (must be 1)'; SBML FQlu = 1

    mr_icg <- fixed(774.96493); label("Molecular weight of indocyanine green (g/mol)")  # Table 2 Mr_icg = 774.96 (SBML 774.96493)
    ti_icg <- fixed(5);         label("Injection time of the intravenous dose (s)")     # Table 2 ti_icg = 5 s

    icgim_vmax   <- fixed(0.0369598840327503);  label("Vmax, hepatic ICG uptake (mmol/min/L liver tissue)")  # Table 2 LI__ICGIM_Vmax = 0.037 (fitted; SBML full precision)
    icgim_km     <- fixed(0.021659178617926);   label("Km, hepatic ICG uptake (mmol/L)")                     # Table 2 LI__ICGIM_Km = 0.0217 (fitted; SBML full precision)
    icgim_ki_bil <- fixed(0.02);                label("Bilirubin inhibition constant of hepatic uptake (mmol/L)")  # Table 2 LI__ICGIM_ki_bil = 0.02
    bil_ext      <- fixed(0.01);                label("Liver-plasma bilirubin concentration in the reference individual (mmol/L)")  # SBML species LI__bil_ext initial concentration 0.01 mM; Table 2 ki_bil comment gives reference 0.005-0.015 mmol/L
    wt_ref       <- fixed(75);                  label("Reference body weight at which bil_ext applies (kg)")       # Table 2 BW = 75 kg
    f_oatp1b3    <- fixed(1);                   label("Hepatic uptake transporter amount scaling factor (unitless)")  # Table 2 scan parameter LI__f_oatp1b3 = 1

    icgli2ca_vmax <- fixed(0.000943672769975891); label("Vmax, biliary export of ICG (mmol/min/L liver tissue)")  # Table 2 LI__ICGLI2CA_Vmax = 0.000944 (fitted; SBML full precision)
    icgli2ca_km   <- fixed(0.012388659243625);    label("Km, biliary export of ICG (mmol/L)")                      # Table 2 LI__ICGLI2CA_km = 0.0124 (fitted; SBML full precision)
    icgli2bi_k    <- fixed(0.000114596604507925); label("Rate constant, bile-to-feces ICG transport (1/min)")      # Table 2 LI__ICGLI2BI_Vmax = 0.000114 1/min (fitted; SBML full precision)

    resection_rate <- fixed(0); label("Fraction of liver volume resected (unitless)")  # Table 2 scan parameter resection_rate = 0; Section 3.1.5 varies it up to 0.9

    propSd <- fixed(0); label("Proportional residual error (fraction; ZERO - not reported by the source)")  # Section 4.2: one identical deterministic parameter set
  })

  model({
    # ================= Liver disease (Section 3.1.4) =====================
    # f_shunts and f_tissue_loss move in lockstep with f_cirrhosis.
    f_cirrhosis   <- 1 - HEPFUNC_REL
    f_shunts      <- f_cirrhosis
    f_tissue_loss <- f_cirrhosis

    # ================= Organ volumes (L) =================================
    fvre <- 1 - (fvbi + fvgi + fvli + fvlu + fvve + fvar + fvpo + fvhv)

    vgi <- WT * fvgi
    vbi <- WT * fvbi
    # Partial hepatectomy removes a fraction of the liver volume (Section 3.1.5).
    vli <- WT * fvli * (1 - resection_rate)
    vlu <- WT * fvlu
    vre <- WT * fvre

    # Large-vessel plasma volumes: each vessel's share of total blood volume,
    # net of the blood already counted inside organs, times (1 - hematocrit).
    fvvessel <- fvar + fvve + fvpo + fvhv
    vblood_organ <- WT * fblood * (1 - fvvessel)
    vve <- (1 - hct) * (WT * fvve - (fvve / fvvessel) * vblood_organ)
    v_ar <- (1 - hct) * (WT * fvar - (fvar / fvvessel) * vblood_organ)
    vpo <- (1 - hct) * (WT * fvpo - (fvpo / fvvessel) * vblood_organ)
    vhv <- (1 - hct) * (WT * fvhv - (fvhv / fvvessel) * vblood_organ)

    vre_plasma <- vre * fblood * (1 - hct)
    vgi_plasma <- vgi * fblood * (1 - hct)
    vli_plasma <- vli * fblood * (1 - hct)
    vlu_plasma <- vlu * fblood * (1 - hct)
    # Functional (parenchymal) liver tissue shrinks with cirrhosis.
    vli_tissue <- vli * (1 - f_tissue_loss) * (1 - fblood)

    # ================= Blood flows (L/min) ===============================
    co <- WT * covbw * f_cardiac_output
    qc <- (co / 1000) * 60
    qlu <- qc * fqlu
    qh <- qc * fqh * f_bloodflow
    qgi <- qc * fqgi * f_bloodflow
    qpo <- qgi
    qha <- qh - qpo
    fqre <- 1 - qh / qlu
    qre <- qc * fqre

    # ================= Concentrations (mg/L) =============================
    c_ve <- venous / vve
    c_ar <- arterial / v_ar
    c_lu <- lung_plasma / vlu_plasma
    c_re <- rest_plasma / vre_plasma
    c_gi <- gut_plasma / vgi_plasma
    c_po <- portal / vpo
    c_li <- liver_plasma / vli_plasma
    c_hv <- hepatic_vein / vhv

    # ================= Intravenous input (mg/min) ========================
    # First-order release with a half-life equal to the injection time.
    ki_icg <- (0.693 / ti_icg) * 60
    iv_in <- ki_icg * depot_iv

    # ================= Blood-flow transport (mg/min) =====================
    flow_ar_re <- qre * c_ar
    flow_re_ve <- qre * c_re
    flow_ar_gi <- qgi * c_ar
    flow_gi_po <- qgi * c_gi
    flow_ar_li <- (1 - f_shunts) * qha * c_ar
    flow_ar_hv <- f_shunts * qha * c_ar
    flow_po_li <- (1 - f_shunts) * qpo * c_po
    flow_po_hv <- f_shunts * qpo * c_po
    flow_li_hv <- (1 - f_shunts) * (qpo + qha) * c_li
    flow_hv_ve <- qh * c_hv
    flow_ve_lu <- qlu * c_ve
    flow_lu_ar <- qlu * c_lu

    # ================= Hepatic transport (mg/min) ========================
    # Rate laws are written in mmol/L and mmol/min as in the source SBML;
    # dividing a mg/L concentration by mr_icg gives mmol/L, and multiplying
    # a mmol/min rate by mr_icg gives mg/min.
    # Bilirubin is a reaction-free species in liver plasma. In the paper's
    # simulations its AMOUNT stayed at the reference individual's value
    # (0.01 mmol/L x liver-plasma volume at 75 kg, no resection) when body
    # weight or resection rate changed, so its concentration rises as the
    # liver-plasma volume shrinks. This as-run behaviour reproduces Figure
    # 7A-D; see the vignette Errata.
    vli_plasma_ref <- wt_ref * fvli * fblood * (1 - hct)
    bil_li <- bil_ext * vli_plasma_ref / vli_plasma

    c_li_mm <- c_li / mr_icg
    c_liver_mm <- liver / vli_tissue / mr_icg
    c_bile_mm <- bile / vbi / mr_icg

    uptake <- mr_icg * f_oatp1b3 * icgim_vmax * vli_tissue * c_li_mm /
      (icgim_km * (1 + bil_li / icgim_ki_bil) + c_li_mm)
    biliary_export <- mr_icg * icgli2ca_vmax * vli_tissue * c_liver_mm /
      (icgli2ca_km + c_liver_mm)
    bile_to_feces <- mr_icg * icgli2bi_k * vli_tissue * c_bile_mm

    # ================= ODEs (amounts, mg) ================================
    d/dt(depot_iv) <- -iv_in
    d/dt(venous) <- iv_in + flow_hv_ve + flow_re_ve - flow_ve_lu
    d/dt(lung_plasma) <- flow_ve_lu - flow_lu_ar
    d/dt(arterial) <- flow_lu_ar - flow_ar_re - flow_ar_gi - flow_ar_li - flow_ar_hv
    d/dt(rest_plasma) <- flow_ar_re - flow_re_ve
    d/dt(gut_plasma) <- flow_ar_gi - flow_gi_po
    d/dt(portal) <- flow_gi_po - flow_po_li - flow_po_hv
    d/dt(liver_plasma) <- flow_ar_li + flow_po_li - flow_li_hv - uptake
    d/dt(hepatic_vein) <- flow_ar_hv + flow_po_hv + flow_li_hv - flow_hv_ve
    d/dt(liver) <- uptake - biliary_export
    d/dt(bile) <- biliary_export - bile_to_feces
    d/dt(a_feces) <- bile_to_feces

    # ================= Observations ======================================
    # Venous plasma ICG (mg/L), the quantity from which the paper computes
    # ICG-R15, ICG-PDR, clearance and half-life (Section 3.1).
    Cc <- c_ve
    C_arterial <- c_ar
    C_hepatic_vein <- c_hv
    # Hepatic extraction ratio with the source's 1e-7 mmol/L guard (SBML ER_icg).
    er_guard <- 1e-7 * mr_icg
    ER <- ((c_ar + er_guard) - (c_hv + er_guard)) / (c_ar + er_guard)

    Cc ~ prop(propSd)
  })
}
