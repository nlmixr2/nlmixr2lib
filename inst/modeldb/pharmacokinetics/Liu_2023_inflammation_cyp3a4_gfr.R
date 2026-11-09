Liu_2023_inflammation_cyp3a4_gfr <- function() {
  description <- paste(
    "Physiologically based clearance model (hypothetical drug, one",
    "compartment). Effect of acute inflammation, represented by the plasma",
    "interleukin-6 (IL-6) concentration, on the clearance and exposure of a",
    "hypothetical small-molecule drug eliminated by hepatic CYP3A4",
    "metabolism, by renal glomerular filtration, or by both. IL-6 lowers",
    "hepatic CYP3A4 activity through a sigmoid Imax function fitted to",
    "digitised human-hepatocyte data (Imax 73.26 percent, IC50 43.65 pg/mL,",
    "Hill 0.59) and raises GFR linearly (0.01556 mL/min/1.73 m^2 per pg/mL)",
    "as fitted to a healthy-volunteer endotoxaemia study. Hepatic clearance",
    "follows the axial-dispersion liver model (dispersion number 0.27) with",
    "hepatic blood flow 25 percent of a 5500 mL/min cardiac output; renal",
    "clearance is fu times GFR. As in the authors' code, the drug's",
    "clearance in the absence of inflammation is a user input (unit value",
    "by default) and inflammation scales it by the ratio of the",
    "physiological clearance at the current IL-6 to that at IL-6 = 0, which",
    "drives a one-compartment model with optional first-order absorption.",
    "The drug-specific properties (unbound fraction, intrinsic clearance,",
    "elimination pathway) are placeholders that the user replaces.",
    "Deterministic typical-value model: no between-subject variability and",
    "no residual error.",
    sep = " "
  )
  reference <- paste(
    "Liu F, Aulin LBS, Manson ML, Krekels EHJ, van Hasselt JGC. Unraveling",
    "the Effects of Acute Inflammation on Pharmacokinetics: A Model-Based",
    "Analysis Focusing on Renal Glomerular Filtration Rate and Cytochrome",
    "P450 3A4-Mediated Metabolism. Eur J Drug Metab Pharmacokinet.",
    "2023;48(6):623-631. doi:10.1007/s13318-023-00852-6. PMID 37715056.",
    "PMCID PMC10624742. Equations 1-13 and Tables 1-2 of the main text; the",
    "supplementary R code (Supplementary Information) gives the unrounded",
    "GFR slope, the fraction metabolised by CYP3A4 and the clearance",
    "normalisation used for every published figure.",
    sep = " "
  )
  vignette <- "Liu_2023_inflammation_cyp3a4_gfr"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    IL6 = list(
      description = "Plasma interleukin-6 concentration, the inflammation marker",
      units = "pg/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "IL-6 = 0 is the no-inflammation reference condition against which the clearance ratio is formed. The paper evaluated constant IL-6 levels of 0, 5, 10, 15, 50, 100, 500, 1000, 1500 and 2500 pg/mL (supplementary R code), spanning the clinically relevant range 0-2500 pg/mL (Table 1). The effect is instantaneous (no enzyme-turnover delay), so a time-varying IL-6 column is applied at the time it is recorded.",
      source_name = "IL6"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Scales the per-kg in vitro intrinsic clearance to a whole-liver value, CLint = WT * clint (supplementary R code `CLint0 <- BW * CLint_bw`). The paper used 70 kg (Table 1). Hepatic blood flow and GFR are not weight-scaled.",
      source_name = "BW"
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "De-normalises the BSA-indexed GFR, GFR = GFR_bsa * BSA / 1.73 (supplementary R code). The paper used 1.73 m^2, at which the factor is 1.",
      source_name = "BSA"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "hypothetical drug",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "hypothetical drug",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA,
    n_studies = NA,
    age_range = "not applicable (single 70 kg, 1.73 m^2 reference adult)",
    weight_range = "70 kg (reference value, Table 1)",
    disease_state = "acute inflammation (infection, sepsis) represented by plasma IL-6 0-2500 pg/mL",
    dose_range = "unit dose of a hypothetical drug, intravenous bolus or oral",
    notes = "No patients were studied. The IL-6 to CYP3A4 relationship was fitted to in vitro human hepatocyte data digitised from two published studies (refs 19, 20 of the paper), and the IL-6 to GFR relationship to a healthy-volunteer endotoxaemia study (ref 12). Simulations covered a panel of hypothetical drugs with in vitro intrinsic clearance 0.5-189 mL/min/kg (12 equally spaced values), unbound fraction 0.01, 0.25, 0.5, 0.75 and 1, and renal, hybrid or hepatic elimination (Table 1)."
  )

  ini({
    # ------------------------------------------------------------------
    # IL-6 effect on CYP3A4 activity (Eq. 1) and on GFR (Eq. 2). These four
    # values were estimated by the authors (Table 2, estimate and SE).
    # ------------------------------------------------------------------
    imax <- 0.7326; label("Maximal fractional inhibition of CYP3A4 activity by IL-6 (fraction)") # Table 2: Imax 73.26 (SE 0.30) percent; supplementary R code Imax <- 73.26 used as Imax/100
    ic50 <- 43.65; label("IL-6 concentration giving half-maximal CYP3A4 inhibition (pg/mL)") # Table 2: IC50 43.65 (SE 1.11) pg/mL
    hill <- 0.59; label("Hill coefficient of the IL-6 effect on CYP3A4 activity (unitless)") # Table 2: gamma 0.59 (SE 0.01)
    e_il6_gfr <- 0.01556; label("Increase in BSA-normalised GFR per unit IL-6 ((mL/min/1.73 m^2)/(pg/mL))") # supplementary R code slope_dGFR <- 0.01556; Table 2 prints alpha 0.016 (SE 0.001) to 2 significant figures

    # ------------------------------------------------------------------
    # System-specific physiology (Table 1), all literature values.
    # ------------------------------------------------------------------
    qc <- fixed(5500); label("Cardiac output (mL/min)") # Table 1: Qco 5500 mL/min; supplementary R code Qco <- 5.5 * 1000
    fq_liver <- fixed(0.25); label("Hepatic blood flow as a fraction of cardiac output (fraction)") # Methods: QH nearly 25 percent of Qco; supplementary R code Qh <- 0.25 * Qco; Table 1: QH 1375 mL/min
    gfr <- fixed(113); label("GFR without inflammation (mL/min/1.73 m^2)") # Table 1: GFRbase 113 mL/min/1.73 m2 (iohexol-derived)
    bpr <- fixed(1); label("Blood-to-plasma concentration ratio (unitless)") # Table 1: BP 1
    dn <- fixed(0.27); label("Axial dispersion number of the hepatic dispersion model (unitless)") # Table 1 and Methods: DN 0.27 for microsomal data

    # ------------------------------------------------------------------
    # Drug-specific properties of the hypothetical drug. The paper varies
    # these across simulation scenarios (Table 1); the defaults are one
    # scenario and are the values a user replaces for a real compound.
    # ------------------------------------------------------------------
    fu <- fixed(0.5); label("Fraction unbound in plasma (fraction)") # Table 1: fu 0.01, 0.25, 0.5, 0.75 or 1; 0.5 is the Fig. 2 scenario
    clint <- fixed(54.4); label("In vitro unbound intrinsic clearance per kg body weight (mL/min/kg)") # Table 1: CLint 0.5-189 mL/min/kg; 54.4 is the Fig. 3 intermediate-ER drug
    fm_cyp3a4 <- fixed(1); label("Fraction of the intrinsic clearance mediated by CYP3A4 (fraction)") # supplementary R code fm_A <- 1
    fcl_h <- fixed(0.5); label("Elimination-pathway index: 0 renal only, 0.5 hepatic plus renal, 1 hepatic only (unitless)") # Eq. 7 and Table 1: fcl_h 0, 0.5 or 1

    # ------------------------------------------------------------------
    # One-compartment PK of the hypothetical drug (Eqs. 10-13, Table 1).
    # The authors normalised the no-inflammation total clearance to 1
    # (supplementary R code: 'total body clearance assumed to be 1',
    # AUC0 <- 1, theta tCL = CL / CL0), so lcl is that reference clearance.
    # ------------------------------------------------------------------
    lcl <- fixed(log(1)); label("Total clearance in the absence of inflammation (L/h)") # supplementary R code AUC0 <- 1 with Dose 1 and F 1, so CL0 = 1
    lvc <- fixed(log(1)); label("Volume of distribution (L)") # Table 1: Vd 1
    lka <- fixed(log(1)); label("First-order absorption rate constant (1/h)") # Table 1: ka 1
    lfdepot <- fixed(log(1)); label("Oral bioavailability (fraction)") # Table 1 and Eq. 9: F 1
  })

  model({
    # CYP3A4 activity relative to no inflammation, Eq. 1 written as in the
    # supplementary R code, (100 - Imax * IL6^h / (IC50^h + IL6^h)) / 100,
    # which is algebraically Eq. 1 and stays finite at IL-6 = 0.
    cyp3a4_act <- 1 - imax * IL6^hill / (ic50^hill + IL6^hill)

    # Hepatic plasma clearance, dispersion model (Eqs. 3-6). Flows and
    # clearances are in mL/min. Only the CYP3A4 share of the intrinsic
    # clearance responds to IL-6. The `0` suffix denotes IL-6 = 0.
    qh <- fq_liver * qc
    clint0 <- clint * WT
    clint_il6 <- clint0 * (fm_cyp3a4 * cyp3a4_act + 1 - fm_cyp3a4)

    rn <- fu / bpr * clint_il6 / qh
    a_dn <- sqrt(1 + 4 * dn * rn)
    eh <- 1 - 4 * a_dn / ((1 + a_dn)^2 * exp((a_dn - 1) / (2 * dn)) - (1 - a_dn)^2 * exp(-(a_dn + 1) / (2 * dn)))
    cl_hep <- qh * eh * bpr

    rn0 <- fu / bpr * clint0 / qh
    a_dn0 <- sqrt(1 + 4 * dn * rn0)
    eh0 <- 1 - 4 * a_dn0 / ((1 + a_dn0)^2 * exp((a_dn0 - 1) / (2 * dn)) - (1 - a_dn0)^2 * exp(-(a_dn0 + 1) / (2 * dn)))
    cl_hep0 <- qh * eh0 * bpr

    # Renal plasma clearance of unbound drug by glomerular filtration
    # (Eqs. 2 and 8), GFR de-normalised from 1.73 m^2 to the actual BSA.
    gfr_il6 <- (gfr + e_il6_gfr * IL6) * BSA / 1.73
    cl_renal <- fu * gfr_il6
    cl_renal0 <- fu * gfr * BSA / 1.73

    # Eq. 7: the hepatic arm is present for fcl_h > 0 and the renal arm for
    # fcl_h < 1, so 1 = hepatic only, 0 = renal only, 0.5 = their sum.
    hep_on <- 0
    if (fcl_h > 0) hep_on <- 1
    renal_on <- 0
    if (fcl_h < 1) renal_on <- 1
    cl_ratio <- (hep_on * cl_hep + renal_on * cl_renal) / (hep_on * cl_hep0 + renal_on * cl_renal0)

    # Eqs. 9 and 14: AUC = F * Dose / CL, so the AUC under inflammation
    # relative to no inflammation is the inverse clearance ratio.
    auc_ratio <- 1 / cl_ratio

    cl <- exp(lcl) * cl_ratio
    vc <- exp(lvc)
    ka <- exp(lka)
    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    f(depot) <- exp(lfdepot)

    Cc <- central / vc
  })
}
