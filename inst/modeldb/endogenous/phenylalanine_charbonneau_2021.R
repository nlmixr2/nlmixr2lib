phenylalanine_charbonneau_2021 <- function() {
  description <- "Phenylalanine model for absorption and metabolism in healthy subjects and patients with PKU"
  reference <- paste(
    "Charbonneau MR, Denney WS, Horvath NG, Cantarella P, Castillo MJ,",
    "Puurunen MK, Brennan AM. Development of a mechanistic model to predict",
    "synthetic biotic activity in healthy volunteers and patients with",
    "phenylketonuria. Commun Biol. 2021;4(1):898.",
    "doi:10.1038/s42003-021-02183-1.",
    "Extends the Kaufman 1999 phenylalanine metabolism model",
    "(Proc Natl Acad Sci USA. 1999;96(6):3160-3164. doi:10.1073/pnas.96.6.3160;",
    "correction Proc Natl Acad Sci USA. 2019;116:8632),",
    "available in this package as Kaufman_1999_phenylalanine.",
    sep = " "
  )
  vignette <- "phenylalanine_charbonneau_2021"
  # Time in hours and plasma Phe in mmol/L (Table 4); meal Phe is dosed into
  # the gut compartment in mg (Table 4, Gut = gut Phe content, mg).
  units <- list(time = "h", dosing = "mg", concentration = "mmol/L")

  # Dietary (meal) phenylalanine enters the gut compartment (Eq. 18).
  dosing <- "gut"

  population <- list(
    species = "human",
    disease_state = paste(
      "Healthy subjects and patients with phenylketonuria (PKU). Residual",
      "phenylalanine hydroxylase (PAH) activity is set with f_pah: 0, 0.01 and",
      "0.02 for classical PKU, 0.5 for a PKU heterozygote and 1 for a healthy",
      "subject (Charbonneau 2021 Results, p. 4-5). bl_phe is the initial plasma",
      "Phe (Table 4: 1.18 mmol/L, the classical-PKU steady state; Fig. 3d uses",
      "0.06 mmol/L for healthy subjects)."
    ),
    dose_range = paste(
      "Dietary Phe as meals into the gut compartment: 2.5 g/day (50 g protein)",
      "split evenly over meals at 0, 4 and 10 h of each day (Methods, p. 10);",
      "a single 1.0 g Phe load in healthy subjects (Results, p. 6; Fig. 3d)."
    ),
    weight_median = "70 kg (Table 4, W; fixed body weight of the simulations)",
    notes = paste(
      "Deterministic mechanistic model (no IIV, no residual error) that fits no",
      "data. The PAH, transaminase and net-protein-breakdown kinetics are",
      "carried from Kaufman 1999 (refs 29-30, with errata corrections); renal",
      "clearance is estimated from Kitagawa 1974 (ref 37); gut absorption is",
      "assumed complete and rapid. Fig. 3d compares the model against the",
      "placebo arm (n = 8 per study day) of the first-in-human SYNB1618 trial",
      "(SYNB1618-CP-001, NCT03516487). The authors' Python implementation is",
      "deposited on Zenodo (ref 38)."
    )
  )

  # What each ODE state holds, in what units, in what biological matrix.
  # Checked against Charbonneau 2021 Table 4 (p. 10) and Eqs. 18 and 20.
  compartmentData <- list(
    gut = list(analyte = "phenylalanine", units = "mg", specimen = "administration site", verified = TRUE),
    phe = list(analyte = "phenylalanine", units = "mmol/L", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Scales the volume of distribution (vd * WT) in the gut-to-plasma unit",
        "conversion f_gut_plasma (Eq. 19). Table 4 fixes W = 70 kg for every",
        "published simulation. Not a fitted covariate effect."
      ),
      source_name = "W"
    )
  )
  depends <- c("f_pah", "bl_phe")
  ini({
    # Initial values of the states (Table 4, Initial value column, p. 10)
    bl_phe <- 1.18; label("Baseline Phenylalanine (Phe) concentration (mmol/L)") # Table 4 (p. 10): Phe initial value = 1.18 mmol/L
    bl_gut <- 0; label("Baseline Phe in the gut (mg)") # Table 4 (p. 10): Gut initial value = 0 mg

    ka_gut <- 0.25; label("Absorption rate from gut to plasma (1/h)") # Table 4 (p. 10): K_a,gut = 0.25 1/hr
    v_npd <- 0.012; label("Rate of net protein breakdown ((mmol/L)/h)") # Table 4 (p. 10): V_npd = 0.012 (mmol/L)/hr

    vmax_pah <- 0.9; label("Maximum rate of Phe breakdown by PAH in a healthy subject ((mmol/L)/h)") # Table 4 (p. 10): V_max,PAH = 0.9 (mmol/L)/hr
    f_pah <- 0; label("Fraction of healthy PAH activity (PKU patient = 0 to 0.02)") # Table 4 (p. 10): F_PAH = 0.0-1.0, unitless; 0 = classical PKU (Results, p. 4)
    km_pah <- 0.51; label("Michaelis-Menten constant for Phe with PAH (mmol/L)") # Table 4 (p. 10): K_m,PAH = 0.51 mmol/L
    kact_pah <- 0.54; label("Phe activation constant for PAH (mmol/L)") # Table 4 (p. 10): K_a,PAH = 0.54 mmol/L

    vmax_trans <- 0.063; label("Maximum rate of Phe breakdown by transaminase ((mmol/L)/h)") # Table 4 (p. 10): V_max,trans = 0.063 (mmol/L)/hr
    km_trans <- 1.37; label("Michaelis-Menten constant for Phe with transaminase (mmol/L)") # Table 4 (p. 10): K_m,trans = 1.37 mmol/L

    lcl_renal <- log(5.696e-4); label("Renal clearance of Phe per body weight ((L/kg)/h)") # Table 4 (p. 10): CL_renal = 5.696e-4 (L/kg)/hr

    vd <- 0.5; label("Body-weight normalized volume distribution of Phe (L/kg)") # Table 4 (p. 10): V_d = 0.5 L/kg
  })
  model({
    # Molecular weight of Phe (g/mol); Table 4 (p. 10): MW_Phe = 165.19
    mw_phe <- 165.19
    # Unit conversion adjustment from Gut to Plasma concentrations (mmol/L)/mg; Eq. 19
    f_gut_plasma <- 1/(mw_phe * vd * WT)

    v_pah <- vmax_pah*f_pah / (1 + km_pah/phe + km_pah*kact_pah/(phe^2)) # Eq. 15; units: (mmol/L)/h
    v_trans <- vmax_trans / (1 + km_trans/phe) # Eq. 16; units: (mmol/L)/h
    cl_renal <- exp(lcl_renal)
    # Eq. 17 as printed: Phe x CL_renal x V_d. Dimensionally this carries an
    # extra (L/kg)^2 relative to (mmol/L)/h; it is kept as published, and the
    # paper's renal share of ~1% of classical-PKU elimination (Methods, p. 10)
    # is reproduced only with this form.
    v_renal <- phe * cl_renal * vd # units: (mmol/L)/h

    d/dt(gut) <- -ka_gut*gut # Eq. 18
    d/dt(phe) <- ka_gut*gut*f_gut_plasma + v_npd - v_pah - v_trans - v_renal # Eq. 20
    gut(0) <- bl_gut
    phe(0) <- bl_phe
    phe_umol <- phe * 1000 # units: umol/L (more commonly used in clinical laboratories)

    # The following is an augmentation of the model reported in the paper.  It
    # indicates the approximate daily Phe intake (in mg) to achieve
    # steady-state.
    daily_phe_intake <- 24 * (v_pah + v_trans + v_renal - v_npd) / f_gut_plasma
  })
}
