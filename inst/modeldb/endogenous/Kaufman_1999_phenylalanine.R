Kaufman_1999_phenylalanine <- function() {
  description <- "Mechanistic model of human phenylalanine metabolism (PAH hydroxylation with substrate activation, transaminase, net protein degradation) in normal subjects, PKU heterozygotes, and phenylketonuric patients"
  reference   <- "Kaufman S. A model of human phenylalanine metabolism in normal subjects and in phenylketonuric patients. Proc Natl Acad Sci USA. 1999 Mar 16;96(6):3160-3164. doi:10.1073/pnas.96.6.3160"
  vignette    <- "Kaufman_1999_phenylalanine"
  units       <- list(time = "h", dosing = "none", concentration = "mmol/L")

  population <- list(
    species        = "human",
    disease_state  = "Applicable to (i) normal control subjects, (ii) obligate PKU heterozygotes (parents of PKU patients), and (iii) classical, moderate, and mild PKU patients. Residual PAH activity is encoded via f_pah (0 = classical PKU with zero PAH; 1 = wild-type PAH; ~0.4 typical for obligate PKU heterozygotes per Kaufman 1999 Table 1). Baseline blood phenylalanine bl_phe (default 0.058 mM from Scriver 1985 mean control level) is the initial condition; set higher (e.g., 1 mM) to reproduce the Kaufman 1999 Fig. 2 phenylalanine-load clearance curves.",
    notes          = "Deterministic mechanistic model; Kaufman 1999 fits no data. The six kinetic constants (vmax_pah, km_pah, kact_pah, vmax_trans, km_trans, v_npd) are derived from external literature or from Kowlessur and Kaufman unpublished measurements as cited in the paper's Methods; f_pah and bl_phe are scenario knobs varied to represent different physiological states (see Kaufman 1999 Table 1). No IIV; no residual error; no dosing events (the endogenous Phe pool starts at bl_phe). Body weight enters only the paper's Methods calculations (a 500 mL/kg volume-of-distribution assumption for translating urinary metabolite excretion into a plasma-rate; see Kaufman 1999 p.3162) and is NOT a covariate on the run-time model. Km_TRANS = 1.37 mM is Kaufman's calculation from Guldberg 1995 (ref 21) load-test data with reported mean +- SD = 1.37 +- 0.14 mM (n = 3); the point estimate is used here."
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    phe = list(analyte = "phenylalanine", units = NA_character_, specimen = "plasma", verified = FALSE)
  )

  ini({
    vmax_pah <- fixed(0.9);   label("Maximum rate of Phe hydroxylation by PAH (mmol/L/h)")   # Kaufman 1999 p.3161 (ref 20: Guettler & Hansen 1977, Scand J Clin Lab Invest 37, 717-722); Kaufman notes this value is approximate ("probably an underestimate")
    km_pah   <- fixed(0.51);  label("Michaelis constant for Phe on PAH (mM)")                # Kaufman 1999 p.3161 (Kowlessur & Kaufman, unpublished data)
    kact_pah <- fixed(0.54);  label("PAH activation constant Ka for Phe (mM)")               # Kaufman 1999 p.3161 (Kowlessur & Kaufman, unpublished data)

    vmax_trans <- fixed(0.063); label("Maximum rate of Phe transamination (mmol/L/h)")       # Kaufman 1999 p.3162 (0.043 from Langenbeck 1992 ref 22 transaminate-metabolite excretion + 0.020 mmol/L/h PAG contribution from Knox 1972 ref 28)
    km_trans   <- fixed(1.37);  label("Michaelis constant for Phe on transaminase (mM)")     # Kaufman 1999 p.3162 (Kaufman calculated 1.37 +- 0.14 mM, n=3, from Guldberg 1995 ref 21 load-test data)

    v_npd <- fixed(0.012); label("Zero-order rate of net protein degradation releasing Phe (mmol/L/h)") # Kaufman 1999 p.3162 (Waterlow & Jackson 1981 ref 31 net protein breakdown + Clowes 1980 ref 32 muscle Phe content)

    f_pah  <- 1.0;   label("Fraction of wild-type PAH activity (0 = classical PKU; 1 = healthy control; ~0.4 typical for obligate PKU heterozygotes per Kaufman 1999 Table 1)")  # Kaufman 1999 Table 1 (scenario knob)
    bl_phe <- 0.058; label("Baseline blood phenylalanine (mM)")                              # Kaufman 1999 p.3162, ref 35 (Scriver 1985 accepted control level 0.058 +- 0.015 mM)
  })

  model({
    # Eq. 2 (Kaufman 1999 p.3161): PAH rate law -- two-site ordered-binding with substrate activation.
    #   v_PAH = Vmax_PAH / (1 + Km/[Phe] + Km*Ka/[Phe]^2)
    # f_pah in [0, 1] scales Vmax_PAH to represent residual PAH activity.
    v_pah <- vmax_pah * f_pah / (1 + km_pah/phe + km_pah*kact_pah/(phe^2))

    # Michaelis-Menten transaminase rate: v_TRANS = Vmax_TRANS / (1 + Km_TRANS/[Phe])
    v_trans <- vmax_trans / (1 + km_trans/phe)

    # Eq. 1 (Kaufman 1999 p.3161): d[Phe]/dt = -v_PAH - v_TRANS + v_NPD
    d/dt(phe) <- -v_pah - v_trans + v_npd
    phe(0)    <- bl_phe

    # Convenience augmentation: micromolar concentration (units commonly used in clinical laboratories)
    phe_umol <- phe * 1000
  })
}
