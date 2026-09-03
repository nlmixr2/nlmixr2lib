McKnite_2026_midazolam_crrt <- function() {
  description <- "Ex vivo (extracorporeal continuous renal replacement therapy circuit). Mechanistic six-state model of midazolam disposition in a Baxter PrisMax CRRT circuit fitted to closed-loop ex vivo data: a blood reservoir split into plasma and red-cell sub-compartments feeds a hemofilter that is likewise split into plasma and red-cell sub-compartments, which exchanges drug with the filter dialysate by diffusion across the membrane and by convection with the filtrate, and the dialysate drains to effluent waste. Circuit adsorption (plus any filter clearance and degradation not otherwise accounted for) is a saturable Michaelis-Menten sink acting on the filter plasma and red-cell states. Flow and hematocrit derivations are Equations 1-6 and the mass transfers Equations 7-11 of Supplement S1. The paper's whole-body child CRRT-PBPK model, which grafts this circuit onto a 17-compartment PK-Sim/MoBi pediatric model, is NOT reproduced here; only the transferable circuit module is."
  reference <- paste(
    "McKnite AM, Hamadeh A, Hunt JP, Green DJ, Imburgia C, Whelan A, Hudson R,",
    "Chevalier A, Mathis CL, Yang MJ, Dwyer JP, Edginton A, Watt KM.",
    "Midazolam Dosing During CRRT: A Combined Ex Vivo and Physiologically-Based",
    "Pharmacokinetic Approach. CPT Pharmacometrics Syst Pharmacol. 2026;15:e70240.",
    "doi:10.1002/psp4.70240. PMCID: PMC13274749.",
    "Circuit structure and the ex vivo prescription: Methods 2.5.",
    "Optimized Michaelis-Menten constants: Results 3.3.",
    "Flow / hematocrit derivations (Equations 1-6) and mass transfers",
    "(Equations 7-11): Supporting Information Appendix S1.",
    "Membrane diffusion coefficient, the unrounded Vmax, the plasma / red-cell",
    "exchange terms and the red-cell partition coefficient come from the authors'",
    "deposited MoBi 10.0 project (Appendix S1, 'PSP-2025-0082-s02'); see the",
    "per-parameter comments and the vignette source-trace table.",
    sep = " "
  )
  vignette <- "McKnite_2026_midazolam_crrt"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # The dosing target is neither `depot` nor `central`: in the ex vivo
  # experiment midazolam is spiked directly into the blood reservoir, so a dose
  # lands in the reservoir plasma and red-cell states. When the module is
  # attached to a body model instead, the reservoir states are replaced by that
  # model's venous blood and the dose goes wherever that model doses.
  dosing <- c("reservoir_plasma", "reservoir_rbc")

  compartmentData <- list(
    reservoir_plasma = list(analyte = "midazolam", units = "mg", specimen = "plasma", verified = TRUE),
    reservoir_rbc = list(analyte = "midazolam", units = "mg", specimen = "blood cell", verified = TRUE),
    circuit_plasma = list(analyte = "midazolam", units = "mg", specimen = "plasma", verified = TRUE),
    circuit_rbc = list(analyte = "midazolam", units = "mg", specimen = "blood cell", verified = TRUE),
    circuit_dialysate = list(analyte = "midazolam", units = "mg", specimen = "dialysate", verified = TRUE),
    circuit_effluent = list(analyte = "midazolam", units = "mg", specimen = "dialysate", verified = TRUE)
  )

  # The ex vivo circuit model carries no subject-level covariates: every
  # quantity that would be a covariate in the patient model (hematocrit,
  # fraction unbound, the CRRT prescription, the filter geometry) is an
  # experimental setting of the circuit and is exposed as a fixed ini() entry
  # so a user can re-point the module with rxSolve(params = ...).
  covariatesDataExcluded <- list(
    HCT = list(
      description = "Hematocrit of the blood circulating through the circuit",
      units = "fraction",
      type = "continuous",
      notes = paste(
        "A per-patient covariate in the paper's whole-body CRRT-PBPK model (critical illness",
        "reduced hematocrit by 28% relative to healthy age-matched children, Methods 2.9), but a",
        "fixed property of the ex vivo apparatus here: the three replicate ex vivo circuits were",
        "run at hematocrit 0.30 (Methods 2.5). Carried as the `hct` ini() parameter."
      )
    ),
    ALB = list(
      description = "Serum albumin, which sets the midazolam fraction unbound through the PK-Sim plasma protein scaling factor",
      units = "g/dL",
      type = "continuous",
      notes = paste(
        "A per-patient covariate in the whole-body model (Equation 1, plasma protein scaling",
        "factor; critical illness reduced albumin by 37%, Methods 2.9). The ex vivo circuits were",
        "run at a deliberately lowered albumin concentration chosen to reproduce the bound /",
        "unbound ratio of the target population (Discussion), so for this module albumin enters",
        "only through the fixed `fu` ini() parameter."
      )
    )
  )

  population <- list(
    species = "ex vivo (closed-loop extracorporeal CRRT circuit primed with whole blood)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    system = paste(
      "Baxter PrisMax continuous renal replacement therapy platform running an HF1000 hemofilter",
      "(Baxter Healthcare, Deerfield, IL) in a closed loop against a 600 mL blood reservoir."
    ),
    n_replicates = 3L,
    hematocrit = "0.30 (mean of three ex vivo replicates)",
    modality = "continuous venovenous hemodiafiltration (CVVHDF)",
    prescription = paste(
      "Blood flow 80 mL/min; dialysate flow 400 mL/h; pre-blood pump 300 mL/h;",
      "replacement fluid 100 mL/h; patient fluid removal 0 mL/h."
    ),
    dose_range = "single spike to a circuit concentration of 0.15 mg/L",
    disease_state = "not applicable (ex vivo)",
    notes = paste(
      "Methods 2.5 and Figure 1A. The ex vivo concentration data the Michaelis-Menten constants",
      "were optimized against were generated in a previously published ex vivo study by the same",
      "group (reference 13 of the source paper) and are not reproduced in this paper; the blood",
      "source and albumin concentration of those circuits are described there, not here.",
      "The paper's downstream patient work -- three children aged 12, 22 and 72 months on CVVHD or",
      "CVVHDF, and five 1000-subject age-stratified virtual populations from neonates to",
      "adolescents (Tables S2, S4, S5) -- exercises the whole-body PK-Sim model that this module",
      "was grafted onto, and is out of scope for this file."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # CRRT prescription (Methods 2.5, the ex vivo circuit). Flows are
    # carried in L/h to match units$time; the paper prints blood flow in
    # mL/min and the rest in mL/h.
    # ------------------------------------------------------------------
    qbfr <- fixed(4.8); label("Blood flow rate through the circuit, Q_BFR (L/h)")               # Methods 2.5: "blood flow rate, 80 mL/min" = 4.8 L/h
    qdia <- fixed(0.4); label("Dialysate flow rate, Q_DIA (L/h)")                               # Methods 2.5: "dialysate flow rate, 400 mL/min". The printed "mL/min" is a unit typo -- 400 mL/min = 24 L/h exceeds the PrisMax maximum dialysate rate and would be 5x the blood flow, every other flow in the same sentence is mL/h, and the deposited MoBi project stores the patient Q_DIA at mL/h scale. Read as 400 mL/h = 0.4 L/h; see the vignette Errata
    qpbp <- fixed(0.3); label("Pre-blood pump flow rate, Q_PBP (L/h)")                          # Methods 2.5: "pre-blood pump, 300 mL/h"
    qrep <- fixed(0.1); label("Total replacement fluid flow rate, Q_REP (L/h)")                 # Methods 2.5: "replacement fluids, 100 mL/h"
    qpfr <- fixed(0); label("Patient fluid removal rate, Q_PFR (L/h)")                          # Methods 2.5: "patient fluid removal, 0 mL/h"
    hct <- fixed(0.30); label("Hematocrit of the circulating blood (fraction)")                 # Methods 2.5: "a reservoir volume of 600 mL and hematocrit of 0.30 were used, representing the average values obtained across three replicates of ex vivo experiments"

    # ------------------------------------------------------------------
    # Circuit volumes (Methods 2.5, HF1000 hemofilter)
    # ------------------------------------------------------------------
    vreservoir <- fixed(0.600); label("Blood reservoir volume (L)")                             # Methods 2.5: "a reservoir volume of 600 mL"
    vfilter <- fixed(0.165); label("Hemofilter plus circuit blood volume (L)")                  # Methods 2.5: "filter volume 165 mL ... The filter volume represented both the HF1000 filter and circuit prime volumes to account for the additional volume of the circuit". Consistent with Table S3, which lists 165 mL as the HF1000 circuit blood volume
    vdialysate <- fixed(0.0591); label("Hemofilter dialysate-side volume (L)")                  # Methods 2.5: "dialysate volume 59.1 mL ... The volume of dialysate was equivalent to HF1000 filter prime volume"

    # ------------------------------------------------------------------
    # Hemofilter membrane geometry and transport (Methods 2.5; the
    # diffusion coefficient is an optimized value the paper does not print
    # and is recovered from the deposited MoBi project). Lengths are in dm
    # and areas in dm^2 so that dfilter * safilter / thfilter is L/h.
    # ------------------------------------------------------------------
    safilter <- fixed(120); label("Hemofilter membrane surface area, SA (dm^2)")                # Methods 2.5: "surface area 1.2 m2" = 120 dm^2. Table S3 lists the HF1000 membrane surface area as 1.1 m2; the Methods value for the ex vivo build is used here and the discrepancy is recorded in the vignette Errata. The term it enters is ~5 orders of magnitude faster than every other circuit process, so neither reading is distinguishable in a simulation
    thfilter <- fixed(5e-4); label("Hemofilter membrane thickness, h (dm)")                     # Methods 2.5: "membrane thickness 50 um" = 5e-4 dm. Corroborated by the deposited MoBi project (MembraneThickness stored as 0.0005 in the dm base unit); Table S3 lists 50 um fiber wall thickness for every filter type
    dfilter <- fixed(0.594); label("Membrane diffusion coefficient, D (dm^2/h)")                # NOT PRINTED IN THE PAPER. Methods 2.5 states the membrane diffusion coefficient was optimized alongside the elimination rate constant but never reports it. Recovered from the authors' deposited MoBi project (Appendix S1, 'PSP-2025-0082-s02'), HF1000 container MembraneDiffusionCoefficient = 0.0099 in the dm^2/min base unit = 0.99 cm^2/min = 0.594 dm^2/h

    # ------------------------------------------------------------------
    # Midazolam properties (Table 1 and the deposited MoBi project)
    # ------------------------------------------------------------------
    fu <- fixed(0.02); label("Fraction of midazolam unbound in plasma (unitless)")              # Table 1: optimized fraction unbound 2.0% (literature values 3.0 and 1.6%). Corroborated by the deposited MoBi project, Midazolam "Fraction unbound (plasma, reference value)" = 0.02
    krbc <- fixed(0.2492); label("Red-cell to plasma midazolam concentration ratio (unitless)") # NOT PRINTED IN THE PAPER. PK-Sim computes this internally; evaluated here from the deposited MoBi project's own stored formula K_rbc = (f_water_rbc + f_lipids_rbc * 10^LogMA + f_proteins_rbc * KProt) * fu with KProt = (0.81 + 0.11 * 10^LogMA) / 24.92 * 5, LogMA = Lipophilicity = 2.98577 (Table 1 rounds it to 2.99), f_water_rbc 0.625, f_lipids_rbc 0.005, f_proteins_rbc 0.325 and fu 0.02, all stored in the project. Independent check: the implied blood-to-plasma ratio krbc * hct + 1 - hct = 0.775 at hematocrit 0.30 matches the published midazolam value
    permrbc <- fixed(0.3004); label("Plasma to red-cell membrane permeability, P (dm/h)")       # NOT PRINTED IN THE PAPER. Evaluated from the deposited MoBi project's stored formula P = (MWEff * 1e9 / 336)^-6 * 10^LogMA / 5 * 1e-5 with effective molecular weight 286.78 g/mol (Table 1) and LogMA 2.98577, giving 5.0066e-3 dm/min = 0.3004 dm/h
    satovrbc <- fixed(1.002e6); label("Red-cell exchange area per unit red-cell volume (1/dm)") # NOT PRINTED IN THE PAPER. The deposited MoBi project computes the plasma / red-cell exchange area as HCT * A_to_V_bc * 0.6 * V with A_to_V_bc = 167000 1/cm = 1.67e6 1/dm (Organism "Surface/Volume ratio (blood cells)"); the constant carried here is A_to_V_bc * 0.6 = 1.002e6 1/dm and the hematocrit and volume factors are applied in model()

    # ------------------------------------------------------------------
    # Saturable circuit drug loss (Results 3.3). The paper describes this
    # as an "elimination rate constant [that] included drug adsorption,
    # clearance by the dialysis filter and, if present, drug degradation"
    # (Appendix S1), optimized to the ex vivo concentration data.
    # ------------------------------------------------------------------
    vmaxcrrt <- fixed(0.047888); label("Maximum circuit adsorption rate, V_max (mg/h)")         # Results 3.3 prints "V max of 0.002 L/min" at one significant figure. The deposited MoBi project's optimized reaction block CRRT_Midazolam_optimized stores V_max = 0.00245 (Parameter Identification 2, 2022-08-15), which is the printed value unrounded; the executable artifact is taken as authoritative over the rounded print. Converted with the Table 1 molecular weight 325.77 g/mol: 0.00245 umol/min * 325.77 ug/umol * 60 min/h = 0.047888 mg/h. The paper's "L/min" unit label is not dimensionally consistent with the V_max * C / (C + K_m) rate law it parameterizes; see the vignette Errata
    kmcrrt <- fixed(0.358347); label("Circuit adsorption half-saturation concentration, K_m (mg/L)") # Results 3.3: "a final K m of 1100 pmol/mL". Reproduced exactly by the deposited MoBi project (K_m stored as 1.1 in the umol/L base unit = 1100 pmol/mL). Converted with the Table 1 molecular weight 325.77 g/mol: 1.1 umol/L * 325.77 ug/umol = 0.358347 mg/L

    # ------------------------------------------------------------------
    # Residual error
    # ------------------------------------------------------------------
    # The circuit module is a deterministic mechanistic fit: the two
    # Michaelis-Menten constants were optimized by Monte Carlo parameter
    # identification in PK-Sim against the ex vivo concentration data
    # (Methods 2.5) and no residual-error model, standard deviation or
    # between-replicate variance is reported anywhere in the paper, the
    # supplement or the deposited project. Held at 0 rather than invented;
    # see the vignette Errata.
    propSd <- fixed(0); label("Proportional residual error (fraction); magnitude not reported")
  })

  model({
    # ---- 1. Derived circuit flows (Appendix S1, Equations 1-4) --------
    qfil <- qpbp + qrep + qpfr                    # Eq 1: blood-to-dialysate filtration flow
    qeff <- qdia + qfil                           # Eq 2: effluent waste flow off the machine
    qret <- qbfr - qrep - qpfr                    # Eq 4: return-line blood flow
    # Eq 3 (Q_PR, the total fluid removal rate) is a reporting quantity that
    # does not enter any mass transfer and is therefore not carried here.

    # ---- 2. Hematocrit across the filter (Equations 5-6) -------------
    # Dilution by the pre-blood pump lowers the hematocrit entering the
    # filter; removal of plasma water as filtrate raises it on the way out.
    hctin <- hct * qbfr / (qbfr + qpbp)           # Eq 5: pre-filter hematocrit
    hctout <- hct * qbfr / qret                   # Eq 6: post-filter hematocrit
    # The deposited MoBi project splits the filter's plasma and red-cell
    # sub-compartment volumes on the mean of the two.
    hctfil <- (hctin + hctout) / 2

    # ---- 3. Sub-compartment volumes ----------------------------------
    vrespls <- vreservoir * (1 - hct)
    vresrbc <- vreservoir * hct
    vfilpls <- vfilter * (1 - hctfil)
    vfilrbc <- vfilter * hctfil

    # ---- 4. Concentrations -------------------------------------------
    crespls <- reservoir_plasma / vrespls
    cresrbc <- reservoir_rbc / vresrbc
    cfilpls <- circuit_plasma / vfilpls
    cfilrbc <- circuit_rbc / vfilrbc
    cdia <- circuit_dialysate / vdialysate

    # ---- 5. Plasma / red-cell exchange -------------------------------
    # Permeability-limited passive diffusion toward the red-cell to plasma
    # partition ratio, following the deposited MoBi project's
    # PassiveDiffusionPl2RBC transport, SA_bc * fu * (C_pls - C_bc / K_rbc).
    # The permeability-surface products are two to three orders of magnitude
    # larger than the blood flow, so in practice the two phases stay at
    # equilibrium and the red-cell states act as a partition-scaled
    # extension of their plasma partners.
    psres <- satovrbc * vreservoir * hct * fu * permrbc
    psfil <- satovrbc * vfilter * hctfil * fu * permrbc
    jresrbc <- psres * (crespls - cresrbc / krbc)
    jfilrbc <- psfil * (cfilpls - cfilrbc / krbc)

    # ---- 6. Blood delivery to and return from the circuit ------------
    # Equations 7-8. The forward legs carry the whole blood flow split on
    # the systemic hematocrit; the return legs carry the return-line flow
    # split on the post-filter hematocrit, so the plasma volume lost as
    # filtrate is not returned while the red-cell flux balances exactly.
    jinpls <- qbfr * (1 - hct) * crespls          # Eq 8, reservoir -> circuit
    jinrbc <- qbfr * hct * cresrbc                # Eq 7, reservoir -> circuit
    joutpls <- qret * (1 - hctout) * cfilpls      # Eq 8, circuit -> reservoir
    joutrbc <- qret * hctout * cfilrbc            # Eq 7, circuit -> reservoir

    # ---- 7. Transfer across the hemofilter membrane ------------------
    jdiff <- (dfilter * safilter / thfilter) * (fu * cfilpls - cdia)  # Eq 9, diffusion
    jconv <- fu * qfil * cfilpls                                     # Eq 10, convection
    jeff <- cdia * qeff                                              # Eq 11, dialysate -> waste

    # ---- 8. Saturable circuit drug loss ------------------------------
    # V_max * C / (C + K_m) on the filter plasma and red-cell states. The
    # deposited MoBi project applies the MidazolamDegradation reaction to
    # every PrisMax container that is not the filter dialysate, i.e. to
    # both filter sub-compartments (and, as an artifact of the container
    # tagging, to the terminal dialysate bag, which cannot influence any
    # upstream concentration and is omitted here; see the vignette Errata).
    adspls <- vmaxcrrt * cfilpls / (cfilpls + kmcrrt)
    adsrbc <- vmaxcrrt * cfilrbc / (cfilrbc + kmcrrt)

    # ---- 9. ODE system -----------------------------------------------
    d/dt(reservoir_plasma) <- -jinpls + joutpls - jresrbc
    d/dt(reservoir_rbc) <- -jinrbc + joutrbc + jresrbc
    d/dt(circuit_plasma) <- jinpls - joutpls - jfilrbc - jdiff - jconv - adspls
    d/dt(circuit_rbc) <- jinrbc - joutrbc + jfilrbc - adsrbc
    d/dt(circuit_dialysate) <- jdiff + jconv - jeff
    d/dt(circuit_effluent) <- jeff

    # ---- 10. Observation ---------------------------------------------
    # The ex vivo assay sampled the circulating blood reservoir, so the
    # model's readout is the reservoir plasma concentration.
    Cc <- crespls
    Cc ~ prop(propSd)
  })
}
