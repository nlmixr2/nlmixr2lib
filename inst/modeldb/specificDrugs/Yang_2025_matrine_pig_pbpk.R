Yang_2025_matrine_pig_pbpk <- function() {
  description <- "PBPK (minimal, flow-limited; pig). Four-compartment physiologically based pharmacokinetic model for the quinolizidine alkaloid matrine in the intestinal lumen of pigs after oral administration, comprising the intestinal lumen (the sampled matrix, reached from a first-order gastric-emptying depot), liver, a lumped 'other organs' compartment carrying renal excretion, and a single well-mixed blood pool. Liver and other organs are perfusion (flow) limited with tissue-to-blood partition coefficients driven by the unbound blood concentration; drug leaves the intestinal lumen by first-order absorption into blood in competition with first-order faecal excretion, and re-enters it by first-order biliary excretion from the liver, which is what produces the observed two-phase luminal decay. Built to support dosage-regimen design against enterotoxigenic Escherichia coli colonising the small intestine, the site of infection in porcine colibacillosis (Yang 2025)."
  reference   <- "Yang B, Jia Y, Wang F, Lv X, Ma S, Tan Y, Zhang W, Wan D, Li R, Zhou D, Yu D. The kinetic behavior of matrine in pig intestinal lumen after oral administration and its physiologically based pharmacokinetic modeling. Front Vet Sci. 2025;12:1620161. doi:10.3389/fvets.2025.1620161"
  vignette    <- "Yang_2025_matrine_pig_pbpk"
  units       <- list(time = "h", dosing = "ug", concentration = "ug/L")
  # Oral dose lands on the reconstructed gastric depot, not on `depot` or
  # `central`, so buildModelDb() cannot infer the route - state it explicitly
  # or the registry records dosing = NA.
  dosing      <- "stomach"

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Every state integrates a drug AMOUNT in ug; the
  # per-compartment concentrations are formed algebraically in model().
  compartmentData <- list(
    # `stomach` is the reconstructed gastric depot - see the Errata note in the
    # vignette and the lfdepot comment below. It holds swallowed drug that has
    # not yet reached a sampled matrix.
    stomach   = list(analyte = "matrine", units = "ug", specimen = "administration site", verified = TRUE),
    # The intestinal lumen IS the sampled matrix in this paper (ileal digesta
    # drawn through a T-cannula). The specimen vocabulary has no
    # 'intestinal contents' token; 'faeces' is its closest luminal-digesta
    # member and is the matrix this state drains into.
    gut_lumen = list(analyte = "matrine", units = "ug", specimen = "faeces", verified = TRUE),
    liver     = list(analyte = "matrine", units = "ug", specimen = "tissue", verified = TRUE),
    other     = list(analyte = "matrine", units = "ug", specimen = "tissue", verified = TRUE),
    blood     = list(analyte = "matrine", units = "ug", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Scales every compartment volume (Vx = Vcx * WT), the cardiac output (Qtot = Qcar * WT) and the renal clearance (CLrenal = Clrenal * WT), all of which Yang 2025 Table 5 reports as fractions of body weight or per-kilogram rates. Because the dose is also given per kilogram, luminal concentrations are body-weight invariant - which is exactly what Yang 2025 Section 4 found by sensitivity analysis (the normalised sensitivity coefficient for body weight was far below 0.1, and the model built on 24.3-31.1 kg pigs predicted the 9.8-10.3 kg experiment 2 cohort without adjustment).",
      source_name        = "BW"
    )
  )

  population <- list(
    species        = "pig (crossbred Landrace x Large White x Duroc for experiment 1; Landrace x Large White for experiment 2)",
    n_subjects     = 37L,
    n_studies      = 1L,
    age_range      = NA_character_,
    weight_range   = "24.3-31.1 kg (experiment 1, n = 12); 9.8-10.3 kg (experiment 2, n = 25) (Yang 2025 Section 2.2 and Table 2)",
    sex_female_pct = NA_real_,
    disease_state  = "Healthy pigs; matrine given as a potential resistance-reversal agent for porcine colibacillosis rather than as treatment of an established infection",
    dose_range     = "Experiment 1: single oral 40 mg/kg and, after a 2-week washout, 70 mg/kg, each given alone (group A) or with amoxicillin 40 mg/kg (group B), n = 6 per group. Experiment 2: 50 mg/kg/day by oral gavage for 5 consecutive days, n = 25 (Yang 2025 Section 2.3 and Table 2)",
    regions        = "China",
    n_observations = NA_integer_,
    notes          = "Model parameterised against a single arm - matrine alone, 40 mg/kg, experiment 1 (Yang 2025 Table 2, 'Model parameterization') - and evaluated against the other four arms. Experiment 1 pigs carried a sterile T-cannula implanted in the terminal ileum about 15 cm cranial to the ileocecal valve, through which roughly 2 g of intestinal contents were drawn at 0.25, 0.5, 1, 2, 3, 4, 5, 8, 12, 16, 24, 36, 48, 72 and 120 h post-dose; the 0.25 h samples could not be collected because the pigs had been fasted for 12 h, so the reported absorption phase begins at 0.5 h (Yang 2025 Section 4). Experiment 2 sacrificed five pigs at each of 0.5, 1, 3, 6 and 12 days after the last dose. Physiological parameters are pig population means from the literature (Yang 2025 refs 22-24); compound-specific parameters were seeded from rat studies (refs 14, 15) and then optimised in acslXtreme 2.5.0.6 by Nelder-Mead maximum likelihood against the 40 mg/kg matrine-alone arm. No individual-level fit was performed, so the model carries no inter-individual random effects and no residual-error model; Yang 2025 evaluated it by linear regression and Pearson correlation of predicted against observed concentrations (Table 6) plus a one-sample t-test on the derived PK parameters."
  )

  ini({
    # ---------------------------------------------------------------------
    # Compound-specific parameters, all from Yang 2025 Table 5 'Optimized
    # value'. Yang 2025 reports point estimates with no standard errors,
    # confidence intervals or random effects - the model is deterministic -
    # so every parameter is FIXED.
    #
    # Sourcing note (Yang 2025 Section 2.6): ka, Pot, ke and kbi are marked
    # 'By optimization' in Table 5; kst, F, Pl and Clrenal carry a literature
    # starting value (refs 14, 15, 21-23) that was then adjusted by the same
    # Nelder-Mead fit, so the OPTIMIZED column - not the starting value - is
    # what this file encodes throughout.
    # ---------------------------------------------------------------------
    lkst      <- fixed(log(0.8544925));   label("First-order gastric emptying rate constant, stomach to intestinal lumen (kst, 1/h)")     # Yang 2025 Table 5 (starting 0.1 from refs 21, 22; optimized 0.8544925)
    lka       <- fixed(log(0.8555538));   label("First-order absorption rate constant from intestinal lumen to blood (ka, 1/h)")          # Yang 2025 Table 5 (starting 0.3; optimized 0.8555538, by optimization)
    lkfec     <- fixed(log(0.007358172)); label("First-order faecal excretion rate constant for unabsorbed drug (ke, 1/h)")               # Yang 2025 Table 5 (starting 0.01; optimized 0.007358172, by optimization)
    lkbile    <- fixed(log(0.05834594));  label("First-order biliary excretion rate constant, liver to intestinal lumen (kbi, 1/h)")      # Yang 2025 Table 5 (starting 0.01; optimized 0.05834594, by optimization)
    lfdepot   <- fixed(log(0.7925891));   label("Oral bioavailability (F, fraction)")                                                     # Yang 2025 Table 5 (starting 0.171 from ref 14; optimized 0.7925891)
    lkp_liver <- fixed(log(2.936615));    label("Liver-to-blood partition coefficient (Pl, unitless)")                                    # Yang 2025 Table 5 (starting 5.5 from ref 15; optimized 2.936615)
    lkp_other <- fixed(log(11.33514));    label("Other organs-to-blood partition coefficient (Pot, unitless)")                            # Yang 2025 Table 5 (starting 5; optimized 11.33514, by optimization)
    lcl_renal <- fixed(log(0.2805897));   label("Renal clearance per unit body weight (Clrenal, L/(h*kg))")                               # Yang 2025 Table 5 (starting 1.182 from ref 14; optimized 0.2805897)

    # Yang 2025 reports the plasma protein BINDING Pb = 0.319, i.e. the BOUND
    # fraction, and uses it in the Table 1 blood equation as the factor
    # (1 - Pb). The canonical nlmixr2lib parameter is the UNBOUND fraction, so
    # this line carries fu = 1 - Pb = 1 - 0.319 = 0.681. See the vignette
    # source-trace table.
    fu        <- fixed(0.681);            label("Fraction of matrine unbound in blood (fu = 1 - Pb, fraction)")                           # Yang 2025 Table 5 (Pb starting 0.3; optimized 0.319, by optimization); fu = 1 - 0.319

    # ---------------------------------------------------------------------
    # Residual error. Yang 2025 is a deterministic PBPK model fitted to
    # cohort mean concentration-time data; it reports no residual-error model
    # (performance was judged by predicted-vs-observed linear regression and
    # Pearson correlation, Table 6). The placeholder below exists only so the
    # model is a syntactically complete nlmixr2 object for forward
    # simulation; it is NOT paper-derived. Same convention as
    # Yang_2023_diclazuril_chicken_pbpk.R and Ai_2024_ractopamine_goat_pbpk.R.
    # ---------------------------------------------------------------------
    propSd_Cgut_lumen <- fixed(0.10); label("Proportional residual error placeholder, intestinal lumen (fraction)") # not reported in Yang 2025; placeholder for syntactic completeness only
  })

  model({
    # ==================================================================
    # Pig physiology. Yang 2025 Table 5 rows sourced to refs 22-24
    # (Guerin 2001 gastric emptying, Kokue 1988 pig sulfa-drug PK and
    # gastric emptying, Upton 2008 organ weights and blood flows). These
    # are literature physiology rather than fitted quantities, so they
    # are carried as traceable literals here rather than as ini()
    # parameters - the same convention as
    # Yang_2023_diclazuril_chicken_pbpk.R and
    # Ai_2024_ractopamine_goat_pbpk.R.
    #
    # Transcription check on Table 5: the two 'By calculation' rows
    # reproduce exactly. Qc_other = 1 - Qc_liver = 1 - 0.3 = 0.7, and
    # Vc_other = 1 - (Vc_gut_lumen + Vc_liver + Vc_blood)
    #          = 1 - (0.0136733 + 0.0294 + 0.06) = 0.8969267, matching
    # Table 5 to all seven printed digits.
    #
    # Yang 2025 Table 5 also reports a haematocrit PCV = 0.333 (starting
    # value = optimized value, ref 23). It appears in none of the four
    # Table 1 mass-balance equations and has no role in this model
    # structure, so it is recorded here for completeness only and is
    # deliberately not used below.
    # ==================================================================
    qcar         <- 4.944                 # cardiac output, L/(h*kg) BW;    Yang 2025 Table 5 (Qcar, ref 23)
    qc_liver     <- 0.3                   # liver blood flow / Qcar;        Yang 2025 Table 5 (Qcl, ref 23)
    qc_other     <- 0.7                   # other-organ blood flow / Qcar;  Yang 2025 Table 5 (Qcot, by calculation)

    vc_gut_lumen <- 0.0136733             # intestinal-content volume / BW; Yang 2025 Table 5 (Vcint, starting 0.01 from ref 23; optimized)
    vc_liver     <- 0.0294                # liver volume / BW;              Yang 2025 Table 5 (Vcl, ref 23)
    vc_blood     <- 0.06                  # blood volume / BW;              Yang 2025 Table 5 (Vcb, ref 23)
    vc_other     <- 0.8969267             # other-organ volume / BW;        Yang 2025 Table 5 (Vcot, by calculation)

    # ------------------------------------------------------------------
    # Compound-specific parameters, back-transformed from the log scale.
    # ------------------------------------------------------------------
    kst      <- exp(lkst)
    ka       <- exp(lka)
    kfec     <- exp(lkfec)
    kbile    <- exp(lkbile)
    kp_liver <- exp(lkp_liver)
    kp_other <- exp(lkp_other)

    # ------------------------------------------------------------------
    # Compartment volumes (L) and plasma flows (L/h). Yang 2025 gives
    # every physiological quantity as a fraction of body weight or a
    # per-kilogram rate, so all of these scale linearly with WT.
    # ------------------------------------------------------------------
    v_gut_lumen <- vc_gut_lumen * WT
    v_liver     <- vc_liver     * WT
    v_blood     <- vc_blood     * WT
    v_other     <- vc_other     * WT

    q_total     <- qcar * WT
    q_liver     <- q_total * qc_liver
    q_other     <- q_total * qc_other
    cl_renal    <- exp(lcl_renal) * WT

    # ------------------------------------------------------------------
    # Compartment concentrations (ug/L; taken as ug/kg for intestinal
    # contents at unit density, per Yang 2025 assumption (iii): the MT
    # concentration in the intestinal lumen equals the value measured in
    # the intestinal contents). States are drug AMOUNTS in ug.
    # ------------------------------------------------------------------
    Cgut_lumen <- gut_lumen / v_gut_lumen   # the sampled matrix: Yang 2025 Figures 3, 4 and 5
    Cc         <- blood     / v_blood       # whole-blood concentration; not measured in this study
    Cliver     <- liver     / v_liver
    Cother     <- other     / v_other

    # Only unbound drug leaves the blood into the perfused tissues. Yang
    # 2025 Table 1 writes this factor explicitly in the blood equation
    # (- Qtot * Ca * (1 - Pb)) but omits it from the two tissue
    # equations, which would leave the printed system creating mass at
    # + Qtot * Ca * Pb per hour - roughly 1.58 * Ca L/(h*kg) against a
    # total elimination of only 0.28 * Ca L/(h*kg), i.e. divergent. The
    # unbound driving concentration below is the mass-conserving reading;
    # see the vignette Errata.
    Cfree      <- fu * Cc

    # ------------------------------------------------------------------
    # Mass-balance ODEs (Yang 2025 Table 1 and Figure 1).
    #
    # Dosing: the oral dose (ug) lands on `stomach` and is multiplied by
    # the bioavailability F through f(stomach) below, then drains into the
    # intestinal lumen at the first-order gastric-emptying rate kst. Yang
    # 2025 Table 1 instead prints the luminal input as a bare
    # `DOSE * kst` with DOSE defined in the footnote as the (constant)
    # oral dose, which as typeset is an infusion that never switches off,
    # and F - a Table 5 optimized parameter with a reported sensitivity
    # result - appears in none of the four printed equations. The depot
    # form below is the reconstruction; it is confirmed against the
    # parameterisation arm by three independent closed forms (Tmax, Cmax,
    # AUC) and by the signs of the paper's own sensitivity analysis, all
    # of which are tabulated in the vignette.
    # ------------------------------------------------------------------
    d/dt(stomach)   <- -kst * stomach
    d/dt(gut_lumen) <- kst * stomach + kbile * liver - ka * gut_lumen - kfec * gut_lumen   # Table 1, intestinal lumen
    d/dt(liver)     <- q_liver * (Cfree - Cliver / kp_liver) - kbile * liver               # Table 1, liver
    d/dt(other)     <- q_other * (Cfree - Cother / kp_other) - cl_renal * Cc               # Table 1, other organs (renal excretion is incorporated here; Section 2.6)
    d/dt(blood)     <- ka * gut_lumen +
                       q_liver * Cliver / kp_liver +
                       q_other * Cother / kp_other -
                       q_total * Cfree                                                     # Table 1, blood

    f(stomach) <- exp(lfdepot)   # Yang 2025 Table 5 bioavailability F, applied to the oral dose

    # ------------------------------------------------------------------
    # Observation. Yang 2025 measured matrine in pig intestinal contents
    # only; blood, liver and other-organ concentrations are model
    # internals available as algebraic outputs but were never sampled.
    # ------------------------------------------------------------------
    Cgut_lumen ~ prop(propSd_Cgut_lumen)
  })
}
