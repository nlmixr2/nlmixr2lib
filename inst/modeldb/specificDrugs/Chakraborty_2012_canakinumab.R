Chakraborty_2012_canakinumab <- function() {
  description <- "Population pharmacokinetic-binding model for canakinumab (anti-IL-1b IgG1/k monoclonal antibody) and its endogenous target IL-1b in adult cryopyrin-associated periodic syndromes (CAPS) patients (Chakraborty 2012). Two physical compartments (central and peripheral) each carry three species: free canakinumab, free IL-1b, and the canakinumab-IL-1b complex. Drug, ligand, and complex share the same volumes of distribution; complex clearance is set equal to free-drug clearance (CLX = CLD). Distribution between compartments uses two permeability-surface-area coefficients: PSD for free drug and complex, PSL for free ligand. Endogenous IL-1b production RLI enters the peripheral compartment. Drug-ligand binding is solved algebraically under a quasi-steady-state assumption with dissociation constant KD (Hayashi 2007 form). Subcutaneous bioavailability F1 was estimated on the logit scale; this file uses the Sp2/0 cell-line value (commercial Ilaris). Body weight modifies CLD, VC, VP via centred power covariates; serum albumin modifies CLD; age modifies SC ka. Two observed analytes: total canakinumab (free + complex) in ug/mL and total IL-1b (free + complex) in pg/mL."
  reference <- "Chakraborty A, Tannenbaum S, Rordorf C, Lowe PJ, Floch D, Gram H, Roy S. Pharmacokinetic and pharmacodynamic properties of canakinumab, a human anti-interleukin-1b monoclonal antibody. Clin Pharmacokinet. 2012;51(6):e1-e18. doi:10.2165/11599820-000000000-00000."
  vignette <- "Chakraborty_2012_canakinumab"
  units <- list(
    time          = "day",
    dosing        = "mg",
    concentration = "ug/mL (total canakinumab); pg/mL (total IL-1b)"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Centred power covariate on CLD (exponent 0.695), VC (0.684), and VP (0.798). All effects centred at 70 kg (Chakraborty 2012 Table IV; equation 10 page e15).",
      source_name        = "Body weight (WT)"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Centred power covariate on subcutaneous ka (exponent -0.555); reference 34 years (Chakraborty 2012 Table IV; described page e14-e15).",
      source_name        = "Age"
    ),
    ALB = list(
      description        = "Serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Centred power covariate on CLD (exponent -0.916); reference 43 g/L (Chakraborty 2012 Table IV; equation 11 page e15).",
      source_name        = "Serum albumin"
    )
  )

  population <- list(
    n_subjects        = 233L,
    n_studies         = 6L,
    study_names       = c(
      "CACZ885B2101 (healthy + mild asthma)",
      "CACZ885A1101 / NCT00421226 (Japanese healthy volunteers)",
      "CACZ885A2101 / NCT00619905 (rheumatoid arthritis)",
      "CAZC885A2202 / EudraCT 2005-004119-31 (psoriasis)",
      "CACZ885A2102 / NCT00487708 (CAPS phase II, adult and paediatric)",
      "CACZ885D2304 / NCT00465985 (CAPS phase III)"
    ),
    age_range         = "4-74 years (paediatric and adult); reference 34 years for ka",
    age_median        = "34 years (typical CAPS patient)",
    weight_range      = "Includes paediatric to adult subjects; reference 70 kg",
    weight_median     = "70 kg (typical CAPS patient)",
    sex_female_pct    = 42.1,
    sex_notes         = "135 male and 98 female subjects in the population analysis (Chakraborty 2012 Section 5.3, page e15).",
    race_ethnicity    = "Mixed: predominantly White (Caucasian healthy volunteers and patients) and East Asian (Japanese healthy volunteers); other categories not separately reported in the model.",
    disease_state     = "Cryopyrin-associated periodic syndromes (CAPS, primary indication; NALP3 mutations); also healthy volunteers, rheumatoid arthritis, mild asthma, and psoriasis cohorts contributing to the model fit.",
    dose_range        = "Single and multiple intravenous (0.3-10 mg/kg) and subcutaneous (75-600 mg, 2 mg/kg) administrations across the six studies (Chakraborty 2012 Table I and Table III).",
    regions           = "Multinational (Europe, North America, Japan).",
    albumin_reference = "43 g/L (median in CAPS analysis cohort)",
    notes             = "CAPS is the primary regulatory indication; the model parameters in this file are the typical-CAPS-patient estimates from Chakraborty 2012 Table IV. Per-population multipliers for Caucasian HV, Japanese HV, RA, asthma, and psoriasis cohorts are reported in Table IV but not encoded as covariates here -- the validation vignette tabulates them for reference. The cell-line covariate (NS0 vs Sp2/0) acts on ka and F1 only; this file uses the Sp2/0 (commercial Ilaris) values."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters -- typical CAPS patient (70 kg, 43 g/L
    # albumin, 34 years), Sp2/0 cell-line product. All values from
    # Chakraborty 2012 Table IV unless otherwise noted.
    # ------------------------------------------------------------------
    lcl   <- log(0.174);  label("Free-canakinumab clearance CLD (L/day) at 70 kg, 43 g/L albumin, CAPS")          # Chakraborty 2012 Table IV
    lvc   <- log(3.30);   label("Central volume VC (L) at 70 kg")                                                  # Chakraborty 2012 Table IV
    lvp   <- log(2.71);   label("Peripheral volume VP (L) at 70 kg")                                               # Chakraborty 2012 Table IV
    lpsd  <- log(0.429);  label("Drug + complex permeability-surface area coefficient PSD (L/day)")                # Chakraborty 2012 Table IV
    lka   <- log(0.269);  label("Subcutaneous absorption rate ka (1/day) at 34 years, Sp2/0 cell line")            # Chakraborty 2012 Table IV
    lcll  <- log(14.2);   label("Free IL-1b clearance CLL (L/day), CAPS")                                          # Chakraborty 2012 Table IV
    lrli  <- log(9.57);   label("Endogenous IL-1b production rate RLI (ng/day), CAPS")                             # Chakraborty 2012 Table IV
    lkd   <- log(1.07);   label("Apparent IL-1b dissociation constant KD (nmol/L), CAPS")                          # Chakraborty 2012 Table IV
    lpsl  <- log(0.386);  label("Free-ligand permeability-surface area coefficient PSL (L/day)")                   # Chakraborty 2012 Table IV

    # Subcutaneous bioavailability F1 (Sp2/0 cell-line, commercial product).
    # Reported as 70.0% (back-calculated from a logit-scale estimate;
    # Table IV footnote a). Encoded on the logit scale here so the model
    # remains numerically stable if the value is later perturbed.
    lfdepot <- log(0.700 / (1 - 0.700)); label("Logit-scale subcutaneous bioavailability (Sp2/0 cell line; F1 = 0.70)")  # Chakraborty 2012 Table IV

    # ------------------------------------------------------------------
    # Covariate exponents -- centred-power model
    # (Chakraborty 2012 Table IV "Covariates" block; equations 10 and 11).
    # ------------------------------------------------------------------
    e_wt_cl  <-  0.695;  label("Power exponent of body weight on CLD (unitless)")                                  # Chakraborty 2012 Table IV
    e_alb_cl <- -0.916;  label("Power exponent of serum albumin on CLD (unitless)")                                # Chakraborty 2012 Table IV
    e_wt_vc  <-  0.684;  label("Power exponent of body weight on VC (unitless)")                                   # Chakraborty 2012 Table IV
    e_wt_vp  <-  0.798;  label("Power exponent of body weight on VP (unitless)")                                   # Chakraborty 2012 Table IV
    e_age_ka <- -0.555;  label("Power exponent of age on ka (unitless)")                                           # Chakraborty 2012 Table IV

    # ------------------------------------------------------------------
    # IIV -- diagonal (no inter-parameter correlations reported).
    # Variance values are taken directly from the "Interindividual
    # variance" column of Table IV. Note: Chakraborty 2012 reports
    # CV (%) as approximately sqrt(variance) (e.g. variance 0.0859 ->
    # 29% CV; variance 0.406 -> 64% CV), so the variance entries below
    # match the column directly without a log-normal CV-to-variance
    # conversion.
    # ------------------------------------------------------------------
    etalcl  ~ 0.0859    # Chakraborty 2012 Table IV (CLD CV 29%)
    etalvc  ~ 0.0589    # Chakraborty 2012 Table IV (VC  CV 24%)
    etalvp  ~ 0.0817    # Chakraborty 2012 Table IV (VP  CV 29%)
    etalpsd ~ 0.280     # Chakraborty 2012 Table IV (PSD CV 53%)
    etalka  ~ 0.406     # Chakraborty 2012 Table IV (ka  CV 64%; reported on the NS0 row but applies to ka)
    etalcll ~ 0.371     # Chakraborty 2012 Table IV (CLL CV 61%)
    etalrli ~ 0.261     # Chakraborty 2012 Table IV (RLI CV 51%)
    etalkd  ~ 0.395     # Chakraborty 2012 Table IV (KD  CV 63%)
    etalpsl ~ 0.254     # Chakraborty 2012 Table IV (PSL CV 50%)

    # ------------------------------------------------------------------
    # Residual error -- the source paper used Y = log(IPRED) + e
    # with intra-subject variance sigma^2 = 0.0527 (canakinumab) and
    # 0.0840 (IL-1b). Encoded as a proportional model with
    # propSd = sqrt(sigma^2), the standard small-sigma approximation
    # of the additive-on-log residual.
    # ------------------------------------------------------------------
    propSd          <- sqrt(0.0527);  label("Proportional residual error on total canakinumab (fraction)")        # Chakraborty 2012 Table IV (canakinumab sigma^2 = 0.0527)
    propSd_totalIL1b <- sqrt(0.0840); label("Proportional residual error on total IL-1b (fraction)")               # Chakraborty 2012 Table IV (IL-1b      sigma^2 = 0.0840)
  })

  model({
    # ------------------------------------------------------------------
    # Constants -- molecular weights (Chakraborty 2012 equation 5
    # footnote: 17000 ng/nmol for IL-1b; canakinumab is reported as
    # ~150 kDa in section 4.2.2 page e10).
    # ------------------------------------------------------------------
    MWX <- 149.117  # canakinumab molecular weight (kDa = ng/nmol; IgG1/k ~149 kDa, exact value from sequence-derived mass)
    MWE <- 17.000   # IL-1b mature protein molecular weight (kDa); paper uses 17000 ng/nmol

    # ------------------------------------------------------------------
    # 1. Individual parameters (centred-power covariate model;
    #    Chakraborty 2012 equation 10 / equation 11).
    # ------------------------------------------------------------------
    cl    <- exp(lcl  + etalcl)  * (WT / 70)^e_wt_cl * (ALB / 43)^e_alb_cl
    vc    <- exp(lvc  + etalvc)  * (WT / 70)^e_wt_vc
    vp    <- exp(lvp  + etalvp)  * (WT / 70)^e_wt_vp
    psd   <- exp(lpsd + etalpsd)
    ka    <- exp(lka  + etalka)  * (AGE / 34)^e_age_ka
    cll   <- exp(lcll + etalcll)
    rli   <- exp(lrli + etalrli)
    kd    <- exp(lkd  + etalkd)
    psl   <- exp(lpsl + etalpsl)
    cl_x  <- cl                    # Complex clearance assumed equal to free-drug clearance (Chakraborty 2012 page e7)
    fdepot <- 1 / (1 + exp(-lfdepot))

    # ------------------------------------------------------------------
    # 2. Quasi-steady-state binding (Chakraborty 2012 equation 7):
    #    KD * V + TD + TL = sum, then
    #    X = [sum - sqrt(sum^2 - 4 * TD * TL)] / 2
    #    applied independently in central (subscript C) and peripheral
    #    (subscript P). All quantities in nmol; KD in nmol/L; V in L.
    # ------------------------------------------------------------------
    sumC <- kd * vc + central + central_il1b
    xC   <- 0.5 * (sumC - sqrt(sumC * sumC - 4 * central * central_il1b))
    sumP <- kd * vp + peripheral1 + peripheral1_il1b
    xP   <- 0.5 * (sumP - sqrt(sumP * sumP - 4 * peripheral1 * peripheral1_il1b))

    # Free amounts (nmol).
    fdC <- central     - xC
    fdP <- peripheral1 - xP
    flC <- central_il1b    - xC
    flP <- peripheral1_il1b - xP

    # ------------------------------------------------------------------
    # 3. ODE system (Chakraborty 2012 equations 1-5).
    #    State units: depot (mg), central / peripheral1 (nmol total drug),
    #    central_il1b / peripheral1_il1b (nmol total ligand).
    #    The depot-to-central transfer converts mg to nmol via 1000/MWX.
    #    Endogenous IL-1b production RLI (ng/day) is converted to nmol/day
    #    via /17000 = /(MWE * 1000) consistent with the paper's eq. 5.
    # ------------------------------------------------------------------
    d/dt(depot)             <- -ka * depot
    d/dt(central)           <-  ka * depot * (1000 / MWX) - fdC * cl / vc - xC * cl_x / vc + (peripheral1 / vp - central / vc) * psd
    d/dt(peripheral1)       <-  (central / vc - peripheral1 / vp) * psd
    d/dt(central_il1b)    <- -flC * cll / vc - xC * cl_x / vc + (flP / vp - flC / vc) * psl + (xP / vp - xC / vc) * psd
    d/dt(peripheral1_il1b) <-  rli / (MWE * 1000) + (flC / vc - flP / vp) * psl + (xC / vc - xP / vp) * psd

    # ------------------------------------------------------------------
    # 4. Steady-state initial conditions for endogenous IL-1b.
    #    With no drug at t = 0: complex amounts are zero, so total
    #    ligand = free ligand. Solving dTLC/dt = 0 and dTLP/dt = 0
    #    gives FLC = RLI * VC / (MWE * 1000 * CLL) and
    #    FLP = (VP/VC) * FLC * (1 + CLL/PSL).
    # ------------------------------------------------------------------
    central_il1b(0)    <- rli * vc / (MWE * 1000 * cll)
    peripheral1_il1b(0) <- (vp / vc) * (rli * vc / (MWE * 1000 * cll)) * (1 + cll / psl)

    # ------------------------------------------------------------------
    # 5. Bioavailability for SC depot (Chakraborty 2012 page e9: F
    #    coded via logit transformation; Sp2/0 typical = 70%).
    # ------------------------------------------------------------------
    f(depot) <- fdepot

    # ------------------------------------------------------------------
    # 6. Observation outputs in assay units.
    #    Dimensional rule: 1 nmol/L * 1 kDa = 1 ng/mL.
    #    Cc (ug/mL)        = (central / VC) [nmol/L] * MWX [kDa] / 1000
    #    totalIL1b (pg/mL) = (central_il1b / VC) [nmol/L] * MWE [kDa] * 1000
    # ------------------------------------------------------------------
    Cc        <- (central        / vc) * MWX / 1000
    totalIL1b <- (central_il1b / vc) * MWE * 1000

    Cc        ~ prop(propSd)
    totalIL1b ~ prop(propSd_totalIL1b)
  })
}
