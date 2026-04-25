Wang_2020_ontamalimab <- function() {
  description <- "Two-compartment population PK model for ontamalimab (SHP647), a fully human IgG2 anti-MAdCAM-1 monoclonal antibody, in adults with moderate-to-severe ulcerative colitis or Crohn's disease (Wang 2020), with first-order SC absorption, absorption lag time, parallel linear and Michaelis-Menten elimination from the central compartment, and allometric weight scaling on CL, Vc, CLd, Vp, and Vmax."
  reference <- "Wang Y, Marier J-F, Kassir N, Chabot JR, Smith B, Cao C, Lewis L, Dorner AJ, Padula SJ, Banfield C. Population Pharmacokinetics and Pharmacodynamics of Ontamalimab (SHP647), a Fully Human Monoclonal Antibody Against Mucosal Addressin Cell Adhesion Molecule-1 (MAdCAM-1), in Patients With Ulcerative Colitis or Crohn's Disease. J Clin Pharmacol. 2020 Jul;60(7):903-914. doi:10.1002/jcph.1590"
  vignette <- "Wang_2020_ontamalimab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight (baseline)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 70 kg. Power effects: CL exp 0.0034, Vc exp 0.635, CLd exp 0.0034, Vp exp 0.635, Vmax exp 1.89 per Wang 2020 Table 2.",
      source_name        = "WT"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 39 g/L (overall median). Power effect on CL with exponent -0.889 per Wang 2020 Table 2 and table footnote (ALB unit g/L).",
      source_name        = "ALB"
    ),
    CRP = list(
      description        = "Baseline C-reactive protein",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 0.837 mg/dL (overall median). Power effect on CL with exponent 0.147 per Wang 2020 Table 2 and table footnote (CRP unit mg/dL). Note: source unit is mg/dL, not the canonical mg/L; the reference value and effect coefficient are inseparable from the unit. 1 mg/dL = 10 mg/L; the equivalent reference in mg/L would be 8.37 mg/L.",
      source_name        = "CRP"
    )
  )

  population <- list(
    n_subjects     = 440L,
    n_observations = 2138L,
    n_studies      = 2L,
    age_range      = "18-68 years",
    age_median     = "37 years",
    weight_range   = "35.6-155 kg",
    weight_median  = "68.8 kg",
    sex_female_pct = 51.1,
    race_ethnicity = c(White = 86.1, Asian = 10.0, Other = 3.9),
    disease_state  = "Moderate-to-severe ulcerative colitis (UC, n = 249; 56.6%) or Crohn's disease (CD, n = 191; 43.4%) with prior failure or intolerance to immunosuppressants and/or anti-TNF agents.",
    dose_range     = "Subcutaneous ontamalimab 7.5, 22.5, 75, or 225 mg every 4 weeks (3 doses; days 1, 28, 56) over a 12-week treatment / induction period.",
    regions        = "Multi-regional (phase 2 OPERA NCT01276509 in CD and phase 2b TURANDOT NCT01620255 in UC).",
    reference_subject = "70-kg patient with UC or CD, baseline albumin 39 g/L, baseline CRP 0.837 mg/dL (overall median values per Wang 2020 Table 2 footnote).",
    notes          = "Pooled phase 2 dataset; 9.4% (202/2138) of ontamalimab samples were below the 10 ng/mL LLOQ and were treated as missing for the population PK analysis. Disease status (CD vs UC), sex, time-varying ADA status, AST, ALT, bilirubin, and age were tested as covariates and were not retained in the final model; only baseline albumin and baseline CRP entered the final CL/F covariate equation alongside allometric weight scaling."
  )

  ini({
    # Structural PK parameters - Wang 2020 Table 2 final-model estimates. Paper
    # reports rate constants in 1/h, clearances in L/h, lag in h, and Km in
    # ng/mL; values below convert to nlmixr2lib defaults (time in day, mass in
    # mg, concentration in mg/L). Conversions:
    #   1/h -> 1/day:  multiply by 24
    #   L/h -> L/day:  multiply by 24
    #   h -> day:      divide by 24
    #   ug/h -> mg/day: multiply by 24/1000
    #   ng/mL -> mg/L: divide by 1000
    lka     <- log(0.0187 * 24);              label("Absorption rate Ka (1/day)")                                # Wang 2020 Table 2: Ka = 0.0187 1/h
    ltlag   <- log(2.61 / 24);                label("Absorption lag time (day)")                                  # Wang 2020 Table 2: Lag = 2.61 h
    lcl     <- log(0.0127 * 24);              label("Apparent linear clearance CL/F (L/day) at reference covariates") # Wang 2020 Table 2: CL/F = 0.0127 L/h
    lvc     <- log(6.53);                     label("Apparent central volume Vc/F (L) at reference covariates")  # Wang 2020 Table 2: Vc/F = 6.53 L
    lcld    <- log(0.000345 * 24);            label("Apparent inter-compartmental clearance CLd/F (L/day) at reference covariates") # Wang 2020 Table 2: CLd/F = 0.000345 L/h
    lvp     <- log(0.0216);                   label("Apparent peripheral volume Vp/F (L) at reference covariates") # Wang 2020 Table 2: Vp/F = 0.0216 L
    lvmax   <- log(5.87 * 24 / 1000);         label("Apparent Michaelis-Menten Vmax/F (mg/day) at reference WT")  # Wang 2020 Table 2: Vmax/F = 5.87 ug/h
    lkm     <- log(19.0 / 1000);              label("Michaelis-Menten constant Km (mg/L)")                        # Wang 2020 Table 2: Km = 19.0 ng/mL

    # Covariate exponents - Wang 2020 Table 2 final-model equations.
    e_wt_cl    <-  0.0034; label("Power exponent of (WT/70 kg) on CL/F (unitless)")                              # Wang 2020 Table 2: CL/F ~ WT
    e_wt_vc    <-  0.635;  label("Power exponent of (WT/70 kg) on Vc/F (unitless)")                              # Wang 2020 Table 2: Vc/F ~ WT
    e_wt_cld   <-  0.0034; label("Power exponent of (WT/70 kg) on CLd/F (unitless)")                             # Wang 2020 Table 2: CLd/F ~ WT
    e_wt_vp    <-  0.635;  label("Power exponent of (WT/70 kg) on Vp/F (unitless)")                              # Wang 2020 Table 2: Vp/F ~ WT
    e_wt_vmax  <-  1.89;   label("Power exponent of (WT/70 kg) on Vmax/F (unitless)")                            # Wang 2020 Table 2: Vmax/F ~ WT
    e_alb_cl   <- -0.889;  label("Power exponent of (ALB/39 g/L) on CL/F (unitless)")                            # Wang 2020 Table 2: CL/F ~ ALB
    e_crp_cl   <-  0.147;  label("Power exponent of (CRP/0.837 mg/dL) on CL/F (unitless)")                       # Wang 2020 Table 2: CL/F ~ CRP

    # Inter-individual variability - Wang 2020 Table 2 reports BSV as CV%.
    # Convert log-normal variance via omega^2 = log(CV^2 + 1):
    #   Ka   CV 61.8% -> log(0.618^2 + 1) = 0.3232
    #   CL/F CV 54.6% -> log(0.546^2 + 1) = 0.2611
    #   Vc/F CV 41.0% -> log(0.410^2 + 1) = 0.1554
    # The paper does not report off-diagonal correlations, so IIVs are taken
    # as independent.
    etalka ~ 0.3232; label("IIV variance on log Ka (Wang 2020 Table 2: 61.8% CV)")
    etalcl ~ 0.2611; label("IIV variance on log CL/F (Wang 2020 Table 2: 54.6% CV)")
    etalvc ~ 0.1554; label("IIV variance on log Vc/F (Wang 2020 Table 2: 41.0% CV)")

    # Residual error - Wang 2020 Table 2 combined additive (166 ng/mL) and
    # proportional (19.6%) error.
    addSd  <- 0.166; label("Additive residual error (mg/L)")     # Wang 2020 Table 2: additive 166 ng/mL = 0.166 mg/L
    propSd <- 0.196; label("Proportional residual error (fraction)") # Wang 2020 Table 2: proportional 19.6%
  })

  model({
    # Individual PK parameters with Wang 2020 covariate models. Reference
    # subject: WT 70 kg, ALB 39 g/L, CRP 0.837 mg/dL.
    cl   <- exp(lcl + etalcl) *
            (WT / 70)^e_wt_cl *
            (ALB / 39)^e_alb_cl *
            (CRP / 0.837)^e_crp_cl
    vc   <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    cld  <- exp(lcld) * (WT / 70)^e_wt_cld
    vp   <- exp(lvp) * (WT / 70)^e_wt_vp
    vmax <- exp(lvmax) * (WT / 70)^e_wt_vmax
    km   <- exp(lkm)
    ka   <- exp(lka + etalka)
    tlag <- exp(ltlag)

    Cc <- central / vc

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot -
                          (cl / vc) * central -
                          vmax * Cc / (km + Cc) -
                          (cld / vc) * central +
                          (cld / vp) * peripheral1
    d/dt(peripheral1) <-  (cld / vc) * central - (cld / vp) * peripheral1

    alag(depot) <- tlag

    Cc ~ add(addSd) + prop(propSd)
  })
}
