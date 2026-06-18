Rosario_2008_maraviroc <- function() {
  description <- "Two-compartment population PK model with first-order absorption and lag time for maraviroc (CCR5 antagonist) coupled with a direct Emax CCR5 receptor occupancy model in healthy adults and HIV-1-positive patients (Rosario 2008). PK is parameterised with dose-dependent relative bioavailability F1 and dose-dependent elimination rate constant K across six dose groups (3, 10, 25, 100 (reference), 300 mg b.i.d. and 600 mg q.d.); receptor occupancy on CD4 T cells is modelled as Occ = E0 + Emax * Cp / (KD + Cp) with a background binding baseline."
  reference <- paste(
    "Rosario MC, Jacqmin P, Dorr P, James I, Jenkins TM, Abel S, van der Ryst E.",
    "Population pharmacokinetic/pharmacodynamic analysis of CCR5 receptor occupancy",
    "by maraviroc in healthy subjects and HIV-positive patients.",
    "Br J Clin Pharmacol. 2008;65 Suppl 1:86-94.",
    "doi:10.1111/j.1365-2125.2008.03140.x"
  )
  vignette <- "Rosario_2008_maraviroc"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    DOSE = list(
      description        = "Per-subject assigned maraviroc dose level (mg) carried as a continuous covariate column. Each subject was randomised to a single dose-cohort across the 12-day (healthy volunteers) or 10-day (HIV patients) treatment period.",
      units              = "mg",
      type               = "continuous",
      reference_category = "100 (b.i.d.) -- the reference dose group anchored at F1 = 1.0 and kel = 0.288 1/h",
      notes              = "Used as a multiplicative-indicator covariate on the depot-compartment bioavailability F and on the central-compartment elimination rate constant kel. Six dose levels exist in the source design (3, 10, 25, 100, 300 mg b.i.d. and 600 mg q.d.); the model() block expresses F and kel as (DOSE == d) sums of the five non-reference fixed-effect ratios, so for any other DOSE value both sums collapse to zero, implying an inactive dose and producing a zero plasma concentration. Set DOSE to one of {3, 10, 25, 100, 300, 600} per subject for valid simulation; downstream users wanting to extrapolate to off-grid dose levels need to add an interpolation rule in their own driver code. Subject-level (time-fixed) within the source study.",
      source_name        = "Dose (Table 2 dose-group label; the underlying NONMEM dataset column name is not reported by the paper)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 88L,
    n_studies      = 2L,
    study_names    = c("A4001002 (healthy volunteers, n = 65 enrolled; n = 55 with PK and n = 36 with receptor-occupancy data)",
                       "A4001007 (asymptomatic HIV-1-positive patients, n = 25 enrolled; n = 17 with PK and n = 24 with receptor-occupancy data; n = 23 with viral-load data)"),
    n_observations = 3621L,
    n_pk_observations = 2770L,
    n_pd_observations = 851L,
    age_range      = "Not reported in the paper; healthy volunteers and asymptomatic HIV-1-positive adult patients (study A4001007 inclusion required asymptomatic HIV-1 infection with CD4 > 250 cells/mm^3, viral load >= 5000 copies/mL, and CCR5-tropic-only virus)",
    sex_female_pct = NA_real_,
    disease_state  = "Mixed: healthy adult volunteers (study A4001002) and asymptomatic HIV-1-positive adults (study A4001007). HIV patients were antiretroviral-drug-naive or off antiretroviral therapy for >= 8 weeks prior to enrolment and had CCR5-tropic-only virus.",
    dose_range     = "Multiple-dose oral maraviroc 3, 10, 25, 100, or 300 mg b.i.d. or 600 mg q.d. for 12 days (study A4001002) or 25 mg q.d. or 100 mg b.i.d. for 10 days (study A4001007). Healthy volunteer 300 mg b.i.d. and 600 mg q.d. cohorts contributed PK only (no receptor-occupancy sampling).",
    regions        = "Pfizer Research Clinic (Hospital Erasme, Brussels, Belgium) for the healthy volunteer study; multicentre design for the HIV patient study (sites not enumerated by the paper).",
    notes          = "Subject distribution (Table 1): 8 placebo + 5 (3 mg) + 5 (10 mg) + 9 (25 mg) + 9 (100 mg) + 9 (300 mg) + 18 (600 mg) healthy volunteers; 8 placebo + 9 (25 mg q.d.) + 8 (100 mg b.i.d.) HIV patients. The 25 mg q.d. patient cohort had one subject who withdrew on day 1 with PK only (no receptor-occupancy or viral-load data). The 100 mg b.i.d. patient cohort had one viral-load exclusion (subject with dual/mixed-tropic virus). Maraviroc plasma LLOQ was 0.5 ng/mL except in the 3 mg and 10 mg b.i.d. cohorts where LLOQ was 0.1 ng/mL (centralised LC-MS/MS, Maxxam Analytics). Receptor occupancy was measured by an ex-vivo MIP-1-beta internalization assay using a fluorescently labelled anti-CCR5 2D7 monoclonal antibody (Esoterix Inc., centralised flow cytometry). The PK model is acknowledged by the authors as exploratory (Discussion: 'The development of a definitive population PK model for maraviroc was not the aim of this analysis'); a comprehensive popPK from Phase 1/2a is described in their reference [9] (Abel et al., Br J Clin Pharmacol 2008 Suppl 1)."
  )

  ini({
    # =========================================================================
    # Structural PK -- Rosario 2008 Table 2. NONMEM ADVAN4 / two-compartment
    # with first-order absorption and lag time, rate-constant parameterisation.
    # Paper's K (elimination) and F1 (relative bioavailability) are reported
    # per dose group; 100 mg b.i.d. is the F1 = 1.00 reference, and the same
    # dose group anchors the kel value used as the structural intercept.
    # NONMEM compartment numbering in the paper (depot=1, central=2,
    # peripheral=3) maps to nlmixr2lib canonical names (depot, central,
    # peripheral1) so paper K23 -> canonical k12 and paper K32 -> canonical k21.
    # =========================================================================

    lka   <- log(1.14);   label("Absorption rate constant ka (1/h)")                       # Table 2: Ka = 1.14 /h (SE 7.6%)
    ltlag <- log(0.01);   label("Absorption lag time ALAG1 (h)")                            # Table 2: ALAG1 = 0.01 h (SE 21%)
    lvc   <- log(754);    label("Central volume of distribution V2 (L)")                    # Table 2: V2 = 754 L (SE 9.5%)
    lkel  <- log(0.288);  label("Elimination rate constant at 100 mg reference (kel, 1/h)") # Table 2: K(100 mg) = 0.288 /h (SE 5.7%)
    lk12  <- log(0.074);  label("Distribution rate central -> peripheral1 k12 (paper K23, 1/h)")  # Table 2: K23 = 0.074 /h (SE 6.7%)
    lk21  <- log(0.051);  label("Distribution rate peripheral1 -> central k21 (paper K32, 1/h)")  # Table 2: K32 = 0.051 /h (SE 2.2%)

    # F1 reference anchor at 100 mg b.i.d. The relative bioavailability for the
    # other dose groups (3, 10, 25, 300, 600 mg) enters in model() via the
    # e_dose_<level>_fdepot multipliers below. Held fixed at log(1) = 0 because
    # the paper anchors the 100 mg dose group F1 = 1.00 (Table 2, no SE).
    lfdepot <- fixed(log(1)); label("Relative bioavailability F1 at 100 mg b.i.d. reference (fraction)") # Table 2: F1 (100 mg) = 1.00 anchor

    # ---- Dose-dependent F1 multipliers (Table 2 F1 rows other than 100 mg) ----
    # Wrapped in fixed() because each value is a single point estimate from the
    # Table 2 dose-stratified F1 fit; reproducing the source model requires
    # holding these at their published values rather than re-estimating.
    e_dose_3mg_fdepot   <- fixed(0.139); label("F1 at 3 mg b.i.d. (fraction; multiplier on lfdepot)")     # Table 2: F1 (3 mg)   = 0.139 (SE 18%)
    e_dose_10mg_fdepot  <- fixed(0.166); label("F1 at 10 mg b.i.d. (fraction)")                            # Table 2: F1 (10 mg)  = 0.166 (SE 12%)
    e_dose_25mg_fdepot  <- fixed(0.265); label("F1 at 25 mg b.i.d. (fraction)")                            # Table 2: F1 (25 mg)  = 0.265 (SE 16%)
    e_dose_300mg_fdepot <- fixed(2.27);  label("F1 at 300 mg b.i.d. (fraction)")                           # Table 2: F1 (300 mg) = 2.27 (SE 17%)
    e_dose_600mg_fdepot <- fixed(2.50);  label("F1 at 600 mg q.d. (fraction)")                             # Table 2: F1 (600 mg) = 2.50 (SE 12%)

    # ---- Dose-dependent kel multipliers (Table 2 K rows; K(dose) / K(100 mg)) ----
    # Five fixed ratios that scale lkel per dose group. The 100 mg reference
    # implies multiplier = 1 and is not encoded as a parameter.
    e_dose_3mg_kel   <- fixed(0.104 / 0.288); label("kel ratio at 3 mg vs 100 mg reference (unitless)")   # Table 2: K(3 mg)   = 0.104 /h
    e_dose_10mg_kel  <- fixed(0.117 / 0.288); label("kel ratio at 10 mg vs 100 mg reference (unitless)")  # Table 2: K(10 mg)  = 0.117 /h
    e_dose_25mg_kel  <- fixed(0.129 / 0.288); label("kel ratio at 25 mg vs 100 mg reference (unitless)")  # Table 2: K(25 mg)  = 0.129 /h
    e_dose_300mg_kel <- fixed(0.358 / 0.288); label("kel ratio at 300 mg vs 100 mg reference (unitless)") # Table 2: K(300 mg) = 0.358 /h
    e_dose_600mg_kel <- fixed(0.376 / 0.288); label("kel ratio at 600 mg vs 100 mg reference (unitless)") # Table 2: K(600 mg) = 0.376 /h

    # =========================================================================
    # PD -- CCR5 receptor occupancy direct-effect Emax with baseline (Equation 1)
    # Occ(%) = E0 + Emax * Cp / (KD + Cp)
    # where Cp is the maraviroc plasma concentration (ng/mL), E0 is the
    # background receptor occupancy at baseline (assay floor ~25%), Emax is the
    # maximum drug-induced occupancy above baseline, and KD is the concentration
    # for 50% of Emax (NOT 50% of total occupancy).
    # =========================================================================
    lemax <- log(66.8);   label("Maximum drug-induced CCR5 receptor occupancy above baseline Emax (%)")   # Table 2: Emax = 66.8% (SE 2.4%)
    lkd   <- log(0.0894); label("Maraviroc plasma concentration for half-Emax receptor occupancy KD (ng/mL)") # Table 2: KD = 0.0894 ng/mL (SE 13%)
    le0   <- log(24.9);   label("Background CCR5 receptor occupancy at baseline E0 (%)")                   # Table 2: E0 = 24.9% (SE 6%)

    # =========================================================================
    # Inter-individual variability -- Table 2 IIV (CV%) column.
    # omega^2 = log(CV^2 + 1) for log-normal IIV.
    #   Ka  CV  89%        -> omega^2 = log(1 + 0.89^2) = 0.5831
    #   V2  CV  31%        -> omega^2 = log(1 + 0.31^2) = 0.0917
    #   K23 CV  24%        -> omega^2 = log(1 + 0.24^2) = 0.0560
    #   KD  CV  21%        -> omega^2 = log(1 + 0.21^2) = 0.0431
    #   E0  CV  28%        -> omega^2 = log(1 + 0.28^2) = 0.0754
    #   ALAG1 CV > 100%    -> approximated as 100% here (omega^2 = log(2) = 0.6931); see vignette Errata.
    # IOV on Ka (60% CV) reported in Table 2 is NOT encoded (no occasion column in the
    # nlmixr2 model interface); see vignette Errata for the dropped IOV.
    # F1 IIV is mentioned in the paper Results but no magnitude is reported in
    # Table 2; not encoded here (set to fixed(0) implicitly by omission); see
    # vignette Errata.
    # =========================================================================
    etalka   ~ 0.5831
    etalvc   ~ 0.0917
    etalk12  ~ 0.0560
    etaltlag ~ 0.6931
    etalkd   ~ 0.0431
    etale0   ~ 0.0754

    # =========================================================================
    # Residual error -- Table 2.
    # EPS1 ("multiplicative error model") on plasma maraviroc concentration:
    # 38% CV. Encoded as proportional residual on Cc (ng/mL).
    # EPS2 ("additive error model") on receptor occupancy: 11%. Encoded as
    # additive residual on Occ (% receptor occupancy).
    # =========================================================================
    propSd     <- 0.38; label("Maraviroc plasma proportional residual SD (fraction)")          # Table 2: EPS1 = 38% CV (SE 6.5%)
    addSd_Occ  <- 11;   label("CCR5 receptor occupancy additive residual SD (% receptor occupancy)") # Table 2: EPS2 = 11% (SE 9.3%)
  })

  model({
    # 1. Dose-group lookup for F1 and kel. (DOSE == d) evaluates to 1 if the
    #    subject is in dose-cohort d and 0 otherwise; the sum picks the
    #    matching multiplier. The 100 mg reference contributes the literal 1.
    fdepot_dose <- (DOSE ==   3) * e_dose_3mg_fdepot +
                   (DOSE ==  10) * e_dose_10mg_fdepot +
                   (DOSE ==  25) * e_dose_25mg_fdepot +
                   (DOSE == 100) * 1 +
                   (DOSE == 300) * e_dose_300mg_fdepot +
                   (DOSE == 600) * e_dose_600mg_fdepot

    kel_dose <- (DOSE ==   3) * e_dose_3mg_kel +
                (DOSE ==  10) * e_dose_10mg_kel +
                (DOSE ==  25) * e_dose_25mg_kel +
                (DOSE == 100) * 1 +
                (DOSE == 300) * e_dose_300mg_kel +
                (DOSE == 600) * e_dose_600mg_kel

    # 2. Individual PK parameters.
    ka   <- exp(lka  + etalka)
    tlag <- exp(ltlag + etaltlag)
    vc   <- exp(lvc  + etalvc)
    kel  <- exp(lkel)              * kel_dose          # no IIV on kel per Table 2
    k12  <- exp(lk12 + etalk12)
    k21  <- exp(lk21)                                  # no IIV on K32 per Table 2

    # 3. Individual PD parameters.
    emax <- exp(lemax)                                 # no IIV on Emax per Methods (Bmax expressed as % of maximum binding)
    kd   <- exp(lkd  + etalkd)
    e0   <- exp(le0  + etale0)

    # 4. ODE system -- two-compartment first-order absorption with lag time.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 5. Bioavailability and lag time.
    f(depot)    <- exp(lfdepot) * fdepot_dose
    alag(depot) <- tlag

    # 6. Observations.
    # central is in mg, vc is in L -> central/vc is mg/L = ug/mL; multiply by
    # 1000 to express Cc in ng/mL (consistent with the paper's reporting units
    # and the KD ng/mL parameterisation).
    Cc  <- 1000 * central / vc                                   # plasma maraviroc, ng/mL
    Occ <- e0 + emax * Cc / (kd + Cc)                            # CCR5 receptor occupancy, % (Equation 1)

    Cc  ~ prop(propSd)
    Occ ~ add(addSd_Occ)
  })
}
