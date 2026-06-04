Belldina_2003_cysteamine <- function() {
  description <- "Two-compartment population PK model with first-order oral absorption and an absorption lag, sequentially linked to a one-compartment effect-site PD model with fractional inhibitory Emax (Hill = 1) for white-blood-cell cystine content reduction by cysteamine in 11 paediatric and young-adult patients (age 3-15 y, weight 14.3-60.2 kg) with nephropathic cystinosis at steady state on cysteamine bitartrate (Cystagon) approximately every 6 hours. PK and PD parameters in the source paper were estimated as individual NONMEM fits per subject and summarised as arithmetic mean / geometric mean / median / min / max across the 11 patients (Tables 2 and 3); this package encodes the arithmetic means as the typical values, with linear allometric weight scaling fixed at exponent 1.0 to reflect the paper's per-kg parameterisation of all clearance and volume terms. Dose is in mg cysteamine bitartrate salt (MW 227.24 g/mol); the model converts internally to plasma cysteamine in micromolar (free-base moiety, MW 77.15 g/mol, the measured analyte). PD output cystine is white-blood-cell cystine content in nmol cystine per mg protein."
  reference <- "Belldina EB, Huang MY, Schneider JA, Brundage RC, Tracy TS. Steady-state pharmacokinetics and pharmacodynamics of cysteamine bitartrate in paediatric nephropathic cystinosis patients. Br J Clin Pharmacol. 2003 Nov;56(5):520-525. doi:10.1046/j.1365-2125.2003.01927.x"
  vignette <- "Belldina_2003_cysteamine"
  units <- list(time = "h", dosing = "mg", concentration = "umol/L (uM)")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline. Reference 36.5 kg is the arithmetic mean weight of the 11 patients in Table 1 (43.8 + 39.0 + 57.2 + 22.5 + 18.4 + 14.3 + 29.1 + 47.0 + 60.2 + 38.7 + 31.3 = 401.5 kg / 11 = 36.5 kg). Allometric exponents are fixed at 1.0 on CL, Q, Vc, Vp to reflect the source paper's per-kg parameterisation of all PK parameters (CL/F mL/min/kg, Vc L/kg, Q mL/min/kg, Vss/F L/kg).",
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 11L,
    n_studies      = 1L,
    age_range      = "3-15 years",
    weight_range   = "14.3-60.2 kg",
    weight_median  = "38.7 kg (median of Table 1; arithmetic mean 36.5 kg)",
    height_range   = "93.0-162.0 cm",
    sex_female_pct = 36.36,
    race_ethnicity = "10 Caucasian + 1 Caucasian/Hispanic per Table 1 (single-centre UCSD cohort)",
    disease_state  = "Nephropathic (infantile) cystinosis without renal transplant; documented evidence of disease and chronic cysteamine therapy >= 12 months prior to study.",
    dose_range     = "225-550 mg oral cysteamine bitartrate (Cystagon) at steady state on approximately every-6-hours dosing.",
    administration = "Oral, with 100 mL of ambient-temperature water; a standardised low-fat / low-protein breakfast was served 30 min before dosing and a standardised lunch 4 h after dosing.",
    regions        = "United States (single centre: University of California, San Diego).",
    notes          = "Single-dose, open-label, steady-state study (Methods). Demographics in Table 1; PK estimates in Table 2; PD estimates in Table 3. Plasma cysteamine measured by HPLC-UV (LLOQ 1.3 uM); WBC cystine measured by 14C-cystine binding-protein assay per Smith et al. The 11 patients were aged 3-15 y; one parent or guardian signed informed consent for subjects under 18, and subjects over 7 also signed an informed assent form."
  )

  ini({
    # ----------------------------------------------------------------
    # Structural PK -- Table 2 arithmetic means at reference WT = 36.5 kg
    # Paper reports per-kg parameters; absolute typical values are computed
    # by multiplying per-kg means from Table 2 by the reference weight
    # 36.5 kg (arithmetic mean across N = 11 patients in Table 1).
    # ----------------------------------------------------------------
    lka   <- log(1.7);    label("Absorption rate constant ka (1/h)")                                                      # Table 2: Ka mean = 1.7 1/h (range 0.7-2.5; geomean 1.6)
    lcl   <- log(70.74);  label("Apparent oral clearance CL/F (L/h)")                                                     # Table 2: CL/F mean = 32.3 mL/min/kg * 36.5 kg * 60/1000 = 70.74 L/h
    lvc   <- log(73.0);   label("Apparent oral central volume Vc/F (L)")                                                  # Table 2: Vc mean = 2.0 L/kg * 36.5 kg = 73.0 L
    lq    <- log(65.26);  label("Apparent oral inter-compartmental clearance Q/F (L/h)")                                  # Table 2: Q mean = 29.8 mL/min/kg * 36.5 kg * 60/1000 = 65.26 L/h
    lvp   <- log(478.15); label("Apparent oral peripheral volume Vp/F (L)")                                               # Table 2: Vss/F mean = 15.1 L/kg * 36.5 kg = 551.15 L; Vp = Vss - Vc = 478.15 L
    ltlag <- log(0.44);   label("Absorption lag time tlag (h)")                                                           # Table 2: Alag mean = 0.44 h (range 0.22-0.92; geomean 0.41)

    # Allometric exponents -- fixed at 1.0 on all clearances and volumes
    # to reflect the paper's per-kg parameterisation (linear weight effect).
    e_wt_cl_q  <- fixed(1.0); label("Shared allometric exponent on CL and Q (paper per-kg = linear)")    # Table 2 per-kg parameterisation
    e_wt_vc_vp <- fixed(1.0); label("Shared allometric exponent on Vc and Vp (paper per-kg = linear)")   # Table 2 per-kg parameterisation

    # ----------------------------------------------------------------
    # PD -- Table 3 arithmetic means; effect-site (Ce) drives a
    # fractional inhibitory Emax (Hill = 1) toward zero WBC cystine.
    # ----------------------------------------------------------------
    lec50 <- log(15.3);  label("Effect-site cysteamine concentration for 50% maximal cystine reduction EC50 (uM)")        # Table 3: EC50 mean = 15.3 uM (range 0.6-61.1; geomean 5.6)
    lke0  <- log(2.2);   label("Effect-compartment equilibration rate constant ke0 (1/h)")                                # Table 3: Ke0 mean = 2.2 1/h (range 0.2-8.9; geomean 1.3)
    lbl   <- log(0.91);  label("Baseline WBC cystine content BL (nmol cystine / mg protein)")                              # Table 3: BL mean = 0.91 nmol/mg (range 0.13-1.9; geomean 0.76)

    # ----------------------------------------------------------------
    # IIV -- the source paper did NOT estimate population OMEGAs.
    # The variances below approximate the cross-individual spread of
    # the 11 individual NONMEM fits via the log-normal identity
    # omega^2 = 2 * log(arithmetic mean / geometric mean). These are
    # descriptive proxies, not formally estimated population variances;
    # see vignette Assumptions and deviations for the full rationale.
    # ----------------------------------------------------------------
    etalka   ~ 0.121  # 2*log(1.7/1.6); Table 2 Ka cross-individual spread
    etalcl   ~ 0.108  # 2*log(32.3/30.6); Table 2 CL/F cross-individual spread
    etalvc   ~ 0.446  # 2*log(2.0/1.6); Table 2 Vc cross-individual spread
    etalq    ~ 0.384  # 2*log(29.8/24.6); Table 2 Q cross-individual spread
    etalvp   ~ 0.634  # 2*log(15.1/11.0); Table 2 Vss/F cross-individual spread, used as Vp/F proxy
    etaltlag ~ 0.141  # 2*log(0.44/0.41); Table 2 Alag cross-individual spread
    etalec50 ~ 2.011  # 2*log(15.3/5.6); Table 3 EC50 cross-individual spread (very large; range spans 0.6-61.1 uM)
    etalke0  ~ 1.052  # 2*log(2.2/1.3); Table 3 Ke0 cross-individual spread
    etalbl   ~ 0.360  # 2*log(0.91/0.76); Table 3 BL cross-individual spread

    # ----------------------------------------------------------------
    # Residual error -- the source paper does not report the residual
    # SD values themselves (PK Methods specify a proportional structure
    # with no magnitude; PD Methods do not specify any error structure).
    # Values below are operator-assumed placeholders documented in the
    # vignette Assumptions and deviations section; a refit can re-estimate.
    # ----------------------------------------------------------------
    propSd        <- 0.15; label("Proportional residual error on plasma cysteamine Cc (fraction; assumed; paper Methods specify proportional structure but SD value is not reported)")
    addSd_cystine <- 0.10; label("Additive residual error on WBC cystine output (nmol cystine / mg protein; assumed; paper does not specify a residual error structure for cystine)")
  })

  model({
    # Reference body weight 36.5 kg (Table 1 arithmetic mean across N = 11).
    # Cysteamine bitartrate (Cystagon) salt MW = 227.24 g/mol; cysteamine
    # free base MW = 77.15 g/mol. The drug is administered as the bitartrate
    # salt (the dose-record amt is mg of salt), but plasma cysteamine is
    # measured as the free-base moiety. The molar ratio is 1:1 (one mol of
    # bitartrate salt liberates one mol of cysteamine), so the unit
    # conversion applied at observation is 1000 / MW_salt = 1000 / 227.24
    # = 4.402 umol cysteamine per mg salt per L.

    # Individual PK parameters with linear (per-kg) weight scaling.
    ka   <- exp(lka   + etalka)
    cl   <- exp(lcl   + etalcl)   * (WT / 36.5)^e_wt_cl_q
    vc   <- exp(lvc   + etalvc)   * (WT / 36.5)^e_wt_vc_vp
    q    <- exp(lq    + etalq)    * (WT / 36.5)^e_wt_cl_q
    vp   <- exp(lvp   + etalvp)   * (WT / 36.5)^e_wt_vc_vp
    tlag <- exp(ltlag + etaltlag)

    # Effect-site parameters (no body-size scaling; concentrations and rates).
    ke0  <- exp(lke0  + etalke0)
    ec50 <- exp(lec50 + etalec50)
    bl   <- exp(lbl   + etalbl)

    # Micro rate constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment PK with first-order oral absorption and absorption lag.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    alag(depot) <- tlag

    # Plasma cysteamine in uM: central (mg of bitartrate salt) / vc (L)
    # yields mg salt / L; convert to umol cysteamine / L via 1000 / MW_salt.
    Cc <- central / vc * 1000 / 227.24

    # Effect compartment driven by plasma Cc with first-order equilibration ke0.
    d/dt(effect) <- ke0 * (Cc - effect)

    # PD model (Belldina 2003 Eq. 1, fractional inhibitory Emax with Hill = 1):
    # E = BL * (1 - Ce / (EC50 + Ce)) = BL * EC50 / (EC50 + Ce). Effect approaches
    # zero as effect-site concentration rises and is constrained from negative.
    cystine <- bl * ec50 / (ec50 + effect)

    # Observations.
    Cc      ~ prop(propSd)
    cystine ~ add(addSd_cystine)
  })
}
