Standing_2008_diclofenac <- function() {
  description <- "One-compartment population PK model for oral diclofenac suspension in children and adult volunteers (Standing 2008): two parallel transit-absorption arms (Savic 2007 analytical input form) feeding two depot compartments that each absorb into a single central disposition compartment with linear elimination. Allometric weight scaling on clearance and volume to a 70 kg reference. Captures the double-peak absorption profile common to immediate-release diclofenac. Separate proportional residual error for paediatric and adult cohorts. Source paper additionally fits between-occasion variability (BOV) on CL/F (20%) and Vd/F (93%) -- BOV is not implemented in this nlmixr2 model file because it requires an OCC column in the user dataset; users who want BOV can add an etalcl_bov / etalvc_bov occasion-level random effect themselves."
  reference <- paste(
    "Standing JF, Howard RF, Johnson A, Savage I, Wong ICK.",
    "Population pharmacokinetics of oral diclofenac for acute pain in children.",
    "Br J Clin Pharmacol. 2008;66(6):846-853.",
    "doi:10.1111/j.1365-2125.2008.03289.x.",
    sep = " "
  )
  vignette <- "Standing_2008_diclofenac"
  units <- list(time = "hour", dosing = "nmol", concentration = "nmol/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used for allometric scaling with 70 kg reference: CL ~ (WT/70)^0.75, V ~ (WT/70)^1. Standing 2008 Methods page 847, equations under 'Pharmacokinetic model building'. Cohort weight range 9-94 kg (Table 1).",
      source_name        = "WT"
    ),
    CHILD = list(
      description        = "Paediatric cohort indicator (1 = paediatric day-surgery patient, 0 = adult volunteer).",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Used only to switch proportional residual error between cohorts. Standing 2008 estimated separate residual error for the two datasets because the two cohorts were analysed with different bioanalytical methods (HPLC/MS with ketoprofen internal standard on paediatric serum vs naproxen internal standard on adult plasma; assay LLOQs 10.1 ng/mL vs 10 ng/mL). Paper Methods page 848 'Pharmacokinetic model building'.",
      source_name        = "CHILD"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 100L,
    n_studies      = 2L,
    age_range      = "1-28 years (children 1-12 years pooled with adults 18-28 years)",
    age_median     = "3 years (paediatric), 21 years (adult); 9 years pooled (Table 1)",
    weight_range   = "9-94 kg",
    weight_median  = "17 kg (paediatric), 72 kg (adult); 34 kg pooled (Table 1)",
    sex_female_pct = 45,
    race_ethnicity = "Not reported separately for the modelling cohort",
    disease_state  = "Paediatric day-surgery patients (dermatology, general, plastic surgery) plus healthy adult volunteers",
    dose_range     = "Single oral dose: 1 mg/kg diclofenac sodium suspension (paediatric, rounded to nearest 5 mg); 50 mg diclofenac sodium suspension (adult); formulation 50 mg/5 mL Rosemont Pharmaceuticals oral suspension",
    regions        = "United Kingdom (Great Ormond Street Hospital, London) and Ireland (Shandon Clinic, Cork; adult bioequivalence cohort)",
    notes          = "Pooled analysis of 558 serum/plasma diclofenac concentrations: 206 from 70 paediatric patients (sparse, 3 samples per dose) plus 352 from 30 adult volunteers (rich, 14 samples per dose). Demographics per Standing 2008 Table 1. Paediatric concentrations are serum, adult concentrations are plasma. NONMEM v6 FOCEI with interaction. Mass units expressed as nanomoles (MW diclofenac sodium = 318.13 g/mol for dose, MW diclofenac free acid = 296.15 g/mol for measured concentration)."
  )

  ini({
    # Structural parameters. All from Standing 2008 Table 2 ('NONMEM
    # parameter estimates from final model'), with values standardised
    # to a 70 kg reference subject by the a priori allometric weight
    # scaling described on page 848 of Methods.

    # Arm 1 (fast absorption): transit chain into depot1, then ka1
    # into central. Paper Table 2 row order.
    lmtt1 <- log(0.68);  label("Arm-1 mean transit time (h)")            # Table 2: MTT1 = 0.68 h (RSE 11.8%)
    ln1   <- log(1.03);  label("Arm-1 number of transit compartments (unitless)") # Table 2: N1 = 1.03 (RSE 28.6%); estimated as a continuous quantity via the gamma-density analytical input
    # F1 = 0.70 = fraction of dose routed via arm 1 (= fdepot here).
    # F2 = 1 - F1 = 0.30 is held to the F1 + F2 = 1 constraint per
    # Methods page 848 ('Bioavailability ... with limits forcing the
    # combined bioavailability to equal 100%'); F2 is therefore not a
    # separately estimated parameter.
    lfdepot <- log(0.70); label("Arm-1 bioavailability fraction (F1, unitless)") # Table 2: F1 = 0.70 (RSE 7.6%); F2 fixed at 1 - F1
    # Paper reports t1/2,A1 = 0.09 h; ka1 = log(2)/t1/2,A1 = 7.70 1/h.
    lka1  <- log(log(2) / 0.09); label("Arm-1 absorption rate constant ka1 (1/h)")  # Table 2: t1/2A1 = 0.09 h (RSE 50.1%); ka1 = ln(2)/t1/2A1

    # Arm 2 (slow absorption): transit chain into depot2, then ka2
    # into central.
    lmtt2 <- log(1.37);  label("Arm-2 mean transit time (h)")            # Table 2: MTT2 = 1.37 h (RSE 6.97%)
    ln2   <- log(41.60); label("Arm-2 number of transit compartments (unitless)") # Table 2: N2 = 41.60 (RSE 73.6%)
    # Paper reports t1/2,A2 = 1.06 h; ka2 = log(2)/t1/2,A2 = 0.654 1/h.
    lka2  <- log(log(2) / 1.06); label("Arm-2 absorption rate constant ka2 (1/h)")  # Table 2: t1/2A2 = 1.06 h (RSE 12.2%); ka2 = ln(2)/t1/2A2

    # Disposition. Per 70 kg standardisation per the allometric size
    # model in Methods (page 848).
    lvc <- log(4.84);  label("Apparent central volume Vd/F per 70 kg (L)") # Table 2: Vd/F = 4.84 L/70 kg (RSE 86.2%)
    lcl <- log(53.98); label("Apparent clearance CL/F per 70 kg (L/h)")    # Table 2: CL/F = 53.98 L/h/70 kg (RSE 4.8%)

    # Allometric exponents fixed a priori per Methods page 848
    # (Anderson and Holford size-and-maturity framework; Kleiber 3/4
    # power rule for CL). Reported without uncertainty so encoded as
    # fixed.
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL (unitless)") # Methods page 848: CL = thetaCL * (W/70)^0.75
    e_wt_vc <- fixed(1.00); label("Allometric exponent on Vc (unitless)") # Methods page 848: V  = thetaV  * (W/70)^1

    # Inter-individual variability. Standing 2008 Table 2 reports
    # IIV as CV%; the internal variance is omega^2 = log(1 + CV^2).
    # No $OMEGA BLOCK is mentioned in the paper, so IIVs are
    # independent (diagonal omega).
    etalmtt1   ~ 0.5142   # Table 2 IIV(MTT1)   = 82%  CV; log(1 + 0.82^2) = log(1.6724) = 0.5142
    etaln1     ~ 0.7131   # Table 2 IIV(N1)     = 102% CV; log(1 + 1.02^2) = log(2.0404) = 0.7131
    etalfdepot ~ 0.0560   # Table 2 IIV(F1)     = 24%  CV; log(1 + 0.24^2) = log(1.0576) = 0.0560
    etalka1    ~ 0.0917   # Table 2 IIV(t1/2A1) = 31%  CV; equivalent on log(ka1) scale
    etalmtt2   ~ 0.8624   # Table 2 IIV(MTT2)   = 117% CV; log(1 + 1.17^2) = log(2.3689) = 0.8624
    etaln2     ~ 1.1509   # Table 2 IIV(N2)     = 147% CV; log(1 + 1.47^2) = log(3.1609) = 1.1509
    etalka2    ~ 0.2152   # Table 2 IIV(t1/2A2) = 49%  CV; equivalent on log(ka2) scale
    etalvc     ~ 0.2562   # Table 2 IIV(Vd/F)   = 54%  CV; log(1 + 0.54^2) = log(1.2916) = 0.2562
    etalcl     ~ 0.0654   # Table 2 IIV(CL/F)   = 26%  CV; log(1 + 0.26^2) = log(1.0676) = 0.0654

    # Residual proportional error stratified by cohort (Methods page 848:
    # 'estimation of residual variability was undertaken separately for
    # the two groups' due to different assay methods).
    propSdAdult <- 0.29; label("Adult proportional residual SD (fraction)")       # Table 2 residual variability (adult)        = 29% (RSE 5.9%)
    propSdPaed  <- 0.18; label("Paediatric proportional residual SD (fraction)")  # Table 2 residual variability (paediatric)   = 18% (RSE 19.9%)
  })

  model({
    # Molecular-weight constant for context (not used in computation;
    # the user-supplied dose is already in nmol). Diclofenac sodium
    # MW = 318.13 g/mol; diclofenac free acid MW = 296.15 g/mol
    # (Methods page 848). Vignette demonstrates the mg -> nmol
    # conversion at the data-construction step.

    # Individual structural parameters (back-transformed from the log
    # scale, with allometric weight scaling on CL and V centred at
    # 70 kg per Methods page 848).
    mtt1   <- exp(lmtt1   + etalmtt1)
    n1     <- exp(ln1     + etaln1)
    fdepot <- exp(lfdepot + etalfdepot)
    ka1    <- exp(lka1    + etalka1)
    mtt2   <- exp(lmtt2   + etalmtt2)
    n2     <- exp(ln2     + etaln2)
    ka2    <- exp(lka2    + etalka2)
    vc     <- exp(lvc     + etalvc) * (WT / 70)^e_wt_vc
    cl     <- exp(lcl     + etalcl) * (WT / 70)^e_wt_cl

    kel <- cl / vc

    # Dual-arm Savic transit absorption (paper Figure 2 schematic):
    # arm 1 (transit chain into depot1, ka1 -> central, fraction
    # F1 = fdepot) and arm 2 (transit chain into depot2, ka2 ->
    # central, fraction F2 = 1 - F1). The user dataset must include
    # two dose events per administration -- one targeting depot1 and
    # one targeting depot2 with the SAME AMT -- so each transit()
    # call sees its own podo() / tad() context. The fdepot and
    # 1 - fdepot fractions enter via the bio argument of transit() so
    # the cumulative gamma-kernel input into depot1 is F1 * AMT and
    # into depot2 is F2 * AMT (sum = AMT, the F1 + F2 = 1 constraint
    # of Methods page 848). f(depot1) = f(depot2) = 0 suppress the
    # default bolus deposition so all mass enters via the transit
    # chains (matches Lee 2015 supplement convention; see
    # Lee_2015_sumatriptan.R for the dual-arm dosing pattern this
    # model inherits). See vignette for the canonical event-table
    # construction.
    d/dt(depot1)  <- transit(n1, mtt1, fdepot)     - ka1 * depot1
    d/dt(depot2)  <- transit(n2, mtt2, 1 - fdepot) - ka2 * depot2
    d/dt(central) <- ka1 * depot1 + ka2 * depot2 - kel * central

    f(depot1) <- 0  # suppress bolus deposition; transit() carries F1 * AMT
    f(depot2) <- 0  # suppress bolus deposition; transit() carries F2 * AMT

    # Concentration in central (nmol/L). The paper expresses
    # concentrations in nmol/L using the diclofenac free-acid MW =
    # 296.15 g/mol (Methods page 848); the user supplies AMT
    # converted to nmol so V_D is already in compatible units.
    Cc <- central / vc

    # Cohort-stratified proportional residual error (Methods page 848).
    # CHILD = 1 for paediatric records, CHILD = 0 for adult records.
    propSd <- propSdPaed * CHILD + propSdAdult * (1 - CHILD)
    Cc ~ prop(propSd)
  })
}
