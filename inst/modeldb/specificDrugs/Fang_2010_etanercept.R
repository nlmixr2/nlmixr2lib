Fang_2010_etanercept <- function() {
  description <- "One-compartment population PK model for rhTNFR-Fc (recombinant human TNF receptor-Fc fusion protein; etanercept-class molecule from Celgen Bio-Pharmaceutical) with first-order subcutaneous absorption, absorption lag time, and linear elimination in healthy Chinese volunteers (single SC doses 12.5-50 mg) and Chinese male patients with ankylosing spondylitis (multiple SC doses 25 mg BIW or 50 mg QW) (Fang 2010). Female sex is the typical-value reference: males have 0.655x lower CL/F. Single-dose administration is the typical-value reference: multi-dose administration in AS patients has 0.674x lower apparent bioavailability F."
  reference   <- "Fang Y, Li LJ, Wang R, Huang F, Song HF, Tang ZM, Li YZ, Guan HS, Zheng QS. Population pharmacokinetics of rhTNFR-Fc in healthy Chinese volunteers and in Chinese patients with Ankylosing spondylitis. Acta Pharmacologica Sinica. 2010;31(11):1500-1507. doi:10.1038/aps.2010.113"
  vignette    <- "Fang_2010_etanercept"
  units       <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (female; Fang 2010 reports CL/F = 0.168 L/h as the female-typical value)",
      notes              = "Fang 2010 encodes sex as a male-indicator (Gender = 1 male, 0 female) in the final-model equation CL/F * theta_Gender^Gender. The canonical SEXF (1 = female, 0 = male) inverts the encoding: the model applies e_sexf_cl^(1 - SEXF) so SEXF = 1 (female) yields the reference factor 1 and SEXF = 0 (male) yields the paper's male-vs-female ratio 0.655.",
      source_name        = "Gender"
    ),
    MULTI_DOSE_PT = list(
      description        = "Multiple-dose-phase-in-patients indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (single-dose / healthy-volunteer record)",
      notes              = "Fang 2010 encodes this as M (M = 1 multiple dosage, M = 0 single dosage) in the final-model equation F * theta_M^M. Single-dose records come from the 32 healthy Chinese volunteers (12.5/25/37.5/50 mg ascending single-dose cohorts); multi-dose records come from the 19 AS patients (25 mg BIW or 50 mg QW for seven consecutive doses). In Fang's design M is effectively subject-level (no subject crosses from single-dose to multi-dose), so MULTI_DOSE_PT can be set per subject from the cohort assignment.",
      source_name        = "M"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 51L,                                            # Fang 2010 Results: 32 healthy + 19 AS patients (1 of 20 enrolled AS patients excluded for missed dose)
    n_observations   = 1187L,                                          # Fang 2010 Results: 1187 plasma samples collected
    n_studies        = 1L,
    age_range        = "20-35 years (healthy 25-35, AS patients 20-31)", # Fang 2010 Patients and healthy subjects + Table 1
    weight_range     = "approximately 50-80 kg (Table 1 mean+/-SD per cohort range 60.0+/-5.3 to 66.1+/-7.7)",
    sex_female_pct   = 31.4,                                            # 16 of 51 subjects (healthy volunteers only, balanced 16M/16F; all AS patients were male)
    race_ethnicity   = "100% Chinese (Han, single-centre Chinese PLA General Hospital, Beijing)",
    disease_state    = "Pooled cohort of 32 healthy Chinese adults (16 male + 16 female) receiving single SC doses and 19 Chinese male adults with moderate-and-active ankylosing spondylitis (AS, New York criteria 1968) receiving 25 mg BIW (n = 10) or 50 mg QW (n = 9) for seven consecutive SC injections.",
    dose_range       = "Healthy volunteers: 12.5, 25, 37.5, 50 mg single SC dose (n = 8 each, 4M/4F per dose level). AS patients: 25 mg BIW (twice weekly) or 50 mg QW (once weekly) for 7 doses.",
    regions          = "China (Beijing)",
    notes            = "All injections delivered subcutaneously in the abdomen at 08:00 (before breakfast). The rhTNFR-Fc product was provided by Celgen Bio-Pharmaceutical Co Ltd (Shanghai) as a lyophilized powder; the paper notes it is the same class of molecule as etanercept (Enbrel, Amgen/Wyeth) but the studied product is a distinct Chinese manufacturer's rhTNFR-Fc, not Enbrel. Plasma rhTNFR-Fc was measured by Quantikine human sTNF RII ELISA (R&D Systems); validated linear range 48.8-3125 pg/mL on the diluted-sample standard curve."
  )

  ini({
    # ----------------------------------------------------------------------
    # Structural parameters - Fang 2010 Table 3 final-model estimates.
    # Typical values represent the reference covariate set:
    #   SEXF = 1 (female), MULTI_DOSE_PT = 0 (single dose).
    # The paper reports CL/F = 0.168 L/h as the female-typical value
    # (Discussion: "CL/F was 0.110 L/h for Chinese males and 0.168 L/h
    # for Chinese females"). Male CL/F = 0.168 * 0.655 = 0.110 L/h is
    # recovered via the e_sexf_cl^(1 - SEXF) covariate term in model().
    # ----------------------------------------------------------------------
    lka   <- log(0.0605); label("First-order SC absorption rate ka (1/h)")                    # Fang 2010 Table 3: Ka = 0.0605 1/h
    lcl   <- log(0.168);  label("Apparent clearance CL/F at SEXF = 1, MULTI_DOSE_PT = 0 (L/h)") # Fang 2010 Table 3: CL/F = 0.168 L/h (female-typical, single-dose reference)
    lvc   <- log(15.5);   label("Apparent central volume V/F (L)")                            # Fang 2010 Table 3: V/F = 15.5 L
    ltlag <- log(1.03);   label("Absorption lag time Tlag (h)")                               # Fang 2010 Table 3: Tlag = 1.03 h
    lfdepot <- fixed(log(1)); label("Bioavailability F at MULTI_DOSE_PT = 0 (fixed at 1 as SC reference)") # SC route - apparent CL/F and V/F absorb the unidentifiable absolute F; the multi-dose effect on F is captured by e_multi_dose_pt_f below.

    # Covariate effects (Fang 2010 Table 3 final-model thetas).
    e_sexf_cl       <- 0.655; label("Multiplicative male-vs-female CL/F ratio (applied as ratio^(1 - SEXF))")             # Fang 2010 Table 3: theta_Gender for CL/F = 0.655 (paper Gender = 1 - SEXF)
    e_multi_dose_pt_f <- 0.674; label("Multiplicative multi-dose-phase F ratio (applied as ratio^MULTI_DOSE_PT)")           # Fang 2010 Table 3: theta_M for F = 0.674

    # ----------------------------------------------------------------------
    # Inter-individual variability - Fang 2010 Table 3 reports IIV as %CV.
    # Translate to log-normal variance: omega^2 = log(CV^2 + 1).
    #   CL/F  33.3% CV -> log(0.333^2 + 1) = 0.1052
    #   V/F   42.7% CV -> log(0.427^2 + 1) = 0.1676
    #   Ka    55.6% CV -> log(0.556^2 + 1) = 0.2694
    #   Tlag  81.8% CV -> log(0.818^2 + 1) = 0.5121
    # No block correlations are reported in Fang 2010 Table 3.
    # ----------------------------------------------------------------------
    etalka  ~ 0.2694   # Fang 2010 Table 3: IIV Ka 55.6% CV
    etalcl  ~ 0.1052   # Fang 2010 Table 3: IIV CL/F 33.3% CV
    etalvc  ~ 0.1676   # Fang 2010 Table 3: IIV V/F 42.7% CV
    etaltlag ~ 0.5121  # Fang 2010 Table 3: IIV Tlag 81.8% CV (SE 74.4% per Table 3)

    # ----------------------------------------------------------------------
    # Residual error - Fang 2010 Table 3 combined proportional + additive.
    # ----------------------------------------------------------------------
    propSd <- 0.203;  label("Proportional residual error (fraction)")        # Fang 2010 Table 3: Proportional error = 20.3%
    addSd  <- 12.6;   label("Additive residual error (ng/mL)")                # Fang 2010 Table 3: Additive error = 12.6 ug/L (= ng/mL)
  })

  model({
    # Individual PK parameters - Fang 2010 final-model exponential
    # covariate form. The paper writes CL/F * theta_Gender^Gender; we use
    # the canonical SEXF (1 = female) storage convention with the effect
    # applied as e_sexf_cl^(1 - SEXF) so SEXF = 1 yields factor 1 and
    # SEXF = 0 yields the paper's male-vs-female ratio 0.655.
    ka  <- exp(lka  + etalka)
    cl  <- exp(lcl  + etalcl) * e_sexf_cl^(1 - SEXF)
    vc  <- exp(lvc  + etalvc)
    tlag <- exp(ltlag + etaltlag)
    fdepot <- exp(lfdepot) * e_multi_dose_pt_f^MULTI_DOSE_PT

    # Elimination rate constant for the one-compartment system.
    kel <- cl / vc

    # ODE system: depot -> central, first-order absorption + elimination.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Bioavailability and absorption-lag time apply to the depot compartment.
    f(depot)    <- fdepot
    alag(depot) <- tlag

    # Plasma concentration: doses in mg, vc in L -> central / vc in mg/L.
    # Multiply by 1000 to express Cc in ng/mL (= ug/L) to match the
    # additive-error unit reported in Fang 2010 Table 3.
    Cc <- central / vc * 1000
    Cc ~ add(addSd) + prop(propSd)
  })
}
