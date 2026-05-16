Woillard_2011_tacrolimus <- function() {
  description <- "Two-compartment population PK model with Erlang-distributed transit absorption (3 transit compartments) for oral tacrolimus in adult renal transplant recipients pooled across the twice-daily immediate-release Prograf formulation and the once-daily prolonged-release Advagraf formulation (Woillard 2011), with multiplicative CYP3A5*1-carrier (expresser) and power-scaled haematocrit effects on apparent clearance, multiplicative formulation effects on the Erlang transit rate constant and on apparent central volume, and combined additive plus proportional residual error."
  reference <- paste(
    "Woillard JB, de Winter BCM, Kamar N, Marquet P, Rostaing L, Rousseau A.",
    "Population pharmacokinetic model and Bayesian estimator for two tacrolimus",
    "formulations -- twice daily Prograf and once daily Advagraf.",
    "Br J Clin Pharmacol 2011; 71(3):391-402.",
    "doi:10.1111/j.1365-2125.2010.03837.x.",
    sep = " "
  )
  vignette <- "Woillard_2011_tacrolimus"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    HCT = list(
      description        = "Haematocrit, expressed as a percentage (packed red-blood-cell volume fraction times 100).",
      units              = "%",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying within subject across the 1-week / 2-week / 1-month / 3-month / 6-month post-transplant Prograf profiles; reported once per subject for the Advagraf cohort (each profile collected >= 12 months post-transplant). Power-law effect on CL/F with reference HCT = 35 % per Woillard 2011 Table 4 covariate equation (`CL/F = theta3 * (HCT / 35)^theta4 * theta5^CYP`). Study-population median HCT differed between formulations: 32.3 % (range 20.9-46.6) for Prograf and 38.5 % (range 26.5-45.1) for Advagraf (Table 1).",
      source_name        = "HT"
    ),
    CYP3A5_EXPR = list(
      description        = "CYP3A5 expresser indicator: 1 if the patient carries at least one functional CYP3A5*1 allele (genotype *1/*1 or *1/*3), 0 if homozygous *3/*3.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP3A5 *3/*3 nonexpresser)",
      notes              = "Time-fixed germline genotype derived from the rs776746 (CYP3A5 6986A>G) polymorphism by TaqMan allelic discrimination assay. The *1 (A) allele encodes functional CYP3A5; the *3 (G) allele creates a cryptic splice site and yields nonfunctional protein. In the Woillard 2011 pooled cohort the genotype distribution was *1/*1 = 1 (1 Advagraf), *1/*3 = 5 (1 Prograf + 4 Advagraf), *3/*3 = 67 (31 Prograf + 36 Advagraf) per Table 1, giving 6 expressers (CYP3A5_EXPR = 1) and 67 nonexpressers (CYP3A5_EXPR = 0). Multiplicative effect on CL/F as `theta5^CYP3A5_EXPR` with `theta5 = 2.00` (Table 4), so expressers have a 2-fold higher apparent clearance; the Discussion confirms `CL/F = 42 L/h` for *1 carriers versus `CL/F = 21 L/h` for *3/*3 carriers.",
      source_name        = "CYP"
    ),
    FORM_TAC_IR = list(
      description        = "Tacrolimus formulation indicator: 1 if the subject received the twice-daily immediate-release Prograf capsule formulation, 0 if the subject received the once-daily prolonged-release Advagraf capsule formulation.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Advagraf prolonged-release, once-daily; the Vc and Ktr typical-value reference in Woillard 2011 Table 4)",
      notes              = "Per-subject (regimen-fixed) categorical covariate; matches the Woillard 2011 `study` covariate (study = 1 for the Prograf cohort, study = 0 for the Advagraf cohort) used as a surrogate for drug formulation. Prograf is the older immediate-release tacrolimus formulation (administered every 12 hours) and Advagraf is the prolonged-release formulation (administered once daily). The Woillard 2011 covariate model applies the formulation indicator multiplicatively to two structural parameters per Table 4: `Ktr = theta1 * theta2^FORM_TAC_IR` with `theta1 = 3.34/h` and `theta2 = 1.53` (Prograf absorption is ~53 % faster than Advagraf) and `Vc/F = theta6 * theta7^FORM_TAC_IR` with `theta6 = 486 L` and `theta7 = 0.29` (Prograf apparent central volume is 29 % of the Advagraf reference). The Discussion notes the formulation-effect surrogate also partly captures the time-post-transplant difference between the two cohorts (Prograf cohort observed within the first 6 months post-transplant, Advagraf cohort > 12 months post-transplant).",
      source_name        = "study"
    )
  )

  population <- list(
    species                = "human",
    n_subjects             = 73L,
    n_studies              = 2L,
    n_profiles             = 186L,
    age_range              = "18-77 years (overall range across both cohorts)",
    age_median             = "Prograf 55 years (18-69); Advagraf 53 years (28-77)",
    weight_range           = "45-116 kg (overall range across both cohorts)",
    weight_median          = "Prograf 65 kg (46-97); Advagraf 69 kg (45-116)",
    sex_female_pct         = 47.9,
    race_ethnicity         = "Not reported in source paper (single-country French cohort).",
    disease_state          = "Adult renal transplant recipients. Prograf cohort: de novo recipients sampled at weeks 1 and 2 and months 1, 3, 6 post-transplantation on a standardized triple immunosuppressive regimen of tacrolimus + mycophenolate mofetil + tapered oral prednisolone; tacrolimus initial dose 0.1 mg/kg/day adjusted by trough TDM (10-15 ng/mL first 6 weeks, then 5-10 ng/mL). Advagraf cohort: stable recipients (> 12 months post-transplant) converted from ciclosporin to once-daily Advagraf > 6 months prior, on mycophenolate mofetil + low-dose prednisolone; tacrolimus initial dose 0.2 mg/kg/day adjusted by trough TDM.",
    dose_range             = "Tacrolimus 0.5-10 mg/day (Prograf median 4 mg, Advagraf median 4 mg per Table 1). Doses titrated to trough target concentrations.",
    regions                = "France (Limoges INSERM U850 and Toulouse Department of Nephrology-Dialysis).",
    cyp3a5_distribution    = "*1/*1 n = 1 (1 Advagraf), *1/*3 n = 5 (1 Prograf + 4 Advagraf), *3/*3 n = 67 (31 Prograf + 36 Advagraf) per Table 1; 6 expressers (CYP3A5_EXPR = 1) and 67 nonexpressers (CYP3A5_EXPR = 0). Hardy-Weinberg equilibrium confirmed.",
    formulation_breakdown  = "Prograf (twice-daily immediate-release): n = 32 subjects, 145 PK profiles. Advagraf (once-daily prolonged-release): n = 41 subjects, 41 PK profiles.",
    sampling_schedule      = "Prograf profiles: 11 blood samples per profile at pre-dose, 0.33, 0.66, 1, 1.5, 2, 3, 4, 6, 9 hours post-dose, plus a 12-hour post-dose sample at W1 and W2 (5 profiles per patient). Advagraf profiles: 12 blood samples per profile at pre-dose and 0.33, 0.66, 1, 1.5, 2, 3, 4, 6, 9, 12, 24 hours post-dose (1 profile per patient).",
    baseline_demographics  = "Median (range) by formulation per Table 1. Haemoglobin: Prograf 10.6 g/dL (6.5-15.7), Advagraf 12.9 g/dL (10.5-15.1). Serum creatinine: Prograf 119 umol/L (63-928), Advagraf 114 umol/L (82-907). Daily prednisolone dose: Prograf 20 mg (0-94), Advagraf 2.5 mg (0-10).",
    bioanalytical          = "Whole-blood tacrolimus quantified by turbulent-flow LC-MS/MS (Cyclone P online extraction, Propel C18 analytical column, TSQ Quantum Discovery MS/MS) with a 1 ng/mL lower limit of quantification.",
    notes                  = "Two patients in the Advagraf cohort had missing haematocrit and haemoglobin values replaced by the cohort median (33.7 % and 11.2 g/dL respectively). For 7 Prograf patients (21 profiles) the pre-dose value was used as a surrogate for the 12-hour post-dose concentration when computing the reference trapezoidal AUC(0,12h)."
  )

  ini({
    # Final non-mixture model estimates from Woillard 2011 Table 4 "Final model
    # obtained in the whole dataset" column (n = 73 subjects, 186 PK profiles).
    # The mixture-model alternative (Table 3) was rejected on the basis of
    # higher OFV, higher shrinkage, and similar IPV / IOV / residual error;
    # the formulation covariate `study` adequately explained the bimodal Ktr
    # distribution without invoking subpopulations (Final model and validation
    # section). The reference subject is an Advagraf-treated (FORM_TAC_IR = 0)
    # CYP3A5 nonexpresser (CYP3A5_EXPR = 0) with HCT = 35 %.

    # Erlang absorption rate constant (transit-chain rate). Multiplicative
    # formulation effect: Prograf (immediate release) absorbs ~53 % faster
    # than Advagraf (prolonged release).
    lktr <- log(3.34); label("Erlang transit absorption rate constant Ktr for the Advagraf reference (1/h)")          # Table 4 whole-dataset theta1 = 3.34 1/h

    # Formulation multiplier on Ktr: Prograf Ktr = 3.34 * 1.53 = 5.11 1/h.
    e_form_tac_ir_ktr <- 1.53; label("Formulation multiplier on Ktr for Prograf vs Advagraf (theta2)")                  # Table 4 whole-dataset theta2 = 1.53

    # Apparent oral clearance reference (CYP3A5 nonexpresser at HCT = 35 %).
    lcl <- log(21.2); label("Apparent oral clearance CL/F at HCT = 35 %, CYP3A5 nonexpresser (L/h)")                    # Table 4 whole-dataset theta3 = 21.2 L/h

    # Haematocrit power-law exponent on CL/F (negative -> CL/F drops as HCT
    # rises; low HCT lowers the red-blood-cell-bound tacrolimus fraction and
    # raises the free plasma fraction available to the liver).
    e_hct_cl <- -1.14; label("Haematocrit power-law exponent on CL/F (theta4; reference HCT = 35 %)")                   # Table 4 whole-dataset theta4 = -1.14

    # CYP3A5 expresser multiplicative factor on CL/F. Expressers (*1 carriers)
    # have ~2-fold higher CL/F than nonexpressers (*3/*3): 42 vs 21 L/h.
    e_cyp3a5_expr_cl <- 2.00; label("CYP3A5 expresser multiplicative factor on CL/F (theta5)")                          # Table 4 whole-dataset theta5 = 2.00

    # Apparent inter-compartmental clearance.
    lq <- log(79); label("Apparent inter-compartmental clearance Q/F (L/h)")                                            # Table 4 whole-dataset Q/F = 79 L/h

    # Apparent central volume of distribution reference (Advagraf). Prograf
    # subjects have Vc/F = 486 * 0.29 = 141 L. The Discussion text reports
    # 205 L vs 527 L for Prograf vs Advagraf -- these alternative numbers do
    # not match the Table 4 equation `Vc/F = theta6 * theta7^study` and are
    # not used in the model file; the vignette flags the inconsistency. The
    # values in `ini()` are the Table 4 estimates.
    lvc <- log(486); label("Apparent central volume Vc/F for the Advagraf reference (L)")                               # Table 4 whole-dataset theta6 = 486 L

    # Formulation multiplier on Vc/F (Prograf relative to Advagraf reference).
    e_form_tac_ir_vc <- 0.29; label("Formulation multiplier on Vc/F for Prograf vs Advagraf (theta7)")                  # Table 4 whole-dataset theta7 = 0.29

    # Apparent peripheral volume of distribution. Note: in the covariate-free
    # model (Table 2) Vp/F was held fixed at 500 L, but in the final
    # covariate model (Table 4) Vp/F is reported with a standard error (SE 7)
    # and an IPV estimate, so it is estimated and is not wrapped in fixed().
    lvp <- log(271); label("Apparent peripheral volume Vp/F (L)")                                                       # Table 4 whole-dataset Vp/F = 271 L

    # Inter-patient variability (IPV) reported as %CV. Conversion to
    # log-scale variance via omega^2 = log(1 + CV^2):
    #   Ktr   CV 24 % -> log(1 + 0.24^2) = 0.05605
    #   CL/F  CV 28 % -> log(1 + 0.28^2) = 0.07546
    #   Vc/F  CV 31 % -> log(1 + 0.31^2) = 0.09161
    #   Q/F   CV 54 % -> log(1 + 0.54^2) = 0.25593
    #   Vp/F  CV 60 % -> log(1 + 0.60^2) = 0.30748
    # The source paper also reported inter-occasion variability (IOV) for
    # Ktr, Vc/F, and CL/F in the Prograf cohort (which had multiple
    # occasions per patient). nlmixr2lib popPK conventions encode only IPV;
    # IOV is documented in the vignette's Assumptions and deviations
    # section and is not carried in this `ini()`.
    etalktr ~ 0.05605                                                                                                   # Table 4 whole-dataset IPV Ktr = 24 % CV
    etalcl  ~ 0.07546                                                                                                   # Table 4 whole-dataset IPV CL/F = 28 % CV
    etalvc  ~ 0.09161                                                                                                   # Table 4 whole-dataset IPV Vc/F = 31 % CV
    etalq   ~ 0.25593                                                                                                   # Table 4 whole-dataset IPV Q/F = 54 % CV
    etalvp  ~ 0.30748                                                                                                   # Table 4 whole-dataset IPV Vp/F = 60 % CV

    # Combined additive + proportional residual error (Table 4 footer,
    # whole-dataset column).
    propSd <- 0.113; label("Proportional residual error (fraction)")                                                    # Table 4 whole-dataset proportional error = 11.3 %
    addSd  <- 0.71;  label("Additive residual error (ng/mL)")                                                           # Table 4 whole-dataset additive error = 0.71 ng/mL
  })

  model({
    # Individual PK parameters. Formulation enters Ktr and Vc/F as a
    # multiplicative factor (theta_X)^FORM_TAC_IR per Woillard 2011 Table 4;
    # haematocrit enters CL/F as a power-law deviation from 35 %; CYP3A5
    # expresser status enters CL/F as a multiplicative factor (theta5)^CYP.
    ktr <- exp(lktr + etalktr) * e_form_tac_ir_ktr ^ FORM_TAC_IR
    cl  <- exp(lcl  + etalcl)  * (HCT / 35) ^ e_hct_cl * e_cyp3a5_expr_cl ^ CYP3A5_EXPR
    vc  <- exp(lvc  + etalvc)  * e_form_tac_ir_vc ^ FORM_TAC_IR
    q   <- exp(lq   + etalq)
    vp  <- exp(lvp  + etalvp)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment disposition with three-compartment Erlang transit
    # absorption (Methods, "Erlang distribution ADVAN5 SS5 ... three
    # transit compartments"). The dose enters `depot`; depot, transit1, and
    # transit2 form the three sequential transit compartments connected by
    # the common rate constant ktr; transit2 empties into central at rate
    # ktr. Bioavailability is implicit in CL/F, Q/F, Vc/F, Vp/F.
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot     - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1  - ktr * transit2
    d/dt(central)     <-  ktr * transit2  - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central   - k21 * peripheral1

    # Tacrolimus is reported in whole blood in ng/mL. Dose units mg and
    # volumes L give central / vc in mg/L; multiply by 1000 to express
    # the prediction in ng/mL.
    Cc <- central / vc * 1000
    Cc ~ add(addSd) + prop(propSd)
  })
}
