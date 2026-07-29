Schaefer_1996_ciprofloxacin <- function() {
  description <- paste(
    "Two-compartment population PK model for intravenous and oral",
    "ciprofloxacin in 10 pediatric cystic fibrosis patients aged 6-16",
    "years (Schaefer 1996). First-order absorption from a depot, two",
    "disposition compartments, plus a cumulative urine compartment driven",
    "by an independently estimated renal clearance. Total clearance is a",
    "linear-with-intercept function of body weight (CL = 8.8 + 0.396 *",
    "WT, L/h), and central / peripheral volumes are directly proportional",
    "to body weight with slopes 0.698 and 1.3 L/kg respectively.",
    "Intercompartmental clearance, absorption rate, renal clearance, and",
    "oral bioavailability are weight-independent. IIV is retained on CL,",
    "Vc, and Vp only; residual error is proportional. Calibrated to the",
    "study weight range (15-42 kg); extrapolation beyond is not",
    "appropriate per the source authors."
  )
  reference <- paste(
    "Schaefer HG, Stass H, Wedgwood J, Hampel B, Fischer C, Kuhlmann J,",
    "Schaad UB. Pharmacokinetics of ciprofloxacin in pediatric cystic",
    "fibrosis patients. Antimicrob Agents Chemother. 1996 Jan;40(1):29-34.",
    "doi:10.1128/aac.40.1.29."
  )
  vignette <- "Schaefer_1996_ciprofloxacin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight at study entry",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at baseline in Schaefer 1996 (study window 24-48 h).",
        "Schaefer 1996 weight range 14.9-42.0 kg (mean 27.5 kg).",
        "Enters CL as a linear additive slope (CL = 8.8 + 0.396 * WT)",
        "and enters Vc and Vp as a proportional scaling at reference",
        "weight 70 kg with fixed exponent 1 (Vc = 0.698 * WT, Vp = 1.3",
        "* WT). The source authors note that extrapolation outside",
        "15-42 kg is not appropriate."
      ),
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 10L,
    n_studies      = 1L,
    age_range      = "5.9-15.7 years",
    age_median     = "10.1 years (mean; Table 1 gives 9.3-15.7 for the 8 oral-dose patients)",
    weight_range   = "14.9-42.0 kg",
    weight_median  = "27.5 kg (mean across N=10)",
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state  = paste(
      "Pediatric cystic fibrosis patients who had completed a standard",
      "2-week induction course of intravenous ceftazidime + amikacin",
      "(with inhaled amikacin and physiotherapy)."
    ),
    dose_range     = paste(
      "Each patient received two 30-min IV infusions of 10 mg/kg",
      "ciprofloxacin (maximum 400 mg) 12 h apart, followed by oral",
      "ciprofloxacin 15 mg/kg every 12 h. Two of the 10 patients did",
      "not receive the oral dose (Patients 7 and 9, Table 1)."
    ),
    regions        = "Switzerland (University of Berne, Pediatric Infectious Diseases) and Germany (Bayer AG sponsor).",
    notes          = paste(
      "Sparse-rich popPK study; 232 ciprofloxacin concentrations (203",
      "plasma + 29 urine) analysed in NONMEM IV with first-order",
      "conditional estimation. Plasma protein binding approximately 34%",
      "(Table 2; mean across 1-20 mg/L). Demographics and dosing in",
      "Table 1; final parameter estimates in Table 4."
    )
  )

  ini({
    # Structural parameters - reference values from Schaefer 1996 Table 4
    # (final model). Volumes encoded at a 70 kg adult reference for
    # uniformity with other paediatric models in the registry; the
    # paper's directly-proportional Vc = 0.698 * WT and Vp = 1.3 * WT are
    # recovered by exp(lvc) * (WT / 70)^1 and exp(lvp) * (WT / 70)^1.

    lka       <- log(0.644);        label("Absorption rate constant (Ka, 1/h)")                                              # Table 4 theta5 = 0.644 /h (SE 25.3%)
    lcl       <- log(8.8);          label("Weight-independent intercept term of total CL (L/h)")                             # Table 4 theta1 = u1 = 8.8 L/h (SE 10.0%); CL = u1 + u8 * WT
    lvc       <- log(0.698 * 70);   label("Central volume Vc at 70 kg reference (L); Vc = 0.698 * WT")                       # Table 4 theta2 = u2 = 0.698 L/kg (SE 12.7%); V2 = u2 * WT
    lvp       <- log(1.3   * 70);   label("Peripheral volume Vp at 70 kg reference (L); Vp = 1.3 * WT")                      # Table 4 theta3 = u3 = 1.3 L/kg (SE 9.5%); V3 = u3 * WT
    lq        <- log(21.0);         label("Intercompartmental clearance Q (L/h)")                                            # Table 4 theta4 = u4 = 21.0 L/h (SE 35.4%)
    lcl_renal <- log(11.4);         label("Renal clearance (L/h)")                                                           # Table 4 theta6 = u6 = CLR = 11.4 L/h (SE 15.1%)
    lfdepot   <- log(0.618);        label("Oral bioavailability (fraction)")                                                 # Table 4 theta7 = u7 = F1 = 0.618 (SE 5.2%)

    # Weight-effect parameters. e_wt_cl is the linear additive slope of
    # CL on WT (L/h/kg), not a fractional / power exponent: this matches
    # the paper's CL = u1 + u8 * WT regression. e_wt_vc and e_wt_vp are
    # fixed allometric exponents = 1, reproducing the paper's Vc = u2 *
    # WT and Vp = u3 * WT proportional scaling.
    e_wt_cl   <- 0.396;             label("Linear WT slope on CL (L/h/kg)")                                                  # Table 4 theta8 = u8 = 0.396 L/h/kg (SE 10.9%)
    e_wt_vc   <- fixed(1);          label("Allometric exponent on Vc (unitless)")                                            # Schaefer 1996 final model: V2 proportional to WT (exponent = 1, fixed)
    e_wt_vp   <- fixed(1);          label("Allometric exponent on Vp (unitless)")                                            # Schaefer 1996 final model: V3 proportional to WT (exponent = 1, fixed)

    # IIV (between-subject variability). Paper text reports the source
    # form CL_i = TVCL * (1 + h_CL_i) with vCL = 7.8% CV (variance
    # approx CV^2 for the small CVs reported). Translated to canonical
    # log-normal form exp(lX + etaX) using omega^2 = log(CV^2 + 1).
    # Diagonal-only: Table 3 row "Final model" notes "Off-diagonal
    # elements in covariance matrix allowed" but Table 4 reports no
    # numeric off-diagonal values; documented as a deviation in the
    # vignette Assumptions section.
    etalcl ~ 0.006061  # Table 4 vCL = 7.8% CV  -> log(1 + 0.078^2) = 0.006061
    etalvc ~ 0.049837  # Table 4 vV2 = 22.6% CV -> log(1 + 0.226^2) = 0.049837
    etalvp ~ 0.021957  # Table 4 vV3 = 14.9% CV -> log(1 + 0.149^2) = 0.021957

    # Residual error. Paper form Y = F * (1 + e), epsilon ~ N(0,
    # sigma^2) with sigma = 0.319 (Table 4 d = 31.9% CV).
    propSd <- 0.319;                label("Proportional residual error (fraction)")                                          # Table 4 d = 31.9% (SE 16.8%)
  })

  model({
    # Individual structural parameters.
    # Total CL: linear additive WT term (Schaefer u1 + u8 * WT).
    # Log-normal IIV applied multiplicatively on total CL via exp(etalcl).
    cl_typ <- exp(lcl) + e_wt_cl * WT
    cl     <- cl_typ * exp(etalcl)

    # Vc, Vp scale proportionally with WT (e_wt_vc = e_wt_vp = 1, fixed)
    # at 70 kg reference. Log-normal IIV via exp(eta).
    vc       <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    vp       <- exp(lvp + etalvp) * (WT / 70)^e_wt_vp

    # Weight-independent parameters: Ka, Q, CL_renal, F1.
    ka       <- exp(lka)
    q        <- exp(lq)
    cl_renal <- exp(lcl_renal)

    # Micro-constants.
    kel       <- cl       / vc
    kel_renal <- cl_renal / vc
    k12       <- q        / vc
    k21       <- q        / vp

    # ODE system: 2-compartment disposition with first-order absorption
    # and a cumulative urine compartment driven by the renal clearance.
    # Total elimination from the central compartment uses kel = cl / vc;
    # the urine compartment additionally accumulates the renal portion
    # (CLR * Cc) without subtracting it again from central (CLR is a
    # subset of total CL, not an additive parallel arm).
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1
    d/dt(urine)       <-  kel_renal * central

    # Oral bioavailability applied to the depot.
    f(depot) <- exp(lfdepot)

    # Observation: plasma concentration in mg/L (dose in mg, vc in L).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
