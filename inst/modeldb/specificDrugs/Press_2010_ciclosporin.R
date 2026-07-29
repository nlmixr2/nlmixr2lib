Press_2010_ciclosporin <- function() {
  description <- paste(
    "Two-compartment population pharmacokinetic model for oral ciclosporin",
    "A (Neoral) in adult kidney transplant recipients (Press 2010).",
    "Delayed absorption is described by one transit compartment with the",
    "first-order transit rate constant set equal to the absorption rate",
    "constant ka (chain: depot -> transit1 -> central at common rate ka;",
    "mean absorption time = (n+1)/ka with n = 1 transit compartment).",
    "Oral bioavailability is FIXED at 0.5 (Methods 'Structural model').",
    "Apparent clearance CL and apparent central volume of distribution Vc",
    "are allometrically scaled to body weight at a 76 kg median reference",
    "with theory-based exponents 0.75 on CL and 1.0 on Vc; the peripheral",
    "volume Vp and intercompartmental clearance Q are not weight-scaled.",
    "Concomitant high-dose oral prednisolone (PRED_DOSE >= 20 mg/day) is",
    "associated with a 55% reduction in the absorption rate constant and a",
    "22% reduction in bioavailability (binary threshold-form covariate).",
    "Inter-occasion variability on bioavailability is encoded here as IIV",
    "on lfdepot because the source does not specify a per-subject occasion",
    "count for downstream simulation (see vignette Assumptions and",
    "deviations)."
  )
  reference <- paste(
    "Press RR, Ploeger BA, den Hartigh J, van der Straaten T, van Pelt H,",
    "Danhof M, de Fijter H, Guchelaar HJ.",
    "Explaining variability in ciclosporin exposure in adult kidney",
    "transplant recipients.",
    "Eur J Clin Pharmacol. 2010;66(6):579-590.",
    "doi:10.1007/s00228-010-0810-9.",
    sep = " "
  )
  vignette <- "Press_2010_ciclosporin"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight (baseline or time-varying)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling: CL = exp(lcl) * (WT / 76)^0.75 and",
        "Vc = exp(lvc) * (WT / 76)^1, with median body weight 76 kg as the",
        "reference and theory-based exponents 0.75 (clearance) and 1.0",
        "(volume) (Press 2010 Results 'CsA pharmacokinetic model'). Cohort",
        "range 49 - 140 kg, median 76 kg. The paper does not weight-scale",
        "the peripheral compartment volume Vp or the intercompartmental",
        "clearance Q. Body weight in the source data was treated as a",
        "subject baseline characteristic at the single visit timepoint for",
        "this covariate model; rxode2 simulations supplied at any time",
        "should populate this column at every observation row."
      ),
      source_name        = "WT"
    ),
    PRED_DOSE = list(
      description        = paste(
        "Concomitant oral prednisolone daily dose (mg/day). Used as a",
        "threshold-form covariate at >= 20 mg/day (binary high-dose",
        "indicator) on the absorption rate constant ka and on oral",
        "bioavailability F."
      ),
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Press 2010 Methods 'Covariate analysis' tested the daily",
        "prednisolone dose against the base model with body weight already",
        "incorporated. Final-model effects (Table 4 / Results 'CsA",
        "pharmacokinetic model'): when PRED_DOSE >= 20 mg/day the",
        "absorption rate constant ka is reduced by 55% (delta-OFV = +233",
        "on deletion) and oral bioavailability F is reduced by 22%",
        "(delta-OFV = +51 on deletion). The binary indicator is derived",
        "inside model() from the continuous PRED_DOSE column.",
        "Prednisolone was tapered from 50 mg twice daily (day 0) to 10 mg",
        "once daily by day 22 in the source cohort; the >= 20 mg/day high-",
        "dose stratum therefore corresponds to the earliest weeks post-",
        "transplant."
      ),
      source_name        = "Prednisolone daily dose"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 33L,
    n_studies      = 1L,
    age_range      = "18 - 70 years",
    age_median     = "43.8 +/- 14.5 years (once-daily arm, n = 17); 48.9 +/- 10.5 years (twice-daily arm, n = 16)",
    weight_range   = "49 - 140 kg",
    weight_median  = "76 kg",
    sex_female_pct = 21.2,
    race_ethnicity = c(Caucasian = 78.8, Other = 21.2),
    disease_state  = paste(
      "De novo adult kidney transplant recipients followed for 1 year",
      "after transplantation. First-graft recipients only, from deceased",
      "or living (non-HLA-identical) donors."
    ),
    dose_range     = paste(
      "Ciclosporin A (Neoral) oral 8 mg/kg/day, randomized to once-daily",
      "or twice-daily regimens. Dose adjusted by therapeutic drug",
      "monitoring against AUC0-24h targets of 10,800 ug*h/L (once daily)",
      "or AUC0-12h target 5,400 ug*h/L (twice daily) in the first 6",
      "weeks, then 6,500 ug*h/L or 3,250 ug*h/L respectively thereafter."
    ),
    regions        = "The Netherlands (Leiden University Medical Center)",
    co_medication  = paste(
      "Quadruple immunosuppression: basiliximab induction (day 0 + 4),",
      "fixed-dose mycophenolate mofetil (1,000 mg twice daily), tapering",
      "prednisolone (50 mg twice daily on day 0, tapered to 10 mg once",
      "daily by day 22), and ciclosporin A."
    ),
    notes          = paste(
      "Rich PK sampling: dense sampling (t = 0, 1, 2, 3, 4, 6, 24 h) on",
      "weeks 2, 6, 12, 26, 52 plus limited TDM sampling (t = 0, 2, 3 h) on",
      "weeks 4, 8, 10, 17, 21, 39. Most patients (22 of 33) provided data",
      "on eleven sampling occasions across the first year post-transplant",
      "(6 on ten occasions, 2 on twelve, 3 on three-to-six).",
      "Ciclosporin A concentrations were measured in whole blood by FPIA",
      "(Abbott TDx); assay linear up to 800 ug/L with inter-day variation",
      "10.4% / 7.8% / 7.5% at 70 / 300 / 600 ug/L respectively.",
      "Model fitted in NONMEM v VI release 1.2 (Icon Development",
      "Solutions) using first-order conditional estimation with",
      "interaction (FOCE-I) and a delta-OFV threshold of 6.63 (chi-square,",
      "1 df, p = 0.01) for covariate inclusion / deletion.",
      "Final-model precision assessed by a 500-replicate bootstrap;",
      "Table 4 reports both the mean and the 2.5 - 97.5 percentile",
      "bootstrap intervals."
    )
  )

  ini({
    # ============================================================
    # Structural PK parameters -- Press 2010 Table 4 'Final model' /
    # bootstrap mean column. Apparent CL/F, Vc/F, Vp/F, Q/F values are
    # reported in the source; here CL, Vc, Vp, Q denote the apparent
    # (oral) parameters that combine with the FIXED bioavailability
    # F = 0.5 to yield the modelled concentrations.
    # ============================================================
    lka <- log(2.0)
    label("Absorption rate constant ka (1/h)")
    # Table 4: ka = 2.0 1/h (95% CI 1.6 - 2.5)

    lcl <- log(15)
    label("Apparent clearance CL/F at WT = 76 kg (L/h)")
    # Table 4: CL = 15 L/h (95% CI 14 - 16); reported on the bootstrap
    # mean line; allometric scaling defined in Results 'CsA
    # pharmacokinetic model' as CL = 15 * (WT/76)^0.75.

    lvc <- log(56)
    label("Apparent central volume of distribution Vc/F at WT = 76 kg (L)")
    # Table 4: Vc = 56 L (95% CI 49 - 64); allometric scaling on Vc with
    # exponent 1 stated in Results 'CsA pharmacokinetic model'.

    lvp <- log(125)
    label("Apparent peripheral volume of distribution Vp/F (L)")
    # Table 4: Vp = 125 L (95% CI 100 - 149); not weight-scaled in the
    # source paper.

    lq <- log(14)
    label("Apparent inter-compartmental clearance Q/F (L/h)")
    # Table 4: Q = 14 L/h (95% CI 12 - 16); not weight-scaled in the
    # source paper.

    lfdepot <- fixed(log(0.5))
    label("Oral bioavailability F (unitless; FIXED)")
    # Methods 'Structural model': "The value for the oral bioavailability
    # was fixed at 50%, as previously described [refs 3, 23] and used in
    # the clinically applied TDM model [ref 4]."

    # ============================================================
    # Allometric exponents -- theory-based and held constant.
    # Results 'CsA pharmacokinetic model': "described allometrically ...
    # CL = 15 * (body weight/76)^0.75 ... typically with a value of 0.75
    # for clearance and 1 for volume of distribution [ref 24]."
    # ============================================================
    allo_cl <- fixed(0.75)
    label("Allometric exponent on CL (unitless; FIXED)")
    allo_vc <- fixed(1.0)
    label("Allometric exponent on Vc (unitless; FIXED)")

    # ============================================================
    # Covariate effects -- concomitant high-dose prednisolone
    # (PRED_DOSE >= 20 mg/day), Press 2010 Table 4. Encoded as
    # multiplicative fractional changes: (1 + e * pred_high) with the
    # tabulated fractional effects.
    # ============================================================
    e_pred_dose_high_ka <- -0.55
    label("Fractional change in ka for PRED_DOSE >= 20 mg/day (unitless)")
    # Table 4: "-55%" for the ka with concomitant DDPR >= 20 mg
    # (delta-OFV = +233 on deletion).

    e_pred_dose_high_f <- -0.22
    label("Fractional change in F for PRED_DOSE >= 20 mg/day (unitless)")
    # Table 4: "-22%" for the bioavailability F with concomitant DDPR
    # >= 20 mg (delta-OFV = +51 on deletion).

    # ============================================================
    # Inter-individual variability -- Press 2010 Table 4. Reported as
    # omega^2 (variance of the eta on the log-parameter scale) in the
    # 'Mean value' column. Quoted "Variability (%)" in the table is
    # sqrt(exp(omega^2) - 1) per the log-normal definition.
    # ============================================================
    etalka ~ 0.09
    # Table 4 IIV absorption rate omega^2 = 0.09 (Variability 30%,
    # CV 31%); 95% bootstrap CI 0.04 - 0.16 on omega^2.

    etalcl ~ 0.03
    # Table 4 IIV clearance omega^2 = 0.03 (Variability 17%, CV 24%);
    # 95% bootstrap CI 0.02 - 0.05.

    etalvc ~ 0.12
    # Table 4 IIV central volume omega^2 = 0.12 (Variability 35%,
    # CV 40%); 95% bootstrap CI 0.05 - 0.24.

    etalfdepot ~ 0.02
    # Table 4 IOV bioavailability omega^2 = 0.02 (Variability 14%,
    # CV 17%); 95% bootstrap CI 0.01 - 0.02. The source paper reports
    # this as inter-occasion variability across PK sampling visits;
    # encoded here as IIV on lfdepot because rxode2 simulations in the
    # downstream library API do not preserve a per-subject occasion
    # structure (see vignette Assumptions and deviations).

    # ============================================================
    # Residual error -- Press 2010 Methods 'Random effects':
    # "log(Cij) = log(Cpred_ij) + epsilon_ij ... epsilon_ij ... mean of
    # zero and variance sigma^2". Table 4 sigma^2 = 0.07 (Variability
    # 26%, CV 10%); 95% bootstrap CI 0.06 - 0.09. Additive error on the
    # log scale is equivalent to a proportional error on the linear
    # scale, encoded here as Cc ~ prop(propSd) with propSd = sqrt(0.07).
    # ============================================================
    propSd <- 0.2646
    label("Proportional residual error (fraction)")
    # sqrt(0.07) = 0.26458; reported as Table 4 sigma^2 = 0.07
    # (Variability 26%).
  })

  model({
    # ------------------------------------------------------------
    # 1. Derived covariate term -- binary high-dose prednisolone
    # indicator (Press 2010 Table 4: applies to ka and F when the daily
    # prednisolone dose is at or above 20 mg/day).
    # ------------------------------------------------------------
    pred_high <- (PRED_DOSE >= 20)

    # ------------------------------------------------------------
    # 2. Individual structural parameters with body-weight allometric
    # scaling on CL and Vc (Press 2010 Results 'CsA pharmacokinetic
    # model'). The peripheral volume Vp and intercompartmental clearance
    # Q are not weight-scaled in the source paper.
    # The high-dose prednisolone covariate enters multiplicatively on
    # ka and on f(depot).
    # ------------------------------------------------------------
    ka <- exp(lka + etalka) * (1 + e_pred_dose_high_ka * pred_high)
    cl <- exp(lcl + etalcl) * (WT / 76)^allo_cl
    vc <- exp(lvc + etalvc) * (WT / 76)^allo_vc
    vp <- exp(lvp)
    q  <- exp(lq)
    fdepot <- exp(lfdepot + etalfdepot) *
              (1 + e_pred_dose_high_f * pred_high)

    # ------------------------------------------------------------
    # 3. Micro-constants for the two-compartment disposition.
    # ------------------------------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ------------------------------------------------------------
    # 4. ODE system: dose -> depot -> transit1 -> central, with the
    # first-order transit rate constant set equal to ka (Press 2010
    # Results 'CsA pharmacokinetic model'). Two-compartment disposition
    # with first-order elimination from central.
    # ------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(transit1)    <-  ka * depot    - ka * transit1
    d/dt(central)     <-  ka * transit1 - kel * central -
                          k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # ------------------------------------------------------------
    # 5. Bioavailability on the dose compartment.
    # ------------------------------------------------------------
    f(depot) <- fdepot

    # ------------------------------------------------------------
    # 6. Observation: ciclosporin A whole-blood concentration. Dose
    # units mg; vc in L means central amount is in mg, so central / vc
    # is in mg/L. Multiply by 1000 to convert to ug/L (the paper's
    # reporting unit and the assay's calibration unit).
    # ------------------------------------------------------------
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
