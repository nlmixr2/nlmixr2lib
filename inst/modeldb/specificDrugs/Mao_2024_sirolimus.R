Mao_2024_sirolimus <- function() {
  description <- "One-compartment population PK model with first-order absorption and first-order elimination for sirolimus in adult liver transplant recipients, developed from routine therapeutic-drug-monitoring whole-blood trough concentrations (Mao 2024). The absorption rate constant was fixed at 0.75 1/h because only trough samples were collected. Hematocrit is the single retained covariate: it acts on apparent clearance through a power form normalized to the cohort median hematocrit of 38 percent with exponent -0.901, so higher hematocrit (more sirolimus sequestered in red blood cells) gives lower CL/F. Between-subject variability on CL/F and Vc/F is a correlated 2x2 block."
  reference <- paste(
    "Mao J, Cheng Y, Liu D, Zhang B, Li X. Dosing Regimen Recommendations",
    "for Sirolimus in Adult Liver Transplant Recipients: Insights from",
    "a Population Pharmacokinetic Model. Drug Des Devel Ther.",
    "2024;18:6379-6388. doi:10.2147/DDDT.S503463.",
    sep = " "
  )
  vignette <- "Mao_2024_sirolimus"
  units    <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Sirolimus C0 was measured in WHOLE BLOOD (Methods,
  # "Concentration Measurement": whole blood sirolimus concentrations by
  # enzyme multiplied immunoassay), so `central` is a whole-blood
  # compartment, not a plasma compartment.
  compartmentData <- list(
    depot   = list(analyte = "sirolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "sirolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    HCT = list(
      description        = "Hematocrit (packed red blood cell volume fraction).",
      units              = "%",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The only covariate retained in the final model (Mao 2024 Results:",
        "adding HCT dropped the objective function value by 32.9, and nothing",
        "was removed during backward elimination). Power form on CL/F",
        "normalized to the cohort median of 38 percent (Table 1 median",
        "hematocrit 38.0 percent; range 17.5-49.5 percent), per Eq. 8:",
        "CL/F = 7.09 * (HCT/38)^-0.901. The exponent is negative because",
        "sirolimus is roughly 95 percent bound to red blood cells, so a",
        "higher hematocrit lowers the unbound fraction and hence CL/F",
        "(Discussion). Values are percentages (e.g. 38), not fractions",
        "(0.38); using a fraction would rescale CL/F by 38^0.901.",
        "Recorded per subject at the time of therapeutic drug monitoring;",
        "the source does not state whether it was carried as baseline-only",
        "or time-varying, so treat it as a baseline covariate."
      ),
      source_name        = "HCT"
    )
  )

  # Covariates that Mao 2024 screened in the stepwise covariate search but
  # did NOT retain in the final model. Documentation only -- these are not
  # referenced in model(). Mycophenolic acid and Wuzhi capsule entered
  # during forward inclusion but were dropped in backward elimination
  # ("they failed to achieve statistical significance and poor estimation
  # of the influence coefficients", Discussion); no point estimate is
  # published for any of them, so none can be encoded.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age.", units = "years", type = "continuous",
      notes = "Screened as a continuous covariate (Methods, Covariate Model); not retained. Cohort median 51.0 years [28.0-74.0]."
    ),
    WT = list(
      description = "Total body weight.", units = "kg", type = "continuous",
      notes = "Screened as a continuous covariate; not retained. Cohort median 68.0 kg [39.0-90.0]. No allometric scaling appears in the final model."
    ),
    POD = list(
      description = "Postoperative day (days since liver transplantation).", units = "days", type = "continuous",
      notes = "Screened as a continuous covariate; not retained. Cohort median 92 days [14-699]."
    ),
    ALT = list(
      description = "Serum alanine aminotransferase activity.", units = "U/L", type = "continuous",
      notes = "Screened as a continuous covariate; not retained. Cohort median 21.0 U/L [5.0-187.0]."
    ),
    CRCL = list(
      description = "Creatinine clearance (Cockcroft-Gault).", units = "mL/min", type = "continuous",
      notes = "Screened as a continuous covariate; not retained. Cohort median 89.4 mL/min [27.0-170.0]. Table 1 footnote a gives the formula used: CRCL = [140 - age] * weight / [0.818 * SCR (umol/L)] * k, with k = 1 for male and 0.85 for female."
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male).", units = "unitless", type = "categorical",
      notes = "Screened as a categorical covariate; not retained. Cohort was 99 male / 4 female, so the female stratum carried almost no information."
    ),
    DOSE_SIROLIMUS_MGD = list(
      description = "Sirolimus daily dose.", units = "mg/day", type = "continuous",
      notes = "Screened as a continuous covariate to test for dose-dependent (nonlinear) PK; not retained. Cohort median 1.0 mg/day [0.5-2.0]."
    ),
    CONMED_MPA = list(
      description = "Concomitant mycophenolic acid.", units = "unitless", type = "categorical",
      notes = "Entered during forward inclusion but removed in backward elimination (Discussion); 24 of 103 patients. No point estimate published."
    ),
    CONMED_WUZHI = list(
      description = "Concomitant Wuzhi capsule (a Schisandra-derived traditional Chinese herbal preparation and potent CYP3A4 inhibitor).", units = "unitless", type = "categorical",
      notes = "Entered during forward inclusion but removed in backward elimination (Discussion); 26 of 103 patients. No point estimate published."
    ),
    CONMED_GCV = list(
      description = "Concomitant ganciclovir.", units = "unitless", type = "categorical",
      notes = "Screened as a categorical covariate; not retained. 5 of 103 patients."
    ),
    CONMED_STEROID = list(
      description = "Concomitant prednisone acetate.", units = "unitless", type = "categorical",
      notes = "Screened as a categorical covariate; not retained. 10 of 103 patients."
    ),
    CONMED_PPI = list(
      description = "Concomitant esomeprazole (proton-pump inhibitor).", units = "unitless", type = "categorical",
      notes = "Screened as a categorical covariate; not retained. 10 of 103 patients."
    ),
    CONMED_CCB = list(
      description = "Concomitant amlodipine (calcium-channel blocker).", units = "unitless", type = "categorical",
      notes = "Screened as a categorical covariate; not retained. 10 of 103 patients."
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 103L,
    n_studies       = 1L,
    n_observations  = 216L,
    age_range       = "28.0-74.0 years",
    age_median      = "51.0 years",
    weight_range    = "39.0-90.0 kg",
    weight_median   = "68.0 kg",
    sex_female_pct  = 3.9,
    race_ethnicity  = "Not reported; single-centre Chinese cohort (Wuhan), so predominantly Han Chinese.",
    disease_state   = "Adult liver transplant recipients receiving oral sirolimus tablets to prevent graft rejection.",
    dose_range      = "0.5-2.0 mg/day. Regimens present in the analysis dataset: 0.5 mg qd, 1 mg qd, 2 mg qd, 1 mg qod, 2 mg qod, 0.5 mg qd alternating with 1 mg qd, and 1 mg qd alternating with 2 mg qd.",
    regions         = "China (single centre: Tongji Hospital, Tongji Medical College, Huazhong University of Science and Technology, Wuhan).",
    hematocrit      = "Median 38.0 percent [17.5-49.5]; the 10th, 50th, and 90th percentiles used for the dosing simulations were 28, 38, and 46 percent.",
    renal_function  = "Creatinine clearance median 89.4 mL/min [27.0-170.0] (Cockcroft-Gault).",
    hepatic_function = "Post-transplant graft function still recovering; ALT median 21.0 U/L [5.0-187.0], AST median 25.0 U/L [12.0-190.0], albumin median 44.2 g/L [27.2-50.9].",
    co_medication   = "Ganciclovir 5/103, prednisone acetate 10/103, esomeprazole 10/103, amlodipine 10/103, mycophenolic acid 24/103, Wuzhi capsule 26/103.",
    assay           = "Enzyme multiplied immunoassay technique on an Architect i1000 analyser with the ARCHITECT Sirolimus Reagent Kit; assay range 2-30 ng/mL, detection limit 0.3 ng/mL, intra- and inter-day CV below 10 percent.",
    therapeutic_window = "4-8 ng/mL whole-blood trough (C0), the target used for the dosing-regimen simulations.",
    notes           = "Retrospective single-centre therapeutic-drug-monitoring cohort collected between January 2018 and August 2024; demographics from Table 1. ALL samples are pre-dose trough (C0) concentrations -- there is no information about the absorption or distribution phase in the data, which is why ka was fixed and why the authors flag Vc/F as less reliably estimated (Discussion, Limitations)."
  )

  ini({
    # Final-model estimates from Mao 2024 Table 2 ("Final Model Estimates"
    # column) and the printed final-model equations Eqs. 7-9 (Results,
    # Model Development). Table 2 and Eqs. 7-9 agree on every value.

    # ka was NOT estimated: "As only C0 data were collected, the value of
    # absorption rate constant (Ka) was fixed at 0.75 h-1, as reported in
    # previous studies" (Methods, Base Model). Table 2 lists the RSE for
    # ka as the literal word "Fixed" and gives no bootstrap interval.
    lka <- fixed(log(0.75)); label("First-order absorption rate constant, ka (1/h)")  # Eq. 7 / Table 2: ka = 0.75 1/h, FIXED

    # CL/F at the reference hematocrit of 38 percent (the cohort median).
    lcl <- log(7.09); label("Apparent clearance, CL/F, at HCT = 38 percent (L/h)")     # Eq. 8 / Table 2: CL/F = 7.09 L/h (RSE 4 percent)
    lvc <- log(496);  label("Apparent central volume of distribution, Vc/F (L)")       # Eq. 9 / Table 2: V/F = 496 L (RSE 15 percent)

    # Hematocrit effect on CL/F. Continuous covariates were entered with
    # the median-normalized power form of Eq. 5,
    #   P_i = TV(P) * (COV_i / COV_median)^theta,
    # instantiated in Eq. 8 as CL/F = 7.09 * (HCT/38)^-0.901. HCT is in
    # percent and 38 is the Table 1 median hematocrit.
    e_hct_cl <- -0.901; label("Power exponent for HCT (percent) on CL/F, normalized to 38 percent (unitless)")  # Eq. 8 / Table 2: HCT on CL/F = -0.901 (RSE 17 percent)

    # Between-subject variability. Eq. 1 is the exponential BSV model
    #   P_i = TV(P) * exp(eta_i),  eta ~ N(0, omega^2),
    # so each nlmixr2 variance is omega^2. Table 2 reports omega itself as
    # a percentage (rows "omega CL/F (%)" = 32.4 and "omega V/F (%)" =
    # 42.7), i.e. the standard deviation of eta on the natural-log scale,
    # so the variances are 0.324^2 and 0.427^2. The off-diagonal is
    # reported directly on the covariance scale (row "omega cov CL/F V/F"
    # = 0.0665, with no percent sign), which is consistent with this
    # reading: the implied correlation is
    #   0.0665 / (0.324 * 0.427) = 0.481,
    # a well-behaved value, and the resulting 2x2 block is positive
    # definite (determinant 0.0147).
    # A "#" comment placed INSIDE this c(...) breaks readModelDb(), so the
    # per-element provenance is recorded here instead:
    #   0.104976 = 0.324^2   Table 2, omega CL/F = 32.4 percent
    #   0.0665               Table 2, omega cov CL/F,V/F (covariance scale)
    #   0.182329 = 0.427^2   Table 2, omega V/F = 42.7 percent
    etalcl + etalvc ~ c(
      0.104976,
      0.0665, 0.182329
    )

    # Residual unexplained variability. Table 2 reports a single residual
    # row, "epsilon_por (%)" = 3.3 (RSE 10 percent, bootstrap 2.41-3.81),
    # and the Table 2 abbreviation list defines epsilon_por as
    # "proportional residual variability" -- so the proportional term is
    # 3.3 percent, i.e. 0.033 on the fraction scale used by prop().
    #
    # The Results text says "A combined error model was used to
    # characterize the RUV" (Eq. 4: Y = F * (1 + eps1) + eps2), but Table 2
    # publishes NO additive term. The combined structure is therefore
    # retained here with the additive standard deviation fixed at zero
    # rather than invented; see the vignette "Assumptions and deviations"
    # section. With addSd = 0 this reduces exactly to Eq. 3 (proportional).
    #
    # SCALE CAVEAT (unresolved in the source; see the vignette). The value
    # is read as a standard deviation because Table 2 tags the row "(%)",
    # and a variance cannot carry a percent unit. The authors appear to
    # have used that tag deliberately: the covariance row in the same block
    # ("omega cov CL/F V/F" = 0.0665) carries NO percent sign, exactly as a
    # covariance should not. Two independent signals nevertheless point the
    # other way -- that the printed 3.3 could be a VARIANCE (sigma^2 =
    # 0.033, i.e. sigma = 18.2 percent):
    #   (a) its RSE of 10 percent matches the asymptotic variance bound
    #       sqrt(2/n_obs) = sqrt(2/216) = 9.6 percent, not the SD bound
    #       sqrt(1/(2*n_obs)) = 4.8 percent;
    #   (b) the scatter of observations about the INDIVIDUAL predictions in
    #       Figure 1A is visibly wider than a 3.3 percent SD would allow --
    #       and eta shrinkage biases that scatter DOWNWARD, so it is a
    #       lower bound on sigma.
    # The published label is followed here because it is the paper's own
    # explicit statement, but the choice is flagged rather than buried.
    # Practical impact is modest: between-subject variability dominates the
    # steady-state trough, so the total trough CV moves only from about
    # 30 percent to about 35 percent between the two readings.
    propSd <- 0.033;   label("Proportional residual error (fraction)")       # Table 2: epsilon_por = 3.3 percent
    addSd  <- fixed(0); label("Additive residual error (ng/mL)")             # Results text says combined (Eq. 4); Table 2 publishes no additive value
  })

  model({
    # Reference hematocrit: the cohort median, Table 1 (38.0 percent) and
    # the denominator printed in Eq. 8.
    ref_hct <- 38

    # Individual parameters. ka carries no eta -- Table 2 reports no
    # between-subject variability for ka, consistent with it being fixed.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (HCT / ref_hct)^e_hct_cl   # Eq. 8
    vc <- exp(lvc + etalvc)                              # Eq. 9

    kel <- cl / vc

    # One-compartment disposition with first-order oral absorption.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Observation: whole-blood sirolimus concentration in ng/mL. Doses are
    # in mg and vc is in L, so central / vc is mg/L; multiply by 1000 to
    # convert to ng/mL (the units of every concentration in the paper).
    Cc <- central / vc * 1000
    Cc ~ add(addSd) + prop(propSd)
  })
}
