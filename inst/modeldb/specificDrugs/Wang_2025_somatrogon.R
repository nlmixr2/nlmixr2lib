Wang_2025_somatrogon <- function() {
  description <- "Two-compartment population PK model with delayed first-order absorption for somatrogon in children with growth hormone deficiency, estimated by fully Bayesian MCMC (Stan/Torsten) (Wang 2025)"
  reference <- "Wang Y, Pei X, Niu T, Korth-Bradley J, Fostvedt L. Implementing a Bayesian approach using Stan with Torsten: Population pharmacokinetics analysis of somatrogon. CPT Pharmacometrics Syst Pharmacol. 2025;14(2):351-364. doi:10.1002/psp4.13279"
  vignette <- "Wang_2025_somatrogon"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Wang 2025 Data S5 (Stan model file):
  # `pmx_solve_twocpt` with a first-order absorption compartment, and
  # `cHat[i] = x[2,i]/theta[i,3]` (central amount / Vc), i.e. compartment 1 =
  # depot (subcutaneous injection site), 2 = central (plasma), 3 = peripheral.
  compartmentData <- list(
    depot       = list(analyte = "somatrogon", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "somatrogon", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "somatrogon", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling on CL/F, Q/F, Vc/F and Vp/F with reference weight 15 kg.",
        "The reference value is not printed in the article; it is read from the",
        "NONMEM control stream shipped as Wang 2025 Data S3, which defines",
        "CLBWT = ((WT/15)**THETA(9)) with the inline comment '15 is the median body",
        "weight for study 004'. Table 3 of the article reports the Study 004 median",
        "body weight as 14.8 kg, i.e. 15 kg is that median rounded. In the Stan model",
        "(Data S5) the same quantity enters pre-computed as the data column WTrt,",
        "so the article and Stan code alone do not disclose the normalizer.",
        "Both exponents are estimated, not fixed at the 0.75/1 allometric defaults."
      ),
      source_name        = "WT"
    ),
    ADA_POS = list(
      description        = "Anti-drug (anti-somatrogon) antibody positive status at the time of the observation",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ADA-negative)",
      notes              = paste(
        "Time-varying: the article's ADAT column is the ADA status at each",
        "observation, dichotomized positive (1) / negative (0) (Methods, 'Modeling",
        "approach'). ADA-positive occasions carry both a proportional shift in CL/F",
        "(e_ada_cl) and an additional gated between-subject random effect on CL/F",
        "(etalcl_ada). In Data S5 both factors sit inside a single",
        "`if (ADAT[i] == 1 && ADAS[i] == 1)` branch, where ADAS is the subject-level",
        "'ever ADA-positive' flag; because ADAT can only be 1 for a subject who is",
        "ever positive, the compound condition reduces to ADAT == 1 for any",
        "internally consistent data set, and this model encodes it as ADA_POS alone.",
        "42.9% of Study 004 participants were ADA-positive overall and 23.8% within",
        "the first year of dosing (article Table 3); only the first year of",
        "observations was used for estimation."
      ),
      source_name        = "ADAT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 42,
    n_studies      = 1,
    n_observations = 560,
    age_range      = "3-11 years",
    age_median     = "5.5 years",
    weight_range   = "10-26.3 kg",
    weight_median  = "14.8 kg",
    sex_female_pct = 33.3,
    race_ethnicity = c(White = 95.2, `African American` = 2.4, `Other or Missing` = 2.4),
    disease_state  = "Pediatric growth hormone deficiency (GHD)",
    dose_range     = "3.05-17.5 mg/week (0.228-0.711 mg/kg/week) by once-weekly subcutaneous injection; median 6.74 mg/week (0.482 mg/kg/week)",
    regions        = "Not reported in the article",
    notes          = paste(
      "Baseline demographics from Wang 2025 Table 3, column 'Phase II (004)'.",
      "Posterior sampling used the 560 observations from the 42 pediatric",
      "participants of the Phase II study CP-4-004 (Results, first paragraph).",
      "Only PK samples collected within 1 year of study start were used, to reduce",
      "MCMC run time (Methods, 'Data inclusion criteria'), even though Study 004",
      "followed participants for up to 5.5 years. The 109 pediatric participants of",
      "the Phase III study CP-4-006 were NOT used for estimation; they served only",
      "as an external set for posterior prediction (Table 2), so their demographics",
      "are not reflected in the parameters encoded here."
    )
  )

  ini({
    # =========================================================================
    # Source of the parameter values: Wang 2025 Table 5, column "Uniform prior
    # set" (posterior means from the semi-centered parameterization, which is
    # the parameterization Table 5 reports).
    #
    # The article reports three prior sets (uniform / weakly informative,
    # moderate-informative, very informative) for one and the same structural
    # model, and does not designate any of them as "the" final model. The
    # uniform set is encoded here because it is the only one whose posterior is
    # driven by the Study 004 data alone: the two informative sets center their
    # priors on estimates from the Phase II PopPK model, which itself was fit to
    # Study 004, and the authors explicitly flag this as non-ideal ("While this
    # is not an ideal source of 'external' information ... It is recommended to
    # avoid 'double-dipping' approaches when constructing prior distributions",
    # Discussion). The article's own conclusion is that the three prior sets are
    # predictively equivalent ("there is no noticeable difference between the
    # three prior sets in the 90% posterior predicted intervals nor the
    # posterior predicted median values", Results). All three parameter sets,
    # and the previously published NONMEM FOCEI estimates, are tabulated and
    # simulated side by side in the validation vignette.
    #
    # Structural parameters are apparent (/F) throughout: bioavailability was
    # fixed to 1 for the pediatric participants (Data S3 NONMEM stream,
    # F1 = THETA(7)*FPROT with THETA(7) = 1 FIX and FPROT = 1 for PROT = 4;
    # Data S5 Stan model, F[1] = 1), so no lfdepot parameter is estimable and
    # none is carried here.
    #
    # Reference subject for the structural values: 15 kg, ADA-negative.
    # =========================================================================
    lcl <- log(0.478);  label("Apparent clearance at 15 kg, ADA-negative (L/h)")          # Wang 2025 Table 5 'CL/F (L/h)', uniform prior set posterior mean 0.478 (90% CrI 0.416, 0.545)
    lq  <- log(0.065);  label("Apparent intercompartmental clearance at 15 kg (L/h)")      # Wang 2025 Table 5 'Q/F (L/h)', uniform prior set posterior mean 0.065 (90% CrI 0.039, 0.098)
    lvc <- log(6.805);  label("Apparent central volume of distribution at 15 kg (L)")     # Wang 2025 Table 5 'Vc/F (L)', uniform prior set posterior mean 6.805 (90% CrI 4.775, 9.503)
    lvp <- log(2.303);  label("Apparent peripheral volume of distribution at 15 kg (L)")  # Wang 2025 Table 5 'Vp/R (L)' [sic; Vp/F per the table abbreviations], uniform prior set posterior mean 2.303 (90% CrI 1.568, 3.125)
    lka <- log(0.178);  label("First-order absorption rate constant (1/h)")                 # Wang 2025 Table 5 'Ka (1/h)', uniform prior set posterior mean 0.178 (90% CrI 0.11, 0.309)

    ltlag <- log(1.116); label("Absorption lag time (h)")                                       # Wang 2025 Table 5 'Lag time (h)', uniform prior set posterior mean 1.116 (90% CrI 0.282, 1.625); applied to the depot only (Data S5: tlag[1] = lag0, tlag[2] = tlag[3] = 0)

    # Allometric exponents. Both are ESTIMATED (not fixed at 0.75/1), and each
    # is a single value shared by two parameters, exactly as the article's
    # Table 5 rows are labelled ("Weight effect on CL/F and Q/F", "Weight effect
    # on Vc/F and Vp/F") and as Data S5 applies them (thetaWT1 to theta[,1] and
    # theta[,2]; thetaWT2 to theta[,3] and theta[,4]).
    e_wt_cl_q  <- 1.258; label("Allometric (WT) exponent shared across CL/F and Q/F (unitless)")   # Wang 2025 Table 5 'Weight effect on CL/F and Q/F', uniform prior set posterior mean 1.258 (90% CrI 0.83, 1.715)
    e_wt_vc_vp <- 1.341; label("Allometric (WT) exponent shared across Vc/F and Vp/F (unitless)")  # Wang 2025 Table 5 'Weight effect on Vc/F and Vp/F', uniform prior set posterior mean 1.341 (90% CrI 0.722, 2.004)

    # Proportional change in CL/F on ADA-positive occasions. Parameterization is
    # P_i = P_pop * (1 + theta_ADAT * ADAT), the displayed (unnumbered) equation
    # in Methods 'Modeling approach'. The negative sign means ADA-positive
    # occasions have LOWER apparent clearance. The 90% CrI spans zero, which the
    # authors attribute to the small Study 004 sample; the effect is retained
    # here because it is a structural component of the fitted model.
    e_ada_cl <- -0.111; label("Proportional change in CL/F when ADA-positive (fraction)")  # Wang 2025 Table 5 'ADAT effect on CL/F', uniform prior set posterior mean -0.111 (90% CrI -0.267, 0.063)

    # =========================================================================
    # Inter-individual variability: a full 4x4 covariance matrix on CL/F, Vc/F,
    # Ka and the ADA-gated CL/F effect, all as multiplicative exponential random
    # effects (Methods, 'Modeling approach'). Unlike the previous NONMEM
    # analyses, which assumed independence, this analysis estimated every
    # off-diagonal element (Wishart prior, identical in all three prior sets).
    # Q/F and Vp/F carry no IIV.
    #
    # The values below are variances and covariances on the log scale, taken
    # directly from Wang 2025 Table 5 (uniform prior set); they are NOT CV%, so
    # no omega^2 = log(CV^2 + 1) conversion applies. Ordering of the block is
    # (CL/F, Vc/F, Ka, ADA-gated CL/F), matching the Omega index order in
    # Data S5 (logeta[j,1] on CL, [j,2] on Vc, [j,3] on Ka, [j,4] on the ADA
    # branch of CL). The matrix as printed is positive definite (eigenvalues
    # 0.386, 0.309, 0.0965, 0.0210), so it is used verbatim with no repair.
    #
    # etalcl_ada is the "additional between-patient variance on the apparent
    # clearance CL/F ... if the observation was with positive anti-drug antibody
    # status" (Methods). It is a SECOND random effect on CL/F that is switched
    # on only when ADA_POS = 1; see model(). Naming follows the gated-eta
    # precedent in Svensson_2012_nevirapine.R (etalfdepot_tb, gated by TB_POS).
    #
    # Element-by-element source trace, in the row order written below. Every
    # value is from Wang 2025 Table 5, uniform prior set:
    #   0.058  'Variance on CL/F'                     (90% CrI 0.021, 0.113)
    #   0.106  'Covariance on CL/F&Vc/F'              (90% CrI 0.02, 0.224)
    #   0.352  'Variance on Vc/F'                     (90% CrI 0.113, 0.685)
    #  -0.024  'Covariance on CL/F&Ka'                (90% CrI -0.115, 0.06)
    #   0.006  'Covariance on Vc/F&Ka'                (90% CrI -0.191, 0.235)
    #   0.307  'Variance on Ka'                       (90% CrI 0.069, 0.711)
    #   0.005  'Covariance on CL/F&ADAT'              (90% CrI -0.04, 0.049)
    #  -0.007  'Covariance on Vc/F&ADAT'              (90% CrI -0.136, 0.108)
    #   0.000  'Covariance on Ka&ADAT'                (90% CrI -0.121, 0.122)
    #   0.096  'Variance of ADAT effect on CL/F'      (90% CrI 0.017, 0.266)
    # The implied CL/F-Vc/F correlation is 0.742, the only off-diagonal whose
    # 90% CrI excludes zero (Results, final paragraph of the Table 5 discussion).
    #
    # NOTE: keep this block free of inline comments. Trailing comments inside an
    # ini() eta block are replaced by a bare ';' when the conventions linter
    # strips comments and re-parses, which turns c(0.058, # ...) into a syntax
    # error. The per-element trace therefore lives above, not on the lines.
    # =========================================================================
    etalcl + etalvc + etalka + etalcl_ada ~
      c(0.058,
        0.106,  0.352,
       -0.024,  0.006, 0.307,
        0.005, -0.007, 0.000, 0.096)

    # Residual error: additive on log-transformed concentrations,
    # log(Y_ij) = log(F_ij) + eps_ij with eps ~ N(0, sigma^2) (Methods,
    # displayed equation). On the linear scale that is a log-normal residual,
    # i.e. Cc ~ lnorm(sigma). The tabulated "Residual Deviance sigma" is the
    # standard deviation, not the variance: Data S5 writes
    # `logCObs ~ normal(log(cHatObs), sigma)` and Stan's normal() takes an SD as
    # its second argument; the corresponding NONMEM stream (Data S3) likewise
    # fixes $SIGMA to 1 and estimates the scale as W in Y = IPRED + W*EPS(1).
    # So no sqrt() is applied here.
    expSd <- 0.691; label("Log-scale (exponential) residual error SD")  # Wang 2025 Table 5 'Residual Deviance sigma', uniform prior set posterior mean 0.691 (90% CrI 0.655, 0.728)
  })

  model({
    # Individual parameters. Allometric scaling is on WT/15 (see covariateData).
    ka  <- exp(lka + etalka)
    vc  <- exp(lvc + etalvc) * (WT / 15)^e_wt_vc_vp
    vp  <- exp(lvp)          * (WT / 15)^e_wt_vc_vp
    q   <- exp(lq)           * (WT / 15)^e_wt_cl_q

    # CL/F carries three multiplicative terms beyond the reference value:
    # allometric weight scaling, the proportional ADA effect (1 + e_ada_cl) and
    # the ADA-gated extra random effect. The latter two are active only on
    # ADA-positive occasions, so with ADA_POS = 0 both collapse to 1 and CL/F
    # reduces to exp(lcl + etalcl) * (WT/15)^e_wt_cl_q -- exactly the two
    # branches of the `if (ADAT == 1 && ADAS == 1)` conditional in Data S5.
    cl  <- exp(lcl + etalcl) * (WT / 15)^e_wt_cl_q *
      (1 + e_ada_cl * ADA_POS) * exp(etalcl_ada * ADA_POS)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Delayed first-order absorption: the lag applies to the depot only.
    tlag <- exp(ltlag)
    alag(depot) <- tlag

    # Concentration: dose in mg, volume in L -> mg/L = ug/mL
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
