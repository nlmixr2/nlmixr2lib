Sood_2026_lu177dotatate_adult <- function() {
  description <- paste(
    "Two-compartment population PK model for the somatostatin-receptor-",
    "targeted radiopharmaceutical [177Lu]Lu-DOTATATE (177Lu-DOTATATE,",
    "lutetium Lu 177 dotatate) in ADULTS with advanced midgut",
    "neuroendocrine tumors, fitted to the 20 NETTER-1 patients who had",
    "radioactivity-blood pharmacokinetic sampling. Dosing is intravenous",
    "with a zero-order input (the clinical infusion), so no absorption",
    "compartment is present; the infusion duration is supplied by the",
    "event table. Radioactivity in blood was converted to a PEPTIDE MASS",
    "concentration before fitting, so the disposition parameters are mass",
    "clearances and volumes (CL 4.95 L/h, Vc 21.59 L, Q 4.78 L/h,",
    "Vp 202.15 L) and the observation Cc is a mass concentration in ng/mL,",
    "NOT a radioactivity concentration. Interindividual variability is",
    "carried on CL and Vc only, as a correlated pair (r = 0.70); all other",
    "disposition parameters are population-only fixed effects. NO",
    "COVARIATE IS RETAINED and that absence is the paper's central result,",
    "not a transcription gap: a stepwise covariate search over age, body",
    "weight, body surface area, creatinine clearance and kidney mass",
    "returned nothing significant, which is what licensed flat 7.4 GBq",
    "dosing to be carried from adults to adolescents. Companion models",
    "from the same paper: Sood_2026_lu177dotatate_adolescent.R (the",
    "NETTER-P adolescent refit), Sood_2026_lu177dotatate_dosimetry_adult.R",
    "and Sood_2026_lu177dotatate_dosimetry_pooled.R (the empirical",
    "kidney / bone-marrow absorbed-dose regressions).",
    sep = " "
  )
  reference <- paste(
    "Sood M, Lachi Silva L, Ho YY, Blumenstein L, Cherfi A, Xu L,",
    "Khanshan F. [177Lu]Lu-DOTATATE population pharmacokinetics and",
    "dosimetry modeling for adolescent and adult patients with",
    "somatostatin receptor-positive gastroenteropancreatic neuroendocrine",
    "tumors. J Nucl Med. 2026;67(6):887-894.",
    "doi:10.2967/jnumed.125.270202",
    sep = " "
  )
  vignette <- "Sood_2026_lu177dotatate"
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL")

  # No covariate was retained in the final model. Everything the paper
  # screened is documented below rather than in covariateData, because a
  # covariateData entry that model() never references is flagged as unused.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline. Table 1 reports the adult DOSIMETRY set (n = 47) as median 56 years, range 29-83; the paper does not tabulate the 20-patient popPK subset separately.",
      units = "years",
      type = "continuous",
      notes = paste(
        "Screened by the MonolixSuite stepwise covariate model (SCM) over",
        "'baseline characteristics (e.g., age, weight, CrCL, and kidney",
        "size) on 177Lu-DOTATATE plasma pharmacokinetics'. Not retained:",
        "'No covariates were statistically significant in the final model'",
        "(supplemental materials, 'Population pharmacokinetics for",
        "adults'). No coefficient is reported, so no effect can be",
        "encoded.",
        sep = " "
      ),
      source_name = "Age (y)"
    ),
    WT = list(
      description = "Total body weight at baseline. Table 1, adult dosimetry set (n = 47): median 75 kg, range 48-145.",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Screened by the SCM and not retained. The Results state it",
        "explicitly: 'None of the baseline characteristics (e.g., age,",
        "weight, body surface area) were significant covariates on",
        "pharmacokinetics parameters (clearance or Vc).' This is the",
        "finding that licensed flat (non-weight-adjusted) dosing:",
        "'Weight-adjusted dosing should not apply to 177Lu-DOTATATE in",
        "adolescents, since no effects of weight, age, or body surface",
        "area were observed on pharmacokinetics or biodistribution in",
        "adults.'",
        sep = " "
      ),
      source_name = "Body weight (kg)"
    ),
    BSA = list(
      description = "Body surface area at baseline. Table 1, adult dosimetry set (n = 47): median 1.9 m^2, range 1.48-2.68.",
      units = "m^2",
      type = "continuous",
      notes = "Screened by the SCM and not retained (Results, 'None of the baseline characteristics ... were significant covariates'). No coefficient is reported.",
      source_name = "BSA (m2)"
    ),
    CRCL = list(
      description = paste(
        "Creatinine clearance, raw (NOT body-surface-area normalized)",
        "mL/min. Table 1, adult dosimetry set (n = 47): median",
        "98.82 mL/min, range 46.97-189.77. The paper does not name the",
        "estimating equation for the adult cohort; the adolescent",
        "companion cohort is anchored to Piepsz 2008, whose whole subject",
        "is escaping the body-surface-area correction, so the column is",
        "absolute mL/min throughout.",
        sep = " "
      ),
      units = "mL/min",
      type = "continuous",
      notes = paste(
        "Screened by the SCM on the plasma PK parameters and not",
        "retained. Creatinine clearance IS retained, with a large effect,",
        "in the companion exposure-dosimetry models",
        "(Sood_2026_lu177dotatate_dosimetry_adult.R and",
        "Sood_2026_lu177dotatate_dosimetry_pooled.R) -- renal function",
        "drives the kidney and bone-marrow absorbed dose without",
        "detectably driving plasma disposition.",
        sep = " "
      ),
      source_name = "CrCL (mL/min)"
    )
  )

  compartmentData <- list(
    central = list(analyte = "[177Lu]Lu-DOTATATE peptide (mass)", units = "ug", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "[177Lu]Lu-DOTATATE peptide (mass)", units = "ug", specimen = "whole blood", verified = FALSE)
  )

  population <- list(
    species = "human",
    n_subjects = 20,
    n_studies = 1,
    age_median = "56 years (adult dosimetry set, n = 47)",
    age_range = "29-83 years (adult dosimetry set, n = 47)",
    weight_median = "75 kg (adult dosimetry set, n = 47)",
    weight_range = "48-145 kg (adult dosimetry set, n = 47)",
    sex_female_pct = 48.9,
    disease_state = "Advanced, progressive, well-differentiated grade 1 or 2 somatostatin-receptor-positive midgut neuroendocrine tumors.",
    renal_function = "Creatinine clearance median 98.82 mL/min, range 46.97-189.77 (Table 1, adult dosimetry set n = 47); kidney mass median 339 g, range 201-575.",
    dose_range = "Intravenous [177Lu]Lu-DOTATATE 7.4 GBq per cycle, 4 cycles 8 weeks apart (cumulative 29.6 GBq). The corresponding PEPTIDE MASS dose -- the quantity this model is actually dosed with -- is described only as the 'actual mass dose' and is never given a number anywhere in the paper or supplement.",
    regions = "NETTER-1 (NCT01578239), multinational phase 3",
    notes = paste(
      "The popPK analysis set is the 20 NETTER-1 patients with",
      "radioactivity-blood sampling; Table 1 tabulates the larger n = 47",
      "adult DOSIMETRY set (20 NETTER-1 plus 27 ERASMUS), so every",
      "demographic above describes that superset and not the 20-patient",
      "popPK subset, which the paper never tabulates separately. Blood",
      "radioactivity was converted to peptide mass before fitting.",
      "Estimation used MonolixSuite 2021R2 with an empirical-Bayes",
      "individual-parameter step; exposure metrics were generated with",
      "Simulx 2021R2. Kidney mass was among the screened covariates",
      "(supplement: 'age, weight, CrCL, and kidney size') and is not",
      "carried as a canonical column here because it was neither",
      "retained nor is it registered.",
      sep = " "
    )
  )

  ini({
    # ==================================================================
    # Structural disposition -- Table 2 ('Adult PopPK Parameter
    # Estimates', NETTER-1 estimate column; parentheses are %RSE). The
    # structural model is supplemental Eqs. S1 and S2, a two-compartment
    # system parameterised in clearance and volume with zero-order input
    # and first-order elimination.
    # ==================================================================
    lcl <- log(4.95)    ; label("Clearance from the central compartment, CL (L/h)")                   # Table 2 CL = 4.95 L/h (%RSE 9.81); Results repeat 'The population clearance ... estimated to be 4.95 L/h'
    lvc <- log(21.59)   ; label("Central volume of distribution, Vc (L)")                             # Table 2 Vc = 21.59 L (%RSE 12.70); Results repeat 'volume of distribution of the central compartment (Vc) ... 21.59 L'
    lq  <- log(4.78)    ; label("Intercompartmental clearance, Q (L/h)")                              # Table 2 Q = 4.78 L/h (%RSE 6.69)
    lvp <- log(202.15)  ; label("Peripheral volume of distribution, Vp (L)")                          # Table 2 Vp = 202.15 L (%RSE 10.12)

    # ==================================================================
    # Interindividual variability. MonolixSuite reports omega_<param> as
    # the STANDARD DEVIATION of the log-scale random effect and
    # corr_<p1>_<p2> as a correlation, so the variances below are the
    # squares of the printed values and the covariance is
    # r * sd_CL * sd_Vc. The supplement fixes which parameters carry a
    # random effect: 'The final popPK model allowed individual random
    # effects on CL and Vc and other PK parameters were fixed effects'.
    # ==================================================================
    # Block entries below, in nlmixr2's lower-triangular order
    # c(var(etalcl), cov, var(etalvc)):
    #   var(etalcl) = 0.41^2 = 0.1681   from Table 2 'IIV on CL' = 0.41 (RSE 18.29)
    #   cov         = 0.70 * 0.41 * 0.52 = 0.149240
    #                                   from Table 2 'Cor CL ~ Vc' = 0.70 (RSE 20.79)
    #   var(etalvc) = 0.52^2 = 0.2704   from Table 2 'IIV on Vc' = 0.52 (RSE 19.37)
    # No trailing comment sits on the eta line itself: rxode2 promotes a
    # trailing comment on an unlabelled ini() line into a label(), and a
    # comment inside the c() breaks the re-parse outright.
    etalcl + etalvc ~ c(0.1681, 0.149240, 0.2704)

    # ==================================================================
    # Residual unexplained variability. Both the main text ('Constant
    # residual error ... 0.40') and the supplement ('additive (constant)
    # error model was assumed') name a CONSTANT error model, and that
    # constant is applied on the LOG-TRANSFORMED observation, i.e. it is
    # an exponential (log-normal) residual with SD 0.40 rather than an
    # additive 0.40 ng/mL. Figure 1A settles it: the prediction-corrected
    # VPC spans 0.01-100 ng/mL with a roughly constant relative spread,
    # and its 90% prediction intervals track the decline to 0.05 ng/mL.
    # An additive SD of 0.40 ng/mL would be eight-fold larger than the
    # median concentration over the whole terminal phase and would drive
    # most simulated observations negative, which no log-axis VPC of that
    # shape can be produced from. See the vignette 'Assumptions and
    # deviations' section, which reproduces the falsification.
    # ==================================================================
    expSd <- 0.40       ; label("Exponential (log-scale constant) residual SD for the plasma concentration (fraction)")  # Table 2 'Constant residual error' = 0.40 (%RSE 4.85)
  })

  model({
    # Individual disposition parameters. Only CL and Vc carry a random
    # effect (supplement, 'Population pharmacokinetics for adults').
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    # Micro-constants for supplemental Eqs. S1 and S2, which are written
    # in amounts:
    #   dAc/dt = (Q/Vp)*Ap - (Q/Vc)*Ac - (CL/Vc)*Ac
    #   dAp/dt = (Q/Vc)*Ac - (Q/Vp)*Ap
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Two-compartment disposition. The zero-order input of the paper's
    # structural model is the clinical intravenous infusion and is
    # supplied by the event table (rate / duration on the dose record),
    # so there is no absorption state.
    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Peptide MASS concentration in blood; radioactivity measurements were
    # converted to mass before fitting (Methods, 'PopPK for Adults'). No unit
    # scaling is needed despite units$dosing being ug and
    # units$concentration ng/mL: the volumes are in L, and ug / L IS ng/mL.
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
