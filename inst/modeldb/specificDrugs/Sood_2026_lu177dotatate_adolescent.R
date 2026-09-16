Sood_2026_lu177dotatate_adolescent <- function() {
  description <- paste(
    "Two-compartment population PK model for the somatostatin-receptor-",
    "targeted radiopharmaceutical [177Lu]Lu-DOTATATE (177Lu-DOTATATE,",
    "lutetium Lu 177 dotatate) in ADOLESCENTS aged 12 to under 18 years",
    "with somatostatin-receptor-positive gastroenteropancreatic",
    "neuroendocrine tumors or pheochromocytomas and paragangliomas,",
    "fitted to the 11 patients of the phase 2 NETTER-P trial. The",
    "structural model is identical to the adult NETTER-1 fit",
    "(Sood_2026_lu177dotatate_adult.R): intravenous zero-order input (the",
    "clinical infusion, supplied by the event table), first-order",
    "elimination, no absorption compartment, and blood radioactivity",
    "converted to a PEPTIDE MASS concentration before fitting, so CL",
    "5.92 L/h, Vc 17.86 L, Q 2.63 L/h and Vp 99.42 L are mass clearances",
    "and volumes and the observation Cc is a mass concentration in ng/mL.",
    "Because of the small sample size only CLEARANCE carries",
    "interindividual variability; every other disposition parameter is a",
    "population-only fixed effect. NO COVARIATE IS RETAINED. That is a",
    "result and not a gap, and it survived two independent searches: a",
    "stepwise covariate model did flag kidney mass and body surface area",
    "on clearance, but at %RSE 251.2, and a horseshoe-prior re-analysis",
    "then returned 95% credible intervals containing zero for every",
    "tested covariate. Together with the adult model this is what",
    "supports flat 7.4 GBq dosing in adolescents. Companion models from",
    "the same paper: Sood_2026_lu177dotatate_adult.R,",
    "Sood_2026_lu177dotatate_dosimetry_adult.R and",
    "Sood_2026_lu177dotatate_dosimetry_pooled.R.",
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

  # No covariate was retained. Everything screened is documented here
  # rather than in covariateData, because a covariateData entry that
  # model() never references is flagged as unused.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline. Table 1, adolescent set (n = 11): median 15 years, range 13-17. Trial eligibility was 12 to under 18 years.",
      units = "years",
      type = "continuous",
      notes = "Screened and not retained; the horseshoe-prior posterior for every tested covariate effect on clearance had a 95% credible interval containing zero (Supplemental Fig. 5). No coefficient is reported.",
      source_name = "Age (y)"
    ),
    WT = list(
      description = "Total body weight at baseline. Table 1, adolescent set (n = 11): median 55 kg, range 39.5-71.",
      units = "kg",
      type = "continuous",
      notes = "Screened and not retained. The Discussion states the conclusion for both cohorts: 'Neither 177Lu-DOTATATE exposure nor organ dosimetry was impacted by age or body weight in a clinically relevant manner'.",
      source_name = "Body weight (kg)"
    ),
    BSA = list(
      description = "Body surface area at baseline. Table 1, adolescent set (n = 11): median 1.58 m^2, range 1.3-1.85.",
      units = "m^2",
      type = "continuous",
      notes = paste(
        "The one covariate pair that the stepwise covariate model did",
        "flag: 'in NETTER-P, the SCM approach indicated kidney mass and",
        "BSA effects on CL'. Including both dropped -2*log-likelihood",
        "from 17.53 to -0.3, but the effects were 'estimated with very low",
        "precision (%RSE 251.2)' and the horseshoe-prior re-analysis found",
        "no covariate whose 95% credible interval excluded zero, so the",
        "covariate-free model was selected as final. No coefficient is",
        "printed, so the flagged effect cannot be encoded even as an",
        "alternative.",
        sep = " "
      ),
      source_name = "BSA (m2)"
    ),
    CRCL = list(
      description = paste(
        "Creatinine clearance, raw (NOT body-surface-area normalized)",
        "mL/min. Table 1, adolescent set (n = 11): median 122.1 mL/min,",
        "range 86-160. NETTER-P required a creatinine clearance of at",
        "least 70 mL/min at entry, so the cohort carries no information",
        "about renal impairment. The absolute (un-normalized) scale is",
        "anchored to Piepsz 2008, the paediatric reference the paper uses",
        "for its simulated creatinine-clearance distribution, whose",
        "subject is escaping the body-surface-area correction.",
        sep = " "
      ),
      units = "mL/min",
      type = "continuous",
      notes = "Screened and not retained on the plasma PK parameters. Creatinine clearance IS retained in the companion exposure-dosimetry models, where it is the dominant covariate on organ absorbed dose.",
      source_name = "CrCL (mL/min)"
    )
  )

  compartmentData <- list(
    central = list(analyte = "[177Lu]Lu-DOTATATE peptide (mass)", units = "ug", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "[177Lu]Lu-DOTATATE peptide (mass)", units = "ug", specimen = "whole blood", verified = FALSE)
  )

  population <- list(
    species = "human",
    n_subjects = 11,
    n_studies = 1,
    age_median = "15 years",
    age_range = "13-17 years (eligibility 12 to under 18 years)",
    weight_median = "55 kg",
    weight_range = "39.5-71 kg",
    sex_female_pct = 54.5,
    disease_state = "Somatostatin-receptor-positive gastroenteropancreatic neuroendocrine tumors (4 of 11, 36.4%) or somatostatin-receptor-positive pheochromocytomas and paragangliomas (7 of 11, 63.6%), metastasized or locally advanced and inoperable, grade 1 or 2 with a Ki-67 index of 20% or less for the neuroendocrine tumors.",
    renal_function = "Creatinine clearance median 122.1 mL/min, range 86-160; entry required at least 70 mL/min. Kidney mass median 273.7 g, range 169.3-346.7.",
    dose_range = "Intravenous [177Lu]Lu-DOTATATE 7.4 GBq per cycle for 4 cycles, every 8 (plus or minus 1) weeks, co-infused with a 2.5% lysine-arginine amino acid solution for renal protection. The corresponding PEPTIDE MASS dose -- the quantity this model is dosed with -- is never given a number in the paper or supplement.",
    regions = "NETTER-P (NCT04711135), multicenter open-label single-arm phase 2",
    notes = paste(
      "Radioactivity in blood was measured with a gamma counter during",
      "treatment cycle 1 (cycle 2 for one patient) at pre-administration,",
      "end of infusion, 2 h, 6 h, 24 h and 72 h after infusion, and was",
      "converted to peptide mass before fitting. Estimation used",
      "MonolixSuite 2021R1. Karnofsky or Lansky performance score at",
      "least 50 was required at entry. Data cutoff 12 March 2024. Ten of",
      "these 11 patients also enter the pooled exposure-dosimetry model",
      "(Sood_2026_lu177dotatate_dosimetry_pooled.R); one patient with a",
      "pheochromocytoma or paraganglioma was excluded there because the",
      "planar and SPECT/CT images could not be used for dosimetry.",
      sep = " "
    )
  )

  ini({
    # ==================================================================
    # Structural disposition -- Table 4 ('Adolescent PopPK Parameter
    # Estimates', NETTER-P estimate column; parentheses are %RSE). The
    # structural model is supplemental Eqs. S1 and S2, the same
    # two-compartment clearance / volume system used for the adults.
    # ==================================================================
    lcl <- log(5.92)   ; label("Clearance from the central compartment, CL (L/h)")   # Table 4 CL = 5.92 L/h (%RSE 5.43); Results repeat 'The clearance for an adolescent patient was estimated at 5.92 L/h (%RSE, 5.43)'
    lvc <- log(17.86)  ; label("Central volume of distribution, Vc (L)")             # Table 4 Vc = 17.86 L (%RSE 7.22); Results repeat 'Vc was 17.86 L (%RSE, 7.22)'
    lq  <- log(2.63)   ; label("Intercompartmental clearance, Q (L/h)")              # Table 4 Q = 2.63 L/h (%RSE 8.89)
    lvp <- log(99.42)  ; label("Peripheral volume of distribution, Vp (L)")          # Table 4 Vp = 99.42 L (%RSE 15.04)

    # ==================================================================
    # Interindividual variability. Only clearance carries a random
    # effect: 'due to limited data, the popPK model only allowed
    # individual random effect on CL and other PK parameters were fixed
    # effects' (supplement, 'Population pharmacokinetics for
    # adolescents'). MonolixSuite reports omega_CL as the STANDARD
    # DEVIATION of the log-scale random effect, so the variance below is
    # the square of the printed 0.14.
    # ==================================================================
    # var(etalcl) = 0.14^2 = 0.0196, from Table 4 'IIV on CL' = 0.14
    # (RSE 30.72); 14.1 percent CV on the linear scale. The source note
    # sits above the line rather than after it because rxode2 promotes a
    # trailing comment on an unlabelled ini() line into a label().
    etalcl ~ 0.0196

    # ==================================================================
    # Residual unexplained variability. As in the adult model the paper
    # names a 'constant' error model and the supplement an 'additive
    # (constant) error model', applied on the LOG-TRANSFORMED
    # observation: an exponential (log-normal) residual with SD 0.25,
    # not an additive 0.25 ng/mL. Figure 1B settles it -- the
    # prediction-corrected VPC spans 0.1-10 ng/mL with a constant
    # relative spread and its 90% prediction intervals still bracket a
    # median near 0.06 ng/mL at 72 h, which an additive SD four times
    # larger than that median cannot produce. See the vignette
    # 'Assumptions and deviations' section.
    # ==================================================================
    expSd <- 0.25      ; label("Exponential (log-scale constant) residual SD for the plasma concentration (fraction)")  # Table 4 'Constant residual error' = 0.25 (%RSE 10.62); Results, 'The constant residual error was low, estimated at 0.25'
  })

  model({
    # Individual disposition parameters. Only CL carries a random effect.
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    # Micro-constants for supplemental Eqs. S1 and S2 (written in
    # amounts):
    #   dAc/dt = (Q/Vp)*Ap - (Q/Vc)*Ac - (CL/Vc)*Ac
    #   dAp/dt = (Q/Vc)*Ac - (Q/Vp)*Ap
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Two-compartment disposition; the zero-order input is the clinical
    # intravenous infusion and is supplied by the event table.
    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Peptide MASS concentration in blood. No unit scaling is needed despite
    # units$dosing being ug and units$concentration ng/mL: the volumes are in
    # L, and ug / L IS ng/mL.
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
