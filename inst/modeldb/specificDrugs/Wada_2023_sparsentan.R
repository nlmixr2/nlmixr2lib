Wada_2023_sparsentan <- function() {
  description <- paste(
    "Two-compartment population PK model with lagged first-order absorption for",
    "oral sparsentan (a single-molecule dual endothelin / angiotensin II receptor",
    "antagonist, DEARA) pooled over nine phase I-III studies in 446 subjects:",
    "236 healthy volunteers, 16 subjects with hepatic impairment, and 194 patients",
    "with primary or genetic focal segmental glomerulosclerosis (FSGS).",
    "Apparent oral clearance carries a first-order CYP3A auto-induction term that",
    "steps CL/F from 3.88 L/h on the first dosing day to 5.11 L/h thereafter",
    "(induction half-life fixed at 0.001 day, i.e. effectively instantaneous),",
    "multiplied by power effects of alkaline phosphatase and creatinine clearance",
    "and by log-additive effects of male sex and of moderate / strong CYP3A4",
    "inhibitor coadministration. Relative bioavailability is a power function of",
    "dose clamped at 200 mg, producing less-than-dose-proportional exposure",
    "(Frel = 1.41, 1.00, 0.71 at 200, 400, 800 mg). Race shifts the apparent",
    "central volume, and formulation (whole tablet or crushed tablet against the",
    "capsule reference) shifts both the absorption rate constant and the",
    "absorption lag time. Residual error is combined proportional plus additive,",
    "with the additive SD fixed at the 2 ng/mL assay lower limit of quantitation."
  )

  reference <- paste(
    "Wada R, Kleijn HJ, Zhang L, Chen S-C.",
    "Population pharmacokinetic analysis of sparsentan in healthy volunteers and",
    "patients with focal segmental glomerulosclerosis.",
    "CPT Pharmacometrics Syst Pharmacol. 2023;12(8):1080-1092.",
    "doi:10.1002/psp4.12996"
  )

  vignette <- "Wada_2023_sparsentan"

  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    DOSE = list(
      description        = "Administered sparsentan dose at the dose record",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Use case (a) of the DOSE canonical: the per-record administered dose",
        "level drives the dose-dependent relative bioavailability",
        "Frel = (max(DOSE, 200) / 400)^-0.495 (Wada 2023 Table 2 note and",
        "Table S3). The 400 mg reference gives Frel = 1, so the reported CL/F,",
        "Vc/F, Q/F and Vp/F values are the apparent parameters at 400 mg.",
        "The clamp at 200 mg reproduces the paper's lower branch",
        "Frel = (200 / 400)^-0.495 for doses at or below 200 mg. Doses studied",
        "ranged from 50 to 1600 mg (Table S1)."
      ),
      source_name        = "Dose"
    ),
    ALP = list(
      description        = "Serum alkaline phosphatase",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power effect on CL/F normalised to the overall population median of",
        "68 U/L (Wada 2023 Table 1 overall median; Table 2 footnote b names",
        "68 U/L as the reference-subject value). Overall range 25-269 U/L.",
        "Elevated ALKP is associated with lower CL/F; the paper tested ALKP as",
        "one of four hepatic markers and retained only this one."
      ),
      source_name        = "ALKP"
    ),
    CRCL = list(
      description        = paste(
        "Creatinine clearance estimated by the Cockcroft-Gault equation,",
        "NOT BSA-normalised"
      ),
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Raw (absolute) Cockcroft-Gault creatinine clearance in mL/min, not",
        "normalised to 1.73 m^2 - same convention as Delattre_2010_amikacin.R,",
        "Georges_2009_ceftazidime.R and Chen_2023_nemonoxacin.R. Power effect",
        "on CL/F normalised to the overall population median of 112 mL/min",
        "(Wada 2023 Table 1 overall median; Table 2 footnote b). Overall range",
        "26-361 mL/min. Because Cockcroft-Gault takes sex as an input and sex",
        "is also a covariate on CL/F, the paper notes the CrCL association may",
        "be partly confounded by sex (Discussion)."
      ),
      source_name        = "CrCL"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (female) in this paper - see notes",
      notes              = paste(
        "The canonical SEXF orientation is kept (1 = female), but Wada 2023",
        "reports the effect for MALES against a FEMALE reference subject",
        "(Table 2 footnote b: 'the reference subject is a white female').",
        "The published coefficient +0.139 is therefore applied as",
        "exp(e_sexf_cl * (1 - SEXF)) so the verbatim source value is preserved",
        "and CL/F is unchanged for females - the same construction used in",
        "Bajaj_2017_nivolumab.R. Men had ~15% higher CL/F than women.",
        "Cohort 33.2% female (Wada 2023 Table 1)."
      ),
      source_name        = "Sex"
    ),
    RACE_BLACK = list(
      description        = "Black or African American race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (White; the most frequent level and the paper's reference)",
      notes              = paste(
        "Log-additive effect on the apparent central volume Vc/F",
        "(Wada 2023 Table 2 'Effect on Vc'); Black subjects had ~36% higher",
        "Vc/F than White subjects. Cohort 23.5% Black or African American",
        "(Table 1). Paired with RACE_ASIAN; both 0 selects the White",
        "reference. Subjects recorded as 'Multiple' (0.7%) or 'Other' (2.2%)",
        "in Table 1 have no published coefficient and fall in the reference",
        "group here."
      ),
      source_name        = "Race"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (White; the most frequent level and the paper's reference)",
      notes              = paste(
        "Log-additive effect on the apparent central volume Vc/F",
        "(Wada 2023 Table 2 'Effect on Vc'); Asian subjects had ~30% higher",
        "Vc/F than White subjects. Cohort 6.1% Asian (Table 1). Paired with",
        "RACE_BLACK; both 0 selects the White reference."
      ),
      source_name        = "Race"
    ),
    FORM_TABLET = list(
      description        = "Whole (uncrushed) tablet formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (capsule, when FORM_CRUSHED_TABLET is also 0)",
      notes              = paste(
        "Three-level formulation stratification {capsule reference, whole",
        "tablet, crushed tablet} encoded as the two binary indicators",
        "FORM_TABLET and FORM_CRUSHED_TABLET; both 0 selects the 100 mg",
        "capsule that carries the published typical Ka = 0.740 /h and",
        "Tlag = 0.32 h (Wada 2023 Table 2 footnote b). Note that the",
        "comparator here is a capsule, not the non-tablet oral liquid named in",
        "the FORM_TABLET register default. The whole tablet is absorbed more",
        "slowly than the capsule (Ka x exp(-0.306)) with a shorter lag time",
        "(Tlag x exp(-0.269)). Formulation had no effect on Frel, so AUC is",
        "unaffected and Cmax differs by <10% (Wada 2023 Results / Figure 3).",
        "Cohort 52.5% tablet (Table 1)."
      ),
      source_name        = "Formulation"
    ),
    FORM_CRUSHED_TABLET = list(
      description        = "Crushed-tablet (dry, not resuspended) formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (capsule, when FORM_TABLET is also 0)",
      notes              = paste(
        "Second indicator of the three-level formulation stratification (see",
        "FORM_TABLET). The crushed tablet is absorbed faster than either the",
        "whole tablet or the capsule (Ka x exp(0.080)) and has a markedly",
        "shorter lag time (Tlag x exp(-1.175), i.e. 0.099 h). Only study",
        "RTRXRE021103 used the crushed tablet; cohort 8.1% (Table 1). The",
        "crushed-tablet Ka coefficient is imprecise (RSE 159.1%)."
      ),
      source_name        = "Formulation"
    ),
    CONMED_CYP3A4_INH_MOD = list(
      description        = "Concomitant moderate CYP3A4 inhibitor coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no moderate CYP3A4 inhibitor)",
      notes              = paste(
        "Log-additive effect on CL/F (Wada 2023 Table 2 'Effect on CL'):",
        "CL/F x exp(-0.273), which the paper translates to a 31.4% increase in",
        "steady-state AUC and a 16.0% increase in Cmax. Cyclosporine is the",
        "moderate inhibitor named in the Discussion (70% increase in observed",
        "sparsentan AUCinf in the dedicated DDI study 021HVOL16006). Cohort",
        "10.5% (Table 1). Weak CYP3A4 inhibitors (12.8% of the cohort) were",
        "tested and NOT retained, so they fall in the 0 reference group",
        "together with subjects on no inhibitor."
      ),
      source_name        = "CYP3A4 inhibitor (moderate)"
    ),
    CONMED_CYP3A4_INH_STRONG = list(
      description        = "Concomitant strong CYP3A4 inhibitor coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no strong CYP3A4 inhibitor)",
      notes              = paste(
        "Log-additive effect on CL/F (Wada 2023 Table 2 'Effect on CL'):",
        "CL/F x exp(-1.069), a 66% reduction in CL/F, which the paper",
        "translates to a 191.3% increase in steady-state AUC and a 99.0%",
        "increase in Cmax - the largest covariate effect in the model and the",
        "only one for which the paper suggests dose adjustment may be",
        "warranted. Itraconazole is the strong inhibitor named in the",
        "Discussion (174% increase in observed sparsentan AUCinf in study",
        "021HVOL16006). Cohort 7.2% (Table 1)."
      ),
      source_name        = "CYP3A4 inhibitor (strong)"
    )
  )

  # Covariates that Wada 2023 screened (Table S2) but did not retain in the
  # final model (Results 'Final model'; Discussion). Documented here so the
  # provenance of the covariate screen is preserved without carrying unused
  # covariateData entries.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "One of four body-size descriptors tested on CL/F and Vc/F",
        "(Table S2 footnote a); no significant effect (Discussion).",
        "Cohort median 78.6 kg, range 21.1-154.0 kg (Table 1)."
      )
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Body-size descriptor tested on CL/F and Vc/F (Table S2); not retained."
    ),
    BSA = list(
      description = "Body surface area",
      units       = "m^2",
      type        = "continuous",
      notes       = "Body-size descriptor tested on CL/F and Vc/F (Table S2); not retained."
    ),
    LBM = list(
      description = "Lean body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Body-size descriptor tested on CL/F and Vc/F (Table S2); not retained."
    ),
    AGE = list(
      description = "Subject age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Tested on CL/F and Vc/F (Table S2); not retained. Steady-state",
        "exposure did rise across the 8-17 / 18-64 / >=65 year brackets",
        "(Figure S2), which the paper attributes to age being an input to the",
        "retained CrCL covariate rather than to an independent age effect.",
        "Cohort median 40 years, range 8-74 years (Table 1)."
      )
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/dL",
      type        = "continuous",
      notes       = "Hepatic marker tested on CL/F (Table S2); not retained."
    ),
    TPRO = list(
      description = "Total serum protein",
      units       = "g/dL",
      type        = "continuous",
      notes       = "Tested on CL/F (Discussion); not retained."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = paste(
        "Hepatic marker tested on CL/F (Table S2); not retained. Hepatic",
        "impairment was not tested directly as a covariate - ALKP, ALT, AST",
        "and total bilirubin were tested in its place and only ALKP survived",
        "(Discussion)."
      )
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Hepatic marker tested on CL/F (Table S2); not retained."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Hepatic marker tested on CL/F (Table S2); not retained."
    ),
    FED = list(
      description = "Fed-versus-fasted dose-record indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Tested on Frel, Ka and Tlag (Table S2); not retained. Phase I",
        "subjects generally dosed after an overnight fast; FSGS patients were",
        "instructed to dose before the first meal of the day, and the analysis",
        "found no significant difference in absorption between the two",
        "(Discussion). Food status was recorded as 'unknown' for all 194 FSGS",
        "patients (Table 1)."
      )
    ),
    CONMED_CYP3A4_IND = list(
      description = "Concomitant CYP3A4 inducer coadministration indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Tested on CL/F (Table S2); not retained. No subject in the cohort",
        "received a weak, moderate or strong CYP3A4 inducer - the 12.1%",
        "flagged in Table 1 are 'unknown' - so the covariate carried no",
        "information."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 446L,
    n_studies      = 9L,
    n_observations = 10957L,
    age_range      = "8-74 years",
    age_median     = "40 years",
    weight_range   = "21.1-154.0 kg",
    weight_median  = "78.6 kg",
    sex_female_pct = 33.2,
    race_ethnicity = c(White = 67.5, Black = 23.5, Asian = 6.1, Multiple = 0.7, Other = 2.2),
    disease_state  = paste(
      "Pooled: 236 healthy volunteers (six phase I studies), 16 subjects with",
      "mild or moderate hepatic impairment (one phase I study), and 194",
      "patients with primary or genetic focal segmental glomerulosclerosis",
      "(phase II DUET and phase III DUPLEX)"
    ),
    renal_function = paste(
      "Normal 74.2%, mild impairment 13.5%, moderate 11.9%, severe 0.4%",
      "(n = 2 severe, so no conclusions were drawn for severe impairment)"
    ),
    hepatic_function = paste(
      "Recorded for 28 subjects only: normal 12, mild 8, moderate 8, severe 0;",
      "missing for 93.7% of the pooled cohort"
    ),
    co_medication  = paste(
      "CYP3A4 inhibitor: none 83.0%, weak 12.8%, moderate 10.5%, strong 7.2%;",
      "P-glycoprotein inhibitor 8.1%; acid-reducing agent 13.2%;",
      "CYP3A4 inducer none 87.9% (remainder unknown)"
    ),
    dose_range     = paste(
      "Single doses 50-1600 mg (capsule, tablet, or crushed tablet) in phase I;",
      "200, 400 or 800 mg once daily in DUET; 400 mg once daily titrating to",
      "800 mg in DUPLEX (Table S1)"
    ),
    notes          = paste(
      "Baseline demographics are Wada 2023 Table 1; the study list is",
      "Table S1. Bioanalysis was validated LC-MS/MS with a 2-4000 ng/mL",
      "calibration range and a 2 ng/mL LLOQ; all BLQ samples were excluded",
      "from the analysis. Three further covariates in Table S2 were screened",
      "and not retained but are not listed in covariatesDataExcluded because",
      "the register has no canonical column for them: population",
      "(FSGS versus non-FSGS, tested on CL/F and Vc/F), acid-reducing agent",
      "(tested on Frel, Ka, Tlag) and P-glycoprotein inhibitor",
      "(tested on CL/F)."
    )
  )

  ini({
    # Structural parameters. Every value below is the estimate for the
    # reference subject defined in Wada 2023 Table 2 footnote b: a White
    # FEMALE receiving a 400 mg CAPSULE, no CYP3A4 inhibitor, CrCL
    # 112 mL/min, ALKP 68 U/L. Frel = 1 at 400 mg, so these are the apparent
    # (per-unit-bioavailability) parameters at the 400 mg dose level.

    # Apparent oral clearance before CYP3A auto-induction. This is the value
    # the paper reports "after a single 400-mg dose".
    lcl <- log(3.88); label("Apparent oral clearance CL/F before CYP3A auto-induction (L/h)")  # Wada 2023 Table 2: CL/F = 3.88 L/h (RSE 4.6%); Table S3 first bracketed term

    # Auto-induction arm of clearance. Wada 2023 Methods writes the
    # time-dependent clearance as
    #   CL(t) = CL0 + dCL * (1 - exp(-k_induction * Day))
    # and Table S3 repeats it verbatim inside the final-model equation. The
    # increment dCL is registered here under the canonical time-varying CL
    # component name lcl_time (the same canonical pair used by
    # Wang_2014_vatalanib.R and Hussein_1997_lamotrigine.R for a first-order
    # auto-induction rise). CL0 + dCL = 3.88 + 1.23 = 5.11 L/h, which the
    # paper reports rounded as the 5.12 L/h steady-state CL/F.
    lcl_time <- log(1.23); label("CYP3A auto-induction increment on CL/F at steady state (L/h)")  # Wada 2023 Table 2: Induction change in CL = 1.23 L/h (RSE 13.6%)

    # First-order onset rate of the auto-induction. Wada 2023 fixed the
    # induction half-life at 0.001 day (Table 2, footnote a: "Induction t1/2
    # and SD of additive error are fixed and not estimated") because the time
    # course was not apparent in Cmin; the corresponding rate constant is
    # log(2) / 0.001 = 693.147 /day. Results notes that when the half-life was
    # estimated instead it came out at 0.3 day (Discussion says 0.2 day),
    # which the authors took as justifying the arbitrary fixed setting.
    lkdes <- fixed(log(693.147)); label("CYP3A auto-induction first-order onset rate constant (1/day)")  # Wada 2023 Table 2: Induction t1/2 = 0.001 day (fixed) -> log(2)/0.001 = 693.147 /day

    lvc <- log(49.3); label("Apparent central volume of distribution Vc/F (L)")            # Wada 2023 Table 2: Vc/F = 49.3 L (RSE 4.3%)
    lq  <- log(2.03); label("Apparent intercompartmental clearance Q/F (L/h)")             # Wada 2023 Table 2: Q/F = 2.03 L/h (RSE 12.0%)
    lvp <- log(12.1); label("Apparent peripheral volume of distribution Vp/F (L)")         # Wada 2023 Table 2: Vp/F = 12.1 L (RSE 10.5%)
    lka <- log(0.740); label("First-order absorption rate constant Ka, capsule (1/h)")     # Wada 2023 Table 2: Ka = 0.740 /h (RSE 6.9%)
    ltlag <- log(0.32); label("Absorption lag time Tlag, capsule (h)")                     # Wada 2023 Table 2: Tlag = 0.32 h (RSE 4.0%)

    # Dose effect on relative bioavailability. Table 2 note:
    #   Frel = (dose/400)^-0.495 if dose >= 200 mg
    #   Frel = (200/400)^-0.495  if dose <  200 mg
    # (Table S3 states the branch boundary as > 200 / <= 200 mg; the two
    # readings agree because both branches give the same value at exactly
    # 200 mg.) The negative exponent means Frel FALLS as dose rises, i.e.
    # less-than-dose-proportional exposure: 1.41, 1.00 and 0.71 at 200, 400
    # and 800 mg respectively (Results 'Final model').
    e_dose_fdepot <- -0.495; label("Power exponent of dose on relative bioavailability Frel (unitless)")  # Wada 2023 Table 2: Dose on Frel = -0.495 (RSE 5.1%)

    # Covariate effects on CL/F. Wada 2023 Methods: "Covariates were added in
    # the log-domain. Continuous covariates were incorporated using a scaled
    # structure based on the median (population) or standard value of the
    # covariate" -> power form (COV / median)^theta. "Categorical covariates
    # were incorporated using the most frequent level of the covariate as the
    # reference" -> log-additive form exp(theta * indicator).
    e_alp_cl  <- -0.208; label("Power exponent of alkaline phosphatase (/68 U/L) on CL/F (unitless)")      # Wada 2023 Table 2 'Effect on CL' ALKP = -0.208 (RSE 27.5%); Table S3 (ALKP/68)^-0.208
    e_crcl_cl <-  0.222; label("Power exponent of creatinine clearance (/112 mL/min) on CL/F (unitless)")  # Wada 2023 Table 2 'Effect on CL' CrCL = 0.222 (RSE 26.5%); Table S3 (CrCL/112)^0.222

    # Sex effect. The published coefficient is for MALES against a FEMALE
    # reference subject (Table 2 footnote b), so it enters as
    # exp(e_sexf_cl * (1 - SEXF)) to preserve the verbatim +0.139 while
    # keeping the canonical SEXF orientation (1 = female).
    e_sexf_cl <- 0.139; label("Log-scale effect of MALE sex on CL/F, female reference (unitless)")  # Wada 2023 Table 2 'Effect on CL' Male = 0.139 (RSE 32.8%); Table S3 (x exp[0.139] if Male)

    e_conmed_cyp3a4_inh_mod_cl    <- -0.273; label("Log-scale effect of moderate CYP3A4 inhibitor coadministration on CL/F (unitless)")  # Wada 2023 Table 2 'Effect on CL' Moderate CYP3A4 = -0.273 (RSE 18.8%)
    e_conmed_cyp3a4_inh_strong_cl <- -1.069; label("Log-scale effect of strong CYP3A4 inhibitor coadministration on CL/F (unitless)")    # Wada 2023 Table 2 'Effect on CL' Strong CYP3A4 = -1.069 (RSE 10.0%)

    # Covariate effects on Vc/F. White is the reference race.
    e_race_black_vc <- 0.309; label("Log-scale effect of Black or African American race on Vc/F (unitless)")  # Wada 2023 Table 2 'Effect on Vc' Black or African American = 0.309 (RSE 18.4%)
    e_race_asian_vc <- 0.265; label("Log-scale effect of Asian race on Vc/F (unitless)")                      # Wada 2023 Table 2 'Effect on Vc' Asian = 0.265 (RSE 48.4%)

    # Formulation effects on absorption. Capsule is the reference.
    e_form_tablet_ka          <- -0.306; label("Log-scale effect of whole tablet on Ka, capsule reference (unitless)")    # Wada 2023 Table 2 'Effect on Ka' Tablet = -0.306 (RSE 34.8%)
    e_form_crushed_tablet_ka  <-  0.080; label("Log-scale effect of crushed tablet on Ka, capsule reference (unitless)")  # Wada 2023 Table 2 'Effect on Ka' Crushed tablet = 0.080 (RSE 159.1%)
    e_form_tablet_tlag         <- -0.269; label("Log-scale effect of whole tablet on Tlag, capsule reference (unitless)")   # Wada 2023 Table 2 'Effect on Tlag' Tablet = -0.269 (RSE 29.1%)
    e_form_crushed_tablet_tlag <- -1.175; label("Log-scale effect of crushed tablet on Tlag, capsule reference (unitless)") # Wada 2023 Table 2 'Effect on Tlag' Crushed tablet = -1.175 (RSE 28.7%)

    # Inter-individual variability. Wada 2023 Methods: log-normal IIV,
    # theta_i = theta_typical * exp(eta_i) with eta ~ N(0, omega^2). Table 2
    # reports both the omega^2 estimate (the "Variance" rows) and its square
    # root as "IIV (%)": sqrt(0.156) = 39.5%, sqrt(0.234) = 48.4%,
    # sqrt(0.474) = 68.9% - so the reported IIV% is sqrt(omega^2), NOT the
    # log-normal CV. The variances are therefore used verbatim, with no
    # omega^2 = log(1 + CV^2) conversion. No IIV was retained on Q/F, Vp/F or
    # Tlag (Table 2 reports "-" or "NA" for those rows), and no correlations
    # between etas are reported, so OMEGA is diagonal.
    etalcl ~ 0.156  # Wada 2023 Table 2: Variance CL = 0.156 (RSE 8.7%); IIV 39.5%, shrinkage 4.7%
    etalvc ~ 0.234  # Wada 2023 Table 2: Variance Vc = 0.234 (RSE 11.3%); IIV 48.4%, shrinkage 17.2%
    etalka ~ 0.474  # Wada 2023 Table 2: Variance Ka = 0.474 (RSE 10.2%); IIV 68.9%, shrinkage 22.2%

    # Residual error: proportional plus additive, with the additive SD fixed
    # at the assay LLOQ. Table 2 labels both rows explicitly as SDs.
    propSd <- 0.365; label("Proportional residual error SD (fraction)")  # Wada 2023 Table 2: SD of proportional error = 0.365 (RSE 1.9%)
    addSd  <- fixed(0.002); label("Additive residual error SD (ug/mL)")  # Wada 2023 Table 2: SD of additive error = 2 ng/mL (fixed at the LLOQ) = 0.002 ug/mL in this model's mg/L scale
  })

  model({
    # ---- 1. Derived covariate terms ---------------------------------------
    # Dose-dependent relative bioavailability (Wada 2023 Table 2 note,
    # Table S3). Clamping the dose at 200 mg reproduces both published
    # branches with a single expression.
    dose_frel <- max(DOSE, 200)
    frel <- (dose_frel / 400)^e_dose_fdepot

    # CYP3A auto-induction of clearance (Wada 2023 Methods; Table S3):
    #   CL(t) = CL0 + dCL * (1 - exp(-k_induction * Day))
    # `Day` in the published equation is the whole-day index since the first
    # dose, taken as 0 on the first dosing day. With the induction half-life
    # fixed at 0.001 day this makes the term a clean step: CL/F = 3.88 L/h
    # throughout the first dosing day and 5.11 L/h from the second dose
    # onward - exactly the paper's "rapid increase to steady-state occurring
    # after the first dose during multiple-dose regimens" (Results, base
    # model) and "occurring instantly after the first dose" (Discussion), and
    # the only reading under which the reported single-dose CL/F of 3.88 L/h
    # is distinguishable from the 5.12 L/h steady-state value. Here it is
    # derived from the simulation clock, which assumes time = 0 is the first
    # dose; see the vignette's Assumptions section.
    day <- floor(time / 24)   # 24 h per day; model time unit is hours
    kdes <- exp(lkdes)        # 1/day
    induction <- 1 - exp(-kdes * day)

    # ---- 2. Individual PK parameters --------------------------------------
    # The eta on clearance multiplies the whole typical-value expression
    # (induction bracket times covariate effects), matching the paper's
    # theta_i = theta_typical * exp(eta_i) parameterisation of Table S3.
    cl_base <- exp(lcl)
    cl_time <- exp(lcl_time) * induction
    cl <-
      (cl_base + cl_time) *
      (ALP / 68)^e_alp_cl *
      (CRCL / 112)^e_crcl_cl *
      exp(e_sexf_cl * (1 - SEXF)) *
      exp(e_conmed_cyp3a4_inh_mod_cl * CONMED_CYP3A4_INH_MOD) *
      exp(e_conmed_cyp3a4_inh_strong_cl * CONMED_CYP3A4_INH_STRONG) *
      exp(etalcl)

    vc <-
      exp(lvc + etalvc) *
      exp(e_race_black_vc * RACE_BLACK) *
      exp(e_race_asian_vc * RACE_ASIAN)

    q  <- exp(lq)
    vp <- exp(lvp)

    ka <-
      exp(lka + etalka) *
      exp(e_form_tablet_ka * FORM_TABLET) *
      exp(e_form_crushed_tablet_ka * FORM_CRUSHED_TABLET)

    tlag <-
      exp(ltlag) *
      exp(e_form_tablet_tlag * FORM_TABLET) *
      exp(e_form_crushed_tablet_tlag * FORM_CRUSHED_TABLET)

    # ---- 3. Micro-constants ----------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- 4. ODE system ---------------------------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(central) <-
      ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ---- 5. Bioavailability and absorption lag ---------------------------
    f(depot) <- frel
    alag(depot) <- tlag

    # ---- 6. Observation and error model ----------------------------------
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
