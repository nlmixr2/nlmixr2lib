Wang_2019a_tacrolimus <- function() {
  description <- paste0(
    "One-compartment population PK model with first-order absorption and ",
    "first-order elimination for twice-daily oral tacrolimus trough ",
    "concentrations in Chinese children with refractory nephrotic syndrome ",
    "(Wang 2019). The absorption rate constant ka is held at 4.48 1/h, ",
    "carried over from earlier paediatric tacrolimus literature, because ",
    "every concentration in the dataset is a pre-dose trough and the ",
    "absorption phase is therefore uninformed; the authors confirmed by ",
    "sensitivity analysis that varying ka five-fold in either direction ",
    "(0.896 to 22.4 1/h) left CL/F, V/F and the objective function nearly ",
    "unchanged. Apparent oral clearance CL/F carries three exponential ",
    "covariate effects entered on UNCENTERED covariates: age in years ",
    "(coefficient +0.0323 per year), serum cystatin C in mg/L (coefficient ",
    "-0.359 per mg/L) and the patient's own total daily tacrolimus dose in ",
    "mg/day (coefficient +0.148 per mg/day). Apparent volume of ",
    "distribution V/F carries no covariate. Inter-individual variability is ",
    "exponential and diagonal on CL/F and V/F; the V/F variance is ",
    "essentially zero (0.002), so V/F is effectively a population constant. ",
    "Residual error is a mixed proportional-plus-additive model. Wang 2019 ",
    "also tabulates eleven previously published paediatric tacrolimus popPK ",
    "models from the literature (its Table V); those are other authors' ",
    "models reproduced in a summary table and are not part of this model ",
    "file."
  )
  reference <- paste0(
    "Wang D, Lu J, Li Q, Li Z. Population pharmacokinetics of tacrolimus in ",
    "pediatric refractory nephrotic syndrome and a summary of other ",
    "pediatric disease models. Exp Ther Med. 2019;17(5):4023-4031. ",
    "doi:10.3892/etm.2019.7446"
  )
  vignette <- "Wang_2019a_tacrolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = TRUE: Wang 2019 Methods 'Analytical method'
  # states "Whole blood concentrations of TAC were measured using the Emit
  # 2000 Tacrolimus assay".
  compartmentData <- list(
    depot = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    AGE = list(
      description = "Patient age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Time-fixed at the value recorded with each therapeutic-drug-",
        "monitoring episode. Enters CL/F as an UNCENTERED exponential term ",
        "exp(0.0323 * AGE), so CL/F rises 3.3% per additional year of age; ",
        "over the observed 2.4-16.4 year span that is a 1.6-fold range. ",
        "Because the term is uncentered, the typical-value theta 5.46 L/h is ",
        "the clearance at AGE = 0, which is outside the studied range and is ",
        "not a physically meaningful value on its own. Wang 2019 Table I: ",
        "mean +/- SD 7.61 +/- 3.92 years, median 6.8 years, range 2.4-16.4 ",
        "years. Wang 2019 Table III: including age on CL/F dropped the ",
        "objective function by 17.160 points, the single largest covariate ",
        "effect in the forward step; removing it raised the objective ",
        "function by 7.988 points in backward elimination. The Discussion ",
        "attributes the effect to developmental maturation of tacrolimus ",
        "clearance and cites a concordant age dependence reported in ",
        "paediatric haematopoietic stem cell transplant recipients."
      ),
      source_name = "age"
    ),
    CYSC = list(
      description = "Serum cystatin C",
      units = "mg/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Renal-function biomarker measured as part of routine care. Enters ",
        "CL/F as an UNCENTERED exponential term exp(-0.359 * CYSC), so CL/F ",
        "falls 30.1% per additional 1 mg/L of cystatin C; over the observed ",
        "0.4-2.3 mg/L span that is a 2.0-fold range. Wang 2019 Table I: mean ",
        "+/- SD 0.85 +/- 0.25 mg/L, median 0.8 mg/L, range 0.4-2.3 mg/L. ",
        "Wang 2019 Table III: including cystatin C on CL/F dropped the ",
        "objective function by 11.282 points; removing it raised the ",
        "objective function by 13.180 points, the largest backward-",
        "elimination penalty of the three retained covariates. Note the sign ",
        "convention: higher cystatin C means WORSE glomerular filtration, so ",
        "a negative coefficient means poorer renal function is associated ",
        "with lower apparent tacrolimus clearance. Tacrolimus is not ",
        "appreciably renally cleared, so Wang 2019 interprets cystatin C ",
        "here as a marker of nephrotic-syndrome disease progression rather ",
        "than as a direct renal-elimination covariate (Discussion: cystatin ",
        "C 'could predict the disease progress' in nephrotic syndrome, and ",
        "'the progression of disease had an impact on CL/F'). This is a ",
        "different mechanistic role from the CYSC uses registered for ",
        "renally-eliminated drugs such as vancomycin and cefuroxime."
      ),
      source_name = "CYSC"
    ),
    DOSE_TAC_MGD = list(
      description = "Patient's own total daily tacrolimus dose",
      units = "mg/day",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Total daily tacrolimus dose summed across the twice-daily ",
        "administrations, updated whenever the prescriber adjusts the dose. ",
        "Enters CL/F as an UNCENTERED exponential term ",
        "exp(0.148 * DOSE_TAC_MGD), so CL/F rises 16.0% per additional ",
        "1 mg/day; over the observed 1.0-4.0 mg/day span that is a 1.6-fold ",
        "range. Wang 2019 Table I: mean +/- SD 1.62 +/- 0.75 mg/day, median ",
        "1.5 mg/day, range 1.0-4.0 mg/day. Wang 2019 Table III: including ",
        "the daily dose on CL/F dropped the objective function by 7.255 ",
        "points; removing it raised the objective function by the same ",
        "7.255 points. IMPORTANT interpretive caveat, stated by the authors ",
        "themselves: this is a dose-on-clearance covariate in a therapeutic-",
        "drug-monitoring dataset where doses were titrated to a trough ",
        "target, so it is confounded by dose individualisation rather than ",
        "being a demonstrated pharmacokinetic nonlinearity. Wang 2019 ",
        "Discussion attributes it to CYP3A5 genotype: patients carrying the ",
        "functional CYP3A5*1 allele clear tacrolimus faster and are ",
        "therefore titrated to higher daily doses, so daily dose acts as a ",
        "surrogate for an unmeasured genotype ('the effect of TAMT on CL/F ",
        "may be primarily derived from CYP3A5 gene polymorphisms ... at ",
        "present, CYP3A5 genotyping is not routinely performed in Chinese ",
        "patients with PRNS'). When simulating, this column must be kept ",
        "consistent with the event table: for the twice-daily regimen the ",
        "per-administration amt is half of DOSE_TAC_MGD. Simulations that ",
        "vary the dose without updating this column will silently lose the ",
        "covariate effect."
      ),
      source_name = "TAMT"
    )
  )

  # Covariates that Wang 2019 collected and screened in the stepwise
  # covariate search (Methods 'Covariate model') but did NOT retain in the
  # final model: only age, cystatin C and daily tacrolimus dose, all on CL/F,
  # met the inclusion and elimination criteria, and no covariate at all was
  # retained on V/F. Documented here so the provenance of the covariate
  # screen is preserved without declaring covariates that model() never
  # references. Globulin, albumin/globulin ratio, gamma-glutamyl
  # transpeptidase, uric acid, mean corpuscular haemoglobin, mean corpuscular
  # haemoglobin concentration and the 21 screened concomitant medications
  # have no canonical entry in inst/references/covariate-columns.md and are
  # recorded in population$notes instead.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened but not retained on either CL/F or V/F. Wang 2019 Table I: mean +/- SD 30.85 +/- 17.12 kg, median 25.0 kg, range 13.5-86.5 kg. Note that no allometric size scaling appears in the final model at all, which is unusual for a paediatric popPK model spanning a 6.4-fold weight range; age was retained instead and the two are strongly collinear in children. This is called out in the vignette's Assumptions and deviations section."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened but not retained. Wang 2019 Table I: mean +/- SD 25.41 +/- 8.87 g/L, median 24.1 g/L, range 12.3-45.3 g/L. The cohort is markedly hypoalbuminaemic, as expected in nephrotic syndrome; tacrolimus is extensively bound to erythrocytes and to alpha-1-acid glycoprotein rather than to albumin, which is consistent with albumin not being retained."
    ),
    TPRO = list(
      description = "Total serum protein",
      units = "g/L",
      type = "continuous",
      notes = "Screened but not retained. Wang 2019 Table I: mean +/- SD 47.51 +/- 10.22 g/L, median 46.9 g/L, range 29.5-69.1 g/L."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened but not retained. Wang 2019 Table I reports the unit as IU/L, which is used interchangeably with U/L: mean +/- SD 15.93 +/- 6.49, median 14.0, range 5.0-35.0."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened but not retained. Wang 2019 Table I reports the unit as IU/L, which is used interchangeably with U/L: mean +/- SD 9.91 +/- 6.48, median 8.0, range 2.0-35.0."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Screened but not retained, in contrast to cystatin C which was retained. Wang 2019 Table I: mean +/- SD 30.49 +/- 12.67 umol/L, median 27.0 umol/L, range 14.0-69.0 umol/L. The Discussion argues explicitly that cystatin C is the better renal and disease-progression marker in primary nephrotic syndrome, being more sensitive than serum creatinine for predicting renal dysfunction in this population."
    ),
    BUN = list(
      description = "Blood urea",
      units = "mmol/L",
      type = "continuous",
      notes = "Screened but not retained. Wang 2019 Table I reports urea (UR): mean +/- SD 4.41 +/- 2.59 mmol/L, median 4.0 mmol/L, range 1.9-18.1 mmol/L. Reported as urea rather than as blood urea nitrogen; divide by 2.14 to convert mmol/L urea to mg/dL BUN if a BUN-scaled dataset is required."
    ),
    HCT = list(
      description = "Haematocrit",
      units = "%",
      type = "continuous",
      notes = "Screened but not retained, which is notable because haematocrit is a retained covariate on tacrolimus apparent clearance in several other registered tacrolimus models (tacrolimus is extensively erythrocyte-bound). Wang 2019 Table I: mean +/- SD 42.62 +/- 4.94%, median 42.6%, range 27.4-55.3%."
    ),
    HGB = list(
      description = "Haemoglobin",
      units = "g/L",
      type = "continuous",
      notes = "Screened but not retained. Wang 2019 Table I: mean +/- SD 144.79 +/- 17.34 g/L, median 146.0 g/L, range 90.0-180.1 g/L."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 65L,
    n_studies = 1L,
    n_observations = 147L,
    age_range = "2.4-16.4 years",
    age_median = "6.8 years",
    age_mean_sd = "7.61 (3.92) years",
    weight_range = "13.5-86.5 kg",
    weight_median = "25.0 kg",
    weight_mean_sd = "30.85 (17.12) kg",
    sex_female_pct = 32.3,
    race_ethnicity = "Chinese (single-centre cohort, Children's Hospital of Fudan University, Shanghai). Wang 2019 describes the patients as Chinese throughout; no further ancestry breakdown is reported.",
    disease_state = paste0(
      "Children under 18 years of age with refractory nephrotic syndrome ",
      "receiving oral tacrolimus, treated at the Children's Hospital of ",
      "Fudan University (Shanghai) between January 2014 and October 2017 ",
      "and analysed retrospectively. Patients with other serious diseases, ",
      "including kidney transplantation, were excluded -- this is ",
      "specifically NOT a transplant cohort, which distinguishes it from ",
      "the majority of paediatric tacrolimus popPK models. 44 male, 21 ",
      "female."
    ),
    dose_range = paste0(
      "Oral tacrolimus capsules (1 mg and 0.5 mg strengths). Initial dose ",
      "0.5-2.0 mg twice daily; total daily dose across the cohort 1.0-4.0 ",
      "mg/day (Table I mean +/- SD 1.62 +/- 0.75 mg/day, median 1.5 ",
      "mg/day). Doses were adjusted on efficacy, adverse effects and the ",
      "therapeutic-drug-monitoring trough concentration."
    ),
    regions = "China (Children's Hospital of Fudan University, Shanghai).",
    sampling_window = paste0(
      "Routine therapeutic drug monitoring, retrospective. 147 whole-blood ",
      "tacrolimus concentrations from 65 patients (a mean of 2.3 samples per ",
      "patient). ALL concentrations are pre-dose troughs: Wang 2019 Methods ",
      "'Drug administration' states 'All blood concentrations were collected ",
      "prior to the subsequent administration. The TAC concentrations used ",
      "in the current research were trough concentrations.' This is why ka ",
      "could not be estimated and why bioavailability and an absorption lag ",
      "time were not estimable."
    ),
    assay = paste0(
      "Whole-blood tacrolimus by Emit 2000 Tacrolimus assay (Siemens ",
      "Healthineers, Erlangen, Germany), linear over 2.0-30.0 ng/mL; ",
      "samples above 30.0 ng/mL were diluted per the manufacturer's ",
      "protocol. Wang 2019 does not report how, or whether, concentrations ",
      "below the lower limit of 2.0 ng/mL were handled."
    ),
    notes = paste0(
      "Single-centre retrospective therapeutic-drug-monitoring study. ",
      "Estimation was by NONMEM 7 FOCE with interaction; internal ",
      "validation was a 2000-replicate bootstrap, of which 1791 minimised ",
      "successfully with an acceptable covariance step. No external ",
      "validation and no visual predictive check were performed. ",
      "Additional covariates screened but not retained, and having no ",
      "canonical entry in inst/references/covariate-columns.md, were: ",
      "globulin (Table I 22.16 +/- 3.32 g/L), the albumin/globulin ratio ",
      "(1.16 +/- 0.44), gamma-glutamyl transpeptidase (32.85 +/- 54.52 ",
      "IU/L, range 9.0-446.0), uric acid (343.42 +/- 117.00 umol/L), mean ",
      "corpuscular haemoglobin (28.91 +/- 1.46 pg) and mean corpuscular ",
      "haemoglobin concentration (340.12 +/- 14.91 g/L). Twenty-one ",
      "concomitant medications were also screened and none was retained ",
      "(Table II); the most frequently co-administered were corticosteroids ",
      "(64/65, 98.5%), spironolactone (22/65, 33.8%), ",
      "dihydrochlorothiazide (23/65, 35.4%), fosinopril (13/65, 20.0%) and ",
      "omeprazole (10/65, 15.4%). Wang 2019 Table II is internally ",
      "inconsistent in its PDF rendering -- several rows have their counts ",
      "and percentages transposed or merged across category lines -- so the ",
      "individual conmed frequencies should be treated as approximate. The ",
      "co-administration frequencies do not affect the model, since no ",
      "conmed was retained."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Fixed-effect estimates: Wang 2019 Table IV 'Parameter estimates of
    # final model and bootstrap validation'. The final-model equations are
    # printed identically in three places -- the Abstract, the Results text
    # immediately below Table IV, and the 'Refractory nephrotic syndrome /
    # Current study' row of Table V:
    #
    #   CL/F (l/h) = 5.46 * EXP(0.0323 * AGE) * EXP(-0.359 * CYSC)
    #                     * EXP(0.148 * TAMT)
    #   V/F (l)    = 57.1
    #   Ka (1/h)   = 4.48 (fixed)
    #
    # NOTE that the three covariate terms are UNCENTERED. 5.46 L/h is
    # therefore the clearance of a hypothetical subject with AGE = 0,
    # CYSC = 0 and TAMT = 0, none of which occurs in the data; it is an
    # equation intercept, not a typical value for this cohort. Substituting
    # the Table I medians (AGE 6.8 y, CYSC 0.8 mg/L, TAMT 1.5 mg/day) gives
    # CL/F = 6.37 L/h for the median child.
    # ---------------------------------------------------------------------
    lcl <- log(5.46) ; label("Apparent oral clearance CL/F equation intercept, at AGE = 0, CYSC = 0, DOSE_TAC_MGD = 0 (L/h)")  # Wang 2019 Table IV CL/F = 5.4600 L/h (SE 22.7%; bootstrap median 5.640)
    lvc <- log(57.1) ; label("Apparent volume of distribution V/F (L)")                                                        # Wang 2019 Table IV V/F = 57.1000 L (SE 46.8%; bootstrap median 59.500)

    # ka was NOT estimated. Wang 2019 Methods: "The absorption rate constant
    # (Ka) of the model was set as 4.48 h-1, according to what was
    # previously set in the literature (28,40,41)", and Table IV lists it as
    # "4.4800 (fixed)" with N/A in every uncertainty column. Results adds a
    # sensitivity analysis: ka was varied five-fold either way, 0.896 to
    # 22.4 1/h, and "the results of CL/F, V/F and the OFV exhibited minimal
    # changes" -- as expected, because every observation is a pre-dose
    # trough and carries almost no information about absorption.
    lka <- fixed(log(4.48)) ; label("Absorption rate constant ka, carried over from earlier paediatric tacrolimus literature (1/h)")  # Wang 2019 Table IV Ka = 4.4800 (fixed)

    # ---------------------------------------------------------------------
    # Covariate coefficients on CL/F, all exponential and all on uncentered
    # covariates. Wang 2019 Table IV rows 'theta AGE', 'theta CYSC',
    # 'theta TAMT'.
    # ---------------------------------------------------------------------
    e_age_cl <- 0.0323       ; label("Exponential coefficient of AGE on CL/F (per year)")                  # Wang 2019 Table IV theta AGE = 0.0323 (SE 35.0%; bootstrap 95% CI 0.007 to 0.062)
    e_cysc_cl <- -0.359      ; label("Exponential coefficient of CYSC on CL/F (per mg/L)")                 # Wang 2019 Table IV theta CYSC = -0.3590 (SE 26.1%; bootstrap 95% CI -0.719 to -0.087)
    e_dose_tac_cl <- 0.148   ; label("Exponential coefficient of DOSE_TAC_MGD on CL/F (per mg/day)")       # Wang 2019 Table IV theta TAMT = 0.1480 (SE 47.9%; bootstrap 95% CI 0.012 to 0.350)

    # ---------------------------------------------------------------------
    # Inter-individual variability. Wang 2019 Results: inter-individual
    # variability "best described by exponential" models, i.e.
    # CL_i = CL_typ * exp(eta_i), so the Table IV omega rows are NONMEM
    # $OMEGA diagonal elements -- VARIANCES on the log scale. They are
    # carried through to ini() unchanged.
    #
    # SCALE DECISION (see the vignette's Assumptions and deviations section
    # for the full argument). Wang 2019's Abstract and its Table V row for
    # the current study restate these same numbers multiplied by 100 and
    # labelled as percentages -- 'The inter-individual variability of CL/F
    # and V/F were 22.2 and 0.2%'. Read literally that would make omega a
    # coefficient of variation of 22.2%, i.e. a variance of 0.0481. It is
    # instead the common reporting slip of printing a NONMEM variance with
    # a percent sign. The decisive evidence is Wang 2019's own Table V,
    # which places this study side by side with eleven other paediatric
    # tacrolimus models under a 'BSV CL (%)' column: every other study in
    # that column reports 24.3, 33.0, 33.5, 40.0, 41.9, 48.7, 50.0, 52.1,
    # 54.8 and 55.6%, while the current study's cell reads 22.2. Taking
    # 0.222 as a variance gives sqrt(0.222) = 47.1% CV, which sits in the
    # middle of that literature range; taking it as a CV makes this study a
    # solitary low outlier. Fig. 1B of the paper corroborates: the spread
    # of observations about the population predictions is far wider than a
    # 22% between-subject CV plus a 36% proportional residual could
    # generate.
    etalcl ~ 0.222   # Wang 2019 Table IV omega CL/F = 0.2220 (SE 18.5%; bootstrap median 0.216); = 47.1% CV
    etalvc ~ 0.002   # Wang 2019 Table IV omega V/F = 0.0020 (SE 48.5%; bootstrap median 0.001); = 4.5% CV, effectively no IIV on V/F

    # ---------------------------------------------------------------------
    # Residual unexplained variability. Wang 2019 Results: residual
    # variability was "best described by ... mixed error models", and the
    # Table IV footnote names the two components explicitly -- 'sigma 1,
    # residual variability, proportional error; sigma 2, residual
    # variability, additive error'. As NONMEM $SIGMA diagonal elements
    # these are variances, so nlmixr2's standard-deviation-scaled prop()
    # and add() take their square roots:
    #   propSd = sqrt(0.359) = 0.5992  (59.9% proportional)
    #   addSd  = sqrt(0.804) = 0.8967  ng/mL
    # The same variance-versus-CV argument given above for the omegas
    # applies here; the two sigma rows are printed in the identical column
    # of the identical table.
    # ---------------------------------------------------------------------
    propSd <- 0.5992 ; label("Proportional residual error (fraction)")  # Wang 2019 Table IV sigma 1 = 0.3590 (SE 8.2%), a variance; sqrt(0.3590) = 0.59917
    addSd <- 0.8967  ; label("Additive residual error (ng/mL)")         # Wang 2019 Table IV sigma 2 = 0.8040 (SE 31.5%), a variance; sqrt(0.8040) = 0.89666
  })

  model({
    # Absorption rate constant, fixed and without inter-individual
    # variability (Wang 2019 estimates no eta on ka -- Table IV lists omega
    # rows for CL/F and V/F only).
    ka <- exp(lka)

    # Wang 2019 final-model clearance equation, reproduced exactly as
    # printed (Abstract, Results below Table IV, and Table V current-study
    # row):
    #   CL/F = 5.46 * EXP(0.0323 * AGE) * EXP(-0.359 * CYSC)
    #              * EXP(0.148 * TAMT)
    # All three covariates enter uncentered, so no reference values appear.
    cl <- exp(lcl + etalcl) *
      exp(e_age_cl * AGE) *
      exp(e_cysc_cl * CYSC) *
      exp(e_dose_tac_cl * DOSE_TAC_MGD)

    # Wang 2019: "No significant effects of covariates on V/F were
    # observed", so V/F carries its typical value and its (near-zero) eta
    # only.
    vc <- exp(lvc + etalvc)

    # One-compartment oral disposition with first-order absorption and
    # first-order elimination. Dose lands in `depot`; bioavailability is
    # not separately identifiable and is absorbed into the apparent CL/F
    # and V/F, per Wang 2019 Methods ("The bioavailability (F) and
    # absorption with a lag time could not be estimated ... Thus, the PK
    # parameters were comprised of apparent oral clearance (CL/F) and
    # apparent volume of distribution (V/F)").
    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Tacrolimus whole-blood concentrations are reported in ng/mL. Doses are
    # in mg and vc is in L, so central/vc is mg/L = ug/mL; multiply by 1000
    # to reach ng/mL.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd) + add(addSd)
  })
}
