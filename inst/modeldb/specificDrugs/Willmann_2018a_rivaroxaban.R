Willmann_2018a_rivaroxaban <- function() {
  description <- "Integrated population PK model for oral rivaroxaban pooled across all approved adult indications (Willmann 2018a; 4,918 patients, 22,843 observations, six phase II and one phase III trial). One-compartment disposition with first-order absorption and first-order elimination parameterized as CL/F, V/F and ka. Relative bioavailability declines with the administered dose as an exponential saturation function between Fmax and Fmin. CL/F carries Tietz-truncated creatinine clearance and body weight in power form, five multiplicative comedication indicators (P-glycoprotein inhibitor, strong / medium / weak CYP3A4 inhibitor, CYP3A4 inducer) and a study / indication factor (atrial fibrillation, acute coronary syndrome, and VTE prevention split at 72 h after the first dose, all relative to the VTE-treatment reference). V/F carries body weight and age in power form and a proportional female shift."
  reference <- paste(
    "Willmann S, Zhang L, Frede M, Kubitza D, Mueck W, Schmidt S, Solms A,",
    "Yan X, Garmann D. Integrated Population Pharmacokinetic Analysis of",
    "Rivaroxaban Across Multiple Patient Populations.",
    "CPT Pharmacometrics Syst Pharmacol. 2018;7(5):309-320.",
    "doi:10.1002/psp4.12288"
  )
  vignette <- "Willmann_2018a_rivaroxaban"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    depot = list(analyte = "rivaroxaban", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "rivaroxaban", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance estimated with the Cockcroft-Gault equation, then Tietz-truncated inside model()",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "RAW Cockcroft-Gault creatinine clearance in mL/min -- NOT",
        "BSA-normalized to mL/min/1.73 m^2. Supply the untruncated",
        "Cockcroft-Gault value; model() applies the Tietz truncation itself",
        "(Willmann 2018a Methods, 'Covariate model development'), capping CRCL",
        "at 140 * BSA / 1.73 mL/min so that implausibly high clearances at",
        "extremes of body surface area cannot drive CL/F. Enters CL/F in power",
        "form as (CRCL_truncated / 93)^0.406 with reference 93 mL/min",
        "(Willmann 2018a Supplementary Eq. 1 and the run12 $PK block in",
        "Supplementary Data S9). Table 1 reports a pooled mean untruncated CrCL",
        "of 97.74 mL/min and a pooled mean truncated CrCL of 97.17 mL/min.",
        "The paper's own discussion notes that weight is carried separately on",
        "CL/F precisely because weight already enters the Cockcroft-Gault",
        "numerator, so the two effects are intentionally confounded and must be",
        "applied together."
      ),
      source_name = "CrCL / CrCLtruncated"
    ),
    BSA = list(
      description = "Body surface area, used only to set the Tietz truncation ceiling on creatinine clearance",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "BSA is not a covariate on any PK parameter. It appears solely in the",
        "Tietz creatinine-clearance ceiling 140 * BSA / 1.73 mL/min",
        "(Willmann 2018a Methods). The run12 $PK block imputes a missing BSA as",
        "1.92 m^2 ('IF (BSA.LT.0) BSA2=1.92'), which matches the pooled mean of",
        "1.93 m^2 in Table 1; supply a real BSA where available."
      ),
      source_name = "BSA"
    ),
    WT = list(
      description = "Total body weight at baseline",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters both CL/F and V/F in power form with the SAME reference weight",
        "of 81 kg: (WT / 81)^-0.278 on CL/F and (WT / 81)^0.216 on V/F",
        "(Willmann 2018a Table 3 and Supplementary Eq. 1). The CL/F exponent is",
        "NEGATIVE, which is not allometry -- Willmann 2018a Discussion explains",
        "that the weight term on CL/F refines the weight dependence already",
        "introduced through the Cockcroft-Gault creatinine clearance, and notes",
        "that weight on CL/F gave no OFV improvement when tested univariately.",
        "Pooled mean 82.48 kg (SD 16.87) in Table 1."
      ),
      source_name = "WT"
    ),
    AGE = list(
      description = "Age at baseline",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters V/F only, in power form as (AGE / 61)^-0.189 (Willmann 2018a",
        "Table 3 and Supplementary Eq. 1). Age was ALSO tested on CL/F: it was",
        "carried in the base model, never improved the fit (Table 2 run 6,",
        "dOFV 0), and was dropped at backward elimination once age and sex on",
        "V/F and weight on CL/F were in (Table 2 run 10, dOFV 1.3). AGE is",
        "therefore a retained covariate -- it belongs in covariateData, not in",
        "covariatesDataExcluded -- but only on V/F. Pooled mean 60.53 years",
        "(SD 11.82) in Table 1; the analysis population is adults only."
      ),
      source_name = "AGE"
    ),
    SEXF = list(
      description = "Sex indicator: 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male; V/F multiplier 1 by definition).",
      notes = paste(
        "Multiplicative proportional shift on V/F: 0.889^SEXF, i.e. apparent",
        "volume of distribution is 11.1% lower in women than in men of the same",
        "weight and age (Willmann 2018a Table 3, 'hV/F, Sex' = 0.889).",
        "The run12 $PK block writes this as THETA(7)**(SEX-1) on a 1 = male /",
        "2 = female source coding, and annotates THETA(7) '; SEX_V (female)',",
        "so the shift applies to women; SEXF = 1 - (SEX - 1) is NOT the",
        "transformation -- SEXF equals SEX - 1 directly. 39.3% of the pooled",
        "cohort were female (Table 1), ranging 22.2-62.7% across studies."
      ),
      source_name = "SEX (1 = male, 2 = female)"
    ),
    CONMED_PGP_INH = list(
      description = "Concomitant P-glycoprotein inhibitor at the record",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant P-glycoprotein inhibitor).",
      notes = paste(
        "Time-varying per record (Willmann 2018a Supplementary Table S1 lists",
        "co-medication as time-varying). Multiplicative effect on CL/F:",
        "0.966^CONMED_PGP_INH, a 3.4% reduction (Table 3, 'hCL/F, PGP'; RSE",
        "1.73%, bootstrap 95% CI 0.933-1.00). The five comedication indicators",
        "multiply together (run12 $PK: COMM is the product of the five terms),",
        "so a patient on both a P-gp inhibitor and a CYP3A4 inhibitor receives",
        "the product of both factors. Willmann 2018a Discussion and",
        "Supplementary Table S2 caution that these population estimates are far",
        "smaller than the dedicated phase I drug-drug-interaction results",
        "(ketoconazole raises AUC 2.6-fold in a DDI study but only 1.02-fold in",
        "this model), because strong inhibitors were a protocol exclusion and",
        "exposure was brief and sometimes topical."
      ),
      source_name = "PGPINH4"
    ),
    CONMED_CYP3A4_INH_STRONG = list(
      description = "Concomitant strong CYP3A4 inhibitor at the record",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant strong CYP3A4 inhibitor).",
      notes = paste(
        "Time-varying per record. Multiplicative effect on CL/F:",
        "0.978^CONMED_CYP3A4_INH_STRONG (Willmann 2018a Table 3). This is the",
        "least reliable coefficient in the model and should be treated as",
        "effectively null: only 6 of 5,041 patients in the pooled phase II/III",
        "program received a concomitant strong CYP3A4 inhibitor (Discussion,",
        "'Strengths and limitations'), strong inhibitors were excluded by",
        "protocol, and the bootstrap 95% CI spans 0.902-1.97. Supplementary",
        "Table S2 contrasts the 1.02-fold model AUC increase with the 2.6-fold",
        "(ketoconazole) and 2.5-fold (ritonavir) increases seen in dedicated",
        "drug-drug-interaction studies. Do not use this coefficient to predict",
        "a ketoconazole or ritonavir interaction."
      ),
      source_name = "SCYP3A4INH3"
    ),
    CONMED_CYP3A4_INH_MOD = list(
      description = "Concomitant moderate ('medium') CYP3A4 inhibitor at the record",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant moderate CYP3A4 inhibitor).",
      notes = paste(
        "Time-varying per record. Multiplicative effect on CL/F:",
        "0.863^CONMED_CYP3A4_INH_MOD, a 13.7% reduction (Willmann 2018a Table 3,",
        "'hCL/F, Medium CYP3A4 inhibitor'; RSE 3.79%, bootstrap 95% CI",
        "0.793-0.920). The paper writes 'Medium' where the FDA/EMA vocabulary",
        "writes 'moderate'; Supplementary Table S2 names erythromycin and",
        "fluconazole as the moderate-inhibitor exemplars, both of which the",
        "model raises AUC 1.16-fold. This is the largest and best-determined of",
        "the five comedication effects."
      ),
      source_name = "MCYP3A4INH2"
    ),
    CONMED_CYP3A4_INH_WEAK = list(
      description = "Concomitant weak CYP3A4 inhibitor at the record",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant weak CYP3A4 inhibitor).",
      notes = paste(
        "Time-varying per record. Multiplicative effect on CL/F:",
        "0.939^CONMED_CYP3A4_INH_WEAK, a 6.1% reduction (Willmann 2018a Table 3,",
        "'hCL/F, Weak CYP3A4 inhibitor'; RSE 2.17%, bootstrap 95% CI",
        "0.900-0.975). Willmann 2018a is the founding model for this canonical:",
        "it is one of the few popPK analyses that estimates all three CYP3A4",
        "inhibitor strengths as separate simultaneous indicators rather than",
        "pooling the weak stratum into the reference. The paper does not name",
        "the specific agents classified as weak."
      ),
      source_name = "WCYP3A4INH1"
    ),
    CONMED_CYP3A4_IND = list(
      description = "Concomitant CYP3A4 inducer at the record",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant CYP3A4 inducer).",
      notes = paste(
        "Time-varying per record. Multiplicative effect on CL/F:",
        "1.30^CONMED_CYP3A4_IND, a 30% INCREASE in apparent clearance",
        "(Willmann 2018a Table 3, 'hCL/F, CYP3A4 inducer'; RSE 6.30%, bootstrap",
        "95% CI 1.16-1.48) -- the largest comedication effect and the only one",
        "above 1. Inducer strengths are pooled into a single indicator;",
        "Supplementary Table S2 names rifampicin (strong CYP3A4 and P-gp",
        "inducer) as the exemplar, with a 0.5-fold AUC change in the dedicated",
        "drug-drug-interaction study versus 0.77-fold in this model."
      ),
      source_name = "CYP3A4IND5"
    ),
    DIS_AF = list(
      description = "Nonvalvular atrial fibrillation cohort indicator (ROCKET AF, study 3001)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the VTE-treatment reference cohort, studies 11223 ODIXa-DVT and 11528 EINSTEIN DVT, whose CL/F factor is 1 by definition).",
      notes = paste(
        "Multiplicative study / indication effect on CL/F: 0.849^DIS_AF, i.e.",
        "15.1% lower apparent clearance in the AF cohort than in VTE treatment",
        "(Willmann 2018a Table 3, 'hCL/F, AF'; RSE 3.48%, bootstrap 95% CI",
        "0.793-0.900). Exactly one of DIS_AF, DIS_ACS and DIS_VTE_P may be 1 on",
        "a given record; all three 0 selects the VTE-treatment reference.",
        "n = 161 patients / 800 observations from the ROCKET AF phase III",
        "substudy (Table 1). The paper is explicit that it cannot fully explain",
        "the indication effects (Discussion): they are 'likely multifactorial',",
        "mixing genuine health-status differences against study-specific peak /",
        "trough sampling windows that emphasise different parts of the profile."
      ),
      source_name = "STU = 3001"
    ),
    DIS_ACS = list(
      description = "Acute coronary syndrome cohort indicator (ATLAS ACS TIMI-46, study 2001)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the VTE-treatment reference cohort, CL/F factor 1).",
      notes = paste(
        "Multiplicative study / indication effect on CL/F: 1.14^DIS_ACS, i.e.",
        "14% higher apparent clearance than in VTE treatment (Willmann 2018a",
        "Table 3, 'hCL/F, ACS'; RSE 1.93%, bootstrap 95% CI 1.10-1.18).",
        "n = 2,251 patients / 9,376 observations, the single largest cohort in",
        "the pool (Table 1). Mutually exclusive with DIS_AF and DIS_VTE_P."
      ),
      source_name = "STU = 2001"
    ),
    DIS_VTE_P = list(
      description = "VTE-prevention (post elective hip or knee replacement thromboprophylaxis) cohort indicator (ODIXa-Hip2 10944, ODIXa-Knee 10945, ODIXa-OD-Hip 11527)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the VTE-treatment reference cohort, CL/F factor 1).",
      notes = paste(
        "Multiplicative study / indication effect on CL/F that is TIME-VARYING",
        "within the cohort: 1.04 at or before 72 h after the first dose and 1.29",
        "thereafter (Willmann 2018a Table 3, 'hCL/F, VTE <=72 h' and",
        "'hCL/F, VTE >72 h'). The run12 $PK block keys the split on the NONMEM",
        "TIME variable ('IF (TIME.LE.72) INDEFF = THETA(16)'), which in this",
        "dataset is time after the first dose, so model() uses the solver time",
        "t and the event table must place the first dose at t = 0. The",
        "Discussion attributes the rise to clearance of rivaroxaban increasing",
        "over the first three days after major orthopaedic surgery. Note that",
        "the immediate postsurgical phase was EXCLUDED from the analysis dataset",
        "(Methods, 'Studies included in the analyses') because patients there",
        "separate into slow and fast absorbers, so this factor describes the",
        "post-exclusion window only. Supplementary Eq. 1 mistypes the first",
        "study number as 10933; the run12 control stream and Table 1 both give",
        "10944 (ODIXa-Hip2)."
      ),
      source_name = "STU in (10944, 10945, 11527)"
    ),
    DOSE_RIVAROXABAN_MG = list(
      description = "Administered rivaroxaban dose per administration",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Supplied as a data column so the dose-dependent relative",
        "bioavailability can be evaluated inside model(), exactly as the run12",
        "control stream does (a DOSE column in $INPUT, separate from AMT, used",
        "in the F1 expression). This is the dose of a SINGLE administration, not",
        "the daily total: Table 1 lists the twice-daily arms as, for example,",
        "'2.5, 5, 10, 20, 30 b.i.d.', and the saturation being modelled is the",
        "solubility-limited absorption of one dose. Set it equal to the amt of",
        "the corresponding dose record. Doses in the pool range from 2.5 mg o.d.",
        "to 30 mg b.i.d."
      ),
      source_name = "DOSE"
    )
  )

  # Screened per Willmann 2018a Supplementary Table S1 but not retained in the
  # final model, so they are documented rather than carried in covariateData.
  covariatesDataExcluded <- list(
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = paste(
        "Listed in Supplementary Table S1 as a derived size covariate on CL and",
        "V, offered as an alternative to body weight; body weight was retained",
        "instead. BMI is still used OUTSIDE the model, to define the virtual",
        "subpopulations of the exposure simulation (<18.5, 18.5 to <25, 25 to",
        "<30, 30 to <40, >=40 kg/m^2)."
      )
    ),
    LBM = list(
      description = "Lean body mass",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Listed in Supplementary Table S1 as a derived size covariate on CL and",
        "V and carried in the analysis dataset (an LBM column appears in the",
        "run12 $INPUT and $TABLE lists), but it does not enter any PK parameter",
        "in the final model. Pooled mean 57.05 kg (SD 10.00) in Table 1."
      )
    ),
    SCR = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      notes = paste(
        "Listed in Supplementary Table S1 as a laboratory covariate on CL.",
        "Renal function entered the final model through Tietz-truncated",
        "creatinine clearance instead, which is itself a function of serum",
        "creatinine via Cockcroft-Gault. Pooled mean 0.93 mg/dL (SD 0.25) in",
        "Table 1."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 4918,
    n_studies = 7,
    n_observations = 22843,
    age_range = "adults; pooled mean 60.53 years (SD 11.82), study means 57.23-66.64 years",
    weight_range = "pooled mean 82.48 kg (SD 16.87), study means 76.75-88.31 kg",
    sex_female_pct = 39.3,
    race_ethnicity = "not reported",
    disease_state = paste(
      "Four indications pooled: VTE prevention after elective hip or knee",
      "replacement (ODIXa-Hip2, ODIXa-OD-Hip, ODIXa-Knee), VTE treatment of",
      "acute symptomatic DVT (ODIXa-DVT, EINSTEIN DVT), acute coronary syndrome",
      "(ATLAS ACS TIMI-46) and nonvalvular atrial fibrillation (ROCKET AF)."
    ),
    dose_range = "2.5 mg once daily to 30 mg twice daily",
    renal_function = "pooled mean Tietz-truncated creatinine clearance 97.17 mL/min (SD 32.34); study means 81.76-108.4 mL/min",
    regions = "global phase II / III program (six phase II studies, one phase III substudy)",
    notes = paste(
      "4,918 of 5,041 enrolled patients contributed PK. Sparse sampling;",
      "median 3-8 observations per patient depending on study. Lower limit of",
      "quantification 0.5 ug/L, with below-LLOQ data excluded. Data from the",
      "immediate postsurgical phase of the VTE-prevention studies were excluded",
      "because patients there separate into slow and fast absorbers.",
      "Estimation by NONMEM 7.3, FOCE with interaction, on non-log-transformed",
      "data."
    )
  )

  ini({
    # Structural parameters -- Willmann 2018a Table 3. The typical values apply
    # at the covariate reference point (CRCL 93 mL/min, WT 81 kg, AGE 61 years,
    # male, no comedication, VTE-treatment indication) defined in
    # Supplementary Eq. 1.
    lka <- log(0.821)
    label("Typical first-order absorption rate constant (1/h)")
    lcl <- log(6.58)
    label("Typical apparent clearance CL/F (L/h)")
    lvc <- log(62.5)
    label("Typical apparent central volume V/F (L)")

    # Dose-dependent relative bioavailability. Willmann 2018a Supplementary
    # Eq. 1 gives F = Fmin + (Fmax - Fmin) * exp(-ln(2)/D50 * DOSE), an
    # exponential decay from Fmax at zero dose to the Fmin asymptote, halfway
    # down at DOSE = D50. Fmax is FIXED, and the anchor is that the function
    # returns 1 at the 10 mg reference dose (Table 3 footnote c): evaluating it
    # gives 0.590 + 0.660 * exp(-0.693147/14.4 * 10) = 0.998.
    fdepot_min <- 0.590
    label("Asymptotic minimum relative bioavailability at high dose (fraction)")
    fdepot_max <- fixed(1.25)
    label("Asymptotic maximum relative bioavailability as dose approaches zero (fraction)")
    led50 <- log(14.4)
    label("Dose at which relative bioavailability has fallen halfway from Fmax to Fmin, D50 (mg)")

    # Covariate effects on CL/F -- Willmann 2018a Table 3, applied in the power
    # / proportional forms of Supplementary Eq. 1.
    e_crcl_cl <- 0.406
    label("Power exponent on (Tietz-truncated CRCL / 93 mL/min) for CL/F (unitless)")
    e_wt_cl <- -0.278
    label("Power exponent on (WT / 81 kg) for CL/F (unitless)")
    e_pgp_inh_cl <- 0.966
    label("Multiplicative CL/F factor for concomitant P-glycoprotein inhibitor (fraction)")
    e_cyp3a4_inh_strong_cl <- 0.978
    label("Multiplicative CL/F factor for concomitant strong CYP3A4 inhibitor (fraction)")
    e_cyp3a4_inh_mod_cl <- 0.863
    label("Multiplicative CL/F factor for concomitant moderate CYP3A4 inhibitor (fraction)")
    e_cyp3a4_inh_weak_cl <- 0.939
    label("Multiplicative CL/F factor for concomitant weak CYP3A4 inhibitor (fraction)")
    e_cyp3a4_ind_cl <- 1.30
    label("Multiplicative CL/F factor for concomitant CYP3A4 inducer (fraction)")
    e_af_cl <- 0.849
    label("Multiplicative CL/F factor for the atrial fibrillation cohort vs VTE treatment (fraction)")
    e_acs_cl <- 1.14
    label("Multiplicative CL/F factor for the acute coronary syndrome cohort vs VTE treatment (fraction)")
    e_vte_p_le72_cl <- 1.04
    label("Multiplicative CL/F factor for VTE prevention at or before 72 h after first dose, vs VTE treatment (fraction)")
    e_vte_p_gt72_cl <- 1.29
    label("Multiplicative CL/F factor for VTE prevention after 72 h from first dose, vs VTE treatment (fraction)")

    # Covariate effects on V/F -- Willmann 2018a Table 3.
    e_wt_vc <- 0.216
    label("Power exponent on (WT / 81 kg) for V/F (unitless)")
    e_age_vc <- -0.189
    label("Power exponent on (AGE / 61 years) for V/F (unitless)")
    e_sexf_vc <- 0.889
    label("Multiplicative V/F factor for female sex vs male (fraction)")

    # Interindividual variability -- Willmann 2018a Table 3. All three etas are
    # exponential ('KA = TVKA*EXP(ETA(1))' and so on in the run12 $PK block).
    # ka is diagonal; CL/F and V/F are an estimated 2x2 block. The parenthesised
    # numbers beside the variances in Table 3 are SHRINKAGE percentages per
    # footnote d, not coefficients of variation.
    etalka ~ 0.628 # Table 3 'x2 ka' variance 0.628 (RSE 5.39%, bootstrap 95% CI 0.544-0.703; shrinkage 32.8%)
    etalcl + etalvc ~ c(
      0.167,
      0.0674, 0.0391
    ) # Table 3 'x2 CL/F' 0.167, 'x2 CL/F, V/F' covariance 0.0674, 'x2 V/F' 0.0391; implied correlation 0.834

    # Residual error -- run12 $ERROR is 'Y = IPRED + W*EPS(1)' with W = IPRED,
    # i.e. purely proportional on the untransformed scale (the paper notes a
    # proportional residual on non-log-transformed data fitted better than
    # log-transformed data).
    propSd <- 0.450555
    label("Proportional residual SD (fraction)")
  })

  model({
    # Tietz truncation of creatinine clearance (Willmann 2018a Methods and the
    # run12 $PK block): cap CRCL at 140 mL/min scaled to body surface area.
    crclCap <- 140 * BSA / 1.73
    crclTrunc <- CRCL * (CRCL <= crclCap) + crclCap * (CRCL > crclCap)

    # Comedication factor: the five indicators multiply, reproducing
    # COMM = COMM_PGPINH4 * COMM_SCYP3A4INH3 * COMM_MCYP3A4INH2 *
    #        COMM_WCYP3A4INH1 * COMM_CYP3A4IND5
    # in the run12 $PK block. Each term is 1 when its indicator is 0.
    comed <- e_pgp_inh_cl^CONMED_PGP_INH *
      e_cyp3a4_inh_strong_cl^CONMED_CYP3A4_INH_STRONG *
      e_cyp3a4_inh_mod_cl^CONMED_CYP3A4_INH_MOD *
      e_cyp3a4_inh_weak_cl^CONMED_CYP3A4_INH_WEAK *
      e_cyp3a4_ind_cl^CONMED_CYP3A4_IND

    # Study / indication factor (INDEFF in the run12 $PK block). VTE treatment
    # is the reference at 1; the VTE-prevention arm switches at 72 h after the
    # first dose, which the control stream keys on NONMEM TIME, so t here must
    # be referenced to the first dose. Written additively so that exactly one
    # active indicator selects its own factor and all-zero gives 1.
    vteFactor <- e_vte_p_le72_cl +
      (e_vte_p_gt72_cl - e_vte_p_le72_cl) * (t > 72)
    study <- 1 +
      DIS_AF * (e_af_cl - 1) +
      DIS_ACS * (e_acs_cl - 1) +
      DIS_VTE_P * (vteFactor - 1)

    ka <- exp(lka + etalka)
    ed50 <- exp(led50)
    cl <- exp(lcl + etalcl) *
      (crclTrunc / 93)^e_crcl_cl *
      (WT / 81)^e_wt_cl *
      comed *
      study
    vc <- exp(lvc + etalvc) *
      (WT / 81)^e_wt_vc *
      (AGE / 61)^e_age_vc *
      e_sexf_vc^SEXF

    kel <- cl / vc

    d / dt(depot) <- -ka * depot
    d / dt(central) <- ka * depot - kel * central

    # Relative bioavailability as a function of the administered dose
    # (Willmann 2018a Supplementary Eq. 1, F1 in the run12 $PK block).
    f(depot) <- fdepot_min +
      (fdepot_max - fdepot_min) *
        exp(-log(2) / ed50 * DOSE_RIVAROXABAN_MG)

    # Dose in mg over V/F in L gives mg/L; x1000 converts to the ug/L used
    # throughout the paper. The run12 $PK block folds the same factor into F1
    # ('; transform mg -> mcg') rather than into the observation.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
