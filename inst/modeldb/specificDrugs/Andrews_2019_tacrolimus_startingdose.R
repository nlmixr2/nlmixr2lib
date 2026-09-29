Andrews_2019_tacrolimus_startingdose <- function() {
  description <- "Two-compartment population PK model with first-order absorption and an absorption lag time for twice-daily oral immediate-release tacrolimus (Prograft) in adult kidney transplant recipients (Andrews 2019 STARTING-DOSE model). This is the companion reduced-covariate model to modellib('Andrews_2019_tacrolimus'), refit on the same 4527 concentrations from 337 adults but restricted to covariates that are known BEFORE transplantation and that do not change substantially afterwards. Albumin, serum creatinine and haematocrit are therefore dropped (Andrews 2019 states their last pre-transplant measurements did not significantly influence CL/F), as is the lean-body-mass effect on V1/F, leaving four covariate effects on apparent oral clearance CL/F: CYP3A5 expresser status (1.62 multiplier for *1/*1 or *1/*3 carriers vs the *3/*3 non-expresser reference), CYP3A4*22 carriage (0.814 multiplier vs CYP3A4*1 or unknown), age (power exponent -0.50 centred at 55.72 years) and body surface area (power exponent 0.72 centred at 1.93 m^2). Bioavailability was fixed to 1, so all disposition parameters are apparent. Inter-individual variability is diagonal on CL/F, V1/F, V2/F and Q/F; residual variability carries separate immunoassay and LC-MS/MS magnitudes selected per sample by the IMMUNOASSAY indicator. The paper pairs this model with a closed-form starting-dose algorithm, Dose = CL/F * AUC, targeting an AUC0-12h of 222 ng*h/mL (a pre-dose concentration of 10 ng/mL); see the validation vignette."
  reference <- "Andrews LM, Hesselink DA, van Schaik RHN, van Gelder T, de Fijter JW, Lloberas N, Elens L, Moes DJAR, de Winter BCM. A population pharmacokinetic model to predict the individual starting dose of tacrolimus in adult renal transplant recipients. Br J Clin Pharmacol. 2019;85(3):601-615. doi:10.1111/bcp.13838"
  vignette <- "Andrews_2019_tacrolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    AGE = list(
      description = "Recipient age at transplantation",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed per subject. Enters CL/F as the power term (AGE/55.72)^-0.50. The centring value 55.72 years is the model-building cohort median hardcoded in the Supporting Information Data S1 starting-dose control stream ('CLAGE = ((AGE/55.72)**THETA(9))'); the published Equation (2) prints it rounded to 56 years. The exponent is steeper than the final model's -0.43 because the starting-dose model must absorb, through the covariates it retains, part of the clearance variability that the dropped post-transplant covariates explained in the final model. Per-cohort medians (Table 1): Rotterdam 58.5 years (range 19.4-79.4), Leiden 54.0 years (range 15.0-77.0).",
      source_name = "AGE"
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters CL/F as the power term (BSA/1.93)^0.72; the centring value 1.93 m^2 appears both in the published Equation (2) and in the Supporting Information Data S1 starting-dose control stream ('CLBSA = ((BSA/1.93)**THETA(12))'). For the starting-dose use case BSA is evaluated at the pre-transplant measurement. Andrews 2019 Discussion argues BSA is a better indicator of metabolic mass than total bodyweight because it is less affected by abnormal adipose mass, and that basing the tacrolimus starting dose on bodyweight alone overexposes a considerable proportion of patients. Per-cohort medians (Table 1): Rotterdam 2.03 m^2 (range 1.24-2.66), Leiden 1.90 m^2 (range 1.33-2.48). The paper does not name the BSA formula used.",
      source_name = "BSA"
    ),
    CYP3A5_EXPR = list(
      description = "CYP3A5 expresser indicator: 1 if the patient carries at least one functional CYP3A5*1 allele (genotype *1/*1 or *1/*3), 0 if a non-expresser (*3/*3 or *3/*6).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CYP3A5 non-expresser: *3/*3 or *3/*6)",
      notes = "Time-fixed germline genotype (CYP3A5*3, rs776746), available before transplantation, which is why it is retained in the starting-dose model. Andrews 2019 Results states that graphical analysis showed no difference in the effect on CL/F between *1/*1 and *1/*3, so the two expresser genotypes are pooled; Equation (2) writes this as '1.62 if CYP3A5*1/*3 or CYP3A5*1/*1'. The Data S1 starting-dose control stream parameterises it as CLCYP3A5 = 1 + THETA(10), so the non-expresser reference multiplier is structurally exactly 1.0. Andrews 2019 Conclusions states the tacrolimus starting dose should be increased to 160 percent in CYP3A5*1 carriers. Genotype distribution is in Table 1; there was no deviation from Hardy-Weinberg equilibrium.",
      source_name = "CYP3A5"
    ),
    SNP_CYP3A4_RS35599367 = list(
      description = "CYP3A4*22 carrier indicator (rs35599367): 1 if the patient carries at least one *22 allele, 0 if CYP3A4*1 homozygous wild-type or if the genotype was not determined.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CYP3A4*1 wild-type, or genotype unknown -- Andrews 2019 Equation (2) pools these two groups)",
      notes = "Time-fixed germline genotype, available before transplantation. Equation (2) pools unknown-genotype subjects with the CYP3A4*1 wild-type reference, writing the reference stratum as '1.0 if CYP3A4*1 or unknown'; the Data S1 starting-dose control stream confirms this with an explicit missing-data branch ('IF(CYP3A4.EQ.-99) CLCYP3A4 = 1 ; Missing data') and parameterises the effect as CLCYP3A4 = 1 + THETA(11), so the wild-type reference multiplier is structurally exactly 1.0. Andrews 2019 Conclusions states the tacrolimus starting dose should be reduced to 80 percent in CYP3A4*22 carriers. Genotype distribution is in Table 1.",
      source_name = "CYP3A4"
    ),
    IMMUNOASSAY = list(
      description = "Per-sample bioanalytical assay indicator: 1 if the tacrolimus whole-blood concentration was measured by immunoassay (ACMIA or EMIT); 0 if measured by LC-MS/MS.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (LC-MS/MS reference method)",
      notes = "Per-sample (per-row) indicator, not a subject-level covariate. The starting-dose model was fit to the same pooled dataset as the final model and carries the same dual residual-error structure: the Data S1 starting-dose control stream branches on DVID exactly as the final-model stream does, giving DVID 1 (immunoassay) a combined proportional plus additive error and DVID 2 (LC-MS/MS) a proportional-only error. Rotterdam samples were assayed by ACMIA (LLOQ 1.5 ng/mL) or EMIT (LLOQ 2.0 ng/mL), pooled under a single immunoassay magnitude; Leiden samples by validated LC-MS/MS (LLOQ 1.0 ng/mL). In a pure-LC-MS/MS prospective dataset, set IMMUNOASSAY = 0 for every row.",
      source_name = "FLAG / DVID"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 337L,
    n_studies = 2L,
    age_range = "15.0-79.4 years (Rotterdam 19.4-79.4; Leiden 15.0-77.0)",
    age_median = "58.5 years (Rotterdam) and 54.0 years (Leiden); combined-cohort median 55.72 years per the Data S1 control-stream centring value",
    weight_range = "37.6-132.0 kg (Rotterdam) and 40.0-114.0 kg (Leiden)",
    weight_median = "79.4 kg (Rotterdam) and 74.0 kg (Leiden)",
    sex_female_pct = 39.5,
    race_ethnicity = c(Caucasian = 78.3, Asian = 9.2, `African descent` = 7.1, Other = 5.3),
    disease_state = "Adult renal transplant recipients during the first 3 months following kidney transplantation, all receiving oral twice-daily immediate-release tacrolimus (Prograft, Astellas) with mycophenolic acid and therapeutic drug monitoring.",
    dose_range = "The starting-dose algorithm of Equation (3) targets a pre-dose concentration of 10 ng/mL, corresponding to an AUC0-12h of 222 ng*h/mL, on a twice-daily schedule. Andrews 2019 also reports the AUC0-12h corresponding to other pre-dose targets: 12.5 ng/mL maps to 277 ng*h/mL and 15 ng/mL to 332 ng*h/mL. For a reference-covariate patient (CYP3A5 non-expresser, CYP3A4*1, age 55.72 years, BSA 1.93 m^2) the algorithm gives 222 * 22.5 / 1000 = 5.0 mg twice daily. The simulation trial compared this against a standard bodyweight-based 0.2 mg/kg/day dose.",
    regions = "The Netherlands (Erasmus MC, Rotterdam -- model building group 1, n = 237; Leiden University Medical Center -- model building group 2, n = 100).",
    n_observations = 4527L,
    model_purpose = "Andrews 2019 developed this reduced model specifically to select the tacrolimus dose at the moment of transplantation, when no post-transplant laboratory values yet exist. Each significant covariate from the final model was assessed for clinical relevance, feasibility of use, and whether it significantly influenced the starting dose. The authors report that the last measured albumin, serum creatinine and haematocrit BEFORE transplantation did not significantly influence CL/F, and that because these parameters change substantially after transplantation they were not incorporated. Time after transplantation was not a significant covariate, so the starting-dose model was fit to the same data as the final model rather than to a restricted early-post-transplant subset.",
    iov_structure = "Inter-occasion variability was estimated on CL/F (14.6 percent CV, Table 2 starting-dose model), with an occasion defined as the measurement of a pre-dose concentration. This model file does NOT encode IOV structurally, following the nlmixr2lib convention for source models with no operational occasion column. A downstream user who wants to simulate IOV can add an occasion indicator and a per-occasion eta on CL/F with variance log(1 + 0.146^2) = 0.021158.",
    external_validation = "Validated internally and externally with prediction-corrected VPCs on the same cohorts used for the final model (Andrews 2019 Figure 4). The authors report that in the external validation of the starting-dose model both the median and the variability were adequately described -- a better outcome than for the final model, whose external validation was limited by the absence of albumin measurements in that cohort.",
    simulation_trial = "A simulated clinical trial using the model-building cohort's own patient characteristics compared the standard bodyweight-based dose against the starting-dose algorithm, simulating C0 and AUC 1000 times per patient at day 10 post-transplant. Model-based dosing put 33.0 percent of patients on target (10-15 ng/mL) versus 26.1 percent for bodyweight-based dosing, reduced the fraction above target from 44.5 to 36.8 percent, reduced markedly supratherapeutic exposure (>20 ng/mL) from 24.6 to 15.6 percent, and reduced markedly subtherapeutic exposure (<5 ng/mL) from 7.2 to 5.2 percent. Median C0 was 13.9 ng/mL (bodyweight-based) versus 12.9 ng/mL (model-based); median AUC 298.5 versus 277.9 ng*h/mL.",
    notes = "Pooled analysis of two Dutch cohorts (Erasmus MC Rotterdam RCT, n = 237; LUMC Leiden routine care, n = 100). Baseline demographics are in Table 1, which reports each cohort separately; the combined-cohort covariate medians used as centring values are recovered from the Supporting Information Data S1 control stream. Parameter estimates are Table 2 column 'Starting dose model (RSE %) [shrinkage]'. The model is intended for IMMEDIATE-RELEASE twice-daily oral tacrolimus in ADULTS at the start of therapy; it is not validated for once-daily extended-release formulations or for children. Andrews 2019 notes the next step is a prospective pilot study of the starting-dose algorithm, for which ethics approval had been obtained."
  )

  ini({
    # Starting-dose-model fixed-effect estimates from Andrews 2019 Table 2,
    # column 'Starting dose model (RSE %) [shrinkage]'. Bioavailability was
    # fixed to 1 (Methods, Base model development), so every clearance and
    # volume below is an APPARENT value. The reference subject for the typical
    # values is a CYP3A5 non-expresser (*3/*3), CYP3A4*1-or-unknown patient of
    # age 55.72 years with a BSA of 1.93 m^2.
    ltlag <- log(0.39); label("Absorption lag time tlag (h)") # Andrews 2019 Table 2 starting dose model tlag = 0.39 h (RSE 12%)
    lka <- log(3.70); label("Absorption rate constant ka (1/h)") # Andrews 2019 Table 2 starting dose model ka = 3.70 (RSE 13%); the Table 2 row label reads 'l h-1' which is a typo -- a first-order absorption rate constant has units 1/h
    lcl <- log(22.5); label("Apparent oral clearance CL/F at the reference subject (L/h)") # Andrews 2019 Table 2 starting dose model CL/F = 22.5 L/h (RSE 3%); also the leading coefficient of Equations (2) and (3)
    lvc <- log(685); label("Apparent central volume V1/F (L)") # Andrews 2019 Table 2 starting dose model V1/F = 685 L (RSE 5%)
    lq <- log(10.6); label("Apparent inter-compartmental clearance Q/F (L/h)") # Andrews 2019 Table 2 starting dose model Q/F = 10.6 L/h (RSE 6%)
    lvp <- log(6590); label("Apparent peripheral volume V2/F (L)") # Andrews 2019 Table 2 starting dose model V2/F = 6590 L (RSE 14%)

    # Covariate effects on CL/F. Andrews 2019 Equation (2):
    #   CL/F = 22.5 * [(1.0 if CYP3A5*3/*3) or (1.62 if CYP3A5*1/*3 or *1/*1)]
    #               * [(1.0 if CYP3A4*1 or unknown) or (0.814 if CYP3A4*22)]
    #               * (Age/56)^-0.50 * (BSA/1.93)^0.72
    # The two exponent signs are stated in Table 2 and confirmed by the
    # Conclusions statement that the starting dose should be higher in younger
    # patients and in those with a higher BSA. Equation (2) prints the age
    # centring value rounded to 56; this file uses the unrounded control-stream
    # constant 55.72 (Data S1), the value the model was actually fitted with.
    # Where Table 2 and Equation (2) differ in precision for the genotype
    # multipliers (Table 2 gives 0.81, Equation (2) gives 0.814) this file
    # takes the equation's value.
    e_cyp3a5_expr_cl <- 1.62; label("CYP3A5 expresser (*1/*1 or *1/*3) multiplier on CL/F") # Andrews 2019 Equation (2) = 1.62; Table 2 starting dose model CYP3A5*1 = 1.62 (RSE 14%)
    e_cyp3a4_22_cl <- 0.814; label("CYP3A4*22 carrier multiplier on CL/F") # Andrews 2019 Equation (2) = 0.814; Table 2 starting dose model CYP3A4*22 = 0.81 (RSE 36%)
    e_age_cl <- -0.50; label("Age power exponent on CL/F, centred at 55.72 years (unitless)") # Andrews 2019 Table 2 starting dose model Age = -0.50 (RSE 15%); sign confirmed by the Conclusions statement that the starting dose should be higher in younger patients
    e_bsa_cl <- 0.72; label("Body surface area power exponent on CL/F, centred at 1.93 m^2 (unitless)") # Andrews 2019 Table 2 starting dose model BSA = 0.72 (RSE 29%); sign confirmed by the Conclusions statement that the starting dose should be higher in patients with a higher BSA

    # Diagonal inter-individual variability. Andrews 2019 Table 2 reports IIV
    # as %CV with no inter-eta correlations (the Data S1 starting-dose control
    # stream uses a diagonal $OMEGA for the four disposition etas). Variances
    # on the internal log scale are omega^2 = log(1 + CV^2):
    #   CL/F CV 39.4% -> log(1 + 0.394^2) = 0.144305
    #   V1/F CV 54.0% -> log(1 + 0.540^2) = 0.255882
    #   V2/F CV 53.7% -> log(1 + 0.537^2) = 0.253377
    #   Q/F  CV 79.6% -> log(1 + 0.796^2) = 0.490796
    etalcl ~ 0.144305 # Andrews 2019 Table 2 starting dose model IIV CL/F 39.4% CV [10% shrinkage]
    etalvc ~ 0.255882 # Andrews 2019 Table 2 starting dose model IIV V1/F 54.0% CV [19% shrinkage]
    etalvp ~ 0.253377 # Andrews 2019 Table 2 starting dose model IIV V2/F 53.7% CV [40% shrinkage]
    etalq ~ 0.490796 # Andrews 2019 Table 2 starting dose model IIV Q/F 79.6% CV [29% shrinkage]

    # Residual error, switched per sample by the IMMUNOASSAY indicator, with
    # the same DVID-branched structure as the final model: immunoassay rows get
    # a combined proportional plus additive error, LC-MS/MS rows a
    # proportional-only error, so the LC-MS/MS additive SD is fixed at 0.
    propSd_immuno <- 0.169; label("Proportional residual SD for immunoassay samples (fraction)") # Andrews 2019 Table 2 starting dose model Proportional Immunoassay = 16.9% (RSE 6%) [25% shrinkage]
    addSd_immuno <- 1.02; label("Additive residual SD for immunoassay samples (ng/mL)") # Andrews 2019 Table 2 starting dose model Additive Immunoassay = 1.02 ug/L (RSE 10%) [25% shrinkage]; 1 ug/L = 1 ng/mL
    propSd_lcms <- 0.244; label("Proportional residual SD for LC-MS/MS samples (fraction)") # Andrews 2019 Table 2 starting dose model Proportional LC-MS/MS = 24.4% (RSE 5%) [12% shrinkage]
    addSd_lcms <- fixed(0); label("Additive residual SD for LC-MS/MS samples (ng/mL)") # Structurally absent: the Data S1 starting-dose control stream gives DVID 2 (LC-MS/MS) a proportional-only error, and Table 2 reports an additive term only for the immunoassay arm
  })

  model({
    # CYP3A5 expresser multiplier on CL/F, parameterised in the Data S1
    # starting-dose control stream as CLCYP3A5 = 1 + THETA(10), so the
    # non-expresser reference multiplier is structurally exactly 1.0.
    f_cyp3a5 <- 1 + (e_cyp3a5_expr_cl - 1) * CYP3A5_EXPR

    # CYP3A4*22 carrier multiplier on CL/F, likewise CLCYP3A4 = 1 + THETA(11).
    # Unknown-genotype subjects are pooled with the wild-type reference.
    f_cyp3a4 <- 1 + (e_cyp3a4_22_cl - 1) * SNP_CYP3A4_RS35599367

    # The two pre-transplant continuous covariates, centred at the
    # combined-cohort medians from the Data S1 control stream.
    f_age <- (AGE / 55.72)^e_age_cl
    f_bsa <- (BSA / 1.93)^e_bsa_cl

    # Individual PK parameters. Unlike the final model, V1/F carries no
    # covariate here -- the Data S1 starting-dose control stream sets
    # 'V2 = TV2 * EXP(ETA(2))' with no covariate factor.
    tlag <- exp(ltlag)
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * f_cyp3a5 * f_cyp3a4 * f_age * f_bsa
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    # Two-compartment oral disposition (NONMEM ADVAN4 TRANS4). The dose lands
    # in `depot`; bioavailability was fixed to 1 and is implicit in the
    # apparent CL/F and V/F parameterisation.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    alag(depot) <- tlag

    # Whole-blood tacrolimus in ng/mL. Dose in mg and vc in L give mg/L =
    # ug/mL, so multiply by 1000; this reproduces the control stream's
    # 'S2 = V2 / 1000'.
    Cc <- central / vc * 1000

    # Per-sample assay-conditional residual error (IMMUNOASSAY: 1 =
    # immunoassay, 0 = LC-MS/MS reference). addSd_lcms is fixed at 0, so
    # LC-MS/MS rows reduce to the proportional-only error the paper specifies.
    addSd <- addSd_immuno * IMMUNOASSAY + addSd_lcms * (1 - IMMUNOASSAY)
    propSd <- propSd_immuno * IMMUNOASSAY + propSd_lcms * (1 - IMMUNOASSAY)
    Cc ~ add(addSd) + prop(propSd)
  })
}
