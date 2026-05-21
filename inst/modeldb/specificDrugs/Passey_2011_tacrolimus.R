Passey_2011_tacrolimus <- function() {
  description <- paste0(
    "Steady-state apparent-clearance regression model for oral tacrolimus ",
    "trough concentrations in adult kidney-transplant recipients (Passey 2011). ",
    "Encoded as a 1-compartment IV continuous-infusion model with a nominal ",
    "fixed central volume of distribution: at steady state, Cc = dose-rate / ",
    "CL/F is independent of V, so the rxode2 simulation reproduces the paper's ",
    "regression-style trough prediction. Apparent clearance CL/F is multiplied ",
    "by five covariate factors: an ordered-categorical days-post-transplant ",
    "effect (3-5 = reference, 6-10 = 0.86, 11-180 = 0.71), three-level ",
    "CYP3A5 genotype (CYP3A5*3/*3 = reference, CYP3A5*1/*3 = 1.70, ",
    "CYP3A5*1/*1 = 2.00), steroid-sparing immunosuppression protocol (0.70), ",
    "a power-form age effect ((Age/50)^-0.40), and concomitant calcium channel ",
    "blocker coadministration (0.94)."
  )
  reference <- paste0(
    "Passey C, Birnbaum AK, Brundage RC, Oetting WS, Israni AK, Jacobson PA. ",
    "Dosing equation for tacrolimus using genetic variants and clinical ",
    "factors. Br J Clin Pharmacol. 2011;72(6):948-957. ",
    "doi:10.1111/j.1365-2125.2011.04039.x."
  )
  vignette <- "Passey_2011_tacrolimus"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    POD = list(
      description        = "Post-operative day (days post-transplant)",
      units              = "days",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Time-varying. Days elapsed since kidney transplantation. Passey 2011 ",
        "Methods / Population modelling: 'days post transplant were converted ",
        "to an ordered categorical covariate and classified as: immediate ",
        "post transplant (days 3-5), early post transplant (days 6-10) and ",
        "late post transplant (days 11-180). Categorization was done as the ",
        "model failed to converge when days post transplant was modelled as ",
        "a continuous function.' Inside model() the continuous POD column ",
        "is decomposed into two binary indicators ((POD>=6 & POD<=10) and ",
        "(POD>=11)) so the typical-value factor reproduces the paper's ",
        "categorical effect. The Passey 2011 model was fitted only to ",
        "concentrations at POD >= 3 (the paper drops POD < 3 because troughs ",
        "before then are not yet near steady state) and POD <= 180 (the ",
        "first 6 months post-transplant). Mean 17 trough concentrations per ",
        "patient (range 1-24); Table 1 cohort summary."
      ),
      source_name        = "days post transplant"
    ),
    AGE = list(
      description        = "Recipient age at the time of transplantation",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Time-fixed (baseline) per subject. Passey 2011 Table 1: ",
        "mean +/- SD 50.2 +/- 12.2 years across all 681 subjects (non-AA ",
        "50.1 +/- 12.2; AA 46.9 +/- 11.5). Centred at the cohort median ",
        "(50 years) for the power effect (Age/50)^-0.40 on CL/F (Passey ",
        "2011 Results: 'CL/F increased until the median age of 50 years and ",
        "then decreased thereafter')."
      ),
      source_name        = "age"
    ),
    CYP3A5_STAR1_HET = list(
      description        = "CYP3A5*1/*3 heterozygote indicator (one functional CYP3A5*1 allele)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP3A5*3/*3 nonexpresser, when paired with CYP3A5_STAR1_HOM = 0)",
      notes              = paste0(
        "Time-fixed per subject (germline genotype). 1 = subject is ",
        "CYP3A5*1/*3 heterozygote (one functional *1 allele at rs776746); ",
        "0 = otherwise (the union of *3/*3 nonexpressers and *1/*1 ",
        "homozygotes; the paired indicator CYP3A5_STAR1_HOM flags the ",
        "homozygous group). Passey 2011 cohort: 129/681 (19%) were *1/*3 ",
        "heterozygotes (non-AA 13%, AA 42%). Reference category for the ",
        "paired indicators (CYP3A5_STAR1_HET = 0 AND CYP3A5_STAR1_HOM = 0) ",
        "is the *3/*3 nonexpresser group (476/681 subjects, 70%)."
      ),
      source_name        = "CYP3A5*1/*3"
    ),
    CYP3A5_STAR1_HOM = list(
      description        = "CYP3A5*1/*1 homozygote indicator (two functional CYP3A5*1 alleles)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP3A5*3/*3 nonexpresser, when paired with CYP3A5_STAR1_HET = 0)",
      notes              = paste0(
        "Time-fixed per subject (germline genotype). 1 = subject is ",
        "CYP3A5*1/*1 homozygote (two functional *1 alleles at rs776746); ",
        "0 = otherwise (the union of *3/*3 nonexpressers and *1/*3 ",
        "heterozygotes; the paired indicator CYP3A5_STAR1_HET flags the ",
        "heterozygous group). Passey 2011 cohort: 72/681 (11%) were *1/*1 ",
        "homozygotes (non-AA 2%, AA 45%); large sample size in the AA ",
        "subgroup motivated the three-level decomposition (Passey 2011 ",
        "Discussion: 'these data demonstrate for the first time that the ",
        "three CYP3A5 genotypes have distinctive CL/F estimates')."
      ),
      source_name        = "CYP3A5*1/*1"
    ),
    CONMED_STEROID_SPARING = list(
      description        = "Steroid-sparing immunosuppressive protocol indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-steroid-sparing centre; continuous corticosteroid use)",
      notes              = paste0(
        "Time-fixed per subject (assigned by transplantation centre). 1 = ",
        "subject was transplanted at a centre using a steroid-sparing ",
        "immunosuppressive protocol (corticosteroids administered for <= 7 ",
        "days post-transplant, per Passey 2011 Methods: 'Centres were ",
        "designated as using a steroid sparing immunosuppressive regimen if ",
        "they administered steroids for <= 7 days post transplant'); 0 = ",
        "subject was at a non-sparing centre (continuous steroid therapy). ",
        "Passey 2011 cohort: 205/681 (30%) at steroid-sparing centres ",
        "(non-AA 36.5%, AA 5.5%). The 0.70 multiplicative factor on CL/F is ",
        "hypothesized to reflect reduced CYP3A induction in the absence of ",
        "ongoing corticosteroid therapy."
      ),
      source_name        = "steroid sparing centre"
    ),
    CONMED_CCB = list(
      description        = "Concomitant calcium channel blocker indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant CCB)",
      notes              = paste0(
        "Time-varying. 1 = subject coadministered any calcium channel ",
        "blocker (CCB) at the time of the trough measurement; 0 = no CCB ",
        "coadministered. Passey 2011 cohort: 5082/11823 (43%) of trough ",
        "samples carried CCB = 1 (non-AA 41%, AA 51%). The CCB class was ",
        "pooled across specific agents (e.g., diltiazem, amlodipine); ",
        "diltiazem is a potent CYP3A inhibitor and amlodipine is not, so ",
        "the 0.94 multiplicative factor on CL/F is a pooled-class estimate ",
        "that under-represents diltiazem-specific inhibition and ",
        "over-represents amlodipine-specific inhibition (Passey 2011 ",
        "Discussion)."
      ),
      source_name        = "calcium channel blocker use"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 681L,
    n_studies        = 1L,
    n_concentrations = 11823L,
    age_range        = ">= 18 years (Table 1 mean +/- SD: 50.2 +/- 12.2 years; non-AA 50.1 +/- 12.2; AA 46.9 +/- 11.5)",
    age_median       = "50 years (used as the reference age in the (Age/50)^-0.40 covariate effect)",
    weight_range     = "Table 1 mean +/- SD 81.3 +/- 18.7 kg (non-AA 81.1 +/- 18.8; AA 81.9 +/- 17.9). Weight was tested as a covariate on CL/F but was NOT significant and is NOT in the final model.",
    sex_female_pct   = 37.0,
    race_ethnicity   = c(`non-African American` = 79.3, `African American` = 20.7),
    disease_state    = paste0(
      "Adult kidney transplant recipients enrolled in the multicentre ",
      "Deterioration of Kidney Allograft Function (DeKAF) Genomics study ",
      "(NCT00270712) between 2006 and 2008. All recipients had end-stage ",
      "renal dysfunction undergoing kidney or kidney-pancreas transplant. ",
      "Living donor 59% (non-AA 66%, AA 31%); first transplant 81% (non-AA ",
      "77%, AA 94%). 30% transplanted at a steroid-sparing centre."
    ),
    dose_range       = paste0(
      "Oral tacrolimus (Prograf), administered once or twice daily during ",
      "the first 6 months post-transplant. Initial dose based on body weight; ",
      "subsequent doses titrated to target trough concentrations (institution-",
      "specific targets, generally 8-12 ng/mL in months 1-3 and 6-10 ng/mL ",
      "in months 3-6 post-transplant). Mean +/- SD daily dose 0.08 +/- 0.05 ",
      "mg/kg/day (non-AA 0.07 +/- 0.05; AA 0.09 +/- 0.04)."
    ),
    regions          = "United States and Canada (multicentre observational trial)",
    sampling_design  = paste0(
      "Sparse trough sampling during routine clinical care: twice weekly in ",
      "weeks 1-8 post-transplant and twice in each of months 3, 4, 5, and 6 ",
      "post-transplant. Mean 17 troughs per patient (range 1-24). Only ",
      "troughs measured after day 2 post-transplant (POD >= 3) were ",
      "retained for modelling, to ensure tacrolimus was at or near steady ",
      "state. 97.1% of trough assays were liquid chromatography-mass ",
      "spectrometry; all assays were CLIA-certified or CLIA-quality."
    ),
    cyp3a5_genotype  = c(`*3/*3` = 70.0, `*1/*3` = 19.0, `*1/*1` = 11.0),
    notes            = paste0(
      "Software: NONMEM v7.1 (FOCE-I); R 2.4.1 for diagnostics; PdxPop ",
      "v4.0. Baseline demographics and CYP3A5 genotype counts per Passey ",
      "2011 Tables 1 and 2; final-model parameter estimates per Passey ",
      "2011 Table 3. The model is a steady-state regression: tacrolimus ",
      "trough Cmin was approximated by the average steady-state ",
      "concentration Css,avg (justified by the long elimination half-life ",
      "of ~12 h relative to the q12h or qd dosing interval), and Cobs = ",
      "(dose-rate) / (CL/F) was the structural prediction. No absorption ",
      "rate constant ka, no volume of distribution V, and no peripheral ",
      "compartments are estimated by Passey 2011. Race (African American ",
      "vs non-AA), recipient weight, gender, donor type, donor gender, ",
      "pre-emptive transplant, number of prior transplants, and the SNPs ",
      "rs12114000 (CYP3A4), rs3734354 (SIM1), rs4926 (SERPING1), rs3135506 ",
      "(APOA5), and rs2608555 (GAN) were tested and were not significant."
    )
  )

  ini({
    # ----- Structural parameter (Passey 2011 Table 3 final estimate) -----
    # Typical apparent oral clearance CL/F for the reference subject:
    # CYP3A5*3/*3 nonexpresser, transplanted at a non-steroid-sparing
    # (i.e., continuous steroid) centre, age 50 years, NOT taking a CCB,
    # POD = 3-5 days. Reported as 38.4 L/h with RSE 4.14% (95% CI 35.3,
    # 41.5; bootstrap median 38.3).
    lcl <- log(38.4); label("Apparent oral clearance CL/F at the reference subject (L/h)")  # Passey 2011 Table 3, typical value of CL/F

    # The Passey 2011 model is a steady-state regression: trough = (dose-
    # rate) / CL/F. The paper does not estimate a volume of distribution
    # ("the concentrations were analysed according to a steady-state
    # infusion model... the trough concentrations were well approximated
    # by the average steady-state plasma concentration"). To run the model
    # as a rxode2 1-compartment IV continuous-infusion solver, a nominal
    # central volume is required for the ODE to be well-defined; vc is
    # fixed at 1000 L purely for that purpose. At steady state the
    # predicted concentration Cc = R/CL is independent of V, so the choice
    # of the nominal value does not affect the model's predictions; the
    # value 1000 L was chosen because it places kel = CL/V near the
    # tacrolimus reported half-life (~12 h at CL/F = 38.4 L/h, V/F = 1000
    # L gives t_{1/2} = ln(2)*1000/38.4 = 18 h, close enough for stable
    # numerical integration over typical observation windows). See the
    # vignette's Assumptions and deviations section for the full
    # rationale.
    lvc <- fixed(log(1000)); label("Nominal central volume of distribution (L) - not estimated by Passey 2011; placeholder for the rxode2 ODE")

    # ----- Categorical days-post-transplant effects (Passey 2011 Table 3) -----
    # Multiplicative factors on CL/F relative to the days 3-5 reference
    # (implicit factor of 1.0). Days 6-10: 0.86 (RSE 4.13%, 95% CI 0.80-
    # 0.93). Days 11-180: 0.71 (RSE 4.14%, 95% CI 0.66-0.77). Bootstrap
    # medians are 0.86 and 0.72, respectively.
    e_dpt_6_10_cl   <- 0.86; label("CL/F multiplicative factor for days 6-10 post-transplant (unitless)")     # Passey 2011 Table 3, 6-10 DPT
    e_dpt_11_180_cl <- 0.71; label("CL/F multiplicative factor for days 11-180 post-transplant (unitless)")   # Passey 2011 Table 3, 11-180 DPT

    # ----- Three-level CYP3A5 genotype effects (Passey 2011 Table 3) -----
    # Multiplicative factors on CL/F relative to the *3/*3 nonexpresser
    # reference (implicit factor of 1.0). Encoded with two paired binary
    # indicators (CYP3A5_STAR1_HET, CYP3A5_STAR1_HOM) following the
    # SLCO1B1_HAP15_HET / SLCO1B1_HAP15_HOM precedent in
    # inst/references/covariate-columns.md. Table 3 reports the original
    # NONMEM point estimate 1.70 with RSE 3.99% for *1/*3 (bootstrap
    # median 1.69); the published final-model equation in the abstract /
    # Results uses the bootstrap-median 1.69. The original-estimate 1.70
    # is used here to match the "Study population Estimate" column and
    # the same convention is applied to the other Table 3 factors (e.g.,
    # 0.71 for 11-180 DPT, bootstrap median 0.72). The 1.69 vs 1.70
    # rounding discrepancy is noted in the vignette Assumptions and
    # deviations section.
    e_cyp3a5_het_cl <- 1.70; label("CL/F multiplicative factor for CYP3A5*1/*3 heterozygote (unitless)")  # Passey 2011 Table 3, CYP3A5*1/*3
    e_cyp3a5_hom_cl <- 2.00; label("CL/F multiplicative factor for CYP3A5*1/*1 homozygote (unitless)")    # Passey 2011 Table 3, CYP3A5*1/*1

    # ----- Steroid-sparing centre effect (Passey 2011 Table 3) -----
    # Multiplicative factor on CL/F when the patient was transplanted at
    # a steroid-sparing centre. 0.70 (RSE 3.50%, 95% CI 0.65-0.75;
    # bootstrap median 0.70). Reference: continuous-steroid (non-sparing)
    # centre, implicit factor of 1.0.
    e_steroid_spare_cl <- 0.70; label("CL/F multiplicative factor for steroid-sparing immunosuppression protocol (unitless)")  # Passey 2011 Table 3, Steroid sparing centre

    # ----- Power-form age effect (Passey 2011 Table 3) -----
    # Exponent on (Age / 50) for the typical-value covariate equation:
    # CL/F *= (Age / 50)^e_age_cl. -0.40 with RSE 13.5% (95% CI -0.50,
    # -0.30; bootstrap median -0.39). The exponent is negative so CL/F
    # decreases as age moves above 50 years; Passey 2011 Results notes
    # that CL/F also increases with age up to 50 years.
    e_age_cl <- -0.40; label("Power exponent of (Age / 50) on CL/F (unitless)")  # Passey 2011 Table 3, (Age/50) exponent

    # ----- Concomitant CCB effect (Passey 2011 Table 3) -----
    # Multiplicative factor on CL/F when the patient was coadministered a
    # calcium channel blocker at the time of the trough measurement.
    # 0.94 (RSE 2.43%, 95% CI 0.89-0.98; bootstrap median 0.94).
    e_ccb_cl <- 0.94; label("CL/F multiplicative factor for concomitant calcium channel blocker (unitless)")  # Passey 2011 Table 3, CCB use

    # ----- Inter-individual variability (Passey 2011 Table 3) -----
    # IIV on CL/F is exponential (log-normal): CL/F_i = TVCL * exp(eta_i)
    # with var(eta) = omega^2. Table 3 reports CV% = 40.1 (95% CI
    # 37.4-43.6); omega^2 = log(1 + CV^2) = log(1 + 0.401^2) =
    # log(1.160801) = 0.149. The base model (no covariates) had IIV =
    # 52.3% CV; the final model explains 23.2% of the relative IIV.
    etalcl ~ 0.149   # Passey 2011 Table 3, IIV for CL/F = 40.1% CV; var = log(1 + 0.401^2) = 0.149

    # ----- Residual error (Passey 2011 Table 3) -----
    # Additive residual error model on the ng/mL concentration scale:
    # Cobs = Cpred + epsilon, epsilon ~ N(0, sigma^2). Table 3 reports
    # SD = 3.19 (95% CI 3.07-3.32; bootstrap median 3.19).
    addSd <- 3.19; label("Additive residual error SD (ng/mL)")  # Passey 2011 Table 3, RUV additive SD
  })

  model({
    # ----- 1. Derived categorical covariate indicators from POD -----
    # The paper categorises POD into three ordered groups (3-5, 6-10,
    # 11-180); the model receives the continuous canonical POD column
    # and derives the two binary indicators here so the typical-value
    # equation reproduces the paper's categorical effect exactly. The
    # day 3-5 reference category corresponds to both indicators = 0.
    dpt_6_10   <- (POD >= 6) * (POD <= 10)
    dpt_11_180 <- (POD >= 11)

    # ----- 2. Multiplicative covariate factors on CL/F -----
    # Days-post-transplant categorical factor (power-of-binary form;
    # the indicators are mutually exclusive so at most one exponent is 1).
    dpt_factor <- e_dpt_6_10_cl^dpt_6_10 * e_dpt_11_180_cl^dpt_11_180

    # CYP3A5 genotype factor (power-of-binary form; the paired
    # indicators are mutually exclusive so at most one exponent is 1).
    # *3/*3 reference: both indicators are 0 so the factor collapses to 1.
    cyp3a5_factor <- e_cyp3a5_het_cl^CYP3A5_STAR1_HET *
                     e_cyp3a5_hom_cl^CYP3A5_STAR1_HOM

    # Steroid-sparing protocol factor (power-of-binary form).
    steroid_factor <- e_steroid_spare_cl^CONMED_STEROID_SPARING

    # Concomitant CCB factor (power-of-binary form).
    ccb_factor <- e_ccb_cl^CONMED_CCB

    # Power-form age effect, centred at the cohort median age of 50 years.
    age_factor <- (AGE / 50)^e_age_cl

    # ----- 3. Individual CL/F -----
    # Exponential (log-normal) IIV on CL/F; multiplicative covariate
    # factors per the Passey 2011 final-model equation. The expanded
    # equation matches the paper's typed-out form:
    #   CL/F = 38.4 * [(0.86 if days 6-10) or (0.71 if days 11-180)] *
    #              [(1.70 if *1/*3) or (2.00 if *1/*1)] *
    #              (0.70 if steroid sparing centre) *
    #              (Age / 50)^-0.40 *
    #              (0.94 if CCB) *
    #              exp(eta_CL)
    cl <- exp(lcl + etalcl) * dpt_factor * cyp3a5_factor *
          steroid_factor * age_factor * ccb_factor

    # Nominal central volume (placeholder for the rxode2 ODE; cancels at
    # steady state).
    vc <- exp(lvc)

    # ----- 4. ODE system (1-compartment IV continuous-infusion) -----
    # Steady-state IV continuous-infusion encoding of the Passey 2011
    # regression model. Dosing events specify a continuous infusion at
    # rate R (mg/h) with SS = 1 so the central compartment reaches Css
    # = R / CL (independent of V). Simulating non-steady-state scenarios
    # with this model file is not supported by Passey 2011 - see the
    # vignette's Assumptions and deviations section.
    kel <- cl / vc
    d/dt(central) <- -kel * central

    # ----- 5. Observation and error -----
    # Convert mg/L (central / vc) to ng/mL (the units reported by Passey
    # 2011): 1 mg/L = 1 ug/mL = 1000 ng/mL. The factor of 1000 below
    # carries the unit conversion only; V cancels at steady state.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd)
  })
}
