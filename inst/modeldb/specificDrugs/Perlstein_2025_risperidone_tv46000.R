# Sequential parent-metabolite population PK model for TV-46000, a
# copolymer-based long-acting subcutaneous risperidone suspension (LASCA), in
# 692 healthy volunteers and adults with schizophrenia or schizoaffective
# disorder pooled from three phase 1 and two phase 3 studies
# (Perlstein 2025, Neurol Ther 14(3):829-848; doi:10.1007/s40120-025-00723-z).

Perlstein_2025_risperidone_tv46000 <- function() {
  description <- paste(
    "Sequential parent-metabolite population PK model for TV-46000, a",
    "long-acting subcutaneous antipsychotic (LASCA) that combines risperidone",
    "with a copolymer-based delivery technology in a suspension given",
    "subcutaneously once monthly (q1m) or once every 2 months (q2m), in 692",
    "healthy volunteers and adults with schizophrenia or schizoaffective",
    "disorder pooled from three phase 1 and two phase 3 studies",
    "(Perlstein 2025). The single subcutaneous depot empties by two parallel",
    "first-order release routes: a fast direct route with rate constant ka1",
    "straight into the risperidone central compartment (the initial release",
    "from the subcutaneous depot, which brings the total active moiety above",
    "10 ng/mL within 24 h), and a slow indirect route with rate constant ka2",
    "into a five-compartment transit chain with rate constant ktr that",
    "delivers the remainder into central over the 28-day or 56-day dosing",
    "interval (the in situ depot). Risperidone is described by a",
    "one-compartment model with first-order elimination CL/F, and because the",
    "fraction metabolized FRMET is fixed to 1 the whole of that elimination",
    "flux forms the equipotent metabolite 9-hydroxyrisperidone (paliperidone),",
    "which is itself described by a one-compartment model with first-order",
    "elimination CLMO. The two analytes are summed into the total active",
    "moiety in risperidone-equivalent units, TAM = [risperidone] +",
    "[9-OH-risperidone] * 410/426, which is the definition the paper gives in",
    "its Introduction and the same molecular-weight correction the upstream",
    "Ivaturi 2017 model applies. Body mass index shifts the balance between the two",
    "release routes (lower BMI gives a faster direct route and a slower",
    "indirect route), a larger injection volume slows the direct route, and",
    "upper-arm rather than abdominal injection raises ka1 by 33 percent.",
    "Interindividual variability is exponential on ka1, ka2 (correlated with",
    "ka1), CL and CLMO, and the log-transformed-both-sides residual error is",
    "proportional on each analyte. The total active moiety TAM = risperidone +",
    "9-hydroxyrisperidone is the clinically reported quantity, and dopamine",
    "D2-receptor occupancy D2RO is derived from it with the literature Emax",
    "model the authors applied in their simulations (Emax 100 percent, kd 10.1",
    "ng/mL). The structural starting point was the RBP-7000 dual-absorption",
    "model of Ivaturi 2017; see modellib('Ivaturi_2017_RBP_7000').",
    sep = " "
  )
  reference <- paste(
    "Perlstein I, Merenlender Wagner A, Elgart A, Zandvliet AS, Hellmann F,",
    "Lin Y, van Maanen E, Plock N, Fauchet F, Singh R (2025).",
    "Population pharmacokinetic modeling of TV-46000, a risperidone",
    "long-acting subcutaneous antipsychotic for the treatment of patients",
    "with schizophrenia. Neurol Ther 14(3):829-848.",
    "doi:10.1007/s40120-025-00723-z.",
    "Structural starting point from Ivaturi V, Gopalakrishnan M, Gobburu JVS,",
    "et al. (2017) Br J Clin Pharmacol 83(7):1476-1498; doi:10.1111/bcp.13246;",
    "see modellib('Ivaturi_2017_RBP_7000').",
    sep = " "
  )
  vignette <- "Perlstein_2025_risperidone_tv46000"
  units    <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    BMI = list(
      description        = "Body mass index",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters as a median-normalized power term on BOTH release rate",
        "constants, with opposite signs: ka1 * (BMI / 28.7)^(-1.1) on the fast",
        "direct route and ka2 * (BMI / 28.7)^1.7 on the slow indirect route",
        "(Perlstein 2025 Table 3, KA1BMI1 and KA2BMI1; continuous-covariate",
        "form Pi = P * (COVi / COVmedian)^theta_i from Methods, 'Population PK",
        "Model Development'). Perlstein 2025 Results states the direction",
        "explicitly: 'Lower BMI was associated with higher KA1 values (faster",
        "release via the direct route) and decreased KA2 (slower release via",
        "the indirect route)', which the authors attribute to the slower route",
        "being distribution into subcutaneous fat. The centring value is the",
        "cohort MEAN of 28.7 kg/m^2 (Table 2, Total Overall column); the paper",
        "does not print the median that its own equation centres on. See the",
        "vignette Assumptions and deviations section.",
        sep = " "
      ),
      source_name        = "BMI"
    ),
    DOSE_TV46000_ML = list(
      description        = "Volume of TV-46000 suspension delivered at this injection",
      units              = "mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Per-dose-record covariate. Enters as a median-normalized power term",
        "on the fast direct-release rate constant only:",
        "ka1 * (DOSE_TV46000_ML / 0.303)^(-0.384) (Perlstein 2025 Table 3,",
        "KA1INJV1 = -0.384). TV-46000 is a fixed-strength suspension of",
        "approximately 360 mg risperidone per mL, so the volume is",
        "proportional to the milligram dose: the paper prints 12.5 mg =",
        "0.035 mL, 25 mg = 0.07 mL, 50 mg = 0.139 mL (Results, 'Risperidone",
        "(Parent) PK Model') and 100 mg = 0.278 mL (Table 2, BA-10148 column,",
        "SD 0 in a single-dose 100 mg study). A downstream user dosing in mg",
        "can derive the column as dose_mg / 360. The centring value is the",
        "cohort MEAN of 0.303 mL over 3287 injections (Table 2, Overall",
        "column); the paper does not print the median. See the vignette",
        "Assumptions and deviations section.",
        sep = " "
      ),
      source_name        = "INJV"
    ),
    INJSITE_ARM = list(
      description        = "Subcutaneous injection site is the upper arm rather than the abdomen",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (abdomen)",
      notes              = paste(
        "Per-dose-record covariate; a subject may alternate sites between",
        "injections, and the paper's conclusion is precisely that the two",
        "sites are interchangeable. Enters through the paper's categorical",
        "covariate form Pi = P * (1 + theta_i)^COVi (Methods), so an arm",
        "injection multiplies ka1 by 1 + 0.331 = 1.331 (Perlstein 2025",
        "Table 3, KA1ADMSITE1 = 0.331; Results: 'Upper-arm administration was",
        "associated with a 33% higher KA1 (i.e., faster release) than",
        "injection in the abdomen'). Abdomen (68.5% of the 3287 injections) is",
        "the reference category, matching the canonical register.",
        sep = " "
      ),
      source_name        = "ADMSITE"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened in the stepwise covariate analysis and reached statistical",
        "significance, but 'Age and sex did not have clinically relevant",
        "effects on apparent clearance (CL/F) compared with the base model,",
        "and both were removed' (Perlstein 2025 Results, 'Risperidone (Parent)",
        "PK Model and the Effect of Covariates'). No point estimate is",
        "reported for the age effect, so it cannot be reconstructed. Patients",
        "older than 65 years were excluded from the source studies.",
        sep = " "
      )
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened and statistically significant on both the parent and the",
        "metabolite, but removed from both. On the parent it was not",
        "clinically relevant on CL/F; on the metabolite, 'The effect of sex on",
        "metabolite clearance (CLMO) was the only covariate identified to have",
        "a statistically significant impact ... however, the impact on CLMO",
        "was small (-1.7%) compared with the base model, and it was removed'",
        "(Perlstein 2025 Results). The -1.7% is the change in objective",
        "function relative to the base model, not a covariate coefficient; no",
        "point estimate for the sex effect is reported.",
        sep = " "
      )
    ),
    WT = list(
      description = "Total body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Screened as an intrinsic factor in the stepwise covariate analysis",
        "(Perlstein 2025 Methods) but not retained; the size descriptor that",
        "survived selection was BMI, on the two release rate constants.",
        sep = " "
      )
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units       = "mL/min",
      type        = "continuous",
      notes       = paste(
        "Screened as an intrinsic factor (Perlstein 2025 Methods) but not",
        "retained. Cohort mean 121 mL/min (SD 34.3; Table 2).",
        sep = " "
      )
    ),
    RACE_BLACK = list(
      description = "Black or African American race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Race was screened as an intrinsic factor (Perlstein 2025 Methods) but",
        "not retained. The cohort was 64.7% Black or African American, 31.9%",
        "White and 1.4% Asian (Table 2); the Discussion flags this skew, plus",
        "the male skew and the exclusion of patients older than 65 years, as",
        "the analysis's main generalizability limitation.",
        sep = " "
      )
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "risperidone", units = "mg", specimen = "administration site", verified = TRUE),
    transit1    = list(analyte = "risperidone", units = "mg", specimen = "administration site", verified = TRUE),
    transit2    = list(analyte = "risperidone", units = "mg", specimen = "administration site", verified = TRUE),
    transit3    = list(analyte = "risperidone", units = "mg", specimen = "administration site", verified = TRUE),
    transit4    = list(analyte = "risperidone", units = "mg", specimen = "administration site", verified = TRUE),
    transit5    = list(analyte = "risperidone", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "risperidone", units = "mg", specimen = "plasma", verified = TRUE),
    central_9oh = list(analyte = "9-hydroxyrisperidone (paliperidone)", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 692L,
    n_studies      = 5L,
    n_observations = paste(
      "Not reported as a record count. 692 participants contributed at least",
      "one measurable post-TV-46000 plasma concentration: 267 from the three",
      "phase 1 studies (50 RISPE1ZG15EU, 97 SAD-10055, 120 BA-10148) and 425",
      "from the two phase 3 studies (352 RISE, 73 SHINE). Two further",
      "participants were excluded for self-administration of risperidone and",
      "41 because all their PK samples were below the limit of quantification",
      "(Perlstein 2025 Results, 'Dataset', and Table S1).",
      sep = " "
    ),
    age_range      = "16-65 years; mean 47.4 (SD 10.9) years (Perlstein 2025 Table 2, Total Overall column). Patients older than 65 years were excluded from the source studies.",
    weight_range   = "42-132 kg; mean 86.3 (SD 16.3) kg (Perlstein 2025 Table 2)",
    bmi_range      = "18-38 kg/m^2; mean 28.7 (SD 4.8) kg/m^2 (Perlstein 2025 Table 2)",
    sex_female_pct = 29.6,
    race_ethnicity = c(Black = 64.7, White = 31.9, Asian = 1.4, Missing = 2.0),
    renal_function = "Creatinine clearance mean 121 mL/min (SD 34.3); not reported for the healthy-volunteer study RISPE1ZG15EU (Perlstein 2025 Table 2)",
    disease_state  = paste(
      "Adults with schizophrenia or schizoaffective disorder, plus 50 healthy",
      "volunteers (7.2% of the analysis population) from the phase 1 study",
      "RISPE1ZG15EU. Phase 3 patients were stabilized on oral risperidone for",
      "12 weeks before the first TV-46000 injection; the two phase 1 patient",
      "studies used a 7-day oral run-in followed by a 7-day washout, so most",
      "participants had undetectable plasma levels at the first injection.",
      sep = " "
    ),
    dose_range     = paste(
      "TV-46000 subcutaneous. Phase 1: single doses 12.5-25 mg",
      "(RISPE1ZG15EU, subtherapeutic, healthy volunteers), single doses",
      "50-225 mg and three consecutive monthly doses 75-150 mg (SAD-10055),",
      "single dose 100 mg from a vial or a prefilled syringe (BA-10148).",
      "Phase 3: q1m 50, 75, 100 or 125 mg and q2m 100, 150, 200 or 250 mg",
      "(RISE and SHINE). Table 4 simulations assume a 28-day month.",
      sep = " "
    ),
    regions        = "Not reported by region; the phase 3 trials RISE (NCT03503318) and SHINE (NCT03893825) were multicenter",
    notes          = paste(
      "Perlstein 2025 Table 2 tabulates demographics for the 733 participants",
      "with baseline data, while the popPK analysis used the 692 with at least",
      "one measurable concentration; the summary statistics quoted here are",
      "the Table 2 values. Injections were given in the abdomen (68.5%) or the",
      "upper arm (31.5%) and delivered from vials (77.8%) or prefilled",
      "syringes (22.2%) across 3287 injections. Product presentation (vial vs",
      "prefilled syringe) was screened as an extrinsic factor and had no",
      "effect on exposure (Conclusion), so it is neither a model covariate nor",
      "a covariatesDataExcluded entry -- no canonical vial-vs-prefilled-syringe",
      "column exists in the register, and the closest entry, DEVICE_AI,",
      "contrasts an autoinjector against a prefilled syringe rather than a",
      "prefilled syringe against a vial.",
      sep = " "
    )
  )

  ini({
    # ======================================================================
    # All structural estimates are Perlstein 2025 Table 3, "Final TV-46000 PK
    # model parameter estimates for risperidone (parent) and 9-OH
    # risperidone (metabolite)", Estimate (95% CI) column. Rate constants are
    # in 1/h; the q1m dosing interval is 28 days = 672 h and q2m is 56 days
    # = 1344 h (Methods, "Model-Based PK Simulations ... for Different Dosing
    # Regimens": "with the assumption that each month had 28 days").
    #
    # CL/F and V/F are apparent (bioavailability is not separately
    # identifiable from a subcutaneous-only dataset), so no lfdepot is
    # estimated and F is left at its rxode2 default of 1.
    # ======================================================================

    # ----------------------- Risperidone (parent) -------------------------
    lcl  <- log(14.3);      label("Apparent clearance of risperidone CL/F (L/h)")                                    # Table 3: CL/F = 14.3 L/h (95% CI 13.4-15.1; %RSE 2.9)
    lvc  <- log(66.3);      label("Apparent central volume of distribution of risperidone V/F (L)")                  # Table 3: V/F = 66.3 L (95% CI 60.6-71.9; %RSE 4.4)
    lka1 <- log(0.000632);  label("Fast direct-route release rate constant KA1, depot to central (1/h)")             # Table 3: KA1 = 0.000632 1/h (95% CI 0.000588-0.000676; %RSE 3.6)
    lka2 <- log(0.000408);  label("Slow indirect-route release rate constant KA2, depot to transit chain (1/h)")     # Table 3: KA2 = 0.000408 1/h (95% CI 0.000329-0.000486; %RSE 9.9)
    lktr <- log(0.0252);    label("Transit rate constant KTR through the five-compartment chain (1/h)")              # Table 3: KTR = 0.0252 1/h (95% CI 0.0228-0.0275; %RSE 4.8)

    # ------------------- Covariate effects on the release rates -----------
    # Continuous covariates enter as Pi = P * (COVi / COVmedian)^theta_i and
    # categorical covariates as Pi = P * (1 + theta_i)^COVi (Perlstein 2025
    # Methods, "Population PK Model Development"). The two centring values
    # (BMI, injection volume) are the cohort MEANS from Table 2, because the
    # paper never prints the medians its own equation refers to -- see the
    # vignette Assumptions and deviations section.
    e_bmi_ka1          <- -1.1;    label("Power exponent on (BMI / 28.7) for KA1 (unitless)")                        # Table 3: KA1BMI1 = -1.1 (95% CI -1.460 to -0.745 as printed "-1.460 to 0.745"; %RSE 16.6)
    e_dose_ka1         <- -0.384;  label("Power exponent on (DOSE_TV46000_ML / 0.303) for KA1 (unitless)")           # Table 3: KA1INJV1 = -0.384 (95% CI -0.469 to -0.298; %RSE 11.3)
    e_injsite_arm_ka1  <- 0.331;   label("Proportional upper-arm-vs-abdomen shift on KA1 (unitless)")                # Table 3: KA1ADMSITE1 = 0.331 (95% CI 0.109-0.553; %RSE 34.2); Results: 33% higher KA1 in the upper arm
    e_bmi_ka2          <- 1.7;     label("Power exponent on (BMI / 28.7) for KA2 (unitless)")                        # Table 3: KA2BMI1 = 1.7 (95% CI 0.992-2.410; %RSE 21.3)

    # ---------------- 9-hydroxyrisperidone (metabolite) -------------------
    # FRMET is fixed to 1 (Table 3 footnote: "VMO (apparent) volume of
    # distribution of metabolite, which can also be interpreted as V/FRMET
    # (fraction metabolized) fixed to 1"), so the whole parent elimination
    # flux forms metabolite and CLMO / VMO are apparent values conditioned on
    # that. No molar-mass correction is applied between the two analytes --
    # the paper prints none, and any such factor is absorbed into the
    # apparent CLMO and VMO. See the vignette Assumptions and deviations
    # section.
    lcl_9oh <- log(5.78);   label("Apparent clearance of 9-hydroxyrisperidone CLMO (L/h)")                           # Table 3: CLMO = 5.78 L/h (95% CI 5.52-6.04; %RSE 2.322)
    lvc_9oh <- log(95.7);   label("Apparent volume of distribution of 9-hydroxyrisperidone VMO (L)")                 # Table 3: VMO = 95.7 L (95% CI 88.5-103.0; %RSE 3.822)

    # ------------------ Interindividual variability -----------------------
    # Perlstein 2025 Methods: "Interindividual variability was evaluated for
    # all parameters and described via exponential error models", and Table 3
    # reports each IIV as a percentage. Those percentages are %CV on the
    # linear scale, so omega^2 = log(1 + CV^2). That reading -- rather than
    # omega = CV/100 -- is what reproduces the Table 3 confidence intervals:
    # propagating the reported %RSE as an RSE on omega through
    # CV = sqrt(exp(omega^2) - 1) gives 95% CIs of 200.6-322.1 (KA2),
    # 44.7-58.2 (KA1), 75.8-89.4 (CL) and 59.3-71.5 (CLMO) against the
    # printed 198.8-320.1, 44.0-57.5, 75.3-89.3 and 58.9-71.1; the omega =
    # CV/100 reading gives 229.5-281.5 for KA2, which is far too narrow.
    # The same log(1 + CV^2) convention is used by the upstream Ivaturi 2017
    # RBP-7000 model.
    #
    # KA2 and KA1 are correlated (Table 3, "Correlation IIV KA2 and KA1, %" =
    # 42.8), so they form a 2x2 block; the covariance is
    # 0.428 * sqrt(2.010 * 0.231) = 0.292. Table 3 reports no IIV on V/F, on
    # KTR or on VMO.
    etalka2 + etalka1 ~ c(2.010,
                          0.292, 0.231)  # Table 3: IIV KA2 = 254.2% CV -> omega^2 = log(1 + 2.542^2) = 2.010 (shrinkage 25.8%); IIV KA1 = 51.0% CV -> log(1 + 0.510^2) = 0.231 (shrinkage 32.7%); correlation 42.8% -> cov = 0.428 * sqrt(2.010 * 0.231) = 0.292
    etalcl            ~ 0.517            # Table 3: IIV on CL   = 82.3% CV -> omega^2 = log(1 + 0.823^2) = 0.517 (shrinkage 6.6%)
    etalcl_9oh        ~ 0.353            # Table 3: IIV on CLMO = 65.1% CV -> omega^2 = log(1 + 0.651^2) = 0.353 (shrinkage 5.5%)

    # ------------------------- Residual error -----------------------------
    # Perlstein 2025 Methods: "The 'log-transformed, both-sides' approach was
    # applied, using an additive residual error to describe the residual
    # variability in the data" -- an additive error on the log scale is a
    # proportional error in linear space, and Table 3's footnote names the
    # reported quantity "EP proportional residual error", so the printed
    # percentages are linear-space CVs and map directly onto propSd.
    propSd     <- 0.405;  label("Proportional residual error on risperidone plasma concentration (fraction)")            # Table 3, parent: EP = 40.5% (%RSE 2.3, shrinkage 5.8%)
    propSd_9oh <- 0.383;  label("Proportional residual error on 9-hydroxyrisperidone plasma concentration (fraction)")   # Table 3, metabolite: EP = 38.3% (%RSE 2.0, shrinkage 2.7%)

    # ------------- Dopamine D2-receptor occupancy (D2RO) ------------------
    # Perlstein 2025 Methods, "Prediction of Individual Exposure Parameters
    # from TAM Concentration-Time Profiles": "the D2RO profile was simulated
    # from the TAM profile with a simple maximum effect (Emax) PK/PD model,
    # where Emax = 100% and apparent equilibrium dissociation constant (kd)
    # was 10.1 ng/mL [22, 23]". Both values are literature constants the
    # authors applied to their simulated TAM profiles rather than parameters
    # estimated in this analysis, hence fixed().
    lemax <- fixed(log(100));  label("Maximum dopamine D2-receptor occupancy (percent); literature value")   # Methods: Emax = 100%
    kd    <- fixed(10.1);      label("Apparent equilibrium dissociation constant for D2RO (ng/mL); literature value")  # Methods: kd = 10.1 ng/mL
  })

  model({
    # 0. Constants. The molecular-weight ratio that converts a
    #    9-hydroxyrisperidone mass concentration into risperidone-equivalent
    #    units when the two analytes are summed. Perlstein 2025 Introduction
    #    defines the reported endpoint as "the total active moiety (TAM; sum
    #    of concentrations of risperidone and its active metabolite [9-OH
    #    risperidone] corrected by molecular weight)" but does not print the
    #    two molecular weights; 410 (risperidone, 410.49 g/mol) and 426
    #    (9-hydroxyrisperidone, 426.49 g/mol) are the rounded values used
    #    throughout the risperidone LAI literature and by the upstream
    #    Ivaturi 2017 model of the same two molecules. See the vignette
    #    Assumptions and deviations section.
    mw_ratio_9oh <- 410 / 426   # unitless

    # 1. Individual PK parameters.
    #
    # The two release rate constants carry the covariate effects. BMI moves
    # them in opposite directions: the exponent is negative on the fast
    # direct route (ka1) and positive on the slow indirect route (ka2), so a
    # leaner participant releases more drug by the direct path. Injection
    # volume and injection site act on ka1 only.
    ka1 <- exp(lka1 + etalka1) *
      (BMI / 28.7)^e_bmi_ka1 *
      (DOSE_TV46000_ML / 0.303)^e_dose_ka1 *
      (1 + e_injsite_arm_ka1)^INJSITE_ARM
    ka2 <- exp(lka2 + etalka2) * (BMI / 28.7)^e_bmi_ka2
    ktr <- exp(lktr)
    cl  <- exp(lcl + etalcl)
    vc  <- exp(lvc)

    cl_9oh <- exp(lcl_9oh + etalcl_9oh)
    vc_9oh <- exp(lvc_9oh)

    # 2. Micro-constants. FRMET is fixed to 1, so the parent's entire
    #    elimination flux kel * central becomes metabolite.
    kel     <- cl / vc
    kel_9oh <- cl_9oh / vc_9oh

    # 3. ODE system, in declaration order so the rxode2 compartment slots are
    #    stable (depot 1, transit1-5 2-6, central 7, central_9oh 8). This is
    #    Perlstein 2025 Figure 1: the single subcutaneous depot has two exit
    #    arrows, KA1 straight to the risperidone compartment ("Faster
    #    Absorption") and KA2 into transit compartment 1 ("Slower
    #    Absorption"); the chain then moves at KTR per stage, including the
    #    final KTR arrow from transit compartment 5 into risperidone. The
    #    implicit fast fraction is ka1 / (ka1 + ka2) = 60.8% at the covariate
    #    reference; the paper reports no separate bioavailability split.
    d/dt(depot)       <- -(ka1 + ka2) * depot
    d/dt(transit1)    <-  ka2 * depot    - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1 - ktr * transit2
    d/dt(transit3)    <-  ktr * transit2 - ktr * transit3
    d/dt(transit4)    <-  ktr * transit3 - ktr * transit4
    d/dt(transit5)    <-  ktr * transit4 - ktr * transit5
    d/dt(central)     <-  ka1 * depot + ktr * transit5 - kel * central
    d/dt(central_9oh) <-  kel * central - kel_9oh * central_9oh

    # 4. Observations. Amounts are in mg and volumes in L, so amount / volume
    #    is mg/L = ug/mL; the factor 1000 converts to the ng/mL the paper
    #    reports.
    Cc     <- 1000 * central     / vc
    Cc_9oh <- 1000 * central_9oh / vc_9oh

    # Total active moiety in risperidone-equivalent ng/mL -- the clinically
    # reported quantity (Perlstein 2025 Table 4 and Figures 2-4) -- and the
    # derived dopamine D2-receptor occupancy in percent (Figure 5).
    TAM  <- Cc + Cc_9oh * mw_ratio_9oh
    D2RO <- exp(lemax) * TAM / (kd + TAM)

    Cc     ~ prop(propSd)
    Cc_9oh ~ prop(propSd_9oh)
  })
}
