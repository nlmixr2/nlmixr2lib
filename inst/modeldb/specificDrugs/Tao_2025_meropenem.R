Tao_2025_meropenem <- function() {
  description <- "One-compartment IV population PK model for meropenem in 144 critically ill adult and elderly surgical-ICU patients with Pseudomonas aeruginosa infections (Tao 2025). Clearance is 4.68 L/h with no retained covariate; central volume is 4.47 L at the population-median age of 63.5 years and scales with age by a power exponent of 0.19, so distribution volume rises with advancing age. Residual variability is additive. Body weight, sex, serum creatinine, serum glucose, and renal function were screened but not retained. Developed to support age-stratified PK/PD target-attainment dosing (40% and 100% fT > 4x MIC)."
  reference <- "Tao R, Chan S, Wang J, Wang X, Shi L, Lan X, Chen L, Mu X. Pharmacokinetic-pharmacodynamic target attainment analyses as support for meropenem dosing regimens in critically ill adult and elderly patients with Pseudomonas aeruginosa infections. Front Pharmacol. 2025;16:1643553. doi:10.3389/fphar.2025.1643553"
  vignette <- "Tao_2025_meropenem"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Checked against the source: Tao 2025 fits a
  # one-compartment model to meropenem PLASMA concentrations measured by
  # LC-MS/MS (Methods 2.2, linear 0.1-100 mg/L); doses are in mg and vc is in L,
  # so `central` holds mg and Cc = central/vc is in mg/L.
  compartmentData <- list(
    central = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    AGE = list(
      description        = "Patient age at ICU admission",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tao 2025 Table 1: median (IQR) 63.50 (46.50, 76.0) years; Results section 3.1 additionally reports mean (SD) 60.63 (18.37) years. Patients under 18 years were excluded. The centering constant 63.5 in the structural equation V = 4.47 * (AGE/63.5)^0.19 is exactly the Table 1 median age, which confirms the covariate is median-normalized rather than referenced to a rounded standard. Age was the ONLY covariate retained in the final model (Results section 3.3) and acts on volume of distribution only; it entered the base model on both CL and V for a 10.36-point OFV drop but was retained on V alone. Dosing simulations were run at 20, 40, 60, and 90 years (Methods section 2.5).",
      source_name        = "Age"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened on CL and V; not retained. Tao 2025 Results section 3.3: 'potential clinical and demographic factors, including age, body weight, sex, serum creatinine, and renal function, were systematically evaluated for their impact on meropenem clearance and volume of distribution.' No weight summary is reported anywhere in the paper (Table 1 omits it), so no distribution is available."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units       = "(binary)",
      type        = "categorical",
      notes       = "Screened on CL and V; not retained. Tao 2025 Results section 3.3. Table 1 reports 86 males (59.72%) and 58 females (40.28%) of n = 144."
    ),
    CREAT = list(
      description = "Serum creatinine concentration",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened on CL and V; not retained. Tao 2025 Results section 3.3. Table 1: median (IQR) 40.00 (22.00, 62.75) umol/L."
    ),
    CRCL = list(
      description = "Renal function (creatinine clearance / estimated glomerular filtration rate)",
      units       = "mL/min/1.73 m^2",
      type        = "continuous",
      notes       = "Screened and eliminated. Tao 2025 Results section 3.3: eGFR entered forward selection but 'during the backward elimination step, estimated glomerular filtration rate (eGFR) was excluded from the model.' Discussion notes the absence of a CrCL effect is attributable to the confounded age trend ('A notable trend of decreasing CrCL was observed as the age of patients increased'). Clinical dosing in the study was nevertheless adjusted by measured CrCL (Methods section 2.1), so CrCL influenced the observed data through the dosing regimen rather than through a model parameter."
    ),
    GLU = list(
      description = "Serum glucose concentration",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Screened and eliminated. Tao 2025 Results section 3.3: 'The inclusion of estimated glucose (GLU) also resulted in a modest, but statistically significant, decrease of 2.16 points in the OFV', but age was 'the only covariate retained'. No glucose summary statistics are reported in the paper."
    ),
    ALB = list(
      description = "Serum albumin concentration",
      units       = "g/L",
      type        = "continuous",
      notes       = "Recorded at baseline and at sampling (Methods section 2.1) but not reported as retained on any parameter. Tao 2025 Table 1 gives median (IQR) 32.00 (17.90, 39.80) with the unit printed as 'mg/L'; 32 mg/L is not a physiologic serum albumin concentration and the values are consistent with g/L, so the printed unit is taken to be a typographical error (see the vignette Errata)."
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 144L,
    n_studies        = 1L,
    age_range        = "18 years and older (lower bound set by the exclusion criterion; no upper bound reported)",
    age_median       = "63.50 years (IQR 46.50-76.0); mean (SD) 60.63 (18.37) years",
    weight_range     = "Not reported (body weight was screened as a covariate but no summary statistics are given)",
    sex_female_pct   = 40.28,
    race_ethnicity   = "Not reported (single-center Chinese surgical ICU cohort)",
    disease_state    = "Critically ill adults and elderly patients admitted to a surgical intensive care unit with Pseudomonas aeruginosa infections. SOFA and APACHE II scores were recorded at admission and at sampling but are not reported numerically. Patients receiving renal replacement therapy or extracorporeal membrane oxygenation during meropenem administration were excluded, as were patients under 18 years.",
    dose_range       = "Meropenem 2000 mg IV loading dose over 30 min in all 144 patients (Table 1), immediately followed by 24-hour continuous infusion reconstituted three times per day. Maintenance dose 0.5 g in 45 patients (31.25%) and 1.0 g in 99 patients (68.75%); dosing interval q8h in 85 patients (59.03%) and q12h in 59 patients (40.97%). Clinical regimens were adjusted according to measured creatinine clearance.",
    regions          = "China (single center: The Second Affiliated Hospital of Guizhou Medical University, Kaili, Guizhou). Surgical ICU admissions between 1 March 2023 and 1 October 2024.",
    renal_function   = "Renal replacement therapy was an exclusion criterion. Serum creatinine median (IQR) 40.00 (22.00, 62.75) umol/L. eGFR was screened as a covariate and eliminated during backward elimination.",
    n_concentrations = 144L,
    notes            = "Limited PK sampling plus therapeutic drug monitoring data; 144 blood samples from 144 patients (approximately one sample per patient), drawn after a minimum of 6 h of meropenem therapy. Observed meropenem concentration median (IQR) 10.40 (3.39, 22.53) mg/L. Five samples (3.47%) were below the 0.1 mg/L LOQ and were handled by the M3 method. Assay LC-MS/MS, linear 0.1-100 mg/L. Estimation in Phoenix NLME 7.0 using FOCE with interaction; final-model OFV 873.81 for the one-compartment model versus 928.18 for the two-compartment model. Shrinkage 0.35 (Vd) and 0.12 (CL); condition number 10.93. Validated by 1000-replicate nonparametric bootstrap (Table 2) and prediction-corrected VPC (Figure 3)."
  )

  ini({
    # Structural parameters (Tao 2025 Table 2, final-model "Estimate" column).
    # The structural equations are printed as a display equation in Results
    # section 3.3 (p. 4):
    #   CL (L/h) = 4.68 * exp(eta_CL)
    #   Vd (L)   = 4.47 * (Age / 63.5)^theta_Age,V * exp(eta_Vd)
    # Reference subject: AGE = 63.5 years (the Table 1 median age).
    lcl <- log(4.68); label("Clearance (L/h)")                  # Tao 2025 Table 2: CL = 4.68 L/h (RSE 8.38%; bootstrap median 4.71, 95% CI 3.36-6.11)
    lvc <- log(4.47); label("Central volume at AGE = 63.5 years (L)")  # Tao 2025 Table 2: Vd = 4.47 L (RSE 13.23%; bootstrap median 4.55, 95% CI 3.89-5.64)

    # Covariate effect: power model of age on central volume,
    #   vc = exp(lvc) * (AGE / 63.5)^e_age_vc
    #
    # VALUE CONFLICT IN THE SOURCE (resolved in favour of Table 2; see the
    # vignette Errata for the full falsification). Tao 2025 reports this
    # exponent twice with different values:
    #   - Table 2 row "theta Age,V": 0.19 (RSE 35.62%; bootstrap median 0.18,
    #     95% CI -0.04 to 0.35)
    #   - the display equation in Results section 3.3 (verified against the
    #     typeset PDF, not just the text extraction): (Age/63.5)^0.29
    #
    # The display equation is DEMONSTRABLY CORRUPTED, which is what breaks the
    # tie. As typeset it reads
    #     CL(L/h) = 4.68 * e^0.045
    #     Vd(L)   = 4.47 * (Age/63.5)^0.29 * e^0.016
    # i.e. the authors mechanically pasted Table 2 numbers into the symbol
    # slots and put the VARIANCES omega^2_CL = 0.045 and omega^2_Vd = 0.016
    # into exp() where the random effects eta_CL and eta_Vd belong. exp(0.045)
    # is a constant, not a random effect, so two of the equation's three
    # substituted values are provably wrong; there is no reason to trust the
    # third over the parameter table it was copied from.
    #
    # Three INDEPENDENT lines of evidence favour 0.19:
    #   (1) Table 2 itself is internally corroborated: the point estimate
    #       (0.19) and the bootstrap median (0.18) agree with each other, and
    #       section 3.4 states the bootstrap medians were "within 5% of the
    #       population parameter estimates ... for all parameters" -- which
    #       holds at 0.19 (5.3%) and fails at 0.29 (37.9%). The bootstrap 95%
    #       CI (-0.04, 0.35) is centred near 0.19.
    #   (2) ARITHMETIC ANSWER KEY. Results section 3.3, immediately below the
    #       display equation, reports "the estimated age-normalized Vd at
    #       steady state were 4.41 (4.21-4.62) L". Read as the typical-value
    #       curve evaluated across the Table 1 age IQR (46.50, 76.0), each
    #       endpoint pins the exponent independently. (That reading is an
    #       inference; the sentence could instead summarise the DISTRIBUTION of
    #       the 144 individual values Figure 2 plots. The exact two-sided match
    #       argues for the curve reading, but nothing here depends on it --
    #       under the other reading these numbers carry no information about
    #       the exponent, leaving (1) and (3) untouched. Neither reading
    #       supports 0.29.) Back-solving from each endpoint SEPARATELY,
    #       assuming nothing:
    #         from 4.21 at age 46.5 -> theta = 0.1923
    #         from 4.62 at age 76.0 -> theta = 0.1837
    #       Two independent constraints both land on ~0.19. Forward-checked:
    #         theta = 0.19: 4.213 (pub 4.21), 4.625 (pub 4.62)  -> exact to 2 dp
    #         theta = 0.29: 4.084 (off 0.13),  4.709 (off 0.09)
    #   (3) FIGURE 2 (p. 7) is a third, wholly independent answer key: it plots
    #       individual age-normalized Vd against four age bands with group
    #       central bars at approximately 4.01, 4.27, 4.57 and 4.86 L (values
    #       digitized off the rendered figure -- see the vignette). Scored at
    #       the band midpoints (30, 50, 75, 92 y):
    #         theta = 0.19: 3.88 4.27 4.61 4.80  RMSE 0.077 L
    #         theta = 0.29: 3.60 4.17 4.69 4.98  RMSE 0.229 L
    #       Least-squares fitting the exponent to those four band centres alone
    #       gives theta_hat = 0.167, near the bootstrap median 0.18.
    #
    # COUNTER-EVIDENCE, disclosed for completeness. The published median of
    # 4.41 L does NOT reproduce at 0.19 evaluated at the median age (any
    # exponent gives exactly 4.47 at Age = 63.5), whereas 4.47*(60.63/63.5)^
    # 0.29 = 4.4104 reproduces it to four decimal places using the section-3.1
    # MEAN age of 60.63 y. This needs no exponent change to explain. Figure 2
    # plots the individual age-normalized values spanning roughly 3.15-5.15 L,
    # i.e. they carry the random effect, so 4.41 is an EMPIRICAL median over
    # 144 patients, not a curve evaluation. Its 1.3% shortfall from 4.47 L is
    # 1.0 standard errors of a sample median at n = 144 with the 12.7% CV
    # implied by omega^2_Vd = 0.016 (asserted in the vignette). Under 0.29 the
    # two interval endpoints would instead correspond to ages 51.7 and 71.2 y,
    # which match no quantity reported anywhere in the paper.
    e_age_vc <- 0.19; label("Power exponent on (AGE / 63.5) for central volume")  # Tao 2025 Table 2: theta Age,V = 0.19 (RSE 35.62%)

    # Inter-individual variability (Tao 2025 Table 2, "Between-subject
    # variation" block, reported directly as VARIANCES omega^2 on the
    # log-normal scale: theta_i = theta_TV * exp(eta_i)). The variance scale is
    # settled by the Table 2 row labels themselves, which read "omega^2 CL" and
    # "omega^2 Vd" rather than "omega", and is corroborated by Results section
    # 3.2 restating them as "4.50% and 1.60% in the final model" -- the
    # variances written as percentages. The RSEs (15.90% and 43.66%) both
    # exceed the sqrt(2/N) = 11.8% lower bound that a variance estimate must
    # satisfy at N = 144, so they are consistent with the variance reading; note
    # this direction is only a consistency check, since an SD's RSE is bounded
    # by roughly half that and would clear 11.8% too (the falsifier runs the
    # other way -- an RSE BELOW the bound would prove an SD).
    etalcl ~ 0.045  # Tao 2025 Table 2: omega^2 CL = 0.045 (RSE 15.90%) -> 21.4% CV
    etalvc ~ 0.016  # Tao 2025 Table 2: omega^2 Vd = 0.016 (RSE 43.66%) -> 12.7% CV

    # Residual error. Tao 2025 Results section 3.2: "Residual variability was
    # adequately accounted for using an additive error model, as the combined
    # proportional and additive error model did not offer a significant
    # improvement."
    addSd <- 0.74; label("Additive residual error (mg/L)")  # Tao 2025 Table 2: sigma additive = 0.74 mg/L (RSE 24.20%; bootstrap median 0.72, 95% CI 0.01-1.07)
  })
  model({
    # Individual PK parameters. Clearance carries no covariate; central volume
    # scales with age by a power function centred at the population-median
    # 63.5 years (Tao 2025 Results section 3.3 display equation).
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc) * (AGE / 63.5)^e_age_vc

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose in mg, vc in L -> central/vc has units mg/L, matching the paper's
    # observed concentrations (LC-MS/MS linear over 0.1-100 mg/L).
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
