Osawa_2018_asunaprevir <- function() {
  description <- "One-compartment population PK model with first-order absorption and linear apparent clearance for oral asunaprevir (ASV, hepatitis C virus NS3/4A protease inhibitor) in Japanese adults with chronic genotype-1 HCV infection (Osawa 2018). Fit by NONMEM 7.2 FOCE to 2626 plasma concentrations from 265 subjects across two trials (AI447017, AI447026) receiving asunaprevir 200 or 600 mg twice daily as a phase-2 film-coated tablet or 100 mg twice daily as a phase-3 soft-gel capsule, always in combination with daclatasvir. Typical apparent clearance is 52.1 L/h and apparent central volume 75.1 L; because Ka (0.228 1/h) is well below the elimination rate constant CL/F divided by V/F (0.694 1/h), the disposition is flip-flop and the terminal slope reflects absorption. Inter-individual variability is carried on CL/F and V/F, the random effect on Ka was fixed to zero by the authors for lack of peak sampling, and the residual error is additive on the natural-log scale (log-transform-both-sides) with inter-individual variability on its magnitude. Four covariates survived backward elimination: baseline AST (power -0.598) and time-varying AST (power -0.382), both normalised to 52 U/L, plus compensated cirrhosis (exponential) on CL/F, and formulation acting on bioavailability. Unlike the companion daclatasvir model, every one of these effects exceeds the 80-125 percent boundaries: the soft-gel capsule has 1.37-fold higher bioavailability than the tablet, and asunaprevir exposure rises with cirrhosis and with worsening AST."

  reference <- paste(
    "Osawa M, Ueno T, Ishikawa H, Imai Y, Garimella T. (2018).",
    "Population Pharmacokinetic Analysis for Daclatasvir and Asunaprevir",
    "in Japanese Subjects With Chronic Hepatitis C Virus Infection.",
    "The Journal of Clinical Pharmacology 58(11):1468-1478.",
    "doi:10.1002/jcph.1274.",
    sep = " "
  )

  vignette <- "Osawa_2018_daclatasvir_asunaprevir"

  # CL/F is reported in L/h and V/F in L (Table 2), so amounts in mg give
  # concentrations in mg/L = ug/mL.
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    AST_BASE = list(
      description = "Baseline serum aspartate aminotransferase, time-fixed per subject. Enters CL/F as the power function (AST_BASE / 52)^-0.598.",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Reference 52 U/L (BASTref), stated in the text immediately below the asunaprevir covariate equations and equal to the Table 1 cohort median. Paired with the time-varying AST column in the same CL/F expression -- the source paper explicitly separates the baseline snapshot from the time-varying value, finding that adding the time-varying term dropped the objective function by a further 139.192 units. Both terms are normalised to the SAME constant 52 U/L; the printed equation divides the time-varying AST by BASTref, not by a separate reference. The negative exponent reproduces the paper's stated forest-plot values exactly: at the 5th and 95th percentiles of baseline AST (22 and 123.6 U/L), CL/F is 1.67-fold higher and 0.60-fold lower than at the 52 U/L reference, matching the paper's '1.67-fold higher' and '0.60-fold lower'. AST was retained in preference to the highly-correlated ALT, which gave a slightly smaller objective-function drop.",
      source_name = "BAST (baseline aspartate aminotransferase)"
    ),
    AST = list(
      description = "Time-varying serum aspartate aminotransferase at each observation time. Enters CL/F as the power function (AST / 52)^-0.382, alongside the separate baseline-AST term.",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Normalised to 52 U/L -- the SAME BASTref constant used by the baseline-AST term, per the printed equation. The Table 2 footnote defines this covariate as 'AST, aspartate aminotransferase at each time'. The encoded exponent reproduces the paper's stated forest-plot values exactly: at the 5th and 95th percentiles of AST at the last sampling time (16 and 115 U/L), CL/F is 1.57-fold higher and 0.74-fold lower than at the 52 U/L reference, matching the paper's '1.57-fold higher' and '0.74-fold lower'. A subject whose AST never moves from its baseline therefore carries a combined AST exponent of -0.598 + -0.382 = -0.980 on (AST/52). NOTE: the Figure 3 caption describes the forest plot's reference subject as having a 'time-varying AST of 26 U/L'; that is a display choice for the forest plot only and is not the model's reference, which the equation text fixes at 52 U/L for both AST terms.",
      source_name = "AST (aspartate aminotransferase at each time)"
    ),
    DIS_CIRRHOSIS = list(
      description = "Compensated cirrhosis indicator. 1 = subject has cirrhosis; 0 = no cirrhosis. Multiplies CL/F by exp(-0.428), i.e. asunaprevir apparent clearance is 0.65-fold that of a non-cirrhotic subject, raising exposure correspondingly.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no cirrhosis; CIRRHOSISref is stated as NO in the text below the covariate equations).",
      notes = "Time-fixed per subject. 8.3 percent of the asunaprevir cohort had cirrhosis (22 of 265, Table 1). Only compensated cirrhosis was enrolled -- the trials excluded decompensated disease -- so this indicator should not be read as spanning Child-Pugh B/C, where the dedicated hepatic-impairment study found 5- to 10-fold and 20- to 30-fold exposure increases. exp(-0.428) = 0.652, matching the paper's stated '0.65-fold lower'. The paper notes that cirrhotic subjects also had higher baseline and time-varying AST, so the cirrhosis and AST terms are partly collinear markers of the same underlying loss of hepatic metabolic capacity.",
      source_name = "CIRRHOSIS (cirrhosis, yes vs no)"
    ),
    FORM_ASV_SOFTGEL = list(
      description = "Asunaprevir formulation indicator. 1 = the phase-3 100 mg twice-daily soft-gel capsule; 0 = the phase-2 film-coated tablet (200 or 600 mg twice daily). Acts on relative bioavailability, which in an oral-only model appears as a common multiplier exp(-0.314) on both CL/F and V/F.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (phase-2 film-coated tablet; FORMref is stated as the phase-2 tablet formulation in the text below the covariate equations).",
      notes = "Time-fixed per subject in this analysis -- 16.2 percent of the cohort received the phase-2 formulation (43 of 265) and 83.8 percent the phase-3 formulation (222 of 265, Table 1), each subject receiving only one. Table 2 tabulates a SINGLE parameter for this effect, labelled 'F ~ FORM' (theta10 = -0.314), while the printed equations apply it to CL/F and V/F separately; those are the same parameter, because a change in bioavailability F is not separately identifiable from oral data and manifests identically in the apparent quantities CL/F and V/F. Scaling both by exp(-0.314) leaves the elimination rate constant CL/V unchanged and raises concentrations by 1 / exp(-0.314) = 1.369, which is precisely the paper's stated 'bioavailability of the phase 3 soft-gel capsule formulation was 1.37-fold higher than the tablet formulation'. The Figure 3 caption confirms the reading: 'Relative bioavailability was computed from CL/F (phase 3 formulation, soft-gel capsule)/CL/F (phase 2 formulation, film-coated tablet), or V/F (phase 3 formulation)/V/F (phase 2 formulation).' Distinct from the registered FORM_ASV_LIQUID, which contrasts liquid against solid asunaprevir formulations in the Wang 2018 meta-analysis; here both arms are solid oral forms and the contrast is soft-gel capsule versus film-coated tablet.",
      source_name = "FORM (formulation; phase 2 film-coated tablet vs phase 3 soft-gel capsule)"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "asunaprevir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "asunaprevir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 265L,
    n_studies = 2L,
    n_observations = 2626L,
    age_range = "24-75 years",
    age_median = "62 years",
    weight_range = "36-93 kg",
    weight_median = "55 kg",
    sex_female_pct = 65.7,
    race_ethnicity = "Japanese (all subjects; eligibility required Japanese men or women aged 20 years or older)",
    disease_state = "Chronic hepatitis C virus genotype-1 infection. 8.3 percent had compensated cirrhosis. By prior-treatment history: 40.8 percent nonresponder / null responder / partial responder, 59.2 percent standard-of-care-ineligible naive or intolerant.",
    hepatic_function = "Baseline AST median 52 U/L; baseline ALT median 55 U/L. Both AST terms and cirrhosis status were retained in the final model, and the paper interprets asunaprevir CL/F as a marker of hepatic function.",
    renal_function = "Baseline creatinine clearance median 84.6 mL/min (range 39.56-171.95); not retained in the final model",
    genotype = "OATP1B1 haplotype was screened (*1B/*1B 19.3 percent, *1B/*1A 32.5 percent, *1A/*1A 10.6 percent, other 24.2 percent, missing 13.6 percent) and showed no statistically significant effect on CL/F.",
    dose_range = "Asunaprevir 200 or 600 mg twice daily as a phase-2 film-coated tablet, or 100 mg twice daily as a phase-3 soft-gel capsule, orally",
    regimens = "Always co-administered with daclatasvir (the all-oral DUAL regimen)",
    regions = "Japan",
    notes = "Baseline demographics from Osawa 2018 Table 1 (asunaprevir panel). Data pooled from two trials -- AI447017 and AI447026 -- listed in Supplementary Table 1, which was not available when this model was built; the per-trial breakdown is therefore not reproduced here. Of 2676 collected pharmacokinetic records, 50 were excluded as below the lower limit of quantification, leaving 2626 for model development. NOTE: the Table 1 asunaprevir panel is headed 'Covariate (n = 256)', but every count in that panel sums to 265 (91 + 174 male/female; 22 + 243 cirrhosis; 43 + 222 formulation) and the Results text states 265 subjects, so the header is a transposition typo for 265. A second apparent transposition in that panel: the asunaprevir AST and ALT rows carry maxima of 595 and 377 U/L respectively, the reverse of the daclatasvir panel's 377 (AST) and 595 (ALT); only the AST median of 52 U/L is load-bearing for the model and it is consistent across the table, the equation text and the forest-plot arithmetic."
  )

  ini({
    # =========================================================================
    # Structural fixed effects -- Osawa 2018 Table 2, asunaprevir panel.
    # The typical values are those of the reference subject: phase-2
    # film-coated tablet, baseline AST 52 U/L, time-varying AST 52 U/L, no
    # cirrhosis. (The Results text additionally describes the typical subject
    # as a 62-year-old 55 kg woman, but age, weight and sex did not survive
    # backward elimination for asunaprevir and are not model inputs.)
    # =========================================================================
    lcl <- log(52.1); label("Apparent clearance CL/F at the reference covariate values (L/h)") # Table 2 asunaprevir theta1 = 52.1 L/h (RSE 14.1%, bootstrap 95% CI 42.1-66.8)
    lvc <- log(75.1); label("Apparent central volume of distribution V/F at the reference covariate values (L)") # Table 2 asunaprevir theta2 = 75.1 L (RSE 18.0%, bootstrap 95% CI 56.7-101.0)
    lka <- log(0.228); label("First-order absorption rate constant Ka (1/h)") # Table 2 asunaprevir theta3 = 0.228 1/h (RSE 1.73%, bootstrap 95% CI 0.221-0.234)

    # =========================================================================
    # Covariate effects -- Osawa 2018 Table 2 and the two covariate equations
    # printed in the Results 'Asunaprevir' subsection:
    #
    #   CL/F_TV = CL/F_TV,ref * (BASTb / BASTref)^CL/F_BAST
    #                         * (ASTb  / BASTref)^CL/F_AST
    #                         * exp(FORM * CL/F_FORM
    #                               + CIRRHOSIS * CL/F_CIRRHOSIS)
    #   V/F_TV  = V/F_TV,ref  * exp(FORM * V/F_FORM)
    #
    # with BASTref 52 U/L, FORMref the phase-2 tablet, CIRRHOSISref NO. Note
    # that BOTH AST terms divide by BASTref -- there is no separate reference
    # for the time-varying term.
    #
    # Table 2 tabulates ONE formulation parameter (theta10, labelled 'F ~
    # FORM'), so CL/F_FORM and V/F_FORM above are the same number: a
    # bioavailability effect is not separately identifiable from oral data and
    # enters CL/F and V/F identically. That is why it is named with the shared-
    # exponent form `e_<cov>_<param1>_<param2>`.
    # =========================================================================
    e_ast_base_cl <- -0.598; label("Power exponent of baseline AST on CL/F, normalised to 52 U/L (unitless)") # Table 2 asunaprevir theta7 = -0.598 (RSE 9.45%, bootstrap 95% CI -0.707 to -0.486)
    e_ast_cl <- -0.382; label("Power exponent of time-varying AST on CL/F, normalised to 52 U/L (unitless)") # Table 2 asunaprevir theta8 = -0.382 (RSE 11.6%, bootstrap 95% CI -0.458 to -0.303)
    e_dis_cirrhosis_cl <- -0.428; label("Exponential effect of compensated cirrhosis on CL/F, no-cirrhosis reference (unitless)") # Table 2 asunaprevir theta11 = -0.428 (RSE 32.7%, bootstrap 95% CI -0.753 to -0.146)
    e_form_asv_softgel_cl_vc <- -0.314; label("Exponential effect of the phase-3 soft-gel capsule on both CL/F and V/F, phase-2 tablet reference (unitless)") # Table 2 asunaprevir theta10 = -0.314 (RSE 45.9%, bootstrap 95% CI -0.551 to -0.0922); single tabulated parameter labelled 'F ~ FORM'

    # =========================================================================
    # Inter-individual variability -- Osawa 2018 Table 2 'Random effects'.
    # Table 2 footnote b: diagonal elements are printed as variance (standard
    # deviation), so the first number in each cell is the value nlmixr2 wants.
    # No CL/F-V/F covariance row is tabulated for asunaprevir, so the block is
    # diagonal (unlike the companion daclatasvir model).
    #
    # Ka carries NO eta: 'The random effect of Ka was fixed as zero because the
    # sampling points around the peak concentrations were not enough for all
    # subjects' (Results, Asunaprevir). Encoding a zero-variance eta would add
    # a singular row to OMEGA, so the eta is omitted entirely and the authors'
    # decision is recorded here and in the vignette.
    # =========================================================================
    etalcl ~ 0.172 # Table 2 asunaprevir omega1,1 = 0.172 (SD 0.415)
    etalvc ~ 0.872 # Table 2 asunaprevir omega2,2 = 0.872 (SD 0.934)

    # =========================================================================
    # Residual error -- Osawa 2018 Table 2 'Residual error' plus the omega4,4
    # row that Table 2 places in the 'Random effects' block under the symbol
    # sigma. The Methods residual model is log-transform-both-sides:
    #
    #   ln(y_ij) = ln(yhat_ij) + theta_ADD * epsilon_ij
    #
    # so the residual is additive on the natural-log scale, which is exactly
    # nlmixr2's lnorm() error structure and the canonical `expSd` name.
    # theta4 is that log-scale SD; omega4,4 is inter-individual variability on
    # its magnitude, i.e. the NONMEM W = THETA(4) * EXP(ETA(4)) construct, so
    # the per-subject residual SD is expSd * exp(etaexpSd). The two rows carry
    # distinct standard errors (RSE 2.84% for theta4, 26.2% for omega4,4),
    # confirming they are separately estimated parameters rather than one value
    # printed twice.
    # =========================================================================
    expSd <- 0.68; label("Residual error SD, additive on the natural-log scale (log-scale SD)") # Table 2 asunaprevir theta4 = 0.68 (RSE 2.84%, bootstrap 95% CI 0.648-0.713)
    etaexpSd ~ 0.0672 # Table 2 asunaprevir omega4,4 = 0.0672 (SD 0.259), listed under 'Random effects' with the symbol sigma; IIV on the residual-error magnitude
  })

  model({
    # -----------------------------------------------------------------------
    # 1. Formulation effect on relative bioavailability.
    #
    # One estimated parameter multiplies BOTH CL/F and V/F (see the ini() note
    # on theta10). Because the same factor divides out of CL/V, the elimination
    # rate constant is formulation-independent and only the concentration level
    # moves: concentrations scale by 1 / exp(-0.314) = 1.369 for the soft-gel
    # capsule relative to the tablet.
    # -----------------------------------------------------------------------
    frel_form <- exp(e_form_asv_softgel_cl_vc * FORM_ASV_SOFTGEL)

    # -----------------------------------------------------------------------
    # 2. Individual PK parameters. The covariate equations are transcribed
    # verbatim from the Results section; the reference constant 52 U/L is the
    # BASTref value stated in the sentence immediately following them and is
    # used by BOTH AST terms.
    # -----------------------------------------------------------------------
    cl <- exp(lcl + etalcl) *
      (AST_BASE / 52)^e_ast_base_cl *
      (AST / 52)^e_ast_cl *
      exp(e_dis_cirrhosis_cl * DIS_CIRRHOSIS) *
      frel_form
    vc <- exp(lvc + etalvc) * frel_form
    ka <- exp(lka)

    # -----------------------------------------------------------------------
    # 3. Micro-constant and ODE system. One compartment with first-order
    # absorption from an oral depot and linear elimination. Bioavailability is
    # not identifiable from oral-only data and is absorbed into CL/F and V/F,
    # so the formulation effect above is applied to those apparent quantities
    # rather than through an f(depot) term.
    #
    # Note that ka (0.228 1/h) is well below kel (52.1/75.1 = 0.694 1/h at the
    # reference covariates), so the system is flip-flop: the terminal slope of
    # the concentration-time profile reflects absorption, not elimination.
    # -----------------------------------------------------------------------
    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # -----------------------------------------------------------------------
    # 4. Observation and residual error. The per-subject residual SD carries
    # its own eta (see the ini() note on omega4,4).
    # -----------------------------------------------------------------------
    Cc <- central / vc
    expSdInd <- expSd * exp(etaexpSd)
    Cc ~ lnorm(expSdInd)
  })
}
