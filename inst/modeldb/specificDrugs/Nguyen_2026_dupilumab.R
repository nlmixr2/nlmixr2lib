Nguyen_2026_dupilumab <- function() {
  description <- "Two-compartment population PK model with a three-compartment transit absorption chain and parallel linear plus Michaelis-Menten elimination for subcutaneous dupilumab in healthy adults and in adults, adolescents and children with eosinophilic esophagitis"
  reference <- "Nguyen JH, Chehade M, Dellon ES, Radin A, Chittenden J, Kamal MA, Louisias M, Xu C, Kosloski MP. Population Pharmacokinetics of Dupilumab in Adults, Adolescents, and Children With Eosinophilic Esophagitis. Clin Pharmacol Ther. 2026. doi:10.1002/cpt.70233"
  vignette <- "Nguyen_2026_dupilumab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Body weight. Baseline weight in adults; time-varying weight in children and adolescents, obtained by linear interpolation between observed measurements (Nguyen 2026 Methods, 'Time-varying body weight').",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric power scaling with a 70 kg reference (WGTREF = 70 in the supplementary control stream).",
        "Three separate exponents: e_wt_cl = 1.08 on linear CL, e_wt_vc_vp = 0.710 shared by Vc and Vp (the",
        "paper reports this as the effect on Vss = Vc + Vp), and e_wt_q = 0.75 fixed on Q. Replacing baseline",
        "weight with time-varying weight was a retained forward step (Table S5 step 2), so simulations of growing",
        "children should supply WT as a time-varying column; rxode2's default LOCF covariate interpolation should",
        "be changed to covsInterpolation = 'linear' to match the source's linear interpolation between measurements."
      ),
      source_name = "WGT / WGTBL"
    ),
    ALB = list(
      description = "Baseline serum albumin concentration.",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power scaling on linear CL with a 45 g/L reference (ALBBLREF = 45 in the supplementary control stream),",
        "exponent -1.16. Reported by the source in g/L, which is already the canonical SI unit, so no unit",
        "conversion is applied inside model(). Cohort mean was 46.4 g/L (Table S3)."
      ),
      source_name = "ALBBL"
    ),
    DIS_EOE = list(
      description = "Eosinophilic esophagitis patient indicator: 1 = patient with EoE, 0 = healthy volunteer.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy adult volunteer; the reference cohort is the 202 healthy adults from the six single-dose phase I studies of Table S1)",
      notes = paste(
        "The structural parameters were estimated from healthy-volunteer data alone and then held fixed; the",
        "EoE effects on CL, Vmax and Vss were estimated when the patient data were added (Nguyen 2026 Methods,",
        "'Base model construction'). Setting DIS_EOE = 0 therefore recovers the published healthy-volunteer",
        "structural model exactly. Note that the residual-error magnitudes also differed between the two",
        "populations; see the model's population$notes."
      ),
      source_name = "EOE"
    )
  )

  # Covariates that the source screened but did not retain in the final model.
  # Documented here so the provenance of the covariate search is preserved
  # without carrying "declared but not referenced" convention warnings.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Baseline age.",
      units = "years",
      type = "continuous",
      notes = paste(
        "Screened on CL, Ka and F1 (Table S2). An age effect on Ka entered in forward selection (Table S5 step 1)",
        "but was removed in backward elimination (Table S5 step 5, dOFV +1.56, P = 0.212); the corresponding",
        "THETA(22) KAAGE is 0.0 FIX in the supplementary control stream. A residual eta-vs-age trend on CL was",
        "seen (Figure S3) but judged physiologically implausible and not pursued. The abstract states plainly",
        "that 'age was not a significant covariate'."
      )
    ),
    EOS = list(
      description = "Baseline peak esophageal intraepithelial eosinophil count (EEOS).",
      units = "cells/uL",
      type = "continuous",
      notes = paste(
        "Screened on CL, Vss and Vmax (Table S2) as a disease-severity marker; not retained. THETA(19) CLEEOS and",
        "THETA(21) VMAXEOSS are both 0.0 FIX in the supplementary control stream. The source reports EEOS in",
        "eos/hpf at 400x magnification (cohort mean 84.6, Table S3), which is a microscopy field count rather",
        "than the volumetric blood count the canonical EOS column carries; this mismatch is a further reason the",
        "column is documented rather than used."
      )
    ),
    ADA_POS = list(
      description = "Maximum anti-drug-antibody titer category.",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Screened on CL as two indicator levels (Table S4 run 4204); the step did not minimize and was rejected",
        "(dOFV +6.83, P = 1.00). THETA(23) CLADA1 and THETA(24) CLADA2 are both 0.0 FIX in the supplementary",
        "control stream. The source states that 'development of anti-dupilumab antibodies was rare and did not",
        "impact functional dupilumab concentrations in these studies'."
      )
    )
  )

  compartmentData <- list(
    depot = list(analyte = "dupilumab", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "dupilumab", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "dupilumab", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "dupilumab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "dupilumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "dupilumab", units = "mg", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 632,
    n_studies = 9,
    n_observations = 4459,
    age_range = "1 to 18+ years (children >= 1 to < 12 years, n = 98; adolescents >= 12 to < 18 years, n = 97; adults >= 18 years, n = 235, in the EoE cohort)",
    age_median = "24.3 years (mean, EoE cohort overall; children 7.11, adolescents 15.1, adults 35.2)",
    weight_range = ">= 5 kg (EoE KIDS entry criterion); >= 40 kg in LIBERTY EoE TREET",
    weight_median = "65.7 kg (mean, EoE cohort overall; children 27.2, adolescents 63.8, adults 82.6)",
    sex_female_pct = 34,
    race_ethnicity = c(White = 90, `Black or African American` = 5, Asian = 2, `Other or not reported` = 3),
    disease_state = "eosinophilic esophagitis (430 patients), plus 202 healthy adult volunteers contributing the dense single-dose PK that identifies the structural model",
    dose_range = "IV 1-12 mg/kg and SC 75-600 mg single doses in healthy adults; SC 100/200/300 mg qw, q2w or q4w weight-tiered multiple doses in EoE",
    regions = "not reported by the source",
    notes = paste(
      "Baseline characteristics of the EoE cohort are Table S3; the healthy-volunteer characteristics are",
      "reported in the cited Li 2020 reference rather than in this paper. Nine EoE patients with implausible",
      "concentration rises more than 40 days after the last recorded dose were excluded (Figure S1). 22% of",
      "samples were below the limit of quantification and were handled by the M3 likelihood method.",
      "RESIDUAL ERROR IS POPULATION- AND ROUTE-SPECIFIC IN THE SOURCE and only the EoE terms are encoded here:",
      "the model file carries the estimated EoE residual (proportional 0.252, additive 11.6 mg/L). The source",
      "additionally fixed a healthy-volunteer IV residual (proportional 0.136, no additive term) and a",
      "healthy-volunteer SC residual (proportional 0.179, additive 0.0284 mg/L) -- Table 2 and the $ERROR block",
      "of the supplementary control stream, which switches on POP and ROUTN. rxode2 carries one error model per",
      "endpoint, so a user simulating healthy volunteers (DIS_EOE = 0) should substitute those magnitudes.",
      "Estimation was FOCEI with a preconditioning step in NONMEM 7.5.0; the final model condition number was 2.05."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Structural parameters. Every one of these was estimated from the
    # healthy-volunteer data alone and then held fixed while the patient
    # covariate effects were estimated -- Table 2 footnote a, 'Parameters
    # were fixed to the values obtained from the model with information from
    # healthy volunteers only'. Hence fixed() throughout this block.
    # Values are typical values for a 70 kg healthy volunteer.
    # ---------------------------------------------------------------------
    lcl     <- fixed(log(0.145))  ; label("Linear clearance in a 70 kg healthy volunteer (L/day)")                  # Table 2, row 'CL, L/day' = 0.145
    lvc     <- fixed(log(2.39))   ; label("Central volume of distribution in a 70 kg healthy volunteer (L)")        # Table 2, row 'Vc, L' = 2.39
    lq      <- fixed(log(0.511))  ; label("Intercompartmental clearance in a 70 kg healthy volunteer (L/day)")      # Table 2, row 'Q, L/day' = 0.511
    lvp     <- fixed(log(1.47))   ; label("Peripheral volume of distribution in a 70 kg healthy volunteer (L)")     # Table 2, row 'Vp, L' = 1.47
    lvmax   <- fixed(log(1.07))   ; label("Maximum nonlinear elimination rate in a healthy volunteer (mg/L/day)")   # Table 2, row 'Vmax, mg/L/day' = 1.07
    lkm     <- fixed(log(0.134))  ; label("Concentration of half-maximum nonlinear clearance (mg/L)")               # Table 2, row 'KM, mg/L' = 0.134
    lka     <- fixed(log(0.284))  ; label("First-order absorption rate constant out of the last transit compartment (1/day)")  # Table 2, row 'Ka, 1/day' = 0.284
    lfdepot <- fixed(log(0.659))  ; label("Subcutaneous bioavailability (fraction)")                                # Table 2, row 'F1' = 0.659; the source estimated this on a logit scale (THETA(8) = 0.658249, FABS = exp(x)/(1+exp(x))), which is numerically identical here because the final model carries no IIV on F1
    lmtt    <- fixed(log(0.0726)) ; label("Mean transit time through the absorption transit chain (day)")           # Table 2, row 'MTT, days' = 0.0726

    # ---------------------------------------------------------------------
    # Covariate effects, estimated in the final model (Table 2, lower block).
    # Continuous covariates enter as power models centred on a reference
    # value; the categorical EoE indicator enters as an exponential shift.
    # Both forms are given in Supplementary Methods, 'Description of
    # covariate and pharmacokinetic parameter relationships'.
    # ---------------------------------------------------------------------
    e_wt_cl         <- 1.08        ; label("Body-weight exponent on linear clearance (unitless)")                                        # Table 2, row 'Time-varying weight on CL (REF: 70 kg)' = 1.08 (1.07, 1.08)
    e_wt_vc_vp      <- 0.710       ; label("Body-weight exponent shared by central and peripheral volume (unitless)")                    # Table 2, row 'Time-varying weight on Vss (REF: 70 kg)' = 0.710 (0.706, 0.714); Vss = Vc + Vp with the effect shared
    e_wt_q          <- fixed(0.75) ; label("Body-weight exponent on intercompartmental clearance (unitless)")                            # Table 2, row 'Time-varying weight on Q (REF: 70 kg)' = 0.75 (fixed); estimated allometry on Q was removed in backward elimination, Table S5 step 4
    e_alb_cl        <- -1.16       ; label("Baseline serum albumin exponent on linear clearance (unitless)")                             # Table 2, row 'Baseline albumin on CL (REF: 45 g/L)' = -1.16 (-1.17, -1.15)
    e_dis_eoe_cl    <- log(0.944)  ; label("log ratio of linear clearance for patients with EoE vs healthy volunteers (unitless)")       # Table 2, row 'Patient with EoE on CL (REF: healthy volunteer)' = 0.944; the table prints the back-transformed ratio
    e_dis_eoe_vmax  <- log(0.782)  ; label("log ratio of maximum nonlinear elimination rate for patients with EoE vs healthy volunteers (unitless)")  # Table 2, row 'Patient with EoE on Vmax (REF: healthy volunteer)' = 0.782
    e_dis_eoe_vc_vp <- log(1.26)   ; label("log ratio of central and peripheral volume for patients with EoE vs healthy volunteers (unitless)")       # Table 2, row 'Patient with EoE on Vss (REF: healthy volunteer)' = 1.26; shared by Vc and Vp

    # ---------------------------------------------------------------------
    # Interindividual variability. Table 2 reports these on the NONMEM
    # variance scale; its footnote c gives the back-transform
    # CV(%) = sqrt(exp(IIV) - 1) * 100, which reproduces the quoted 31.9%
    # and 16.2% exactly and so confirms the values are log-scale variances.
    # A single random effect is shared by Vc and Vp, matching the source's
    # ETA_VSS being added to both TVVC and TVVP.
    # ---------------------------------------------------------------------
    etalcl    ~ 0.0970   # Table 2, row 'IIV on CL' = 0.0970 (0.0969, 0.0971); footnote c gives 31.9% CV, eta-shrinkage 8.62%
    etalvc_vp ~ 0.0260   # Table 2, row 'IIV on Vss' = 0.0260 (0.0260, 0.0261); footnote c gives 16.2% CV, eta-shrinkage 64.9%

    # ---------------------------------------------------------------------
    # Residual error in patients with EoE. See population$notes for the two
    # additional healthy-volunteer residual-error regimes that the source's
    # $ERROR block switches to on POP and ROUTN and that rxode2 cannot carry
    # simultaneously.
    # ---------------------------------------------------------------------
    propSd <- 0.252 ; label("Proportional residual error in patients with EoE (fraction)")  # Table 2, row 'Proportional error in patients with EoE' = 0.252 (0.250, 0.254)
    addSd  <- 11.6  ; label("Additive residual error in patients with EoE (mg/L)")           # Table 2, row 'Additive error in patients with EoE, mg/L' = 11.6 (11.5, 11.8)
  })

  model({
    # -- Individual parameters -------------------------------------------
    # Power covariate model on the continuous covariates (Supplementary
    # Methods eq. for continuous covariates) and exponential shift for the
    # categorical EoE indicator (Supplementary Methods eq. for categorical
    # covariates). Reference values WGTREF = 70 kg and ALBBLREF = 45 g/L are
    # listed in the 'COV REFERENCE VALUES' block of the supplementary
    # control stream and restated in the Table 2 row headers.
    cl <- exp(lcl + e_dis_eoe_cl * DIS_EOE + etalcl) * (WT / 70)^e_wt_cl * (ALB / 45)^e_alb_cl
    q  <- exp(lq) * (WT / 70)^e_wt_q
    # Vc and Vp share both the EoE shift and the weight exponent, and share a
    # single random effect: the source adds ETA_VSS to TVVC and to TVVP and
    # sets VCWT = VPWT = VSSWT. This is what the paper means by Vss "with
    # covariate and random effects shared".
    vc <- exp(lvc + e_dis_eoe_vc_vp * DIS_EOE + etalvc_vp) * (WT / 70)^e_wt_vc_vp
    vp <- exp(lvp + e_dis_eoe_vc_vp * DIS_EOE + etalvc_vp) * (WT / 70)^e_wt_vc_vp

    ka   <- exp(lka)
    mtt  <- exp(lmtt)
    vmax <- exp(lvmax + e_dis_eoe_vmax * DIS_EOE)
    km   <- exp(lkm)

    # Three transit transfers at rate ktr carry the subcutaneous dose from
    # the injection-site depot to the absorption compartment, which empties
    # into central at rate ka. NN = 3 and KTR = NN / MTT in the supplementary
    # control stream, so MTT is the mean time spent traversing the chain.
    ktr <- 3 / mtt

    # -- Micro-constants --------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # -- ODE system -------------------------------------------------------
    # Matches DADT(1)-DADT(6) of the supplementary $DES block. Compartment
    # map: depot = DOSESC (A1), transit1-3 = TR1-TR3 (A4-A6), central = A2,
    # peripheral1 = PERIPH (A3).
    Cc <- central / vc

    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * (depot - transit1)
    d/dt(transit2)    <-  ktr * (transit1 - transit2)
    d/dt(transit3)    <-  ktr * transit2 - ka * transit3
    d/dt(central)     <-  ka * transit3 - kel * central + k21 * peripheral1 - k12 * central -
      central * vmax / (km + Cc)
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # -- Bioavailability --------------------------------------------------
    # F1 applies to the subcutaneous injection-site compartment (F1 = FABS on
    # compartment 1 in the supplementary $PK block). An intravenous dose is
    # given directly into central and is therefore unaffected.
    f(depot) <- exp(lfdepot)

    # -- Observation and error --------------------------------------------
    Cc ~ add(addSd) + prop(propSd)
  })
}
