JacoboCabral_2015_tacrolimus <- function() {
  description <- paste0(
    "Two-compartment population PK model for oral tacrolimus in Mexican ",
    "paediatric renal-transplant recipients (Jacobo-Cabral 2015): first-order ",
    "absorption with a lag time, no allometric scaling, three-level CYP3A5 ",
    "genotype effect on apparent oral clearance (*3/*3 reference, *1/*3 ",
    "+50%, *1/*1 +93%), formulation-type effects on Ka and on relative ",
    "bioavailability F (pooled Prograf + Framebin + Tenacrine reference ",
    "vs Limustin generic vs unrecorded), an exponential per-dose effect on F ",
    "centred at 2 mg, exponential inter-patient variability on Ka, V/F and F, ",
    "and a residual error described in the paper as additive on the natural-log ",
    "concentration scale (encoded as proportional residual error in linear ",
    "space, which is the standard nlmixr2 equivalent for SD <= 0.15)."
  )
  reference <- paste0(
    "Jacobo-Cabral CO, Garcia-Roca P, Romero-Tejeda EM, Reyes H, Medeiros M, ",
    "Castaneda-Hernandez G, Troconiz IF. Population pharmacokinetic analysis ",
    "of tacrolimus in Mexican paediatric renal transplant patients: role of ",
    "CYP3A5 genotype and formulation. Br J Clin Pharmacol. 2015;80(4):630-641. ",
    "doi:10.1111/bcp.12649"
  )
  vignette <- "JacoboCabral_2015_tacrolimus"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    CYP3A5_STAR1_HET = list(
      description        = "CYP3A5*1/*3 heterozygote indicator (one functional CYP3A5*1 allele)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP3A5*3/*3 nonexpresser, when paired with CYP3A5_STAR1_HOM = 0)",
      notes              = paste0(
        "Time-fixed per subject (germline genotype at rs776746). ",
        "1 = CYP3A5*1/*3 heterozygote (one functional *1 allele); ",
        "0 = otherwise (the union of *3/*3 nonexpressers and *1/*1 homozygotes; ",
        "the paired indicator CYP3A5_STAR1_HOM flags the homozygous group). ",
        "Jacobo-Cabral 2015 cohort (Table 1): 21/53 (40%) were *1/*3 heterozygotes, ",
        "3/53 (5.7%) were *1/*1 homozygotes, 29/53 (54%) were *3/*3 nonexpressers. ",
        "The three-level decomposition is encoded with two binary indicators ",
        "(SLCO1B1_HAP15 / Passey 2011 precedent); reference category is the ",
        "*3/*3 nonexpresser stratum with both indicators = 0."
      ),
      source_name        = "CYP3A5*1/*3"
    ),
    CYP3A5_STAR1_HOM = list(
      description        = "CYP3A5*1/*1 homozygote indicator (two functional CYP3A5*1 alleles)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP3A5*3/*3 nonexpresser, when paired with CYP3A5_STAR1_HET = 0)",
      notes              = paste0(
        "Time-fixed per subject (germline genotype at rs776746). ",
        "1 = CYP3A5*1/*1 homozygote (two functional *1 alleles); ",
        "0 = otherwise (the union of *3/*3 nonexpressers and *1/*3 heterozygotes; ",
        "the paired indicator CYP3A5_STAR1_HET flags the heterozygous group). ",
        "Jacobo-Cabral 2015 cohort (Table 1): 3/53 (5.7%) were *1/*1 homozygotes. ",
        "Paired with CYP3A5_STAR1_HET to encode the three-level genotype factor."
      ),
      source_name        = "CYP3A5*1/*1"
    ),
    FORM_TAC_LIMUSTIN = list(
      description        = "Limustin Mexican generic tacrolimus formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (pooled Prograf + Framebin + Tenacrine; reference Ka and F)",
      notes              = paste0(
        "Time-fixed per subject (formulation does not switch within subject in this cohort). ",
        "1 = subject received Limustin (a Mexican generic IR tacrolimus product); ",
        "0 = subject received any of the three reference-pool formulations (Prograf ",
        "innovator, Framebin generic, Tenacrine generic) OR formulation was unknown ",
        "(the unknown stratum is captured separately by FORM_TAC_UNK). The paper ",
        "groups the Prograf + Framebin + Tenacrine formulations as the reference ",
        "after backward elimination found no significant difference between them. ",
        "Jacobo-Cabral 2015 cohort (Methods): 9/53 (17%) received Limustin; the ",
        "reference pool covers 37/53 (70%; Prograf 29, Framebin 5, Tenacrine 3); ",
        "7/53 (13%) had unknown formulation (-> FORM_TAC_UNK = 1)."
      ),
      source_name        = "FOR"
    ),
    FORM_TAC_UNK = list(
      description        = "Tacrolimus unrecorded-formulation indicator (Jacobo-Cabral 2015 cohort missingness stratum)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (formulation documented)",
      notes              = paste0(
        "Time-fixed per subject. 1 = subject's tacrolimus formulation was not ",
        "recorded in the source cohort (n = 7 of 53; 13%); ",
        "0 = subject's formulation was documented (any of Prograf, Limustin, ",
        "Framebin or Tenacrine). The Jacobo-Cabral 2015 final model retained ",
        "the unknown stratum as a separate categorical level with its own Ka ",
        "and F coefficients rather than imputing or excluding. Downstream ",
        "simulators normally set FORM_TAC_UNK = 0 for every virtual subject ",
        "(formulation is known by design); the indicator is included for ",
        "faithful reproduction of the published equations only."
      ),
      source_name        = "FOR"
    ),
    DOSE = list(
      description        = "Current per-dose administered tacrolimus amount",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Per-dose administered amount in mg. Used inside the dose-dependent ",
        "bioavailability term `F_DTOT = exp(e_dtot_fdepot * (DOSE - 2))` (Jacobo-",
        "Cabral 2015 Table 2 'F = e^[theta10 * (Dose - 2)]'; theta10 = -0.30), so ",
        "the median 2 mg dose maps to F_DTOT = 1 and higher doses give F_DTOT < 1. ",
        "Jacobo-Cabral 2015 cohort (Table 1): per-dose median 2 mg, range 0.5-6 mg ",
        "(weighted dose median 0.047 mg/kg, range 0.009-0.268 mg/kg in 11.2-75.5 kg ",
        "paediatric recipients). The covariate column is set on each observation row ",
        "to the per-dose amount that produced the most recent dosing event; for ",
        "steady-state q12h dosing in the original cohort, DOSE is constant across ",
        "each subject. This is the same canonical 'DOSE' covariate used elsewhere ",
        "in nlmixr2lib (Zheng 2016 sifalimumab, Hansson 2013 sunitinib, etc.); the ",
        "paper labels it DTOT."
      ),
      source_name        = "DTOT"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 53L,
    n_studies       = 1L,
    n_observations  = 405L,
    age_range       = "2-19 years (median 16)",
    age_median      = "16 years",
    weight_range    = "11.2-75.5 kg (mean 48.2, SD 15.2; median 48)",
    weight_median   = "48 kg",
    sex_female_pct  = 35.8,
    race_ethnicity  = "Mexican paediatric renal-transplant recipients (single-centre, Federico Gomez Children's Hospital, Mexico City)",
    disease_state   = paste0(
      "Paediatric renal-transplant recipients on standardized maintenance ",
      "immunosuppression (tacrolimus + mycophenolate mofetil +/- prednisone). ",
      "Standardized triple therapy in 44 of 53 subjects (tacrolimus + MMF + ",
      "prednisone); 9 received tacrolimus + MMF only. Patients were in a stable ",
      "post-transplant condition (mean post-operative day 391.6, SD 327.2; median ",
      "244 days, range 50-1230)."
    ),
    dose_range      = paste0(
      "Oral tacrolimus 0.5-6 mg twice daily (every 12 h) titrated to target ",
      "trough concentrations of 5-10 ng/mL. Per-dose median 2 mg; weighted dose ",
      "median 0.047 mg/kg/dose. Patients received one of Prograf (innovator, ",
      "n = 29), Limustin (generic, n = 9), Framebin (generic, n = 5), Tenacrine ",
      "(generic, n = 3); 7 patients had unknown formulation."
    ),
    regions         = "Mexico (Federico Gomez Children's Hospital of Mexico, Mexico City)",
    sampling_design = paste0(
      "One full PK profile per patient at steady state (after morning dose): ",
      "pre-dose plus 0.5, 1, 2, 3, 4, 6, 8 and 12 h postdose (nine timepoints per ",
      "subject; 405 evaluable concentrations total, two LLOQ samples excluded). ",
      "Whole-blood tacrolimus measured by chemiluminescent microparticle ",
      "immunoassay (ARCHITECT, Abbott Laboratories)."
    ),
    cyp3a5_genotype = c(`*3/*3` = 54.7, `*1/*3` = 39.6, `*1/*1` = 5.7),
    formulation_distribution = c(
      Prograf = 54.7, Limustin = 17.0, Framebin = 9.4, Tenacrine = 5.7, Unknown = 13.2
    ),
    notes           = paste0(
      "Software: NONMEM v7.2 (FOCE + INTERACTION); PsN v3.4.2 for stepwise ",
      "covariate modelling (SCM) and bootstrap; Xpose / R 3.0.1 for diagnostics. ",
      "Steady-state q12h dosing modelled using the NONMEM SS = 1 option. ",
      "Baseline demographics: Jacobo-Cabral 2015 Table 1. Final-model parameters: ",
      "Jacobo-Cabral 2015 Table 2 (bootstrap n = 2000). The CYP3A5 effect on CL/F ",
      "explained almost all of the inter-patient variability in CL/F (base IPV 30% ",
      "-> final model IPV negligible). The formulation effects on Ka and F reduced ",
      "the Ka IPV from 79% to 37% and the F IPV from 72% to 38%. ABCB1 genotype ",
      "and the other covariates (age, weight, BSA, Cr, Hb, Hct, AST, ALT, total ",
      "protein, albumin, post-operative days, prednisone / verapamil co-medication) ",
      "were tested but were NOT significant in the final model. Bioavailability F ",
      "was fixed at 1 (anchored at the reference subject / per-dose median); no IV ",
      "data were available to estimate absolute F."
    )
  )

  ini({
    # Structural PK parameters -- Jacobo-Cabral 2015 Table 2 final estimates.
    # Reference subject: CYP3A5 *3/*3 nonexpresser, on the pooled
    # Prograf + Framebin + Tenacrine formulation, receiving the median 2 mg
    # per-dose amount. All values are population means; bootstrap medians
    # (n = 2000) are nearly identical (see vignette source-trace table).
    lka  <- log(0.52);  label("Absorption rate constant Ka at the reference formulation pool (1/h)")  # Jacobo-Cabral 2015 Table 2, theta1 = 0.52 (RSE 27%, IPV 37%)
    lcl  <- log(11.98); label("Apparent oral clearance CL/F at CYP3A5*3/*3 reference (L/h)")           # Jacobo-Cabral 2015 Table 2, theta3 = 11.98 (RSE 8%)
    lvc  <- log(24.16); label("Apparent central volume of distribution V/F (L)")                       # Jacobo-Cabral 2015 Table 2, V/F = 24.16 (RSE 39%, IPV 66%)
    lvp  <- log(383.5); label("Apparent peripheral volume of distribution V_T/F (L)")                  # Jacobo-Cabral 2015 Table 2, V_T/F = 383.5 (RSE 34%)
    lq   <- log(32.49); label("Apparent inter-compartmental clearance Q/F (L/h)")                      # Jacobo-Cabral 2015 Table 2, Q/F = 32.49 (RSE 20%)
    ltlag <- log(0.39); label("Absorption lag time (h)")                                               # Jacobo-Cabral 2015 Table 2, t_lag = 0.39 (RSE 6%)
    lfdepot <- fixed(log(1)); label("Relative bioavailability F at the reference subject (anchor, unitless)")  # Jacobo-Cabral 2015 Table 2 footnote: F = 100*, was not estimated (fixed at 1 due to absence of IV data)

    # CYP3A5 genotype effects on CL/F -- Jacobo-Cabral 2015 Table 2 covariate
    # equation 'CL/F = theta3 * INF_CYP3A5' with INF_CYP3A5 = 1 + theta8 (if
    # *1/*3) or 1 + theta9 (if *1/*1) or 1 (if *3/*3 reference). Both
    # indicators are 0 for the *3/*3 reference; the paired indicators are
    # mutually exclusive so at most one contributes a nonzero deviation.
    e_cyp3a5_het_cl <- 0.50; label("CL/F linear-deviation effect of CYP3A5*1/*3 heterozygote (+50%; unitless)")  # Jacobo-Cabral 2015 Table 2, theta8 = 0.50 (RSE 38%)
    e_cyp3a5_hom_cl <- 0.93; label("CL/F linear-deviation effect of CYP3A5*1/*1 homozygote (+93%; unitless)")    # Jacobo-Cabral 2015 Table 2, theta9 = 0.93 (RSE 33%)

    # Formulation effects on Ka -- Jacobo-Cabral 2015 Table 2 covariate
    # equation 'Ka (h^-1) = theta1 * Ka_FOR' with Ka_FOR = 1 + theta13
    # (Limustin) or 1 + theta14 (Unknown) or 1 (reference pool). The
    # indicators are mutually exclusive (a subject is either in the
    # documented-non-Limustin reference pool, the Limustin stratum, or the
    # unknown stratum).
    e_form_limustin_ka <- -0.76; label("Ka linear-deviation effect of Limustin formulation (-76%; unitless)")  # Jacobo-Cabral 2015 Table 2, theta13 = -0.76 (RSE 6%)
    e_form_unk_ka      <- -0.51; label("Ka linear-deviation effect of unknown formulation (-51%; unitless)")   # Jacobo-Cabral 2015 Table 2, theta14 = -0.51 (RSE 23%)

    # Formulation effects on relative bioavailability F -- Jacobo-Cabral
    # 2015 Table 2 covariate equation 'F = 100 * F_DTOT * F_FOR' with
    # F_FOR = 1 + theta11 (Limustin) or 1 + theta12 (Unknown) or 1
    # (reference pool). The reference subject's F is anchored at 100% so
    # the relative-bioavailability multipliers are unitless.
    e_form_limustin_fdepot <- -0.53; label("F linear-deviation effect of Limustin formulation (-53%; unitless)")  # Jacobo-Cabral 2015 Table 2, theta11 = -0.53 (RSE 22%)
    e_form_unk_fdepot      <- -0.53; label("F linear-deviation effect of unknown formulation (-53%; unitless)")   # Jacobo-Cabral 2015 Table 2, theta12 = -0.53 (RSE 16%)

    # Per-dose-amount effect on F -- Jacobo-Cabral 2015 Table 2 covariate
    # equation 'F_DTOT = exp(theta10 * (Dose - 2))', where Dose is the
    # per-dose tacrolimus amount in mg (median 2 mg) and theta10 = -0.30.
    # Higher per-dose amounts reduce relative bioavailability (consistent
    # with saturable first-pass metabolism / dose-dependent absorption).
    e_dose_fdepot <- -0.30; label("Exponential per-dose effect on F, centred at 2 mg (1/mg)")  # Jacobo-Cabral 2015 Table 2, theta10 = -0.30 (RSE 19%)

    # Inter-patient variability (IPV) -- Jacobo-Cabral 2015 Table 2 reports
    # IPV as %CV. Convert to log-scale variance via omega^2 = log(1 + CV^2):
    #   Ka  CV 37% -> log(1 + 0.37^2) = log(1.1369) = 0.1283
    #   V/F CV 66% -> log(1 + 0.66^2) = log(1.4356) = 0.3613
    #   F   CV 38% -> log(1 + 0.38^2) = log(1.1444) = 0.1349
    # IPV on CL/F was reduced to negligible by the CYP3A5 covariate effect
    # (Jacobo-Cabral 2015 Results: 'The inclusion of the effect of CYP3A5
    # genotype on CL/F made the estimate of the IPV on CL/F negligible')
    # and is NOT carried in the final model -- omitted here. The
    # off-diagonal elements of the Omega matrix were tested and were not
    # significant (Jacobo-Cabral 2015 Methods, Base population model).
    etalka     ~ 0.1283  # Jacobo-Cabral 2015 Table 2 IPV Ka = 37% CV
    etalvc     ~ 0.3613  # Jacobo-Cabral 2015 Table 2 IPV V/F = 66% CV
    etalfdepot ~ 0.1349  # Jacobo-Cabral 2015 Table 2 IPV F = 38% CV

    # Residual error -- Jacobo-Cabral 2015 Table 2 reports a residual SD of
    # 0.12 on the natural-log concentration scale ('Residual error
    # [ln(ng/mL)] 0.12, RSE 8%'). The paper's Methods state 'residual error
    # was described using an additive error model on the logarithmic
    # scale'. For sigma <= 0.15 the additive-on-log form is numerically
    # equivalent to a proportional residual error in linear space with
    # propSd ~ sigma; using prop() here matches the convention used by the
    # other tacrolimus models in nlmixr2lib (Storset 2014 propSd = 0.149;
    # Bergmann 2014 propSd = 0.183) and reproduces the observed-vs-
    # predicted spread in Jacobo-Cabral 2015 Figure 2.
    propSd <- 0.12; label("Proportional residual error (fraction; encoded from the paper's additive-on-log SD = 0.12)")  # Jacobo-Cabral 2015 Table 2, residual SD = 0.12 on ln scale
  })

  model({
    # ----- 1. Individual structural PK parameters -----
    # Ka and F absorb the formulation effects via additive-multiplier
    # factors; CL/F absorbs the CYP3A5 effects similarly. Reference
    # subject: CYP3A5 *3/*3 nonexpresser, on pooled Prograf + Framebin +
    # Tenacrine formulation, receiving the 2 mg per-dose amount.

    # Formulation factor on Ka: Ka_FOR = 1 for the reference pool
    # (both indicators 0), 1 - 0.76 = 0.24 for Limustin, and 1 - 0.51 =
    # 0.49 for the unknown stratum.
    ka_for <- 1 +
              e_form_limustin_ka * FORM_TAC_LIMUSTIN +
              e_form_unk_ka      * FORM_TAC_UNK

    # Formulation factor on F: F_FOR = 1 for the reference pool,
    # 1 - 0.53 = 0.47 for Limustin, and 1 - 0.53 = 0.47 for the unknown
    # stratum.
    f_for <- 1 +
             e_form_limustin_fdepot * FORM_TAC_LIMUSTIN +
             e_form_unk_fdepot      * FORM_TAC_UNK

    # Dose-dependent F factor: F_DTOT = exp(theta10 * (DOSE - 2)). At the
    # 2 mg median dose, F_DTOT = 1; at 0.5 mg, F_DTOT = exp(0.30 * 1.5) =
    # 1.57; at 6 mg, F_DTOT = exp(-0.30 * 4) = 0.30.
    f_dose <- exp(e_dose_fdepot * (DOSE - 2))

    # CYP3A5 factor on CL/F: INF_CYP3A5 = 1 for *3/*3 (both indicators
    # 0); 1 + 0.50 = 1.50 for *1/*3; 1 + 0.93 = 1.93 for *1/*1.
    inf_cyp3a5 <- 1 +
                  e_cyp3a5_het_cl * CYP3A5_STAR1_HET +
                  e_cyp3a5_hom_cl * CYP3A5_STAR1_HOM

    # Individual PK parameters with covariate effects and log-normal IIV.
    ka   <- exp(lka   + etalka)     * ka_for
    cl   <- exp(lcl)                * inf_cyp3a5
    vc   <- exp(lvc   + etalvc)
    q    <- exp(lq)
    vp   <- exp(lvp)
    tlag <- exp(ltlag)

    # Relative bioavailability F: anchored at 1 for the reference subject
    # (lfdepot fixed at log(1) = 0); carries the formulation factor
    # (Limustin / unknown) and the dose-dependent factor. IIV on F is
    # exponential, applied as exp(etalfdepot) on the linear scale.
    fdepot <- exp(lfdepot + etalfdepot) * f_for * f_dose

    # ----- 2. Micro-rate constants -----
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ----- 3. ODE system (two-compartment with first-order absorption + lag) -----
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    # Bioavailability and lag-time applied to the depot compartment.
    f(depot)    <- fdepot
    alag(depot) <- tlag

    # ----- 4. Observation -----
    # Tacrolimus whole-blood concentration in ng/mL. central / vc gives
    # mg/L; multiply by 1000 to convert to ng/mL (= ug/L).
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
