# Population pharmacokinetic model for oral lumefantrine in the WWARN
# pooled meta-analysis of 4,122 patients with uncomplicated Plasmodium
# falciparum malaria (Kloprogge 2018, PLOS Medicine 15(6):e1002579;
# doi:10.1371/journal.pmed.1002579).

Kloprogge_2018_lumefantrine <- function() {
  description <- paste(
    "Population PK model for oral lumefantrine in 1,347 patients",
    "(children, non-pregnant adults, and second-/third-trimester",
    "pregnant women) from 26 studies in 12 African, Oceanian, and",
    "Southeast Asian countries with uncomplicated Plasmodium falciparum",
    "malaria treated with the standard fixed-dose artemether-lumefantrine",
    "regimen (Kloprogge 2018 PLOS Medicine). Two-compartment disposition",
    "with first-order absorption; F fixed at 1 with log-normal IIV",
    "(Box-Cox shape -0.343 on the F IIV departure from log-normal not",
    "reproduced here -- see Errata); allometric scaling of CL/F and Q/F",
    "(power 3/4) and of Vc/F and Vp/F (power 1) on body weight centered",
    "at the model-building median 42 kg; dose-saturable absorption on F",
    "with Dose50 = 3.86 mg/kg; exponential effect of log10 admission",
    "parasitaemia on F centered at log10(15,800/uL) = 4.2 (coefficient",
    "-0.643 per log10 unit); proportional pregnancy effect on ka",
    "(+35.2% in second and third trimester). IIV on Vc/F (CV 144%) and F",
    "(CV 70.3%); additive log-scale residual SD 0.323.",
    sep = " "
  )
  reference <- paste(
    "Kloprogge F, Workman L, Borrmann S, Tekete M, Lefevre G, Hamed K,",
    "et al. (2018). Artemether-lumefantrine dosing for malaria treatment",
    "in young children and pregnant women: A pharmacokinetic-",
    "pharmacodynamic meta-analysis. PLOS Medicine 15(6):e1002579.",
    "doi:10.1371/journal.pmed.1002579.",
    sep = " "
  )
  vignette <- "Kloprogge_2018_lumefantrine"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject at admission. Kloprogge 2018 pooled the",
        "PK model-building dataset across children, non-pregnant adults,",
        "and pregnant women (n = 1,347; Table 1 'Model building data'",
        "row, body weight range 6-150 kg, median 42 kg). Body weight is",
        "applied as an allometric scaler with exponent 3/4 on CL/F and",
        "Q/F and exponent 1 on Vc/F and Vp/F, centered at the model-",
        "building median 42 kg (Table 2 footnote: 'Clearance and volume",
        "parameters were centred on the median body weight (WT) and",
        "scaled allometrically: CL and Q = theta(n) * (WT/42)^(3/4),",
        "V = theta(n) * (WT/42)').",
        sep = " "
      ),
      source_name        = "WT"
    ),
    PREG = list(
      description        = "Pregnancy status indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "1 = pregnant (second or third trimester), 0 = non-pregnant.",
        "Time-fixed per subject. Kloprogge 2018 enrolled 3.1% pregnant",
        "women (n = 42 of 1,347 in the PK model-building dataset; Table",
        "1) in their second or third trimester (median gestational age",
        "23.0 weeks, range 13.1-38.0). The paper applies pregnancy as a",
        "proportional relative effect on the absorption rate ka",
        "(Table 2 footnote: 'the categorical pregnancy effect was",
        "implemented as a proportional effect:",
        "ka = theta(n) * (1 + theta_pregnancy)') with theta_pregnancy =",
        "+0.352, i.e. ka in pregnant women is 35.2% higher than in",
        "non-pregnant adults. The +35.2% increase in ka shortens the",
        "absorption phase and is associated with a 20.2% reduction in",
        "day-7 venous plasma lumefantrine concentration (Results section",
        "'Pregnant women') because the distribution kinetics is",
        "substantially altered. Reference category 0 = non-pregnant.",
        sep = " "
      ),
      source_name        = "PREG"
    ),
    PARA = list(
      description        = "Plasmodium falciparum parasitaemia at admission (asexual parasites/uL)",
      units              = "parasites/uL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Admission-only (time-fixed) parasitaemia. Kloprogge 2018",
        "Table 1 'Model building data' row reports median 9,450 (range",
        "13-450,000; interquartile range 2,240-38,700) parasites/uL.",
        "Table 2 footnote: 'All parameters were centred on a non-",
        "pregnant patient weighing 42 kg with an admission parasitaemia",
        "of 15,800 parasites/uL', and the published exponential effect",
        "centering value 4.2 = log10(15,800) is consistent with 15,800",
        "being the geometric-mean (i.e. mean log10) parasitaemia in the",
        "model-building dataset rather than the arithmetic median (the",
        "two diverge because the parasitaemia distribution is right-",
        "skewed). Applied as an exponential",
        "effect on relative bioavailability F using a log10 transform",
        "inside the model:",
        "F_para = exp(e_lnpc_f * (log10(max(PARA, 1)) - 4.2))",
        "with e_lnpc_f = -0.643 per log10 unit. The centering value 4.2",
        "in the published equation matches log10(15,800) = 4.20 (the",
        "model-building median parasitaemia from Table 1, also stated as",
        "the centering reference in the Table 2 footnote); the paper's",
        "wording 'centred on the median natural logarithm transformed",
        "value' is treated as a terminology slip given the unambiguous",
        "numeric match 4.2 = log10(15,800), see vignette Errata. PARA",
        "values below 1 (effectively zero or below detection) are gated",
        "to 1 via max(PARA, 1) so that the exponential effect collapses",
        "to F_para = exp(-0.643 * (0 - 4.2)) = exp(2.7) at PARA = 1; for",
        "downstream simulations where the absence of parasites should",
        "give F_para = 1 (i.e. no covariate effect), set PARA = 15,800",
        "(the centering value) instead of 0. Distinct from the LNPC",
        "canonical (natural-log admission parasitaemia used by Birgersson",
        "2019); PARA carries the raw count and the log10 transform is",
        "performed inside model() to match the source paper's",
        "convention.",
        sep = " "
      ),
      source_name        = "PARA"
    ),
    DOSE = list(
      description        = "Per-dose lumefantrine amount administered (mg)",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Per-dose lumefantrine amount in milligrams, supplied as a",
        "per-dose-record covariate aligned with the corresponding event-",
        "table dosing row (use case (b) of the canonical DOSE entry,",
        "applied per dose record). Required by the model to compute the",
        "per-dose milligram-per-kilogram exposure DOSEMGKG = DOSE / WT,",
        "which then drives the dose-saturable bioavailability term",
        "F_dose = 1 - DOSEMGKG / (dose50 + DOSEMGKG) with dose50 = 3.86",
        "mg/kg (Kloprogge 2018 Table 2). Standard fixed-dose Coartem",
        "(Novartis) provides 120 mg lumefantrine/tablet, so for a 42-kg",
        "adult on the standard 4-tablet/dose regimen, DOSE = 480 mg per",
        "twice-daily dose; pediatric tablet counts scale by weight band",
        "(1 tablet for 5-14 kg, 2 for 15-24 kg, 3 for 25-34 kg, 4 for",
        "ge 35 kg). Set DOSE to the milligram amount administered at",
        "each dose event in the rxode2 event table, alongside amt (mg).",
        sep = " "
      ),
      source_name        = "DOSE"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 1347L,
    n_studies       = 26L,
    n_pregnant      = 42L,
    age_range       = "0.5-78.0 years (Table 1, PK LF model building data)",
    age_median      = "16.6 years (Table 1)",
    weight_range    = "6-150 kg (Table 1, PK LF model building data)",
    weight_median   = "42 kg (Table 1; centering reference in Table 2 footnote)",
    sex_female_pct  = 44.1,
    disease_state   = paste(
      "Uncomplicated Plasmodium falciparum malaria. Pooled WWARN meta-",
      "analysis of 26 published clinical studies across 12 countries",
      "(Benin, Guinea-Bissau, Tanzania, Uganda, Kenya, Mali, Mozambique,",
      "Liberia, Papua New Guinea, Laos, Thailand, Cambodia). Approximately",
      "37% of patients were below 10 years of age and 3.1% were pregnant",
      "(median gestational age 23.0 weeks, range 13.1-38.0).",
      sep = " "
    ),
    dose_range      = paste(
      "Standard fixed-dose Coartem (Novartis): 20 mg artemether + 120 mg",
      "lumefantrine per tablet; weight-band tablet count administered",
      "twice daily for 3 days at 0, 8, 24, 36, 48, and 60 hours, with",
      "fat to optimise lumefantrine bioavailability. Weight-band tablet",
      "counts: 1 tablet (5-14 kg), 2 tablets (15-24 kg), 3 tablets",
      "(25-34 kg), 4 tablets (>=35 kg). Per-dose lumefantrine mg/kg",
      "range 3.2-20.9 mg/kg (median 10.2; Table 1).",
      sep = " "
    ),
    regions         = "Africa (Benin, Guinea-Bissau, Tanzania, Uganda, Kenya, Mali, Mozambique, Liberia), Oceania (Papua New Guinea), Southeast Asia (Laos, Thailand, Cambodia)",
    notes           = paste(
      "Demographics from Kloprogge 2018 Table 1 row 'PK LF model data,",
      "Model building data'. Concentration-time data from patients",
      "contributing 2 or more venous plasma samples per patient were",
      "used to build the PK model; sparser-sampled patients (single",
      "venous-plasma sample, n = 400) were reserved for external",
      "validation. The model-building dataset is multi-matrix-restricted",
      "to venous plasma after intact tablets; sampling-matrix and",
      "formulation correction factors for venous blood (-11.3%),",
      "capillary plasma (-18.3%), capillary blood (-23.8%), dispersible",
      "tablets (+6.31%), and crushed tablets (+26.0%) are reported",
      "post-hoc on the residual-error scale (Methods, Results 'Matrix",
      "and formulation effects') and are NOT encoded here. Eta-",
      "shrinkage 47.5% on Vc/F and 15.2% on F (Table 2 footnote). Of",
      "the 4,122 PK patients in the WWARN repository, 1,347 contributed",
      "to model building, 400 to intact-tablet external validation, 278",
      "(crushed) + 287 (dispersible) to formulation post-hoc, and 595 +",
      "191 + 840 to venous-blood / capillary-plasma / capillary-blood",
      "post-hoc; 154 of the 4,276 were excluded for missing dosing",
      "(n=71) or repeated dosing (n=83).",
      sep = " "
    )
  )

  ini({
    # Structural population mean parameters come from Kloprogge 2018
    # Table 2 'Population estimate' column. The paper reports back-
    # transformed (linear-scale) values; log() is applied here for the
    # nlmixr2 internal scale. Population estimates are centered on a
    # non-pregnant 42-kg patient with admission parasitaemia of 15,800
    # parasites/uL (Table 2 footnote).
    lka  <- log(0.0386) ; label("Absorption rate constant ka (1/h, non-pregnant)")               # Kloprogge 2018 Table 2: ka     = 0.0386 (RSE 2.72%; 95% CI 0.0368-0.0410)
    lcl  <- log(1.35)   ; label("Apparent elimination clearance CL/F (L/h, typical 42-kg)")        # Kloprogge 2018 Table 2: CL/F   = 1.35   (RSE 29.7%; 95% CI 0.538-2.19)
    lvc  <- log(11.2)   ; label("Apparent central volume of distribution Vc/F (L, typical 42-kg)") # Kloprogge 2018 Table 2: Vc/F   = 11.2   (RSE 30.3%; 95% CI 4.65-18.8)
    lq   <- log(0.344)  ; label("Apparent intercompartmental clearance Q/F (L/h, typical 42-kg)")  # Kloprogge 2018 Table 2: Q/F    = 0.344  (RSE 29.8%; 95% CI 0.137-0.566)
    lvp  <- log(59.0)   ; label("Apparent peripheral volume of distribution Vp/F (L, typical 42-kg)") # Kloprogge 2018 Table 2: Vp/F = 59.0   (RSE 29.7%; 95% CI 23.6-96.4)

    # Relative bioavailability anchored at 1 (structural, fixed by the
    # source paper; Table 2 reports F as '1 (fixed)'). All variability
    # in F is captured via the etalfdepot IIV (CV 70.3%) plus the dose-
    # saturable and parasitaemia covariate effects below. The paper
    # additionally reports a Box-Cox shape parameter on F (-0.343,
    # RSE 19.5%) that mildly distorts the F IIV distribution away from
    # strict log-normality; this Box-Cox departure is NOT reproduced in
    # the encoded log-normal IIV here (see vignette Errata).
    lfdepot <- fixed(log(1)) ; label("Relative bioavailability F (unitless, fixed)")             # Kloprogge 2018 Table 2: F = 1 (fixed)

    # Allometric exponents. Fixed at the canonical Mahidol-Oxford
    # malaria-popPK values used by the source paper: 3/4 on clearance
    # (CL/F, Q/F) and 1 on volume (Vc/F, Vp/F). Reported in Table 2
    # footnote: 'Clearance and volume parameters were centred on the
    # median body weight (WT) and scaled allometrically (CL and Q =
    # theta(n) * (WT/42)^(3/4); V = theta(n) * (WT/42))'.
    allo_cl <- fixed(3/4) ; label("Allometric exponent on CL/F and Q/F (unitless, fixed)")        # Kloprogge 2018 Table 2 footnote (allometric scaling)
    allo_vc <- fixed(1)   ; label("Allometric exponent on Vc/F and Vp/F (unitless, fixed)")        # Kloprogge 2018 Table 2 footnote (allometric scaling)

    # Dose-saturable absorption (saturation of relative bioavailability).
    # Table 2 footnote: 'dose-dependent absorption was implemented as a
    # saturation model: F = theta(n) * (1 - Dosage / (Dose50 + Dosage))'.
    # Encoded here so that DOSEMGKG = DOSE / WT (computed inside model())
    # drives F_dose = 1 - DOSEMGKG / (dose50 + DOSEMGKG); at the
    # population-median 10.2 mg/kg dose F_dose = 0.275, at Dose50
    # F_dose = 0.5, and at ~36 mg/kg (~9x Dose50) F_dose ~= 0.10 ('dose90
    # of absorption saturation' per the Results 'Disease- and dosage-
    # related covariate model' paragraph). The 3.86 mg/kg value comes
    # from Table 2; the Results text mentions 3.42 mg/kg from an earlier
    # fit -- the final Table 2 value is used here.
    dose50 <- 3.86 ; label("Dose at 50% saturation of absorption (mg/kg)")                       # Kloprogge 2018 Table 2: Dose50 = 3.86 (RSE 41.5%; 95% CI 1.25-8.04)

    # Pregnancy proportional effect on ka. Table 2 footnote: 'the
    # categorical pregnancy effect was implemented as a proportional
    # effect: ka = theta(n) * (1 + theta_pregnancy)' with
    # theta_pregnancy = +0.352, so ka in second-/third-trimester
    # pregnant women is 35.2% higher than in non-pregnant patients.
    e_preg_ka <- 0.352 ; label("Pregnancy effect on ka: ka_pregnant / ka_nonpregnant - 1 = +0.352") # Kloprogge 2018 Table 2: Pregnancy on ka = 0.352 (RSE 21.2%; 95% CI 0.212-0.510)

    # Parasitaemia exponential effect on F. Results text 'Disease- and
    # dosage-related covariate model':
    #   F = theta(n) * exp(theta_parasitaemia * (parasitaemia - 4.2))
    # with theta_parasitaemia = -0.643 per log10 unit, centered at
    # 4.2 = log10(15,800) (the median admission parasitaemia in the
    # model-building dataset, reported as the centering reference in
    # Table 2 footnote). The covariate name LNPC is used here for the
    # effect parameter even though the active covariate column is the
    # raw PARA count, because the log10 transform is applied inside
    # model() (see model() block). See vignette Errata regarding the
    # paper's 'natural logarithm' wording -- the 4.2 centering value
    # numerically requires log10, not natural log (ln(15,800) = 9.67).
    e_lnpc_f <- -0.643 ; label("Exponential effect of log10 admission parasitaemia on F (per log10 parasites/uL, centered at 4.2)") # Kloprogge 2018 Table 2: Parasitaemia on F = -0.643 (RSE 13.0%; 95% CI -0.793 to -0.473)

    # IIV. Kloprogge 2018 Table 2 reports CV% (footnote: 100 *
    # sqrt(exp(omega^2) - 1)). The internal log-normal variances are
    # recovered by omega^2 = log((CV/100)^2 + 1). Only Vc/F and F carry
    # IIV in the source model; the other structural parameters (ka, CL,
    # Q, Vp) have dashes in the IIV columns of Table 2.
    #   Vc CV 144%  -> log(1.44^2 + 1)  = 1.12286
    #   F  CV 70.3% -> log(0.703^2 + 1) = 0.40158
    etalvc     ~ 1.12286  # Kloprogge 2018 Table 2: BSV on Vc/F = 144% CV (RSE 10.8%; 95% CI 115-165)
    etalfdepot ~ 0.40158  # Kloprogge 2018 Table 2: BSV on F   = 70.3% CV (RSE 5.95%; 95% CI 65.3-75.3). Box-Cox shape parameter -0.343 NOT applied here -- see model description and vignette Errata.

    # Residual error. Kloprogge 2018 modelled the natural logarithm of
    # the lumefantrine plasma concentration with an additive error on
    # the log scale, which is equivalent to proportional error in
    # nlmixr2's linear-concentration space (Methods: 'residual
    # variability was described using an additive error model on
    # logarithmic data'). Table 2 reports the additive log-scale SD as
    # sigma = 0.323; encoded here as propSd.
    propSd <- 0.323 ; label("Proportional residual SD on linear concentration scale (= SD on log scale)") # Kloprogge 2018 Table 2: sigma = 0.323 (RSE 4.88%; 95% CI 0.293-0.357)
  })

  model({
    # Individual structural parameters with allometric WT scaling.
    # Pregnancy applies a proportional shift to ka (no log additivity);
    # the relevant covariates (WT, PREG, PARA, DOSE) are time-fixed
    # except DOSE which is per-dose-record (see covariateData notes).
    ka  <- exp(lka)              * (1 + e_preg_ka * PREG)
    cl  <- exp(lcl)              * (WT / 42)^allo_cl
    vc  <- exp(lvc  + etalvc)    * (WT / 42)^allo_vc
    q   <- exp(lq)               * (WT / 42)^allo_cl
    vp  <- exp(lvp)              * (WT / 42)^allo_vc

    # Two-compartment disposition micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODE system: first-order absorption from depot into a
    # two-compartment disposition model.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                 k12 * central - k21 * peripheral1

    # Relative bioavailability F applied to the depot compartment.
    # Composed of three multiplicative factors:
    #   1. Log-normal IIV: exp(lfdepot + etalfdepot), centered at F = 1.
    #   2. Dose-saturable absorption: F_dose = 1 - DOSEMGKG / (dose50 +
    #      DOSEMGKG), reducing F as the per-dose mg/kg increases.
    #      DOSEMGKG = DOSE / WT (per-dose-event covariate).
    #   3. Parasitaemia effect: F_para = exp(e_lnpc_f *
    #      (log10(max(PARA, 1)) - 4.2)), centered at log10(15,800) = 4.2.
    #      The max(PARA, 1) gate keeps F_para finite for PARA <= 0.
    dose_mgkg <- DOSE / WT
    fdose     <- 1 - dose_mgkg / (dose50 + dose_mgkg)
    fpara     <- exp(e_lnpc_f * (log10(max(PARA, 1)) - 4.2))
    f(depot)  <- exp(lfdepot + etalfdepot) * fdose * fpara

    # Lumefantrine plasma concentration in ug/mL. Dose units are mg,
    # Vc is L, so central/vc has units mg/L = ug/mL. The paper reports
    # day-7 concentrations in ng/mL (e.g. 596 ng/ml target; 175 / 200
    # ng/ml efficacy thresholds); multiply Cc by 1000 in the vignette
    # to compare to those thresholds.
    Cc <- central / vc

    # Proportional residual error on the linear-concentration scale
    # (NONMEM additive-on-log-scale maps to proportional in nlmixr2;
    # see ini() comment on propSd).
    Cc ~ prop(propSd)
  })
}
