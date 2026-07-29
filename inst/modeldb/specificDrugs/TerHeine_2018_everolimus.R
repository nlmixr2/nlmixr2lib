TerHeine_2018_everolimus <- function() {
  description <- paste(
    "Semi-mechanistic two-compartment population PK model for everolimus",
    "in pooled adult oncology (metastatic thyroid or breast cancer) and",
    "renal transplant patients (ter Heine 2018). Oral absorption is",
    "modelled with a chain of four transit compartments parameterised by",
    "the mean absorption time MAT and the Savic 2007 convention",
    "ktr = (n + 1) / MAT (n = 4 transit compartments). Hepatic disposition",
    "uses a well-stirred liver model: hepatic plasma flow QHP = QH *",
    "(1 - HCT); hepatic extraction EH = fu * CLint / (QHP + fu * CLint)",
    "with FIXED unbound fraction fu = 0.27; oral bioavailability",
    "F = 1 - EH and systemic plasma clearance CLH = QHP * EH. Volume",
    "parameters (VC, VP) and flow parameters (QH = 90 L/h FIXED, Q) are",
    "allometrically scaled to fat-free mass FFM at a 57.2 kg reference",
    "(equivalent to a 70 kg, 1.80 m adult male) with theory-based",
    "exponents 0.75 on flows and 1.0 on volumes (Anderson and Holford).",
    "Concomitant high-dose oral prednisolone (PRED_DOSE >= 20 mg/day, a",
    "CYP3A4 inducer) increases apparent CLint by 31%. Modelled plasma",
    "concentrations were derived externally from observed whole-blood",
    "concentrations and HCT via a Langmuir-plus-linear erythrocyte",
    "binding model (Bmax = 0.964 mg/L, Kd = 0.0920 mg/L, Kns = 0.153);",
    "the vignette uses the same back-calculation to compare simulated",
    "plasma concentrations against the paper's whole-blood trough",
    "targets."
  )
  reference <- paste(
    "ter Heine R, van Erp NP, Guchelaar HJ, de Fijter JW, Reinders MEJ,",
    "van Herpen CM, Burger DM, Moes DJAR.",
    "A pharmacological rationale for improved everolimus dosing in",
    "oncology and transplant patients.",
    "Br J Clin Pharmacol. 2018;84(9):1575-1586.",
    "doi:10.1111/bcp.13591.",
    sep = " "
  )
  vignette <- "TerHeine_2018_everolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  covariateData <- list(
    FFM = list(
      description        = paste(
        "Fat-free mass derived from total body weight, height, and sex",
        "via the Janmahasatian (2005) formula as referenced by Holford",
        "et al. (paper reference 36). Time-fixed at baseline."
      ),
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Reference FFM = 57.2 kg, corresponding to a 70 kg, 1.80 m adult",
        "male (BMI = 21.6, male Janmahasatian formula yields exactly 57.2",
        "kg). Per Methods 'Structural model development', all flow",
        "parameters (QH = 90 L/h FIXED, intercompartmental Q) and volume",
        "parameters (VC, VP) are allometrically scaled to FFM using",
        "theory-based exponents 0.75 on flows and 1.0 on volumes",
        "(Anderson and Holford). The paper compared total body weight,",
        "normal fat mass, and FFM as size descriptors and found FFM best",
        "described the PK / body-size relationship (Results 'Base model",
        "development')."
      ),
      source_name        = "Fat-free mass (derived from weight, length, sex per Janmahasatian / Holford et al.)"
    ),
    HCT = list(
      description        = paste(
        "Hematocrit expressed as a volume fraction (0 - 1). Enters the",
        "hepatic plasma flow via QHP = QH * (1 - HCT) in the well-stirred",
        "liver model and the whole-blood-to-plasma back-conversion in the",
        "vignette."
      ),
      units              = "fraction (0 - 1)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Hematocrit drives two distinct effects in the source paper: (1)",
        "hepatic plasma flow via QHP = QH * (1 - HCT) inside the",
        "well-stirred liver model (Methods Equation 3), and (2) the",
        "whole-blood-to-plasma concentration conversion used by the",
        "authors to derive the modelled plasma concentrations from the",
        "assayed whole-blood concentrations (Methods 'Structural model",
        "development', binding equation; Bmax = 0.964 mg/L, Kd = 0.0920",
        "mg/L, Kns = 0.153). Cohort range 0.28 - 0.50 (Table 1; transplant",
        "median 0.36, cancer median 0.38)."
      ),
      source_name        = "Ht (fraction)"
    ),
    PRED_DOSE = list(
      description        = paste(
        "Concomitant oral prednisolone daily dose (mg/day). Used as a",
        "threshold-form covariate at >= 20 mg/day (binary high-dose",
        "indicator) on apparent intrinsic clearance."
      ),
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Methods 'Covariate analysis' tested concomitant high-dose",
        "prednisolone (defined as >= 20 mg/day, a known CYP3A4 inducer)",
        "as a binary covariate on apparent intrinsic clearance.",
        "Final-model effect (Table 2): +31% on CLint when PRED_DOSE >=",
        "20 mg/day (95% CI 9.10 - 55.8%). The binary indicator is",
        "derived inside model() from the continuous PRED_DOSE column.",
        "Oncology cohort: 0 mg/day in all 71 subjects (no concomitant",
        "prednisolone). Transplant cohort: 7.5 - 40 mg/day (mean 14.5,",
        "median 10; Table 1)."
      ),
      source_name        = "Prednisolone dose (total daily dose, mg/day)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 126L,
    n_studies      = 5L,
    age_range      = "19 - 80 years",
    age_median     = "53 years (transplant), 62 years (cancer)",
    weight_range   = "45 - 110.3 kg",
    weight_median  = "79 kg (transplant), 73 kg (cancer)",
    sex_female_pct = 54.0,
    race_ethnicity = NULL,
    disease_state  = paste(
      "Two pooled adult subpopulations: (1) 71 oncology patients with",
      "metastatic thyroid or breast cancer on everolimus (Afinitor) 10",
      "mg once daily (clinicaltrials.gov NCT01118065 and NCT01948960);",
      "(2) 55 renal transplant recipients on everolimus (Certican) 0.75 -",
      "3 mg twice daily for prophylaxis of allograft rejection (Dutch",
      "Trial Register NTR567 and NTR1615; clinicaltrials.gov NCT02387151)."
    ),
    dose_range     = paste(
      "Oncology cohort: everolimus 10 mg orally once daily (Afinitor;",
      "all 71 subjects).",
      "Transplant cohort: everolimus 1.5 - 3 mg orally twice daily",
      "(Certican; mean 2.4 mg, median 3 mg, range 1.5 - 3 mg per",
      "Table 1).",
      "Twenty-two of the transplant subjects also received tacrolimus",
      "after alemtuzumab induction; the remainder were on",
      "everolimus + prednisolone."
    ),
    regions        = "The Netherlands",
    notes          = paste(
      "Rich PK sampling: approximately 8 or more samples per dosing",
      "interval (transplant; mean 6.3, median 7 samples per patient) and",
      "12.6 mean samples per cancer patient (range 3 - 21).",
      "Total 1239 PK observations across 126 subjects (893 oncology +",
      "347 transplant samples per Table 1).",
      "NONMEM v7.3 with Stochastic Approximation Expectation",
      "Maximization (SAEM) followed by Importance Sampling (SAEM-IMP);",
      "95% CIs computed by Sampling Importance Resampling in PsN 4.7.0.",
      "Whole-blood everolimus concentrations were converted to plasma",
      "concentrations before fitting using paired whole-blood and HCT",
      "observations and a Langmuir-plus-linear erythrocyte-binding model",
      "(Bmax = 0.964 mg/L, Kd = 0.0920 mg/L, Kns = 0.153) fit externally",
      "to in vitro [3H]everolimus blood distribution data from the",
      "manufacturer (Methods 'Structural model development', Equation 1).",
      "All reported parameter estimates are plasma PK parameters and",
      "apparent (oral / F) -- absolute bioavailability of everolimus is",
      "unknown. NONMEM control stream is provided as supplemental S2 in",
      "the source publication."
    )
  )

  ini({
    # ============================================================
    # Structural PK parameters -- ter Heine 2018 Table 2 'Final model'
    # column. All volume and flow parameters are allometrically scaled
    # to a reference fat-free mass of 57.2 kg (Methods 'Structural model
    # development'; theory-based Anderson and Holford exponents 0.75 on
    # flows and 1.0 on volumes -- the paper states the parameters were
    # allometrically scaled to FFM but does not restate the exponents,
    # so the canonical theory-based values are used here and the
    # assumption is documented in the vignette's Assumptions section).
    # ============================================================
    lmat <- log(0.404)
    label("Mean absorption time MAT (h)")
    # Table 2 Final: MAT = 0.404 h (95% CI 0.347 - 0.457)
    lcl_int <- log(340)
    label("Apparent intrinsic clearance CLint (L/h)")
    # Table 2 Final: CLint = 340 L/h (95% CI 319 - 362)
    lvc <- log(175)
    label("Apparent central volume VC at FFM = 57.2 kg (L)")
    # Table 2 Final: VC = 175 L (95% CI 157 - 193)
    lvp <- log(577)
    label("Apparent peripheral volume VP at FFM = 57.2 kg (L)")
    # Table 2 Final: VP = 577 L (95% CI 542 - 609)
    lq <- log(85.7)
    label("Apparent inter-compartmental clearance Q at FFM = 57.2 kg (L/h)")
    # Table 2 Final: Q = 85.7 L/h (95% CI 80.3 - 91.5)
    lqh <- fixed(log(90))
    label("Hepatic blood flow QH at FFM = 57.2 kg (L/h; FIXED)")
    # Table 2 Final: QH = 90 L/h (FIXED; Methods 'Structural model development')

    # Reference unbound fraction (Methods 'Structural model development',
    # Equation 4). fu = 0.27, FIXED from the literature (paper reference
    # 35). The paper notes an interindividual variability of 3% on fu
    # but this is not retained in the typical-value model encoded here.
    lfu <- fixed(log(0.27))
    label("Unbound fraction fu (unitless; FIXED)")

    # Allometric exponents fixed at theory-based Anderson and Holford
    # values (Methods 'Structural model development' cites the Holford
    # et al. paper for the allometric form but does not restate the
    # exponents; canonical 0.75 on flows and 1.0 on volumes is the
    # convention used in the well-stirred liver / FFM scaling literature
    # the paper references).
    e_ffm_flow <- fixed(0.75)
    label("Allometric exponent on Q, QH (unitless; FIXED)")
    e_ffm_volume <- fixed(1.0)
    label("Allometric exponent on VC, VP (unitless; FIXED)")

    # Prednisolone covariate effect on CLint (Methods 'Covariate
    # analysis'; Table 2 Final 'Increase in CLint due high dose
    # prednisolone'). Applied multiplicatively as
    # (1 + e_pred_dose_high_cl_int * pred_high) where pred_high = 1
    # when PRED_DOSE >= 20 mg/day, else 0.
    e_pred_dose_high_cl_int <- 0.31
    label("Fractional increase in CLint for PRED_DOSE >= 20 mg/day (unitless)")
    # Table 2 Final: +31% (95% CI 9.10 - 55.8)

    # ============================================================
    # Inter-individual variability -- ter Heine 2018 Table 2 'Final
    # model'. The paper reports CV%; converted to log-normal omega^2 =
    # log(1 + CV^2).
    # ============================================================
    etalcl_int ~ log(1 + 0.339^2)
    # Table 2 Final IIV CLint = 33.9% CV (95% CI 29.4 - 38.2%)
    etalvc ~ log(1 + 0.406^2)
    # Table 2 Final IIV VC = 40.6% CV (95% CI 29.3 - 50.0%)
    etalmat ~ log(1 + 1.10^2)
    # Table 2 Final intra-individual (inter-occasion) variability on
    # MAT = 110% CV (95% CI 90.9 - 127%). The source describes this as
    # an intra-individual / IOV term across dosing occasions; we encode
    # it as IIV on lmat here because the source does not specify the
    # per-subject number of dosing occasions explicitly. The simulated
    # marginal magnitude of MAT variability is preserved; the inter-
    # occasion structure is documented as a deviation in the vignette
    # Assumptions and deviations section.

    # ============================================================
    # Residual error -- Table 2 Final 'Residual variability', CV 17.9%.
    # Methods 'Statistics' state 'The residual error was modelled with
    # an additive error on the log scale, thus approximating a
    # proportional error' -- this maps directly onto nlmixr2's prop()
    # error on the linear scale.
    # ============================================================
    propSd <- 0.179
    label("Proportional residual error (fraction)")
    # Table 2 Final residual variability = 17.9% (95% CI 17.0 - 18.7%)
  })

  model({
    # ------------------------------------------------------------
    # Individual structural parameters with FFM-based allometric
    # scaling (Methods 'Structural model development'). Reference FFM
    # = 57.2 kg corresponds to a 70 kg, 1.80 m adult male via the
    # Janmahasatian male formula.
    # ------------------------------------------------------------
    pred_high <- (PRED_DOSE >= 20)
    mat <- exp(lmat + etalmat)
    cl_int <- exp(lcl_int + etalcl_int) *
              (1 + e_pred_dose_high_cl_int * pred_high)
    vc <- exp(lvc + etalvc) * (FFM / 57.2)^e_ffm_volume
    vp <- exp(lvp) * (FFM / 57.2)^e_ffm_volume
    q  <- exp(lq) * (FFM / 57.2)^e_ffm_flow
    qh <- exp(lqh) * (FFM / 57.2)^e_ffm_flow
    fu <- exp(lfu)

    # Transit-compartment rate constant (Methods Equation 2). The
    # source paper states n = 4 transit compartments; with the Savic
    # 2007 convention ktr = (n + 1) / MAT, four intermediate transit
    # compartments give ktr = 5 / MAT.
    ktr <- 5 / mat

    # Well-stirred liver model (Methods Equations 3 - 5):
    #   QHP = QH * (1 - HCT)                    hepatic plasma flow
    #   EH  = fu * CLint / (QHP + fu * CLint)   hepatic extraction
    #   CLH = QHP * EH                          systemic plasma CL
    #   F   = 1 - EH                            oral bioavailability
    qhp <- qh * (1 - HCT)
    eh  <- fu * cl_int / (qhp + fu * cl_int)
    cl_h <- qhp * eh
    fdepot <- 1 - eh

    # Micro-constants
    kel <- cl_h / vc
    k12 <- q / vc
    k21 <- q / vp

    # ------------------------------------------------------------
    # Four transit compartments feeding into a two-compartment
    # disposition model (Methods 'Structural model development',
    # Figure 1). The chain depot -> transit1 -> transit2 -> transit3
    # -> transit4 -> central propagates at the common rate ktr.
    # ------------------------------------------------------------
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot    - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1 - ktr * transit2
    d/dt(transit3)    <-  ktr * transit2 - ktr * transit3
    d/dt(transit4)    <-  ktr * transit3 - ktr * transit4
    d/dt(central)     <-  ktr * transit4 - kel * central -
                          k12 * central  + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central  - k21 * peripheral1

    # Oral bioavailability via hepatic first pass (well-stirred liver).
    f(depot) <- fdepot

    # ------------------------------------------------------------
    # Observation: plasma concentration. Dose units mg; vc in L means
    # central amount is in mg, so central / vc is in mg/L. Multiply by
    # 1000 to convert to ug/L (the paper's reporting unit).
    # ------------------------------------------------------------
    Cc <- 1000 * central / vc

    # ------------------------------------------------------------
    # Whole-blood concentration via the source paper's Langmuir-plus-
    # linear erythrocyte binding model (Methods 'Structural model
    # development' Equation 1; Results 'Base model development'):
    #   Crb = Bmax * Cp / (Kd + Cp) + Kns * Cp
    #   Cwb = HCT * Crb + (1 - HCT) * Cp
    # Bmax = 964 ug/L, Kd = 92 ug/L, Kns = 0.153 (unitless partition
    # coefficient; the paper text reads 'Bmax, Kd and Kns ... 0.964,
    # 0.0920 and 0.153 mg l-1, respectively' but Equation 1
    # dimensional analysis requires Kns to be dimensionless so the
    # 'mg l-1' qualifier applies to Bmax and Kd only). These constants
    # are not adjustable parameters of the PK fit; they were fit
    # externally to in vitro [3H]everolimus blood distribution data
    # from the manufacturer.
    # ------------------------------------------------------------
    Crb <- 964 * Cc / (92 + Cc) + 0.153 * Cc
    Cwb <- HCT * Crb + (1 - HCT) * Cc

    Cc ~ prop(propSd)
  })
}
