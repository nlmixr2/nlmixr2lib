Petain_2008_imatinib <- function() {
  description <- paste0(
    "Joint parent-metabolite population PK model for oral imatinib and its ",
    "main active metabolite N-desmethyl imatinib (CGP74588) in children ",
    "with solid malignancies and adults with gastrointestinal stromal ",
    "tumor, aged 2-84 years (Petain 2008). Imatinib is described by a ",
    "one-compartment model with first-order absorption and first-order ",
    "elimination; the metabolite is described by a one-compartment model ",
    "whose disposition parameters are apparent with respect to the ",
    "fraction metabolized (CLm/fm and V1m/fm). Imatinib CL/F carries ",
    "power-form effects of total body weight (reference 54 kg), alpha-1 ",
    "acid glycoprotein (reference 1.13 g/L) and albumin (reference 38 ",
    "g/L); imatinib Vc/F carries weight and AAG effects; metabolite CLm/fm ",
    "carries weight and AAG effects plus a 0.7-fold shift on the ",
    "day-29-or-later occasion. This is one of only two models in the ",
    "source evaluation developed on data from both children and adults. ",
    "TRANSCRIBED FROM A SECONDARY SOURCE: the parameter values come from ",
    "Table 1 of Yang 2025, an external evaluation of 15 published imatinib ",
    "population PK models, not from the primary publication. Re-extract ",
    "from Petain 2008 when that paper is obtained."
  )
  reference <- paste0(
    "Petain A, Kattygnarath D, Azard J, Chatelut E, Delbaldo C, Geoerger ",
    "B, Barrois M, Seronie-Vivien S, LeCesne A, Vassal G. Population ",
    "pharmacokinetics and pharmacogenetics of imatinib in children and ",
    "adults. Clin Cancer Res. 2008;14(21):7102-7109. ",
    "doi:10.1158/1078-0432.CCR-08-0950. ",
    "PARAMETER SOURCE (secondary): Yang T, Rasmussen ASB, Weimann A, ",
    "Thastrup M, Rank CU, Als-Nielsen B, Malmros J, Wik HS, Lohi O, ",
    "Overgaard U, Johannsdottir IMR, Vaitkeviciene G, Dalhoff K, ",
    "Schmiegelow K, Lund TM. Published population pharmacokinetic models ",
    "of imatinib perform poorly on TDM data from pediatric patients. ",
    "Target Oncol. 2025;20(5):871-886. Table 1, row 'Petain et al. (2008)' ",
    "and Table 1 footnote b. doi:10.1007/s11523-025-01172-2."
  )
  vignette <- "Yang_2025_imatinib_external_evaluation"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot          = list(analyte = "imatinib", units = "mg", specimen = "administration site", verified = FALSE),
    central        = list(analyte = "imatinib", units = "mg", specimen = "plasma", verified = FALSE),
    central_ndmima = list(
      analyte  = "N-desmethyl imatinib (CGP74588)",
      units    = "mg divided by the fraction of imatinib metabolized to N-desmethyl imatinib (fm)",
      specimen = "plasma",
      verified = FALSE
    )
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Enters imatinib CL/F as (WT/54)^0.56, imatinib Vc/F as ",
        "(WT/54)^0.79 and metabolite CLm/fm as (WT/54)^-0.62 (Yang 2025 ",
        "Table 1). The reference 54 kg is the centring constant printed ",
        "inside every covariate term and is plausible as the cohort median ",
        "for a pooled 2-84 year cohort of children and adults. Note the ",
        "NEGATIVE exponent on metabolite clearance, which is unusual and ",
        "should be confirmed against the primary. Because this model was ",
        "developed on both children and adults and already carries a body ",
        "weight effect, Yang 2025 did NOT apply its own allometric scaling ",
        "to it (Yang 2025 Results: 'Allometric scaling was not used in the ",
        "models by Menon-Andersen et al. and Petain et al., as their models ",
        "were based on a population of both adults and children, and they ",
        "used body weight scaling in their original models')."
      ),
      source_name        = "TBW"
    ),
    AAG = list(
      description        = "Alpha-1 acid glycoprotein concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Enters imatinib CL/F as (AAG/1.13)^-0.65, imatinib Vc/F as ",
        "(AAG/1.13)^-1.01 and metabolite CLm/fm as (AAG/1.13)^-0.81 (Yang ",
        "2025 Table 1). Yang 2025 Table 1 abbreviation list gives the unit ",
        "explicitly: 'AGP plasma alpha1-acid glycoprotein level (g/L)', ",
        "matching the canonical AAG unit, so no conversion is required. ",
        "The reference 1.13 g/L is the centring constant printed inside ",
        "the covariate terms. The negative exponents are mechanistically ",
        "expected: imatinib is extensively bound to AAG, so a higher AAG ",
        "concentration lowers the unbound fraction and therefore both the ",
        "apparent clearance and the apparent volume measured on total ",
        "plasma concentrations. The exponent on Vc/F is close to -1, which ",
        "is what a pure binding-driven effect on a total-concentration ",
        "volume would predict."
      ),
      source_name        = "AGP"
    ),
    ALB = list(
      description        = "Serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Enters imatinib CL/F as (ALB/38)^0.66 (Yang 2025 Table 1). Yang ",
        "2025 Table 1 abbreviation list: 'ALB albumin (g/L)', matching the ",
        "canonical ALB unit. The reference 38 g/L is the centring constant ",
        "printed inside the covariate term. Albumin appears on imatinib ",
        "CL/F only, not on Vc/F and not on the metabolite."
      ),
      source_name        = "ALB"
    ),
    OCC = list(
      description        = "Occasion indicator: 1 = the day-1 (first-dose) occasion, 2 = the day-29-or-later occasion",
      units              = "(count)",
      type               = "categorical",
      reference_category = "1 (day 1, first dose)",
      notes              = paste0(
        "Yang 2025 Table 1 attaches footnote b to this model's metabolite ",
        "clearance term 'CLm/fm: 52.2 x ... x 0.7^occ'. Footnote b reads ",
        "'Occasion (OCC) with OCC = 0 for day 1 and OCC = 1 (used in this ",
        "current validation study) for day >= 29'. The canonical OCC ",
        "column is 1-based, so the paper's OCC = 0 maps to OCC = 1 here ",
        "and the paper's OCC = 1 maps to OCC = 2. The effect is ",
        "MULTIPLICATIVE and applies to the METABOLITE clearance only: a ",
        "0.7-fold (30% lower) CLm/fm on the day-29-or-later occasion. No ",
        "eta is attached to the occasion, so this is a structural ",
        "typical-value shift, not inter-occasion variability. Same column ",
        "and coding as Schmidli_2005_imatinib.R and Demetri_2009_imatinib.R."
      ),
      source_name        = "occ"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 67L,
    n_studies      = 1L,
    n_observations = "505 concentrations (Yang 2025 Table 1)",
    age_range      = "2-84 years",
    disease_state  = paste0(
      "Children with solid malignancies and adults with gastrointestinal ",
      "stromal tumor (GIST). Yang 2025 Results notes that the primary ",
      "'included 33 patients with solid malignancies aged 2-22 years' in a ",
      "phase II study."
    ),
    dose_range     = "Oral imatinib 340 mg/m^2 daily (children); 400 mg daily (adults)",
    regions        = "France",
    bioanalytical  = "HPLC; limit of quantification 10 ng/mL for imatinib and 10 ng/mL for N-desmethyl imatinib (Yang 2025 Table 1)",
    notes          = paste0(
      "One of only two models among the 15 evaluated by Yang 2025 that was ",
      "developed on data from both children and adults (the other is ",
      "MenonAndersen_2009_imatinib.R). Despite this, Yang 2025 found that ",
      "it 'did not show superior performance compared to the other ",
      "models' on their pediatric-inclusive external dataset (Yang 2025 ",
      "Results 3.3): median prediction error 15.03%, median absolute ",
      "prediction error 53.83% for imatinib (Table 3). Demographic detail ",
      "beyond the row above (weight range, sex split, race) is not ",
      "reported by the secondary source and must be read from the primary ",
      "publication."
    )
  )

  ini({
    # ----- Imatinib (parent) structural parameters -----
    # Typical values are for the REFERENCE subject: WT = 54 kg,
    # AAG = 1.13 g/L, ALB = 38 g/L, on the day-1 occasion.
    lka <- log(1.03); label("Imatinib first-order absorption rate constant ka (1/h)")  # Yang 2025 Table 1: Ka = 1.03
    lcl <- log(7.29); label("Imatinib apparent oral clearance CL/F at the reference subject (L/h)")  # Yang 2025 Table 1: CL/F = 7.29 x (TBW/54)^0.56 x (AGP/1.13)^-0.65 x (ALB/38)^0.66
    lvc <- log(202); label("Imatinib apparent central volume Vc/F at the reference subject (L)")  # Yang 2025 Table 1: Vc/F = 202 x (TBW/54)^0.79 x (AGP/1.13)^-1.01

    # ----- Imatinib covariate exponents -----
    e_wt_cl <- 0.56; label("Power exponent of (WT / 54) on imatinib CL/F (unitless)")  # Yang 2025 Table 1
    e_aag_cl <- -0.65; label("Power exponent of (AAG / 1.13) on imatinib CL/F (unitless)")  # Yang 2025 Table 1
    e_alb_cl <- 0.66; label("Power exponent of (ALB / 38) on imatinib CL/F (unitless)")  # Yang 2025 Table 1
    e_wt_vc <- 0.79; label("Power exponent of (WT / 54) on imatinib Vc/F (unitless)")  # Yang 2025 Table 1
    e_aag_vc <- -1.01; label("Power exponent of (AAG / 1.13) on imatinib Vc/F (unitless)")  # Yang 2025 Table 1

    # ----- N-desmethyl imatinib (metabolite) structural parameters -----
    # Both metabolite disposition parameters are APPARENT with respect to
    # the fraction of imatinib metabolized to N-desmethyl imatinib (fm),
    # which is not separately identifiable from plasma data alone. Yang
    # 2025 Table 1 labels them 'CLm/fm' and 'V1m/fm' and its abbreviation
    # list defines them as the apparent clearance and apparent central
    # volume of the metabolite. Because BOTH are divided by the same fm,
    # the ratio CLm/V1m -- and therefore the predicted metabolite
    # concentration -- is unaffected; see the ODE comment in model().
    lcl_ndmima <- log(52.2); label("N-desmethyl imatinib apparent clearance CLm/fm at the reference subject on the day-1 occasion (L/h)")  # Yang 2025 Table 1: CLm/fm = 52.2 x (TBW/54)^-0.62 x (AGP/1.13)^-0.81 x 0.7^occ
    lvc_ndmima <- log(32.2); label("N-desmethyl imatinib apparent central volume V1m/fm (L)")  # Yang 2025 Table 1: V1m/fm = 32.2

    # ----- Metabolite covariate effects -----
    e_wt_cl_ndmima <- -0.62; label("Power exponent of (WT / 54) on N-desmethyl imatinib CLm/fm (unitless)")  # Yang 2025 Table 1
    e_aag_cl_ndmima <- -0.81; label("Power exponent of (AAG / 1.13) on N-desmethyl imatinib CLm/fm (unitless)")  # Yang 2025 Table 1
    # Multiplicative factor raised to the occasion indicator, exactly as the
    # published term 0.7^occ is written.
    e_occ_cl_ndmima <- 0.7; label("Multiplicative factor on N-desmethyl imatinib CLm/fm for the day-29-or-later occasion (unitless)")  # Yang 2025 Table 1: x 0.7^occ

    # ----- Inter-individual variability -----
    # Yang 2025 Table 1 reports, for this row, 'Imatinib: CV%(CL): 19%,
    # CV%(Ka): 82%' and 'Metabolite: CV%(CLm): 28%, CV%(V1m): 81%'. No
    # covariances are tabulated, so all four etas are carried as
    # independent diagonal elements. The tabulated CV% is taken as omega
    # (the log-scale standard deviation), so the variance is (CV/100)^2 --
    # the convention used throughout this transcription and in the shipped
    # Jiang_2023_imatinib.R. Note that NO IIV is reported on imatinib Vc/F.
    etalcl ~ 0.0361  # Yang 2025 Table 1: imatinib CV%(CL) 19% -> omega^2 = 0.19^2
    etalka ~ 0.6724  # Yang 2025 Table 1: imatinib CV%(Ka) 82% -> omega^2 = 0.82^2
    etalcl_ndmima ~ 0.0784  # Yang 2025 Table 1: metabolite CV%(CLm) 28% -> omega^2 = 0.28^2
    etalvc_ndmima ~ 0.6561  # Yang 2025 Table 1: metabolite CV%(V1m) 81% -> omega^2 = 0.81^2

    # ----- Residual unexplained variability -----
    # Separate proportional error models for the two analytes.
    propSd <- 0.336; label("Imatinib proportional residual error (fraction)")  # Yang 2025 Table 1: imatinib Prop 33.6%
    propSd_ndmima <- 0.348; label("N-desmethyl imatinib proportional residual error (fraction)")  # Yang 2025 Table 1: metabolite Prop 34.8%
  })

  model({
    # ----- 1. Occasion indicator -----
    # Canonical OCC is 1-based; occ_late is 1 on the day-29-or-later
    # occasion and 0 on the day-1 occasion, matching the paper's own 0/1
    # occ coding term for term.
    occ_late <- (OCC == 2)

    # ----- 2. Individual parameters -----
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) *
      (WT / 54)^e_wt_cl * (AAG / 1.13)^e_aag_cl * (ALB / 38)^e_alb_cl
    vc <- exp(lvc) * (WT / 54)^e_wt_vc * (AAG / 1.13)^e_aag_vc

    cl_ndmima <- exp(lcl_ndmima + etalcl_ndmima) *
      (WT / 54)^e_wt_cl_ndmima * (AAG / 1.13)^e_aag_cl_ndmima *
      e_occ_cl_ndmima^occ_late
    vc_ndmima <- exp(lvc_ndmima + etalvc_ndmima)

    # ----- 3. ODE system -----
    # Imatinib: one compartment, first-order absorption from depot.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - cl * central / vc

    # N-desmethyl imatinib: one compartment, formed from imatinib.
    # The metabolite state is carried in fm-SCALED amount units, i.e. the
    # state holds A_m / fm rather than A_m. Writing the true metabolite
    # balance as
    #     dA_m/dt = fm * CL * Cc - CLm * A_m / Vm
    # and dividing throughout by fm gives
    #     d(A_m/fm)/dt = CL * Cc - (CLm/fm) * (A_m/fm) / (Vm/fm),
    # which is the equation below with cl_ndmima = CLm/fm and
    # vc_ndmima = Vm/fm. The observed metabolite concentration is
    # unaffected by the scaling because (A_m/fm) / (Vm/fm) = A_m / Vm, so
    # Cc_ndmima below is the true plasma concentration even though fm is
    # not identifiable. The formation flux is the FULL imatinib
    # elimination flux cl * central / vc: any fraction of imatinib
    # clearance that does not form N-desmethyl imatinib is absorbed into
    # the unidentifiable fm and cancels.
    d/dt(central_ndmima) <-
      cl * central / vc - cl_ndmima * central_ndmima / vc_ndmima

    # ----- 4. Observations and error -----
    # States are mg (imatinib) or mg/fm (metabolite) and volumes are L, so
    # each ratio is mg/L; x 1000 gives ng/mL, the unit of the tabulated
    # limits of quantification and of imatinib TDM targets.
    Cc <- 1000 * central / vc
    Cc_ndmima <- 1000 * central_ndmima / vc_ndmima

    Cc ~ prop(propSd)
    Cc_ndmima ~ prop(propSd_ndmima)
  })
}
