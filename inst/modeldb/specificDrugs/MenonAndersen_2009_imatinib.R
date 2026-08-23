MenonAndersen_2009_imatinib <- function() {
  description <- paste0(
    "Joint parent-metabolite population PK model for oral imatinib and its ",
    "main active metabolite N-desmethyl imatinib (CGP74588) in children ",
    "and young adults aged 6-24 years with Philadelphia chromosome-positive ",
    "leukemia or solid tumors (Menon-Andersen 2009). Imatinib is described ",
    "by a one-compartment model with zero-order absorption of duration D1 ",
    "into the central compartment; the metabolite is described by a ",
    "TWO-compartment model whose disposition parameters are apparent with ",
    "respect to the fraction metabolized (CLm/fm, V1m/fm, Qm/fm, V2m/fm). ",
    "Every disposition parameter of both analytes is scaled allometrically ",
    "by total body weight referenced to 70 kg, with the standard fixed ",
    "exponents 0.75 on clearances and 1 on volumes. Residual error is ",
    "log-scale additive (equivalent to proportional in linear space) for ",
    "both analytes. TRANSCRIBED FROM A SECONDARY SOURCE: the parameter ",
    "values come from Table 1 of Yang 2025, an external evaluation of 15 ",
    "published imatinib population PK models, not from the primary ",
    "publication. Re-extract from Menon-Andersen 2009 when that paper is ",
    "obtained."
  )
  reference <- paste0(
    "Menon-Andersen D, Mondick JT, Jayaraman B, Thompson PA, Blaney SM, ",
    "Bernstein M, Bond M, Champagne M, Fossler MJ, Barrett JS. Population ",
    "pharmacokinetics of imatinib mesylate and its metabolite in children ",
    "and young adults. Cancer Chemother Pharmacol. 2009;63(2):229-238. ",
    "doi:10.1007/s00280-008-0730-x. ",
    "PARAMETER SOURCE (secondary): Yang T, Rasmussen ASB, Weimann A, ",
    "Thastrup M, Rank CU, Als-Nielsen B, Malmros J, Wik HS, Lohi O, ",
    "Overgaard U, Johannsdottir IMR, Vaitkeviciene G, Dalhoff K, ",
    "Schmiegelow K, Lund TM. Published population pharmacokinetic models ",
    "of imatinib perform poorly on TDM data from pediatric patients. ",
    "Target Oncol. 2025;20(5):871-886. Table 1, row 'Menon-Andersen et al. ",
    "(2009)'. doi:10.1007/s11523-025-01172-2."
  )
  vignette <- "Yang_2025_imatinib_external_evaluation"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    central            = list(analyte = "imatinib", units = "mg", specimen = "plasma", verified = FALSE),
    central_ndmima     = list(
      analyte  = "N-desmethyl imatinib (CGP74588)",
      units    = "mg divided by the fraction of imatinib metabolized to N-desmethyl imatinib (fm)",
      specimen = "plasma",
      verified = FALSE
    ),
    peripheral1_ndmima = list(
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
        "Scales every disposition parameter of both analytes ",
        "allometrically, referenced to 70 kg (Yang 2025 Table 1): imatinib ",
        "CL/F x (WT/70)^0.75 and Vc/F x (WT/70)^1; metabolite CLm/fm x ",
        "(WT/70)^0.75, V1m/fm x (WT/70)^1, Qm/fm x (WT/70)^0.75 and ",
        "V2m/fm x (WT/70)^1. The exponents are the standard allometric ",
        "pair (0.75 on flows, 1 on volumes) and are encoded as fixed(), ",
        "since a table that prints exactly 0.75 and an implicit 1 with no ",
        "uncertainty is reporting theory-imposed constants rather than ",
        "estimates. Because this model was developed on children and young ",
        "adults and already carries weight-based scaling, Yang 2025 did NOT ",
        "apply its own allometric scaling to it (Yang 2025 Results: ",
        "'Allometric scaling was not used in the models by Menon-Andersen ",
        "et al. and Petain et al., as their models were based on a ",
        "population of both adults and children, and they used body weight ",
        "scaling in their original models')."
      ),
      source_name        = "TBW"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 41L,
    n_studies      = 1L,
    n_observations = "842 imatinib and 424 N-desmethyl imatinib plasma concentrations (Yang 2025 Table 1)",
    age_range      = "6-24 years",
    disease_state  = "Children and young adults with Philadelphia chromosome-positive (Ph+) leukemia or solid tumors",
    dose_range     = "Oral imatinib 260-570 mg/m^2 total daily dose",
    regions        = "USA",
    bioanalytical  = "HPLC; limit of quantification 1 or 10 ng/mL for imatinib and 2 ng/mL for N-desmethyl imatinib (Yang 2025 Table 1)",
    notes          = paste0(
      "One of only two models among the 15 evaluated by Yang 2025 that was ",
      "developed on data including children (the other is ",
      "Petain_2008_imatinib.R), and the only one developed exclusively in a ",
      "pediatric / young-adult cohort. Despite this, Yang 2025 found that ",
      "it 'did not show superior performance compared to the other models' ",
      "on their pediatric-inclusive external dataset (Yang 2025 Results ",
      "3.3): median prediction error -28.75%, median absolute prediction ",
      "error 41.06% for imatinib (Table 3). Demographic detail beyond the ",
      "row above (weight range, sex split, race) is not reported by the ",
      "secondary source and must be read from the primary publication."
    )
  )

  ini({
    # ----- Imatinib (parent) structural parameters -----
    # Typical values are for a 70 kg subject, at which every allometric
    # factor equals 1.
    lcl <- log(10.8); label("Imatinib apparent oral clearance CL/F at WT = 70 kg (L/h)")  # Yang 2025 Table 1: CL/F = 10.8 x (TBW/70)^0.75
    lvc <- log(284); label("Imatinib apparent central volume Vc/F at WT = 70 kg (L)")  # Yang 2025 Table 1: Vc/F = 284 x (TBW/70)
    ld1 <- log(1.67); label("Zero-order absorption duration D1 (h)")  # Yang 2025 Table 1: D1 = 1.67

    # ----- N-desmethyl imatinib (metabolite) structural parameters -----
    # All four are APPARENT with respect to the fraction of imatinib
    # metabolized to N-desmethyl imatinib (fm), which is not separately
    # identifiable from plasma data. Because the SAME fm divides every
    # metabolite clearance and every metabolite volume, the predicted
    # metabolite concentration is unaffected; see the ODE comment in
    # model().
    lcl_ndmima <- log(9.65); label("N-desmethyl imatinib apparent clearance CLm/fm at WT = 70 kg (L/h)")  # Yang 2025 Table 1: CLm/fm = 9.65 x (TBW/70)^0.75
    lvc_ndmima <- log(11.6); label("N-desmethyl imatinib apparent central volume V1m/fm at WT = 70 kg (L)")  # Yang 2025 Table 1: V1m/fm = 11.6 x (TBW/70)
    lq_ndmima <- log(2.9); label("N-desmethyl imatinib apparent intercompartmental clearance Qm/fm at WT = 70 kg (L/h)")  # Yang 2025 Table 1: Qm/fm = 2.9 x (TBW/70)^0.75
    lvp_ndmima <- log(256); label("N-desmethyl imatinib apparent peripheral volume V2m/fm at WT = 70 kg (L)")  # Yang 2025 Table 1: V2m/fm = 256 x (TBW/70)

    # ----- Allometric exponents (theory-imposed, not estimated) -----
    # Yang 2025 Table 1 prints 0.75 on every clearance and an implicit
    # exponent of 1 on every volume, with no uncertainty. Two parameters
    # are enough because the same pair is reused for both analytes.
    e_wt_cl <- fixed(0.75); label("Allometric exponent of (WT / 70) on every clearance (unitless)")  # Yang 2025 Table 1: (TBW/70)^0.75 on CL/F, CLm/fm and Qm/fm
    e_wt_vc <- fixed(1); label("Allometric exponent of (WT / 70) on every volume (unitless)")  # Yang 2025 Table 1: (TBW/70) on Vc/F, V1m/fm and V2m/fm

    # ----- Inter-individual variability -----
    # Yang 2025 Table 1 reports, for this row, 'Imatinib: CV%(CL): 31.5%,
    # CV%(D1): 92.6%' and 'Metabolite: CV%(CLm): 33.4%'. No covariances are
    # tabulated, so the three etas are carried as independent diagonal
    # elements. The tabulated CV% is taken as omega (the log-scale standard
    # deviation), so the variance is (CV/100)^2 -- the convention used
    # throughout this transcription and in the shipped
    # Jiang_2023_imatinib.R. Note that NO IIV is reported on imatinib Vc/F
    # or on any metabolite volume.
    etalcl ~ 0.099225  # Yang 2025 Table 1: imatinib CV%(CL) 31.5% -> omega^2 = 0.315^2
    etald1 ~ 0.857476  # Yang 2025 Table 1: imatinib CV%(D1) 92.6% -> omega^2 = 0.926^2
    etalcl_ndmima ~ 0.111556  # Yang 2025 Table 1: metabolite CV%(CLm) 33.4% -> omega^2 = 0.334^2

    # ----- Residual unexplained variability -----
    # Yang 2025 Table 1 reports 'Log-scale add: 40.8%' (imatinib) and
    # 'Log-scale add: 34.7%' (metabolite). An additive residual on the
    # log-transformed observation is equivalent to a proportional residual
    # in linear space for small-to-moderate values, and the secondary
    # source reports the magnitude as a PERCENT, i.e. already converted to
    # the coefficient-of-variation scale. It is therefore encoded as
    # prop().
    propSd <- 0.408; label("Imatinib proportional residual error (fraction; log-scale additive in the source)")  # Yang 2025 Table 1: imatinib log-scale add 40.8%
    propSd_ndmima <- 0.347; label("N-desmethyl imatinib proportional residual error (fraction; log-scale additive in the source)")  # Yang 2025 Table 1: metabolite log-scale add 34.7%
  })

  model({
    # ----- 1. Allometric factors -----
    wt_cl <- (WT / 70)^e_wt_cl
    wt_vc <- (WT / 70)^e_wt_vc

    # ----- 2. Individual parameters -----
    cl <- exp(lcl + etalcl) * wt_cl
    vc <- exp(lvc) * wt_vc
    d1 <- exp(ld1 + etald1)

    cl_ndmima <- exp(lcl_ndmima + etalcl_ndmima) * wt_cl
    vc_ndmima <- exp(lvc_ndmima) * wt_vc
    q_ndmima <- exp(lq_ndmima) * wt_cl
    vp_ndmima <- exp(lvp_ndmima) * wt_vc

    # ----- 3. ODE system -----
    # Imatinib: one compartment, zero-order input. Dose records must carry
    # rate = -2 so rxode2 uses the modelled duration dur(central) = d1.
    d/dt(central) <- -cl * central / vc
    dur(central) <- d1

    # N-desmethyl imatinib: two compartments, formed from imatinib. The
    # metabolite states are carried in fm-SCALED amount units, i.e. they
    # hold A_m / fm rather than A_m. Writing the true central-metabolite
    # balance as
    #     dA_m/dt = fm * CL * Cc - CLm * A_m / Vm - Qm * (...)
    # and dividing throughout by fm gives the equations below with
    # cl_ndmima = CLm/fm, vc_ndmima = V1m/fm, q_ndmima = Qm/fm and
    # vp_ndmima = V2m/fm. The observed metabolite concentration is
    # unaffected because (A_m/fm) / (Vm/fm) = A_m / Vm, so Cc_ndmima below
    # is the true plasma concentration even though fm is not identifiable.
    # The formation flux is the FULL imatinib elimination flux
    # cl * central / vc: any fraction of imatinib clearance that does not
    # form N-desmethyl imatinib is absorbed into the unidentifiable fm and
    # cancels.
    d/dt(central_ndmima) <-
      cl * central / vc -
      cl_ndmima * central_ndmima / vc_ndmima -
      q_ndmima * central_ndmima / vc_ndmima +
      q_ndmima * peripheral1_ndmima / vp_ndmima
    d/dt(peripheral1_ndmima) <-
      q_ndmima * central_ndmima / vc_ndmima -
      q_ndmima * peripheral1_ndmima / vp_ndmima

    # ----- 4. Observations and error -----
    # States are mg (imatinib) or mg/fm (metabolite) and volumes are L, so
    # each ratio is mg/L; x 1000 gives ng/mL.
    Cc <- 1000 * central / vc
    Cc_ndmima <- 1000 * central_ndmima / vc_ndmima

    Cc ~ prop(propSd)
    Cc_ndmima ~ prop(propSd_ndmima)
  })
}
