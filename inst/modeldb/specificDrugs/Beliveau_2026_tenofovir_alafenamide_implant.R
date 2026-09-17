Beliveau_2026_tenofovir_alafenamide_implant <- function() {
  description <- paste(
    "Three-analyte population PK model for tenofovir alafenamide (TAF)",
    "delivered by an ultralong-acting subdermal implant in South African",
    "cisgender women (CAPRISA 018 first-in-human trial). Zero-order in vivo",
    "TAF release feeds a one-compartment plasma TAF model whose apparent",
    "central volume carries a power effect of the total TAF dose delivered by",
    "the implant(s); the whole TAF elimination flux forms plasma tenofovir",
    "(TFV) in a two-compartment model allometrically scaled on body weight;",
    "and peripheral blood mononuclear cell (PBMC) TFV diphosphate (TFV-DP) is",
    "driven from plasma TFV by Michaelis-Menten formation with first-order",
    "loss. Every systemic disposition parameter is fixed from the literature;",
    "only the dose effect on apparent TAF volume and the relative",
    "bioavailability of implant-derived TFV versus oral dosing were estimated."
  )
  reference <- paste(
    "Beliveau M, Chang C, Lewis L, Letsoalo MP, Abdool Karim Q,",
    "Abdool Karim SS, Marzinke MA, Moss JA, Gengiah TN, Baum MM. Population",
    "pharmacokinetics of tenofovir alafenamide delivered via an annual",
    "subdermal implant in South African women. Sci Rep. 2026;16:18424.",
    "doi:10.1038/s41598-026-48746-2"
  )
  vignette <- "Beliveau_2026_tenofovir_alafenamide_implant"

  # Beliveau 2026 ran the model on a molar basis (Methods, "Implementation and
  # evaluation of implant pharmacokinetic model": "TAF doses and analyte
  # concentrations were converted to molar amounts"). This file instead carries
  # TAF and TFV as masses in ug, so that concentrations come out directly in
  # ug/L = ng/mL -- the units in which the paper reports every TAF and TFV
  # concentration and in which Km is tabulated (Table 3, Km = 29.3 ug/L). The
  # two parameterisations are algebraically identical; the only price is the
  # single TFV/TAF molecular-weight ratio applied at the conversion step in
  # model(). PBMC TFV-DP is carried in its reported unit, fmol per 10^6 cells.
  units <- list(time = "h", dosing = "ug", concentration = "ug/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Tenofovir alafenamide is the dosed parent and so keeps
  # the bare canonical `central` / `Cc` names; tenofovir takes the registered
  # `_tfv` metabolite suffix and tenofovir diphosphate the `_tfvdp` suffix
  # (compartment-names.md, parent-wins rule; precedents
  # Thoueille_2023_tenofovir_alafenamide.R and Yu_2026_tenofovir.R).
  # `pbmc_tfvdp` is a BIOPHASE state carried in CONCENTRATION units, not an
  # amount: Vmax is tabulated per 10^6 cells, so no volume divides it at the
  # observation step.
  compartmentData <- list(
    central = list(analyte = "tenofovir alafenamide", units = "ug", specimen = "plasma", verified = TRUE),
    central_tfv = list(analyte = "tenofovir", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1_tfv = list(analyte = "tenofovir", units = "ug", specimen = "plasma", verified = TRUE),
    pbmc_tfvdp = list(
      analyte = "tenofovir diphosphate",
      units = "fmol/10^6 cells",
      specimen = "blood cell",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Body weight at baseline.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters every plasma TFV disposition parameter as a power function",
        "centred on 70.8 kg, carried over unchanged from the oral TAF",
        "literature model of Ji et al. (Beliveau 2026 Table 3): exponent 1 on",
        "(V/F)TFV and -0.25 on K12, K21 and Ke(TFV). Plasma TAF and PBMC",
        "TFV-DP carry no weight effect. The 70.8 kg centring value is the",
        "reference of the source literature model, NOT the CAPRISA 018 median",
        "-- study weights were 49.1-90.9 kg with group means 68.2-72.5 kg",
        "(Table 1), so the cohort happens to sit close to the reference.",
        "Beliveau 2026 screened baseline body weight as a covariate on the",
        "pre-systemic (implant release / relative availability) parameters and",
        "did not retain it there (Methods, third modelling assumption)."
      ),
      source_name = "Weight"
    ),
    DOSE_TAF_MG = list(
      description = paste(
        "Total mass of tenofovir alafenamide delivered in vivo by the",
        "implant(s) over the whole insertion period, estimated from the",
        "residual drug assayed in the used implants."
      ),
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters apparent central TAF volume (V/F)TAF as a power function",
        "centred on the cohort median 17.3 mg: (V/F)TAF x (Dose/17.3)^0.618",
        "(Table 4). This is the ONLY covariate Beliveau 2026 retained, and it",
        "was significant at P < 0.001 (Results, 'Model performance'). It is a",
        "per-subject constant describing the whole implant course and CANNOT",
        "be derived from the dosing amount on a record: the model is dosed as",
        "a zero-order release RATE in ug/h, whereas this covariate is the",
        "cumulative mg delivered. Set it on every record of a subject.",
        "Observed values: Group 1 median 0.2 mg, Group 2A median 22.0 mg,",
        "Group 2C median 36.8 mg, overall range 0.05-109 mg (Table 1). The",
        "authors read the effect as saturable absorption from the",
        "subcutaneous space: doubling the dose from 17.3 to 34.6 mg raises",
        "V/F by 54% (2^0.618 = 1.535), i.e. a 35% fall in apparent",
        "bioavailability (Results)."
      ),
      source_name = "Dose"
    )
  )

  # Covariates Beliveau 2026 screened graphically against the pre-systemic
  # (implant release and relative availability) parameters and did NOT retain
  # (Methods, third modelling assumption). Only the dose effect survived.
  # Documented here for provenance; none is referenced in model().
  covariatesDataExcluded <- list(
    CRCL = list(
      description = "Baseline creatinine clearance.",
      units = "mL/min",
      type = "continuous",
      notes = paste(
        "Screened on the pre-systemic parameters, not retained. Group means",
        "147-166 mL/min, range 98.0-203 mL/min (Table 1)."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 30L,
    n_studies = 1L,
    n_observations = paste(
      "493 plasma samples assayed for TAF and for TFV (TAF: 355 quantifiable,",
      "138 post-dose below the limit of quantification; TFV: only 46",
      "quantifiable and 447 post-dose BLQ) and 172 PBMC samples assayed for",
      "TFV-DP (91 quantifiable, 81 post-dose BLQ), all from Table 2"
    ),
    age_range = "18-38 years",
    weight_range = "49.1-90.9 kg",
    sex_female_pct = 100,
    race_ethnicity = c(Black = 100),
    disease_state = paste(
      "Healthy, HIV-negative cisgender women enrolled for HIV-1 pre-exposure",
      "prophylaxis; not a disease population"
    ),
    dose_range = paste(
      "One or two subdermal implants each containing 110 +/- 10 mg tenofovir",
      "alafenamide. Group 1 (n = 6, lead-in) one implant for 4 weeks; Group 2A",
      "(n = 12) one implant and Group 2C (n = 12) two implants, both for up to",
      "48 weeks. Only the three active arms entered this analysis; the placebo",
      "arms (Groups 2B and 2D) did not. The mass actually released in vivo was",
      "far smaller than the payload: estimated delivered doses were 0.05-109 mg"
    ),
    regions = "South Africa (CAPRISA 018; PACTR201809520959443)",
    notes = paste(
      "Baseline characteristics from Beliveau 2026 Table 1. All 30",
      "participants were Black South African women. Implant insertion",
      "duration differed sharply by group (Group 1 median 677 h by design;",
      "Group 2A median 7930 h; Group 2C median 7170 h), with the medicated",
      "groups shortened by early removals for implant-site reactions.",
      "Additional covariates screened against the pre-systemic parameters and",
      "not retained: implant arm (left/right), compliance (early versus",
      "scheduled removal), implant dosing duration, contraception method",
      "(contraceptive implant, injectable, or intrauterine device),",
      "mid-upper arm circumference (means 29.4-30.4 cm), and the presence of",
      "an implant-site reaction, hyperpigmentation, pruritis or induration.",
      "Only creatinine clearance among these has a canonical covariate name",
      "and so is the only one carried in covariatesDataExcluded.",
      "Concentrations below the assay limits of quantification (TAF 0.03",
      "ng/mL; TFV 1 ng/mL; TFV-DP 1.7 fmol/10^6 cells) were set to missing",
      "rather than censored: the authors tried an M3-style censoring model and",
      "abandoned it as 'highly unstable with high %CV values on parameters'",
      "(Results). That matters when comparing simulations with the published",
      "observed summaries, because 26% of post-dose TAF and 86% of post-dose",
      "TFV samples were BLQ and the surviving quantifiable values are",
      "therefore biased upward."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Plasma TAF. Both values are literature constants carried over from the
    # oral TAF model of Ji et al. and were NOT re-estimated here (Results,
    # 'Model performance': 'Systemic model parameters ... were taken from the
    # literature (Table 3) and only the relative bioavailability from oral to
    # implant dosing needed to be estimated').
    # ---------------------------------------------------------------------
    lkel <- fixed(log(0.924))
    label("TAF elimination rate constant Ke(TAF), literature value (1/h)")
    # Table 3, row 'K e (TAF)' = 0.924 /h. Footnote a builds it from an assumed
    # 0.5 h plasma half-life and the sex effect of Ji et al.: ln(2)/0.5 h x 0.66
    # = 1.4 x 0.66 = 0.924. The footnote's trailing unit 'L' is a typo; this is
    # a first-order rate constant. The 0.66 factor is a FEMALE-specific
    # adjustment and is already baked into this number -- every CAPRISA 018
    # participant was female, so the model carries no sex covariate.

    lvc <- fixed(log(62.6))
    label("TAF apparent central volume of distribution (V/F)TAF at the reference 17.3 mg delivered dose, literature value (L)")
    # Table 3, row '( V/F ) TAF' = 62.6 L.

    # ---------------------------------------------------------------------
    # The two parameters Beliveau 2026 actually fitted (Table 4, both marked
    # 'Model fit'). Neither is wrapped in fixed().
    # ---------------------------------------------------------------------
    e_dose_taf_mg_vc <- 0.618
    label("Power of total delivered TAF dose on (V/F)TAF, centred on 17.3 mg (unitless)")
    # Table 4, row 'Dose effect on TAF volume' = (V/F)TAF x (Dose/17.3)^0.618;
    # the exponent is also quoted in Results ('The exponent for the dose effect
    # on the typical (V/F)TAF value was 0.618'). Cross-check: the Results
    # paragraph on saturable absorption states that doubling 17.3 -> 34.6 mg
    # gives a 54% increase in V/F, and 2^0.618 = 1.535, i.e. +54%.

    fm_tfv <- 9.24
    label("Apparent relative bioavailability of implant-derived TFV versus oral dosing, multiplying TAF-to-TFV formation (unitless)")
    # Table 4, row 'F rel (TAF)' = 9.24. NOTE the value exceeds 1 and is an
    # APPARENT quantity: (V/F)TFV below is an oral-apparent volume that already
    # carries the oral bioavailability of the literature model, so this factor
    # absorbs the implant-versus-oral difference in TFV availability rather
    # than being a true mass fraction. See the extended note in the vignette's
    # Errata: Table 4 labels the row 'Frel(TAF)', but Fig. S1's legend defines
    # it as 'Frel, TFV bioavailability' and the supplementary simulation
    # figures place it unambiguously on the TFV arm.

    # ---------------------------------------------------------------------
    # Plasma TFV. All literature constants (Table 3), all allometrically
    # scaled on body weight centred at 70.8 kg.
    # ---------------------------------------------------------------------
    lvc_tfv <- fixed(log(1360))
    label("TFV apparent central volume of distribution (V/F)TFV at 70.8 kg, literature value (L)")
    # Table 3, row '( V/F ) TFV' = 1360 x (Weight/70.8)^1.

    e_wt_vc_tfv <- fixed(1)
    label("Power of body weight on (V/F)TFV, centred on 70.8 kg (unitless)")
    # Table 3, row '( V/F ) TFV' exponent = 1.

    lk12_tfv <- fixed(log(0.2257))
    label("TFV central-to-peripheral transfer rate constant K12 at 70.8 kg, literature value (1/h)")
    # Table 3, row 'K 12' = 0.2257 x (Weight/70.8)^-0.25.

    e_wt_k12_tfv <- fixed(-0.25)
    label("Power of body weight on K12, centred on 70.8 kg (unitless)")
    # Table 3, row 'K 12' exponent = -0.25.

    lk21_tfv <- fixed(log(0.2981))
    label("TFV peripheral-to-central transfer rate constant K21 at 70.8 kg, literature value (1/h)")
    # Table 3, row 'K 21' = 0.2981 x (Weight/70.8)^-0.25.

    e_wt_k21_tfv <- fixed(-0.25)
    label("Power of body weight on K21, centred on 70.8 kg (unitless)")
    # Table 3, row 'K 21' exponent = -0.25.

    lkel_tfv <- fixed(log(0.039))
    label("TFV elimination rate constant Ke(TFV) at 70.8 kg, literature value (1/h)")
    # Table 3, row 'K e (TFV)' = 0.039 x (Weight/70.8)^-0.25.

    e_wt_kel_tfv <- fixed(-0.25)
    label("Power of body weight on Ke(TFV), centred on 70.8 kg (unitless)")
    # Table 3, row 'K e (TFV)' exponent = -0.25.

    # ---------------------------------------------------------------------
    # PBMC TFV-DP. Literature constants from Duwal et al. (Table 3), carried
    # over unchanged.
    # ---------------------------------------------------------------------
    km_tfvdp <- fixed(29.3)
    label("Michaelis-Menten constant for TFV-to-TFV-DP conversion, literature value (ug/L)")
    # Table 3, row 'K m' = 29.3 ug/L. Expressed against the plasma TFV
    # concentration, which is why this file keeps TFV in mass units.

    lvmax_tfvdp <- fixed(log(1.44))
    label("Maximum velocity of TFV-to-TFV-DP conversion, literature value (fmol per 10^6 cells per h)")
    # Table 3, row 'V max' = 1.44 fmol/10^6 cells/h.

    lkel_tfvdp <- fixed(log(0.006))
    label("PBMC TFV-DP elimination rate constant Ke(TFV-DP), literature value (1/h)")
    # Table 3, row 'K e (TFV-DP)' = 0.006 /h.

    # ---------------------------------------------------------------------
    # Inter-individual variability is DELIBERATELY ABSENT. Beliveau 2026 did
    # carry a random effect on F -- the simulation Methods use 'individual
    # random effect values of F' -- but reports no variance, %CV or shrinkage
    # for it anywhere in the paper or the supplement, and reports no residual
    # error model at all. No variance is invented here. The etas are omitted
    # rather than written as `~ fixed(0)`, because a zero-variance diagonal
    # makes OMEGA singular and rxSolve then fails in chol(). Downstream users
    # who want a stochastic cohort must supply their own omega.
    # ---------------------------------------------------------------------
    propSd <- fixed(0)
    label("Proportional residual SD, plasma TAF (fraction; 0 -- not reported in the source)")
    propSd_tfv <- fixed(0)
    label("Proportional residual SD, plasma TFV (fraction; 0 -- not reported in the source)")
    propSd_Cpbmc_tfvdp <- fixed(0)
    label("Proportional residual SD, PBMC TFV-DP (fraction; 0 -- not reported in the source)")
  })

  model({
    # Molecular weights (g/mol) of tenofovir alafenamide free base
    # (C21H29N6O5P) and tenofovir (C9H14N5O4P). These are standard chemical
    # constants and are NOT reported by Beliveau 2026, which states only that
    # the analysis was run on a molar basis. They appear here solely as the
    # ratio that converts the molar 1:1 TAF-to-TFV conversion into the mass
    # units this file uses; in the paper's own molar parameterisation the
    # ratio is absent because it equals 1.
    mwTaf <- 476.47
    mwTfv <- 287.21

    # Plasma TAF. The apparent central volume carries the only covariate the
    # authors retained: a power function of the total TAF mass the implant(s)
    # released in vivo, centred on the cohort median 17.3 mg.
    kel <- exp(lkel)
    vc <- exp(lvc) * (DOSE_TAF_MG / 17.3)^e_dose_taf_mg_vc

    # Plasma TFV. Every disposition parameter is allometrically scaled on body
    # weight centred at the literature model's 70.8 kg reference.
    wtNorm <- WT / 70.8
    vc_tfv <- exp(lvc_tfv) * wtNorm^e_wt_vc_tfv
    k12_tfv <- exp(lk12_tfv) * wtNorm^e_wt_k12_tfv
    k21_tfv <- exp(lk21_tfv) * wtNorm^e_wt_k21_tfv
    kel_tfv <- exp(lkel_tfv) * wtNorm^e_wt_kel_tfv

    kel_tfvdp <- exp(lkel_tfvdp)
    vmax_tfvdp <- exp(lvmax_tfvdp)

    # Zero-order in vivo release from the implant is supplied by the event
    # table as an infusion into `central` (Results: 'The current model used a
    # zero-order absorption rate rather than the first-order rates observed
    # following oral dosing'). TAF declines mono-exponentially thereafter
    # (Methods, fourth modelling assumption).
    d/dt(central) <- -kel * central

    # The WHOLE TAF elimination flux forms TFV, 1:1 on a molar basis, scaled
    # by the apparent relative bioavailability fm_tfv. Following the registered
    # `fm_<metabolite>` convention the formation flux is NOT subtracted from
    # the parent: `central` retains its full -kel*central elimination, so the
    # TAF arm is notionally split rather than extended. The system is
    # therefore deliberately NOT mass-conserving -- (V/F)TFV is an
    # oral-apparent volume, so fm_tfv exceeds 1 by construction.
    d/dt(central_tfv) <-
      fm_tfv * kel * central * (mwTfv / mwTaf) -
      kel_tfv * central_tfv -
      k12_tfv * central_tfv + k21_tfv * peripheral1_tfv
    d/dt(peripheral1_tfv) <- k12_tfv * central_tfv - k21_tfv * peripheral1_tfv

    Cc <- central / vc
    Cc_tfv <- central_tfv / vc_tfv

    # PBMC TFV-DP is linked to the plasma TFV CONCENTRATION by Michaelis-Menten
    # formation with first-order loss (Methods, sixth modelling assumption).
    # The state is itself a concentration in fmol per 10^6 cells -- Vmax is
    # already tabulated per 10^6 cells per hour -- so no volume divides it. The
    # formation flux is a read-out and does not feed back on plasma TFV; the
    # amounts involved are negligible and the units are not commensurable.
    d/dt(pbmc_tfvdp) <-
      vmax_tfvdp * Cc_tfv / (km_tfvdp + Cc_tfv) - kel_tfvdp * pbmc_tfvdp
    Cpbmc_tfvdp <- pbmc_tfvdp

    Cc ~ prop(propSd)
    Cc_tfv ~ prop(propSd_tfv)
    Cpbmc_tfvdp ~ prop(propSd_Cpbmc_tfvdp)
  })
}
