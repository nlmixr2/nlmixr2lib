Beliveau_2026_tenofovir_alafenamide_implant_weibull <- function() {
  description <- paste(
    "Weibull in vivo release variant of the CAPRISA 018 tenofovir alafenamide",
    "(TAF) subdermal implant model in South African cisgender women. The",
    "systemic half is identical to",
    "Beliveau_2026_tenofovir_alafenamide_implant: one-compartment plasma TAF",
    "with a dose effect on apparent volume, two-compartment plasma tenofovir",
    "(TFV) formed from the whole TAF elimination flux and scaled",
    "allometrically on body weight, and PBMC TFV diphosphate driven from",
    "plasma TFV by Michaelis-Menten formation with first-order loss. It",
    "replaces the zero-order release input with the Weibull function the",
    "authors fitted to per-participant Wagner-Nelson deconvolutions of plasma",
    "TAF, so the implant empties along a shape-controlled release profile",
    "instead of at a constant rate. Release parameters are the cohort-median",
    "descriptive statistics of those individual fits, not population estimates."
  )
  reference <- paste(
    "Beliveau M, Chang C, Lewis L, Letsoalo MP, Abdool Karim Q,",
    "Abdool Karim SS, Marzinke MA, Moss JA, Gengiah TN, Baum MM. Population",
    "pharmacokinetics of tenofovir alafenamide delivered via an annual",
    "subdermal implant in South African women. Sci Rep. 2026;16:18424.",
    "doi:10.1038/s41598-026-48746-2"
  )
  vignette <- "Beliveau_2026_tenofovir_alafenamide_implant"

  # Units, and the reason this file uses mass rather than the paper's molar
  # basis, are identical to the zero-order sibling; see the extended note in
  # Beliveau_2026_tenofovir_alafenamide_implant.R.
  units <- list(time = "h", dosing = "ug", concentration = "ug/L")

  # Issue #482. Same states as the zero-order sibling plus `depot`, the
  # undelivered TAF still inside the implant reservoir, which the Weibull
  # hazard empties.
  compartmentData <- list(
    depot = list(analyte = "tenofovir alafenamide", units = "ug", specimen = "administration site", verified = TRUE),
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
        "(V/F)TFV and -0.25 on K12, K21 and Ke(TFV). Plasma TAF, the implant",
        "release profile and PBMC TFV-DP carry no weight effect."
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
        "(Table 4). In this Weibull variant the SAME quantity is also the",
        "reference mass the release profile is expressed against, because the",
        "paper defines its relative-fraction-absorbed y-axis as modelled mass",
        "released divided by measured mass released (Results, 'Assessment of",
        "implant performance'). Dose the implant payload into `depot` as this",
        "mass converted to ug and let lfdepot carry the F_inf scaling; do not",
        "pre-multiply the two. Observed values: Group 2A median 22.0 mg,",
        "Group 2C median 36.8 mg, overall range 0.05-109 mg (Table 1)."
      ),
      source_name = "Dose"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 30L,
    n_studies = 1L,
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
      "alafenamide, in place for 4 weeks (Group 1) or up to 48 weeks (Groups",
      "2A and 2C); estimated delivered doses 0.05-109 mg"
    ),
    regions = "South Africa (CAPRISA 018; PACTR201809520959443)",
    notes = paste(
      "Same cohort as Beliveau_2026_tenofovir_alafenamide_implant; see that",
      "file for the full baseline description and the BLQ handling. What",
      "differs here is the STATUS of the release parameters. Table 5 reports",
      "descriptive statistics -- mean, SD, median, CV%, min and max -- of b,",
      "F_inf and MDT across per-participant Weibull fits to Wagner-Nelson",
      "deconvolutions, NOT a population model with typical values and an",
      "estimated omega. This file therefore carries the OVERALL median of",
      "each (n = 30) as its typical value and encodes no IIV, even though the",
      "spread is large and real: CV% across individuals was 65.8 for b, 102.6",
      "for F_inf and 131.5 for MDT, with MDT ranging 97.2-73,800 h. Group 1",
      "differs systematically from Groups 2A and 2C because its implants were",
      "removed on schedule at about 28 days, truncating both the fraction",
      "released and the apparent MDT (Group 1 median MDT 323 h versus 10,600",
      "h and 7810 h); a user simulating a 4-week course should prefer the",
      "Group 1 column of Table 5 to the pooled medians used here. The authors",
      "classify 57% of implants as showing linear (zero-order) release and",
      "43% as nonlinear (Results), so the zero-order sibling model is not a",
      "worse description of the cohort -- it is the better description of the",
      "majority of it."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # In vivo implant release. Beliveau 2026 deconvolved each participant's
    # plasma TAF profile by the Wagner-Nelson method and fitted the resulting
    # fraction absorbed with (Methods, 'Evaluation of implant performance'):
    #
    #     y(t) = F_inf * {1 - exp[-(t / MDT)^b]}
    #
    # which is the canonical nlmixr2lib Weibull release input with ra = 1/MDT
    # and gam1 = b, so the fraction still unreleased is exp(-(ra*t)^gam1).
    # All three values are the Overall (n = 30) MEDIAN column of Table 5 and
    # are wrapped in fixed() because they are descriptive statistics of
    # individual fits rather than estimated population parameters.
    # ---------------------------------------------------------------------
    lra <- fixed(log(1 / 6390))
    label("Weibull release rate-scaling parameter, reciprocal of the mean dissolution time MDT (1/h)")
    # Table 5, block 'MDT', row 'Median (CV%)', column 'Overall ( N =30)'
    # = 6390 h (CV 131.5%); ra = 1/MDT. NOTE the paper's prose calls MDT 'the
    # mid-point of release (i.e., time at which half of the measured total
    # implant dose has been delivered)', but the printed equation makes
    # y(MDT) = F_inf * (1 - 1/e) = 0.632 * F_inf, not half. The equation is
    # encoded here and the prose is recorded as a deviation in the vignette.

    lgam1 <- fixed(log(1.36))
    label("Weibull release shape (slope) parameter b (unitless)")
    # Table 5, block 'b', row 'Median (CV%)', column 'Overall ( N =30)' = 1.36
    # (CV 65.8%). Individual values spanned 0.653-4.97; b > 1 gives the
    # sigmoidal, initially-accelerating release the model assumes, and values
    # below 1 would make the hazard below infinite at t = 0.

    lfdepot <- fixed(log(1.46))
    label("Weibull total fraction absorbed F_inf, relative to the measured mass released from the implant (unitless)")
    # Table 5, block 'F inf', row 'Median (CV%)', column 'Overall ( N =30)'
    # = 1.46 (CV 102.6%). This is NOT a bioavailability bounded by 1: the
    # paper defines it as modelled mass released divided by MEASURED mass
    # released, so 1.46 means the deconvolution implies 46% more TAF reached
    # the circulation than the residual-drug assay of the used implants
    # accounted for (Results, 'Assessment of implant performance': 'a value of
    # 2 indicates that the model-predicted mass released is double what was
    # measured in the used implant'). Individual values spanned 1.01-15.9.

    # ---------------------------------------------------------------------
    # Everything below is identical to
    # Beliveau_2026_tenofovir_alafenamide_implant.R; see that file for the
    # full source-trace commentary on each value.
    # ---------------------------------------------------------------------
    lkel <- fixed(log(0.924))
    label("TAF elimination rate constant Ke(TAF), literature value (1/h)")
    # Table 3, row 'K e (TAF)' = 0.924 /h (footnote a: ln(2)/0.5 h x 0.66).

    lvc <- fixed(log(62.6))
    label("TAF apparent central volume of distribution (V/F)TAF at the reference 17.3 mg delivered dose, literature value (L)")
    # Table 3, row '( V/F ) TAF' = 62.6 L.

    e_dose_taf_mg_vc <- 0.618
    label("Power of total delivered TAF dose on (V/F)TAF, centred on 17.3 mg (unitless)")
    # Table 4, row 'Dose effect on TAF volume' = (Dose/17.3)^0.618.

    fm_tfv <- 9.24
    label("Apparent relative bioavailability of implant-derived TFV versus oral dosing, multiplying TAF-to-TFV formation (unitless)")
    # Table 4, row 'F rel (TAF)' = 9.24; placed on the TFV arm per Fig. S1's
    # legend ('Frel, TFV bioavailability'). See the vignette Errata.

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

    km_tfvdp <- fixed(29.3)
    label("Michaelis-Menten constant for TFV-to-TFV-DP conversion, literature value (ug/L)")
    # Table 3, row 'K m' = 29.3 ug/L.

    lvmax_tfvdp <- fixed(log(1.44))
    label("Maximum velocity of TFV-to-TFV-DP conversion, literature value (fmol per 10^6 cells per h)")
    # Table 3, row 'V max' = 1.44 fmol/10^6 cells/h.

    lkel_tfvdp <- fixed(log(0.006))
    label("PBMC TFV-DP elimination rate constant Ke(TFV-DP), literature value (1/h)")
    # Table 3, row 'K e (TFV-DP)' = 0.006 /h.

    # No IIV: Table 5's CV% columns describe the spread of INDIVIDUAL Weibull
    # fits, which also carries their estimation error, and are not an
    # estimated omega. See the sibling file for the same reasoning about the
    # unreported random effect on F.
    propSd <- fixed(0)
    label("Proportional residual SD, plasma TAF (fraction; 0 -- not reported in the source)")
    propSd_tfv <- fixed(0)
    label("Proportional residual SD, plasma TFV (fraction; 0 -- not reported in the source)")
    propSd_Cpbmc_tfvdp <- fixed(0)
    label("Proportional residual SD, PBMC TFV-DP (fraction; 0 -- not reported in the source)")
  })

  model({
    # Molecular weights (g/mol) of tenofovir alafenamide free base and
    # tenofovir. Standard chemical constants, NOT reported by Beliveau 2026;
    # they appear only as the ratio converting the molar 1:1 TAF-to-TFV
    # conversion into this file's mass units.
    mwTaf <- 476.47
    mwTfv <- 287.21

    ra <- exp(lra)
    gam1 <- exp(lgam1)

    kel <- exp(lkel)
    vc <- exp(lvc) * (DOSE_TAF_MG / 17.3)^e_dose_taf_mg_vc

    wtNorm <- WT / 70.8
    vc_tfv <- exp(lvc_tfv) * wtNorm^e_wt_vc_tfv
    k12_tfv <- exp(lk12_tfv) * wtNorm^e_wt_k12_tfv
    k21_tfv <- exp(lk21_tfv) * wtNorm^e_wt_k21_tfv
    kel_tfv <- exp(lkel_tfv) * wtNorm^e_wt_kel_tfv

    kel_tfvdp <- exp(lkel_tfvdp)
    vmax_tfvdp <- exp(lvmax_tfvdp)

    # Weibull in vivo release. The fraction of the implant payload still
    # unreleased is S(t) = exp(-(ra*t)^gam1), so the amount remaining is
    # emptied by the Weibull hazard h(t) = gam1 * ra * (ra*t)^(gam1 - 1).
    # Integrating that against a depot loaded with F_inf * Dose reproduces the
    # paper's y(t) = F_inf * {1 - exp[-(t/MDT)^b]} exactly, with ra = 1/MDT
    # and gam1 = b. Time is measured from implant insertion, so tad() is read
    # once into a local; with the typical gam1 = 1.36 > 1 the hazard is 0 at
    # insertion and rises, which is the accelerating release the shape implies.
    tRel <- tad()
    krel <- gam1 * ra * (ra * tRel)^(gam1 - 1)

    d/dt(depot) <- -krel * depot
    d/dt(central) <- krel * depot - kel * central

    # F_inf scales the payload dosed into the implant reservoir. It is a
    # fraction of the MEASURED released mass and legitimately exceeds 1.
    f(depot) <- exp(lfdepot)

    # As in the zero-order sibling, the whole TAF elimination flux forms TFV
    # 1:1 on a molar basis scaled by fm_tfv, and is NOT subtracted from the
    # parent, so the system is deliberately not mass-conserving.
    d/dt(central_tfv) <-
      fm_tfv * kel * central * (mwTfv / mwTaf) -
      kel_tfv * central_tfv -
      k12_tfv * central_tfv + k21_tfv * peripheral1_tfv
    d/dt(peripheral1_tfv) <- k12_tfv * central_tfv - k21_tfv * peripheral1_tfv

    Cc <- central / vc
    Cc_tfv <- central_tfv / vc_tfv

    d/dt(pbmc_tfvdp) <-
      vmax_tfvdp * Cc_tfv / (km_tfvdp + Cc_tfv) - kel_tfvdp * pbmc_tfvdp
    Cpbmc_tfvdp <- pbmc_tfvdp

    Cc ~ prop(propSd)
    Cc_tfv ~ prop(propSd_tfv)
    Cpbmc_tfvdp ~ prop(propSd_Cpbmc_tfvdp)
  })
}
