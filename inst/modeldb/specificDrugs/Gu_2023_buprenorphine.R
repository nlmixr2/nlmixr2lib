Gu_2023_buprenorphine <- function() {
  description <- paste(
    "Parent-metabolite population PK model for sublingual buprenorphine and",
    "its metabolite norbuprenorphine after the buprenorphine/naloxone (4:1)",
    "sublingual tablet in Chinese male healthy volunteers and patients with",
    "opioid use disorder (Gu 2023). Buprenorphine: two-compartment",
    "disposition with Savic transit absorption (8 transit compartments fixed,",
    "estimated mean transit time) into an absorption compartment emptying at",
    "first-order ka; relative bioavailability falls linearly with the number",
    "of swallows after dosing (reference 3 swallows). Norbuprenorphine:",
    "two-compartment disposition formed from the whole buprenorphine",
    "clearance (fraction metabolised fixed to 1, multiplied by 1.98 in the",
    "16 mg dose group). The naloxone component is a separately fitted model",
    "(Gu_2023_naloxone)."
  )
  reference <- paste(
    "Gu M, Li A, Mak W, Dong F, Xu N, Zhang J, Shi Y, Zheng N, Tang Z, He Q,",
    "Ruan C, Guo W, Xiang X, Wang C, Han B, Zhu X (2023). Population",
    "pharmacokinetics of buprenorphine and naloxone sublingual combination",
    "in Chinese healthy volunteers and patients with opioid use disorder:",
    "Model-based dose optimization. Front. Pharmacol. 14:1089862.",
    "doi:10.3389/fphar.2023.1089862."
  )
  vignette <- "Gu_2023_buprenorphine_naloxone"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Amounts are mg of buprenorphine; norbuprenorphine is carried in
  # buprenorphine-mass equivalents (no molecular-weight conversion in the
  # supplementary control stream, CM = A(4)/V4*1000000).
  compartmentData <- list(
    depot = list(analyte = "buprenorphine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "buprenorphine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "buprenorphine", units = "mg", specimen = "tissue", verified = TRUE),
    central_norbup = list(analyte = "norbuprenorphine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_norbup = list(analyte = "norbuprenorphine", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    NSWALLOWS = list(
      description = "Number of swallows performed after sublingual administration of the dose",
      units = "(count)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Linear effect on buprenorphine relative bioavailability,",
        "F = 1 + (-0.035) * (NSWALLOWS - 3) (supplementary control stream",
        "'TVSLDO=DOSE*(1+THETA(8)*(NOS-3))'; Table 4 'The effect of NOS on",
        "F' = -0.035). The control stream centres on 3 swallows, not on the",
        "cohort median of 2 printed in Table 1. Observed range 0-12 in the",
        "single-dose groups; it was not recorded in the multiple-dose group",
        "and was imputed there with the 8 mg single-dose median (2)."
      ),
      source_name = "NOS"
    ),
    DOSE = list(
      description = "Labelled buprenorphine dose of the administered sublingual buprenorphine/naloxone dose",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Use case (a) of the DOSE canonical: the dose-group level (4, 8 or",
        "16 mg buprenorphine). Enters only as the 16 mg indicator on the",
        "fraction metabolised to norbuprenorphine, Fm = 1 * 1.98^(DOSE == 16)",
        "(control stream 'IF(DOSE.EQ.16) FFACTOR = THETA(18)'). Any other",
        "dose level takes Fm = 1, so doses other than 4, 8 and 16 mg are",
        "extrapolations for norbuprenorphine."
      ),
      source_name = "DOSE"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Tested in the stepwise covariate search (Section 4.3) and not retained."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Tested in the stepwise covariate search (Section 4.3) and not retained; height was dropped beforehand as highly correlated with weight."
    ),
    DIS_OUD = list(
      description = "Opioid use disorder patient (1) vs healthy volunteer (0)",
      units = "(binary)",
      type = "binary",
      notes = "Disease state was tested in the stepwise covariate search and no difference between healthy volunteers and OUD subjects was found (Section 4.3)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 46L,
    n_studies = 1L,
    age_range = "20-44 years",
    age_median = "27.5 years",
    weight_range = "53-88 kg",
    weight_median = "65.0 kg",
    sex_female_pct = 0,
    race_ethnicity = "Han Chinese (100%)",
    disease_state = "34 healthy volunteers (single dose) and 12 patients with opioid use disorder (multiple dose)",
    dose_range = paste(
      "Sublingual buprenorphine/naloxone tablets (2 mg buprenorphine +",
      "0.5 mg naloxone per tablet): single 4, 8 or 16 mg buprenorphine doses",
      "in healthy volunteers (n = 12, 12, 10) and 8 mg once daily for 8 days",
      "in OUD subjects (n = 12)"
    ),
    regions = "China (Beijing Anding Hospital)",
    notes = paste(
      "Table 1. 699 buprenorphine and 691 norbuprenorphine plasma",
      "concentrations; below-LLOQ (50 ng/L) data were under 10% and were",
      "discarded (M1). Healthy volunteers received naltrexone blockade around",
      "the single dose. Number of swallows after dosing median 2 (range",
      "0-12) in the single-dose groups."
    )
  )

  ini({
    # Table 4 (final buprenorphine-norbuprenorphine model). All structural
    # parameters are apparent (relative to the unknown absolute sublingual
    # bioavailability and fraction metabolised).
    lcl <- log(270)
    label("Apparent buprenorphine clearance CL/F (L/h)") # Table 4 'CL BUP /F, L/h' = 270
    lvc <- log(377)
    label("Apparent buprenorphine central volume Vc/F (L)") # Table 4 'V C BUP /F, L' = 377
    lvp <- log(5879)
    label("Apparent buprenorphine peripheral volume Vp/F (L)") # Table 4 'V P BUP /F, L' = 5879
    lq <- log(404)
    label("Apparent buprenorphine intercompartmental clearance Q/F (L/h)") # Table 4 'Q BUP /F, L/h' = 404
    lka <- log(0.397)
    label("Absorption rate constant from the absorption compartment (1/h)") # Table 4 'ka,/h' = 0.397
    lmtt <- log(0.234)
    label("Mean transit time from the first transit compartment to the absorption compartment (h)") # Table 4 'MTT, h' = 0.234
    lntr <- fixed(log(8))
    label("Number of transit compartments (unitless)") # Figure 1A 'Transit n=8'; control stream THETA(6) = 8 FIX
    e_nswallows_fdepot <- -0.035
    label("Linear change in relative bioavailability per additional swallow from 3 (1/swallow)") # Table 4 'The effect of NOS on F' = -0.035

    # Norbuprenorphine. The control stream assigns 264 to the clearance
    # (TVCLM = THETA(10)) and 22.2 to the central volume (TVV4 = THETA(9));
    # Table 4 prints the two row labels the other way round. The control
    # stream assignment is used (see the vignette for the NCA check).
    fm <- fixed(1)
    label("Fraction of buprenorphine clearance forming norbuprenorphine, 4 and 8 mg (unitless)") # control stream THETA(17) = 1 FIX; Section 3.1.1 'Fm was fixed to 1'
    e_dose16_fm <- 1.98
    label("Multiplier on Fm for the 16 mg dose group (unitless)") # Table 4 'Fm increment of 16 mg' = 1.98
    lcl_norbup <- log(264)
    label("Apparent norbuprenorphine clearance CL/Fm (L/h)") # Table 4 row printed 'V C norBUP /Fm, L' = 264; control stream TVCLM = THETA(10) = 264
    lvc_norbup <- log(22.2)
    label("Apparent norbuprenorphine central volume Vc/Fm (L)") # Table 4 row printed 'CL norBUP /Fm, L/h' = 22.2; control stream TVV4 = THETA(9) = 22.2
    lvp_norbup <- log(5170)
    label("Apparent norbuprenorphine peripheral volume Vp/Fm (L)") # Table 4 'V P norBUP /Fm, L' = 5170
    lq_norbup <- log(705)
    label("Apparent norbuprenorphine intercompartmental clearance Q/Fm (L/h)") # Table 4 'Q norBUP /Fm, L/h' = 705

    # IIV: Table 4 prints sqrt(omega^2) * 100 as '%CV' (V C BUP 135.6% =
    # sqrt(1.84), the control-stream $OMEGA value), so omega^2 = (CV/100)^2.
    # No IIV on ka or on the number of transit compartments (control stream
    # $OMEGA 0 FIX).
    etalcl ~ 0.110224 # Table 4 IIV 'CL BUP /F' = 33.2% -> 0.332^2
    etalvc ~ 1.838736 # Table 4 IIV 'V C BUP /F' = 135.6% -> 1.356^2
    etalvp ~ 0.171396 # Table 4 IIV 'V P BUP /F' = 41.4% -> 0.414^2
    etalq ~ 0.244036 # Table 4 IIV 'Q BUP /F' = 49.4% -> 0.494^2
    etalmtt ~ 0.140625 # Table 4 IIV 'MTT' = 37.5% -> 0.375^2
    etalcl_norbup ~ 0.255025 # Table 4 IIV row printed 'V C norBUP /Fm' = 50.5% -> 0.505^2; control stream ETA(9) on CLM = 0.255
    etalvc_norbup ~ 1.430416 # Table 4 IIV row printed 'CL norBUP /Fm' = 119.6% -> 1.196^2; control stream ETA(8) on V4 = 1.43
    etalvp_norbup ~ 0.203401 # Table 4 IIV 'V P norBUP /Fm' = 45.1% -> 0.451^2
    etalq_norbup ~ 0.321489 # Table 4 IIV 'Q norBUP /Fm' = 56.7% -> 0.567^2

    # Proportional residual error per analyte; the additive components are
    # fixed to 0 in the control stream.
    propSd <- 0.281
    label("Proportional residual error, buprenorphine (fraction)") # Table 4 'sigma prop (BUP)' = 28.1%
    propSd_norbup <- 0.22
    label("Proportional residual error, norbuprenorphine (fraction)") # Table 4 'sigma prop (norBUP)' = 22%
  })
  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    q <- exp(lq + etalq)
    ka <- exp(lka)
    mtt <- exp(lmtt + etalmtt)
    ntr <- exp(lntr)

    cl_norbup <- exp(lcl_norbup + etalcl_norbup)
    vc_norbup <- exp(lvc_norbup + etalvc_norbup)
    vp_norbup <- exp(lvp_norbup + etalvp_norbup)
    q_norbup <- exp(lq_norbup + etalq_norbup)

    # Relative bioavailability per swallow, centred on 3 swallows as in the
    # control stream (TVSLDO = DOSE * (1 + THETA(8) * (NOS - 3))).
    fdepot <- 1 + e_nswallows_fdepot * (NSWALLOWS - 3)
    fm_i <- fm * e_dose16_fm^(DOSE == 16)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    kel_norbup <- cl_norbup / vc_norbup
    k12_norbup <- q_norbup / vc_norbup
    k21_norbup <- q_norbup / vp_norbup

    # Savic analytical transit input (ktr = (ntr + 1) / mtt inside
    # transit()) into the absorption compartment; the control stream sets
    # F1 = 0 so the dose itself enters only through the transit term.
    d / dt(depot) <- transit(ntr, mtt, fdepot) - ka * depot
    d / dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d / dt(peripheral1) <- k12 * central - k21 * peripheral1
    d / dt(central_norbup) <- fm_i * kel * central - kel_norbup * central_norbup -
      k12_norbup * central_norbup + k21_norbup * peripheral1_norbup
    d / dt(peripheral1_norbup) <- k12_norbup * central_norbup - k21_norbup * peripheral1_norbup
    f(depot) <- 0

    # Dose in mg, volumes in L -> mg/L; x 1000 gives ng/mL.
    Cc <- 1000 * central / vc
    Cc_norbup <- 1000 * central_norbup / vc_norbup
    Cc ~ prop(propSd)
    Cc_norbup ~ prop(propSd_norbup)
  })
}
