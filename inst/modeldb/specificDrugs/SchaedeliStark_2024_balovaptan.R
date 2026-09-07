SchaedeliStark_2024_balovaptan <- function() {
  description <- "One-compartment population PK model of balovaptan with transit-compartment absorption, a dose-dependent central volume from empirical saturable binding, a turnover-gated gut extraction process, and brain V1a receptor occupancy, in adults and children with autism spectrum disorder (Schaedeli Stark 2024)"
  reference <- "Schaedeli Stark F, Chavanne C, Derks M, Jolling K, Lagraauw HM, Lindbom L, Prins K, Silber Baumann HE. A population pharmacokinetics model of balovaptan to support dose selection in adult and pediatric populations. J Pharmacokinet Pharmacodyn. 2024;51(3):227-242. doi:10.1007/s10928-023-09898-0"
  vignette <- "SchaedeliStark_2024_balovaptan"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # `propSd_early` / `propSd_late` / `propSd_exp_kdes` are the three parameters
  # of the paper's time-varying proportional RUV (Table 2 footnote):
  #   RUVprop = RUVearly - (RUVearly - RUVlate) * (1 - exp(-RUVrate * TAD))
  # The two SDs are not the canonical single `propSd`, and the third is the
  # decay rate constant of that time course, so all three are declared here.
  paper_specific_residual_sds <- c("propSd_early", "propSd_late",
                                   "propSd_exp_kdes")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Fig. 3 of the source.
  compartmentData <- list(
    depot      = list(analyte = "balovaptan", units = "mg", specimen = "administration site", verified = TRUE),
    transit1   = list(analyte = "balovaptan", units = "mg", specimen = "administration site", verified = TRUE),
    transit2   = list(analyte = "balovaptan", units = "mg", specimen = "administration site", verified = TRUE),
    transit3   = list(analyte = "balovaptan", units = "mg", specimen = "administration site", verified = TRUE),
    transit4   = list(analyte = "balovaptan", units = "mg", specimen = "administration site", verified = TRUE),
    transit5   = list(analyte = "balovaptan", units = "mg", specimen = "administration site", verified = TRUE),
    transit6   = list(analyte = "balovaptan", units = "mg", specimen = "administration site", verified = TRUE),
    central    = list(analyte = "balovaptan", units = "mg", specimen = "plasma", verified = TRUE),
    moderator1 = list(analyte = NA_character_, units = "(unitless fraction)", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Directly proportional (exponent fixed to 1.00) scaling of the baseline central volume, normalised to a 76 kg reference. 76 kg was the median weight of the neurotypical adult volunteers in the pooled dataset. Weight was NOT retained on CL/F: the estimated exponent (0.256) lost significance once the age maturation function was added.",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Drives the asymptotic CL/F maturation function. Reaches 50% of adult CL/F at AGE50 = 5.34 years and ~90% at 14 years; effectively at plateau by 20 years, and no further age effect up to 65 years.",
      source_name        = "AGE"
    ),
    FED = list(
      description        = "Fed state at dosing",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted)",
      notes              = "Multiplicative power-form effect on the mean transit time: MTT = MTT_fasted * 3.39^FED. The factor was fixed (not estimated) because only 4.3% of the participants with PK data were fasted. Both phase II ASD studies dosed with food, so the paper's own typical-participant simulations (Figs. 5 and 6) correspond to FED = 1.",
      source_name        = "FOOD"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 370,
    n_studies      = 5,
    age_range      = "5-64 years (median 19)",
    weight_range   = "18.8-152.0 kg (median 72.8)",
    sex_female_pct = 8.1,
    disease_state  = "autism spectrum disorder with IQ >= 70 (n = 315) and neurotypical adults (n = 55)",
    dose_range     = "1.5-52 mg once daily, oral",
    regions        = "not reported",
    notes          = "Pooled from three phase I studies in neurotypical adults (NCT01418963 n = 24, NCT03579719 n = 15, NCT03586726 n = 16; rich sampling, 5-52 mg) and two phase II ASD studies (VANILLA NCT01793441 n = 146 adults; aV1ation NCT02901431 n = 169 children and adolescents 5-17 years; sparse sampling, 1.5-10 mg). 3985 PK observations. The model was subsequently used to simulate IV infusion dosing and brain V1a receptor occupancy for a phase II malignant-cerebral-edema trial (NCT05399550); no IV data were fitted."
  )

  ini({
    # Table 2: apparent clearance CL = 8.52 L/h (95% CI 7.70-9.64)
    lcl <- log(8.52)
    label("Apparent clearance CL/F (L/hr)")
    # Table 2: apparent volume of distribution at baseline V0 = 565 L
    # (95% CI 468-686). This is the volume at zero central amount, before the
    # saturable-binding reduction; the dynamic Vc is derived in model().
    lvc <- log(565)
    label("Apparent central volume of distribution at baseline V0 (L)")
    # Table 2: mean transit time MTT = 0.367 h (95% CI 0.34-0.39). This is the
    # FASTED value; the fed value is 3.39 x larger (e_fed_mtt below).
    lmtt <- log(0.367)
    label("Mean transit time in the fasted state (hr)")
    # Table 2: number of transit compartments Ntr = 5.81 (95% CI 4.29-7.26).
    # Fixed here because the chain must be realised at an integer length; see
    # the model() note on integerisation.
    ntr <- fixed(5.81)
    label("Number of transit compartments")
    # Table 2: VLmax = 0.805 (95% CI 0.771-0.839), the maximum fraction by
    # which the central volume is reduced at full saturation of the
    # hypothetical binding sites.
    lvbmax <- log(0.805)
    label("Maximum fractional reduction of the central volume from saturable binding (fraction)")
    # Table 2: VLA50 = 3.36 (95% CI 2.16-5.21). Table 2 prints the unit as
    # "ug", but the value is in mg -- see the vignette Errata; the ug reading
    # saturates the binding term at every clinical dose and abolishes the
    # non-linearity the model exists to describe.
    lvba50 <- log(3.36)
    label("Central amount giving half the maximum volume reduction (mg)")
    # Table 2: Kgut = 1.89 1/h (95% CI 1.16-2.88)
    lkgutex <- log(1.89)
    label("Gut extraction rate constant from the absorption compartment (1/hr)")
    # Table 2: Kout = 0.00205 1/h (Fixed). Fixing Kout was reported to improve
    # numerical stability. Kin = Kout by the printed micro-constant relation.
    lkout <- fixed(log(0.00205))
    label("Gut-extraction turnover pool loss rate constant (1/hr)")
    # Table 2: S = 22.0 (95% CI 13.8-30.5), the factor by which the amount in
    # the absorption compartment stimulates Kout. Units 1/mg (depot is in mg).
    lsmod <- log(22.0)
    label("Scaling factor for the depot-amount stimulation of the turnover loss rate (1/mg)")
    # Table 2: WT on V0 = 1.00 (Fixed) -- "Scaling of Vc/F was directly
    # proportional to total body weight and fixed as such".
    e_wt_vc <- fixed(1.00)
    label("Body-weight exponent on the baseline central volume (76 kg reference)")
    # Table 2: FOOD on MTT = 3.39 (Fixed). Multiplicative, matching the sibling
    # "WT V0 = 1.00 (Fixed)" exponent row; the Results text calls it "an
    # estimate of 3.39 h", which is a unit slip -- see the vignette Errata.
    e_fed_mtt <- fixed(3.39)
    label("Multiplicative effect of the fed state on the mean transit time")
    # Table 2: AGEslope = 0.267 (95% CI 0.106-0.447)
    e_age_cl <- 0.267
    label("Slope of the CL/F age maturation function (1/year)")
    # Table 2: AGE50 = 5.34 years (95% CI 1.36-6.89), "age where 50% of adult
    # CL/F is reached".
    age50_cl <- 5.34
    label("Age at which 50 percent of adult CL/F is reached (years)")
    # Back-solved from Table 3: the paper prints the receptor-occupancy
    # equation RO = 100 * Cp * fu_plasma / (Kb + Cp * fu_plasma) but reports
    # neither Kb nor fu_plasma anywhere on disk. Only the ratio Kb/fu_plasma is
    # identifiable, and it is over-determined by Table 3's seven paired
    # (concentration, occupancy) medians; the rounding-interval intersection
    # across all seven rows is [3.076, 3.446] ng/mL. Value NOT from the paper
    # text or tables -- back-solved from printed medians; see vignette Errata.
    lkd_fu <- log(3.16)
    label("Brain V1a dissociation constant per unit total plasma concentration, Kb/fu_plasma (ng/mL)")

    # Table 2 footnote b: IIV reported as variances of the log-normal
    # distribution (39%, 31% and 21% CV respectively). Diagonal OMEGA.
    etalcl ~ 0.150
    etalvc ~ 0.0944
    etalmtt ~ 0.0423

    # Table 2: the time-varying proportional RUV decays mono-exponentially from
    # RUVearly = 1.20 (95% CI 0.811-1.77) immediately after dosing to
    # RUVlate = 0.216 (95% CI 0.201-0.232), at RUVrate = 1.32 1/h (95% CI
    # 0.992-1.670; a 0.53 h half-life, matching the Results text).
    propSd_early <- 1.20
    label("Proportional residual SD immediately after dosing (fraction)")
    propSd_late <- 0.216
    label("Proportional residual SD after the absorption phase (fraction)")
    propSd_exp_kdes <- 1.32
    label("Decay rate constant of the time-varying proportional residual SD (1/hr)")
    # Table 2 footnote a: additive RUV fixed to 0.025 ng/mL SD, half the assay
    # lower limit of quantitation. (Table 2 prints the variance, 0.000625;
    # 0.025^2 = 0.000625, and nlmixr2 parameterises the SD.)
    addSd <- fixed(0.025)
    label("Additive residual SD (ng/mL)")
    # The paper reports no residual error on receptor occupancy: RO is a
    # deterministic transform of plasma concentration, simulated rather than
    # fitted. Fixed to zero rather than invented.
    addSd_RO <- fixed(0)
    label("Additive residual SD on brain receptor occupancy (percentage points; not reported)")
  })

  model({
    # Apparent clearance with the asymptotic age-maturation function.
    #
    # The printed equation (p. 234) is
    #   CLi = CL * (1 - 1 / exp(-AGEslope * (AGE50 - AGE)))
    # which returns -143.9% of adult CL at age 2 and exactly 0% at AGE50 --
    # contradicting Table 2's own definition of AGE50 as the "age where 50% of
    # adult CL/F is reached". A single dropped "1 +" in the denominator repairs
    # it to the logistic below, which reproduces four independent printed
    # anchors (0.29 at age 2 and the 0.5 crossing per Fig. 7; ~90% at 14 y and
    # the ~20 y plateau per the Results text). See the vignette Errata.
    cl <- exp(lcl + etalcl) / (1 + exp(e_age_cl * (age50_cl - AGE)))

    # Dynamic central volume (Fig. 3 and p. 234):
    #   Vc = V0 * (1 - VLmax * Ac / (Ac + VLA50)),  Vc,i = Vc * (WT/76)^1.00
    # `central` is in mg, so VLA50 is in mg (see the ini() note).
    vbmax <- exp(lvbmax)
    vba50 <- exp(lvba50)
    vc <- exp(lvc + etalvc) * (WT / 76)^e_wt_vc *
      (1 - vbmax * central / (central + vba50))

    # Transit absorption. MTT = MTT_fasted * 3.39^FED; Ktr = (Ntr + 1) / MTT.
    # Ntr = 5.81 is not an integer, so the chain is realised with 6 transit
    # compartments (7 transfers: depot -> transit1..6 -> central). The printed
    # Ktr relation is kept intact, so the realised mean transit time is
    # 7/(Ntr+1) = 1.028 x MTT (+2.8%). See the vignette Errata.
    mtt <- exp(lmtt + etalmtt) * e_fed_mtt^FED
    ktr <- (ntr + 1) / mtt

    ke <- cl / vc

    # Gut extraction out of the absorption compartment only (Fig. 3 puts the
    # "Kgut x At" arrow on the Aa box, and the TCAM label on the Aa -> Ac
    # arrow). The turnover pool `moderator1` (the paper's At) starts at 1 and
    # is depleted by the amount in the depot, so extraction is large on the
    # first dose and negligible once the pool is suppressed by repeated dosing.
    kgutex <- exp(lkgutex)
    kout <- exp(lkout)
    kin <- kout
    smod <- exp(lsmod)

    moderator1(0) <- 1

    d/dt(depot) <- -ktr * depot - kgutex * moderator1 * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(transit4) <- ktr * transit3 - ktr * transit4
    d/dt(transit5) <- ktr * transit4 - ktr * transit5
    d/dt(transit6) <- ktr * transit5 - ktr * transit6
    d/dt(central) <- ktr * transit6 - ke * central
    d/dt(moderator1) <- kin - kout * moderator1 * (1 + smod * depot)

    # `central` is in mg and `vc` in L, so central/vc is mg/L = ug/mL; the
    # factor 1000 converts to the ng/mL used throughout the paper.
    Cc <- 1000 * central / vc

    # Brain V1a receptor occupancy (Methods, p. 230):
    #   RO = 100 * Cp * fu_plasma / (Kb + Cp * fu_plasma)
    #      = 100 * Cp / (Kb/fu_plasma + Cp)
    # Only the ratio Kb/fu_plasma is identifiable; see the ini() note.
    kd_fu <- exp(lkd_fu)
    RO <- 100 * Cc / (kd_fu + Cc)

    # Time-varying proportional RUV (Table 2 footnote). The printed expression
    # omits the minus sign in the exponent, which would send the SD to
    # infinity; the Results text ("mono-exponential decay from a high value
    # immediately after dosing to a lower level ... half-life of 0.53 h")
    # fixes the sign unambiguously. TAD is measured from the oral dose: every
    # study contributing to the fit was oral.
    propSdT <- propSd_late + (propSd_early - propSd_late) *
      exp(-propSd_exp_kdes * tad(depot))

    Cc ~ prop(propSdT) + add(addSd)
    RO ~ add(addSd_RO)
  })
}
