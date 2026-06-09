Sanhajariya_2018_snake_venom <- function() {
  description <- "Exploratory population PK meta-analysis of snake venom in humans (Sanhajariya 2018): one-compartment model with zero-order input (duration D1 = 1 h, fixed) and first-order elimination, fit in NONMEM 7.2 to 218 timed venom concentrations from 145 snakebite patients pooled across 24 published case reports / series. Snake family (Elapidae vs Viperidae) modifies F1; Viperidae is the reference (F1 = 1, fixed). Authors describe the model as a preliminary prior for future snake-envenoming PK modelling; F1 also absorbs the large bite-to-bite variability in injected venom mass."
  reference <- "Sanhajariya S, Duffull SB, Isbister GK. Pharmacokinetics of snake venom. Toxins (Basel). 2018;10(2):73. doi:10.3390/toxins10020073"
  vignette <- "Sanhajariya_2018_snake_venom"
  units <- list(time = "hour", dosing = "mcg", concentration = "mcg/L")

  covariateData <- list(
    SNAKEFAMILY_ELAPID = list(
      description        = "Snake-family categorical indicator for venom-source classification: 1 = Elapidae family (front-fanged elapids: cobras, kraits, mambas, sea snakes, Australian terrestrial elapids); 0 = Viperidae family (true vipers and pit vipers).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Viperidae)",
      notes              = "Per-bite covariate (property of the snake that delivered the dose, not of the patient). Sanhajariya 2018 Table A1 covariate model: F1(Viperidae) = 1 (fixed reference), F1(Elapidae) = 0.569 (RSE 43%). The two-family split was the deepest stratification the sparse meta-analysis data could support; per-species or per-genus effects could not be estimated (Section 2.2.3 prose).",
      source_name        = "Elapidae / Viperidae family label (no NONMEM column name disclosed in Appendix Table A1)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 145,
    n_studies      = 24,
    age_range      = "Not reported in the main text; the underlying 24 primary studies are case reports / series of snake-envenoming patients (Table 6).",
    weight_range   = "Not reported in the main text.",
    sex_female_pct = "Not reported in the main text.",
    race_ethnicity = "Not reported in the main text.",
    disease_state  = "Snake envenoming, pre-antivenom samples only. Elapid bites (cobras, kraits, sea snakes, Australian elapids, taipans) and viperid bites (Bothrops, Crotalus, Daboia / Vipera russelli, Vipera aspis / berus / ammodytes, Bitis, Hypnale, Cerastes).",
    dose_range     = "Single snakebite; the venom mass per bite is not directly measured. The model treats the dose as a nominal unit and F1 absorbs the per-bite relative venom amount (CV 275%) on top of the snake-family multiplier (0.569 for Elapidae vs 1 for Viperidae). For simulation, supply a representative venom mass in mcg.",
    regions        = "Australia, Brazil, France, Slovenia, Thailand, Taiwan, Sri Lanka, Myanmar, UK, Switzerland and UK (Cerastes), Martinique, Papua New Guinea.",
    notes          = "Data extracted from text or digitised from concentration-time figures (WebPlotDigitizer v3.12) of 24 published reports listed in Sanhajariya 2018 Table 6; pre-antivenom observations only. 218 timed concentrations total. Most subjects contributed a single sample; only five primary studies reported serial concentrations. The authors describe the model as a preliminary prior for future, richer snake-venom PK datasets (Section 2.2.3 and Section 4)."
  )

  ini({
    # Sanhajariya 2018 Table A1 (Appendix A), 'Covariate Model' column.
    lcl <- log(13.3);  label("Clearance (L/h)")                                                  # Table A1: CL = 13.3 L/h (RSE 14 %)
    lvc <- log(184);   label("Volume of distribution (L)")                                       # Table A1: V  = 184 L (RSE 10 %)

    # Zero-order input duration into the central compartment. Held fixed
    # at 1 h by the authors because the bite event time is approximate
    # and the data do not constrain D1.
    ld1 <- fixed(log(1)); label("Duration of zero-order input D1 into central (h, fixed)")        # Table A1: D1 = 1 h (FIX)

    # Bioavailability anchor for the reference (Viperidae) family is fixed
    # at 1; SNAKEFAMILY_ELAPID supplies a log-multiplicative shift to
    # 0.569 for elapid bites. F1 here is the bioavailability of the
    # input route into the central compartment (no depot compartment
    # exists in the structural model). See vignette 'Assumptions and
    # deviations' for the semantic-stretch note on naming lfdepot.
    lfdepot                     <- fixed(log(1));    label("Bioavailability anchor for Viperidae bites (unitless, fixed reference)")   # Table A1: F1 (Viperidae) = 1 FIX
    e_snakefamily_elapid_fdepot <- log(0.569);       label("Log-multiplicative effect on F1 for Elapidae vs Viperidae (unitless)")     # Table A1: F1 (Elapidae) = 0.569 (RSE 43 %)

    # IIV (CV %) reported in Table A1 'Covariate Model' column. Converted
    # to log-normal variance via omega^2 = log(1 + CV^2).
    etalcl     ~ log(1 + 0.437^2)   # Table A1: BSV CL CV  43.7 % (RSE 52 %)  -> omega^2 = log(1 + 0.437^2) = 0.1762
    etalvc     ~ log(1 + 0.298^2)   # Table A1: BSV V  CV  29.8 % (RSE 101 %) -> omega^2 = log(1 + 0.298^2) = 0.0853
    etald1     ~ log(1 + 0.441^2)   # Table A1: BSV D1 CV  44.1 % (RSE 17 %)  -> omega^2 = log(1 + 0.441^2) = 0.1791
    etalfdepot ~ log(1 + 2.754^2)   # Table A1: BSV F1 CV 275.4 % (RSE 7 %)   -> omega^2 = log(1 + 2.754^2) = 2.1501  (very large; absorbs unknown per-bite venom mass)

    # Residual error -- Table A1 'Proportional error' = 0.047 (RSE 25 %).
    # Interpreted as the NONMEM $SIGMA variance (standard convention --
    # the surrounding BSV rows carry explicit '%CV' formatting while
    # this row does not); the proportional SD is therefore sqrt(0.047)
    # = 0.2168, i.e. ~21.7 % CV. This sidecar interpretation was
    # confirmed by the operator before drafting.
    propSd <- sqrt(0.047); label("Proportional residual error (fraction)")                       # Table A1: SIGMA = 0.047 (variance) -> propSd = sqrt(0.047) = 0.2168
  })

  model({
    # Individual PK parameters
    cl     <- exp(lcl + etalcl)
    vc     <- exp(lvc + etalvc)
    d1     <- exp(ld1 + etald1)
    fdepot <- exp(lfdepot + e_snakefamily_elapid_fdepot * SNAKEFAMILY_ELAPID + etalfdepot)

    kel <- cl / vc

    # One-compartment central. Each snakebite is encoded by the user as a
    # single dose record into central with rate = -2 (modeled duration),
    # so that f(central) and dur(central) apply.
    d/dt(central) <- -kel * central

    f(central)   <- fdepot
    dur(central) <- d1

    # Plasma venom concentration. amt is in mcg, vc in L -> central / vc
    # is in mcg/L (the ELISA reporting unit; Figure 5 axis label).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
