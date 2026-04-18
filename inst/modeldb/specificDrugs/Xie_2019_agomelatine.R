# Load the xie_2019_agomelatine model from the paper supplement

Xie_2019_agomelatine <- function () {
  description <- "A semiphysiological population pharmacokinetic model of agomelatine and its metabolites in Chinese healthy volunteers"
  reference <- "Xie F, Vermeulen A, Colin P, Cheng Z. A semiphysiological population pharmacokinetic model of agomelatine and its metabolites in Chinese healthy volunteers. Br J Clin Pharmacol. 2019 May;85(5):1003-1014. doi: 10.1111/bcp.13902. Epub 2019 Mar 21. PMID: 30761579; PMCID: PMC6475681."
  units <-
    list(
      time = "hr",
      dosing = "mg",
      concentration = "ng/mL" # applied to all three plasma outputs: calmt (agomelatine), c3oh (3-hydroxy-agomelatine), c7dm (7-desmethyl-agomelatine)
    )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used inside liver-volume allometry: lv = 0.05012 * WT^0.78. No explicit reference weight reported; the allometric form uses the raw WT value directly.",
      source_name        = "WT"
    ),
    ooc1 = list(
      description        = "Occasion indicator for period 1 of the four-period crossover study",
      units              = "(binary)",
      type               = "binary",
      reference_category = "Not applicable; ooc1..ooc4 are a mutually exclusive set (exactly one is 1 per observation)",
      notes              = "Lower-case name preserved from source per covariate-columns.md register. Used to select the period-specific IOV eta across k13, alag2, k23, clint, and the logit-fraction partitioning absorption between DEPOT1 and DEPOT2.",
      source_name        = "ooc1"
    ),
    ooc2 = list(
      description        = "Occasion indicator for period 2 of the four-period crossover study",
      units              = "(binary)",
      type               = "binary",
      reference_category = "Not applicable; ooc1..ooc4 are a mutually exclusive set (exactly one is 1 per observation)",
      notes              = "Lower-case name preserved from source per covariate-columns.md register.",
      source_name        = "ooc2"
    ),
    ooc3 = list(
      description        = "Occasion indicator for period 3 of the four-period crossover study",
      units              = "(binary)",
      type               = "binary",
      reference_category = "Not applicable; ooc1..ooc4 are a mutually exclusive set (exactly one is 1 per observation)",
      notes              = "Lower-case name preserved from source per covariate-columns.md register.",
      source_name        = "ooc3"
    ),
    ooc4 = list(
      description        = "Occasion indicator for period 4 of the four-period crossover study",
      units              = "(binary)",
      type               = "binary",
      reference_category = "Not applicable; ooc1..ooc4 are a mutually exclusive set (exactly one is 1 per observation)",
      notes              = "Lower-case name preserved from source per covariate-columns.md register.",
      source_name        = "ooc4"
    )
  )

  population <- list(
    n_subjects     = "TODO: from source paper",
    n_studies      = 1,
    age_range      = "TODO: from source paper",
    age_median     = "TODO: from source paper",
    weight_range   = "TODO: from source paper",
    weight_median  = "60 kg (used as representative weight in vignette simulations per Table 1 of source)",
    sex_female_pct = "TODO: from source paper",
    race_ethnicity = c(Asian = 100),
    disease_state  = "Healthy Chinese volunteers",
    dose_range     = "25 mg single oral dose (vignette simulation; confirm full design in source)",
    regions        = "China",
    notes          = "Four-period crossover design (occasions ooc1..ooc4) with IOV on multiple PK parameters. TODO: fill exact demographics from Table 1 of Xie 2019."
  )

  ini({
    ltvk13 <- log(4.54); label("K13 (1/h)")
    ltvv4 <- log(64.6); label("V4 (L)")
    ltvclint <- log(111000); label("CLint (L/h)")
    lBA3 <- log(exp(-1.78)); label("BA3") # log(exp()) so that the original parameter value is maintained in the ini block
    lBA4 <- log(exp(-3.95)); label("BA4") # log(exp()) so that the original parameter value is maintained in the ini block
    ltvcl3oh <- log(44.9); label("CL3OH (L/h)")
    ltvcl7dm <- log(52.9); label("CL7DM (L/h)")
    ltvalag1 <- log(0.185); label("ALAG1 (h)")
    ltvk23 <- log(4.23); label("K23 (1/h)")
    ltvalag2 <- log(0.305); label("ALAG2 (h)")
    F1 <- 0.681; label("F1")
    lvq7dm <- log(28.1); label("Q7DM (L/h)")
    lvv7dm <- log(536); label("V7DM (L)")
    lvqalmt <- log(4.36); label("QALMT (L/h)")
    lvvalmt <- log(157); label("VALMT (L)")
    sdalmt <- 0.39; label("Residual standard deviation of agomelatine (log-scale, additive error)")
    sd3oh <- 0.228; label("Residual standard deviation of 3-hydroxy-agomelatine (log-scale, additive error)")
    sd7dm <- 0.297; label("SD7DM")
    etaltvk13 ~ 0.101
    etaltvv4 ~ 0.0374
    etaltvclint ~ 0.998
    etalBA3 ~ 0.175
    etalBA4 ~ c(0.153, 0.231)
    etaltvcl3oh ~ 0.0414
    etaltvcl7dm ~ c(0.0492, 0.0774)
    etaltvalag1 ~ 0.0628
    etaltvk23 ~ 3.78
    etaltvalag2 ~ 3.3
    etaF1 ~ 0.526
    etalvq7dm ~ fix(0.001)
    etalvv7dm ~ fix(0.001)
    etalvqalmt ~ 1.46
    etalvvalmt ~ fix(0.001)
    e.IOV1 ~ 1.52
    eta17 ~ fix(1.52)
    eta18 ~ fix(1.52)
    eta19 ~ fix(1.52)
    e.IOV2 ~ 4.32
    eta21 ~ fix(4.32)
    eta22 ~ fix(4.32)
    eta23 ~ fix(4.32)
    e.IOV3 ~ 5.01
    eta25 ~ fix(5.01)
    eta26 ~ fix(5.01)
    eta27 ~ fix(5.01)
    e.IOV4 ~ 0.0779
    eta29 ~ fix(0.0779)
    eta30 ~ fix(0.0779)
    eta31 ~ fix(0.0779)
    e.IOV5 ~ 2.32
    eta33 ~ fix(2.32)
    eta34 ~ fix(2.32)
    eta35 ~ fix(2.32)
  })
  model({
    iov1 <- ooc1 * e.IOV1 + ooc2 * eta17 + ooc3 * eta18 + ooc4 * eta19
    iov2 <- ooc1 * e.IOV2 + ooc2 * eta21 + ooc3 * eta22 + ooc4 * eta23
    iov3 <- ooc1 * e.IOV3 + ooc2 * eta25 + ooc3 * eta26 + ooc4 * eta27
    iov4 <- ooc1 * e.IOV4 + ooc2 * eta29 + ooc3 * eta30 + ooc4 * eta31
    iov5 <- ooc1 * e.IOV5 + ooc2 * eta33 + ooc3 * eta34 + ooc4 * eta35

    k13 <- exp(ltvk13 + etaltvk13) * exp(iov1)
    v4 <- exp(ltvv4 + etaltvv4)
    clint <- exp(ltvclint + etaltvclint + iov4)
    alag1 <- exp(ltvalag1 + etaltvalag1)
    k23 <- exp(ltvk23 + etaltvk23) * exp(iov3)
    alag2 <- exp(ltvalag2 + etaltvalag2 + iov2)

    expp <- log(F1/(1 - F1)) + etaF1
    fDepot1 <- exp(expp + iov5)/(1 + exp(expp + iov5))
    fDepot2 <- 1 - fDepot1
    lv <- 0.05012 * WT^0.78
    v3 <- lv
    pbr <- 1/0.69
    qh <- 50.4 * lv
    FU <- 0.05
    eh <- FU * pbr * clint/(qh + FU * pbr * clint)
    fh <- 1 - eh
    clh <- qh * eh/pbr
    cl <- clh
    ba3 <- exp(lBA3 + etalBA3)
    ba4 <- exp(lBA4 + etalBA4)
    fm3oh <- ba3/(1 + ba3 + ba4)
    fm7dm <- ba4/(1 + ba3 + ba4)
    cl3oh <- exp(ltvcl3oh + etaltvcl3oh)
    cl7dm <- exp(ltvcl7dm + etaltvcl7dm)
    q7dm <- exp(lvq7dm + etalvq7dm)
    v7dm <- exp(lvv7dm + etalvv7dm)
    qalmt <- exp(lvqalmt + etalvqalmt)
    valmt <- exp(lvvalmt + etalvvalmt)
    mpr1 <- 259/243
    mpr2 <- 229/243
    v5 <- v4
    v6 <- v4

    d/dt(DEPOT1) <- -k13 * DEPOT1
    d/dt(DEPOT2) <- -k23 * DEPOT2
    d/dt(LIVER) <- k13 * DEPOT1 + k23 * DEPOT2 - qh * fh * LIVER/v3 + qh * CENTPRNT/v4 - clh * LIVER/v3
    d/dt(CENTPRNT) <- qh * fh * LIVER/v3 - qh * CENTPRNT/v4 - CENTPRNT * qalmt/v4 + ALMTPERI * qalmt/valmt
    d/dt(METB3OH) <- fm3oh * clh * LIVER/v3 * mpr1 - METB3OH * cl3oh/v5
    d/dt(METB7DM) <- fm7dm * clh * LIVER/v3 * mpr2 - METB7DM * cl7dm/v6 - METB7DM * q7dm/v6 + METB7DMPERI * q7dm/v7dm
    d/dt(METB7DMPERI) <- METB7DM * q7dm/v6 - METB7DMPERI * q7dm/v7dm
    d/dt(ALMTPERI) <- CENTPRNT * qalmt/v4 - ALMTPERI * qalmt/valmt

    alag(DEPOT1) <- alag1
    alag(DEPOT2) <- alag2
    f(DEPOT1) <- fDepot1
    f(DEPOT2) <- fDepot2

    # Concentration of agomelatine in plasma. The unit conversion of *1000 was
    # not in the original model. It is used to allow dosing units to be in mg and
    calmt <- CENTPRNT/v4 * 1000
    lcalmt <- log(calmt)
    # Concentration of 3-hydroxy-agomelatine in plasma. The unit conversion of *1000 was
    # not in the original model. It is used to allow dosing units to be in mg and
    c3oh <- METB3OH/v5 * 1000
    lc3oh <- log(c3oh)
    # Concentration of 7-desmethyl-agomelatine in plasma. The unit conversion of *1000 was
    # not in the original model. It is used to allow dosing units to be in mg and
    c7dm <- METB7DM/v6 * 1000
    lc7dm <- log(c7dm)

    lloqalmt <- log(0.0457)
    lloq7dm <- log(0.1372)

    lcalmt ~ add(sdalmt)
    lc3oh ~ add(sd3oh)
    lc7dm ~ add(sd7dm)
  })
}
