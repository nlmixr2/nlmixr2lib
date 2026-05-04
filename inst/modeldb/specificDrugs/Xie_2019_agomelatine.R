# Load the xie_2019_agomelatine model from the paper supplement

Xie_2019_agomelatine <- function () {
  description <- "A semiphysiological population pharmacokinetic model of agomelatine and its metabolites in Chinese healthy volunteers"
  reference <- "Xie F, Vermeulen A, Colin P, Cheng Z. A semiphysiological population pharmacokinetic model of agomelatine and its metabolites in Chinese healthy volunteers. Br J Clin Pharmacol. 2019 May;85(5):1003-1014. doi: 10.1111/bcp.13902. Epub 2019 Mar 21. PMID: 30761579; PMCID: PMC6475681."
  vignette <- "Xie_2019_agomelatine"
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
      notes              = "Lower-case name preserved from source per covariate-columns.md register. Used to select the period-specific IOV eta across k13, alag2, k23, clint, and the logit-fraction partitioning absorption between depot and depot2.",
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
    ltvalag <- log(0.185); label("ALAG1 (h)")
    ltvk23 <- log(4.23); label("K23 (1/h)")
    ltvalag2 <- log(0.305); label("ALAG2 (h)")
    fpop <- 0.681; label("Typical population fraction of dose absorbed via the primary depot (F1; unitless)")
    lvq7dm <- log(28.1); label("Q7DM (L/h)")
    lvv7dm <- log(536); label("V7DM (L)")
    lvqalmt <- log(4.36); label("QALMT (L/h)")
    lvvalmt <- log(157); label("VALMT (L)")
    addSd_lcalmt <- 0.39;  label("Additive residual SD of log-agomelatine plasma concentration (log-ng/mL)")
    addSd_lc3oh  <- 0.228; label("Additive residual SD of log-3-hydroxy-agomelatine plasma concentration (log-ng/mL)")
    addSd_lc7dm  <- 0.297; label("Additive residual SD of log-7-desmethyl-agomelatine plasma concentration (log-ng/mL)")
    etaltvk13 ~ 0.101
    etaltvv4 ~ 0.0374
    etaltvclint ~ 0.998
    etalBA3 ~ 0.175
    etalBA4 ~ c(0.153, 0.231)
    etaltvcl3oh ~ 0.0414
    etaltvcl7dm ~ c(0.0492, 0.0774)
    etaltvalag ~ 0.0628
    etaltvk23 ~ 3.78
    etaltvalag2 ~ 3.3
    etafpop ~ 0.526
    etalvq7dm ~ fix(0.001)
    etalvv7dm ~ fix(0.001)
    etalvqalmt ~ 1.46
    etalvvalmt ~ fix(0.001)
    # Inter-occasion variability (IOV). Each chain is a single estimated
    # variance (occasion 1) propagated as fix(...) through occasions 2-4 so
    # the four occasions share a common between-occasion variance. The eta
    # name encodes which structural parameter the IOV applies to and the
    # occasion index: etaiov_<param>_<occ>.
    etaiov_k13_1 ~ 1.52
    etaiov_k13_2 ~ fix(1.52)
    etaiov_k13_3 ~ fix(1.52)
    etaiov_k13_4 ~ fix(1.52)
    etaiov_alag2_1 ~ 4.32
    etaiov_alag2_2 ~ fix(4.32)
    etaiov_alag2_3 ~ fix(4.32)
    etaiov_alag2_4 ~ fix(4.32)
    etaiov_k23_1 ~ 5.01
    etaiov_k23_2 ~ fix(5.01)
    etaiov_k23_3 ~ fix(5.01)
    etaiov_k23_4 ~ fix(5.01)
    etaiov_clint_1 ~ 0.0779
    etaiov_clint_2 ~ fix(0.0779)
    etaiov_clint_3 ~ fix(0.0779)
    etaiov_clint_4 ~ fix(0.0779)
    etaiov_fpop_1 ~ 2.32
    etaiov_fpop_2 ~ fix(2.32)
    etaiov_fpop_3 ~ fix(2.32)
    etaiov_fpop_4 ~ fix(2.32)
  })
  model({
    iov_k13   <- ooc1 * etaiov_k13_1   + ooc2 * etaiov_k13_2   + ooc3 * etaiov_k13_3   + ooc4 * etaiov_k13_4
    iov_alag2 <- ooc1 * etaiov_alag2_1 + ooc2 * etaiov_alag2_2 + ooc3 * etaiov_alag2_3 + ooc4 * etaiov_alag2_4
    iov_k23   <- ooc1 * etaiov_k23_1   + ooc2 * etaiov_k23_2   + ooc3 * etaiov_k23_3   + ooc4 * etaiov_k23_4
    iov_clint <- ooc1 * etaiov_clint_1 + ooc2 * etaiov_clint_2 + ooc3 * etaiov_clint_3 + ooc4 * etaiov_clint_4
    iov_fpop    <- ooc1 * etaiov_fpop_1    + ooc2 * etaiov_fpop_2    + ooc3 * etaiov_fpop_3    + ooc4 * etaiov_fpop_4

    k13 <- exp(ltvk13 + etaltvk13) * exp(iov_k13)
    v4 <- exp(ltvv4 + etaltvv4)
    clint <- exp(ltvclint + etaltvclint + iov_clint)
    alag <- exp(ltvalag + etaltvalag)
    k23 <- exp(ltvk23 + etaltvk23) * exp(iov_k23)
    alag2 <- exp(ltvalag2 + etaltvalag2 + iov_alag2)

    expp <- log(fpop/(1 - fpop)) + etafpop
    fDepot <- exp(expp + iov_fpop)/(1 + exp(expp + iov_fpop))
    fDepot2 <- 1 - fDepot
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
    mpr_3oh <- 259/243
    mpr_7dm <- 229/243
    v5 <- v4
    v6 <- v4

    d/dt(depot)           <- -k13 * depot
    d/dt(depot2)          <- -k23 * depot2
    d/dt(liver)           <- k13 * depot + k23 * depot2 - qh * fh * liver/v3 + qh * central/v4 - clh * liver/v3
    d/dt(central)         <- qh * fh * liver/v3 - qh * central/v4 - central * qalmt/v4 + peripheral1 * qalmt/valmt
    d/dt(central_3oh)     <- fm3oh * clh * liver/v3 * mpr_3oh - central_3oh * cl3oh/v5
    d/dt(central_7dm)     <- fm7dm * clh * liver/v3 * mpr_7dm - central_7dm * cl7dm/v6 - central_7dm * q7dm/v6 + peripheral1_7dm * q7dm/v7dm
    d/dt(peripheral1_7dm) <- central_7dm * q7dm/v6 - peripheral1_7dm * q7dm/v7dm
    d/dt(peripheral1)     <- central * qalmt/v4 - peripheral1 * qalmt/valmt

    alag(depot)  <- alag
    alag(depot2) <- alag2
    f(depot)     <- fDepot
    f(depot2)    <- fDepot2

    # Concentration of agomelatine in plasma. The unit conversion of *1000 was
    # not in the original model. It is used to allow dosing units to be in mg and
    calmt <- central/v4 * 1000
    lcalmt <- log(calmt)
    # Concentration of 3-hydroxy-agomelatine in plasma. The unit conversion of *1000 was
    # not in the original model. It is used to allow dosing units to be in mg and
    c3oh <- central_3oh/v5 * 1000
    lc3oh <- log(c3oh)
    # Concentration of 7-desmethyl-agomelatine in plasma. The unit conversion of *1000 was
    # not in the original model. It is used to allow dosing units to be in mg and
    c7dm <- central_7dm/v6 * 1000
    lc7dm <- log(c7dm)

    lloqalmt <- log(0.0457)
    lloq7dm <- log(0.1372)

    lcalmt ~ add(addSd_lcalmt)
    lc3oh  ~ add(addSd_lc3oh)
    lc7dm  ~ add(addSd_lc7dm)
  })
}
