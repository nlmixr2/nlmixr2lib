Wu_2023_SPI_62 <- function() {
  description <- paste(
    "Simultaneous population target-mediated drug disposition (TMDD) PK / PD model for the",
    "small-molecule 11-beta-hydroxysteroid dehydrogenase type 1 (HSD-1) inhibitor SPI-62",
    "(formerly ASP3662) in healthy adults (Wu 2023). Oral absorption through a chain of four",
    "identical first-order transit steps (depot -> transit1 -> transit2 -> transit3 -> central,",
    "all governed by a single Ktr), two-compartment linear disposition (central + peripheral1)",
    "with clearance CL and distribution flow Q, and explicit second-order binding of SPI-62 in",
    "the central compartment to its pharmacological target HSD-1 (association rate constant Kon,",
    "dissociation rate constant Koff) to form a drug-target complex. The total target amount",
    "Rtotal is held constant -- no target synthesis or degradation is included -- so free target",
    "is Rtotal minus the complex amount. This saturable, high-affinity / low-capacity binding is",
    "what produces SPI-62's striking nonlinear PK: extremely low plasma exposure after the first",
    "low doses (most of the dose is trapped on target), turning into linear PK with unusually",
    "high accumulation ratios once the target is saturated by repeated dosing. Hepatic HSD-1",
    "activity, measured as the urinary ratio of (tetrahydrocortisol + allotetrahydrocortisol) to",
    "tetrahydrocortisone expressed as a percentage of each subject's own baseline, is driven",
    "directly (no effect delay) by the central-compartment SPI-62 concentration through an",
    "inhibitory sigmoid Imax function with half-maximal inhibitory concentration IC50 and power",
    "coefficient gamma. Fit by simultaneous PK-PD estimation in NONMEM 7.4.3 (FOCEI, ADVAN13) to",
    "pooled data from the SPI-62 first-in-human single-ascending-dose (1-10 mg) and",
    "multiple-ascending-dose (0.2-2 mg once daily) phase 1 trials. All PK parameters are apparent",
    "(per unit bioavailability) because F is unknown. Exponential inter-individual variability on",
    "Vcentral, CL, Ktr, Koff, Rtotal, and IC50; proportional residual error on both the plasma",
    "concentration and the HSD-1 activity endpoint. No covariates were retained: age, sex, body",
    "weight, and race showed no meaningful relationship with any parameter, so the final model",
    "equals the base model. Amounts are carried internally in nmol so that A / V is directly in",
    "nmol/L (= nM), the units of Kon, Koff, and IC50; convert an mg dose to nmol by dividing by",
    "the SPI-62 molecular weight of 424.4 g/mol and multiplying by 1e6.",
    sep = " "
  )
  reference <- paste(
    "Wu N, Katz DA, An G. Population target-mediated pharmacokinetic/pharmacodynamic modeling to",
    "evaluate SPI-62 exposure and hepatic 11beta-hydroxysteroid dehydrogenase type 1 (HSD-1)",
    "inhibition in healthy adults. Clin Pharmacokinet. 2023;62(9):1275-1288.",
    "doi:10.1007/s40262-023-01278-8. PMID:37452986. PMCID:PMC10449972.",
    sep = " "
  )
  vignette <- "Wu_2023_SPI_62"

  units <- list(
    time          = "hour",
    dosing        = "nmol",
    concentration = "nmol/L",
    dosing_notes  = paste(
      "Amounts are carried internally in nmol and volumes in L, so C = A / V is in nmol/L,",
      "which is identical to nM -- the unit in which Wu 2023 reports Kon (nM^-1 h^-1), Koff",
      "(h^-1), and IC50 (nM). Rtotal is an amount (nmol), matching Table 2. Convert an mg dose",
      "to nmol with nmol = mg * 1e6 / 424.4, using the SPI-62 molecular weight of 424.4 g/mol.",
      "That molecular weight is derived from the paper's own explicit unit conversion in Sect.",
      "3.3, which states IC50 = 0.0787 nM = 0.0334 ng/mL (0.0334 / (0.0787e-3) = 424.4 g/mol);",
      "the paper does not tabulate a molecular weight directly. Example doses: 0.2 mg = 471.3",
      "nmol, 0.4 mg = 942.5 nmol, 0.7 mg = 1649.4 nmol, 2 mg = 4712.6 nmol, 3 mg = 7068.9 nmol,",
      "10 mg = 23562.9 nmol. The HSD-1 activity output is a percentage of each subject's own",
      "pre-dose baseline urinary steroid-metabolite ratio (100% = uninhibited), per Eq. 8.",
      sep = " "
    )
  )

  # No covariates are used by this model. Wu 2023 Sect. 3.2 reports that
  # exploratory analysis found no meaningful impact of age, sex, body weight, or
  # race on any of the six parameters carrying IIV (Vcentral, CL, Ktr, Koff,
  # Rtotal, IC50), so the formal stepwise covariate search was never run and the
  # final PK/PD model is identical to the base PK/PD model.
  covariateData <- list()

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at study entry.",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened graphically against the individual random effects of Vcentral, CL, and IC50",
        "(Wu 2023 ESM Fig. S2; distribution in ESM Fig. S1) but NOT retained: Sect. 3.2 reports",
        "p > 0.05 for all plots, so formal forward-addition / backward-elimination covariate",
        "testing was not performed. Cohort mean +/- SD 35.9 +/- 9.65 years.",
        sep = " "
      ),
      source_name = "age"
    ),
    WT = list(
      description = "Body weight at study entry.",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Screened graphically (Wu 2023 ESM Figs. S1-S2) but NOT retained (Sect. 3.2, p > 0.05).",
        "No allometric scaling is applied in this model. Cohort mean +/- SD 76.5 +/- 12.1 kg.",
        sep = " "
      ),
      source_name = "body weight"
    ),
    SEXF = list(
      description        = "Female sex indicator (1 = female, 0 = male).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = male",
      notes              = paste(
        "Screened graphically as the categorical covariate 'gender' (Wu 2023 ESM Fig. S2) but",
        "NOT retained (Sect. 3.2, p > 0.05). Cohort composition 33 males and 11 females",
        "(25.0% female).",
        sep = " "
      ),
      source_name        = "gender"
    ),
    RACE_WHITE = list(
      description        = "White race indicator (1 = White, 0 = non-White).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = non-White",
      notes              = paste(
        "Screened graphically as the categorical covariate 'race' (Wu 2023 ESM Fig. S2) but NOT",
        "retained (Sect. 3.2, p > 0.05). The paper reports the race composition of the 44",
        "analysis subjects as 26 White, 14 Black, 3 Asian, and 1 Native American; it does not",
        "state how the four levels were grouped for the exploratory boxplots, and no effect",
        "coefficient is published for any level. The single Native American subject has no",
        "canonical indicator column in inst/references/covariate-columns.md and none is",
        "introduced here, because no race term enters model().",
        sep = " "
      ),
      source_name        = "race"
    ),
    RACE_BLACK = list(
      description        = "Black / African American race indicator (1 = Black, 0 = other).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = other",
      notes              = "Screened but NOT retained; see the RACE_WHITE note. 14 of 44 subjects.",
      source_name        = "race"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator (1 = Asian, 0 = other).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = other",
      notes              = "Screened but NOT retained; see the RACE_WHITE note. 3 of 44 subjects.",
      source_name        = "race"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 44L,
    n_studies      = 2L,
    study_names    = c(
      "SAD -- SPI-62 (ASP3662) first-in-human single ascending dose in healthy adults; 24 of 48 randomised subjects contributed (the 1, 3, 6, and 10 mg cohorts; n = 6 active per cohort)",
      "MAD -- SPI-62 (ASP3662) multiple ascending dose in healthy adults; 20 subjects from the low-dose cohorts only (0.2, 0.4, 0.7, and 2 mg)"
    ),
    age_range      = "adults; mean +/- SD 35.9 +/- 9.65 years (individual min-max not tabulated; distribution in ESM Fig. S1).",
    weight_range   = "mean +/- SD 76.5 +/- 12.1 kg (individual min-max not tabulated; distribution in ESM Fig. S1).",
    sex_female_pct = 25.0,
    sex_notes      = "33 males and 11 females among the 44 analysis subjects (Wu 2023 Sect. 2.2.1).",
    race_ethnicity = c(White = 59.1, Black = 31.8, Asian = 6.8, `Native American` = 2.3),
    race_notes     = "Counts as reported in Wu 2023 Sect. 2.2.1: 26 White, 14 Black, 3 Asian, 1 Native American (percentages of 44).",
    disease_state  = "Healthy adult volunteers.",
    dose_range     = paste(
      "SAD: 1, 3, 6, and 10 mg single oral dose (the 30 and 60 mg cohorts exist in the trial but",
      "were excluded from this analysis). MAD low-dose arms: 3 mg loading dose on day 1 then 0.2",
      "mg once daily on days 2-14; 0.4 mg once daily on days 1-14; and 0.7 or 2 mg single dose on",
      "day 1 followed by a 6-day washout then once-daily dosing on days 7-20. The MAD high-dose",
      "arms (10, 20, 50 mg) were deliberately excluded -- Sect. 2.1 reports that including the",
      ">= 20 mg groups degraded the fit at the low doses that are the intended clinical range.",
      sep = " "
    ),
    regions        = "Not reported.",
    data_records   = paste(
      "996 SPI-62 plasma concentrations (222 below the limit of quantification, imputed as",
      "LLOQ/2) and 279 observed baseline-corrected hepatic HSD-1 activity values (one outlier",
      "excluded). LLOQ was 0.1 ng/mL in the SAD trial and 4 pg/mL in the MAD trial.",
      sep = " "
    ),
    notes          = paste(
      "Secondary analysis of the two SPI-62 phase 1 trials; the underlying datasets are",
      "proprietary and were not released. PK and PD data were fit simultaneously. The PK model",
      "structure is the same two-compartment-plus-three-transit TMDD model the same group",
      "previously published from the PK data alone; the parameters here are re-estimated from the",
      "joint PK + PD fit and Sect. 3.3 notes they are very close to the PK-only values",
      "(e.g. Rtotal 5460 nmol here vs 6070 nmol PK-only).",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Absorption -- four identical sequential first-order transit steps
    # (depot -> transit1 -> transit2 -> transit3 -> central), all sharing
    # a single rate constant Ktr. Wu 2023 Eqs. 1-5 and Fig. 3.
    # ------------------------------------------------------------------
    lktr <- log(8.82); label("Transit absorption rate constant Ktr (1/h)")            # Wu 2023 Table 2 (Ktr = 8.82 1/h, RSE 11%)

    # ------------------------------------------------------------------
    # Two-compartment linear disposition. Table 2 footnote a: these are
    # apparent parameters because bioavailability F is unknown.
    # ------------------------------------------------------------------
    lcl  <- log(10.1); label("Apparent clearance CL/F (L/h)")                          # Wu 2023 Table 2 (CL = 10.1 L/h, RSE 6%; apparent per Table 2 footnote a)
    lvc  <- log(152);  label("Apparent central volume of distribution Vcen/F (L)")      # Wu 2023 Table 2 (V_cen = 152 L, RSE 16%; apparent per Table 2 footnote a)
    lq   <- log(2.38); label("Apparent distribution flow Q/F (L/h)")                    # Wu 2023 Table 2 (Q = 2.38 L/h, RSE 12%; apparent per Table 2 footnote a)
    lvp  <- log(116);  label("Apparent peripheral volume of distribution Vperi/F (L)")  # Wu 2023 Table 2 (V_peri = 116 L, RSE 6%; apparent per Table 2 footnote a)

    # Bioavailability is not identifiable (all L and L/h above are apparent),
    # so F is anchored at 1 rather than estimated. Wu 2023 Eq. 1 writes the
    # depot initial condition as Dose * F.
    lfdepot <- fixed(log(1)); label("Bioavailability F of SPI-62 (fraction, fixed at 1)")  # Wu 2023 Eq. 1 and Table 2 footnote a (F unknown; only apparent parameters estimated, so F is anchored at 1)

    # ------------------------------------------------------------------
    # Target binding (TMDD). Concentrations are nM and Rtotal is an
    # amount in nmol, so Kon * C * (Rtotal - RC) has units nmol/h.
    # Wu 2023 Eqs. 5 and 7. The total target amount is constant: no
    # synthesis or degradation of HSD-1 is included (Sect. 2.2).
    # ------------------------------------------------------------------
    lkon  <- log(8.43);  label("Second-order association rate constant Kon (1/(nM*h))")   # Wu 2023 Table 2 (Kon = 8.43 nM^-1 h^-1, RSE 6%)
    lkoff <- log(0.229); label("First-order dissociation rate constant Koff (1/h)")       # Wu 2023 Table 2 (Koff = 0.229 1/h, RSE 31%); implied Kd = Koff/Kon = 0.0272 nM (Sect. 3.3)
    lrtot <- log(5460);  label("Total amount of HSD-1 target Rtotal (nmol)")              # Wu 2023 Table 2 (R_total = 5460 nmol, RSE 7%); Sect. 3.3 notes this corresponds to roughly 2.2 mg of SPI-62

    # ------------------------------------------------------------------
    # PD -- inhibitory sigmoid Imax model on hepatic HSD-1 activity,
    # driven directly by the central-compartment concentration with no
    # effect delay (Sect. 3.1: no lag between Cmax and maximal
    # inhibition). Wu 2023 Eq. 10.
    #
    # Table 2 reports Imax in percent (99.9%); Eq. 10 consumes it as a
    # fraction inside 100 * (1 - Imax * C^gamma / (IC50^gamma + C^gamma)),
    # so the stored value is 0.999. Imax was estimated, not fixed --
    # ESM Table S1 row "13-final" lists "Imax = theta" (no IIV, not
    # fixed) and Table 2 gives it an RSE of 2%.
    # ------------------------------------------------------------------
    limax <- log(0.999);  label("Maximum fractional inhibition Imax of hepatic HSD-1 activity (fraction)")  # Wu 2023 Table 2 (Imax = 99.9%, RSE 2%), expressed as the fraction Eq. 10 consumes
    lec50 <- log(0.0787); label("SPI-62 concentration at 50% of Imax, IC50 (nmol/L)")                        # Wu 2023 Table 2 (IC50 = 0.0787 nM, RSE 16%); Sect. 3.3 gives the equivalent 0.0334 ng/mL and notes it is close to Kd = 0.0272 nM
    lhill <- log(0.441);  label("Sigmoidicity (power) coefficient gamma of the Imax model (unitless)")       # Wu 2023 Table 2 (gamma = 0.441, RSE 7%); the paper's "power coefficient r" in Eq. 10

    # ------------------------------------------------------------------
    # Inter-individual variability. Wu 2023 Sect. 2.2.1: "IIV was
    # estimated using an exponential model ... with a mean of 0 and a
    # variance of omega^2". Table 2 reports each IIV as a percent
    # coefficient of variation, so the internal variance is
    # omega^2 = log(CV^2 + 1) (the log-normal identity). Diagonal only
    # -- Wu 2023 reports no off-diagonal covariances.
    #   Vcentral CV 51.9% -> log(1 + 0.519^2) = 0.238514
    #   CL       CV 19.2% -> log(1 + 0.192^2) = 0.036201
    #   Ktr      CV 50.2% -> log(1 + 0.502^2) = 0.224745
    #   Koff     CV  119% -> log(1 + 1.19^2)  = 0.882155
    #   Rtotal   CV 27.5% -> log(1 + 0.275^2) = 0.072902
    #   IC50     CV 39.5% -> log(1 + 0.395^2) = 0.144987
    # ------------------------------------------------------------------
    etalvc   ~ 0.238514   # Wu 2023 Table 2 (IIV Vcentral = 51.9%, RSE 33%, shrinkage 14%); omega^2 = log(CV^2 + 1)
    etalcl   ~ 0.036201   # Wu 2023 Table 2 (IIV CL = 19.2%, RSE 84%, shrinkage 21%); omega^2 = log(CV^2 + 1)
    etalktr  ~ 0.224745   # Wu 2023 Table 2 (IIV Ktr = 50.2%, RSE 38%, shrinkage 7%); omega^2 = log(CV^2 + 1)
    etalkoff ~ 0.882155   # Wu 2023 Table 2 (IIV Koff = 119%, RSE 55%, shrinkage 15%); omega^2 = log(CV^2 + 1)
    etalrtot ~ 0.072902   # Wu 2023 Table 2 (IIV Rtotal = 27.5%, RSE 31%, shrinkage 9%); omega^2 = log(CV^2 + 1)
    etalec50 ~ 0.144987   # Wu 2023 Table 2 (IIV IC50 = 39.5%, RSE 63%, shrinkage 37%); omega^2 = log(CV^2 + 1)

    # ------------------------------------------------------------------
    # Residual error -- proportional on both endpoints. A combined
    # proportional + additive and a pure additive model were both tested
    # (Sect. 2.2.1; ESM Table S1 row 14 tried additive residual error and
    # was 80.782 OFV units worse), and proportional was retained.
    # ------------------------------------------------------------------
    propSd                 <- 0.274; label("Proportional residual error on SPI-62 plasma concentration (fraction)")  # Wu 2023 Table 2 (sigma_PK = 27.4%, RSE 3%, shrinkage 7%)
    propSd_hsd1activity    <- 0.193; label("Proportional residual error on hepatic HSD-1 activity (fraction)")        # Wu 2023 Table 2 (sigma_PD = 19.3%, RSE 12%, shrinkage 5%)
  })

  model({
    # ------------------------------------------------------------------
    # Unit convention (see units$dosing_notes):
    #   amounts   nmol      volumes   L      time   h
    #   => A / V is nmol/L == nM, matching Kon, Koff, and IC50.
    # Doses passed to rxode2::et() must therefore be in nmol; convert
    # from mg with nmol = mg * 1e6 / 424.4, where 424.4 g/mol is the
    # SPI-62 molecular weight implied by the paper's own conversion
    # IC50 = 0.0787 nM = 0.0334 ng/mL (Wu 2023 Sect. 3.3).
    # ------------------------------------------------------------------

    # 1. Individual parameters. Exponential (log-normal) IIV on Vcentral,
    #    CL, Ktr, Koff, Rtotal, and IC50; none on Q, Vperi, Kon, Imax, or
    #    gamma (Wu 2023 Table 2).
    ktr    <- exp(lktr + etalktr)
    cl     <- exp(lcl + etalcl)
    vc     <- exp(lvc + etalvc)
    q      <- exp(lq)
    vp     <- exp(lvp)
    fdepot <- exp(lfdepot)

    kon    <- exp(lkon)
    koff   <- exp(lkoff + etalkoff)
    rtot   <- exp(lrtot + etalrtot)

    imax   <- exp(limax)
    ec50   <- exp(lec50 + etalec50)
    hill   <- exp(lhill)

    # 2. Central-compartment SPI-62 concentration (nM). This single
    #    quantity drives the target binding (Eqs. 5, 7) and the PD
    #    (Eq. 10).
    Cc <- central / vc

    # 3. Target binding flux (nmol/h). Free target is Rtotal - complex
    #    because Rtotal is constant: Wu 2023 Sect. 2.2 states "The total
    #    amount of HSD-1 (R_total) was assumed to be a constant. No
    #    target synthesis or degradation processes were included".
    bind_rate   <- kon * Cc * (rtot - complex)
    unbind_rate <- koff * complex

    # 4. ODE system. Wu 2023 Eqs. 1-7. All initial conditions are 0 (the
    #    rxode2 default); the dose enters `depot` as Dose * F per Eq. 1.
    d/dt(depot)       <- -ktr * depot                                        # Eq. 1
    d/dt(transit1)    <-  ktr * depot    - ktr * transit1                    # Eq. 2
    d/dt(transit2)    <-  ktr * transit1 - ktr * transit2                    # Eq. 3
    d/dt(transit3)    <-  ktr * transit2 - ktr * transit3                    # Eq. 4
    d/dt(central)     <-  ktr * transit3 - bind_rate + unbind_rate -
                          (cl / vc) * central -
                          (q / vc) * central + (q / vp) * peripheral1        # Eq. 5
    d/dt(peripheral1) <-  (q / vc) * central - (q / vp) * peripheral1        # Eq. 6
    d/dt(complex)     <-  bind_rate - unbind_rate                            # Eq. 7

    # 5. Bioavailability (F = 1 fixed; Eq. 1 initial condition Dose * F).
    f(depot) <- fdepot

    # 6. Observations.
    #
    #    Cc            -- SPI-62 plasma concentration (nM).
    #    hsd1activity  -- hepatic HSD-1 activity as a percentage of the
    #                     subject's own pre-dose baseline urinary
    #                     (tetrahydrocortisol + allotetrahydrocortisol) /
    #                     tetrahydrocortisone ratio (Eq. 8), predicted by
    #                     the inhibitory sigmoid Imax model of Eq. 10:
    #                       100 * (1 - Imax * C^gamma / (IC50^gamma + C^gamma))
    #                     100% is uninhibited; hepatic HSD-1 inhibition
    #                     in percent is 100 - hsd1activity (Eq. 9).
    #
    #    `max(Cc, 0.0)` is a numerical guard only, not a change to Eq. 10:
    #    gamma = 0.441 is a fractional power, so a transient slightly
    #    negative Cc from ODE solver error -- a real risk here because the
    #    concentrations of interest are in the pg/mL range -- would make
    #    C^gamma NaN. It has no effect on any non-negative Cc.
    cpos         <- max(Cc, 0.0)
    hsd1activity <- 100 * (1 - imax * cpos^hill / (ec50^hill + cpos^hill))   # Eq. 10

    Cc           ~ prop(propSd)
    hsd1activity ~ prop(propSd_hsd1activity)
  })
}
