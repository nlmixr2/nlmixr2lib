Wang_2023_zoledronicAcid_mbma <- function() {
  description <- paste0(
    "MBMA. Kinetic-pharmacodynamic (K-PD) model of lumbar-spine bone mineral ",
    "density (BMD) for intravenous zoledronic acid in primary osteoporosis, ",
    "built on aggregate BMD-time data digitised from 10 published randomised ",
    "trials (6,014 patients, 1-5 mg single and once-yearly dosing, follow-up ",
    "12-72 months). The administered dose lands in a virtual K-PD amount ",
    "compartment (depot_kpd, ng) that is eliminated first-order at rate kel ",
    "(the paper's KDE); the resulting virtual dose-driving rate ir = kel * ",
    "depot_kpd (ng/month) stimulates the zero-order BMD synthesis rate ",
    "kin = BMD_BL * kout through a unit-Emax factor (1 + ir / (ekd50 + ir)), ",
    "so a saturating dose can at most double BMD synthesis. The stimulated ",
    "synthesis rate is multiplied by a first-order desensitisation (tolerance) ",
    "factor exp(-kdes * t) on absolute study time, which the authors added to ",
    "reproduce the attenuation of BMD gain observed over repeated annual doses ",
    "(attributed to loss of osteoclast sensitivity); BMD is lost first-order at ",
    "kout. The single random effect is between-study variability on ekd50; ",
    "residual error is additive on BMD. Because the analysis units are ",
    "published trial arms rather than individuals, this model simulates ",
    "arm-level mean BMD trajectories, not individual patient BMD. The paper's ",
    "companion exposure-response analysis of the acute phase reaction is ",
    "descriptive (observed incidence by dose group) and contributes no ",
    "estimated parameters, so it is not encoded here."
  )
  reference <- paste(
    "Wang H, Liu Q, Jiang M, Song C, Liu D (2023).",
    "Optimization of the dosage regimen of zoledronic acid with a",
    "kinetic-pharmacodynamic model and exposure-response analysis.",
    "Front Pharmacol 14:1089774.",
    "doi:10.3389/fphar.2023.1089774.",
    sep = " "
  )
  vignette <- "Wang_2023_zoledronicAcid"
  units <- list(
    time          = "month",
    dosing        = "ng",
    concentration = "g/cm^2"
    # units$dosing is ng, NOT mg. The virtual K-PD amount compartment must
    # carry the dose in nanograms for the published EDK50 of 41300 ng/month
    # (Wang 2023 Table 2) to be dimensionally consistent with
    # ir = kel * depot_kpd: a 5 mg infusion is amt = 5e6, giving
    # ir(0) = 0.08150 * 5e6 = 4.075e5 ng/month and a stimulation fraction
    # ir / (ekd50 + ir) = 0.908. Entering the dose in mg instead would make
    # the drug effect vanish (fraction ~1e-5). units$concentration documents
    # the single model output: lumbar-spine areal BMD in g/cm^2; this is a
    # K-PD model, so there is no observed plasma zoledronic-acid
    # concentration anywhere in it.
  )

  compartmentData <- list(
    depot_kpd = list(
      analyte  = "zoledronic acid",
      units    = "ng",
      specimen = "administration site",
      verified = TRUE
    ),
    BMD_LS = list(
      analyte  = "bone mineral density (lumbar spine)",
      units    = "g/cm^2",
      specimen = "tissue",
      verified = TRUE
    )
  )

  covariateData <- list(
    BMD_BL = list(
      description        = paste(
        "Pre-treatment lumbar-spine (vertebral) bone mineral density measured by",
        "dual-energy X-ray absorptiometry (DXA). Time-fixed per study arm. Sets",
        "the BMD state initial condition (`BMD_LS(0) <- BMD_BL`) and, through",
        "Wang 2023 Equation 4 (KS = BASE * KD), the zero-order BMD synthesis",
        "rate, so that an undosed arm holds exactly at its own baseline (up to",
        "the desensitisation term)."
      ),
      units              = "g/cm^2",
      type               = "continuous",
      reference_category = "n/a -- per-arm anchor; no covariate coefficient is estimated on it",
      notes              = paste(
        "Wang 2023 Table 1 reports the per-study baseline vertebral BMD for each",
        "of the 10 included trials: 1.03, 1.03, 0.81, 0.79, 0.64, 0.93, 0.75,",
        "0.66, 1.06 and 1.03 g/cm^2 (range 0.64-1.06). The paper calls this",
        "quantity BASE (= R(0) = KS / KD) in Equations 3 and 4 but does not list",
        "it as an estimated parameter in Table 2, because it is supplied per arm",
        "from the digitised source trial. The model is exactly scale-invariant in",
        "BMD_BL for percent-change-from-baseline output -- BMD_LS(t) / BMD_BL",
        "does not depend on BMD_BL -- so the choice of baseline affects only the",
        "absolute g/cm^2 scale, which is why the paper's own figures are all",
        "plotted as percent change from baseline."
      ),
      source_name        = "Baseline of vertebral BMD (Wang 2023 Table 1); BASE (Wang 2023 Equations 3 and 4)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age in years. Wang 2023 Table 1 reports per-study mean ages from 57.2 to 85.4 years.",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Wang 2023 Methods 2.4 lists age among the candidate covariates screened",
        "by NONMEM stepwise SCM (forward dOFV > 6.63, backward dOFV > 10.83).",
        "Results 3.2: 'At present, no covariates associated with it were filtered",
        "out.' The Discussion attributes the null screen to the limited amount of",
        "aggregate data. Not retained in the final model."
      )
    ),
    SEXF = list(
      description = "Sex indicator (1 = female). Wang 2023 Table 1 reports 100% female in 9 of the 10 included trials and 93.6% female in the tenth (Nakamura 2017), i.e. 99.7% female overall.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened per Wang 2023 Methods 2.4; not retained (Results 3.2). The Discussion notes 'the number of male subjects included in the model was too few' to support a sex effect."
    ),
    BMI = list(
      description = "Body mass index. Wang 2023 Table 1 reports per-study mean BMI from 22.7 to 28.2 kg/m^2 for the 6 trials that reported it.",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened per Wang 2023 Methods 2.4; not retained (Results 3.2). BMI and body weight were both candidates and Methods 2.4 states that for covariate pairs with correlation > 0.8 only one was carried forward."
    ),
    WT = list(
      description = "Body weight in kg. Wang 2023 Table 1 reports per-study mean weights from 52.4 to 68.0 kg for the 6 trials that reported it.",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened per Wang 2023 Methods 2.4; not retained (Results 3.2). The Introduction motivates weight as a candidate because zoledronic-acid dose requirement should track skeletal size, but the covariate search returned nothing significant."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 6014L,
    n_studies      = 10L,
    age_range      = "57.2-85.4 years (per-study means; Wang 2023 Table 1)",
    weight_range   = "52.4-68.0 kg (per-study means, 6 of 10 trials reported weight; Wang 2023 Table 1)",
    sex_female_pct = 99.7,
    race_ethnicity = c(
      `Chinese or Japanese (Liang 2017, Li 2022, Nakamura 2017)` = 16.0,
      `predominantly White (Grey, Black, Greenspan trials)`      = 84.0
    ),
    disease_state  = "Primary osteoporosis or osteopenia; baseline vertebral BMD 0.64-1.06 g/cm^2 (Wang 2023 Table 1)",
    dose_range     = "Zoledronic acid 1, 2.5 or 5 mg intravenously as a single dose, or 5 mg intravenously once yearly for up to 6 years (Wang 2023 Table 1)",
    regions        = "New Zealand, United States, multinational, China, Japan",
    followup_range = "12-72 months (Wang 2023 Results 3.1)",
    data_source    = paste(
      "Aggregate (arm-level) data digitised from published figures and tables",
      "with GetData Graph Digitizer 1.9; no individual patient data. 454 trials",
      "were screened and 10 met the inclusion criteria (randomised controlled",
      "trials of zoledronic acid in primary osteoporosis reporting both lumbar",
      "spine and total hip DXA BMD)."
    ),
    notes          = paste(
      "Wang 2023 Table 1 lists the 10 contributing trials: Grey 2014 (n = 180),",
      "Grey 2012a (n = 180), Black 2012 (n = 616), Black 2007 (n = 3875),",
      "Liang 2017 (n = 175), Greenspan 2015 (n = 89), Li 2022 (n = 458),",
      "Nakamura 2017 (n = 330), Grey 2012b (n = 20) and Grey 2022 (n = 91).",
      "The race/ethnicity split above is the patient-weighted share of the",
      "explicitly East-Asian cohorts (Liang 2017, Li 2022, Nakamura 2017;",
      "963 of 6014 patients) versus the remainder, which the source trials",
      "describe as New Zealand / United States / multinational postmenopausal",
      "cohorts without a tabulated race breakdown -- Wang 2023 itself reports",
      "no race column, so this is an inference from the contributing trials and",
      "is flagged in the vignette Assumptions and deviations. The separate",
      "exposure-response dataset for the acute phase reaction (Wang 2023",
      "Table 3) draws on 6 further trial arms (Grey 2012a, Shiraki 2017,",
      "Grey 2012b, Popp 2017, Sieber 2013, Li 2022) and is not part of the",
      "K-PD fit."
    )
  )

  ini({
    # -----------------------------------------------------------------
    # Structural parameters -- Wang 2023 Table 2, "Estimate" column
    # (final model, FOCE-I, NONMEM 7.2.0). The bootstrap median and 95%
    # CI from the same table are recorded per line for transparency; the
    # point estimates are what is encoded. No structural parameter was
    # held fixed, so none is wrapped in fixed().
    #
    # All four are strictly positive rate / potency constants and are
    # therefore log-transformed.
    # -----------------------------------------------------------------

    # Eq 1-2: virtual K-PD amount compartment. `kel` is the canonical
    # K-PD elimination rate (registered source alias `kde`), which is
    # the paper's KDE.
    lkel <- log(0.08150)
    label("K-PD virtual-compartment first-order elimination rate KDE (1/month)")
    # Wang 2023 Table 2: KDE = 0.08150 1/month (RSE 71.5%; bootstrap median 0.12844, 95% CI 0.01812-1.08467)

    # Eq 3-4: first-order BMD degradation. `kout` is the canonical
    # indirect-response elimination rate; the paper writes it KD.
    lkout <- log(0.00474)
    label("First-order lumbar-spine BMD degradation rate KD (1/month)")
    # Wang 2023 Table 2: KD = 0.00474 1/month (RSE 17.9%; bootstrap median 0.00516, 95% CI 0.00385-0.00697)

    # Eq 5: dose-driving rate at half-maximal stimulation of BMD
    # synthesis. `ekd50` is the K-PD drug-signal potency canonical
    # already used by the sibling zoledronic-acid extraction
    # Mori_2018_zoledronicAcid.R; Wang 2023 writes the same quantity
    # EDK50 (letters transposed).
    lekd50 <- log(41300.0)
    label("Virtual dose-driving rate giving 50% stimulation of BMD synthesis EDK50 (ng/month)")
    # Wang 2023 Table 2: EDK50 = 41300.0 ng/month (RSE 40.9%; bootstrap median 28925.7, 95% CI 335.7-108289.8)

    # Eq 6: first-order tolerance / desensitisation rate on the
    # stimulated synthesis rate. `kdes` is the registered
    # desensitisation-rate canonical; Wang 2023 writes it k (Table 2
    # row "K") and Discussion attributes it to declining osteoclast
    # receptor sensitivity under repeated dosing.
    lkdes <- log(0.00754)
    label("First-order desensitisation (tolerance) rate k on stimulated BMD synthesis (1/month)")
    # Wang 2023 Table 2: K = 0.00754 1/month (RSE 20.3%; bootstrap median 0.00674, 95% CI 0.00001-0.01569)

    # -----------------------------------------------------------------
    # Random effects -- Wang 2023 Table 2.
    #
    # The units of analysis are published trial ARMS (arm-level median
    # BMD-time curves digitised from 10 trials), so the reported "ETA"
    # is between-STUDY variability, not between-subject variability.
    # It is therefore named eta_study_* per the MBMA convention rather
    # than etalekd50, and the model must not be read as an
    # individual-level popPK/PD model.
    #
    # ETA of EDK50 is the only estimated random effect. ETA of KDE,
    # ETA of KD and ETA of K are all reported as "FIXED" in Table 2,
    # i.e. fixed to zero (no variability); they are OMITTED here rather
    # than written as `~ fixed(0)`, because a zero-variance diagonal
    # makes OMEGA singular and breaks the Cholesky sampler used by
    # rxSolve.
    #
    # The variance is applied on the LOG scale (exponential random
    # effect), not additively as Wang 2023 Equation 8 (`Kei = TVKe +
    # etai`) literally prints. Equation 8 is generic boilerplate -- it
    # is written about "Ke", "the observed concentration value" and
    # "the individual PK parameters", none of which exist in this
    # paper's K-PD BMD model, and the promised second (residual)
    # equation is missing entirely. Read additively, omega = sqrt(0.4541)
    # = 0.674 ng/month on an EDK50 of 41300 ng/month is a coefficient of
    # variation of 0.0016%, i.e. no variability at all, which cannot
    # produce the visibly wide simulated 95% bands in Wang 2023 Figures
    # 4A/4B. Read exponentially it is 75.8% CV, and it reproduces the
    # published Figure 4A upper 95% bound (7.1% BMD gain) to two
    # significant figures. See vignette Errata.
    eta_study_lekd50 ~ 0.45410
    # Wang 2023 Table 2: ETA of EDK50 = 0.45410 (RSE 66.4%; bootstrap median 0.40691, 95% CI 0.00001-2.54713; shrinkage 18.7%)

    # -----------------------------------------------------------------
    # Residual error -- Wang 2023 Table 2 row "EPS".
    #
    # Wang 2023 Methods 2.3 states that the residual epsilon "follows a
    # (0, sigma^2) normal distribution", i.e. additive. NONMEM $SIGMA is
    # reported as a VARIANCE, so the additive SD is sqrt(0.00005) =
    # 0.00707 g/cm^2. That is ~0.9% of the 0.64-1.06 g/cm^2 baseline
    # range in Table 1, consistent with DXA precision. (Reading 0.00005
    # as an SD instead would imply 0.006% precision, which no DXA
    # instrument achieves.) The expression is left as sqrt(0.00005) so
    # the published number stays visible in the source trace.
    # -----------------------------------------------------------------
    addSd <- sqrt(0.00005)
    label("Additive residual SD on lumbar-spine BMD (g/cm^2); sqrt of the reported EPS variance 5e-5")
    # Wang 2023 Table 2: EPS = 0.00005 (variance; RSE 25.8%; bootstrap median 0.00003, 95% CI 0.00003-0.00006; shrinkage 7.1%)
  })

  model({
    # -----------------------------------------------------------------
    # 1. Arm-level parameters. Only ekd50 carries a random effect.
    # -----------------------------------------------------------------
    kel   <- exp(lkel)
    kout  <- exp(lkout)
    ekd50 <- exp(lekd50 + eta_study_lekd50)
    kdes  <- exp(lkdes)

    # -----------------------------------------------------------------
    # 2. Zero-order BMD synthesis rate (Wang 2023 Equation 4):
    #      KS = BASE * KD
    # With BASE = BMD_BL this makes the drug-free, tolerance-free
    # steady state exactly BMD_BL (Equation 3: R(0) = KS / KD = BASE).
    # Units: (g/cm^2) * (1/month) = (g/cm^2)/month, matching
    # d/dt(BMD_LS).
    # -----------------------------------------------------------------
    kin <- BMD_BL * kout

    # -----------------------------------------------------------------
    # 3. Virtual K-PD compartment (Wang 2023 Equations 1 and 2):
    #      dA/dt = -KDE * A
    #      IR    =  KDE * A
    # The dose (ng) is administered into depot_kpd. `ir` is the virtual
    # dose-driving rate in ng/month, the same units as ekd50.
    # -----------------------------------------------------------------
    d/dt(depot_kpd) <- -kel * depot_kpd
    ir <- kel * depot_kpd

    # -----------------------------------------------------------------
    # 4. Unit-Emax stimulation of BMD synthesis (Wang 2023 Equation 5):
    #      KS' = KS * (1 + IR / (EDK50 + IR))
    # The paper prints no separate Emax parameter and Table 2 lists
    # none, so the maximum achievable stimulation is exactly 2-fold.
    # Dimensionless because ir and ekd50 share units.
    # -----------------------------------------------------------------
    kin_stim <- kin * (1 + ir / (ekd50 + ir))

    # -----------------------------------------------------------------
    # 5. Tolerance / desensitisation factor (Wang 2023 Equation 6):
    #      Tol = exp(-k * time)
    # `time` is absolute time since the start of the regimen (rxode2's
    # `t`), not time since the most recent dose: Wang 2023 introduces
    # the term specifically to reproduce the attenuation of BMD gain
    # across successive annual doses (Results 3.2, Discussion), which a
    # per-dose reset would not produce. Dimensionless.
    # -----------------------------------------------------------------
    tol <- exp(-kdes * t)

    # -----------------------------------------------------------------
    # 6. Lumbar-spine BMD turnover (Wang 2023 Equation 7):
    #      dR/dt = KS' * Tol - KD * R,   R(0) = BASE
    # Note that Tol multiplies the WHOLE stimulated synthesis rate,
    # including its baseline part, exactly as printed -- not only the
    # drug-driven increment. A consequence is that a long simulation
    # with no further dosing drifts below baseline (the drug-free
    # attractor is BMD_BL * exp(-kdes * t) rather than BMD_BL). This is
    # the literal published equation and is retained; the behaviour is
    # quantified in the vignette Errata.
    # Units: (g/cm^2)/month on both terms.
    # -----------------------------------------------------------------
    d/dt(BMD_LS) <- kin_stim * tol - kout * BMD_LS
    BMD_LS(0)    <- BMD_BL

    # -----------------------------------------------------------------
    # 7. Observation. BMD_LS is itself the ODE state, so event-table
    # observation rows use cmt = "BMD_LS" and dose rows use
    # cmt = "depot_kpd"; there is no algebraic observable to alias.
    # -----------------------------------------------------------------
    BMD_LS ~ add(addSd)
  })
}
