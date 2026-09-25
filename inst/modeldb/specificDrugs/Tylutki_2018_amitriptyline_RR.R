Tylutki_2018_amitriptyline_RR <- function() {
  description <- paste(
    "Sigmoid Emax model relating total plasma amitriptyline concentration",
    "to the absolute electrocardiographic R-R interval length (observation",
    "variable rr, ms) in humans, from Tylutki 2018's PBPK-QSTS cardiac-",
    "safety analysis. The model is",
    "rr = e0 + (rmax - e0) * C^hill / (ec50^hill + C^hill),",
    "with e0 = 995.3 ms (drug-free baseline R-R, i.e. 60 beats/min),",
    "rmax = 500.8 ms (asymptotic R-R at maximal effect, i.e. 120",
    "beats/min), ec50 = 0.4 umol/L and hill = 1.5. Amitriptyline SHORTENS",
    "the R-R interval (speeds the heart) so rmax lies BELOW e0 and the",
    "curve is monotonically decreasing; the maximal attainable effect is",
    "e0 - rmax = 494.5 ms of R-R shortening. The fit is dominated by",
    "acute-overdose case reports, which supply almost all of the",
    "concentration range above 1 umol/L.",
    "PD-ONLY model: amitriptyline total plasma concentration is supplied",
    "as the time-varying covariate CP_AMITRIPTYLINE_UM (umol/L). The",
    "source paper's PK layer is a full-PBPK model for amitriptyline linked",
    "to a minimal-PBPK model for nortriptyline that Tylutki 2018 imports",
    "UNCHANGED from Tylutki 2018 J Pharm Sci (doi:10.1016/j.xphs.2017.11.012)",
    "and whose organ volumes, blood flows, partition coefficients and",
    "clearances appear nowhere in this paper; that layer is therefore not",
    "reproduced here and users must supply their own concentration",
    "trajectory. Amitriptyline molecular weight is 277.4 g/mol, so",
    "1 umol/L = 277.4 ng/mL.",
    "NOTE: the R-R equation is printed incorrectly in the source (Eq. 4",
    "omits its leading baseline term and evaluates to 0 at zero",
    "concentration). The form encoded here is the corrected one; it is the",
    "unique reconstruction from the printed parameter symbols that",
    "reproduces the published Figure 2 curve. See the ini() block and the",
    "vignette Errata for the digitisation evidence.",
    sep = " "
  )

  reference <- paste(
    "Tylutki Z, Mendyk A, Polak S. (2018). Physiologically based",
    "pharmacokinetic-quantitative systems toxicology and safety",
    "(PBPK-QSTS) modeling approach applied to predict the variability of",
    "amitriptyline pharmacokinetics and cardiac safety in populations and",
    "in individuals. Journal of Pharmacokinetics and Pharmacodynamics",
    "45(5):663-677. doi:10.1007/s10928-018-9597-6. PMCID PMC6182726.",
    "Parameter values are from the Results section 'Emax model'; the",
    "structural equation is Eq. 4 of the Methods section 'Pharmacodynamic",
    "models', corrected against the published Figure 2 (see ini()).",
    sep = " "
  )

  vignette <- "Tylutki_2018_amitriptyline"

  units <- list(
    time = "h",
    dosing = "(none; PD-only model driven by an external amitriptyline plasma-concentration covariate)",
    concentration = "(observation rr is the absolute R-R interval length in ms; driving covariate CP_AMITRIPTYLINE_UM is in umol/L)"
  )

  covariateData <- list(
    CP_AMITRIPTYLINE_UM = list(
      description = "Instantaneous TOTAL (not unbound) amitriptyline plasma concentration at the time of each R-R observation, supplied as a time-varying covariate from observed plasma samples or an upstream PK source.",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying per event row. Drives the sigmoid Emax expression",
        "rr = e0 + (rmax - e0) * CP_AMITRIPTYLINE_UM^hill /",
        "(ec50^hill + CP_AMITRIPTYLINE_UM^hill).",
        "Tylutki 2018 Eq. 4 defines C as 'AT total plasma concentration",
        "[uM]', so the column is TOTAL drug in plasma and is NOT the free",
        "cardiac-tissue concentration that the same paper feeds to the",
        "Cardiac Safety Simulator for its QT endpoint. Those two",
        "quantities differ by roughly four orders of magnitude: Tylutki",
        "2018 Eq. 5 gives Cfree,cardiac = Ctotal,plasma * Kp_ht * fu_ht",
        "with Kp_ht = 11.77 and fu_ht = 0.0012 for amitriptyline (and",
        "Kp_ht = 35.63, fu_ht = 0.001 for nortriptyline). Supplying a free",
        "cardiac concentration to this model would understate the driver",
        "by a factor of about 71.",
        "Amitriptyline molecular weight is 277.4 g/mol (Tylutki 2018",
        "Methods, citing PubChem CID 2160), so 1 umol/L = 277.4 ng/mL and",
        "a concentration reported in ng/mL must be divided by 277.4 before",
        "being supplied on this column.",
        "Reference values observed: the Figure 2 dataset spans roughly 0.1",
        "to 17.8 umol/L (28 to 4900 ng/mL). Therapeutic amitriptyline",
        "trough concentrations of 50-300 ng/mL correspond to 0.18-1.08",
        "umol/L, i.e. they straddle the fitted ec50 of 0.4 umol/L; the",
        "concentrations above about 2 umol/L come almost entirely from the",
        "acute-overdose case reports in the fitting set.",
        "Set to 0 for the drug-free reference, at which the expression",
        "collapses to rr = e0 = 995.3 ms.",
        "The natural input source is the paper's own full-PBPK",
        "amitriptyline model, which is imported unchanged from Tylutki",
        "2018 J Pharm Sci 107:1167-1177 and is not reproducible from this",
        "publication (no organ volume, blood flow, tissue partition",
        "coefficient or clearance is printed here); no amitriptyline popPK",
        "model is in the nlmixr2lib registry, so users must supply their",
        "own trajectory."
      ),
      source_name = "C (AT total plasma concentration, Tylutki 2018 Eq. 4)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 18L,
    age_range = "not reported for the pooled R-R dataset; the constituent sources range from a paediatric intoxication series to elderly overdose case reports",
    weight_range = "not reported",
    sex_female_pct = NA_real_,
    race_ethnicity = "not reported",
    disease_state = paste(
      "Pooled literature dataset of paired (total plasma amitriptyline",
      "concentration, R-R interval length) observations assembled by",
      "Tylutki 2018 from 18 published sources (their references 29-46).",
      "The sources are heterogeneous by design and span three settings:",
      "(1) healthy-volunteer cardiovascular studies of therapeutic",
      "amitriptyline dosing (Warrington 1989 Br J Clin Pharmacol 27:343;",
      "Wester 1980 Curr Med Res Opin 6:29; Stern 1985 Pharmacopsychiatry",
      "78:272), (2) a heart-rate analysis in 24 depressed patients taking",
      "150 mg/day (Rechlin 1994 Psychopharmacology 116:110), and (3)",
      "acute-amitriptyline-overdose case reports and case series",
      "(Diaz-Buxo 1978, Kansal 2017, Yang 1991, Amitai 1993, Rudorfer",
      "1982, Schmidt 2015, Ozayar 2012, Karaci 2013, Baysal 2007, Sein",
      "Anand 2005, Zakynthinos 2000, Gomolin 1983, Huge 2011, Spiker",
      "1975). The overdose sources supply essentially the whole",
      "concentration range above about 2 umol/L, so the plateau parameter",
      "rmax is informed by intoxicated rather than therapeutically dosed",
      "subjects."
    ),
    dose_range = "not applicable -- the fitting dataset is indexed by measured plasma concentration, not by dose",
    regions = "not reported (the constituent publications are European, North American, Middle Eastern and East Asian)",
    notes = paste(
      "The number of paired observations is not stated numerically in the",
      "paper; Figure 2 displays approximately 60 points, annotated by sex",
      "as female, male or not available. Sex is shown as a plotting",
      "aesthetic only -- it is not a covariate in the fitted model.",
      "Fitting was by simulated annealing ('SANN') in the FME package",
      "under R 3.4.0 (Tylutki 2018 Methods 'Pharmacodynamic models').",
      "The reported model RMSE is 120.98 ms; see the ini() block for why",
      "that is encoded as the additive residual standard deviation.",
      "This file encodes ONLY the R-R sub-model of Tylutki 2018. The",
      "paper's other two components are not extractable from it: the",
      "amitriptyline full-PBPK / nortriptyline minimal-PBPK layer is",
      "imported unchanged from a prior publication whose parameters are",
      "not restated here, and the QT / pseudo-ECG layer is the ten",
      "Tusscher & Panfilov 2006 ventricular cardiomyocyte model as",
      "implemented in the proprietary Cardiac Safety Simulator v2.1",
      "(Simcyp / Certara) platform. The oral-absorption parameters this",
      "paper does contribute (ka = 0.24 1/h and a mean lag time of 1.33 h",
      "with 30% CV, plus F drawn from N(0.459, 0.093) truncated to",
      "[0.33, 0.62] and fa*Fg drawn from a log-normal with mean 0.832 and",
      "CV 0.131) are an input function to that imported PBPK model and",
      "are not runnable on their own; they are recorded in the vignette",
      "rather than encoded here."
    )
  )

  ini({
    # ==================================================================
    # Sigmoid Emax model for the R-R interval (Tylutki 2018 Eq. 4 and
    # Results section 'Emax model'). All four estimates are printed in a
    # single Results sentence:
    #
    #   'The estimates of Emax model parameters were as follows:
    #    RR0 = 995.3 [ms], RRmax = 500.8 [ms], EC50 = 0.4 [uM], and
    #    n = 1.5. The RMSE of established Emax model equaled 120.98.'
    #
    # THE PUBLISHED EQUATION IS WRONG AND HAS BEEN CORRECTED HERE.
    # Eq. 4 is printed as
    #
    #   RR = (RR0 - RRmax) * C^n / (EC50^n + C^n)
    #
    # which is zero at C = 0 and INCREASES to 494.5 ms - the opposite of
    # the paper's own Figure 2, which starts at ~995 ms at C = 0 and
    # DECREASES to a plateau just above 500 ms. The printed form also
    # contradicts the paper's stated mechanism (amitriptyline raises
    # heart rate, i.e. shortens R-R) and its own definition of RR0 as
    # 'the baseline R-R interval length'. The equation is missing its
    # leading baseline term.
    #
    # Two reconstructions restore RR0 at C = 0 using only the printed
    # symbols:
    #   (A) RR = RR0 - RRmax        * C^n/(EC50^n + C^n)  -> plateau 494.5
    #   (B) RR = RR0 - (RR0-RRmax)  * C^n/(EC50^n + C^n)  -> plateau 500.8
    # (B) is the printed numerator with the dropped 'RR0 -' restored, and
    # it is the reading under which RRmax is what the paper says it is,
    # an R-R interval LENGTH in ms rather than a change in one.
    #
    # Figure 2 discriminates between them. The figure was rendered at
    # 600 dpi, the axes calibrated on the panel gridlines (residuals of
    # the linear y-axis fit < 0.4 ms over 300-1000 ms), and the fitted
    # curve traced at 399 concentrations from 0.11 to 17.8 uM:
    #
    #   measured plateau (median over 68 clean columns) = 503.1 ms
    #     (B) predicts 500.8 ms  -> +2.3 ms, within the 3 px line width
    #     (A) predicts 494.5 ms  -> +8.6 ms
    #   mean residual over all 399 traced points
    #     (B) -0.63 ms     (A) +5.47 ms
    #
    # (B) is therefore adopted. It is written below in the equivalent
    # baseline-plus-plateau form used by Dings_2026_cafedrine_*.R,
    #   rr = e0 + (rmax - e0) * C^hill / (ec50^hill + C^hill),
    # so that both printed values (995.3 and 500.8) appear verbatim as
    # ini() parameters and the 494.5 ms amplitude is derived rather than
    # transcribed.
    # ==================================================================

    le0 <- log(995.3)
    label("Baseline drug-free R-R interval length (ms)")
    # Tylutki 2018 Results 'Emax model': RR0 = 995.3 ms. Equivalent to a
    # drug-free heart rate of 60000/995.3 = 60.3 beats/min. Anchors the
    # curve at C = 0.

    lrmax <- log(500.8)
    label("Asymptotic R-R interval length at maximal amitriptyline effect (ms)")
    # Tylutki 2018 Results 'Emax model': RRmax = 500.8 ms. This is the
    # PLATEAU of the curve, not an amplitude, and it lies BELOW the
    # baseline because amitriptyline shortens R-R; 500.8 ms corresponds
    # to 119.8 beats/min. The maximal attainable effect is therefore
    # e0 - rmax = 494.5 ms of shortening. Named lrmax after the
    # 'maximum attainable response' canonical used by
    # Dings_2026_cafedrine_theodrenaline_ephedrine.R (lrmax_hr /
    # lrmax_map / lrmax_sbp), which has the identical
    # baseline-plus-plateau structure, rather than lemax (an increment).
    # The paper's own gloss 'RRmax is the maximum R-R interval length
    # [ms]' is the plateau reading; its other gloss, 'EC50 is the AT
    # concentration that produces 50% of RRmax', would make RRmax an
    # effect magnitude instead. The two glosses are mutually
    # inconsistent and Figure 2 selects the plateau reading (see above).

    lec50 <- log(0.4)
    label("Amitriptyline concentration producing half the maximal R-R effect (umol/L)")
    # Tylutki 2018 Results 'Emax model': EC50 = 0.4 uM. Printed to a
    # single significant figure, which is the limiting precision of this
    # model. 0.4 umol/L = 111 ng/mL, i.e. inside the usual therapeutic
    # plasma range, so the curve's steepest region coincides with
    # therapeutic exposure.

    lhill <- log(1.5)
    label("Sigmoidicity (Hill) exponent n")
    # Tylutki 2018 Results 'Emax model': n = 1.5. The paper calls this
    # 'the sigmoidicity factor' (Methods, Eq. 4 glossary); canonical
    # name lhill per the register's Hill-coefficient convention.

    addSd <- 120.98
    label("Additive residual error standard deviation on the R-R interval (ms)")
    # Tylutki 2018 Results 'Emax model': 'The RMSE of established Emax
    # model equaled 120.98.' The fit was a least-squares fit by the
    # 'SANN' method of the FME package (Methods 'Pharmacodynamic
    # models'), for which the root-mean-square error IS the estimated
    # residual standard deviation on the modelled scale, so it is
    # encoded directly as an additive residual SD in ms rather than as
    # fixed(0). It is large relative to the 494.5 ms signal, which is
    # consistent with the scatter visible in Figure 2 and with the
    # pooled, heterogeneous provenance of the fitting data.

    # ==================================================================
    # Inter-individual variability: none. The Emax model was fitted to a
    # pooled literature dataset in which each source contributes one or
    # a few observations, and Tylutki 2018 reports no variance component
    # for it (the inter-individual variability the paper does discuss
    # belongs to the imported PBPK layer, not to this PD model). Per the
    # standing operator policy on unreported IIV this file ships a
    # typical-value-only model; the vignette Errata documents the gap.
    # ==================================================================
  })

  model({
    e0 <- exp(le0)
    rmax <- exp(lrmax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # ==================================================================
    # Corrected Tylutki 2018 Eq. 4 (see the ini() block for the
    # reconstruction and its validation against Figure 2). Written in
    # baseline-plus-plateau form; because rmax < e0 the second term is
    # negative and rr decreases monotonically from e0 towards rmax.
    #
    # CP_AMITRIPTYLINE_UM is the TOTAL plasma concentration in umol/L,
    # supplied as a time-varying covariate; ec50 is already on that
    # scale so no unit rescaling is applied.
    #
    # There are no ODE states: this is a purely algebraic PD model
    # driven by an external concentration column, so rxSolve() must NOT
    # be given omega = NA (there are no etas to suppress).
    # ==================================================================
    rr <- e0 + (rmax - e0) * CP_AMITRIPTYLINE_UM^hill /
      (ec50^hill + CP_AMITRIPTYLINE_UM^hill)

    rr ~ add(addSd)
  })
}
