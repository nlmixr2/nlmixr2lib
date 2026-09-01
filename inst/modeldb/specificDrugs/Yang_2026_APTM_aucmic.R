Yang_2026_APTM_aucmic <- function() {
  description <- "Preclinical (chicken, specific-pathogen-free). Inhibitory sigmoid Emax PK/PD-index model for the in vivo anti-mycoplasma effect of APTM (14-O-[(4-amino-6-hydroxy-pyrimidine-2-yl) thioacetyl] mutilin), a novel semi-synthetic pleuromutilin derivative, against Mycoplasma gallisepticum strain S6 in an intratracheal chicken infection model, driven by the AUC/MIC index. Yang 2026 Materials and methods parameterises the effect as E = E0 - Imax * Ce^gamma / (IC50^gamma + Ce^gamma), where E is the SIGNED change in lung mycoplasma load in log10 CFU/mL accrued over the 72 h treatment course (negative = bacterial reduction), E0 is the corresponding change in the untreated control, Imax is the maximum attainable reduction, Ce is the PK/PD index and gamma is the Hill coefficient. Parameters from Yang 2026 Table 3, AUC0-24h/MIC column: Imax = 3.986 log10 CFU/mL, IC50 = 490.449, E0 = 0.179 log10 CFU/mL, gamma = 1.525. TWO READINGS OF THE PRINTED EQUATION WERE ADJUDICATED NUMERICALLY. The printed denominator reads 'IC50 + Ce^gamma', omitting the exponent on IC50 that the companion time-kill equation on the same page does carry; the definition in the surrounding text ('IC50 is the index value that produces 50% of the maximum inhibitory effect') requires IC50^gamma. Substituting Table 3 into the IC50^gamma form returns an inhibitory term of exactly 1.0001 at the paper's stated 1-log10 target of 239.39 and exactly 2.0001 at its stated 2-log10 target of 492.75, whereas the literal printed form returns 3.84 at the 2-log10 target. The IC50^gamma reading is therefore used, and the same substitution establishes that the 'Log10CFU/mL drop' rows of Table 3 are the value of the INHIBITORY TERM rather than of the net change E (at the 1-log target the model predicts E = -0.821, not -1). There is NO PK component: exposure enters as the externally supplied AUC covariate AUC_APTM divided by the parameter mic, because Yang 2026 analysed the plasma concentrations non-compartmentally in Phoenix WinNonlin (Table 2) and published no structural PK model. The model predicts the CHANGE in lung load directly rather than integrating a bacterial density, because Yang 2026's in vivo pharmacodynamic readout is a single cross-sectional count per dose group at 72 h expressed relative to a baseline group whose absolute lung load, although measured, is never reported; predicting the change is exactly the quantity the paper fitted and requires no unreported baseline. Yang 2026 reports neither between-subject variability nor a residual error magnitude, so no eta parameters are present and addSd is FIXED at 0 for deterministic typical-value simulation. Sibling models from the same paper: Yang_2026_APTM_cmaxmic (the same structure re-fitted against the Cmax/MIC index) and Yang_2026_APTM_timekill (the in vitro time-kill model)."
  reference <- paste(
    "Yang W, Ding H, Ma X, Lv T, Wang L. (2026).",
    "Pharmacokinetic/pharmacodynamic relationship of a novel pleuromutilin",
    "derivative APTM against Mycoplasma gallisepticum.",
    "Poultry Science 105:106560.",
    "doi:10.1016/j.psj.2026.106560. PMCID: PMC12919259.",
    "Model equation: Materials and methods, 'Pharmacokinetic, pharmacodynamic,",
    "and statistical analysis'. Parameter estimates and PK/PD-index targets:",
    "Table 3, AUC 0-24h/MIC column. Non-compartmental plasma exposures:",
    "Table 2. MIC: Results, 'In vitro susceptibility of MG to APTM'.",
    "In vivo dose-response: Results, 'In vivo efficacy of APTM against",
    "M. gallisepticum' and Figure 4A-B.",
    sep = " "
  )
  vignette <- "Yang_2026_APTM"
  units <- list(
    time = "h",
    dosing = "ug*h/mL (APTM plasma AUC over the treatment course, supplied as a covariate)",
    concentration = "log10 CFU/mL (observation; the signed change in lung mycoplasma load over 72 h)"
  )

  depends <- c("AUC_APTM")

  covariateData <- list(
    AUC_APTM = list(
      description        = "APTM area under the plasma concentration-time curve, used as the numerator of the AUC/MIC PK/PD index driving the in vivo effect",
      units              = "ug*h/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Yang 2026 calls this exposure 'AUC0-24h' throughout, but the values it actually regressed against dose are the AUCINF column of Table 2. Two independent checks establish this: least-squares regression of the Table 2 AUCINF values (2382.53, 8551.66 and 30671.54 h*ng/mL at 5, 15 and 40 mg/kg) on dose returns R^2 = 0.99480, matching the paper's stated 0.9948 to five significant figures; and inverting that same regression at the Table 3 index targets reproduces the paper's own back-calculated daily doses of 12.22 and 21.83 mg/kg to within 0.2% (12.24 and 21.86). Set to 0 for the untreated control so the sigmoid term vanishes and the predicted change equals E0. Because the paper published no structural PK model, this covariate is the model's only route for drug exposure. Note the unit change relative to Table 2: the paper tabulates h*ng/mL, so divide the tabulated value by 1000 to obtain the ug*h/mL this column expects, which keeps the AUC/MIC ratio consistent with mic in ug/mL.",
      source_name        = "AUC0-24h (Yang 2026 Materials and methods, 'Pharmacokinetic, pharmacodynamic, and statistical analysis'; Figure 3B; index targets in Table 3). Numerically the AUCINF column of Table 2."
    )
  )

  population <- list(
    species        = "chicken (specific-pathogen-free, one-day-old at purchase, 35-45 g)",
    n_subjects     = 60L,
    n_studies      = 1L,
    organism       = "Mycoplasma gallisepticum standard strain S6 (ATCC 15302; China Institute of Veterinary Drug Control). APTM MIC = 0.03125 ug/mL by broth microdilution and 0.0625 ug/mL by broth macrodilution; the microdilution value is the one used to form the PK/PD indices and, as a surrogate for MIC90, in the dose calculation",
    system         = "Intratracheal M. gallisepticum infection model: chickens inoculated intratracheally with 0.2 mL of a 1 x 10^9 CFU/mL exponential-phase suspension once daily for three consecutive days",
    disease_state  = "Experimental M. gallisepticum respiratory infection; efficacy read as the change in lung mycoplasma load relative to a baseline group euthanised before the first dose",
    dose_range     = "0 (vehicle control), 5, 10, 15, 20, 25, 30, 35 and 40 mg/kg APTM by oral gavage, once daily for three consecutive days (45% APTM soluble powder dissolved in water)",
    design         = "60 infected chickens: eight treatment groups plus a vehicle control group (n = 6 each) and an additional baseline group (n = 6) euthanised before the first administration to establish the initial lung load",
    sampling       = "Lungs aseptically collected and homogenised 24 h after the final dose (72 h after the first), homogenate adjusted to 1 mL with sterile saline, ten-fold serially diluted and plated for viable mycoplasma counts",
    regions        = "China (South China Agricultural University, Guangzhou)",
    notes          = "Ethics approval 2025C037 (Animal Ethics Committee of South China Agricultural University). Observed mean changes in lung load: +0.20 log10 CFU/mL in the untreated control and reductions of 0.15, 0.31, 1.13, 1.85, 1.63, 2.67, 2.48 and 2.80 log10 CFU/mL at 5, 10, 15, 20, 25, 30, 35 and 40 mg/kg. The companion Cmax/MIC fit is packaged separately as Yang_2026_APTM_cmaxmic; the two indices are statistically indistinguishable in this study (R^2 = 0.9424 for AUC0-24h/MIC versus 0.9428 for Cmax/MIC) and the paper prefers AUC0-24h/MIC on mechanistic grounds (concentration-dependent killing plus a long half-life), which is why both parameterisations are carried rather than only the better-correlating one. The paper's printed dose equation, Dose = (AUC/MIC breakpoint x MIC90 x Cl) / (fu x F) with fu set to 1 and F not factored in, yields daily oral doses of 12.22 mg/kg (1-log10 target) and 21.83 mg/kg (2-log10 target), the latter rounded to the paper's headline 22 mg/kg recommendation. The 264 chickens of the pharmacokinetic arm are not counted in n_subjects, which refers to the pharmacodynamic cohort this model was fitted to; the whole study used 324 birds."
  )

  ini({
    # =================================================================
    # Yang 2026 Table 3, AUC 0-24h/MIC column
    # =================================================================
    # Yang 2026 Materials and methods, "Pharmacokinetic, pharmacodynamic,
    # and statistical analysis", printed equation:
    #
    #   E = E0 - Imax * Ce^gamma / (IC50 + Ce^gamma)
    #
    # The printed denominator omits the exponent on IC50 that the
    # time-kill equation three paragraphs earlier does carry, and the
    # surrounding definition ("IC50 is the index value that produces 50%
    # of the maximum inhibitory effect") requires IC50^gamma. The
    # IC50^gamma reading is used, and is confirmed by the paper's own
    # Table 3 targets: the inhibitory term evaluates to 1.0001 at the
    # stated 1-log10 target of 239.39 and to 2.0001 at the stated
    # 2-log10 target of 492.75, whereas the literal printed form gives
    # 3.84 at the 2-log10 target. That substitution also establishes
    # that the "Log10CFU/mL drop" rows of Table 3 are the value of the
    # inhibitory TERM, not of the net change E. See the vignette.
    #
    # E is the SIGNED change in lung log10 CFU/mL over the 72 h course,
    # so E0 is a small POSITIVE number (the untreated control grew) and
    # stays on the natural scale because it is not sign-constrained.
    e0 <- 0.179
    label("Change in lung mycoplasma load in the untreated control over 72 h E0 (log10 CFU/mL)")  # Yang 2026 Table 3, AUC 0-24h/MIC column, E 0 = 0.179
    limax <- log(3.986)
    label("Log maximum attainable reduction in lung mycoplasma load Imax (log10 CFU/mL)")  # Yang 2026 Table 3, AUC 0-24h/MIC column, I max = 3.986
    lec50 <- log(490.449)
    label("Log AUC/MIC index value producing 50% of the maximum inhibitory effect IC50 (h)")  # Yang 2026 Table 3, AUC 0-24h/MIC column, IC 50 = 490.449
    lhill <- log(1.525)
    label("Log Hill coefficient gamma defining the steepness of the index-effect curve (unitless)")  # Yang 2026 Table 3, AUC 0-24h/MIC column, gamma = 1.525

    # =================================================================
    # Minimum inhibitory concentration
    # =================================================================
    # The PK/PD index driving the sigmoid is AUC/MIC, so the MIC is
    # required numerically. FIXED because it is a measured property of
    # the challenge strain, not an estimated parameter. Change it to
    # apply the model to an isolate with a different susceptibility.
    # The microdilution value is the one the index used: reproducing the
    # paper's back-calculated doses from the Table 2 exposure regression
    # requires 0.03125, not the macrodilution 0.0625.
    mic <- fixed(0.03125)
    label("APTM MIC against M. gallisepticum S6 (ug/mL; broth microdilution)")  # Yang 2026 Results, "In vitro susceptibility of MG to APTM": 0.03125 ug/mL by microdilution (0.0625 ug/mL by macrodilution), constant across 10^5, 10^6 and 10^7 CFU/mL inocula

    # =================================================================
    # Residual error
    # =================================================================
    # Yang 2026 Table 3 reports the coefficient of determination of the
    # sigmoid fit (R^2 = 0.942) and no residual standard deviation, so
    # the residual SD is held at zero for deterministic typical-value
    # simulation. See the vignette Assumptions and deviations section.
    addSd <- fixed(0)
    label("Additive residual SD on the 72 h change in log10 CFU/mL (0; not reported in Yang 2026)")  # Yang 2026 Table 3 reports R^2 only, no residual SD
  })

  model({
    imax <- exp(limax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # PK/PD index: APTM plasma AUC divided by the MIC of the challenge
    # strain (Yang 2026 Materials and methods). Units: h.
    aucmic <- AUC_APTM / mic

    # Yang 2026 inhibitory sigmoid Emax equation, in the IC50^gamma
    # reading adjudicated above. Cc is the SIGNED change in lung
    # mycoplasma load over the 72 h treatment course in log10 CFU/mL
    # (negative = reduction); it equals e0 at zero exposure and
    # approaches e0 - imax at saturating exposure.
    Cc <- e0 - imax * aucmic^hill / (ec50^hill + aucmic^hill)
    Cc ~ add(addSd)
  })
}
