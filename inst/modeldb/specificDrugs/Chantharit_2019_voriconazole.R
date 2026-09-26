Chantharit_2019_voriconazole <- function() {
  description <- "One-compartment population pharmacokinetic model with first-order absorption and linear elimination for oral voriconazole in Thai adults treated for invasive aspergillosis (Chantharit 2019); apparent clearance carries serum albumin and log10-transformed gamma-glutamyltransferase, apparent volume of distribution scales linearly with actual body weight, and the first-order absorption rate constant is fixed at 1 per hour"
  reference <- "Chantharit P, Tantasawat M, Kasai H, Tanigawara Y. 1566. Population pharmacokinetics of voriconazole: serum albumin status as a novel marker of clearance and dosage optimization. Open Forum Infect Dis. 2019;6(Suppl 2):S572. doi:10.1093/ofid/ofz360.1430"
  vignette <- "Chantharit_2019_voriconazole"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. The cohort received ORAL voriconazole exclusively
  # ("One hundred and six patients using oral VRCZ were included", Results),
  # so the model carries a depot and the disposition parameters are apparent
  # oral values (CL/F, V/F) - see the note above `lcl` in `ini()`.
  compartmentData <- list(
    depot = list(analyte = "voriconazole", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "voriconazole", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters clearance as (ALB/28)^0.93, from the final-model equation printed in the Results: 'CL (L/hr) = theta CL x (albumin/28) theta2 x (logGGT/2.4) theta3 x exp (eta CL)'. The 28 g/L divisor is the cohort centring value; the abstract prints no demographics table, so it cannot be cross-checked against a tabulated median (the 2019 abstract carries only Table 1, parameters, and Table 2, dosing simulations). The POSITIVE exponent means clearance FALLS as albumin falls, which is the direction the abstract states in prose and on which the whole paper turns: 'Patient with serum albumin 30 g/L had CL lower than patient having serum albumin > 30 g/L, P = 0.0007' and 'Serum albumin is a novel marker influencing VRCZ CL'. The centring value of 28 g/L sits just below the 30 g/L threshold the abstract uses to stratify its dosing simulations in Table 2, so the typical patient in this cohort is hypoalbuminaemic. Voriconazole is roughly 58% protein bound, and the abstract's own reading (echoed by Jiang 2022, which cites this work) is that hypoalbuminaemia raises the unbound fraction; note that the direction fitted here is on TOTAL-concentration apparent clearance, which is what the assay measured.",
      source_name = "albumin"
    ),
    GGT = list(
      description = "Serum gamma-glutamyltransferase activity",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters clearance as (log10(GGT)/2.4)^-0.58, from the final-model equation printed in the Results. The NEGATIVE exponent means clearance falls as the cholestatic marker rises, the expected direction for a hepatically metabolised triazole, and the same qualitative role GGT plays in Wang_2025_voriconazole.R (which cites this paper as precedent for retaining hepatic markers on voriconazole clearance). TWO readings of the printed term had to be settled and only one is defensible; see the extended comment above `e_ggt_cl` in `ini()` for the argument. In brief: the abstract writes 'logGGT' without naming the base, and the centring constant 2.4 is the cohort median of that log. Base-10 gives a median GGT of 10^2.4 = 251 U/L, roughly 4-5x the upper limit of normal and entirely ordinary for patients with invasive aspergillosis on prolonged voriconazole, a recognised hepatotoxin. Natural log would give e^2.4 = 11 U/L, which is BELOW the normal reference range and not credible as the median of a cohort in which the authors retained GGT as a statistically significant covariate. Base-10 is therefore used. Because the covariate enters through log10(GGT), the model is defined only for GGT > 1 U/L, and the term is a ratio of logs rather than a log of ratios, so it is NOT the usual power-of-median form.",
      source_name = "GGT"
    ),
    WT = list(
      description = "Actual body weight",
      units = "kg",
      type = "continuous",
      notes = "Enters the apparent volume of distribution as (WT/55)^1, from the final-model equation printed in the Results: 'V (Liter) = theta V x (Actual body weight/55) theta1 x exp (eta V)'. Table 1 reports theta1 as '1 (fixed)' in both the final-model and the bootstrap columns, so the weight exponent was HELD at unity rather than estimated - linear, not allometric, scaling. The 55 kg divisor is the cohort centring value and is consistent with an adult Thai population; it is noticeably lower than the 70 kg default common in Western cohorts, which matters here because the exponent is 1 and the term therefore rescales volume in direct proportion. The abstract prints no demographics table, so 55 kg cannot be cross-checked against a tabulated median.",
      source_name = "Actual body weight"
    )
  )

  # Screened and explicitly NOT retained. CYP2C19 is the headline negative
  # result of this abstract, so its provenance is worth preserving even
  # though it appears nowhere in `model()`.
  covariatesDataExcluded <- list(
    CYP2C19_PHENOTYPE = list(
      description = "CYP2C19 metabolizer phenotype",
      units = "categorical",
      type = "categorical",
      notes = "Phenotyped in 88 of the 106 patients - 43 extensive metabolizers, 37 intermediate metabolizers and 8 poor metabolizers (Results). Screened and rejected: 'CYP2C19 phenotypes did not influence any PK parameter during the covariate model building, then all 106 patients were included in the model construction.' The abstract returns to this when defending serum albumin as the operative marker, noting that the albumin effect held 'irrespective of PM status because of there were all phenotypes which distributed across two groups; %EM: %IM: %PM for 55: 38.3: 6.6 and 48:36:16'. Voriconazole is a CYP2C19 substrate and many published models DO retain the phenotype, so the negative finding is a genuine result of this cohort rather than an omission."
    )
  )

  population <- list(
    n_subjects = 106,
    n_studies = 1,
    disease_state = "invasive aspergillosis",
    species = "human",
    regions = "Thailand",
    route = "oral",
    genotype = "CYP2C19 phenotype available in 88 of 106 (43 extensive, 37 intermediate, 8 poor metabolizers); not retained in the final model",
    notes = "Combined dataset of intensive pharmacokinetic sampling and routine therapeutic-drug-monitoring trough concentrations. Estimation by FOCE ELS in Phoenix NLME; the final model was checked by bootstrap (1000 runs, 96% successful), visual predictive check and goodness-of-fit plots. The source is a one-page IDWeek 2019 poster abstract and prints NO baseline-demographics table, so age, sex and the covariate distributions are not reported; the covariate centring constants in the final-model equation (albumin 28 g/L, actual body weight 55 kg, log10 gamma-glutamyltransferase 2.4, i.e. about 251 U/L) are the only window onto the cohort's central tendency and are used as such in the validation vignette. The therapeutic trough range targeted by the paper's dosing simulations is 2.0-5.0 mg/L."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters: Chantharit 2019 Table 1, "Final model /
    # Estimate value" column. Table 1 is embedded in the abstract page as
    # a raster image and does NOT appear in the PDF text layer, so it is
    # absent from the preprocessed `_trimmed.md`; the values below were
    # read from the extracted image (`pdfimages -png`). The Results prose
    # independently prints the same two point estimates - "Estimated
    # clearance (CL) and volume of distribution (V) values were 7.33 L
    # hour-1 and 439.69 L, respectively" - which corroborates the table
    # read for the two parameters it covers.
    #
    # These are APPARENT ORAL values (CL/F, V/F). The cohort received oral
    # voriconazole exclusively, so bioavailability is not identifiable and
    # no F term is carried; the abstract's equations likewise write plain
    # "CL" and "V". Treating them as apparent is what makes the numbers
    # cohere: V/F = 439.69 L against a typical voriconazole V of roughly
    # 4.6 L/kg (about 250 L at 55 kg) is the expected inflation for F < 1.
    # ------------------------------------------------------------------
    lka <- fixed(log(1)) # Table 1, row 'Ka (fixed)' = 1 in BOTH the final-model and bootstrap columns, 95% CI 'NA'. Held at 1 /h, not estimated; the abstract states 'The linear one-compartment model with first-order absorption and elimination by fixing Ka value at 1.0 well described the data.' Fixing absorption is routine when the data are dominated by TDM troughs, which carry almost no information about the absorption phase.
    label("Absorption rate constant (1/h)")
    lcl <- log(7.33) # Table 1, row 'thetaCL (L/h)' = 7.33 (95% CI 6.78, 7.89; %CV 3.83; bootstrap median 7.27, 95% CI 6.26, 8.22). Also printed in the Results prose.
    label("Apparent oral clearance (L/h)")
    lvc <- log(439.69) # Table 1, row 'thetaV (L)' = 439.69 (95% CI 430.62, 448.77; %CV 1.05; bootstrap median 436.26, 95% CI 150.35, 947.82). Also printed in the Results prose. Note the bootstrap interval is far wider than the asymptotic one, so the volume is much less well determined than the 1.05% CV suggests - unsurprising for a trough-dominated dataset.
    label("Apparent oral volume of distribution (L)")

    # ------------------------------------------------------------------
    # Covariate effects, all from the single final-model equation printed
    # in the Results:
    #
    #   V  (Liter) = thetaV  * (Actual body weight/55)^theta1 * exp(etaV)
    #   CL (L/hr)  = thetaCL * (albumin/28)^theta2 * (logGGT/2.4)^theta3
    #                        * exp(etaCL)
    #
    # with theta1, theta2 and theta3 given numerically in Table 1.
    # ------------------------------------------------------------------
    e_alb_cl <- 0.93 # Table 1, row 'theta2' = 0.93 (95% CI 0.75, 1.12; %CV 10.18; bootstrap median 0.93, 95% CI 0.57, 1.31). Exponent on (albumin / 28 g/L) in clearance. The confidence interval spans 1, so the data are consistent with clearance being simply PROPORTIONAL to albumin; the estimated value is retained here rather than rounded to 1 because Table 1 reports it as estimated, not fixed.
    label("Exponent of (ALB / 28 g/L) on apparent clearance (unitless)")

    # The base of "logGGT" is not stated anywhere in the abstract, and the
    # choice is load-bearing: it decides what GGT value a user must supply
    # to reproduce the typical patient, and it changes clearance for every
    # other patient. It is settled here by the centring constant, which
    # must be the cohort median of the transformed covariate.
    #
    #   base 10 -> median GGT = 10^2.4 = 251 U/L
    #   base e  -> median GGT = e^2.4  =  11 U/L
    #
    # 11 U/L is below the lower bound of the normal adult reference range
    # (roughly 9-48 U/L in women, 12-73 U/L in men). A cohort whose MEDIAN
    # gamma-glutamyltransferase was subnormal could not plausibly have
    # yielded GGT as a retained covariate with an interval as tight as
    # (-0.67, -0.49). 251 U/L - about 4-5x the upper limit of normal - is
    # the ordinary picture for patients with invasive aspergillosis on
    # prolonged voriconazole, itself a recognised hepatotoxin, and makes
    # the covariate's significance intelligible. Base 10 is therefore used.
    # This is an INFERENCE, not a printed fact: it is the one encoding
    # choice in this file that the source does not settle directly, it is
    # flagged in the vignette Errata, and it should be confirmed against
    # the full journal version (Chantharit 2020, Ther Drug Monit
    # 42(6):872-879, doi:10.1097/FTD.0000000000000799), which is not open
    # access and could not be retrieved.
    #
    # Note the form is a ratio OF LOGS, log10(GGT)/2.4, not a log of a
    # ratio - the abstract writes "(logGGT/2.4)", with the division inside
    # the parentheses and applied to the transformed value.
    e_ggt_cl <- -0.58 # Table 1, row 'theta3' = -0.58 (95% CI -0.67, -0.49; %CV -7.95; bootstrap median -0.57, 95% CI -0.97, -0.22). Exponent on (log10(GGT) / 2.4) in clearance. The sign is printed with the estimate and agrees with the bootstrap column.
    label("Exponent of (log10(GGT) / 2.4) on apparent clearance (unitless)")

    e_wt_vc <- fixed(1) # Table 1, row 'theta1 (fixed)' = 1 in BOTH the final-model and bootstrap columns, 95% CI 'NA'. Exponent on (actual body weight / 55 kg) in the apparent volume. Held at unity, so volume scales LINEARLY with weight; this is not an allometric 0.75/1.0 pair because no weight term was retained on clearance at all.
    label("Exponent of (WT / 55 kg) on apparent volume of distribution (unitless)")

    # ------------------------------------------------------------------
    # Interindividual variability. Table 1 labels these rows 'omega^2CL'
    # and 'omega^2V' and annotates each with '(omega x 100)', so the first
    # number is the VARIANCE and the parenthesised one is 100 x the
    # standard deviation. The arithmetic confirms the reading and fixes
    # the scale: sqrt(0.190) = 0.4359 -> 43.58, and sqrt(0.080) = 0.28284
    # -> 28.28, both matching the printed parentheses exactly. The
    # variances are therefore what go here, and the exponential form is
    # explicit in the printed equations as exp(etaCL) and exp(etaV).
    # ------------------------------------------------------------------
    etalcl ~ 0.190 # Table 1, row 'omega^2CL (omegaCL x 100)' = 0.190 (43.58); 95% CI 0.147, 0.234; %CV 11.58; bootstrap median 0.181 (42.650), 95% CI 0.120, 0.253. Reported eta-shrinkage 0.19.
    etalvc ~ 0.080 # Table 1, row 'omega^2V (omegaV x 100)' = 0.080 (28.28); 95% CI 0.063, 0.098; %CV 11.083; bootstrap median 0.080 (28.28), 95% CI 0.002, 0.503. Reported eta-shrinkage 0.92 - very high, so the individual volume estimates are largely driven back to the typical value and the volume IIV should be treated as weakly identified.

    # ------------------------------------------------------------------
    # Residual error. Table 1's last row is 'Residual variability(sigma)'
    # = 0.59 (95% CI 0.54, 0.64; %CV 4.21; bootstrap median 0.58, 95% CI
    # 0.51, 0.66). The symbol is sigma, not sigma^2, and the two
    # interindividual rows in the same table are explicitly labelled
    # omega^2, so 0.59 is a STANDARD DEVIATION on its own scale.
    #
    # The abstract never names the error MODEL. Figure 1, the visual
    # predictive check, settles it. Its blue 5th-percentile prediction
    # line lies on the x-axis (about 0.03 mg/L) while the blue median sits
    # near 2.4 mg/L - a 5%/50% ratio of roughly 0.012. Simulating this
    # model's own IIV at 200 mg q12h and applying each candidate error
    # form gives:
    #
    #   proportional (CV 0.59)  5%/50% = 0.029   4.5% of draws negative
    #   additive     (0.59 mg/L) 5%/50% = 0.291
    #   exponential  (sigma 0.59) 5%/50% = 0.290
    #
    # Only the proportional form drives the lower prediction interval onto
    # the axis, and it does so for a concrete reason: a 59% proportional
    # CV makes about 4.5% of simulated observations negative, so the 5th
    # percentile of the predicted distribution sits essentially at zero -
    # exactly the flat-at-zero blue line Figure 1 shows. Additive and
    # exponential errors both place the 5th percentile near 0.64 mg/L,
    # which would be plainly visible well above the axis and is not what
    # the figure shows. The vignette re-runs this comparison as a gate.
    # ------------------------------------------------------------------
    propSd <- 0.59 # Table 1, row 'Residual variability(sigma)' = 0.59 (95% CI 0.54, 0.64; %CV 4.21; bootstrap median 0.58, 95% CI 0.51, 0.66); proportional form adjudicated against the Figure 1 visual predictive check, see above.
    label("Proportional residual error (fraction)")
  })

  model({
    ka <- exp(lka)

    # Clearance: the Results final-model equation. Both covariates are
    # median-normalised, albumin directly and gamma-glutamyltransferase
    # through its base-10 logarithm.
    cl <- exp(lcl + etalcl) *
      (ALB / 28)^e_alb_cl *
      (log10(GGT) / 2.4)^e_ggt_cl

    # Apparent volume: linear in actual body weight (exponent fixed at 1).
    vc <- exp(lvc + etalvc) * (WT / 55)^e_wt_vc

    kel <- cl / vc

    # One-compartment disposition with a first-order absorption depot.
    # Dosing is oral into `depot`; no bioavailability term is carried
    # because the cohort received no intravenous voriconazole and F is
    # therefore not identifiable - `cl` and `vc` are apparent oral values.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Dose in mg divided by apparent volume in L gives mg/L, the unit of
    # the reported concentrations (equivalently ug/mL).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
