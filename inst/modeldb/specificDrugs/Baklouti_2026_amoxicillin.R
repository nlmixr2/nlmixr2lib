Baklouti_2026_amoxicillin <- function() {
  description <- paste(
    "Population PK model of oral amoxicillin in breastfeeding women,",
    "describing maternal plasma and breast-milk concentrations",
    "simultaneously. A gastrointestinal depot feeds a central (plasma)",
    "compartment by first-order absorption; the central compartment is",
    "drained both by first-order systemic elimination and by a",
    "UNIDIRECTIONAL first-order transfer into a breast-milk compartment,",
    "and milk is cleared by its own first-order loss constant. Nothing",
    "returns from milk to plasma: the authors tested a bidirectional",
    "peripheral-milk parameterisation and it failed to converge. The milk",
    "volume is not identifiable from the data and was fixed equal to the",
    "central volume, so the milk state is scaled by V/F. Because the milk",
    "loss constant runs continuously rather than only at feeds, the model",
    "understates milk amounts between feeds and is therefore a",
    "conservative (low) estimator of the relative infant dose; the authors",
    "flag this as its main structural limitation. Inter-individual",
    "variability on all five structural parameters, separate proportional",
    "residual errors for plasma and for milk, and no retained covariates.",
    sep = " "
  )
  reference <- paste(
    "Baklouti S, Rigourd V, Panchaud A, Nordeng H, Allegaert K, Annaert P,",
    "Huang MC, Monfort A, Guidi M, Gandia P.",
    "Population pharmacokinetic modelling of amoxicillin in human breast",
    "milk - A contribution from the ConcePTION project.",
    "Br J Clin Pharmacol. 2026;92(6):1833-1844. doi:10.1002/bcp.70434.",
    "All parameter values are the final estimates of Table 3; the ODE",
    "system and the Vmilk = Vc assumption are the display equations and",
    "prose of Results section 3.2 (pp. 1836-1837), and the model topology",
    "is Figure 1.",
    sep = " "
  )
  vignette <- "Baklouti_2026_amoxicillin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # No covariate was retained in the final model (Results section 3.2,
  # "No socio-demographic or clinico-biological covariates were included
  # in the final model."). Everything the authors screened is documented
  # in covariatesDataExcluded below.
  covariateData <- list()

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Maternal age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened by the forward-addition / backward-elimination covariate search and not retained (Baklouti 2026 Methods section 2.3 and Results section 3.2). Cohort mean 35 years, SD 4, CV 11% (Table 1)."
    ),
    WT = list(
      description = "Maternal body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened and not retained (Methods section 2.3; Discussion section 4.1, 'The final model did not retain any of the documented variables'). Cohort mean 64.9 kg, SD 12.1, CV 18.6% (Table 1). Still needed OUTSIDE the PK model to compute the relative infant dose, which normalises the maternal dose per kg of maternal body weight; the simulations of Table 4 assume 65 kg (SD 12)."
    ),
    HT = list(
      description = "Maternal body height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened and not retained (Methods section 2.3). Cohort mean 164.5 cm, SD 5.9, CV 3.6% (Table 1)."
    ),
    CRCL = list(
      description = "Maternal glomerular filtration rate, estimated with the CKD-EPI equation and BSA-normalised",
      units       = "mL/min/1.73 m^2",
      type        = "continuous",
      notes       = "Screened and not retained (Methods section 2.3, 'maternal glomerular filtration rate (estimated GFR using CKD-EPI formula)'). Cohort mean 114.7 mL/min/1.73 m^2, SD 8.5, CV 7.4% (Table 1). The Discussion (section 4.1) attributes the null result to insufficient covariate spread rather than to absence of a relationship: 'This is typically the case for amoxicillin elimination clearance, which is primarily governed by renal function, given the low eGFR IIV (CV = 7.5%) in our study.' Source column eGFR."
    ),
    ALB = list(
      description = "Maternal serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened and not retained (Methods section 2.3, 'maternal albuminemia (g/L, continuous variable)'). No cohort summary is reported in Table 1."
    ),
    BREASTFEED_EXCLUSIVE = list(
      description = "Whether the mother was breastfeeding exclusively",
      units       = "(binary)",
      type        = "categorical",
      notes       = "Screened and not retained (Methods section 2.3, 'breastfeeding exclusivity (yes/no, categorical variable)'). The per-category counts are not reported. Reference category 0 (not exclusive)."
    ),
    AGE_INFANT = list(
      description = "Postnatal age of the breastfed infant, screened as a proxy for milk maturity",
      units       = "months",
      type        = "continuous",
      notes       = "Screened and not retained (Methods section 2.3: 'Since the infant's age reflects milk maturity - which can influence its composition and consequently amoxicillin concentrations - the infant's age (in months), body weight (kg) and height (cm) were recorded on the day of the study'). Cohort mean 5.5 months, SD 4.9, CV 89% (Table 2). Enrolment required the infant to be over 4 weeks old and breastfeeding to have started at least 2 weeks after birth, so all milk sampled was mature milk."
    ),
    WT_INFANT = list(
      description = "Body weight of the breastfed infant",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened and not retained (Methods section 2.3). Cohort mean 6.7 kg, SD 2.7, CV 40.2% (Table 2). Used outside the PK model in the relative-infant-dose calculation via the assumed daily milk intake of 150 mL per kg of infant body weight (Methods section 2.4). Infant HEIGHT (mean 62.0 cm, SD 8.9, Table 2) was screened alongside it and likewise not retained; it is recorded here rather than under a canonical of its own because no model in the library references it."
    )
  )

  compartmentData <- list(
    depot   = list(analyte = "amoxicillin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "amoxicillin", units = "mg", specimen = "plasma", verified = TRUE),
    milk    = list(analyte = "amoxicillin", units = "mg", specimen = "milk", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 25,
    n_studies      = 1,
    n_observations = 150,
    age_mean       = "35 years",
    weight_mean    = "64.9 kg",
    height_mean    = "164.5 cm",
    sex_female_pct = 100,
    renal_function = "Normal; CKD-EPI eGFR mean 114.7 mL/min/1.73 m^2 (SD 8.5, CV 7.4%). The narrow spread is why no renal covariate could be identified.",
    disease_state  = "Breastfeeding women treated with oral immediate-release amoxicillin, with or without clavulanic acid, for at least 2 days for an approved indication prescribed in routine community care. Mothers of preterm infants (gestational age under 34 weeks), of twins, or requiring specialised medical supervision were excluded.",
    dose_range     = "Oral immediate-release amoxicillin 1 g twice daily (n = 6) or 1 g three times daily (n = 19); mean total daily dose 2760 mg (SD 436). No specific brand was mandated. Sampling was at steady state after at least 2 days of treatment.",
    regions        = "France (recruited from the lactarium of Necker-Enfants Malades Hospital, Paris). EudraCT 2021-002247-30; French ethics approval 21.01446.000014.",
    infant_partner = "One breastfed infant per mother, all over 4 weeks old: 12 male / 13 female, mean postnatal age 5.5 months (SD 4.9), mean weight 6.7 kg (SD 2.7), mean height 62.0 cm (SD 8.9). Infant characteristics enter the relative-infant-dose calculation, not the PK model. No adverse events were observed in any breastfed infant.",
    feeding_pattern = "Mothers reported a mean of 9 breastfeeds per day (SD 4, CV 47%; Table 1). The steady-state Monte Carlo simulations of Table 4 instead draw a feeding frequency per woman with mean 11 per day over the range 6-18, holding total daily milk intake constant, and the relative infant dose assumes a daily milk intake of 150 mL per kg of infant body weight (Methods section 2.4).",
    notes          = "Baseline demographics from Baklouti 2026 Tables 1 and 2. Three blood and three milk samples per mother (75 + 75 = 150 observations), all drawn at home within one dosing interval: nominally 15-30 min, 1-2 h and 3-4 h post-dose. Blood and milk could not be drawn at exactly the same minute, so pairs are near-simultaneous rather than synchronous. Milk was collected with an electric breast pump, about 20 mL of FOREMILK per sample without emptying the breast; because amoxicillin is water-soluble and concentrates in the low-fat foremilk, the authors describe this as a deliberate worst-case sampling scheme for milk exposure. Assay LC-MS/MS with LLOQ 0.2 mg/L (plasma) and 0.01 mg/L (milk); no sample fell below either limit and no data were missing. Three mothers took other drugs (one antihypertensive, one anti-inflammatory, one another antibiotic); one had a history of hypertension; none smoked or drank alcohol."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters. Every value is a FINAL estimate from
    # Baklouti 2026 Table 3 ("Final parameter estimates of amoxicillin
    # breast milk popPK model, n = 25 breastfeeding mothers, n = 75 blood
    # samples and n = 75 milk samples"), fitted in MonolixSuite 2024R1.
    #
    # The model is parameterised on RATE CONSTANTS plus one volume, not
    # on clearances: the paper estimates ke and V/F separately with an
    # independent eta on each. Reparameterising to lcl + lvc would force
    # etalkel = etalcl - etalvc, i.e. a correlated eta block the authors
    # did not fit, so it would change the model. lkel is therefore used
    # alongside an explicit lvc (operator sidecar oare_PMC13206287
    # request-001 / response-001, question q2, option A).
    #
    # V/F is apparent: amoxicillin was given orally and F was not
    # estimated, so bioavailability is folded into the volume and no
    # lfdepot appears below. The paper's G(0) = F * D initial condition
    # is therefore realised as a plain dose into depot.
    # ------------------------------------------------------------------
    lka <- log(0.17)
    label("First-order absorption rate constant from the gastrointestinal depot, ka (1/h)")                          # Table 3 'k a (h-1) 0.17, RSE 39.9%, 95%CI 0.085-0.35'
    lvc <- log(82.38)
    label("Apparent central (plasma) volume of distribution, V/F (L)")                                               # Table 3 'V/F (L) 82.38, RSE 11.0%, 95%CI 66.46-102.1'
    lkel <- log(0.31)
    label("First-order systemic elimination rate constant from the central compartment, ke (1/h)")                   # Table 3 'K e (h-1) 0.31, RSE 23.1%, 95%CI 0.2-0.48'

    # Unidirectional plasma-to-milk transfer. Member of the registered
    # k_<from>_<to> directional-transfer family (see the k_central_elf
    # entry in inst/references/parameter-names.md, whose Notes cite
    # k_central_milk / k_milk_central in Wattanakul_2024_primaquine.R).
    # Here it is ESTIMATED rather than derived from a q / vc pair, so it
    # takes the log-transformed ini() form.
    lk_central_milk <- log(0.028)
    label("First-order transfer rate constant from the central compartment to milk, kmilk (1/h)")                    # Table 3 'K milk (h-1) 0.028, RSE 12.5%, 95%CI 0.022-0.036'

    # First-order loss out of the milk compartment. Registered as a
    # member of the lkeff_<pool> / keff_<pool> efflux family alongside
    # lkeff_rbc and lkeff_pbmc (operator sidecar oare_PMC13206287
    # request-001 / response-001, question q1, option A). Lumped in the
    # same sense as its siblings: it combines removal of milk from the
    # breast by feeding with milk turnover. Nothing returns to plasma,
    # so this is NOT k_milk_central.
    lkeff_milk <- log(0.33)
    label("First-order elimination rate constant of amoxicillin out of the milk compartment, kmilk_e (1/h)")         # Table 3 'K milk_e (h-1) 0.33, RSE 19.4%, 95%CI 0.23-0.48'

    # ------------------------------------------------------------------
    # Inter-individual variability, lognormal on every structural
    # parameter (Methods section 2.3, "IIV was evaluated assuming
    # lognormally distributed individual parameters").
    #
    # SCALE. Monolix reports omega as a STANDARD DEVIATION on the log
    # scale, and Table 3's IIV rows are those omegas, not variances.
    # Confirmed against an independently printed number: the Discussion
    # (section 4.1) quotes "a mean ka of 0.17 h-1 (CV = 37%)" against
    # Table 3's IIV ka of 0.37, i.e. the paper quotes CV% as omega * 100.
    # Reading 0.37 as a variance would give omega = 0.608 and CV = 66%,
    # which contradicts the printed 37%. rxode2 wants VARIANCES, so each
    # value below is the tabulated omega SQUARED.
    #
    # COVERAGE. Methods section 2.3 says "IIV was included in ka, ke,
    # kmilk and kmilk_e" -- omitting V/F -- but Table 3 tabulates
    # "IIV V/F (se) 0.086" with RSE 47.7% and a 95% CI of 0.038-0.2. A
    # row carrying an RSE and a confidence interval is an estimated
    # component, so the table is taken as authoritative and the eta on
    # V/F is retained; the Methods sentence is simply incomplete.
    # ------------------------------------------------------------------
    etalka ~ 0.1369                # Table 3 'IIV k a (se) 0.37, RSE 44.2%, 95%CI 0.17-0.8'; 0.37^2
    etalvc ~ 0.007396              # Table 3 'IIV V/F (se) 0.086, RSE 47.7%, 95%CI 0.038-0.2'; 0.086^2
    etalkel ~ 0.2401               # Table 3 'IIV K e (se) 0.49, RSE 24.9%, 95%CI 0.31-0.78'; 0.49^2
    etalk_central_milk ~ 0.0121    # Table 3 'IIV K milk (se) 0.11, RSE 64.8%, 95%CI 0.039-0.31'; 0.11^2
    etalkeff_milk ~ 0.1444         # Table 3 'IIV K milk_e (se) 0.38, RSE 28.3%, 95%CI 0.22-0.64'; 0.38^2

    # ------------------------------------------------------------------
    # Residual error: a separate proportional model for each matrix
    # (Results section 3.2, "A proportional residual error model was used
    # for plasma and milk data"). Monolix's b coefficient multiplies the
    # prediction, y = f + b * f * e with e ~ N(0, 1), which is exactly
    # nlmixr2's prop().
    # ------------------------------------------------------------------
    propSd <- 0.48
    label("Proportional residual error for maternal plasma amoxicillin (fraction)")                                  # Table 3 'b1 (se) 0.48, RSE 10.9%, 95%CI 0.39-0.59'; Note: 'b1 and b2 are the proportional errors of plasma and milk concentrations, respectively'
    propSd_Cmilk <- 0.26
    label("Proportional residual error for breast-milk amoxicillin (fraction)")                                      # Table 3 'b2 (se) 0.26, RSE 11.3%, 95%CI 0.21-0.33'
  })

  model({
    # ---- 1. Individual parameters -------------------------------------
    ka <- exp(lka + etalka)
    vc <- exp(lvc + etalvc)
    kel <- exp(lkel + etalkel)
    k_central_milk <- exp(lk_central_milk + etalk_central_milk)
    keff_milk <- exp(lkeff_milk + etalkeff_milk)

    # The milk volume is not identifiable from these data and was fixed
    # equal to the central volume (Results section 3.2, p. 1837: "The
    # actual milk volume could not be determined ... Consequently, the
    # milk volume (Vmilk) was arbitrarily set to a value identical to
    # that of the systemic compartment (Vc)."). It is an identity rather
    # than a parameter, so it carries its own eta only through vc. The
    # Discussion (section 4.1) acknowledges this is not physiological -
    # real daily milk output averages about 2 L - and defends it as
    # avoiding over-parameterisation.
    vmilk <- vc

    # ---- 2. ODE system ------------------------------------------------
    # Baklouti 2026 Results section 3.2, display equations on p. 1836,
    # topology in Figure 1. Transcribed from the PDF with
    # `pdftotext -layout`; the trimmed markdown renders them as
    # <!-- formula-not-decoded -->, and the publisher's symbol font
    # encodes the minus and multiplication operators as C0 control bytes
    # (0x03 and 0x04), so the signs below were decoded with `cat -A`:
    #
    #   dG(t)/dt = -ka * G(t)                            G(0) = F * D
    #   dX(t)/dt =  ka * G(t) - ke * X(t) - kmilk * X(t) X(0) = 0
    #   dM(t)/dt =  kmilk * X(t) - kmilk_e * M(t)        M(0) = 0
    #
    # Note that kmilk is subtracted from the central compartment as well
    # as added to milk, so transfer into milk is a genuine elimination
    # pathway for plasma; total plasma loss is ke + kmilk = 0.338 1/h.
    # Nothing returns from milk (Results section 3.2: "incorporating
    # transfer from plasma to milk but without a transfer from milk to
    # plasma. Instead, an elimination rate constant was applied to the
    # milk compartment.").
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k_central_milk * central
    d/dt(milk) <- k_central_milk * central - keff_milk * milk

    # ---- 3. Observations and residual error ---------------------------
    # Amounts are in mg and volumes in L, so both concentrations are in
    # mg/L, the reporting unit of the assay and of Table 4.
    Cc <- central / vc
    Cmilk <- milk / vmilk

    Cc ~ prop(propSd)
    Cmilk ~ prop(propSd_Cmilk)
  })
}
