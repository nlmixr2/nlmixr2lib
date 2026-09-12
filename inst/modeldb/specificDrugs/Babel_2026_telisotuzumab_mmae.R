Babel_2026_telisotuzumab_mmae <- function() {
  description <- paste0(
    "Population pharmacokinetic model of the UNCONJUGATED monomethyl ",
    "auristatin E (MMAE) payload released from telisotuzumab vedotin ",
    "(Teliso-V) in adults with c-Met protein overexpressing advanced ",
    "solid tumours (Babel 2026, n = 304 pooled from a phase 1 study ",
    "and the LUMINOSITY phase 2 study). One-compartment linear ",
    "disposition fed by a first-order deconjugation rate ka from a ",
    "depot that represents the administered antibody-drug conjugate ",
    "(Babel 2026 Figure S1B). Interindividual variability on CL, Vc ",
    "and ka; combined proportional plus additive residual error. ",
    "Baseline albumin and renal function act on CL; age, baseline ",
    "albumin and race act on ka; body weight acts on Vc. Because ka ",
    "is much smaller than CL/Vc the payload profile is flip-flop and ",
    "its apparent terminal half-life is set by deconjugation, not by ",
    "MMAE elimination. Babel 2026 developed this model INDEPENDENTLY ",
    "of the conjugate model rather than as a coupled system, so the ",
    "depot is dosed directly and the conjugate model is not a ",
    "prerequisite; the companion conjugate model is ",
    "Babel_2026_telisotuzumab. NOTE: the source does not report the ",
    "drug-antibody ratio, the molecular weights, or the systemically ",
    "available payload fraction needed to convert a milligram dose of ",
    "conjugate into the payload mass this depot receives, so the ",
    "depot must be dosed in MMAE-equivalent mass; see the vignette ",
    "Assumptions and deviations section."
  )
  reference <- paste(
    "Babel H, Brunsdon P, Engelhardt B, Schmitt V, Ratajczak C, Mensing S,",
    "Menon RM, Parikh A. Population pharmacokinetics and exposure-response",
    "analyses for telisotuzumab vedotin in patients with c-Met protein",
    "overexpressing tumors.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15(1):e70219.",
    "doi:10.1002/psp4.70219. PMCID PMC12945708.",
    "All parameter values are from Data S1 (Supporting Information) Table S7,",
    "'Final Model Parameter Estimates and Variability of Teliso-V Conjugate",
    "and Unconjugated MMAE Payload Pharmacokinetics',",
    "Unconjugated MMAE Payload block; the structure is Figure S1B.",
    sep = " "
  )
  vignette <- "Babel_2026_telisotuzumab"

  units <- list(
    time          = "day",
    dosing        = "mg (MMAE equivalents entering the depot, NOT milligrams of conjugate)",
    concentration = "ug/mL"
  )

  # CL is in L/day and Vc in L (Table S7), so a depot amount in mg gives
  # a concentration in mg/L == ug/mL. Babel 2026 plots the payload on a
  # ng/mL axis (Figure 5, CavgMMAE 0-4 ng/mL), which is 1000 * Cc. The
  # same ug/mL-internal / ng/mL-reported convention is used by the
  # sibling vedotin payload model Choules_2024_enfortumab.
  compartmentData <- list(
    depot   = list(analyte = "telisotuzumab vedotin conjugate, as MMAE equivalents awaiting deconjugation", units = "mg", specimen = "not applicable", verified = FALSE),
    central = list(analyte = "unconjugated monomethyl auristatin E (MMAE)", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Normalised to the overall-population median of 68.9 kg inside model() (Babel 2026 Table S4, All Participants; corroborated by the Table S3 forest-plot reference group 'Body Weight <= 68.9 kg'). Acts on Vc only in the payload model. Babel 2026 Results: 'body weight was significantly correlated with Vc in the MMAE population PK model'.",
      source_name        = "Body Weight"
    ),
    AGE = list(
      description        = "Age at baseline.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Normalised to the overall-population median of 65 years inside model() (Babel 2026 Table S4, All Participants median 65, range 30-87). Acts on the deconjugation rate ka only, with a negative exponent, so older patients deconjugate more slowly.",
      source_name        = "Age"
    ),
    ALB = list(
      description        = "Baseline serum albumin concentration.",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Normalised to the overall-population median of 41.6 g/L inside model() (Babel 2026 Table S4, All Participants median 41.6 g/L, range 29.0-52.0). Acts on both CL and ka, and with opposite signs: the CL exponent is +1.99 and the ka exponent is -1.01, so a higher albumin both clears MMAE faster and releases it more slowly. Both push payload exposure down, which is the direction Babel 2026 Figure 1B shows (albumin > 41.6 g/L versus <= 41.6 g/L gives Cmax 0.663 and AUCtau 0.652) and which the Discussion states in words: 'unconjugated MMAE payload exposure was lower in patients with higher baseline albumin'.",
      source_name        = "Baseline Albumin"
    ),
    RACE_BLACK = list(
      description        = "Black or African American race indicator; 1 = Black or African American, 0 otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (White; Babel 2026 Table S3 reference group 'Race: White')",
      notes              = "Multiplicative factor on the deconjugation rate ka. Babel 2026 Table S4: 8 of 304 (3%) Black or African American. Figure 1B reports n = 8 versus n = 198 White for this comparison, and the authors caution that the small n limits conclusions.",
      source_name        = "Race: Black or African American"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator; 1 = Asian, 0 otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (White; Babel 2026 Table S3 reference group 'Race: White')",
      notes              = "Multiplicative factor on the deconjugation rate ka. Babel 2026 Table S4: 98 of 304 (32%) Asian. The Discussion quantifies the consequence: Asian patients had 'about 25% lower unconjugated MMAE exposures' than White patients, still inside the White exposure range. RACE_BLACK and RACE_ASIAN are mutually exclusive; a White patient carries 0 for both.",
      source_name        = "Race: Asian"
    ),
    RENALIMP_MILD = list(
      description        = "Mild renal impairment indicator; 1 = mild, 0 otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal renal function; Babel 2026 Table S7 row 'Mild vs. Normal Renal Impairment on CL')",
      notes              = "Multiplicative factor on MMAE CL. Babel 2026 Table S2 footnote a defines the strata on estimated glomerular filtration rate computed by Cockcroft-Gault: normal at least 90, mild 60-90, moderate 30-60 and severe below 30 mL/min/1.73 m2. Table S4: 132 of 304 (43%) mild versus 110 of 304 (36%) normal. Mutually exclusive with RENALIMP_MOD and RENALIMP_SEV; a patient with normal renal function carries 0 for all three.",
      source_name        = "Baseline Renal Function: Mild Impairment"
    ),
    RENALIMP_MOD = list(
      description        = "Moderate renal impairment indicator; 1 = moderate, 0 otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal renal function)",
      notes              = "Babel 2026 POOLED moderate and severe renal impairment into a single stratum for estimation (Table S7 row 'Moderate/Severe vs. Normal Renal Impairment on CL'), because only 2 of 304 patients had severe impairment. The single published factor is therefore applied to RENALIMP_MOD and RENALIMP_SEV alike inside model(); the two indicators are kept separate here so that a user's data can retain the clinical distinction the source could not estimate. Table S4: 59 of 304 (19%) moderate. Mutually exclusive with RENALIMP_MILD and RENALIMP_SEV.",
      source_name        = "Baseline Renal Function: Moderate Impairment"
    ),
    RENALIMP_SEV = list(
      description        = "Severe renal impairment indicator; 1 = severe, 0 otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal renal function)",
      notes              = "Carries the SAME published factor as RENALIMP_MOD because Babel 2026 estimated a single pooled 'Moderate/Severe vs. Normal' effect; see the RENALIMP_MOD notes. Table S4: 2 of 304 (1%) severe. The Discussion is explicit that this stratum is not informative on its own: 'conclusions regarding severe renal impairment on unconjugated MMAE payload exposure are limited, given the small number of patients with severe renal impairment (n = 2)'. Mutually exclusive with RENALIMP_MILD and RENALIMP_MOD.",
      source_name        = "Baseline Renal Function: Severe Impairment"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex indicator; 1 = female, 0 = male.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on payload CL, Vc and ka in the Babel 2026 Table S2 covariate sweep but not retained in the final payload model (Table S7 lists no sex term for the payload). Sex WAS retained on conjugate Vc; see Babel_2026_telisotuzumab."
    ),
    ADA_POS = list(
      description = "Treatment-emergent anti-drug antibody status; 1 = positive, 0 = negative.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Babel 2026 Table S2 lists treatment-emergent ADA status as a covariate of interest for conjugate clearance only, not for the payload, and the Discussion confirms that treatment-emergent ADAs were 'significant covariates on Teliso-V conjugate CL, but not MMAE CL'. Documented here to preserve the covariate screen."
    ),
    CONMED_CYP3A_INHIB = list(
      description = "Concomitant strong CYP3A inhibitor indicator; 1 = coadministered, 0 = not.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on payload CL (Babel 2026 Table S2, 'Concomitant medications (strong CYP3A inhibitors, strong CYP3A inducers)') but not retained in the final model, so no point estimate exists anywhere on disk. MMAE is a CYP3A4 substrate, which is why the sweep included it. The sibling vedotin payload model Choules_2024_enfortumab does carry a CYP3A perpetrator effect, back-calculated from a dedicated drug-interaction simulation rather than estimated from patient data."
    ),
    CONMED_CYP3A_IND = list(
      description = "Concomitant strong CYP3A inducer indicator; 1 = coadministered, 0 = not.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on payload CL (Babel 2026 Table S2) but not retained; see CONMED_CYP3A_INHIB."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 304L,
    n_studies      = 2L,
    age_range      = "median 65 years, range 30-87 (Babel 2026 Table S4, All Participants)",
    weight_range   = "median 68.9 kg, range 36.0-144 (Babel 2026 Table S4, All Participants)",
    sex_female_pct = 37.8,
    race_ethnicity = c(White = 65, Asian = 32, `Black or African American` = 3),
    disease_state  = "Advanced solid tumours likely to express c-Met (phase 1, NCT02099058, n = 35) and locally advanced or metastatic c-Met protein overexpressing non-small cell lung cancer (LUMINOSITY phase 2, NCT03539536, n = 269).",
    dose_range     = "Phase 1: 0.15-3.3 mg/kg every 3 weeks and 1.6-2.2 mg/kg every 2 weeks. LUMINOSITY: 1.6 or 1.9 mg/kg every 2 weeks of telisotuzumab vedotin conjugate.",
    regions        = "Europe 26%, North America 26%, Asia 28%, rest of world 20% (Babel 2026 Table S4)",
    renal_function = "normal 110 of 304 (36%), mild 132 (43%), moderate 59 (19%), severe 2 (1%), missing 1 (Babel 2026 Table S4)",
    notes          = paste0(
      "The payload analysis uses the same 304 patients as the ",
      "conjugate analysis. 4.69% of payload records were below the ",
      "limit of quantitation (Babel 2026 Data S1 Methods). Payload ",
      "parameters were estimated with relative standard errors of at ",
      "most 29.3%. Shrinkage was 4.64% on CL, 14.2% on Vc and 18.2% ",
      "on ka."
    )
  )

  ini({
    # ==================================================================
    # Structural disposition. Values are the population estimates in
    # Babel 2026 Data S1 Table S7, Unconjugated MMAE Payload block,
    # reported on the natural scale with %RSE and 95% CI and therefore
    # wrapped in log() here.
    #
    # ka is the DECONJUGATION rate, not an absorption rate: Babel 2026
    # Methods describes the payload model as 'a first-order release of
    # MMAE from the ADC structure', and Figure S1B draws it as a single
    # MMAE compartment with an inbound arrow labelled Ka and an
    # outbound clearance arrow. Structurally this is identical to a
    # first-order-input one-compartment model, which is how it is
    # encoded here.
    #
    # ka (0.154 /day, half-life 4.50 days) is well below kel
    # (76.3 / 96.0 = 0.795 /day, half-life 0.87 days), so the system is
    # flip-flop and the observed payload terminal slope reflects
    # deconjugation.
    # ==================================================================
    lka <- log(0.154); label("First-order deconjugation rate of MMAE from the conjugate (1/day)")  # Table S7 'Ka (1/day)' 0.154, %RSE 2.58, 95% CI (0.147, 0.162)
    lcl <- log(76.3);  label("Unconjugated MMAE clearance (L/day)")                                # Table S7 'CL (L/day)' 76.3, %RSE 4.83, 95% CI (69.4, 83.9)
    lvc <- log(96.0);  label("Unconjugated MMAE volume of distribution (L)")                       # Table S7 'Vc (L)' 96.0, %RSE 5.62, 95% CI (86.0, 107)

    # ==================================================================
    # Covariate effects. Continuous covariates are power functions of
    # covariate / overall-population median and categorical covariates
    # are multiplicative factors relative to the reference group
    # (Babel 2026 Methods).
    # ==================================================================
    e_alb_cl         <- 1.99;   label("Power exponent on (ALB/41.6 g/L) for MMAE CL (unitless)")   # Table S7 'Albumin on CL' 1.99, %RSE 10.5, 95% CI (1.58, 2.39)
    e_renalmild_cl   <- 0.842;  label("Multiplicative factor on MMAE CL for mild renal impairment versus normal (unitless)")            # Table S7 'Mild vs. Normal Renal Impairment on CL' 0.842, %RSE 5.58, 95% CI (0.755, 0.939)
    e_renalmodsev_cl <- 0.755;  label("Multiplicative factor on MMAE CL for moderate or severe renal impairment versus normal (unitless)") # Table S7 'Moderate/Severe vs. Normal Renal Impairment on CL' 0.755, %RSE 7.30, 95% CI (0.655, 0.871)
    e_age_ka         <- -0.454; label("Power exponent on (AGE/65 years) for the deconjugation rate (unitless)")   # Table S7 'Age on Ka' -0.454, %RSE 25.4, 95% CI (-0.680, -0.228)
    e_alb_ka         <- -1.01;  label("Power exponent on (ALB/41.6 g/L) for the deconjugation rate (unitless)")   # Table S7 'Albumin on Ka' -1.01, %RSE 18.4, 95% CI (-1.38, -0.648)
    e_black_ka       <- 0.853;  label("Multiplicative factor on the deconjugation rate for Black or African American versus White (unitless)") # Table S7 'Black or African American vs. White on Ka' 0.853, %RSE 15.5, 95% CI (0.631, 1.15)
    e_asian_ka       <- 0.874;  label("Multiplicative factor on the deconjugation rate for Asian versus White (unitless)")                     # Table S7 'Asian vs. White on Ka' 0.874, %RSE 4.44, 95% CI (0.802, 0.954)
    e_wt_vc          <- 0.612;  label("Power exponent on (WT/68.9 kg) for MMAE Vc (unitless)")     # Table S7 'Body Weight on Vc' 0.612, %RSE 29.3, 95% CI (0.261, 0.964)

    # ==================================================================
    # Interindividual variability. As in the conjugate block, the
    # 'Population Estimate' column for the IIV rows is the VARIANCE:
    # the Table S7 footnote defines %CV as SQRT(exp(omega2)-1)*100, and
    # sqrt(exp(0.240) - 1) * 100 = 52.1%, sqrt(exp(0.486) - 1) * 100 =
    # 79.1% and sqrt(exp(0.0735) - 1) * 100 = 27.6%, reproducing the
    # printed %CV column exactly.
    #
    # Babel 2026 Results states that the payload model included a
    # CORRELATION between the CL and Vc random effects, but Table S7
    # tabulates only the three diagonal variances and no covariance or
    # correlation coefficient. The correlation magnitude is therefore
    # unreported and the block is left diagonal here rather than
    # invented; see the vignette Assumptions and deviations section.
    # ==================================================================
    etalcl ~ 0.240   # Table S7 'IIV on CL' variance 0.240, 52.1 %CV, 4.64% shrinkage
    etalvc ~ 0.486   # Table S7 'IIV on Vc' variance 0.486, 79.1 %CV, 14.2% shrinkage
    etalka ~ 0.0735  # Table S7 'IIV on Ka' variance 0.0735, 27.6 %CV, 18.2% shrinkage

    # ==================================================================
    # Residual unexplained variability, entered as square roots because
    # Table S7 reports variances and nlmixr2 error terms take standard
    # deviations. The additive term is effectively zero
    # (sqrt(7.53e-10) = 2.744e-5 ug/mL, i.e. 0.027 ng/mL), so the
    # payload residual model is proportional in all but name; it is
    # retained because the source estimated it with a 4.72% RSE.
    # ==================================================================
    propSd <- sqrt(0.0833);    label("Proportional residual error on Cc (fraction)")   # Table S7 'Proportional Error (Variance)' 0.0833, %RSE 0.849, 95% CI (0.0819, 0.0847)
    addSd  <- sqrt(7.53e-10);  label("Additive residual error on Cc (ug/mL)")          # Table S7 'Additive Error (Variance)' 7.53e-10, %RSE 4.72, 95% CI (6.84e-10, 8.23e-10)
  })

  model({
    # ----- Individual parameters -----
    # The published 'Moderate/Severe vs. Normal' factor is applied to
    # the moderate and the severe indicator alike, because Babel 2026
    # estimated the two strata as one pooled group. The indicators are
    # mutually exclusive, so the exponent is 0 or 1 for any subject.
    cl <- exp(lcl + etalcl) *
      (ALB / 41.6)^e_alb_cl *
      e_renalmild_cl^RENALIMP_MILD *
      e_renalmodsev_cl^(RENALIMP_MOD + RENALIMP_SEV)
    vc <- exp(lvc + etalvc) * (WT / 68.9)^e_wt_vc
    ka <- exp(lka + etalka) *
      (AGE / 65)^e_age_ka *
      (ALB / 41.6)^e_alb_ka *
      e_black_ka^RACE_BLACK *
      e_asian_ka^RACE_ASIAN

    # ----- Micro-constant -----
    kel <- cl / vc

    # ----- One-compartment disposition with first-order deconjugation input -----
    # (Babel 2026 Figure S1B.) depot holds the MMAE still bound to the
    # circulating conjugate; dose it with the MMAE-equivalent mass, not
    # with the milligram dose of conjugate.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ----- Observation -----
    # Cc is the unconjugated MMAE plasma concentration in ug/mL;
    # multiply by 1000 for the ng/mL scale of Babel 2026 Figure 5.
    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
