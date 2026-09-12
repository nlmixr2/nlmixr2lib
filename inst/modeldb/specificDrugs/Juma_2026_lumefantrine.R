# Population pharmacokinetic model for oral lumefantrine in pregnant and
# non-pregnant Kenyan women with uncomplicated Plasmodium falciparum malaria
# (Juma 2026, Br J Clin Pharmacol 92(3):911-921; doi:10.1111/bcp.70318).

Juma_2026_lumefantrine <- function() {
  description <- paste(
    "Population PK model for oral lumefantrine in 50 pregnant (second and",
    "third trimester) and 25 non-pregnant Kenyan women with uncomplicated",
    "Plasmodium falciparum malaria treated with the standard six-dose",
    "artemether-lumefantrine regimen (Juma 2026 Br J Clin Pharmacol).",
    "Absorption is a three-transit chain (ktr = 3 / MTT) feeding a separately",
    "estimated first-order absorption step (ka) into a two-compartment",
    "disposition model with first-order elimination. Relative bioavailability",
    "F is fixed at 1 for the population and carries a Box-Cox transformed",
    "inter-occasion random effect (lambda = -0.690) plus an a priori",
    "dose-saturable absorption term with Dose50 = 3.84 mg/kg. Allometric body",
    "weight scaling on all clearance (exponent 0.75) and volume (exponent 1)",
    "parameters centred at 70 kg. Pregnancy is a proportional covariate that",
    "increases CL/F by 23.2% and Vc/F by 28.1%, giving roughly 30% lower",
    "lumefantrine exposure in pregnant women. Inter-individual variability on",
    "MTT (CV 74.4%) and Vc/F (CV 14.1%); inter-occasion variability across two",
    "pooled dose occasions on ka (CV 109%) and F (CV 64.0%). Combined",
    "proportional (33.0%) and additive (12.0 ng/mL) residual error.",
    sep = " "
  )
  reference <- paste(
    "Juma E, Ding J, Ongas M, Koskei N, Onyango K, Oloo F, Aman R,",
    "Kokwaro G, Tarning J, Ogutu B (2026). Population pharmacokinetics of",
    "lumefantrine in pregnant and non-pregnant women with uncomplicated",
    "Plasmodium falciparum malaria in Western Kenya. British Journal of",
    "Clinical Pharmacology 92(3):911-921. doi:10.1111/bcp.70318.",
    sep = " "
  )
  vignette <- "Juma_2026_lumefantrine"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Juma 2026 Figure 1 (dose compartment
  # -> three transit-absorption compartments -> central -> peripheral) and
  # Methods 2.6 (venous plasma sampling).
  compartmentData <- list(
    depot       = list(analyte = "lumefantrine", units = "mg", specimen = "administration site", verified = TRUE),
    transit1    = list(analyte = "lumefantrine", units = "mg", specimen = "administration site", verified = TRUE),
    transit2    = list(analyte = "lumefantrine", units = "mg", specimen = "administration site", verified = TRUE),
    transit3    = list(analyte = "lumefantrine", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "lumefantrine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "lumefantrine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject at enrolment. Juma 2026 Methods 2.11:",
        "'Body weight was added as an allometric function on all clearance",
        "and volume parameters a priori, considering the strong biological",
        "prior of this relationship. The allometric function was centred on",
        "a typical body weight of 70 kg and scaled with an exponent of 0.75",
        "and 1 for clearance and volume parameters, respectively",
        "(Equations 2 and 3)'. The Table 2 note confirms the centring value:",
        "'Population estimates are given for a typical non-pregnant woman",
        "weighting 70 kg'. Note that the printed Equations 2 and 3 write the",
        "denominator symbolically as BW_median while the surrounding text and",
        "the Table 2 note both give 70 kg; the observed cohort medians are",
        "62.0 kg (pregnant) and 59.5 kg (non-pregnant), so 70 kg is a rounded",
        "reference rather than the study median (see vignette Errata).",
        "Table 1 body-weight ranges: pregnant 62.0 kg (40.0-86.6),",
        "non-pregnant 59.5 kg (45.0-89.0). Body weight also enters the",
        "dose-saturable bioavailability term through DOSE / WT.",
        sep = " "
      ),
      source_name        = "BW"
    ),
    PREG = list(
      description        = "Pregnancy status indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "1 = pregnant (second or third trimester), 0 = non-pregnant.",
        "Time-fixed per subject. Juma 2026 enrolled 50 pregnant women (30",
        "second trimester, 20 third trimester; gestational age median 26",
        "weeks, range 13-40) and 25 non-pregnant women (Table 1). Pregnancy",
        "is applied as a proportional covariate on apparent elimination",
        "clearance and apparent central volume, using the Table 2 note form",
        "theta = theta_TV * (1 + theta_cov) where theta_TV is the typical",
        "value and theta_cov the categorical pregnancy effect relative to a",
        "non-pregnant reference. Table 2: Pregnancy on CL = 23.2% and",
        "Pregnancy on Vc = 28.1%, i.e. theta_cov = 0.232 and 0.281.",
        "Trimester was evaluated as a separate categorical covariate and was",
        "NOT significant on any PK parameter, so a single pregnant-versus-",
        "non-pregnant indicator captures the whole published effect",
        "(Results 3.1.2). Reference category 0 = non-pregnant.",
        sep = " "
      ),
      source_name        = "PREG"
    ),
    DOSE = list(
      description        = "Per-dose lumefantrine amount administered (mg)",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Per-dose lumefantrine amount in milligrams, supplied as a",
        "per-dose-record covariate aligned with the corresponding event-table",
        "dosing row (use case (b) of the canonical DOSE entry, applied per",
        "dose record). Required by the model to compute the per-dose",
        "milligram-per-kilogram exposure DOSEMGKG = DOSE / WT, which drives",
        "the a priori dose-saturable bioavailability term of Juma 2026",
        "Equation 4, F = 100% * (1 - Dose / (Dose + Dose50)) with Dose50 =",
        "3.84 mg/kg fixed. Every participant in the study received four",
        "tablets of artemether/lumefantrine 20/120 mg (Coartem, Novartis) at",
        "hours 0, 8, 24, 36, 48 and 60 (Methods 2.3), so DOSE = 480 mg at",
        "each of the six dose records; a 70 kg woman therefore takes 6.86",
        "mg/kg per dose and the saturation term evaluates to 0.359. Set DOSE",
        "to the milligram amount administered at each dose event in the",
        "rxode2 event table, alongside amt (mg). Matches the usage in the",
        "sibling lumefantrine models Kloprogge_2018_lumefantrine.R (the",
        "pooled analysis this paper takes Dose50 from) and",
        "Chotsiri_2019_lumefantrine.R.",
        sep = " "
      ),
      source_name        = "DOSE"
    ),
    OCC = list(
      description        = "Dose-occasion index for the inter-occasion random effects",
      units              = "(count)",
      type               = "categorical",
      reference_category = 1,
      notes              = paste(
        "1 or 2. Juma 2026 Methods 2.10: 'Few PK samples were available from",
        "second to fifth dose, and dose occasions were therefore pooled",
        "together and evaluated as two separate dose occasions (i.e. first to",
        "third dose, and fourth to sixth dose)'. So OCC = 1 covers the doses",
        "at 0, 8 and 24 h and OCC = 2 covers the doses at 36, 48 and 60 h,",
        "and OCC must be carried as a time-varying column that switches from",
        "1 to 2 at 36 h. Decomposed inside model() into binary indicators oc1",
        "and oc2 that multiplex the inter-occasion etas on ka and on relative",
        "bioavailability, exactly as a NONMEM $PK block would",
        "(rxode2 parses but cannot simulate the eta ~ var | occ syntax).",
        "Any OCC value other than 1 or 2 sets both indicators to zero and so",
        "removes the inter-occasion random effects entirely; that is the",
        "convenient way to simulate an IOV-free typical profile.",
        sep = " "
      ),
      source_name        = "OCC"
    )
  )

  # Covariates the paper screened but did NOT retain in the final model.
  # Documentation only -- checkModelConventions() does not require these to
  # appear in model(), and they must not.
  covariatesDataExcluded <- list(
    GA = list(
      description = "Gestational age",
      units       = "weeks",
      type        = "continuous",
      notes       = paste(
        "Median 26 weeks (range 13-40) in the 50 pregnant participants",
        "(Table 1), determined from the date of the last monthly period",
        "(Methods 2.4). Evaluated as a continuous covariate on all PK",
        "parameters and WAS significant on clearance (p < .01) and on central",
        "volume (p < .05), but 'The model with pregnancy status (categorical",
        "covariate) had the lowest OFV and was chosen as the final model'",
        "(Results 3.1.2), and no gestational-age coefficient is reported in",
        "Table 2. The final model therefore carries the binary PREG effect",
        "only; the gestational-age parameterisation cannot be reconstructed",
        "from the published values.",
        sep = " "
      )
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Median 20.2 years (18.0-35.3) in pregnant and 25.2 years",
        "(18.1-35.6) in non-pregnant participants (Table 1). Screened on all",
        "PK parameters by forward selection / backward elimination and not",
        "retained: 'Other covariates, such as baseline parasitaemia and age,",
        "did not have significant impact on the PK parameters of",
        "lumefantrine' (Results 3.1.2).",
        sep = " "
      )
    ),
    BODYTEMP = list(
      description = "Baseline body temperature",
      units       = "degC",
      type        = "continuous",
      notes       = paste(
        "Median 36.8 degC (35.6-39.4) in pregnant and 37.6 degC (36.3-39.4)",
        "in non-pregnant participants (Table 1). Evaluated on all PK",
        "parameters in the stepwise covariate screen (Methods 2.11) and not",
        "retained in the final model.",
        sep = " "
      )
    ),
    PARA = list(
      description = "Plasmodium falciparum parasitaemia at enrolment",
      units       = "parasites/uL",
      type        = "continuous",
      notes       = paste(
        "Median 10 960/uL (1000-199 360) in pregnant and 24 320/uL",
        "(2560-152 960) in non-pregnant participants (Table 1). Screened and",
        "not retained. The Discussion notes this is inconsistent with the",
        "pooled analysis of reference 8 (Kloprogge 2018), in which relative",
        "bioavailability decreased with increasing enrolment parasitaemia,",
        "'This discrepancy may be attributed to the narrow distribution of",
        "baseline parasite density in our study'. The sibling model",
        "Kloprogge_2018_lumefantrine.R does encode that parasitaemia effect.",
        sep = " "
      )
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 75L,
    n_studies       = 1L,
    n_pregnant      = 50L,
    age_range       = "18.0-35.3 years (pregnant) and 18.1-35.6 years (non-pregnant) (Table 1)",
    age_median      = "20.2 years (pregnant) and 25.2 years (non-pregnant) (Table 1)",
    weight_range    = "40.0-86.6 kg (pregnant) and 45.0-89.0 kg (non-pregnant) (Table 1)",
    weight_median   = "62.0 kg (pregnant) and 59.5 kg (non-pregnant) (Table 1); model centred at 70 kg",
    sex_female_pct  = 100,
    ga_range        = "13-40 weeks, median 26 weeks (Table 1); 30 women in the second and 20 in the third trimester",
    disease_state   = paste(
      "Uncomplicated Plasmodium falciparum malaria. Eligibility required",
      "fever >= 37.5 degC or a history of fever in the preceding 24 h or",
      "other malaria symptoms, P. falciparum mono-infection with a parasite",
      "count of 1000-200 000 parasites/uL, haemoglobin >= 8 g/dL, and no",
      "antimalarial or antimalarial-active antimicrobial in the preceding",
      "month (Methods 2.2). The day-28 PCR-corrected cure rate was 100% in",
      "both groups; all three recurrences during follow-up were PCR-confirmed",
      "re-infections.",
      sep = " "
    ),
    dose_range      = paste(
      "Directly observed oral artemether/lumefantrine 20/120 mg (Coartem,",
      "Novartis AG), four tablets per dose (480 mg lumefantrine) at hours 0,",
      "8, 24, 36, 48 and 60, for a total of 2880 mg lumefantrine. Doses were",
      "given with water followed by 250 mL of milk or a light meal",
      "(Methods 2.3).",
      sep = " "
    ),
    regions         = "Kenya (Ahero County Hospital, Kisumu County, western Kenya; area of high perennial malaria transmission)",
    notes           = paste(
      "Enrolled August 2013 to April 2014; Pan African Clinical Trial",
      "Registry PACTR201211000451437. A total of 1151 venous plasma",
      "lumefantrine concentrations were available from all 75 participants",
      "with at least one post-dose sample. Rich sampling over the first 72 h",
      "was split between two sampling arms within each of the three groups",
      "(non-pregnant, second-trimester, third-trimester), plus single samples",
      "on days 4, 5, 6, 7 and 14 for everyone (Methods 2.6). Lumefantrine was",
      "quantified in venous plasma by LC-MS/MS with an LLOQ of 10 ng/mL;",
      "86 samples (7.4%) were below the LLOQ and were omitted (M1 approach),",
      "which the authors justified with categorical visual predictive checks.",
      "Model building used NONMEM 7.4 with FOCE-I; parameter uncertainty came",
      "from sampling importance resampling (SIR). Eta shrinkage was 16.7% on",
      "MTT, 60.4% on ka, 16.5% on F and 54.8% on Vc/F (Table 2).",
      sep = " "
    )
  )

  ini({
    # Structural population parameters are the NONMEM population estimates in
    # Juma 2026 Table 2, which are centred on a 'typical' non-pregnant woman
    # weighing 70 kg with acute P. falciparum malaria (Table 2 note). The
    # paper reports linear-scale values; log() is applied here for the
    # nlmixr2 internal scale.
    lmtt <- log(4.70)  ; label("Mean transit time of the three-transit absorption chain MTT (h)")      # Juma 2026 Table 2: MTT = 4.70 h (RSE 6.9%; SIR median 4.62, 95% CI 4.01-5.28)
    lka  <- log(0.411) ; label("First-order absorption rate constant ka from the last transit compartment into central (1/h)") # Juma 2026 Table 2: Ka = 0.411 1/h (RSE 12.7%; SIR median 0.401, 95% CI 0.322-0.533)
    lcl  <- log(4.03)  ; label("Apparent elimination clearance CL/F (L/h, typical non-pregnant 70-kg woman)")  # Juma 2026 Table 2: CL/F = 4.03 L/h (RSE 5.1%; SIR median 4.03, 95% CI 3.66-4.45)
    lvc  <- log(108)   ; label("Apparent central volume of distribution Vc/F (L, typical non-pregnant 70-kg woman)") # Juma 2026 Table 2: VC/F = 108 L (RSE 8.6%; SIR median 108, 95% CI 90.3-127)
    lq   <- log(1.67)  ; label("Apparent inter-compartmental clearance Q/F (L/h, typical non-pregnant 70-kg woman)")  # Juma 2026 Table 2: Q/F = 1.67 L/h (RSE 7.1%; SIR median 1.67, 95% CI 1.46-1.93)
    lvp  <- log(184)   ; label("Apparent peripheral volume of distribution Vp/F (L, typical non-pregnant 70-kg woman)") # Juma 2026 Table 2: Vp/F = 184 L (RSE 6.1%; SIR median 185, 95% CI 168-211)

    # Relative bioavailability anchored at 1. Methods 2.9: 'Relative
    # bioavailability (F) was set to unity for the population but included to
    # allow variability to be estimated between patients in the absorption of
    # lumefantrine'. Table 2 reports F (%) = 100 (fixed). All F variability is
    # carried by the Box-Cox-transformed inter-occasion etas below plus the
    # dose-saturable term in model().
    lfdepot <- fixed(log(1)) ; label("Relative oral bioavailability F (unitless)")  # Juma 2026 Table 2: F (%) = 100 (fixed)

    # Box-Cox shape parameter for the random effect on relative
    # bioavailability (Petersson 2009 form, as used by
    # GonzalezSales_2015_testosterone.R and Friberg_2012_voriconazole.R):
    # eta_bc = (exp(eta)^lambda - 1) / lambda, which is coded below in the
    # numerically equivalent form (exp(lambda * eta) - 1) / lambda. lambda = 0
    # recovers a plain log-normal random effect, and eta = 0 gives eta_bc = 0
    # for any lambda, so the typical value of F is unchanged. Results 3.1.2:
    # 'Applying a Box-Cox transformation of inter-individual variability of
    # relative bioavailability improved model fit further (delta AIC = -14.9)'.
    boxcox_fdepot <- -0.690 ; label("Box-Cox shape parameter lambda for the random effect on relative bioavailability F (unitless)")  # Juma 2026 Table 2: Box-cox on F = -0.690 (RSE 13.2%; SIR median -0.664, 95% CI -0.830 to -0.463)

    # Allometric exponents, fixed a priori at the canonical values stated in
    # Methods 2.11 (0.75 on clearance parameters, 1 on volume parameters) and
    # written out in Equations 2 and 3.
    e_wt_cl <- fixed(0.75) ; label("Allometric exponent on CL/F and Q/F (unitless)")   # Juma 2026 Methods 2.11 and Equation 2: (BW_i / BW_median)^0.75
    e_wt_vc <- fixed(1)    ; label("Allometric exponent on Vc/F and Vp/F (unitless)")  # Juma 2026 Methods 2.11 and Equation 3: (BW_i / BW_median)^1.0

    # Dose-saturable absorption on relative bioavailability, added a priori
    # and not estimated. Methods 2.11: 'The saturation parameter, dose50
    # (3.84 mg/kg), was not estimated but taken from a large pooled analysis
    # of lumefantrine and refers to the dosage at which a typical patient
    # reaches 50% saturation of the absorption'. The pooled analysis is
    # reference 8 = Kloprogge 2018, whose own Table 2 value is 3.86 mg/kg;
    # this model uses the 3.84 mg/kg value printed by Juma 2026.
    dose50 <- fixed(3.84) ; label("Per-dose amount at 50% saturation of absorption (mg/kg)")  # Juma 2026 Table 2: Dose50 (mg/kg) on F = 3.84 (fixed); Methods 2.11 and Equation 4

    # Pregnancy proportional covariate effects. Table 2 note: 'Pregnancy
    # status was implemented using a proportional covariate model
    # theta = theta_TV * (1 + theta_cov), where theta_TV is the typical value
    # of a given parameter, and theta_cov is the categorical pregnancy effect
    # compared to a non-pregnant reference population'. Results 3.1.2:
    # 'Pregnant women exhibited a 23.2% (95%CI: 10.8-34.8%) and 28.1%
    # (95%CI: 7.3-51.8%) higher CL and V, respectively'.
    e_preg_cl <- 0.232 ; label("Pregnancy effect on CL/F: CL_pregnant / CL_non-pregnant - 1 = +0.232")  # Juma 2026 Table 2: Pregnancy on CL (%) = 23.2 (RSE 25.8%; SIR median 23.3, 95% CI 10.8-34.8)
    e_preg_vc <- 0.281 ; label("Pregnancy effect on Vc/F: Vc_pregnant / Vc_non-pregnant - 1 = +0.281")  # Juma 2026 Table 2: Pregnancy on Vc (%) = 28.1 (RSE 38.8%; SIR median 27.0, 95% CI 7.3-51.8)

    # Random effects. Table 2 note: 'Coefficients of variation for inter-
    # individual variability (IIV) and inter-occasion variability (IOV),
    # presented in the table, were calculated as 100 * (e^variance - 1)^1/2',
    # so the internal variances are recovered by
    # omega^2 = log((CV/100)^2 + 1):
    #   MTT IIV 74.4% -> log(0.744^2 + 1) = 0.4405336
    #   Vc  IIV 14.1% -> log(0.141^2 + 1) = 0.0196860
    #   ka  IOV 109%  -> log(1.09^2  + 1) = 0.7830336
    #   F   IOV 64.0% -> log(0.640^2 + 1) = 0.3433060
    # Only these four terms appear in Table 2; in particular the final model
    # carries NO inter-individual variability on CL/F, Q/F or Vp/F, despite
    # the generic exponential-IIV statement of Equation 1 (see vignette
    # Errata).
    etalmtt ~ 0.4405336  # Juma 2026 Table 2: MTT (IIV) = 74.4% CV (RSE 9.5%; SIR median 75.5, 95% CI 63.1-89.5; shrinkage 16.7%)
    etalvc  ~ 0.0196860  # Juma 2026 Table 2: VC/F (IIV) = 14.1% CV (RSE 41.9%; SIR median 14.6, 95% CI 3.7-25.4; shrinkage 54.8%)

    # Inter-occasion variability across the two pooled dose occasions
    # (Methods 2.10). One eta per occasion, all occasions sharing a single
    # variance -- the NONMEM '$OMEGA BLOCK(1) SAME' idiom -- so the second
    # occasion's variance is fixed to the first. rxode2 parses but cannot
    # simulate the native 'eta ~ var | occ' multi-level syntax, so the
    # occasion-indicator expansion in model() is used instead.
    etaiov_ka_1 ~ 0.7830336       # Juma 2026 Table 2: Ka (IOV) = 109% CV (RSE 11.6%; SIR median 111, 95% CI 85.0-135; shrinkage 60.4%)
    etaiov_ka_2 ~ fix(0.7830336)  # same variance as occasion 1 by construction ('$OMEGA BLOCK(1) SAME')

    etaiov_fdepot_1 ~ 0.3433060       # Juma 2026 Table 2: F (IOV) = 64.0% CV (RSE 6.9%; SIR median 64.1, 95% CI 54.6-72.8; shrinkage 16.5%). Enters F through the Box-Cox transformation.
    etaiov_fdepot_2 ~ fix(0.3433060)  # same variance as occasion 1 by construction ('$OMEGA BLOCK(1) SAME')

    # Residual error. Results 3.1.2: 'A combination of a proportional and
    # additive error showed the best residual diagnostics and was implemented
    # throughout in the model development'. Table 2 reports the two terms on
    # the linear concentration scale: proportional 33.0% and additive
    # 12.0 ng/mL. The additive term is converted to the ug/mL concentration
    # unit used by this model (12.0 ng/mL = 0.0120 ug/mL).
    propSd <- 0.330  ; label("Proportional residual SD (unitless fraction of the prediction)")  # Juma 2026 Table 2: Proportional (%) = 33.0 (RSE 3.4%; SIR median 32.9, 95% CI 31.0-35.2)
    addSd  <- 0.0120 ; label("Additive residual SD (ug/mL)")                                     # Juma 2026 Table 2: Additive (ng/mL) = 12.0 (RSE 26.3%; SIR median 12.5, 95% CI 8.6-18.9) = 0.0120 ug/mL
  })

  model({
    # Occasion indicators for the two pooled dose occasions (Methods 2.10:
    # doses 1-3 and doses 4-6). This expansion is what a NONMEM $PK block
    # writes out, and it is required because rxode2 cannot simulate the
    # native multi-level IOV syntax.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_ka     <- oc1 * etaiov_ka_1     + oc2 * etaiov_ka_2
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2

    # Absorption chain rate constants. Juma 2026 Figure 1 draws the dose
    # compartment transferring at ktr into a block of n = 3 transit
    # compartments, which then transfers at ka into the central compartment,
    # and Results 3.1.2 states the transit and absorption rate constants were
    # 'separately estimated'. MTT therefore spans only the ktr-governed part
    # of the chain -- three first-order transfers, depot -> transit1 ->
    # transit2 -> transit3 -- so ktr = 3 / MTT, and the final transit3 ->
    # central step is governed by ka. This is the same convention as the
    # co-author's Ding_2024_amodiaquine.R (n = 2 transits, ktr = 2 / MTT,
    # separate ka on the last step).
    mtt <- exp(lmtt + etalmtt)
    ktr <- 3 / mtt
    ka  <- exp(lka + iov_ka)

    # Individual disposition parameters. Allometric body-weight scaling on
    # all clearance and volume parameters centred at 70 kg (Equations 2 and
    # 3), then the proportional pregnancy effect on CL/F and Vc/F (Table 2
    # note). Q/F and Vp/F carry allometry only -- pregnancy was eliminated
    # from them in the stepwise backward elimination (Results 3.1.2).
    cl <- exp(lcl)           * (WT / 70)^e_wt_cl * (1 + e_preg_cl * PREG)
    vc <- exp(lvc + etalvc)  * (WT / 70)^e_wt_vc * (1 + e_preg_vc * PREG)
    q  <- exp(lq)            * (WT / 70)^e_wt_cl
    vp <- exp(lvp)           * (WT / 70)^e_wt_vc

    # Two-compartment disposition micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODE system (Juma 2026 Figure 1). Amounts are mg, volumes are L.
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot    - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1 - ktr * transit2
    d/dt(transit3)    <-  ktr * transit2 - ka  * transit3
    d/dt(central)     <-  ka  * transit3 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                   k12 * central - k21 * peripheral1

    # Relative bioavailability applied to the dose compartment, composed of
    # two multiplicative factors:
    #   1. The population value F = 1 (fixed) carrying a Box-Cox transformed
    #      inter-occasion random effect. Petersson 2009 form
    #      eta_bc = (exp(eta)^lambda - 1) / lambda, written equivalently as
    #      (exp(lambda * eta) - 1) / lambda. At eta = 0 this is 0 for any
    #      lambda, so the typical F is unaffected.
    #   2. The a priori dose-saturable absorption term of Equation 4,
    #      F = 100% * (1 - Dose / (Dose + Dose50)) with the per-dose amount
    #      expressed in mg/kg as DOSE / WT and Dose50 = 3.84 mg/kg. For the
    #      study's 480 mg dose in a 70 kg woman this evaluates to 0.359.
    dose_mgkg     <- DOSE / WT
    fdose         <- 1 - dose_mgkg / (dose50 + dose_mgkg)
    eta_fdepot_bc <- (exp(boxcox_fdepot * iov_fdepot) - 1) / boxcox_fdepot
    f(depot)      <- exp(lfdepot + eta_fdepot_bc) * fdose

    # Lumefantrine venous plasma concentration. Dose units are mg and vc is
    # in L, so central / vc is mg/L = ug/mL. The paper reports Cmax in ug/mL
    # (Table 3) and day-7 concentrations in ng/mL; multiply Cc by 1000 to
    # compare against the 200 ng/mL day-7 efficacy threshold.
    Cc <- central / vc

    # Combined proportional and additive residual error on the linear
    # concentration scale (Table 2).
    Cc ~ prop(propSd) + add(addSd)
  })
}
