Oosten_2016_fentanyl <- function() {
  description <- "One-compartment population PK model for fentanyl administered by continuous subcutaneous infusion and transdermal matrix patch in adult cancer patients, with separate first-order absorption for each route, transdermal lag time, allometric body-weight scaling on CL/F and V/F (V/F fixed at 280 L), IIV on Ka (sc and td), F (td), and CL/F, IOV on transdermal Ka multiplexed by occasion, and proportional residual error (Oosten 2016)."
  reference <- paste(
    "Oosten AW, Abrantes JA, Jonsson S, de Bruijn P, Kuip EJM, Falcao A,",
    "van der Rijt CCD, Mathijssen RHJ.",
    "Treatment with subcutaneous and transdermal fentanyl: results from a",
    "population pharmacokinetic study in cancer patients.",
    "Eur J Clin Pharmacol. 2016 Apr;72(4):459-467.",
    "doi:10.1007/s00228-015-2005-x.",
    sep = " "
  )
  vignette <- "Oosten_2016_fentanyl"
  units <- list(time = "hr", dosing = "ug", concentration = "ng/mL") # dose ug + Vc L -> Cc ug/L = ng/mL; matches Oosten 2016 plasma units (Methods + Table 2)

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling per Oosten 2016 Methods + Results: CL/F = theta * (WT/70)^0.75 and V/F = theta * (WT/70)^1; reference weight 70 kg implied by Table 2 'CL 70kg/F' and 'V 70kg/F' headings. The paper notes 'allometrically scaled body weight on CL/F and V/F was found to explain some variability and was kept to increase model stability' (Results, Fentanyl pharmacokinetics) but does not numerically list the exponents; the conventional theoretical-allometric values 0.75 / 1.0 are inferred and encoded as fixed.",
      source_name        = "WT"
    ),
    OCC = list(
      description        = "Integer-valued transdermal-occasion indicator for IOV multiplexing on transdermal Ka.",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Oosten 2016 Patients/Materials/Methods defines an occasion as 'a transdermal dose followed by at least one observation'; IOV on td Ka was retained in the final model (Table 2: 32.8% CV). Values 1..10 identify successive transdermal occasions within subject; up to ten patch occasions are multiplexed via binary indicators oc1..oc10 inside model(). Non-transdermal records or transdermal occasions outside 1..10 set every indicator to 0 and yield the typical-value ka_td (no IOV applied). Extending beyond ten occasions requires adding etaiov_ka_td_<n> blocks; the cap reflects typical hospital-stay simulations (patches every 72 h).",
      source_name        = "OCC"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 52L,
    n_studies       = 1L,
    age_range       = "23-80 years",
    age_median      = "63 years",
    weight_range    = "Not tabulated as kg; BMI 18-40 kg/m^2 (median 25)",
    weight_median   = "Not directly reported; computed back from BMI median 25 kg/m^2 yields ~75 kg at a typical height of 1.73 m -- treated as approximate.",
    sex_female_pct  = 37,
    race_ethnicity  = c(Caucasian = 90, Other = 2, Unknown = 8),
    disease_state   = "Adults admitted to the Erasmus MC Cancer Institute (Rotterdam, The Netherlands) for moderate-to-severe cancer-related nociceptive pain receiving subcutaneous and/or transdermal fentanyl for pain control. Primary tumor sites included breast (15%), urinary tract incl. kidney (15%), prostate (13%), soft-tissue sarcoma/GIST (12%), colorectal (10%), and other (35%) (Oosten 2016 Table 1). No patient met the combined liver-failure criterion (AST or ALT > ULN, bilirubin > ULN, and albumin < LLN). WHO performance status 1 in 37%, 2 in 33%, 3 in 8%, unknown in 23%. Median NRS pain-at-rest 5 (range 2-10).",
    dose_range      = "Subcutaneous: continuous infusion 10-300 ug/h (median 75); transdermal: matrix patch (Fentanyl Sandoz Matrix) 12-400 ug/h (median 50), replaced every 72 h. Rotations from sc to td used a 1:1 dose-conversion ratio with the sc infusion continued at the same dose for 6 h after patch application then tapered 50% for another 6 h.",
    regions         = "The Netherlands (single tertiary cancer centre, January 2010 - November 2013).",
    n_observations  = 942L,
    sampling        = "Sparse opportunistic sampling: 942 fentanyl plasma samples from 52 patients (median 15 per patient, range 1-86) over up to 72 h after each change in opioid regimen; protocol prescribed twice-daily samples (around 8 am and 8 pm), a baseline sample before every regimen change, and a series around any extra sc bolus (baseline, 5, 15, 30, 60 min). Median observed concentration 1.33 ng/mL (range 0.122-10.7 ng/mL). 32 patients had semi-simultaneous sc and td exposure.",
    co_medication   = "One patient used the strong CYP3A4 inducer carbamazepine 200 mg during the study period; all others were screened to be free of strong CYP3A4 inhibitors or inducers.",
    notes           = "Three patients participated in the study twice. Covariate analysis beyond a priori allometric weight scaling was not pursued due to limited sample size (Discussion). Trial registration NTR4369 (Dutch Trial Register)."
  )

  ini({
    # Structural PK parameters reported at the reference 70 kg adult subject
    # (Oosten 2016 Table 2). Bootstrap mean and 95% CI given alongside the
    # NONMEM estimate.
    lka_sc     <- log(0.0358);     label("Subcutaneous absorption rate constant (Ka_sc, 1/h)")              # Oosten 2016 Table 2: ka_sc = 0.0358 (RSE 24.4%); bootstrap mean 0.0374 (95% CI 0.0248-0.0555)
    lka_td     <- log(0.0135);     label("Transdermal absorption rate constant (Ka_td, 1/h)")               # Oosten 2016 Table 2: ka_td = 0.0135 (RSE 16.8%); bootstrap mean 0.0140 (95% CI 0.0105-0.0188)
    llag_td    <- log(4.73);       label("Transdermal absorption lag time (Tlag_td, h)")                    # Oosten 2016 Table 2: t_lag_td = 4.73 (RSE 21.2%); bootstrap mean 4.65 (95% CI 2.25-6.98)
    lcl        <- log(49.6);       label("Apparent clearance CL/F at 70 kg (L/h)")                          # Oosten 2016 Table 2: CL/F at 70 kg = 49.6 (RSE 9.36%); bootstrap mean 50.4 (95% CI 40.9-61.6)
    lvc        <- fixed(log(280)); label("Apparent volume of distribution V/F at 70 kg (L)")                # Oosten 2016 Table 2 + Results: V70kg/F fixed to 280 L (citation [25] of Oosten 2016); sensitivity analysis with V/F +/-50% showed insensitivity of other parameters
    lfdepot_td <- fixed(log(1));   label("Transdermal bioavailability (typical fraction; F_sc = 1 anchor)") # Oosten 2016 Table 2: the typical F_td is not separately reported (only the IIV on F_td is listed); per the standard popPK convention for a reference route, F_td is anchored at 1 at the population level with between-subject variability allowed (etalfdepot_td below). F for the subcutaneous route is fixed at 1 implicitly by leaving f(depot) unset.

    # Allometric body-weight scaling. The exponents are encoded fixed at the
    # theoretical allometric values (0.75 on CL, 1.0 on V); the paper's narrative
    # records "allometrically scaled body weight on CL/F and V/F" without listing
    # numeric exponents.
    e_wt_cl    <- fixed(0.75);     label("Body-weight allometric exponent on CL/F (unitless)")              # Oosten 2016 Methods + Results: allometric scaling on CL/F; exponent inferred at the theoretical 0.75 because the paper does not list a custom value and the "70 kg" normalization in Table 2 implies the canonical form
    e_wt_vc    <- fixed(1);        label("Body-weight allometric exponent on V/F (unitless)")               # Oosten 2016 Methods + Results: allometric scaling on V/F; exponent inferred at the theoretical 1.0 (linear weight on volume) under the same reasoning as e_wt_cl

    # Inter-individual variability. Oosten 2016 Table 2 reports IIV as CV%; the
    # log-normal variance is omega^2 = log(1 + CV^2).
    etalka_sc     ~ log(1 + 0.935^2)  # Oosten 2016 Table 2: IIV Ka_sc 93.5% CV (RSE 15.2%); bootstrap 91.1% (95% CI 59.6-119)
    etalka_td     ~ log(1 + 0.424^2)  # Oosten 2016 Table 2: IIV Ka_td 42.4% CV (RSE 23.9%); bootstrap 41.4% (95% CI 10.5-59.2)
    etalfdepot_td ~ log(1 + 0.423^2)  # Oosten 2016 Table 2: IIV F_td  42.3% CV (RSE 30.0%); bootstrap 45.7% (95% CI 19.7-67.8)
    etalcl        ~ log(1 + 0.432^2)  # Oosten 2016 Table 2: IIV CL/F  43.2% CV (RSE 15.2%); bootstrap 41.6% (95% CI 27.1-53.9)

    # Inter-occasion variability on transdermal Ka. Oosten 2016 Methods defines
    # an occasion as a transdermal dose followed by at least one observation;
    # Table 2 reports IOV Ka_td 32.8% CV (RSE 51.1%); bootstrap 39.2% (95% CI
    # 12.0-77.0). A single variance is estimated and shared across occasions via
    # the NONMEM $OMEGA BLOCK(1) SAME idiom; nlmixr2 has no SAME shortcut so
    # occasion 1 carries the estimated variance and occasions 2..10 are fix()-ed
    # to the same value.
    etaiov_ka_td_1  ~ log(1 + 0.328^2)
    etaiov_ka_td_2  ~ fixed(log(1 + 0.328^2))
    etaiov_ka_td_3  ~ fixed(log(1 + 0.328^2))
    etaiov_ka_td_4  ~ fixed(log(1 + 0.328^2))
    etaiov_ka_td_5  ~ fixed(log(1 + 0.328^2))
    etaiov_ka_td_6  ~ fixed(log(1 + 0.328^2))
    etaiov_ka_td_7  ~ fixed(log(1 + 0.328^2))
    etaiov_ka_td_8  ~ fixed(log(1 + 0.328^2))
    etaiov_ka_td_9  ~ fixed(log(1 + 0.328^2))
    etaiov_ka_td_10 ~ fixed(log(1 + 0.328^2))

    # Residual error. Oosten 2016 Methods describes the residual error as
    # "additive on the log-scale", which is equivalent to a proportional error
    # model in nlmixr2's linear-concentration space for small-to-moderate CV.
    # Table 2 reports a proportional residual of 23.4% CV (RSE 5.17%).
    propSd <- 0.234; label("Proportional residual error (fraction)") # Oosten 2016 Table 2: proportional residual 23.4% CV; bootstrap 23.2% (95% CI 20.6-25.6)
  })

  model({
    # Decompose the integer-valued OCC column into binary occasion indicators
    # for IOV multiplexing on transdermal Ka. OCC = 1..10 selects the matching
    # per-occasion eta; OCC = 0 or any value outside 1..10 zeros every
    # indicator and yields the typical-value ka_td (no IOV applied). This is
    # appropriate for non-transdermal records (subcutaneous infusion or
    # plasma observations between patches) and for simulations beyond ten
    # transdermal occasions.
    oc1  <- (OCC == 1)
    oc2  <- (OCC == 2)
    oc3  <- (OCC == 3)
    oc4  <- (OCC == 4)
    oc5  <- (OCC == 5)
    oc6  <- (OCC == 6)
    oc7  <- (OCC == 7)
    oc8  <- (OCC == 8)
    oc9  <- (OCC == 9)
    oc10 <- (OCC == 10)
    iov_ka_td <- oc1 * etaiov_ka_td_1 + oc2 * etaiov_ka_td_2 + oc3 * etaiov_ka_td_3 +
                 oc4 * etaiov_ka_td_4 + oc5 * etaiov_ka_td_5 + oc6 * etaiov_ka_td_6 +
                 oc7 * etaiov_ka_td_7 + oc8 * etaiov_ka_td_8 + oc9 * etaiov_ka_td_9 +
                 oc10 * etaiov_ka_td_10

    # Individual PK parameters with allometric body-weight scaling on CL/F and
    # V/F (Oosten 2016 Methods). The transdermal Ka carries both IIV and the
    # per-occasion IOV; the transdermal bioavailability carries IIV around the
    # population anchor of 1; V/F is fixed at 280 L per 70 kg and carries no
    # IIV.
    ka_sc      <- exp(lka_sc + etalka_sc)
    ka_td      <- exp(lka_td + etalka_td + iov_ka_td)
    tlag_td    <- exp(llag_td)
    fdepot_td  <- exp(lfdepot_td + etalfdepot_td)
    cl         <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc         <- exp(lvc)          * (WT / 70)^e_wt_vc
    kel        <- cl / vc

    # Parallel first-order absorption from two depots into a single central
    # compartment. cmt column on each dose row selects the route:
    #   depot  (cmt = "depot"  / cmt 1): subcutaneous infusion. f(depot)
    #                                    unset -> F_sc = 1 (reference); no lag.
    #   depot2 (cmt = "depot2" / cmt 2): transdermal matrix patch loaded as a
    #                                    bolus (amt = patch_rate_ug_per_h *
    #                                    wear_time_h). f(depot2) = fdepot_td
    #                                    with IIV; lag(depot2) = tlag_td.
    d/dt(depot)   <- -ka_sc * depot
    d/dt(depot2)  <- -ka_td * depot2
    d/dt(central) <-  ka_sc * depot + ka_td * depot2 - kel * central

    lag(depot2)   <- tlag_td
    f(depot2)     <- fdepot_td

    # Plasma concentration. Dose units ug, Vc units L -> concentration ug/L
    # (= ng/mL), matching the Oosten 2016 plasma units.
    Cc <- central / vc

    # Proportional residual error on the linear concentration scale. The
    # "additive on the log-scale" parameterization in the paper maps to a
    # proportional error model in nlmixr2 for small-to-moderate CV.
    Cc ~ prop(propSd)
  })
}
