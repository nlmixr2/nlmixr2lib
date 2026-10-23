Han_2024_cabotegravir <- function() {
  description <- "Two-compartment population PK model for cabotegravir after oral tablets and long-acting (LA) intramuscular injections into either the gluteal or the lateral thigh (vastus lateralis) muscle, in adults with and without HIV-1 (Han 2024, phase 1 study 208832 and the ATLAS-2M thigh PK substudy). Three parallel first-order depots (oral, IM gluteal, IM thigh) feed one central compartment, so a gluteal depot and a thigh depot can absorb at the same time when the injection site is switched. Only the thigh absorption rate constant and the thigh bioavailability were estimated; every other parameter is fixed to the Han 2022 oral + gluteal model. The typical thigh absorption rate constant is the typical gluteal value plus 0.0002527 1/h, so it inherits the gluteal effects of sex, BMI, split injection and needle length; thigh bioavailability is 89.9% of gluteal."
  reference <- "Han K, D'Amico RD, Spreen WR, Ford SL. Population pharmacokinetics of cabotegravir following intramuscular thigh injections in adults with and without HIV. Antimicrob Agents Chemother. 2024;68(12):e00880-24. doi:10.1128/aac.00880-24. Fixed parameters from: Han K, Baker M, Lovern M, et al. Br J Clin Pharmacol. 2022;88(10):4607-4622. doi:10.1111/bcp.15439"
  vignette <- "Han_2024_cabotegravir"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Fixed from Han 2022: power functions centred at the Han 2022 model-building median of 76.6 kg (not the Han 2024 median of 75.75 kg), exponent 0.618 shared by CL/F and Q/F and 0.702 shared by V2/F and V3/F. Supply a per-subject constant. Han 2024 cohort 50.1-129.7 kg (Table 1).",
      source_name = "BWT"
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Fixed from Han 2022: the typical gluteal absorption rate constant is 50.9% lower in women (Han 2024 Table 2 'Female vs male on KA gluteal -50.9% (fixed)'). It reaches the thigh rate through TV_KA3 = TV_KA2 + THETA1, so KA_thigh is 34.5% lower in women at the typical values (Han 2024 Results). 37.88% of the Han 2024 cohort was female, assigned sex at birth (Table 1).",
      source_name = "female"
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Fixed from Han 2022: power function on the gluteal absorption rate constant centred at 25.4 kg/m^2 with exponent -0.766 (Han 2024 Table 2 'BMI exponent on KA gluteal -0.766 (fixed)' and the Results equation '(BMI/25.4)^-0.766'). Reaches the thigh rate through the typical gluteal value. Han 2024 cohort 17.87-50.9 kg/m^2 (Table 1).",
      source_name = "BMI"
    ),
    SMOKE = list(
      description = "Current smoker, 1 = current smoker, 0 = never or former smoker",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not a current smoker)",
      notes = "Fixed from Han 2022: CL/F is 17.4% higher in current smokers. Not printed in Han 2024 Table 2, but Table 2 footnote a fixes every parameter not related to the thigh injection to the Han 2022 final model, and the Han 2024 Discussion and Figure 5 simulate smoking status as a retained covariate. 27.3% of the Han 2024 cohort were current smokers (Table 1).",
      source_name = "Smoke"
    ),
    INJ_SPLIT = list(
      description = "Long-acting dose given as two split intramuscular injections, 1 = split, 0 = single (unsplit) injection",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (unsplit injection)",
      notes = "Fixed from Han 2022: the gluteal absorption rate constant is 47.8% higher for a split dose (Han 2024 Table 2 'Split injection on KA gluteal 47.8% (fixed)'). Because KA_thigh is built on the typical gluteal value, the same factor enters the thigh rate. Han 2024 Discussion: the dataset 'only included non-split injections', so the effect on a thigh dose was not estimated; leave INJ_SPLIT = 0 for thigh-only simulations. A single column serves both IM depots, so when a gluteal and a thigh depot absorb simultaneously they share the value on the current record.",
      source_name = "split"
    ),
    NEEDLE_LEN_IN = list(
      description = "Length of the needle used for the long-acting intramuscular injection",
      units = "inch",
      type = "continuous",
      reference_category = NULL,
      notes = "Fixed from Han 2022: power function on the gluteal absorption rate constant centred at 1.5 inch with exponent 0.478 (Han 2024 Table 2 'NDL exponent on KA gluteal 0.478 (fixed)' and the Results equation '(NDL/1.5)^0.478'). Supply INCHES. Reaches the thigh rate through the typical gluteal value; only 3 of 366 thigh injections (0.8%) used a needle other than 1.5 inch (Han 2024 Discussion), so the thigh-specific effect was not evaluated. The Figure 5 simulations assigned 1.5 inch to every injection.",
      source_name = "NDL"
    )
  )

  covariatesDataExcluded <- list(
    HIV_POS = list(
      description = "HIV-1 infection status, 1 = person living with HIV-1",
      units = "(binary)",
      type = "binary",
      notes = "Screened on the thigh parameters and not retained: 'KA thigh was not associated with HIV status' and the thigh relative bioavailability 'was not associated with any covariate evaluated' (Han 2024 Results). 89.39% of the cohort had HIV-1 (Table 1).",
      source_name = "HIV status"
    )
  )

  compartmentData <- list(
    depot_oral = list(analyte = "cabotegravir", units = "mg", specimen = "administration site", verified = TRUE),
    depot_im_gluteal = list(analyte = "cabotegravir", units = "mg", specimen = "administration site", verified = TRUE),
    depot_im_thigh = list(analyte = "cabotegravir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "cabotegravir", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cabotegravir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 132L,
    n_studies = 2L,
    age_range = "20-67 years",
    age_median = "42.5 years",
    weight_range = "50.1-129.7 kg",
    weight_median = "75.75 kg",
    bmi_range = "17.87-50.9 kg/m^2",
    bmi_median = "25.35 kg/m^2",
    sex_female_pct = 37.88,
    race_ethnicity = c(White = 78.79, Black = 18.18, Asian = 1.52, Other = 1.52),
    disease_state = "Adults living with HIV-1 who were virologically suppressed on cabotegravir + rilpivirine LA gluteal injections for at least 3 years (ATLAS-2M thigh PK substudy, n = 118, 89.39%) and HIV-1-negative adults (phase 1 study 208832, n = 14).",
    dose_range = "Study 208832: 4 weeks of oral cabotegravir 30 mg daily, a 10-14 day washout, then one 600 mg (3 mL) IM thigh injection. ATLAS-2M: last IM gluteal injection at Week -4, then four 400 mg (2 mL) monthly thigh injections (n = 64) or two 600 mg (3 mL) thigh injections 8 weeks apart (n = 54), then back to gluteal injections.",
    regions = "Multinational (ViiV Healthcare / GSK phase 1 and phase 3b studies).",
    smoking = "Never 56.8%, former 15.9%, current 27.3%",
    n_concentrations = "3481 plasma concentrations: 1249 after 366 thigh injections, 1998 after 1618 gluteal injections and 234 after oral dosing; BQL after LA injection handled with the M3 method (LLOQ 0.025 ug/mL)",
    notes = "Han 2024 Table 1 (baseline characteristics) and Table 3 (study designs). NONMEM 7.3, FOCE with interaction. Gluteal, oral and disposition parameters were fixed to the Han 2022 final model, which was fitted to 1647 adults in 16 studies."
  )

  ini({
    # Han 2024 estimated only the two thigh parameters and their IIV (Table 2).
    # Every other parameter was fixed to the Han 2022 final oral + gluteal
    # model (Han 2024 Table 2 footnote a and Methods), so it is wrapped in
    # fixed(). Han 2024 Table 2 reprints most of those values; KAoral, the
    # body-weight exponents, the smoking effect, the gluteal / oral / CL / V2
    # IIV and the residual error are read from Han 2022 Table 3 ('Final model
    # Estimate' column).
    #
    # IIV SCALE. Both papers print IIV as a percentage; it is read as a
    # log-normal CV, omega^2 = log(1 + CV^2), the same convention used for
    # Han_2022_cabotegravir.

    # --- Thigh injection (estimated in Han 2024) -------------------------
    dka_thigh <- 0.0002527; label("Additive increment of the typical thigh over the typical gluteal absorption rate constant, KA thigh - KA gluteal (1/h)") # Han 2024 Table 2 'KA thigh - KA gluteal (h-1) 0.0002527', 90% CI 0.0001559-0.0003573, RSE 20.7%
    lfdepot_thigh <- log(0.899); label("Thigh IM injection bioavailability relative to gluteal IM injection, F thigh/gluteal (fraction)") # Han 2024 Table 2 'F thigh/gluteal 89.9%', 90% CI 85.7%-94.7%, RSE 2.9%

    # --- Gluteal injection and oral tablet (fixed from Han 2022) ---------
    lka_gluteal <- fixed(log(0.000733)); label("Gluteal IM absorption rate constant KA gluteal for a 25.4 kg/m^2 male, unsplit dose, 1.5-inch needle (1/h)") # Han 2024 Table 2 'KA gluteal (h-1) 0.000733 (fixed)'; the Results equation rounds it to 0.00073
    lka_oral <- fixed(log(1.41)); label("Oral tablet absorption rate constant KAoral (1/h)") # Han 2022 Table 3 'KAoral (h-1) 1.41' (final model)
    lfdepot_oral <- fixed(log(0.756)); label("Oral tablet bioavailability relative to gluteal IM injection, F oral/gluteal (fraction)") # Han 2024 Table 2 'F oral/gluteal 75.6% (fixed)'

    # --- Disposition, typical subject 76.6 kg (fixed from Han 2022) ------
    lcl <- fixed(log(0.151)); label("Apparent clearance CL/F at 76.6 kg (L/h)") # Han 2024 Table 2 'CL/F (L/h) 0.151 (fixed)'
    lvc <- fixed(log(5.27)); label("Apparent central volume V2/F at 76.6 kg (L)") # Han 2024 Table 2 'V2/F (L) 5.27 (fixed)'
    lq <- fixed(log(0.507)); label("Apparent intercompartmental clearance Q/F at 76.6 kg (L/h)") # Han 2024 Table 2 'Q/F (L/h) 0.507 (fixed)'
    lvp <- fixed(log(2.43)); label("Apparent peripheral volume V3/F at 76.6 kg (L)") # Han 2024 Table 2 'V3/F (L) 2.43 (fixed)'

    # --- Covariate effects (fixed from Han 2022) -------------------------
    e_wt_cl_q <- fixed(0.618); label("Power exponent of baseline weight on CL/F and Q/F (unitless)") # Han 2022 Table 3 'BWT on CL/F and Q/F 0.618'
    e_wt_vc_vp <- fixed(0.702); label("Power exponent of baseline weight on V2/F and V3/F (unitless)") # Han 2022 Table 3 'BWT on V2/F and V3/F 0.702'
    e_smoke_cl <- fixed(0.174); label("Fractional change in CL/F for current smokers (unitless)") # Han 2022 Table 3 'Smoke on CL/F 17.4%'
    e_sexf_ka_gluteal <- fixed(-0.509); label("Fractional change in KA gluteal for women (unitless)") # Han 2024 Table 2 'Female vs male on KA gluteal -50.9% (fixed)'
    e_bmi_ka_gluteal <- fixed(-0.766); label("Power exponent of BMI on KA gluteal, reference 25.4 kg/m^2 (unitless)") # Han 2024 Table 2 'BMI exponent on KA gluteal -0.766 (fixed)'
    e_inj_split_ka_gluteal <- fixed(0.478); label("Fractional change in KA gluteal for a split dose (unitless)") # Han 2024 Table 2 'Split injection on KA gluteal 47.8% (fixed)'
    e_needle_len_in_ka_gluteal <- fixed(0.478); label("Power exponent of needle length on KA gluteal, reference 1.5 inch (unitless)") # Han 2024 Table 2 'NDL exponent on KA gluteal 0.478 (fixed)'

    # --- Between-subject variability --------------------------------------
    # Han 2024 Methods: the random effects for KA gluteal, F oral/gluteal,
    # KA thigh and F thigh/gluteal 'were estimated separately' (no
    # covariances). Han 2022 estimated covariances among KA gluteal, CL/F and
    # V2/F but did not publish them, so those etas are diagonal as well.
    etalka_thigh ~ 0.15961 # Han 2024 Table 2 'KA thigh IIV 41.6%': log(1 + 0.416^2) = 0.15961; shrinkage 24.6%
    etalfdepot_thigh ~ 0.051112 # Han 2024 Table 2 'F thigh/gluteal IIV 22.9%': log(1 + 0.229^2) = 0.051112; shrinkage 27.5%
    etalka_gluteal ~ fixed(0.28911) # Han 2022 Table 3 'KALA IIV 57.9%': log(1 + 0.579^2) = 0.28911
    etalka_oral ~ fixed(0.58736) # Han 2022 Table 3 'KAoral IIV 89.4%': log(1 + 0.894^2) = 0.58736
    etalfdepot_oral ~ fixed(0.029827) # Han 2022 Table 3 'F1 IIV 17.4%': log(1 + 0.174^2) = 0.029827
    etalcl ~ fixed(0.052867) # Han 2022 Table 3 'CL/F IIV 23.3%': log(1 + 0.233^2) = 0.052867
    etalvc ~ fixed(0.040383) # Han 2022 Table 3 'V2/F IIV 20.3%': log(1 + 0.203^2) = 0.040383

    # --- Residual error (fixed from Han 2022) ------------------------------
    addSd <- fixed(0.0319); label("Additive residual error (ug/mL)") # Han 2022 Table 3 'Add err (ug/mL) 0.0319'
    propSd <- fixed(0.273); label("Proportional residual error (fraction)") # Han 2022 Table 3 'Prop err 27.3%'
  })

  model({
    # Typical gluteal absorption rate constant with the Han 2022 covariate
    # effects (Han 2024 Results equation; power functions for continuous and
    # (1 + theta * IND) for categorical covariates, Han 2024 Methods).
    tv_ka_gluteal <- exp(lka_gluteal) *
      (1 + e_sexf_ka_gluteal * SEXF) *
      (1 + e_inj_split_ka_gluteal * INJ_SPLIT) *
      (BMI / 25.4)^e_bmi_ka_gluteal *
      (NEEDLE_LEN_IN / 1.5)^e_needle_len_in_ka_gluteal
    ka_gluteal <- tv_ka_gluteal * exp(etalka_gluteal)
    # Supplemental Table, final model: TV_KA3 = TV_KA2 + THETA1 and
    # KA3 = TV_KA3 * EXP(ETA). The increment is added to the TYPICAL gluteal
    # value, so the gluteal eta does not enter the thigh rate.
    ka_thigh <- (tv_ka_gluteal + dka_thigh) * exp(etalka_thigh)
    ka_oral <- exp(lka_oral + etalka_oral)

    cl <- exp(lcl + etalcl) * (WT / 76.6)^e_wt_cl_q * (1 + e_smoke_cl * SMOKE)
    vc <- exp(lvc + etalvc) * (WT / 76.6)^e_wt_vc_vp
    q <- exp(lq) * (WT / 76.6)^e_wt_cl_q
    vp <- exp(lvp) * (WT / 76.6)^e_wt_vc_vp

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot_oral) <- -ka_oral * depot_oral
    d/dt(depot_im_gluteal) <- -ka_gluteal * depot_im_gluteal
    d/dt(depot_im_thigh) <- -ka_thigh * depot_im_thigh
    d/dt(central) <- ka_oral * depot_oral + ka_gluteal * depot_im_gluteal + ka_thigh * depot_im_thigh -
      kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # The gluteal injection is the bioavailability reference of the apparent
    # parameters; oral and thigh doses carry relative bioavailabilities.
    f(depot_oral) <- exp(lfdepot_oral + etalfdepot_oral)
    f(depot_im_thigh) <- exp(lfdepot_thigh + etalfdepot_thigh)

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
