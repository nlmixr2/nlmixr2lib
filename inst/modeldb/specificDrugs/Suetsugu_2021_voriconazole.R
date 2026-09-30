Suetsugu_2021_voriconazole <- function() {
  description <- paste0(
    "Steady-state Michaelis-Menten population PK model relating the ",
    "voriconazole steady-state trough concentration to the daily ",
    "maintenance dose in Japanese adult allogeneic haematopoietic stem cell ",
    "transplant recipients (Suetsugu 2021; n = 47, 216 trough samples). ",
    "The model is algebraic -- Css_trough = Km * F * dose / (Vmax - F * ",
    "dose), with the daily dose supplied as the DOSE_VORI_MGD covariate and ",
    "F fixed to 1 -- because every observation was a steady-state trough and ",
    "no volume or absorption parameter was estimated. Vmax is 1.72-fold ",
    "higher with concomitant letermovir and 1.30-fold higher with ",
    "concomitant methylprednisolone (multiplicative), with log-normal IIV on ",
    "Vmax and an additive residual error. No steady state exists when F * ",
    "dose >= Vmax; the equation then returns a negative or infinite value."
  )
  reference <- paste(
    "Suetsugu K, Muraki S, Fukumoto J, Matsukane R, Mori Y, Hirota T,",
    "Miyamoto T, Egashira N, Akashi K, Ieiri I. Effects of Letermovir and/or",
    "Methylprednisolone Coadministration on Voriconazole Pharmacokinetics in",
    "Hematopoietic Stem Cell Transplantation: A Population Pharmacokinetic",
    "Study. Drugs R D. 2021;21(4):419-429. doi:10.1007/s40268-021-00365-0"
  )
  vignette <- "Suetsugu_2021_voriconazole"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    DOSE_VORI_MGD = list(
      description = "Patient's own total daily maintenance voriconazole dose",
      units = "mg/day",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "The dose-rate input of the steady-state Michaelis-Menten equation ",
        "(Suetsugu 2021 Eq. 1, 'Daily dose'). Sum of all maintenance doses ",
        "given in a day (e.g. 200 mg twice daily -> 400). Oral and ",
        "intravenous doses are treated identically because F is fixed to 1 ",
        "(Section 2.5). The model is only defined for DOSE_VORI_MGD < Vmax ",
        "(670 mg/day typical without interacting co-medication)."
      ),
      source_name = "Daily dose"
    ),
    CONMED_LETERMOVIR = list(
      description = "Concomitant letermovir administration indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant letermovir)",
      notes = paste0(
        "Suetsugu 2021 Eq. 4: 'LMVi is concomitant letermovir use in ",
        "subject i (used = 1, not used = 0)'. Multiplies Vmax by 1.72 when 1 ",
        "(Table 2 'Effect of LMV on Vmax'). Letermovir 480 mg/day was ",
        "generally started from the day of stem cell infusion (Section 2.1); ",
        "19 of 47 patients and 71 of 216 samples were on letermovir (Table 1, ",
        "Fig. 1). Time-varying at the sample level."
      ),
      source_name = "LMV"
    ),
    CONMED_METHYLPREDNISOLONE = list(
      description = "Concomitant methylprednisolone administration indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant methylprednisolone)",
      notes = paste0(
        "Suetsugu 2021 Eq. 4: 'mPSLi is concomitant methylprednisolone use ",
        "in subject i (used = 1, not used = 0)'. Multiplies Vmax by 1.30 ",
        "when 1 (Table 2 'Effect of mPSL on Vmax'). 18 patients received ",
        "methylprednisolone, median 40.0 mg/day (range 25.0-80.0) at the ",
        "first voriconazole measurement (Section 3.1). Prednisolone and ",
        "dexamethasone were screened separately and not retained."
      ),
      source_name = "mPSL"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Biological sex, 1 = female",
      units = "(binary)",
      type = "binary",
      notes = "Screened (Section 2.6); not retained. Table 1: 23 male, 24 female."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened (Section 2.6); not retained. Table 1 median 51 years (range 22-69)."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened (Section 2.6); not retained. Table 1 median 55.0 kg (range 31.3-90.3)."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened (Section 2.6); not retained. Table 1 median 3.8 g/dL (range 2.5-5.7), i.e. 38 g/L (25-57 g/L)."
    ),
    CRP = list(
      description = "C-reactive protein",
      units = "mg/L",
      type = "continuous",
      notes = "Screened (Section 2.6); not retained. Table 1 median 0.32 mg/dL (range 0.01-19.60), i.e. 3.2 mg/L (0.1-196 mg/L)."
    ),
    CONMED_PPI = list(
      description = "Concomitant proton-pump inhibitor use",
      units = "(binary)",
      type = "binary",
      notes = "Screened (Section 2.6: rabeprazole, esomeprazole, lansoprazole, omeprazole, vonoprazan); not retained. Table 1 counts 16, 13, 7, 6 and 5 patients respectively."
    ),
    CONMED_STEROID = list(
      description = "Concomitant systemic corticosteroid",
      units = "(binary)",
      type = "binary",
      notes = "Prednisolone (n = 19) and dexamethasone (n = 2) were screened (Sections 2.6 and 3.1) and not retained; of the three corticosteroids only methylprednisolone was retained (Discussion)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 47L,
    n_studies = 1L,
    n_observations = 216L,
    age_range = "22-69 years",
    age_median = "51 years",
    weight_range = "31.3-90.3 kg",
    weight_median = "55.0 kg",
    sex_female_pct = 51.1,
    race_ethnicity = "Japanese (single centre, Kyushu University Hospital, Fukuoka)",
    disease_state = paste0(
      "Adult allogeneic haematopoietic stem cell transplant recipients ",
      "(April 2016 - March 2020) on tacrolimus graft-versus-host disease ",
      "prophylaxis, receiving voriconazole for prophylaxis (n = 14) or ",
      "treatment (n = 33) of invasive fungal infection with therapeutic drug ",
      "monitoring. Diagnoses AML/MDS 23, ALL 12, lymphoma 6, other 6. ",
      "Excluded: voriconazole started > 100 days after HSCT, age < 20 years, ",
      "and liver dysfunction (AST/ALT > 5 x ULN or total bilirubin > 2.0 ",
      "mg/dL)."
    ),
    dose_range = paste0(
      "Voriconazole maintenance doses from routine care; 38 patients oral ",
      "only, 1 intravenous only, 8 both (28 of 216 concentrations, about 13%, ",
      "during intravenous administration). The paper does not tabulate the ",
      "observed daily-dose range."
    ),
    co_medication = paste0(
      "Letermovir 19 patients; methylprednisolone 18, prednisolone 19, ",
      "dexamethasone 2; PPIs rabeprazole 16, esomeprazole 13, lansoprazole ",
      "7, omeprazole 6, vonoprazan 5 (Table 1, Section 3.1)."
    ),
    regions = "Japan (Kyushu University Hospital, Fukuoka)",
    sampling_window = paste0(
      "Steady-state trough concentrations only; concentrations within 5 ",
      "days of voriconazole initiation were excluded (Section 2.2)."
    ),
    assay = paste0(
      "UPLC or UPLC-MS/MS outsourced to SRL (LLOQ 0.1 mg/L) or LSI Medience ",
      "(LLOQ 0.3 mg/L); CV < 15% (Section 2.2)."
    ),
    notes = paste0(
      "Retrospective single-centre TDM analysis; NONMEM 7.4.3 FOCE-I. ",
      "Evaluation by GOF, case-deletion diagnostics and a 1000-replicate ",
      "bootstrap (Table 2)."
    )
  )

  ini({
    # Final-model estimates - Suetsugu 2021 Table 2 'Parameter estimates
    # from final model and bootstrap evaluation'.
    lkm <- log(1.97); label("Michaelis-Menten constant Km (mg/L)") # Table 2: Km = 1.97 ug/mL (RSE 23.5%); 1 ug/mL = 1 mg/L
    lvmax <- log(670); label("Maximum elimination rate Vmax without letermovir or methylprednisolone (mg/day)") # Table 2: Vmax = 670 mg/day (RSE 10.2%); Eq. 4
    # Section 2.5: 'Because the bioavailability of voriconazole is known to be
    # essentially 100% [19], the F-value was set to 1.'
    lfdepot <- fixed(log(1)); label("Bioavailability F applied to the daily dose (fraction)") # Section 2.5: F set to 1

    # Categorical covariate effects on Vmax, Eq. 3 form P = Ppop * theta^X
    # (Eq. 4: Vmax_i = 670 x 1.72^LMV x 1.30^mPSL).
    e_conmed_letermovir_vmax <- 1.72; label("Multiplicative factor on Vmax with concomitant letermovir (unitless)") # Table 2: Effect of LMV on Vmax = 1.72 (RSE 2.3%)
    e_conmed_methylprednisolone_vmax <- 1.30; label("Multiplicative factor on Vmax with concomitant methylprednisolone (unitless)") # Table 2: Effect of mPSL on Vmax = 1.30 (RSE 7.5%)

    # IIV on Vmax, exponential model (Section 2.5). Table 2 reports 23.6 CV%;
    # omega^2 = log(1 + 0.236^2) = 0.05421.
    etalvmax ~ 0.05421 # Table 2: IIV Vmax 23.6 CV% (RSE 13.9%, shrinkage 11.2%)

    # Additive residual error (Section 2.5), in concentration units.
    addSd <- 0.77; label("Additive residual error SD (mg/L)") # Table 2: Additive error = 0.77 ug/mL (RSE 10.8%)
  })

  model({
    # Eq. 4: Vmax_i = 670 x 1.72^LMV_i x 1.30^mPSL_i, with exponential IIV.
    vmax <- exp(lvmax + etalvmax) *
      e_conmed_letermovir_vmax^CONMED_LETERMOVIR *
      e_conmed_methylprednisolone_vmax^CONMED_METHYLPREDNISOLONE
    km <- exp(lkm)
    fdepot <- exp(lfdepot)

    # Eq. 1: Css_trough = Km x F x Daily dose / (Vmax - F x Daily dose).
    # Daily dose (mg/day) and Vmax (mg/day) share units, so Cc carries the
    # units of Km (mg/L). The equation has no steady-state solution when
    # F x Daily dose >= Vmax (elimination saturates below the input rate); it
    # then returns a negative or infinite value, which is not a concentration.
    Cc <- km * fdepot * DOSE_VORI_MGD / (vmax - fdepot * DOSE_VORI_MGD)

    Cc ~ add(addSd)
  })
}
