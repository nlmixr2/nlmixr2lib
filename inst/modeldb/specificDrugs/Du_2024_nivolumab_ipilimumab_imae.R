Du_2024_nivolumab_ipilimumab_imae <- function() {
  description <- paste0(
    "Stratified Cox proportional-hazards exposure-SAFETY model for the ",
    "time to first GRADE 2 OR HIGHER IMMUNE-MEDIATED ADVERSE EVENT ",
    "(gr2+ IMAE; immune-related adverse events substituted for the ",
    "ipilimumab-monotherapy studies) with nivolumab, ipilimumab or ",
    "nivolumab plus ipilimumab, pooled across 15 adult and pediatric ",
    "studies (Du 2024, n = 3507, mostly melanoma; 42 patients < 12 years ",
    "and 55 adolescents 12 to < 18 years). The model returns the ",
    "RELATIVE hazard against the paper's Figure 1 reference patient (a ",
    "60-year-old 75 kg white male with first-line advanced melanoma, ",
    "ECOG 0, tumour PD-L1 < 5%, LDH at the upper limit of normal and no ",
    "drug exposure). Exposure enters as the TIME-VARYING daily average ",
    "serum concentrations of nivolumab and ipilimumab (both linear, in ",
    "ug/mL) plus their product: the ipilimumab term (HR 1.008 per ug/mL) ",
    "and the interaction (HR 0.9995 per (ug/mL)^2) are significant, the ",
    "nivolumab term is not. Baseline covariates: age, body weight, line ",
    "of therapy, tumour type / treatment setting (adjuvant melanoma, ",
    "other tumours), PD-L1 status, ECOG, sex, race and log(LDH/ULN). NO ",
    "BASELINE HAZARD IS ENCODED: a Cox regression is semiparametric, and ",
    "this one was additionally STRATIFIED by treatment (nivolumab ",
    "monotherapy, ipilimumab monotherapy, nivolumab plus ipilimumab), so ",
    "three separate h0(t) were left unspecified by the fit; the ",
    "covariate coefficients are shared across the strata, and only ",
    "RATIOS of hr within one treatment stratum are meaningful. No PK ",
    "layer and no ODE: the daily Cavg columns come from individual ",
    "empirical-Bayes PK predictions in the source analysis. No ",
    "between-subject random effect and no residual error are estimated."
  )
  reference <- paste(
    "Du S, Zhao Y, Hu Z, Liu S, Roy A, Shen J, Zhu L, Hamuro L.",
    "Pediatric model-based dose optimization using a pooled",
    "exposure-response safety analysis for nivolumab and nivolumab plus",
    "ipilimumab combination in melanoma.",
    "CPT Pharmacometrics Syst Pharmacol. 2024;13(1):168.",
    "doi:10.1002/psp4.13070.",
    "Coefficients are transcribed from Table 1; the model formula is",
    "Appendix S1 Table S4/S5 (bold row). The daily average exposures",
    "derive from the same group's adult and pediatric population PK",
    "analysis, published as Hu 2024 (doi:10.1002/psp4.13098; see",
    "modellib('Hu_2024_nivolumab') and modellib('Hu_2024_ipilimumab')).",
    sep = " "
  )
  vignette <- "Du_2024_nivolumab_ipilimumab"
  units <- list(
    time = "n/a (semiparametric stratified Cox regression; the baseline hazards and hence the time scale are left unspecified by the method)",
    dosing = "n/a (no dose events; exposure enters as the time-varying CAV_NIVO and CAV_IPI data columns)",
    concentration = "hr (relative hazard of a first gr2+ IMAE, unitless; also lhr, the linear predictor)"
  )

  covariateData <- list(
    CAV_NIVO = list(
      description = "Time-varying daily average serum nivolumab concentration (daily Cavg)",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters LINEARLY, uncentred; 0 for patients and days with no",
        "nivolumab (ipilimumab monotherapy). Du 2024 Methods: time-varying",
        "daily Cavg simulated from the individual empirical-Bayes PK",
        "parameters of the pooled adult/pediatric popPK model, chosen",
        "over a per-patient summary because nivolumab and ipilimumab",
        "treatment durations differ (Figure S1 shows the profile for",
        "N1I3 followed by nivolumab 3 mg/kg q2w). The daily average is",
        "AUC over the calendar day divided by one day; supply one value",
        "per day of follow-up. Appendix S1 Table S4 source column CAVGN."
      ),
      source_name = "CAVGN"
    ),
    CAV_IPI = list(
      description = "Time-varying daily average serum ipilimumab concentration (daily Cavg)",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters LINEARLY, uncentred; 0 for patients and days with no",
        "ipilimumab (nivolumab monotherapy, as in the Figure 2b",
        "simulation). Same daily-average construction as CAV_NIVO.",
        "Ipilimumab is the significant exposure driver (HR 1.008 per",
        "ug/mL, 95% CI 1.001-1.014). Appendix S1 Table S4 source column",
        "CAVGI."
      ),
      source_name = "CAVGI"
    ),
    AGE = list(
      description = "Age at baseline",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear on the log-hazard, centred here at the Figure 1 reference of 60 years (Cox coefficients do not depend on centring; the centring only fixes what hr = 1 means). Analysis range 1-90 years, P5-P95 19.1-78 years (Figure 1).",
      source_name = "AGE"
    ),
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear on the log-hazard, centred here at the Figure 1 reference of 75 kg. Analysis range 9.3-183 kg, P5-P95 52-111 kg (Figure 1). A baseline value, constant within a patient.",
      source_name = "WTB"
    ),
    LINE_1L = list(
      description = "First-line-therapy indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (first line) is the paper's reference; the canonical 0 = second line or later carries the coefficient",
      notes = "Du 2024 Table 1 row 'Line of therapy [>=2L:1L]' estimates the effect of second-or-later line against first line, so the model applies it as e_line_2l_haz * (1 - LINE_1L). Line was missing for 23 patients (Table S1); the paper does not say how they were coded.",
      source_name = "LINEM"
    ),
    TUMTP_MEL_ADJUVANT = list(
      description = "Adjuvant (resected) melanoma indicator within the tumour-type / treatment-setting factor",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (advanced, unresectable or metastatic, melanoma when TUMTP_OTHER is also 0)",
      notes = "Du 2024 Table 1 row 'Treatment setting [adjMEL:MEL]'; Figure 1 labels the same factor 'Tumor Type [adjMEL:MEL]' (n = 1352 vs 2043). Mutually exclusive with TUMTP_OTHER; both 0 selects advanced melanoma.",
      source_name = "TTYPEFM"
    ),
    TUMTP_OTHER = list(
      description = "Other-tumour-type / other-treatment-setting indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (advanced melanoma when TUMTP_MEL_ADJUVANT is also 0)",
      notes = "Du 2024 Table 1 row 'Treatment setting [others:MEL]'; Figure 1 'Tumor Type [others:MEL]'. The 'Other' level of Table S1 (n = 112: 39 patients < 12 years, 43 adolescents and 30 adults) is not itemised by histology; it is dominated by the pediatric patients enrolled in non-melanoma pediatric studies. Not interchangeable with another paper's TUMTP_OTHER.",
      source_name = "TTYPEFM"
    ),
    PDL1_TUM_POS = list(
      description = "Tumour PD-L1 expression >= 5% indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (tumour PD-L1 < 5%)",
      notes = "Positivity threshold is 5% of tumour cells (Du 2024 Table 1 row 'PD-L1 status [>=5%:<5%]'). Set to 0 when PDL1_TUM_MISSING = 1; the model also gates the term off for missing PD-L1.",
      source_name = "PDL15M"
    ),
    PDL1_TUM_MISSING = list(
      description = "Tumour PD-L1 result missing or unknown indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (PD-L1 result available)",
      notes = "Du 2024 Table 1 row 'PD-L1 status [missing:<5%]'; 936 of 3507 patients, including all but two of the 97 pediatric patients (Table S1).",
      source_name = "PDL15M"
    ),
    ECOG_GE1 = list(
      description = "Baseline ECOG performance status >= 1 indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ECOG 0)",
      notes = "Du 2024 Table 1 row 'Performance status [>=1:0]'; Figure 1 'PS [>0:0]'.",
      source_name = "PSM"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Du 2024 Table 1 row 'Sex [female:male]'; HR 1.408, the paper's headline 41% higher risk in women.",
      source_name = "SEXM"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White when RACE_BLACK and RACE_OTHER are also 0)",
      notes = "Du 2024 Table 1 row 'Race [Asian:white]'; n = 52.",
      source_name = "RACEM"
    ),
    RACE_BLACK = list(
      description = "Black / African American race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White when RACE_ASIAN and RACE_OTHER are also 0)",
      notes = "Du 2024 Table 1 row 'Race [Black/African American:White]'; n = 15.",
      source_name = "RACEM"
    ),
    RACE_OTHER = list(
      description = "Other or unknown race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White when RACE_ASIAN and RACE_BLACK are also 0)",
      notes = "Du 2024 Table 1 row 'Race [others/unknown:White]' pools 'Other' with 'unknown'; n = 52.",
      source_name = "RACEM"
    ),
    LDHR = list(
      description = "Baseline lactate dehydrogenase divided by the laboratory upper limit of normal",
      units = "(unitless ratio)",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters as the NATURAL log, log(LDHR), so the term is zero at the ULN (Figure 1 reference 'median normalized LDH value of 1'). The natural log is pinned by Figure 1: with -0.0209 per log unit the P5 (0.57) and P95 (2.94) give 1.012 and 0.978, as printed (1.01, 0.978); log10 would give 0.990 at the P95, not 0.978 (at the P05 it gives 1.005, which also rounds to 1.01). Source column LIBLDHR. LDH was missing for 146 patients, 79 of them pediatric (Table S1).",
      source_name = "LIBLDHR"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 3507L,
    n_studies = 15L,
    age_range = "1-90 years (median 58); 42 patients < 12 years, 55 adolescents 12 to < 18 years, 3410 adults (Appendix S1 Table S1)",
    weight_range = "9.3-183 kg (median 79.4)",
    sex_female_pct = 39.2,
    race_ethnicity = c(White = 96.6, `Black/African American` = 0.4, Asian = 1.5, `Other/unknown` = 1.5),
    disease_state = "Primarily melanoma: advanced (58.3%) or adjuvant resected (38.6%); 3.2% other tumour types / settings, mostly pediatric. 81.7% first line.",
    dose_range = "Nivolumab, ipilimumab or nivolumab plus ipilimumab at a wide variety of weight-based and flat IV regimens; treatment split 57.2% nivolumab, 23.9% ipilimumab, 18.9% nivolumab plus ipilimumab",
    regions = "Multinational (Bristol Myers Squibb clinical program)",
    notes = paste(
      "Baseline demographics from Appendix S1 Table S1. ECOG 0 in 78.3%;",
      "PD-L1 >= 5% in 23.4%, < 5% in 49.9%, missing in 26.7%. The",
      "endpoint is time to first gr2+ IMAE regardless of immune-modulating",
      "treatment; ipilimumab-monotherapy studies reported irAEs, which",
      "were used in their place (Table S2)."
    )
  )

  ini({
    # ==================================================================
    # Du 2024 Table 1, "Parameter estimates of the E-R of gr2+ IMAEs
    # (full model)". Estimates are on the log-hazard scale, the scale a
    # Cox model estimates on; the Hazard ratio column is exp(Estimate).
    # Nothing is back-solved or digitised. The fitted formula (Appendix
    # S1 Table S4/S5, bold row) is
    #   Surv(START, STOP, EVENT) ~ CAVGN * CAVGI + AGE + WTB + LINEM +
    #     TTYPEFM + PDL15M + PSM + SEXM + RACEM + LIBLDHR + strata(TRTP)
    # fitted with the R survival package (version 3.3-1). The full model
    # is the final model: Du 2024 reports and applies it without a
    # backward-elimination step.
    #
    # NO BASELINE HAZARD IS ENCODED. A Cox regression is semiparametric,
    # and strata(TRTP) gives each treatment (nivolumab monotherapy,
    # ipilimumab monotherapy, nivolumab plus ipilimumab) its own
    # unspecified h0(t); none is a parameter the fit produced. This
    # model returns the RELATIVE hazard only, matching the pattern of
    # Liu_2024_saf189s_pfs.R and Chan_2025_atezolizumab_os.R.
    #
    # No variance components are reported, so there is no IIV and no
    # residual error, and no observation endpoint is declared.
    # ==================================================================
    e_cav_nivo_haz <- -0.0004655; label("Log hazard ratio per 1 ug/mL of nivolumab daily Cavg (log scale; HR 0.9995, not significant)")                  # Du 2024 Table 1, 'Nivolumab daily C avg (ug/mL)': -0.0004655 (SE 0.0009231), HR 0.9995 (0.9977-1.001)
    e_cav_ipi_haz <- 0.007693; label("Log hazard ratio per 1 ug/mL of ipilimumab daily Cavg (log scale; HR 1.008)")                                          # Du 2024 Table 1, 'Ipilimumab daily C avg (ug/mL)': 0.007693 (SE 0.003228), HR 1.008 (1.001-1.014)
    e_cav_nivo_ipi_haz <- -0.000549; label("Log hazard ratio per (ug/mL)^2 of the nivolumab x ipilimumab daily Cavg product (log scale; HR 0.9995)")         # Du 2024 Table 1, 'C avg nivolumab:C avg ipilimumab': -0.000549 (SE 0.000159), HR 0.9995 (0.9991-0.9998)
    e_age_haz <- 0.00414; label("Log hazard ratio per year of age (log scale; HR 1.004)")                                                                    # Du 2024 Table 1, 'Age (years)': 0.00414 (SE 0.001987), HR 1.004 (1-1.008)
    e_wt_haz <- 0.006033; label("Log hazard ratio per kg of baseline body weight (log scale; HR 1.006)")                                                     # Du 2024 Table 1, 'Body weight (kg)': 0.006033 (SE 0.00156), HR 1.006 (1.003-1.009)
    e_line_2l_haz <- -0.2079; label("Log hazard ratio for second-or-later line versus first line (log scale; HR 0.8123)")                                   # Du 2024 Table 1, 'Line of therapy [>=2L:1L]': -0.2079 (SE 0.09439), HR 0.8123 (0.6751-0.9774)
    e_tumtp_mel_adjuvant_haz <- -0.2972; label("Log hazard ratio for adjuvant versus advanced melanoma (log scale; HR 0.7429)")                             # Du 2024 Table 1, 'Treatment setting [adjMEL:MEL]': -0.2972 (SE 0.0889), HR 0.7429 (0.6241-0.8843)
    e_tumtp_other_haz <- 0.3119; label("Log hazard ratio for other tumour type / setting versus advanced melanoma (log scale; HR 1.366)")                   # Du 2024 Table 1, 'Treatment setting [others:MEL]': 0.3119 (SE 0.207), HR 1.366 (0.9105-2.05)
    e_pdl1_pos_haz <- -0.02175; label("Log hazard ratio for tumour PD-L1 >= 5% versus < 5% (log scale; HR 0.9785)")                                         # Du 2024 Table 1, 'PD-L1 status [>=5%:<5%]': -0.02175 (SE 0.06456), HR 0.9785 (0.8622-1.11)
    e_pdl1_missing_haz <- -0.1229; label("Log hazard ratio for missing tumour PD-L1 versus < 5% (log scale; HR 0.8844)")                                    # Du 2024 Table 1, 'PD-L1 status [missing:<5%]': -0.1229 (SE 0.08553), HR 0.8844 (0.7479-1.046)
    e_ecog_ge1_haz <- 0.0877; label("Log hazard ratio for ECOG performance status >= 1 versus 0 (log scale; HR 1.092)")                                     # Du 2024 Table 1, 'Performance status [>=1:0]': 0.0877 (SE 0.06695), HR 1.092 (0.9574-1.245)
    e_sexf_haz <- 0.3423; label("Log hazard ratio for female versus male (log scale; HR 1.408)")                                                             # Du 2024 Table 1, 'Sex [female:male]': 0.3423 (SE 0.05882), HR 1.408 (1.255-1.58)
    e_race_asian_haz <- 0.2104; label("Log hazard ratio for Asian versus White race (log scale; HR 1.234)")                                                 # Du 2024 Table 1, 'Race [Asian:white]': 0.2104 (SE 0.2139), HR 1.234 (0.8115-1.877)
    e_race_black_haz <- -0.3134; label("Log hazard ratio for Black / African American versus White race (log scale; HR 0.7309)")                           # Du 2024 Table 1, 'Race [Black/African American:White]': -0.3134 (SE 0.4529), HR 0.7309 (0.3008-1.776)
    e_race_other_haz <- -0.3504; label("Log hazard ratio for other / unknown versus White race (log scale; HR 0.7044)")                                     # Du 2024 Table 1, 'Race [others/unknown:White]': -0.3504 (SE 0.2276), HR 0.7044 (0.4509-1.101)
    e_log_ldhr_haz <- -0.0209; label("Log hazard ratio per natural-log unit of baseline LDH / ULN (log scale; HR 0.9793)")                                  # Du 2024 Table 1, 'Log LDH (xULN)': -0.0209 (SE 0.04906), HR 0.9793 (0.8895-1.078)
  })

  model({
    # ------------------------------------------------------------------
    # Linear predictor against the Du 2024 Figure 1 reference patient:
    # 60 years, 75 kg, white, male, first-line advanced melanoma,
    # ECOG 0, tumour PD-L1 < 5%, LDH at the ULN, no drug exposure.
    # Centring age and weight changes only what hr = 1 refers to, not
    # any coefficient. The subject hazard is h(t) = h0_s(t) * hr(t),
    # where h0_s is the unspecified baseline hazard of the patient's
    # treatment stratum s; hr varies in time only through the daily
    # Cavg columns.
    #
    # The interaction is the plain product of the two daily Cavg
    # values, so it is zero in both monotherapy strata and acts only
    # while both drugs are present. The missing-PD-L1 indicator
    # switches the >= 5% term off so that the two PD-L1 levels stay
    # mutually exclusive.
    # ------------------------------------------------------------------
    lhr_exposure <- e_cav_nivo_haz * CAV_NIVO +
      e_cav_ipi_haz * CAV_IPI +
      e_cav_nivo_ipi_haz * CAV_NIVO * CAV_IPI

    lhr_baseline <- e_age_haz * (AGE - 60) +
      e_wt_haz * (WT - 75) +
      e_line_2l_haz * (1 - LINE_1L) +
      e_tumtp_mel_adjuvant_haz * TUMTP_MEL_ADJUVANT +
      e_tumtp_other_haz * TUMTP_OTHER +
      e_pdl1_pos_haz * PDL1_TUM_POS * (1 - PDL1_TUM_MISSING) +
      e_pdl1_missing_haz * PDL1_TUM_MISSING +
      e_ecog_ge1_haz * ECOG_GE1 +
      e_sexf_haz * SEXF +
      e_race_asian_haz * RACE_ASIAN +
      e_race_black_haz * RACE_BLACK +
      e_race_other_haz * RACE_OTHER +
      e_log_ldhr_haz * log(LDHR)

    lhr <- lhr_exposure + lhr_baseline

    hr <- exp(lhr)
  })
}
