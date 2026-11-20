Daud_2011_peginterferon_alfa_2b_anc_auc <- function() {
  description <- paste(
    "Static exposure-response (inhibitory Imax) model relating the",
    "steady-state weekly AUCtau of subcutaneous peginterferon alfa-2b to the",
    "absolute neutrophil count expressed as percent of baseline (Daud 2011),",
    "fitted to dose-group mean data pooled from adults with resected",
    "high-risk melanoma and adults with advanced solid tumours. ANC declines",
    "hyperbolically from E0 towards the reported Imax as AUCtau increases."
  )
  reference <- paste(
    "Daud AI, Xu C, Hwu WJ, Urbas P, Andrews S, Papadopoulos NE, Floren LC,",
    "Yver A, DeConti RC, Sondak VK.",
    "Pharmacokinetic/pharmacodynamic analysis of adjuvant pegylated",
    "interferon alpha-2b in patients with resected high-risk melanoma.",
    "Cancer Chemother Pharmacol. 2011;67(3):657-666.",
    "doi:10.1007/s00280-010-1326-9"
  )
  vignette <- "Daud_2011_peginterferon_alfa_2b"
  units <- list(
    time = "n/a (static steady-state exposure-response relationship; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the AUC_PEGINTERFERON_ALFA_2B covariate column)",
    concentration = "ANC as percent of baseline (ANC)"
  )

  covariateData <- list(
    AUC_PEGINTERFERON_ALFA_2B = list(
      description = "Steady-state area under the serum peginterferon alfa-2b concentration-time curve over the once-weekly (168 h) dosing interval, AUCtau.",
      units = "pg*h/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Daud 2011 used the observed dose-group mean AUCtau at week 8 (melanoma, 6 ug/kg/week) and week 12 (melanoma, 3 ug/kg/week), and previously published week-4 AUCtau values from patients with solid tumours given 0.75-7.5 ug/kg/week (Bukowski 2002). Supply the steady-state weekly AUCtau, e.g. dose / (CL/F) from modellib('Daud_2011_peginterferon_alfa_2b') where CL/F = kel * V/F. Set to 0 for no drug; the drug term then vanishes and ANC = E0.",
      source_name = "AUCtau"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 66L,
    n_studies = 2L,
    disease_state = "Pooled: 32 adults with resected high-risk (stage IIB-IIIC) melanoma from Daud 2011 and 34 adults with advanced solid tumours (6 of them melanoma) from the phase I/II study of Bukowski 2002.",
    dose_range = "Peginterferon alfa-2b SC once weekly: 0.75, 1.5, 3, 4.5, 6 and 7.5 ug/kg/week (solid tumours, 3-12 patients per dose group); 6 and 3 ug/kg/week (melanoma).",
    notes = "Fitted to dose-group mean AUCtau and mean ANC change (percent of baseline) in WinNonlin (Daud 2011 Fig. 2; supplementary Table 1 'Imax model parameters (pooled data, n = 3 to 32)', R = 0.986, AIC 46.0). The solid-tumour ANC data were unpublished data provided by Schering-Plough. Demographics of the pooled cohort are not reported; for the melanoma cohort see Daud_2011_peginterferon_alfa_2b."
  )

  ini({
    # Inhibitory Imax model on AUCtau (Daud 2011 supplementary Methods;
    # supplementary Table 1 'Imax model parameters (pooled data, n=3 to 32)').
    # The supplementary Table %CV column is the precision of each estimate,
    # not between-subject variability: the model was fitted to dose-group
    # means, so it carries no IIV.
    lrbase <- log(98.7); label("Baseline effect E0, ANC at zero exposure (% of baseline)") # Supp Table 1: E0 = 98.7% (CV 3.8%)
    limax <- log(33.3); label("Reported Imax, the ANC approached at saturating exposure (% of baseline)") # Supp Table 1: Imax = 33.3% (CV 11.5%); role confirmed against Fig. 2 (see vignette)
    lauc50 <- log(80900); label("AUCtau giving half of the maximal ANC decrease, IC50 (pg*h/mL)") # Supp Table 1: IC50 = 80,900 pg*hr/mL (CV 28.2)

    # Residual error -- not reported; fixed to zero.
    addSd_ANC <- fixed(0); label("Additive residual error on ANC (% of baseline); not reported") # not reported in Daud 2011
  })

  model({
    rbase <- exp(lrbase)
    imax <- exp(limax)
    auc50 <- exp(lauc50)

    # The supplementary Methods print E = E0 - Imax*C/(IC50 + C), but the
    # model-fit curve in Daud 2011 Fig. 2 and the dose-group means level off
    # near 40% of baseline, which only the form below reproduces with the
    # tabulated values: Imax is the response at saturating AUCtau.
    ANC <- rbase - (rbase - imax) * AUC_PEGINTERFERON_ALFA_2B / (auc50 + AUC_PEGINTERFERON_ALFA_2B)

    ANC ~ add(addSd_ANC)
  })
}
