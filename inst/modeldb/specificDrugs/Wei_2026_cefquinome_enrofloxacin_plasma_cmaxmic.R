Wei_2026_cefquinome_enrofloxacin_plasma_cmaxmic <- function() {
  description <- paste0(
    "Preclinical (chick, yellow-feathered broiler). Inhibitory sigmoid Emax PK/PD-index model for ",
    "the antibacterial effect of cefquinome (CEQ) given together with a fixed 20 mg/kg background ",
    "dose of enrofloxacin (ENR) against the multidrug-resistant Klebsiella quasipneumoniae subsp. ",
    "similipneumoniae clinical isolate CLS2, in an intratracheal chick pneumonia model. This file ",
    "carries the fit driven by the Cmax/MIC (combine) index computed from plasma free drug PK data. ",
    "Wei 2026 'Analysis of antibacterial effects and fitting of PK/PD' parameterises the effect as ",
    "E = Emax - (Emax - E0) * Ce^N / (EC50^N + Ce^N), where E is the SIGNED difference in lung ",
    "K. pneumoniae load in log10 CFU/mL between the treated group and the untreated blank control ",
    "at 72 h (negative = bacterial reduction), and Ce is the PK/PD index. ",
    "NOTE THE INVERTED NAMING, identical to the sibling model Gao_2025_cefquinome_pkpd_index: ",
    "Wei 2026 defines its Emax as the change in the BLANK CONTROL group and its E0 as the maximum ",
    "antibacterial effect, which is the reverse of the usual convention. This file therefore maps ",
    "Wei 2026's Emax = 0.001 onto the canonical `e0` (effect at zero exposure) and Wei 2026's ",
    "E0 = -4.039 onto the canonical `emax` (maximal drug effect), and writes the algebraically ",
    "identical E = e0 - (e0 - emax) * Ce^N / (Ce^N + EC50^N). No value is altered by the remapping. ",
    "Parameters from Wei 2026 Table 6, plasma free drug / Cmax/MIC (combine) column: e0 = 0.001, emax = -4.039 ",
    "log10 CFU/mL, EC50 = 57.921, Hill N = 1.556 (R^2 = 0.803, AIC = 20.905). ",
    "Wei 2026 reports no 3 log10 CFU/mL target index for this fit (the Table 6 row is a dash); it calculated targets only for the two better-correlating plasma indices. ",
    "There is NO PK component: Wei 2026 analysed the plasma and lung concentrations ",
    "non-compartmentally in Phoenix 8.4 (Tables 3 and 4 report only Tmax, Cmax, T1/2beta and ",
    "AUC72h) and published no structural PK model, so exposure enters as the externally supplied ",
    "covariate CMAXMIC_CEFQ. That covariate carries the PK/PD index as an ALREADY-FORMED RATIO rather ",
    "than an absolute exposure divided by a model `mic` parameter, because the relevant ",
    "susceptibility here is not a fixed strain property but MIC(combine) -- the MIC of CEQ ",
    "measured in the presence of the average ENR concentration achieved in that matrix -- and it ",
    "is NOT constant across the arms Wei 2026 fitted: reconstructing Table 5 from Tables 3 and 4 ",
    "gives plasma MIC(combine) = 0.031 ug/mL for all six arms, exactly as Wei 2026 states, and that single value reproduces every plasma Cmax/MIC entry of Table 5 to five significant figures. ",
    "The model predicts the control-referenced CHANGE directly rather than integrating a bacterial ",
    "density, because Wei 2026's readout is a cross-sectional lung count per dose group at 72 h ",
    "expressed relative to the blank control, whose absolute load is never reported. ",
    "Wei 2026 reports neither between-subject variability nor a residual error magnitude for the ",
    "PK/PD fit, so there are no eta parameters and addSd is FIXED at 0; the model is intended for ",
    "typical-value simulation. Sibling models from the same paper: Wei_2026_cefquinome_enrofloxacin_plasma_aucmic, Wei_2026_cefquinome_enrofloxacin_plasma_tmic, Wei_2026_cefquinome_enrofloxacin_lung_cmaxmic, Wei_2026_cefquinome_enrofloxacin_lung_aucmic and Wei_2026_cefquinome_enrofloxacin_lung_tmic."
  )
  reference <- paste(
    "Wei Y, Zhang F, Liu X, Li X, Li T, Li Y, Ding H. (2026).",
    "Pharmacokinetic/pharmacodynamic relationships and development of resistance of",
    "enrofloxacin and cefquinome in combination therapy against Klebsiella pneumoniae in chicks.",
    "BMC Veterinary Research 22:258.",
    "doi:10.1186/s12917-026-05301-5. PMCID: PMC13134116.",
    "Model equation: Materials and methods, 'Analysis of antibacterial effects and fitting of",
    "PK/PD'. Parameter estimates: Table 6, plasma free drug / Cmax/MIC (combine) column.",
    "PK/PD index and effect data: Table 5. Non-compartmental PK: Tables 3 and 4.",
    sep = " "
  )
  vignette <- "Wei_2026_enrofloxacin_cefquinome"

  units <- list(
    time = "h",
    dosing = "unitless (cefquinome Cmax/MIC (combine) index, supplied as a covariate)",
    concentration = "log10 CFU/mL (observation)"
  )

  depends <- c("CMAXMIC_CEFQ")

  covariateData <- list(
    CMAXMIC_CEFQ = list(
      description        = "Cefquinome PK/PD index: the peak free plasma concentration of cefquinome divided by MIC(combine), the MIC of cefquinome measured in the presence of the average free enrofloxacin plasma concentration achieved over 72 h",
      units              = "unitless",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Wei 2026 tabulates this index for all six combination-therapy arms in Table 5, separately ",
        "for plasma free drug and for lung. It is carried as an already-formed ratio rather than ",
        "as an absolute exposure divided by a model `mic` parameter because the denominator is ",
        "MIC(combine), the MIC of cefquinome measured by CAMH broth microdilution in the presence ",
        "of the average free enrofloxacin concentration achieved over 72 h in that matrix ",
        "(Materials and methods, 'Determination of the MIC combine of cefquinome against ",
        "K. pneumoniae'). MIC(combine) is a property of the ENR co-exposure and dosing schedule, ",
        "not of the isolate alone, and it is not constant across the fitted arms: ",
        "plasma MIC(combine) = 0.031 ug/mL for all six arms, exactly as Wei 2026 states, and that single value reproduces every plasma Cmax/MIC entry of Table 5 to five significant figures. ",
        "Set to 0 for the untreated blank control so the sigmoid term vanishes and the predicted ",
        "72 h change reduces to e0. See the vignette Assumptions and deviations section."
      ),
      source_name        = "Cmax/MIC (combine) (Wei 2026 Table 5, plasma free drug block; Materials and methods, 'Analysis of antibacterial effects and fitting of PK/PD')"
    )
  )

  population <- list(
    species        = "chick (yellow-feathered, Guangdong Academy of Agricultural Sciences), infected with Klebsiella pneumoniae",
    n_subjects     = 352L,
    n_studies      = 1L,
    age_range      = "purchased at one day old, observed for one week, challenged intratracheally at seven days old; treatment began 24 h after infection",
    organism       = "Klebsiella quasipneumoniae subsp. similipneumoniae clinical isolate CLS2 (identified by 16S rDNA sequencing), provided by the Veterinary Pharmacology Laboratory, South China Agricultural University. MIC 2.00 ug/mL for enrofloxacin and 0.13 ug/mL (0.125) for cefquinome. Whole-genome sequencing identified 18 antibiotic-resistance genes (Table 2), including the efflux pumps acrB, mdtQ, kpnE, oqxB, smeE and tet(A), the beta-lactamase shv-209 and the porins ompA and ompK37; no quinolone-resistance-determining-region chromosomal mutations were found. Escherichia coli ATCC 25922 was the quality-control strain",
    system         = "Intratracheal pneumonia model: 0.4 mL of a 10^8 CFU/mL CLS2 saline suspension injected into the respiratory tract. Treatment started 24 h after infection and continued for three days; lungs from eight chicks per group were harvested at 24, 48 and 72 h after the first dose, homogenised (1.0 g in 1.0 mL sterile saline) and plated by tenfold serial dilution on MacConkey agar with and without 1x MIC of ENR or CEQ. Plate-count detection limit 50 CFU/mL. Group allocation and counting were blinded by letter labels",
    disease_state  = "experimentally induced Klebsiella pneumoniae pneumonia",
    dose_range     = "All treatment arms received a fixed 20 mg/kg intramuscular enrofloxacin background. Single-dose combination arms gave ENR 20 mg/kg plus CEQ 2, 5 or 20 mg/kg once daily for three days; split-dose combination arms gave the same daily totals of ENR 20 mg/kg plus CEQ 4, 10 or 20 mg/kg divided into two equal doses 12 h apart for three days. Monotherapy control arms received ENR 20 mg/kg or CEQ 20 mg/kg, single or split",
    regions        = "China (South China Agricultural University, Guangzhou)",
    notes          = "The 352 chicks are the pharmacodynamic cohort (11 groups, n = 8 per group per time point); a separate 600-chick cohort in six groups supplied the non-compartmental pharmacokinetics of Tables 3 and 4 (n = 6 per group per time point, single sample per bird by cardiac puncture after cervical dislocation). Plasma protein binding was 20.18% for enrofloxacin (from previously published data) and 13.24% for cefquinome by equilibrium dialysis, so the plasma PK of Tables 3 and 4 is reported as FREE drug. Enrofloxacin concentrations are the sum of enrofloxacin and its active metabolite ciprofloxacin. Wei 2026 fitted the same inhibitory sigmoid Emax structure independently against three PK/PD indices (Cmax/MIC, AUC72h/MIC and %T>MIC) in each of two matrices (plasma free drug and lung), and reports all six parameter sets side by side in Table 6, so all six are packaged. The paper's own model selection favours %T>MIC in plasma (R^2 = 0.953, AIC = 14.695), concluding that cefquinome behaves as a time-dependent agent when combined with enrofloxacin, and its recommended regimen is ENR 10 mg/kg plus CEQ 2.6 mg/kg intramuscularly twice daily for three days"
  )

  ini({
    # =================================================================
    # Wei 2026 Table 6 -- inhibitory sigmoid Emax PK/PD-index model
    # matrix: plasma free drug   index: Cmax/MIC (combine)
    # =================================================================
    # Wei 2026 equation as printed (Materials and methods, 'Analysis of
    # antibacterial effects and fitting of PK/PD'):
    #   E = Emax - (Emax - E0) * Ce^N / (EC50^N + Ce^N)
    # with Wei 2026's own definitions, which INVERT the usual roles:
    #   "Emax: Change in bacterial count in the blank control group"
    #   "E0: Maximum antibacterial effect produced by the drug"
    # Table 6 confirms the inversion: the value that the curve takes at
    # Ce = 0 is small and non-negative (the control row of Table 5 is
    # E = 0.00 by construction, because E is measured RELATIVE to the
    # blank control), while the saturating value is a large negative
    # number (net kill).
    #
    # This file uses the canonical roles, so:
    #   e0   <- Wei 2026's Emax = 0.001   (effect at zero exposure)
    #   emax <- Wei 2026's E0   = -4.039  (maximal drug effect)
    # and model() writes the algebraically identical
    #   E = e0 - (e0 - emax) * Ce^N / (Ce^N + EC50^N)
    # No value is changed by the remapping -- only the label each
    # number carries. The same remapping is applied in the library's
    # other cefquinome PK/PD-index model, Gao_2025_cefquinome_pkpd_index.
    #
    # E is a SIGNED log10 CFU/mL difference from the blank control, so
    # emax is negative and neither e0 nor emax can be log-transformed;
    # both stay on the natural scale. EC50 and the Hill coefficient are
    # strictly positive and are carried as logs.
    e0 <- 0.001
    label("Difference in lung K. pneumoniae load from the blank control at zero exposure (log10 CFU/mL)")  # Wei 2026 Table 6, plasma free drug block, Cmax/MIC (combine) column, "E max" = 0.001 (defined in Materials and methods as the blank-control change)
    emax <- -4.039
    label("Maximal difference in lung K. pneumoniae load from the blank control (log10 CFU/mL; negative = reduction)")  # Wei 2026 Table 6, plasma free drug block, Cmax/MIC (combine) column, "E 0" = -4.039 (defined in Materials and methods as the maximum antibacterial effect)
    lec50 <- log(57.921)
    label("Log Cmax/MIC (combine) index producing 50% of the maximal antibacterial effect EC50 (unitless)")  # Wei 2026 Table 6, plasma free drug block, Cmax/MIC (combine) column, EC 50 = 57.921
    lhill <- log(1.556)
    label("Log Hill coefficient N defining the steepness of the index-effect curve (unitless)")  # Wei 2026 Table 6, plasma free drug block, Cmax/MIC (combine) column, Hill's slope = 1.556

    # =================================================================
    # Residual error
    # =================================================================
    # Wei 2026 Table 6 reports the correlation coefficient (R^2 = 0.803)
    # and the AIC (20.905) of this fit and gives no residual standard
    # deviation and no standard errors on the estimates, so the residual
    # SD is held at zero for deterministic typical-value simulation.
    # See the vignette Assumptions and deviations section.
    addSd <- fixed(0)
    label("Additive residual SD on the control-referenced 72 h change in log10 CFU/mL (0; not reported in Wei 2026)")  # Wei 2026 Table 6 reports R^2 and AIC only, no residual SD and no standard errors
  })

  model({
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # PK/PD index driving the sigmoid, supplied as an already-formed
    # ratio (see the covariateData notes for why it is not split into
    # an absolute exposure divided by a mic parameter).
    ce <- CMAXMIC_CEFQ

    # Wei 2026 inhibitory sigmoid Emax equation, written with canonical
    # parameter roles (see the ini() note on the inverted naming). Cc is
    # the SIGNED difference in lung K. pneumoniae load from the blank
    # control at 72 h in log10 CFU/mL (negative = reduction); it equals
    # e0 at zero exposure and approaches emax at saturating exposure.
    Cc <- e0 - (e0 - emax) * ce^hill / (ce^hill + ec50^hill)
    Cc ~ add(addSd)
  })
}
