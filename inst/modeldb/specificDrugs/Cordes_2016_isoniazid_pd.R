Cordes_2016_isoniazid_pd <- function() {
  description <- paste0(
    "QSP. Semi-mechanistic pharmacodynamic model of Mycobacterium tuberculosis killing in the ",
    "human lung during isoniazid (INH) monotherapy, giving the early bactericidal activity (EBA) ",
    "as the net of mycobacterial growth, host immune killing, and a sigmoid Emax INH effect. ",
    "Cordes 2016 equation 3: dN/dt = mu - beta0 - Emax * (C/MIC)^h / (Km^h + (C/MIC)^h), where N ",
    "is the mycobacterial burden in log10 CFU/liter and C is the UNBOUND INH concentration in the ",
    "interstitial space of the lung. All three terms are rates in log10 CFU/day, so the state is ",
    "integrated here as d/dt(bact) = ln(10) * knet * bact, which makes log10(bact) change by ",
    "exactly knet per day while keeping the state on the linear CFU/liter scale used by the ",
    "sibling bacterial-dynamics models in this library. Immune status enters as Cordes 2016's ",
    "factor d (here the canonical f_immune, in [0, 1]): beta* = beta0 * f_immune, with 1 = fully ",
    "immunocompetent and 0 = fully immune-deficient, so the user can sweep the severity of an ",
    "HIV- or immunosuppressant-driven deficiency. ",
    "THERE IS NO PK COMPONENT IN THIS FILE, AND THAT IS DELIBERATE. Cordes 2016's drug input is a ",
    "whole-body PBPK model of INH and six metabolites built in the proprietary platform PK-Sim ",
    "6.0.3; its tissue:plasma partition coefficients and its tissue-specific NAT2 / NAAA / GLYAT / ",
    "SLC16A10 abundance profiles are generated inside the platform from a named method and a ",
    "built-in gene-expression database, and are not published in the paper, the supplement, or a ",
    "deposited project file. The paper prints no clearance, no volume, no AUC and no Cmax ",
    "anywhere, so no reduced compartmental PK model can be recovered from it either. Exposure ",
    "therefore enters this file as the externally supplied covariate CEFFECT, on the same pattern ",
    "as Gao_2025_cefquinome_pkpd_index and Crass_2025_pegcetacoplan_ga_exposureresponse. Pair it ",
    "with any INH PK model that can report an unbound lung-interstitial concentration; ",
    "Vinnard_2017_isoniazid supplies a NAT2-aware apparent-oral-clearance popPK model for INH. ",
    "TWO PUBLISHED-VALUE CONFLICTS ARE RESOLVED HERE AND DOCUMENTED IN THE VIGNETTE ERRATA. (1) ",
    "The growth rate mu is taken as 0.0428 log10 CFU/day from the supplement, not the 0.048 ",
    "printed in Table 5, because Table 5's own beta0 = 0.0219 is exactly the supplement's ",
    "mu_ID(0.0428) - mu_IC(0.0209) and only mu = 0.0428 reproduces the paper's literature-derived ",
    "untreated immunocompetent growth rate of 0.0209 log10 CFU/day. (2) Cordes 2016's Km is ",
    "carried as the canonical ec50 on the DIMENSIONLESS C/MIC scale (half-maximal effect at ",
    "C = 25.19 * MIC = 36.8 umol/L), because that is the only reading under which the printed ",
    "equation 2 is dimensionally valid; Table 5's 'umol/liter' annotation on Km and the ",
    "accompanying prose definition imply instead a half-effect at 25.19 umol/L. ",
    "Deterministic typical-value model: Cordes 2016 reports no between-subject variability and no ",
    "residual error for the PD parameters, so there are no eta terms and addSd is FIXED at 0."
  )
  reference <- paste(
    "Cordes H, Thiel C, Aschmann HE, Baier V, Blank LM, Kuepfer L. (2016).",
    "A physiologically based pharmacokinetic model of isoniazid and its application in",
    "individualizing tuberculosis chemotherapy.",
    "Antimicrobial Agents and Chemotherapy 60(10):6134-6145.",
    "doi:10.1128/AAC.00508-16.",
    sep = " "
  )
  vignette <- "Cordes_2016_isoniazid"

  units <- list(
    time = "day",
    dosing = "umol/L (unbound isoniazid in lung interstitial fluid, supplied as the covariate CEFFECT)",
    concentration = "log10 CFU/L (observation)"
  )

  depends <- c("CEFFECT")
  paper_specific_compartments <- c("bact")

  compartmentData <- list(
    bact = list(
      analyte = "Mycobacterium tuberculosis",
      units = "CFU/L",
      specimen = "tissue",
      verified = TRUE,
      description = paste0(
        "Mycobacterial burden at the site of infection, the interstitial space of the human ",
        "lung. Cordes 2016 calibrated the model against sputum CFU counts (the EBA assay) and ",
        "treats them as reporting on the lung burden."
      )
    )
  )

  covariateData <- list(
    CEFFECT = list(
      description = paste0(
        "Unbound isoniazid concentration in the interstitial space of the lung, the on-target ",
        "site of infection. Cordes 2016 generated this trajectory with its PK-Sim whole-body ",
        "PBPK model and fed it to the PD model as the effective drug input (Results, 'PBPK/PD ",
        "model development'); the paper argues it is a reasonable proxy for on-target ",
        "availability because INH concentrations in plasma, epithelial lining fluid and alveolar ",
        "cells do not differ significantly (Discussion, citing Conte 2002). The user supplies ",
        "the trajectory; set it to 0 for the untreated reference."
      ),
      units = "umol/L",
      type = "continuous",
      source_name = "C",
      reference_category = NA_character_,
      notes = paste0(
        "Time-varying. Cordes 2016 equation 2 divides C by the MIC before it enters the sigmoid, ",
        "so CEFFECT must be a molar concentration on the same scale as the mic parameter ",
        "(umol/L). To convert from mass units, INH molecular weight is 137.14 g/mol (Table 1), ",
        "so 1 mg/L = 7.292 umol/L."
      )
    )
  )

  population <- list(
    species = "human",
    disease_state = "active pulmonary tuberculosis (PD calibration); healthy volunteers and tuberculosis patients (PBPK layer, not ported here)",
    notes = paste0(
      "The PD parameters Emax, Km and h were fitted to early-bactericidal-activity data from ",
      "sputum of NAT2-phenotype-specific pulmonary tuberculosis patients over the first two days ",
      "of INH monotherapy (Donald 1997 AJRCCM 156:895-900, main-text reference 23; Donald 2004 ",
      "Clin Infect Dis 39:1425-1430, main-text reference 45). Supplement, 'PD Model Construction ",
      "& Validation': only patient subgroups with more than three individuals were used, to limit ",
      "the effect of outliers, and the original studies' sampling patterns were reproduced in the ",
      "simulations. QD INH doses spanning 9 mg to 600 mg are represented in the calibration set ",
      "(Fig. 5B). The untreated growth rate mu and the immune killing rate beta0 were not fitted: ",
      "they were derived from a literature review of M. tuberculosis growth rates in untreated ",
      "humans and in immunocompetent versus immune-deficient mice (Supplement equations A2-A3). ",
      "Cordes 2016's population work simulated 1,000 virtual individuals per NAT2 acetylator ",
      "phenotype, but that variability lives entirely in the PBPK layer (anatomy and physiology, ",
      "supplement Table S2) and so is not represented in this PD-only file; the per-individual ",
      "immune-competence factors sampled for the immune-deficient populations are supplement ",
      "Table S3."
    )
  )

  ini({
    # =====================================================================
    # Cordes 2016 equation 3, as printed:
    #
    #   dN(t)/dt = N0 * [mu - beta0 - Emax*(C/MIC)^h / (Km^h + (C/MIC)^h)]
    #
    # THE PRINTED 'N0 *' FACTOR IS NOT APPLIED HERE. Three independent
    # lines of evidence say the rate balance carries no N0 multiplier:
    #   1. The supplement writes the same model as equation A1 with no
    #      such factor at all:
    #        mu(MT)_IC^human = mu(MT)_ID^human - beta^human - gamma(INH)
    #      and gives every term units of log10 CFU/day.
    #   2. Only without the factor does the paper's own arithmetic close:
    #      mu - beta0 = 0.0428 - 0.0219 = 0.0209 log10 CFU/day, exactly
    #      the literature-averaged untreated immunocompetent growth rate
    #      the supplement reports. With N0 = 10 it would be 0.209.
    #   3. Only without the factor is the saturated kill rate
    #      (Emax = 0.534 log10 CFU/day) the right size for isoniazid,
    #      whose measured 2-day EBA is about 0.5 log10 CFU/day. With
    #      N0 = 10 the model would predict 5.3 log10 CFU/day.
    # Because dN/dt does not depend on N, N0 is a pure initial condition
    # and Table 5 itself calls it 'Arbitrary'.
    # =====================================================================

    # ---- Mycobacterial growth and host immune killing -------------------
    # Both are literature-derived rather than fitted, hence fixed().
    #
    # mu: DELIBERATELY 0.0428, NOT Table 5's printed 0.048. Supplement
    # equation A2 derives mu(MT)_ID^human = 0.0428 log10 CFU/day (the
    # growth rate with no immune contribution) and equation A3 then gives
    # beta^human = 0.0428 - 0.0209 = 0.0219, which is exactly the beta0
    # Table 5 prints. Table 5's mu = 0.048 is inconsistent with Table 5's
    # own beta0: it would put the untreated immunocompetent growth rate at
    # 0.048 - 0.0219 = 0.0261 rather than the 0.0209 log10 CFU/day the
    # supplement derives from the literature. See the vignette Errata.
    lkgrowth <- fixed(log(0.0428))
    label("Log M. tuberculosis growth rate with no immune contribution (log10 CFU/day)")

    lkimm <- fixed(log(0.0219))
    label("Log immune-system-mediated M. tuberculosis killing rate (log10 CFU/day)")

    # Cordes 2016 Materials and Methods: 'beta0 is replaced by
    # beta* = beta0 * d, where d accounts for the strength of the immune
    # response in an individual (the value of d ranges from 1 for a fully
    # immunocompetent individual to 0 for a fully immune-deficient
    # individual)'. Carried on the natural scale as a [0, 1] fraction, on
    # the Dogra_2023_covid19vaccine pattern. Override per simulation with
    # rxode2::rxSolve(mod, params = c(f_immune = 0.3), ...); the
    # per-individual values Cordes 2016 sampled for its immune-deficient
    # populations are supplement Table S3.
    f_immune <- fixed(1)
    label("Fractional strength of the host immune response (1 = immunocompetent, 0 = fully immune-deficient)")

    # ---- Isoniazid effect: inhibitory sigmoid Emax (Table 5, 'Fitted') --
    lemax <- log(0.534)
    label("Log maximal isoniazid-induced M. tuberculosis killing rate (log10 CFU/day)")  # Table 5, row 'Emax' = 0.534 /day, Fitted

    # Cordes 2016's 'Km', carried under its canonical role name ec50.
    #
    # UNITS ARE THE ONE PLACE THIS FILE DEPARTS FROM Table 5. Equation 2
    # reads Emax*(C/MIC)^h / (Km^h + (C/MIC)^h). The second term of that
    # denominator is dimensionless, so Km must be dimensionless too, i.e.
    # a MULTIPLE OF THE MIC -- half-maximal effect at C = 25.19 * MIC =
    # 36.8 umol/L. Table 5 instead annotates Km as 'umol/liter' and the
    # surrounding prose calls it 'the INH concentration at which half the
    # maximal antimicrobial effect is reached', which would put the
    # half-effect at 25.19 umol/L and require equation 2's denominator to
    # have been (Km/MIC)^h. The equation is followed here because it is
    # the only reading that is dimensionally valid as printed, and because
    # the C/MIC normalisation is inherited from Czock & Keller 2007
    # (reference 47), whose antimicrobial sigmoid is written in MIC
    # multiples. The two readings differ by a factor of MIC = 1.46, and
    # the paper reports no absolute lung concentration against which they
    # could be discriminated. See the vignette Errata.
    lec50 <- log(25.19)
    label("Log C/MIC ratio giving half-maximal isoniazid killing (dimensionless multiple of MIC)")  # Table 5, row 'Km' = 25.19, Fitted

    lhill <- log(0.56)
    label("Log Hill coefficient of the isoniazid killing sigmoid (unitless)")  # Table 5, row 'h' = 0.56, Fitted

    # ---- MIC ------------------------------------------------------------
    # Measured susceptibility of the target organism, not an estimated
    # parameter, hence fixed(). Cordes 2016 Discussion: 'we chose a
    # conservative estimate of the MIC of 0.2 mg/liter to include most
    # INH-susceptible M. tuberculosis strains (0.05 mg/liter < MIC <
    # 0.1 mg/liter for most strains)'. 0.2 mg/L / 137.14 g/mol =
    # 1.458 umol/L, the 1.46 umol/L of Table 5. Change it to apply the
    # model to an isolate of different susceptibility; the paper notes the
    # conservative choice makes the model underestimate EBA.
    mic <- fixed(1.46)
    label("Isoniazid MIC against M. tuberculosis (umol/L)")  # Table 5, row 'MIC' = 1.46 umol/liter, from reference 41 (Schon 2009)

    # ---- Initial mycobacterial burden ------------------------------------
    # dN/dt does not depend on N, so this only sets where the trajectory
    # starts; Table 5 marks it 'Arbitrary'.
    log10_cfu0 <- fixed(10)
    label("Initial M. tuberculosis burden (log10 CFU/L)")  # Table 5, row 'N0' = 10 log10 CFU/liter, Arbitrary

    # ---- Residual error --------------------------------------------------
    # Cordes 2016 reports goodness of fit for the PD layer only as a
    # correlation across the EBA data set (R^2 = 0.6, P < 0.001; Results),
    # with no residual standard deviation, so the residual SD is held at
    # zero for deterministic typical-value simulation.
    addSd <- fixed(0)
    label("Additive residual SD on log10 CFU/L (not reported in Cordes 2016)")  # Results: R^2 = 0.6 only, no residual SD reported
  })

  model({
    kgrowth <- exp(lkgrowth)
    kimm <- exp(lkimm)
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # Cordes 2016 equation 2. The drug input is normalised by the MIC
    # before it enters the sigmoid, and ec50 is on that same C/MIC scale
    # (see the ini() note).
    ce_mic <- CEFFECT / mic
    gamma <- emax * ce_mic^hill / (ec50^hill + ce_mic^hill)

    # Cordes 2016 equation 3, with beta* = beta0 * d (Materials and
    # Methods). knet is the net rate of change of the burden in
    # log10 CFU/day; its negative is the early bactericidal activity.
    knet <- kgrowth - kimm * f_immune - gamma

    # Integrating the linear state at ln(10) * knet * bact makes
    # log10(bact) move by exactly knet per day, which is the paper's
    # equation written on the scale its parameters are reported in.
    d/dt(bact) <- log(10) * knet * bact
    bact(0) <- 10^log10_cfu0

    # log10 CFU/L observation, with a 1-CFU/L floor so the log10 stays
    # finite if the burden is driven below 1 CFU/L.
    log_cfu <- log10(bact + 1)
    log_cfu ~ add(addSd)
  })
}
