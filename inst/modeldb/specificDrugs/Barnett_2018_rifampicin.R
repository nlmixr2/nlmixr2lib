Barnett_2018_rifampicin <- function() {
  description <- "One-compartment population PK model with Wilkins/Savic transit-compartment absorption for a single 600 mg oral dose of rifampicin in healthy adult males (Barnett 2018), refit from the Wilkins 2008 structural form. The rifampicin model is one of three popPK models developed jointly in Barnett 2018 to support OATP1B drug-drug-interaction modeling with coproporphyrin I and rosuvastatin; the rifampicin compartmental output is the time-varying CRIF input that drives the competitive OATP1B inhibition term in the sibling coproporphyrin I and rosuvastatin models."
  reference <- paste(
    "Barnett S, Ogungbenro K, Menochet K, Shen H, Lai Y, Humphreys WG, Galetin A.",
    "Gaining Mechanistic Insight Into Coproporphyrin I as Endogenous Biomarker",
    "for OATP1B-Mediated Drug-Drug Interactions Using Population Pharmacokinetic",
    "Modeling and Simulation.",
    "Clin Pharmacol Ther. 2018;104(3):564-574.",
    "doi:10.1002/cpt.983.",
    "Structural rifampicin popPK model adapted from",
    "Wilkins JJ et al. Antimicrob Agents Chemother. 2008;52(6):2138-2148;",
    "see modellib('Wilkins_2008_rifampicin') for the original transit-chain form.",
    sep = " "
  )
  vignette <- "Barnett_2018_rifampicin"
  units <- list(time = "hour", dosing = "mg", concentration = "umol/L")

  covariateData <- list(
    OCC = list(
      description        = "Integer-valued occasion / period indicator (Barnett 2018 study design: OCC1 = rifampicin-only period, OCC2 = rosuvastatin-only period, OCC3 = combined rifampicin + rosuvastatin period).",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Time-varying within subject; constant within an occasion. Rifampicin was dosed on OCC1 and OCC3 in the source clinical study (Lai et al. 2016, the n=12 healthy-male SLCO1B1-wildtype cohort that supplied the data for the Barnett 2018 fit). Decomposed inside model() into binary indicators oc1, oc2, oc3 that multiplex the per-occasion IOV etas on log-Ka, log-V, log-MTT (Barnett 2018 Table 1 reports a single shared IOV variance per parameter across occasions; the encoding mirrors the Wilkins_2008_rifampicin OMEGA BLOCK(1) SAME pattern). OCC2 (RSV-only period) has no rifampicin data in the source fit, so its eta is unused for in-paper simulations; it is retained for users who want to simulate at all three study occasions.",
      source_name        = "OCC"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 12L,
    n_studies        = 1L,
    n_observations   = 276L,
    age_range        = "Healthy adult males; demographic distribution not tabulated by Barnett 2018 (source dataset Lai et al. 2016, 12 healthy male subjects, SLCO1B1 c.521 T>C wildtype only; no OATP1B1*5 / *15 carriers).",
    weight_range     = "(not extracted; Barnett 2018 Methods do not tabulate per-subject weights for the n=12 cohort.)",
    sex_female_pct   = 0,
    disease_state    = "Healthy adult male volunteers in a three-occasion (7-day washout between OCC1-2 and OCC2-3) drug-drug-interaction crossover study.",
    dose_range       = "Single 600 mg oral rifampicin dose at the start of OCC1 and at the start of OCC3 (co-administered with 5 mg rosuvastatin on OCC3).",
    regions          = "(not extracted; the underlying clinical study Lai et al. 2016 region was not explicitly stated in Barnett 2018 Methods.)",
    notes            = "Demographics inferred from Barnett 2018 Methods (Clinical data section) and the cited source clinical study Lai Y et al., Pharmacol Res Perspect 2016;4(3):e00207. The 276 rifampicin plasma samples for modeling purposes were collected over the 24 h post-dose window on OCC1 (n=144) and OCC3 (n=132). The NONMEM covariance step failed for this model (Table 1 footnote a), so the parameter standard errors are not available."
  )

  ini({
    # Structural parameters -- Barnett 2018 Table 1, rifampicin rows
    # (column 'Estimate (SE%) Population'). The NONMEM covariance step
    # for this run failed (Table 1 footnote a), so the SE% column is
    # not reported in the source. Point estimates are taken verbatim
    # from Table 1. The transit-chain form follows the Wilkins 2008 /
    # Savic 2007 analytical input rate -- see the model() block below
    # for the rxode2 transit() implementation.
    lka  <- log(2.09);  label("Absorption rate constant Ka (1/h)")                        # Table 1, RIF row 'ka (1/h)' = 2.09
    lcl  <- log(3.97);  label("Apparent oral clearance CL/F (L/h)")                       # Table 1, RIF row 'CL (L/h)' = 3.97
    lvc  <- log(24.7);  label("Apparent central volume of distribution V/F (L)")          # Table 1, RIF row 'V (L)' = 24.7
    lmtt <- log(0.74);  label("Mean transit time MTT through the absorption chain (h)")   # Table 1, RIF row 'MTT (h)' = 0.74
    lnn  <- log(8.63);  label("Number of transit compartments NN (continuous, unitless)") # Table 1, RIF row 'n' = 8.63

    # IIV -- log-normal variances computed from Table 1 IIV CV% column
    # via the standard omega^2 = log(1 + CV^2) conversion.
    etalka  ~ 0.19199   # Table 1, RIF IIV ka  = 46.0% CV; log(1 + 0.46^2)  = 0.19199
    etalcl  ~ 0.08991   # Table 1, RIF IIV CL  = 30.7% CV; log(1 + 0.307^2) = 0.08991
    etalvc  ~ 0.03734   # Table 1, RIF IIV V   = 19.5% CV; log(1 + 0.195^2) = 0.03734
    etalmtt ~ 0.39017   # Table 1, RIF IIV MTT = 69.1% CV; log(1 + 0.691^2) = 0.39017
    etalnn  ~ 0.05731   # Table 1, RIF IIV n   = 24.3% CV; log(1 + 0.243^2) = 0.05731

    # IOV -- Barnett 2018 reports a single shared IOV variance per
    # parameter (ka, V, MTT) across study occasions. Encoded with the
    # first occasion's eta estimated and the others fixed to the same
    # value (matching the Wilkins_2008_rifampicin OMEGA BLOCK(1) SAME
    # idiom; nlmixr2 has no SAME shortcut).
    etaiov_ka_1  ~ 0.21674          # Table 1, RIF IOV ka  = 49.2% CV; log(1 + 0.492^2) = 0.21674
    etaiov_ka_2  ~ fix(0.21674)     # same IOV variance across occasions
    etaiov_ka_3  ~ fix(0.21674)
    etaiov_vc_1  ~ 0.00434          # Table 1, RIF IOV V   = 6.6%  CV; log(1 + 0.066^2) = 0.00434
    etaiov_vc_2  ~ fix(0.00434)
    etaiov_vc_3  ~ fix(0.00434)
    etaiov_mtt_1 ~ 0.20428          # Table 1, RIF IOV MTT = 47.6% CV; log(1 + 0.476^2) = 0.20428
    etaiov_mtt_2 ~ fix(0.20428)
    etaiov_mtt_3 ~ fix(0.20428)

    # Residual error -- combined additive + proportional. The additive
    # component is paper-reported in umol/L and fixed at estimation.
    propSd <- 0.313;       label("Proportional residual error (fraction)")             # Table 1, RIF row 'r prop (%)' = 31.3 (interpreted as 31.3% CV -> fraction 0.313)
    addSd  <- fixed(0.01); label("Additive residual error (umol/L, fixed)")             # Table 1, RIF row 'r add (uM)' = 0.01 FIXED
  })

  model({
    # Decompose the integer-valued OCC column into binary indicators
    # for IOV multiplexing on log-Ka, log-V, log-MTT.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)

    iov_ka  <- oc1 * etaiov_ka_1  + oc2 * etaiov_ka_2  + oc3 * etaiov_ka_3
    iov_vc  <- oc1 * etaiov_vc_1  + oc2 * etaiov_vc_2  + oc3 * etaiov_vc_3
    iov_mtt <- oc1 * etaiov_mtt_1 + oc2 * etaiov_mtt_2 + oc3 * etaiov_mtt_3

    # Individual PK parameters
    ka  <- exp(lka  + etalka  + iov_ka)
    cl  <- exp(lcl  + etalcl)
    vc  <- exp(lvc  + etalvc  + iov_vc)
    mtt <- exp(lmtt + etalmtt + iov_mtt)
    nn  <- exp(lnn  + etalnn)

    kel <- cl / vc

    # 1-compartment PK with the Savic 2007 / Wilkins 2008 analytical
    # transit-chain input rate. rxode2's built-in transit(n, mtt)
    # function delivers the closed-form gamma-PDF input to depot;
    # bioavailability on depot is set to 0 so that the bolus content
    # does not enter and only the transit-chain rate drives absorption.
    d/dt(depot)   <- transit(nn, mtt) - ka * depot
    d/dt(central) <-                    ka * depot - kel * central

    f(depot) <- 0

    # Plasma concentration. Dose in mg, V in L -> central/vc has units
    # of mg/L; convert to umol/L by multiplying by 1000 / MW_rifampicin
    # (rifampicin free-acid MW = 822.94 g/mol).
    MW_rifampicin <- 822.94
    Cc <- (central / vc) * 1000 / MW_rifampicin

    Cc ~ add(addSd) + prop(propSd)
  })
}
