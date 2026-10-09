Lin_2022_PF05251749_750mg <- function() {
  description <- paste(
    "Two-compartment oral pharmacokinetic reduction of the Simcyp",
    "minimal-PBPK-with-single-adjusting-compartment (SAC) base model for",
    "the casein kinase 1 delta/epsilon inhibitor PF-05251749 at 750 mg",
    "q.d. in healthy adults (Lin 2022). The authors built separate",
    "dose-specific Simcyp models at 400 and 750 mg q.d. because exposure",
    "was less than dose-proportional above 200 mg; this file is the",
    "750 mg model and Lin_2022_PF05251749_400mg is its sibling. The",
    "Simcyp whole-body equations are not published, but the compound",
    "layer (Table 2: ka, fa, fu, B:P, observed oral clearance CLpo, Vss,",
    "Vsac and the SAC inter-compartmental clearance Q) is enough to",
    "rebuild the disposition as first-order absorption into a systemic",
    "compartment exchanging with the SAC (peripheral1). Systemic",
    "clearance and bioavailability follow from CLpo through the",
    "well-stirred liver model; the systemic volume is",
    "(Vss - Vsac) x 70 kg minus a liver volume. Hepatic blood flow",
    "(90 L/h) and liver volume (1.648 L) are not printed in the paper and",
    "are standard Simcyp healthy-adult values. Nothing is fitted. The",
    "reduction reproduces the paper's predicted day-14 Tmax to within 3%",
    "and gives Cmax and AUCtau 7% below the paper's geometric-mean",
    "predictions (see the validation vignette).",
    "Typical-value model: the source reports no variance components and",
    "no residual-error model, so there are no etas and propSd is zero.",
    "The CYP3A-induction (midazolam DDI) layer depends on the",
    "proprietary Simcyp midazolam compound file and enzyme-turnover",
    "defaults and is not encoded.",
    sep = " "
  )
  reference <- paste(
    "Lin J, Gaudreault F, Johnson N, Lin Z, Nouri P, Goosen TC,",
    "Sawant-Basak A. (2022). Investigation of CYP3A induction by",
    "PF-05251749 in early clinical development: comparison of linear",
    "slope physiologically based pharmacokinetic prediction and",
    "biomarker response. Clin Transl Sci 15(9):2184-2194.",
    "doi:10.1111/cts.13352.",
    sep = " "
  )
  vignette <- "Lin_2022_PF05251749"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(
      analyte = "PF-05251749",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "PF-05251749",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "PF-05251749",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    )
  )

  # The Simcyp simulations draw body weight per virtual subject, and the
  # source states its volumes in L/kg, so the platform model scales Vss and
  # Vsac with weight. That relationship is not a reported covariate model:
  # the reduction fixes a 70 kg reference instead.
  covariatesDataExcluded <- list(
    WT = list(
      description = paste(
        "Body weight. Lin 2022 Table 2 states Vss and Vsac in L/kg, so",
        "the Simcyp model scales distribution volume with each virtual",
        "subject's weight. This reduction evaluates the volumes at a",
        "70 kg reference instead of carrying a weight term, because the",
        "systemic-compartment volume also subtracts a liver volume whose",
        "weight scaling Simcyp does not publish."
      ),
      units = "kg",
      type = "continuous",
      notes = "Implicit in the L/kg volume inputs; not carried in this reduction."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 8L,
    n_studies = 1L,
    age_range = "18-55 years (B8001002 eligibility); Simcyp virtual population 20-50 years",
    age_mean = "42.5 years (SD 10.5) across all 61 Part A participants",
    weight_mean = "76.4 kg (SD 9.5) across all 61 Part A participants; volumes evaluated at a 70 kg reference",
    sex_female_pct = 16.4,
    race_ethnicity = c(White = 85.2, Black = 14.8),
    disease_state = "Healthy adults (male, or female of non-childbearing potential)",
    dose_range = "750 mg PF-05251749 orally q.d. for 14 days, fasted",
    regions = "United States (single clinical research unit)",
    studies = paste(
      "B8001002 Part A (NCT02691702): randomized, double-blind,",
      "sequential, placebo-controlled multiple-ascending-dose study of",
      "PF-05251749 50-750 mg q.d. for 14 days. The 750 mg cohort (n = 8",
      "on PF-05251749) supplied the observed CL/F of 22.5 L/h that is",
      "the Simcyp clearance input, and the day-14 profile Vsac and Q",
      "were fitted to."
    ),
    notes = paste(
      "Demographics (Table S3) are for all 61 Part A participants",
      "(PF-05251749, placebo and melatonin arms). The Simcyp",
      "simulations used a virtual healthy-volunteer population of 10",
      "trials x 10 subjects, aged 20-50 years, 50% female, fasted. This",
      "is a PBPK analysis rather than a population-PK fit, so there is",
      "no pooled analysis dataset and no estimated variance components.",
      "Elimination is mainly metabolic (in vitro CYP1A2 34%, CYP3A4/5",
      "26%, CYP2B6 22%, CYP2C19 18%; < 0.2% of the dose excreted",
      "unchanged in urine); the Simcyp model used the in vivo clearance",
      "without an enzyme-pathway split."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Every parameter is held constant: nothing is estimated in this
    # reduction. Values are Lin 2022 Table 2 inputs (750 mg column where
    # the table gives 400 mg / 750 mg pairs) or arithmetic consequences of
    # them. Shared quantities:
    #   fu,p  = 0.231     Table 2, fraction unbound in plasma
    #   B:P   = 0.7       Table 2, blood-to-plasma ratio
    #   fa    = 1         Table 2, fraction absorbed (assumed by the authors)
    #   CLpo  = 22.5 L/h  Table 2 'CL po (L/h)' 29.3/22.5, clinical data
    #                     (= B8001002 Table S5 steady-state CL/F at 750 mg)
    #   QH    = 90 L/h    hepatic blood flow; not printed in Lin 2022.
    #                     Standard Simcyp healthy-adult value, as printed
    #                     in Hanley 2024 (doi:10.1002/psp4.13106)
    #                     Appendix S1 Equation S2.
    # ------------------------------------------------------------------

    lka <- fixed(log(2))
    label("First-order absorption rate constant (1/h)")
    # Table 2 'K a (1/h)' = 2, chosen by sensitivity analysis (Figure S4)
    # to recover the observed Tmax.

    # Well-stirred liver with fa = 1 and fG = 1. For an oral dose the
    # apparent plasma clearance is fu,p * CLint, so
    #   fu,b * CLint = CLpo / (B:P) = 22.5 / 0.7      = 32.142857 L/h
    #   fH = QH / (QH + fu,b * CLint) = 90 / 122.142857 = 0.736842
    #   CL (systemic, plasma) = fH * CLpo = 0.736842 * 22.5 = 16.578947 L/h
    lcl <- fixed(log(16.578947))
    label("Systemic plasma clearance (L/h)")
    # Derived from Table 2 (CLpo 22.5 L/h, B:P 0.7, fa 1) and QH 90 L/h.

    # Systemic-compartment volume at the 70 kg reference weight:
    #   (Vss - Vsac) * 70 = (3.1 - 1.3) L/kg * 70 kg = 126 L
    #   minus the liver   = 126 - 1.648            = 124.352 L
    # Vss 3.1 L/kg and Vsac 1.3 L/kg are Table 2 ('V ss (L/kg)' 3.1;
    # 'V SAC (L/kg)' 1.6/1.3). The 1.648 L liver volume is not printed in
    # Lin 2022; it is the Simcyp default liver weight of 1648 g printed in
    # Hanley 2024 Appendix S1 Equation S2.
    lvc <- fixed(log(124.352))
    label("Systemic compartment volume at the 70 kg reference weight (L)")

    # Single adjusting compartment: Vsac 1.3 L/kg * 70 kg = 91 L.
    lvp <- fixed(log(91))
    label("Single adjusting compartment volume at the 70 kg reference weight (L)")
    # Table 2 'V SAC (L/kg)' = 1.3 (750 mg), fitted to the day-14 profile.

    lq <- fixed(log(11))
    label("Inter-compartmental clearance between systemic compartment and SAC (L/h)")
    # Table 2 'Q (L/h)' = 11 (750 mg), fitted to the day-14 profile.

    # Oral bioavailability F = fa * fG * fH = 1 * 1 * 0.736842. fa = 1 is
    # Table 2; fG = 1 because the in vivo clearance input carries no gut
    # intrinsic clearance (Table 2 'F u, gut' = 1, no enzyme pathways).
    lfdepot <- fixed(log(0.736842))
    label("Oral bioavailability (fraction)")
    # Derived from Table 2 (fa 1, CLpo 22.5 L/h, B:P 0.7) and QH 90 L/h.

    # Lin 2022 is a PBPK simulation analysis: the %CV values in Table 3 are
    # the spread of the Simcyp virtual population, not estimated variances,
    # and no residual-error model is reported. Residual error is zero.
    propSd <- fixed(0)
    label("Proportional residual error (fraction)")
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    vp <- exp(lvp)
    q <- exp(lq)
    fdepot <- exp(lfdepot)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Lin 2022 Table 2 'Model: Minimal PBPK' with first-order absorption;
    # the liver and portal-vein compartments are lumped into the systemic
    # compartment and the SAC is peripheral1.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fdepot

    # Doses in mg and volumes in L give mg/L; x1000 reports ng/mL as in
    # Lin 2022 Table 3.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
