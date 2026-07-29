Varatharajan_2016_daunorubicin <- function() {
  description <- paste(
    "Population PK model for IV daunorubicin (Dnr) and its primary",
    "carbonyl-reductase metabolite daunorubicinol (DOL) in adult de novo",
    "acute myeloid leukaemia (AML) patients (Varatharajan 2016). Each",
    "component (parent and metabolite) is described by an independent",
    "two-compartment disposition parameterised on apparent clearance,",
    "central volume, and the inter-compartmental rate constants K12 and",
    "K21. Daunorubicin is converted to daunorubicinol via parent",
    "elimination (the model assumes the fraction metabolised fm = 1, so",
    "the published DOL CL and V are 'apparent' values that absorb fm).",
    "No covariates were retained in the final structural model;",
    "demographic / pharmacogenetic associations in the paper are reported",
    "on post hoc empirical-Bayes estimates rather than as fixed-effects",
    "covariate parameters."
  )
  reference <- "Varatharajan S, Panetta JC, Abraham A, Karathedath S, Mohanan E, Lakshmi KM, Arthur N, Srivastava VM, Nemani S, George B, Srivastava A, Mathews V, Balasubramanian P. Population pharmacokinetics of Daunorubicin in adult patients with acute myeloid leukemia. Cancer Chemother Pharmacol. 2016;78(5):1051-1058. doi:10.1007/s00280-016-3166-8"
  vignette <- "Varatharajan_2016_daunorubicin"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ng/mL"
  )

  covariateData <- list()

  population <- list(
    n_subjects     = 70,
    n_studies      = 1,
    age_range      = "16-60 years (median 38)",
    age_median     = "38 years",
    sex_female_pct = 47.1,
    species        = "Human (adult AML)",
    disease_state  = "Adult de novo acute myeloid leukaemia (excluding AML-M3) treated with standard induction chemotherapy comprising cytarabine + daunorubicin.",
    dose_range     = "Daunorubicin 60 mg/m^2/day as 1-hour IV infusion on days 1-3 (PK sampling on day 1 only).",
    regions        = "Christian Medical College, Vellore, Tamil Nadu, India.",
    enrollment_period = "2009-2014",
    cytogenetic_risk_pct = c(Favorable = 12, Intermediate = 68, Adverse = 20),
    notes          = paste(
      "Patient demographics from Table 1. Plasma sampling on day 1 of",
      "induction at 0, 0.25, 1, 2, 4, 6, and 24 h. Plasma Dnr and DOL were",
      "quantified in all 70 patients by HPLC with fluorescence detection",
      "(LOD 1 ng/mL, LOQ 10 ng/mL). Sex split per Table 1 is reported as",
      "42 males / 33 females (sums to 75) although the cohort size is",
      "consistently stated as n = 70 throughout the abstract, methods, and",
      "Table 1 header; sex_female_pct here uses 33 / 70 = 47.1% to match",
      "the n = 70 denominator. Population is a single-centre Indian AML",
      "cohort; the published PK parameters carry that ethnic-geographic",
      "context."
    )
  )

  ini({
    # Daunorubicin (parent) -- two-compartment disposition.
    # Population means and IIV CV% from Table 2 (Monolix SAEM final
    # estimates). The paper estimates K12 and K21 directly (Monolix
    # micro-constant parameterisation) rather than Q and Vp; the
    # paper-reported values are kept here so the source-trace is
    # one-to-one with Table 2.
    lcl       <- log(269.8)            ; label("Daunorubicin total clearance CL (L/h)")                      # Table 2 row Dnr Clearance
    lvc       <- log(15.0)             ; label("Daunorubicin central volume V (L)")                          # Table 2 row Dnr Volume
    lk12      <- log(22.4)             ; label("Daunorubicin central->peripheral rate constant K12 (1/h)")   # Table 2 row Dnr K12
    lk21      <- log(0.6)              ; label("Daunorubicin peripheral->central rate constant K21 (1/h)")   # Table 2 row Dnr K21

    # Daunorubicinol (metabolite, DOL) -- two-compartment disposition.
    # CL and V are apparent values (CL_DOL,real / fm and V_DOL,real / fm
    # respectively) because the fraction of Dnr metabolised to DOL
    # (fm) is not separately identifiable from plasma Dnr / DOL data;
    # the metabolite ODE below assumes fm = 1 so these apparent values
    # reproduce the observed DOL concentrations.
    lcl_dol   <- log(23.6)             ; label("Daunorubicinol apparent clearance CL/fm (L/h)")              # Table 2 row DOL Clearance
    lvc_dol   <- log(7.6)              ; label("Daunorubicinol apparent central volume V/fm (L)")            # Table 2 row DOL Volume
    lk12_dol  <- log(32.0)             ; label("Daunorubicinol central->peripheral rate constant K12 (1/h)") # Table 2 row DOL K12
    lk21_dol  <- log(0.4)              ; label("Daunorubicinol peripheral->central rate constant K21 (1/h)") # Table 2 row DOL K21

    # IIV. Diagonal omega blocks (paper does not report off-diagonal
    # correlations). CV% from Table 2 converted to log-scale variance
    # via omega^2 = log(1 + CV^2).
    etalcl       ~ 0.5341  # Dnr CL CV% = 84% -> log(1 + 0.84^2)
    etalvc       ~ 1.1135  # Dnr V  CV% = 143% -> log(1 + 1.43^2)
    etalk12      ~ 0.3895  # Dnr K12 CV% = 69% -> log(1 + 0.69^2)
    etalk21      ~ 0.4178  # Dnr K21 CV% = 72% -> log(1 + 0.72^2)
    etalcl_dol   ~ 0.3434  # DOL CL CV% = 64% -> log(1 + 0.64^2)
    etalvc_dol   ~ 0.3525  # DOL V  CV% = 65% -> log(1 + 0.65^2)
    etalk12_dol  ~ 0.4273  # DOL K12 CV% = 73% -> log(1 + 0.73^2)
    etalk21_dol  ~ 0.3801  # DOL K21 CV% = 68% -> log(1 + 0.68^2)

    # Residual error -- proportional only on each observed analyte
    # ("relative error CV %" Table 2 footer; Monolix proportional
    # residual model assumed normally distributed residuals on the
    # parameter scale).
    propSd     <- 0.47  ; label("Daunorubicin proportional residual SD (fraction)")     # Table 2 Dnr relative error 47%
    propSd_dol <- 0.34  ; label("Daunorubicinol proportional residual SD (fraction)")   # Table 2 DOL relative error 34%
  })

  model({
    # Individual parameters.
    cl       <- exp(lcl      + etalcl)
    vc       <- exp(lvc      + etalvc)
    k12      <- exp(lk12     + etalk12)
    k21      <- exp(lk21     + etalk21)

    cl_dol   <- exp(lcl_dol  + etalcl_dol)
    vc_dol   <- exp(lvc_dol  + etalvc_dol)
    k12_dol  <- exp(lk12_dol + etalk12_dol)
    k21_dol  <- exp(lk21_dol + etalk21_dol)

    # Elimination rate constants (kel = CL / Vc); kel_dol uses the
    # apparent DOL CL and Vc (the ratio is the real DOL elimination
    # rate constant; fm cancels in the quotient).
    kel     <- cl     / vc
    kel_dol <- cl_dol / vc_dol

    # Daunorubicin: 2-compartment IV. Parent elimination drives DOL
    # formation via kel * central (mass equivalent of the parent
    # eliminated; assumes fm = 1 for the structural model). Monolix
    # supplemental Fig 1 schematic.
    d/dt(central)         <- -kel * central -
                              k12 * central + k21 * peripheral1
    d/dt(peripheral1)     <-  k12 * central - k21 * peripheral1

    # Daunorubicinol: 2-compartment with first-order input from parent
    # central elimination.
    d/dt(central_dol)     <-  kel * central -
                              kel_dol * central_dol -
                              k12_dol * central_dol +
                              k21_dol * peripheral1_dol
    d/dt(peripheral1_dol) <-  k12_dol * central_dol -
                              k21_dol * peripheral1_dol

    # Plasma concentrations. Compartment amounts are in mg (matching
    # the dose unit), volumes in L -> mg/L; multiplied by 1000 to
    # report in the paper's unit ng/mL (= ug/L).
    Cc     <- (central     / vc)     * 1000
    Cc_dol <- (central_dol / vc_dol) * 1000

    Cc     ~ prop(propSd)
    Cc_dol ~ prop(propSd_dol)
  })
}
