deJonge_2005_paclitaxel <- function() {
  description <- "Semi-mechanistic population pharmacokinetic model for orally administered paclitaxel formulated in Cremophor EL (CrEL) and coadministered with cyclosporin A in adult cancer patients. Free paclitaxel in the gastrointestinal tract (depot) absorbs first-order (kabs) into a two-compartment plasma disposition (central + peripheral1; linear elimination CL/F, volume V/F, intercompartmental clearance Q derived from the paper's k23 = Q/Vc and k32 = Q/Vp). A second GI-tract paclitaxel pool (`bound`) holds drug encapsulated in CrEL micelles; the depot <-> bound equilibrium is governed by a single rate constant keq whose forward binding rate scales with the GI-tract CrEL amount (`cremophor`), which itself decays first-order with rate kcrem. Bioavailability F1 is fixed at 1 with log-normal between-subject variability; the paper found no dose-dependence in F. Inter-occasion variability on CL collapses to between-subject variability in this packaged form because the source dataset's occasion column is not encoded -- see vignette Assumptions and deviations."
  reference <- paste(
    "de Jonge ME, Huitema ADR, Schellens JHM, Rodenhuis S, Beijnen JH. (2005).",
    "Population pharmacokinetics of orally administered paclitaxel formulated in Cremophor EL.",
    "British Journal of Clinical Pharmacology 59(3):325-334.",
    "doi:10.1111/j.1365-2125.2004.02325.x.",
    sep = " "
  )
  vignette <- "deJonge_2005_paclitaxel"
  paper_specific_compartments <- c("bound", "cremophor")

  units <- list(
    time          = "hour",
    dosing        = "mg",
    concentration = "mg/L"
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age in years (median 53, range 26-69).",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested in covariate screen on CL and V; not retained -- the covariate analysis did not reveal significant relationships with any pharmacokinetic parameter (de Jonge 2005 Results, final paragraph; Table 2)."
    ),
    WT = list(
      description = "Body weight in kg (median 71, range 50-103).",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested in covariate screen using the centered power form Cl_pop = theta1 * (WT/71)^theta2 (Covariate analysis section); not retained. The paper never assigns retained allometric or covariate effects. See vignette Assumptions and deviations."
    ),
    BSA = list(
      description = "Body surface area in m^2 (median 1.81, range 1.52-2.22).",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested; not retained. Doses in the source studies were prescribed per BSA (60-360 mg/m^2) but BSA itself did not enter the final structural model."
    ),
    SCR = list(
      description = "Serum creatinine in umol/L (median 82, range 61-127).",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested; not retained. All patients had normal renal function per the inclusion criteria."
    ),
    AST = list(
      description = "Aspartate amino transferase in U/L (median 20, range 6-57).",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested; not retained."
    ),
    ALT = list(
      description = "Alanine amino transferase in U/L (median 23, range 4-86).",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested; not retained."
    ),
    ALP = list(
      description = "Alkaline phosphatase in U/L (median 134, range 36-439).",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested; not retained."
    ),
    GGT = list(
      description = "Gamma-glutamyl transferase in U/L (median 98, range 10-873).",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested; not retained."
    ),
    ALB = list(
      description = "Serum albumin in g/L (median 44, range 29-54).",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested; not retained."
    ),
    BILT = list(
      description = "Total bilirubin in umol/L (median 7, range 4-21).",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested; not retained."
    ),
    LDH = list(
      description = "Lactate dehydrogenase in U/L (median 344, range 116-1129).",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested; not retained."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 55L,
    n_studies      = 4L,
    age_range      = "26-69 years",
    age_median     = "53 years",
    weight_range   = "50-103 kg",
    weight_median  = "71 kg",
    bsa_range      = "1.52-2.22 m^2",
    bsa_median     = "1.81 m^2",
    sex_female_pct = 69,
    disease_state  = "Adult cancer patients receiving oral paclitaxel for advanced solid tumours (six patients had advanced breast cancer continuing on oral paclitaxel after high-dose chemotherapy with cyclophosphamide / thiotepa / carboplatin and peripheral blood progenitor cell transplantation; the remaining 49 patients were from three previously published studies in solid tumours). All patients had normal cardiac, renal, hepatic, haematopoietic and pulmonary function.",
    dose_range     = "Oral paclitaxel 60-360 mg/m^2 once or twice daily; 67 courses; 797 samples. Each oral paclitaxel administration was preceded by oral cyclosporin A 12-15 mg/kg (Neoral solution 10 min before, or capsules 30 min before, paclitaxel) which inhibits intestinal P-gp and CYP3A4 to permit measurable paclitaxel absorption. Patients fasted overnight and received a standard breakfast 2 h post-dose. Patients who vomited within 1 h of paclitaxel were excluded.",
    regions        = "The Netherlands (Netherlands Cancer Institute / Slotervaart Hospital, Amsterdam).",
    notes          = paste(
      "Pooled across four treatment protocols (de Jonge 2005 Table 1):",
      "(1) weekly 200 mg, n=6 (not previously published);",
      "(2) one administration of 60-360 mg/m^2, n=33 (ref 23 -- Meerum Terwogt 1999);",
      "(3) two administrations 60-160 mg/m^2 with 7 h interval, n=10 (ref 24 -- Malingre 2000);",
      "(4) two administrations of 60 mg/m^2 one week apart, with and without added CrEL (5 or 15 mL/m^2), n=6 (ref 21 -- Malingre 2001).",
      "NONMEM V level 1.1 double precision, FOCE-INTER. Goodness-of-fit assessed via Xpose 2.0 in S-Plus 2000.",
      "Cremophor EL is the CrEL formulation vehicle in the i.v. paclitaxel solution Paxene (6 mg paclitaxel per mL of CrEL/ethanol 1:1 w/v); the oral formulation studied here is the i.v. solution swallowed with 100 mL of tap water.",
      "Multi-dose simulations: the structural model assumes the GI-tract CrEL pool `cremophor` is replenished to its normalised baseline at each oral paclitaxel administration.",
      "For single-dose simulations the model file sets `cremophor(0) <- 1` via the model() block.",
      "For multi-dose simulations the user must add explicit dosing events for the `cremophor` compartment with amt = 1 (normalised) at each oral paclitaxel administration time.",
      "Covariate analysis (Table 2): age, body weight, body surface area, serum creatinine, AST, ALT, alkaline phosphatase, GGT, serum albumin, total bilirubin, lactate dehydrogenase were tested on CL and V; none reached significance and none were retained.",
      sep = " "
    )
  )

  ini({
    # ----- Structural PK parameters (de Jonge 2005 Table 3, final estimates) -----
    # All values are paper-reported point estimates (RSE in parentheses) from
    # Table 3 of de Jonge 2005; reparameterised to canonical CL/Vc/Q/Vp from the
    # paper's k_abs / k_23 / k_32 / CL / V parameterisation: Q = k23 * Vc and
    # Vp = Q / k32. The reparameterisation is exact at the typical-value level
    # and preserves the variance of the IIV on k23 when translated onto Q
    # (see "Inter-individual variability" below).

    lka <- log(0.62)
    label("Absorption rate constant kabs from free GI-tract pool into central (1/h)") # Table 3: k_abs = 0.62 (RSE 37.7%)

    lcl <- log(127)
    label("Apparent clearance CL/F from central (L/h)") # Table 3: CL = 127 (RSE 9.61%); paper's notation is CL/F with F = 1 FIX

    lvc <- log(409)
    label("Apparent central volume of distribution V/F (L)") # Table 3: V = 409 (RSE 16.8%); paper's notation is V/F with F = 1 FIX

    lq <- log(98.569)
    label("Intercompartmental clearance Q = k23 * Vc (L/h)") # Derived from Table 3: k23 = 0.241 (RSE 17.5%) and Vc = 409 -> Q = 0.241 * 409

    lvp <- log(1041.96)
    label("Peripheral volume of distribution Vp = Q / k32 (L)") # Derived from Table 3: k32 = 0.0946 (RSE 12.9%) and Q = 98.569 -> Vp = 98.569 / 0.0946

    lkcrem <- log(1.73)
    label("First-order rate constant for CrEL elimination from the GI-tract pool kcrem (1/h)") # Table 3: k_crem = 1.73 (RSE 12.4%)

    lkeq <- log(0.334)
    label("Equilibrium rate constant for paclitaxel free <-> bound in the GI tract keq (1/h); paper reports keq = k15/k51, encoded here as k_forward = keq * cremophor and k_backward = keq") # Table 3: k_eq = 0.334 (RSE 43.4%)

    lfdepot <- fixed(log(1))
    label("Oral bioavailability anchor F1 (FIXED at 1; paper found F was independent of paclitaxel dose and of CrEL amount over the studied range)") # Table 3: Bioavailability = 1 FIX

    # ----- Between-subject variability (de Jonge 2005 Table 3, % IIV column) -----
    # Paper reports CV% on an exponential (log-normal) error model; converted to
    # eta variance via omega^2 = log(1 + CV^2). Diagonal block; the paper does
    # not report any cross-parameter correlation.

    # k_abs IIV 61.1% (RSE 33.5%) -> log(1 + 0.611^2) = 0.31712
    etalka ~ 0.31712

    # k_23 IIV 31.3% (RSE 66.5%) translated onto Q (canonical CL/Vc/Q/Vp form):
    # log(k_23) = log(Q) - log(Vc); with no IIV on Vc the IIV on log(k_23)
    # equals the IIV on log(Q), so omega^2_Q = log(1 + 0.313^2) = 0.09349.
    # The reparameterisation moves the IIV onto Q without changing its variance.
    etalq ~ 0.09349

    # F1 IIV 55.0% (RSE 22.5%) -> log(1 + 0.550^2) = 0.26429.
    # F1 is FIXED at 1 (typical value) so this IIV is the between-subject spread
    # in bioavailability around the population anchor of 1; the standard NONMEM
    # idiom F = 1 FIX with eta on F = 0 + omega is reproduced here as
    # f_depot_i = exp(lfdepot + etalfdepot) = exp(etalfdepot) in model().
    etalfdepot ~ 0.26429

    # CL inter-occasion variability 31.8% (RSE 38.3%) encoded as ordinary IIV.
    # nlmixr2lib has no canonical pattern for occasion-keyed IOV without an
    # OCC column, so the IOV variance is carried as a between-subject random
    # effect on the clearance log-typical-value. The packaged model therefore
    # loses the within-subject between-occasion drift the paper estimates but
    # preserves the population-level CL spread. See vignette Assumptions and
    # deviations. omega^2 = log(1 + 0.318^2) = 0.09633.
    etalcl ~ 0.09633

    # ----- Residual error (de Jonge 2005 Table 3, last row) -----
    # Table 3: residual proportional error = 45.1% (RSE 10.2%).
    # Proportional CV in % maps directly to propSd in nlmixr2's prop() error
    # form: propSd = CV/100 = 0.451.
    propSd <- 0.451
    label("Proportional residual error on plasma paclitaxel concentration (fraction)") # Table 3: residual proportional error = 45.1%
  })

  model({
    # ----- Individual parameters -----
    # Log-normal typical values + log-normal IIV.
    ka     <- exp(lka + etalka)
    cl     <- exp(lcl + etalcl)
    vc     <- exp(lvc)
    q      <- exp(lq + etalq)
    vp     <- exp(lvp)
    kcrem  <- exp(lkcrem)
    keq    <- exp(lkeq)
    fdepot <- exp(lfdepot + etalfdepot)

    # ----- Two-compartment central / peripheral micro-constants -----
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ----- Mechanistic GI-tract CrEL <-> paclitaxel binding -----
    # The paper's keq is the equilibrium rate constant for free <-> bound
    # paclitaxel in the GI tract, written as keq = k15/k51 with the
    # rapid-equilibrium narrative "the amount of CrEL in compartment 4
    # influences the equilibrium between unbound and bound paclitaxel"
    # (Results, Model development paragraph 2). The published value
    # 0.334 (1/h) is one rate constant; the paper does not write out the
    # explicit functional form of the Acrem-modulation. The most
    # parsimonious mass-action mechanism that reproduces the narrative is:
    #
    #     forward (free -> bound):  k15 * depot = keq * cremophor * depot
    #     backward (bound -> free): k51 * bound = keq * bound
    #
    # At rapid equilibrium this gives bound/depot = cremophor; when the
    # normalised cremophor pool is at its baseline (1) about half the
    # GI-tract paclitaxel is bound, and as cremophor decays with kcrem
    # the equilibrium shifts to free, releasing paclitaxel for absorption.
    # This matches the paper's narrative and is faithful to the rate-constant
    # value in Table 3. See vignette Assumptions and deviations for a fuller
    # discussion of this mechanistic choice in the absence of a control
    # stream.
    bind_fwd <- keq * cremophor * depot
    bind_bwd <- keq * bound

    # ----- ODE system (de Jonge 2005 Figure 2) -----
    d/dt(depot)       <- -ka * depot - bind_fwd + bind_bwd
    d/dt(bound)       <-  bind_fwd - bind_bwd
    d/dt(cremophor)   <- -kcrem * cremophor
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Single-dose initial condition for the GI-tract CrEL pool. The pool is
    # normalised so that one oral paclitaxel administration is paired with a
    # unit amount of CrEL (the paper found bioavailability and the absorption
    # dynamics were independent of the absolute CrEL amount over the studied
    # 5-15 mL/m^2 range). For multi-dose simulations the user adds explicit
    # dosing events into the `cremophor` compartment with amt = 1 at each
    # oral paclitaxel administration time so each dose is paired with a
    # fresh CrEL replenishment.
    cremophor(0) <- 1

    # ----- Bioavailability -----
    f(depot) <- fdepot

    # ----- Observation and error model -----
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
