Weatherley_2018_fosdagrocorat <- function() {
  description <- paste(
    "Simultaneous parent-metabolite population PK model for oral",
    "fosdagrocorat (PF-04171327, a phosphate-ester prodrug of the",
    "dissociated glucocorticoid receptor agonist PF-00251802) in adult",
    "rheumatoid-arthritis patients receiving stable background",
    "methotrexate (Weatherley 2018). The prodrug is fully cleaved by",
    "alkaline phosphatase in the gut wall before absorption; only the",
    "active Metabolite-1 (PF-00251802) and its circulating N-oxide",
    "Metabolite-2 (PF-04015475) are modelled. Metabolite-1 is described",
    "by a two-compartment disposition (apparent CL, V2, Q, V4 fixed at",
    "209 L) with first-order absorption (K12) and bioavailability F1",
    "fixed at 1 (apparent F absorbed via the prodrug-to-Metabolite-1",
    "conversion). Metabolite-2 is described by a one-compartment",
    "disposition (apparent Vm, CLm) with Fm fixed at 1 (assumed 100",
    "percent molar conversion of Metabolite-1 to Metabolite-2).",
    "Standard allometric weight scaling is fixed on Metabolite-1",
    "disposition (exponent 0.75 on CL and Q; exponent 1.00 on V2 and",
    "V4); body-weight scaling on Metabolite-2 CLm is estimated as a",
    "power-form covariate (exponent 0.450) and no weight effect is",
    "applied to Vm (rejected in stepwise covariate testing). Retained",
    "covariates on Metabolite-1 CL are female-vs-male (-26.8 percent)",
    "and a small linear age effect (-0.00633 L/h per year above 40).",
    "Retained covariates on Metabolite-2 CLm are female-vs-male",
    "(-34.1 percent) and body weight. Inter-individual variability is",
    "reported on Metabolite-1 CL (33 percent CV) and absorption rate",
    "K12 (249 percent CV), and on Metabolite-2 Vm (44 percent CV) and",
    "CLm (26 percent CV). The publication's interoccasion variability",
    "on Metabolite-1 F1 (23.8 percent CV across dosing occasions) is",
    "encoded here as IIV on F1 because nlmixr2lib simulation does not",
    "carry an occasion column; this approximation is documented in the",
    "vignette Assumptions section. Residual error is combined additive",
    "plus proportional on the linear-concentration scale separately",
    "for each analyte (Metabolite-1 proportional 19.9 percent +",
    "additive 0.305 ng/mL; Metabolite-2 proportional 7.8 percent +",
    "additive 0.10 ng/mL).")
  reference <- "Weatherley B, McFadyen L, Tammara B. Population pharmacokinetics of fosdagrocorat (PF-04171327), a dissociated glucocorticoid receptor agonist, in patients with rheumatoid arthritis. Clin Transl Sci. 2018;11(1):54-62. doi:10.1111/cts.12515"
  vignette <- "Weatherley_2018_fosdagrocorat"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Total body weight at baseline.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Fixed allometric scaling on Metabolite-1 disposition: exponent 0.75 on CL and Q, exponent 1.00 on V2 and V4. Estimated power exponent 0.450 on Metabolite-2 CLm (Weatherley 2018 Table 2). No weight effect on Metabolite-2 Vm (rejected in stepwise testing because the estimated exponent 0.241 was far from the canonical allometric value 1.0; Weatherley 2018 Results page 58). Reference weight 70 kg per Weatherley 2018 Methods.",
      source_name        = "BWT (body weight)"
    ),
    AGE = list(
      description        = "Subject age at baseline.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear-additive effect on Metabolite-1 CL: change of -0.00633 L/h per year above the reference age of 40 (Weatherley 2018 Table 2 'Without BLQ or taper' column). The publication characterizes this effect as approximately 3 percent CL reduction per 30-year increment (Discussion page 60) and notes it could likely have been dropped without clinical relevance.",
      source_name        = "Age"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male).",
      notes              = "Linear multiplicative effects: -26.8 percent on Metabolite-1 CL and -34.1 percent on Metabolite-2 CLm (Weatherley 2018 Table 2). After accounting for allometric body-weight scaling, female subjects therefore have higher exposure to both Metabolite-1 and Metabolite-2 than male subjects of the same weight. The mechanistic basis of the female effect was not identified; in vitro CYP3A4 metabolism does not show a documented sex difference (Weatherley 2018 Discussion page 60).",
      source_name        = "Sex (1 = male, 2 = female in NONMEM dataset; canonical SEXF = as.integer(SEX == 2))"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 179L,
    n_studies      = 1L,
    age_range      = "18-84 years (combined male and female)",
    weight_range   = "~40-140 kg (paper Figure 5 illustrates exposure across this range; per-arm SDs reported in Table 1)",
    sex_female_pct = 76,
    race_ethnicity = "Reported but not stratified in the public summary; covariate testing on race did not retain race in the final model.",
    disease_state  = "Patients with moderate-to-severe rheumatoid arthritis (RA) receiving stable background methotrexate, enrolled in study A9391010 (NCT01393639), a 12-week phase II randomized double-blind dose-ranging study.",
    dose_range     = "Oral fosdagrocorat 1, 5, 10, or 15 mg once daily for 8 weeks, followed by a blinded 1 mg tapered regimen (Q48 h for 2 weeks, then Q72 h for 2 weeks) and a 1-week off-drug washout. Only the 8-week active-dosing period contributed to the final model after below-quantitation and taper observations were excluded.",
    regions        = "Multinational; specific country composition not summarized in the public publication.",
    notes          = "Demographics (Weatherley 2018 Table 1, per-arm means and SDs): 76 percent female overall. Mean weight 80.4 kg (male) and 72.9 kg (female). Mean baseline creatinine clearance ~100-130 mL/min. NONMEM 7.2 (ICON Development Solutions, Ellicott City, MD) with first-order conditional estimation was used; observations below the limit of quantitation (1.00 ng/mL for Metabolite-1 and 0.50 ng/mL for Metabolite-2) and concentrations measured during the taper period were excluded from the final model. The final model was selected from the 'Without BLQ or taper' simultaneous Metabolite-1 / Metabolite-2 fit (Table 2 right-most column)."
  )

  ini({
    # ------------------------------------------------------------------
    # METABOLITE-1 (PF-00251802) STRUCTURAL PARAMETERS
    # ------------------------------------------------------------------
    # Reference subject: male, 70 kg, 40 years old. All structural
    # values from Weatherley 2018 Table 2, 'Without BLQ or taper' column
    # (the selected final model).

    lka <- log(0.377)
    label("Absorption rate constant K12 from prodrug depot to Metabolite-1 central (1/h)")  # Weatherley 2018 Table 2: K12 = 0.377 1/h

    lcl <- log(7.29)
    label("Apparent Metabolite-1 clearance CL/F at male / 70 kg / 40 yr reference (L/h)")  # Weatherley 2018 Table 2: CL = 7.29 L/h

    lvc <- log(85.9)
    label("Apparent Metabolite-1 central volume V2/F at 70 kg reference (L)")  # Weatherley 2018 Table 2: V2 = 85.9 L

    lq <- log(11.6)
    label("Apparent Metabolite-1 inter-compartmental clearance Q/F at 70 kg reference (L/h)")  # Weatherley 2018 Table 2: Q = 11.6 L/h

    lvp <- fixed(log(209))
    label("Apparent Metabolite-1 peripheral volume V4/F at 70 kg reference (L; fixed)")  # Weatherley 2018 Table 2 + Results page 58: V4 was re-estimated to 209 L, then re-fixed at this value because re-estimation produced only a 0.464 OFV drop

    lfdepot <- fixed(log(1))
    label("Metabolite-1 fraction-of-dose absorbed F1 (fixed at 1; apparent F absorbed via prodrug conversion)")  # Weatherley 2018 Methods + Table 2: F1 FIXED = 1; all parameters are apparent values F-scaled by the prodrug-to-Metabolite-1 conversion

    # ------------------------------------------------------------------
    # METABOLITE-1 COVARIATE EFFECTS
    # ------------------------------------------------------------------

    e_wt_cl <- fixed(0.75)
    label("Allometric exponent for body weight on Metabolite-1 CL (unitless; fixed)")  # Weatherley 2018 Results page 58: 'fixed standard allometric scaling ... power 0.75 for CL and Q'

    e_wt_vc <- fixed(1.00)
    label("Allometric exponent for body weight on Metabolite-1 V2 (unitless; fixed)")  # Weatherley 2018 Results page 58: 'fixed standard allometric scaling ... power 1.0 for V2 and V4'

    e_wt_q <- fixed(0.75)
    label("Allometric exponent for body weight on Metabolite-1 Q (unitless; fixed)")  # Weatherley 2018 Results page 58

    e_wt_vp <- fixed(1.00)
    label("Allometric exponent for body weight on Metabolite-1 V4 (unitless; fixed)")  # Weatherley 2018 Results page 58

    e_sexf_cl <- -0.268
    label("Linear coefficient for female-vs-male on Metabolite-1 CL (unitless; -26.8% change)")  # Weatherley 2018 Table 2: -32.9 / -26.8 / -31.2 / -26.8% across the four data subsets; final 'Without BLQ or taper' column uses -26.8

    e_age_cl <- -0.00633
    label("Additive change in Metabolite-1 CL per year above the reference age of 40 (L/h per year)")  # Weatherley 2018 Table 2: -0.00627 / -0.00589 / -0.00675 / -0.00633 across subsets; final column = -0.00633

    # ------------------------------------------------------------------
    # METABOLITE-1 INTER-INDIVIDUAL VARIABILITY
    # ------------------------------------------------------------------
    # Table 2 reports IIV as approximate lognormal CV percent. Variance
    # used by ini() is the log-scale variance omega^2 = log(CV^2 + 1).
    # Per the source: 'IIV could only be estimated on CL and K12' for
    # Metabolite-1 (Results page 58); no IIV on V2, V4, Q. The
    # interoccasion variability on F1 (23.8% CV across dosing occasions
    # within a subject) is encoded here as IIV on F1 because nlmixr2lib
    # simulation does not carry an occasion column; see vignette
    # Assumptions for the approximation rationale.

    etalcl ~ 0.1034
    # omega^2 = log(0.33^2 + 1) = 0.1034; Weatherley 2018 Table 2 'Without BLQ or taper' column: IIV CL = 33%

    etalka ~ 1.974
    # omega^2 = log(2.49^2 + 1) = 1.974; Weatherley 2018 Table 2 'Without BLQ or taper': IIV K12 = 249%

    etalfdepot ~ 0.0550
    # omega^2 = log(0.238^2 + 1) = 0.0550; Weatherley 2018 Table 2: IOV on F1 = 23.8% (modelled as IIV; see vignette Assumptions)

    # ------------------------------------------------------------------
    # METABOLITE-1 RESIDUAL ERROR
    # ------------------------------------------------------------------

    propSd <- 0.199
    label("Metabolite-1 proportional residual SD on linear concentration (fraction)")  # Weatherley 2018 Table 2 'Without BLQ or taper': proportional = 19.9%

    addSd <- 0.305
    label("Metabolite-1 additive residual SD on linear concentration (ng/mL)")  # Weatherley 2018 Table 2 'Without BLQ or taper': additive = 0.305 ng/mL (note: paper reports RSE 172%, see vignette Errata)

    # ------------------------------------------------------------------
    # METABOLITE-2 (PF-04015475 N-oxide) STRUCTURAL PARAMETERS
    # ------------------------------------------------------------------

    lcl_noxide <- log(17.2)
    label("Apparent Metabolite-2 clearance CLm/F/Fm at male / 70 kg reference (L/h)")  # Weatherley 2018 Table 2: CLm = 17.2 L/h

    lvc_noxide <- log(62.8)
    label("Apparent Metabolite-2 volume of distribution Vm/F/Fm (L)")  # Weatherley 2018 Table 2: Vm = 62.8 L

    # ------------------------------------------------------------------
    # METABOLITE-2 COVARIATE EFFECTS
    # ------------------------------------------------------------------

    e_wt_cl_noxide <- 0.450
    label("Estimated power exponent for body weight on Metabolite-2 CLm (unitless)")  # Weatherley 2018 Table 2 'Without BLQ or taper': BWT power on CLm = 0.450; Results page 59 explicitly notes 'standard allometric scaling was not utilized for Metabolite-2'

    e_sexf_cl_noxide <- -0.341
    label("Linear coefficient for female-vs-male on Metabolite-2 CLm (unitless; -34.1% change)")  # Weatherley 2018 Table 2 'Without BLQ or taper': -34.1%

    # ------------------------------------------------------------------
    # METABOLITE-2 INTER-INDIVIDUAL VARIABILITY
    # ------------------------------------------------------------------

    etalvc_noxide ~ 0.1771
    # omega^2 = log(0.44^2 + 1) = 0.1771; Weatherley 2018 Table 2 'Without BLQ or taper': IIV Vm = 44%

    etalcl_noxide ~ 0.0654
    # omega^2 = log(0.26^2 + 1) = 0.0654; Weatherley 2018 Table 2 'Without BLQ or taper': IIV CLm = 26%

    # ------------------------------------------------------------------
    # METABOLITE-2 RESIDUAL ERROR
    # ------------------------------------------------------------------

    propSd_noxide <- 0.078
    label("Metabolite-2 proportional residual SD on linear concentration (fraction)")  # Weatherley 2018 Table 2 'Without BLQ or taper': proportional = 7.8%

    addSd_noxide <- 0.10
    label("Metabolite-2 additive residual SD on linear concentration (ng/mL)")  # Weatherley 2018 Table 2 'Without BLQ or taper': additive = 0.10 ng/mL (footnote c flags very-small reported uncertainty)
  })

  model({
    # Reference covariate values (Weatherley 2018 Methods and Table 2
    # footer: reference subject = male, 70 kg, 40 years old).
    ref_wt  <- 70
    ref_age <- 40

    # ------------------------------------------------------------------
    # METABOLITE-1 individual parameters
    # ------------------------------------------------------------------
    # Age enters as an additive offset on the typical CL in L/h
    # (Weatherley 2018 Table 2: 'Age on CL (change in L/h every year
    # above 40)'). Sex enters as a linear-fractional multiplier and
    # body weight enters as fixed-exponent allometric scaling on top
    # of the age-adjusted baseline. The lognormal eta on Metabolite-1
    # CL multiplies the full deterministic typical value.

    cl <- (exp(lcl) + e_age_cl * (AGE - ref_age)) *
          (1 + e_sexf_cl * SEXF) *
          (WT / ref_wt)^e_wt_cl *
          exp(etalcl)

    vc <- exp(lvc) * (WT / ref_wt)^e_wt_vc

    q  <- exp(lq)  * (WT / ref_wt)^e_wt_q

    vp <- exp(lvp) * (WT / ref_wt)^e_wt_vp

    ka <- exp(lka + etalka)

    # ------------------------------------------------------------------
    # METABOLITE-2 individual parameters
    # ------------------------------------------------------------------

    cl_noxide <- exp(lcl_noxide + etalcl_noxide) *
                 (WT / ref_wt)^e_wt_cl_noxide *
                 (1 + e_sexf_cl_noxide * SEXF)

    vc_noxide <- exp(lvc_noxide + etalvc_noxide)

    # ODE system. Doses are administered to the depot compartment as
    # fosdagrocorat (apparent F1 = 1). Elimination of Metabolite-1
    # from the central compartment is assumed to be 100 percent
    # conversion to Metabolite-2 (Fm = 1) so the molar flux out of the
    # parent compartment equals the molar flux into the N-oxide
    # compartment. No stoichiometric mass correction is applied: the
    # reported apparent volumes Vm/F/Fm already absorb the molecular-
    # weight difference between PF-00251802 and PF-04015475 (the
    # +16 g/mol N-oxide adduct, a ~3.5 percent mass increment).
    d/dt(depot)          <- -ka * depot
    d/dt(central)        <-  ka * depot - (cl / vc) * central - (q / vc) * central + (q / vp) * peripheral1
    d/dt(peripheral1)    <-  (q / vc) * central - (q / vp) * peripheral1
    d/dt(central_noxide) <-  (cl / vc) * central - (cl_noxide / vc_noxide) * central_noxide

    # Bioavailability into depot. With lfdepot fixed at log(1) and
    # etalfdepot supplying the IIV approximation of the paper's IOV on
    # F1, this evaluates to a lognormal F1 ~ exp(N(0, omega^2)) per
    # subject. No absorption lag is applied (Weatherley 2018 Table 2:
    # ALAG1 = 0 FIXED).
    f(depot) <- exp(lfdepot + etalfdepot)

    # Observations. Compartment amounts are in mg (matching the dose
    # unit) and apparent volumes are in L, giving a raw concentration
    # in mg/L. Concentrations are scaled by 1000 to express output in
    # ng/mL (mg/L * 1000 = ng/mL), matching the unit Weatherley 2018
    # reports across all figures and tables.
    Cc        <- (central        / vc)        * 1000
    Cc_noxide <- (central_noxide / vc_noxide) * 1000

    Cc        ~ prop(propSd)        + add(addSd)
    Cc_noxide ~ prop(propSd_noxide) + add(addSd_noxide)
  })
}
