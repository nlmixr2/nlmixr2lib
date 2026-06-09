Goggin_2004_emfilermin <- function() {
  description <- paste(
    "One-compartment population PK model for subcutaneous emfilermin",
    "(recombinant human leukaemia inhibitory factor, r-hLIF) in healthy",
    "postmenopausal women and in infertile women undergoing in vitro",
    "fertilization and embryo transfer (IVF-ET) (Goggin 2004).",
    "Absorption is zero-order (D1 = 0.84 h, invariant, no IIV) directly",
    "into the central compartment, followed by first-order elimination.",
    "Apparent clearance CL/F is decreased by 35% in IVF-ET patients",
    "(typical 37 L/h) relative to healthy postmenopausal women (typical",
    "57 L/h). Apparent volume V/F is linear in body weight on the natural",
    "scale: V/F = 235 L at the median 62 kg, increasing or decreasing by",
    "6.7 L/kg (~29% per 10 kg) -- an absolute-linear covariate form, not",
    "log-multiplicative. Inter-individual variability is log-normal on",
    "CL/F (17% CV) and V/F (28% CV); inter-occasion variability is",
    "log-normal on V/F (23% CV) across three protocol-defined occasions",
    "(first dosing day = 1, intermediate dosing days = 2, last dosing",
    "day = 3). Residual error is proportional (20% CV). Studied weight",
    "range was 48-83 kg; the linear V/F-WT term is extrapolation-unsafe",
    "below ~27 kg where the typical V/F would become negative."
  )
  reference <- paste(
    "Goggin T, Nguyen QTX, Munafo A.",
    "Population pharmacokinetic modelling of Emfilermin (recombinant",
    "human leukaemia inhibitory factor, r-hLIF) in healthy",
    "postmenopausal women and in infertile patients undergoing in vitro",
    "fertilization and embryo transfer.",
    "Br J Clin Pharmacol. 2004;57(4):412-418.",
    "doi:10.1111/j.1365-2125.2003.02064.x"
  )
  vignette <- "Goggin_2004_emfilermin"
  units <- list(time = "h", dosing = "ug", concentration = "pg/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear absolute effect on V/F centered at the cohort-median 62 kg (Goggin 2004 Methods Step 2 and Table 5): TVP_V = 235 + 6.7 * (WT - 62). Studied weight range 48-83 kg; the typical-V/F term becomes negative below ~27 kg, so the model is extrapolation-unsafe outside the studied range.",
      source_name        = "WGT"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-participant cohort indicator: 1 = healthy postmenopausal woman on hormone replacement therapy (phase I studies 1 and 2), 0 = infertile premenopausal woman undergoing IVF-ET (proof-of-concept study 3).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (IVF-ET patient)",
      notes              = "Time-fixed per subject. Goggin 2004 Methods Step 2 encodes the type-of-population covariate with TYPE = 1 for IVF-ET patients and reports the typical CL/F = 57 L/h standardised for healthy postmenopausal women (Table 4, row 'CL/F'). The canonical DIS_HEALTHY convention uses 0 = patient as reference, so the structural typical lcl is shifted to the IVF-ET (patient) state: lcl = log(57 * 0.649) = log(37.0); the e_dis_healthy_cl coefficient (= -log(0.649)) restores the paper's PM-typical 57 L/h at DIS_HEALTHY = 1. Mathematically identical to the paper's TVP_CL = Ppop * 0.649^TYPE with TYPE = 1 - DIS_HEALTHY.",
      source_name        = "TYPE"
    ),
    OCC = list(
      description        = "Integer-valued occasion / period indicator for inter-occasion-variability multiplexing.",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Three protocol-defined occasions per Goggin 2004 Methods Step 1 ('In this analysis three occasions were defined: occasion 1, the first day of treatment; occasion 3, the last day of treatment (either day 6 or 7); and occasion 2, all intermediate days (mostly day 4).'). Decomposed inside model() into binary indicators oc1 / oc2 / oc3 that multiplex the three IOV etas on log-V/F. Single estimated variance shared across occasions (NONMEM $OMEGA BLOCK(1) SAME pattern; nlmixr2 has no SAME shortcut so occasions 2 and 3 are fix()'d to the occasion-1 variance).",
      source_name        = "OCC"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 64L,
    n_studies      = 3L,
    age_range      = "30-64 years (across cohorts)",
    age_median     = "58 years (postmenopausal cohorts) / 34 years (IVF-ET cohort)",
    weight_range   = "48-83 kg (across cohorts)",
    weight_median  = "62 kg (overall cohort median; PM = 63 kg, IVF-ET = 61 kg)",
    sex_female_pct = 100,
    race_ethnicity = "Not reported in the paper. Studies conducted in the UK (Cambridge, Nottingham).",
    disease_state  = "Two cohorts: (1) healthy postmenopausal women on Cyclo-progynova/Progynova/Utrogestan hormone replacement therapy (studies 1 and 2; n = 25 after exclusion of two subjects with bioanalytical or dosing problems); (2) premenopausal women with recurrent implantation failure undergoing in vitro fertilization or intracytoplasmic sperm injection and embryo transfer (study 3; n = 39), pretreated with nafarelin pituitary down-regulation, recombinant FSH ovarian stimulation, and recombinant hCG triggering.",
    dose_range     = "Subcutaneous emfilermin 100 ug or 250 ug single dose (study 1), 150 ug BID for 7 days (study 2), or 150 ug BID for 7 days starting on the day of embryo transfer (study 3).",
    regions        = "United Kingdom (Cambridge: Addenbrooke's Hospital, Bourn Hall Clinic, Papworth Hospital; Nottingham: The Park Hospital).",
    notes          = "Subject characteristics (Table 1): study 1 (n = 14) median age 58 years (range 50-64), median WT 63 kg (range 55-70); study 2 (n = 11) median age 58 years (range 49-63), median WT 63 kg (range 56-80); study 3 (n = 39) median age 34 years (range 30-37), median WT 61 kg (range 48-83). Analysis dataset: 342 samples from 64 subjects (226 from postmenopausal women, 116 from IVF-ET patients) after exclusion of subject 113 (suspected dosing error) and subject 116 (haemolyzed samples). Below-LoQ samples (LoQ = 50 pg/mL) treated as missing. ALT, AST, creatinine, and bilirubin were screened as univariate covariates but not retained in the final model (Methods Step 3, Table 3). Population type and age were confounded; the analysis retains population type rather than age because the OBJ drop with age alone was only 3.5 points vs 33.4 for population type."
  )

  ini({
    # -------------------------------------------------------------------
    # Structural PK parameters - Goggin 2004 Table 4 (final population
    # pharmacokinetics model estimates). Paper standardises CL/F at the
    # healthy postmenopausal reference; the canonical DIS_HEALTHY register
    # uses 0 = patient as reference, so lcl is shifted to the IVF-ET
    # typical 37 L/h (= 57 * 0.649) and the healthy effect restores 57
    # L/h at DIS_HEALTHY = 1. See covariateData$DIS_HEALTHY notes.
    # -------------------------------------------------------------------
    lcl <- log(57) + log(0.649)  ; label("Apparent clearance CL/F at IVF-ET (DIS_HEALTHY = 0) reference (L/h); = log(57.0 * 0.649) = log(37.0)")  # Goggin 2004 Table 4: CL/F = 57.0 L/h (PM-typical) and Type-on-CL/F = 0.649 multiplier for IVF-ET
    lvc <- log(235)              ; label("Apparent central volume of distribution V/F at WT = 62 kg (L)")                                       # Goggin 2004 Table 4: V/F = 235 L (typical at the cohort-median 62 kg)
    ld1 <- log(0.84)             ; label("Duration of zero-order subcutaneous absorption D1 (h)")                                               # Goggin 2004 Table 4: D1 = 0.84 h (no IIV estimated)

    # -------------------------------------------------------------------
    # Covariate effects.
    # -------------------------------------------------------------------
    # Healthy-vs-patient (DIS_HEALTHY) multiplicative shift on log-CL.
    # Paper form: TVP_CL = Ppop * 0.649^TYPE with TYPE = 1 for IVF-ET.
    # With DIS_HEALTHY = 1 - TYPE, the equivalent log-additive form is
    # log(TVP_CL) = lcl + log(57/37) * DIS_HEALTHY = lcl + 0.4325 * DIS_HEALTHY.
    e_dis_healthy_cl <- -log(0.649)  ; label("Effect of DIS_HEALTHY = 1 (healthy PM) on log-CL/F (unitless); = -log(0.649) = log(57/37) = +0.4325")  # Goggin 2004 Table 4: Type on CL/F = 0.649 (multiplier for IVF-ET = 1 - DIS_HEALTHY)

    # Body-weight linear additive effect on V/F on the natural scale.
    # Paper form: TVP_V = Ppop + theta_WT * (WT - 62), with reference WT
    # = 62 kg (cohort median, Methods Step 2). Effect is NOT
    # log-multiplicative; do not exponentiate.
    e_wt_vc <- 6.7  ; label("Linear additive effect of WT on V/F (L per kg WT above 62 kg reference)")  # Goggin 2004 Table 4: Weight on V/F = 6.7 (L per kg, linear absolute slope)

    # -------------------------------------------------------------------
    # Inter-individual variability - Goggin 2004 Table 4. Paper Methods
    # Step 1 reports IIV as a log-normal multiplicative random effect
    # Pi = TVP * exp(eta_i) and Table 4 lists sqrt(omega^2) as %CV. The
    # internal log-scale variance follows omega^2 = log(CV^2 + 1).
    # -------------------------------------------------------------------
    etalcl ~ 0.02849  # Goggin 2004 Table 4: ISV(CL/F) = 17% CV  -> omega^2 = log(0.17^2 + 1) = 0.02849
    etalvc ~ 0.07551  # Goggin 2004 Table 4: ISV(V/F)  = 28% CV  -> omega^2 = log(0.28^2 + 1) = 0.07551

    # -------------------------------------------------------------------
    # Inter-occasion variability on V/F. Three protocol-defined occasions
    # (Methods Step 1) share a single estimated variance (p^2). nlmixr2
    # has no $OMEGA BLOCK(1) SAME shortcut so the second and third
    # occasions are fix()'d to the occasion-1 variance, matching the
    # Jonsson 2011 ethambutol pattern.
    # -------------------------------------------------------------------
    etaiov_vc_1 ~ 0.05155            # Goggin 2004 Table 4: IOV(V/F) = 23% CV -> p^2 = log(0.23^2 + 1) = 0.05155 (estimated, occasion 1)
    etaiov_vc_2 ~ fix(0.05155)       # occasion-2 variance fixed equal to occasion-1 ($OMEGA BLOCK(1) SAME pattern)
    etaiov_vc_3 ~ fix(0.05155)       # occasion-3 variance fixed equal to occasion-1 ($OMEGA BLOCK(1) SAME pattern)

    # -------------------------------------------------------------------
    # Residual error - proportional (Goggin 2004 Methods Step 1; combined
    # additive + proportional explored but additive component's 95% CI
    # spanned zero, so dropped in the final model, Results paragraph 2).
    # -------------------------------------------------------------------
    propSd <- 0.20  ; label("Proportional residual error on Cc (fraction)")  # Goggin 2004 Table 4: RV = 20% CV
  })

  model({
    # Decompose the integer-valued occasion column into binary indicators
    # for IOV multiplexing on log-V/F. Three occasions defined in the
    # paper Methods Step 1: occasion 1 = first dosing day; occasion 2 =
    # all intermediate days (mostly day 4); occasion 3 = last dosing day
    # (day 6 or 7).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    iov_vc <- oc1 * etaiov_vc_1 + oc2 * etaiov_vc_2 + oc3 * etaiov_vc_3

    # Individual PK parameters.
    # CL/F: multiplicative log-normal IIV plus DIS_HEALTHY shift. lcl is
    # parameterised at the IVF-ET (DIS_HEALTHY = 0) reference; the
    # +e_dis_healthy_cl term restores the PM typical 57 L/h at
    # DIS_HEALTHY = 1.
    cl <- exp(lcl + etalcl + e_dis_healthy_cl * DIS_HEALTHY)

    # V/F: linear-additive WT effect on the natural scale, then
    # multiplicative log-normal IIV + IOV. Paper Methods Step 2:
    #   TVP_V(WT) = Ppop + theta_WT * (WT - 62)
    #   V_i      = TVP_V(WT) * exp(eta_i + iov_i)
    tvp_vc <- exp(lvc) + e_wt_vc * (WT - 62)
    vc <- tvp_vc * exp(etalvc + iov_vc)

    # D1: zero-order absorption duration, no IIV.
    d1 <- exp(ld1)

    # Elimination rate constant for the one-compartment system.
    kel <- cl / vc

    # One-compartment SC PK with zero-order input directly into central.
    # Dose targets cmt = "central" with rxode2's dur() control parameter
    # setting the infusion duration to D1; the dose event must carry
    # rate = -2 (dataset RATE column) to request model-defined duration.
    d/dt(central) <- -kel * central

    # Zero-order input duration on the central compartment.
    dur(central) <- d1

    # Plasma concentration. Dose units ug, central in ug, vc in L:
    # central / vc gives ug/L. The Goggin 2004 ELISA reports r-hLIF in
    # pg/mL (paper Methods 'Analysis of r-hLIF', LoQ = 50 pg/mL,
    # standard range 31.2-2000 pg/mL), and 1 ug/L = 1 ng/mL = 1000 pg/mL,
    # so multiply by 1000 to land in pg/mL.
    Cc <- central / vc * 1000

    Cc ~ prop(propSd)
  })
}
