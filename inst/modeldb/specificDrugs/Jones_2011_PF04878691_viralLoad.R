Jones_2011_PF04878691_viralLoad <- function() {
  description <- paste(
    "Combined PK + OAS + HCV viral RNA pharmacodynamic chain for oral",
    "PF-04878691 (TLR7 agonist) used to predict the antiviral efficacy of",
    "PF-04878691 in chronic hepatitis C (HCV) patients (Jones 2011 BJCP",
    "Figure 10 simulation). PK is the two-compartment time-varying clearance",
    "model from Jones_2011_PF04878691.R (Table 1; all PK structural",
    "parameters fixed at the published Table 1 values plus the IIVs on CL_SS",
    "and ka). The PF-04878691 OAS indirect-response sub-model is the same as",
    "Jones_2011_PF04878691_oas.R (Table 2; OAS in fold-change units, baseline",
    "rbase_oas = 0.96, drug stimulates production through slope * Cc^gamma).",
    "The OAS-viral-load relationship was fit on TLR9-agonist (CPG-10101)",
    "data in HCV patients (Jones 2011 Table 4) and is assumed transferable to",
    "PF-04878691 under the paper's explicit translation assumption that",
    "'both TLR7 and TLR9 work through the same pathway'. The viral-load",
    "model is an inhibitory sigmoid Imax driven by the change in OAS from",
    "baseline expressed as a fold change",
    "oas_fc_above = oas(t) / rbase_oas - 1, so at the OAS baseline the",
    "viral-load deviation from BASE is zero. The viral RNA observation",
    "(vload) is in log10 copies/mL."
  )
  reference <- paste(
    "Jones HM, Chan PLS, van der Graaf PH, Webster R. Use of modelling and",
    "simulation techniques to support decision making on the progression of",
    "PF-04878691, a TLR7 agonist being developed for hepatitis C.",
    "Br J Clin Pharmacol. 2012;73(1):77-92.",
    "doi:10.1111/j.1365-2125.2011.04047.x.",
    "PK + OAS structure inherited from companion files Jones_2011_PF04878691.R",
    "and Jones_2011_PF04878691_oas.R (Table 1, Table 2). The OAS-viral-load",
    "layer parameters were fit on the McHutchison 2007 CPG-10101 HCV-patient",
    "OAS / HCV-viral-RNA data set referenced as Jones 2011 reference [16]",
    "(McHutchison JG, et al. Hepatology 2007;46(4):1341-9)."
  )
  vignette <- "Jones_2011_PF04878691_HCV"
  paper_specific_compartments <- c("oas")
  paper_specific_etas <- c(
    "etalrbase_oas", "etalkout_oas",
    "etalimax_vl", "etalvo50_vl"
  )
  paper_specific_residual_sds <- c("addSd_vload")
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required scaling covariate for the inherited per-kg PK parameters; see Jones_2011_PF04878691.R. Population median 79 kg (PF-04878691 healthy volunteer cohort); HCV patient cohort body-weight distribution is not tabulated in Jones 2011 -- vignette users simulate at the PF-04878691 cohort median.",
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 39L,
    n_studies      = 2L,
    age_range      = "21-55 years (PF-04878691 cohort); HCV cohort age not separately tabulated in Jones 2011",
    age_median     = "34 years (PF-04878691 cohort)",
    weight_range   = "57-97 kg (PF-04878691 cohort); HCV cohort weight not separately tabulated in Jones 2011",
    weight_median  = "79 kg (PF-04878691 cohort)",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not tabulated in Jones 2011.",
    disease_state  = "PF-04878691 PK + OAS parameters fit on healthy adult volunteers (NCT00810758, n = 24); OAS-viral-load parameters fit on chronic-hepatitis-C patients receiving subcutaneous CPG-10101 (a TLR9 agonist, Jones 2011 reference [16], n = 39 contributing patients of 60 total HCV randomised to placebo or CPG-10101 0.25 / 1 / 4 / 10 / 20 mg twice weekly for 4 weeks or 0.5 / 0.75 mg/kg once weekly for 4 weeks).",
    dose_range     = "PF-04878691 in healthy volunteers: 3, 6, 9 mg orally twice weekly. CPG-10101 in HCV patients: 0.25-20 mg subcutaneously twice weekly or 0.5-0.75 mg/kg once weekly (only used to fit the OAS-viral-load layer; not the PF-04878691 PK or OAS layers).",
    regions        = "Not specified.",
    notes          = "Chain assumes 'no PK or biomarker difference between the healthy volunteer and HCV patient populations' and 'the subsequent intracellular signalling pathways are similar after binding of both TLR7 and TLR9 agonists to their receptors; thus the magnitude of OAS and IFN response required for an antiviral effect by each of the pathways would be comparable' (Jones 2011 Discussion). These are explicit paper assumptions and were the basis for the Figure 10 simulations that recommended discontinuation of PF-04878691."
  )

  ini({
    # ----------------------------------------------------------------------
    # PK structural parameters from Jones_2011_PF04878691.R (Table 1). All
    # fixed because the OAS / viral-load layers were fit sequentially using
    # the upstream popPK EBE PK parameters.
    # ----------------------------------------------------------------------
    lcl       <- fixed(log(1.7));    label("Steady-state apparent clearance per kg body weight (CL_SS = paper CLF, L/h/kg)")              # Table 1 (CLF = 1.7 L/h/kg)
    lcl_time  <- fixed(log(1.8));    label("Initial offset of the time-varying clearance component per kg (CL_TIME0 = CL0 - CLF, L/h/kg)") # Derived from Table 1 (CL0 = 3.5, CLF = 1.7)
    lkdeg     <- fixed(log(0.24));   label("Exponential decay rate of the time-varying clearance component (paper DEG, 1/h)")              # Table 1 (DEG = 0.24 1/h)
    lvc       <- fixed(log(3.3));    label("Apparent central volume of distribution per kg body weight (Vc, L/kg)")                        # Table 1 (Vc = 3.3 L/kg)
    lq        <- fixed(log(0.74));   label("Apparent intercompartmental clearance per kg body weight (Q, L/h/kg)")                         # Table 1 (Q = 0.74 L/h/kg)
    lvp       <- fixed(log(21));     label("Apparent peripheral volume of distribution per kg body weight (Vp, L/kg)")                     # Table 1 (Vp = 21 L/kg)
    lka       <- fixed(log(0.078));  label("First-order absorption rate constant (ka, 1/h)")                                               # Table 1 (ka = 0.078 1/h)

    etalcl ~ fixed(0.067)                                                                                                                  # Table 1 (IIV CLF = 0.067)
    etalka ~ fixed(0.19)                                                                                                                   # Table 1 (IIV ka  = 0.19)

    # ----------------------------------------------------------------------
    # OAS indirect-response PD layer (Jones 2011 Table 2). Same structure
    # and parameter values as Jones_2011_PF04878691_oas.R; OAS parameters
    # carry the _oas suffix here to disambiguate from the viral-load layer.
    # ----------------------------------------------------------------------
    lkout_oas   <- log(0.034); label("OAS first-order elimination rate constant kout (1/h)")                                                # Table 2 (kout = 0.034 1/h)
    lslope_oas  <- log(3.5);   label("Slope of the drug stimulation of OAS production per (ng/mL)^gamma (paper SLP)")                       # Table 2 (SLP = 3.5)
    lrbase_oas  <- log(0.96);  label("OAS baseline fold change (rbase = kin/kout, unitless fold change)")                                   # Table 2 (BASE = 0.96)
    lgamma_oas  <- log(1.6);   label("Sigmoidicity exponent on Cc in the OAS stimulation power function (gamma, unitless)")                 # Table 2 (gamma = 1.6)

    etalkout_oas  ~ 1.7                                                                                                                    # Table 2 (IIV kout = 1.7)
    etalrbase_oas ~ 0.18                                                                                                                   # Table 2 (IIV BASE = 0.18)

    # OAS residual error is not modelled here because the OAS state is an
    # intermediate variable in the simulation chain (its only role is to
    # drive the viral-load response). Users who want the OAS observation
    # error should use the standalone Jones_2011_PF04878691_oas.R model.

    # ----------------------------------------------------------------------
    # OAS -> HCV viral RNA inhibitory sigmoid Imax model (Jones 2011 Table 4,
    # CPG-10101 cohort). The viral-load model takes the OAS deviation from
    # baseline expressed as fold change above 1 (oas_fc_above = oas / rbase
    # - 1) as its input, so the no-drug intercept is BASE = 7.3 log10
    # copies/mL. Imax is reported as a negative log10 copies/mL value
    # because increasing OAS reduces viral load; it is encoded here as
    # imax_abs <- exp(limax_vl) and then sign-flipped inside model() as
    # imax_signed = -imax_abs, so the log-transform applies to the
    # positive-only magnitude.
    # ----------------------------------------------------------------------
    lbase_vl <- log(7.3);     label("Baseline HCV viral RNA at OAS fold change above baseline = 0 (log10 copies/mL)")                       # Table 4 (BASE = 7.3 log10 copies/mL, %CV 1.4)
    limax_vl <- log(2.7);     label("Magnitude of the maximum HCV viral RNA reduction (|Imax|, log10 copies/mL; applied as -Imax_abs inside model())") # Table 4 (Imax = -2.7 log10 copies/mL, %CV 15)
    lvo50_vl <- log(3.6);     label("OAS fold change above baseline producing 50 percent of |Imax| (VO50, unitless fold change)")           # Table 4 (VO50 = 3.6, %CV 11)
    lgamma_vl <- log(0.68);   label("Sigmoidicity exponent on the OAS fold change in the viral-load inhibition function (gamma, unitless)") # Table 4 (gamma = 0.68, %CV 25)

    # Viral-load IIVs (Table 4). The Jones 2011 Results text describes the
    # IIV block as "variance covariance matrix for Imax and gamma" and
    # later attributes the high 64.3 percent IIV to "g" (gamma). Table 4
    # however lists OM2 = IIV Imax (0.29) and OM3 = IIV VO50 (0.25) -- the
    # table indexing therefore places IIV on Imax and VO50, not on Imax
    # and gamma. This is a paper inconsistency; the table is treated as
    # authoritative for the encoded variances (the variances 0.29 and 0.25
    # are the reported point estimates regardless of which parameter they
    # belong to). The text-vs-table discrepancy is recorded in the vignette
    # Errata. No off-diagonal covariance is reported in Table 4, so the
    # IIVs are encoded as independent diagonal omegas.
    etalimax_vl ~ 0.29                                                                                                                     # Table 4 (IIV Imax = 0.29, %CV 36)
    etalvo50_vl ~ 0.25                                                                                                                     # Table 4 (IIV VO50 = 0.25, %CV 64)

    # Viral-load residual error: the only model in this paper that uses an
    # additive residual on the log10 domain (Methods: "the TLR9 OAS-viral
    # load model where residual error was described by an additive model in
    # the log10 domain"). Table 4 reports residual error = 0.41 (%CV 6.3).
    addSd_vload <- 0.41; label("HCV viral RNA additive residual SD on the log10 scale (log10 copies/mL)")                                   # Table 4 (additive residual = 0.41, %CV 6.3)
  })

  model({
    # ----------------------------------------------------------------------
    # Inherited PK (Table 1). See Jones_2011_PF04878691.R for the standalone
    # PK model rationale.
    # ----------------------------------------------------------------------
    cl_ss    <- exp(lcl + etalcl) * WT
    cl_time0 <- exp(lcl_time) * WT
    kdeg_cl  <- exp(lkdeg)
    vc       <- exp(lvc) * WT
    q        <- exp(lq)  * WT
    vp       <- exp(lvp) * WT
    ka       <- exp(lka + etalka)

    cl  <- cl_ss + cl_time0 * exp(-kdeg_cl * time)
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    Cc <- 1000 * central / vc   # PF-04878691 plasma concentration in ng/mL

    # ----------------------------------------------------------------------
    # OAS indirect-response PD (Jones 2011 Table 2). Identical to the
    # standalone OAS file; see Jones_2011_PF04878691_oas.R.
    # ----------------------------------------------------------------------
    rbase_oas <- exp(lrbase_oas + etalrbase_oas)
    kout_oas  <- exp(lkout_oas  + etalkout_oas)
    slope_oas <- exp(lslope_oas)
    gamma_oas <- exp(lgamma_oas)

    kin_oas <- rbase_oas * kout_oas

    oas(0) <- rbase_oas

    d/dt(oas) <- kin_oas * (1 + slope_oas * Cc^gamma_oas) - kout_oas * oas

    # ----------------------------------------------------------------------
    # OAS-driven HCV viral RNA inhibition (Jones 2011 Table 4). The input
    # to the sigmoid Imax function is the OAS fold change ABOVE baseline:
    #
    #     oas_fc_above = oas(t) / rbase_oas - 1
    #
    # so at the OAS baseline (oas = rbase_oas) the deviation is zero and
    # vload = BASE_vl. Above baseline, vload decreases by Imax * sigmoid.
    # max(., 0) guards the sigmoid against transient negative inputs that
    # the ODE could produce in early simulation time-steps (the power
    # gamma_vl = 0.68 is non-integer so a negative base would be invalid).
    # ----------------------------------------------------------------------
    base_vl  <- exp(lbase_vl)
    imax_vl  <- -exp(limax_vl + etalimax_vl)   # signed Imax: negative log10 copies/mL
    vo50_vl  <- exp(lvo50_vl + etalvo50_vl)
    gamma_vl <- exp(lgamma_vl)

    oas_fc_above <- max(oas / rbase_oas - 1, 0)

    vload <- base_vl + imax_vl * oas_fc_above^gamma_vl /
                       (vo50_vl^gamma_vl + oas_fc_above^gamma_vl)

    vload ~ add(addSd_vload)
  })
}
