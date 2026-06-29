Geldof_2008_fluvoxamine_rat <- function() {
  description <- paste(
    "Preclinical (rat, male Wistar). Non-linear pharmacokinetic brain",
    "distribution model for fluvoxamine in plasma, brain extracellular",
    "fluid (ECF) and total brain tissue, fit by Geldof et al. (2008, Pharm",
    "Res 25(4):792-804) using simultaneous analysis of microdialysate ECF",
    "(n = 26 rats, frontal-cortex CMA/12 probe) and total brain tissue",
    "(n = 35 rats, destructive brain sampling) after a single 30 min IV",
    "infusion of 1, 3.7 or 7.3 mg/kg fluvoxamine. The structural model is",
    "a three-compartment plasma disposition (central + 2 peripherals, with",
    "PK parameters fixed at the mean post-hoc estimates from the upstream",
    "Geldof 2007 rat population PK model, Table I 'Microdialysis + brain",
    "sampling' row) coupled to a single-state lumped brain compartment",
    "(brain_total) whose dynamics follow dCT/dt = kin*Cp - kout*CSP (paper",
    "Eq 10), with the shallow perfusion-limited CSP and deep brain CDB",
    "(= ECF) concentrations recovered algebraically at every time step",
    "from CT via the rapid-equilibrium saturable-efflux quadratic (paper",
    "Appendix Eq 47). The single lumped efflux parameter N***max",
    "(NstarMax in this file) and C50 govern the saturation of the active",
    "(Pgp / MRP-mediated) removal flux from the deep brain back to the",
    "shallow brain. Inter-individual variability is on kin and kout only,",
    "with the correlation reported in Table II. The proportional residual",
    "error sigma^2 = 0.042 is shared between the ECF (Cecf) and total-",
    "brain (Cbrain) observations per Table II.",
    sep = " "
  )
  reference <- paste(
    "Geldof M, Freijer J, van Beijsterveldt L, Danhof M.",
    "Pharmacokinetic modeling of non-linear brain distribution of",
    "fluvoxamine in the rat.",
    "Pharm Res. 2008;25(4):792-804.",
    "doi:10.1007/s11095-007-9390-5.",
    "Plasma PK parameters fixed from the upstream Geldof 2007 rat",
    "popPK model (Eur J Pharm Sci 30(1):45-55;",
    "doi:10.1016/j.ejps.2006.10.001) per Table I",
    "'Microdialysis + brain sampling' row of the 2008 paper.",
    sep = " "
  )
  vignette <- "Geldof_2008_fluvoxamine_rat"

  paper_specific_compartments <- c("brain_total")
  paper_specific_etas <- c("etalkin", "etalkout")
  paper_specific_residual_sds <- c("propSd_Cecf", "propSd_Cbrain")

  units <- list(
    time          = "min",
    dosing        = "ng",
    concentration = "ng/mL"
  )

  covariateData <- list()

  population <- list(
    species        = "rat (male Wistar, Charles River Wiga GmbH, Sulzfeld, Germany)",
    n_subjects     = 61L,
    n_studies      = 1L,
    age_range      = "adult (specific age not reported; housed 6-10 days post-arrival under standard conditions, then 1 week post-cannulation recovery for microdialysis rats or 2 days for brain-sampling rats)",
    weight_range   = "226-250 g body weight at the start of the experiments",
    sex_female_pct = 0,
    race_ethnicity = NA,
    disease_state  = paste(
      "Healthy rats prepared with chronic right-jugular-vein and",
      "left-femoral-artery cannulae for fluvoxamine administration and",
      "arterial blood sampling. The 26 microdialysis-study rats also",
      "carried an indwelling CMA/12 microdialysis probe in the right",
      "frontal cortex (AP +3.2, L -3.0, V -1.5 mm from bregma, Paxinos",
      "and Watson atlas) for sampling brain ECF. The 35 brain-sampling",
      "rats received only the vascular cannulae and were sacrificed at",
      "predetermined times for destructive total-brain tissue assays."
    ),
    dose_range     = paste(
      "Single 30 min IV infusion of fluvoxamine free base into the right",
      "jugular vein at 1, 3.7 or 7.3 mg/kg (flow rate 20 uL/min via a",
      "BAS BeeHive pump). Group sizes in the microdialysis study: 8 / 8 /",
      "10 rats at 1 / 3.7 / 7.3 mg/kg; in the brain-sampling study:",
      "19 / 0 / 16 rats at 1 / 3.7 / 7.3 mg/kg."
    ),
    regions        = "preclinical (in-vivo rat); Leiden University, The Netherlands",
    notes          = paste(
      "Inter-animal weight range is narrow (226-250 g) and the model",
      "does not estimate body weight as a covariate; doses are reported",
      "per kg but the structural plasma PK parameters are absolute (mL,",
      "mL/min) and reflect the typical mid-range Wistar rat. The plasma",
      "concentration-time profile for each animal in the brain submodel",
      "is the post-hoc empirical-Bayes prediction from the upstream",
      "Geldof 2007 popPK model (n = 187 rats from prior plasma-only",
      "studies); inter-individual variability could not be re-estimated",
      "in the current brain submodel and was inherited only as the",
      "central + peripheral structural mean. Table I of Geldof 2008",
      "reports the mean post-hoc plasma PK parameter values pooled",
      "across the 'Microdialysis + brain sampling' combined cohort",
      "(used here), as well as the separate microdialysis-only and",
      "brain-sampling-only cohorts (not used here); a covariate",
      "analysis found no significant difference between the three dose",
      "groups or between the two studies (paper Results, p798).",
      "Brain ECF concentrations measured in the 26 microdialysis rats",
      "ranged 1-214 ng/mL and total brain tissue concentrations were",
      "measured in the 35 brain-sampling rats up to 750 min post-dose;",
      "the LOQ for fluvoxamine was 1 ng/mL in plasma, ECF and brain",
      "(paper Drug Analysis section, p797). An in vivo recovery factor",
      "of 0.27 (mean across the 6 rats without an individual",
      "retrodialysis recovery determination) was used to back-correct",
      "microdialysate concentrations to true ECF concentrations; the",
      "remaining 20 rats used their own retrodialysis-measured",
      "individual recovery values."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Plasma PK -- fixed at the mean post-hoc estimates of the upstream
    # Geldof 2007 rat popPK model (Eur J Pharm Sci 30:45-55), reported in
    # Geldof 2008 Table I row 1 ("Microdialysis + brain sampling" combined
    # cohort, n = 187 rats in the upstream study). The 2008 paper does
    # NOT re-estimate these in the brain submodel; the individual plasma
    # concentration trace per animal was supplied as the input function
    # to the brain ODEs via the post-hoc empirical-Bayes predictions.
    # Inter-individual variability on the plasma parameters is therefore
    # not carried in this model -- only the typical population means
    # are reproduced here.
    # ------------------------------------------------------------------
    lcl  <- fixed(log(31.6)); label("Plasma systemic clearance CL (mL/min)")                        # Geldof 2008 Table I 'Microdialysis + brain sampling' row: CL = 31.6 mL/min
    lvc  <- fixed(log(321));  label("Plasma central volume of distribution V1 (mL)")                # Geldof 2008 Table I row 1: V1 = 321 mL
    lvp  <- fixed(log(949));  label("Plasma peripheral 1 volume of distribution V2 (mL)")           # Geldof 2008 Table I row 1: V2 = 949 mL
    lq   <- fixed(log(33.7)); label("Inter-compartmental clearance central <-> peripheral1 Q2 (mL/min)") # Geldof 2008 Table I row 1: Q2 = 33.7 mL/min
    lvp2 <- fixed(log(136));  label("Plasma peripheral 2 volume of distribution V3 (mL)")           # Geldof 2008 Table I row 1: V3 = 136 mL (V3 had no estimable IIV in the upstream model and uses the n = 187 population estimate; identical across all three Table I rows)
    lq2  <- fixed(log(1.0));  label("Inter-compartmental clearance central <-> peripheral2 Q3 (mL/min)") # Geldof 2008 Table I row 1: Q3 = 1.0 mL/min (also no estimable IIV in the upstream model; identical across rows)

    # ------------------------------------------------------------------
    # Brain distribution structural parameters (Geldof 2008 Table II).
    # kin and kout are the brain-perfusion-derived influx / efflux rate
    # constants on the total brain compartment (paper Eqs 8-10):
    #   kin  = QB / VT       (brain perfusion / brain volume)
    #   kout = QB / (VT * P) (perfusion / (brain volume * partition coef))
    # NstarMax is the lumped Pgp / MRP active-removal-flux saturation
    # parameter N***max from paper Appendix Eq 56 -- formally
    #   N***max = Nmax * VSP / (kdiff * VT)
    # with [conc] units (this skill's interpretation; see Errata block in
    # the validation vignette regarding the 'ng/h' units printed in
    # Table II, which is inconsistent with the appearance of N***max as
    # an additive term against CT and C50 in the partition-coefficient
    # quadratic and is treated here as a paper typo for ng/mL).
    # ------------------------------------------------------------------
    lkin       <- log(0.16);    label("Brain influx rate constant kin (1/min): QB / VT")            # Geldof 2008 Table II: kin = 0.16 /min (CV 13.6%)
    lkout      <- log(0.019);   label("Brain efflux rate constant kout (1/min): QB / (VT * P)")     # Geldof 2008 Table II: kout = 0.019 /min (CV 8.1%)
    lNstarMax  <- log(30700);   label("Lumped saturable active removal flux capacity N***max (ng/mL): Nmax * VSP / (kdiff * VT)") # Geldof 2008 Table II: N***max = 30,700 (CV 92.5%); Table II prints units as ng.h^-1 but dimensional analysis of paper Eq 57 requires N***max to be a concentration -- this is treated as a paper typo for ng/mL (see vignette Errata)
    lc50       <- log(710);     label("Fluvoxamine concentration in deep brain at which the active removal flux is 50% saturated, C50 (ng/mL)") # Geldof 2008 Table II: C50 = 710 ng/mL (CV 96.8%)

    # ------------------------------------------------------------------
    # Inter-individual variability (Table II "w^2*kin", "w^2*kout",
    # "w*kin / w*kout"). The paper reports variance estimates on the
    # log-scale (NONMEM omega^2): IIV(kin) = 0.50, IIV(kout) = 0.17.
    # The third row "w*kin / w*kout = 0.24" is the correlation
    # coefficient r between the two etas.
    #   cov(eta_kin, eta_kout) = r * sqrt(var_kin) * sqrt(var_kout)
    #                          = 0.24 * sqrt(0.50) * sqrt(0.17)
    #                          = 0.0700
    # IIV could not be adequately estimated on N***max or C50 and was
    # fixed to zero in the paper.
    # ------------------------------------------------------------------
    etalkin + etalkout ~ c(0.50, 0.0700, 0.17)  # Geldof 2008 Table II: omega^2(kin) = 0.50 (CV 25.7%), omega^2(kout) = 0.17 (CV 28.5%), corr(eta_kin, eta_kout) = 0.24 (CV 32.9%) -- covariance = 0.24 * sqrt(0.50 * 0.17) = 0.0700

    # ------------------------------------------------------------------
    # Residual error. Geldof 2008 Eq 19:
    #   Cmij = Cij * (1 + eps_ij),  eps_ij ~ N(0, sigma^2)
    # where Cmij is the measured concentration (ECF or total brain) and
    # Cij the model prediction. The same proportional-error term sigma^2
    # = 0.042 (Table II 'sigma^2 (eps1ij)') is applied to BOTH the ECF
    # and the total-brain observations -- the paper estimates a single
    # residual variance pooled across the two output streams. The
    # corresponding proportional SD is sqrt(0.042) = 0.2049 (fraction).
    # Encoded here as two separate parameters with the same numerical
    # value so each multi-output residual line ('Cecf ~ prop(...)',
    # 'Cbrain ~ prop(...)') can declare its own per-output residual SD;
    # the shared-value nature is documented in the vignette Assumptions
    # and deviations section.
    # ------------------------------------------------------------------
    propSd_Cecf   <- 0.2049; label("Proportional residual SD on brain ECF observations Cecf (fraction)")          # Geldof 2008 Table II: sigma^2(eps1ij) = 0.042 (CV 17.3%); SD = sqrt(0.042) = 0.2049; shared with Cbrain per paper Eq 19
    propSd_Cbrain <- 0.2049; label("Proportional residual SD on total brain tissue observations Cbrain (fraction)") # Geldof 2008 Table II: sigma^2(eps1ij) = 0.042 (CV 17.3%); SD = sqrt(0.042) = 0.2049; shared with Cecf per paper Eq 19
  })

  model({
    # ------------------------------------------------------------------
    # Plasma 3-compartment disposition (fixed structural means; no IIV).
    # The dose is administered as a 30 min IV infusion to the central
    # compartment via the AMT + RATE columns of the event table.
    # ------------------------------------------------------------------
    cl   <- exp(lcl)
    vc   <- exp(lvc)
    vp   <- exp(lvp)
    q    <- exp(lq)
    vp2  <- exp(lvp2)
    q2   <- exp(lq2)

    kel  <- cl  / vc
    k12  <- q   / vc
    k21  <- q   / vp
    k13  <- q2  / vc
    k31  <- q2  / vp2

    # ------------------------------------------------------------------
    # Brain catenary distribution (Geldof 2008 main text Eqs 8-10 and
    # Appendix Eqs 36, 47, 55-57). After eliminating the diffusion rate
    # constant kdiff via the rapid-equilibrium assumption (paper Methods
    # paragraph 4 and Appendix Eqs 25-27), the shallow and deep brain
    # compartments collapse to a single ODE on the total brain
    # concentration CT:
    #   dCT/dt = kin*Cp - kout*CSP                                  (Eq 10)
    # where CSP is the concentration in the shallow perfusion-limited
    # brain area, computed at every t from CT via the partition
    # coefficients fDB = CDB/CT and fSP = CSP/CT.
    # The partition coefficient fDB solves the rapid-equilibrium
    # quadratic obtained from the saturable Michaelis-Menten efflux
    # from the deep brain back to the shallow brain (paper Appendix
    # Eq 47, recast in concentration units and assuming VSP = VDB so
    # that the lumped N***max enters as a concentration):
    #   CDB^2 + (C50 + N***max - CT) * CDB - CT * C50 = 0
    # which has physical root
    #   CDB = ( (CT - C50 - N***max)
    #          + sqrt( (C50 + N***max - CT)^2 + 4*CT*C50 ) ) / 2
    # i.e.
    #   fDB = CDB / CT
    #       = ( (CT - C50 - N***max) + sqrt(disc) ) / (2 * CT)
    # NOTE: The published formula in Eq 57 of Geldof 2008 prints
    # (+N***max + CT - C50 + sqrt(...)) in the numerator, which gives
    # fDB -> infinity as CT -> 0 and is therefore inconsistent with
    # the limit fDB -> C50/(C50 + N***max) implied by the low-CT
    # linearization of Eq 47. The sign of N***max in Eq 57 of Geldof
    # 2008 is treated here as a publisher typesetting error; the
    # implementation uses the derivation-consistent (CT - C50 - N***max)
    # form. See vignette Errata for the full provenance.
    # The corresponding partition coefficient fSP for the shallow
    # compartment is read off the brain mass balance VT*CT = VSP*CSP +
    # VDB*CDB and Eq 16 of the main text:
    #   fSP = 1 + (VDB / VSP) * (1 - fDB)
    # Setting VSP = VDB (the simplest paper-consistent assumption,
    # under which N***max enters as a concentration in the CT quadratic
    # and CT = (CSP + CDB) / 2 ), the relation reduces to
    #   fSP = 2 - fDB
    # ------------------------------------------------------------------
    kin      <- exp(lkin + etalkin)
    kout     <- exp(lkout + etalkout)
    NstarMax <- exp(lNstarMax)
    c50      <- exp(lc50)

    cp <- central / vc

    # Algebraic partition coefficients from the brain quadratic
    # (use a small epsilon in the denominator to keep the limit at
    # brain_total = 0 well behaved during numerical integration).
    bt    <- brain_total
    diff_ <- bt - c50 - NstarMax
    disc  <- diff_ * diff_ + 4 * bt * c50
    cdb   <- (diff_ + sqrt(disc)) / 2
    csp   <- 2 * bt - cdb

    # ODE system
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2
    d/dt(brain_total) <-  kin * cp - kout * csp

    # ------------------------------------------------------------------
    # Observations
    # ------------------------------------------------------------------
    Cc     <- cp           # plasma fluvoxamine concentration (ng/mL); diagnostic output, not fit in this paper (no residual SD reported)
    Cecf   <- cdb          # brain ECF (deep brain) fluvoxamine concentration measured via microdialysis
    Cbrain <- bt           # total brain tissue fluvoxamine concentration measured by destructive sampling

    Cecf   ~ prop(propSd_Cecf)
    Cbrain ~ prop(propSd_Cbrain)
  })
}
