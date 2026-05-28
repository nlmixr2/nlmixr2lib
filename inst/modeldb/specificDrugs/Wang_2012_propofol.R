Wang_2012_propofol <- function() {
  description <- "Three-compartment intravenous population PK model for propofol across the human life-span (Wang 2012; 174 subjects pooled across seven previously published studies covering preterm and term neonates, infants, toddlers, children, adolescents, and adults; body weight 0.68-122.7 kg, age 1 day-81 years). Final 'bodyweight-dependent exponent (BDE)' model (Model IV / Final PK model, Table IV): clearance is scaled by total body weight via a power function whose exponent k changes sigmoidally with body weight from k0 = 1.34 at a theoretical 0 kg to k0 - kmax = 0.55 at large body weights, with k50 = 3.78 kg and a Hill coefficient gamma = 5.24 governing the steepness of the decline. The slow inter-compartmental clearance Q3 and the second peripheral volume V3 scale linearly with body weight (BW/70); the first peripheral volume V2 scales as (BW/70)^0.55; the fast inter-compartmental clearance Q2 is independent of body weight. The central volume V1 = 7.58 L is constant for subjects with postnatal age >= 100 days and scales linearly as V1 * (BW/70) for younger subjects. Inter-individual variability (log-normal) was retained on CL, V1, V2, V3, and Q3; no IIV on Q2. Additive residual error on log-transformed concentrations was used, equivalent to a proportional error on the linear concentration scale."
  reference <- paste(
    "Wang C, Peeters MYM, Allegaert K, Blusse van Oud-Alblas HJ,",
    "Krekels EHJ, Tibboel D, Danhof M, Knibbe CAJ. (2012).",
    "A bodyweight-dependent allometric exponent for scaling clearance",
    "across the human life-span. Pharmaceutical Research 29(6):1570-1581.",
    "doi:10.1007/s11095-012-0668-x.",
    sep = " "
  )
  vignette <- "Wang_2012_propofol"
  units    <- list(time = "minute", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight (baseline; time-fixed within each study).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used (i) as the allometric scaler on CL via the bodyweight-dependent exponent k = k0 - kmax * WT^hill / (k50^hill + WT^hill) with reference 70 kg, (ii) linearly on Q3 and V3 (reference 70 kg), (iii) as a power on V2 with estimated exponent m = 0.55 (reference 70 kg), and (iv) linearly on V1 only when postnatal age is less than 100 days (reference 70 kg); V1 is constant otherwise. Source range 0.68-122.7 kg pooled across the seven studies (Table I).",
      source_name        = "WT"
    ),
    PNA = list(
      description        = "Postnatal age (chronological time since birth).",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used as a binary switch on the V1 covariate model: V1 = V1p * (WT/70) when PNA < 100 days (~3.286 months), V1 = V1p otherwise (Wang 2012 Results paragraph identifying the V1 linear-body-weight relationship for children younger than 100 days, Table IV V1 row). The canonical PNA register stores PNA in months; the source paper's 100-day threshold is converted internally as 100 / 30.4375. Source range encompasses neonates with PNA 1-25 days (median 8 days) through adults at 81 years (Table I; Methods 'Neonates (24)' subsection).",
      source_name        = "PNA"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 174L,
    n_studies      = 7L,
    age_range      = "1 day to 81 years",
    weight_range   = "0.68-122.7 kg",
    sex_female_pct = NA_real_,
    disease_state  = "Pooled paediatric and adult cohorts spanning preterm and term neonates, infants, toddlers, children, adolescents, and adults receiving propofol for procedural sedation, post-operative sedation, anaesthesia induction, or as an investigational pharmacokinetic study drug in healthy volunteers.",
    dose_range     = "Cohort-specific regimens (Methods): neonates 3 mg/kg IV bolus for chest-tube removal / placement or endotracheal intubation; infants 2-4 mg/kg/h IV infusion for post-craniofacial-surgery sedation (median 12.5 h); toddlers 4 mg/kg IV bolus before bathing of minor burns; children 3 mg/kg IV bolus, or 3.5 mg/kg loading dose followed by maintenance infusions of 0.15 or 0.20 then 0.125 mg/kg/min; adolescents 2-10 mg/kg/h IV infusion for scoliosis surgery (median 6.8 h); adult women 2.5 mg/kg IV bolus over 60 s for gynaecological-surgery induction; adult healthy volunteers IV bolus plus 60-min infusion at 25, 50, 100, or 200 mg/kg/min.",
    regions        = "Netherlands (pooled studies from Leiden University, Erasmus MC Sophia Children's Hospital, St. Antonius Hospital Nieuwegein, and other Dutch and Belgian centres).",
    notes          = "4,396 propofol plasma concentration observations contributed across the seven cohorts (Table I 'Samples c' column gives per-subject sample counts 4-21). Cohort sizes after exclusions (Table I): neonates 25, infants 20, toddlers 12, children 53, adolescents 14, female adults 24, healthy adult volunteers 24. Concentrations were logarithmically transformed and fitted simultaneously in NONMEM VI (FOCEI) because the concentration range across cohorts spanned more than 1,000 fold (Methods, 'Pharmacokinetic Modeling')."
  )

  ini({
    # Final PK model estimates (Wang 2012 Table IV). Time in minutes, dose in
    # mg, central-compartment concentration in mg/L (= ug/mL). The paper
    # writes the three-compartment parameterisation in terms of total
    # clearance CL, central volume V1, fast inter-compartmental clearance
    # Q2, fast peripheral volume V2, slow inter-compartmental clearance Q3,
    # and slow peripheral volume V3. The nlmixr2lib convention maps:
    #   V1  (paper, central)            -> vc
    #   V2  (paper, fast peripheral)    -> vp
    #   V3  (paper, slow peripheral)    -> vp2
    #   Q2  (paper, central <-> fast)   -> q
    #   Q3  (paper, central <-> slow)   -> q2
    lcl  <- log(2.02);  label("Clearance at 70 kg reference body weight (L/min)")              # Wang 2012 Table IV Final model: Cl_p = 2.02 (CV% 2.6)
    lvc  <- log(7.58);  label("Central volume of distribution V1 (L)")                         # Wang 2012 Table IV Final model: V1_p = 7.58 (CV% 12.4)
    lq   <- log(1.77);  label("Fast inter-compartmental clearance Q2 (L/min)")                 # Wang 2012 Table IV Final model: Q2 = 1.77 (CV% 6.3)
    lvp  <- log(15.5);  label("Fast peripheral volume of distribution V2 at 70 kg (L)")        # Wang 2012 Table IV Final model: V2_p = 15.5 (CV% 14.7)
    lq2  <- log(1.65);  label("Slow inter-compartmental clearance Q3 at 70 kg (L/min)")        # Wang 2012 Table IV Final model: Q3_p = 1.65 (CV% 5.0)
    lvp2 <- log(221);   label("Slow peripheral volume of distribution V3 at 70 kg (L)")        # Wang 2012 Table IV Final model: V3_p = 221 (CV% 8.9)

    # Bodyweight-dependent allometric exponent (BDE) parameters for the WT
    # effect on CL (Wang 2012 Eq. 5; Table IV Final model). The exponent
    # k = k0 - kmax * WT^hill / (k50^hill + WT^hill) declines sigmoidally
    # from k0 at WT = 0 kg to k0 - kmax at large WT, with the half-decrease
    # body weight k50 and the Hill coefficient governing steepness.
    # k0 and the Hill coefficient were fixed in the Final model because
    # the covariance step did not succeed when they were estimated; the
    # values come from the successful run without covariance (Wang 2012
    # Results paragraph on Model IV).
    k0_cl   <- fixed(1.34);     label("BDE exponent of WT/70 on CL at WT = 0 kg (unitless)")   # Wang 2012 Table IV Final model: k0 = 1.34 FIX (bootstrap mean 1.35, CV% 6.2)
    kmax_cl <- 0.79;            label("BDE maximum decrease in WT/70 exponent on CL (unitless)") # Wang 2012 Table IV Final model: kmax = 0.79 (CV% 12.2)
    k50_cl  <- 3.78;            label("BDE body weight at which half of maximum decrement is reached (kg)") # Wang 2012 Table IV Final model: k50 = 3.78 (CV% 15.1)
    hill_cl <- fixed(5.24);     label("BDE Hill coefficient (steepness of sigmoidal decline; unitless)") # Wang 2012 Table IV Final model: gamma = 5.24 FIX (bootstrap mean 5.25, CV% 41.6)

    # Body-weight exponent on the fast peripheral volume V2 (Wang 2012
    # Table IV Final model footnote: V2 = V2_p * (BW/70)^m). The slow
    # peripheral volume V3 and the slow inter-compartmental clearance Q3
    # both use a linear (BW/70) scaling -- exponent fixed at 1 -- so they
    # do not carry a separate estimated exponent parameter here; their
    # scaling is encoded directly in model().
    e_wt_vp <- 0.55;            label("Allometric exponent of WT/70 on V2 (unitless)")          # Wang 2012 Table IV Final model: m = 0.55 (CV% 17.5)

    # Inter-individual variability (log-normal eta on log-scale parameters).
    # Wang 2012 Table IV Final model reports IIV directly as the variance
    # omega^2 of the eta normal distribution (Methods, Statistical Model:
    # theta_i = theta_TV * exp(eta_i), eta_i ~ N(0, omega^2)). IIV was
    # retained on CL, V1, V2, V3, and Q3 only; no IIV on Q2 (Table IV
    # leaves the Q2 row blank under 'Inter-individual variability').
    etalcl  ~ 0.09        # Wang 2012 Table IV Final model: omega^2 CL = 0.09 (CV% 18.0)
    etalvc  ~ 1.19        # Wang 2012 Table IV Final model: omega^2 V1 = 1.19 (CV% 41.3)
    etalvp  ~ 0.52        # Wang 2012 Table IV Final model: omega^2 V2 = 0.52 (CV% 40.8)
    etalvp2 ~ 0.71        # Wang 2012 Table IV Final model: omega^2 V3 = 0.71 (CV% 44.0)
    etalq2  ~ 0.25        # Wang 2012 Table IV Final model: omega^2 Q3 = 0.25 (CV% 17.0)

    # Residual error. Wang 2012 fits the residual error as additive on
    # log-transformed concentrations (Methods, Statistical Model Eq.
    # following Eq. for log-normal IIV): log(C_ij) = log(C_pred,ij) +
    # eps_ij with eps ~ N(0, sigma^2). The text explicitly notes that
    # this 'corresponds to proportional error on untransformed data',
    # so the equivalent linear-scale proportional SD is sqrt(sigma^2) =
    # sqrt(0.06) = 0.2449.
    propSd  <- sqrt(0.06);  label("Proportional residual error on linear concentration (fraction)")  # Wang 2012 Table IV Final model: sigma^2 = 0.06 (CV% 10.3); propSd = sqrt(0.06) = 0.2449
  })

  model({
    # Bodyweight-dependent allometric exponent for CL (Wang 2012 Eq. 5).
    # The exponent declines sigmoidally with body weight from k0_cl at
    # WT = 0 kg to (k0_cl - kmax_cl) at large WT, passing through the
    # half-decrement at WT = k50_cl.
    k_cl <- k0_cl - kmax_cl * WT^hill_cl / (k50_cl^hill_cl + WT^hill_cl)

    # Central-volume covariate model (Wang 2012 Table IV V1 row): V1 scales
    # linearly with (WT/70) for subjects with postnatal age below 100 days
    # and is constant otherwise. The canonical PNA register carries PNA
    # in months, so the source paper's 100-day threshold is converted
    # internally as 100 / 30.4375 (mean days per month).
    v1_scale <- ifelse(PNA < 100 / 30.4375, WT / 70, 1)

    # Individual PK parameters. Allometric (WT/70)^k_cl on CL with the
    # bodyweight-dependent exponent; linear (WT/70) on Q3 and V3; power
    # (WT/70)^m on V2; conditional (WT/70) on V1 via v1_scale; no
    # body-weight scaling on Q2 (Wang 2012 Table IV).
    cl  <- exp(lcl  + etalcl)  * (WT / 70)^k_cl
    vc  <- exp(lvc  + etalvc)  * v1_scale
    q   <- exp(lq)
    vp  <- exp(lvp  + etalvp)  * (WT / 70)^e_wt_vp
    q2  <- exp(lq2  + etalq2)  * (WT / 70)
    vp2 <- exp(lvp2 + etalvp2) * (WT / 70)

    # Three-compartment IV PK with first-order elimination from the central
    # compartment (Wang 2012 Methods 'Structural Model': three-compartment
    # model with central CL, central V1, fast Q2/V2, and slow Q3/V3).
    # Dose enters `central` directly (intravenous bolus or infusion) via
    # the cmt column of the user data set.
    d/dt(central)     <-  q  / vp  * peripheral1 + q2 / vp2 * peripheral2 -
                          (cl + q + q2) / vc * central
    d/dt(peripheral1) <-  q  / vc  * central     - q  / vp  * peripheral1
    d/dt(peripheral2) <-  q2 / vc  * central     - q2 / vp2 * peripheral2

    # Propofol plasma / whole-blood concentration in the central
    # compartment. Dose units mg, vc units L -> Cc units mg/L (= ug/mL).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
