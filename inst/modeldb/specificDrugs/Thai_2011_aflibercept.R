Thai_2011_aflibercept <- function() {
  description <- "Mechanism-based population PK model for free and bound aflibercept (anti-VEGF Fc-fusion 'trap' protein; VEGF-Trap) in healthy adult male subjects (Thai 2011 BJCP). Two-compartment disposition of free aflibercept with linear elimination from the central compartment plus Michaelis-Menten binding of free aflibercept to VEGF occurring in the peripheral (tissue) compartment, producing a one-compartment bound-aflibercept species that is eliminated by first-order internalisation (kint). This is the second Michaelis-Menten approximation of the TMDD model of Gibiansky et al. (irreversible binding), with the bound complex carried as an explicit state. The bound-aflibercept volume of distribution Vb is fixed equal to the central volume Vc for identifiability. Pooled data from two phase 1 single-dose IV-infusion studies in healthy males (1, 2, 4 mg/kg over 1 h). No covariates were tested or retained in the final model."
  reference <- "Thai HT, Veyrat-Follet C, Vivier N, Dubruc C, Sanderink G, Mentre F, Comets E. A mechanism-based model for the population pharmacokinetics of free and bound aflibercept in healthy subjects. Br J Clin Pharmacol. 2011;72(3):402-414. doi:10.1111/j.1365-2125.2011.04015.x"
  vignette <- "Thai_2011_aflibercept"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 56L,
    n_studies      = 2L,
    age_range      = NA_character_,
    weight_range   = NA_character_,
    sex_female_pct = 0,
    race_ethnicity = NA_character_,
    disease_state  = "Healthy male volunteers from two phase 1 single-dose IV-infusion studies.",
    dose_range     = "1, 2 or 4 mg/kg as a single 1-h IV infusion. The crossover SC arm (2 mg/kg) of study 2 was excluded from the population analysis.",
    regions        = "South Africa (study 1; Pharma-Ethics IRB) and Germany (study 2; Ethik-Kommission der Landesarztekammer Baden-Wuerttemberg).",
    notes          = "Total of 1476 concentrations were analysed (732 free aflibercept + 744 bound aflibercept), of which 242 (32.5%) bound concentrations were below the bound LOQ of 0.0314 ug eq/mL (handled in MONOLIX 3.1 with the SAEM censored-data extension). Bound aflibercept assay concentrations were converted by the authors to free-aflibercept-equivalent concentrations by multiplication by 0.717 (the ratio of MW of free to bound aflibercept), so all bound observations are in ug eq/mL. No covariates were tested because the population is uniformly healthy males. Demographics (precise age, weight, race) are not enumerated in the paper text or Table 1."
  )

  ini({
    # Structural parameters - typical values for the healthy adult subject
    # (Thai 2011 Table 1). All values are population point estimates from
    # the final mechanism-based MM-TMDD model.
    lcl    <- log(0.88);   label("Linear clearance of free aflibercept from central compartment (CL, L/day)")    # Table 1: CL = 0.88 L/day (RSE 4%)
    lvc    <- log(4.94);   label("Central volume of distribution of free aflibercept (Vp in paper notation, L)") # Table 1: Vp = 4.94 L (RSE 4%)
    lq     <- log(1.39);   label("Inter-compartmental clearance of free aflibercept (Q, L/day)")                 # Table 1: Q = 1.39 L/day (RSE 9%)
    lvp    <- log(2.33);   label("Peripheral / tissue volume of distribution of free aflibercept (Vt, L)")       # Table 1: Vt = 2.33 L (RSE 7%)
    lvmax  <- log(0.99);   label("Michaelis-Menten maximum binding capacity Vmax (mg/day)")                      # Table 1: Vmax = 0.99 mg/day (RSE 5%)
    lkm    <- log(2.91);   label("Michaelis-Menten half-saturation constant Km (ug/mL = mg/L)")                  # Table 1: Km = 2.91 ug/mL (RSE 11%)
    lkint  <- log(0.028);  label("First-order internalisation rate constant of bound aflibercept (kint, 1/day)") # Table 1: kint = 0.028 /day (RSE 5%)
    # Vb (volume of distribution of bound aflibercept) is structurally fixed
    # equal to Vc for identifiability (Thai 2011 Table 1 row 'Vb (l): 4.94
    # ( = Vp)' and main-text page 405: 'Vb was fixed to the population value
    # of Vp with no interindividual variability, instead of being estimated'),
    # and is derived inside model() as vb <- vc rather than carried as an
    # estimated/fixed ini() parameter.

    # IIV - exponential, diagonal (Thai 2011 Methods 'Statistical model').
    # The Table 1 column 'Interindividual variability w (%)' reports the
    # standard deviation of the random effect eta on the log scale (the
    # paper's 'w' symbol per Methods 'variance w^2 ... e.g. for CL, w^2_CL').
    # The internal variance is therefore omega^2 = (w/100)^2.
    etalcl   ~ 0.0784    # Table 1: w(CL) = 28.0%   -> omega^2 = 0.28^2  = 0.0784
    etalvc   ~ 0.07453   # Table 1: w(Vp) = 27.3%   -> omega^2 = 0.273^2 = 0.07453
    etalq    ~ 0.24800   # Table 1: w(Q)  = 49.8%   -> omega^2 = 0.498^2 = 0.24800
    etalvp   ~ 0.15840   # Table 1: w(Vt) = 39.8%   -> omega^2 = 0.398^2 = 0.15840
    etalvmax ~ 0.01850   # Table 1: w(Vmax) = 13.6% -> omega^2 = 0.136^2 = 0.01850
    etalkm   ~ 0.20794   # Table 1: w(Km) = 45.6%   -> omega^2 = 0.456^2 = 0.20794
    # No IIV retained on kint (Thai 2011 Table 1: w(kint) = '-'; main text
    # page 405: 'interindividual variability on the internalization rate
    # constant (kint) was found to be small (11.9%) and badly estimated.
    # The likelihood ratio test demonstrated that removing this variability
    # ... did not significantly change the fit').

    # Residual error - free aflibercept uses combined additive + proportional;
    # bound aflibercept uses proportional only (Thai 2011 Table 1 right block
    # and page 405: 'best residual error model of bound aflibercept was the
    # proportional model').
    propSd          <- 0.171;  label("Proportional residual error on free aflibercept (fraction)") # Table 1: sigma_p free = 17.1% (RSE 3%)
    addSd           <- 0.05;   label("Additive residual error on free aflibercept (ug/mL)")        # Table 1: sigma_a free = 0.05 ug/mL (RSE 9%)
    propSd_complex  <- 0.126;  label("Proportional residual error on bound aflibercept (fraction)") # Table 1: sigma_p bound = 12.6% (RSE 4%)
  })

  model({
    # Individual parameters
    cl   <- exp(lcl   + etalcl)
    vc   <- exp(lvc   + etalvc)
    q    <- exp(lq    + etalq)
    vp   <- exp(lvp   + etalvp)
    vmax <- exp(lvmax + etalvmax)
    km   <- exp(lkm   + etalkm)
    kint <- exp(lkint)             # no IIV; not retained
    vb   <- vc                     # paper: Vb fixed = Vp (= vc in canonical naming) for identifiability

    # Concentrations
    cfree_central <- central     / vc   # ug/mL (free aflibercept in plasma / central)
    cfree_tissue  <- peripheral1 / vp   # ug/mL (free aflibercept in peripheral / tissue)
    cbound        <- complex     / vb   # ug eq/mL (bound aflibercept; vb = vc by constraint)

    # Michaelis-Menten binding occurs in the peripheral (tissue) compartment
    # (Thai 2011 page 404, 'The nonlinear peripheral binding model was found
    # to describe the data better than the central binding one'). Units:
    # vmax in mg/day, km in mg/L, cfree_tissue in mg/L, so binding_rate in mg/day.
    binding_rate <- vmax * cfree_tissue / (km + cfree_tissue)

    # ODE system - Thai 2011 main text Eq. (final model) / Appendix Eq. A14.
    # Free aflibercept distributes between central and peripheral by ordinary
    # 2-compartment kinetics with linear elimination from central; the
    # peripheral state additionally feeds the bound-aflibercept compartment
    # through saturable MM binding; bound aflibercept is eliminated by
    # first-order internalisation kint.
    d/dt(central)     <-  q * (cfree_tissue - cfree_central) - cl * cfree_central
    d/dt(peripheral1) <-  q * (cfree_central - cfree_tissue) - binding_rate
    d/dt(complex)     <-  binding_rate - kint * complex

    # Observations
    Cc         <- cfree_central
    Cc_complex <- cbound

    Cc         ~ add(addSd) + prop(propSd)
    Cc_complex ~ prop(propSd_complex)
  })
}
