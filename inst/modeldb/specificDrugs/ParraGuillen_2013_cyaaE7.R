ParraGuillen_2013_cyaaE7 <- function() {
  description <- paste(
    "Preclinical (mouse, female C57BL/6 with subcutaneous TC1 tumor expressing HPV E7).",
    "Semi-mechanistic K-PD tumor-growth-dynamics model of single-dose CyaA-E7 cancer",
    "vaccine: a virtual vaccine compartment feeds a two-compartment transit chain to a",
    "vaccine-elicited inhibitory signal SVAC that reduces tumor size via a second-order",
    "k3 * SVAC * tumor_size term, inhibited by a Hill-function regulator REG driven by",
    "tumor size; a binary mixture covariate MIX_VAC_RELAPSE gates the SVAC degradation",
    "rate (k2 = 0 for cure, k2 = k1 for relapse)."
  )
  reference <- paste(
    "Parra-Guillen ZP, Berraondo P, Grenier E, Ribba B, Troconiz IF.",
    "Mathematical Model Approach to Describe Tumour Response in Mice After Vaccine",
    "Administration and its Applicability to Immune-Stimulatory Cytokine-Based",
    "Strategies. The AAPS Journal 2013;15(3):797-807. doi:10.1208/s12248-013-9483-5."
  )
  vignette <- "ParraGuillen_2013_tumor_immunotherapy"
  paper_specific_compartments <- c("vac", "tran", "svac", "reg", "tumor_size")
  units <- list(time = "day", dosing = "(arbitrary unit, set to 1 at vaccine injection)", concentration = "(K-PD, no PK)")

  covariateData <- list(
    MIX_VAC_RELAPSE = list(
      description        = paste(
        "Per-subject binary mixture-model class indicator. 1 = relapser",
        "subpopulation (transient vaccine-elicited inhibitory signal, SVAC",
        "degradation rate k2 = k1; tumor regrowth observed after initial",
        "response; 15.6 % of mice); 0 = responder / cure subpopulation",
        "(permanent vaccine-elicited inhibitory signal, k2 = 0 FIX; complete",
        "tumor regression maintained; 84.4 % of mice). Time-fixed per subject."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (responder / cure subpopulation; majority class at 84.4 %)",
      notes              = paste(
        "Latent class from a NONMEM $MIXTURE block in Parra-Guillen 2013",
        "(Methods Data Evaluation bullet (e); Discussion p. 803). Population",
        "probability of relapse 1 - P(1) = 0.156 was obtained from a",
        "sub-analysis of only the 4 / 7 / 11-day vaccine-administration",
        "groups (where tumor-size-induced resistance had not yet emerged)",
        "and was held fixed during the full-cohort fit. For typical-value",
        "simulation set MIX_VAC_RELAPSE = 0 (cure phenotype, Figure 3 light",
        "grey individual predictions) or 1 (relapser, Figure 3 dark grey).",
        "For population simulation, draw MIX_VAC_RELAPSE ~ Bernoulli(0.156)",
        "per subject."
      ),
      source_name        = "MIXTURE (NONMEM $MIXTURE assignment; component 2 = relapser)"
    )
  )

  population <- list(
    species        = "mouse (female C57BL/6, 5 weeks old; subcutaneous TC1 tumor expressing HPV E7)",
    n_subjects     = 93L,
    n_studies      = 1L,
    age_range      = "5 weeks at tumor cell inoculation",
    weight_range   = NA_character_,
    sex_female_pct = 100,
    disease_state  = paste(
      "Subcutaneous TC1 tumor model (5 x 10^5 cells injected into the shaved",
      "back on day 0; TC1 = mouse lung-epithelial cell line transformed with",
      "HPV-16 E6 / E7 and activated H-Ras). Tumor size measured as the mean",
      "of two perpendicular diameters; limit of quantification 2 mm (61.6 %",
      "of all tumor-size observations were BQL)."
    ),
    dose_range     = paste(
      "Single IV dose of 50 ug CyaA-E7 vaccine on one of day 4 (n = 15), 7",
      "(n = 16), 11 (n = 14), 18 (n = 15), 25 (n = 14), or 30 (n = 15) after",
      "tumor cell inoculation. Vehicle (PBS) control n = 19 on day 4. CyaA-E7",
      "= fusion of the adenylate cyclase of Bordetella pertussis to the HPV E7",
      "oncoprotein, targeting antigen to dendritic cell CD11b receptors."
    ),
    regions        = "Preclinical (Berraondo et al. 2007 dataset; Centre for Applied Medical Research / CIMA, Pamplona, Spain)",
    notes          = paste(
      "Training (model-building) dataset: 93 mice (6 vaccine-arm groups + 1",
      "control). External validation dataset (Parra-Guillen 2013 Figure 1",
      "but not used in this model fit): 34 additional mice (n = 23 PBS day 4,",
      "n = 11 CyaA-E7 day 25). Tumor sizes > 20 mm triggered euthanasia per",
      "the institutional animal-care protocol. Tumor-size data were",
      "log-transformed before fitting and BQL observations (Ts <= 2 mm) were",
      "treated as censored using the NONMEM Laplacian M3 likelihood."
    )
  )

  ini({
    # Initial tumor size at day 0 (Table I row 'Ts0')
    lrbase  <- log(0.324);   label("Initial tumor size Ts0 at cell inoculation (mm)")                # Table I row 'Ts0' (CyaA-E7 column)

    # Linear (zero-order) tumor growth rate (Table I row 'lambda')
    llam    <- log(0.354);   label("Zero-order tumor growth rate lambda (mm / day)")                  # Table I row 'lambda' (CyaA-E7 column)

    # Shared first-order rate constant for vaccine elimination, transit, and
    # relapser-population SVAC turnover (Methods last paragraph p. 799:
    # 'processes reflecting vaccine elimination and turnover of the vaccine
    # elicited inhibitory signal (SVAC) share the same parameter k1').
    lk1     <- log(0.0907);  label("First-order rate constant k1 for vaccine kinetics (1 / day)")    # Table I row 'k1' (CyaA-E7 column)

    # Vaccine efficacy second-order rate constant (Table I row 'k3')
    lk3     <- log(1.08);    label("Vaccine efficacy second-order rate constant k3 (1 / day)")        # Table I row 'k3' (CyaA-E7 column)

    # Regulator-compartment dynamics rate (Table I row 'k4')
    lk4     <- log(0.0390);  label("First-order regulator-compartment dynamics rate k4 (1 / day)")    # Table I row 'k4' (CyaA-E7 column)

    # Regulator amount at half-maximal inhibition of vaccine efficiency (Hill
    # EC50 analogue; units = mm because REG accumulates in tumor-size units
    # via Eq. 5 dREG / dt = k4 * (Ts - REG)).
    lreg50  <- log(3.18);    label("Regulator amount at half inhibition REG50 (mm)")                  # Table I row 'REG50' (CyaA-E7 column)

    # Hill steepness for the REG-driven inhibition of k3 (Eq. 4 sigmoidal
    # inhibition form: k3 * REG50^gamma / (REG50^gamma + REG^gamma)).
    lhill   <- log(5.24);    label("Hill steepness gamma for REG inhibition of k3 (unitless)")        # Table I row 'gamma' (CyaA-E7 column)

    # Inter-animal variability (IAV%) was modelled exponentially (Methods
    # 'Model Building and Selection'); the log-scale variance is therefore
    # (IAV / 100)^2. Only lambda and REG50 carry IAV in the final model
    # (Discussion 'inclusion of IAV in lambda and REG50 was supported, AIC
    # = -356 vs no IAV').
    etallam   ~ 0.010201   # IAV 10.1 % (Table I row 'lambda'):   0.101^2
    etalreg50 ~ 0.114244   # IAV 33.8 % (Table I row 'REG50'):    0.338^2

    # Residual error: additive on log-transformed tumor diameter (Methods
    # 'residual variability was modelled considering an additive error in the
    # logarithmic domain of the transformed data'). In nlmixr2 this maps to
    # lnorm() with expSd interpreted as the SD on the log scale.
    expSd  <- 0.206;   label("Log-additive residual SD on log(tumor_size) (Log(mm))")                  # Table I row 'Residual error [Log (mm)]' (CyaA-E7 column)
  })

  model({
    # 1. Individual typical-value parameters
    rbase  <- exp(lrbase)
    lam    <- exp(llam + etallam)
    k1     <- exp(lk1)
    k3     <- exp(lk3)
    k4     <- exp(lk4)
    reg50  <- exp(lreg50 + etalreg50)
    hill   <- exp(lhill)

    # 2. Mixture-gated SVAC degradation rate (Eq. 3 in the paper). For
    # responders (MIX_VAC_RELAPSE = 0) k2 = 0 so SVAC remains constant once
    # established -> cure. For relapsers (MIX_VAC_RELAPSE = 1) k2 = k1 so
    # SVAC decays at rate k1 -> tumor regrowth after vaccine clearance.
    k2 <- MIX_VAC_RELAPSE * k1

    # 3. Hill-function inhibition of vaccine efficacy by the regulator (Eq. 4
    # bracketed factor): k3_eff = k3 * REG50^gamma / (REG50^gamma + REG^gamma).
    inh    <- reg50^hill / (reg50^hill + reg^hill)

    # 4. ODE system (Eqs. 1-5 of Parra-Guillen 2013)
    d/dt(vac)        <- -k1 * vac                              # Eq. 1
    d/dt(tran)       <-  k1 * vac  - k1 * tran                 # Eq. 2
    d/dt(svac)       <-  k1 * tran - k2 * svac                 # Eq. 3
    d/dt(tumor_size) <-  lam - k3 * inh * tumor_size * svac    # Eq. 4
    d/dt(reg)        <-  k4 * tumor_size - k4 * reg            # Eq. 5

    # 5. Initial conditions at tumor cell inoculation (Methods p. 799 last
    # paragraph). vac is reset to 1 at vaccine injection time via an event-
    # table dose row (amt = 1, cmt = vac); tran / svac / reg start at 0 by
    # rxode2 default.
    tumor_size(0) <- rbase

    # 6. Observation: tumor size in mm. The paper transformed observations to
    # log(mm) and fit with an additive residual; nlmixr2's lnorm() error model
    # is equivalent (SD applies on the log scale).
    tumor_size ~ lnorm(expSd)
  })
}
attr(ParraGuillen_2013_cyaaE7, "message") <- paste(
  "Vaccine dose enters the vac compartment as a unit bolus (amt = 1, cmt = vac).",
  "MIX_VAC_RELAPSE (0 = cure / responder, 1 = relapser; population probability of 1",
  "is 1 - P(1) = 0.156, fixed) must be supplied per subject in the simulation data.",
  "Observation tumor_size is the mean of two perpendicular tumor diameters (mm)."
)
ParraGuillen_2013_cyaaE7
