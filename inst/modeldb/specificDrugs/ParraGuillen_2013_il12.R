ParraGuillen_2013_il12 <- function() {
  description <- paste(
    "Preclinical (mouse, female C57BL/6 with subcutaneous MC38 tumor).",
    "Applicability re-fit of the Parra-Guillen 2013 semi-mechanistic K-PD",
    "tumor-growth-dynamics model to a single dose of murine IL-12 delivered",
    "by hydrodynamic plasmid injection. Structural equations are identical",
    "to the CyaA-E7 build (ParraGuillen_2013_cyaaE7); cell-line and",
    "immunotherapeutic-kinetics parameters (Ts0, lambda, k1, REG50) are",
    "re-estimated against the Medina-Echeverz 2014 dataset, while vaccine-",
    "efficacy and regulator-shape parameters (k3, k4, gamma) and the mixture",
    "probability P(1) = 0.844 are carried fixed from the CyaA-E7 fit."
  )
  reference <- paste(
    "Parra-Guillen ZP, Berraondo P, Grenier E, Ribba B, Troconiz IF.",
    "Mathematical Model Approach to Describe Tumour Response in Mice After",
    "Vaccine Administration and its Applicability to Immune-Stimulatory",
    "Cytokine-Based Strategies. The AAPS Journal 2013;15(3):797-807.",
    "doi:10.1208/s12248-013-9483-5. IL-12 cohort: Medina-Echeverz J,",
    "Berraondo P et al. (2014, J Immunol; cited as reference 25 of the",
    "Parra-Guillen 2013 paper)."
  )
  vignette <- "ParraGuillen_2013_tumor_immunotherapy"
  paper_specific_compartments <- c("vac", "tran", "svac", "reg", "tumor_size")
  units <- list(time = "day", dosing = "(arbitrary unit, set to 1 at plasmid injection)", concentration = "(K-PD, no PK)")

  covariateData <- list(
    MIX_VAC_RELAPSE = list(
      description        = paste(
        "Per-subject binary mixture-model class indicator (carried from the",
        "CyaA-E7 fit). 1 = relapser subpopulation (transient inhibitory",
        "signal, k2 = k1); 0 = responder / cure subpopulation (permanent",
        "inhibitory signal, k2 = 0 FIX). Time-fixed per subject."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (responder / cure subpopulation)",
      notes              = paste(
        "Latent-class assignment from the NONMEM $MIXTURE block (Parra-Guillen",
        "2013 Methods Data Evaluation bullet (e); Discussion p. 803). The",
        "mixture probability P(1) = 0.844 FIX was carried over from the",
        "CyaA-E7 fit when re-fitting IL-12 because the paper did not re-",
        "estimate P(1) for the Medina-Echeverz cohort (Table I row 'P(1)',",
        "IL-12 column). For typical-value simulation set MIX_VAC_RELAPSE = 0",
        "(cure) or 1 (relapser); for population simulation draw",
        "MIX_VAC_RELAPSE ~ Bernoulli(0.156) per subject."
      ),
      source_name        = "MIXTURE (NONMEM $MIXTURE assignment; component 2 = relapser)"
    )
  )

  population <- list(
    species        = "mouse (female C57BL/6, 5 weeks old; subcutaneous MC38 tumor)",
    n_subjects     = 33L,
    n_studies      = 1L,
    age_range      = "5 weeks at tumor cell inoculation",
    weight_range   = NA_character_,
    sex_female_pct = 100,
    disease_state  = paste(
      "Subcutaneous MC38 tumor model (5 x 10^5 cells subcutaneously injected",
      "on day 0; MC38 = murine colon adenocarcinoma). Tumor size measured as",
      "the mean of two perpendicular diameters; limit of quantification 2 mm",
      "(treated as censored M3 likelihood)."
    ),
    dose_range     = paste(
      "Single hydrodynamic injection of 10 ug of plasmid encoding murine",
      "IL-12 on day 23 after tumor cell inoculation (n = 21). PBS control",
      "group n = 12. Hydrodynamic injection delivers the plasmid to",
      "hepatocytes where the cytokine is expressed and released systemically."
    ),
    regions        = "Preclinical (Medina-Echeverz 2014 dataset; CIMA, Pamplona, Spain)",
    notes          = paste(
      "Cohort from Medina-Echeverz et al. (2014, cited as reference 25 of",
      "Parra-Guillen 2013). Tumor-size data were log-transformed before",
      "fitting; BQL observations (<= 2 mm) treated as censored. Parameters",
      "k3 (vaccine-efficacy second-order rate), k4 (regulator dynamics),",
      "gamma (Hill steepness), and the mixture probability P(1) were held",
      "fixed at the CyaA-E7 estimates because the IL-12 cohort (n = 33,",
      "one administration time) does not support their re-estimation."
    )
  )

  ini({
    # Initial tumor size at day 0 (Table I row 'Ts0', IL-12 column). The MC38
    # cell-line model is fit with a very small Ts0 (~1e-6 mm), i.e. tumor is
    # essentially absent at cell inoculation and grows linearly from zero.
    lrbase  <- log(1.16e-6); label("Initial tumor size Ts0 at cell inoculation (mm)")                 # Table I row 'Ts0' (IL-12 column)

    # Linear (zero-order) tumor growth rate -- re-estimated for MC38
    llam    <- log(0.335);   label("Zero-order tumor growth rate lambda (mm / day)")                  # Table I row 'lambda' (IL-12 column)

    # Shared first-order rate constant for plasmid kinetics, transit, and
    # relapser-population SVAC turnover -- re-estimated for IL-12.
    lk1     <- log(0.189);   label("First-order rate constant k1 for plasmid kinetics (1 / day)")    # Table I row 'k1' (IL-12 column)

    # k3, k4, gamma, and P(1) were held FIXED at the CyaA-E7 estimates when
    # fitting the IL-12 cohort because the MC38 dataset does not support
    # their identification (Methods 'Applicability of the Mathematical Model'
    # subsection p. 800: 'estimating those parameters dependent upon the
    # tumour cell line (Ts0, lambda and REG50) or upon the immunotherapeutic
    # agent (k1) together with the corresponding inter-animal variability').
    lk3     <- fixed(log(1.08));    label("Vaccine efficacy second-order rate constant k3 (1 / day) -- fixed from CyaA-E7 fit")   # Table I row 'k3' (IL-12 column = '1.08 FIX')
    lk4     <- fixed(log(0.0390));  label("First-order regulator dynamics rate k4 (1 / day) -- fixed from CyaA-E7 fit")           # Table I row 'k4' (IL-12 column = '0.0390 FIX')
    lhill   <- fixed(log(5.24));    label("Hill steepness gamma for REG inhibition (unitless) -- fixed from CyaA-E7 fit")         # Table I row 'gamma' (IL-12 column = '5.24 FIX')

    # Regulator half-inhibition amount -- re-estimated for MC38
    lreg50  <- log(2.08);    label("Regulator amount at half inhibition REG50 (mm)")                  # Table I row 'REG50' (IL-12 column)

    # IAV -- re-estimated for IL-12 cohort
    etallam   ~ 0.037249   # IAV 19.3 % (Table I row 'lambda', IL-12):   0.193^2
    etalreg50 ~ 0.130321   # IAV 36.1 % (Table I row 'REG50',  IL-12):   0.361^2

    # Residual error: log-additive (NONMEM 'additive on the logarithmic
    # domain of the transformed data' -> nlmixr2 lnorm())
    expSd  <- 0.168;   label("Log-additive residual SD on log(tumor_size) (Log(mm))")                  # Table I row 'Residual error [Log (mm)]' (IL-12 column)
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

    # 2. Mixture-gated SVAC degradation rate (Eq. 3). For responders
    # (MIX_VAC_RELAPSE = 0) k2 = 0; for relapsers (MIX_VAC_RELAPSE = 1)
    # k2 = k1.
    k2 <- MIX_VAC_RELAPSE * k1

    # 3. Hill-function inhibition of vaccine efficacy by the regulator (Eq. 4
    # bracketed factor)
    inh    <- reg50^hill / (reg50^hill + reg^hill)

    # 4. ODE system (Eqs. 1-5 of Parra-Guillen 2013, IL-12 applicability fit)
    d/dt(vac)        <- -k1 * vac                              # Eq. 1
    d/dt(tran)       <-  k1 * vac  - k1 * tran                 # Eq. 2
    d/dt(svac)       <-  k1 * tran - k2 * svac                 # Eq. 3
    d/dt(tumor_size) <-  lam - k3 * inh * tumor_size * svac    # Eq. 4
    d/dt(reg)        <-  k4 * tumor_size - k4 * reg            # Eq. 5

    # 5. Initial conditions at tumor cell inoculation
    tumor_size(0) <- rbase

    # 6. Observation: tumor size in mm (log-additive residual)
    tumor_size ~ lnorm(expSd)
  })
}
attr(ParraGuillen_2013_il12, "message") <- paste(
  "IL-12 plasmid dose enters the vac compartment as a unit bolus (amt = 1, cmt = vac).",
  "MIX_VAC_RELAPSE (0 = cure / responder, 1 = relapser; population probability of 1",
  "is 1 - P(1) = 0.156, fixed from CyaA-E7) must be supplied per subject. Parameters",
  "k3 / k4 / gamma / P(1) are FIXED at the CyaA-E7 fit values; only Ts0 / lambda /",
  "k1 / REG50 (+ IAV) are re-estimated for the Medina-Echeverz MC38 cohort."
)
ParraGuillen_2013_il12
