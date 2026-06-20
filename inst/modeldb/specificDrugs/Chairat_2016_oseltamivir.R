Chairat_2016_oseltamivir <- function() {
  description <- paste(
    "Joint population pharmacokinetic model for oral oseltamivir (parent)",
    "and its active antiviral metabolite oseltamivir carboxylate in 12 obese",
    "(BMI >= 30 kg/m^2) and 12 non-obese (BMI < 30 kg/m^2) healthy Thai adult",
    "volunteers (Chairat 2016 BJCP). First-order absorption (ka) into a",
    "one-compartment parent (OS) disposition, an intermediate metabolism",
    "compartment delaying carboxylate appearance (rate km), and a",
    "one-compartment oseltamivir carboxylate (OC) disposition. Relative oral",
    "bioavailability F is fixed to unity with interindividual variability on",
    "F absorbing absorption differences. Creatinine clearance computed using",
    "Janmahasatian fat-free mass (CLCR(FFM); raw Cockcroft-Gault with FFM",
    "substituted for total body weight) is a linear covariate on CL/FOC,",
    "centred at the population median CLCR(FFM) of 73 mL/min (3.84% increase",
    "per 10 mL/min increase). Obesity itself was not a retained covariate in",
    "the formal model. Residual error is additive on log-transformed",
    "concentrations of OS and OC (encoded here as a log-normal residual on",
    "Cc and Cc_oc).")
  reference <- "Chairat K, Jittamala P, Hanpithakpong W, Day NPJ, White NJ, Pukrittayakamee S, Tarning J. Population pharmacokinetics of oseltamivir and oseltamivir carboxylate in obese and non-obese volunteers. Br J Clin Pharmacol. 2016;81(6):1103-1112. doi:10.1111/bcp.12892"
  vignette <- "Chairat_2016_oseltamivir"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  # Intermediate metabolism delay compartment between parent OS central and
  # metabolite OC central; not an absorption-chain transit (Savic) and not a
  # well-stirred-hepatic compartment (Standing 2012 transit_oc). Declared as
  # paper-specific so checkModelConventions() does not flag the name.
  paper_specific_compartments <- "metabolism"

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance computed by the Cockcroft-Gault formula with fat-free mass (FFM, Janmahasatian 2005) substituted for total body weight; NOT BSA-normalised. Chairat 2016 Methods Eq. 9 (FFM) and the CG variant used as CLCR(FFM).",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear covariate on apparent oseltamivir carboxylate clearance CL/FOC: CL/FOC = exp(lcl_oc) * (1 + e_crcl_cl_oc * (CRCL - 73)). The 0.00384 per mL/min coefficient is 3.84% per 10 mL/min, i.e. Chairat 2016 Table 1 covariate-effect row and Eq. 11. Centering 73 mL/min is the population median CLCR(FFM) used by the paper; the cohort range was 48.0-114 mL/min (Chairat 2016 Discussion paragraph 3).",
      source_name        = "CLCR(FFM)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 24L,
    n_studies      = 1L,
    age_range      = "18-60 years",
    weight_range   = "Non-obese median BMI 22.2 kg/m^2 (range 18.8-24.2); obese median BMI 33.8 kg/m^2 (range 30.8-43.2); pooled BMI range 18.8-43.2 kg/m^2.",
    sex_female_pct = NA_real_,
    race_ethnicity = c(Asian = 100),
    disease_state  = "Healthy adult Thai volunteers; 12 obese (BMI >= 30 kg/m^2) and 12 non-obese (BMI < 30 kg/m^2). Volunteers eligible if healthy male or non-pregnant female aged 18-60 years.",
    dose_range     = "Single oral 75 mg and 150 mg oseltamivir (fasted) in a randomised cross-over design with a 7-day washout between visits.",
    regions        = "Thailand (Hospital for Tropical Diseases, Faculty of Tropical Medicine, Mahidol University, Bangkok).",
    n_observations = "624 venous plasma samples total; 103 (16.5%) oseltamivir concentrations and 15 (2.40%) oseltamivir carboxylate concentrations were below the LLOQ (1 ng/mL OS, 10 ng/mL OC). Beal M3 method was evaluated but the final model used omitted-LLOQ data.",
    sampling       = "Pre-dose and 0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 10, 12 and 24 h post-dose.",
    notes          = "Open-label, crossover, randomised PK study. ClinicalTrials.gov NCT01049763. Full clinical and NCA details in Jittamala et al. 2014 (Antimicrob Agents Chemother 58:1615-1621). Demographic baselines from Chairat 2016 Methods and Results."
  )

  ini({
    # Structural parameters (Chairat 2016 Table 1 population estimates).
    lka     <- log(2.81);  label("Absorption rate constant ka (1/h)")                              # Chairat 2016 Table 1: ka = 2.81 1/h
    lcl     <- log(585);   label("Oseltamivir apparent clearance CL/FOS (L/h)")                    # Chairat 2016 Table 1: CL/FOS = 585 L/h
    lvc     <- log(1110);  label("Oseltamivir apparent volume of distribution V/FOS (L)")          # Chairat 2016 Table 1: V/FOS = 1110 L
    lkm     <- log(2.13);  label("Metabolism rate constant km for OC formation (1/h)")             # Chairat 2016 Table 1: km = 2.13 1/h
    lcl_oc  <- log(20.6);  label("Oseltamivir carboxylate apparent clearance CL/FOC (L/h)")        # Chairat 2016 Table 1: CL/FOC = 20.6 L/h
    lvc_oc  <- log(159);   label("Oseltamivir carboxylate apparent volume of distribution V/FOC (L)") # Chairat 2016 Table 1: V/FOC = 159 L
    lfdepot <- fixed(log(1)); label("Relative oral bioavailability F (fixed to unity)")            # Chairat 2016 Table 1: F = 100% (fixed)

    # Covariate effect on CL/FOC (Chairat 2016 Eq. 11 and Table 1 row "Effect
    # of CLCR on CL/FOC (% change per 10 units of CLCR)" = 3.84). The
    # coefficient stored here is in /mL/min units (= 0.0384 per 10 mL/min).
    e_crcl_cl_oc <- 0.00384; label("Linear effect of CLCR on CL/FOC (per 1 mL/min, around 73 mL/min)")  # Chairat 2016 Table 1 / Eq. 11: 3.84% per 10 mL/min around 73 mL/min

    # Random effects (Chairat 2016 Table 1). The paper reports BSV (IIV) and
    # BOV (interoccasion variability, IOV) as %CV computed by NONMEM as
    # [exp(omega^2)-1]^(1/2) * 100. Convert back via omega^2 = log(1+CV^2):
    #   BSV F  = 17.6% -> 0.030505
    #   IOV ka = 98.7% -> 0.680216
    #   BSV CL/FOS = 16.6% -> 0.027183
    #   IOV V/FOS = 18.6% -> 0.034007
    #   IOV km = 43.2% -> 0.171042
    #   BSV V/FOC = 18.7% -> 0.034375
    # IOV components are encoded here as IIV because rxode2 / nlmixr2 popPK
    # validation simulations typically use one occasion per subject; the
    # deviation is called out in the vignette Assumptions and deviations.
    # No IIV or IOV was retained for CL/FOC in the final model (Table 1 shows
    # "-" for the CL/FOC IIV column).
    etalfdepot ~ 0.030505  # Chairat 2016 Table 1: BSV(F) = 17.6%CV
    etalka     ~ 0.680216  # Chairat 2016 Table 1: IOV(ka) = 98.7%CV; encoded as IIV
    etalcl     ~ 0.027183  # Chairat 2016 Table 1: BSV(CL/FOS) = 16.6%CV
    etalvc     ~ 0.034007  # Chairat 2016 Table 1: IOV(V/FOS) = 18.6%CV; encoded as IIV
    etalkm     ~ 0.171042  # Chairat 2016 Table 1: IOV(km) = 43.2%CV; encoded as IIV
    etalvc_oc  ~ 0.034375  # Chairat 2016 Table 1: BSV(V/FOC) = 18.7%CV

    # Residual error. The paper used additive errors on log-transformed
    # concentrations, which corresponds to a log-normal residual on the
    # natural concentration scale. Encoded here as ~ lnorm(<sd>) with the
    # log-scale SDs from Table 1.
    expSd    <- 0.431; label("Log-scale residual SD on oseltamivir Cc")                        # Chairat 2016 Table 1: additive on log(OS) = 0.431
    expSd_oc <- 0.161; label("Log-scale residual SD on oseltamivir carboxylate Cc_oc")         # Chairat 2016 Table 1: additive on log(OC) = 0.161
  })

  model({
    # Centering value for CLCR(FFM) -- population median per Chairat 2016
    # Eq. 11 reads CL/FOC = 20.6 * [1 + 0.0384 * (CLCR(FFM) - 73)/10].
    crcl_ref <- 73

    # Individual structural parameters (typical * IIV * covariate scaling).
    ka     <- exp(lka     + etalka)
    cl     <- exp(lcl     + etalcl)
    vc     <- exp(lvc     + etalvc)
    km     <- exp(lkm     + etalkm)
    cl_oc  <- exp(lcl_oc)                * (1 + e_crcl_cl_oc * (CRCL - crcl_ref))
    vc_oc  <- exp(lvc_oc  + etalvc_oc)
    fbio   <- exp(lfdepot + etalfdepot)

    # ODE system (Chairat 2016 Figure 1):
    #   depot --ka--> central (OS) --(cl/vc)--> metabolism --km--> central_oc
    #   (OC) --(cl_oc/vc_oc)--> eliminated
    # All states are amounts (mg). Volumes are L; clearances L/h; rates 1/h.
    d/dt(depot)      <- -ka * depot
    d/dt(central)    <-  ka * depot - (cl / vc) * central
    d/dt(metabolism) <-  (cl / vc) * central - km * metabolism
    d/dt(central_oc) <-  km * metabolism - (cl_oc / vc_oc) * central_oc
    f(depot)         <-  fbio

    # Concentrations in ng/mL (dose in mg, V in L gives mg/L = 1000 ng/mL).
    Cc    <- 1000 * central    / vc
    Cc_oc <- 1000 * central_oc / vc_oc

    Cc    ~ lnorm(expSd)
    Cc_oc ~ lnorm(expSd_oc)
  })
}
