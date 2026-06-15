Landersdorfer_2012_vildagliptin <- function() {
  description <- "Mechanism-based population PK plus DPP-4 activity model for vildagliptin in patients with type 2 diabetes. Target-mediated drug disposition with capacity-limited slow-tight binding of vildagliptin to DPP-4 in plasma and tissue and partial hydrolysis of vildagliptin by DPP-4."
  reference <- "Landersdorfer CB, He YL, Jusko WJ. Mechanism-based population pharmacokinetic modelling in diabetes: vildagliptin as a tight binding inhibitor and substrate of dipeptidyl peptidase IV. Br J Clin Pharmacol. 2012;73(3):391-401. doi:10.1111/j.1365-2125.2011.04108.x"
  vignette <- "Landersdorfer_2012_vildagliptin"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  paper_specific_compartments <- c("complex_peripheral")

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 13L,
    n_studies      = 1L,
    age_range      = "37-64 years",
    age_median     = "53.5 years",
    weight_range   = "65-116 kg",
    weight_median  = "91 kg (mean)",
    height_range   = "148-183 cm",
    sex_female_pct = 53.8,
    disease_state  = "type 2 diabetes mellitus (washout of hypoglycaemic drugs for up to 4 weeks)",
    dose_range     = "10, 25, 100 mg oral vildagliptin BID for 28 days, four-way crossover with placebo",
    notes          = paste(
      "Demographics from Results section: 12 subjects completed all four periods,",
      "one subject completed only the 10 and 25 mg treatments; 6 male, 7 female.",
      "Bioavailability F (77.2%) was fixed from co-modelling with a separate",
      "intravenous-dose study (He YL 2007, ref [6] in the paper)."
    )
  )

  ini({
    # Structural PK parameters (typical values, log-transformed)
    lka1    <- log(1.26)         ; label("Gut1 -> Gut2 absorption rate constant k_a1 (1/h)")          # Table 1: k_a1 = 1.26 h^-1 (BSV 46%, SE 15%)
    lka2    <- log(1.05)         ; label("Gut2 -> central absorption rate constant k_a2 (1/h)")      # Table 1: k_a2 = 1.05 h^-1 (BSV 14%, SE 4%)
    lcl     <- log(36.4)         ; label("Non-saturable vildagliptin clearance CL (L/h)")            # Table 1: CL = 36.4 L/h (BSV 25%, SE 9%)
    lvc     <- log(22.2)         ; label("Central volume of distribution V_C (L)")                   # Table 1: V_C = 22.2 L (BSV 37%, SE 11%)
    lvp     <- log(97.3)         ; label("Peripheral (tissue) volume V_P (L)")                       # Table 1: V_P = 97.3 L (BSV 37%, SE 13%)
    lq      <- log(40.1)         ; label("Inter-compartmental clearance CL_ic (L/h)")                # Table 1: CL_ic = 40.1 L/h (BSV 34%, SE 11%)
    lfdepot <- fixed(log(0.772)) ; label("Oral bioavailability F (fraction)")                        # Table 1 footnote *: F fixed at 0.772 from co-modelling with i.v. data (He 2007)

    # TMDD slow-tight-binding parameters
    lkd     <- log(71.9)         ; label("Equilibrium dissociation constant K_d (nmol/L)")           # Table 1: K_d = 71.9 nmol/L (BSV 54%, SE 16%)
    lk2     <- log(23.4)         ; label("Conversion rate weak -> high-affinity complex k_2 (1/h)")  # Table 1: k_2 = 23.4 h^-1 (BSV 70%, SE 22%)
    lkoff   <- log(0.612)        ; label("Dissociation rate of vildagliptin from DPP-4 k_off (1/h)") # Table 1: k_off = 0.612 h^-1 (BSV 94%, SE 27%)
    lkdeg   <- log(0.110)        ; label("Hydrolysis rate of vildagliptin by DPP-4 k_deg (1/h)")     # Table 1: k_deg = 0.110 h^-1 (BSV 81%, SE 26%)

    # Total DPP-4 amounts and activity conversion factor
    lrmaxc  <- log(5.0)          ; label("Total DPP-4 in central compartment R_maxC (nmol)")         # Table 1: R_maxC = 5.0 nmol (BSV 12%, SE 4%)
    lrmaxp  <- log(13000)        ; label("Total DPP-4 in peripheral compartment R_maxP (nmol)")     # Table 1: R_maxP = 13.0 umol = 13000 nmol (BSV 64%, SE 23%)
    lcf1    <- log(2.80)         ; label("DPP-4-amount to activity conversion cf1 (mU/mL/min per nmol)") # Table 1: cf1 = 2.80 (BSV 17%, SE 5%)

    # Between-subject variability (diagonal only).
    # Paper "Between-subject variability model" subsection: full variance-covariance matrix on PK
    # parameters was estimated in S-ADAPT, but Table 1 reports only the diagonal as
    # BSV (%) = sqrt(omega^2) * 100 (paper convention). Off-diagonal covariances were not reported;
    # we use diagonal-only IIV here. See vignette Assumptions and deviations.
    etalka1   ~ 0.2116                                                                                # BSV k_a1 46% -> var = 0.46^2
    etalka2   ~ 0.0196                                                                                # BSV k_a2 14% -> var = 0.14^2
    etalcl    ~ 0.0625                                                                                # BSV CL 25% -> var = 0.25^2
    etalvc    ~ 0.1369                                                                                # BSV V_C 37% -> var = 0.37^2
    etalvp    ~ 0.1369                                                                                # BSV V_P 37%
    etalq     ~ 0.1156                                                                                # BSV CL_ic 34%
    etalkd    ~ 0.2916                                                                                # BSV K_d 54%
    etalk2    ~ 0.4900                                                                                # BSV k_2 70%
    etalkoff  ~ 0.8836                                                                                # BSV k_off 94%
    etalkdeg  ~ 0.6561                                                                                # BSV k_deg 81%
    etalrmaxc ~ 0.0144                                                                                # BSV R_maxC 12%
    etalrmaxp ~ 0.4096                                                                                # BSV R_maxP 64%
    etalcf1   ~ 0.0289                                                                                # BSV cf1 17%

    # Residual error: combined additive + proportional for both outputs (paper "Residual error model").
    propSd      <- 0.487 ; label("Vildagliptin proportional residual error (fraction)")               # Table 1: CV_Vilda = 48.7%
    addSd       <- 0.99  ; label("Vildagliptin additive residual error (ng/mL)")                      # Table 1: SD_Vilda = 0.99 ng/mL
    propSd_DPP4 <- 0.196 ; label("DPP-4 activity proportional residual error (fraction)")             # Table 1: CV_DPP-4 = 19.6%
    addSd_DPP4  <- 0.061 ; label("DPP-4 activity additive residual error (mU/mL/min)")                # Table 1: SD_DPP-4 = 0.061
  })

  model({
    # Vildagliptin molecular weight (PubChem CID 6918537, C17H25N3O2 = 303.40 g/mol).
    # Used to convert the user-facing mg dosing to internal nmol amounts and to convert
    # internal nM concentrations back to ng/mL for the bioanalytical readout.
    mw_vilda <- 303.40

    # Individual parameters
    ka1   <- exp(lka1   + etalka1)
    ka2   <- exp(lka2   + etalka2)
    cl    <- exp(lcl    + etalcl)
    vc    <- exp(lvc    + etalvc)
    vp    <- exp(lvp    + etalvp)
    q     <- exp(lq     + etalq)
    fdep  <- exp(lfdepot)
    kd    <- exp(lkd    + etalkd)
    k2    <- exp(lk2    + etalk2)
    koff  <- exp(lkoff  + etalkoff)
    kdeg  <- exp(lkdeg  + etalkdeg)
    rmaxc <- exp(lrmaxc + etalrmaxc)
    rmaxp <- exp(lrmaxp + etalrmaxp)
    cf1   <- exp(lcf1   + etalcf1)

    # Vildagliptin concentrations (nmol/L = nM) in the two disposition compartments
    Cc_nM <- central     / vc
    Cp_nM <- peripheral1 / vp

    # Free DPP-4 amounts (nmol): total minus the drug-bound complex in each compartment
    freeR_c <- rmaxc - complex
    freeR_p <- rmaxp - complex_peripheral

    # Slow-tight binding rate (Michaelis-Menten form, paper Methods Eq for V_maxC / V_maxP):
    #   V_max = k_2 * C_drug * (R_max - DR) / (K_d + C_drug)
    # Units: (1/h) * (nmol/L) * (nmol) / (nmol/L) = nmol/h, i.e. amount of complex formed per hour.
    v_maxc <- k2 * Cc_nM * freeR_c / (kd + Cc_nM)
    v_maxp <- k2 * Cp_nM * freeR_p / (kd + Cp_nM)

    # ODE system (all amounts in nmol; volumes in L). Paper Methods "Structural models".
    d/dt(depot)              <- -ka1 * depot
    d/dt(transit1)           <-  ka1 * depot - ka2 * transit1
    d/dt(central)            <-  ka2 * transit1 - (cl / vc) * central - (q / vc) * central + (q / vp) * peripheral1 - v_maxc + koff * complex
    d/dt(peripheral1)        <-  (q / vc) * central - (q / vp) * peripheral1 - v_maxp + koff * complex_peripheral
    d/dt(complex)            <-  v_maxc - koff * complex - kdeg * complex
    d/dt(complex_peripheral) <-  v_maxp - koff * complex_peripheral - kdeg * complex_peripheral

    # Bioavailability and mg -> nmol dose-unit conversion applied at the depot input.
    # User passes amt in mg; bioavailable amount (nmol) = amt * F * 1e6 / MW (g/mol).
    # Worked example: 100 mg * 0.772 * 1e6 / 303.40 = 254452 nmol enters depot.
    f(depot) <- fdep * 1e6 / mw_vilda

    # Initial conditions: at baseline (before dosing) no drug is bound, so free DPP-4 = R_max
    # in each compartment. Baseline DPP-4 activity in plasma = cf1 * R_maxC (typical value
    # cf1 * R_maxC = 2.80 * 5.0 = 14.0 mU/mL/min, consistent with placebo baseline in Figure 3).

    # Observations
    #   Cc   - vildagliptin plasma concentration (ng/mL); convert nM by * MW / 1000
    #   DPP4 - DPP-4 activity in plasma (mU/mL/min); cf1 acts on free DPP-4 amount in central
    Cc   <- Cc_nM * mw_vilda / 1000
    DPP4 <- cf1 * freeR_c

    Cc   ~ prop(propSd)      + add(addSd)
    DPP4 ~ prop(propSd_DPP4) + add(addSd_DPP4)
  })
}
