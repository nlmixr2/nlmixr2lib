Winchell_2024_posaconazole <- function() {
  description <- "One-compartment population PK model for posaconazole intravenous solution and powder for oral suspension (PFS) in pediatric patients aged 2 to 17 years with documented or expected neutropenia (Winchell 2024). First-order absorption from the depot compartment for the PFS formulation; IV doses enter the central compartment directly. Clearance and central volume are allometrically scaled by body weight with estimated (not fixed) exponents. Relative bioavailability of the PFS formulation is estimated on the logit scale so it is constrained to (0, 1), and carries a large logit-domain interindividual variability. No demographic or clinical covariate (age, weight beyond allometry, eGFR, sex, ethnicity) and no food effect was retained, so the base structural model is the final model. Residual variability is additive on the log scale (lnorm)."
  reference   <- "Winchell G, de Greef R, Ouerdani A, Fauchet F, Wrishko RE, Mangin E, Bruno C, Waskin H. A population pharmacokinetic model for posaconazole intravenous solution and oral powder for suspension formulations in pediatric patients with neutropenia. Antimicrob Agents Chemother. 2024;68(4):e01197-23. doi:10.1128/aac.01197-23."
  vignette    <- "Winchell_2024_posaconazole"
  units       <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot   = list(analyte = "posaconazole", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "posaconazole", units = "mg", specimen = "plasma",                   verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject in the source analysis. Enters both CL and Vc as an allometric power effect on (WT / 28.6), with both exponents estimated rather than fixed (Winchell 2024 Table 2: 0.624 for CL, 0.971 for Vc). The normalisation constant 28.6 kg is the overall median body weight of the analysis population (Winchell 2024 Table 1, All N = 114 column); the paper does not print the reference weight used in the NONMEM covariate model, so it was back-solved against the paper's own Supplementary Tables 2 and 3 -- see the vignette Assumptions and deviations section. Body weight was also screened as an ordinary covariate in the stepwise analysis (beyond the allometric terms) and was not retained.",
      source_name        = "Weight"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 114L,
    n_studies      = 1L,
    n_observations = "1,236 plasma posaconazole PK observations (Winchell 2024 Results 'Participant characteristics'). 80 samples were excluded: 33 from 2 participants with incomplete oral PFS dose intakes, 31 with no available sampling time, and 13 with a duplicate sampling time.",
    age_range      = "2-17 years (median 8; 2-<7 years subgroup median 3, 7-17 years subgroup median 13)",
    weight_range   = "10.2-102 kg (median 28.6; 2-<7 years subgroup median 16, range 10.2-41.7; 7-17 years subgroup median 45.4, range 18.2-102)",
    sex_female_pct = 41,
    race_ethnicity = "83% White; 87% not Hispanic (Winchell 2024 Results 'Participant characteristics')",
    disease_state  = "Immunocompromised pediatric patients aged 2 to 17 years with documented or anticipated neutropenia in the setting of acute leukemia, myelodysplasia, severe aplastic anemia, autologous HSCT, high-risk neuroblastoma, advanced-stage non-Hodgkin lymphoma, allogeneic HSCT during the pre-engraftment (neutropenic) period, or hemophagocytic lymphohistiocytosis. Posaconazole was given for prophylaxis of invasive fungal disease.",
    dose_range     = "3.5, 4.5, or 6 mg/kg posaconazole, IV BID on day 1 then IV QD on days 2-10, followed by either oral PFS or IV QD on days 10-28; absolute dose capped at 300 mg. Cohort sizes: 3.5 mg/kg n = 35, 4.5 mg/kg n = 31, 6 mg/kg n = 48; IV only n = 54, IV and PFS n = 60.",
    regions        = "Not reported (multicenter phase 1b study P097 / MK-5592-097, ClinicalTrials.gov NCT02452034)",
    renal_function = "Median eGFR 146 mL/min/1.73 m2 (range 12.2-314); eGFR was screened as a covariate and not retained.",
    notes          = "Fit in NONMEM 7.2 with FOCE and an additive residual-error model on log-transformed concentrations. Observations below the 5.00 ng/mL LLOQ were excluded rather than modelled. Shrinkage was low for CL (5%) but substantial for Vc (34%) and for the logit-scale bioavailability random effect (42%); the paper's own simulations (Supplementary Tables 2 and 3) behave as if the bioavailability random effect were considerably smaller than the reported omega, which is consistent with that 42% shrinkage -- see the vignette Assumptions and deviations section. Food intake within 2 h before to 1 h after a PFS dose (as a binary meal indicator and as light / medium / heavy categories) was tested on bioavailability and showed no significant effect, supporting administration of PFS with or without food."
  )

  ini({
    # ==================================================================
    # Structural PK parameters (Winchell 2024 Table 2, "Parameter
    # estimates for the final PK model"). Reference subject: body
    # weight 28.6 kg (the overall median of the analysis population,
    # Winchell 2024 Table 1 "All (N = 114)" column). The paper does not
    # state the allometric reference weight; see the model() block and
    # the vignette Assumptions and deviations section for the
    # back-solve that establishes it.
    # ==================================================================
    lcl <- log(4.71);  label("Clearance CL at the reference body weight (L/h)")                  # Winchell 2024 Table 2: CL = 4.71 L/h, RSE 3.86%
    lvc <- log(112);   label("Central volume of distribution Vc at the reference body weight (L)")  # Winchell 2024 Table 2: Vc = 112 L, RSE 5.18%
    lka <- log(0.212); label("First-order absorption rate constant ka for the PFS formulation (1/h)")  # Winchell 2024 Table 2: KA = 0.212 1/h, RSE 17.9%

    # Relative bioavailability of the oral PFS formulation, estimated on
    # the logit scale so it is constrained to (0, 1) -- Winchell 2024
    # Results "Final model": "estimated bioavailability using a logit
    # function (in order to constrain its value between 0 and 1)".
    # qlogis(0.826) = 1.55754.
    logitfdepot <- qlogis(0.826); label("Logit of the relative bioavailability F1 of the oral PFS formulation (unitless)")  # Winchell 2024 Table 2: F1 = 0.826, RSE 5.58%

    # Allometric exponents on body weight. Both are ESTIMATED in the
    # source analysis (RSE reported), not fixed at the theoretical
    # 0.75 / 1 values -- so they are deliberately not wrapped in
    # fixed(). Winchell 2024 Results: "The estimate of allometric
    # exponents for body weight on CL (0.624) and Vc (0.971)
    # approximated theoretical values (0.75 and 1, respectively)."
    e_wt_cl <- 0.624; label("Allometric exponent of body weight on CL (unitless; estimated)")  # Winchell 2024 Table 2: alpha for CL = 0.624, RSE 9.86%
    e_wt_vc <- 0.971; label("Allometric exponent of body weight on Vc (unitless; estimated)")  # Winchell 2024 Table 2: alpha for Vc = 0.971, RSE 7.86%

    # ==================================================================
    # Interindividual variability.
    #
    # Winchell 2024 Table 2 footnote b defines the IIV column for CL and
    # Vc as "a CV%, calculated by the square root of Omega, multiplied by
    # 100" -- i.e. the tabulated 37.1 and 27.7 are 100 * omega (the
    # log-scale SD), NOT lognormal CV% values requiring the
    # omega^2 = log(1 + CV^2) conversion. Variances are therefore
    # omega^2 = 0.371^2 and 0.277^2 directly.
    #
    # Winchell 2024 Results "Final model" states there is a correlation
    # between the CL and Vc random effects, but neither Table 2 nor the
    # supplement reports its magnitude. The block structure is retained
    # with the off-diagonal set to zero rather than inventing a
    # covariance; this is recorded in the vignette Assumptions and
    # deviations section.
    # ==================================================================
    etalcl + etalvc ~ c(0.137641,
                        0, 0.076729)

    # Random effect on bioavailability, applied on the LOGIT scale.
    # Winchell 2024 Table 2 footnote c: "The variability in F1 is shown
    # as a standard deviation in the logit domain" = 2.02, "This
    # corresponds to 95% of patients having F1 comprising between 8% and
    # 99%". Taking 2.02 as the logit-domain SD reproduces that interval
    # exactly: expit(1.55754 -/+ 1.96 * 2.02) = (0.083, 0.996) = 8% to
    # 99%, which confirms the SD (not variance) reading. Variance is
    # therefore 2.02^2 = 4.0804.
    etalogitfdepot ~ 4.0804  # Winchell 2024 Table 2: IIV (F1) = 2.02 as a logit-domain SD, RSE 18.9%

    # ==================================================================
    # Residual error -- additive on log-transformed concentrations
    # (Winchell 2024 Results "Final model": "an additive error model in
    # the logarithmic scale"; Methods: "an additive model for residual
    # variability on log-transformed data"). Encoded as ~ lnorm(expSd)
    # so the SD applies directly on log(Cc), matching the sibling
    # posaconazole model vanIersel_2018_posaconazole.R.
    # ==================================================================
    expSd <- 0.331; label("Log-scale (additive-on-log) residual standard deviation on Cc")  # Winchell 2024 Table 2: Residual error, SD = 0.331, RSE 4.71%
  })

  model({
    # ==================================================================
    # 1. Reference constant
    # ==================================================================
    # Allometric normalisation weight (kg). Winchell 2024 does not print
    # the reference weight used in the NONMEM covariate model. 28.6 kg is
    # the overall median body weight of the analysis population
    # (Winchell 2024 Table 1, "All (N = 114)" column) and is the value
    # supported by back-solving the paper's own simulation outputs:
    # Supplementary Tables 2 and 3 bracket the reference weight at
    # roughly 29-35 kg and exclude the conventional 70 kg by a wide
    # margin (70 kg overpredicts the Supplementary Table 2 IV Cavg
    # geometric means by 54-78%). Normalising to the analysis-population
    # median also matches the approach used by the same authors in the
    # adult posaconazole tablet model that this analysis was based on
    # (vanIersel_2018_posaconazole.R, Winchell 2024 reference 18).
    wt_ref <- 28.6

    # ==================================================================
    # 2. Individual PK parameters -- allometric scaling on body weight
    # ==================================================================
    wt_ratio <- WT / wt_ref

    cl <- exp(lcl + etalcl) * wt_ratio^e_wt_cl
    vc <- exp(lvc + etalvc) * wt_ratio^e_wt_vc
    ka <- exp(lka)

    # Relative bioavailability of the oral PFS formulation: the random
    # effect is added on the logit scale, then transformed back to (0, 1).
    fdepot <- expit(logitfdepot + etalogitfdepot)

    # ==================================================================
    # 3. Micro-constants
    # ==================================================================
    kel <- cl / vc

    # ==================================================================
    # 4. ODE system -- one compartment with first-order absorption.
    # Oral PFS dose records target cmt = depot; IV solution dose records
    # target cmt = central directly ("with IV administration being
    # assumed to be directed into the central compartment", Winchell 2024
    # Methods "Population PK model development").
    # ==================================================================
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ==================================================================
    # 5. Bioavailability -- applies to the oral PFS depot dose only; the
    # IV solution enters central with an implicit F of 1, so F1 here is
    # the bioavailability of PFS relative to the IV solution.
    # ==================================================================
    f(depot) <- fdepot

    # ==================================================================
    # 6. Observation and residual error
    # ==================================================================
    # Cc (ng/mL) = central (mg) / vc (L) * 1e6 (mg to ng) / 1e3 (L to mL)
    #            = central / vc * 1000.
    Cc <- (central / vc) * 1000
    Cc ~ lnorm(expSd)
  })
}
