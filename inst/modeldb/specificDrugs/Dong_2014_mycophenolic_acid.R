Dong_2014_mycophenolic_acid <- function() {
  description <- "Population PK-PD model for oral mycophenolic acid (MPA, the active moiety of mycophenolate mofetil MMF) in paediatric renal transplant recipients in the early post-transplant period (Dong 2014). Two-compartment disposition with a Savic 2007-style 8-transit-compartment absorption chain followed by a first-order absorption step from depot to central; dose-dependent relative bioavailability described by a power function of dose per body surface area (DBSA) with reference 450 mg/m^2; estimated body-weight exponent of 0.31 on CL/F (not the canonical allometric 0.75). The PD layer links MPA plasma concentration to inosine monophosphate dehydrogenase (IMPDH) activity in peripheral blood mononuclear cells via a simplified inhibitory Emax model with Emax fixed at 0 (i.e., complete inhibition achievable in the limit of high MPA concentration)."
  reference <- paste(
    "Dong M, Fukuda T, Cox S, de Vries MT, Hooper DK, Goebel J, Vinks AA.",
    "Population pharmacokinetic-pharmacodynamic modelling of mycophenolic",
    "acid in paediatric renal transplant recipients in the early",
    "post-transplant period.",
    "Br J Clin Pharmacol. 2014;78(5):1102-1112.",
    "doi:10.1111/bcp.12426.",
    sep = " "
  )
  vignette <- "Dong_2014_mycophenolic_acid"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at the per-subject value listed in Dong 2014 Table 1 (paediatric renal transplant cohort; mean 39.8 kg, median 38.2 kg, range 10.3-106.4 kg). Used for size-scaling on apparent clearance with reference 70 kg and an ESTIMATED exponent of 0.31 (not the canonical allometric 0.75; Dong 2014 Discussion notes that the 0.75 exponent could not be confirmed in this small paediatric cohort and a data-driven exponent of 0.31 was retained). Weight was tested but not retained on Vc/F, Vp/F, or Q/F.",
      source_name        = "WT"
    ),
    BSA = list(
      description        = "Body surface area (assumed DuBois or Mosteller formula; not specified in source)",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at the per-subject value listed in Dong 2014 Table 1 (PK cohort mean 1.21 m^2, median 1.26 m^2, range 0.49-2.21 m^2). Used as the denominator of dose per body surface area (DBSA = DOSE / BSA, mg/m^2), which enters the relative bioavailability power function with reference 450 mg/m^2 (the lower of the two starting-dose-per-BSA protocols). Source paper does not state which BSA formula was used; assume DuBois or Mosteller for downstream simulations.",
      source_name        = "BSA"
    ),
    DOSE = list(
      description        = "Current administered mycophenolate mofetil (MMF) dose at the most recent dosing event",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-dose-record value carried forward between doses. Used together with BSA to derive the dose per body surface area DBSA = DOSE / BSA in mg/m^2, which drives the dose-dependent relative bioavailability factor BIO = (DBSA / 450)^e_dbsa_f with reference 450 mg/m^2 (the lower of the two starting-dose protocols in Dong 2014). Subjects were started pre-surgery on 450 or 600 mg/m^2 MMF twice daily according to each institutional protocol and subsequently titrated; the population mean DBSA on the PK study day was 444 mg/m^2 (range 244.6-589.6 mg/m^2; Dong 2014 Table 1). At simulation time, set DOSE equal to the prescribed MMF amount in mg at each dose record (the rxode2 amt event-column value is the same number on dose rows).",
      source_name        = "DOSE"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age in years",
      units       = "years",
      type        = "continuous",
      notes       = "Tested in the covariate stepwise selection but did not significantly improve the model (Dong 2014 Results 'Population PK modelling' and Discussion). PK study cohort age range 2.1-20.2 years, median 14.2 years (Dong 2014 Table 1)."
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Tested as a popPK covariate but not retained (Dong 2014 Results)."
    ),
    CRCL = list(
      description = "Creatinine clearance calculated by the Schwartz formula (paediatric; BSA-normalized)",
      units       = "mL/min/1.73 m^2",
      type        = "continuous",
      notes       = "Tested as a popPK covariate but not retained (Dong 2014 Methods 'Patients' and Results 'Population PK modelling'). PK study cohort range 20.5-228.3 mL/min/1.73 m^2, mean 118.1 (Dong 2014 Table 1)."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/dL",
      type        = "continuous",
      notes       = "Tested on CL/F (not significant) and on PD baseline IMPDH activity E0 (delta-OFV = -4.26, P > 0.01) but not retained because it did not meet the required significance threshold (Dong 2014 Results 'Population PK modelling' and 'Population PK-PD modelling')."
    ),
    HGB = list(
      description = "Hemoglobin",
      units       = "g/dL",
      type        = "continuous",
      notes       = "Reported in Table 1 baseline demographics; not retained in the final popPK model."
    ),
    RACE_BLACK = list(
      description        = "African-American race indicator (1 = African-American, 0 = otherwise; Caucasian in this cohort)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Caucasian; 20 of 24 subjects in the PK cohort, 15 of 17 in the PK-PD cohort)",
      notes              = "Tested on PD baseline IMPDH activity E0 with a statistically significant OFV drop (delta-OFV = -8.84) consistent with lower IMPDH baseline activity in the two African-American patients (mean 2.13 vs 3.86 nmol h-1 mg-1 protein in Caucasians). However, the bootstrap 95% confidence interval (-1.09, +0.50) indicated an unreliable covariate effect given only n = 2 African-American patients, and race was therefore not included in the final PK-PD model (Dong 2014 Results 'Population PK-PD modelling')."
    ),
    CONMED_THYMOGLOBULIN = list(
      description = "Co-medication with rabbit antithymocyte globulin (Thymoglobulin) for induction immunosuppression",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested but not retained (Dong 2014 Results 'Population PK modelling')."
    ),
    CONMED_BASILIXIMAB = list(
      description = "Co-medication with basiliximab (Simulect) for induction immunosuppression",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested but not retained (Dong 2014 Results 'Population PK modelling')."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 24L,
    n_studies      = 1L,
    age_range      = "2.1-20.2 years (PK cohort; PK-PD subset 4.1-20.2 years)",
    age_median     = "14.2 years (PK cohort); 14.7 years (PK-PD subset)",
    weight_range   = "10.3-106.4 kg (PK cohort; PK-PD subset 10.3-106.0 kg)",
    weight_median  = "38.2 kg (PK cohort); 43.3 kg (PK-PD subset)",
    sex_female_pct = 37.5,
    race_ethnicity = c(Caucasian = 83, African_American = 17),
    disease_state  = "Paediatric kidney transplant recipients in the early post-transplant period (days 4-9 post-transplantation). All patients were on tacrolimus + prednisone backbone immunosuppression with basiliximab (Simulect) or rabbit antithymocyte globulin (Thymoglobulin) as induction therapy; patients on concomitant cyclosporine were excluded.",
    dose_range     = "Oral MMF (CellCept) at 450 or 600 mg/m^2 twice daily at study entry per institutional protocol, with subsequent dose adjustments at the prescribing physician's discretion. Mean dose per BSA on the PK study day was 444 mg/m^2 (range 244.6-589.6 mg/m^2; Dong 2014 Table 1). For patients weighing less than 15 kg a sparser PK sampling schedule was used due to maximum blood draw volume restrictions.",
    regions        = "Cincinnati, OH, USA (Cincinnati Children's Hospital Medical Center; single-site).",
    notes          = "Pooled PK cohort n = 24 with 214 MPA plasma concentrations (none below the 0.25 mg/L LLOQ). PK-PD subset n = 17 with 97 IMPDH activity measurements available from peripheral blood mononuclear cells; the other 7 patients had insufficient sample volume for the IMPDH assay. A pre-transplant baseline IMPDH activity was collected once before transplantation. Demographics from Dong 2014 Table 1."
  )

  ini({
    # ========================================================================
    # Structural PK parameters for the two-compartment disposition with a
    # Savic 2007-style 8-transit-compartment absorption chain feeding a depot
    # that empties at rate ka into central. Apparent values (scaled by oral
    # bioavailability F). Reference body weight 70 kg.
    # All structural values from Dong 2014 Table 2 final estimates.
    # ========================================================================
    lcl  <- log(22.0);     label("Apparent oral clearance CL/F at 70 kg (L/h)")                                # Dong 2014 Table 2 row 'CL/F (L/h/70 kg)' theta1 = 22.0
    lvc  <- log(45.4);     label("Apparent central volume of distribution Vc/F (L)")                           # Dong 2014 Table 2 row 'Vc/F (L) = 45.4'
    lvp  <- log(411);      label("Apparent peripheral volume of distribution Vp/F (L)")                        # Dong 2014 Table 2 row 'Vp/F (L) = 411'
    lq   <- log(22.4);     label("Apparent inter-compartmental clearance Q/F (L/h)")                           # Dong 2014 Table 2 row 'Q/F (L/h) = 22.4'
    lka  <- log(2.5);      label("First-order absorption rate constant from depot to central (1/h)")           # Dong 2014 Table 2 row 'Ka (h-1) = 2.5'
    lmtt <- log(0.25);     label("Mean absorption transit-chain time MTT (h)")                                 # Dong 2014 Table 2 row 'MTT (h) = 0.25'
    lnn  <- fixed(log(8)); label("Number of transit compartments n (FIXED at 8)")                              # Dong 2014 Table 2 row 'Number of transit compartment (n) = 8 fix'

    # ========================================================================
    # Allometric exponent on CL/F. The estimated value 0.31 is data-driven
    # rather than the canonical theoretical 0.75; the source authors note in
    # the Discussion that the 0.75 law could not be confirmed in this small
    # paediatric cohort and other paediatric MPA studies have reported
    # similarly low exponents.
    # ========================================================================
    e_wt_cl  <- 0.31;      label("Body weight exponent on CL/F (unitless)")                                    # Dong 2014 Table 2 row 'theta2 (BW exponent on CL/F) = 0.31'

    # ========================================================================
    # Dose-dependent relative bioavailability: BIO = theta3 * (DBSA/450)^theta4
    # with theta3 fixed to 1, so the bioavailability factor is the power
    # function (DBSA/450)^e_dbsa_f. The negative exponent -0.43 captures the
    # decrease in F with increasing MMF dose per BSA, consistent with the
    # adult value of -0.41 reported by de Winter et al. (Dong 2014 Discussion).
    # ========================================================================
    e_dbsa_f <- -0.43;     label("Exponent of dose per body surface area on relative bioavailability (unitless; BIO = (DBSA/450)^e_dbsa_f)")  # Dong 2014 Table 2 row 'theta4 (DBSA exponent on F) = -0.43'

    # ========================================================================
    # Inter-individual variability. The paper estimated IIV on CL/F, Ka, and
    # MTT; IIV on Vc/F, Vp/F, Q/F and on the number of transit compartments n
    # was fixed to zero during model building (Dong 2014 Results, 'Removal of
    # inter-individual variability on V_c/F, V_p/F and Q/F provided more
    # stable models and did not significantly compromise the model fit'). The
    # source reports %CV; convert to log-scale variance via
    #   omega^2 = log(1 + CV^2)
    # consistent with Dong 2014 Table 2 footnote c.
    # ========================================================================
    etalcl  ~ 0.0649  # Dong 2014 Table 2 omega(CL/F) = 25.9% CV; omega^2 = log(1 + 0.259^2) = 0.0649
    etalka  ~ 2.300   # Dong 2014 Table 2 omega(Ka)   = 299.6% CV; omega^2 = log(1 + 2.996^2) = 2.300
    etalmtt ~ 1.130   # Dong 2014 Table 2 omega(MTT)  = 144.8% CV; omega^2 = log(1 + 1.448^2) = 1.130

    # ========================================================================
    # PD parameters of the simplified inhibitory Emax model linking MPA
    # plasma concentration Cc to IMPDH activity:
    #     impdh = E0 * EC50 / (EC50 + Cc)
    # The asymptotic minimum activity Emax was fixed to 0 in the paper (its
    # estimated value was unreliable, RSE > 50%), so complete inhibition is
    # achievable in the limit of high MPA concentration (Dong 2014 Results
    # 'Population PK-PD modelling' and Table 3).
    # ========================================================================
    lrbase  <- log(3.45);  label("Baseline IMPDH activity E0 (nmol h-1 mg-1 protein)")                         # Dong 2014 Table 3 row 'E0 (nmol h-1 mg-1 protein) = 3.45'
    lec50   <- log(1.73);  label("MPA concentration for 50% inhibition of IMPDH activity EC50 (mg/L)")        # Dong 2014 Table 3 row 'EC50 (mg/L) = 1.73'

    # IIV on PD parameters
    etalrbase ~ 0.1459  # Dong 2014 Table 3 omega(E0)   = 39.6% CV; omega^2 = log(1 + 0.396^2) = 0.1459
    etalec50  ~ 0.4226  # Dong 2014 Table 3 omega(EC50) = 72.5% CV; omega^2 = log(1 + 0.725^2) = 0.4226

    # ========================================================================
    # Residual error.
    # - PK: the paper fit log-transformed concentrations with an additive
    #   error model (Methods 'Population PK-PD modelling'); the NONMEM
    #   variance sigma corresponds to the SD on the log scale, which is the
    #   nlmixr2 `expSd` for `lnorm`. Per Table 2 footnote c, the reported
    #   "Additive residual error (%CV) = 51.0" is calculated as
    #   100*sqrt(sigma) where sigma is the NONMEM variance, giving
    #   sigma = (51.0/100)^2 = 0.2601 and expSd = sqrt(0.2601) = 0.51 (log
    #   scale; equivalent linear-space lognormal CV approximately
    #   sqrt(exp(0.2601) - 1) = 54.5%).
    # - PD: proportional residual only (Dong 2014 Table 3 row 'Proportional
    #   residual error (%CV) = 42.2'). Per Table 3 footnote c, %CV =
    #   100*sqrt(sigma) with sigma = NONMEM variance of the proportional
    #   error, so propSd_impdh = sqrt(0.422^2) = 0.422 on the linear scale.
    # ========================================================================
    expSd        <- 0.51;   label("MPA PK residual SD on log-transformed concentration (lognormal)")          # Dong 2014 Table 2 row 'Additive residual error (%CV) = 51.0'; per footnote c sigma = (51.0/100)^2 = 0.2601 and expSd = sqrt(sigma) = 0.51
    propSd_impdh <- 0.422;  label("IMPDH activity proportional residual SD (fraction)")                       # Dong 2014 Table 3 row 'Proportional residual error (%CV) = 42.2'; per footnote c sigma = (42.2/100)^2 = 0.1781 and propSd = sqrt(sigma) = 0.422
  })

  model({
    # --------------------------------------------------------------------------
    # 1. Derived covariate-effect terms
    # --------------------------------------------------------------------------
    # Dose per body surface area at the current dose event (mg/m^2). The DOSE
    # covariate column carries the current MMF dose amount in mg forward
    # between dose events, matching the per-record value used in the Dong 2014
    # NONMEM control stream.
    dbsa <- DOSE / BSA

    # Relative bioavailability power function: BIO = (DBSA/450)^e_dbsa_f.
    # The lead coefficient theta3 was fixed to 1 in the paper so the
    # bioavailability factor is just the power function. Reference DBSA is
    # 450 mg/m^2, the lower of the two starting-dose-per-BSA protocols.
    bio_factor <- (dbsa / 450) ^ e_dbsa_f

    # --------------------------------------------------------------------------
    # 2. Individual PK parameters (allometric size scaling on CL only)
    # --------------------------------------------------------------------------
    cl  <- exp(lcl + etalcl) * (WT / 70) ^ e_wt_cl
    vc  <- exp(lvc)
    vp  <- exp(lvp)
    q   <- exp(lq)
    ka  <- exp(lka + etalka)
    mtt <- exp(lmtt + etalmtt)
    nn  <- exp(lnn)

    # Individual PD parameters
    e0   <- exp(lrbase + etalrbase)
    ec50 <- exp(lec50  + etalec50)

    # --------------------------------------------------------------------------
    # 3. Micro-constants for the two-compartment disposition
    # --------------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # --------------------------------------------------------------------------
    # 4. ODE system. Savic 2007-style analytical transit-compartment chain
    # feeding the depot; depot empties at rate ka into the two-compartment
    # disposition. rxode2's transit() built-in evaluates the analytical
    # (Stirling-approximated) gamma input function; f(depot) <- 0 suppresses
    # the bolus into depot so transit() delivers the full dose via the chain.
    # The bioavailability factor enters transit() as the bio argument because
    # nominal F is fixed at 1 in the paper.
    # --------------------------------------------------------------------------
    d/dt(depot)       <- transit(nn, mtt, bio_factor) - ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    f(depot) <- 0

    # --------------------------------------------------------------------------
    # 5. Observations and residual error
    # --------------------------------------------------------------------------
    # Plasma MPA concentration (mg / L = ug / mL = paper units mg/L).
    Cc <- central / vc
    Cc ~ lnorm(expSd)

    # IMPDH activity (nmol h-1 mg-1 protein); simplified inhibitory Emax with
    # Emax (the asymptotic minimum at C = infinity) fixed at 0, so the
    # paper's E = E0 * (1 - C/(EC50 + C)) simplifies algebraically to
    # E0 * EC50 / (EC50 + Cc).
    impdh <- e0 * ec50 / (ec50 + Cc)
    impdh ~ prop(propSd_impdh)
  })
}
