Ngo_2023_rivaroxaban_human <- function() {
  description <- "Two-compartment population PK model for rivaroxaban (RIV) in humans with a carbamazepine (CBZ) drug-drug-interaction covariate (Ngo 2023). Parameters were NOT fitted to human data: the structural model and all estimates come from a rat popPK fit (Ngo 2020) and were extrapolated to a 60 kg human by simple allometry with exponent 1.00 for the volumes (Vc, Vp) and by liver-blood-flow (LBF) scaling for the clearances (CL, Q) following Ward and Smith 2004; the first-order absorption rate constant Ka was taken from the Mueck 2007 human rivaroxaban popPK model, and D2, Alag2, F1, the IIV magnitudes and the residual error were assumed unchanged between rat and human. Absorption is a parallel mixed zero- and first-order process: a fraction F1 of the dose enters a depot and is absorbed first-order with rate Ka (lymphatic uptake of the lipophilic drug), while the complementary fraction 1 - F1 enters the central compartment directly as a zero-order input of duration D2 after a lag Alag2 (solubility-limited transfer across the enterocyte capillary network and hepatic portal vein). Concomitant carbamazepine, a strong CYP3A4 / CYP2J2 and P-gp / BCRP inducer, acts through two fractional-change covariate terms: it increases apparent clearance CL/F by 211 % and prolongs the zero-order absorption time D2 by 33.9 %. Simulating 20 mg once-daily rivaroxaban with and without 900 mg/day carbamazepine reproduces the paper's predicted interaction: AUC falls by 52.3 % after the first dose and by 64.1 % at steady state, and Cmax falls by 41.0 % and 49.8 % respectively. (Ngo 2023 prints 68.5 % for the steady-state AUC decrease, but that does not follow from the AUC values in its own Table 2, which give 64.1 %; see the vignette Errata.)"
  reference   <- "Ngo LT, Yun H-y, Chae J-w. Application of the Population Pharmacokinetics Model-Based Approach to the Prediction of Drug-Drug Interaction between Rivaroxaban and Carbamazepine in Humans. Pharmaceuticals. 2023;16(5):684. doi:10.3390/ph16050684. The underlying rat popPK estimates (Table 1, column 'Estimated in Rats') were first published in Ngo LT, Yang S, Tran QT, Kim SK, Yun H, Chae J. Effects of Carbamazepine and Phenytoin on Pharmacokinetics and Pharmacodynamics of Rivaroxaban. Pharmaceutics. 2020;12(11):1040. doi:10.3390/pharmaceutics12111040. The human first-order absorption rate constant Ka = 0.97 1/h is inherited from Mueck W, Becka M, Kubitza D, Voith B, Zuehlsdorf M. Population model of the pharmacokinetics and pharmacodynamics of rivaroxaban - an oral, direct Factor Xa inhibitor - in healthy subjects. Int J Clin Pharmacol Ther. 2007;45(6):335-344."
  vignette    <- "Ngo_2023_rivaroxaban"
  units       <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    CONMED_CBZ = list(
      description        = "Concomitant carbamazepine indicator: 1 = subject pretreated with carbamazepine (900 mg/day in the simulated human scenario; 45 mg/kg twice daily for 6 consecutive days in the underlying rat study), 0 = rivaroxaban administered alone",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (rivaroxaban alone; the 'control group' of Ngo 2023)",
      notes              = "Time-fixed in this model: the underlying rat design pretreated the test group for 6 days before the rivaroxaban dose, and the human simulation assumes carbamazepine autoinduction has already reached steady state, so the indicator does not change within a subject. Enters the model twice, both as fractional-change (1 + coefficient * CONMED_CBZ) terms rather than log-multiplicative terms: on apparent clearance CL/F (Ngo 2023 Equation 1, coefficient CL_CBZ = 2.11) and on the zero-order absorption duration D2 (Ngo 2023 Equation 2, coefficient D2_CBZ = 0.339). The two effects have opposite consequences for exposure magnitude but reinforce each other in lowering Cmax: faster elimination plus slower absorption. Ngo 2023 also tested carbamazepine on the first-order absorption rate constant Ka but did not retain it (Ngo 2023 Section 4.2).",
      source_name        = "Pretreatment with CBZ (index variable; 'control group' vs 'test group')"
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "rivaroxaban", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "rivaroxaban", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "rivaroxaban", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 1000L,
    n_studies      = 1L,
    weight_median  = "60 kg (assumed reference human body weight for the allometric extrapolation; Ngo 2023 Section 4.3)",
    disease_state  = "Simulated healthy adults. No human subjects were studied: the 1000 replicates are a Monte-Carlo cohort generated from the extrapolated typical values and the rat-derived IIV (Ngo 2023 Section 2.3).",
    dose_range     = "20 mg oral rivaroxaban once daily, with or without 900 mg/day oral carbamazepine (Ngo 2023 Section 2.3). 20 mg/day is the FDA-labelled rivaroxaban dose for stroke risk reduction in nonvalvular atrial fibrillation; 900 mg/day sits inside the 800-1200 mg/day labelled carbamazepine maintenance range for epilepsy (Ngo 2023 Discussion).",
    regions        = "Not applicable (simulation). The underlying rat study was conducted at Chungnam National University, Daejeon, Republic of Korea.",
    notes          = "PARAMETER PROVENANCE: this model was never fitted to human concentration data. The source popPK fit used 12 Sprague-Dawley-type rats (n = 6 per group, control and carbamazepine-pretreated) dosed 3 mg/kg oral rivaroxaban with sampling at 0.25, 0.5, 1, 2, 4, 8, 10 and 24 h (Ngo 2023 Section 4.1, originally Ngo 2020); see the sibling model Ngo_2023_rivaroxaban_rat for that fit in its native per-kilogram form. Ngo 2023 verified the extrapolation by checking that the predicted human AUC and Cmax fall within two-fold of pooled clinical observations: predicted AUC 1291.7 vs observed mean 1847.5 ng*h/mL (range 1559.6-2244.9) and predicted Cmax 133.2 vs observed mean 219.7 ng/mL (range 189.9-291.6) after a single 20 mg fed dose; at steady state predicted 2157.5 vs observed 2619.2 ng*h/mL (2329.9-2908.4) and predicted 172.2 vs observed 341.5 ng/mL (273.5-409.5) (Ngo 2023 Discussion). The extrapolated CL/F of 9.03 L/h and Vc/F of 42.1 L also agree with the Mueck 2007 human popPK values of 9.2 L/h and 55.3 L. The unbound fraction differs markedly between species (1.3 % in rat vs about 6.5 % in human) but was deliberately NOT used as an allometric correction factor because no drug-specific guidance was available (Ngo 2023 Discussion, 'limitations')."
  )

  ini({
    # ====================================================================
    # Structural PK parameters -- Ngo 2023 Table 1, column 'Extrapolated
    # in Humans'. Every value in that column is the rat estimate carried
    # through one of three routes, all of which reproduce the printed
    # human number exactly from the printed rat number:
    #
    #   Volumes (Vc, Vp), Ngo 2023 Equation 4, simple allometry with
    #   exponent 1.00 and BW_human / BW_rat = 60 / 0.25:
    #       Vc/F : 0.701 L/kg * 60 kg                        =  42.06 L
    #       Vp/F : 5.60  L/kg * 60 kg                        = 336    L
    #
    #   Clearances (CL, Q), Ngo 2023 Equation 5, liver-blood-flow scaling
    #   of Ward and Smith 2004 with LBF_human / LBF_rat = 21 / 85
    #   mL/min/kg, then multiplied by the 60 kg human body weight:
    #       CL/F : 0.609 L/h/kg * (21/85) * 60 kg            =   9.03 L/h
    #       Q/F  : 0.665 L/h/kg * (21/85) * 60 kg            =   9.86 L/h
    #
    #   Everything else ("Same" in Table 1) is assumed species-invariant,
    #   except Ka which is replaced by a published human value.
    #
    # The RSE (%) column of Table 1 reports the precision of the ORIGINAL
    # RAT estimates (Ngo 2020); it does not describe uncertainty in the
    # extrapolated human values, which carry the additional and
    # unquantified error of the interspecies scaling itself.
    # ====================================================================
    lcl <- log(9.03);   label("Apparent clearance CL/F from the central compartment in the absence of carbamazepine (L/h)")  # Ngo 2023 Table 1: 0.609 L/h/kg in rats (RSE 36.9 %) -> 9.03 L/h in humans by LBF scaling
    lvc <- log(42.06);  label("Apparent central volume of distribution Vc/F (L)")                                            # Ngo 2023 Table 1: 0.701 L/kg in rats (RSE 44.2 %) -> 42.06 L in humans by simple allometry
    lq  <- log(9.86);   label("Apparent inter-compartmental clearance Q/F (L/h)")                                            # Ngo 2023 Table 1: 0.665 L/h/kg in rats (RSE 21.4 %) -> 9.86 L/h in humans by LBF scaling
    lvp <- log(336);    label("Apparent peripheral volume of distribution Vp/F (L)")                                         # Ngo 2023 Table 1: 5.60 L/kg in rats (RSE 38.1 %) -> 336 L in humans by simple allometry

    # Absorption. Ka is NOT the rat estimate (2.31 1/h): Ngo 2023 replaced
    # it with the human value from the Mueck 2007 rivaroxaban popPK model,
    # so it is an inherited literature constant rather than something this
    # analysis estimated -- hence fixed().
    lka     <- fixed(log(0.97)); label("First-order absorption rate constant Ka of the depot arm (1/h)")                     # Ngo 2023 Table 1, human column: Ka = 0.97 1/h taken from Mueck 2007 (reference [25]); the rat estimate was 2.31 1/h
    ld2     <- log(6.62);        label("Zero-order absorption duration D2 into the central compartment without carbamazepine (h)")  # Ngo 2023 Table 1 and Equation 2: D2 = 6.62 h in the control group, assumed identical in humans
    ltlag   <- log(0.501);       label("Lag time Alag2 before the zero-order absorption arm begins (h)")                     # Ngo 2023 Table 1: Alag2 = 0.501 h, assumed identical in humans
    lfdepot <- log(0.260);       label("Fraction F1 of the dose absorbed by the first-order (depot) arm (unitless)")         # Ngo 2023 Table 1: F1 = 0.260 (RSE 9.10 %), assumed identical in humans; the complementary 1 - F1 = 0.740 takes the zero-order route

    # ====================================================================
    # Carbamazepine drug-drug-interaction coefficients. Ngo 2023 codes
    # carbamazepine pretreatment as a 0/1 index variable entering as a
    # FRACTIONAL CHANGE, not as a log-multiplicative shift:
    #     Equation 1 / 3:  CL/F = TVCL * (1 + CL_CBZ)
    #     Equation 2:      D2   = TVD2 * (1 + D2_CBZ)
    # Ngo 2023 Section 4.3 states the carbamazepine effects in humans were
    # assumed identical to those measured in rats, so the same two
    # coefficients apply to both this model and its rat sibling.
    # ====================================================================
    e_conmed_cbz_cl <- 2.11;  label("Fractional increase in apparent clearance CL/F with concomitant carbamazepine (unitless)")            # Ngo 2023 Equation 1 and Section 2.1: CL_CBZ = 2.11 (a 211 % increase); check 9.03 * (1 + 2.11) = 28.1 L/h, the Table 1 test-group value
    e_conmed_cbz_d2 <- 0.339; label("Fractional increase in the zero-order absorption duration D2 with concomitant carbamazepine (unitless)")  # Ngo 2023 Equation 2 and Section 2.1: D2_CBZ = 0.339 (a 33.9 % increase); 6.62 * (1 + 0.339) = 8.86 h vs the 8.84 h printed in Table 1 -- see vignette Errata

    # ====================================================================
    # Inter-individual variability. Ngo 2023 Table 1 reports IIV as a
    # percentage and Section 2.3 states "the variability of each PK
    # parameter was assumed to be the same in humans and rats". The
    # percentages are read as coefficients of variation and converted to
    # the log-scale variance by the lognormal identity
    #     omega^2 = log(CV^2 + 1).
    # ====================================================================
    etalcl ~ log(0.490^2 + 1)  # Ngo 2023 Table 1: IIV for CL/F = 49.0 % CV (RSE 35.0 %) -> omega^2 = 0.2151
    etalvc ~ log(0.470^2 + 1)  # Ngo 2023 Table 1: IIV for Vc/F = 47.0 % CV (RSE 33.0 %) -> omega^2 = 0.1996

    # ====================================================================
    # Residual unexplained variability -- combined additive and
    # proportional, on the ng/mL scale of the observed plasma
    # concentrations. Assumed identical in humans and rats.
    # ====================================================================
    addSd  <- 13.6;  label("Additive residual error on plasma rivaroxaban concentration (ng/mL)")   # Ngo 2023 Table 1: additive error 13.6 ng/mL (RSE 35.8 %)
    propSd <- 0.232; label("Proportional residual error on plasma rivaroxaban concentration (fraction)")  # Ngo 2023 Table 1: proportional error 23.2 % (RSE 17.2 %)
  })

  model({
    # ====================================================================
    # Individual PK parameters. The carbamazepine terms are applied on the
    # NATURAL scale as (1 + coefficient * indicator), matching Ngo 2023
    # Equations 1-3 exactly; writing them as exp(l... + coef * CONMED_CBZ)
    # would be a log-multiplicative form the paper did not use and would
    # not reproduce the printed test-group values.
    # ====================================================================
    cl     <- exp(lcl + etalcl) * (1 + e_conmed_cbz_cl * CONMED_CBZ)
    vc     <- exp(lvc + etalvc)
    q      <- exp(lq)
    vp     <- exp(lvp)
    ka     <- exp(lka)
    d2     <- exp(ld2) * (1 + e_conmed_cbz_d2 * CONMED_CBZ)
    tlag   <- exp(ltlag)
    fdepot <- exp(lfdepot)

    # ====================================================================
    # Two-compartment disposition with parallel mixed-order absorption
    # (Ngo 2023 Figure 2). The oral dose leaves the GI tract by two routes
    # that both terminate in the central compartment:
    #
    #   depot   : fraction F1 = 0.260, absorbed first-order at rate Ka,
    #             representing lymphatic uptake of rivaroxaban carried in
    #             the lipid vehicle. No lag on this arm.
    #   central : fraction 1 - F1 = 0.740, delivered as a zero-order input
    #             of duration D2 beginning after the lag Alag2,
    #             representing solubility-limited transfer through the
    #             capillary network around the enterocytes and the hepatic
    #             portal vein.
    #
    # A user therefore encodes each oral administration as TWO dose
    # records at the same time: one with cmt = "depot" and one with
    # cmt = "central". The f() multipliers split the dose; dur(central)
    # imposes the zero-order duration and alag(central) the lag. Because
    # every disposition parameter is apparent (divided by the true
    # bioavailability F), F1 and 1 - F1 sum to 1 and describe only how the
    # absorbed dose is partitioned between the two routes, not how much of
    # the dose is absorbed.
    # ====================================================================
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - (cl / vc) * central -
                          (q / vc) * central + (q / vp) * peripheral1
    d/dt(peripheral1) <-  (q / vc) * central - (q / vp) * peripheral1

    f(depot)      <- fdepot
    f(central)    <- 1 - fdepot
    dur(central)  <- d2
    alag(central) <- tlag

    # Plasma rivaroxaban concentration. Doses are in mg and vc is in L, so
    # central / vc is mg/L = ug/mL; the factor 1000 converts to the ng/mL
    # units in which Ngo 2023 reports every concentration, Cmax and AUC.
    Cc <- (central / vc) * 1000

    Cc ~ add(addSd) + prop(propSd)
  })
}
