Hanley_2024_brigatinib <- function() {
  description <- paste(
    "Two-compartment oral population pharmacokinetic reduction of the",
    "Simcyp minimal-PBPK-with-single-adjusting-compartment (SAC) model",
    "for the ALK inhibitor brigatinib in healthy adults (Hanley 2024).",
    "The source model was built in the Simcyp Population-based Simulator",
    "(versions 15 and 17) and its whole-body mass-balance equations are",
    "not published, so the platform model itself cannot be encoded here.",
    "What IS fully reported is the brigatinib compound layer, and it is",
    "sufficient to reconstruct the disposition as an ordinary",
    "compartmental model: first-order absorption into a depot,",
    "distribution between a systemic compartment and the SAC (the",
    "paper's k_in / k_out, encoded as the canonical k12 / k21), and",
    "first-order elimination from the systemic compartment. Systemic",
    "clearance is obtained from the reported unbound hepatic intrinsic",
    "clearance through the well-stirred liver model plus the reported",
    "renal clearance; bioavailability is the reported f_a times f_G",
    "times the f_H implied by that same well-stirred calculation. No",
    "parameter is fitted here and none is imported from a Simcyp",
    "population file: every value is either a Table 1 / Appendix S1",
    "input or an arithmetic consequence of one. The reduction",
    "reproduces the paper's own predicted Cmax, Tmax and AUC after",
    "single 90 mg and 180 mg oral doses to within 6.4%, and Cmax at",
    "both doses to within 1% (see the validation vignette).",
    "This is a typical-value simulation model:",
    "the source reports no inter-individual variance components and no",
    "residual-error model, so there are no etas and propSd is fixed at",
    "zero. The drug-drug-interaction predictions that are the paper's",
    "main contribution (itraconazole, rifampin, diltiazem, verapamil,",
    "efavirenz, and the transporter substrates) depend on proprietary",
    "Simcyp perpetrator compound files and are NOT reproducible from",
    "this model; see the vignette for the full list of deviations.",
    sep = " "
  )
  reference <- paste(
    "Hanley MJ, Rowland Yeo K, Tugnait M, Iwasaki S, Narasimhan N,",
    "Zhang P, Venkatakrishnan K, Gupta N. (2024).",
    "Evaluation of the drug-drug interaction potential of brigatinib",
    "using a physiologically-based pharmacokinetic modeling approach.",
    "CPT Pharmacometrics Syst Pharmacol 13(4):624-637.",
    "doi:10.1002/psp4.13106.",
    sep = " "
  )
  vignette <- "Hanley_2024_brigatinib"
  units    <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Hanley 2024 Figure 1b (the minimal
  # PBPK + SAC schematic) and the Figure 1 legend.
  compartmentData <- list(
    depot = list(
      analyte = "brigatinib", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "brigatinib", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "brigatinib", units = "mg",
      specimen = "tissue", verified = TRUE
    )
  )

  # No covariates are carried in this reduction. The Simcyp simulations
  # behind Hanley 2024 vary body weight, age, sex, albumin, haematocrit
  # and serum creatinine across the virtual population, but those act
  # through Simcyp population files that are not reported, so no covariate
  # relationship in the published paper is reproducible here. Recorded as
  # screened-but-not-carried rather than silently dropped.
  covariatesDataExcluded <- list(
    WT = list(
      description = paste(
        "Body weight. The source expresses Vss and Vsac in L/kg",
        "(Hanley 2024 Table 1), so the Simcyp model does scale",
        "distribution volume with weight. This reduction fixes the",
        "70 kg reference value instead of carrying a weight term,",
        "because the systemic-compartment volume also subtracts a",
        "liver volume whose weight-scaling Simcyp does not publish;",
        "a linear weight term on vc would therefore be a guess about",
        "the platform rather than a reported relationship."
      ),
      units = "kg",
      type  = "continuous",
      notes = "Implicit in the L/kg volume inputs; not carried in this reduction."
    ),
    ALB = list(
      description = paste(
        "Serum albumin. Drives the fraction unbound in the Simcyp",
        "Sim-Cancer population via Appendix S1 Equation S6, which",
        "calibrates a healthy-participant fu to a disease-population",
        "fu using the ratio of albumin concentrations. Appendix S1",
        "Table S5 reports albumin of 50.3 / 49.38 g/L (male / female",
        "healthy) versus 38.1 / 39.25 g/L in patients with cancer,",
        "but the resulting cancer fu is never printed, so the",
        "relationship cannot be reproduced quantitatively."
      ),
      units = "g/L",
      type  = "continuous",
      notes = "Reported qualitatively only; the calibrated cancer fu is not printed."
    ),
    CRCL = list(
      description = paste(
        "Creatinine clearance. Appendix S1 Equation S7 generates a",
        "creatinine clearance per virtual subject and divides it by",
        "120 (female) or 130 (male) to produce a renal-function",
        "scalar applied directly to the renal clearance input.",
        "Equation S7 itself is one of the equations the PDF renders",
        "as an undecodable formula image, so the scalar cannot be",
        "reconstructed."
      ),
      units = "mL/min",
      type  = "continuous",
      notes = "Scales CLR in the Simcyp Sim-Cancer population; equation not decodable."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 36L,
    n_studies      = 6L,
    age_range      = "23-62 years across the simulated healthy-participant scenarios",
    weight_median  = "70 kg (reference weight used for the L/kg volume inputs)",
    sex_female_pct = 30,
    disease_state  = paste(
      "Healthy adult participants. The paper also simulates patients",
      "with ALK-positive advanced malignancies using the Simcyp",
      "Sim-Cancer population, but that arm depends on an unpublished",
      "albumin-calibrated fraction unbound and is not reproduced by",
      "this model."
    ),
    dose_range     = "Single oral doses of 90 mg and 180 mg; multiple-dose 90 mg q.d. for 7 days then 180 mg q.d.",
    regions        = "North European Caucasian virtual population (Simcyp default) for the healthy-participant simulations.",
    studies        = paste(
      "AP26113-15-106 bioequivalence study (n = 36 healthy",
      "participants) is the study the disposition parameters ka, Vss,",
      "Vsac, kin and kout were optimized against and the study whose",
      "observed CL/F of 13.47 L/h anchors the retrograde clearance",
      "calculation. AP26113-13-104 human ADME study (n = 6 healthy",
      "males, single 180 mg / 100 uCi [14C]-brigatinib) supplies the",
      "mass-balance split that fixes f_a and the renal fraction.",
      "AP26113-15-105 drug-drug-interaction study (three arms of",
      "n = 20) supplies the itraconazole and rifampin observations.",
      "AP26113-15-107 (hepatic impairment) and AP26113-15-108 (renal",
      "impairment) contributed the ex vivo plasma protein binding",
      "(n = 17 pooled) that gives fu = 0.0883. AP26113-11-101",
      "(NCT01449461) is the phase I/II cancer study used for the",
      "multiple-dose verification arm."
    ),
    notes          = paste(
      "n_subjects records the 36 participants of AP26113-15-106, the",
      "single study the disposition parameters were optimized against;",
      "n_studies counts the six clinical studies listed above that",
      "contributed inputs. This is a PBPK analysis rather than a",
      "population-PK fit, so there is no single pooled analysis dataset",
      "and no estimated variance components."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Every parameter below is fixed: nothing was estimated in building
    # this reduction. Values are either verbatim Hanley 2024 Table 1 /
    # Appendix S1 inputs, or arithmetic consequences of them. The two
    # derivations (clearance and systemic volume) are spelled out in
    # full so they can be checked line by line.
    #
    # Shared quantities, all Table 1 / Appendix S1:
    #   fu  = 0.0883      fraction unbound in plasma
    #   B:P = 0.69        blood-to-plasma ratio
    #   fub = fu / (B:P) = 0.0883 / 0.69 = 0.127971
    #   QH  = 90 L/h      hepatic blood flow (Appendix S1 Equation S2)
    #   CLuH,int = 12.21 uL/min/mg microsomal protein (Table 1)
    #
    # Rescaling CLuH,int to a whole-liver flow with the two scalars
    # Appendix S1 states for exactly this purpose (liver weight 1648 g;
    # 39.8 mg microsomal protein per g liver):
    #   12.21 * 39.8 * 1648 uL/min = 800859 uL/min
    #                              = 0.800859 L/min * 60 = 48.05 L/h
    #
    # NOTE ON A SUPPLEMENT TYPO. Appendix S1 states this whole-liver
    # value as "43.02 L/h". That is inconsistent with its own two
    # scalars, which give 48.05 L/h, and it is 43.02 that is wrong:
    # substituting 48.05 L/h back into the retrograde Equation S2
    # returns the observed CL/F of 13.47 L/h to three significant
    # figures, whereas 43.02 L/h returns 12.71 L/h (5.6% low). The
    # 12.21 uL/min/mg figure in Table 1 -- the value actually entered
    # into Simcyp, and the value apportioned across pathways in
    # Appendix S1 Table S3 -- corresponds to 48.05 L/h, not 43.02 L/h.
    # This file therefore uses 48.05 L/h and treats 43.02 as a
    # transcription error in the supplement.
    # ------------------------------------------------------------------

    lka <- fixed(log(1.8))
    label("First-order absorption rate constant ka (1/h)")
    # Hanley 2024 Table 1: ka = 1.8 1/h, optimized against study
    # AP26113-15-106. Appendix S1 records that the permeability-predicted
    # value of 0.13 1/h gave a Tmax about 3-fold too high and that 1.8
    # was chosen by sensitivity analysis to recover the observed ~2 h Tmax.

    # Systemic (plasma) clearance, well-stirred liver plus renal:
    #   fub * CLuH,int          = 0.127971 * 48.05  = 6.1492 L/h
    #   CL_H,blood = QH * 6.1492 / (QH + 6.1492)
    #              = 90 * 6.1492 / 96.1492           = 5.7559 L/h
    #   CL_H,plasma = CL_H,blood * (B:P)
    #              = 5.7559 * 0.69                   = 3.9716 L/h
    #   CL = CL_H,plasma + CLR = 3.9716 + 3.21       = 7.1816 L/h
    # CLR = 3.21 L/h is Hanley 2024 Table 1, derived in Appendix S1 as
    # 21.38 / 89.75% x 13.47 L/h from the AP26113-13-104 mass balance.
    lcl <- fixed(log(7.1816))
    label("Systemic plasma clearance CL (L/h)")
    # Derived from Hanley 2024 Table 1 (CLuH,int 12.21 uL/min/mg;
    # CLR 3.21 L/h; fu 0.0883; B:P 0.69) and Appendix S1 Equation S2
    # (QH 90 L/h; liver weight 1648 g; 39.8 mg protein/g liver).

    # Systemic-compartment volume. Table 1 gives the whole-body
    # steady-state volume Vss = 6.91 L/kg and the SAC volume
    # Vsac = 6.0 L/kg. In the Simcyp minimal-PBPK layout of Figure 1b
    # the remaining volume is split between the systemic compartment and
    # the liver, so at the 70 kg reference weight:
    #   (6.91 - 6.0) L/kg * 70 kg = 63.70 L
    #   minus the liver           = 63.70 - 1.648 = 62.05 L
    # The 1648 g liver weight is the same Appendix S1 default used in
    # the clearance rescaling above. The portal-vein volume of Figure 1b
    # is not reported and is not subtracted; see the vignette Errata.
    lvc <- fixed(log(62.05))
    label("Systemic compartment volume vc (L) at the 70 kg reference weight")
    # Derived from Hanley 2024 Table 1 (Vss 6.91 L/kg; Vsac 6.0 L/kg)
    # and the Appendix S1 default liver weight of 1648 g.

    # The paper's k_in and k_out are first-order rate constants acting on
    # the MASSES of drug in the systemic compartment and in the SAC
    # respectively (Hanley 2024 Figure 1 legend). That is exactly the
    # definition of the canonical central-to-peripheral1 and
    # peripheral1-to-central micro-constants, so they are encoded as
    # k12 / k21. Note that the SAC volume Vsac = 6.0 L/kg does not enter
    # the plasma prediction at all under a mass-based parameterisation
    # and is therefore not carried as a parameter; a peripheral volume
    # implied by vc * k12 / k21 would NOT equal Vsac (see vignette).
    lk12 <- fixed(log(28.13))
    label("Systemic-to-SAC transfer rate constant k12 (paper k_in, 1/h)")
    # Hanley 2024 Table 1: kin = 28.13 1/h, optimized against AP26113-15-106.

    lk21 <- fixed(log(20.00))
    label("SAC-to-systemic transfer rate constant k21 (paper k_out, 1/h)")
    # Hanley 2024 Table 1: kout = 20.00 1/h, optimized against AP26113-15-106.

    # Oral bioavailability is the product of the three sequential
    # fractions the source reports or implies:
    #   f_a = 0.63   fraction absorbed (Table 1; Appendix S1 derives
    #                63.26% of dose absorbed from the mass balance)
    #   f_G = 0.9    fraction escaping gut metabolism (Appendix S1
    #                Equation S2 text, "predicted to be approximately
    #                0.9 based on preliminary simulations")
    #   f_H = 0.936  fraction escaping hepatic first pass, which is
    #                QH / (QH + fub * CLuH,int) = 90 / 96.1492 from the
    #                same well-stirred calculation used for lcl above
    #   F = 0.63 * 0.9 * 0.936046 = 0.5307
    lfdepot <- fixed(log(0.5307))
    label("Oral bioavailability F (fraction)")
    # Derived from Hanley 2024 Table 1 (fa 0.63) and Appendix S1
    # Equation S2 (fG 0.9; fH from QH, fu, B:P and CLuH,int).

    # Hanley 2024 is a PBPK simulation analysis, not a population-PK
    # fit. It reports no residual-error model and no inter-individual
    # variance components; the percent-CV figures in Tables 2-4 are the
    # spread of a Simcyp virtual population driven by unpublished
    # population files, not estimated omegas. Rather than invent a
    # variance, the residual error is fixed at zero, which makes this a
    # deterministic typical-value simulation model.
    propSd <- fixed(0)
    label("Proportional residual error SD (fraction; zero, no error model reported by the source)")
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual parameters. No covariates and no random effects.
    # ------------------------------------------------------------------
    ka     <- exp(lka)
    cl     <- exp(lcl)
    vc     <- exp(lvc)
    k12    <- exp(lk12)
    k21    <- exp(lk21)
    fdepot <- exp(lfdepot)

    # ------------------------------------------------------------------
    # 2. Micro-constant for elimination.
    # ------------------------------------------------------------------
    kel <- cl / vc

    # ------------------------------------------------------------------
    # 3. ODE system, corresponding to Hanley 2024 Figure 1b once the
    # liver and portal-vein compartments are lumped into the systemic
    # compartment. `peripheral1` is the paper's single adjusting
    # compartment (SAC): k12 acts on the mass in the systemic
    # compartment and moves drug into the SAC, k21 acts on the mass in
    # the SAC and returns it, per the Figure 1 legend definitions of
    # k_in and k_out. Amounts are in mg.
    # ------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central -
                          k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # ------------------------------------------------------------------
    # 4. Bioavailability on the oral depot.
    # ------------------------------------------------------------------
    f(depot) <- fdepot

    # ------------------------------------------------------------------
    # 5. Observation. Doses are in mg and vc is in L, so central / vc is
    # in mg/L = ug/mL; multiply by 1000 to report ng/mL, the units used
    # throughout Hanley 2024 Tables 2-4.
    # ------------------------------------------------------------------
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
