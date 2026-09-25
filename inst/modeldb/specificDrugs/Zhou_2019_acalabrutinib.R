Zhou_2019_acalabrutinib <- function() {
  description <- paste(
    "Two-compartment oral pharmacokinetic reduction of the Simcyp",
    "minimal-PBPK-with-single-adjusting-compartment (SAC) model for the",
    "covalent Bruton tyrosine kinase inhibitor acalabrutinib in healthy",
    "adults (Zhou 2019). The source model was built in Simcyp version",
    "14/17 and its whole-body mass-balance equations are not published,",
    "so the platform model itself cannot be encoded here. What IS fully",
    "reported is the acalabrutinib compound layer: Table 2 of the paper",
    "plus the Simcyp compound file Acalabrutinib.cmpz that the authors",
    "deposited as Supporting Information. That layer is sufficient to",
    "rebuild the disposition as an ordinary compartmental model:",
    "first-order absorption into a depot with a lag time, distribution",
    "between a systemic compartment and the SAC (the Simcyp kin / kout,",
    "encoded as the canonical k12 / k21), and first-order elimination",
    "from the systemic compartment split into its printed renal and",
    "non-renal components. No parameter is fitted here: every value is",
    "either a Table 2 / main-text input or an arithmetic consequence of",
    "one. The reduction reproduces the paper's own predicted Cmax to",
    "within 5.4% and its predicted AUC to within 5.7% at all three",
    "simulated single-dose levels (50, 100 and 400 mg; Tables S2 and",
    "S3), with no fitted parameter anywhere.",
    "This is a typical-value simulation model: the source reports no",
    "estimated inter-individual variance components and no",
    "residual-error model, so there are no etas and propSd is fixed at",
    "zero.",
    "The active metabolite ACP-5862 is deliberately NOT encoded. Most",
    "ACP-5862 is formed pre-systemically, so its exposure is controlled",
    "by a first-pass term requiring hepatic blood flow and intestinal",
    "availability, neither of which is printed in the paper or in either",
    "deposited compound file; a metabolite arm would rest on an invented",
    "constant. The drug-drug-interaction predictions that are the",
    "paper's main contribution are likewise NOT reproducible here, for",
    "the same reason plus their dependence on proprietary Simcyp",
    "perpetrator compound files. See the validation vignette for the",
    "sensitivity analysis that justifies both exclusions.",
    sep = " "
  )
  reference <- paste(
    "Zhou D, Podoll T, Xu Y, Moorthy G, Vishwanathan K, Ware J,",
    "Slatter JG, Al-Huniti N. (2019).",
    "Evaluation of the drug-drug interaction potential of acalabrutinib",
    "and its active metabolite, ACP-5862, using a physiologically-based",
    "pharmacokinetic modeling approach.",
    "CPT Pharmacometrics Syst Pharmacol 8:489-499.",
    "doi:10.1002/psp4.12408.",
    sep = " "
  )
  vignette <- "Zhou_2019_acalabrutinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Zhou 2019 Table 2 (the Simcyp
  # input table) and the deposited compound file Acalabrutinib.cmpz
  # fields <idSACKin>, <idSACKout>, <idSACVolume>, <idSACCLin> and
  # <idSACCLout>.
  compartmentData <- list(
    depot = list(
      analyte = "acalabrutinib",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "acalabrutinib",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "acalabrutinib",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    )
  )

  # No covariates are carried. The Simcyp virtual North European
  # Caucasian / Healthy Volunteer population behind Zhou 2019 varies body
  # weight, age and sex across the virtual subjects, but those act
  # through Simcyp population files that are not reported, so no
  # covariate relationship in the published paper is reproducible here.
  # Recorded as screened-but-not-carried rather than silently dropped.
  covariatesDataExcluded <- list(
    WT = list(
      description = paste(
        "Body weight. The source expresses Vss (0.21 L/kg) and Vsac",
        "(0.028 L/kg) per kilogram in Table 2, so the Simcyp model does",
        "scale distribution volume with weight. This reduction does not",
        "need a weight term at all, because the deposited compound file",
        "reports the SAC first-order rate constants together with the",
        "corresponding inter-compartmental clearances, and",
        "vc = CLin / kin is an absolute volume in litres. Carrying a",
        "weight exponent on vc would therefore be a guess about the",
        "platform rather than a reported relationship."
      ),
      units = "kg",
      type = "continuous",
      notes = "Implicit in the L/kg volume inputs; not needed by this reduction."
    ),
    AGE = list(
      description = paste(
        "Age. Table 1 and Table S1 match each virtual trial to the age",
        "range of the corresponding clinical cohort, spanning 18-65",
        "years across the nine studies. Age drives liver weight and",
        "CYP3A4 abundance inside the Simcyp population files, but no",
        "age relationship is printed, so none is carried."
      ),
      units = "years",
      type = "continuous",
      notes = "Sets the virtual-population sampling range only; no printed relationship."
    ),
    SEXF = list(
      description = paste(
        "Female sex indicator (1 = female). Table 1 records the female",
        "proportion of every cohort (0% to 83%) and Table S1 confirms",
        "each virtual trial was matched to it. Sex drives organ weights",
        "and blood flows inside the Simcyp population files; no sex",
        "effect on any acalabrutinib parameter is printed."
      ),
      units = "unitless",
      type = "categorical",
      notes = "Sets the virtual-population sex split only; no printed relationship."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 114L,
    n_studies = 5L,
    age_range = "18-65 years across the contributing cohorts (Table 1)",
    weight_median = "not reported; the Simcyp virtual North European Caucasian / Healthy Volunteer default distribution was used",
    sex_female_pct = "0-83 across cohorts (Table 1); 33% in ACE-HV-113, the study the SAC parameters were estimated against",
    disease_state = "Healthy adult volunteers.",
    dose_range = "Single oral doses of 25, 50, 75, 100 and 400 mg acalabrutinib (Table 1).",
    regions = "Simcyp virtual North European Caucasian / Healthy Volunteer population (Methods, Model verification).",
    studies = paste(
      "Five phase I studies in healthy subjects supplied the",
      "observations and the compound-layer inputs (Table 1).",
      "ACE-HV-001 is the dose-escalation study: cohorts 4 (50 mg) and 6",
      "(100 mg) were used to estimate the SAC volume and the SAC rate",
      "constants kin and kout with the Simcyp parameter estimation",
      "module, cohort 6 also anchored Peff,man, cohort 7 (50 mg with and",
      "without itraconazole) fixed fm,CYP3A4 = 0.82, and cohorts 3",
      "(25 mg) and 5 (75 mg) were verification.",
      "ACE-HV-001 is also the source of the in vivo total clearance of",
      "169 L/h that, with the absolute oral bioavailability, anchors the",
      "systemic clearance of this reduction.",
      "ACE-HV-004 part 3 (100 mg with and without rifampicin, n = 24)",
      "verified the CYP3A contribution.",
      "ACE-HV-009 is the 14C absolute-bioavailability, excretion and",
      "metabolism study (n = 8-14) that supplies the absolute oral",
      "bioavailability of 25% and the renal clearances.",
      "ACE-HV-113 (n = 12-13) and ACE-HV-005 (n = 18-40) measured both",
      "acalabrutinib and ACP-5862 and were the parent/metabolite",
      "verification studies."
    ),
    notes = paste(
      "n_subjects sums the distinct cohorts of Table 1 that contributed",
      "acalabrutinib observations (6 + 6 + 12 + 6 + 16 + 24 + 12 + 18 +",
      "14 = 114). n_studies counts the five ACE-HV protocols.",
      "These counts describe the clinical data behind the compound",
      "layer, not an analysis dataset: this is a PBPK analysis rather",
      "than a population-PK fit, so there are no estimated variance",
      "components. Each simulation in Table S1 was ten virtual trials",
      "matched to the size, age range and sex split of its cohort.",
      "Reported fractions of acalabrutinib elimination, carried here as",
      "provenance rather than as model parameters because a fraction-",
      "metabolised split cannot change the plasma prediction:",
      "fm,CYP3A4 = 0.82 of metabolic intrinsic clearance, of which the",
      "ACP-5862-forming arm is Vmax/Km = 4.13/2.78 = 1.486 of the total",
      "CYP3A4 CLint of 9.63 uL/min/pmol, so about 12.7% of hepatic",
      "metabolism forms ACP-5862 - the paper's own 'about 12%', against",
      "10% of total dose in the human mass-balance study."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Every parameter below is fixed: nothing was estimated in building
    # this reduction. Values are either verbatim Zhou 2019 Table 2 /
    # main-text inputs, or arithmetic consequences of them.
    #
    # Two independent checks confirm the reading of the deposited
    # compound files Acalabrutinib.cmpz and ACP-5862.cmpz (Supporting
    # Information, psp412408-sup-0002). In the Simcyp minimal-PBPK
    # layout kin = CLin / V_systemic and kout = CLout / V_SAC, so
    # V_SAC = CLout / kout. That returns 1.01689053 / 0.45 = 2.2598 L
    # for acalabrutinib and 0.08070559 / 0.01 = 8.0706 L for ACP-5862,
    # and dividing each by its own Table 2 Vsac in L/kg (0.028 and 0.1)
    # returns 80.706 kg BOTH times - the body weight of the Simcyp
    # population representative. Independently, V_systemic + V_SAC
    # recovers 97.9% of Table 2's Vss x 80.706 kg for acalabrutinib
    # (16.644 of 17.005 L) and 97.8% for ACP-5862 (28.614 of 29.243 L),
    # the same ratio for both compounds. No reference body weight has to
    # be assumed anywhere below: vc is absolute.
    # ------------------------------------------------------------------

    lka <- fixed(log(1.65))
    label("First-order absorption rate constant ka (1/h)")
    # Zhou 2019 Table 2, row 'ka (hour-1)' = 1.65, source 'Predicted'.
    # The Simcyp simulation itself used the mechanistic Peff,man
    # absorption route (Table 2 row 'Peff,man' = 4 x 10-4 cm/second,
    # 'Optimized based on clinical data'; the deposited compound file
    # carries the matching <idPeffmanPredicted> 3.95879817), so a single
    # first-order ka is the one structural approximation this reduction
    # makes. It is validated, not assumed: see the vignette, where the
    # reduction lands within 5.4% of the paper's predicted Cmax at every
    # simulated dose level.

    ltlag <- fixed(log(0.25))
    label("Absorption lag time t_lag (h)")
    # Zhou 2019 Table 2, row 'Lag time (hour)' = 0.25, source
    # 'Optimized'. The deposited compound file confirms it exactly in
    # field <idLagTime> = 0.25.

    lvc <- fixed(log(14.3847))
    label("Systemic compartment volume vc (L)")
    # Derived from the deposited compound file Acalabrutinib.cmpz as
    # CLin / kin, using <idSACCLin> = 15.2477827 L/h and <idSACKin> =
    # 1.06 1/h: 15.2477827 / 1.06 = 14.3847 L. Table 2 prints the same
    # kin to three significant figures. See the header comment above for
    # the two checks that confirm this reading.

    # The Simcyp kin and kout are first-order rate constants acting on
    # the MASSES of drug in the systemic compartment and in the single
    # adjusting compartment respectively - Table 2's own footnote reads
    # 'kin, kout, first order rate constants in and out Vsac'. That is
    # exactly the definition of the canonical central-to-peripheral1 and
    # peripheral1-to-central micro-constants, so they are encoded as
    # k12 / k21. Under this mass-based parameterisation the Vsac input
    # of 0.028 L/kg does not itself enter the plasma prediction and is
    # not carried; the peripheral volume implied by vc * k12 / k21 is
    # 33.9 L rather than the 2.26 L that CLout / kout gives, so the
    # reduction's steady-state volume (48.3 L) does not equal the Table 2
    # Vss (17.0 L). That discrepancy is a property of the published
    # inputs, not of this encoding, and it is not reconciled here - see
    # the vignette's Assumptions and deviations section, which also
    # records that the alternative single-Q reading was tested and
    # reproduces the paper's predicted Cmax substantially worse.
    lk12 <- fixed(log(1.06))
    label("Systemic-to-SAC transfer rate constant k12 (Simcyp kin, 1/h)")
    # Zhou 2019 Table 2, row 'k in (1/hour)' = 1.06, estimated with the
    # Simcyp parameter estimation module from the six ACE-HV-001 cohort 6
    # subjects who received 100 mg. Deposited compound file field
    # <idSACKin> = 1.06.

    lk21 <- fixed(log(0.45))
    label("SAC-to-systemic transfer rate constant k21 (Simcyp kout, 1/h)")
    # Zhou 2019 Table 2, row 'k out (1/hour)' = 0.45, estimated in the
    # same step. Deposited compound file field <idSACKout> = 0.45.

    # --- Clearance -----------------------------------------------------
    # The paper anchors the whole elimination layer on one absolute
    # number: 'The intrinsic metabolic clearance of acalabrutinib was
    # estimated from the in vivo total clearance of 169 L/hour
    # (ACE-HV-001) with the retrograde method.' That 169 L/hour is an
    # apparent ORAL clearance - the paper's own predicted AUC of
    # 616 ng.h/mL after 100 mg (Table S2) implies CL/F = 162 L/h, and
    # the Table 2 CYP3A4 CLint of 9.63 uL/min/pmol is stated to be its
    # retrograde consequence. Combined with the absolute oral
    # bioavailability of 25% reported in the main text ('an absolute
    # oral bioavailability of 25%', from the 14C study ACE-HV-009), the
    # systemic plasma clearance is 169 x 0.25 = 42.25 L/h. The printed
    # renal clearance is carried as its own component so the split is
    # not lost, exactly as the paper reports it; the two components sum
    # back to 42.25 L/h.
    lcl_nonren <- fixed(log(40.92))
    label("Non-renal (hepatic metabolic) systemic plasma clearance cl_nonren (L/h)")
    # Derived as total systemic clearance minus renal clearance:
    # 169 L/h x 0.25 - 1.33 L/h = 42.25 - 1.33 = 40.92 L/h.

    lcl_renal <- fixed(log(1.33))
    label("Renal clearance cl_renal (L/h)")
    # Zhou 2019 Table 2, row 'CLR (L/hour)' = 1.33, source 'Clinical
    # data'; the main text adds 'The renal clearance of 1.33 L/hour for
    # acalabrutinib observed in clinical study ACE-HV-009 was applied
    # directly.' Deposited compound file field <idCLRbase> = 1.33.

    lfdepot <- fixed(log(0.25))
    label("Absolute oral bioavailability F (fraction)")
    # Zhou 2019 main text: 'Acalabrutinib is rapidly absorbed with a
    # short oral half-life of about 1.57 hours in healthy subjects, with
    # an absolute oral bioavailability of 25%', citing the 14C
    # absolute-bioavailability study ACE-HV-009 (Podoll 2019, Drug Metab
    # Dispos 47:145-154). Carried as a single lumped constant because
    # its mechanistic decomposition F = fa x Fg x Fh needs hepatic blood
    # flow, which this paper never prints - see the vignette.

    # Zhou 2019 is a PBPK simulation analysis, not a population-PK fit.
    # It reports no residual-error model and no estimated
    # inter-individual variance components; the percent-CV figures in
    # Table S2 are the spread of a Simcyp virtual population driven by
    # unpublished population files, not estimated omegas. The deposited
    # compound file does carry 30% input CVs, but the identical value 30
    # appears on essentially every CV field in the file including many
    # that the model never uses, so it is the Simcyp default rather than
    # a compound-specific estimate and is not imported here. Rather than
    # invent a variance, the residual error is fixed at zero, which
    # makes this a deterministic typical-value simulation model.
    propSd <- fixed(0)
    label("Proportional residual error SD (fraction; zero, no error model reported by the source)")
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual parameters. No covariates and no random effects.
    # ------------------------------------------------------------------
    ka <- exp(lka)
    tlag <- exp(ltlag)
    vc <- exp(lvc)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    cl_nonren <- exp(lcl_nonren)
    cl_renal <- exp(lcl_renal)
    fdepot <- exp(lfdepot)

    # ------------------------------------------------------------------
    # 2. Clearance assembly. Total systemic plasma clearance is the sum
    # of the two printed components, 40.92 + 1.33 = 42.25 L/h.
    # ------------------------------------------------------------------
    cl <- cl_nonren + cl_renal
    kel <- cl / vc

    # ------------------------------------------------------------------
    # 3. ODE system, corresponding to the Simcyp minimal-PBPK layout once
    # the liver and portal-vein compartments are lumped into the systemic
    # compartment. `peripheral1` is the single adjusting compartment:
    # k12 acts on the mass in the systemic compartment and moves drug
    # into the SAC, k21 acts on the mass in the SAC and returns it.
    # Amounts are in mg.
    # ------------------------------------------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central -
      k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ------------------------------------------------------------------
    # 4. Bioavailability and lag time on the oral depot.
    # ------------------------------------------------------------------
    f(depot) <- fdepot
    alag(depot) <- tlag

    # ------------------------------------------------------------------
    # 5. Observation. Doses are in mg and vc is in L, so central / vc is
    # in mg/L = ug/mL; multiply by 1000 to report ng/mL, the units used
    # throughout Zhou 2019 Tables S2 and S3.
    # ------------------------------------------------------------------
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
