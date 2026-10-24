Pan_2019_fentanyl_lumped_b <- function() {
  description <- paste(
    "PBPK (proper-lumped, 4 states). Lumped fentanyl PBPK model 'b' of Pan",
    "and Duffull 2019 (Figure 4b), obtained by proper lumping of the 16-state",
    "Bjorkman 1994 whole-body system (Pan_2019_fentanyl_pbpk) with arterial",
    "blood constrained to stay unlumped. States: `lumped1` (Figure 4b 'L1' =",
    "venous blood, carcass, gut, gut_deep), `arterial`, `lumped2` (Figure 4b",
    "'L2' = heart, kidney, spleen, skin, adipose, liver, liver_deep) and",
    "`lumped3` (Figure 4b 'L3' = brain, muscle, gut_deep2, spleen_deep). A",
    "lumped state is the SUM of the concentrations (ug/L) of the original",
    "states it pools. Rate constants are the 4 x 4 lumped",
    "micro-rate-constant matrix K_hat = M K M+ (Supplementary Section 3); no between-subject variability and",
    "no residual error (deterministic system)."
  )
  reference <- paste(
    "Pan S, Duffull SB (2019). Automated proper lumping for simplification of",
    "linear physiologically based pharmacokinetic systems. J Pharmacokinet",
    "Pharmacodyn 46(4):361-370. doi:10.1007/s10928-019-09644-5.",
    "Lumped matrix from Supplementary Material 1 Section 3; lumping",
    "assignments from Figure 4b; venous volume and dose scaling from",
    "Section 1A and the deposited Runfile_Autolumping_FenPBPK.m."
  )
  vignette <- "Pan_2019_fentanyl_lumping"
  units <- list(time = "min", dosing = "ug", concentration = "ug/L")

  compartmentData <- list(
    lumped1 = list(analyte = "fentanyl", units = "ug/L", specimen = "whole blood", verified = TRUE),
    arterial = list(analyte = "fentanyl", units = "ug/L", specimen = "whole blood", verified = TRUE),
    lumped2 = list(analyte = "fentanyl", units = "ug/L", specimen = "tissue", verified = TRUE),
    lumped3 = list(analyte = "fentanyl", units = "ug/L", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 1L,
    age_range = "adult (reference man)",
    weight_range = "73 kg (reference man)",
    disease_state = "Not applicable: a reduction of a physiology-based reference individual, not a fitted cohort.",
    dose_range = "750 ug intravenous bolus",
    notes = paste(
      "Same 73 kg reference man as Pan_2019_fentanyl_pbpk (Supplementary",
      "Material 1 Section 1A). One of five simulated-annealing solutions",
      "meeting the acceptance criterion of an absolute relative difference",
      "of at most 0.002% in arterial AUC over 0-210 min against the full",
      "system. `lumped1` pools venous blood with tissues; its specimen is",
      "recorded as whole blood because it carries the dose and the venous",
      "return."
    )
  )

  ini({
    # Pan 2019 Eq 2, K_hat = M K M+, with K the Section 1B / Model_parameter_values.m
    # 16 x 16 matrix (see Pan_2019_fentanyl_pbpk) and M the Figure 4b lumping
    # matrix: {venous, carcass, gut, gut_deep} | arterial | {heart, kidney,
    # spleen, skin, adipose, liver, liver_deep} | {brain, muscle, gut_deep2,
    # spleen_deep}. Supplementary Material 1 Section 3 prints the result to two
    # decimals, 'K = [-0.61 0.87 0.03 0.03; 0.84 -3.35 0 0; 0.002 0.60 -0.06
    # 0.01; 0.02 15.9 6.74 -5.53]', with rows and columns ordered L1, arterial,
    # L3, L2 -- the two tissue lumps the other way round from the Figure 4b
    # labels. The arterial-inflow column fixes that order: 0.60 = brain 0.56 +
    # muscle 0.038 (Figure 4b L3) and 15.9 = heart + kidney + pancreas/spleen +
    # skin + fat (Figure 4b L2). Here lumped2 = Figure 4b L2 (printed row and
    # column 4) and lumped3 = Figure 4b L3 (printed row and column 3). Every
    # value below rounds to the print; the unrounded values are carried
    # because the two-decimal print moves the terminal phase off the paper's
    # own Figure 6b. k_<from>_<to> is the off-diagonal entry K_hat[to, from];
    # kloss_<state> is the negated diagonal entry.
    kloss_lumped1 <- fixed(0.608415); label("Loss rate constant of lumped1 (1/min)") # S3 K[1,1] printed -0.61
    k_arterial_lumped1 <- fixed(0.87); label("Arterial-to-lumped1 rate constant (1/min)") # S3 K[1,2] printed 0.87; = 0.04 + 0.83 (S1A per-kg flows of carcass, gut)
    k_lumped3_lumped1 <- fixed(0.0346151); label("Lumped3-to-lumped1 rate constant (1/min)") # S3 K[1,3] printed 0.03
    k_lumped2_lumped1 <- fixed(0.0266616); label("Lumped2-to-lumped1 rate constant (1/min)") # S3 K[1,4] printed 0.03
    k_lumped1_arterial <- fixed(0.8375); label("Lumped1-to-arterial rate constant (1/min)") # S3 K[2,1] printed 0.84; = (Qcad/Vart) / 4
    kloss_arterial <- fixed(3.35); label("Loss rate constant of arterial (1/min)") # S3 K[2,2] printed -3.35
    k_lumped1_lumped3 <- fixed(0.00230061); label("Lumped1-to-lumped3 rate constant (1/min)") # S3 K[3,1] printed 0.002
    k_arterial_lumped3 <- fixed(0.598); label("Arterial-to-lumped3 rate constant (1/min)") # S3 K[3,2] printed 0.60; = 0.56 + 0.038 (S1A per-kg flows of brain, muscle)
    kloss_lumped3 <- fixed(0.0647907); label("Loss rate constant of lumped3 (1/min)") # S3 K[3,3] printed -0.06
    k_lumped2_lumped3 <- fixed(0.0141576); label("Lumped2-to-lumped3 rate constant (1/min)") # S3 K[3,4] printed 0.01
    k_lumped1_lumped2 <- fixed(0.0153704); label("Lumped1-to-lumped2 rate constant (1/min)") # S3 K[4,1] printed 0.02
    k_arterial_lumped2 <- fixed(15.9236); label("Arterial-to-lumped2 rate constant (1/min)") # S3 K[4,2] printed 15.9
    k_lumped3_lumped2 <- fixed(6.7378); label("Lumped3-to-lumped2 rate constant (1/min)") # S3 K[4,3] printed 6.74
    kloss_lumped2 <- fixed(5.53493); label("Loss rate constant of lumped2 (1/min)") # S3 K[4,4] printed -5.53

    v_venous <- fixed(3.9); label("Venous blood volume, scales the bolus into a concentration (L)") # S1A 'Vven = 3.9'; Runfile 'IC = [Dose/Vven; ...]'
  })

  model({
    d/dt(lumped1) <- -kloss_lumped1 * lumped1 + k_arterial_lumped1 * arterial +
      k_lumped2_lumped1 * lumped2 + k_lumped3_lumped1 * lumped3
    d/dt(arterial) <- k_lumped1_arterial * lumped1 - kloss_arterial * arterial
    d/dt(lumped2) <- k_lumped1_lumped2 * lumped1 + k_arterial_lumped2 * arterial -
      kloss_lumped2 * lumped2 + k_lumped3_lumped2 * lumped3
    d/dt(lumped3) <- k_lumped1_lumped3 * lumped1 + k_arterial_lumped3 * arterial +
      k_lumped2_lumped3 * lumped2 - kloss_lumped3 * lumped3

    # Lumped initial condition IC_hat = M IC: the bolus (ug) enters the
    # venous-containing lump as concentration Dose / Vven.
    f(lumped1) <- 1 / v_venous

    # Observation: arterial blood concentration (unlumped output state).
    Cc <- arterial
  })
}
