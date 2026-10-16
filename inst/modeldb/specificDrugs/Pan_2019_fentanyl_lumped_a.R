Pan_2019_fentanyl_lumped_a <- function() {
  description <- paste(
    "PBPK (proper-lumped, 4 states). Lumped fentanyl PBPK model 'a' of Pan",
    "and Duffull 2019 (Figure 4a), obtained by proper lumping of the 16-state",
    "Bjorkman 1994 whole-body system (Pan_2019_fentanyl_pbpk) with arterial",
    "blood constrained to stay unlumped. States: `venous` (venous blood,",
    "unlumped; Figure 4a 'L1'), `arterial`, `lumped1` (Figure 4a 'L2' = brain,",
    "heart, kidney, gut, gut_deep, gut_deep2, spleen, spleen_deep, liver,",
    "liver_deep, skin) and `lumped2` (Figure 4a 'L3' = muscle, adipose,",
    "carcass). A lumped state is the SUM of the concentrations (ug/L) of the",
    "original states it pools. Rate constants are the 4 x 4 lumped",
    "micro-rate-constant matrix K_hat = M K M+ (Supplementary Section 2); no between-subject variability",
    "and no residual error (deterministic system)."
  )
  reference <- paste(
    "Pan S, Duffull SB (2019). Automated proper lumping for simplification of",
    "linear physiologically based pharmacokinetic systems. J Pharmacokinet",
    "Pharmacodyn 46(4):361-370. doi:10.1007/s10928-019-09644-5.",
    "Lumped matrix from Supplementary Material 1 Section 2; lumping",
    "assignments from Figure 4a; venous volume and dose scaling from",
    "Section 1A and the deposited Runfile_Autolumping_FenPBPK.m."
  )
  vignette <- "Pan_2019_fentanyl_lumping"
  units <- list(time = "min", dosing = "ug", concentration = "ug/L")

  compartmentData <- list(
    venous = list(analyte = "fentanyl", units = "ug/L", specimen = "whole blood", verified = TRUE),
    arterial = list(analyte = "fentanyl", units = "ug/L", specimen = "whole blood", verified = TRUE),
    lumped1 = list(analyte = "fentanyl", units = "ug/L", specimen = "tissue", verified = TRUE),
    lumped2 = list(analyte = "fentanyl", units = "ug/L", specimen = "tissue", verified = TRUE)
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
      "Material 1 Section 1A). Found by simulated annealing (and by full",
      "enumeration) as a 4-state model meeting the acceptance criterion of",
      "an absolute relative difference of at most 0.002% in arterial AUC over",
      "0-210 min against the full system."
    )
  )

  ini({
    # Pan 2019 Eq 2, K_hat = M K M+, with K the Section 1B / Model_parameter_values.m
    # 16 x 16 matrix (see Pan_2019_fentanyl_pbpk) and M the Figure 4a lumping
    # matrix: venous | arterial | {brain, heart, kidney, liver, skin, gut,
    # gut_deep, gut_deep2, spleen, spleen_deep, liver_deep} | {muscle, adipose,
    # carcass}. Supplementary Material 1 Section 2 prints the result to two
    # decimals, 'K = [-1.72 0 0.02 0.08; 3.35 -3.35 0 0; 0 17.3 -1.15 0; 0 0.11
    # 0 -0.02]'; every value below rounds to that print. The unrounded values
    # are carried because the 0.02 print of the slowest constant (0.01726) moves
    # the terminal phase off the paper's own Figure 6a; the vignette shows both.
    # k_<from>_<to> is the off-diagonal entry K_hat[to, from]; kloss_<state> is
    # the negated diagonal entry. Zero entries are structural and omitted.
    kloss_venous <- fixed(1.71795); label("Loss rate constant of venous (1/min)") # S2 K[1,1] printed -1.72; = Qcad/Vven = 6.7/3.9
    k_lumped1_venous <- fixed(0.021266); label("Lumped1-to-venous rate constant (1/min)") # S2 K[1,3] printed 0.02
    k_lumped2_venous <- fixed(0.0818176); label("Lumped2-to-venous rate constant (1/min)") # S2 K[1,4] printed 0.08
    k_venous_arterial <- fixed(3.35); label("Venous-to-arterial rate constant (1/min)") # S2 K[2,1] printed 3.35; = Qcad/Vart
    kloss_arterial <- fixed(3.35); label("Loss rate constant of arterial (1/min)") # S2 K[2,2] printed -3.35
    k_arterial_lumped1 <- fixed(17.2856); label("Arterial-to-lumped1 rate constant (1/min)") # S2 K[3,2] printed 17.3
    kloss_lumped1 <- fixed(1.15462); label("Loss rate constant of lumped1 (1/min)") # S2 K[3,3] printed -1.15
    k_arterial_lumped2 <- fixed(0.106); label("Arterial-to-lumped2 rate constant (1/min)") # S2 K[4,2] printed 0.11; = 0.038 + 0.028 + 0.04 (S1A per-kg flows of muscle, fat, carcass)
    kloss_lumped2 <- fixed(0.0172604); label("Loss rate constant of lumped2 (1/min)") # S2 K[4,4] printed -0.02

    v_venous <- fixed(3.9); label("Venous blood volume, scales the bolus into a concentration (L)") # S1A 'Vven = 3.9'; Runfile 'IC = [Dose/Vven; ...]'
  })

  model({
    d/dt(venous) <- -kloss_venous * venous + k_lumped1_venous * lumped1 + k_lumped2_venous * lumped2
    d/dt(arterial) <- k_venous_arterial * venous - kloss_arterial * arterial
    d/dt(lumped1) <- k_arterial_lumped1 * arterial - kloss_lumped1 * lumped1
    d/dt(lumped2) <- k_arterial_lumped2 * arterial - kloss_lumped2 * lumped2

    # Lumped initial condition IC_hat = M IC: the bolus (ug) enters venous as
    # concentration Dose / Vven.
    f(venous) <- 1 / v_venous

    # Observation: arterial blood concentration (unlumped output state).
    Cc <- arterial
  })
}
