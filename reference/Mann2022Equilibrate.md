# Pre-equilibrate the Mann 2022 respiratory-physiology states

The Mann 2022 / FDA delaymymod.c physiology layer is sensitive to
initial conditions because the chemoreflex feedback sees a 48.6 mm Hg
error in \`p_b_co2 - p_b_co2_0\` at t = 0 when the FDA delaystates.R
nominal values (\`palv_co2 = 40.28\`, \`cb_co2 = 0.645\`) are used
verbatim. The transient over-ventilation that follows leaves the system
in a state that subsequent opioid doses cannot push into sustained PaO2
\< 15 mm Hg, severely under-estimating cardiac arrest incidence.

## Usage

``` r
Mann2022Equilibrate(model, params = NULL, duration_min = 90)
```

## Arguments

- model:

  An \`rxode2\` model containing at least the Mann 2022 respiratory
  physiology states (\`palv_co2\`, \`palv_o2\`, \`cb_co2\`, \`cb_o2\`,
  \`ct_co2\`, \`ct_o2\`, \`yco2\`, \`yo2\`, \`dp_state\`, \`dc_state\`,
  \`alpha_h\`). Works with the standalone
  \`Mann_2022_respiratory_physiology\` model, with composed chains (PK +
  binding + physiology), or with any user-extended variant.

- params:

  Named vector of model parameters. For a composed chain, supply the
  full PK + binding + physiology parameter set (the opioid PK states
  stay at zero throughout the pre-equilibration because no dose is
  given).

- duration_min:

  Pre-equilibration duration in minutes. The default \`90\` is
  sufficient for the Mann 2022 physiology to settle to four-decimal
  precision; the slowest state is \`alpha_h\` with a 5-minute time
  constant.

## Value

Named list of physiology equilibrium state values, suitable for passing
as the \`inits\` argument to \`rxode2::rxSolve\` on subsequent
dose-bearing simulations.

## Details

The FDA \`simulateToGetOD_IM.R\` reference script handles this by
running a no-drug pre-simulation via \`fundede\` per subject and using
the FINAL state as the dose-time initial state. This function does the
equivalent: it runs the supplied \`model\` for \`duration_min\` minutes
with zero opioid dose and zero antagonist input, then returns the
final-time values of the standard Mann 2022 physiology states as a named
list suitable for passing as \`inits\` to \`rxode2::rxSolve\`.

The equilibrium does NOT depend on opioid PK parameters (no drug is
present during pre-equilibration). The standalone
\`Mann_2022_respiratory_physiology\` model and the Laffont 2025 vignette
inline chain ship with the FDA nominal initial conditions as \`ini()\`
defaults (matching \`delaystates.R\`), so any code that uses these
models for downstream overdose simulation MUST call this function and
pass the result via the \`inits\` argument of \`rxode2::rxSolve\`;
without it the FDA-published cardiac-arrest rates cannot be reproduced.
