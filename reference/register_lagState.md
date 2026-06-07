# Register the lagState C user function with rxode2

lagState() implements a true transport delay on an algebraic or state
expression evaluated during ODE integration. It is required by the
Mann_2022_respiratory_physiology model to encode the FDA delaymymod.c
peripheral (~7 s) and central (~11 s) chemoreflex transport delays
exactly, without resorting to first-order-lag approximations (which damp
amplitude as well as delaying, distorting the overdose trajectory).

## Usage

``` r
register_lagState()
```

## Details

The function signature is: \`lagState(t, val, lag, channel, init_val)\`

Arguments inside an rxode2 \`model()\` block: \* \`t\` - current
integration time (use the model symbol \`t\`) \* \`val\` - the algebraic
/ state expression to be delayed \* \`lag\` - delay duration, in the
same time units as the model \* \`channel\` - integer channel index in
\`\[0, 8)\` so multiple delays may coexist in one model without buffer
collision \* \`init_val\` - value returned when \`t - lag \< 0\` (the
pre-history window) - typically the steady-state of \`val\`

Behaviour: \* Per-thread, per-channel ring buffer of recent (time,
value) pairs. \* Detects new-subject start (t returning to ~0) and
resets the buffer. \* Handles solver rollback (entries with t greater
than the current integration time are discarded). \* Linear
interpolation between bracketing buffer entries. \* If the buffer
overflows (the requested lagged time has fallen off the oldest
still-kept entry), \`Rf_error\` is raised with a clear diagnostic
message - never silently return a wrong value.

The C code uses GCC thread-local storage (\`\_\_thread\`), so each
rxode2 worker thread maintains independent buffers. This is safe because
rxode2 parallelises across subjects (one subject per thread at a time).

This function is called automatically at package load via \`.onLoad\`,
so end users do not need to invoke it. Registration is idempotent.
