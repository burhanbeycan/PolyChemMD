# Density sampling review

`analysis/sampling.py` evaluates an explicitly selected, uniformly sampled
window from one uninterrupted state and ensemble. It executes no simulation.
Inputs are rows containing exactly `step`, `temperature_k`, `pressure_bar` and
`density_g_cm3`; pressure and density must already use canonical units. Duplicate
restart-boundary frames, uneven sampling, unknown fields and nonfinite values
are rejected. A continuation of one cell remains one cell.

```python
from polychemmd.analysis.sampling import (
    DensitySamplingSpec, analyze_density_sampling,
)
spec = DensitySamplingSpec(timestep_fs=0.25, discard_steps=200000,
                           min_duration_ps=100.0, block_size=50)
review = analyze_density_sampling(rows, spec)
```

This example declares analysis settings; it supplies neither a trajectory nor a
validated polymer protocol. Record the specification before examining results.
The package does not select a discard or tune thresholds to make a series pass.
`examples/density_sampling.json` uses a 10 ps illustrative review threshold to
match its campaign example; the module default remains 100 ps. Neither duration
is a universal physical threshold or evidence that the example has been run.
Post-hoc adjustment invalidates a predeclared qualification comparison.

## Estimators

For retained density values, finite autocovariances use the common divisor
`N`, including only the available `N-k` products at lag `k`:

\[
\widehat\gamma_k=\frac1N\sum_{i=1}^{N-k}(x_i-\bar x)(x_{i+k}-\bar x),
\qquad \widehat\rho_k=\widehat\gamma_k/\widehat\gamma_0.
\]

Adjacent pairs start at lags `(0,1)`, `(2,3)`, etc. They are truncated before the
first nonpositive pair and constrained to be nonincreasing. With the resulting
positive pair sums `P_j`, the estimate is

\[
\widehat g=\max(1,-1+2\sum_jP_j),\qquad
\widehat N_{\rm eff}=N/\widehat g.
\]

This follows the initial monotone sequence construction of Geyer, using
normalized finite autocovariances and an additional conservative floor at one.
Geyer's theoretical setting is a reversible stationary Markov chain; those
assumptions are **not established** for a supplied MD series. Reaching the
requested maximum lag before observing a nonpositive pair fails the closure
gate. Observing a negative pair does not establish that slower modes are absent.
See [Geyer (1992), DOI 10.1214/ss/1177011137](https://doi.org/10.1214/ss/1177011137)
and the [author's method documentation](https://www.stat.umn.edu/geyer/mcmc/library/mcmc/html/initseq.html).

The link between integrated autocorrelation, statistical inefficiency and
effective samples in molecular simulation is described by
[Chodera (2016), DOI 10.1021/acs.jctc.5b00784](https://doi.org/10.1021/acs.jctc.5b00784),
[equations 9–13 and discussion](https://pmc.ncbi.nlm.nih.gov/articles/PMC4945107/).
Here the finite autocovariance already incorporates loss of products with lag;
another `(1-k/N)` weight is not applied. This implementation uses a predeclared
discard rather than Chodera's automatic discard optimization.

Nonoverlapping complete blocks give `SE_block = s(block means)/sqrt(B)`.
Pairing consecutive blocks gives a second SE at twice the block length.
A third estimate is `s(density)*sqrt(g/N)`. The reported conditional normal
half-width is `1.959963984540054` times the largest of those three estimates.
It is an approximate nominal 95% quantity, **not a calibrated confidence
interval or coverage guarantee**. Finite correlation estimates and biased
preparation can both make uncertainty estimates unreliable.

## Predeclared review gates

Defaults are transparent development review settings, not published universal
physical thresholds or a polymer equilibration method.

| Check | Default |
| --- | --- |
| Physical span of analyzed complete-block frames | At least 100 ps |
| Complete blocks and coarsened blocks | At least 16 and 8 |
| Estimated effective samples | At least 50 |
| Block length relative to estimated `g` | At least `5*g` samples |
| Absolute lag-one block autocorrelation | At most 0.20 |
| Coarsened/original SE ratio | Between `1/1.25` and `1.25` |
| Conditional density half-width / mean | At most 0.02 |
| Relative split-half and four-window density changes | At most 0.02 each |
| Temperature sample-mean deviation from target | At most 15 K |
| Pressure sample-mean deviation from target | At most 100 bar |

The temperature and pressure gates concern sample means, never one instantaneous
pressure excursion. They remain heuristic state checks. Tolerances must be
justified for a declared chemistry, cell size, ensemble and planned duration;
these defaults are not a substitute for that study. The analyzed physical span
excludes any incomplete trailing block. Both the full retained span and the
analyzed span are reported. Constant densities or identical block means are
inconclusive and never receive a zero-uncertainty certificate.

`diagnostics_passed` means only that these checks pass. Every return keeps
`equilibrium_status="not_assessed"`, `production_ready=false` and
`density_prediction_validation="not_performed"`. Phase, chain relaxation,
preparation sensitivity, cell-size/cutoff convergence, force-field accuracy and
state-matched external comparison remain separate requirements. The existing
50 fs preparation pilot fails the physical-duration and block-count gates.

## Independent cells

`summarize_density_replicates(records, expected_cell_ids=[...])` requires every
predeclared cell, including failures. Each record contains exactly
`independent_cell_id`, `state_id`, `protocol_id`, `status` and
`sampling_analysis`. Execution status is `completed`, `failed` or `inconclusive`;
failed/inconclusive execution has `sampling_analysis=None`. Completed execution
can still have insufficient sampling. State, protocol and sampling settings must
match. Duplicate identities, missing cells and inconsistent passed summaries
are rejected.

Only when at least three declared cells all pass the diagnostics and their
means differ is an equal-cell mean and between-cell standard error returned.
Frame counts never increase the replicate count; no density interval is
reported from failed or incomplete replicate sets. Cell identity and independent
preparation remain declared provenance rather than cryptographically verified
facts. Sampling-review status cannot authenticate those claims.

The tests use hand-computed four-point correlations, an arithmetic block-mean
example and explicitly synthetic independent/AR(1) series. They verify
estimators and failure handling; these series are not simulated polymer data.
