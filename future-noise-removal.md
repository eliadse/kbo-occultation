# Future work: noise removal without hurting detection

Deferred from `noise-implementation-plan.md`. Depends on its Part 3.5 (FAR curves) and Part 1 (signal band).

- **Figure of merit:** efficiency ε(R, mag) at a fixed real-data FAR, built from the Part 3.5 curve.
- **Methods to compare:**
  - Butterworth cutoff scan (0.3–3 Hz) with the template always reshaped;
  - PSD-whitened matched filter;
  - DC detrend window scan (2, 5, 10 s);
  - fixed DC FFT swap: crossover ≤ 0.8 × DC Nyquist, raised-cosine transition, tapered record ends;
  - DC despike used as a veto;
  - no cross-telescope common-mode subtraction, because a real event is itself common-mode within ≤ 8 ms.
- **Fixes:**
  - rewrite `check_highpass_signal_safety.py` with the embedded-template method;
  - unify the scattered cutoffs (1 / 3 / 3.5 Hz) into one config constant;
  - fix `DC_CUT_HZ`.
- **Model consistency:**
  - add the variance-estimator floor √(2k/N) and a red-noise term to `noise_budget`;
  - use `scint_coeff = 1.5`;
  - inject the template on the star part only (exclude NSB and pedestal).

