# Statistical foundations for the diffusion-driven Hawkes project

This note summarizes the statistical ideas used in the project. The goal is to understand the inference pipeline behind the scripts, not only to run the code.

The model combines three ingredients:

1. a continuously evolving covariate $X(t)$;
2. a point process of events $N(t)$;
3. an intensity $\lambda(t)$ describing the instantaneous event rate.

The basic model used in the first scripts is

$$
\lambda(t)
=
\exp(\mu_0+\mu_1 X(t))
+
\sum_{t_j<t} \alpha e^{-\beta(t-t_j)}.
$$

The first term is the covariate-driven baseline. The second term is the Hawkes self-excitation term.

---

## 1. Point process and event counts

A point process is a random collection of event times. In this project, events could represent neural spikes, failures, messages, or any sequence of time-localized events.

The counting process $N(t)$ counts how many events have occurred up to time $t$.

On a time grid, we write $dN_k$ for the number of events in the small time bin around time $t_k$.

In the scripts, `dN[k]` stores the number of events in bin $k$.

If `dN[k] = 0`, no event occurred in that bin.

If `dN[k] = 1`, one event occurred in that bin.

For small time steps, most bins have 0 events.

---

## 2. Intensity

The intensity $\lambda(t)$ is the instantaneous event rate.

Informally,

$$
\lambda(t) dt
$$

is approximately the probability of observing an event in a small interval of length $dt$.

So if $\lambda(t)$ is high, events are more likely.

If $\lambda(t)$ is low, events are less likely.

In the binned approximation, we use

$$
dN_k \sim \mathrm{Poisson}(\lambda_k dt).
$$

This means that the expected number of events in bin $k$ is $\lambda_k dt$.

---

## 3. Covariate-driven baseline

The baseline intensity is

$$
\exp(\mu_0+\mu_1 X(t)).
$$

The exponential function is used because it guarantees positivity.

The parameter $\mu_0$ controls the general baseline level.

The parameter $\mu_1$ controls how strongly the covariate affects the event rate.

If $\mu_1 = 0$, then $X(t)$ has no effect on the event rate.

If $\mu_1 > 0$, larger values of $X(t)$ increase the event rate.

If $\mu_1 < 0$, larger values of $X(t)$ decrease the event rate.

In our simulations, we often used $\mu_1 = 0.7$, so the covariate effect was positive.

---

## 4. Hawkes self-excitation

The Hawkes term is

$$
\sum_{t_j<t} \alpha e^{-\beta(t-t_j)}.
$$

Each past event increases the current intensity.

The parameter $\alpha$ controls the size of the excitation jump after an event.

The parameter $\beta$ controls how quickly the excitation decays.

Large $\alpha$ means events strongly excite future events.

Large $\beta$ means the excitation decays quickly.

The ratio $\alpha / \beta$ is important for stability. In the one-dimensional exponential Hawkes model, one typically wants

$$
\frac{\alpha}{\beta} < 1.
$$

In the scripts, we often used

$$
\alpha = 0.7,
\qquad
\beta = 2.5,
$$

so

$$
\frac{\alpha}{\beta}=0.28.
$$

---

## 5. Likelihood

A likelihood measures how plausible the observed data are under a given parameter value.

In the binned approximation, we use

$$
dN_k \sim \mathrm{Poisson}(\lambda_k dt).
$$

For each bin, the model assigns a probability to the observed count $dN_k$.

The likelihood is the product of these probabilities over all time bins.

A high likelihood means the parameter values make the observed data plausible.

A low likelihood means the parameter values make the observed data implausible.

---

## 6. Log-likelihood

Products of many probabilities are numerically inconvenient.

Therefore we use the log-likelihood.

Taking logs turns products into sums.

For the binned point-process model, after dropping constants independent of the parameters, the log-likelihood is

$$
\ell
=
\sum_k
\left[
dN_k \log(\lambda_k)
-
\lambda_k dt
\right].
$$

The first term,

$$
dN_k \log(\lambda_k),
$$

rewards the model for assigning high intensity when events occur.

The second term,

$$
-\lambda_k dt,
$$

penalizes the model for assigning high intensity everywhere.

This balance is essential.

A model cannot simply make $\lambda$ huge all the time, because it would pay a large penalty in bins without events.

---

## 7. Maximum likelihood estimation

Maximum likelihood estimation chooses the parameters that maximize the log-likelihood.

In the scripts, the unknown parameters are

$$
\mu_0,\quad \mu_1,\quad \alpha,\quad \beta.
$$

The estimator is the parameter vector that makes the observed data most plausible according to the model.

Since R's `optim()` function minimizes, the scripts minimize the negative log-likelihood.

Thus maximizing

$$
\ell
$$

is equivalent to minimizing

$$
-\ell.
$$

---

## 8. Why optimize log_alpha and log_beta?

The parameters $\alpha$ and $\beta$ must be positive.

To guarantee this, the scripts optimize over

$$
\log(\alpha),
\qquad
\log(\beta).
$$

Then the code recovers

```r
alpha <- exp(log_alpha)
beta <- exp(log_beta)
```

This ensures that the estimated $\alpha$ and $\beta$ are always positive.

---

## 9. BFGS

BFGS is the numerical optimization method used by `optim()`.

It is a deterministic quasi-Newton method.

It starts from an initial guess, evaluates the objective function, estimates a useful descent direction, updates an approximation of curvature, and iterates.

In our scripts, BFGS is used to find parameter estimates.

BFGS does not compute p-values directly.

It only provides fitted parameters and maximized log-likelihoods.

---

## 10. Likelihood-ratio statistic

A likelihood-ratio test compares two nested models:

1. a null model;
2. a larger full model.

The likelihood-ratio statistic is

$$
\mathrm{LR}
=
2
\left(
\ell_{\mathrm{full}}
-
\ell_{\mathrm{null}}
\right).
$$

If the full model fits much better than the null model, then LR is large.

If the full model does not improve the fit much, then LR is small.

The word "ratio" comes from the original likelihood ratio

$$
\frac{L_{\mathrm{full}}}{L_{\mathrm{null}}}.
$$

Since we work with log-likelihoods, the ratio becomes a difference:

$$
\log
\left(
\frac{L_{\mathrm{full}}}{L_{\mathrm{null}}}
\right)
=
\ell_{\mathrm{full}}
-
\ell_{\mathrm{null}}.
$$

---

## 11. Covariate likelihood-ratio test

For the covariate effect, the null hypothesis is

$$
H_0:\mu_1=0.
$$

The full model estimates $\mu_1$ freely.

The null model fixes $\mu_1=0$.

The test asks:

> Does allowing the covariate effect improve the fit enough to reject the hypothesis that the covariate has no effect?

In the simulations where $\mu_1=0.7$, the covariate test strongly rejected $H_0$.

In the null calibration experiment where $\mu_1=0$, the false rejection rate was close to 5%.

This suggests that the standard likelihood-ratio p-value is reasonably well calibrated for the covariate effect.

---

## 12. Hawkes likelihood-ratio test

For Hawkes self-excitation, the null hypothesis is

$$
H_0:\alpha=0.
$$

The full model estimates $\alpha>0$.

The null model fixes $\alpha=0$.

This test asks:

> Does allowing Hawkes self-excitation improve the fit enough to reject the hypothesis that past events do not excite future events?

This test is more delicate than the covariate test.

The reason is that $\alpha=0$ is a boundary point. Since $\alpha$ is constrained to be nonnegative, the null is at the edge of the parameter space.

Also, when $\alpha=0$, the Hawkes term disappears, so $\beta$ is not identifiable.

Indeed,

$$
\alpha e^{-\beta(t-t_j)} = 0
$$

for every value of $\beta$ when $\alpha=0$.

Therefore $\beta$ becomes meaningless under the null hypothesis.

This explains why the naive chi-square p-value is only approximate for the Hawkes test.

---

## 13. P-value

A p-value answers the question:

> If the null hypothesis were true, how often would we see a test statistic at least as extreme as the one observed?

For the likelihood-ratio test, this becomes:

> If the null model were true, how often would LR be at least as large as the observed LR?

A small p-value means the observed statistic would be surprising under the null hypothesis.

A large p-value means the observed statistic is not surprising under the null hypothesis.

A p-value is not the probability that the null hypothesis is true.

---

## 14. Critical value

A critical value gives a threshold for rejection.

At the 5% level, the critical value is chosen so that only 5% of null test statistics exceed it.

For a chi-square approximation with 1 degree of freedom, the 95% critical value is

```r
qchisq(0.95, df = 1)
```

which is approximately

```text
3.841459
```

So the chi-square test rejects at 5% if

$$
\mathrm{LR} > 3.841459.
$$

The p-value and the critical-value approach are equivalent ways of performing the same test.

---

## 15. Parametric bootstrap

A parametric bootstrap computes a p-value by simulation rather than by relying on a theoretical chi-square approximation.

For the Hawkes test, the procedure is:

1. fit the null model $\alpha=0$ to the observed data;
2. simulate many fake datasets under this fitted null model;
3. for each fake dataset, fit the null model and the full model;
4. compute a bootstrap LR statistic for each fake dataset;
5. compare the observed LR to the bootstrap LR distribution.

The bootstrap p-value is

$$
\frac{
1+\#\{\mathrm{LR}_{\mathrm{boot}}\geq \mathrm{LR}_{\mathrm{obs}}\}
}{
1+B
},
$$

where $B$ is the number of bootstrap repetitions.

The `+1` correction avoids returning exactly zero with a finite number of bootstrap samples.

The parametric bootstrap is useful here because it learns the null distribution of LR from simulated null data.

This is especially helpful when the chi-square approximation is questionable.

---

## 16. Bootstrap critical value

The bootstrap critical value is the empirical 95% quantile of the bootstrap LR values.

In R, it is computed by

```r
crit_boot_95 <- quantile(boot_LR, probs = 0.95)
```

This means:

> 95% of the bootstrap LR values are below this value.

So the bootstrap test rejects at 5% if

$$
\mathrm{LR}_{\mathrm{obs}} > \mathrm{crit}_{\mathrm{boot},95}.
$$

In the no-Hawkes bootstrap experiment, the bootstrap critical value was larger than the chi-square critical value.

This showed that the chi-square test was too permissive for the Hawkes boundary case.

---

## 17. What script 01 showed

`01_basic_simulation_and_fit.R` simulated one dataset and fitted the full model.

It showed that the model could recover the parameters reasonably well from one trajectory.

It also computed two likelihood-ratio tests:

1. covariate relevance test: $H_0:\mu_1=0$;
2. Hawkes self-excitation test: $H_0:\alpha=0$.

In the example run, both effects were present and both tests strongly rejected the null hypotheses.

---

## 18. What script 02 showed

`02_monte_carlo_consistency.R` repeated the simulation and fitting many times for several observation lengths $T$.

The goal was to check whether the estimates improve as $T$ increases.

This is an empirical check of consistency.

The estimates became more concentrated around the true values for larger $T$.

The parameter $\beta$ was harder to estimate than $\mu_1$ and $\alpha$, but it also improved with more data.

---

## 19. What script 03 showed

`03_null_calibration.R` studied the tests under true null scenarios.

In the `no_covariate` scenario, the true value was $\mu_1=0$.

The covariate false rejection rate was close to 5% when the number of repetitions was increased to 100.

This suggests that the covariate likelihood-ratio test is reasonably calibrated.

In the `no_hawkes` scenario, the true value was $\alpha=0$.

The Hawkes false rejection rate was larger than expected under the naive chi-square p-value.

This confirmed that the Hawkes test is nonstandard under $\alpha=0$.

---

## 20. What script 04 showed

`04_hawkes_bootstrap_pvalue.R` computed a bootstrap p-value in a no-Hawkes scenario with $\alpha=0$.

The observed LR was zero.

Both the chi-square p-value and the bootstrap p-value were equal to 1.

This correctly indicated no evidence of Hawkes self-excitation.

The bootstrap 95% critical value was larger than the chi-square 95% critical value.

This showed that the bootstrap test is stricter and better adapted to the nonregular null hypothesis.

---

## 21. What script 05 showed

`05_hawkes_bootstrap_positive_control.R` computed a bootstrap p-value in a Hawkes-present scenario with $\alpha=0.7$.

The observed LR was very large.

The bootstrap p-value was equal to the smallest possible value with 100 bootstrap repetitions,

$$
\frac{1}{101}
\approx 0.0099.
$$

This means that none of the 100 bootstrap null datasets produced an LR as large as the observed one.

Therefore the bootstrap test correctly detected Hawkes self-excitation.

Together, scripts 04 and 05 validate the bootstrap test in both directions:

```text
alpha = 0:
    do not reject

alpha = 0.7:
    reject
```

---

## 22. Current statistical conclusions

At this point, the statistical pipeline supports the following conclusions.

First, the maximum-likelihood estimator behaves reasonably in repeated simulations.

Second, the covariate effect can be tested using the standard likelihood-ratio chi-square approximation.

Third, the Hawkes self-excitation test is nonregular at $\alpha=0$.

Fourth, a parametric bootstrap provides a more reliable calibration of the Hawkes test.

Fifth, the parameter $\beta$ should not be interpreted when $\alpha=0$, because it is not identifiable under the null.

---

## 23. Glossary

**Point process:** A random sequence of event times.

**Counting process:** The process $N(t)$ that counts how many events have occurred up to time $t$.

**Intensity:** The instantaneous event rate $\lambda(t)$.

**Covariate:** An external variable $X(t)$ that influences the event rate.

**Hawkes process:** A point process where past events increase the probability of future events.

**Likelihood:** A function measuring how plausible the observed data are under a parameter value.

**Log-likelihood:** The logarithm of the likelihood, easier to optimize.

**Maximum likelihood estimator:** The parameter value that maximizes the likelihood.

**BFGS:** A numerical optimization method used to minimize the negative log-likelihood.

**Likelihood-ratio statistic:** A statistic comparing the fit of a full model to a null model.

**P-value:** The probability, under the null hypothesis, of seeing a statistic at least as extreme as the observed one.

**Critical value:** A threshold above which the null hypothesis is rejected.

**Parametric bootstrap:** A simulation-based method that uses a fitted parametric model to approximate the null distribution of a statistic.

**Boundary case:** A hypothesis where the null value lies at the edge of the allowed parameter space.

**Identifiability:** A parameter is identifiable if different values lead to distinguishable distributions of the data.

**Power:** The probability that a test rejects the null when the alternative is true.

**False rejection rate:** The probability that a test rejects the null when the null is true.