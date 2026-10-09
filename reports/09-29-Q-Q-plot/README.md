## 09/29 (Q-Q plot)

## 1. Review

- Maaslin3’s method is actually slightly less powerful than our method in the positive part when the mean sequencing depth is low.
- After reviewing several methods, I found that most of them handle the zero part and the positive part separately and then apply certain transformations (such as RINT, analyzing residuals after fitting a simple model, or considering principal components via PCA). In particular, the zero part is consistently modeled using logistic regression.
- Specifically, further inspection revealed that three studies explicitly stated in the main text that they performed rarefaction (two of them clearly indicated re-sampling to 10,000 reads and directly discarding samples with fewer than 10,000 reads), while one study included sequencing depth as a covariate.
- This week, my main work focused on examining sampling zeros and continuing to complete the code.

## 2. Q-Q plot

The only remaining part of the GWAS code is the addition of parallelization; therefore, the numerical tests are still being run locally and have not yet been applied to real data.

We first evaluate the algorithm under the null hypothesis. We estimate the true abundance data (1000 samples, 300 taxa) based on spike-in data, then generate a count matrix through sampling. Next, we fix a case/control ratio and simulate 10,000 SNPs before drawing QQ plots. We compare different methods, namely:

- ZIHPN: zero-inflated hurdle Poisson–Normal prior model
- LN-Wald: log-normal + logistic, using the Wald test (no results for the zero part, as it could not be simulated)
- LN-Score: log-normal + logistic, using the score test
- PALM

We consider tests for both beta and beta0 separately, all with the simplest median correction applied. The results are as follows.

![](assets/1.png)
## 3. Future
- code for real data on server
