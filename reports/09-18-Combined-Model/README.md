## 09/18 (Combined Model)

## 1. Review

- We found that the estimate for $b$ may be inadequate, so we tried raising the Poisson truncation point. However, after several days of experimentation, it turns out there was actually an error in the original normal generation process; therefore, increasing the truncation point does not in fact change the FDR control.

## 2. Combined Model

### 2.1. Model Definition

$$\pi_i = P(Y_i>0) = \mathrm{logistic}(x_i^\top \alpha_0)=\frac{1}{1+e^{-x_i^\top\alpha_0}}$$
$$Y_i\mid Y_i>0, b_i \sim \mathrm{\text{Zero-Truncated-Poisson}}(\mu_i)$$
$$\begin{aligned}
\ell(\alpha_0,\alpha,b,\tau) &=
\sum_{i\in S_0}\log ((1-\pi_i)+\pi_i P_{0i})\\
&\quad + \sum_{i\in S_+}\Big[Y_i(x_i^\top\alpha+\log N_i+b_i)-e^{x_i^\top\alpha+\log N_i+b_i} + \log\pi_i\Big]\\
&\quad - \frac{1}{2\tau}\sum_{i\in S_+}b_i^2 - \frac{|S_+|}{2}\log(2\pi\tau).
\end{aligned}$$
### 2.2. Simulation

![](assets/1.png)
![](assets/2.png)
![](assets/3.png)

## 3. Future
- data generation with annotation
- GWAS
- search mbGWAS paper
