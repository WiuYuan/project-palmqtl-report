## 09/08 (Model in AA)

## 1. Review

- We built a pipeline that simulates based on spike-in data.
- We found that the Maaslin3 model actually performs very well, even though it directly models TSS.
- We tried using information from the positive part when modeling the zero component, but found that when the case/control ratio is relatively small, the FDR cannot be controlled.

## 2. Explore

- Maaslin3 is able to avoid the outlier problem for a simple reason: on the log scale, outliers no longer exert such a large influence on the score.
- Under the meta-test in MaAsLin3, the heterogeneity (het) remains well controlled.

![](assets/robust.png)
![](assets/zero.png)
![](assets/positive.png)
![](assets/joint.png)

- GWAS Simulation
![](assets/gwas_zero_example.png)

## 3. Future
- feed positive data only, maaslin2/3 + palm + palmqtl
