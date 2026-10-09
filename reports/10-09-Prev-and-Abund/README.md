## 10/09 (Prev and Abund)

## 1. Review
- 上周我们讨论算法在真实数据上的结果, 具体而言是
  - ZIHPN算法和maaslin3的算法在真实数据上看不出区别
  - abundance部分对于某些taxon出现一定的inflation, 具体表现是qqplot上翘, 结果是只有一个taxon被检测出来有效应
  - prevalence部分出现较为严重的inflation, 导致多个假阳性, 猜测是未做SPA的原因, 简单采用SAIGE测试, 发现新的qqplot贴合对角线, 但大规模测试未做
- 同时讨论了一个奇怪的现象, seq depth取很高, 无sample zero的时候, 效应被加上abundance下, 且depth在case control之间不平衡, prevalence部分校正depth会出现问题, abundance部分则不会

## 2. Prevalence Part
