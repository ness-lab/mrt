## Estimating effective population size from diversity and population-wide recombination

### Catherine Burzynski BIO481 project
#### Supervised by: James S. Santangelo, Ahmed Hasan, and Joel Levine

### Project

One summary statistic that is commonly used to describe populations is the effective population size (Ne), which roughly speaking describes the amount of genetic drift in populations. The lower the Ne, the more drift, and hence the greater the random fluctuations in allele frequencies from one generation to the next. This is important since it directly impacts the efficacy of natural selection: more drift = lower efficacy of selection, which can be especially problematic for small populations trying to recover from disturbances (think conservation planning). Therefore, accurately estimating Ne is important for understanding the relative roles of drift vs. selection in influencing the evolutionary potential and trajectory of natural populations.

Ne can be estimated from genomic sequence data in two ways: the first relies on patterns of genetic diversity while the second relies on patterns of recombination. Using the first, if we know the mutation rate and the level of genetic diversity, we can estimate Ne. Using the second, if we know the population recombination rate (think frequency of sex), we can estimate Ne. Under ideal conditions in an unchanging population (i.e., one at equilibrium), the estimates of Ne using both approaches will be the same. However, perturbations (e.g., population bottlenecks) can influence recombination and diversity to varying extents, resulting in disagreements between the estimates of Ne obtained from both approaches. This is the crux of the project, which aims to address the following question: **How do non-equilibrium conditions influence estimates of Ne obtained from genetic diversity and recombination?**

To answer this, we’ll simulate genomic data in SLiM (a common population genetic simulation software) using a number of demographic scenarios (e.g., constant size, reductions in size, etc.). Because we are simulating, mutation and recombination rates will be known, allowing us to estimate Ne on the simulated data using both approaches without issue.

## Resources

- [Youtube playlist by Dr. Mohamed Noor](https://www.youtube.com/watch?v=UeWU1yOz8lQ&list=PL4n6Uk3aii8iottl_J2OWrn-8RVveZWri)
- [Beginner programmign lessons](https://utm-coders.github.io/studyGroup/lessons/)
    - [Youtube videos accompanying above lessons](https://www.youtube.com/channel/UCTTUketB568idCSUsXUug6w)
