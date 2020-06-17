# 2020-06-17 14:10:07

Quick detour on some of the models I'm running:

```
So we can think of ~0+x or equally ~x+0 as an equation of the form: y=ax+b. By adding 0 we are forcing b to be zero, that means that we are looking for a line passing the origin (no intercept). If we indicated a model like ~x+1 or just ~x, there fitted equation could possibily contain a non-zero term b. Equally we may restrict b by a formula ~x-1 or ~-1+x that both mean: no intercept (the same way we exclude a row or column in R by negative index). However something like ~x-2 or ~x+3 is meaningless.
```

So, in R, removing the intercept terms makes it a more ANOVA-like model. In the
intercept model below, the vs coefficient becomes difference. Without the
intercept, there are two vs coefficients modeling the conditional expectations for both groups.

```
> dat <- mtcars
> dat$vs <- factor(dat$vs)
> lm(mpg ~ vs + hp + vs:hp, data = dat)

Call:
lm(formula = mpg ~ vs + hp + vs:hp, data = dat)

Coefficients:
(Intercept)          vs1           hp       vs1:hp  
   24.49637     14.50418     -0.04153     -0.11657  

> lm(mpg ~ 0 + vs + hp + vs:hp, data = dat)

Call:
lm(formula = mpg ~ 0 + vs + hp + vs:hp, data = dat)

Coefficients:
     vs0       vs1        hp    vs1:hp  
24.49637  39.00055  -0.04153  -0.11657  
```

The RSE for both models is identical, as is the significance of all terms but
the vs term.

But this also makes a nice compelling case for keeping the intercept:

https://blog.minitab.com/blog/adventures-in-statistics-2/regression-analysis-how-to-interpret-the-constant-y-intercept

But most limma and edgeR tutorials set the intercept to zero, I bet because they
always test a contrast. For example:

https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4937821/pdf/f1000research-5-18347.pdf

The latter even says:

```
For a given experiment, there are usually several equivalent ways to set up an appropriate design matrix. For example, ~0+group+lane removes the intercept from the first factor, group, but an intercept remains in the second factor lane. Alternatively, ~group+lane could be used to keep the intercepts in both group and lane. Understanding how to interpret the coefficients estimated in a given model is key here. We choose the first model for our analysis, as setting up model contrasts is more straight forward in the absence of an inter- cept for group.
```


Also helpful in understanding contrasts:

https://genomicsclass.github.io/book/pages/interactions_and_contrasts.html

Especially this part:

```
The question of whether the push vs. pull difference is different in L2 compared to L1, is answered by a single term in the model: the typepush:legL2 estimated coefficient corresponding to the yellow arrow in the plot. A p-value for whether this coefficient is actually equal to zero can be read off from the table printed with summary(fitX) above.
```

In our case, looking at the Region:DX term in the model only tests whether the
Caudate vs ACC difference is different in ADHD compared to normals. That's not
what we want (at least not just that). We want to know where ADHD differs from
NV regardless of region. Then we could partition it by region. 

