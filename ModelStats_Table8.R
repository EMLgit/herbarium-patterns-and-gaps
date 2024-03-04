---
title: "TableS8_statsAnalyses"
author: "emlombardi"
date: "2024-01-18"
output: html_document
---

## Running statistical analyses on supplemental table 8

Honestly we just need some better quantitative measures to talk about. I'm going to use GLMM to analyze how different fixed and random variables impact the number of collections made (across environmental variables)

Internet resources that are very useful:
1. Basically all resources available from Dr. Ben Bolker
https://bbolker.github.io/goettingen_2019/notes/modeling_inference.html
https://bbolker.github.io/goettingen_2019/
https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf
https://bbolker.github.io/mixedmodels-misc/ecostats_chap.html

Maybe: https://cran.r-project.org/web/packages/glmm/vignettes/vignettes.pdf

```{r libraries}
library(dplyr)
library(glmm)
```

## Read data
#This is a de-formated version of Supplemental Table S8 in the manuscript ['citation']


tab8 <- read.csv("~/Desktop/Research/UNM/Herbarium collections patterns and gaps project/Tables/TableS8.csv")
summary(tab8)


#Check the variance etc. 
It's bascially count data, even if we decide to also model Shannon diversity response. Thus, Poisson distribution is most accurate. 
First we need to explore the data generally though
```{r}
ggplot(tab8,aes(no.records,no.species, color=category))+
  geom_point() +
  theme_bw()

ggplot(data=tab8, aes(x=environment, y=no.records, group=category, color=category))+
  geom_point()+
  facet_grid(category ~ .)

ggplot(tab8,aes(x=no.records,y=shannon,colour=layer))+
    stat_sum(aes(size=after_stat(n)),alpha=0.7)+
    scale_y_log10()+
    geom_smooth(method="glm",method.args=list(family=quasipoisson)) +
    theme_bw()



#VERY BASIC
order <- c("Biotic", "Abiotic", "Sociopolitical")

# Convert the 'category' variable to a factor with the desired order
tab8$category <- factor(tab8$category, levels = order)


catcolors <- c( "#009E73", "#56B4E9", "#CC79A7")
s5<-ggplot(tab8,aes(no.records, log(no.species), color=category))+
  geom_point() +
  geom_smooth(method="loess") +
  scale_color_manual(values=catcolors) +
  theme_bw() +
  facet_wrap(~category) +
  labs(x = "Number of records", y = "Log-transformed Number of species")

ggsave("/Users/elizabethlombardi/Desktop/Research/UNM/Herbarium collections patterns and gaps project/Figures/supplementSppAccumCurve.png", s5, width = 10, height = 6, dpi = 300)


ggplot(tab8,aes(no.records, no.species))+
  geom_point() +
  geom_smooth(method="loess", color="darkslategrey") +
  theme_bw() +
  labs(x = "Number of records", y = "Number of species")


library(broom)

###BASIC model
mod1 <- glm(no.species ~ no.records + environment, data = tab8)
tidymod1 <- as.data.frame(tidy(mod1))
summary(mod1)

pvals.mod1 <- tidymod1 %>%
  select(term, p.value)

mod1.sigEnv <- pvals.mod1 %>% 
  filter(p.value < 0.05)


#try by category; I don't think the results here are sensible...? Intercept is a bit strange. 
mod2 <- glm(no.species ~ no.records + category, data = tab8)
tidymod2 <- as.data.frame(tidy(mod2))
summary(mod2)


#look at relative diversity given clusters
plot(tab8.rev$rel.diversity ~ tab8.rev$clustering)
mod3 <- glm(rel.diversity ~ clustering + environment, data = tab8.rev)
tidymod3 <- as.data.frame(tidy(mod3))
summary(mod3)

#remove the clustering outlier, which is the rock house-nestor mountain homestead acec
tab8.rev2 <- tab8.rev %>% filter(!environment == "rock house-nestor martin homestead acec")
plot(tab8.rev2$rel.diversity ~ tab8.rev2$clustering)
mod4 <- glm(rel.diversity ~ clustering + environment, data = tab8.rev2)
tidymod4 <- as.data.frame(tidy(mod4))
summary(mod4)


#I think we should probably keep the rock house-nestor martin homestead acec out of the previous model, too...
#same as mod1 but without nestor martin acec
mod1b <- glm(no.species ~ no.records + environment, data = tab8.rev2)
tidymod1b <- as.data.frame(tidy(mod1b))
summary(mod1b)

pvals.mod1b <- tidymod1b %>%
  select(term, p.value)

mod1.sigEnv <- pvals.mod1b %>% 
  filter(p.value < 0.05)


```
#Jans suggestions
-saturation
-rate of saturation

Instead of using Shannon diversity, look at Anne Chao and/or Jost (Hill numbers)
(r package iNext)
-completion based indices to check out gaps (but probs still don't want to actually talk about gaps)


start with online resource on iNEXT and hill numbers based on Chao's work:
https://bookdown.org/c_w_beirne/wildCo-Data-Analysis/composition.html

#Basic model
```{r}
m1 <- lmer(log(shannon) ~ log(no.records) + (1|category), data=tab8)
summary(m1)
hist(residuals(m1))
anova(m1)

#diagnostic plots
qqmath(ranef(m1), distribution = qnorm)

```

#Next chunk
```{r}
m2 <- lmer(shannon ~ no.records+(1|category), data=tab8)
summary(m2)
plot(m2)

#should I rescale no.records? 
tab8.rev <- transform(tab8, C.no.records=no.records - mean(no.records))
m3 <- lmer(shannon ~ no.records+(1|category), data=tab8.rev)
#didn't fix it, so I'll take a different approach

#i. z-score scaling
tab8.rev$z.shannon <- scale(tab8.rev$shannon)
tab8.rev$z.recs <- scale(tab8.rev$no.records)

m4 <- lmer(z.shannon ~ z.recs * (1|category), data=tab8.rev)
summary(m4)

#ii. log-log rescaling
tab8.rev$log.shannon <- log(tab8.rev$shannon)
tab8.rev$log.recs <- log(tab8.rev$no.records)

m5 <- lmer(log.shannon ~ log.recs + (1 | category), data = tab8.rev)
summary(m5)


###CHECK
hist(residuals(m1))
plot(predict(m1), residuals(m1))
summary(m1)
anova(m1)


```


###Try glmm 

```{r}
m6 <- glmer(no.species ~ log(no.records) + (1|environment), data=tab8.rev, family="poisson")


#estimate theta
fit <- glm.nb(no.species ~ no.records, data=tab8.rev)
theta <- fit$theta

m7 <- glmer(log(no.species) ~ log(no.records) + (1|environment), data=tab8.rev, family=negative.binomial(theta))
hist(residuals(m7))
plot(predict(m7))

```