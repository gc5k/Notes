# https://en.wikipedia.org/wiki/Ordinal_regression
# https://stats.idre.ucla.edu/r/dae/ordinal-logistic-regression/
library(foreign)
library(ggplot2)
require(MASS)
require(Hmisc)
require(reshape2)
dat <- read.dta("https://stats.idre.ucla.edu/stat/data/ologit.dta")

lapply(dat[, c("apply", "pared", "public")], table)
ftable(xtabs(~ public + apply + pared, data = dat))
summary(dat$gpa)
sd(dat$gpa)

ggplot(dat, aes(x = apply, y = gpa)) +
  geom_boxplot(size = .75) +
  geom_jitter(alpha = .5) +
  facet_grid(pared ~ public, margins = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

m <- polr(apply ~ pared, data = dat)
summary(m)

m <- polr(apply ~ pared + public + gpa, data = dat, Hess=TRUE)
summary(m)

(ctable <- coef(summary(m)))

p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2

## combined table
(ctable <- cbind(ctable, "p value" = p))
(ci <- confint(m)) # default method gives profiled CIs
confint.default(m) # CIs assuming normality

## https://stats.idre.ucla.edu/r/faq/ologit-coefficients/
exp(coef(m)) ## odds ratios
exp(cbind(OR = coef(m), ci)) ## OR and CI
