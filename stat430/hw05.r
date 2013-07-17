# 1.a Boxplot
rev <- read.table("RevMutant.Rtxt", header=T)
R1 <- subset(rev, Mutant == "R1")
png(filename="R1_boxplot.png")
boxplot(pg ~ Experiment, data=R1)
dev.off()

# 1.b Test homoscedasticity
bartlett.test(x=R1$pg, g=R1$Experiment)
library("car")
leveneTest(y=R1$pg, group=R1$Experiment)

# 1.c Perform ANOVA
R1.anova <- aov(pg ~ Experiment, data=R1)
summary(R1.anova)
kruskal.test(R1$pg, g=R1$Experiment)

# 2.a Estimate mean activity
rev <- read.table("RevMutant.Rtxt", header=T)
mean.pg <- as.vector(tapply(rev$pg, rev$Mutant, mean))

# 2.b Compute CIs with Tukey method
anova <- aov(pg ~ Mutant, data=rev)
summary(anova)
names <- sort(as.vector(unique(rev$Mutant)))
n.pg <- as.vector(tapply(rev$pg, rev$Mutant, length))
df <- 77
k <- 10
s.pool <- sqrt(158.9)
alpha <- 0.05
results <- matrix(ncol=6, nrow=45)
row = 1
for(i in 1:(k-1))
{
  for(j in (i+1):k)
  {
    results[row, 1] <- names[i]
    results[row, 2] <- names[j]
    results[row, 3] <- round((mean.pg[i] - mean.pg[j]) - qtukey(1-alpha, nmeans=k, df=df)*s.pool*sqrt((1/n.pg[i])+(1/n.pg[j]))/sqrt(2), 4)
    results[row, 4] <- round((mean.pg[i] - mean.pg[j]) + qtukey(1-alpha, nmeans=k, df=df)*s.pool*sqrt((1/n.pg[i])+(1/n.pg[j]))/sqrt(2), 4)
    t.confint <- t.test(rev$pg[rev$Mutant == names[i]], rev$pg[rev$Mutant == names[j]])$conf.int
    results[row, 5] <- round(t.confint[1], 4)
    results[row, 6] <- round(t.confint[2], 4)
    row <- row+1
  }
}
head(results)

# 2.d
# subset(results, results[,1] == "R1" | results[,2] == "R1")
names
R1.index <- 5
signif <- qtukey(1-alpha, nmeans=k, df=df)
sig <- round(signif, 2)
m.R1 <- mean.pg[R1.index]
m.R1.rnd <- round(m.R1,2)
for(i in 1:k)
{
  m.i <- mean.pg[i]
  meandiff <- m.R1 - m.i
  if(meandiff > signif)
  {
    m.i.rnd <- round(m.1,2)
    m.diff <- round(meandiff, 2)
    cat(names[i]," (", m.R1.rnd, "-", m.i.rnd, "=", m.diff, " > ", sig, ")\n")
  }
}

