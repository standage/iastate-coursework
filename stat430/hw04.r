# 1. Load the data
hiv <- read.table(file("hiv_status_data.Rtxt"))
c1 <- subset(hiv, Clinic.ID == 1)
c2 <- subset(hiv, Clinic.ID == 2)
c1.Y <- sort(c1$Y)
c2.Y <- sort(c2$Y)
c1.n <- length(c1.Y)
c2.n <- length(c2.Y)

# 1.a - Assess normality of data
#       Create normal probability plots manually
c1.seq <- seq(from=1, to=c1.n, by=1)/(c1.n+1)
c2.seq <- seq(from=1, to=c2.n, by=1)/(c2.n+1)
c1.q <- qnorm(c1.seq)
c2.q <- qnorm(c2.seq)
plot(c1.q, c1.Y)
plot(c2.q, c2.Y)

#       Create plots with built-in functions
png(file="qqplot-c1-Y.png")
qqnorm(c1.Y, main="QQ-Plot of Y Values (Clinic 1)")
qqline(c1.Y, col="red")
dev.off()
png(file="qqplot-c2-Y.png")
qqnorm(c2.Y, main="QQ-Plot of Y Values (Clinic 2)")
qqline(c2.Y, col="red")
dev.off()

#       1.b - Demonstrate relationship between X and Var(X)
lower <- quantile(hiv$Y, probs=.25, names=FALSE)
higher <- quantile(hiv$Y, probs=.75, names=FALSE)
Y.small <- hiv$Y[hiv$Y > lower & hiv$Y < higher]
Y.big <- hiv$Y[hiv$Y <= lower | hiv$Y >= higher]
length(Y.small)
length(Y.big)
var(Y.small)
var(Y.big)

#       1.c - Perform log transformation
c1.Y.transformed <- log(abs(c1.Y))
png(file="c1Y-vs-c1Ytransformed.png")
plot(c1.Y, c1.Y.transformed, main="Clinic 1 Response Values", xlab="Y", ylab="log(Y)")
dev.off()
c2.Y.transformed <- log(abs(c2.Y))
png(file="c2Y-vs-c2Ytransformed.png")
plot(c2.Y, c2.Y.transformed, main="Clinic 2 Response Values", xlab="Y", ylab="log(Y)")
dev.off()

#       1.d - Test if clinic impacts response (transformed data)
t.test(c1.Y.transformed, c2.Y.transformed)
var.test(c1.Y.transformed, c2.Y.transformed)

#       1.e - Non-parametric test of clinic impact on response (non-transformed data)
wilcox.test(c1.Y, c2.Y, paired=F)

# 2. Load the data
hiv <- read.table(file("hiv_status_data.Rtxt"))
c1 <- subset(hiv, Clinic.ID == 1)
c2 <- subset(hiv, Clinic.ID == 2)
c1.trimmed <- subset(c1, Z != 0 & !is.na(Z))
c1.trimmed.n <- length(c1.trimmed$Z)
c2.trimmed <- subset(c2, Z != 0 & !is.na(Z))
c2.trimmed.n <- length(c2.trimmed$Z)

#       2.a - Power calculation
n <- c1.trimmed.n + c2.trimmed.n
delta <- 0.0003
sigma <- 0.003
alpha <- 0.05
1 - pnorm( qnorm(1-(alpha/2)) - delta/(sigma*sqrt(2/n)) ) + pnorm( -qnorm(1-(alpha/2)) - delta/(sigma*sqrt(2/n)) )

#       2.b - Test clinic effect on adherence
#             Test normality of adherence (Z) values
c1.trimmed.data <- sort(c1.trimmed$Z)
c2.trimmed.data <- sort(c2.trimmed$Z)
png(file="qqplot-c1-Z.png")
qqnorm(c1.trimmed.data, main="QQ-Plot of Z Values (Clinic 1)")
qqline(c1.trimmed.data, col="red")
dev.off()
png(file="qqplot-c2-Z.png")
qqnorm(c2.trimmed.data, main="QQ-Plot of Z Values (Clinic 2)")
qqline(c2.trimmed.data, col="red")
dev.off()

#             Test normality of log(Z) values
c1.trimmed.data.log <- log(c1.trimmed.data)
c2.trimmed.data.log <- log(c2.trimmed.data)
png(file="qqplot-c1-Z-log.png")
qqnorm(c1.trimmed.data.log, main="QQ-Plot of Log-Z Values (Clinic 1)")
qqline(c1.trimmed.data.log, col="red")
dev.off()
png(file="qqplot-c2-Z-log.png")
qqnorm(c2.trimmed.data.log, main="QQ-Plot of Log-Z Values (Clinic 2)")
qqline(c2.trimmed.data.log, col="red")
dev.off()

#             Test clinic effect on log(Z) values
t.test(c1.trimmed.data.log, c2.trimmed.data.log)
var.test(c1.trimmed.data.log, c2.trimmed.data.log)

#       2.c - Non-parametric test of clinic effect on adherence
T.hat <- mean(c1.trimmed.data) - mean(c2.trimmed.data)
data.combined <- c(c1.trimmed.data, c2.trimmed.data)
perm <- function(B=1000) {
  T <- c()
  for(i in 1:B) {
    pseudo.c1.positions <- sample(1:length(data.combined), size=length(c1.trimmed.data))
    pseudo.c1.data <- data.combined[pseudo.c1.positions]
    pseudo.c2.data <- data.combined[-pseudo.c1.positions]
    T[i] <- mean(pseudo.c1.data) - mean(pseudo.c2.data)
  }
  p.value <- sum(abs(T) > abs(T.hat))/B
  return(p.value)
}
perm(10000)

# 3 - Find confidence interval for difference in location
hiv <- read.table(file("hiv_status_data.Rtxt"))
c1 <- subset(hiv, Clinic.ID == 1)
c2 <- subset(hiv, Clinic.ID == 2)
c1.Y <- sort(c1$Y)
c2.Y <- sort(c2$Y)
central.delta <- optimize(function(d) wilcox.test(c1.Y+d, c2.Y)$p.value, c(-.0005, .0005), maximum=T)$maximum
uniroot(function(d) wilcox.test(c1.Y+d, c2.Y)$p.value - 0.05, c(-.5, central.delta), tol=0.00000001)$root
uniroot(function(d) wilcox.test(c1.Y+d, c2.Y)$p.value - 0.05, c(central.delta, .5), tol=0.00000001)$root

# 4 - Bootstrap confidence intervals for clinic effect on response and adherence
#     Load data
hiv <- read.table(file("hiv_status_data.Rtxt"))
c1 <- subset(hiv, Clinic.ID == 1)
c2 <- subset(hiv, Clinic.ID == 2)
c1.Y <- sort(c1$Y)
c2.Y <- sort(c2$Y)
c1.n <- length(c1.Y)
c2.n <- length(c2.Y)

#     Boostrap CI for response (Y) values
B <- 1000
sample.meandiff <- mean(c1.Y) - mean(c2.Y)
sample.pooled.sd <- sqrt( ((c1.n-1)*var(c1.Y) + (c2.n)*var(c2.Y))/(c1.n + c2.n - 2) )
c1.rep <- matrix(c1.Y, nrow=B, ncol=c1.n, byrow=T)
c1.boot <- apply(c1.rep, 1, sample, replace=T)
c1.boot.mean <- apply(c1.boot, 2, mean)
c1.boot.sd <- apply(c1.boot, 2, sd)
c2.rep <- matrix(c2.Y, nrow=B, ncol=c2.n, byrow=T)
c2.boot <- apply(c2.rep, 1, sample, replace=T)
c2.boot.mean <- apply(c2.boot, 2, mean)
c2.boot.sd <- apply(c2.boot, 2, sd)
boot.meandiff <- c1.boot.mean - c2.boot.mean
boot.pooled.sd <- sqrt( ((c1.n-1)*c1.boot.sd^2 + (c2.n-1)*c2.boot.sd^2)/(c1.n + c2.n - 2) )
z.boot <- (boot.meandiff - sample.meandiff)/(boot.pooled.sd * sqrt((1/c1.n)+(1/c2.n)))
mean(boot.meandiff) + quantile(z.boot, prob=0.025)*sample.pooled.sd/sqrt((1/c1.n)+(1/c2.n))
mean(boot.meandiff) + quantile(z.boot, prob=0.975)*sample.pooled.sd/sqrt((1/c1.n)+(1/c2.n))

