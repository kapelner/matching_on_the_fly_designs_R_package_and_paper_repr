
library(survival)
set.seed(1)
n = 100
time = rexp(n)
status = rbinom(n, 1, 0.7)
arm = rbinom(n, 1, 0.5)
# Ensure max times are different
time[arm == 0][which.max(time[arm == 0])] = 10
time[arm == 1][which.max(time[arm == 1])] = 5

obj = survfit(Surv(time, status) ~ arm)
res = summary(obj)$table
print(res)
print(max(time[arm==0]))
print(max(time[arm==1]))

# Check if we can extract the cutoff used?
# print(obj) gives some info
print(obj)
