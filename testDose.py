"General script to check physical dose versus BED, changes often depending on model"""

import math
# input information
fractions = 1.0
duration = 1210.4 / (3600.0 * 24.0)
PhysicalDose = 15
doses = [8.5, 6.5]

time = (1210.4/2) / (3600.0 * 24.0)
t1 = (330.93/2) / (3600.0 * 24.0)
t2 = (909.47/2) / (3600.0 * 24.0)
times = [t1, t2]

# hardcoded biological factors
alpha = 0.15
beta = 0.048
mu = 11.09
Tp = 61.6
Tk = 0.0

# simple version of BED for fractions of equal dose; does not account for proliferation
BEDsimple = (PhysicalDose*fractions)*(1+(PhysicalDose/(alpha/beta)))

# formula from Chen's poster; full BED that accounts for proliferation
g_factor = 2 * (math.exp(-mu*time)+(mu*time)-1) / ((mu*time)**2)
RE = 1 + (((PhysicalDose/2.0)/(alpha/beta))*g_factor)
cellProliferation = (0.693/(alpha*Tp))*(duration-Tk)
BEDfull = ((PhysicalDose*fractions)*RE) - cellProliferation

# formula that accounts for different fractional doses and cell proliferation
g = 0.0
for j in range(0,N):
    factor_1 = (doses[j]**2)*(times[j] - ((1/mu)*(1-math.exp(-mu*times[j]))))/(times[j]**2)
    factor_2 = 0.0
    for k in range(0,j):
        factor_2 += (doses[k]*doses[j]*((math.exp(-mu*times[0]))*(math.exp(mu*times[k]) - 1)*(math.exp(-mu*times[j])-1)))/(times[k]*times[j])
    g += factor_1 - ((1/mu)*factor_2)
    g_factor = (2/mu)*g
    prolif = ((0.693/(alpha*Tp))*(duration-Tk))
    BED = (D + (1/((alpha/beta)) * g_factor)) - prolif

# difference between the full model and simple model
difference = abs(BEDtest - BEDfull)

print "Test BED calculation: ", BEDtest
print "Full BED calculation: ", BEDfull
