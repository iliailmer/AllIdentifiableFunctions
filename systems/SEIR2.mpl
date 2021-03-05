# Taken from
# N. Tuncer, T. Le
# "Structural and practical identifiability analysis of outbreak models"
# https://doi.org/10.1016/j.mbs.2018.02.004
# Equation (2.2) with cumulative incidence observations
read "../ComputeIdentifiableFunctionsRG.mpl"; #"../IdentifiabilityODE.mpl":
kernelopts(printbytes=false, assertlevel=1):
interface(echo=0, prettyprint=0):
model :=[
  diff(S(t), t) = -b * S(t) * In(t) / N(t),
  diff(E(t), t) = b * S(t) * In(t) / N(t) - nu * E(t),
  diff(In(t), t) = nu * E(t) - a * In(t),
  diff(N(t), t) = 0,
  diff(Cu(t), t) = nu * E(t),
  y1(t) = Cu(t),
  y2(t) = N(t)
]:

me:= MultiExperimentIdentifiableFunctions(model, simplified_generators=true, no_bound=true): print(me[3]);#IdentifiabilityODE(sigma, GetParameters(sigma)):
