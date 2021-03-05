# Example (with initial conditions assumed being unknown) from Section IV of "DAISY: an Efficient Tool to Test Global Identifiability. Some Case Studies"
# by G. Bellu, M.P. Saccomani

read "../ComputeIdentifiableFunctionsRG.mpl"; #"../IdentifiabilityODE.mpl";
kernelopts(printbytes=false, assertlevel=1):
interface(echo=0, prettyprint=0):
model :=[
  diff(x1(t), t) = -beta * x1(t) * x4(t) - d * x1(t) + s,
  diff(x2(t), t) = beta * q1 * x1(t) * x4(t) - k1 * x2(t) - mu1 * x2(t),
  diff(x3(t), t) = beta * q2 * x1(t) * x4(t) + k1 * x2(t) - mu2 * x3(t),
  diff(x4(t), t) = -c * x4(t) + k2 * x3(t),
  y1(t) = x1(t),
  y2(t) = x4(t)
]:

me:= MultiExperimentIdentifiableFunctions(model, simplified_generators=true, no_bound=true): print(me[3]);#IdentifiabilityODE(sigma, GetParameters(sigma));
