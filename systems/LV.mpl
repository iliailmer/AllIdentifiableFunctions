read "../ComputeIdentifiableFunctionsRG.mpl"; #"../IdentifiabilityODE.mpl";
kernelopts(printbytes=false, assertlevel=1):
interface(echo=0, prettyprint=0):
model :=[
  diff(x1(t), t) = a * x1(t) - b * x1(t) * x2(t),
  diff(x2(t), t) = -c * x2(t) + d * x1(t) * x2(t),
  y(t) = x1(t) + u(t)
]:

me:= MultiExperimentIdentifiableFunctions(model, simplified_generators=true, no_bound=true): print(me[3]);#IdentifiabilityODE(sigma, [a, b, c, d, x1(0), x2(0)]);
