read "../ComputeIdentifiableFunctionsRG.mpl"; #"../IdentifiabilityODE.mpl";
kernelopts(printbytes=false, assertlevel=1):
interface(echo=0, prettyprint=0):
model :=[
  diff(S(t), t) = -b * In(t) * S(t),
  diff(In(t), t) = b * In(t) * S(t) - g * In(t),
  diff(R(t), t) = g * In(t),
  diff(aux(t), t) = 0,
  y1(t) = In(t),
  y2(t) = b / g + aux(t)
]:

me:= MultiExperimentIdentifiableFunctions(model, simplified_generators=true, no_bound=true): print(me[3]);#IdentifiabilityODE(sigma, [aux(0)], infolevel = 3):
