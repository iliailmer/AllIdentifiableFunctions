# The system is taken from
# Wodarz, D., Nowak, M.
# Mathematical models of HIV pathogenesis and treatment  
# System (6)
read "../ComputeIdentifiableFunctions.mpl"; #"../IdentifiabilityODE.mpl":
kernelopts(printbytes=false, assertlevel=1):
interface(echo=0, prettyprint=0):
model :=[
  diff(x(t), t) = lm - d * x(t) - beta * x(t) * v(t),
  diff(y(t), t) = beta * x(t) * v(t) - a * y(t),
  diff(v(t), t) = k * y(t) - u * v(t),
  diff(w(t), t) = c * z(t) * y(t) * w(t) - c * q * y(t) * w(t) - b * w(t),
  diff(z(t), t) = c * q * y(t) * w(t) - h * z(t),
  y1(t) = w(t),
  y2(t) = z(t)
]:

me:= MultiExperimentIdentifiableFunctions(model, simplified_generators=true, no_bound=true, infolevel=1): 
print(me[3]);#IdentifiabilityODE(sigma, GetParameters(sigma)):
