kernelopts(printbytes=false, assertlevel=1):
interface(echo=0, prettyprint=0):
pth:= "/Users/iliailmer/Documents/CUNY/Research/AlgebraicGeom/AllIdentifiableFunctions/TDDS/TDDS.lib":
libname := pth, libname;

with(TDDS);
with(DifferentialThomas);
infolevel[TDDS] := 2;
infolevel[DifferentialThomas]:=2;


# sigma := [
#   diff(x(t), t) = lm - d * x(t) - beta * x(t) * v(t),
#   diff(y(t), t) = beta * x(t) * v(t) - a * y(t),
#   diff(v(t), t) = k * y(t) - u * v(t),
#   diff(w(t), t) = c * z(t) * y(t) * w(t) - c * q * y(t) * w(t) - b * w(t),
#   diff(z(t), t) = c * q * y(t) * w(t) - h * z(t),
#   y1(t) = w(t),
#   y2(t) = z(t)
# ];

sigma := [
  diff(S(t), t) = b * S(t) * In(t) / N,
  diff(E(t), t) = b * S(t) * In(t) / N - nu * E(t),
  diff(In(t), t) = nu * E(t) - a * In(t),
  y1(t) = In(t),
  y2(t) = N
]:

ComputeRanking([t], [[S,E,In], [y1, y2]]);
tdds_out := DifferentialThomasDecomposition(map(x->lhs(x) - rhs(x), sigma), [], "stop"=1): # TDDS not using equations

Ranking([t], [[S,E,In], [y1, y2]]);
dt_out := ThomasDecomposition(sigma, [], "stop"=1); # Equations():

# Rosenfeld Groebner Code

Rorig := DifferentialRing(blocks = [[y1, y2], [S,E,In]], derivations = [t], arbitrary = []):
chset_orig := RosenfeldGroebner(model, Rorig)[1]

(*
[
diff(y2(t),t)*eB+k2*(xa(t)*y3(t)+y2(t)*y4(t)-y1(t)),

k2*xb(t)-diff(y2(t),t), 

-y2(t)+xc(t),
((-k2*y4(t)-eB*(k1-k2))*y3(t)+eB^2*k1)*diff(y2(t),t)-(-y3(t)*diff(y1(t),t)+k1*(-y2(t)*y4(t)+y1(t))*(eB-y3(t)))*k2,

diff(y2(t),t $ 2)*y3(t)+(eB*k1+k2*y3(t))*diff(y2(t),t)-k1*k2*(-y2(t)*y4(t)+y1(t))
]


[k1*k2*y1(t)-diff(y2(t),t $ 2)*y3(t)-k2*diff(y2(t),t)*y3(t)-eB*k1*diff(y2(t),t)-k1*k2*y2(t)*y4(t), 
diff(y2(t),t $ 3)+diff(y2(t),t $ 2)*k1+diff(y2(t),t $ 2)*k2+diff(y2(t),t)*k1*k2]*)
