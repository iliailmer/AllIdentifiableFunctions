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