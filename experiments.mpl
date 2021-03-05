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
  diff(xA(t), t) = -k1 * xA(t),
  diff(xB(t), t) = k1 * xA(t) - k2 * xB(t),
  diff(xC(t), t) = k2 * xB(t),
  diff(eA(t), t) = 0,
  diff(eC(t), t) = 0,
  y1(t) = xC(t),
  y2(t) = eA(t) * xA(t) + eB * xB(t) + eC(t) * xC(t),
  y3(t) = eA(t),
  y4(t) = eC(t)
]:

ComputeRanking([t], [[eA, eC, xA, xB, xC], [y1, y2, y3, y4]]);
tdds_out := DifferentialThomasDecomposition(map(x->lhs(x) - rhs(x), sigma), [], "stop"=1): # TDDS not using equations

Ranking([t], [[eA, eC, xA, xB, xC], [y1, y2, y3, y4]]);
dt_out := ThomasDecomposition(sigma, [], "stop"=1); # Equations():

# Rosenfeld Groebner Code

Rorig := DifferentialAlgebra:-DifferentialRing(blocks = [ [y1, y2, y3, y4], [eA, eC, xA, xB, xC]], derivations = [t], arbitrary = [eb, k1, k2]):
chset_orig := DifferentialAlgebra:-RosenfeldGroebner(map(x->lhs(x)-rhs(x), sigma), Rorig, singsol=none)[1];

