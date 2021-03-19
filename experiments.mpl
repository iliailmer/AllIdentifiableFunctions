
pth:= "/Users/iliailmer/Documents/CUNY/Research/AlgebraicGeom/AllIdentifiableFunctions/TDDS/TDDS.lib":
libname := pth, libname;

# with(TDDS);
# with(DifferentialThomas);
# infolevel[TDDS] := 2;
# infolevel[DifferentialThomas]:=2;
read "ComputeIdentifiableFunctions.mpl";

kernelopts(printbytes=false, assertlevel=1):
interface(echo=0, prettyprint=0):
model :=[
  diff(x1(t), t) = -(k3 + k7) * x1(t) + k4 * x2(t),
  diff(x2(t), t) = k3 * x1(t) - (k4 + a(t) * k5 + b(t) * d(t) * k5) * x2(t) + k6 * x3(t) + k6 * x4(t) + k5 * x2(t) * x3(t) + k5 * x2(t) * x4(t),
  diff(x3(t), t) = a(t) * k5 * x2(t) - k6 * x3(t) - k5 * x2(t) * x3(t),
  diff(x4(t), t) = b(t) * d(t) * k5 * x2(t) - k6 * x4(t) - k5 * x2(t) * x4(t),
  diff(x5(t), t) = k7 * x1(t),
  diff(a(t), t) = 0,
  diff(b(t), t) = 0,
  diff(d(t), t) = 0,
  y1(t) = x5(t),
  y2(t) = a(t),
  y3(t) = b(t),
  y4(t) = d(t)
]:
states, inputs, outputs, params, model_eqs := op(ParseInput(model)):
model_denomfree := ExtractDenominator(model_eqs):

io_eqs := GetIOEquations(model_denomfree, states, inputs, outputs, params, 1):
print(io_eqs);
if nops(inputs)>0 then
    state_outs_inputs := [[op(states)], [op(outputs)], [op(inputs)]]:
else
    state_outs_inputs := [[op(states)], [op(outputs)]]:
end if:

io_jets:= map(DifferentialThomas:-Tools:-ToJet, io_eqs):
DifferentialThomas:-Ranking([t], state_outs_inputs):
leaders := map(DifferentialThomas:-Leader, io_jets);
degrees := seq(degree(io_jets[i], leaders[i]), i=1..nops(io_jets));


read "ComputeIdentifiableFunctionsRG.mpl";

io_eqs := GetIOEquations(model_denomfree, states, inputs, outputs, params, 1):


io_jets:=map(x->DifferentialAlgebra:-Tools:-ToJet(x, states), io_eqs):
Rorig := DifferentialRing(blocks = [[op(outputs)], [op(states)], [op(inputs)]], derivations = [t], arbitrary = params):
leaders := map(x->DifferentialAlgebra:-Tools:-LeadingDerivative(x, Rorig), io_jets);
degrees := seq(degree(io_jets[i], leaders[i]), i=1..nops(io_jets));


