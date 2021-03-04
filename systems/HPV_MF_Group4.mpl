read "../IdentifiabilityODE_dev.mpl":

eq[1] := (i,j) -> diff(S[i](t),t) = mu/2 + gammaG[i]*IG[i](t) + gammaO[i]*IO[i](t) - S[i](t)*mu 
                                - S[i](t)*(betaOO[j,i] * (IO[j](t) + IOG[j](t)) + betaGO[j,i] * (IG[j](t) + IOG[j](t)))
                                - S[i](t)*(betaOG[j,i] * (IO[j](t) + IOG[j](t)) + betaGG[j,i] * (IG[j](t) + IOG[j](t))):

eq[2] := (i,j) -> diff(IO[i](t),t) = S[i](t) * (betaOO[j,i] * (IO[j](t) + IOG[j](t)) + betaGO[j,i] * (IG[j](t) + IOG[j](t))) + gammaG[i]*IOG[i](t)
                             - IO[i](t)*(nuOG[i] + gammaO[i] + mu + betaOG[j,i] * (IO[j](t) + IOG[j](t)) + betaGG[j,i] * (IG[j](t) + IOG[j](t))):

eq[3] := (i,j) -> diff(IG[i](t),t) = S[i](t) * (betaOG[j,i] * (IO[j](t) + IOG[j](t)) + betaGG[j,i] * (IG[j](t) + IOG[j](t))) + gammaO[i]*IOG[i](t)
                             - IG[i](t)*(nuGO[i] + gammaG[i] + mu + betaOO[j,i] * (IO[j](t) + IOG[j](t)) + betaGO[j,i] * (IG[j](t) + IOG[j](t))):

eq[4] := (i,j) -> diff(IOG[i](t),t) = IO[i](t)*(nuOG[i] + betaOG[j,i]* (IO[j](t) + IOG[j](t)) + betaGG[j,i]*(IG[j](t) + IOG[j](t)))
                                  + IG[i](t)*(nuGO[i] + betaOO[j,i]* (IO[j](t) + IOG[j](t)) + betaGO[j,i]*(IG[j](t) + IOG[j](t)))
                                  - IOG[i](t)*(gammaO[i] + gammaG[i] + mu):

sigma := [seq(eq[i](M,F),i=1..4), 
          seq(eq[i](F,M),i=1..4),
          y1(t)=IG[M](t) + IOG[M](t), 
          y2(t)=IO[M](t) + IOG[M](t),
          y3(t)=IG[F](t) + IOG[F](t),
          y4(t)=IO[F](t) + IOG[F](t)
 ]:

with(StringTools):
sigmastr := map(convert,sigma,string):
sigmastrsubs := map2(Subs, {"[F]"="F","[M]"="M","[M,F]"="MF","[F,M]"="FM"}, sigmastr):
sigmaconv := map(parse,sigmastrsubs):

IdentifiabilityODE(sigmaconv, GetParameters(sigmaconv));
