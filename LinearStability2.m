clc,clear all,close all
%% Linear stability in this code is written to find stability of the fixed point

syms vpl b drs kappa a k
syms v theta t th_0 

f = @(y) [
    y(1) * ((-k * (y(1) - vpl) - b * ((1 / y(2)) - (y(1) / drs))) / (kappa * y(1) + a));
    1 - (y(1) * y(2)) / drs
];

y=[v,theta];
F=f(y);
J11=diff(F(1),v);
J12=diff(F(1),theta);
J21=diff(F(2),v);
J22=diff(F(2),theta);

J11=subs(J11,[v,theta,kappa],[vpl,drs/vpl,0]);
J12=subs(J12,[v,theta,kappa],[vpl,drs/vpl,0]);
J21=subs(J21,[v,theta,kappa],[vpl,drs/vpl,0]);
J22=subs(J22,[v,theta,kappa],[vpl,drs/vpl,0]);
J=[J11,J12;J21 J22];
eigens=simplify(eig(J));

pretty(eigens(1))
pretty(eigens(2))
