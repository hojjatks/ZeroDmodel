clc,clear all,close all
%% Linear stability in this code is written to find stability with respect to v=0 and theta=theta_0>0

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

J11=subs(J11,[v,theta],[0,t+th_0]);
J12=subs(J12,[v,theta],[0,t+th_0]);
J21=subs(J21,[v,theta],[0,t+th_0]);
J22=subs(J22,[v,theta],[0,t+th_0]);
J=[J11,J12;J21 J22];
pretty(J)
pretty(eig(J))