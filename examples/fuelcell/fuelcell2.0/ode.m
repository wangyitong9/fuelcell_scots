function Prob6_ODE45(t_end, u)
%
refine=4;
RelTol = 10^-7;
AbsTol = 10^-8;
options = odeset('Refine',refine, 'RelTol',RelTol,'AbsTol',AbsTol);
t_start=0;
[A,r]=myfunc(360,370,0.6,0.75);
sz=size(A);
line=sz(1);
[tout, xout] = ode45(@f, [t_start t_end], [300 0.15], options, u)
function [A, r] = myfunc(xl, xu, yl, yu)
A = [];
r = [];
midx = (xl+xu)/2;
midy = (yl+yu)/2;
dx = (xu-xl)/9;
dy = (yu-yl)/9;
for x = xl:dx:xu 
for y = yl:dy:yu
A = [A; x y];
r = [r; abs(midx-x) abs(midy-y)];
end
end
function out = pow(base, exp)
out=base^exp;
function dxdt=f(t,x,u)
x1=x(1);
x2=x(2);
u1=u(1);
u2=u(2);
%dxdt0 = 2*x1+1;
%dxdt0=1/100*(2.26e-4*u1*150+2*83145/96485*asinh((u1 + 0.001)*1e5 / 16)*u1*150-1/26*pow(u2 / 3, 1 / 1.35))*x1+(0.05*u1/100-0.12*log(1 - u1 / 0.7))*u1*150;
dxdt0 = 1 / 100 * ((2.26e-4*x1 + 2 * 8.3145*x1 / 96485 * asinh((u1 + 0.001)*1e5 / 16) + 0.05*u1 - 0.12*log(1 - u1 / 0.7))*u1 * 150 - 1 / 26 * x1 * pow(u2 / 3, 1 / 1.35));
%dxdt0 = 1 / 100 * ((2.26e-4*x1+ 2 * 8.3145*x1 / 96485 * asinh((u1 + 0.001)*1e5 / 16))*u1 * 150);
dxdt1 = (1.2*x2 + 3 - sqrt(pow(1.2*x2 + 3, 2) - 4 * 0.18 / 24 * (3 + 4.65e-3*pow(u2, 2) / 3 / 125 * 1000 / 60 - u1 * 150 *3 * (1.256 - 2.26e-4*x1 - 2 * 8.3145*x1 / 96485 * asinh((u1 + 0.001)*1e5 / 16) - 0.05*u1 + 0.12*log(1 - u1 / 0.7)))))/3600/1.875/2/0.18;
%
%if 0 % Simulate using linear model 
%    dx1=x2;
%    dx2=-x1-x2+u
%    %
%else  % Simulate NL model
%    dx1=x2;
%    dx2=-x1-sin(x2)+u;    
%end
%
dxdt=[dxdt0;dxdt1];
return
function drdt = growthbound(t, x, u)
u0=u(1);
u1=u(2);
r0=x(1);
r1=x(2);
x0=400;
x1=0;
Ncells=3;
Np=3;
Ns=8;
a00 = (2.26e-4 + 2 * 8.3145 / 96485 * asinh((u0 + 0.001) / 16*1e5))*u0 * 150/100- 1 / 26 * pow(u1 / Ncells, 1 / 1.35)/100;
a01 = 0;
a10 = 1 / 4/0.18/3600/1.875 * pow(pow(1.2*x1 + 3, 2) - 4 * 0.18 / Np/Ns * (Ncells + 4.65e-3*pow(u1, 2) / Ncells / 125 * 1000 / 60 - u0 * 150 * Ncells * (1.256 - 2.26e-4*x0 - 2 * 8.3145*x0 / 96485 * asinh((u0 + 0.001)*1e5 / 16) - 0.05*u0 + 0.12*log(1 - u0 / 0.7))), -1 / 2) * 4 * 0.18 / Np/Ns * u0 * 150*Ncells*(2.26e-4 + 2 * 8.3145 / 96485 * asinh((u0 + 0.001) / 16 * 1e5));
a11 = 1.2/3600/1.875/2/0.18; %- 1 / 2/3600/1.875/2/0.18 * pow(pow(1.2*r1 + 3, 2) - 4 * 0.18 / Np/Ns * (Ncells + 4.65e-3*pow(u1, 2) / Ncells / 125 * 1000 / 60 - u0 * 150 * Ncells * (1.256 - 2.26e-4*r0 - 2 * 8.3145*r0 / 96485 * asinh((u0 + 0.001)*1e5 / 16) - 0.05*u0 + 0.012*log(1 - u0 / 0.7))), -1 / 2) * 2 * 1.2*(1.2*r1 + 3);
drdt0 = a00 * r0 + a01 * r1;
drdt1 = a10 * r0 + a11 * r1+0.05;
drdt = [drdt0;drdt1];



