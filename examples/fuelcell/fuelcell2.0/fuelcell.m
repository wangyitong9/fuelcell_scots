
function fuelcell()
clear set
close all

%% simulation
% x0 is a 2d vector representing initial state
x0=[296 0.6];
% load the symbolic set containing the controller
C=SymbolicSet('reach_set5.bdd');
target=SymbolicSet('system_target.bdd');
y=x0;
v=[];
T=100; 
refine=4;
RelTol = 10^-7;
AbsTol = 10^-8; 
% for t=1:T/0.5
% 
%   u=C.getInputs(y(end,:));
%   v=[v; u(end,:)];
%   [t x]=ode45(@fuelcell_ode,[0 .5], y(end,:), odeset('abstol',1e-4,'reltol',1e-4),u(end,:));
% 
%   y=[y; x(end,:)];
% end
options = odeset('Refine',refine, 'RelTol',RelTol,'AbsTol',AbsTol);
%  while(1)
% 
%   if (target.isElement(y(end,:))) 
%     break;
%   end 
% 
%   u=C.getInputs(y(end,:));
%   v=[v; u(1,:)];
%   [t x]=ode45(@fuelcell_ode,[0 .5], y(end,:), [],u(1,:));
% 
%   y=[y; x(end,:)];
% end

%% plot the backward reachable set
% colors
colors=get(groot,'DefaultAxesColorOrder');

% plot the domain of the controller
set=SymbolicSet('reach_set5.bdd','projection',[1 2 3]);
x=set.points;
y=[];
% for i = 1:30506
%     if x(i,3)==170
%         y= [y; x(i,1), x(i,2)];
%     end
% end
% plot(y(:,1),y(:,2),'.','color',[0.8 0.8 0.8])
% axis([273 293 0.99 1])
plot3(x(:,1),x(:,2),x(:,3),'.','color',[0.8 0.8 0.8])
hold on

% plot initial state  and trajectory
% plot(y(:,1),y(:,2),'.-','color',colors(1,:))
% plot(y(1,1),y(1,1),'.','color',colors(5,:),'markersize',20)
% 
% plot target set
% v=[1.15 5.45; 1.55 5.45; 1.15 5.85; 1.55 5.85 ];
% v=[280 0.996; 290 0.996; 280 0.997; 290 0.997 ];
% patch('vertices',v,'faces',[1 2 4 3],'facecolor','none','edgec',colors(2,:),'linew',1)
% hold on
% 
% 
% v=[273 0.993; 293 0.993; 273 1; 293 1 ];
% patch('vertices',v,'faces',[1 2 4 3],'facecolor','none','edgec',colors(5,:),'linew',1)
% hold on
% 
% box on
% axis([273 293 0.993 1])

end

function dxdt = fuelcell_ode(t,x,u)
x1=x(1);
x2=x(2);
u1=u(1);
u2=u(2);
dxdt0 = 1 / 100 * ((2.26e-4*x1 + 2 * 8.3145*x1 / 96485 * asinh((u1 + 0.001)*1e5 / 16) + 0.05*u1 - 0.12*log(1 - u1 / 0.7))*u1 * 150 - 1 / 26 * x1 * pow(u2 / 3, 1 / 1.35));
dxdt1 = (1.2*x2 + 3 - sqrt(pow(1.2*x2 + 3, 2) - 4 * 0.18 / 24 * (3 + 4.65e-3*pow(u2, 2) / 3 / 125 * 1000 / 60 - u1 * 150 *3 * (1.256 - 2.26e-4*x1 - 2 * 8.3145*x1 / 96485 * asinh((u1 + 0.001)*1e5 / 16) - 0.05*u1 + 0.12*log(1 - u1 / 0.7)))))/3600/1.875/2/0.18;

dxdt= [dxdt0;dxdt1];
end
function out = pow(base, exp)
out=base^exp;
end

