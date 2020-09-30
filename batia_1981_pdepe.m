function batia_1981_pdepe
% From:
% Bhatia SK, Perlmutter DD, 
% "A random pore model for fluid-solid reactions: II Diffusion and transport effects", 
% AIChE Journal, 27(2) (1981), 247-254.
psi=1;
Z=1; 
% looks like when Z=1, eps0 is irrelevant
eps0=0.3;
beta=1;
phi=3;
Sh=100;
%Sh=0; % dirichlet constraint

etamesh=linspace(0,1,50);
%tspan=linspace(0,2,10);
tspan=[0 .15 .9 1.65 2.4 3.15 5.4]';
m=2;
icf = @(x) icfun(x);
bcf = @(xl,ul,xr,ur,t) bcfun(xl,ul,xr,ur,t,Sh);
pdef = @(x,t,u,dudx) pdefun(x,t,u,dudx,Z,eps0,psi,beta,phi);
sol=pdepe(m,pdef,icf,bcf,etamesh,tspan);
cstar=sol(:,:,1);
X=sol(:,:,2);

figure; plot(tspan, cstar(:,1)); grid
title 'C at the center as a function of time time';

% Figure 4 from reference
plotResultAtTimes(etamesh, cstar, tspan, ...
  'C at selected time points');
% Figure 5 from reference
plotResultAtTimes(etamesh, X, tspan, ...
  'X at selected time points');

end

function [c,f,s] = pdefun(x,t,u,dudx,Z,eps0,psi,beta,phi)
CStar=u(1,:);
X=u(2,:);
c=[0; 1];
DEffStar=1-((Z-1)*(1-eps0)/eps0)*X;
f=[DEffStar; 0].*dudx;
ss=sqrt(1-psi*log(1-X));
dXDt=CStar.*(1-X).*ss./(1+(beta*Z/psi.*ss-1));
s=[-phi^2*dXDt; dXDt];
end

function [pl,ql,pr,qr] = bcfun(xl,ul,xr,ur,t,Sh)
pl = [0 0]';
ql = [1 1]';
if Sh==0
  pr=[ur(1)-1 0]';
  qr =[0 1]' ;
else
  pr=[Sh*(1-ur(1));0];
  qr = [1 1]';
end
end

function u0 = icfun(x)
nx = length(x);
if 0
  % smooth change of C from outer boundary
  C0 = 1./exp(5e6*(1-x));
  u0 =[C0; zeros(1,nx)];
else
  % sharp change of C at outer boundary
  R=1;
  if abs(x-R) < 100*eps
    C0=1;
  else
    C0=0;
  end
  u0 =[C0 0]'*ones(1,nx);
end
end

function plotResultAtTimes(x, y, tspan, tit)
nt = size(y,1);
leg={};
figure;
hold on;
for i=1:nt
  plot(x, y(i,:));
  leg{i} = sprintf('%f', tspan(i));
end
grid
title(tit);
legend(leg);
hold off;
end


