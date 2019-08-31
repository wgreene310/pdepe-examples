function beamFreeVibrationCantileverBC
% NOTE: This example requires a fix to the R2019a version
% of pdepe. Replace this line in pdepe:
% D( c == 0, 2:nx-1) = 0; 
% (line 256 in R2019a)
% with this line:
% D( c == 0, :) = 0;
E=200e9;
thick = .2;
width = .1;
I = width*thick^3/12;
EI=E*I;
A=thick*width;
L=10;
F=0;
N=0;
rho=7700; % density of steel, kg/m^3
betaL=1.875104068711961;
m = 0;
amp=.1;
numElems=20; % even number of elements gives a node in the middle of the beam
numNodes = numElems + 1;
x = linspace(0,L,numNodes);
t=linspace(0, 1.2, 100);
% use anonymous functions to pass parameters to functions required by pdepe
pde = @(x,t,u,DuDx) beampde(x,t,u,DuDx, EI, rho, A, F, N);
bc = @(xl,ul,xr,ur,t) beambc(xl,ul,xr,ur,t);
ic = @(x) beamic(x, L, amp, betaL);
sol = pdepe(m,pde,ic,bc,x,t);

% plot the results
uTip = sol(:, end, 1);
uAnalVibr = analFreeVibr(L, EI, rho*A, amp, x, t, betaL);
figure; plot(t, uTip, t, uAnalVibr(:,end), 'o'); grid;
title('Displacement at beam tip as a function of time');
figure; plot(x,sol(end,:,1), x, uAnalVibr(end, :), 'o'); grid;
title('Displacement at final time');
s1 = sol(:,:,1);
fprintf('Max displacement error=%g\n', max(abs(s1(:)-uAnalVibr(:))));
end

function [cr,fr,sr] = beampde(x,t,u,DuDx, EI, rho, A, F, N)
cr = [1 rho*A 0]';
fr = [0 -EI*DuDx(3)-N*DuDx(1) -DuDx(1)]';
sr = [u(2) F u(3)]';
end

function u0 = beamic(x, L, amp, betaL)
  % Choose initial condition in the shape of the eigenvector 
  % corresponding to the lowest eigenvalue. This makes comparison
  % with the analytical solution easier.
  nx=length(x);
  beta=betaL/L;
  R=(cosh(betaL) + cos(betaL))/(sinh(betaL) + sin(betaL));
  betaX = beta*x;
  w0 = amp*(cosh(betaX) - cos(betaX) - R*(sinh(betaX) - sin(betaX)));
  d2w0dx2=amp*beta^2*(cosh(betaX) + cos(betaX) - R*(sinh(betaX) + sin(betaX)));
  u0 = [w0; zeros(1,nx); d2w0dx2];
end
	
function [pl,ql,pr,qr] = beambc(xl,ul,xr,ur,t)
% clamped end
pl = [ul(1);ul(2);0];
ql = [0;0;1];
% free end
pr = [0;0;ur(3)];
qr = [1;1;0];
end

function u=analFreeVibr(L, EI, rhoA, amp, x, t, betaL)
 omega=(betaL/L)^2*sqrt(EI/rhoA);
 u0 = beamic(x, L, amp, betaL);
 u=cos(omega*t)'*u0(1,:);
 end
  


