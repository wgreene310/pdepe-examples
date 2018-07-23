function beamFreeVibration
E=200e9;
thick = .2;
width = .1;
I = width*thick^3/12;
EI=E*I;
A=thick*width;
L=10;
xMidPoint = L/2;
rho=7700; % density of steel, kg/m^3
m = 0;
amp=.1;
numElems=20; % need a node in the middle of the beam
elemLength = L/numElems;
numNodes = numElems + 1;

x = linspace(0,L,numNodes);
t=linspace(0, .75, 100);
%figure; plot(x, wa); grid;
pde = @(x,t,u,DuDx) beampde(x,t,u,DuDx, EI, rho, A);
bc = @(xl,ul,xr,ur,t) beambc(xl,ul,xr,ur,t);
ic = @(x) beamic(x, L, amp, EI);
tic
if(1)
options.vectorized = 'on';
sol = pde1d(m,pde,ic,bc,x,t, options);
else
sol = pdepe(m,pde,ic,bc,x,t);
end
toc
% plot the results

midNode = int32(numElems/2 + 1);
uMidPt = sol(:, midNode, 1);
max(uMidPt)
uAnalVibr = analFreeVibr(L, EI, rho*A, amp, x, t);
figure; plot(t, uMidPt, t, uAnalVibr(:,midNode), 'o'); grid;
figure; plot(x,sol(end,:,1), x, uAnalVibr(end, :), 'o'); grid;
s1 = sol(:,:,1);
max(abs(s1(:)-uAnalVibr(:)))
%uAnalVibr-s1
assertElementsAlmostEqual(uAnalVibr, s1, 'absolute', .05);
end

function [cr,fr,sr] = beampde(x,t,u,DuDx, EI, rho, A)
nx = length(x);
o=ones(1,nx);
z=zeros(1,nx);
cr = [o ; rho*A*o; z];
fr = [z; -DuDx(3,:); -EI*DuDx(1,:)];
f=1;
sr = [u(2,:); z; u(3,:)];
end

function [cr,fr,sr] = beampdeX(x,t,u,DuDx, EI, rho, A)
nx = length(x);
o=ones(1,nx);
z=zeros(1,nx);
cr = [o ; rho*A*o; z];
fr = [z; -EI*DuDx(3,:); -DuDx(1,:)];
f=1;
sr = [u(2,:); z; u(3,:)];
end

function u0 = beamic(x, L, amp, EI)
  % half sin wave initial condition
  s = sin(pi*x/L);
  u0 = [amp*s; 0; -amp*EI*(pi/L)^2*s];
end
	
function [pl,ql,pr,qr] = beambc(xl,ul,xr,ur,t)
  pl = [ul(1);ul(2);ul(3)];
  ql = [0;0;0];
  pr = [ur(1);ur(2);ur(3)];
  qr = [0;0;0];
end

function u=analFreeVibr(L, EI, rhoA, amp, x, t)
 omega=(pi/L)^2*sqrt(EI/rhoA);
 u=amp*cos(omega*t)'*sin(pi*x/L);
 end
  


