
x = linspace(0,1,30);
t = linspace(0,1,10);

% define s coefficient in PDE
A = @(x) cos(3*pi*x);
% define initial condition
B = @(x) sin(pi*x);

% calculate analytical solution
ua=uAnal(x,t,A, B);

pdeFunc = @(x,t,u,DuDx) pde(x,t,u,DuDx,A);
icFunc = @(x) ic(x,B);
bcFunc = @(xl,ul,xr,ur,t) bc(xl,ul,xr,ur,t);
m=0;
u = pdepe(m, pdeFunc,icFunc,bcFunc,x,t);

err=ua(:)-u(:);
fprintf('Solution error: max=%12.3e, average=%12.3e\n', ...
  max(abs(err)), norm(err)/length(err));
figure; plot(t, u(:,end), t, ua(:,end), 'o'); grid on;
xlabel('Time'); ylabel('u'); 
legend('pdepe', 'analytical', 'Location','northwest' );
title('Solution at right end as a function of time');

figure; plot(t, u(:,1), t, ua(:,1), 'o'); grid on;
xlabel('Time'); ylabel('u'); 
legend('pdepe', 'analytical', 'Location','northwest' );
title('Solution at left end as a function of time');

figure; plot(x, u(end,:), x, ua(end,:), 'o'); grid on;
xlabel('x'); ylabel('u'); 
legend('pdepe', 'analytical', 'Location','northwest' );
title('Solution along the length at final time');

function [c,f,s] = pde(x,t,u,DuDx,A)
nx=length(x);
c = ones(1,nx);
f = zeros(1,nx);
s = A(x);
end

function u0 = ic(x,B)
u0 = B(x);
end

function [pl,ql,pr,qr] = bc(xl,ul,xr,ur,t)
pl = 0;
ql = 1;
pr = 0;
qr = 1;
end

function u=uAnal(x,t,A, B)
nt = length(t);
u=t'*A(x) + repmat(B(x),nt,1);
end

%% Copyright (C) 2018 William H. Greene
