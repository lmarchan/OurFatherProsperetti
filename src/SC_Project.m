%% Scientific Computing for Mechanical Engineers Project
% Leon Marchant
% Professors Andrea Prosperetti & Amit Amritkar 
% Project A - Poisson Equation APc1-6
% Write a computer code to solve the two-dimensional Poisson Equationin the
% domain of interest with specified boundary conditions

clc
clear
% Domain of interest is a rectangle with:
% ax<x<bx     ay<y<by
% Where ax=ay=-pi and bx=by=pi
ax=-pi;
bx=pi;
ay=-pi;
by=pi;

% Discretizing the domain into nx and ny increments to obtain segment
% lengths hx and hy
nx=20;
ny=20;

hx = (bx-ax)/(nx-1);
hy = (by-ay)/(ny-1);

% Boundary conditions for u(x,ay) and u(x,by)
%phi(x) = (x-ax)*(x-ax)*sin(pi(x-ax)/(2*(bx-ax)));
%psi(x) = (cos(pi*(x-ax))-1)*cosh(bx-x);

% Boudary conditions for du/dx(x=ax)=0 & du/dx(x=bx)=0 

% Outlining the function F(x,y) 
%F(x,y) = sin(pi*(x-ax)/(bx-ax))*cos(pi*(2*(y-ay)/(by-ay)+1)/2);

% Ghost nodes
axneg=ax-hx;
ayneg=ay-hy;
nxplus=nx+2;
nyplus=ny+2;

% Discretizing u(x,y)
for i = 1:nxplus
    x(i) = axneg+hx*(i-1);
end

for j = 1:nyplus
    y(j) = ayneg+hy*(j-1);
end

