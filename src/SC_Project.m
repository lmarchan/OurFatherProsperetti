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

% Boudary conditions for du/dx(x=ax)=0 & du/dx(x=bx)=0 

% Outlining the function F(x,y) 
%F(x,y) = sin(pi*(x-ax)/(bx-ax))*cos(pi*(2*(y-ay)/(by-ay)+1)/2);

% Ghost nodes
axneg=ax-hx;
ayneg=ay-hy;
nxplus=nx+2;
nyplus=ny+2;

% Discretizing x and y, and boundary conditions for u(x,ay),
% u(x,by),(du/dx)@ax = 0 & (du/dx)@bx = 0

%preallocation
x=zeros(1,nxplus);
phi=zeros(1,nxplus);
psi=zeros(1,nxplus);
y=zeros(1,nyplus);
rho=zeros(1,nyplus);
theta=zeros(1,nyplus);

for j = 1:nxplus
    x(j) = axneg+hx*(j-1);
    phi(j) = (x(j)-ax)*(x(j)-ax)*sin(pi*(x(j)-ax)/(2*(bx-ax)));
    psi(j) = (cos(pi*(x(j)-ax))-1)*cosh(bx-x(j));
end

for k = 1:nyplus
    y(k) = ayneg+hy*(k-1);
    rho(k) = phi(2);
    theta(k) = psi(nyplus-1);
end

    


