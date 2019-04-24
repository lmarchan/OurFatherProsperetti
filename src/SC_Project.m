%% Scientific Computing for Mechanical Engineers Project
% Leon Marchant
% Professors Andrea Prosperetti & Amit Amritkar 
% Project A - Poisson Equation APc1-6
% Write a computer code to solve the two-dimensional Poisson Equationin the
% domain of interest with specified boundary conditions

% Domain of interest is a rectangle with:
% ax<x<bx     ay<y<by
% Where ax=ay=-pi and bx=by=pi
ax=-pi;
bx=pi;
ay=-pi;
by=pi;

% Discretizing the domain into nx and ny increments to obtain segmet
% lengths hx and hy
hx = (bx-ax)/(nx+1);
hy = (by-ay)/(ny+1);

% something
phi(x) = (x-ax)*(x-ax)*sin(pi(x-ax)/(2*(bx-ax)));
psi(x) = (cos(pi*(x-ax))-1)*cosh(bx-x);
F(x,y) = sin(pi*(x-ax)/(bx-ax))*cos(pi*(2*(y-ay)/(by-ay)+1)/2);