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

% Ghost nodes
axneg=ax-hx;
% ayneg=ay-hy; This makes ghost nodes on bottom boundary (not needed)
nxplus=nx+2;
nyplus=ny+2; %This makes ghost nodes on top boundary (not needed)

%preallocation
x=zeros(1,nxplus);
phi=zeros(1,nxplus);
psi=zeros(1,nxplus);
y=zeros(1,ny);               %       y=zeros(1,nyplus);    (old)
                             %       rho=zeros(1,nyplus);
                             %       theta=zeros(1,nyplus);

% Discretizing x, and boundary conditions for u(x,ay) and u(x,by),


for k = 1:nxplus
    x(k) = axneg+hx*(k-1);
    phi(k) = (x(k)-ax)*(x(k)-ax)*sin(pi*(x(k)-ax)/(2*(bx-ax))); %Boundary Conditions @y=-pi
    psi(k) = (cos(pi*(x(k)-ax))-1)*cosh(bx-x(k));   %Boundary Conditions @y=pi
end
% Discretizing y
% u(x,by),(du/dx)@ax = 0 & (du/dx)@bx = 0
for j = 1:ny                %      k = 1:nyplus    (old)
    y(j) = ay+hy*(j-1);     %       y(k) = ayneg+hy*(k-1);   (old)   
                            %     rho(k) = phi(2);      (incorrect)
                            %     theta(k) = psi(nyplus-1);  (incorrect)
end

% Discretizing the function F(x,y) **
F=zeros(ny,nxplus);
u=zeros(ny,nxplus);
for j = 1:ny
    for k = 1:nxplus       
        F(j,k) = sin(pi*(x(k)-ax)/(bx-ax)).*cos(pi*(2*(y(j)-ay)/(by-ay)+1)/2);
    end
end
u(1,:) = psi;
u(ny,:) = phi;
u(:,1) = u(:,3);
u(:,nxplus) = u(:,nxplus-2);
for i = 1:100
    for j = 2:ny-1 
        u(:,1) = u(:,3);
        for k = 2:nxplus-1       
            u(j,k) = u(j,k-1)+u(j,k+1)+u(j-1,k)+u(j+1,k)-4*u(j,k)+F(j,k)*hy*hx;
        end
    end
end

