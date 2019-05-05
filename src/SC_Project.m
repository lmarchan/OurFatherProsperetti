%% Scientific Computing for Mechanical Engineers Project
% Leon Marchant
% Professors Andrea Prosperetti & Amit Amritkar 
% Project A - Poisson Equation APc1-6
% Write a computer code to solve the two-dimensional Poisson Equationin the
% domain of interest with specified boundary conditions

clc
clear
tic         %Start timer
% Domain of interest is a rectangle with:
% ax<x<bx     ay<y<by
% Where ax=ay=-pi and bx=by=pi
ax=-pi;
bx=pi;
ay=-pi;
by=pi;

% Discretizing the domain into nx and ny increments to obtain segment
% lengths hx and hy
nx=100;
ny=100;
hx = (bx-ax)/(nx-1);
hy = (by-ay)/(ny-1);

% Preallocation
x=zeros(1,nx);
phi=zeros(1,nx);
psi=zeros(1,nx);
y=zeros(1,ny);               

% Discretizing x, and boundary conditions for u(x,ay) and u(x,by),
for k = 1:nx
    x(k) = ax+hx*(k-1);
    phi(k) = (x(k)-ax)*(x(k)-ax)*sin(pi*(x(k)-ax)/(2*(bx-ax))); %Boundary Conditions @y=-pi
    psi(k) = (cos(pi*(x(k)-ax))-1)*cosh(bx-x(k));   %Boundary Conditions @y=pi
end
% Discretizing y
% u(x,by),(du/dx)@ax = 0 & (du/dx)@bx = 0
for j = 1:ny               
    y(j) = ay+hy*(j-1);    
end
% Preallocation for F and u
F=zeros(ny,nx);
u=zeros(ny,nx);

% Setting boundary conditions in place for upper and lower boundaries
u(1,:) = psi;
u(ny,:) = phi;
E=1;        %Setting initial error 

while E > 10^-7   % Code to remaine running until error is less than set amount
    GSu = u;        % Saving prior u to use for Gauss-Seidel method
    for k = 1:nx 
        for j = 2:ny-1      
            if k==1   % Setting boundary conditions in place for left boundary
                u(j,k) = (2*u(j,k+1)+u(j-1,k)+GSu(j+1,k)+F(j,k)*hy*hx)/4;% Using previous u for Gauss-Seidel method
            elseif k==nx   % Setting boundary conditions in place for right boundary
                u(j,k) = (2*u(j,k-1)+u(j-1,k)+GSu(j+1,k)+F(j,k)*hy*hx)/4;% Using previous u for Gauss-Seidel method
            else
            F(j,k) = sin(pi*(x(k)-ax)/(bx-ax)).*cos(pi*(2*(y(j)-ay)/(by-ay)+1)/2);
            u(j,k) = (u(j,k-1)+GSu(j,k+1)+u(j-1,k)+GSu(j+1,k)+F(j,k)*hy*hx)/4;% Using previous u for Gauss-Seidel method
            end
        end 
    end
    E = max(max(abs(GSu-u))); % Using the L infinite error equation to find the difference between iterations, to decide whether to keep iterating
end
graph = surf(x,y,u); % Graph the surface plot for x, y, and the function u(x,y)
xlabel('x','Fontsize',16); % Label axes and make fonts larger to improve readbility 
ylabel('y','Fontsize',16);
zlabel('u(x,y)','Fontsize',16);
title('Solution to the Poisson Equation Using the Gauss-Seidel Method','Fontsize',16);
colorbar('vertical') % Use color bar to better visualize the value of u(x,y)
colormap('cool') % Use cool colormap to avoid missrepresenting data to people whom are color blind
toc 
