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
% axneg=ax-hx;
% nxplus=nx+2;
% nyplus=ny+2; 

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

% Preallocation
F=zeros(ny,nx);
u=zeros(ny,nx);

% Setting boundary conditions in place
u(1,:) = psi;
u(ny,:) = phi;

tic
for i = 1:100

    for j = 2:ny-1 
%           u(:,1) = u(:,3);
%           u(:,nxplus) = u(:,nxplus-2);
%          F(j,2) = sin(pi*(x(2)-ax)/(bx-ax)).*cos(pi*(2*(y(j)-ay)/(by-ay)+1)/2);
%          
%          
%          F(j,nxplus-1) = sin(pi*(x(nxplus-1)-ax)/(bx-ax)).*cos(pi*(2*(y(j)-ay)/(by-ay)+1)/2);
        for k = 1:nx      
            if k==1
                u(j,k) = (2*u(j,k+1)+u(j-1,k)+u(j+1,k)+F(j,k)*hy*hx)/4;
            elseif k==nx
                u(j,k) = (2*u(j,k-1)+u(j-1,k)+u(j+1,k)+F(j,k)*hy*hx)/4;
                
            else
                    
              
            F(j,k) = sin(pi*(x(k)-ax)/(bx-ax)).*cos(pi*(2*(y(j)-ay)/(by-ay)+1)/2);
            u(j,k) = (u(j,k-1)+u(j,k+1)+u(j-1,k)+u(j+1,k)+F(j,k)*hy*hx)/4;
            end
        end
             
        
    end
    i
    graph = surf(x,y,u);
    drawnow;
    
    refreshdata(graph)
end
    
toc
