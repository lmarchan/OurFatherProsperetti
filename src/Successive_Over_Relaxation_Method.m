%% Scientific Computing for Mechanical Engineers Project
% Leon Marchant
% Professors Andrea Prosperetti & Amit Amritkar 
% Project A - Poisson Equation APc1-6
% Write a computer code to solve the two-dimensional Poisson Equationin the
% domain of interest with specified boundary conditions
%% Successive Over Relaxation Method 
clc
clear
tic         %Start timer
checkpoint='checkpoint.mat';
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
F=zeros(ny,nx);
u=zeros(ny,nx);

% Discretizing x, and boundary conditions for u(x,ay) and u(x,by),
for k = 1:nx
    x(k) = ax+hx*(k-1);
    u(ny,k) = (x(k)-ax)*(x(k)-ax)*sin(pi*(x(k)-ax)/(2*(bx-ax)));   %Boundary Conditions @y=-pi
    u(1,k) = (cos(pi*(x(k)-ax))-1)*cosh(bx-x(k));   %Boundary Conditions @y=pi
end
% Discretizing y
% u(x,by),(du/dx)@ax = 0 & (du/dx)@bx = 0
for j = 1:ny               
    y(j) = ay+hy*(j-1);    
end
% Set up the initial SRO parameter w
wopt = 2/(1+sin(pi*hx)); % optimal parameter w found using equation obtained from: "https://userpages.umbc.edu/~gobbert/papers/YangGobbertAML2007.pdf"
%Setting initial error 
E=1;
%Setting wanted magnitude of error
Ewanted=10^-7; 
% Initialize iteration and checkpointing values
i=0;
t=0;
while E > Ewanted   % Code to remaine running until error is less than wanted amount
    SORu = u;        % Saving prior u to use for Gauss-Seidel method
        
    for k = 1:nx 
        for j = 2:ny-1      
            F(j,k) = sin(pi*(x(k)-ax)/(bx-ax)).*cos(pi*(2*(y(j)-ay)/(by-ay)+1)/2); % Find value for F to be later plugged into u(x,y) equations
            if k==1   % Setting boundary conditions in place for left boundary
                u(j,k) = (2*u(j,k+1)+u(j-1,k)+SORu(j+1,k)+F(j,k)*hy*hx)/4;% Using previous u for Gauss-Seidel method
            elseif k==nx   % Setting boundary conditions in place for right boundary
                u(j,k) = (2*u(j,k-1)+u(j-1,k)+SORu(j+1,k)+F(j,k)*hy*hx)/4;% Using previous u for Gauss-Seidel method
            else
            u(j,k) = (1-wopt)*SORu(j,k)+wopt*(u(j,k-1)+SORu(j,k+1)+u(j-1,k)+SORu(j+1,k)+F(j,k)*hy*hx)/4;% Using previous u for Gauss-Seidel method
            end
        end 
    end
    t=t+1;
    if t == 1000 % Set to save checkpoint every 1000 iterations
        disp('Checkpointing program, 1000 iterations have occured'); % Let user know that their data is being saved, also lets them know 1000 iterations have occured
        save(checkpoint);
        t=0; % Reset checkpoint ticker
    end
    i=i+1;
    E = max(max(abs(SORu-u))); % Using the L infinite error equation to find the difference between iterations, to decide whether to keep iterating
end
r = max(nx,ny); %This variable will decide whether to include gridlines for the surface plot
if r<100        % At around 100 nodes is when the gridline begin to obscure the color, this if statement 
    graph = surf(x,y,u); % Graph the surface plot for x, y, and the function u(x,y)
    % Label axes and make fonts larger to improve readbility 
    xlabel('x','Fontsize',16); 
    ylabel('y','Fontsize',16);
    zlabel('u(x,y)','Fontsize',16);
    title('Solution to the Poisson Equation Using the Successive Over Relaxation Method','Fontsize',16);
    % Use color bar to better visualize the value of u(x,y)
    colorbar('vertical')
    % Use cool colormap to avoid missrepresenting data to people whom are color blind
    colormap('cool') 
else
    graph = surf(x,y,u); % Graph the surface plot for x, y, and the function u(x,y)
    % Label axes and make fonts larger to improve readbility 
    xlabel('x','Fontsize',16); 
    ylabel('y','Fontsize',16);
    zlabel('u(x,y)','Fontsize',16);
    title('Solution to the Poisson Equation Using the Successive Over Relaxation Method','Fontsize',16);
    % Use color bar to better visualize the value of u(x,y)
    colorbar('vertical')
    % Use cool colormap to avoid missrepresenting data to people whom are color blind
    colormap('cool') 
    set(graph,'edgecolor','none') %This turns off grid lines since they obscure the color when there are many nodes
end 
toc % end timer