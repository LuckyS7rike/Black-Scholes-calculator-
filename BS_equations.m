clear all; close all;

% Given table values

sigma = 0.895; % Volatility
r = 0.0403; % 
T = 0.27397260274 % 
K = 20; % 
S0 = 30.94; % S^0

% Parameters for the grid
M = 50; % Time steps
N = 50; % Steps for asset price
dt = T/M; % Size of the time step
dS = 2*S0/N; % Size of the asset price step
S = 0:dS:2*S0; % Asset price values
t = 0:dt:T; % Time 

% Create the option price
option = zeros(M+1, N+1);

% Makeing sure this is a European put option
%option(M+1,:) = max(K-S, 0); % Set as put option
option(M+1,:) = max(S-K, 0); % Set as call option 

% Finite difference (Crank Nicolson) 
for i = M:-1:1
    % Setting up the matrix to solve Crank Nicolson method
    M = zeros(N-1, N-1);
    B = zeros(N-1, 1);
    
    for j = 2:N
        % Coefficients
        c1 = 0.25*dt*(sigma^2*(j-1)^2 - r*(j-1));
        c2 = -dt*0.5*(sigma^2*(j-1)^2 + r);
        c3 = 0.25*dt*(sigma^2*(j-1)^2 + r*(j-1));
        
        % Fill the matrix M
        if j > 2
            M(j-1,j-2) = -c1; % Lower diagonal of the matrix
        end
        M(j-1,j-1) = 1 - c2; % Main diagonal of the matrix
        if j < N
            M(j-1,j) = -c3; % Upper diagonal of the matrix
        end
        
        % Sub for the vector B 
        B(j-1) = c1*option(i+1,j-1) + (1+c2)*option(i+1,j) + c3 *option(i+1,j+1);
    end
    
    % Compute using LU factorization function
    [L, U] = lu(M);
    
    % Solving the system
    y = L \ B;
    option(i,2:N) = U \ y;
    
    % Boundary Conditions
    option(i,1) = 0; 
    option(i,N+1) = 2*option(i,N) - option(i,N-1); 
end

% Fair put option price
option_price = option(1, S0/dS + 1)

%%
[Call, Put] = blsprice(S0, K, r, T, sigma)