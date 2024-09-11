clear; % clear Workspace
clc;   % clear Command Window
close all;

% Parameters
beta = 0.9; % discount factor
A = [0.8, 1.2];
alpha = 0.3;

% Calculation of Steady State
% Euler equation gives 1 = beta(0.3Kss^{-0.7}+0.3)
Kss = ((1 / beta / 0.3) - 1) ^ (-1 / 0.7);

% 500 grids around steady state
Kgrid = linspace(0, 2 * Kss, 500);
N = 500 * 2; % number of states
% Initializing value function
V = zeros(500, 2);
% First input is K and second input is z
% Policy function
Copt = zeros(500, 2);

% Create matrices of consumption and instantaneous utility possibilities
k = Kgrid;
c = (k.^alpha)' * ones(size(k)) - ones(size(k))' * k; % kpoints x kpoints consumption matrix
ctiny = 0.0001; % a very small number
c(c == 0) = ctiny; % set all zero elements of c to a very small number
c(c < 0) = NaN; % set all negative elements of c to "not a number"
u = log(c); % kpoints x kpoints instantaneous utility matrix

%% Iteration
iter = 0;
norm = 1;
oldV = V;
newV = V; % They are used to compute norm
min_consumption = 0; % Define the minimum consumption level

while norm > 1e-6
    iter = iter + 1;
    for i = 1:N
        % Find the current state
        if i <= 500
            kindex = i;
            z = 1;
            kt = Kgrid(kindex);
            At = A(z);
        else
            kindex = i - 500;
            z = 2;
            kt = Kgrid(kindex);
            At = A(z);
        end

        % Find the new value function
        aux = u(:, kindex) + beta * (0.5 * oldV(:, 1) + 0.5 * oldV(:, 2));
        [newV(kindex, z), I] = max(aux);

        % Calculate the consumption based on the optimal policy
        C = At * kt^(0.3) + 0.3 * kt - Kgrid(I);
        
        % Apply minimum consumption level
        Copt(kindex, z) = max(C, min_consumption);
    end

    % Off-grid interpolation for both value function and policy function
    off_grid_indices = find(Kgrid < Kss | Kgrid > Kss);
    for i = 1:length(off_grid_indices)
        index = off_grid_indices(i);
        if index <= 500 % for z = 1
            newV(index, 1) = interp1(Kgrid, newV(:, 1), Kgrid(index), 'linear', 'extrap');
            Copt(index, 1) = interp1(Kgrid, Copt(:, 1), Kgrid(index), 'linear', 'extrap');
        else % for z = 2
            newV(index - 500, 2) = interp1(Kgrid, newV(:, 2), Kgrid(index), 'linear', 'extrap');
            Copt(index - 500, 2) = interp1(Kgrid, Copt(:, 2), Kgrid(index), 'linear', 'extrap');
        end
    end

    % Update the value function using the interpolated values
    norm = max(abs(newV - oldV), [], 'all');
    oldV = newV;
    
    % Update the index I using the interpolated values
    for i = 1:N
        if i <= 500
            z = 1;
        else
            z = 2;
        end
        kt = Kgrid(mod(i - 1, 500) + 1);
        [~, I] = max(u(:, mod(i - 1, 500) + 1) + beta * (0.5 * interp1(Kgrid, newV(:, 1), kt, 'linear', 'extrap') + 0.5 * interp1(Kgrid, newV(:, 2), kt, 'linear', 'extrap')));
        if z == 1
            I1(i) = I;
        else
            I2(i) = I;
        end
    end
    disp([norm, iter]);
end
