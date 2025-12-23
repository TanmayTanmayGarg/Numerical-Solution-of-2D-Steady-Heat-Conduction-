clear;      
clc;        
close all;  

N=[11, 21, 21,41];
M=[11, 21, 41,41];
itr = [0,0,0,0];
run_time_qus_3 = [0,0,0,0];
vec_resi_qus_3 = cell(1,4);    % one cell per grid size

%given
T_up= 423;
T_bottom = 323;
T_right = 473;
h = 250; T_infinity = 573; K = 5;

for(a= 1:size(N,2))
    tic
    n = N(a); % grid is n*m
    m = M(a);
    res_hist = [];

    L = 80e-3; delta_x = L/(m-1); delta_y = L/(n-1);
    p_1 = (delta_y^2 + delta_x^2)/(2*delta_y^2);
    p_2 = delta_x^2/delta_y^2;
    c_1 = -(4*p_1)-(2*delta_x*h)/K;
    c_2 = -(2*delta_x*h*T_infinity)/K;
    
    %initialize (n)*(n) grid points
    T_field = ones(n,m);
    T_field(1,:)= T_up;
    T_field(n,:)=T_bottom;
    T_field(:,m)=T_right;
    
    tolerance = 1e-6;
    residual  = 10;
    row_coeff_mat = make_row_coeff_mat(T_field,c_1,p_1); %pre define one time; because it is same for each row.
    no_of_itr = 0;
    while(residual>tolerance)
        no_of_itr = no_of_itr+1;
        T_old = T_field;
        for(i=2:n-1)
            x = solve_row(T_field,i,c_2,row_coeff_mat,p_2);
            T_field(i,1:m-1) = x(:)';
        end
        residual = norm(T_old-T_field);
        res_hist(end+1) = residual; 
    end
    vec_resi_qus_3{a} = res_hist;
    solver_time = toc;
    itr(a) = no_of_itr;
    run_time_qus_3(a) = solver_time;

    % to plot the countour I have taken the help of AI of which code is written
    % below
    % Define the physical coordinates for the x and y axes
    x_coords = linspace(0, L, m); % Creates n points from 0 to L
    y_coords = linspace(L, 0, n); % Creates n points from L down to 0 for correct mapping with the T_matrix as in T_matrix (x,1) is the point at x,L in real.
    
    % Create a 2D grid of coordinates
    [X, Y] = meshgrid(x_coords, y_coords);
    
    % Create a filled contour plot
    subplot(1,4,a); 
    contourf(X*1000, Y*1000, T_field, 50, 'LineStyle', 'none'); % 20 contour levels
    
    % Add labels and a color bar for clarity
    title("Grid Size n = " + n + ", m = " + m);
    xlabel('X-axis (mm)');
    ylabel('Y-axis (mm)');
    c = colorbar;
    c.Label.String = sprintf('Temperature (^{\\circ}C)');
    axis equal; % Ensures the plot is square, like your body
    
end
sgtitle("2D Temperature Distribution for Different Grid Sizes (line-by-line method ,row sweep)"); % common figure title for subplots
%%
%residual vs number of iteration for different grid size (for only row sweep) 

figure;
loglog(1:length(vec_resi_qus_3{1}), vec_resi_qus_3{1},'LineWidth',1.3);hold on;%both x and y log
for(a= 2:size(N,2))
    loglog(1:length(vec_resi_qus_3{a}), vec_resi_qus_3{a},'LineWidth',1.3);%both x and y log
end

hold off;
legend('n=11,m=11','n=21,m=21','n=21,m=41','n=41,m=41','Location', 'northeast');
xlabel('Iteration(log scale)');
ylabel('Residual (log scale)');
title('Residual vs Iteration (line-by-line method ,row sweep)');
grid on;

%% 
%residual vs number of iteration for 41x41 grid size (compare row sweep and gauss seidel) 


load('qus2_data.mat');
load('qus1_data.mat');
% extract the two residual vectors (41x41 corresponds to index 4)

res41_qus_3 = vec_resi_qus_3{4};   % Row Sweep
runTime41_qus_3 = run_time_qus_3;  

res41_qus_2   = vec_resi_qus_2{4}; % Gauss-Seidel
runTime41_qus_2 = run_time_qus_2;

runTime41_qus_1 = run_time_qus_1;   % Gaussian Elimination
grid_points_qus_1 = N_qus_1.*N_qus_1;


% plot on log log (y log scale, x log scale)
figure;
loglog(1:length(res41_qus_3), res41_qus_3, '-o', 'LineWidth', 1.4, 'MarkerSize',2); hold on;
loglog(1:length(res41_qus_2),   res41_qus_2,   '-s', 'LineWidth', 1.4, 'MarkerSize',2);
hold off;
legend('Line-by-line Row Sweep(41×41)', 'Gauss-Seidel (41×41)', 'Location','northeast');
xlabel('Iteration (log scale)');
ylabel('Residual (log scale)');
title('Residual vs Iteration for grid 41 x 41');
grid on;


%% 
% CPU run times vs grid points for methods:- "Gaussian Elimination", "Gauss Seidel", "Line by Line Row Sweep"  

% Use semilogy to plot the y-axis on a log scale because, gaussian
% elimination take very long time as compared to the gauss seidel and line
% by line approach.

figure;
total_grid_points = N.*M;
semilogy(total_grid_points, runTime41_qus_3, '-o',...
         total_grid_points, runTime41_qus_2, '-o',...
         grid_points_qus_1, runTime41_qus_1,'-o' ,'LineWidth', 1.5);

title("CPU run time (Log Scale) vs total number of grid points");
legend('Line-by-line Row Sweep', 'Gauss-Seidel','Gaussian Elimination', 'Location','northwest');
xlabel('Number of Grid Points');
ylabel('Run Time in Seconds (Log Scale)');
grid on; % Adding a grid can help with readability on log plots
save('qus3_data.mat', 'run_time_qus_3', 'vec_resi_qus_3');
%% 


function [x] = solve_row(T_field,i,c_2,row_coeff_mat,p_2)
    b = row_sweep_col_b(T_field,i,c_2,p_2);
    x = TDMA_solver(row_coeff_mat,b);
end

function [mat] = make_row_coeff_mat(T_field,c_1,p_1)
    m=size(T_field,2);
    mat = zeros(m-1,m-1);
    mat(1,1) = c_1; mat(1,2)=2;
    left = 1; middle = 2; right =3; 
    for(j=2:m-1)
        if(middle<m-1)
            mat(j,left)=1;
            mat(j,middle)=-4*p_1;
            mat(j,right)=1;
        else
            mat(j,left)=1;
            mat(j,middle)=-4*p_1;
        end
        left = left+1; middle = middle+1; right = right+1;
    end
end

function[b]=row_sweep_col_b(T_field,i,c_2,p_2)
    m=size(T_field,2);
    b = zeros(m-1,1);
    for(j = 1:m-1)
        if(j==1)
            b(j,1) = c_2-p_2*T_field(i-1,j)-p_2*T_field(i+1,j);
        elseif(j<m-1)
            b(j,1) = -p_2*T_field(i-1,j)-p_2*T_field(i+1,j);
        else
            b(j,1)=-T_field(i,j+1)-p_2*T_field(i-1,j)-p_2*T_field(i+1,j);
        end
    end
end

function [x] = TDMA_solver(A,b)
    n = size(A,1);
    for(i = 2:n)
        factor = A(i,i-1)/A(i-1,i-1);
        A(i,i) = A(i,i) - A(i-1,i) * factor;
        b(i) = b(i) - b(i-1)*factor;
        A(i,i-1) = 0;
    end
  
    x = zeros(n,1);
    x(n) = b(n)/A(n,n); 
    for(i = n-1:-1:1)
        x(i) =(b(i) - A(i,i+1)*x(i+1) )/A(i,i);
    end
end