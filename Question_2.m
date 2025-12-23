clear;      
clc;        
close all;  


%gauss seidel medthod;
N=[11, 21, 21,41];
M=[11, 21, 41,41];
itr = [0,0,0,0];
run_time_qus_2 = [0,0,0,0];
vec_resi_qus_2 = cell(1,4);   % one cell per grid size


%given
T_up= 423;
T_bottom = 323;
T_right = 473;
h = 250; T_infinity = 573; K = 5;

%
for(a= 1:size(N,2))
    tic
    n = N(a); % grid is n*m
    m = M(a);
    res_hist = [];

    L = 80e-3; delta_x = L/(m-1); delta_y = L/(n-1);
    
    p_1 = (delta_y^2 + delta_x^2)/(2*delta_y^2);
    p_2 = (delta_x^2)/(delta_y^2);
    c_1 = -(4*p_1)-(2*delta_x*h)/K;
    c_2 = -(2*delta_x*h*T_infinity)/K;
    
    %initialize (n)*(n) grid points
    T_field = ones(n,m);
    T_field(1,:)= T_up;
    T_field(n,:)=T_bottom;
    T_field(:,m)=T_right;
    
    tolerance = 1e-6;
    residual  = 10;
    
    no_of_itr = 0;
    while(residual>tolerance)
        no_of_itr = no_of_itr +1;
        T_old = T_field;
        for(i=2:n-1)
            for(j=1:m-1)
                if(j == 1)
                    T_field(i,1) = (c_2- (2*T_field(i,2)) - p_2*T_field(i-1,1) -p_2*T_field(i+1,1) )/c_1;
                else
                    T_field(i,j) =(-p_2*T_field(i-1,j)-p_2*T_field(i+1,j)-T_field(i,j-1)-T_field(i,j+1))/(-4*p_1);
                end 
            end
        end
        residual = norm(T_old-T_field);
        res_hist(end+1) = residual;   
    end
    vec_resi_qus_2{a} = res_hist;
    solver_time = toc;
    itr(a) = no_of_itr;
    run_time_qus_2(a) = solver_time;
    % to plot the countour I have taken the help of AI of which code is written
    % below
    % Define the physical coordinates for the x and y axes
    x_coords = linspace(0, L, m); % Creates n points from 0 to L
    y_coords = linspace(L, 0, n); % Creates n points from L down to 0 for correct mapping with the T_matrix as in T_matrix (x,1) is the point at x,L in real.
    
    % Create a 2D grid of coordinates
    [X, Y] = meshgrid(x_coords, y_coords);
    
    % Create a filled contour plot
    subplot(1,4,a);
    contourf(X*1000, Y*1000, T_field, 20, 'LineStyle', 'none'); % 20 contour levels
    
    % Add labels and a color bar for clarity
    title("Grid Size n = " + n + ", m = " + m);
    xlabel('X-axis (mm)');
    ylabel('Y-axis (mm)');
    c = colorbar;
    c.Label.String = sprintf('Temperature (^{\\circ}C)');
    axis equal; % Ensures the plot is square, like your body
    
    
end

sgtitle('2D Temperature Distribution for Different Grid Sizes (Gauss Seidel Method)'); % common figure title for subplots

%% 

figure;
total_grid_points = N.*M;
plot(total_grid_points, run_time_qus_2,'-o','LineWidth',1.5);
title("CPU run time vs total number of grid points (Gauss Seidel)");
xlabel('Number of Grid Points');
ylabel('Run Time in Seconds');

%% 

figure;
loglog(1:length(vec_resi_qus_2{1}), vec_resi_qus_2{1},'LineWidth',1.3);hold on;
for(a= 2:size(N,2))
    loglog(1:length(vec_resi_qus_2{a}), vec_resi_qus_2{a},'LineWidth',1.3);
end

hold off;
legend('n=11,m=11','n=21,m=21','n=21,m=41','n=41,m=41','Location', 'northeast');
xlabel('Iteration (log scale)');
ylabel('Residual (log scale)');
title('Residual vs Iteration (Gauss-Seidel)');
grid on;

save('qus2_data.mat', 'run_time_qus_2', 'vec_resi_qus_2');
%%

load('qus1_data.mat');
% extract the two residual vectors (41x41 corresponds to index 4)

res41_qus_2   = vec_resi_qus_2{4}; % Gauss-Seidel
runTime41_qus_2 = run_time_qus_2;

runTime41_qus_1 = run_time_qus_1;   % Gaussian Elimination
grid_points_qus_1 = N_qus_1.*N_qus_1;

figure;
total_grid_points = N.*M;
semilogy(total_grid_points, runTime41_qus_2, '-o',...
         grid_points_qus_1, runTime41_qus_1,'-o' ,'LineWidth', 1.5);

title("CPU run time (Log Scale) vs total number of grid points");
legend('Gauss-Seidel','Gaussian Elimination', 'Location','northwest');
xlabel('Number of Grid Points');
ylabel('Run Time in Seconds (Log Scale)');
grid on; % Adding a grid can help with readability on log plots
