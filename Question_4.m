clear;      
clc;        
close all;  
n = 21; % grid is n*m
m = 41;
tic
res2141_ADI = [];  % Column Sweep
runTime21x41_ADI = 0;
%gauss seidel medthod;

%given
T_up= 423;
T_bottom = 323;
T_right = 473;
h = 250; T_infinity = 573; K = 5;
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
col_1_coeff_mat= make_col_1_coeff_mat(T_field,c_1,p_2);
col_i_coeff_mat= make_col_i_coeff_mat(T_field,p_1,p_2);


row_sweep = 1;

while(residual>tolerance)
    T_old = T_field;
    if(row_sweep==1)
        for(i=2:n-1)
            x = solve_row(T_field,i,c_2,row_coeff_mat,p_2);
            T_field(i,1:m-1) = x(:)'; % x(:)' represent transpose.
        end
        row_sweep=0;
    else
        for(j=1:m-1)
            x = solve_col(T_field,j,c_2,col_1_coeff_mat,col_i_coeff_mat,p_2);
            T_field(2:n-1,j) = x;
        end
        row_sweep = 1;
    end
    residual = norm(T_old-T_field);
    res2141_ADI(end+1) = residual;
end
runTime21x41_ADI = toc;

%%
% to plot the countour I have taken the help of AI of which code is written
% below
% Define the physical coordinates for the x and y axes
x_coords = linspace(0, L, m); % Creates n points from 0 to L
y_coords = linspace(L, 0, n); % Creates n points from L down to 0 for correct mapping with the T_matrix as in T_matrix (x,1) is the point at x,L in real.

% Create a 2D grid of coordinates
[X, Y] = meshgrid(x_coords, y_coords);

% Create a filled contour plot
figure; % Creates a new figure window
contourf(X*1000, Y*1000, T_field, 20, 'LineStyle', 'none'); % 20 contour levels

% Add labels and a color bar for clarity
title('2D Temperature Distribution(ADI),Grid 21x41');
xlabel('X-axis (mm)');
ylabel('Y-axis (mm)');
colorbar; % Shows the temperature scale
axis equal; % Ensures the plot is square, like your body
%%

load('qus2_data.mat');
load('qus3_data.mat');
load('col_Sweep_data.mat')

% extract the two residual vectors (41x41 corresponds to index 3)

res3_col_sweep = vec_resi_col_sweep{3};   % Column Sweep
runTime3_col_sweep = run_time_col_sweep(3);

res3_qus_2   = vec_resi_qus_2{3};     % Gauss-Seidel
runTime41_qus_2 = run_time_qus_2(3);

res3_qus_3   = vec_resi_qus_3{3};     % Row Sweep 
runTime41_qus_3 = run_time_qus_3(3);

% plot on loglog (x,y both log scale)
figure;
loglog(1:length(res3_col_sweep), res3_col_sweep, '-o', 'LineWidth', 1, 'MarkerSize',6); hold on;
loglog(1:length(res3_qus_2),   res3_qus_2,   '-s', 'LineWidth', 1, 'MarkerSize',6);
loglog(1:length(res3_qus_3),   res3_qus_3,   '-x', 'LineWidth', 1, 'MarkerSize',1);
loglog(1:length(res2141_ADI),   res2141_ADI,'-d','LineWidth', 1, 'MarkerSize',0.1);
hold off;

legend('Line-by-line Column Sweep', ...
    'Gauss-Seidel', ...
    'Line-by-line Row Sweep' ...
    ,'Alternating Direction Implicit (ADI) method' ...
    ,'Location','southwest');
xlabel('Iteration (log scale)');
ylabel('Residual (log scale)');
title('Residual vs Iteration for grid 41x41');
grid on;
%%

% Extract runtimes for the 21x41 grid
runTime21x41_col = run_time_col_sweep(3);   % Column Sweep
runTime21x41_gs  = run_time_qus_2(3);       % Gauss-Seidel
runTime21x41_row = run_time_qus_3(3);       % Row Sweep

methods = {'Column Sweep','Gauss-Seidel','Row Sweep','ADI'};
runtimes = [runTime21x41_col, runTime21x41_gs, runTime21x41_row, runTime21x41_ADI];

% Plot as bar chart
figure;
bar(runtimes);
set(gca,'XTickLabel',methods,'XTickLabelRotation',20);
ylabel('Runtime (s)');
title('Runtime comparison for grid 21x41');
grid on;


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

%column sweep
function [mat] = make_col_i_coeff_mat(T_field,p_1,p_2)
    n=size(T_field,1);
    mat = zeros(n-2,n-2);
    mat(1,1) = -4*p_1; mat(1,2)=p_2;
    left = 1; middle = 2; right =3; 
    for(i=2:n-2)
        if(middle<n-2)
            mat(i,left)=p_2;
            mat(i,middle)=-4*p_1;
            mat(i,right)=p_2;
        else
            mat(i,left)=p_2;
            mat(i,middle)=-4*p_1;
        end
        left = left+1; middle = middle+1; right = right+1;
    end
end
function [mat] = make_col_1_coeff_mat(T_field,c_1,p_2)
    n=size(T_field,1);
    mat = zeros(n-2,n-2);
    mat(1,1) = c_1; mat(1,2)=p_2;
    left = 1; middle = 2; right =3; 
    for(i=2:n-2)
        if(middle<n-2)
            mat(i,left)=p_2;
            mat(i,middle)=c_1;
            mat(i,right)=p_2;
        else
            mat(i,left)=p_2;
            mat(i,middle)=c_1;
        end
        left = left+1; middle = middle+1; right = right+1;
    end
end

function [x] = solve_col(T_field,j,c_2,col_1_coeff_mat,col_i_coeff_mat,p_2)
    b = col_sweep_col_b(T_field,j,c_2,p_2);
    n=size(T_field,1);
    x = zeros(n-2,1);
    if(j==1)
        x = TDMA_solver(col_1_coeff_mat,b);
    else
        x = TDMA_solver(col_i_coeff_mat,b);
    end
end

function[b]=col_sweep_col_b(T_field,j,c_2, p_2)
    n=size(T_field,1);
    b = zeros(n-2,1);
    if(j==1)
        for(i=2:n-1)
            b(i-1,1)=c_2 - 2*T_field(i,j+1);                
        end
        b(1,1)= b(1,1) -p_2*T_field(1,1);
        b(n-2,1) = b(n-2,1) - p_2*T_field(n,1);
                                           
    else
        for(i=2:n-1)
            b(i-1,1)= -T_field(i,j+1) -T_field(i,j-1);                
        end
        b(1,1)= b(1,1) - p_2*T_field(1,j);
        b(n-2,1) = b(n-2,1) -p_2*T_field(n,j);                              
    end
end