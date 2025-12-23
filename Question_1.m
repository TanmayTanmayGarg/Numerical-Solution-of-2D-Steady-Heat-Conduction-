clear;      
clc;        
close all;  

N_qus_1 = [11,21,41]; % grid is n*n
run_time_qus_1 = [0,0,0];

%gaussian elimination medthod;

T_up= 423;
T_bottom = 323;
T_right = 473;
h = 250; T_infinity = 573; K = 5;
L = 80e-3; 

for(a = 1:3)
    tic
    n = N_qus_1(a); % grid is n*m
    unknown_grid_point= (n-2)*(n-1);
    m = unknown_grid_point;
    temp_coeff_matrix = zeros(m,m);
    delta_x = L/(n-1);
    c_1 = -4-(2*delta_x*h)/K;
    c_2 = -(2*delta_x*h*T_infinity)/K;
    
    
    %set main diagonal
    i = 1;
    while(i <= m)
        temp_coeff_matrix(i,i) = c_1;
        i = i+1;
        for(j = 1:n-2)
            temp_coeff_matrix(i,i)=-4;
            i = i+1;
        end
    end
    
    %set up right diagonal to the main diagonal
    i=1; j=2;
    while(i <= m && j <= m)
        temp_coeff_matrix(i,j) = 2;
        i = i+1; j = j+1;
        for(k = 1:n-3)
            if(i <= m && j <= m)
                temp_coeff_matrix(i,j)= 1;
                i = i+1; j = j+1;
            end
        end
        if(i <= m && j <= m)
            temp_coeff_matrix(i,j) = 0;
            i=i+1;j=j+1;
        end
    end
    
    %set up left diagonal to the main diagonal
    i=2; j=1;
    while(i <= m && j <= m)
        for(k = 1:n-2)
            if(i <= m && j <= m)
                temp_coeff_matrix(i,j)= 1;
                i = i+1; j = j+1;
            end
        end
        if(i <= m && j <= m)
            temp_coeff_matrix(i,j) = 0;
            i=i+1;j=j+1;
        end
    end
    
    %set up most left side diagonal 
    
    i=n; j=1;
    while(i <= m && j <= m)
        temp_coeff_matrix(i,j) = 1;
        i=i+1;j=j+1;
    end
    
    %set up most right side diagonal 
    i=1; j=n;
    while(i <= m && j <= m)
        temp_coeff_matrix(i,j) = 1;
        i=i+1;j=j+1;
        
    end
    
    % AX = Y;
    y = zeros(m,1);
    y_index = 1;
    for(i = 2:n-1)
        for(j=1:n-1)
            %upper BC
            if(i==2)
                if(j == 1), y(y_index,1) = c_2-T_up;
                elseif(j==n-1), y(y_index)=-T_right-T_up;
                else, y(y_index,1)=-T_up;
                end
    
            elseif( i==n-1)
                if(j==1), y(y_index,1) = c_2-T_bottom;
                elseif(j==n-1), y(y_index)=-T_bottom-T_right;
                else, y(y_index)=-T_bottom;
                end
            
            else
                if(j==1), y(y_index)=c_2;
                elseif (j==n-1), y(y_index)=-T_right;
                end
            end
            y_index = y_index+1;
        end
    end
    
    x=myGaussElim(temp_coeff_matrix,y);
    
    
    %assigning 
    T_field = zeros(n,n);
    for(i = 1:n) %applying BCS;
        T_field(1,i)=T_up;
        T_field(n,i)=T_bottom;
        T_field(i,n)=T_right;
    end
    
    x_index = 1;
    for(i = 2:n-1)
        for(j=1:n-1)
            T_field(i,j) = round(x(x_index,1), 2);
            x_index = x_index+1;
        end
    end
    solver_time = toc;
    run_time_qus_1(a) = solver_time;
    % to plot the countour I have taken the help of AI of which code is written
    % below
    % Define the physical coordinates for the x and y axes
    x_coords = linspace(0, L, n); % Creates n points from 0 to L
    y_coords = linspace(L, 0, n); % Creates n points from L down to 0 for correct mapping with the T_matrix as in T_matrix (x,1) is the point at x,L in real.
    
    % Create a 2D grid of coordinates
    [X, Y] = meshgrid(x_coords, y_coords);
    
    % Create a filled contour plot
    %figure; % Creates a new figure window
    subplot(1,3,a);  
    contourf(X*1000, Y*1000, T_field, 25, 'LineStyle', 'none'); % 20 contour levels
    
    % Add labels and a color bar for clarity
    title("Grid Size n = " + n + ", m = " + n);
    xlabel('X-axis (mm)');
    ylabel('Y-axis (mm)');
    c = colorbar;
    c.Label.String = sprintf('Temperature (^{\\circ}C)');

    axis equal; % Ensures the plot is square, like your body
    
    
end
sgtitle('2D Temperature Distribution for Different Grid Sizes (Gaussian Elimination Method)'); % common figure title for subplots


figure;
total_grid_points = N_qus_1.*N_qus_1;
plot(total_grid_points, run_time_qus_1,'-o','LineWidth',1.5);
title("CPU run time vs total number of grid points for (Gaussian Elimination)");
xlabel('Number of Grid Points');
ylabel('Run Time in Seconds');
grid on;

save('qus1_data.mat', 'N_qus_1','run_time_qus_1');
