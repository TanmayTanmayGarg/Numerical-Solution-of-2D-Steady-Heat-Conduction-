function x = myGaussElim(A, b)
    % Solve Ax = b using Gaussian Elimination without pivoting
    % Input:  A - coefficient matrix (n x n),b - right-hand side vector (n x 1)
    
    [n, m] = size(A);
    if n ~= m
        error('Matrix A must be square.');
    end
    if length(b) ~= n
        error('Dimension mismatch between A and b.');
    end

    % Form augmented matrix
    Ab = [A b];

    % Forward elimination
    for k = 1:n-1
        % No pivoting
        for i = k+1:n
            factor = Ab(i, k) / Ab(k, k);
            Ab(i, k:end) = Ab(i, k:end) - factor * Ab(k, k:end);
        end
    end

    % Back substitution
    x = zeros(n, 1);
    for i = n:-1:1
        x(i) = (Ab(i, end) - Ab(i, i+1:n) * x(i+1:n)) / Ab(i, i);
    end
end

