n_val = 50;
a_val = -1;
b_val = 2;
c_val = -1;
w_val = .01;
temp = 0;

matrix = create_matrix(a_val, b_val, c_val, n_val);
b_matrix = ones(n_val, 1);
jacobi_err = [];
GS_err = [];
SOR_err = [];

fprintf("Gaussian Elimination: ");
output_val = gauss_elim(matrix, b_matrix);
disp(output_val);

fprintf("Gaussian Elimination (PP): ");
output_val_pp = gauss_elim_pp(matrix, b_matrix);
disp(output_val_pp);

fprintf("Gaussian Elimination (SPP): ");
output_val_spp = gauss_elim_spp(matrix, b_matrix);
disp(output_val_spp);

%------------------ USED FOR OMEGA
%{
omega = 0.01;
min_iter = inf;
min_omega = 0;
omega_vals = [];
iter_vals = [];
while omega <= 2
    [sol, iter] = SOR_method(a_val, b_val, 10^-10, w_val);
    iter_vals(end+1) = iter;
    omega_vals(end+1) = omega;
    if iter < min_iter
        min_omega = omega;
        min_iter = iter;
    end
    omega = omega + 0.01;
end
if temp == 0
    fprintf("Minimum number of iterations: %f W: %f\n", min(iter_vals), min_omega);
    temp = 1;
end
%}
%------------------

output_val_jacobi = jacobi_method(matrix, b_matrix, 10^-10);
disp(output_val_jacobi);

temp = 0;

output_GS = gauss_seidel_method(matrix, b_matrix, 10^-10);
disp(output_GS);

[SOR_method_output, iter] = SOR_method(matrix, b_matrix, 10^-10, w_val);
disp(SOR_method_output);  

% Making the matrix
function [output_matrix] = create_matrix(a, b, c, n)
    output_matrix = zeros(n, n);
    for x = 1:n
        for y = 1:n
            if x == y
                output_matrix(x, y) = b;
            elseif y+1 == x
                output_matrix(x, y) = a;
            elseif x+1 == y
                output_matrix(x, y) = c;
            else
                output_matrix(x, y) = 0;
            end
        end
    end
end

function [gaussElimOutput] = gauss_elim(A,b)
    n = length(b);
    % Initialize a matrix with all zeros
    gaussElimOutput = zeros(n, 1);
    % Apply forward elimination
    for k = 1: n-1
        for i = k+1:n
            % If the diagonal element is zero, set it to a very small value to avoid division by zero
            % Calculate the multiple for the corresponding row
            multiple = A(i,k)/A(k,k);
            % Apply the multiple to the corresponding row and column
            for j = k+1:n
                A(i,j) = A(i,j) - multiple*A(k,j);
            end
            % Apply the multiple to the b matrix
            b(i) = b(i) - multiple*b(k);
        end
    end
    % If the last diagonal element is zero, set it to a very small value to avoid division by zero
    if A(n,n) == 0
        A(n,n) = 0.0000000000000001;
    end
    % Apply backward substitution on the matrix
    gaussElimOutput(n) = b(n)/A(n,n);
    for i = n-1:-1:1
        sum = b(i);
        for j=i+1:n
            sum = sum-A(i,j)*gaussElimOutput(j);
        end
        % If the diagonal element is zero, set it to a very small value to avoid division by zero
        if A(i,i) == 0
            A(i,i) = 0.0000000000000001;
        end
        % Calculate the solution matrix
        gaussElimOutput(i) = sum/A(i,i);
    end
end

function [gaussElimPP] = gauss_elim_pp(A_matrix, b_matrix)
    %Get the length of the b_matrix, which is the number of equations
    n = length(b_matrix);
    %Initialize the gaussElimPP vector with zeros to hold the solution
    gaussElimPP = zeros(n, 1);
    %Perform forward elimination with partial pivoting
    for k = 1:n-1
        %Partial pivoting to avoid division by zero
        [~, index] = max(abs(A_matrix(k:n, k)));
        index = index + k - 1;
        if A_matrix(index, k) == 0
            error('Singular Matrix');
        end
        if index ~= k
            %Swap rows
            A_matrix([k,index],:) = A_matrix([index,k],:);
            b_matrix([k,index]) = b_matrix([index,k]);
        end
        %Gaussian elimination
        for i = k+1:n
            multiple = A_matrix(i,k)/A_matrix(k,k);
            A_matrix(i,k+1:n) = A_matrix(i,k+1:n) - multiple*A_matrix(k,k+1:n);
            b_matrix(i) = b_matrix(i) - multiple*b_matrix(k);
        end
    end
    %Perform backward substitution
    gaussElimPP(n) = b_matrix(n)/A_matrix(n,n);
    for i = n-1:-1:1
        gaussElimPP(i) = (b_matrix(i) - A_matrix(i,i+1:n)*gaussElimPP(i+1:n))/A_matrix(i,i);
    end
end

function [gaussElimOutput] = gauss_elim_spp(A,b)
    n = length(b);
    idx = (1:n)';
    s = max(abs(A'));
    for k = 1:(n-1)
        r = abs(A(idx(k),k)/s(idx(k)));
        p = k;
        for i = (k+1):n
            % Scaled partial pivoting
            t = abs(A(idx(i),k)/s(idx(i)));
            if t > r
                r = t;
                p = i;
            end
        end
        % Swap indices based on pivot selection
        q = idx(p);
        idx(p) = idx(k);
        idx(k) = q;
        for i = (k+1):n
            % Choosing scaled pivot position for computation
            A(idx(i),k) = A(idx(i),k)/A(idx(k),k);
            for j = (k+1):n
                A(idx(i),j) = A(idx(i),j)-A(idx(i),k)*A(idx(k),j);
            end
        end
    end
    % Update corresponding matrix
    y = zeros(n,1);
    y(1) = b(idx(1));
    for i = 2:n
        y(i) = b(idx(i));
        for j = 1:(i-1)
            y(i) = y(i)-A(idx(i),j)*y(j);
        end
    end
    % Backward substitution and solving
    gaussElimOutput = zeros(n,1);
    gaussElimOutput(n) = y(n)/A(idx(n),n);
    for i = (n-1):-1:1
        gaussElimOutput(i) = y(i);
        for j = (i+1):n
            % Update jth position
            gaussElimOutput(i) = gaussElimOutput(i)-A(idx(i),j)*gaussElimOutput(j);
        end
        % Update ith position
        gaussElimOutput(i) = gaussElimOutput(i)/A(idx(i),i);
    end
end

function [x] = jacobi_method(A, b, tol)
    n = length(b);
    x = zeros(n, 1);
    iters = 0;
    % Loop until desired tolerance
    while true
        iters = iters + 1;
        % Update previous x
        x_prev = x;
        % Compute x values using Jacobi Method
        for i = 1:n
            x(i) = b(i);
            for j = 1:n
                if i ~= j
                    x(i) = x(i) - A(i, j) * x_prev(j);
                end
            end
            x(i) = x(i) / A(i, i);
        end
        % Compute the error and check if it's within the desired tolerance
        err = norm(x - x_prev);
        if err < tol
            break;
        end
    end
    % Display the number of iterations
    fprintf("Jacobi Method Iterations: %d\n", iters);
end

function [outputGS] = gauss_seidel_method(Ain, bin, tol)
    n = length(bin);
    outputGS = zeros(n, 1);
    err = inf;
    iter = 0;
    %Iterate until desired tolerance
    while err>=tol
        %Update previous x values
        xPrev = outputGS;
        %Loop through the rows
        for i = 1:n
            sum = 0;
            %Calculate sum for recently calculated x values
            for j = 1:i-1
                sum = sum + Ain(i,j)*outputGS(j);
            end
            %Calculate sum for previous x values
            for j = i+1:n
                sum = sum + Ain(i,j)*xPrev(j);
            end
            %Update using A(i,i) value
            outputGS(i) = (1/Ain(i,i))*(bin(i)-sum);
        end
        %Update iteration and error
        iter = iter+1;
        err = norm(xPrev-outputGS);
    end
    fprintf("GaussSeidelMethod Iterations: %f\n", iter);
end

%Function to implement SOR Method. Inputs are matrices A and b, desired
%tolerance, and omega value. Output is solution matrix.
function [SORoutput, i] = SOR_method(A, b, tol, w)
    %Create variables
    i = 0;
    currentEr = inf;
    n = length(b);
    SORoutput = zeros(n, 1);
    %Loop until desired tolerance 
    while  tol <= currentEr
        currentEr = 0;
        %Loop through rows
        for i = 1 : n
            sum = 0;
            for j = 1 : n
                %Get to the correct sum by getting rid of extra values
                sum = sum-A(i,j)*SORoutput(j);
            end
            %Update sum using SOR
            sum = w * (sum+b(i)) / A(i,i);
            %Update error
            if abs(sum) > currentEr
                currentEr = abs(sum);
            end
            %Update SOR
            SORoutput(i) = SORoutput(i) + sum;
        end 
        i = i + 1;
    end
    fprintf("SORmethod Iterations: %f\n", i);
end