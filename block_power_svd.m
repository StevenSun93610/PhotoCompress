function [R_U, R_V, U_k, V_k] = block_power_svd(A, r, tol)
% initialize V_0 as the first
% r columns of the n × n identity matrix
n = size(A, 2);
V_0 = eye(n);
V_0 = V_0(:, 1:r);
% variable set up for loop
V_previous = V_0;
res = Inf;
%variable set up for returning
R_U = zeros(n);
R_V = zeros(n);
U_k = [];
V_k = [];
while res > tol
    U_hat = A * V_previous;
    % Get reduced QR factorization U_hat
    [Q_U, R_U] = qrPosDiagR(U_hat);
    U_k = Q_U;
    V_hat = A' * U_k;
    % Get reduced QR factorization V_hat
    [Q_V, R_V] = qrPosDiagR(V_hat);
    V_k = Q_V;
    V_previous = V_k;
    % get the norm
    res = max(norm(triu(R_U,1)), norm(triu(R_V,1)));
end


end

