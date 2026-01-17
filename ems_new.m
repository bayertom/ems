function [xs, ys] = ems(C, Q, n, delta, dist_max, lambda, eps, max_iters, k)

%Get coordinates
xc = C(1:n); yc = C(n+1:2*n);

%Spline initialization
x = xc; y = yc; 
s  = [x; y];

%Difference matrices
E = eye(n, n);
D1 = diff(E, 1);
D2 = diff(E, 2);

%Smoothing matrix
if k==1
    D12 = sqrt(2)/2*[D1; D1];
elseif k==2
    D12 = sqrt(2)/2*[D2; D2];
else
    D12 = [D1; D2];
end

%Discrete Laplacian operator
if k==1  %First derivative
    L =  D1.' * D1;
elseif k==2 %Second derivative
    L = D2.' * D2;
else %First and second derivative
    L = (D1.' * D1) + (D2.' * D2);
end

%Iteration process
for it = 1:max_iters

    %Create stacked x, y matrix
    s2 = [s(1:n), s(n+1:2*n)];

    % %Find nearest obstacle point
    knn = 1;
    [imin, dmin] = knnsearch(Q, s2, 'K', knn, 'Distance','euclidean');
    QN = Q(imin, :);
    qnx = QN(:, 1); qny = QN(:, 2);

    %Compute normalized displacement vectors: u = qk - s
    U = zeros(n, 2);
    UT = QN - s2;
    nz = dmin > 0;
    U(nz,:) = UT(nz,:) ./ dmin(nz);

    %Change weights of points closer than dist_max for E4
    w4 = double(dmin < dist_max);
    closer_idxs = (w4 > 0);

    %Change weights of points closer than dist_max for E5
    mask1 = (D1 ~= 0); mask2 = (D2 ~= 0);
    w51 = double( (mask1 * double(closer_idxs)) > 0 );
    w52 = double( (mask2 * double(closer_idxs)) > 0 );

    %First derivative
    if k==1
        W5 = diag([w51; w51]);
    %Second derivative
    elseif k==2
        W5 = diag([w52; w52]);
    %First and derivatives
    else
        W5 = diag([w51; w52]);
    end

    % **************Internal energy*****************************************

    % E1: residuals
    A1_x = eye(n);
    b1x  = xc; b1y  = yc;
    A1 = blkdiag(A1_x, A1_x);

    %E2: smoothing
    A2_x = L;
    b2x = zeros(n,1); b2y = zeros(n,1);
    A2 = blkdiag(A2_x, A2_x);

    %E3: preserving shape of c(t)
    A3 = A2;
    b3x = L*xc; b3y = L*yc;

    %***************External energy*****************************************

    %E4: squares of residuals + offset
    ux = U(:,1); uy = U(:,2);
    tt = ux.*qnx + uy.*qny - delta;
    A4x = diag((w4.*ux).^2);
    A4y = diag((w4.*uy).^2);
    A4xy = diag((w4.^2).*ux.*uy);
    A4 = [A4x, A4xy; A4xy A4y];

    b4x = (w4.^2).*ux.*tt;
    b4y = (w4.^2).*uy.*tt;

    %E5: preserving shape of qk(t)
    A5_x = D12' * W5 * D12;
    b5x = A5_x * qnx;
    b5y = A5_x * qny;
    A5 = blkdiag(A5_x, A5_x);

    %Assembly matrix Ax, Ay, Axy
    A = lambda(1)*A1 + lambda(2)*A2 + lambda(3)*A3 + lambda(4)*A4 + lambda(5)*A5;

    %Assembly vector bx, by
    bx = lambda(1)*b1x + lambda(2)*b2x + lambda(3)*b3x + lambda(4)*b4x+ lambda(5)*b5x ;
    by = lambda(1)*b1y + lambda(2)*b2y + lambda(3)*b3y + lambda(4)*b4y+ lambda(5)*b5y;

    %Assembly matrix A
    b = [bx; by];

    % Avoid singularity, add diagonal matrix
    A = A + 1e-12*eye(2*n);

    %Find solution (simple inversion)
    s_new = inv(A) * b;

    %Stopping condition
    step = norm(s_new - s, inf);
    fprintf('Outer %2d: ||delta s||_inf = %.3e\n', it, step);
    if step < eps
        break;
    end

    %Update the solution
    s = s_new;

end

%Get coordinates
xs = s(1:n);  ys = s(n+1:2*n);
end

%*************************************************************
clc
clear
format long g
axis equal
hold on

%Testing example
r = 2;
dt1 = 0.01;
t = 0:dt1*pi:2*pi;

%Create obstacles qk(t)
xq1 = r *cos(t);yq1 = r*sin(t);
xq2 = xq1 + 10; yq2 = yq1;
xq3 = xq1 + 10;yq3 = yq1 + 10;
xq4 = xq1; yq4 = yq1 + 10;
xq5 = xq1 + 5; yq5 = yq1 + 5;
xq6 = xq2 + 10; yq6 = yq2 + 0;
xq7 = xq3 + 10;yq7 = yq3 + 0;
xq8 = xq5 + 10;yq8 = yq5 + 0;

%Approximated curve c(t)
dt2 = 0.1;
s =  -2*pi:dt2:8*pi;
xc = (s + 0.95)';
yc = (3.9*sin(s/3) + 5)';

%Plot obstacles
hold on
plot (xq1, yq1, 'k');
plot (xq2, yq2, 'k');
plot (xq3, yq3, 'k');
plot (xq4, yq4, 'k');
plot (xq5, yq5, 'k');
plot (xq6, yq6, 'k');
plot (xq7, yq7, 'k');
plot (xq8, yq8, 'k');

%Plot curve
plot (xc, yc, 'k');

%Stacked vectors
C = [xc; yc];
Q = [[xq1', yq1']; [xq2', yq2']; [xq3', yq3']; [xq4', yq4']; [xq5', yq5']; [xq6', yq6']; [xq7', yq7']; [xq8', yq8']];

%Input parameters of s(t) for estimators
delta = 0.30; dist_max = 3.0;
lambda = [1.0,  10.0,  1.0,  1.0,  1];
max_iters = 1;
eps = 1.0e-9;

%Construct spline
k = 12; %First and second derivative
n = length(xc);
[xs, ys] = ems(C, Q, n, delta, dist_max, lambda, eps, max_iters, k);

%Plot results
hold on

%figure('Color','w'); hold on; grid on;
plot(C(1:n), C(n+1:2*n), '-b', 'LineWidth',1.3, 'DisplayName','c_1(t) input');
plot(xs(1:n), ys(1:n), '-r', 'LineWidth',2, 'MarkerSize',5, 'DisplayName','s_1 nodes');

axis equal; xlabel('x'); ylabel('y');





