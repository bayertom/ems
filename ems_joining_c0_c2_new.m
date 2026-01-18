function [xs, ys] = ems_c0_c2(C, Q, n1, n2, h, delta, dist_max, lambda, eps, max_iters, method, cont, stencil)
% Join two splines s1, s2 with energies E1..E5 and enforce
%   cont = 0 : C0  only (position)
%   cont = 1 : C0 + C1 (position + 1st derivative)
%   cont = 2 : C0 + C1 + C2 (position + 1st + 2nd derivative)
% Parameters:
%   C        : (n1+n2)x2 data points [xc yc]
%   Q        : nox2 obstacle boundary samples
%   n1, n2   : number of nodes of the first/second spline
%   h        : length step
%   delta    : E4 offset
%   dist_max : E4/E5 active distance threshold
%   lambda   : [lambda1..lambda5] scalar parameterd for E1..E5
%   eps      : stopping tolerance on ||delta s||_inf
%   max_iters: maximum amount of iterations
%   method   : 1 (backslash on KKT), 2 (LDL), 3 (null-space)
%   cont     : 0,1,2 continuity mode
%   stencil  : 3, 5 stencil size

%Total nodes
n12 = n1 + n2;

%Get coordinates
xc = C(:,1); yc = C(:,2);

%Spline initialization
x = xc; y = yc;

%Stacked vector (2*n12)
s = [x; y];

%Difference matrices
E1 = eye(n1);
D11 = diff(E1,1);             % (n1-1) x n1
D21 = diff(E1,2);             % (n1-2) x n1

E2 = eye(n2);
D12 = diff(E2,1);             % (n2-1) x n2
D22 = diff(E2,2);             % (n2-2) x n2

%Stencil formed by 3 points
if stencil == 3
    [D11, D21] = diffs3pts(n1);
    [D12, D22] = diffs3pts(n2);

%Stencil formed by 5 points
else
    [D11, D21] = diffs5pts(n1, h);
    [D12, D22] = diffs5pts(n2, h);
end

%Block diagonal matrices D1,D2 for the joined spline
D1 = blkdiag(D11, D12);       % (n12-2) x n12
D2 = blkdiag(D21, D22);       % (n12-4) x n12

%Combined derivatives
DS = [D1; D2];                % (2*n12-6) x n12

%Discrete Laplacian operator
L = D1.'*D1 + D2.'*D2;        % n12 x n12

%Junction indices
i_last  = n1;                 % last index of s1
i_first = n1 + 1;             % first index of s2

% C0: position equality p_{n1} = p_{n1+1}
g0 = zeros(1,n12);
g0(i_last)  =  1;
g0(i_first) = -1;

% C1: equality of derivatives p^{1}_{n1} = p^{1}_{n1+1}
if cont >= 1
    %Last forward diff of s1
    d1L = D1(n1-1,:);

    %First forward diff on s2
    d1R = D1(n1,:);

    %Their difference
    g1  = d1L - d1R;
end

% C2: equality of derivatives p^{2}_{n1} = p^{2}_{n1+1}
if cont >= 2
    %Last forward diff of s1
    d2L = D2(n1-2,:);         % last 2nd-diff on s1

    %First forward diff on s2
    d2R = D2(n1-1,:);

    %Their difference
    g2  = d2L - d2R;
end

%Create matrix G
switch cont
    %C0
    case 0
        m = 2;
        G = zeros(m, 2*n12);

        % C0 in x, y
        G(1, 1:n12) = g0;
        G(2, n12+1:2*n12)= g0;

    % C0 + C1
    case 1
        m = 4;
        G = zeros(m, 2*n12);

        % C0 in x,y
        G(1, 1:n12) = g0;
        G(2, n12+1:2*n12) = g0;

        % C1 in x,y
        G(3, 1:n12) = g1;
        G(4, n12+1:2*n12) = g1;

    % C0 + C1 + C2
    case 2
        m = 6;
        G = zeros(m, 2*n12);

        %C0 in x, y
        G(1, 1:n12) = g0;
        G(2, n12 + 1:2*n12) = g0;

        %C1 in x, y
        G(3, 1:n12)  = g1;
        G(4, n12 + 1:2*n12) = g1;

        %C2 in x, y
        G(5, 1:n12) = g2;
        G(6, n12 + 1:2*n12) = g2;

end

%Iteration process
for it = 1:max_iters

    %Create stacked x, y matrix for knn search
    s2  = [s(1:n12), s(n12+1:2*n12)];   % n12 x 2

    %Find nearest obstacle point
    [imin, dmin] = knnsearch(Q, s2, 'K', 1, 'Distance','euclidean');
    QN  = Q(imin,:);
    qnx = QN(:,1); qny = QN(:,2);

    % Normalized displacement u = (q - s)/||q-s||
    U = zeros(n12,2);
    UT = QN - s2;
    nz  = dmin > 0;
    U(nz,:) = UT(nz,:) ./ dmin(nz);
    ux = U(:,1);  uy = U(:,2);

    %Change weights of points closer than dist_max for E4
    w4 = double(dmin < dist_max);
    closer = (w4 > 0);

    %Change weights of points closer than dist_max for E5
    mask1 = (D1 ~= 0); mask2 = (D2 ~= 0);
    w51 = (mask1 * double(closer) > 0);
    w52 = (mask2 * double(closer) > 0);
    W5 = diag([w51; w52]);

    % **************Internal energy*****************************************

    %E1: residuals
    A1  = eye(n12);
    b1x = xc; b1y = yc;

    %E2: smoothing
    A2  = L;
    b2x = zeros(n12,1); b2y = zeros(n12,1);

    %E3: preserving shape of c(t)
    A3 = A2;
    b3x = L * xc; b3y = L * yc;

    %***************External energy*****************************************

    %E4: squares of residuals + offset
    tt = ux.*qnx + uy.*qny - delta;
    A4x = diag((w4 .* ux).^2);
    A4y = diag((w4 .* uy).^2);
    A4xy = diag((w4.^2) .* (ux .* uy));

    b4x = (w4.^2) .* (ux .* tt);
    b4y = (w4.^2) .* (uy .* tt);

    %E5: preserving shape of qk(t)
    A5 = DS' * (W5 * DS);
    b5x = A5 * qnx;
    b5y = A5 * qny;

    %Assembly matrix Ax, Ay, Axy
    Ax = lambda(1)*A1 + lambda(2)*A2 + lambda(3)*A3 + lambda(4)*A4x + lambda(5)*A5;
    Ay = lambda(1)*A1 + lambda(2)*A2 + lambda(3)*A3 + lambda(4)*A4y + lambda(5)*A5;
    Axy = lambda(4)*A4xy;

    %Assembly vector bx, by
    bx = lambda(1)*b1x + lambda(2)*b2x + lambda(3)*b3x + lambda(4)*b4x + lambda(5)*b5x;
    by = lambda(1)*b1y + lambda(2)*b2y + lambda(3)*b3y + lambda(4)*b4y + lambda(5)*b5y;

    %Full system A s = b
    A2D = [Ax,  Axy; Axy, Ay];
    b2D = [bx; by];

    %Avoid singularity, add diagonal matrix
    A2D = A2D + 1e-12*eye(2*n12);

    %Create KKT: [A2D G'; G 0] [s; lambdaK] = [b2D; 0]
    K = [A2D, G.'; G, zeros(m,m)];
    bs = [b2D; zeros(m,1)];

    % Solve KKT, GE
    if method == 1
        sol = K \ bs;

    % Solve KKT, LDL
    elseif method == 2
        [Lk, Dk, Pk] = ldl(K);
        rhs_perm = Pk.' * bs;
        v = Lk \ rhs_perm;
        w = Dk \ v;
        z = Lk.' \ w;
        sol = Pk * z;

    % Solve KKT, Null-space
    else
        Z = null(G,'r');   
        Ared = Z.' * A2D * Z;
        bred = Z.' * b2D;
        y_ns = Ared \ bred;
        s_ns = Z * y_ns;
        sol  = [s_ns; zeros(m,1)];
    end

    % Extract spline
    s_new = sol(1:2*n12);

    %Stopping condition
    step = norm(s_new - s, inf);

    %Print error
    fprintf('Outer %2d: ||delta s||_inf = %.3e\n', it, step);

    %Update the solution
    s = s_new;

    %Stop iteration process
    if step < eps
        break;
    end
end

%Get coordinates
xs = s(1:n12); ys = s(n12+1:2*n12);

end

clc; clear; close all;
hold on; grid on;

%Input parameters of c(t)
n1 = 30; n2 = 30;
t1 = linspace(0,0.8,n1);
t2 = linspace(0,0.8,n2);
h = 1/(n1-1);

%Curves c1, c2
c1 = [t1.',  zeros(n1)+0.2];
c2 = [1.0 + t2.', 0.8 + zeros(n2)];
xc = [c1(:,1); c2(:,1)];
yc = [c1(:,2); c2(:,2)];

%Two circular obstacles
t = 0:0.01:2*pi;
x1 = 0.6; y1 = 0.8; r1 = 0.2;
x2 = 1.2; y2 = 0.2; r2 = 0.2;
Q1 = [(x1 + r1*cos(t))', (y1 + r1*sin(t))'];
Q2 = [(x2 + r2*cos(t))', (y2 + r2*sin(t))'];

%Stacked vectors
C = [xc, yc];
Q = [Q1; Q2];

%Input parameters of s(t)
delta = 0.05; dist_max  = 0.50;
lambda = [1 1 5 1 0.0];
max_iters = 20;
eps = 1e-9;
method = 1;
stencil = 3;

%C0
[xc0,yc0] = ems_c0_c2(C, Q, n1, n2, h, delta,dist_max, lambda, eps, max_iters, method, 0, stencil);

%C0 + C1
[xc1,yc1] = ems_c0_c2(C, Q, n1, n2, h, delta, dist_max, lambda, eps, max_iters, method, 1, stencil);

%C0 + C1 + C2
[xc2,yc2] = ems_c0_c2(C, Q, n1, n2, h, delta, dist_max, lambda, eps, max_iters, method,2, stencil);

%Plot curve c(t)
plot(c1(:,1), c1(:,2), '-', 'LineWidth',1.3, 'DisplayName','c_1(t) input');
plot(c2(:,1), c2(:,2), '-', 'LineWidth',1.3, 'DisplayName','c_2(t) input');

% Plot C0 spline
plot(xc0(1:n1), yc0(1:n1), '-', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName','C^0 (s_1)');
plot(xc0(n1+1:end), yc0(n1+1:end), '-', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName','C^0 (s_2)');

% Plot C1 spline
plot(xc1(1:n1), yc1(1:n1), '-', 'LineWidth',1.5, 'MarkerSize', 4, 'DisplayName','C^0+C^1 (s_1)');
plot(xc1(n1+1:end), yc1(n1+1:end),'-', 'LineWidth',1.5, 'MarkerSize', 4, 'DisplayName','C^0+C^1 (s_2)');

% Plot C2 spline
plot(xc2(1:n1), yc2(1:n1), '-', 'LineWidth',1.5, 'MarkerSize', 4, 'DisplayName','C^0+C^1+C^2 (s_1)');
plot(xc2(n1+1:end), yc2(n1+1:end),'-', 'LineWidth',1.5, 'MarkerSize', 4, 'DisplayName','C^0+C^1+C^2 (s_2)');

%Plot obstacles
plot(Q1(:,1), Q1(:,2), 'k-', 'LineWidth',1.25, 'DisplayName','Obstacle 1');
plot(Q2(:,1), Q2(:,2), 'k-', 'LineWidth',1.25, 'DisplayName','Obstacle 2');

% Junction points
plot(xc0(n1), yc0(n1), 'ko', 'MarkerFaceColor','r', 'DisplayName','junction C^0');
plot(xc1(n1), yc1(n1), 'ko', 'MarkerFaceColor','g', 'DisplayName','junction C^1');
plot(xc2(n1), yc2(n1), 'ko', 'MarkerFaceColor','b', 'DisplayName','junction C^2');

axis equal; xlabel('x'); ylabel('y');
title('Joined splines with C^0, C^1, C^2 continuity');
legend('Location','bestoutside');

function [D1_5, D2_5] = diffs5pts(n, h)
% 5-point central differences for interior nodes i=3..n-2
% First derivative:  (s_{i-2} - 8*s_{i-1} + 8*s_{i+1} - s_{i+2})/(12*h)
% Second derivative: (-s_{i-2} + 16*s_{i-1} - 30*s_i + 16*s_{i+1} - s_{i+2})/(12*h^2)

n2 = n - 4;          
D1_5 = zeros(ni, n);
D2_5 = zeros(ni, n);

%Evaluate difference matrices
for k = 1:n2
    i = k + 2;

    %First derivative
    D1_5(k, i-2) = 1/(12*h);
    D1_5(k, i-1) = -8/(12*h);
    D1_5(k, i+1) = 8/(12*h);
    D1_5(k, i+2) = -1/(12*h);

    %Second derivative
    D2_5(k, i-2) = -1/(12*h^2);
    D2_5(k, i-1) = 16/(12*h^2);
    D2_5(k, i) = -30/(12*h^2);
    D2_5(k, i+1) = 16/(12*h^2);
    D2_5(k, i+2) = -1/(12*h^2);
end
end

function [D1_3, D2_3] = diffs3pts(n)
% 5-point central differences for interior nodes i=2..n-1
E1 = eye(n);

%Get difference matrices
D1_3 = diff(E1,1);
D2_3 = diff(E1,2);
end
