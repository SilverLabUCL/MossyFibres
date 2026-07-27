function [r_vals, K_vals, L_vals] = ripleyKL2D(xy, letplot, fov_size)
% Computes Ripley's K and L function in 2D with isotropic edge correction
%
% Inputs:
%   xy       : N x 2 matrix of point coordinates (X, Y) in um
%   letplot  : 1 to plot, 0 to suppress
%   fov_size : [width, height] of FOV in um (e.g., [250, 250])
%              If omitted, uses bounding box of point cloud
%
% Outputs:
%   r_vals : vector of radii (um)
%   K_vals : Ripley's K-function
%   L_vals : centered L-function (L(r)-r), CSR baseline = 0

%% Input check
if size(xy, 2) ~= 2
    error('xy must be an N x 2 matrix.');
end
N = size(xy, 1);
if N < 2
    error('Need at least 2 points.');
end

%% Parameters
r_max  = 60;
n_r    = 60;
r_vals = linspace(r_max / n_r, r_max, n_r);

%% Define study area bounding box
if nargin < 3 || isempty(fov_size)
    x0 = min(xy(:,1));  x1 = max(xy(:,1));
    y0 = min(xy(:,2));  y1 = max(xy(:,2));
else
    % FOV origin at (0,0)
    x0 = 0;  x1 = fov_size(1);
    y0 = 0;  y1 = fov_size(2);
end

studyArea = (x1 - x0) * (y1 - y0);
if studyArea <= 0
    error('Study area is zero or negative. Check fov_size or point coordinates.');
end

%% Point density
lambda = N / studyArea;

%% Distance from each point to the four edges (N x 1)
d_left   = xy(:,1) - x0;
d_right  = x1 - xy(:,1);
d_bottom = xy(:,2) - y0;
d_top    = y1 - xy(:,2);

%% Check for points outside FOV
out_of_fov = d_left < 0 | d_right < 0 | d_bottom < 0 | d_top < 0;
if any(out_of_fov)
    warning('%d point(s) lie outside the FOV and will bias the results.', sum(out_of_fov));
end

%% Pairwise distance matrix (N x N)
D = squareform(pdist(xy));

%% Edge correction weight matrix (vectorized)
W = compute_weights(D, d_left, d_right, d_bottom, d_top, N);

%% K function (unbiased denominator: (N-1)*lambda)
K_vals = zeros(size(r_vals));
for k = 1:length(r_vals)
    r = r_vals(k);
    mask = (D <= r) & (D > 0);
    K_vals(k) = sum(W(mask)) / ((N - 1) * lambda);
end

%% Centered L function
L_vals = sqrt(max(K_vals, 0) / pi) - r_vals;

%% Plot
if letplot == 1
    K_csr = pi * r_vals.^2;
    figure('Position', [100 100 750 320]);

    subplot(1,2,1);
    plot(r_vals, K_vals, '-k', 'LineWidth', 2); hold on;
    plot(r_vals, K_csr,  '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
    xlabel('r (um)');  ylabel('K(r)');
    legend('Observed','CSR','Location','best','Box','off');
    set(gca,'LineWidth',1,'FontSize',13);

    subplot(1,2,2);
    plot(r_vals, L_vals,             '-k',  'LineWidth', 2); hold on;
    plot(r_vals, zeros(size(r_vals)), 'r--', 'LineWidth', 2);
    xlabel('r (um)');  ylabel('L(r) - r');
    legend('Observed','CSR','Location','best','Box','off');
    set(gca,'LineWidth',1,'FontSize',13);
end
end


%% ====================================================================
function W = compute_weights(D, d_left, d_right, d_bottom, d_top, N)
% Vectorized isotropic edge correction weight matrix.
% For each point i, computes the proportion of a circle of radius D(i,j)
% that lies within the rectangular FOV, using 360 angular samples.

n_angles  = 360;
theta     = linspace(0, 2*pi, n_angles + 1);
theta     = theta(1:end-1);
cos_theta = cos(theta);
sin_theta = sin(theta);

W = ones(N, N);

for i = 1:N
    dl = d_left(i);   dr = d_right(i);
    db = d_bottom(i); dt = d_top(i);

    r_vec = D(i, :)';       % N x 1

    dx = r_vec * cos_theta; % N x 360
    dy = r_vec * sin_theta;

    inside = (dx >= -dl) & (dx <= dr) & ...
             (dy >= -db) & (dy <= dt);

    prop = sum(inside, 2) / n_angles; % N x 1

    zero_prop = (prop == 0) & (r_vec > 0);
    if any(zero_prop)
        warning('ripleyKL2D:zeroProp', ...
            'Point %d has %d neighbor(s) whose edge-correction circle falls entirely outside FOV; those pairs are excluded.', ...
            i, sum(zero_prop));
    end

    w_row          = zeros(N, 1);
    valid          = (prop > 0) & (r_vec > 0);
    w_row(valid)   = 1 ./ prop(valid);
    W(i, :)        = w_row';
end
end