clc
clear variables
close all

%% Re-weight of particle size distribution
%   Theory
%   If one has a PSD weighted by, say, number but wants to instead view the
%   volume-weighted distribution (or some other weighting) it would be nice
%   to have concise formulas to express this transformation.
%   
%   Here, we assume we start with q_0(d,l) and want to compute 
%   q_2(d,l) --> weighted by surface area of cylinder
%   q_3(d,l) --> weighted by volume of cylinder
%
%   The conversion is done via the generalized moment method
%   q_r(d,l) = k_r(d,l)q_0(d,l)/ int_0^Inf k_r(d,l)q_0(d,l) dd dl
%   In this case
%   q_2(d,l) --> k_2(d,l) = pi*d*l + 1/2 pi*d^2
%   q_3(d,l) --> k_3(d,l) = pi/4 * d^2 * l
%
%   One can then show that
%   q_2(d,l) = w1*LN(mu+Sigma*[1;1], Sigma) + w2*LN(mu+Sigma*[2;0],Sigma)
%   q_3(d,l) = LN(mu+Sigma*[2;1], Sigma)
%   where w1 and w2 are weighting factors dependent on Sigma and mu

% q_0 distribution

avg_d = 20;                                                                 % median diameter in nm
avg_l = 100;                                                                % median length in nm

spread_d = 1.1;                                                            % 68% of prob in [avg_d/spread_d, avg_d*spread_d] and 95% in [avg_d/spread_d^2, avg_d*spread_d^2]
spread_l = 1.1;                                                            % 68% of prob in [avg_l/spread_l, avg_l*spread_l] and 95% in [avg_l/spread_l^2, avg_l*spread_l^2]
corr     = 0.3;                                                             % Correlation between parameters for the 2D distribution

% Log-normal distribution parameters
mu = log([avg_d, avg_l]);
Sigma = [log(spread_d), corr*log(spread_d)*log(spread_l);
    corr*log(spread_d)*log(spread_l), log(spread_l)];

% Function for graphing
q0 = @(d,l) fcn_LN2D(d,l,mu,Sigma);

%% Surface Area re-weight

% Adjustments to mean values
bsa_1 = [1;1];
bsa_2 = [2;0];

mu_s1 = mu' + Sigma*bsa_1;
mu_s2 = mu' + Sigma*bsa_2;

% Compute weighting factors for each log-normal contribution
c1 = exp( 1/2 * bsa_1' * (Sigma*bsa_1 - 2*mu) );
c2 = exp( 1/2 * bsa_2' * (Sigma*bsa_2 - 2*mu) );

norm_factor = pi * c1 + 1/2 * pi * c2;

w1 = pi * c1 / norm_factor;
w2 = 1/2 * pi * c2 / norm_factor;

% Function for graphing
q2 = @(d,l) w1*fcn_LN2D(d,l,mu_s1',Sigma) + w2*fcn_LN2D(d,l,mu_s2',Sigma);

%% Volume re-weight

% Adjustment to mean value
bv = [2;1];
mu_v = mu' + Sigma*bv;

% Function for graphing
q3 = @(d,l) fcn_LN2D(d,l,mu_v',Sigma);



%%  Generate Data for Joint Plot

% Generate grid for plot
dlims = 0:0.1:60;
llims = 0:0.1:250;
[D,L] = meshgrid(dlims,llims);

% Evaluate q0, q2, q3
q0_vals = q0(D,L);
q2_vals = q2(D,L);
q3_vals = q3(D,L);

% Find mode of each
q0_mode = fcn_FindMode(q0,mu);
q2_mode = fcn_FindMode(q2,mu_s1);
q3_mode = fcn_FindMode(q3,mu_v);


% Find 95% probability region for each
confidence_level = 0.95;
q0_Dens95 = fcn_FindConfidenceCountour(q0,confidence_level);
q2_Dens95 = fcn_FindConfidenceCountour(q2,confidence_level);
q3_Dens95 = fcn_FindConfidenceCountour(q3,confidence_level);

%%  Create Joint Plot

% Specify colors for each PSD
color_q0 = [0 0.4470 0.7410];
color_q2 = [0.8500 0.3250 0.0980];
color_q3 = [0.9290 0.6940 0.1250];

% Contour plots
mode_size = 150;
figure
contour(D,L,q0_vals,[q0_Dens95 q0_Dens95],...
    'LineColor',color_q0,...
    'LineWidth',10,...
    'DisplayName','q_0(d,l)')
hold on

contour(D,L,q2_vals,[q2_Dens95 q2_Dens95],...
    'LineColor',color_q2,...
    'LineWidth',10,...
    'DisplayName','q_2(d,l)')

contour(D,L,q3_vals,[q3_Dens95 q3_Dens95],...
    'LineColor',color_q3,...
    'LineWidth',10,...
    'DisplayName','q_3(d,l)')

scatter(q0_mode(1),q0_mode(2),...
    mode_size, color_q0, 'filled','pentagram',...
    'HandleVisibility','off')

scatter(q2_mode(1),q2_mode(2),...
    mode_size, color_q2, 'filled','pentagram',...
    'HandleVisibility','off')

scatter(q3_mode(1),q3_mode(2),...
    mode_size, color_q3, 'filled','pentagram',...
    'HandleVisibility','off')

hold off
xlim([0 60])
xlabel('Diameter / nm')
ylabel('Length / nm')
fontsize(gcf,scale=1.8)
hLegend = legend('Location','northeast');
set(hLegend, 'Color','none');

print(gcf,"figure2","-dsvg")
print(gcf,'-depsc','-image',"figure2");
%%  Create Individual Plots

plot_range = [0 100 0 300];

figure
colormap('parula')
fcontour(q0,...
    plot_range,...
    'MeshDensity',100,...
    'Fill','on')
colorbar
xlabel('Diameter / nm')
ylabel('Length / nm')
title('q_0(d,l)')
fontsize(gcf,scale=1.5)

figure
colormap('parula')
fcontour(q2,...
    plot_range,...
    'MeshDensity',100,...
    'Fill','on')
colorbar
xlabel('Diameter / nm')
ylabel('Length / nm')
title('q_2(d,l)')
fontsize(gcf,scale=1.5)

figure
colormap('parula')
fcontour(q3,...
    plot_range,...
    'MeshDensity',100,...
    'Fill','on')
colorbar
xlabel('Diameter / nm')
ylabel('Length / nm')
title('q_3(d,l)')
fontsize(gcf,scale=1.5)

%%  Helper Functions

%   Probability density function for log-normal
function pdf = fcn_LN2D(D,L,mu,sigma)
    X = log(D);
    Y = log(L);
    x = X(:);
    y = Y(:);
    R = [x,y];
    psd = mvnpdf(R,mu,sigma)./exp(x)./exp(y);
    [r,c] = size(X);
    pdf = reshape(psd,[r,c]);
end

%   Find maximum of PSD
function mode = fcn_FindMode(fcn_psd,guess)
    opt_fcn = @(x) -fcn_psd(x(1),x(2));
    mode = fminsearch(opt_fcn, guess);
end

%   Search for confidence region based on maximum density value
function density_threshold = fcn_FindConfidenceCountour(fcn_pdf, level)
    fcn_integrand = @(x,y,t) (fcn_pdf(x,y)>=t).*fcn_pdf(x,y);
    num_draws = 5e6;
    rng(1)
    mc_draws = unifrnd(1e-1,500,num_draws,2);
    volume = 500*500;
    mc_integral = @(t) volume / num_draws * sum( fcn_integrand(mc_draws(:,1),mc_draws(:,2),t) );
    fcn_FindRoot = @(t) mc_integral(t) - level;
    % fcn_FindRoot = @(t) integral2(@(x,y) fcn_integrand(x,y,t),1e-1,500,1e-1,500) - level;
    tic
    [density_threshold, fval, ~,~] = fzero(fcn_FindRoot,0);
    toc

end