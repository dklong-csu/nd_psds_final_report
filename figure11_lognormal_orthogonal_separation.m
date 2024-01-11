clc
clear variables
close all

%%  Physical Constants

%   Density of H2O is 0.99704702 g/mL @ 25C
rhoH2O = 0.99704702 * 1e-21;    % converted to g/nm^3
%   Bulk density of gold is 19.3 g/cm^3
rhoAu = 19.3 * 1e-21;           % converted to g/nm^3
%   Viscosity of water is 1.0016 mPa*s @ 20C
eta = 1.0016e-7;            % converted to mPa*s = g/nm/s
%   Surface density charge
rhoch = 1.7e-14;


%%  Definition of qfeed

avg_d = 40;                                                                 % median diameter in nm
avg_l = 100;                                                                % median length in nm

spread_d = 1.1;                                                             % 68% of prob in [avg_d/spread_d, avg_d*spread_d] and 95% in [avg_d/spread_d^2, avg_d*spread_d^2]
spread_l = 1.1;                                                             % 68% of prob in [avg_l/spread_l, avg_l*spread_l] and 95% in [avg_l/spread_l^2, avg_l*spread_l^2]
corr     = 0.3;                                                             % Correlation between parameters for the 2D distribution

% Log-normal distribution parameters
mu = log([avg_d, avg_l]);
Sigma = [log(spread_d), corr*log(spread_d)*log(spread_l);
    corr*log(spread_d)*log(spread_l), log(spread_l)];

% Function for graphing
fcn_Qfeed = @(d,l) fcn_LN2D(d,l,mu,Sigma);


%%  Graph visualizing the two cuts on the feed

%   Sample cut parameters
sed_lo = 2.5e-11;
sed_hi = 3.5e-11;
fcn_SedCalc = @(d,l) fcn_SedimentationCoeff(d,l,rhoAu, rhoH2O, eta);
fcn_Tsed = @(d,l) fcn_1DEfficiency(d,l,fcn_SedCalc,sed_lo,sed_hi);

mobil_lo = 2.0e-6;
mobil_hi = 4.0e-6;
fcn_MobCalc = @(d,l) fcn_ElectronMobility(d,l,eta,rhoch);
fcn_Tmob = @(d,l) fcn_1DEfficiency(d,l,fcn_MobCalc,mobil_lo,mobil_hi);

%   Plots of each cut in 1D and the cuts against the feed



%%  Optimization to maximize the area of a particular region

%   Region 1, 2, 3, or 4
maximize_region = 2;

%   Starting parameter values
%   Computations work better when parameters are ~1 so these are normalized
%   within the optimization routine
sed_lo_start = 2.5e-11;
sed_hi_start = 3.5e-11;
sed_mag = 1.0e-11;

mobil_lo_start = 2.0e-6;
mobil_hi_start = 4.0e-6;
mobil_mag = 1.0e-6;
x0 = [ sed_lo_start/sed_mag, sed_hi_start/sed_mag, mobil_lo_start/mobil_mag, mobil_hi_start/mobil_mag];

%   Some (example) bounds on the parameters to be realistic
sed_lb = 2e-11;
sed_ub = 4e-11;
mobil_lb = 1e-6;
mobil_ub = 5e-6;
lb = [ sed_lb/sed_mag, sed_lb/sed_mag, mobil_lb/mobil_mag, mobil_lb/mobil_mag];
ub = [ sed_ub/sed_mag, sed_ub/sed_mag, mobil_ub/mobil_mag, mobil_ub/mobil_mag];

%   Enforce sed_lo <= sed_hi and mobil_lo <= mobil_hi
%   i.e. sed_lo - sed_hi <= 0, and mobil_lo - mobil_hi <= 0
%   This is formulated as a linear constraint
%   For simplicity let x = [sed_lo; sed_hi; mobil_lo; mobil_hi]
%   Then for A=[1,-1,0,0;0,0,1,-1] and b=[0;0]
%   The desired linear constraints follow: Ax <= b
%   which is the format fmincon wants
constrA = [1,-1,0,0;0,0,1,-1];
%   However, realistically, one can only implement a physical process where
%   the separation parameters are a minimum of some fixed amount apart.
%   Here let's pretend this separation can be performed such that the width
%   is at a minimum half a decimal point.
%   Then -- since the parameters are normalized to be roughly 1 -- the b
%   vector would negative of the accessible threshold. So b=[-0.1;-0.1].
constrB = [-0.5;-0.5];

%   Desired area of the specified region
%   The objective function will compute the area A of 'maximize_region'
%   and compare this to the desired area. The objective function to be
%   minimized is therefore
%       Cost = (A - desired_area)^2
%   Desired_area should be in [0,1] (since the area is a probability)
%   Desired_area = 0 or 1 will likely result in the cuts not actually
%   separating any particles, so moderate values between are recommended.
desired_area = 0.25;

%   Define optimization problem
costFcn = @(x) fcn_MaximizeAreaCostFcn(x, [sed_mag, mobil_mag], fcn_Qfeed, maximize_region, rhoAu, rhoH2O, eta, rhoch, 0.25);
useParallel = true;                                                         %   This can help makes things faster
options = optimoptions("fmincon",'Display','iter-detailed', 'UseParallel',true);
optimal_prm = fmincon(costFcn, x0,...
    constrA, constrB,...
    [], [],...
    lb, ub,...
    [],...
    options);

%%  Define the optimized cut

%   Cut based on sedimentation coefficient
sed_lo = optimal_prm(1)*sed_mag;
sed_hi = optimal_prm(2)*sed_mag;
fcn_SedCalc = @(d,l) fcn_SedimentationCoeff(d,l,rhoAu, rhoH2O, eta);
fcn_Tsed = @(d,l) fcn_1DEfficiency(d,l,fcn_SedCalc,sed_lo,sed_hi);

%   Cut based on electron mobility
mobil_lo = optimal_prm(3)*mobil_mag;
mobil_hi = optimal_prm(4)*mobil_mag;
fcn_MobCalc = @(d,l) fcn_ElectronMobility(d,l,eta,rhoch);
fcn_Tmob = @(d,l) fcn_1DEfficiency(d,l,fcn_MobCalc,mobil_lo,mobil_hi);

fprintf("The optimized parameters are:\n")
fprintf("Sedimentation coefficient lower bound: %.6e\n",sed_lo)
fprintf("Sedimentation coefficient upper bound: %.6e\n",sed_hi)
fprintf("Electrophoretic mobility lower bound:  %.6e\n",mobil_lo)
fprintf("Electrophoretic mobility upper bound:  %.6e\n",mobil_hi)

%%  COMPUTATIONAL WARNING
fprintf("Calculating Separated PSDs (this may take a few minutes)...")


%%  PSD for region I
r = 1;
%   Normalizing factor
fcn_Tr1 = fcn_ComputeRegionEfficiency(r,fcn_Tsed, fcn_Tmob);
A1 = fcn_ComputeCutArea(fcn_Qfeed, fcn_Tr1);
%   PSD
fcn_Qr1 = @(d,l) fcn_Tr1(d,l) .* fcn_Qfeed(d,l) / max(A1,1e-12);            %   The max function is for numerical stability

%%  PSD for region II
r = 2;
%   Normalizing factor
fcn_Tr2 = fcn_ComputeRegionEfficiency(r,fcn_Tsed, fcn_Tmob);
A2 = fcn_ComputeCutArea(fcn_Qfeed, fcn_Tr2);
%   PSD
fcn_Qr2 = @(d,l) fcn_Tr2(d,l) .* fcn_Qfeed(d,l) / max(A2,1e-12);            %   The max function is for numerical stability

%%  PSD for region III
r = 3;
%   Normalizing factor
fcn_Tr3 = fcn_ComputeRegionEfficiency(r,fcn_Tsed, fcn_Tmob);
A3 = fcn_ComputeCutArea(fcn_Qfeed, fcn_Tr3);
%   PSD
fcn_Qr3 = @(d,l) fcn_Tr3(d,l) .* fcn_Qfeed(d,l) / max(A3,1e-12);            %   The max function is for numerical stability

%%  PSD for region IV
r = 4;
%   Normalizing factor
fcn_Tr4 = fcn_ComputeRegionEfficiency(r,fcn_Tsed, fcn_Tmob);
A4 = fcn_ComputeCutArea(fcn_Qfeed, fcn_Tr4);
%   PSD
fcn_Qr4 = @(d,l) fcn_Tr4(d,l) .* fcn_Qfeed(d,l) / max(A4,1e-12);            %   The max function is for numerical stability

fprintf("done!\n")
%%  Create Plots
fprintf("Creating plots...")

nIDs = 5;
alphabet = ('a':'z').';
chars = num2cell(alphabet(1:nIDs));
chars = chars.';
charlbl = strcat('(',chars,')'); 

figure('Position',[1,49,1920,955])
t = tiledlayout(2,3,...
    "TileSpacing","tight",...
    "Padding","loose");

dlims = [5 75];
llims = [5 190];
plotdlims = [0 80];
plotllims = [0 200];
show_spread = true;

%   Feed PSD
fcn_MakeSeparationPlot(t, fcn_Qfeed, ...
    fcn_SedCalc, [sed_lo sed_hi],...
    fcn_MobCalc, [mobil_lo mobil_hi],...
    dlims, llims,...
    plotdlims, plotllims,...
    linspace(1e-4,1.5e-3),...
    "Feed",...
    show_spread, ...
    true,...
    [2 1])
hold on
text(0.025,0.1,charlbl{1},'Units','normalized','FontSize',12)
hold off

%   Region I
fcn_MakeSeparationPlot(t, fcn_Qr1, ...
    fcn_SedCalc, [sed_lo sed_hi],...
    fcn_MobCalc, [mobil_lo mobil_hi],...
    dlims, llims,...
    plotdlims, plotllims,...
    linspace(1e-4,1.5e-3),...
    "Region I",...
    show_spread, ...
    false,...
    [1 1])
hold on
text(0.025,0.1,charlbl{2},'Units','normalized','FontSize',12)
hold off

%   Region II
fcn_MakeSeparationPlot(t, fcn_Qr2, ...
    fcn_SedCalc, [sed_lo sed_hi],...
    fcn_MobCalc, [mobil_lo mobil_hi],...
    dlims, llims,...
    plotdlims, plotllims,...
    linspace(1e-4,1.5e-3),...
    "Region II",...
    show_spread, ...
    false,...
    [1 1])
hold on
text(0.025,0.1,charlbl{3},'Units','normalized','FontSize',12)
hold off

%   Region III
fcn_MakeSeparationPlot(t, fcn_Qr3, ...
    fcn_SedCalc, [sed_lo sed_hi],...
    fcn_MobCalc, [mobil_lo mobil_hi],...
    dlims, llims,...
    plotdlims, plotllims,...
    linspace(1e-3,1.5e-2),...
    "Region III",...
    show_spread, ...
    false,...
    [1 1])
hold on
text(0.025,0.1,charlbl{4},'Units','normalized','FontSize',12)
hold off

%   Region IV
fcn_MakeSeparationPlot(t, fcn_Qr4, ...
    fcn_SedCalc, [sed_lo sed_hi],...
    fcn_MobCalc, [mobil_lo mobil_hi],...
    dlims, llims,...
    plotdlims, plotllims,...
    linspace(1e-4,1.5e-3),...
    "Region IV",...
    show_spread, ...
    false,...
    [1 1])
hold on
text(0.025,0.1,charlbl{5},'Units','normalized','FontSize',12)
hold off

% xlabel(t,'Diameter / nm')
% ylabel(t,'Length / nm')

fontsize(gcf,scale=2.25)

fprintf("done!\n")
fprintf("Saving figure to 'figure11.svg' and 'figure11.eps'...")
print("figure11","-dsvg")
print(gcf,'-depsc','-image',"figure11");
fprintf("done!\n")


%%  Helper Functions
function V = fcn_CylinderVolume(d,l)
    V = pi/4 * d.^2 .* l;
end

function A = fcn_CylinderSurfaceArea(d,l)
    A = pi/2 * d.^2 + pi*d.*l;
end

function Xv = fcn_VolumeEquivDiam(d,l)
    V = fcn_CylinderVolume(d,l);
    Xv = nthroot(6/pi * V, 3);
end

function T = fcn_1DEfficiency(d,l, fcn_Conversion, v_lower,v_upper)
    V = fcn_Conversion(d,l);
    m = 1./(v_upper-v_lower);
    b = -v_lower./(v_upper-v_lower);
    T = min(1,max(0,m.*V+b));
end

function f = fcn_FrictionRatio(d,l)
    y = log(l./d);
    f = 1.0304 + ...
        0.0193*y + ...
        0.06229*y.^2 + ...
        0.00476*y.^3 + ...
        0.00166*y.^4 + ...
        2.66e-6*y.^5;
end

function S = fcn_SedimentationCoeff(d,l,rhoParticle,rhoSolvent,viscositySolvent)
    V = fcn_CylinderVolume(d,l);
    mass = V .* rhoParticle;
    fricfactor = (1 - rhoSolvent/rhoParticle);
    numer = mass .* fricfactor;
    denom = 3*pi*viscositySolvent*fcn_FrictionRatio(d,l).*fcn_VolumeEquivDiam(d,l);
    S = numer ./ denom;
end

function mu = fcn_ElectronMobility(d,l,viscositySolvent,surfDensityCharge)
    A = fcn_CylinderSurfaceArea(d,l);
    Q = A*surfDensityCharge;
    denom = 3*pi*viscositySolvent*fcn_FrictionRatio(d,l).*fcn_VolumeEquivDiam(d,l);
    mu = Q./denom;
end

function A = fcn_ComputeCutArea(fcn_qfeed, fcn_Efficiency)
    fcn_integrand = @(d,l) fcn_qfeed(d,l) .* fcn_Efficiency(d,l);
    A = integral2(fcn_integrand,1e-6,Inf,1e-6,Inf);
end

function fcn_Teff = fcn_ComputeRegionEfficiency(region, fcn_T1, fcn_T2)
    switch region
        case 1
            fcn_Teff = @(d,l) fcn_T1(d,l) .* fcn_T2(d,l);
        case 2
            fcn_Teff = @(d,l) fcn_T1(d,l) .* (1-fcn_T2(d,l));
        case 3
            fcn_Teff = @(d,l) (1-fcn_T1(d,l)) .* fcn_T2(d,l);
        case 4
            fcn_Teff = @(d,l) (1-fcn_T1(d,l)) .* (1-fcn_T2(d,l));
        otherwise
            msg = 'Error: Region must be 1, 2, 3, or 4';
            error(msg)
    end
end

function A = fcn_ComputeRegionArea(prm, fcn_qfeed, region, rhoParticle, rhoSolvent, viscositySolvent, surfDensityCharge)
    % prm = [sedLo, sedHi, mobilLo, mobilHi]

    %   Define cut based on sedimentation coefficient
    fcn_Sed = @(d,l) fcn_SedimentationCoeff(d,l,rhoParticle,rhoSolvent,viscositySolvent);
    sedLo = prm(1);
    sedHi = prm(2);
    fcn_TSed = @(d,l) fcn_1DEfficiency(d,l,fcn_Sed,sedLo, sedHi);

    %   Define cut based on electron mobility
    fcn_Mobil = @(d,l) fcn_ElectronMobility(d,l, viscositySolvent, surfDensityCharge);
    mobilLo = prm(3);
    mobilHi = prm(4);
    fcn_TMobil = @(d,l) fcn_1DEfficiency(d,l,fcn_Mobil, mobilLo, mobilHi);

    %   Define separation efficiency function based on region of interest
    fcn_Efficiency = fcn_ComputeRegionEfficiency(region, fcn_TSed, fcn_TMobil);

    A = fcn_ComputeCutArea(fcn_qfeed, fcn_Efficiency);
end

function fcn_MakeSeparationPlot(t, fcn_PSD,fcn_Cut1, cut1_lims, fcn_Cut2, cut2_lims, dlims, llims, plotdlims, plotllims, psd_levels, title_txt, show_spread, plot_legend, fig_size)
    ax1 = nexttile(fig_size);
    box on
    ax2 = copyobj(ax1,t);
    ax3 = copyobj(ax1,t);
    ax4 = copyobj(ax1,t);
    ax5 = copyobj(ax1,t);
    hold on
    mymap = [1,1,1;parula];
    colormap(ax1,mymap)
    fcontour(ax1,fcn_PSD,...
        [dlims, llims],...
        'Fill','on',...
        'MeshDensity',200,...
        'LevelList',psd_levels,...
        'HandleVisibility','off')
    
    fcontour(ax2,fcn_Cut1,...
        [plotdlims, plotllims],...
        'LevelList',mean(cut1_lims),...
        'LineStyle','--',...
        'LineColor',[0.6350 0.0780 0.1840],...
        'LineWidth',2)
    if show_spread
    % fcontour(ax3,fcn_Cut1,...
    %     [plotdlims, plotllims],...
    %     'LevelList',cut1_lims,...
    %     'LineStyle','-',...
    %     'LineColor',[0.6350 0.0780 0.1840],...
    %     'LineWidth',2)
        fc1 = fcontour(ax3,fcn_Cut1,...
                [plotdlims, plotllims],...
                'LevelList',cut1_lims(1),...
                'LineStyle','-',...
                'LineColor',[0.6350 0.0780 0.1840],...
                'LineWidth',2,...
                'Fill','off');
            contourLo = fc1.ContourMatrix;
            x1 = contourLo(1,3:end-2);
            y1 = contourLo(2,3:end-2);

            fc2 = fcontour(ax3,fcn_Cut1,...
                [plotdlims, plotllims],...
                'LevelList',cut1_lims(2),...
                'LineStyle','-',...
                'LineColor',[0.6350 0.0780 0.1840],...
                'LineWidth',2,...
                'Fill','off');
            
            contourHi = fc2.ContourMatrix;
            interp = griddedInterpolant(contourHi(1,3:end-2), contourHi(2,3:end-2),...
                'linear','linear');

            shade(ax3, x1, y1,...
                x1, interp(x1),...
                'FillType',[1 2;2 1],...
                'FillColor',[0.6350 0.0780 0.1840],...
                'FillAlpha',0.2,...
                'LineStyle','none')
    end
    fcontour(ax4,fcn_Cut2,...
        [plotdlims, plotllims],...
        'LevelList',mean(cut2_lims),...
        'LineStyle','--',...
        'LineColor',[0.4940 0.1840 0.5560],...
        'LineWidth',2)
    if show_spread
        fc1 = fcontour(ax5,fcn_Cut2,...
            [plotdlims, plotllims],...
            'LevelList',cut2_lims(1),...
            'LineStyle','-',...
            'LineColor',[0.4940 0.1840 0.5560],...
            'LineWidth',2,...
            'Fill','off');
        contourLo = fc1.ContourMatrix;
        x1 = contourLo(1,3:end-2);
        y1 = contourLo(2,3:end-2);
        fc2 = fcontour(ax5,fcn_Cut2,...
            [plotdlims, plotllims],...
            'LevelList',cut2_lims(2),...
            'LineStyle','-',...
            'LineColor',[0.4940 0.1840 0.5560],...
            'LineWidth',2,...
            'Fill','off');
        
        contourHi = fc2.ContourMatrix;
        interp = griddedInterpolant(contourHi(1,3:end-2), contourHi(2,3:end-2),...
            'linear','linear');
        shade(ax5, x1, y1,...
            x1, interp(x1),...
            'FillType',[1 2;2 1],...
            'FillColor',[0.4940 0.1840 0.5560],...
            'FillAlpha',0.2,...
            'LineStyle','none')
    end

    if plot_legend
        plot(ax1,NaN,NaN,...,
            '--',...,
            'Color',[0.6350 0.0780 0.1840],...
            'HandleVisibility','on',...
            'DisplayName',"Sedimentation",...
            'LineWidth',2)
        plot(ax1,NaN,NaN,...,
            '--',...,
            'Color',[0.4940 0.1840 0.5560],...
            'HandleVisibility','on',...
            'DisplayName',"Mobility",...
            'LineWidth',2)
        legend(ax1,'Location','northwest')
    end

    xlim(plotdlims)
    xticks(linspace(plotdlims(1),plotdlims(2),5))
    ylim(plotllims)
    yticks(linspace(plotllims(1),plotllims(2),5))
    colorbar(ax1)
    xlabel("Diameter / nm")
    ylabel("Length / nm")
    title(title_txt)
    
    ax2.UserData = linkprop([ax1,ax2,ax3,ax4,ax5],...
        {'Position','InnerPosition','DataAspectRatio','xtick','ytick', ...
        'ydir','xdir','xlim','ylim'}); % add more props as needed
    ax2.Visible = 'off';
    ax3.UserData = linkprop([ax1,ax2,ax3,ax4,ax5],...
        {'Position','InnerPosition','DataAspectRatio','xtick','ytick', ...
        'ydir','xdir','xlim','ylim'}); % add more props as needed
    ax3.Visible = 'off';
    ax4.UserData = linkprop([ax1,ax2,ax3,ax4,ax5],...
        {'Position','InnerPosition','DataAspectRatio','xtick','ytick', ...
        'ydir','xdir','xlim','ylim'}); % add more props as needed
    ax4.Visible = 'off';
    ax5.UserData = linkprop([ax1,ax2,ax3,ax4,ax5],...
        {'Position','InnerPosition','DataAspectRatio','xtick','ytick', ...
        'ydir','xdir','xlim','ylim'}); % add more props as needed
    ax5.Visible = 'off';
    hold off
end

function Cost = fcn_MaximizeAreaCostFcn(prm, prm_mag, fcn_qfeed, region, rhoParticle, rhoSolvent, viscositySolvent, surfDensityCharge, desired_area)
    sed_prm = prm(1:2)*prm_mag(1);
    mobil_prm = prm(3:4)*prm_mag(2);
    prm_unnorm = [sed_prm(:);mobil_prm(:)];
    A = fcn_ComputeRegionArea(prm_unnorm, fcn_qfeed, region, rhoParticle, rhoSolvent, viscositySolvent, surfDensityCharge);
    Cost = (A-desired_area).^2;
end
