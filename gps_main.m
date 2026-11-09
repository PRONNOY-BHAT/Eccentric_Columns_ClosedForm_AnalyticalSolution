clear; clc; close all;

%% Input Parameters from the Source Lin and Teng (2019)
par.E_c     = 30e9;            % concrete Young's modulus       [Pa]
par.nu_c    = 0.18;            % concrete Poisson ratio         [-]
par.R       = 0.075;           % column radius                  [m]
par.t_f     = 0.165e-3;        % FRP thickness                  [m]
par.E_f     = 250e9;           % FRP modulus                    [Pa]
par.eps0    = -1.0e-3;         % centroidal axial strain        [-]
factor      = 0.09375; 
par.e       = factor*par.R;
par.kappa   = (4*abs(par.eps0)*par.e)/(par.R^2);


%% FE benchmark results 
FE.q0       = 97.34e3;         % axisymmetric component (n=0) from FE  [Pa]
FE.q1_tilde = -35.09e3;        % first-harmonic compoenent (n=1) from FE  [Pa]


%% Solve closed-form
sol = gps_solver(par);


%% Plotting setup 
theta = linspace(0, 2*pi, 721);  
deg   = rad2deg(theta);


%%  Figure 2: Confinement pressure around circumference (line + polar)
fl_full = sol.q0 + sol.q1_tilde*cos(theta);    % Combined n=0 + n=1
fl_axi  = sol.q0 * ones(size(theta));          % axisymmetric only
fl_onlymode1 = sol.q1_tilde*cos(theta); 

Matrix1(:,1) = deg;
Matrix1(:,2) = fl_full;
Matrix1(:,3) = fl_onlymode1;
Matrix1(:,4) = fl_axi;
save('ConfinementPressureDistribution.txt','Matrix1','-ascii')

figure('Color','w','Position',[80 80 1200 460],'Name','Fig 2: Confinement');
subplot(1,2,1);
plot(deg, fl_axi/1e3, 'k--', 'LineWidth', 1.6); 
hold on;
plot(deg, fl_full/1e3, 'b-',  'LineWidth', 2.0);
xlabel('Angular position \theta (deg)');
ylabel('Confinement pressure f_l(\theta)  [kPa]');
title('Confinement pressure around the circumference');
legend(sprintf('Axisymmetric only: q^{(0)} = %.1f kPa', sol.q0/1e3),'n = 0 + n = 1 combined','Location','south');
xlim([0 360]); xticks(0:45:360); grid on;

subplot(1,2,2);
polarplot(theta, fl_full/1e3, 'b-', 'LineWidth', 2.0); 
hold on;
polarplot(theta, fl_axi/1e3,  'k--', 'LineWidth', 1.6);
title('Polar view: f_l(\theta) in kPa');
legend('n=0 + n=1','Axisymmetric only','Location','southoutside');

sgtitle('Fig. 2  –  Confinement pressure f_l(\theta) at FRP-concrete interface', ...
        'FontWeight','bold','FontSize',13);


%%  Figure 3: Radial stress at r = R around the circumference
sigma_r_R = -fl_full;
figure('Color','w','Position',[100 100 900 460],'Name','Fig 3: sigma_r at R');
plot(deg, -sol.q0/1e3*ones(size(deg)), 'k--', 'LineWidth', 1.6); hold on;
plot(deg, sigma_r_R/1e3, 'b-', 'LineWidth', 2.0);
xlabel('Angular position \theta (deg)');
ylabel('\sigma_r(R,\theta)  [kPa, tension positive]');
title('Fig. 3  –  Radial stress at FRP interface');
legend('Axisymmetric only','n=0 + n=1 combined','Location','best');
xlim([0 360]); xticks(0:45:360); grid on;


%%  Figure 4: Hoop stress at r = R around the circumference
sigma_th_R = -sol.q0 + 6*sol.A*par.R*cos(theta);
figure('Color','w','Position',[120 120 900 460],'Name','Fig 4: sigma_th at R');
plot(deg, -sol.q0/1e3*ones(size(deg)), 'k--', 'LineWidth', 1.6); hold on;
plot(deg, sigma_th_R/1e3, 'r-', 'LineWidth', 2.0);
xlabel('Angular position \theta (deg)');
ylabel('\sigma_\theta(R,\theta)  [kPa]');
title('Fig. 4  –  Hoop stress at FRP interface  (n=1 amplitude is 3× radial)');
legend('Axisymmetric only','n=0 + n=1 combined','Location','best');
xlim([0 360]); xticks(0:45:360); grid on;


%%  Figure 5: Shear stress at r = R around the circumference
Tau_Rtheta = 2*sol.A*par.R*sin(theta);
Matrix2(:,1) = deg;
Matrix2(:,2) = sigma_r_R;
Matrix2(:,3) = sigma_th_R;
Matrix2(:,4) = Tau_Rtheta;
Matrix2(:,5) = -sol.q0*ones(size(deg))
save('SigmaRSigmaThetaDistribution.txt','Matrix2','-ascii')


%%  Figures (6-8): Stress fields inside concrete cross-section
N       = 200;
xv      = linspace(-par.R, par.R, N);
yv      = linspace(-par.R, par.R, N);
[X,Y]   = meshgrid(xv, yv);
[TH,RR] = cart2pol(X,Y);
inside  = RR <= par.R;

f       = gps_field(RR, TH, par, sol);
% Mask outside circular concrete section
sr  = f.sigma_r;      sr(~inside)  = NaN;
sth = f.sigma_th;     sth(~inside) = NaN;
trt = f.tau_rth;      trt(~inside) = NaN;

% Common plotting settings
xplot = xv*1e3;
yplot = yv*1e3;
ticks = -75:25:75;

% Combined Figure: sigma_r, sigma_theta and tau_rtheta
figure('Color','w','Position',[100 150 1500 500],'Name','Stress fields inside concrete');

% (a) Radial stress sigma_r
subplot(1,3,1);
h = imagesc(xplot, yplot, sr/1e3);
set(h,'AlphaData',inside);
set(gca,'YDir','normal');
set(gca,'Color','w');
axis image;
xticks(ticks); yticks(ticks);
hold on;
plotCircle(par.R*1e3);
colormap(redblue(256));
clim([-140 0]);
colorbar;
title('\sigma_r(r,\theta)  [kPa]');
xlabel('x (mm)');
ylabel('y (mm)');

% (b) Hoop stress sigma_theta
subplot(1,3,2);
h = imagesc(xplot, yplot, sth/1e3);
set(h,'AlphaData',inside);
set(gca,'YDir','normal');
set(gca,'Color','w');
axis image;
xticks(ticks);
yticks(ticks);
hold on;
plotCircle(par.R*1e3);
colormap(redblue(256));
clim([-250 0]);
colorbar;
title('\sigma_\theta(r,\theta)  [kPa]');
xlabel('x (mm)');
ylabel('y (mm)');

% (c) Shear stress tau_rtheta
subplot(1,3,3);
h = imagesc(xplot, yplot, trt/1e3);
set(h,'AlphaData',inside);
set(gca,'YDir','normal');
set(gca,'Color','w');
axis image;
xticks(ticks); yticks(ticks);
hold on;
plotCircle(par.R*1e3);
colormap(redblue(256));
clim([-40 40]);
colorbar;
title('\tau_{r\theta}(r,\theta)  [kPa]');
xlabel('x (mm)');
ylabel('y (mm)');

% Overall figure title
sgtitle('Stress fields inside the concrete cross-section','FontWeight','bold');
exportgraphics(gcf,'Fig5_StressFields.tif','Resolution',1200);



% ---------- Fig 8: sigma_z field + FRP hoop strain ----------
sz  = f.sigma_z;     sz(~inside)  = NaN;
figure('Color','w','Position',[200 200 1200 460],'Name','Fig 8: sigma_z & eh');
subplot(1,2,1);
h = imagesc(xv*1e3, yv*1e3, sz/1e3); set(h,'AlphaData',inside); set(gca,'Color','w'); 
set(gca,'YDir','normal'); axis image;
hold on; plotCircle(par.R*1e3); colormap(redblue(256)); colorbar;
title('\sigma_z(r,\theta)  [MPa]');
xlabel('x (mm)'); ylabel('y (mm)');
XYZ_sigma_z = [X(:)*1e3, Y(:)*1e3, sz(:)/1e6];
writematrix(XYZ_sigma_z, 'AxialStressDistribution.txt', 'Delimiter', 'tab');

eh_full = sol.eh_axisym - sol.eh_bend*cos(theta);
subplot(1,2,2);
plot(deg, sol.eh_axisym*ones(size(deg))*1e4, 'k--', 'LineWidth', 1.6); hold on;
plot(deg, eh_full*1e4, 'b-', 'LineWidth', 2.0);
plot([180 180], [0 sol.eh_max*1e4], 'r:', 'LineWidth', 1.0);
plot(180, sol.eh_max*1e4, 'r.', 'MarkerSize', 24);
text(185, sol.eh_max*1e4, sprintf(' ECF: \\epsilon_h^{max} = %.2f\\times10^{-4}', sol.eh_max*1e4),'Color','r','FontSize',10);
xlabel('Angular position \theta (deg)');
ylabel('FRP hoop strain \epsilon_h  (\times 10^{-4})');
title('FRP hoop strain around circumference');
legend(sprintf('Axisymmetric: %.2f\\times10^{-4}', sol.eh_axisym*1e4),'n=0 + n=1 combined','Location','best');
xlim([0 360]); xticks(0:45:360); grid on;
sgtitle('Fig. 8  –  Axial stress \sigma_z and FRP hoop strain','FontWeight','bold');
Matrix3(:,1) = deg;
Matrix3(:,2) = eh_full;
save('FRP hoop strain Distribution.txt','Matrix3','-ascii')


%% =====================================================================
%  Local helper functions
function plotCircle(R)
    th = linspace(0, 2*pi, 400);
    plot(R*cos(th), R*sin(th), 'k', 'LineWidth', 1.6);
end

function c = redblue(n)
    half = floor(n/2);
    R1 = linspace(0.20,1,half)'; G1 = linspace(0.30,1,half)'; B1 = linspace(0.70,1,half)';
    R2 = linspace(1,0.75,n-half)'; G2 = linspace(1,0.10,n-half)'; B2 = linspace(1,0.15,n-half)';
    c = [R1 G1 B1; R2 G2 B2];
end