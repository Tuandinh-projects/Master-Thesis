%% ========================================================================
%  Behavioral TANK Model for Small Open Economy in a Currency Union:
%  (1) Spending Cuts (-1% g)  &  (2) Lump-sum Tax Hike (+1% tau)
%  BEHAVIORAL SOE-TANK (m = 0.85) VS. RATIONAL SOE-TANK (m = 1.00)
%% ========================================================================
clear; close all; clc;

% 1. Parameter 
beta    = 0.99;   % Discount factor (\beta)
gamma   = 5.0;    % Risk aversion (\gamma)
varphi  = 1.0;    % Inverse Frisch elasticity (\varphi)
theta   = 0.875;  % Calvo price stickiness (\theta)
epsilon = 5.3;    % Elasticity of substitution among domestic varieties (\epsilon)
lam     = 0.20;   % Share of Hand-to-Mouth households (\lambda)
alpha   = 0.30;   % Import share / openness (\alpha)
etaX    = 1.50;   % Export price elasticity (\eta_X)
psi     = 0.01;   % Debt-elastic risk premium (\psi)
phib    = 0.33;   % Tax feedback on public debt (\phi_b)
phig    = 0.10;   % Tax feedback on gov spending (\phi_g)
rhog    = 0.80;   % Persistence of gov spending shock (\rho_g)
rhoa    = 0.90;   % Persistence of TFP shock (\rho_a)
rhoiea  = 0.80;   % Persistence of Euro Area interest rate shock (\rho_{iea})
rhoy    = 0.80;   % Persistence of foreign output shock (\rho_{y*})
m       = 0.85;   % Cognitive discounting parameter (m)

% 2. Steady State Values
Wbar  = (epsilon - 1) / epsilon;    % \overline{W} = \overline{MC}
Gbar  = 0.28;                       % \overline{G}
Cbar  = 1 - Gbar;                   % \overline{C}
bbar  = 4 * 0.60;                   % \overline{b} (60% GDP)
Xbar  = alpha * Cbar;               % \overline{X}
Tbar  = Gbar + (1 - beta) * bbar;   % \overline{T}
Dbar  = 1 - Wbar;                   % \overline{D} = 1/\epsilon

netF   = (1 - beta)*(bbar/(1 - lam)) + Dbar/(1 - lam) - Tbar;
nh_min = (Tbar / Wbar) + 1e-4;
nf_max = (1 - lam * nh_min) / (1 - lam);
nf_min = max(0.01, (-netF / Wbar) + 1e-4);

ss_obj = @(nf) nf / ((1 - (1 - lam)*nf)/lam) - ( (Wbar*((1 - (1 - lam)*nf)/lam) - Tbar) / (Wbar*nf + netF) )^(gamma/varphi);
nfbar  = fzero(ss_obj, [nf_min, nf_max]);
nhbar  = (1 - (1 - lam)*nfbar) / lam;
cfbar  = Wbar*nfbar + netF;
chbar  = Wbar*nhbar - Tbar;

Lfbar  = Wbar * nfbar;              % \overline{L}^F = \overline{W}\overline{N}^F
afbar  = bbar / (1 - lam);          % \overline{a}^F = \overline{b}/(1-\lambda)
KF     = cfbar + (gamma / varphi) * Lfbar; % K^F = \overline{C}^F + (\gamma/\varphi)\overline{L}^F
kappa  = (1 - theta) * (1 - beta * theta) / theta; % \kappa_p

% 3. Create b_soe_compare.mod
fid = fopen('b_soe_compare.mod', 'w');
fprintf(fid, 'var cf ch nf nh c n y a w mc p pid pic x g b f tau deltaD i iea ystar;\n');
fprintf(fid, 'varexo eps_g eps_tau eps_a eps_iea eps_ystar;\n');
fprintf(fid, 'parameters beta gamma varphi theta epsilon lam alpha etaX psi phib phig rhog rhoa rhoiea rhoy Wbar Cbar cfbar chbar nfbar nhbar Gbar Xbar Tbar bbar Dbar Lfbar afbar KF kappa m;\n\n');

P = {'beta',beta; 'gamma',gamma; 'varphi',varphi; 'theta',theta; 'epsilon',epsilon; ...
     'lam',lam; 'alpha',alpha; 'etaX',etaX; 'psi',psi; 'phib',phib; 'phig',phig; ...
     'rhog',rhog; 'rhoa',rhoa; 'rhoiea',rhoiea; 'rhoy',rhoy; 'm',m; ...
     'Wbar',Wbar; 'Gbar',Gbar; 'Cbar',Cbar; 'bbar',bbar; 'Xbar',Xbar; 'Tbar',Tbar; 'Dbar',Dbar; ...
     'cfbar',cfbar; 'chbar',chbar; 'nfbar',nfbar; 'nhbar',nhbar; 'Lfbar',Lfbar; 'afbar',afbar; ...
     'KF',KF; 'kappa',kappa};

for k = 1:size(P,1)
    fprintf(fid, '%s = %.16g;\n', P{k,1}, P{k,2});
end

fprintf(fid, '\nmodel(linear);\n');
fprintf(fid, '  w = gamma*cf + varphi*nf - alpha*p;\n');                                % (L1)
fprintf(fid, '  w = gamma*ch + varphi*nh - alpha*p;\n');                                % (L2)
fprintf(fid, '  chbar*ch = Wbar*nhbar*(w + nh + alpha*p) - tau;\n');                    % (L3)
fprintf(fid, '  Cbar*c = (1 - lam)*cfbar*cf + lam*chbar*ch;\n');                        % (L4)
fprintf(fid, '  n = (1 - lam)*nfbar*nf + lam*nhbar*nh;\n');                             % (L5)
fprintf(fid, '  y = a + n;\n');                                                         % (L6)
fprintf(fid, '  mc = w - a;\n');                                                        % (L7)
fprintf(fid, '  pid = beta*pid(+1) + kappa*mc;\n');                                     % (L8)
fprintf(fid, '  pid = p - p(-1);\n');                                                   % (L9)
fprintf(fid, '  pic = (1 - alpha)*pid;\n');                                             % (L10)
fprintf(fid, '  deltaD = y - Wbar*(w + n);\n');                                         % (L11)
fprintf(fid, '  x = ystar - etaX*p;\n');                                                % (L12)
fprintf(fid, '  y = (1 - alpha)*Cbar*(c - alpha*p) + g + Xbar*x;\n');                   % (L13)
fprintf(fid, '  f = (f(-1) + Xbar*(x + alpha*p) - alpha*Cbar*(c + (1 - alpha)*p))/beta;\n'); % (L14)
fprintf(fid, '  i = iea - psi*f(-1);\n');                                               % (L15)
fprintf(fid, '  b = (b(-1) - bbar*pic + g + alpha*Gbar*p - tau)/beta + bbar*i;\n');     % (L16)
fprintf(fid, '  tau = phib*b(-1) + phig*g + eps_tau;\n');                               % (L17)
fprintf(fid, '  cf = m*cf(+1) - (i - m*pic(+1))/gamma + (1 - beta)*(1 - m)/KF * ( (b + f)/(1 - lam) + afbar*pic(+1) );\n'); % (B10)
fprintf(fid, '  g = rhog*g(-1) + eps_g;\n');
fprintf(fid, '  a = rhoa*a(-1) + eps_a;\n');
fprintf(fid, '  iea = rhoiea*iea(-1) + eps_iea;\n');
fprintf(fid, '  ystar = rhoy*ystar(-1) + eps_ystar;\n');
fprintf(fid, 'end;\n\n');
fprintf(fid, 'steady; check;\n');
fprintf(fid, 'shocks; var eps_g; stderr 0.01; var eps_tau; stderr 0.01; var eps_a; stderr 0.01; var eps_iea; stderr 0.001; var eps_ystar; stderr 0.01; end;\n');
fprintf(fid, 'stoch_simul(order=1, irf=24, nograph) y c cf ch nf nh pid pic p b f tau g x i;\n');
fclose(fid);

%% ========================================================================
%  4. Simulating 2 Models: BEHAVIORAL (m = 0.85) VS. RATIONAL (m = 1.00)
%% ========================================================================

% Stage A: behavioral baseline model (m = 0.85)
dynare b_soe_compare noclearall;
irf_behav = oo_.irfs;

% State B: update m = 1.00 (Rational SOE-TANK) and re-run
set_param_value('m', 1.00);
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);
irf_rat = oo_.irfs;

%% ========================================================================
%  5. Create IRFs and Comparing Figures
%% ========================================================================
H = 24;
t = 1:H;

vars_all   = {'y', 'c', 'cf', 'ch', 'g', ...
    'nf', 'nh', 'pid', 'pic', 'p', ...
    'b', 'f', 'tau', 'x', 'i'};

titles_all = {'y (Output)', 'c (Total Cons.)', 'c^F (Forward-looking Cons.)', ...
    'c^H (HtM Cons.)', 'g (Gov. Spending)', ...
    'n^F (Forward-looking Labor)', 'n^H (HtM Labor)', ...
    '\pi^D (Dom. Inflation)', '\pi^C (CPI Inflation)', 'p (Relative Price)', ...
    'b (Public Debt)', 'f (Net Foreign Assets)', '\tau (Lump-sum Tax)', ...
    'x (Exports)', 'i (Nominal Interest Rate)'};

% --- FIGURE 1: Spending Cuts (-1% g) ---
figure('Name', 'Fig 1: All IRFs to -1% Spending Cut', 'Color', 'w', 'Position', [40 80 1400 720]);
for k = 1:length(vars_all)
    subplot(3, 5, k);
    plot(t, -irf_behav.([vars_all{k} '_eps_g']), 'k-',  'LineWidth', 1.6); hold on;
    plot(t, -irf_rat.([vars_all{k} '_eps_g']),   'b--', 'LineWidth', 1.6);
    yline(0, 'r-', 'LineWidth', 0.8);
    title(titles_all{k}, 'FontWeight', 'bold');
    xlim([1 H]); grid on;
    if k == 1
        legend({'Behavioral (m = 0.85)', 'Rational (m = 1.00)'}, 'Location', 'SouthEast', 'Box', 'off');
    end
end

% --- FIGURE 2: Tax Hike (+1% tau) ---
figure('Name', 'Fig 2: All IRFs to +1% Tax Hike', 'Color', 'w', 'Position', [70 50 1400 720]);
for k = 1:length(vars_all)
    subplot(3, 5, k);
    plot(t, irf_behav.([vars_all{k} '_eps_tau']), 'k-',  'LineWidth', 1.6); hold on;
    plot(t, irf_rat.([vars_all{k} '_eps_tau']),   'b--', 'LineWidth', 1.6);
    yline(0, 'r-', 'LineWidth', 0.8);
    title(titles_all{k}, 'FontWeight', 'bold');
    xlim([1 H]); grid on;
    if strcmp(vars_all{k}, 'g')
        ylim([-0.01 0.01]);
    end
    if k == 1
        legend({'Behavioral (m = 0.85)', 'Rational (m = 1.00)'}, 'Location', 'Best', 'Box', 'off');
    end
end