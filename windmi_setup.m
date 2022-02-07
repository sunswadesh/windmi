%% Load all variables needed for the windmi model to run.
R_Earth = 6380e3;  % radius if earth in meters

if nom ==1 
    % Nominal values, used for testing, taken from IEEE paper 2004(Horton, Doxas )
    L = 100;
    M = 2/3;
    C = 48000;
    Sigma = 7.8;
    Omega = 10000*R_Earth^3;
    mu0 = 4.2e-09;
    Ic = 1.78e+8;
    alpha = 7.550246e+11;
    taup = 10*60;
    tauE = 0.5*3600;
    L1 = 19;
    C1 = 800;
    SigmaI = 3;
    beta = 1;
    alpha_F = 100;
    R_prc = 0.1;
    tau_rc = 12*3600;
    L2 = 8;
    R_A2 =0.3;
    B_tr = 5e-9;
    A_eff = 2*(R_Earth^2);
    L_y = 5*R_Earth;

    Alph = 4.3;
%-------------------------
else
    L = WindmiCoeff(1); M = WindmiCoeff(2); C = WindmiCoeff(3); Sigma = WindmiCoeff(4);
    Omega = WindmiCoeff(5); mu0 = WindmiCoeff(6); Ic = WindmiCoeff(7); alpha = WindmiCoeff(8);
    taup = WindmiCoeff(9); tauE = WindmiCoeff(10); L1 = WindmiCoeff(11); C1 = WindmiCoeff(12);
    SigmaI = WindmiCoeff(13); beta = WindmiCoeff(14); alpha_F = WindmiCoeff(15); R_prc = WindmiCoeff(16); 
    tau_rc = WindmiCoeff(17); L2 = WindmiCoeff(18); R_A2 =WindmiCoeff(19); 
    B_tr = WindmiCoeff(20); A_eff = WindmiCoeff(21); L_y = WindmiCoeff(22); Alph = WindmiCoeff(23);
end
%%

Ic = 3.0e7;

% tauE = 5*3600/10;
% C=1e5;

Icsp = ones(1,length(datrng)); Icsp(900:end)=0.1*sin(0:length(Icsp)-900);
% assignin('base','Icsp',Icsp);

% fixed values
% del_I = 1.250000e+05;

del_I = 1.250000e+05;


A = (L1*L-M^2)^(-1);
inv_dI = (del_I)^(-1);
B_E = 3.1e-5;
pf = 3500; 
dst_f = -(2*1e-7)*(1/(B_E*(R_Earth^3)));


if opt ==1
    % Select between nominal values +/- a tolerance(R_factor), or expected
    % ranges.% Set up parameter ranges

    R_factor = 0.5;

    L_range = [L*(1-1e-0*R_factor) L*(1+1e-0*R_factor)];
    M_range = [M*(1-1e-0*R_factor) M*(1+1e-0*R_factor)];
    C_range = [C*(1-1e-0*R_factor) C*(1+1e-0*R_factor)];
    Sigma_range = [Sigma*(1-1e-0*R_factor) Sigma*(1+1e-0*R_factor)];
    Omega_range = [Omega*(1-1e-0*R_factor) Omega*(1+1e-0*R_factor)];
    mu0_range = [mu0*(1-1e-0*R_factor) mu0*(1+1e-0*R_factor)];
    Ic_range = [Ic*(1-1e-0*R_factor) Ic*(1+1e-0*R_factor)];
    alpha_range = [alpha*(1-1e-0*R_factor) alpha*(1+1e-0*R_factor)];
    taup_range = [taup*(1-1e-0*R_factor) taup*(1+1e-0*R_factor)];
    tauE_range = [tauE*(1-1e-0*R_factor) tauE*(1+1e-0*R_factor)];
    L1_range = [L1*(1-1e-0*R_factor) L1*(1+1e-0*R_factor)];
    C1_range = [C1*(1-1e-0*R_factor) C1*(1+1e-0*R_factor)];
    SigmaI_range = [SigmaI*(1-1e-0*R_factor) SigmaI*(1+1e-0*R_factor)];
    beta_range = [beta*(1-0.00000001*R_factor) beta*(1+0.00000001*R_factor)];
    alpha_F_range = [alpha_F*(1-1e-0*R_factor) alpha_F*(1+1e-0*R_factor)];
    R_prc_range = [R_prc*(1-1.0e-0*R_factor) R_prc*(1+1e-0*R_factor)];

    tau_rc_range = [tau_rc*(1-1.e-0*R_factor) tau_rc*(1+8e-0*R_factor)];
    % tau_rc_range = [6*3600 60*3600];


    L2_range = [L2*(1-1e-0*R_factor) L2*(1+1e-0*R_factor)];
    R_A2_range = [R_A2*(1-1e-0*R_factor) R_A2*(1+1e-0*R_factor)];
    B_tr_range = [B_tr*(1-1.0e-0*R_factor) B_tr*(1+1e-0*R_factor)];
    A_eff_range = [A_eff*(1-1.0e-0*R_factor) A_eff*(1+1e-0*R_factor)];
    L_y_range = [L_y*(1-1.0e-0*R_factor) L_y*(1+1e-0*R_factor)];

    Alph_range = [Alph*(1-1e-0*R_factor) Alph*(1+1e-0*R_factor)];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%Set Ranges for all the variables;

    varlim(1,:) = L_range;varlim(2,:) = M_range;varlim(3,:) = C_range;varlim(4,:) = Sigma_range;
    varlim(5,:) = Omega_range;varlim(6,:) = mu0_range;varlim(7,:) = Ic_range;varlim(8,:) = alpha_range;
    varlim(9,:) = taup_range;varlim(10,:) = tauE_range;varlim(11,:) = L1_range;varlim(12,:) = C1_range;
    varlim(13,:) = SigmaI_range;varlim(14,:) = beta_range;varlim(15,:) = alpha_F_range;
    varlim(16,:) = R_prc_range;varlim(17,:) = tau_rc_range;varlim(18,:) = L2_range;varlim(19,:) = R_A2_range;
    varlim(20,:) = B_tr_range;varlim(21,:) = A_eff_range;varlim(22,:) = L_y_range; varlim(23,:) = Alph_range;

    N_var=size(varlim,1); %%Number variables GA has to optimize.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
