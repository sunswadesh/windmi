%% Run the windmi model and extract currents

% t_stamp = stim;

Inputy = ave_filt(VSW(datrng)',1);    %% Choose Vsw only for the storm interval.

%% The time is modified for the windmi model since simulink likes uniformly spaced time points starting from zero.
if (Winddat ==0 )
    % Choose time only for the storm interval (in secs).
    Inputx = linspace(0,vsw_time(datrng(end))/24/3600-stm_int(1),length(Inputy))*24*3600; %
%     Inputyy(ii,:) = Inputy;
%     Inputxx(ii,:) = Inputx;
    Input = [Inputx',Inputy',Icsp'];
else
    sdnan = isnan(Inputy); % remove nan values
    Inputy(sdnan) = [];
%     Icsp(sdnan) = [];
    % make time equally spaced, simulink doesn't like it otherwise
%     Inputx = linspace(min(t_stamp),max(t_stamp), length(Inputy));  
    Inputx = vsw_time(sdnan==0); %real time has gaps
%     Inputx = linspace(min(Inputx), max(Inputx),length(Inputx)); % real time 
%     Inputx = linspace(0,vsw_time(datrng(end))/24/3600,length(Inputy))*24*3600; %
%     Inputyy(ii+1,:) = Inputy;
%     Inputxx(ii+1,:) = Inputx;
    % Simulink needs time starting from 0. Hence modify inputx
    minx = min(Inputx); % save the minimum value from the data with nans
    finx = Inputx-minx; %Fake time starting at 0 but duration is the same as actual time
%     finx = linspace(0, max(Inputx),length(Inputx)); % fake time 
    Input = [finx',Inputy'];
    Inputx = finx;
end

% %% Make ring current time constant tau_rc a function of inputy.
% tau_rc = 

% ref_y = -(ave_filt(AL_data,1));
ref_y = AL_data;
tspan = Inputx(size(Inputx,2));

% sim('windmi_8',tspan);
sim('windmi_8_SP',tspan);

dst_time = Inputx;


                                                                                                                                                                                                                                                                                          YY = yout(:,2)';
y_sim = interp1(tout',YY,dst_time);
dst_sim = dst_f*interp1(tout',yout(:,7),dst_time);
% dst_scale = max(abs(dst_data))/(1e9*max(abs(dst_sim)));
PVterm = interp1(tout',yout(:,11),dst_time);
I = interp1(tout',yout(:,1),dst_time);
I1 = interp1(tout',yout(:,2),dst_time);
I2 = interp1(tout',yout(:,8),dst_time);
Volt = interp1(tout',yout(:,3),dst_time);
VSW_delay = interp1(tout',yout(:,10)',dst_time);
p_cps = interp1(tout',yout(:,5)',dst_time);
I_geotail = interp1(tout',yout(:,1),dst_time);
theta = interp1(tout',yout(:,9),dst_time);
I_ps = interp1(tout',yout(:,12),dst_time);
pVAeff = interp1(tout',yout(:,11),dst_time);


thin = interp1(tout',yout(:,13),dst_time);
icout = interp1(tout',yout(:,14),dst_time);

I2Rprc = R_prc*I2.^2;
% make time the same as vsw_time
if Winddat==1 % make time accurate again.
    dst_time = dst_time + minx;
%     Inputx = Inputx  + Inputxx(ii,1);
end