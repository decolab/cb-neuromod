function FCD_LSD_model

%%%%%%%%%%%%%%%%%%%%%%%%
%   
%  Computes simulations with the Dynamic Mean Field Model (DMF) using 
%  Feedback Inhibitory Control (FIC) and Regional Drug Receptor Modulation (RDRM):
%
%  - the optimal coupling (we=2.1) for fitting the placebo condition 
%  - the optimal neuromodulator gain for fitting the LSD condition (wge=0.2)
%
%
%   Whole-brain multimodal neuroimaging model using serotonin receptor maps explain non-linear functional effects of LSD
%   Deco,G., Cruzat,J., Cabral, J., Knudsen,G.M., Carhart-Harris,R.L., Whybrow,P.C., 
%       Logothetis,N.K. & Kringelbach,M.L. (2018) Current Biology
%
%  Code written by Gustavo Deco gustavo.deco@upf.edu 2017
%  Reviewed by Josephine Cruzat and Joana Cabral 
%  
%%%%%%%%%%%%%%%%%%%%%%%%

% Set global variables that are used in sub-functions
global wgaine wgaini Receptor;

% Load Structural Connectivity Matrix
load all_SC_FC_TC_76_90_116.mat sc90
C=sc90/max(sc90(:))*0.2;
N=size(C,1);

% Load Regional Drug Receptor Map
load mean5HT2A_bindingaal.mat mean5HT2A_aalsymm
Receptor=mean5HT2A_aalsymm(:,1)/max(mean5HT2A_aalsymm(:,1));

% Set General Model Parameters
dtt   = 1e-3;  % Sampling rate of simulated neuronal activity (seconds)
dt    = 0.1;    
taon=100;      
taog=10;
gamma=0.641;
sigma=0.01;
JN=0.15;
I0=0.382;
Jexte=1.;
Jexti=0.7;
w=1.4;

TR    = 2;     % Sampling rate of saved simulated BOLD (seconds)
NSUB=15;       % Number of Subjects in empirical fMRI dataset 
Tmax=220;      % Number of timepoints in each fMRI session
Tmaxneuronal=(Tmax+10)*(TR/dtt); % Numer of simulated time points

% FILTER SETTINGS
fnq=1/(2*TR);                 % Nyquist frequency
flp = .02;                    % lowpass frequency of filter
fhi = 0.1;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter

Isubdiag = find(tril(ones(N),-1)); % Indices of upper triagular part of FC matrix

% Define optimal parameters 
we=2.1;  % Global Coupling parameter
J=Balance_J9(we,C); % This is the Feedback Inhibitory Control
% I is calculated this only once, then saved
save J_Balance J
 
% J can be calculated only once and then load J_Balance J

%% SIMULATION OF OPTIMAL PLACEBO
wge=0; % 0 for placebo, 0.2 for LSD
wgaini=0;
wgaine=wge;

N_windows=length(1:3:190);
cotsampling_pla_s=zeros(NSUB,N_windows*(N_windows-1)/2);

[status,cmdout] = system('od -vAn -N4 -tu4 < /dev/urandom');
rng(str2num(cmdout));

for nsub=1:NSUB
    kk=1;
    disp(['Subject ' num2str(nsub)])
    neuro_act=zeros(Tmaxneuronal,N);
    sn=0.001*ones(N,1);
    sg=0.001*ones(N,1);
    nn=1;
    % Model simulations
    for t=0:dt:Tmaxneuronal-1
        xn=I0*Jexte+w*JN*sn+we*JN*C*sn-J.*sg;
        xg=I0*Jexti+JN*sn-sg;
        rn=phie9(xn);
        rg=phii9(xg);
        sn=sn+dt*(-sn/taon+(1-sn)*gamma.*rn./1000.)+sqrt(dt)*sigma*randn(N,1);
        sn(sn>1) = 1;
        sn(sn<0) = 0;
        sg=sg+dt*(-sg/taog+rg./1000.)+sqrt(dt)*sigma*randn(N,1);
        sg(sg>1) = 1;
        sg(sg<0) = 0;
        if abs(mod(t,1))<0.01
            neuro_act(nn,:)=rn';
            nn=nn+1;
        end
    end
    nn=nn-1;
    
    %%%% BOLD
    % Friston BALLOON-WINDKESSEL MODEL
    T = nn*dtt; % Total time in seconds    
    B = BOLD(T,neuro_act(:,1)'); % B=BOLD activity, bf=Foutrier transform, f=frequency range)
    BOLD_act = zeros(length(B),N);
    BOLD_act(:,1) = B;
    
    for nnew=2:N
        B = BOLD(T,neuro_act(:,nnew));
        BOLD_act(:,nnew) = B;
    end
    
    bds=BOLD_act(TR/dtt:TR/dtt:end,:);
    
    %%%%%%%%%%%%
    
    % Compute the Hilbert phase from BOLD signals
    for seed=1:N
        ts=detrend(demean(bds(1:Tmax,seed)'));
        ts((ts>3*std(ts)))=3*std(ts);
        ts((ts<-3*std(ts)))=-3*std(ts);
        tss(seed,:) = filtfilt(bfilt,afilt,ts);
    end
    
    ii2=1;
    for t=1:3:190
        jj2=1;
        cc=corrcoef((tss(:,t:t+30))');
        for t2=1:3:190            
            if jj2>ii2
                cc2=corrcoef((tss(:,t2:t2+30))');
                ca=corr2(cc(Isubdiag),cc2(Isubdiag));
                cotsampling_pla_s(nsub,kk)=ca;
                kk=kk+1;
            end
            jj2=jj2+1;
        end
        ii2=ii2+1;
    end   
end

figure; hist(cotsampling_pla_s(:))

save FCD_values_placebo cotsampling_pla_s
        

%% SIMULATION OF OPTIMAL LSD fit
wge=0.2; % 0 for placebo, 0.2 for LSD
wgaini=0;
wgaine=wge;

N_windows=length(1:3:190);
cotsampling_lsd_s=zeros(NSUB,N_windows*(N_windows-1)/2);

[status,cmdout] = system('od -vAn -N4 -tu4 < /dev/urandom');
rng(str2num(cmdout));

for nsub=1:NSUB
    kk=1;
    disp(['Subject ' num2str(nsub)])
    neuro_act=zeros(Tmaxneuronal,N);
    sn=0.001*ones(N,1);
    sg=0.001*ones(N,1);
    nn=1;
    for t=0:dt:Tmaxneuronal-1
        xn=I0*Jexte+w*JN*sn+we*JN*C*sn-J.*sg;
        xg=I0*Jexti+JN*sn-sg;
        rn=phie9(xn);
        rg=phii9(xg);
        sn=sn+dt*(-sn/taon+(1-sn)*gamma.*rn./1000.)+sqrt(dt)*sigma*randn(N,1);
        sn(sn>1) = 1;
        sn(sn<0) = 0;
        sg=sg+dt*(-sg/taog+rg./1000.)+sqrt(dt)*sigma*randn(N,1);
        sg(sg>1) = 1;
        sg(sg<0) = 0;
        if abs(mod(t,1))<0.01
            neuro_act(nn,:)=rn';
            nn=nn+1;
        end
    end
    nn=nn-1;
    
    %%%% BOLD
    % Friston BALLOON-WINDKESSEL MODEL
    T = nn*dtt; % Total time in seconds
    
    B = BOLD(T,neuro_act(:,1)'); % B=BOLD activity, bf=Foutrier transform, f=frequency range)
    BOLD_act = zeros(length(B),N);
    BOLD_act(:,1) = B;
    
    for nnew=2:N
        B = BOLD(T,neuro_act(:,nnew));
        BOLD_act(:,nnew) = B;
    end
    
    bds=BOLD_act(TR/dtt:TR/dtt:end,:);
    
    %%%%%%%%%%%%
    
    for seed=1:N
        ts=detrend(demean(bds(1:Tmax,seed)'));
        ts((ts>3*std(ts)))=3*std(ts);
        ts((ts<-3*std(ts)))=-3*std(ts);
        tss(seed,:) = filtfilt(bfilt,afilt,ts);
    end
    
    ii2=1;
    for t=1:3:190
        jj2=1;
        cc=corrcoef((tss(:,t:t+30))');
        for t2=1:3:190            
            if jj2>ii2
                cc2=corrcoef((tss(:,t2:t2+30))');
                ca=corr2(cc(Isubdiag),cc2(Isubdiag));
                cotsampling_lsd_s(nsub,kk)=ca;
                kk=kk+1;
            end
            jj2=jj2+1;
        end
        ii2=ii2+1;
    end   
end

figure; hist(cotsampling_lsd_s(:))

save FCD_values_lsd cotsampling_lsd_s


[h_lsd_s, x]=hist(cotsampling_lsd_s(:),-.1:.025:1);
h_pla_s=hist(cotsampling_pla_s(:),-.1:.025:1);

%figure

subplot(1,2,2)
bar(x,[h_lsd_s' h_pla_s'])
ylabel('Count')
xlabel('FCD values')
legend('LSD','Placebo')
title('Simulated data')
box off

