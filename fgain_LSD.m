clc; clear all;

%%  OPTIMIZATION LSD
%
%   Whole-brain multimodal neuroimaging model using serotonin receptor maps explain non-linear functional effects of LSD
%   Deco,G., Cruzat,J., Cabral, J., Knudsen,G.M., Carhart-Harris,R.L., Whybrow,P.C., 
%       Logothetis,N.K. & Kringelbach,M.L. (2018) Current Biology
%
%   July, 2018, Barcelona-Spain
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% LOAD DATA

load all_SC_FC_TC_76_90_116.mat;
load mean5HT2A_bindingaal.mat
load LSDnew.mat;

global wgaine wgaini Receptor;

C = sc90;                   % Structural connectivity data
C = C/max(max(C))*0.2;
N = 90;                     % Number of brain regions

Receptor = mean5HT2A_aalsymm(:,1)/max(mean5HT2A_aalsymm(:,1));

NSUB     = 15;                   % Number of subjects
Isubdiag = find(tril(ones(N),-1));

%% DATA (There are different cases: 5=PLACEBO, 2=LSD)

CASE = 2;

TC1=LR_version_symm(tc_aal{1,CASE});
TC2=LR_version_symm(tc_aal{2,CASE});
TC3=LR_version_symm(tc_aal{3,CASE});
TC4=LR_version_symm(tc_aal{4,CASE});
TC5=LR_version_symm(tc_aal{5,CASE});
TC6=LR_version_symm(tc_aal{6,CASE});
TC7=LR_version_symm(tc_aal{7,CASE});
TC8=LR_version_symm(tc_aal{8,CASE});
TC9=LR_version_symm(tc_aal{9,CASE});
TC10=LR_version_symm(tc_aal{10,CASE});
TC11=LR_version_symm(tc_aal{11,CASE});
TC12=LR_version_symm(tc_aal{12,CASE});
TC13=LR_version_symm(tc_aal{13,CASE});
TC14=LR_version_symm(tc_aal{14,CASE});
TC15=LR_version_symm(tc_aal{15,CASE});

xs   = eval(sprintf('TC%d',1));
Tmax = size(xs,2);

%%%%%%%%%%%%%%% FILTER SETTINGS  

delt    = 2;                        % sampling interval
k       = 2;                        % 2nd order butterworth filter
fnq     = 1/(2*delt);               % Nyquist frequency
flp     = .01;                      % lowpass frequency of filter
fhi     = .1;                       % highpass
Wn      = [flp/fnq fhi/fnq];        % butterworth bandpass non-dimensional frequency

[bfilt2,afilt2] =  butter(k,Wn);    % construct the filter

%%%%%%%%%%%%%%

kk  = 1;
kk3 = 1;
for nsub = 1:NSUB
    signaldata = eval(sprintf('TC%d', nsub));
    FCe(nsub,:,:) = corrcoef(signaldata');
    
    %Get the BOLD phase using the Hilbert transform
    for seed = 1:N
        x = demean(detrend(signaldata(seed,:)));
        x(find(x>3*std(x)))  = 3*std(x);
        x(find(x<-3*std(x))) = -3*std(x);
        timeseriedata(seed,:)= filtfilt(bfilt2,afilt2,x);    % zero phase filter the data
    end
       
    ii2 = 1;
    for t = 1:3:190
        jj2 = 1;
        cc = corrcoef((timeseriedata(:,t:t+30))');
        for t2 = 1:3:190
            cc2 = corrcoef((timeseriedata(:,t2:t2+30))');
            ca = corrcoef(cc(Isubdiag),cc2(Isubdiag));
            if jj2 >ii2
                cotsampling(kk) = ca(2);
                kk = kk+1;
            end
            jj2 = jj2+1;
        end
        ii2 = ii2+1;
    end
end

FCemp = squeeze(mean(FCe,1));

%%%%%%%%% Set General Model Parameters

dtt    = 1e-3;   % Sampling rate of simulated neuronal activity (seconds)
dt     = 0.1;
taon   = 100;
taog   = 10;
gamma  = 0.641;
sigma  = 0.01;
JN     = 0.15;
I0     = 0.382;
Jexte  = 1.;
Jexti  = 0.7;
w      = 1.4;

% Define optimal parameters 
we     = 2.1000;
wgaine = 0;
wgaini = 0;

J            = Balance_J9(we,C);
Tmaxneuronal = (Tmax+10)*2000;
WG           = 0:0.002:0.4;

% Model Simulations
for wge = WG
   for ii = 1:20
	[status,cmdout] = system('od -vAn -N4 -tu4 < /dev/urandom');
	rng(str2num(cmdout));
        wgaine = wge;
        wgaini = 0;
        kk = 1;
        kk3 = 1;
        for nsub = 1:NSUB
            neuro_act = zeros(Tmaxneuronal,N);
            sn = 0.001*ones(N,1);
            sg = 0.001*ones(N,1);
            nn = 1;
            for t = 0:dt:Tmaxneuronal
                xn = I0*Jexte+w*JN*sn+we*JN*C*sn-J.*sg;
                xg = I0*Jexti+JN*sn-sg;
                rn = phie9(xn);
                rg = phii9(xg);
                sn = sn+dt*(-sn/taon+(1-sn)*gamma.*rn./1000.)+sqrt(dt)*sigma*randn(N,1);
                sn(sn>1) = 1;
                sn(sn<0) = 0;
                sg = sg+dt*(-sg/taog+rg./1000.)+sqrt(dt)*sigma*randn(N,1);
                sg(sg>1) = 1;
                sg(sg<0) = 0;
                if abs(mod(t,1)) < 0.01
                    neuro_act(nn,:) = rn';
                    nn = nn+1;
                end
            end
            nn = nn-1;
            
            %%% BOLD empirical
            %   Friston BALLOON MODEL
            
            T = nn*dtt;                     % Total time in seconds          
            B = BOLD(T,neuro_act(1:nn,1)'); % B=BOLD activity, bf=Foutrier transform, f=frequency range)
            BOLD_act = zeros(length(B),N);
            BOLD_act(:,1) = B;
            
            for nnew = 2:N
                B = BOLD(T,neuro_act(1:nn,nnew));
                BOLD_act(:,nnew) = B;
            end
            
            bds = BOLD_act(2000:2000:end,:);
            
            FCs(nsub,:,:) = corrcoef(bds);

            %%%%%%%%%%%%
            
            for seed = 1:N
                ts = detrend(demean(bds(1:Tmax,seed)'));
                ts(find(ts>3*std(ts))) = 3*std(ts);
                ts(find(ts<-3*std(ts))) = -3*std(ts);
                tss(seed,:) = filtfilt(bfilt2,afilt2,ts);
            end
            
            ii2 = 1;
            for t = 1:3:190
                jj2 = 1;
                cc  = corrcoef((tss(:,t:t+30))');
                for t2 = 1:3:190
                    cc2 = corrcoef((tss(:,t2:t2+30))');
                    ca  = corrcoef(cc(Isubdiag),cc2(Isubdiag));
                    if jj2 > ii2
                        cotsamplingsim(kk) = ca(2);
                        kk = kk+1;
                    end
                    jj2 = jj2+1;
                end
                ii2 = ii2+1;
            end      
        end
        
        [hh pp FCDfitt(ii)] = kstest2(cotsampling,cotsamplingsim);
        
        FCsimul = squeeze(mean(FCs,1));
        
        cc = corrcoef(FCemp(Isubdiag),FCsimul(Isubdiag));
        fitting(ii) = cc(2);
  end
end

FCDfitt_LSD = FCDfitt;
fitting_LSD = fitting;

save('fgain_LSD.mat','WG','FCDfitt_LSD','fitting_LSD');


