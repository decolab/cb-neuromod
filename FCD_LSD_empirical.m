%function FCD_LSD_empiric

%% FUNCTION TO CALCULATE FCD DISTRIBUTION FOR EACH SUBJECT
%
%   Whole-brain multimodal neuroimaging model using serotonin receptor maps explain non-linear functional effects of LSD
%   Deco,G., Cruzat,J., Cabral, J., Knudsen,G.M., Carhart-Harris,R.L., Whybrow,P.C., 
%       Logothetis,N.K. & Kringelbach,M.L. (2018) Current Biology
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Subjects=15;
Conditions=[2 5]; % 2=LSD rest, 5=PLACEBO rest
TR=2;

%load fMRI data
load LSDnew.mat tc_aal
[N, Tmax]=size(tc_aal{1,1}); % N = number of areas; Tmax = total time

% FILTER SETTINGS             
fnq=1/(2*TR);                 % Nyquist frequency
flp = .04;                    % lowpass frequency of filter
fhi = 0.07;                   % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter

N_windows=length(1:3:190);
Isubdiag = find(tril(ones(N),-1)); % Indices of triangular lower part of matrix
cotsampling=zeros(length(Conditions),Subjects,N_windows*(N_windows-1)/2);

% Loop over condtions and subjects

for task=1:length(Conditions)

    for s=1:Subjects
        kk=1;
        disp(['Subject ' num2str(s)])
        signal = tc_aal{s,Conditions(task)};
        signal_filt=zeros(size(signal));
        
        %Get the BOLD phase using the Hilbert transform
        for seed=1:N
            ts=demean(detrend(signal(seed,:)));
            ts(ts>3*std(ts))=3*std(ts); % Remove strong artefacts
            ts(ts<-3*std(ts))=-3*std(ts); % Remove strong artefacts
            signal_filt(seed,:) =filtfilt(bfilt,afilt,ts); % Band pass filter
        end
        
        % For each pair of sliding windows calculate the FC at t and t2 and
        % compute the correlation between the two.
        ii2=1;
        for t=1:3:190
            jj2=1;
            cc=corrcoef((signal_filt(:,t:t+30))');
            for t2=1:3:190
                cc2=corrcoef((signal_filt(:,t2:t2+30))');
                ca=corr2(cc(Isubdiag),cc2(Isubdiag));
                if jj2>ii2
                    cotsampling(task,s,kk)=ca;
                    kk=kk+1;
                end
                jj2=jj2+1;
            end
            ii2=ii2+1;
        end
    end
end

[h_lsd, x]=hist(cotsampling(1,:),-.1:.025:1);
h_pla=hist(cotsampling(2,:),-.1:.025:1);

save FCD_values_Empirical

figure

subplot(1,2,1)
bar(x,[h_lsd' h_pla'])
ylabel('Count')
xlabel('FCD values')
legend('LSD','Placebo')
title('Empirical data')
box off


