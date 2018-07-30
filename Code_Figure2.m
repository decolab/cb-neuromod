clc; clear all;

%%  CODE TO PLOT FIGURE 2
%
%   Whole-brain multimodal neuroimaging model using serotonin receptor maps explain non-linear functional effects of LSD
%   Deco,G., Cruzat,J., Cabral, J., Knudsen,G.M., Carhart-Harris,R.L., Whybrow,P.C., 
%       Logothetis,N.K. & Kringelbach,M.L. (2018) Current Biology
%
%   Illustration of the sliding-window, FCD and histogram.
%   Written by Joana Cabral joana.cabral@psyc.ox.ac.uk
%   October 2017, Barcelona-Spain
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Load BOLD signal from one subject
load LSDnew.mat tc_aal
subject=7;
scan=5; % 5 = Placebo, 2 = LSD
BOLD=tc_aal{subject,scan};
clear tc_aal
[N, Tmax]=size(BOLD);
TR=2;

% Bandpass filter settings              
fnq=1/(2*TR);                 % Nyquist frequency
flp = .02;                    % lowpass frequency of filter
fhi = 0.1;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter

BOLD_filt=zeros(size(BOLD));
        
%Get the BOLD phase using the Hilbert transform
for seed=1:N
    ts=detrend(BOLD(seed,:)-mean(BOLD(seed,:)));
    ts(ts>3*std(ts))=3*std(ts); % Remove strong artefacts
    ts(ts<-3*std(ts))=-3*std(ts); % Remove strong artefacts
    BOLD_filt(seed,:) =filtfilt(bfilt,afilt,ts); % Band pass filter
end

clear BOLD ts 

% Define parameters for Dynamic FC analysis
slide=3;
window=30;
t_windows=1:slide:Tmax-window;
N_windows=length(t_windows);
FC=zeros(N_windows,N,N);
FCD=zeros(N_windows);
Isubdiag = find(tril(ones(N),-1)); % Indices of triangular lower part of matrix
Order=[1:2:N N:-2:2];

% Compute FC for each sliding window
for nw=1:N_windows
    t=t_windows(nw);
    FC(nw,:,:)=corrcoef((BOLD_filt(:,t:t+window))');
end

% Compute the FCD matrix
for nw1=1:N_windows
    FC1=FC(nw1,:,:);
    for nw2=nw1:N_windows
        FC2=FC(nw2,:,:);
        FCD(nw1,nw2)=corr2(FC1(Isubdiag),FC2(Isubdiag));
    end
end
 
% PLOTS

figure
colormap(jet)

% Plot the filtered BOLD signal
subplot(3,5,1:3)
plot(1:TR:Tmax*TR,BOLD_filt')
box off
ylabel('BOLD signal')
xlabel('Time (s)')
set(gca,'XTick',50:50:400)
xlim([0 Tmax*TR])

% Plot some X FC matrices (equally spaced)
X=10;
for u=1:X
    nw=round(u*N_windows/X);
    subplot(3,X,X+u)
    imagesc(squeeze(FC(nw,Order,Order)))
    axis square
    set(gca,'XTick',[])
    set(gca,'YTick',[])
end
    
% Plot the FCD matrix
subplot(3,5,11)
imAlpha=triu(ones(N_windows),1);
imagesc(t_windows*TR,t_windows*TR,FCD,'AlphaData',imAlpha)
xlabel('Time X (s)')
ylabel('Time Y (s)')
title ('FCD')
axis square
colorbar

% Plot the FCD histogram
subplot(3,5,13)
histogram(FCD(imAlpha==1),50,'EdgeColor','w','Facecolor','b')
box off
xlabel('FCD values')
ylabel('Count')


