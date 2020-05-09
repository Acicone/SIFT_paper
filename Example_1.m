clear ; close all ;
clc

addpath('.\Conceft') ;
addpath('.\Conceft\tool') ;
%initstate(1) ;
name_file='Ex1';

%% How we generated the simulated data
% 	% the sampling time (100Hz sampling rate)
% 	% high sampling rate to avoid sampling issue
% Hz = 100 ;
% time = [1/Hz:1/Hz:16]' ;
% N = length(time) ;
% 
%         % Setup parameters
%         % the number of random Multitaper (rMT)
% MT = 100 ;
% DDD = 1200*1000*5 ;
% 
% initstate(1) ;
% dim = 2 ;
% 
% 
% 	%% the amplitude modulation of the simulated signal
% am1 = smooth(cumsum(randn(N,1)) ./ Hz, 200, 'loess') ;
% am1 = 2 + am1 ./ max(abs(am1)) ;
% am2 = smooth(cumsum(randn(N,1)) ./ Hz, 200, 'loess') ;
% am2 = 2 + am2 ./ max(abs(am2)) ;
% 
% % am3 = smooth(cumsum(randn(N,1)) ./ Hz, 200, 'loess') ;
% % am3 = 2 + am3 ./ max(abs(am3)) ;
% % 
% % am4 = smooth(cumsum(randn(N,1)) ./ Hz, 200, 'loess') ;
% % am4 = 2 + am4 ./ max(abs(am4)) ;
% % am1(1:500) = 0 ;
% % am2(end-600:end) = 0 ;
% 
% 
%     %% the instantaneous frequency of the simulated signal
% 
% if1 = 8 + exp(time*0.1) ;
% 
% if2 = 6 + exp(time*0.1) ;
% 
% % if3 = smooth(cumsum(randn(N,1)) ./ Hz, 200, 'loess') ;
% % if3 = 20 - 5 * if3 ./ max(abs(if3)) ;
% % if4 = smooth(cumsum(randn(N,1)) ./ Hz, 300, 'loess') ;
% % if4 = 30 + 8 * if4 ./ max(abs(if4)) ;
% 
% phi1 = cumsum(if1) / Hz ; 
% phi2 = cumsum(if2) / Hz ; 
% % phi3 = cumsum(if3) / Hz ; 
% % phi4 = cumsum(if4) / Hz ; 
% 
% 	% the simulated signal.
% s1 = am1 .* cos(2*pi*phi1) ; 
% s2 = am2 .* cos(2*pi*phi2) ; 
% 
% % s3 = am3 .* cos(2*pi*phi3) ; 
% % s4 = am4 .* cos(2*pi*phi4) ; 
% 
% clean = s1 + s2;% + s3 + s4;
% 
% % if1(1:500) = nan ;
% % if2(end-600:end) = nan ;
% % am1(1:500) = nan ;
% % am2(end-600:end) = nan ;
% 
% 
% 	%% add noise (Gaussian white noise)
% sigma = 1 ;%sqrt( var(clean)*10.^( -snrdb /10 ) );
% noise = random('T',4,N,1) ;
% noise = sigma * noise ; 
% var(noise)
% snrdb = 20 * log10(std(clean)./std(noise)) ;
% fprintf(['snrdb = ',num2str(snrdb),'\n']) ;
% 
% 	%% simulated observed time series
% xm = clean + noise ;
% 
% Smooth = 0 ;
% Hemi = 0 ;

load('Ex1_data')

%% plot Signal

hq = figure;

subplot(2, 2, 1) ;
plot(time, s1, 'color', [.7 .7 .7]) ; axis([1 15 -3 3]) ; set(gca, 'fontsize', 20) ; 
hold on; plot(time, am1, 'k', 'linewidth', 2) ;
set(gca, 'xtick', []) ; ylabel('$s_1(t)$','Interpreter','latex');

subplot(2, 2, 2) ;
plot(time, if1, 'k', 'linewidth', 2) ; set(gca, 'fontsize', 20) ; hold on ;
plot(time, if2, 'k--', 'linewidth', 2) ; %axis([0 12 0 18]) ; 
text(4, if1(400)+1, '$$\varphi''_1(t)$$', 'Interpreter','latex', 'fontsize', 24) ;
text(6, if2(600)+1, '$$\varphi''_2(t)$$', 'Interpreter','latex', 'fontsize', 24) ; set(gca, 'xtick', []) ;
axis([1 15 5 15])

subplot(2, 2, 3) ;
plot(time, s2, 'color', [.7 .7 .7]) ; axis([1 15 -3 3]) ; set(gca, 'fontsize', 20) ; 
hold on; plot(time, am2, 'k', 'linewidth', 2) ; set(gca, 'xtick', []) ; ylabel('$s_2(t)$','Interpreter','latex');

subplot(2, 2, 4) ;
plot(time, clean, 'k', 'linewidth', 2) ; set(gca, 'fontsize', 20) ;
set(gca, 'xtick', []) ; ylabel('$s_1(t)+s_2(t)$','Interpreter','latex'); axis tight ; axis([1 15 -8 10])

% subplot(3, 2, 5) ;
% plot(time, noise, 'color', [.7 .7 .7]) ; axis tight ; set(gca, 'fontsize', 20) ;
% xlabel('Time (s)','Interpreter','latex'); ylabel('$\xi(t)$','Interpreter','latex'); axis([1 15 -8 10])
% 
% subplot(3, 2, 6) ; hold off ;
% plot(time, xm, 'color', [.3 .3 .3]) ; axis tight ; set(gca, 'fontsize', 20) ;
% xlabel('Time (s)','Interpreter','latex'); ylabel('$s_1(t)+s_2(t)+\xi(t)$','Interpreter','latex'); axis([1 15 -8 10])
set(hq,'PaperPositionMode','auto');
%axis xy ; set(gca,'fontsize', 30) ; axis tight ;
set(hq,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(hq,[name_file '_Fig1.eps'],'epsc')
%% TFR calculation

% setup parameters for the SST or ConceFT

% number of chosen orthonormal windows for ConceFT
NoWindowsInConceFT = 2 ;
% number of random linear combinations of chosen windows
NoConceFT = 2 ;
%the window length. Ideally, it should be chosen so that
% roughly 7-10 oscillations (ignore the multiples) are
% included in the window. Must be odd, See hermf.m
WindowLength = 377 ;
% this is the bandwith of the chosen window. See hermf.m
% in the attached code for details.
WindowBandwidth = 10 ;
SamplingRate = Hz ;
% Setup the frequency range for the analysis
% The allowed range is 0-0.5
% This part might be tricky. This 0-0.5 limitation is
% setup under the assumption that the sampling rate is 1Hz
% After the analysis, the sampling rate should be adjusted
% so that the result is associated with the original
% sampling rate.
% In this example, the true range is [0, 0.5]*SamplingRate
HighFrequencyLimit = 0.5 ;
LowFrequencyLimit = 0 ;
% the frequency axis resolution in the final time-frequency representation
FrequencyAxisResolution = 0.0001 ;

% call the main code, which is the ConceFT based on
% synchrosqueezed short time Fourier transform (STFT)
% Output:
% tfr: STFT result
% tfrtic: frequency axis tic for the STFT
% tfrsq: synchrosqueezed STFT (it is equivalent to running ConceFT only one time)
% ConceFT: ConceFT of synchrosqueezed STFT.
% tfrsqtic: frequency axis tic for the tfrsq and ConceFT
[tfr, tfrtic, tfrsq, ConceFT, tfrsqtic] = ConceFT_STFT(clean, LowFrequencyLimit,...
    HighFrequencyLimit, FrequencyAxisResolution, 1, WindowLength, NoWindowsInConceFT, WindowBandwidth, NoConceFT, 0, 0, 0) ;

%% We use the ground truth frequency patterns for the instantaneous frequency curves

int1=[8,20];
int2=[1,11];
freq=tfrsqtic*Hz;

M_n(:,1)=if1;
M_n(:,2)=if2;

for ii=1:length(if1)
    [val,curve_n(ii,1)]=min(abs(freq-if1(ii)));
    [val,curve_n(ii,2)]=min(abs(freq-if2(ii)));
end
%%
hq = figure;
imageSQ(time, tfrsqtic*Hz, abs(tfrsq), .995) ; colormap(1-gray) ;
hold on
axis xy ; set(gca,'fontsize', 30) ;
% plot(time, if1, 'k', 'linewidth', 2) ; set(gca, 'fontsize', 30) ; hold on ;
% plot(time, if2, 'k', 'linewidth', 2) ; %axis([0 12 0 18]) ; 
%text(4, if1(400)+1, '$$\varphi''_1(t)$$', 'Interpreter','latex', 'fontsize', 30) ;
%text(6, if2(600)+1, '$$\varphi''_2(t)$$', 'Interpreter','latex', 'fontsize', 30) ; %set(gca, 'xtick', []) ;
x1=xlabel('time (s)','Interpreter','latex');
y1=ylabel('frequency (Hz)','Interpreter','latex');
plot(time, M_n(:,1), 'b', 'linewidth', 3) ;
plot(time, M_n(:,2), 'r', 'linewidth', 3) ;
axis([1 15 4 16])
lg=legend('$$\varphi''_1(t)$$','$$\varphi''_2(t)$$');
set(lg,'Interpreter','latex');
set(hq,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(hq,[name_file '_if.eps'],'epsc')

%% ALIF decomposition using manualy designed curves build using the observed frequencies

addpath('.\ALIF')

opt = Settings_ALIF('ALIF.NIMTs',1,'plots',0,'saveplots',0,'ALIF.xi',1.4,'ALIF.delta',10^-6);

mask2 = opt.ALIF.xi*1./M_n(:,1)'*Hz;

IMT2=ALIFv5_1(clean,opt,mask2);

IMT=[IMT2];


%% Comparisons with SST decomposition (we devide by the number of elements in frequency domain we added)

freq=tfrsqtic*Hz; % frequencies in Hz
diff(freq) % to check visulally for 0.1 Hz differences
Band=0.3;

[h, Dh, ~] = hermf(WindowLength, 1, WindowBandwidth) ;
[IMT_SST] = Recon_sqSTFT_v2(tfrsq, tfrsqtic, Hz, curve_n(:,1), Band, h((WindowLength+1)/2)) ;

[IMT_SST(2,:)] = Recon_sqSTFT_v2(tfrsq, tfrsqtic, Hz, curve_n(:,2), Band, h((WindowLength+1)/2)) ;

%% Comparisons with band pass filter approach

IMT_BP=bandpass(clean,int1,Hz);

IMT_BP(:,2)=bandpass(clean,int2,Hz);

%% TFR calculation after win tuning

% setup parameters for the SST or ConceFT

% number of chosen orthonormal windows for ConceFT
NoWindowsInConceFT = 1 ;
% number of random linear combinations of chosen windows
NoConceFT = 1 ;
%the window length. Ideally, it should be chosen so that
% roughly 7-10 oscillations (ignore the multiples) are
% included in the window. Must be odd, See hermf.m
WindowLength = 677 ;
% this is the bandwith of the chosen window. See hermf.m
% in the attached code for details.
WindowBandwidth = 10 ;
SamplingRate = Hz ;
% Setup the frequency range for the analysis
% The allowed range is 0-0.5
% This part might be tricky. This 0-0.5 limitation is
% setup under the assumption that the sampling rate is 1Hz
% After the analysis, the sampling rate should be adjusted
% so that the result is associated with the original
% sampling rate.
% In this example, the true range is [0, 0.5]*SamplingRate
HighFrequencyLimit = 0.5 ;
LowFrequencyLimit = 0 ;
% the frequency axis resolution in the final time-frequency representation
FrequencyAxisResolution = 0.0001 ;

% call the main code, which is the ConceFT based on
% synchrosqueezed short time Fourier transform (STFT)
% Output:
% tfr: STFT result
% tfrtic: frequency axis tic for the STFT
% tfrsq: synchrosqueezed STFT (it is equivalent to running ConceFT only one time)
% ConceFT: ConceFT of synchrosqueezed STFT.
% tfrsqtic: frequency axis tic for the tfrsq and ConceFT
[tfr, tfrtic, tfrsq, ConceFT, tfrsqtic] = ConceFT_STFT(clean, LowFrequencyLimit,...
    HighFrequencyLimit, FrequencyAxisResolution, 1, WindowLength, NoWindowsInConceFT, WindowBandwidth, NoConceFT, 0, 0, 0) ;

%% Comparisons with SST decomposition (we devide by the number of elements in frequency domain we added)

freq=tfrsqtic*Hz; % frequencies in Hz
diff(freq) % to check visulally for 0.1 Hz differences
Band=0.3;

[h, Dh, ~] = hermf(WindowLength, 1, WindowBandwidth) ;
[IMT_SST_nWin] = Recon_sqSTFT_v2(tfrsq, tfrsqtic, Hz, curve_n(:,1), Band, h((WindowLength+1)/2)) ;

[IMT_SST_nWin(2,:)] = Recon_sqSTFT_v2(tfrsq, tfrsqtic, Hz, curve_n(:,2), Band, h((WindowLength+1)/2)) ;


%% differences IMTs

hf = figure ;
plot(time,IMT_BP(:,1)'-s1'+5,'b','linewidth',2)
hold on
plot(time,real(IMT_SST(1,:))-s1','g','linewidth',2)
plot(time,real(IMT_SST_nWin(1,:))-s1'-4,'k','linewidth',2)
plot(time,IMT(1,:)-s1'-6,'r','linewidth',2)
lg=legend('BPF','SST','SST after tuning','SIFT');
set(lg,'Interpreter','latex')
axis xy ; set(gca,'fontsize', 30) ;
axis([1 15 -8 8])
set(hf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(hf,[name_file '_IMT1_diff.eps'],'epsc')

err_IMT1_BP=norm(IMT_BP(:,1)'-[s1'],2)/norm(s1,2)
err_IMT1_SST=norm(real(IMT_SST(1,:))-[s1'],2)/norm(s1,2)
err_IMT1_SST_nWin=norm(real(IMT_SST_nWin(1,:))-[s1'],2)/norm(s1,2)
err_IMT1_ALIF=norm(IMT(1,:)-[s1'],2)/norm(s1,2)

hf = figure ;
plot(time,IMT_BP(:,2)'-[s2']+6,'b','linewidth',2)
hold on
plot(time,real(IMT_SST(2,:))-[s2'],'g','linewidth',2)
plot(time,real(IMT_SST_nWin(2,:))-s2'-5,'k','linewidth',2)
plot(time,IMT(2,:)-[s2']-7,'r','linewidth',2)
lg=legend('BPF','SST','SST after tuning','SIFT');
set(lg,'Interpreter','latex')
axis xy ; set(gca,'fontsize', 30) ; 
axis([1 15 -9 9])
set(hf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(hf,[name_file '_IMT2_diff.eps'],'epsc')

err_IMT2_BP=norm(IMT_BP(:,2)'-[s2'],2)/norm(s2,2)
err_IMT2_SST=norm(real(IMT_SST(2,:))-[s2'],2)/norm(s2,2)
err_IMT2_SST_nWin=norm(real(IMT_SST_nWin(2,:))-[s2'],2)/norm(s2,2)
err_IMT2_ALIF=norm(IMT(2,:)-[s2'],2)/norm(s2,2)


%%

hq=plot_imf_v12([IMT_BP(Hz:end-Hz,1)';real(IMT_SST(1,Hz:end-Hz));...
    real(IMT_SST_nWin(1,Hz:end-Hz));IMT(1,Hz:end-Hz)],[IMT_BP(Hz:end-Hz,2)';...
    real(IMT_SST(2,Hz:end-Hz));real(IMT_SST_nWin(2,Hz:end-Hz));IMT(2,Hz:end-Hz)],...
    [ones(4,1)*s1(Hz:end-Hz)'],[ones(4,1)*s2(Hz:end-Hz)'],time(Hz:end-Hz),4,[],[],[],[],[],...
    {'BPF','SST','SST after tuning','SIFT'},{'Ground truth','Ground truth','Ground truth','Ground truth'});
saveas(hq,[name_file '_IMTs.eps'],'epsc')
