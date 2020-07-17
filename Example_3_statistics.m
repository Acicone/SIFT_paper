clear ; close all ;
clc

addpath('.\Conceft') ;
addpath('.\Conceft\tool') ;
%initstate(1) ;

name_file='Ex3_stat';

%% How we generated the simulated data
% Ncases=100;
% load Ex3_data
% for i=1:Ncases
%     % add noise (Gaussian white noise)
%     sigma = 1;%sqrt( var(clean)*10.^( -snrdb /10 ) );
%     noise(:,i) = random('T',4,N,1) ;
%     noise(:,i) = sigma * noise(:,i) ;
%     var(noise(:,i))
%     snrdb = 20 * log10(std(clean)./std(noise(:,i))) ;
%     fprintf(['snrdb = ',num2str(snrdb),'\n']) ;
%     
%     % simulated observed time series
%     xm(:,i) = clean + noise(:,i) ;
% end
% 
% Smooth = 0 ;
% Hemi = 0 ;

load('Ex3_stat_data')

err_IMT1_BP=zeros(1,Ncases);
err_IMT1_SST=zeros(1,Ncases);
err_IMT1_SST_nWin=zeros(1,Ncases);
err_IMT1_ALIF=zeros(1,Ncases);
err_IMT2_BP=zeros(1,Ncases);
err_IMT2_SST=zeros(1,Ncases);
err_IMT2_SST_nWin=zeros(1,Ncases);
err_IMT2_ALIF=zeros(1,Ncases);

for i=1:Ncases
    %% TFR calculation
    
    
    % setup parameters for the SST or ConceFT
    
    % number of chosen orthonormal windows for ConceFT
    NoWindowsInConceFT = 2 ;
    % number of random linear combinations of chosen windows
    NoConceFT = 1 ;
    %the window length. Ideally, it should be chosen so that
    % roughly 7-10 oscillations (ignore the multiples) are
    % included in the window. Must be odd, See hermf.m
    WindowLength = 377; %677 ;
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
    [tfr, tfrtic, tfrsq, ConceFT, tfrsqtic] = ConceFT_STFT(xm(:,i), LowFrequencyLimit,...
        HighFrequencyLimit, FrequencyAxisResolution, 1, WindowLength, NoWindowsInConceFT, WindowBandwidth, NoConceFT, 0, 0, 0) ;

    %% We use the ground truth frequency patterns for the instantaneous frequency curves
    int1=[6.9,20];
    int2=[4,8];
    freq=tfrsqtic*Hz;
    
    M_n(:,1)=if1;
    M_n(:,2)=if2;
    
    for ii=1:length(if1)
        [val,curve_n(ii,1)]=min(abs(freq-if1(ii)));
        [val,curve_n(ii,2)]=min(abs(freq-if2(ii)));
    end
    
    %% ALIF decomposition using manualy designed curves build using the observed frequencies
    
    addpath('.\ALIF')
    
    % First component
    
    opt = Settings_ALIF('ALIF.NIMFs',1,'plots',0,'saveplots',0,'ALIF.xi',1.1,'ALIF.delta',10^-4);
    
    mask1 = opt.ALIF.xi*1./(M_n(:,1)'+1)*Hz;
    
    IMT1=ALIFv5_1(xm(:,i),opt,mask1);
        
    % second IMT
    opt = Settings_ALIF('ALIF.NIMFs',1,'plots',0,'saveplots',0,'ALIF.xi',1.4,'ALIF.delta',10^-5);
    
    mask2 = opt.ALIF.xi*1./M_n(:,1)'*Hz;
    
    IMT2=ALIFv5_1(IMT1(2,:),opt,mask2);
    
    % third IMT
    opt = Settings_ALIF('ALIF.NIMFs',1,'plots',0,'saveplots',0,'ALIF.xi',1.0,'ALIF.delta',10^-6);
    
    mask3 = opt.ALIF.xi*1./(M_n(:,2)')*Hz;
    
    IMT3=ALIFv5_1(IMT2(2,:),opt,mask3);
    
    % fourth IMT
    opt = Settings_ALIF('ALIF.NIMFs',1,'plots',0,'saveplots',0,'ALIF.xi',1.4,'ALIF.delta',10^-4);
    
    mask4 = opt.ALIF.xi*1./(M_n(:,2)')*Hz;
    
    IMT4=ALIFv5_1(IMT3(2,:),opt,mask4);
    
    close all
    IMT=[IMT1(1,:);IMT2(1,:);IMT3(1,:);IMT4];
    
    
    
    %% Comparisons with SST decomposition (we devide by the number of elements in frequency domain we added)
    
    freq=tfrsqtic*Hz; % frequencies in Hz
    diff(freq) % to check visulally for 0.1 Hz differences
    Band=0.3;
    
    [h, Dh, ~] = hermf(WindowLength, 1, WindowBandwidth) ;
    [IMT_SST] = Recon_sqSTFT_v2(tfrsq, tfrsqtic, Hz, curve_n(:,1), Band, h((WindowLength+1)/2)) ;
    
    [IMT_SST(2,:)] = Recon_sqSTFT_v2(tfrsq, tfrsqtic, Hz, curve_n(:,2), Band, h((WindowLength+1)/2)) ;
    
    %% Comparisons with band pass filter approach
    
    IMT_BP=bandpass(xm(:,i),int1,Hz);
    
    IMT_BP(:,2)=bandpass(xm(:,i),int2,Hz);
    
    %% TFR calculation new Win
    
    
    % setup parameters for the SST or ConceFT
    
    % number of chosen orthonormal windows for ConceFT
    NoWindowsInConceFT = 2 ;
    % number of random linear combinations of chosen windows
    NoConceFT = 1 ;
    %the window length. Ideally, it should be chosen so that
    % roughly 7-10 oscillations (ignore the multiples) are
    % included in the window. Must be odd, See hermf.m
    WindowLength = 677 ; %377
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
    [tfr, tfrtic, tfrsq, ConceFT, tfrsqtic] = ConceFT_STFT(xm(:,i), LowFrequencyLimit,...
        HighFrequencyLimit, FrequencyAxisResolution, 1, WindowLength, NoWindowsInConceFT, WindowBandwidth, NoConceFT, 0, 0, 0) ;
    
    %% Comparisons with SST decomposition (we devide by the number of elements in frequency domain we added)
    
    freq=tfrsqtic*Hz; % frequencies in Hz
    diff(freq) % to check visulally for 0.1 Hz differences
    Band=0.3;
    
    [h, Dh, ~] = hermf(WindowLength, 1, WindowBandwidth) ;
    [IMT_SST_nWin] = Recon_sqSTFT_v2(tfrsq, tfrsqtic, Hz, curve_n(:,1), Band, h((WindowLength+1)/2)) ;
    
    
    [IMT_SST_nWin(2,:)] = Recon_sqSTFT_v2(tfrsq, tfrsqtic, Hz, curve_n(:,2), Band, h((WindowLength+1)/2)) ;
    
    
    %% differences IMTs
    
    err_IMT1_BP(i)=norm(IMT_BP(:,1)'-[s1'],2)/norm(s1,2);
    err_IMT1_SST(i)=norm(real(IMT_SST(1,:))-[s1'],2)/norm(s1,2);
    err_IMT1_SST_nWin(i)=norm(real(IMT_SST_nWin(1,:))-[s1'],2)/norm(s1,2);
    err_IMT1_ALIF(i)=norm(IMT(2,:)-[s1'],2)/norm(s1,2);
        
    err_IMT2_BP(i)=norm(IMT_BP(:,2)'-[s2'],2)/norm(s2,2);
    err_IMT2_SST(i)=norm(real(IMT_SST(2,:))-[s2'],2)/norm(s2,2);
    err_IMT2_SST_nWin(i)=norm(real(IMT_SST_nWin(2,:))-[s2'],2)/norm(s2,2);
    err_IMT2_ALIF(i)=norm(IMT(4,:)-[s2'],2)/norm(s2,2);
    
end


%% Figure
figure
subplot(4,2,1)
plot(err_IMT1_BP)
subplot(4,2,3)
plot(err_IMT1_SST)
subplot(4,2,5)
plot(err_IMT1_SST_nWin)
subplot(4,2,7)
plot(err_IMT1_ALIF)
subplot(4,2,2)
plot(err_IMT2_BP)
subplot(4,2,4)
plot(err_IMT2_SST)
subplot(4,2,6)
plot(err_IMT2_SST_nWin)
subplot(4,2,8)
plot(err_IMT2_ALIF)

%% stats
disp('BP 1')
mean(err_IMT1_BP)
std(err_IMT1_BP)

disp('SST 1')
mean(err_IMT1_SST)
std(err_IMT1_SST)

mean(err_IMT1_SST_nWin)
std(err_IMT1_SST_nWin)

disp('ALIF 1')
mean(err_IMT1_ALIF)
std(err_IMT1_ALIF)

disp('BP 2')
mean(err_IMT2_BP)
std(err_IMT2_BP)

disp('SST 2')
mean(err_IMT2_SST)
std(err_IMT2_SST)

mean(err_IMT2_SST_nWin)
std(err_IMT2_SST_nWin)

disp('ALIF 2')
mean(err_IMT2_ALIF)
std(err_IMT2_ALIF)

