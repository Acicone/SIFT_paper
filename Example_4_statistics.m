clear ; close all ;
clc

addpath('.\Conceft') ;
addpath('.\Conceft\tool') ;
%initstate(1) ;

name_file='Ex4_stat';

%% How we generated the simulated data
% Ncases=100;
% load Ex4_data
% for i=1:Ncases
% 	% add noise (Gaussian white noise)
% sigma = 1 ;%sqrt( var(clean)*10.^( -snrdb /10 ) );
% noise(:,i) = random('T',4,N,1) ;
% noise(:,i) = sigma * noise(:,i) ;
% var(noise(:,i))
% snrdb = 20 * log10(std(clean)./std(noise(:,i))) ;
% fprintf(['snrdb = ',num2str(snrdb),'\n']) ;
%
% 	% simulated observed time series
% xm(:,i) = clean + noise(:,i) ;
% end
%
% Smooth = 0 ;
% Hemi = 0 ;

load('Ex4_stat_data')

err_IMT1_ALIF=zeros(1,Ncases);
err_IMT2_ALIF=zeros(1,Ncases);

for i=1:Ncases
  
    load([name_file '_curves'])
    
    % ALIF decomposition using manualy designed curves build using the observed frequencies
    
    addpath('.\ALIF')
    
    % First component
    
    opt = Settings_ALIF('ALIF.NIMFs',1,'plots',0,'saveplots',0,'ALIF.xi',1.4,'ALIF.delta',10^-6);
    
    mask1 = opt.ALIF.xi*1./(M_n(:,1)')*Hz;
    
    IMT1=ALIFv5_1(xm(:,i),opt,mask1);
    
    % second IMT
    opt = Settings_ALIF('ALIF.NIMFs',1,'plots',0,'saveplots',0,'ALIF.xi',1.4,'ALIF.delta',10^-6);
    
    mask2 = opt.ALIF.xi*1./M_n(:,1)'*Hz;
    
    IMT2=ALIFv5_1(IMT1(2,:),opt,mask2);
    
    % third IMT
    % opt = Settings_ALIF('ALIF.NIMFs',1,'plots',0,'saveplots',0,'ALIF.xi',1.1,'ALIF.delta',10^-5);
    %
    % mask3 = opt.ALIF.xi*1./(M_n(:,2)')*Hz;
    %
    % IMT3=ALIFv5_1(IMT2(2,:),opt,mask3);
    
    IMT3=[zeros(size(IMT2(2,:)));IMT2(2,:)];
    
    % fourth IMT
    opt = Settings_ALIF('ALIF.NIMFs',1,'plots',0,'saveplots',0,'ALIF.xi',1.4,'ALIF.delta',10^-7);
    
    mask4 = opt.ALIF.xi*1./(M_n(:,2)')*Hz;
    
    IMT4=ALIFv5_1(IMT3(2,:),opt,mask4);

    % fifth IMT
    opt = Settings_ALIF('ALIF.NIMFs',1,'plots',0,'saveplots',0,'ALIF.xi',1.1,'ALIF.delta',10^-5);
    
    mask5 = opt.ALIF.xi*1./(M_n(:,3)')*Hz;
    
    IMT5=ALIFv5_1(IMT4(2,:),opt,mask5);
    
    % sixth IMT
    opt = Settings_ALIF('ALIF.NIMFs',1,'plots',0,'saveplots',0,'ALIF.xi',1.4,'ALIF.delta',10^-6);
    
    mask6 = opt.ALIF.xi*1./(M_n(:,3)')*Hz;
    
    IMT6=ALIFv5_1(IMT5(2,:),opt,mask6);
    
    % 7th IMT
    opt = Settings_ALIF('ALIF.NIMFs',1,'plots',0,'saveplots',0,'ALIF.xi',1.1,'ALIF.delta',10^-5);
    
    mask7 = opt.ALIF.xi*1./(M_n(:,4)')*Hz;
    
    IMT7=ALIFv5_1(IMT6(2,:),opt,mask7);
    
    % 8th IMT
    opt = Settings_ALIF('ALIF.NIMFs',1,'plots',0,'saveplots',0,'ALIF.xi',1.4,'ALIF.delta',10^-4);
    
    mask8 = opt.ALIF.xi*1./(M_n(:,4)')*Hz;
    
    IMT8=ALIFv5_1(IMT7(2,:),opt,mask8);
    
    close all
    IMT=[IMT1(1,:);IMT2(1,:);IMT3(1,:);IMT4(1,:);IMT5(1,:);IMT6(1,:);IMT7(1,:);IMT8];
    
    
    %% differences IMTs
    
    err_IMT1_ALIF(i)=norm(IMT(6,Hz:end-Hz)-s1(Hz:end-Hz)',2)/norm([s1(Hz:end-Hz)'],2);
   
    err_IMT2_ALIF(i)=norm(sum(IMT([4 8],Hz:end-Hz),1)-s2(Hz:end-Hz)',2)/norm(s2(Hz:end-Hz)',2);
    
end

%% Figure
figure
subplot(1,2,1)
plot(err_IMT1_ALIF)
subplot(1,2,2)
plot(err_IMT2_ALIF)

%% stats
disp('ALIF 1')
mean(err_IMT1_ALIF)
std(err_IMT1_ALIF)

disp('ALIF 2')
mean(err_IMT2_ALIF)
std(err_IMT2_ALIF)
