%%% Demo code for estimating the degree of noncircularity of STFT subbands,
%%% as described in
%%%
%%% [1] S. Wisdom, G. Okopal, L. Atlas, and J. Pitton, "Voice Activity
%%%     Detection Using Subband Noncircularity," in Proc. ICASSP, Brisbane,
%%%     Austrailia, 2015.
%%%
%%% If you find this code useful, please cite the paper [1].
%%%
%%% Written by S. Wisdom and G. Okopal
%%%

clear variables;

% needed for thermal colormap:
addpath('colormaps');

tstart=tic;

%% define parameters

fs=8e3; %sampling frequency
hop = 16;   %STFT hop, determines downsampling of subbands

% SUT sliding snapshot parameters:
M=2048/hop; Mhop=floor(0.01*fs)/hop; Nwin=1024;   %used for ICASSP 2015 submission
% M=512/hop; Mhop=floor(0.01*fs)/hop; Nwin=64;    %reveals speech harmonics
Mds=1;  %don't downsample when averaging over snapshots; let hop take care of this

% STFT filterbank parameters:
Nfft = Nwin; %length of FFT
h = hamming(Nwin,'periodic');   %analysis window
[h,winCheck]=ola(h,hop);    %ensure window is perfect overlap-add

% QUT-NOISE parameters
% scen='CAFE-CAFE-1';
scen='CAR-WINDOWNB-1';  %one of the 20 noise types in QUT-NOISE
% scen='STREET-CITY-1';
% scen='REVERB-POOL-1';
sesh='sA';              %session ('sA' or 'sB')
len='l060';             %duration of files ('l060' or 'l120')
snr='n+00'; SNR = 0;    %SNR in dB
fidx=15;                 %file index within condition (between 1 and 50)

name=[scen '_' sesh '_' len '_' snr];

x=wavread('CAR-WINDOWNB-1_sA_l060_n+00_i74429_x9b53d_mix.wav');
s=wavread('CAR-WINDOWNB-1_sA_l060_n+00_i74429_x9b53d_speech.wav');
v=x-s;

save_dir='';

% optional way to scale SNR:
mult=1;
if mult~=1
    x = s+mult.*v;
    disp(['Scaling SNR by ' num2str(10*log10(1/(mult^2))) ' dB']);
end

% use default mixing type (uses estimated RIR):
mixing_type='default';
% mixing_type='delays';
switch(mixing_type)
    case 'delays'
        x = s + [v(:,1),[v(2:end,1);0]];
end

% flags that control saving and loading of circularity coefficients k:
flag_load_k=1;  %if saved k's exist, load them
flag_save_k=0;  %save off the computed k's, uses about 120MB of space

%% run the algorithm
k_savedir='.';
run_name=[name '_M' num2str(M) '_Nwin' num2str(Nwin) '_hop' num2str(hop)];
if strncmp(mixing_type,'default',length('default'))
    kname=fullfile(k_savedir,['k_singleChan_' run_name '.mat']);
else
    kname=fullfile(k_savedir,['k_singleChan_' run_name '_mixing' mixing_type '.mat']);
end
if exist(kname,'file') && flag_load_k
    load(kname);
else
%     xpe=preemph(x,0.98);
    xpe=x;
    disp('Computing ground truth k of speech...');
    tic; [ks(:,:,1),S(:,:,1)] = circsb_ds(s(:,1),fs,Nfft,Nwin,h,hop,M,Mhop,Mds,0,0,1); toc;
    disp('Computing ground truth k of noise...');
    tic; [kv(:,:,1),V(:,:,1)] = circsb_ds(v(:,1),fs,Nfft,Nwin,h,hop,M,Mhop,Mds,0,0,1); toc;
    disp('Computing k of mix...');
    tic; [kx(:,:,1),X(:,:,1)] = circsb_ds(x(:,1),fs,Nfft,Nwin,h,hop,M,Mhop,Mds,0,0,1); toc;
    if flag_save_k
        save(kname,'ks','kv','kx');
    end
end

t=(0:(size(kx,2)-1))./fs.*Mhop*hop;
load('CAR-WINDOWNB-1_sA_l060_n+00_i74429_x9b53d_labels.mat','ll');

%% compute detection statistic
mx=sum(abs(kx))./size(kx,1);

mxs=sum(abs(kx))./size(kx,1);
mxs1=mxs(ll==1);
mxs0=mxs(ll==0);
edges=linspace(0,1,200);
h1=histc(mxs1,edges);
h0=histc(mxs0,edges);
figure('Position',[360,278,2*560,420]);
plot(linspace(0,1,length(h1)),h1./sum(h1),'r','LineWidth',2); hold on; plot(linspace(0,1,length(h1)),h0./sum(h0),'--','LineWidth',2);
lh=legend('Speech and noise','Noise only','FontSize',28);
ylabel('P(SCC)');
xlabel('SCC');
xlim([0.4 0.8]);
set(findall(findall(gcf,'type','axes'),'type','text'),'fontSize',28);
set(findall(gcf,'type','text'),'fontSize',28,'interpreter','latex');
set(findall(gcf,'type','axes'),'fontSize',28);
grid on;

if ~exist('S')
    S=swsgram(s(:,1),fs,1024,1024);
end
if ~exist('X')
    X=swsgram(x(:,1),fs,1024,1024);
end
if ~exist('V','var')
    V=swsgram(v(:,1),fs,1024,1024);
end
    
%% 2x3 panel figure:
axpos=0.0;
axscale=1.3;

figure('Position',[360   278   2*560   420]);
ah(3)=subplot(233);
imagesc(t,1:size(X,1),db(abs(X(:,1:8:end,1))));
axis xy; title('$\left|\mathrm{STFT}\;\,\right|^2$ in dB of mix'); climdb(60);
ch(1)=colorbar;
set(gca,'XTick',(1200:400:2400).*(Mhop*hop)./fs);

ah(2)=subplot(232);
imagesc(t,1:size(X,1),db(abs(V(:,1:8:end,1))));
axis xy; title('$\left|\mathrm{STFT}\;\,\right|^2$ in dB of isolated noise'); climdb(60);
ch(1)=colorbar;
set(gca,'XTick',(1200:400:2400).*(Mhop*hop)./fs);

ah(1)=subplot(231);
imagesc(t,1:size(X,1),db(abs(S(:,1:8:end,1))));
ylabel('Freq. bin');
axis xy; title('$\left|\mathrm{STFT}\;\,\right|^2$ in dB of isolated speech'); ylabel('Freq. bin'); climdb(60);
ch(1)=colorbar;
set(gca,'XTick',(1200:400:2400).*(Mhop*hop)./fs);

% fiddle with the axes sizes and get all the CLims to be the same:
clim_max=-Inf;
for ii=1:3
    clim=get(gca,'CLim');
    if clim(1)>clim_max
        clim_max=clim(2);
        clims(2)=clim_max;
    end
end
for ii=1:3
    axes(ah(ii));
    set(gca,'Clim',[clim_max-1,clim_max]);
    switch(ii)
        case 1
            axpos=0.05;
        case 2
            axpos=0.025;
        case 3
            axpos=0;
    end
    apos=get(gca,'Position'); set(gca,'Position',[apos(1)-axpos,apos(2),axscale*apos(3),apos(4)])     
    climdb(60);
end

ah(5)=subplot(235); imagesc(t,1:size(S,1),(abs(kv(:,1:4:end,1).^2))); axis xy; title('$\hat{k}^2(\omega,n)$ of isolated noise');
xlabel('Time (s)');
ch(2)=colorbar;
set(gca,'CLim',[0,1]);
set(gca,'XTick',(1200:400:2400).*(Mhop*hop)./fs);

ah(4)=subplot(234); imagesc(t,1:size(S,1),(abs(ks(:,1:4:end,1).^2))); axis xy; title('$\hat{k}^2(\omega,n)$ of isolated speech');
ylabel('Freq. bin'); xlabel('Time (s)');
ch(3)=colorbar;
set(gca,'CLim',[0,1]);
set(gca,'XTick',(1200:400:2400).*(Mhop*hop)./fs);

ah(6)=subplot(236); imagesc(t,1:size(S,1),(abs(kx(:,1:4:end,1).^2))); axis xy; title('$\hat{k}^2(\omega,n)$ of mix'); 
xlabel('Time (s)'); 
ch(4)=colorbar;
set(gca,'CLim',[0,1]);
linkaxes(ah);
xlim([1200 2400].*(Mhop*hop)./fs);
set(gca,'XTick',(1200:400:2400).*(Mhop*hop)./fs);

% fiddle with the axes sizes:
for ii=4:6
    axes(ah(ii));
    switch(ii)
        case 4
            axpos=0.05;
        case 5
            axpos=0.025;
        case 6
            axpos=0;
    end
    apos=get(gca,'Position'); set(gca,'Position',[apos(1)-axpos,apos(2),axscale*apos(3),apos(4)])     
end

cmap=thermal;
colormap(flipud(cmap));
set(findall(findall(gcf,'type','axes'),'type','text'),'fontSize',14);
set(findall(gcf,'type','text'),'fontSize',14,'interpreter','latex');
set(findall(gcf,'type','axes'),'fontSize',14);
