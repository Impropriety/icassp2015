clear variables;

% load 0dB SNR mix of speech and car noise:
[x,fs]=wavread('CAR-WINDOWNB-1_sA_l060_n+00_i74429_x9b53d_mix.wav');

% set options for VAD algorithm:
opts.fs=fs;
opts.flag_verbose=1;

% range of thresholds to try:
threshold=linspace(0.4,0.8,11);

% run the VAD
[stat,labels]=vad_circ(x,opts,threshold);

% load ground truth labels and compare
load('CAR-WINDOWNB-1_sA_l060_n+00_i74429_x9b53d_labels.mat','ll');
llrep=repmat(ll,[length(threshold),1]);
% miss rate:
mr =sum( double(llrep==1 & labels==0), 2)./sum( double(ll==1) );
% false alarm rate:
far=sum( double(llrep==0 & labels==1), 2)./sum( double(ll==0) );

% plot the miss and false alarm rates
figure;
plot(threshold,mr);
hold on;
plot(threshold,far,'r');
legend('MR','FAR');
xlabel('Threshold');
