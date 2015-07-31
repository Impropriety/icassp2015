function [stat,labels]=vad_circ(x,opts,threshold)

    %sampling frequency
    if isfield(opts,'fs')
        fs=opts.fs;
    else
        fs=8e3; 
    end
    
    %STFT hop, determines downsampling of subbands
    if isfield(opts,'hop')
        hop=opts.hop;
    else
        hop = 16;   
    end
    
    % sliding snapshot window length for estimation of circ. coef.s:
    if isfield(opts,'M')
        M=opts.M;
    else
        M=2048/hop;
    end
    
    % hop for sliding snapshot window
    if isfield(opts,'Mhop')
        Mhop=opts.Mhop;
    else
        Mhop=floor(0.01*fs)/hop;
    end
    
    % STFT window length
    if isfield(opts,'Nwin')
        Nwin=opts.Nwin;
    else
        Nwin=1024;   %used for ICASSP 2015 submission
    end
    
    % length of FFT
    Nfft=Nwin;
    % STFT window
    h = hamming(Nwin,'periodic');   %analysis window
    [h,winCheck]=ola(h,hop);    %ensure window is perfect overlap-add
    
    % median filter window length
    if isfield(opts,'Nmedfilt')
        Nmedfilt=opts.Nmedfilt;
    else
        % set to 1 second
        Nmedfilt=floor(1/(Mhop*hop/fs));
    end
    
    % verbose flag
    if isfield(opts,'flag_verbose')
        flag_verbose=opts.flag_verbose;
    else
        flag_verbose=0;
    end
    
    Mds=1;  %don't downsample when averaging over snapshots; let hop take care of this

    % compute single-channel circularity coefficients for each
    % time-frequency:
    [kx,X] = circsb_ds(x(:,1),fs,Nfft,Nwin,h,hop,M,Mhop,Mds,0,flag_verbose,1);
    
    % compute detection statistic, the mean across frequency of the
    % magnitude of the circularity coefficient
    stat=mean(abs(kx));
    
    % apply each threshold, median-filtering the resulting labels:
    labels=zeros(length(threshold),length(stat));
    for ithresh=1:length(threshold)
        labels(ithresh,:)=stat>=threshold(ithresh);
        labels(ithresh,:)=round(medfilt1(double(labels(ithresh,:)),Nmedfilt)-0.1);
    end
