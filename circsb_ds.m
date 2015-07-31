function [k,X,H,C]=circsb_ds(x,fs,Nfft,Nwin,h,hop,M,Mhop,Mds,flag_recon,flag_verbose,flag_keep_transients,flag_exact_length,flag_demean)
% function [k,X]=circsb_ds(x,fs,Nfft,Nwin,h,hop,M,Mhop,Mds,flag_recon,...
%                                                          flag_verbose,...
%                                                          flag_keep_transients,...
%                                                          flag_exact_length,...
%                                                          flag_demean)
% 
% Estimates circularity coefficients of short windows of subbands of a
% signal x sampled at rate fs. 
% Once the STFT X is computed, the subbands have their phase unwrapped,
% resulting in an efficiently-computed filterbank STFT. Then the sample
% circularity coefficient is computed using a "snapshot" length M and snapshot
% hop Mhop, all with a snapshot downsampling factor of Mds.
%
% Inputs are
%  x     - the 2xNdur matrix of noisy mixture observations
%  fs    - sampling rate
%  Nfft  - size of FFT used in STFT
%  Nwin  - duration of STFT window
%  h     - analysis window for STFT
%  hop   - hop in samples for STFT
%  M     - number of samples in snapshot used for estimating Hermitian and
%          complementary covariances
%  Mhop  - hop in samples of snapshot window
%  Mds   - downsampling factor for snapshot
%
% Set any input to [] to use default parameters.
%
% Nfft1 is number of one-sided FFT bins, given by Nfft/2 + 1
% Nframes is number of frames in STFT
% numSnapshots is number of M-length snapshots of STFT, given by
% 1+(Nframes-M)/(M/Mhop)
% 
% Outputs are
%  k     - Nfft1 x numSnapshots matrix of real nonnegative circularity
%          spectrum for each time/frequency point, values between 0 and 1.
%  X     - Nfft1 x numFrames STFT matrices of noisy mixtures
%


    if isempty(hop)
        hop=1;
    end
    
    if ~exist('flag_recon','var')
        flag_recon = 1;
    end
    if ~exist('flag_verbose','var')
        flag_verbose = 0;
    end
    if ~exist('flag_keep_transients','var')
        flag_keep_transients = 0;
    end
    if ~exist('flag_exact_length','var')
%         if flag_keep_transients
%             flag_exact_length = 1;
%         else
%             flag_exact_length = 0;
%         end
        flag_exact_length = 0;
    end
    if ~exist('flag_demean','var')
        flag_demean=1;
    end
    
%     if ~(size(x,1)==2)
%         if size(x,2)==2
%             x = x.';
%         else
%             error('Input matrix x is wrong size! One dimension should be equal to 2.');
%         end
%     end

    if flag_verbose, disp('Computing filterbank STFT of input data...'); tic; end
    X = swSTFT_fb(x(:),fs,Nwin,Nfft,h,hop,flag_keep_transients,flag_exact_length);  %STFTfb of channel 1 of noisy measurements
    if flag_verbose, toc; end
%     if flag_verbose, disp('Computing filterbank STFT of channel 2...'); tic; end
%     X(:,:,2) = swSTFT_fb(x(2,:).',fs,Nwin,Nfft,h,hop,flag_keep_transients,flag_exact_length);  %STFTfb of channel 2 of noisy measurements
%     if flag_verbose, toc; end

    Nfft1   = size(X,1);    %number of one-sided FFT bins
    Nframes = size(X,2);    %number of STFT frames
    
    if flag_recon
        [hrecon winCheck]=ola(hanning(M),Mhop);
    end
    
    remFrames = mod(Nframes-M,Mhop); %get remainder from segmenting STFT into overlapping frames of length M
    if remFrames
        % pad out the remainder with zeros
        Xpad = zeros(Nfft1,Nframes+(Mhop-remFrames));
        Xpad(:,1:Nframes) = X;
        X=Xpad;
        clear Xpad;
        Nframes = size(X,2);
    end
    
    numSnapshots = 1 + (Nframes-M)/Mhop;

    k = zeros(Nfft1,numSnapshots);      %initialize circularity spectrum matrix
    H = zeros(size(k));                 %initialize Hermitian variance matrix
    C = zeros(size(k));                 %initialize complementary variance matrix
%     Aall = zeros(Nfft1,numSnapshots,2,2);
%     kidx = zeros(1,numSnapshots);
%     kest = zeros(Nfft1,numSnapshots,2);    %initialize estimated circularity coeff matrix
%     yr = zeros(size(X));
    idx=1;
%     X=permute(X,[2,3,1]); %permute dimensions for more efficient computation below (gets rid of squeeze call)
    if flag_verbose, tic; end
    for mm=1:numSnapshots
        
        if flag_verbose
            if (mm-1) && mod(mm,floor(numSnapshots/10))==0
                toc;
                disp(['Processed ' num2str(mm) ' of ' num2str(numSnapshots) ' snapshots.']);
                tic;
            end
        end
        
        Xsnap = X( :, idx : Mds : idx+M-1); %grab a snapshot, Nfft x M/Mds
        idx = idx + Mhop;
        
        if flag_demean
            Xsnap = bsxfun(@minus,Xsnap,mean(Xsnap,2));
        end
        C(:,mm) = sum(Xsnap.^2,2);
        H(:,mm) = eps+sum(abs(Xsnap).^2,2);
        k(:,mm) = C(:,mm)./H(:,mm);
        
        
        
    end
    
%     X=permute(X,[3,1,2]); %undo permutation from earlier
    
end
