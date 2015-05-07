function [X F T numSampsOrig winCheck]=swSTFT(x,fs,Nwin,Nfft,Nhop,h)
% function [X F T numSampsOrig]=swSTFT(x,fs,Nwin,Nfft,Nhop,h)
%
% Performs STFT of X using overlap-add.
%
% Inputs
%  x    -    Ndur x Nch input data
%  fs   -    sampling frequency
%  Nwin -    length of time-domain window
%  Nfft -    length of FFT
%  Nhop -    window hop in samples
%  h (optional) - analysis window, should be a OLA window. If not
%  specified, a Hanning window is used.
%
% Outputs
%  X    -    struct containing the (Nfft-1)/2+1 x numFrames STFT of x
%  F    -    frequencies of the spectral bins of X in Hz
%  T    -    time indices of frames of X in seconds
%  numSampsOrig - original length of signal
%

    if ~exist('h')
        h = hamming(Nwin,'periodic');
    end

    numFrames = ceil( (size(x,1)-Nwin)./Nhop );
    numSampsOrig = size(x,1);
    % find number of samples represented by buffer matrix:
    numSamps = numFrames*Nhop+Nwin;
    % pad the end of x with zeros to equal number of samples represented
    % by buffer matrix:
    x = [x; zeros(numSamps-size(x,1),size(x,2))];

    % put an extra window length on both ends to deal with tapering
    % effects:
    x = [zeros(Nwin,1); x; zeros(Nwin,1)];
    numSamps = numSamps + 2*Nwin;               %new number of samples
    numFrames = ceil( (size(x,1)-Nwin)./Nhop ); %compute new number of frames
    
    % initialize buffer matrix
    Bx = zeros(Nwin,numFrames); 
    % break the signal into windowed frames:
    fidx=1; % frame index
    winCheck = zeros(numSamps,1);   %add up windows to check scaling
    for iframe=1:numFrames
        Bx(:,iframe) = h.*x(fidx:fidx+Nwin-1);   %window the data
%         Bx(:,iframe) = x(fidx:fidx+Nwin-1);   %window the data
        winCheck(fidx:fidx+Nwin-1) = winCheck(fidx:fidx+Nwin-1)+h;  %add up windows to check scaling
        fidx = fidx + Nhop;   %hop the frame
    end
%     Bx = bsxfun(@times,Bx,h(:)); %doesn't seem to be faster...
    
    Bx = Bx./max(winCheck); %normalize frames so that reconstruction has same amplitude
    
    % take FFT of frames:
    X = fft(Bx,Nfft);
    X = X(1:Nfft/2+1,:);  %one-sided FFT
    X = X(:,2:end);   %don't include zero-padding at beginning (iswSTFT adds this back)
    F = (0:(Nfft/2+1))./Nfft.*fs;
    T = (0:(numFrames-1)).*Nhop./fs;
    
end
