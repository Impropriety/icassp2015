function X=swSTFT_fb(x,fs,Nwin,Nfft,h,hop,flag_keep_transients,flag_exact_length)
% function X=swSTFT_fb(x,fs,Nwin,Nfft,h,hop)
% 
% Returns STFT with phase unwrapped in subbands and the filter transient at
% the beginning removed.
% 

    if ~exist('hop')
        hop = 1;
    end

    if ~exist('flag_keep_transients') || isempty('flag_keep_transients')
        flag_keep_transients = 0;
    end
    
    if ~exist('flag_exact_length')
        flag_exact_length = 0;
    end

%     Xstft = swSTFT(x,fs,Nwin,Nfft,1,h);
    Xstft = swSTFT(x,fs,Nwin,Nfft,hop,h);
    
    % for each subband, unwrap the phase:
    X = zeros(size(Xstft));
    Nframes = size(Xstft,2);
    for ii=1:size(Xstft,1)
        fsb = (ii-1)/Nfft;
        X(ii,:) = abs(Xstft(ii,:)).*exp(1i.*(unwrap(angle(Xstft(ii,:)))-cumsum(2*pi*fsb*hop.*ones(1,Nframes))));
    end
%     X=X(:,1:hop:end);
    
    if ~flag_keep_transients
        X = X(:,(round(Nwin/hop)+1):(end-round(Nwin/hop))); %remove the filter transients at beginning and end
    elseif flag_exact_length
        X = X(:,1:length(x)/hop); %remove the filter transients at beginning and end
    end

end
