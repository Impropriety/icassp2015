function [Swr F T] = swsgram(X,fs,varargin)
% function [Swr F T] = swsgram(X,fs,(Mstft=256),(Nstft=512),(time_flag=0),(overlap=16))
%
% Wrapper for spectrogram function.
%
%

    Mstft=256;
    Nstft=512;

    if nargin>2
        if ~isempty(varargin{1})
            Mstft = varargin{1};
        end
    end
    if nargin>3
        if ~isempty(varargin{2})
            Nstft = varargin{2};
        end
    end
    if nargin>4
        time_flag=varargin{3};
    else
        time_flag=0;
    end
    if nargin>5
        overlap = varargin{4};
    else
        overlap = 16;
    end
    
    k=(0:(Nstft-1))./Nstft;
    [Swr F T] = spectrogram(X,Mstft,(overlap-1)*Mstft/overlap,Nstft,fs);
    if ~time_flag
        T = 0:(length(X)-1);
    end
    F=k(1:floor((end/2))).*fs;
%     imagesc( T ,F,20*log10(abs(Swr)));
%     climdb(60);
%     set(gca,'YDir','Normal');
end
