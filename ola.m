function [hola winCheck]=ola(h,hop)

    Nwin=length(h);
    Nframes=4*ceil(Nwin/hop);
    
    winCheck=zeros(Nwin+(Nframes-1)*hop,1);
    idx=1;
    for ii=1:Nframes
        winCheck(idx:idx+Nwin-1) = winCheck(idx:idx+Nwin-1)+h;
        idx=idx+hop;
    end

    hola=h./median(winCheck);
    
end