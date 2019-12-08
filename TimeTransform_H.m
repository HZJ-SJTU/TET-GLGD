function [Wx,t,f,Tx,GroupDelay,q,Rep] = TimeTransform_H(x,fs,s,type,th)
%% input:
%x:Signal
%fs:Sampling frequency
%s:Window function width
%type:Type of time-frequency transform(STFT、TSST1、TSST2、TET1、TET2)
%th:Anti-noise threshold
%% output:
%Wx:Spectrum of STFT
%t:Time
%f:Frequency
%Tx:Time-frequency representation
%GroupDelay:GroupDelay matrix
%q,Rep:Reconstruction parameters
%% 参数判断
    if (nargin > 5)
       error('Input has too many parameters.');
    elseif(nargin == 4)
       th = 180;  
    elseif(nargin == 3)
       th = 180; type = 'STFT';
    elseif(nargin == 2)
       th = 180; type = 'STFT';s = 0.0015;
    elseif(nargin == 1 || nargin == 0 )
        error('Input missing parameters.');
    end
    if (~strcmp(type, 'STFT'))&&(~strcmp(type, 'TSST1'))&&(~strcmp(type, 'TSST2'))&&(~strcmp(type, 'TET1'))&&(~strcmp(type, 'TET2'))
        error('type should be STFT,TSST1,TSST2,TET1,TET2');
    end
    %% 判断是否为列向量，是则转置
    [xrow,~] = size(x);
    if (xrow~=1)
        x = x';
    end
    %% 预处理
    N = length(x);
    E = mean(abs(x));
    tao = (0:N-1)/fs;   
    t = (0:N-1)/fs; 
    if mod(N,2)==0
        f = (0:N/2)*fs/N;
    else
        f = (0:(N-1)/2)*fs/N;
    end
    L = length(f);
    dt = 1/fs;
    %% 计算
    %计算Wx
    gt = @(t) s^(-1/2)*pi^(-1/4).*exp(-t.^2/s^2/2);
    Wx = zeros(L, N);
    for ptr = 1:N
        gh = gt(tao-t(ptr));
        gh = conj(gh);
        xcpsi = fft(gh .* x);
        Wx(:, ptr) = xcpsi(1:L);
    end
    if ~strcmp(type, 'STFT')      
        %计算dWx
        gt = @(t) s^(-1/2)*pi^(-1/4).*exp(-t.^2/s^2/2).*(-1i*t);
        dWx = zeros(L, N);
        for ptr = 1:N
            gh = gt(tao-t(ptr));
            gh = conj(gh);
            xcpsi = fft(gh .* x);
            dWx(:, ptr) = xcpsi(1:L);
        end 
        if (strcmp(type, 'TET2'))||(strcmp(type, 'TSST2'))
            %计算ddWx
            gt = @(t) s^(-1/2)*pi^(-1/4).*exp(-t.^2/s^2/2).*(-t.^2);
            ddWx = zeros(L, N);
            for ptr = 1:N
                gh = gt(tao-t(ptr));
                gh = conj(gh);
                xcpsi = fft(gh .* x);
                ddWx(:, ptr) = xcpsi(1:L);
            end
            %计算wWx
            gt = @(t) s^(-1/2)*pi^(-1/4).*exp(-t.^2/s^2/2).*(1i*t/s/s);
            wWx = zeros(L, N);
            for ptr = 1:N
                gh = gt(tao-t(ptr));
                gh = conj(gh);
                xcpsi = fft(gh .* x);
                wWx(:, ptr) = xcpsi(1:L);
            end
            %计算wdWx
            gt = @(t) s^(-1/2)*pi^(-1/4).*exp(-t.^2/s^2/2).*(t.*t/s/s-1);
            wdWx = zeros(L, N);
            for ptr = 1:N
                gh = gt(tao-t(ptr));
                gh = conj(gh);
                xcpsi = fft(gh .* x);
                wdWx(:, ptr) = xcpsi(1:L);
            end
        end      
        %% 求群延迟
        if (strcmp(type, 'TET2'))||(strcmp(type, 'TSST2'))
            %二阶
            Denominator = wdWx.*Wx-dWx.*wWx;
            Numerator = ddWx.*wWx-dWx.*wdWx;
            Numerator2 = dWx.*dWx - ddWx.*Wx;
            p = Numerator./Denominator;
            q = Numerator2./Denominator;
            for ptr = 1:N
                p(:,ptr) = p(:,ptr) - 1i*t(ptr);
            end
            Rep = real(p);
            GroupDelay = -imag(p);
        else
            %一阶
            GroupDelay = imag(dWx./Wx);
            for ptr = 1:N
                GroupDelay(:,ptr) = GroupDelay(:,ptr) + t(ptr);
            end
            Rep = 0;
            q = 0;
        end
        GroupDelay( abs(Wx) < 0.0001 ) = 0;
        
        %% TFR
        Tx = zeros(L,N);
        if (strcmp(type, 'TET1'))||(strcmp(type, 'TET2'))
            %extract
            for prt=1:L
                for b=1:N
                    m = min(max(1 + round((GroupDelay(prt,b)-0)/dt),1),N);
                    if(b == m && abs(Wx(prt,b))>th*E)
                        Tx(prt, b) = Wx(prt, b);
                    end
                end
            end
        else
            %reassign
            for prt=1:L
                for b=1:N
                    if(abs(Wx(prt,b))>th*E)
                        m = min(max(1 + round((GroupDelay(prt,b)-0)/dt),1),N);
                        Tx(prt, m) = Tx(prt, m) + Wx(prt, b)*dt;
                    end
                end
            end
        end
    else
        Rep = 0;
        q = 0;
        GroupDelay = 0;
        Tx = 0;
    end
end


