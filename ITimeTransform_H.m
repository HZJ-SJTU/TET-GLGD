function [Rx,t,multi,multi2] = ITimeTransform_H(Tx,fs,s,type,num,delta_or_Rep,q)
%% input:
%Tx:Time-frequency representation
%fs:Sampling frequency
%s:Window function width
%type:Type of time-frequency transform(STFT、TSST1、TSST2、TET1、TET2)
%num:Number of reconstructed modes
%delta_or_Rep,q:Reconstruction parameters for TET2 or TSST1/TSST2
%% output:
%Rx:the sum of reconstructed modes
%t:Time
    %% 参数判断
    if (nargin < 4)
       error('Input missing parameters.');
    end
    if (~strcmp(type, 'STFT'))&&(~strcmp(type, 'TSST1'))&&(~strcmp(type, 'TSST2'))&&(~strcmp(type, 'TET1'))&&(~strcmp(type, 'TET2'))
        error('type should be STFT,TSST1,TSST2,TET1,TET2');
    end
    if (strcmp(type, 'STFT'))&&(nargin > 4)
        error('Input has too many parameters.');
    end
    if (strcmp(type, 'TSST1')||strcmp(type, 'TSST2'))&&(nargin > 6)
        error('Input has too many parameters.');
    end
    if (strcmp(type, 'TET1'))&&(nargin > 5)
        error('Input has too many parameters.');
    end
    if (strcmp(type, 'TET2'))&&(nargin > 7)
        error('Input has too many parameters.');
    end
    %% 预处理
    [L,N] = size(Tx);
    t = (0:N-1)/fs; 
    %% STFT
     if (strcmp(type, 'STFT'))
        f = (0:N/2)*fs/N;
        gt = @(t) s^(-1/2)*pi^(-1/4).*exp(-t.^2/s^2/2);
        gh = conj(gt(0));
        for i = 1:L
            for j = 1:N
                Tx(i,j) = Tx(i,j)*exp(1i*2*pi*f(i)*t(j));
            end
        end
        b = sum(Tx * (f(2)-f(1)),1);
        Rx = real(b)/gh/fs*2;
        multi = zeros(1,N);
        multi2 = zeros(L,N);
     end
     %% TSST
     if (strcmp(type, 'TSST1')||strcmp(type, 'TSST2'))
        delta = delta_or_Rep;
        %脊线位置提取
        [Cs, ~] = brevridge_mult(abs(Tx'), 0:N-1, num, 0.009, 10);
        %关键参数计算
        gf = @(w) sqrt(2 * s)*pi^(1/4)*exp(-(s * w).^2/2);
        gh = conj(gf(0));
        %各个分量重构
        multi = zeros(N,num);
        multi2 = zeros(L,N);
        for k = 1:num
            %脊线提取
            Ex = zeros(L,N);
            for i =1:L
                minvalue = max(Cs(k,i)-round(delta/2),1);
                maxvalue = min(Cs(k,i)+round(delta/2),N);
                multi2(i,minvalue:maxvalue) = multi2(i,minvalue:maxvalue) + abs(Tx(i,minvalue:maxvalue));
                Ex(i,minvalue:maxvalue) = Tx(i,minvalue:maxvalue)/gh;
            end
            temp = sum(Ex,2);
            %频谱扩展
            if mod(N,2)==0        
                for i = 0:N/2
                    multi(i+1,k) = temp(i+1);
                end
                for i = N/2+1:N-1
                    multi(i+1,k) = conj(temp(N-i+1));
                end
            else
                for i = 0:(N-1)/2
                    multi(i+1,k) = temp(i+1);
                end
                for i = (N+1)/2:N-1
                    multi(i+1,k) = conj(temp(N-i+1));
                end
            end
            %重构计算
            multi(:,k) = ifft(multi(:,k));
        end  
        multi = real(multi);
        %重构之和
        Rx = sum(multi,2);
     end
     %% TET1
     if (strcmp(type, 'TET1'))
        [Cs, ~] = brevridge_mult(abs(Tx'), 0:N-1, num, 0.009, 10);  
        gt = @(t) s^(-1/2)*pi^(-1/4).*exp(-t.^2/s^2/2);
        gh = conj(gt(0));
        multi = zeros(N,num);
        multi2 = zeros(L,N);
        for k = 1:num
            if mod(N,2)==0        
                for i = 0:N/2
                    multi(i+1,k) = Tx(i+1,Cs(k,i+1));
                    multi2(i+1,Cs(k,i+1)) = multi2(i+1,Cs(k,i+1)) + abs(multi(i+1,k));
                end
                for i = N/2+1:N-1
                    multi(i+1,k) = conj(Tx(N-i+1,Cs(k,N-i+1)));
                end
            else
                for i = 0:(N-1)/2
                    multi(i+1,k) = Tx(i+1,Cs(k,i+1));
                    multi2(i+1,Cs(k,i+1)) = multi2(i+1,Cs(k,i+1)) + abs(multi(i+1,k));
                end
                for i = (N+1)/2:N-1
                    multi(i+1,k) = conj(Tx(N-i+1,Cs(k,N-i+1)));
                end
            end
            multi(:,k) = ifft(multi(:,k))/gh;
        end    
        Rx = sum(multi,2);
     end
     if (strcmp(type, 'TET2'))
        multi2 = zeros(L,N);
        Rep = delta_or_Rep;
        %脊线位置提取
        [Cs, ~] = brevridge_mult(abs(Tx'), 0:N-1, num, 0.009, 10);
        %关键参数计算
        M0 = sqrt(-q+s^2);
        N0 = exp(-0.5*(-Rep./M0).^2); 
        %各个分量重构
        multi = zeros(N,num);
        for k = 1:num
            %脊线提取
            Ex = zeros(L,1);
            for i =1:L
                multi2(i,Cs(k,i)) = multi2(i,Cs(k,i)) + abs(Tx(i,Cs(k,i)));
                Ex(i) = Tx(i,Cs(k,i))*pi^(1/4)/sqrt(s)*M0(i,Cs(k,i))*N0(i,Cs(k,i));
            end
            %频谱扩展
            if mod(N,2)==0        
                for i = 0:N/2
                    multi(i+1,k) = Ex(i+1);
                end
                for i = N/2+1:N-1
                    multi(i+1,k) = conj(Ex(N-i+1));
                end
            else
                for i = 0:(N-1)/2
                    multi(i+1,k) = Ex(i+1);
                end
                for i = (N+1)/2:N-1
                    multi(i+1,k) = conj(Ex(N-i+1));
                end
            end
            %重构计算
            multi(:,k) = ifft(multi(:,k));
        end  
        multi = real(multi);%去除计算误差
        %重构之和
        Rx = sum(multi,2);
     end
end