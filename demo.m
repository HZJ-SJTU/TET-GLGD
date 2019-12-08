clc; clear; close all;
%% Parameters
    N = 2048;
    fs = 200;
    f = (0:N/2)*fs/N;
    t = (0:N-1)/fs;
    num = 1;
    load MyColormaps;
    map = mymap;
%% Test Signal    
%Mode1
    A_f1 = exp(0.005*f);
    Phi_f1 = 4*f+0.07/2*f.^2-0.0005/3*f.^3;
    GD_t1 = 4+0.07*f-0.0007*f.^2;
    X1 = A_f1.*exp(-1i*2*pi*Phi_f1);
    X1(end) = -A_f1(end);
    Y1 = [X1  conj(fliplr(X1(2:end-1)))];    
    y1 = ifft(Y1);
%Mode2
    A_f2 = 1;
    Phi_f2 = 4*f+0.01*f.*f;
    GD_t2 = 4+0.02*f;
    X2 = A_f2.*exp(-1i*2*pi*Phi_f2);
    X2(end) = -A_f2(end);
    Y2 = [X2  conj(fliplr(X2(2:end-1)))];    
    y2 = ifft(Y2);
%Mode3
    A_f3 = exp(0.008*f);
    Phi_f3 = 5*f;
    GD_t3 = 5;
    X3 = A_f3.*exp(-1i*2*pi*Phi_f3);
    X3(end) = -A_f3(end);
    Y3 = [X3  conj(fliplr(X3(2:end-1)))];    
    y3 = ifft(Y3);
%% y1
    s = 0.1;
    x = y1;
    figure
    plot(t,x)
    xlabel({'Time (s)','(a1)'},'FontSize',24);set(gca,'XTick',0:1:10);
    ylabel('Amplitude','FontSize',24);set(gca,'YTick',-0.3:0.2:0.3);
    set(gca,'FontSize',24);axis([0 10 -0.3 0.3])
    set(gca,'looseInset',[0 0 0 0]);
    set(gcf,'color','white');
    %TET1
    [~,t,f,Tx,~,~,~] = TimeTransform_H(x,fs,s,'TET1',0);   
    figure()
    imagesc(t,f,abs(Tx));axis xy;axis tight;colormap(map);
    xlabel({'Time (s)','(b1)'},'FontSize',24);set(gca,'XTick',3:1:7);
    ylabel('Frequency (Hz)','FontSize',24); set(gca,'YTick',0:25:100); 
    set(gca,'FontSize',24);axis([3 7 0 100])
    set(gca,'looseInset',[0 0 0 0]);
    set(gcf,'color','white');
    
    [Rx,t,~,~] = ITimeTransform_H(Tx,fs,s,'TET1',num);
    SNRoutput = SNR(x,Rx');        
    figure();hold on;box on;
    plot(t,x,'r')
    plot(t,x-Rx','black','Linewidth',3);  
    xlabel({'Time (s)','(c1)'},'FontSize',24);set(gca,'XTick',3:1:7);
    ylabel('Amplitude','FontSize',24);set(gca,'YTick',-0.3:0.2:0.3);
    set(gca,'FontSize',24);axis([3 7 -0.3 0.3])
    legend('Original','Error','Orientation','horizon')
    set(gca,'looseInset',[0 0 0 0]);
    set(gcf,'color','white'); 
    text(3.5,-0.2,['RQF:',num2str(SNRoutput,3),'dB'],'Color','black','FontWeight','bold','FontSize',24)
    %TET2
    [~,t,f,Tx,~,q,Rep] = TimeTransform_H(x,fs,s,'TET2',0);   
    figure()
    imagesc(t,f,abs(Tx));axis xy;axis tight;colormap(map);
    xlabel({'Time (s)','(d1)'},'FontSize',24);set(gca,'XTick',3:1:7);
    ylabel('Frequency (Hz)','FontSize',24); set(gca,'YTick',0:25:100); 
    set(gca,'FontSize',24);axis([3 7 0 100])
    set(gca,'looseInset',[0 0 0 0]);
    set(gcf,'color','white');
    
    [Rx,t,~,~] = ITimeTransform_H(Tx,fs,s,'TET2',num,Rep,q);
    SNRoutput = SNR(x,Rx');        
    figure();hold on;box on;
    plot(t,x,'r')
    plot(t,x-Rx','black','Linewidth',3);  
    xlabel({'Time (s)','(e1)'},'FontSize',24);set(gca,'XTick',3:1:7);
    ylabel('Amplitude','FontSize',24);set(gca,'YTick',-0.3:0.2:0.3);
    set(gca,'FontSize',24);axis([3 7 -0.3 0.3])
    legend('Original','Error','Orientation','horizon')
    set(gca,'looseInset',[0 0 0 0]);
    set(gcf,'color','white'); 
    text(3.5,-0.2,['RQF:',num2str(SNRoutput,3),'dB'],'Color','black','FontWeight','bold','FontSize',24)
%% y2
    s = 0.1;
    x = y2;
    figure
    plot(t,x)
    xlabel({'Time (s)','(a2)'},'FontSize',24);set(gca,'XTick',0:1:10);
    ylabel('Amplitude','FontSize',24);set(gca,'YTick',-0.12:0.06:0.12);
    set(gca,'FontSize',24);axis([0 10 -0.12 0.12])
    set(gca,'looseInset',[0 0 0 0]);
    set(gcf,'color','white');
    %TET1
    [~,t,f,Tx,~,~,~] = TimeTransform_H(x,fs,s,'TET1',0);   
    figure()
    imagesc(t,f,abs(Tx));axis xy;axis tight;colormap(map);
    xlabel({'Time (s)','(b2)'},'FontSize',24);set(gca,'XTick',3:1:7);
    ylabel('Frequency (Hz)','FontSize',24); set(gca,'YTick',0:25:100); 
    set(gca,'FontSize',24);axis([3 7 0 100])
    set(gca,'looseInset',[0 0 0 0]);
    set(gcf,'color','white');
    
    [Rx,t,~,~] = ITimeTransform_H(Tx,fs,s,'TET1',num);
    SNRoutput = SNR(x,Rx');    
    figure();hold on;box on;
    plot(t,x,'r')
    plot(t,x-Rx','black','Linewidth',3);  
    xlabel({'Time (s)','(c2)'},'FontSize',24);set(gca,'XTick',3:1:7);
    ylabel('Amplitude','FontSize',24);set(gca,'YTick',-0.12:0.06:0.12);
    set(gca,'FontSize',24);axis([3 7 -0.12 0.12])
    legend('Original','Error')
    set(gca,'looseInset',[0 0 0 0]);
    set(gcf,'color','white'); 
    text(5,-0.1,['RQF:',num2str(SNRoutput,3),'dB'],'Color','black','FontWeight','bold','FontSize',24)
    %TET2
    [~,t,f,Tx,~,q,Rep] = TimeTransform_H(x,fs,s,'TET2',0); 
    figure()
    imagesc(t,f,abs(Tx));axis xy;axis tight;colormap(map);
    xlabel({'Time (s)','(d2)'},'FontSize',24);set(gca,'XTick',3:1:7);
    ylabel('Frequency (Hz)','FontSize',24); set(gca,'YTick',0:25:100); 
    set(gca,'FontSize',24);axis([3 7 0 100])
    set(gca,'looseInset',[0 0 0 0]);
    set(gcf,'color','white');
    
    [Rx,t,~,~] = ITimeTransform_H(Tx,fs,s,'TET2',num,Rep,q);
    SNRoutput = SNR(x,Rx'); 
    figure();hold on;box on;
    plot(t,x,'r')
    plot(t,x-Rx','black','Linewidth',3);  
    xlabel({'Time (s)','(e2)'},'FontSize',24);set(gca,'XTick',3:1:7);
    ylabel('Amplitude','FontSize',24);set(gca,'YTick',-0.12:0.06:0.12);
    set(gca,'FontSize',24);axis([3 7 -0.12 0.12])
    legend('Original','Error')
    set(gca,'looseInset',[0 0 0 0]);
    set(gcf,'color','white'); 
    text(5,-0.1,['RQF:',num2str(SNRoutput,3),'dB'],'Color','black','FontWeight','bold','FontSize',24)
%% y3
    s = 0.1;
    x = y3;
    figure
    plot(t,x)
    xlabel({'Time (s)','(a3)'},'FontSize',24);set(gca,'XTick',0:1:10);
    ylabel('Amplitude','FontSize',24);set(gca,'YTick',-0.4:0.4:1.6);
    set(gca,'FontSize',24);axis([0 10 -0.4 1.6])
    set(gca,'looseInset',[0 0 0 0]);
    set(gcf,'color','white');
    %TET1
    [~,t,f,Tx,~,~,~] = TimeTransform_H(x,fs,s,'TET1',0); 
    figure()
    imagesc(t,f,abs(Tx));axis xy;axis tight;colormap(map);
    xlabel({'Time (s)','(b3)'},'FontSize',24);set(gca,'XTick',4.75:0.25:5.25);
    ylabel('Frequency (Hz)','FontSize',24); set(gca,'YTick',0:25:100); 
    set(gca,'FontSize',24);axis([4.75 5.25 0 100])
    set(gca,'looseInset',[0 0 0 0]);
    set(gcf,'color','white');  
    
    [Rx,t,~,~] = ITimeTransform_H(Tx,fs,s,'TET1',num);
    SNRoutput = SNR(x,Rx');  
    figure();hold on;box on;
    plot(t,x,'r')
    plot(t,x-Rx','black','Linewidth',3);  
    xlabel({'Time (s)','(c3)'},'FontSize',24);set(gca,'XTick',475:0.25:5.25);
    ylabel('Amplitude','FontSize',24);set(gca,'YTick',-0.4:0.4:1.6);
    set(gca,'FontSize',24);axis([4.75 5.25 -0.4 1.6])
    legend('Original','Error')
    set(gca,'looseInset',[0 0 0 0]);
    set(gcf,'color','white'); 
    text(4.775,-0.3,['RQF:',num2str(SNRoutput,3),'dB'],'Color','black','FontWeight','bold','FontSize',24)
    %TET2
    [~,t,f,Tx,~,q,Rep] = TimeTransform_H(x,fs,s,'TET2',0);
    figure()
    imagesc(t,f,abs(Tx));axis xy;axis tight;colormap(map);
    xlabel({'Time (s)','(d3)'},'FontSize',24);set(gca,'XTick',4.75:0.25:5.25);
    ylabel('Frequency (Hz)','FontSize',24); set(gca,'YTick',0:25:100); 
    set(gca,'FontSize',24);axis([4.75 5.25 0 100])
    set(gca,'looseInset',[0 0 0 0]);
    set(gcf,'color','white'); 
    
    [Rx,t,~,~] = ITimeTransform_H(Tx,fs,s,'TET2',num,Rep,q);
    SNRoutput = SNR(x,Rx');     
    figure();hold on;box on;
    plot(t,x,'r')
    plot(t,x-Rx','black','Linewidth',3);  
    xlabel({'Time (s)','(e3)'},'FontSize',24);set(gca,'XTick',4.75:0.25:5.25);
    ylabel('Amplitude','FontSize',24);set(gca,'YTick',-0.4:0.4:1.6);
    set(gca,'FontSize',24);axis([4.75 5.25 -0.4 1.6])
    legend('Original','Error')
    set(gca,'looseInset',[0 0 0 0]);
    set(gcf,'color','white'); 
    text(4.775,-0.3,['RQF:',num2str(SNRoutput,3),'dB'],'Color','black','FontWeight','bold','FontSize',24)