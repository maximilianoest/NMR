function [jw] = Get_J_Lipari_Szabo(COR,omega,DeltaT,Steps,picoSecond)

[Np,lags] = size(COR);
Time = (0:lags-1)';

lim1 = 3500*picoSecond/(Steps*DeltaT);
lim2 = 40000*picoSecond/(Steps*DeltaT);
lim2 = 0.8;

Lipari1 = @ (b,x)   b(1)^2*exp(-x/b(3)) +b(2) ;
Lipari2 = @ (b,x)   (1-b(1)^2)*exp(-x/b(2))+b(1)^2*exp(-x/b(3)) ;

opts = optimset('Display','off');
for ii=1:Np
    F=(COR(ii,:)')/COR(ii,1);

    p1=lsqcurvefit( Lipari1, [1,200,0], Time(1:1:round(lim1*lags)) ...
        ,F(1:1:round(lim1*lags)),[0 0 20],[1 1 600],opts); %perform fit
    Fit1= p1(1)^2*exp(-Time(1:round(0.1*lags))/p1(3))+p1(2);
    
    p2=lsqcurvefit( Lipari2, [0.7,5000,p1(3)] ...
        , Time(1:4:round(lim2*lags)),F(1:4:round(lim2*lags)) ...
        ,[0 100 0.8*p1(3)],[1 2000000 1.2*p1(3)],opts); %perform fit
    Fit2=(1-p2(1)^2)*exp(-Time/p2(2))+p2(1)^2*exp(-Time/p2(3));
    
%     plot(abs(real(Fit1)))
%     hold on
%     plot(abs(real(F)),'x')
%     hold off
%     pause(0.2)
    
%     plot(real(Fit2))
%     hold on
%     plot(real(F(1:round(lim2*lags))))
%     hold off
%     pause(0.2)
    
    
    Amp=p2(1);
    Tau_glo=p2(3)*DeltaT;
    Tau_eff=p2(2)*DeltaT;
    
%          jw_fft=@(t) ((1-Amp^2)*exp(-t/Tau_eff) ...
%               +Amp^2*exp(-t/Tau_glo)).*cos(omega*t);
%          jw(ii)=integral(jw,0,inf);

    jw(ii)=COR(ii,1)*( Amp^2*Tau_glo/(1+(Tau_glo*omega)^2) ...
        +  (1-Amp^2)*Tau_eff/(1+(Tau_eff*omega)^2)   );

end

jw=jw';

end

