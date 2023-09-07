%% Theoretical wave spectrum
%f = array of frequency values of equal width
%type           param
%'PM'           [fp]
%'ITTC_1'       [Hs Tp]
%'ITTC_2'       [Hs Te]
%'JONSWAP_1'    [Hs Tp gamma]
%'JONSWAP_2'    [Hs Te gamma]

function [Sf] = waveSpectrum(f, param, specType)

%% Global parameters
g = 9.81;

%% ITTC - 2 parameter Hs Tp
switch specType
    case 'Gaussian'
        Hs = param(1);
        Tp = param(2);
        fp = 1/Tp;
        sigma = 0.08;
        w = f.*(2*pi);
        wp = fp*(2*pi);
        
        Sf = (((Hs/4)^2)/(sigma*sqrt(2*pi)*wp))*exp(-((w-wp).^2)/(2*(sigma^2)*(wp^2)));
        
        Sf = Sf.*(2*pi);
        
    case 'PM'
        
        g=9.81;
        alpha = 0.0081;
        fp = param(1);
        w = f.*(2*pi);
        wp = fp.*(2*pi);
       
        Sf = (alpha.*(g.^2)./(w.^5)).*exp(-1.25.*(wp./w).^4);
               
        Sf = Sf.*(2*pi);
        
    case 'ITTC_1'
        Hs = param(1);
        Tp = param(2);
        K = (Tp/2.492)*sqrt(g/Hs);

        A = (0.0081/K^4)*g^2;
        B = (0.0081/K^4)*(4*g^2/Hs^2);

        Sf = (A./f.^5).*exp(-(B./f.^4));
    
    case 'ITTC_2'
        Hs = param(1);
        Te = param(2);
        K = (Te/2.137)*sqrt(g/Hs);

%         Tp = 7/sqrt(15);

        A = (0.0081/K^4)*g^2;
        B = (0.0081/K^4)*(4*g^2/Hs^2);

        Sf = (A./f.^5).*exp(-(B./f.^4));
        
    case 'JONSWAP_1'
        fp = 1/param(2);
        Tp = param(2);
        Hs = param(1);
        sigma_a = 0.07; 
        sigma_b = 0.09;
        alpha = 0.0081;
        gamma = param(3);
        w = f.*(2*pi);
        wp = fp*(2*pi);
        
        %as defined in the EDL wave synthesiser
        ind = find(w <= wp);
        Sf(ind,1) = ((5/(16/Tp.*(1.15+0.168.*gamma-0.925./(1.909+gamma)).*2.*pi()))...
            .*Hs^2.*gamma.^exp(-((w(ind)./(2*pi()/Tp)-1).^2)/(2.*sigma_a.^2)))...
            ./((w(ind)./(2.*pi()/Tp)).^5).*exp(-1.25./((w(ind)./(2*pi()/Tp)).^4));
        
        ind = find(w > wp);
        Sf(ind,1) = ((5/(16/Tp.*(1.15+0.168.*gamma-0.925./(1.909+gamma)).*2.*pi()))...
            .*Hs^2.*gamma.^exp(-((w(ind)./(2*pi()/Tp)-1).^2)/(2.*sigma_b.^2)))...
            ./((w(ind)./(2.*pi()/Tp)).^5).*exp(-1.25./((w(ind)./(2*pi()/Tp)).^4));
        
        Sf = Sf.*(2*pi);
        
    case 'JONSWAP_2'
        
    case 'Wide Band'
        Hs = param(1);
        fmax = param(3);
        fmin = param(2);
        
        Sf = (((0.25*Hs)^2)/(fmax-fmin))+(f.*0);
        ind = ~logical(f >= fmin & f <= fmax);
        Sf(ind) = 0;
        
    case 'B-S'
        Hs = param(1);
        Tp = param(2);
        
        Sf = 0.08.*((Hs^2.*(Tp./1.4059))./((Tp./1.4059).*f).^5).*exp(-0.318.*(1./((Tp./1.4059).*f)).^4);
        
    case 'ISSC'
       
       Hs = param(1);
       fp = param(2);
       w = f.*(2*pi);
       wp = fp*(2*pi);
        
       Sf = ((5./16).*(Hs.^2).*(wp.^4).*(w.^-5)).*exp((-5./4).*(w./wp).^-4);
       Sf = Sf.*(2*pi);
end
        
%% Plot Spectrum
% figure;
% plot(f, Sf)
% 
% setPlotSize