%% Characteristic permittivity of a covered microstrip line with a dielectric overlay (Barbuto, Alú, Bilotti, Toscano, Vegni)
%inputs:
%   - hc: height of the cover substrate
%   - h: height of the microstrip substrate
%   - w: width of the line
%   - eps_r: relative permittivity of the microstrip substrate
%   - eps_rc: relative permittivity of the cover substrate

function epsr_eff=uStripOverlayModel(hc,h,w,eps_r,eps_rc)
    % regression coefficients provided
    k1=0.52; k2=0.241; k3=0.715; k4=0.446; k5=1.814; k6=1.798; k7=12.52; 
    k8=1.877; k9=0.904; k10=0.367; k11=1.782; k12=0.782; k13=0.214;

    wn=w/h;
	hn=hc/h;

%% effective permittivity of a microstrip line with an infinite cover substrate of any epsilon
    if wn>=1 
        s = (1+12/wn)^(-0.5);
    elseif wn<1 
        s = (1+12/wn)^(-0.5) + 0.04*(1-wn);
    end

    if wn<=2 
        expCorr = k2/(wn^k3+k4);
    elseif wn>2 
        expCorr = k5/(wn^k6+k7);
    end

    if eps_r>=eps_rc
        corr = 1;
    elseif eps_r<eps_rc 
        corr = (eps_rc/eps_r).^expCorr;
    end
    
    epsr_eff.eps_r_eff_inf = 0.5*(eps_r+eps_rc) + k1*(eps_r-eps_rc)*s*corr; %infinite cover

%% effective relative permittivity of the microstrip with a (finite) dielectric overlay
    if wn>=1 
        eps_rcov = 1+2/pi*(eps_rc-1)*atan(k8*hn^k9*wn^k10);
    elseif wn<1 
        eps_rcov = 1+2/pi*(eps_rc-1)*atan(k11*hn^k12*wn^k13); 
    end
   
    if eps_r>=eps_rcov
        corr = 1;
    elseif eps_r<eps_rcov
        corr = (eps_rcov/eps_r).^expCorr;
    end
    
    epsr_eff.eps_r_eff_cov = 0.5*(eps_r+eps_rcov) + k1*(eps_r-eps_rcov)*s*corr; %finite cover
end