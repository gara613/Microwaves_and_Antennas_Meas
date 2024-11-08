%% function plotPattern(E,theta,phi,varargin)
% Radiation pattern plot
% Inputs:
%   * E:                (Array), Supported sizes: N_theta x N_Phi x N_Pats, N_theta x N_Phi, N_theta x N_Pats, N_theta x 1, N_phi x 1
%   * theta:            (Array), 1D or 2D
%   * phi:              (Array), 1D or 2D
%   * varargin:     
%       - norma:        (Logic), true if normalized plot is desired
%       - scale:        (String), indicate the field axis scaling, 'dB', 'lin'
%       - plStyle:      (String), representation of the input field, 'asis' -> E, 'fld' -> abs(E), 'pwr' -> abs(E).^2, 'dir' -> dir*abs(E)/max(abs(E))
%       - coords:       (String), format of the axes, 1D fields: 'rect','polar', 2D fields: 'rect','proj2D','pol3D'
%                           'rect' for 3D patterns creates a projection on a rectangular (theta, phi) grid
%                           'proj2D' for 3D patterns creates a projection on a Direction cosine coordinate system (u=sin(th)cos(ph),v=sin(th)sin(ph)).
%       - limsScale:    (Vector), [minScale,maxScale]
%       - thetaCut:     (scalar), -1 for sweep, other value for an specific angle cut (3D patterns)
%       - phiCut:       (scalar), -1 for sweep, other value for an specific angle cut (3D patterns)
%       - multPatts:	(logic), true to indicate that multiple patterns are provided, false by default.
%       - aggrPatts:    (logic), true to indicate that patterns should be added before plotting, false by default. 
% Outputs:
%   * No variable neither file output is generated, just the pattern plot
%
% ex: 
% N = 8; 
% R = 30; 
% ampt_az = myChebyCoeffs(N,R);
% ampt_el = myChebyCoeffs(N,R);
% 
% pha_az1 = exp(-1i*2*pi*0.6*(1:N)'*sin(0)*cos(0));
% pha_el1 = exp(-1i*2*pi*0.6*(1:N)'*sin(0)*sin(0));
% 
% pha_az2 = exp(-1i*2*pi*0.6*(1:N)'*sin(pi/6)*cos(0));
% pha_el2 = exp(-1i*2*pi*0.6*(1:N)'*sin(pi/6)*sin(0));
% 
% pha_az3 = exp(-1i*2*pi*0.6*(1:N)'*sin(pi/6)*cos(pi/6));
% pha_el3 = exp(-1i*2*pi*0.6*(1:N)'*sin(pi/6)*cos(pi/6));
% 
% B1 = (ampt_az.*pha_az1)*(ampt_el.*pha_el1).';
% B2 = (ampt_az.*pha_az2)*(ampt_el.*pha_el2).';
% B3 = (ampt_az.*pha_az3)*(ampt_el.*pha_el3).';
% 
% AF(:,:,1) = arrayFactCalc(B1, [0.6,0.6], 'RA','N_Samples',[181,361]);
% AF(:,:,2) = arrayFactCalc(B2, [0.6,0.6], 'RA','N_Samples',[181,361]);
% AF(:,:,3) = arrayFactCalc(B3, [0.6,0.6], 'RA','N_Samples',[181,361]);
% 
% theta = (0:180)*pi/180;
% phi = (0:360)*pi/180;
% 
% plotPattern(AF,theta,phi,'norma',false,'scale','dB','plStyle','fld',...
%     'coords','pol3D','limsScale',[-15,25],'cut2Plot',[-1,-1],'multPatts',true,'aggrPatts',false)
% 
% plotPattern(AF,theta,phi,'norma',false,'scale','dB','plStyle','fld',...
%     'coords','polar','limsScale',[-15,25],'cut2Plot',[-1,90],'multPatts',true,'aggrPatts',false)
% 
% plotPattern(AF,theta,phi,'norma',false,'scale','dB','plStyle','pwr',...
%     'coords','proj2D','limsScale',[-15,25],'cut2Plot',[-1,-1],'multPatts',true,'aggrPatts',true)
% 
% For further examples, see  'test_plotPattern.m'
%
% Germán Augusto Ramírez
% EPFL MAG, May 2021


function plotPattern(E,theta,phi,varargin)
    [norma,scale,plStyle,coords,limsScale,thetaCut,phiCut,multPatts,aggrPatts,drawAxes] = ...
        readVarIns(varargin);
       
%% Process angular sampling
    if length(size(theta))==2 & length(size(phi))==2 & any(size(theta)==1 | size(phi)==1)
        theta = theta(:);
        phi = phi(:);        
        [th,ph] = ndgrid(theta,phi);
    elseif length(size(theta))==2 & length(size(phi))==2 & size(theta) == size(phi)
        [th,ph] = deal(theta,phi);
        theta = unique(theta);
        phi = unique(phi);
    end
    
%% Read and validate field sampling
    % Single 1D pattern, theta sweep
	if length(size(E))==2 & any(size(E)==1) & length(E) == length(theta)
        E = E(:);
        Ntheta = length(E);
        Nphi = 1; 
        Npatts = 1;
        axLabel = '\theta';
    % Single 1D pattern, phi sweep
    elseif length(size(E))==2 & any(size(E)==1) & length(E) == length(phi)
        E = E(:);
        theta = phi; 
        phi = 1; 
        Ntheta = length(E);
        Nphi = 1; 
        Npatts = 1;
        axLabel = '\phi';
	elseif length(size(E))==2 
        % Single 3D pattern
        if size(E,1) == size(th,1) & size(E,2) == size(ph,2) & ~multPatts
            % size(E,1) == length(theta) & size(E,2) == length(phi) & ~multPatts
            [Ntheta,Nphi] = size(E);
            Npatts = 1;
            axLabel = '\theta';
        % Single 3D pattern, swapped angular convention
        elseif size(E,1) == length(phi) & size(E,2) == length(theta) & ~multPatts
            E = E.';
            [Ntheta,Nphi] = size(E);
            Npatts = 1;
            axLabel = '\phi';
        % Multiple 1D patterns
        elseif size(E,1) == length(theta) | size(E,2) == length(theta) & multPatts
            axLabel = '\theta';
            if size(E,1) == length(theta) 
                E = E.';
            end
            [Npatts,Ntheta] = size(E);
            Nphi = 1;
        end 
	% Multiple 3D patterns
	elseif length(size(E))== 3
        [Ntheta,Nphi,Npatts] = size(E);
        if phiCut ~= -1
            axLabel = '\theta';
        elseif thetaCut ~= -1
            axLabel = '\phi';
        end
    else
        error('Inconsistent fields sampling');
	end
       
%% Optional variables
%% Field representation
    switch plStyle 
        case 'asis'
            F2plot = E;
        case 'fld'
            F2plot = abs(E);
        case 'pwr'
            F2plot = abs(E).^2;
        case 'dir' 
            if Ntheta == 1 | Nphi == 1 
                error('Cannot calculate directivity from 1D pattern');
            end
            % Assuming provided E is a complex field value
            F2plot = normField(abs(E).^2,Nphi,multPatts);   % Normalize and then scale by directivity
            FDir = zeros(Npatts,1);
            for contPatts = 1:Npatts 
                FDir(contPatts) = calcDir(F2plot(:,:,contPatts),th,ph);	% calcDir output is in dB!
            end  
            F2plot = bsxfun(@times,abs(F2plot),10.^(reshape(FDir,1,1,[])/10));
    end
    
%% Normalize pattern if required
    if norma
        F2plot = normField(F2plot,Nphi,multPatts); 
    end
%% Linear or dB scale, 
    if strcmpi(scale, 'dB')
    	F2plot = 10*log10(F2plot);
    % Polar axes (require nonegative values)
        if any(strcmpi(coords,{'polar','pol3D'}))
            F2plot = normVal(F2plot,min(limsScale),max(limsScale));
            pola_Scale_Labels = linspace(min(limsScale),max(limsScale),5);
        elseif any(strcmpi(coords,{'rect','proj2D'}))
            F2plot = normVal(F2plot,min(limsScale),max(limsScale))*(max(limsScale)-min(limsScale))+min(limsScale);
        end
    % Polar axes (require nonegative values)
    elseif strcmpi(scale, 'lin') & any(strcmpi(coords,{'polar','pol3D'})) 
        F2plot = normVal(F2plot,min(limsScale),max(limsScale));
        pola_Scale_Labels = linspace(min(limsScale),max(limsScale),5);	        
    end   
       
    %% plots
%%  3D patterns
    if (thetaCut == -1) & (phiCut == -1)
	% whole 3D sweep plots
        n = 64; 
        base = [0 1 0; 1 1 0; 1 0 0]; 
        my_map = interp1(linspace(2*n,0,size(base,1)), base, fliplr(0:2*n), 'pchip');
        
        % 'rect' for 3D patterns creates a rectangular projection on a (theta, phi) grid
        if strcmpi(coords,'rect') 
            if aggrPatts
                F2plot = sum(F2plot,3);
                figure,
                imagesc(F2plot);
                colormap(my_map);	
                colorbar; 
                xlabel('\phi (deg)'); ylabel('\theta (deg)');
                title(['Field ','(',plStyle,')']);
                set(gca,'XTick',1:round(Nphi/10):Nphi); set(gca,'XTickLabel',phi(1:round(Nphi/10):end)*180/pi); 
                set(gca,'YTick',1:round(Ntheta/10):Ntheta); set(gca,'YTickLabel',theta(1:round(Ntheta/10):end)*180/pi);
            else
                for cont = 1:Npatts
                    figure,
                    imagesc(F2plot(:,:,cont));
                    colormap(my_map);	
                    colorbar; 
                    xlabel('\phi (deg)'); ylabel('\theta (deg)');
                    title(['Field ','(',plStyle,')']);
                    set(gca,'XTick',1:round(Nphi/10):Nphi); set(gca,'XTickLabel',phi(1:round(Nphi/10):end)*180/pi); 
                    set(gca,'YTick',1:round(Ntheta/10):Ntheta); set(gca,'YTickLabel',theta(1:round(Ntheta/10):end)*180/pi);
                end
            end
            
        % 'proj2D' for 3D patterns creates a Direction cosine coordinate system (u=sin(th)cos(ph),v=sin(th)sin(ph)).
        elseif strcmpi(coords,'proj2D')	% 
            xx = sin(th).*cos(ph);
            yy = sin(th).*sin(ph);
            if aggrPatts
                F2plot = sum(F2plot,3);
                figure,
                contourf(xx, yy, F2plot, n,'edgecolor','none');
                colormap(my_map);	
                colorbar;
                shading flat; axis('equal');
                xlabel('u=sin\theta cos\phi'); ylabel('v=sin\theta sin\phi');
                title('2D Pattern Projection');
            else
                for cont = 1:Npatts
                    figure,
                    zz = F2plot(:,:,cont); %.*cos(thetaGrid) * Multiplying by cos(th) is a sometimes used as weighting for visualization
                    contourf(xx, yy, zz, n,'edgecolor','none');
                    colormap(my_map);	
                    colorbar;
                    shading flat; axis('equal')
                    xlabel('u=sin\theta cos\phi'); ylabel('v=sin\theta sin\phi');
                    title('2D Pattern Projection');
                end
            end

        % 3D polar plot, common format for antenna radiation patterns
        elseif strcmpi(coords,'pol3D')
            if aggrPatts
                F2plot = sum(F2plot,3);
                
                if length(theta) == length(th(:))
                    XX = sin(theta)*cos(phi').*F2plot;        
                    YY = sin(theta)*sin(phi').*F2plot;        
                    ZZ = cos(theta)*ones(1,Nphi).*F2plot;                         
                else
                    XX = sin(th).*cos(ph).*F2plot;        
                    YY = sin(th).*sin(ph).*F2plot;        
                    ZZ = cos(th).*F2plot;     
                end
                CC = abs(sqrt(XX.^2+YY.^2+ZZ.^2));

                figure,
                if drawAxes
                    plot3([0,1.1*max(CC(:))],[0,0],[0,0], 'r', [0,0],[0,1.1*max(CC(:))],[0,0], 'g', [0,0],[0,0],[0,1.1*max(CC(:))], 'b', 'linewidth',3)
                    text(1.1*max(CC(:)),0,0,' x' ,'FontSize',18);
                    text(0,1.1*max(CC(:)),0,' y' ,'FontSize',18);
                    text(0,0,1.1*max(CC(:)),' z' ,'FontSize',18);
                end
                view(3)
                hold on, 
                surf(XX,YY,ZZ,CC, 'EdgeAlpha',0.0);
                colormap(my_map);	
                alpha(0.99); 
                title('Normalized Power Pattern (dB)');    
                colorbar('Ticks', linspace(0,1,5), 'TickLabels', {num2str(pola_Scale_Labels')});
                axis equal; axis auto; box on;
                axis off; 
                camlight; lightangle(0,45), lighting gouraud;
                camproj ('perspective') ;
            else
                for cont = 1:Npatts
                    
                    if length(theta) == length(th(:))
                        XX = sin(theta)*cos(phi').*F2plot(:,:,cont);        
                        YY = sin(theta)*sin(phi').*F2plot(:,:,cont);        
                        ZZ = cos(theta)*ones(1,Nphi).*F2plot(:,:,cont);                         
                    else
                        XX = sin(th).*cos(ph).*F2plot(:,:,cont);        
                        YY = sin(th).*sin(ph).*F2plot(:,:,cont);        
                        ZZ = cos(th).*F2plot(:,:,cont);     
                    end
                    CC = abs(sqrt(XX.^2+YY.^2+ZZ.^2));

                    figure,    
                    if drawAxes
                        plot3([0,1.1*max(CC(:))],[0,0],[0,0], 'r', [0,0],[0,1.1*max(CC(:))],[0,0], 'g', [0,0],[0,0],[0,1.1*max(CC(:))], 'b', 'linewidth',3)
                        text(1.1*max(CC(:)),0,0,' x' ,'FontSize',18);
                        text(0,1.1*max(CC(:)),0,' y' ,'FontSize',18);
                        text(0,0,1.1*max(CC(:)),' z' ,'FontSize',18);
                    end
                    view(3)
                    hold on,                    
                    surf(XX,YY,ZZ,CC, 'EdgeAlpha',0.0);
                    colormap(my_map);	
                    alpha(1);%0.99); 
                    title('Normalized Power Pattern (dB)');                        
                    colorbar('Ticks', linspace(0,1,5), 'TickLabels', {num2str(pola_Scale_Labels')});
                    axis equal; axis auto; box on;
                    axis off; 
                    camlight; lightangle(0,45), lighting gouraud;
                    camproj ('perspective') ;
                end
            end
        end        

%% 2D Patterns, specified as cuts of a complete pattern, or inputs specifying only azimuth or elevation variations
    % phi cut 
    elseif (thetaCut == -1) & any(intersect(size(F2plot), size(phi)))
        if isequal(setdiff(unique(sign(theta)),0), unique(sign(diff(theta))))        
            [~,idx] = min(abs(bsxfun(@minus,phi*180/pi,phiCut.'))); %min(abs(phi*180/pi-phiCut));
            F2plot = [F2plot(end:-1:2,mod((Nphi-1)/2+idx,361),:);F2plot(:,idx,:)];
            %F2plot = [F2plot(:,idx,:); F2plot(end:-1:1,mod((Nphi-1)/2+idx,361),:)];
            
            if strcmpi(coords,'polar')              
                theta = [-theta(end:-1:2);theta];
                %theta = [theta;-theta(end:-1:1)];
            elseif strcmpi(coords,'rect')
                theta = [-theta(end:-1:2);theta];
                %theta = [theta;-theta(end:-1:1)];
            end            
        end
	% theta cut
    elseif (phiCut == -1) & any(intersect(size(F2plot), size(theta)))
        if isequal(setdiff(unique(sign(theta)),0), unique(sign(diff(theta))))
            [~,idx] = min(abs(bsxfun(@minus,theta*180/pi,thetaCut.'))); %min(abs(theta*180/pi-thetaCut));
            F2plot = F2plot(idx,:,:);
            theta = phi;
        end
    elseif ~(length(size(F2plot))==2)
        error('Invalid pattern cut or sweep specified');
    end
 
%% Input is assumed as a Set (can be only one) of azimuth or elevation cuts
	if (thetaCut ~= -1) | (phiCut ~= -1)
        figure,    
        if strcmpi(coords,'rect')
            plot(theta*180/pi,squeeze(F2plot),'linewidth',2);
            grid on; axis([min(theta*180/pi),max(theta*180/pi),min(limsScale),max(limsScale)]);
            xlabel([axLabel, ' (deg)']); ylabel(['(' scale ')']);
        elseif strcmpi(coords,'polar')
            polarplot(theta, squeeze(F2plot),'linewidth',2);
            rticks([0, 0.25, 0.5, 0.75, 1]); rticklabels(pola_Scale_Labels);
            pax = gca; pax.ThetaDir = 'clockwise'; pax.ThetaZeroLocation = 'top';
        end            
        if norma 
            title(['Normalized radiation pattern (',plStyle,')']);         
        else
            title(['Radiation pattern (',plStyle,')']);         
        end
	end    

end

% Returns the normalized field for any of the field formats allowed
function out = normField(F,Nph,multPat)
	if length(size(F))==2 & Nph == 1 & ~multPat
    	out = F/ max(F);
	elseif length(size(F))==2 & Nph == 1 & multPat
        out = bsxfun(@rdivide, F, max(F,[],2));
	elseif length(size(F))==2 & ~multPat
    	out = F/ max(max(F));
	elseif length(size(F))==3
    	out = bsxfun(@times,F,1/max(max(F)));
	end       
end

% Read optional inputs
function [norma,scale,plStyle,coords,limsScale,thetaCut,phiCut,multPatts,aggrPatts,drawAxes] = ...
	readVarIns(optInputsArr)

    found=zeros(length(optInputsArr),1);%binary vector telling which parameter has been processed correctly

    % Create an input vector with the options and values found
    pos = find(strcmpi(optInputsArr,'norma'));
    if pos
        found(pos:pos+1)=1;
        norma = optInputsArr{pos+1};
    else
        norma = true;
    end
        
    pos = find(strcmpi(optInputsArr,'scale'));
    if pos
        found(pos:pos+1)=1;
        scale = optInputsArr{pos+1};
    else
        scale = 'dB';
    end
           
	pos = find(strcmpi(optInputsArr,'plStyle'));
    if pos
        found(pos:pos+1)=1;
        plStyle = optInputsArr{pos+1};
    else
        plStyle = 'pwr';
    end

	pos = find(strcmpi(optInputsArr,'coords'));
    if pos
        found(pos:pos+1)=1;
        coords = optInputsArr{pos+1};
    else
        coords = 'rect';
    end

    pos = find(strcmpi(optInputsArr,'limsScale'));
    if pos
        found(pos:pos+1)=1;
        limsScale = optInputsArr{pos+1};
    else
        limsScale = [-50,0];
    end
    
	pos = find(strcmpi(optInputsArr,'multPatts'));
    if pos
        found(pos:pos+1)=1;
        multPatts = optInputsArr{pos+1};
    else
        multPatts = false;
    end

	pos = find(strcmpi(optInputsArr,'aggrPatts'));
    if pos
        found(pos:pos+1)=1;
        aggrPatts = optInputsArr{pos+1};
    else
        aggrPatts = false;
    end
    
    pos = find(strcmpi(optInputsArr,'cut2Plot'));
    if pos
        found(pos:pos+1)=1;
        cut2Plot = optInputsArr{pos+1};
        if iscell(cut2Plot)
            thetaCut = cut2Plot{1}(:);
            phiCut = cut2Plot{2}(:);
        else            
            thetaCut = cut2Plot(:,1);
            phiCut = cut2Plot(:,2);
        end
    elseif any(strcmpi(coords,{'proj2D','pol3D'}))
        thetaCut = -1;
        phiCut = -1;        
    else
        thetaCut = 0;
        phiCut = 0;
    end

	pos = find(strcmpi(optInputsArr,'drawAxes'));
    if pos
        found(pos:pos+1)=1;
        drawAxes = optInputsArr{pos+1};
    else
        drawAxes = true;
    end
    
    
	if ~all(found)
        UnkPar=find(~found);
        error(['unrecognized optional input par: ',optInputsArr{UnkPar(1)}]);
	end
end