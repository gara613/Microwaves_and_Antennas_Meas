% plotSpars(freqs, Spars, varargin)
% Plots the S parameters contained in an ND-array of Spars (like the one returned by the routine readSpars)
% 
% Inputs:
%	- freqs = array containing frequency sweep    
%	- Spars = S parameters matrix or list of matrices (each with dimensions Nports x Nports x Nfreqs)
%	- varargin options = {linear, fase, cmp, parkind, {names list}} 
%       - linear: string flag indicating if the magnitude plot is in linear or logarithmic scale
%     	- fase: string flag to specify if the phase must be plotted or not, fase='fase'=true, fase=otherwise=false
%       - cmp: string flag to indicate a comparison of S matrices (if enabled, a list of S matrices must be passed) 
%       - parkind: string to select and label one S parameter even for the compared plots e.g: 'S_{11}', 'S_{21}' 
%       - name: string to be used in the title
%       - names: list of the names for the S matrices compared.
%
% Universidad Nacional de Colombia
% Germán Augusto Ramírez Arroyave
% CMUN - 2014

function plotSpars(freqs, Spars, varargin)
%defaults
nVarargs = length(varargin); 
NumMats = 1;
cmp = [];
type = 'dB';%[];
fase = [];
parkind='S_{11}';
name=[];
namelist = [];
[m,n,r] = size(Spars); 
Nfreqs=length(freqs);

% vars={'cmp','type','fase','parkind','name','namelist'};
% vars(ismember(vars,varargin));

% logic order for variables' priority: fase->name->parkind
% something more elaborated as a logical vector can be done to determine if any optional variable is passed, independently of the others
if ~isempty(varargin) 
	if nVarargs == 1
        fase = varargin{1}; % flag indicating that phase plots are required
    elseif nVarargs == 2
        fase = varargin{1}; 
        cmp = varargin{2}; %  flag indicating that comparison of S matrices is required, hence Spars must be a list of S matrices
    elseif nVarargs == 3
        fase = varargin{1}; 
        cmp = varargin{2};
        parkind = varargin{3};
	elseif nVarargs == 4
        fase = varargin{1}; 
        cmp = varargin{2};
        parkind = varargin{3};
        name = varargin{4};
	elseif nVarargs == 5
        fase = varargin{1}; 
        cmp = varargin{2};
        parkind = varargin{3};
        name = varargin{4};
    	namelist = varargin{5:end}; % list of names for the S parameters matrices compared
	end
end

if strcmp(cmp,'cmp')
    % Different S parameters matrices (all of the same size) should be passed as elements of a cell array 
    NumMats = n; % number of S matrices
    [m,n,r] = size(Spars{1}); % Size of each S matrix
end

legendstr = cell(NumMats*(m*n),1); % Cell array for the legends in plots

    %% S matrix correspond to just one S parameter
    if m == 1
        if ~strcmp(cmp,'cmp') % only one device
            if strcmp(type,'dB') 
                garPlot(freqs,squeeze(20*log10(abs(Spars(1,1,:)))),2,['Magnitude ' name],'Frecuency (Hz)',['Magnitude ' parkind ' (dB)'],'');
            elseif strcmp(type,'linear')
                garPlot(freqs,squeeze(abs(Spars(1,1,:))),2,['Magnitude ' name],'Frecuency (Hz)',['Magnitude ' parkind ' (linear)'],'');
            end
            if strcmp(fase,'fase') % Phase plots are enabled (unwrap y/n?)
                garPlot(freqs,squeeze(180/pi*( (angle(Spars(1,1,:))))),2,['Phase ', name],'Frecuency (Hz)',['Phase ' parkind ' (degrees)'],'');
            end
        elseif strcmp(cmp,'cmp') % more than one device
            matindex = 1;
            mags=zeros(NumMats,Nfreqs);
            for cont = 1:NumMats
                mags(matindex,:) = squeeze(abs(Spars{cont}(1,1,:)));%squeeze(20*log10(abs(Spars{cont}(1,1,:))));
                if ~isempty(namelist)
                    legendstr{matindex} = [parkind ' ' namelist{cont}];
            else
                	legendstr{matindex} = [parkind ' ' num2str(matindex)];
                end
                matindex = matindex + 1;
            end
            if strcmp(type,'dB') 
                garPlot(freqs,20*log10(mags),2,['Magnitude ' name],'Frecuency (Hz)','Magnitude (dB)',legendstr);
            elseif strcmp(type,'linear') 
                garPlot(freqs,mags,2,['Magnitude ' name],'Frecuency (Hz)','Magnitude (linear)',legendstr);
            end
            if strcmp(fase,'fase') 
                matindex = 1;
                phases=zeros(NumMats,Nfreqs);
                for cont = 1:NumMats                    %unwrap
                    phases(matindex,:) = squeeze(180/pi*((angle(Spars{cont}(1,1,:))))); 
                    if ~isempty(namelist)
                        legendstr{matindex} = [parkind ' ' namelist{cont}];
                    else
                        legendstr{matindex} = [parkind ' ' num2str(matindex)];
                    end
                    matindex = matindex + 1;
                end
                garPlot(freqs,phases,2,['Phase ' name],'Frecuency (Hz)','Phase (degrees)',legendstr);    
            end
        end       
    %% S matrix of just one device    
    elseif ~strcmp(cmp,'cmp') % there is not a comparison of S matrices
        matindex = 1;
        mags=zeros(m*n,Nfreqs);
        for cont1 = 1:m
            for cont2 = 1:n
                mags(matindex,:)= squeeze(abs(Spars(cont1,cont2,:)));%squeeze(20*log10(abs(Spars(cont1,cont2,:))));
                legendstr{matindex} = ['S_{' num2str(cont1) num2str(cont2) '} '];
                matindex = matindex + 1;
            end
        end
        if strcmp(type,'dB') 
            garPlot(freqs,20*log10(mags),2,['Magnitude ' name],'Frecuency (Hz)','Magnitude (dB)',legendstr);
        elseif strcmp(type,'linear') 
            garPlot(freqs,mags,2,['Magnitude ' name],'Frecuency (Hz)','Magnitude (dB)',legendstr);
        end

        if strcmp(fase,'fase') 	% Phase plots are enabled
            phases=zeros(m*n,Nfreqs);
            matindex = 1;
            for cont1 = 1:m
                for cont2 = 1:n                         %unwrap
                    phases(matindex,:)= squeeze(180/pi*((angle(Spars(cont1,cont2,:)))));
                    legendstr{matindex} = ['S_{' num2str(cont1) num2str(cont2) '} '];
                    matindex = matindex + 1;
                end
            end
            garPlot(freqs,phases,2,['Phase ' name],'Frecuency (Hz)','Phase (degrees)',legendstr);
        end
        
    %% Comparison of S Matrices for many devices    
	elseif strcmp(cmp,'cmp')
        matindex = 1;
        for cont = 1:NumMats
            for cont1 = 1:m
                for cont2 = 1:n
                    mags(matindex,:)=squeeze(abs(Spars{cont}(cont1,cont2,:)));%squeeze(20*log10(abs(Spars{cont}(cont1,cont2,:))));
                    if ~isempty(namelist)
                        legendstr{matindex} = ['S_{' num2str(cont1) num2str(cont2) '} ' namelist{cont}];
                    else
                        legendstr{matindex} = ['S_{' num2str(cont1) num2str(cont2) '} '];
                    end
                    matindex = matindex + 1;
                end
            end
        end
        if strcmp(type,'dB') 
            garPlot(freqs,20*log10(mags),2,['Magnitude ' name],'Frecuency (Hz)','Magnitude (dB)',legendstr);
        elseif strcmp(type,'linear') 
            garPlot(freqs,mags,2,['Magnitude ' name],'Frecuency (Hz)','Magnitude (dB)',legendstr);
        end
        if strcmp(fase,'fase') 
            matindex = 1;
%            phases = zeros(NumMats,length(freqs));
            for cont = 1:NumMats
                for cont1 = 1:m
                    for cont2 = 1:n                         % unwrap
                        phases(matindex,:) = squeeze(180/pi*((angle(Spars{cont}(cont1,cont2,:))))); 
                        if ~isempty(namelist)
                            legendstr{matindex} = ['S_{' num2str(cont1) num2str(cont2) '} ' namelist{cont}];
                        else
                            legendstr{matindex} = ['S_{' num2str(cont1) num2str(cont2) '} '];
                        end
                        matindex = matindex + 1;
                    end
                end
            end
            garPlot(freqs,phases,2,['Phase ' name],'Frecuency (Hz)','Phase (degrees)',legendstr);
        end
    end
end    
function garPlot(x,y,lineW,name,xlab,ylab,legendstr)
    figure    
    set(gca,'ColorOrder',[0.1 0.1 0.8; 0.5 0.5 0; 0.1 0.8 0.1; 0.8 0.1 0.1],...
        'LineStyleOrder','-|-.|--|:|^-|*-|o-',...
        'NextPlot', 'replacechildren');
    plot(x,y,'LineWidth',lineW); 
    title(name); xlabel(xlab); ylabel(ylab); legend(legendstr); grid on;
end