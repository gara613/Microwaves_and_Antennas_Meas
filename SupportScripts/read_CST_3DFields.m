% VisFields
% Format: 
% x [mm]    y [mm]	z [mm]	ExRe [V/m]	EyRe [V/m]	EzRe [V/m]	ExIm [V/m]	EyIm [V/m]	EzIm [V/m]  
% close all, clearvars, clc

%filesPath = 'D:\gramirez\CST\04 AntennaCoupling\TwoMonopolesCoupling\Export\3d';
%filename = 'E_field_f=9.75_pt1_VolEx.txt';
%filename = 'E_field_f=9.75_pt1_SurEx.txt';     % This option exports the fields only over the conductor surfaces
%filename = 'H_field_f=9.75_pt1_VolEx.txt';      % CST default meshing
% filesPath = 'D:\gramirez\CST\04 AntennaCoupling\00 ReactionTh_Based\TwoDipoles_PerforatedScreen_Jonsson\Export\3d';
% filename = 'E-field_f=0.3_2_FS.txt';

function [F,sampGrid] = read_CST_3DFields(filename)    

    fileID = fopen(fullfile(filename));

    hdr_text = textscan(fileID,'%s',18);
    frewind(fileID);    % Set the file position indicator to the beginning of the file.
    F_data = textscan(fileID,[repmat('%f',[1,3]),repmat('%f',[1,6])],'CollectOutput',1, 'HeaderLines',2,'Delimiter','\t'); %10.6f
    fclose(fileID);

    sampGrid = F_data{1}(:,1:3);
    fieldData = F_data{1}(:,4:9);

    Fx = complex(fieldData(:,1),fieldData(:,4));
    Fy = complex(fieldData(:,2),fieldData(:,5));
    Fz = complex(fieldData(:,3),fieldData(:,6));
    F = [Fx, Fy, Fz];

%%
    debug = false
    if debug        
        magField = sqrt(Fx.*conj(Fx) + Fy.*conj(Fy) + Fz.*conj(Fz));

        X = sampGrid(:,1);
        Y = sampGrid(:,2);
        Z = sampGrid(:,3);
    
        nx = length(unique(X));
        ny = length(unique(Y));
        nz = length(unique(Z));
    
        scVal = magField;
        scVal(scVal==0) = nan;

        figure,
        scatter3(X,Y,Z,scVal,10*log10(scVal))
        xlabel('x (mm)');   ylabel('y (mm)');   zlabel('z (mm)');

        x_0_cut = X==0;
        
        figure,
        %scatter(Y(x_0_cut),Z(x_0_cut),magField(x_0_cut),10*log10(magField(x_0_cut))); colorbar
        scatter(Y(x_0_cut),Z(x_0_cut),100,10*log10(magField(x_0_cut))); colorbar
        xlabel('y (mm)');   ylabel('z (mm)'); 

        figure, 
        imagesc(10*log10(reshape(magField(x_0_cut),ny,nz))), colorbar
        view(-90,90);
        %unique(Y(x_0_cut)) 
        %xticks()
        %xticklabels();
        % unique(Z(x_0_cut))
        %yticks()
        %yticklabels();
        xlabel('y (mm)');   ylabel('z (mm)'); 
    end
end