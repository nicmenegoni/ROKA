%% Description of the dxf_planefit
%This code fit to a dxf file (containg the several 3D trace of the
%disconinuities) a plane represented as a circular disc. 
%This vesion of the code works with the DXF library of CC v2.9,
%but accept only polyline (no points).
%
%To read the dxf, the function f_LectDxf, developed by Steven Michel
%(2005) is used.

%% function [c_Poly,c_Poi] = f_LectDxf(nomArch)
% FUNCTION dxf = read_dxf(filename)
%
% Author:  Steven Michael (smichael@ll.mit.edu)
%
% Date:    3/10/2005
%
% Description
%
%   The following compiled MATLAB file reads in an
%   ascii DXF file and loads the information into
%   the "dxf" variable
%
% Inputs:
%
%   Filename   :    The filename of the ASCII DXF File
%
% Outputs:
%
%   dxf        :    A 3-D variable of size (NX3X3).  The first
%                   index N is the number of facets.  The second
%                   index references the 3 vertices of each
%                   facet, and the third index is the
%                   (x,y,z) location of the vertex.
%
% Note:
%
%   This may not work with all DXF files.  The DXF file format
%   is complicated and continuously evolving.  This has worked
%   with every file I happened to use it on, but that will
%   not be true for everyone.
%
%   Also, the function does not load any color or texture
%   information.
%%% Read entities information of dxf file
%author = SebaTM
%Jul 27 2009
%Based in dxf2coord 1.1 matrix of lukas wischounig, but is not dependent of the Code
%Group relative position. That is better way to read dxf format. Don't fail if the
%polyline has arcs (budges), but yet don't read them. Don't read arcs as circles. Read
%properties (see case 'LINE' by examples of modifications). Group Codes and Associated
%Values is read in accurately way (accept associated values with space).
%
%Use cell2mat(cell(:,1)) to acquire geometry data in matrix,
%by example cell2mat(c_Cir(:,1))
clear all
close all
[filename, pathname]=uigetfile({'*.dxf', 'Select a DXF file'}, 'Select a DXF file',...
    'F:\DATI\D_data\dottorato\DATI\Antola\Outcrop_models\St_280116\Stazione 1_new\DXF');% <- MODIFY the PATH
% select path in which DXF of Circular Plane will be saved
uiwait(msgbox('Select folder for CircPlane.dxf'));
pathCircPlane=uigetdir(pathname, 'Select folder for CircPlane.dxf');
disp('########### START OF FITTING PROCCES ##########')
tic
%% Read file
fId = fopen(fullfile(pathname,filename));
c_ValAsoc = textscan(fId,'%d%s','Delimiter','\n');
fclose(fId);
% Code Group Matrix
m_GrCode = c_ValAsoc{1};
% Associated value String Cell
c_ValAsoc = c_ValAsoc{2};
%[m_GrCode,c_ValAsoc] = c_ValAsoc{:};

%% Entities
m_PosCero = find(m_GrCode==0);
%Is searched by (0,SECTION),(2,ENTITIES)
indInSecEnt = strmatch('ENTITIES',c_ValAsoc(m_PosCero(1:end-1)+1),'exact');
%(0,ENDSEC)
m_indFinSecEnt = strmatch('ENDSEC',c_ValAsoc(m_PosCero(indInSecEnt:end)),'exact');
% Entities Position
m_PosCero = m_PosCero(indInSecEnt:indInSecEnt-1+m_indFinSecEnt(1));
% Variable initiation
%accelerate?
     c_Line = cell(sum(strcmp('LINE',c_ValAsoc(m_PosCero))),2);
     c_Poly = cell(sum(strcmp('POLYLINE',c_ValAsoc(m_PosCero))),2);
     c_Cir = cell(sum(strcmp('CIRCLE',c_ValAsoc(m_PosCero))),2);
     c_Arc = cell(sum(strcmp('ARC',c_ValAsoc(m_PosCero))),2);
     c_Poi = cell(sum(strcmp('POINT',c_ValAsoc(m_PosCero))),2);
c_Line = cell(1,2);
c_Poly = cell(1,2);
c_Cir = cell(1,2);
c_Arc = cell(1,2);
c_Poi = cell(1,2);
%
iLine = 1;
iPoly = 1;
iCir = 1;
iArc = 1;
iPoi = 1;
% Loop on the Entities
for iEnt = 1:length(m_PosCero)-2
    m_GrCodeEnt = m_GrCode(m_PosCero(iEnt+1):m_PosCero(iEnt+2)-1);
    c_ValAsocEnt = c_ValAsoc(m_PosCero(iEnt+1):m_PosCero(iEnt+2)-1);
    nomEnt = c_ValAsocEnt{1};  %c_ValAsocEnt{m_PosCero(iEnt+1)}
    %In the entitie's name is assumed uppercase
    switch nomEnt
        case 'VERTEX'
            % (X,Y,Z) Position
            c_Poi{iPoi,1} = [str2double(f_ValGrCode(10,m_GrCodeEnt,c_ValAsocEnt)),...
                str2double(f_ValGrCode(20,m_GrCodeEnt,c_ValAsocEnt)),...
                str2double(f_ValGrCode(30,m_GrCodeEnt,c_ValAsocEnt))];
            % Layer
            c_Poi(iPoi,2) = f_ValGrCode(8,m_GrCodeEnt,c_ValAsocEnt);
            % Add properties
            %
            iPoi = iPoi+1;
            %case Add Entities
        
            % (X,Y,Z) vertexs
            
   
    end
end
%% create XYZNCloud matrix with X, Y and Z coordinate and Number of pointcloud
for i=1:numel(c_Poi)/2
   
    clear temp_name
    temp_name=char(c_Poi(i,2));
    Ncloud_c_Poi (i,1) = str2num(temp_name(10:end));
    
end
XYZNCloud=cell2mat(c_Poi(:,1));
XYZNCloud(:,4)=Ncloud_c_Poi(:,1);
%end
%% Fitting Plane step
plane=zeros(1,13);
for i=1 : max(XYZNCloud(:,4))%max(XYZNCloud(:,4)=number of planes
    
    temp_PC=zeros(1,3);
    
    for j=1:numel(XYZNCloud(:,4))%numel(XYZNCloud(:,4)=number of points (numel(XYZNCloud(:,4)== numel(XYZNCloud(:,1)== numel(XYZNCloud(:,2) == numel(XYZNCloud(:,3)
        % Write a temporary matrix containing points reffered to a plane
        if XYZNCloud(j,4)==i
            if temp_PC(1,1) == 0 && temp_PC(1,2) == 0 && temp_PC(1,3) == 0
                temp_PC(1,:)=XYZNCloud(j,1:3);
            else
                temp_PC(1+numel(temp_PC)/3,:)=XYZNCloud(j,1:3);
            end
        end
        
        
    end
    if (numel(temp_PC)/3)>2
        if plane(1,1) == 0 && plane(1,2) == 0 && plane(1,3) == 0
            plane(1,:)=fit3dplane_2(temp_PC);
        else
            plane(1+numel(plane)/13,:)=fit3dplane_2(temp_PC);
        end
    end
    
end
Tplane = table(plane(:,1),plane(:,2),plane(:,3),plane(:,4),plane(:,5),plane(:,6),plane(:,7),plane(:,8),plane(:,9),plane(:,10),plane(:,11),plane(:,12), plane(:,13));
Tplane.Properties.VariableNames = {'Dip' 'DipDirection' 'Radius' 'Xcenter' 'Ycenter' 'Zcenter' 'Nx' 'Ny' 'Nz' 'Mcoplanarity' 'Kcolinearity' 'MeanAbsError' 'StDevAbsError' };
tablefilenameTXT = (['Fit_',filename(1:end-4),'.txt']);
writetable(Tplane,fullfile(pathname,tablefilenameTXT));
tablefilenameXLSX = (['Fit_',filename(1:end-4),'.xlsx']);
writetable(Tplane,fullfile(pathname,tablefilenameXLSX));
toc
disp('########### END OF FITTING PROCCES ##########')


disp('########### START OF DISCONTINUITY PLOTTING AND SAVING PROCCESSES ##########')
tic
nplane=numel(plane(:,1));
% Plotting and saving dxf of circular plane
% figure (5)
% title(['Discontinuity plane, N = ', num2str(nplane)])
% %%--------------------Plot Circular Discontinuit Planes ------------------

N_set=zeros(nplane, 7); %matrix to storage number of plane for each set, every coloumn represent a set.

% Remember that radius of disconinuity its egual to diagonal of rectagnle
% formed by its values Horizontal and Vertical extent reached by
% CloudCompare
xyz(:,1)=plane(:,4);
xyz(:,2)=plane(:,5);
xyz(:,3)=plane(:,6);
Nxyz(:,1)=plane(:,7);
Nxyz(:,2)=plane(:,8);
Nxyz(:,3)=plane(:,9);
radius(:,1)=plane(:,3);
for i=1:nplane
    
    theta=0:0.1:2*pi; %Setting the spacing in which coordinate of the discontinuity disc are calculated
    v=null(Nxyz(i,:));% calculate vectors needed to plot the discs
    points=repmat(xyz(i,:)',1,size(theta,2))+radius(i)*(v(:,1)*cos(theta)+v(:,2)*sin(theta));%calculate points coordinate of the discs
    
    % Change name to the coordinate in a better and more resonable name
    X=points(1,:)';
    Y=points(2,:)';
    Z=points(3,:)';
    
    % With the below loop matlab could export DXF with different color in base
    % of set
    %for j=1:nclu
        %if idx3(i)==j
            % use DXFLib, Version 0.9.1 (Copyright (c) 2009-2011 Grzegorz Kwiatek)
            % to export DXF for each plane
            filename_mod=filename(1:end-4);
            Color = {'k','b','r','g','y',[.5 .6 .7],[.8 .2 .6]};
            % name is written in this way:
            % (number of set)_CircPlane_(number of row in CSV matrix)_(name of CSV file).dxf
            dxf_name=(['CircPlane',num2str(i),'_',num2str(filename_mod),'.dxf']);
            FID=dxf_open(pathCircPlane,dxf_name);
            FID = dxf_set(FID,'Color',Color{1});
            dxf_polyline(FID, X, Y, Z);
            dxf_close(FID); % end of DXF exportation
        %end
    %end
    
    %     % Uncomment for plot of discs of discontinuities
    %     hold on % hold on, box on and grid are useful for plotting
    %     grid on
    %     box on
    %     fill3(points(1,:),points(2,:),points(3,:),'b', 'FaceAlpha', 0.5);%Plot 3D disc for each plane
    %
    %     xlabel('x-axis (East)')
    %     ylabel('y-axis (Nord)')
    %     zlabel('z-axis (Elev)')
    %     hold on
    %     grid on
    %     box on
    
    %         %Uncomment if you wnat to see the normals of the planes
    %         quiver3(xyz(i,1), xyz(i,2), xyz(i,3), Nxyz(i,1), Nxyz(i,2), Nxyz(i,3),0.4)
    %         hold on
    %         grid on
    %         box on
    %
    
    
    
    
end
%-------------------------------------------------------------------------
toc
disp('########### END OF DISCONTINUITY PLOTTING AND SAVING PROCCESSES ##########')

%% Function to fit a 3D plane
function plane3d = fit3dplane(PC)
% This function is developed following the work of Jones et al. (2015)
% and Seers and Hodgetts (2016);
% REFERENCES:
%
% - Jones, R. R., Pearce, M. A., Jacquemyn, C., & Watson, F. E. (2015). 
%   Robust best-fit planes from geospatial data. Geosphere, 12(1), 196-202.
%
% - Seers, T. D., & Hodgetts, D. (2016). Probabilistic constraints on
%   structural lineament best fit plane precision obtained through
%   numerical analysis. Journal of Structural Geology, 82, 37-47.
% 
% INPUT
% 'PC' (point cloud) is a vector or a matrix N x 3 (N=number of points)
% that contain X Y Z coordinate of the points of the cloud

centroid = mean(PC); % calculate centroid of point cloud

%[U, S, W] = svd(PC);

C = cov(PC); % calculate the covariance matrix 'C'

[V, D] = eig(cov(PC)); % calculate eigenvector 'V' and eigenvalue 'D' of the
% covariance matrix 'C'


M = log(D(3,3)/D(1,1));% calculate vertex coplanarity 'M'
K = log(D(3,3)/D(2,2))/log(D(2,2)/D(1,1));% calculate vertex collinearity 
%'K'


N = V(:,1)';

% Normals have to be oriented upward!!! (Normals have to be horiented in same
% direction)

if N(1,3)<0 %if normal is oriented downward
    N(1,1)=-N(1,1);
    N(1,2)=-N(1,2);
    N(1,3)=-N(1,3);    
end
cosAlpha = N(1,1)/norm(N);
cosBeta = N(1,2)/norm(N);
cosGamma = N(1,3)/norm(N);

dip = 90 + rad2deg(asin(-cosGamma));

if cosAlpha > 0 && cosBeta > 0
dipdir = rad2deg(atan(cosAlpha/cosBeta));
elseif cosAlpha > 0 && cosBeta < 0
dipdir = 180 + rad2deg(atan(cosAlpha/cosBeta));
elseif cosAlpha < 0 && cosBeta < 0
dipdir = 180 + rad2deg(atan(cosAlpha/cosBeta));
elseif cosAlpha < 0 && cosBeta > 0
dipdir = 360 + rad2deg(atan(cosAlpha/cosBeta));
end

%calculate radius of plane = max distance point-centroid of fitted plane
for i=1:numel(PC(:,1))
pdist1 (i) =sqrt((centroid(1,1)-PC(i,1))^2 + (centroid(1,2)-PC(i,2))^2 + (centroid(1,3)-PC(i,3))^2);
end

[pM1, pI1] = max(pdist1);

for i=1:numel(PC(:,1))
    if i==pI1
        pdist2(i)=0;
    else
        pdist2 (i) =sqrt((PC(pI1,1)-PC(i,1))^2 + (PC(pI1,2)-PC(i,2))^2 + (PC(pI1,3)-PC(i,3))^2);
    end
end
[pM2, pI2] = max(pdist2);

centroid=mean([PC(pI1,:); PC(pI2,:)]);
radius =sqrt((PC(pI1,1)-PC(pI2,1))^2 + (PC(pI1,2)-PC(pI2,2))^2 + (PC(pI1,3)-PC(pI2,3))^2)/2;

%calculation of mean error and standard deviation of fitted plane
d = -sum(bsxfun(@times, N, bsxfun(@minus, centroid, PC)), 2);
mean_err=mean(abs(d));
stdev_err=std(abs(d));


plane3d = [dip, dipdir, radius, centroid, N, M, K, mean_err, stdev_err];
end
%% Functions to read and write dxf
% by Grzegorz Kwiatek 2011
function c_Val = f_ValGrCode(grCode,m_GrCode,c_ValAsoc)
c_Val = c_ValAsoc(m_GrCode==grCode);
end
function c_XData = f_XData(grCode,XDatNom,m_GrCode,c_ValAsoc)
m_PosXData = find(m_GrCode==1001);
if ~isempty(m_PosXData)
    indInXData = m_PosXData(strmatch(upper(XDatNom),c_ValAsoc(m_PosXData),'exact'));
    m_indFinXData = find(m_GrCode(indInXData+2:end)==1002)+indInXData+1;
    m_IndXData = indInXData+2:m_indFinXData(1)-1;
    c_XData = f_ValGrCode(grCode,m_GrCode(m_IndXData),c_ValAsoc(m_IndXData));
else
    c_XData = {[]};
end
end
function FID = dxf_open(pathname,fname) 
%DXF_OPEN Open DXF file.
%   Fid = DXF_OPEN(fname) opens DXF file for writing and writes the DXF 
%   file header. The function returns the matlab structure 'Fid' with 
%   various parameters used later by other drawing functions. One must 
%   use this structure in the subsequent calls to drawing routines.
%   
%   REMARKS
%     The assumed units are meters.
%   
%   See also DXF_CLOSE

%   Copyright 2010-2011 Grzegorz Kwiatek
%   $Revision: 1.1.4 $  $Date: 2011.11.15 $

try
  
  fid = fopen(fullfile(pathname,fname),'w');
  

  % Setup default values.
  FID.filename = (fname);
  FID.fid = fid;
  FID.show = false;
  FID.dump = true;
  
  FID.layer = 0;
  FID.color = 255;
  FID.textheight = 1;  
  FID.textrotation = 0;  
	FID.textthickness = 0;  
	FID.textalignment = 0;  
	FID.textvalignment = 0;  
	FID.textextrusion = [0 0 1];  
	FID.textobliqueangle = 0;  

  fprintf(fid,'0\nSECTION\n');
  fprintf(fid,'2\nHEADER\n');
  fprintf(fid,'9\n$ACADVER\n1\nAC1006\n'); % Default units: meters.
  fprintf(fid,'9\n$INSUNITS\n70\n6\n'); % Default units: meters.
  fprintf(fid,'0\nENDSEC\n');

  % Start entities or produce TABLE first.
  %if nargin ~= 1
  %  dxf_layertable(FID, varargin);
  %end
  
  % Dump entities section.
  fprintf(FID.fid,'0\nSECTION\n');
  fprintf(FID.fid,'2\nENTITIES\n');

catch exception
  if FID.fid >= 0
    fclose(FID.fid);
  end
  rethrow(exception);
end

end
function FID = dxf_set(FID, varargin)
%DXF_SET Set property.
%   DXF_SET(FID, 'PropertyName', 'PropertyValue', ...) sets up property
%   in DXF file.
%
%   Properties:
%   'Layer' - Layer number, followed by a scalar value in a range 0-256.
%   'Color' - Current color, followed either by a scalar value in a range
%             1-255 (ACI color) or MATLAB color vector [R G B], where 
%             R,G,B ranges 0-1. [0 0 0] is black, [1 1 1] is white, 
%             [1 0 0] is pure red, [0 1 0] is pure green and [0 0 1] is 
%             pure blue. Also ACI value 0 represent BYBLOCK and 256 represent
%             BYLAYER (currently both values are not used!)
%   
%   Copyright 2010-2011 Grzegorz Kwiatek.
%   $Revision: 1.1.2 $  $Date: 2011.09.06 $

p = inputParser;
p.addRequired('FID', @isstruct);
p.addParamValue('Layer', FID.layer, @(x) ischar(x) || (isnumeric(x) && x > 0) );  
p.addParamValue('Color', FID.color, @(x) (isscalar(x) && x>=0 && x<=256) || (isvector(x) && size(x,1)==1 && size(x,2) ==3) );  
p.addParamValue('TextHeight', FID.textheight, @(x)x >= 0);  
p.addParamValue('TextRotation', FID.textrotation, @(x)x >= 0);  
p.addParamValue('TextThickness', FID.textthickness, @(x)x >= 0);  
p.addParamValue('TextAlignment', FID.textalignment, @(x)x >= 0);  
p.addParamValue('TextVAlignment', FID.textvalignment, @(x)x >= 0);  
p.addParamValue('TextExtrusion', FID.textextrusion, @(x)(isvector(x) && size(x,1)==1 && size(x,2) ==3));  
p.addParamValue('TextObliqueAngle', FID.textobliqueangle, @(x)x >= 0);  
p.parse(FID, varargin{:});

FID.layer = p.Results.Layer;
if size(p.Results.Color,2) == 3
  FID.color = dxf_rgb2aci(p.Results.Color);
else
  FID.color = p.Results.Color;
end
FID.textheight = p.Results.TextHeight;  
FID.textrotation = p.Results.TextRotation;  
FID.textthickness = p.Results.TextThickness;  
FID.textalignment = p.Results.TextAlignment;  
FID.textvalignment = p.Results.TextVAlignment;  
FID.textextrusion = p.Results.TextExtrusion;  
FID.textobliqueangle = p.Results.TextObliqueAngle;  
  
end
function FID = dxf_polyline(FID, X, Y, Z)
%DXF_POLYLINE Store lines in DXF file.
%   DXF_POLYLINE(FID, X, Y, Z) writes polylines into DXF file. 
%   FID is valid DXFLIB structure created with DXF_OPEN routine. 
%   X, Y and Z are matrices of the same size containing coordinates of 
%   points for polylines. One separate line per column of matrix is 
%   drawn. 
%   
%   See also DXF_OPEN, DXF_CLOSE.

%   Copyright 2011 Grzegorz Kwiatek
%   $Revision: 1.2.0 $  $Date: 2011.11.17 $

try

  if FID.show
    line(X,Y,Z,'Color',dxf_aci2rgb(FID.color));
  end
  
  if FID.dump
    for i=1:size(X,2) % for each column!

      fprintf(FID.fid,'0\n');
      fprintf(FID.fid,'POLYLINE\n');

      dxf_print_layer(FID);
      fprintf(FID.fid,'66\n');  % entities follow (not necessary)
      fprintf(FID.fid,'1\n');  
      dxf_print_point(FID,0,0.0,0.0,0.0); % dummy point before vertices

      % 8 = This is a 3D polyline
      % 16 = This is a 3D polygon mesh
      % 32 = The polygon mesh is closed in the N direction
      % 64 = The polyline is a polyface mesh

      fprintf(FID.fid,'70\n');  
      fprintf(FID.fid,'8\n');   
      dxf_print_vertex(FID, [X(:,i) Y(:,i) Z(:,i)],32);
      dxf_print_seqend(FID);
    end % for loop
  end

catch exception
  if FID.fid >= 0
    fclose(FID.fid);
  end
  rethrow(exception);
end
end
function dxf_print_point(FID,pointno,x,y,z)
%DXF_PRINT_POINT Dump entity properties.
%   Internal function: Not usable directly.
%
%   Copyright 2011 Grzegorz Kwiatek
%   $Revision: 1.0.0 $  $Date: 2011.08.25 $%

try
  fprintf(FID.fid,'1%1d\n%1.8g\n2%1d\n%1.8g\n3%1d\n%1.8g\n',pointno,x,pointno,y,pointno,z);
catch exception
  if FID.fid >= 0
    fclose(FID.fid);
  end
  rethrow(exception);
end
end
function dxf_print_layer(FID)
%DXF_PRINT_LAYER Dump entity properties.
%   Internal function: Not usable directly.
%
%   Copyright 2010-2011 Grzegorz Kwiatek.
%   $Revision: 1.0.3 $  $Date: 2011.08.25 $%

try
  fprintf(FID.fid,'8\n');  % layer, default 0
  if isnumeric(FID.layer)
    fprintf(FID.fid,'%d\n', FID.layer);  
  else
    fprintf(FID.fid,'%s\n', FID.layer);  
  end
  
  % Dump color.
  fprintf(FID.fid,'62\n%d\n',FID.color); % color group code.

catch exception
  if FID.fid >= 0
    fclose(FID.fid);
  end
  rethrow(exception);
end
end
function FID = dxf_print_seqend(FID)
%DXF_PRINT_SEQEND Dump entity properties.
%   Internal function: Not usable directly.
%
%   Copyright 2011 Grzegorz Kwiatek.
%   $Revision: 1.0.0 $  $Date: 2011.08.25 $%

try 
  fprintf(FID.fid,'0\nSEQEND\n8\n%d\n',FID.layer);
catch exception
  if FID.fid >= 0
    fclose(FID.fid);
  end
  rethrow(exception);
end  
end
function dxf_print_vertex(FID, VERTICES, vertextype, varargin)
%DXF_PRINT_VERTEX Dump entity properties.
%   Internal function: Not usable directly.
%
%   Copyright 2011 Grzegorz Kwiatek.
%   $Revision: 1.0.0 $  $Date: 2011.08.25 $%

try
  if nargin == 4
    % only per-face colouring is accepted.
    COLOURS = varargin{1};
    % if COLOURS is (N x 3) matrix, reduce it to ACI colour index table (N x 1).
    if size(COLOURS,2) == 3
      COLOURS = dxf_rgb2aci(COLOURS);
    end
  end
  
  for i=1:size(VERTICES,1)
    
    fprintf(FID.fid,'0\n');
    fprintf(FID.fid,'VERTEX\n');

    if vertextype ~= 128
      % = 128+64 dump as polyface mesh vertex (called by dxf_polymesh).
      % = 32 dump as polyline vertex (called by dxf_polyline)
      % Per-vertex colouring does not work for DXF files.
      dxf_print_layer(FID);
      dxf_print_point(FID,0,VERTICES(i,1),VERTICES(i,2),VERTICES(i,3));
    else
      % =128 - dump faces of polyface mesh.
      % Setup different per-face colour depending on COLOURS matrix.
      currentcolor = FID.color;
      if nargin == 4
        FID.color = COLOURS(i);
      end
      dxf_print_layer(FID);
      FID.color = currentcolor;
      dxf_print_point(FID,0,0,0,0);
    end
    
    % 70 Vertex flags:
    % 1 = Extra vertex created by curve-fitting
    % 2 = Curve-fit tangent defined for this vertex. A curve-fit tangent direction of 0 may be 
    %     omitted from DXF output but is significant if this bit is set
    % 4 = Not used
    % 8 = Spline vertex created by spline-fitting
    % 16 = Spline frame control point
    % 32 = 3D polyline vertex
    % 64 = 3D polygon mesh
    % 128 = Polyface mesh vertex

    fprintf(FID.fid,'70\n');  
    fprintf(FID.fid,'%d\n',vertextype); % layer  
    
    % =128, dump faces of polyface mesh.
    if vertextype == 128
      fprintf(FID.fid,'71\n');  
      fprintf(FID.fid,'%d\n',VERTICES(i,1)); % layer  
      fprintf(FID.fid,'72\n');  
      fprintf(FID.fid,'%d\n',VERTICES(i,2)); % layer  
      fprintf(FID.fid,'73\n');  
      fprintf(FID.fid,'%d\n',VERTICES(i,3)); % layer  
      if size(VERTICES,2) == 4
        fprintf(FID.fid,'74\n');  
        fprintf(FID.fid,'%d\n',VERTICES(i,4)); % layer  
      else
        fprintf(FID.fid,'74\n');  
        fprintf(FID.fid,'0\n'); % layer  
      end
    end
  end
  
catch exception
  if FID.fid >= 0
    fclose(FID.fid);
  end
  rethrow(exception);
end  
end
function dxf_close(FID)
%DXF_CLOSE Close DXF file.
%   DXF_CLOSE(FID) closes the DXF file opened with DXF_OPEN command.
%   
%   See also DXF_OPEN

%   Copyright 2010-2011 Grzegorz Kwiatek
%   $Revision: 1.1.2 $  $Date: 2011.08.25 $

try
  fprintf(FID.fid,'0\nENDSEC\n0\nEOF\n');
  fclose(FID.fid);
catch exception
  if FID.fid >= 0
    fclose(FID.fid);
  end
  rethrow(exception);
end
end
