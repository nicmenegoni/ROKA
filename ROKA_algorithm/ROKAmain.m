%% ROck slope Kinematic Analysis (ROKA)
close all; clear variables; tic
disp('########### ROKA algorithm has been launched ############')% Display a message onto the command window
%% 0) DEFINITION OF SOME PARAMETERS
userSPr=0.1;%radius of the scan-volume and scan-plane

disccutoff=2*userSPr*0.9;%define the cutoff value of the intersection between discontinuities and scan-plane used to perform or not the KA (
%the default value is 75% of the diameter of the scan-plane

uselatlimits=1;% Define if you use the lateral limits or not, and give a value (YES=1 and NO=0)

latlimits=20;%Define lateral limits (Godmann report 20°)

frictionangle=30;%Define the friction angle (commonly 30°)

export_AllIntersection=0;%Define if you want to export all the disconinuity intersection (also the non-crtitical)

% % prompt = {['Please, define the scan-radius: ']};
% %     dlg_title = 'Define the scan-radius';
% % userSPr=inputdlg(prompt);%define the size of the scan-radius
% % userSPr=str2num(userSPr{1});

%% 1) IMPORT FRACTURE GEOMETRY DATA
% This Matlab code works using a xlsx file in which xyz (center), Nxyz
% (orientation) and radius (dimension) of the discontinuities are stored.
warning ('off', 'all')%warnings disabled 
uiwait(msgbox('Select fracture geometry data (XLSX file) to load'));
[filename, pathname] = uigetfile({'*.xlsx', 'Select fracture geometry data (XLSX file) to load'},'Select fracture geometry data (XLSX file) to load',...
    'Z:\Menegoni\Temporanea\ROKA_files_test');% <- MODIFY the PATH

Fracdata=readtable(fullfile(pathname, filename));%read  fracture data stored
%in the previously defined XLSX file.

%% 2) IMPORT FRACTURE SET DATA
uiwait(msgbox('Select Fracture set data (XLSX file) to load'));
[filenameset, pathnameset] = uigetfile({'*.xlsx', 'Select Fracture set data (XLSX file) to load'},...
    'Select Fracture set data (XLSX file) to load',...
    pathname);
Fracset=readtable(fullfile(pathnameset, filenameset));%read fracture set
%stored in the previously defined XLSX file.

%% 3) DEFINE DISCONTINUITES VARIABLES FOR CALCULATION
nplane=numel(Fracdata.Dip);%number of fractures/discontinuity planes)
xyz=[Fracdata.Xcenter(:),Fracdata.Ycenter(:),Fracdata.Zcenter(:)];
Nxyz=[Fracdata.Nx(:),Fracdata.Ny(:),Fracdata.Nz(:)];
dipdir =  Fracdata.DipDirection(:);%Fracture dip direction
dip =  Fracdata.Dip(:);%fracture dip angle
radius = Fracdata.Radius(:);
Set = Fracset.Set;% Fracture set (defined in the imported xlsx file)
Set(isnan(Set)) = max(Set)+1;% Definition of the random set (not defined in the imported xslx file)
Color = {'k','b','r','g','y',[.5 .6 .7],[.8 .2 .6]};% define color for set from 1 to 7
%For what concern the color it is ossible to increase the number of sets


%% 6) IMPORT AND READ  THE POINTCLOUD
uiwait(msgbox('Select PointCloud  (TXT file) to load'));
[PCfn,PCpn]=uigetfile({'*.txt', 'Select PointCloud  (TXT file) to load'},'Select PointCloud  (TXT file) to load',...
    pathname);%Import the path point cloud txt file
PC=importdata(fullfile(PCpn,PCfn));%import and read the pont cloud
ptCloud=pointCloud(PC(:,1:3));
PCnormals=PC(:,4:6);
nPC=length(PC(:,1));
%% Calculate Disconunuity intersections
%Creating some matrix containing the dip and dip direction of the
%disconinuity planes and the trend and plung of the disconinuity
%intersections
%Intersection definitions
disp(['Starting to calculate all the disconinuities intersections at ', num2str(toc/60),' minutes'])
[Int_data]=intersectionCalculator(pathname,Fracdata,Fracset,export_AllIntersection);%matrix contating all the information about the disconinuity intersections
disp(['End of the calculation of the disconinuities intersections at ', num2str(toc/60),' minutes'])
num_Int=numel(Int_data.Trend);

Trend=Int_data.Trend;
Plunge=Int_data.Plunge;
Int_Sets(:,1)=Int_data.Set_i;
Int_Sets(:,2)=Int_data.Set_j;
pole_Trend = zeros(num_Int,1);%pre-allocating the intersection trend vector
pole_Plunge = zeros(num_Int,1);%pre-allocating the intersection plunge vector
for i= 1 : num_Int %calculate discontinuity line pole (dip and dip direction)
    %this is an trick to plot plunge and trend line vector
    pole_Plunge(i,1) = 90 - Plunge(i);%discontinuity pole dip
    if Trend(i) < 180%discontinuity pole dip direction
        pole_Trend(i,1) = Trend(i) + 180;
    else
        pole_Trend(i,1) = Trend(i) - 180;
    end
end
%Disconinuity definition
Dip=Fracdata.Dip(:);%read Dip value of discontinuities
DipDir=Fracdata.DipDirection(:);%read Dip Direction value of discontinuities
Set=Fracset.Set(:);%read Set value of discontinuities (random values are defined as NaN as default)
Set(isnan(Set))=0;%change Set value of random discontinuities from NaN to 0
Set_name=unique(Set);
nSet=numel(unique(Set));
pole_Dip=zeros(nplane,1);
pole_DipDir=zeros(nplane,1);
for i= 1 : nplane %calculate discontinuity line pole (dip and dip direction)
    %this is an trick to plot plunge and trend line vector
    pole_Dip(i,1) = 90 - Dip(i);%discontinuity pole dip
    if DipDir(i) < 180%discontinuity pole dip direction
        pole_DipDir(i,1) = DipDir(i) + 180;
    else
        pole_DipDir(i,1) = DipDir(i) - 180;
    end
end

%% Pre calculation variable definition
CDM=zeros(nplane,2);%Critical Disconinuity Matrix (CDM): record the possible
% mode of failure of the disconuities as 0 (no failure) or 1 (possible
% failure). First value Planar Sliding (PS) second Flexural Toppling (FT)
CDMpc=cell(length(PC(:,1)),2);%this variable is neeeded for the parfor 
% version of the code, becase it record onto the PC points the index of the
% critical disconinuity for planar sliding (CDMpc{i,1}) and flexural
% toppling (CDMpc{i,2})

nCritic = zeros(length(PC(:,1)),5);%for every scan plane calulated at each
%point, the number of crticial discontunuity is stored in this order: PS,
%FT, WS, DT, OT

CIM=zeros(num_Int,3);%Critical Intersection Matrix (CIM): record the possible
% mode of failure of the disconuities intersectionas 0 (no failure) or 1
% (possible failure). First value Wedge Sliding (WS) second Oblique Toppling (OT)
% third value Direct Toppling (DT)
CIMpc=cell(length(PC(:,1)),3);%this variable is neeeded for the parfor 
% version of the code, becase it record onto the PC points the index of the
% critical disconinuity for wedge sliding (CIMpc{i,1}) and direct toppling
% (CIMpc{i,2}) and oblique toppling (CIMpc{i3})

PCDicsInters=zeros(length(PC(:,1)),1);% record of the number of disconunity
% that intersect the scan planes;

PCIntInters=zeros(length(PC(:,1)),1);%record of the number of disconunity
% intersections that intersect the scan planes;

slopeDipdirDip=zeros(length(PC(:,1)),2);%record of the orientation in
% Dip direction and Dip of each scan plane
overhanging=zeros(length(PC(:,1)),1);%record if the slope is overhanging
SPr=zeros(length(PC(:,1)),1);%record the 'radius' of the scan plane
%% Kinematic Analysis Loop
disp(['Starting the kinematic analysis at ', num2str(toc/60),' minutes'])
tic
reqKnn=requiredKNNcalc(PC(:,1:3),ptCloud,userSPr);
disp(['Time to calculate the best Knn points number ', num2str(toc/60), ' min.'])
tic
parfor ptloop= 1: nPC %define a loop foe every point of the cloud
    %clearvars idxK ptdist SPpoints A_SP SP_N SP_center v1 v2 v3 v4 nCritic_line cosAlpha cosBeta cosGamma slope_dipdir slope_dip
    
    [idxK, ptdist] = findNearestNeighbors(ptCloud,PC(ptloop,1:3),2*reqKnn);%calculate the index of the Knn points.
    % the idxK indicates the position of the Knn points inside the PC matrix;
    % the ptdist indicates the distance of the points from the 'ptloop' one
    idxK(ptdist>userSPr)=[];
    ptdist(ptdist>userSPr)=[];
    SPpoints = PC(idxK,1:3);%define the points belonging to the scan area (based on Knn algorithm)
    
    %Calculating the center, normals, radius and attitude (dipdir and dip) of the scanplane
    SP_center = mean(PC(idxK,1:3));%scanplane center
    if length(idxK(:,1))>1 % if the analyzed point has one or more nearest neighbour points closer that fall in the scan radius, the code perfomr the consequent analysis
        % This if loop is necessary becasue sometimes, sparse points could 
        % be present and they have a distance from the other point of the 
        % pointcloud longer than the defined scan radius. Therefore, it is 
        % necessary to define exclude these points from calculation
        SP_N = mean(PC(idxK,4:6));%scanplane normals
        SPr(ptloop,1) = max(ptdist);%radius equal to the most distance filtered Knn points (always lower than userSP)
        [sdipdir,sdip]= normal2attitude(SP_N);%attitude of the scanplane
        slopeDipdirDip(ptloop,:)=[sdipdir,sdip];
        if SP_N (1,3)<0
            overhanging(ptloop,1)=1;
        else
            overhanging(ptloop,1)=0;
        end
        %calculating intersection between discontinties and finite scan plane
        [in_out] = intSPdisc(SP_center, SP_N, SPr(ptloop,1), xyz, Nxyz, radius,disccutoff);
        PCDicsInters(ptloop)=sum(in_out>0);
        [int_in_out] = intSPinters([Int_data.x1,Int_data.y1,Int_data.z1],...
            [Int_data.x2,Int_data.y2,Int_data.z2],...
            Int_data.Length_Intersection,SP_center, SP_N, SPr(ptloop,1));
        
        %[int_in_out] = lineinsidewindow(num_Int,v1,v2,v3,v4,A_SP,SP_N,SP_center,Int_data);
        
        PCIntInters(ptloop)=sum(int_in_out>0);
        idxPS=[];idxFT=[];idxWS=[];idxDT=[];idxOT=[];
        if (sum(in_out)>0 || sum(int_in_out)>0) && overhanging(ptloop,1)==0
            [idxPS,idxFT,idxWS,idxDT,idxOT]=KinematicAnalysis_parfor(Fracdata,Int_data,uselatlimits,latlimits,frictionangle,in_out,int_in_out,slopeDipdirDip(ptloop,:));
        elseif (sum(in_out)>0 || sum(int_in_out)>0) && overhanging(ptloop,1)==1
            [idxPS,idxFT,idxWS,idxDT,idxOT]=KinematicAnalysis_parfor(Fracdata,Int_data,uselatlimits,latlimits,frictionangle,in_out,int_in_out,[sdipdir,90]);
        end
        
        
        CDMpc(ptloop,:)={idxPS,idxFT};
        CIMpc(ptloop,:)={idxWS,idxDT,idxOT};
        nCritic(ptloop,:)=[numel(idxPS),numel(idxFT),numel(idxWS),numel(idxDT),numel(idxOT)];
    else %if the considered point has no nearest neighbour that fall inside the scan radius, the code does not perform the KA analysis
        overhanging(ptloop,1)=NaN;
        PCIntInters(ptloop)=0;
        PCDicsInters(ptloop)=0;
        slopeDipdirDip(ptloop,:)=[0,0];
        idxPS=[];idxFT=[];idxWS=[];idxDT=[];idxOT=[];
        CDMpc(ptloop,:)={idxPS,idxFT};
        CIMpc(ptloop,:)={idxWS,idxDT,idxOT};
        nCritic(ptloop,:)=[numel(idxPS),numel(idxFT),numel(idxWS),numel(idxDT),numel(idxOT)];
    end
end
toc
%Previous to the result point cloud exportation, it removes the NaN value
PC(isnan(overhanging),:)=[];
slopeDipdirDip(isnan(overhanging),:)=[];
PCDicsInters(isnan(overhanging),:)=[];
PCIntInters(isnan(overhanging),:)=[];
nCritic(isnan(overhanging),:)=[];
CDMpc(isnan(overhanging),:)=[];
CIMpc(isnan(overhanging),:)=[];
overhanging(isnan(overhanging),:)=[];

%Calculating the number of positive crtical test for each disconinuity plane
%and intersection
PCidPS=vertcat(CDMpc{:,1});
PCidFT=vertcat(CDMpc{:,2});
PCidWS=vertcat(CIMpc{:,1});
PCidDT=vertcat(CIMpc{:,2});
PCidOT=vertcat(CIMpc{:,3});
idPS=unique(PCidPS);
idFT=unique(PCidFT);
idWS=unique(PCidWS);
idDT=unique(PCidDT);
idOT=unique(PCidOT);
for i= 1 : length(idPS)
CDM(idPS(i),1) = sum(PCidPS == idPS(i));
end
for i= 1 : length(idFT)
CDM(idFT(i),2) = sum(PCidFT == idFT(i));
end
for i= 1 : length(idWS)
CIM(idWS(i),1) = sum(PCidWS == idWS(i));
end
for i= 1 : length(idDT)
CIM(idDT(i),2) = sum(PCidDT == idDT(i));
end
for i= 1 : length(idOT)
CIM(idOT(i),3) = sum(PCidOT == idOT(i));
end

%% Calculate the 'percentual' criticity of the discontinuity planes and intersections
CDMperc=zeros(size(CDM));%Initializing the matrix where the criticalk value is expressed in percentual
CIMperc=zeros(size(CIM));%Initializing the matrix where the criticalk value is expressed in percentual

%calculating the percentual of critical value, where the 100% refer to the
%critical trace of the discontinuity with max exposure
CDMperc(:,1) = (CDM(:,1)/max(CDM(:,1)))*100;
CDMperc(:,2) = (CDM(:,2)/max(CDM(:,2)))*100;
CIMperc(:,1) = (CIM(:,1)/max(CIM(:,1)))*100;
CIMperc(:,2) = (CIM(:,2)/max(CIM(:,2)))*100;
CIMperc(:,3) = (CIM(:,3)/max(CIM(:,3)))*100;

%% Calculate the 'crtiticity' onto the pointcloud

PCcritic = zeros(length(PC(:,1)),5);
%for each point gives the max 'critical value' of the interescted and unstable discontinuity planes and interesctions (if present)
for i = 1 : length(PC(:,1))
    %max for Planar Sliding
    if sum(CDMpc{i,1})>0
    PCcritic(i,1)=max(CDMperc(CDMpc{i,1},1));
    end
    
    %max for Flexural Toppling
    if sum(CDMpc{i,2})>0
        PCcritic(i,2)=max(CDMperc(CDMpc{i,2},2));
    end
    %max for Wedge Sliding
    if sum(CIMpc{i,1})>0
    PCcritic(i,3)=max(CIMperc(CIMpc{i,1},1));
    end
    
    %max for Direct Toppling
    if sum(CIMpc{i,2})>0
    PCcritic(i,4)=max(CIMperc(CIMpc{i,2},2));
    end
    
    %max for Oblique Toppling
    if sum(CIMpc{i,3})>0
    PCcritic(i,5)=max(CIMperc(CIMpc{i,3},3));
    end
end


%% Exporting Critical Discontinuites and Intersections
expAnalysis(pathname,Fracdata,Fracset,Int_data,CDM,CIM,PC(:,1:3),PCcritic,uselatlimits,slopeDipdirDip,PCDicsInters,PCIntInters,overhanging,CDMperc,CIMperc)
%expFracture(


disp('#######################')
disp(['END OF TOTAL PROCESSES at ', num2str(toc/60),' minutes'])
disp('#######################')

disp('#####STATISTICS########')
disp(['ScanPlane radius (SPr)------------------'])
disp(['mean)', num2str(mean(SPr)),' meters'])
disp(['max)', num2str(max(SPr)),' meters'])
disp(['min)', num2str(min(SPr)),' meters'])


% END - Audio message
WarnWave = [sin(1:.6:400), sin(1:.7:400), sin(1:.4:400)];
Audio = audioplayer(WarnWave, 22050);
play(Audio);

