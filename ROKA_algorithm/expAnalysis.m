function expAnalysis(pathname,Fracdata,Fracset,Int_data,CDM,CIM,PC,PCcritic,uselatlimits,slopeDipdirDip,PCDicsInters,PCIntInters,overhanging,CDMperc,CIMperc)
%% LOG
%modified 2020/04/27--> an 'exist' function is added before the exportation
%of the dxf in order to avoid the multiple exporting of the crticial
%discontinuity planes and intersetcions
%Modified 2020/05/19--> the matrix CDM and CIM are added in order to
%simplify calculation, and also a process toexport the point cloiud with
%the number of ceritical disconinuities was added;
%CDM matrix contains the critical value of the disconinuity in this order:
%                   1)PS; 2) FT;
%CIM matrix contains the critical value of the intersections in this order:
%                   1)WS; 2) DT; 3) OT
%Modified 2020/11/03--> the crtical value are expressed in percentual
%respect to the max disconinuity trace
%Removed the Set folder dutring the saving


%% FUNCTION CODE
warning('off')
Color = {'k','b','r','g','y',[.5 .6 .7],[.8 .2 .6]};
disp('Strarting to export the critical discontinuities')
for i=1:length(Fracdata{:,1}) %export critical plane
    clearvars FID
    %% export PlanarSliding critical CDM(:,1)
    if CDM(i,1)>0
        theta=0:0.1:2*pi; %Setting the spacing in which coordinate of the discontinuity disc are calculated
        v=null(Fracdata{i,7:9});% calculate vectors needed to plot the discs
        points=repmat(Fracdata{i,4:6}',1,size(theta,2))+Fracdata{i,3}*(v(:,1)*cos(theta)+v(:,2)*sin(theta));%calculate points coordinate of the discs
        
        % Change name to the coordinate in a better and more resonable name
        X=points(1,:)'; Y=points(2,:)'; Z=points(3,:)';
        % use DXFLib, Version 0.9.1 (Copyright (c) 2009-2011 Grzegorz Kwiatek)
        % to export DXF for each plane
        
        if ispc
            if uselatlimits==1
                mkdir (pathname,'KinAnalysis\Disc\PlanarSlidingLatLimits')
                dxf_name=([num2str(CDMperc(i,1),'%03.f'),'_CriticDisc_',num2str(i),'.dxf']);
                FID=dxf_open(fullfile(pathname,'KinAnalysis\Disc\PlanarSlidingLatLimits\'),dxf_name);
                FID = dxf_set(FID,'Color',Color{5});
                dxf_polyline(FID, X, Y, Z);
                dxf_close(FID); % end of DXF exportation
                
            elseif uselatlimits==0
                
                mkdir (pathname,'KinAnalysis\Disc\PlanarSlidingNoLatLimits')
                dxf_name=([num2str(CDMperc(i,1),'%03.f'),'_CriticDisc_',num2str(i),'.dxf']);
                FID=dxf_open(fullfile(pathname,'KinAnalysis\Disc\PlanarSlidingNoLatLimits\'),dxf_name);
                FID = dxf_set(FID,'Color',Color{3});
                dxf_polyline(FID, X, Y, Z);
                dxf_close(FID); % end of DXF exportation
                
            end
        elseif ismac
            if uselatlimits==1
                mkdir (pathname,'KinAnalysis/Disc/PlanarSlidingLatLimits')
                dxf_name=([num2str(CDMperc(i,1),'%03.f'),'_CriticDisc_',num2str(i),'.dxf']);
                FID=dxf_open(fullfile(pathname,'KinAnalysis/Disc/PlanarSlidingLatLimits/'),dxf_name);
                FID = dxf_set(FID,'Color',Color{5});
                dxf_polyline(FID, X, Y, Z);
                dxf_close(FID); % end of DXF exportation
                
            elseif uselatlimits==0
                
                mkdir (pathname,'KinAnalysis/Disc/PlanarSlidingNoLatLimits')
                dxf_name=([num2str(CDMperc(i,1),'%03.f'),'_CriticDisc_',num2str(i),'.dxf']);
                FID=dxf_open(fullfile(pathname,'KinAnalysis/Disc/PlanarSlidingNoLatLimits/'),dxf_name);
                FID = dxf_set(FID,'Color',Color{3});
                dxf_polyline(FID, X, Y, Z);
                dxf_close(FID); % end of DXF exportation
            end
        end
    end
    
    %% Export Flexural Toppling critical discontintuies (CDM(:,2)
    if CDM(i,2)>0
        theta=0:0.1:2*pi; %Setting the spacing in which coordinate of the discontinuity disc are calculated
        v=null(Fracdata{i,7:9});% calculate vectors needed to plot the discs
        points=repmat(Fracdata{i,4:6}',1,size(theta,2))+Fracdata{i,3}*(v(:,1)*cos(theta)+v(:,2)*sin(theta));%calculate points coordinate of the discs
        % Change name to the coordinate in a better and more resonable name
        X=points(1,:)'; Y=points(2,:)'; Z=points(3,:)';
        % use DXFLib, Version 0.9.1 (Copyright (c) 2009-2011 Grzegorz Kwiatek)
        % to export DXF for each plane
        if ispc
            
            mkdir (pathname,'KinAnalysis\Disc\FlexuralToppling')
            dxf_name=([num2str(CDMperc(i,2),'%03.f'),'_CriticDisc_',num2str(i),'.dxf']);
            FID=dxf_open(fullfile(pathname,'KinAnalysis\Disc\FlexuralToppling\'),dxf_name);
            FID = dxf_set(FID,'Color',Color{5});
            dxf_polyline(FID, X, Y, Z);
            dxf_close(FID); % end of DXF exportation'
            
        elseif ismac
            
            mkdir (pathname,'KinAnalysis/Disc/FlexuralToppling')
            dxf_name=([num2str(CDMperc(i,2),'%03.f'),'_CriticDisc_',num2str(i),'.dxf']);
            FID=dxf_open(fullfile(pathname,'KinAnalysis/Disc/FlexuralToppling/'),dxf_name);
            FID = dxf_set(FID,'Color',Color{5});
            dxf_polyline(FID, X, Y, Z);
            dxf_close(FID); % end of DXF exportation
            
        end
    end
    
end
for i=1:length(Int_data{:,1}) %export critical Intersection
    %% export WedgeSliding critical: CIM(:,1)
    if CIM(i,1)>0
        % Change name to the coordinate in a better and more resonable name
        X=[Int_data.x1(i);Int_data.x2(i)]; Y=[Int_data.y1(i);Int_data.y2(i)]; Z=[Int_data.z1(i);Int_data.z2(i)];
        % use DXFLib, Version 0.9.1 (Copyright (c) 2009-2011 Grzegorz Kwiatek)
        % to export DXF for each plane
        if ispc
            if uselatlimits==1
                
                mkdir (pathname,'KinAnalysis\Inters\WedgeSlidingLatLimits')
                dxf_name=([num2str(CIMperc(i,1),'%03.f'),'_Inters_',num2str(i),'_Discs',num2str(Int_data.Nplane_i(i)),'_',num2str(Int_data.Nplane_j(i)),'.dxf']);
                FID=dxf_open(fullfile(pathname,'KinAnalysis\Inters\WedgeSlidingLatLimits\'),dxf_name);
                FID = dxf_set(FID,'Color',Color{5});
                dxf_polyline(FID, X, Y, Z);
                dxf_close(FID); % end of DXF exportation
                
            elseif uselatlimits==0
                
                mkdir (pathname,'KinAnalysis\Inters\WedgeSlidingNoLatLimits')
                dxf_name=([num2str(CIMperc(i,1),'%03.f'),'_Inters_',num2str(i),'_Discs',num2str(Int_data.Nplane_i(i)),'_',num2str(Int_data.Nplane_j(i)),'.dxf']);
                FID=dxf_open(fullfile(pathname,'KinAnalysis\Inters\WedgeSlidingNoLatLimits\'),dxf_name);
                FID = dxf_set(FID,'Color',Color{3});
                dxf_polyline(FID, X, Y, Z);
                dxf_close(FID); % end of DXF exportation
                
            end
        elseif ismac
            if uselatlimits==1
                
                mkdir (pathname,'KinAnalysis/Inters/WedgeSlidingLatLimits')
                dxf_name=([num2str(CIMperc(i,1),'%03.f'),'_Inters_',num2str(i),'_Discs',num2str(Int_data.Nplane_i(i)),'_',num2str(Int_data.Nplane_j(i)),'.dxf']);
                FID=dxf_open(fullfile(pathname,'KinAnalysis/Inters/WedgeSlidingLatLimits/'),dxf_name);
                FID = dxf_set(FID,'Color',Color{5});
                dxf_polyline(FID, X, Y, Z);
                dxf_close(FID); % end of DXF exportation
                
            elseif uselatlimits==0
                
                mkdir (pathname,'KinAnalysis/Inters/WedgeSlidingNoLatLimits')
                dxf_name=([num2str(CIMperc(i,1),'%03.f'),'_Inters_',num2str(i),'_Discs',num2str(Int_data.Nplane_i(i)),'_',num2str(Int_data.Nplane_j(i)),'.dxf']);
                FID=dxf_open(fullfile(pathname,'KinAnalysis/Inters/WedgeSlidingNoLatLimits/'),dxf_name);
                FID = dxf_set(FID,'Color',Color{3});
                dxf_polyline(FID, X, Y, Z);
                dxf_close(FID); % end of DXF exportation
                
            end
        end
    end
    %% export DirectToppling critical Intersections: CIM(:,2)
    if CIM(i,2)>0
        
        % Change name to the coordinate in a better and more resonable name
        X=[Int_data.x1(i);Int_data.x2(i)]; Y=[Int_data.y1(i);Int_data.y2(i)]; Z=[Int_data.z1(i);Int_data.z2(i)];
        % use DXFLib, Version 0.9.1 (Copyright (c) 2009-2011 Grzegorz Kwiatek)
        % to export DXF for each plane
        
        if ispc
            if uselatlimits==1
                
                mkdir (pathname,'KinAnalysis\Inters\DirectTopplingLatLimits\')
                dxf_name=([num2str(CIMperc(i,2),'%03.f'),'_Inters_',num2str(i),'_Discs',num2str(Int_data.Nplane_i(i)),'_',num2str(Int_data.Nplane_j(i)),'.dxf']);
                FID=dxf_open(fullfile(pathname,'KinAnalysis\Inters\DirectTopplingLatLimits\'),dxf_name);
                FID = dxf_set(FID,'Color',Color{5});
                dxf_polyline(FID, X, Y, Z);
                dxf_close(FID); % end of DXF exportation
                
            elseif uselatlimits==0
                
                mkdir (pathname,'KinAnalysis\Inters\DirectTopplingNoLatLimits\')
                dxf_name=([num2str(CIMperc(i,2),'%03.f'),'_Inters_',num2str(i),'_Discs',num2str(Int_data.Nplane_i(i)),'_',num2str(Int_data.Nplane_j(i)),'.dxf']);
                FID=dxf_open(fullfile(pathname,'KinAnalysis\Inters\DirectTopplingNoLatLimits\'),dxf_name);
                FID = dxf_set(FID,'Color',Color{3});
                dxf_polyline(FID, X, Y, Z);
                dxf_close(FID); % end of DXF exportation
                
            end
        elseif ismac
            if uselatlimits==1
                
                mkdir (pathname,'KinAnalysis/Inters/DirectTopplingLatLimits/')
                dxf_name=([num2str(CIMperc(i,2),'%03.f'),'_Inters_',num2str(i),'_Discs',num2str(Int_data.Nplane_i(i)),'_',num2str(Int_data.Nplane_j(i)),'.dxf']);
                FID=dxf_open(fullfile(pathname,'KinAnalysis/Inters/DirectTopplingLatLimits/'),dxf_name);
                FID = dxf_set(FID,'Color',Color{5});
                dxf_polyline(FID, X, Y, Z);
                dxf_close(FID); % end of DXF exportation
                
            elseif uselatlimits==0
                
                mkdir (pathname,'KinAnalysis/Inters/DirectTopplingNoLatLimits/')
                dxf_name=([num2str(CIMperc(i,2),'%03.f'),'_Inters_',num2str(i),'_Discs',num2str(Int_data.Nplane_i(i)),'_',num2str(Int_data.Nplane_j(i)),'.dxf']);
                FID=dxf_open(fullfile(pathname,'KinAnalysis/Inters/DirectTopplingNoLatLimits/'),dxf_name);
                FID = dxf_set(FID,'Color',Color{3});
                dxf_polyline(FID, X, Y, Z);
                dxf_close(FID); % end of DXF exportation
                
            end
        end
    end
    %% export ObliqueToppling critical Intersections: CIM(:,3)
    if CIM(i,3)>0
        
        % Change name to the coordinate in a better and more resonable name
        X=[Int_data.x1(i);Int_data.x2(i)]; Y=[Int_data.y1(i);Int_data.y2(i)]; Z=[Int_data.z1(i);Int_data.z2(i)];
        % use DXFLib, Version 0.9.1 (Copyright (c) 2009-2011 Grzegorz Kwiatek)
        % to export DXF for each plane
        
        if ispc
            
            mkdir (pathname,'KinAnalysis\Inters\ObliqueToppling\')
            dxf_name=([num2str(CIMperc(i,3),'%03.f'),'_Inters_',num2str(i),'_Discs',num2str(Int_data.Nplane_i(i)),'_',num2str(Int_data.Nplane_j(i)),'.dxf']);
            FID=dxf_open(fullfile(pathname,'KinAnalysis\Inters\ObliqueToppling\'),dxf_name);
            FID = dxf_set(FID,'Color',Color{5});
            dxf_polyline(FID, X, Y, Z);
            dxf_close(FID); % end of DXF exportation
            
        elseif ismac
            
            mkdir (pathname,'KinAnalysis/Inters/ObliqueToppling/')
            dxf_name=([num2str(CIMperc(i,3),'%03.f'),'_Inters_',num2str(i),'_Discs',num2str(Int_data.Nplane_i(i)),'_',num2str(Int_data.Nplane_j(i)),'.dxf']);
            FID=dxf_open(fullfile(pathname,'KinAnalysis/Inters/ObliqueToppling/'),dxf_name);
            FID = dxf_set(FID,'Color',Color{5});
            dxf_polyline(FID, X, Y, Z);
            dxf_close(FID); % end of DXF exportation
            
        end
    end
end
disp('End of critical disconinuities export process')
disp('Starting to export the PointCloud')
disp(['->Exporting_results_at_= ',num2str(toc/60,'%.1f'), '(min)'])



%Previous to the result point cloud exportation, it removes the NaN value
PC(isnan(overhanging),:)=[];
slopeDipdirDip(isnan(overhanging),:)=[];
PCDicsInters(isnan(overhanging),:)=[];
PCIntInters(isnan(overhanging),:)=[];
%nCritic(isnan(overhanging),:)=[];
overhanging(isnan(overhanging),:)=[];
if ispc
    fid = fopen(fullfile(pathname,'KinAnalysis\','Roka_pointcloud.txt'),'wt');
elseif ismac
    fid = fopen(fullfile(pathname,'KinAnalysis/','ROKA_poincloud.txt'),'wt');
else
    fid = fopen(fullfile(pathname,'Roka_pointcloud.txt'),'wt');
end
fprintf(fid, 'X Y Z DipDir Dip Overhanging AllDisc_inters AllInt_inters PlanarSliding FlexuralToppling WedgeSliding DirectToppling ObliqueToppling \n');
%fprintf(fid, '%.6f %.6f %.6f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f \n', [PC(i,1:3),slopeDipdirDip(i,:),overhanging(i,:),PCDicsInters(i,:),PCIntInters(i,:),nCritic(i,:)]);
fprintf(fid, '%.6f %.6f %.6f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f \n', [PC(:,1:3),slopeDipdirDip(:,:),overhanging(:,:),PCDicsInters(:,:),PCIntInters(:,:),PCcritic(:,:)]');

fclose(fid);
if ispc
    save(fullfile(pathname,['KinAnalysis\ROKA_workspace.mat']))
elseif ismac
    save(fullfile(pathname,['KinAnalysis/ROKA_workspace.mat']))
else
    save(fullfile(pathname,['ROKA_workspace.mat']))
end
disp('End of the PC export process')
warning('on')
end