function [Tintersections] = intersectionCalculator(pathname,data,dataDips,export_AllIntersection)
%% Read data and write vectors
nplane = numel(data.Dip(:));
radius = data.Radius(:);
xyz(:,1) = data.Xcenter(:);
xyz(:,2) = data.Ycenter(:);
xyz(:,3) = data.Zcenter(:);
Nxyz(:,1) = data.Nx(:);
Nxyz(:,2)= data.Ny(:);
Nxyz(:,3) = data.Nz(:);
Set=dataDips.Set(:);%read Set value of discontinuities (random values are defined as NaN as default)
Set(isnan(Set))=0;%change Set value of random discontinuities from NaN to 0
nSet=numel(unique(Set));
Setval=(unique(Set));
if export_AllIntersection==1
    for i = 1: nSet%Create directory in which the intersection file.dxf are saved
        for j = 1:nSet
            Int_directory (1,1) = Setval(i);
            Int_directory (2,1) = Setval(j);
            Int_directory = sort(Int_directory);
            if ispc
                if Int_directory (1,1) == 0 && Int_directory (2,1)== 0
                    mkdir([pathname,'\intersection\SetRandom_SetRandom']);
                elseif Int_directory (1,1) == 0 && Int_directory (2,1) > 0
                    mkdir([pathname,'\intersection\SetRandom_Set',num2str(Int_directory (2,1))]);
                else
                    mkdir([pathname,'\intersection\Set', num2str(Int_directory (1,1)),'_Set', num2str(Int_directory (2,1))]);
                end
            elseif ismac
                if Int_directory (1,1) == 0 && Int_directory (2,1)== 0
                    mkdir([pathname,'/intersection/SetRandom_SetRandom']);
                elseif Int_directory (1,1) == 0 && Int_directory (2,1) > 0
                    mkdir([pathname,'/intersection/SetRandom_Set',num2str(Int_directory (2,1))]);
                else
                    mkdir([pathname,'/intersection/Set', num2str(Int_directory (1,1)),'_Set', num2str(Int_directory (2,1))]);
                end
            end
            
        end
    end
end
% Normals have to be oriented upward!!! (Normals have to be horiented in same
% direction if you want to calculate spacing)
Nxyz(Nxyz(:,3)<0,:)=- Nxyz(Nxyz(:,3)<0,:);

%-------------------------------------------------------------------------
%% Intersection calclulation
%% ------------------------- INTERSECTION ---------------------------------
countintersection=zeros(nplane,nplane);
countedge=zeros(nplane,nplane);
% figure(5)
hold on
Intersvalues=zeros(1,13);



%this paragrpah works in the same way of the MTL 2D calulation
%(disconituinies are considered insetad CW).
p0x = zeros(nplane-1,nplane-1);
p0y = zeros(nplane-1,nplane-1);
p0z = zeros(nplane-1,nplane-1);
dpx = zeros(nplane-1,nplane-1);
dpy = zeros(nplane-1,nplane-1);
dpz = zeros(nplane-1,nplane-1);
n1 = zeros(nplane-1,3);
n2 = zeros(nplane-1,3);
for i= 1 : nplane-1
    for j= i+1 :nplane
        n1(i,:)=Nxyz(i,:);
        n2(j,:)=Nxyz(j,:);
        tol = 1e-14;
        
        % Uses Hessian form, ie : N.p = d
        % I this case, d can be found as : -N.p0, when N is normalized
        d1 = dot(n1(i,:), xyz(i,:), 2);
        d2 = dot(n2(j,:), xyz(j,:), 2);
        
        % compute dot products
        dot1 = dot(n1(i,:), n1(i,:), 2);
        dot2 = dot(n2(j,:), n2(j,:), 2);
        dot12 = dot(n1(i,:), n2(j,:), 2);
        
        % intermediate computations
        det= dot1*dot2 - dot12*dot12;
        c1= (d1*dot2 - d2*dot12)./det;
        c2= (d2*dot1 - d1*dot12)./det;
        
        % compute line origin and direction
        p0= c1*n1(i,:) + c2*n2(j,:);
        dp = cross(n1(i,:), n2(j,:), 2);
        % test if planes are parallel
        if abs(cross(n1(i,:), n2(j,:), 2)) < tol
            p0x(i,j)= NaN;
            p0y(i,j)= NaN;
            p0z(i,j)= NaN;
            dpx(i,j)= NaN;
            dpy(i,j)= NaN;
            dpz(i,j)= NaN;
            
        else
            
            p0x(i,j)= p0(1,1);
            p0y(i,j)= p0(1,2);
            p0z(i,j)= p0(1,3);
            dpx(i,j)= dp(1,1);
            dpy(i,j)= dp(1,2);
            dpz(i,j)= dp(1,3);
            line=[p0x(i,j), p0y(i,j), p0z(i,j), dpx(i,j), dpy(i,j), dpz(i,j)];
            
            point1= [xyz(i,1), xyz(i,2), xyz(i,3)];
            point2= [xyz(j,1), xyz(j,2), xyz(j,3)];
            dpointL_P1 = norm(cross(line(:,4:6), line(:,1:3)-point1)) / norm(line(:,4:6));
            dpointL_P2 = norm(cross(line(:,4:6), line(:,1:3)-point2)) / norm(line(:,4:6));
            %define radius
            
            radius1=radius(i);
            radius2=radius(j);
            
            if dpointL_P1<radius1 && dpointL_P2<radius2
                
                %Define Sphere for Sphere-Line Intersection
                SPHERE1 = [xyz(i,1), xyz(i,2), xyz(i,3),  radius1];
                SPHERE2 = [xyz(j,1), xyz(j,2), xyz(j,3),  radius2];
                
                %Use Spehere-Line intersection function of geom3d package
                %(INRA)
                PTS1 = intersectLineSphere(line, SPHERE1);
                PTS2 = intersectLineSphere(line, SPHERE2);
                PTS=[PTS1; PTS2];
                
                %calcolo dei punti dalle sfere (es d_point11= distanza punto 1
                %da sfera 1)
                
                %Point1
                d_point11=dist2points(PTS(1,:), xyz(i,:));
                d_point12=dist2points(PTS(1,:), xyz(j,:));
                if (-d_point11+radius1)>-0.0005 && (-d_point12+radius2)>-0.0005
                    PTS(1,1)=PTS(1,1);
                    PTS(1,2)=PTS(1,2);
                    PTS(1,3)=PTS(1,3);
                else
                    PTS(1,1)=NaN;
                    PTS(1,2)=NaN;
                    PTS(1,3)=NaN;
                end
                %Point 2
                d_point21=dist2points(PTS(2,:), xyz(i,:));
                d_point22=dist2points(PTS(2,:), xyz(j,:));
                if (-d_point21+radius1)>-0.0005 && (-d_point22+radius2)>-0.0005
                    PTS(2,1)=PTS(2,1);
                    PTS(2,2)=PTS(2,2);
                    PTS(2,3)=PTS(2,3);
                else
                    PTS(2,1)=NaN;
                    PTS(2,2)=NaN;
                    PTS(2,3)=NaN;
                end
                %Point 3
                d_point31=dist2points(PTS(3,:), xyz(i,:));
                d_point32=dist2points(PTS(3,:), xyz(j,:));
                if (-d_point31+radius1)>-0.0005 && (-d_point32+radius2)>-0.005
                    PTS(3,1)=PTS(3,1);
                    PTS(3,2)=PTS(3,2);
                    PTS(3,3)=PTS(3,3);
                else
                    PTS(3,1)=NaN;
                    PTS(3,2)=NaN;
                    PTS(3,3)=NaN;
                end
                %Point 4
                d_point41=dist2points(PTS(4,:), xyz(i,:));
                d_point42=dist2points(PTS(4,:), xyz(j,:));
                if (-d_point41+radius1)>-0.0005 && (-d_point42+radius2)>-0.0005
                    PTS(4,1)=PTS(4,1);
                    PTS(4,2)=PTS(4,2);
                    PTS(4,3)=PTS(4,3);
                else
                    PTS(4,1)=NaN;
                    PTS(4,2)=NaN;
                    PTS(4,3)=NaN;
                end
                
                
                
                d_pp13 = dist2points(PTS(1,:), PTS(3,:));
                d_pp14 = dist2points(PTS(1,:), PTS(4,:));
                d_pp23 = dist2points(PTS(2,:), PTS(3,:));
                d_pp24 = dist2points(PTS(2,:), PTS(4,:));
                
                % Eliminare i punti 3 e 4 se coindidono con 1 e 2
                if (d_pp13<0.005 && d_pp24<0.0005)
                    PTS(3,1)=NaN;
                    PTS(3,2)=NaN;
                    PTS(3,3)=NaN;
                    PTS(4,1)=NaN;
                    PTS(4,2)=NaN;
                    PTS(4,3)=NaN;
                    countedge(i,j)=1;
                end
                if (d_pp23<0.0005 && d_pp14<0.0005)
                    PTS(3,1)=NaN;
                    PTS(3,2)=NaN;
                    PTS(3,3)=NaN;
                    PTS(4,1)=NaN;
                    PTS(4,2)=NaN;
                    PTS(4,3)=NaN;
                    countedge(i,j)=1;
                end
                
                % selezionare solo i valori non nulli
                select = ~isnan( PTS ) ;
                pointIntersecting=PTS(select);
                
                npointIntersecting= (numel(pointIntersecting))/3;
                
                if  npointIntersecting==1
                    Point_intersection=[pointIntersecting(1), pointIntersecting(2), pointIntersecting(3)];
                    %                     figure(1)
                    %                     drawPoint3d(Point_intersection,'marker', '+', 'markerSize', 10, 'linewidth', 3);
                    %                     hold on
                end
                if npointIntersecting>1
                    Point_intersection=[pointIntersecting(1), pointIntersecting(3), pointIntersecting(5),...
                        pointIntersecting(2), pointIntersecting(4), pointIntersecting(6)];
                    %                     figure(1)
                    %                     drawEdge3d(Point_intersection, 'color', 'r', 'linewidth', 4);
                    %                     hold on
                    countedge(i,j)=1;
                    
                    
                    if pointIntersecting(5)>pointIntersecting(6)
                        %set the point with higher Z coord. as Point1 (not Point i !!)
                        %in the final table Point i could be not equal to
                        %Point1. Point1 has to be the point with higher Z
                        %to obtain the correct plunge and trend values
                        X=[pointIntersecting(1); pointIntersecting(2)];
                        Y=[pointIntersecting(3); pointIntersecting(4)];
                        Z=[pointIntersecting(5); pointIntersecting(6)];
                    else
                        X=[pointIntersecting(2); pointIntersecting(1)];
                        Y=[pointIntersecting(4); pointIntersecting(3)];
                        Z=[pointIntersecting(6); pointIntersecting(5)];
                    end
                    
                    %Calculation of plunge and trend/azimuth
                    Int_Length= sqrt((X(2,1) - X(1,1))^2 + (Y(2,1) - Y(1,1))^2 + (Z(2,1) - Z(1,1))^2);
                    cosAlpha=(X(2,1) - X(1,1))/Int_Length;
                    cosBeta=(Y(2,1) - Y(1,1))/Int_Length;
                    cosGamma=(Z(2,1) - Z(1,1))/Int_Length;
                    plungeRad = asin(-cosGamma);
                    azimuthRad = atan(cosAlpha/cosBeta);
                    plungeDeg=rad2deg(plungeRad);
                    azimuthDeg=rad2deg(azimuthRad);
                    if cosBeta<0
                        azimuthDeg_c= 180 + azimuthDeg;
                    else
                        if cosAlpha>0
                            azimuthDeg_c= azimuthDeg;
                        elseif cosAlpha<0
                            azimuthDeg_c= 360+azimuthDeg;
                        end
                    end
                    azimuthDeg=azimuthDeg_c;
                    
                    
                    
                    %Writing a XLSX file where Point1 and Point2 XYZ, trend
                    %and plunge are saved
                    nwrite= numel(Intersvalues)/13;
                    if nwrite ==1 && sum(Intersvalues(1,:))== 0
                        
                        Intersvalues(1,1)=azimuthDeg;
                        Intersvalues(1,2)=plungeDeg;
                        Intersvalues(1,3)=Set(i);
                        Intersvalues(1,4)=Set(j);
                        Intersvalues(1,5)=(i);
                        Intersvalues(1,6)=(j);
                        Intersvalues(1,7)=Int_Length;
                        Intersvalues(1,8)=pointIntersecting(1);
                        Intersvalues(1,9)=pointIntersecting(3);
                        Intersvalues(1,10)=pointIntersecting(5);
                        Intersvalues(1,11)=pointIntersecting(2);
                        Intersvalues(1,12)=pointIntersecting(4);
                        Intersvalues(1,13)=pointIntersecting(6);
                        
                    else
                        Intersvalues(nwrite+1,1)=azimuthDeg;
                        Intersvalues(nwrite+1,2)=plungeDeg;
                        Intersvalues(nwrite+1,3)=Set(i);
                        Intersvalues(nwrite+1,4)=Set(j);
                        Intersvalues(nwrite+1,5)=(i);
                        Intersvalues(nwrite+1,6)=(j);
                        Intersvalues(nwrite+1,7)=Int_Length;
                        Intersvalues(nwrite+1,8)=pointIntersecting(1);
                        Intersvalues(nwrite+1,9)=pointIntersecting(3);
                        Intersvalues(nwrite+1,10)=pointIntersecting(5);
                        Intersvalues(nwrite+1,11)=pointIntersecting(2);
                        Intersvalues(nwrite+1,12)=pointIntersecting(4);
                        Intersvalues(nwrite+1,13)=pointIntersecting(6);
                        
                    end
                    
                    if export_AllIntersection==1
                        %Save DXF intersection file
                        dxf_name=(['Int',num2str(i),'_',num2str(j),'.dxf']);
                        Int_directory = [Set(i); Set(j)];
                        Int_directory =sort(Int_directory);
                        if ispc
                            if Int_directory(1,1)==0 && Int_directory(2,1)==0
                                FID=dxf_open(([pathname,'\intersection\SetRandom_SetRandom']),dxf_name);
                            elseif Int_directory(1,1)==0 && Int_directory(2,1)>0
                                FID=dxf_open(([pathname,'\intersection\SetRandom_Set',num2str(Int_directory(2,1))]),dxf_name);
                            elseif Int_directory(1,1)>0 && Int_directory(2,1)>0
                                FID=dxf_open(([pathname,'\intersection\Set',num2str(Int_directory(1,1)),'_Set',num2str(Int_directory(2,1))]),dxf_name);
                            end
                        elseif ismac
                            if Int_directory(1,1)==0 && Int_directory(2,1)==0
                                FID=dxf_open(([pathname,'/intersection/SetRandom_SetRandom']),dxf_name);
                            elseif Int_directory(1,1)==0 && Int_directory(2,1)>0
                                FID=dxf_open(([pathname,'/intersection/SetRandom_Set',num2str(Int_directory(2,1))]),dxf_name);
                            elseif Int_directory(1,1)>0 && Int_directory(2,1)>0
                                FID=dxf_open(([pathname,'/intersection/Set',num2str(Int_directory(1,1)),'_Set',num2str(Int_directory(2,1))]),dxf_name);
                            end
                        end
                        FID = dxf_set(FID,'Color',[1 0 0]);
                        dxf_polyline(FID, X, Y, Z);
                        dxf_close(FID);
                    end
                end
                
                
                
                
                
                %                 drawLine3d(line, 'color', 'y', 'linewidth', 2);
                %                 hold on
                
                countintersection(i,j)=1;
                
                
                
            else
                countintersection(i,j)=0;
                countedge(i,j)=0;
            end
            
        end
    end
    
    
end
Tintersections = table(Intersvalues(:,1),Intersvalues(:,2),Intersvalues(:,3),...
    Intersvalues(:,4),Intersvalues(:,5),Intersvalues(:,6),...
    Intersvalues(:,7),Intersvalues(:,8),Intersvalues(:,9),...
    Intersvalues(:,10),Intersvalues(:,11),Intersvalues(:,12),...
    Intersvalues(:,13));
Tintersections.Properties.VariableNames = {'Trend' 'Plunge' 'Set_i'...
    'Set_j' 'Nplane_i' 'Nplane_j'...
    'Length_Intersection' 'x1' 'y1'...
    'z1' 'x2' 'y2'...
    'z2' };
if export_AllIntersection==1
    if ispc
        tablefilenameTXT = (['Intersection_calculated.txt']);
        writetable(Tintersections,fullfile([pathname,'\intersection\'],tablefilenameTXT));
        tablefilenameXLSX = (['Intersection_calculated.xlsx']);
        writetable(Tintersections,fullfile([pathname,'\intersection\'],tablefilenameXLSX));
    elseif ismac
        tablefilenameTXT = (['Intersection_calculated.txt']);
        writetable(Tintersections,fullfile([pathname,'/intersection/'],tablefilenameTXT));
        tablefilenameXLSX = (['Intersection_calculated.xlsx']);
        writetable(Tintersections,fullfile([pathname,'/intersection/'],tablefilenameXLSX));
    end
end
%disp('###EnD oF iNtErSeCtIoN eLaBoRaTiOn PrOcEsS ###')
end