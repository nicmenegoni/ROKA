
%% Read data and write Dip, Dip Direction, Set, Pole Dip Direction and Dip vectors
%Discontinuity planes
function [idxPS,idxFT,idxWS,idxDT,idxOT,nCritic_line]=KinematicAnalysis_parfor(data,Int_data,uselatlimits,latlimit,friction,in_out,int_in_out,planeslope)

%CDM matrix contains the critical value of the disconinuity in this order:
%                   1)PS; 2) FT;
%CIM matrix contains the critical value of the intersections in this order:
%                   1)WS; 2) DT; 3) OT


%% Read data and write Dip, Dip Direction, Set, Pole Dip Direction and Dip vectors
%Discontinuity planes
%calcualte orientation of slope
num_Int=numel(Int_data.Trend);
num_disc = numel(data.Dip);
DipDir=data.DipDirection;
Dip=data.Dip;
Plunge=Int_data.Plunge;
Trend=Int_data.Trend;
% and calculate the plunge and trend of the disconuity poles
pole_Dip = 90 - Dip;% calculate the plunge of the pole
pole_DipDir(DipDir < 181) = DipDir(DipDir < 181) + 180;% calculate the pole trend
pole_DipDir(DipDir > 180) = DipDir(DipDir > 180) - 180;% calculate the pole trend
%% Setting storage variables
nPS=0; nFT=0; nWS=0; nDT=0; nOT=0;%indicators for the possible MoF onto the point cloud
nCritic_line=zeros(1,5);%here, the previous MoF indicators will be stored
criticDisc=zeros(num_Int,2);%array where the the criticallity for the intersection are stored (1=positive to MoFs; 0=negative)
criticInt=zeros(num_disc,3);%array where the the criticallity for the intersection are stored (1=positive to MoFs; 0=negative)

%% Kynematic analysis

%Setting properly planeslope
if planeslope(1,1)<90%calculate strike of the slope
    strikeslope = 360 + planeslope(1,1) - 90;
else
    strikeslope = planeslope(1,1) - 90;
end

if planeslope(1,1) < latlimit
    sxlatlimit = 360 + planeslope(1,1) - latlimit;%sxlatlimit=laterl limit at left of dip direction of slope
    dxlatlimit = planeslope(1,1) + latlimit;%dxlatlimit=laterl limit at rigth of dip direction of slope
elseif planeslope(1,1) >= 360 - latlimit
    sxlatlimit = planeslope(1,1) - latlimit;
    dxlatlimit = planeslope(1,1) + latlimit - 360;
else
    sxlatlimit = planeslope(1,1) - latlimit;
    dxlatlimit = planeslope(1,1) + latlimit;
end

%% Planare sliding (PS) and Flexural Toppling kinematic analysis
%A plane could act as a sliding plane if is dip-slope and dipping with an
%angle lower than apparent angle of the slope along dipdirection of the plane
%and higher than frction angle.
app_angle=zeros(num_disc,1);
slopelimit=[planeslope(1,1),  planeslope(1,2)-friction];
pole_app_angle=zeros(num_disc,1);
for i = 1 : num_disc % calculate if discontinuity plane could be a sliding plane
    
    if in_out(i)>0
        %% Planar Sliding Test
        if (DipDir(i) > strikeslope && DipDir(i) < strikeslope + 180)%Calculate apparent angle of the slope along Dip Direction of the discontinuity
            app_angle(i,1) = atand(tand(planeslope(1,2)) * cosd(DipDir(i)- planeslope(1,1)));
            %tic;app_angle = atand(tand(planeslope(1,2)) * cosd(test- planeslope(1,1)));toc
        elseif (strikeslope+180> 360 && DipDir(i) < strikeslope - 180)
            app_angle(i,1) = atand(tand(planeslope(1,2)) * cosd(360 + DipDir(i)- planeslope(1,1)));
        else
            app_angle(i,1)=NaN;
        end
        if Dip(i) > friction && Dip(i) < app_angle(i)
            if uselatlimits ==1 && (DipDir(i)>sxlatlimit && DipDir(i)<dxlatlimit) || (sxlatlimit >dxlatlimit && DipDir(i)<sxlatlimit && DipDir(i)<dxlatlimit)|| (sxlatlimit >dxlatlimit && DipDir(i)>sxlatlimit && DipDir(i)>dxlatlimit)
                %1st case if sx is < dx latlimit, 2nd and 3rd when dx < sx
                %lat limit with DipDir of discontinuity >340 && <360 (3rdcase) and
                %DipDir >0 && <20 (2nd case).
                criticDisc(i,1) = 1;%value assigned to critical discontinuity inside the lateral limits
                nPS=nPS+1;%save the number of critical disconinuities for the specific localition/point
            elseif uselatlimits == 0
                criticDisc(i,1) = 1;%value assigned to crtitical discontinuity outside the lateral limits
                nPS=nPS+1;
            end
        end
        
        %% Flexural Toppling test
        %Calculate apparent angle of the slope along Dip Direction of the discontinuity
        if (pole_DipDir(i)>sxlatlimit && pole_DipDir(i)<dxlatlimit) || (sxlatlimit >dxlatlimit && pole_DipDir(i)<sxlatlimit && pole_DipDir(i)<dxlatlimit)|| (sxlatlimit >dxlatlimit && pole_DipDir(i)>sxlatlimit && pole_DipDir(i)>dxlatlimit)
            %if cylcle for Pole DipDirections included in the lateral limit
            
            if (pole_DipDir(i) > strikeslope && pole_DipDir(i) < strikeslope + 180)
                pole_app_angle(i,1) = atand(tand(slopelimit(1,2)) * cosd(pole_DipDir(i)- slopelimit(1,1)));
            elseif (strikeslope+180> 360 && pole_DipDir(i) < strikeslope - 180)
                pole_app_angle(i,1) = atand(tand(slopelimit(1,2)) * cosd(360 + pole_DipDir(i)- slopelimit(1,1)));
            else
                pole_app_angle(i,1)=NaN;
            end
        end
        
        %Calculate if discontinuity could be affected by flexural toppling
        if pole_Dip(i) < pole_app_angle(i)
            criticDisc(i,2) = 1;
            nFT=nFT+1;
        end
    end
end
%% WedgeSliding (WS), DirectToppling (DT) and ObliqueToppling (OT) kinematic analysis
%Preconfiguration for wedge sliding
Int_app_angle=zeros(num_Int,1);
Int_app_angle_friction=zeros(num_Int,1);

%Preconfiguration for direct and oblique toppling
upperTopplingLimit = 90 - planeslope(1,2);%angle between slope face and intersection
%(ho qualche dubbio sul considerare questo limite costante e non variabile
%in base all'inclinazione apparente della parete lungo la direzione
%d'immersione dell'intersezione.
if sxlatlimit(1,1) < 180%calculate the opposite of the lateral limit (sx)
    oppsxlatlimit = sxlatlimit+180;
else
    oppsxlatlimit = sxlatlimit-180;
end
if dxlatlimit(1,1) < 180%calculate the opposite of the lateral limits (dx)
    oppdxlatlimit = dxlatlimit+180;
else
    oppdxlatlimit = dxlatlimit-180;
end


for i = 1 : num_Int % calculate if discontinuity plane could be a sliding plane
    
    if int_in_out(i)>0
        %% Wedge Sliding test
        if (Trend(i) > strikeslope && Trend(i) < strikeslope + 180)%Calculate apparent angle of the slope along Dip Direction of the discontinuity
            Int_app_angle(i,1) = atand(tand(planeslope(1,2)) * cosd(Trend(i)- planeslope(1,1)));
            Int_app_angle_friction(i,1) = atand(tand(friction) * cosd(Trend(i)- planeslope(1,1)));
            
        elseif (strikeslope+180> 360 && Trend(i) < strikeslope - 180)
            
            Int_app_angle(i,1) = atand(tand(planeslope(1,2)) * cosd(360 + Trend(i)- planeslope(1,1)));
            Int_app_angle_friction(i,1) = atand(tand(friction) * cosd(Trend(i)- planeslope(1,1)));
        else
            Int_app_angle(i,1)=NaN;
            Int_app_angle_friction(i,1)=NaN;
        end
        
        if Plunge(i) < Int_app_angle(i)
            if uselatlimits==1 && Plunge(i)>friction
                criticInt(i,1) = 1;%value assigned to critical discontinuity inside the lateral limits
                nWS=nWS+1;
            elseif uselatlimits==0 && Plunge(i)>friction
                %elseif uselatlimits==0 && (Plunge(i)>friction|| (Plunge(i)<friction && Plunge(i)>Int_app_angle_friction(i)))
                criticInt(i,1) = 1;%value assigned to critical discontinuity inside the lateral limits
                nWS=nWS+1;
            end
        end
        %% DirectToppling (DT) and ObliqueToppling (OT) test
        if (Trend(i) > oppsxlatlimit && Trend(i) < oppdxlatlimit) || ...
                (oppdxlatlimit < oppsxlatlimit && Trend(i) > oppsxlatlimit && Trend(i) > oppdxlatlimit) || ...
                (oppdxlatlimit < oppsxlatlimit && Trend(i) < oppsxlatlimit && Trend(i) < oppdxlatlimit)
            
            if Plunge(i) < (90 - planeslope(1,2))
                criticInt(i,2)=1;
                nDT=nDT+1;
            end
            
        elseif (Plunge(i) < (90 - planeslope(1,2)) && strikeslope <= 180 && planeslope(1,1)<= 180 && (Trend(i) < strikeslope || Trend(i)> (strikeslope+180))) || ...
                (Plunge(i) < (90 - planeslope(1,2)) && strikeslope>= 180 && (planeslope(1,1)<= 360 || planeslope(1,1)<= 0) && Trend(i)>(strikeslope-180) && Trend(i)<strikeslope)
            criticInt(i,3)=1;
            nOT=nOT+1;
            
        end
    end
end
idxPS=find(criticDisc(:,1));
idxFT=find(criticDisc(:,2));
idxWS=find(criticInt(:,1));
idxDT=find(criticInt(:,2));
idxOT=find(criticInt(:,3));

% for i = 1 : num_disc %Discontinuity Planes
%     if in_out>0
%         if (pole_DipDir(i) > oppsxlatlimit && pole_DipDir(i) < oppdxlatlimit) || ...
%                 (oppdxlatlimit < oppsxlatlimit && pole_DipDir(i) > oppsxlatlimit && pole_DipDir(i) > oppdxlatlimit) || ...
%                 (oppdxlatlimit < oppsxlatlimit && pole_DipDir(i) < oppsxlatlimit && pole_DipDir(i) < oppdxlatlimit)
%             if pole_Dip(i) > upperTopplingLimit && pole_Dip(i) < (90 - friction)
%                 BasalPlane(i) = 2;
%             elseif pole_Dip(i)>(90-friction)
%                 BasalPlane(i) = 1;
%             else
%                 BasalPlane(i) = NaN;
%             end
%         elseif (pole_Dip(i) > (90 - friction) && strikeslope <= 180 && planeslope(1,1)<= 180 && (pole_DipDir(i) < strikeslope || pole_DipDir(i)> (strikeslope+180))) || ...
%                 (pole_Dip(i) > (90 - friction) && strikeslope>= 180 && (planeslope(1,1)<= 360 || planeslope(1,1)<= 0) && pole_DipDir(i)>(strikeslope-180) && pole_DipDir(i)<strikeslope)
%             BasalPlane(i) = 3;
%         else
%             BasalPlane(i)=NaN;
%         end
%
%     end
% end