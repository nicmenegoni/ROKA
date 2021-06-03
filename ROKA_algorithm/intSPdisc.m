function [countCWDisc_inter] = intSPdisc(CWxyz, CWNxyz, CWr, xyz, Nxyz, radius,disccutoff)
%% INFORMATION
% This function count verify if a discontinuity intersects the scan plane
% of the RoKA code.
% It is based onto another script of the DIEapp suite that perform the P21
%
%--------------INPUT parameters:------------------------------------------
%   CWxyz = center of the Circular Window;
%   CWNxyz = normal vector of the Circular Window;
%   CWr = radius of the Circular Window;
%   xyz = center coordinates of the disconutinuity;
%   Nxyz = normal of the discontinuity plane;
%   radius = radius of the disconunuity disc;
%
%-------------OUTPUT parameters:------------------------------------------
%   countCWDisc_inter = vector array where the results of the intersection
%   test is stored (0 = no intersection, 1 = intersection)
%
%% Circular window - discontinuities intersections
nplane = length(radius);

% initializig a matrix for counting the intersection of each discontinuity and set
countCWDisc_inter=zeros(nplane,1);%
% Preallocating memory
lineAB=zeros(nplane,6);

for i=1:nplane
    %%Compute intersection between two plane
    tol = 1e-14;% setting the angle cutoff for which pseudo parallel planes (parallel at the CW plane) are not consider.
    if radius(i)<disccutoff
        
    else
        if abs(cross(Nxyz(i,:), CWNxyz, 2)) > tol
            % The 10 rows below  are based onto a code achived from the website
            % tbirdal.blogspot.com and are based onto a function created
            % according the the solution prioposed by Krumm (2010).
            %References:
            % Krumm J. (2010). Interesection of Two Planes. Microsoft Reseacrh. Available
            % at: https://www.microsoft.com/en-us/research/publication/intersection-of-two-planes/
            M = [2 0 0 CWNxyz(1) Nxyz(i,1)
                0 2 0 CWNxyz(2) Nxyz(i,2)
                0 0 2 CWNxyz(3) Nxyz(i,3)
                CWNxyz(1) CWNxyz(2) CWNxyz(3) 0 0
                Nxyz(i,1) Nxyz(i,2) Nxyz(i,3) 0 0];
            
            b4 = CWxyz(1).*CWNxyz(1) + CWxyz(2).*CWNxyz(2) + CWxyz(3).*CWNxyz(3);
            b5 = xyz(i,1).*Nxyz(i,1) + xyz(i,2).*Nxyz(i,2) + xyz(i,3).*Nxyz(i,3);
            b = [2*CWxyz(1) ; 2*CWxyz(2) ; 2*CWxyz(3); b4 ; b5];
            
            x = M\b;
            Pline = x(1:3)';%position of the vector line
            Vline = cross(CWNxyz, Nxyz(i,:));%direction vector of the intersection line
            
            lineAB(i, 1:6)=[Pline,Vline];%Infinite line representing intersection between two infinte planes of discontinuity and CW
            
            point1AB(i,:)= xyz(i,:);% Center of discontinuity
            point2AB= CWxyz;% center of Circular Window
            
            %Distance between disconintuity center and previous intersection line
            %dpointL_PA(i) = abs(cross((xyz(i,:)-Pline),(xyz(i,:)-(Pline+Vline)))) / abs(Pline - (Pline+Vline));
            dpointL_PA(i) = norm(cross(Vline,(Pline-xyz(i,:))) /norm(Vline));
            
            %Distance between circular window center and previous intersection line
            dpointL_PB(i) = norm(cross(Vline,(Pline-CWxyz))/norm(Vline));
            
            if abs(dpointL_PA(i))<radius(i) && abs(dpointL_PB(i))<CWr %if distance is lower than radius of 2 discs (discont. and CW) intersection is considered
                clear d_point1A d_point1B d_point2A d_point2B d_point3A d_point3B d_point4A d_point4B
                selectAB=zeros(4,4);
                
                %Define Sphere (fitted for the 2 considered discs) for Sphere-Line Intersection
                SPHEREA = [xyz(i,1), xyz(i,2), xyz(i,3),  radius(i)];%definisco la sfera per la discontinuità
                SPHEREB = [CWxyz(1,1), CWxyz(1,2), CWxyz(1,3),  CWr];%definisco la sfera per la CW
                
                %Use Spehere-Line intersection function of geom3d package
                %(INRA) to find points of intersection
                PTSA= intersectLineSphere(lineAB(i,:), SPHEREA);% for line and discontinuity sphere
                P_interA(i,1:6)=[PTSA(1,1:3), PTSA(2,1:3)];
                PTSB = intersectLineSphere(lineAB(i,:), SPHEREB);% for line and CW sphere
                P_interB(i,1:6)=[PTSB(1,1:3), PTSB(2,1:3)];
                
                
                
                % Calculation of l, the cutted trace length, trace visble
                % inside the CW
                
                %calculation distance between sphere-line intersection points
                %and centers of discontinuity and CW
                
                %Point1 = PTSA(1,:)
                d_point1A = sqrt( (xyz(i,1)-PTSA(1,1))^2 + (xyz(i,2)-PTSA(1,2))^2 + (xyz(i,3)-PTSA(1,3))^2);
                d_point1B = sqrt( (CWxyz(1)-PTSA(1,1))^2 + (CWxyz(2)-PTSA(1,2))^2 + (CWxyz(3)-PTSA(1,3))^2);
                radiusB=CWr;
                
                if ((-d_point1A+radius(i))>-0.0005) && ((-d_point1B+CWr)>-0.0005)
                    PTSAB(1,1)=PTSA(1,1);
                    PTSAB(1,2)=PTSA(1,2);
                    PTSAB(1,3)=PTSA(1,3);
                else
                    PTSAB(1,1)=NaN;
                    PTSAB(1,2)=NaN;
                    PTSAB(1,3)=NaN;
                end
                %Same procedure is done for all other 3 spheres-line intersection points calculated
                %Point 2
                d_point2A = sqrt( (xyz(i,1)-PTSA(2,1))^2 + (xyz(i,2)-PTSA(2,2))^2 + (xyz(i,3)-PTSA(2,3))^2);
                d_point2B = sqrt( (CWxyz(1)-PTSA(2,1))^2 + (CWxyz(2)-PTSA(2,2))^2 + (CWxyz(3)-PTSA(2,3))^2);
                
                if ((-d_point2A+radius(i))>-0.0005) && ((-d_point2B+radiusB)>-0.0005)
                    PTSAB(2,1)=PTSA(2,1);
                    PTSAB(2,2)=PTSA(2,2);
                    PTSAB(2,3)=PTSA(2,3);
                else
                    PTSAB(2,1)=NaN;
                    PTSAB(2,2)=NaN;
                    PTSAB(2,3)=NaN;
                end
                
                %Point 3
                d_point3A = sqrt( (xyz(i,1)-PTSB(1,1))^2 + (xyz(i,2)-PTSB(1,2))^2 + (xyz(i,3)-PTSB(1,3))^2);
                d_point3B = sqrt( (CWxyz(1)-PTSB(1,1))^2 + (CWxyz(2)-PTSB(1,2))^2 + (CWxyz(3)-PTSB(1,3))^2);
                
                if ((-d_point3A+radius(i))>-0.0005) && ((-d_point3B+radiusB)>-0.005)
                    PTSAB(3,1)=PTSB(1,1);
                    PTSAB(3,2)=PTSB(1,2);
                    PTSAB(3,3)=PTSB(1,3);
                else
                    PTSAB(3,1)=NaN;
                    PTSAB(3,2)=NaN;
                    PTSAB(3,3)=NaN;
                end
                %Point 4
                d_point4A = sqrt( (xyz(i,1)-PTSB(2,1))^2 + (xyz(i,2)-PTSB(2,2))^2 + (xyz(i,3)-PTSB(2,3))^2);
                d_point4B = sqrt( (CWxyz(1)-PTSB(2,1))^2 + (CWxyz(2)-PTSB(2,2))^2 + (CWxyz(3)-PTSB(2,3))^2);
                
                if ((-d_point4A+radius(i))>-0.0005) && ((-d_point4B+radiusB)>-0.0005)
                    PTSAB(4,1)=PTSB(2,1);
                    PTSAB(4,2)=PTSB(2,2);
                    PTSAB(4,3)=PTSB(2,3);
                else
                    PTSAB(4,1)=NaN;
                    PTSAB(4,2)=NaN;
                    PTSAB(4,3)=NaN;
                end
                
                %--Distances between the points
                d_pp13AB = dist2points(PTSAB(1,:), PTSAB(3,:));
                d_pp14AB = dist2points(PTSAB(1,:), PTSAB(4,:));
                d_pp23AB = dist2points(PTSAB(2,:), PTSAB(3,:));
                d_pp24AB = dist2points(PTSAB(2,:), PTSAB(4,:));
                
                % Erase similar point (save only 1 point for two or more
                % similar points)
                if (d_pp13AB<0.005 && d_pp24AB<0.0005)
                    PTSAB(3,1)=NaN;
                    PTSAB(3,2)=NaN;
                    PTSAB(3,3)=NaN;
                    PTSAB(4,1)=NaN;
                    PTSAB(4,2)=NaN;
                    PTSAB(4,3)=NaN;
                    
                end
                if (d_pp23AB<0.0005 && d_pp14AB<0.0005)
                    PTSAB(3,1)=NaN;
                    PTSAB(3,2)=NaN;
                    PTSAB(3,3)=NaN;
                    PTSAB(4,1)=NaN;
                    PTSAB(4,2)=NaN;
                    PTSAB(4,3)=NaN;
                    
                end
                
                
                
                % Select non-NaN values
                selectAB = ~isnan( PTSAB ) ;
                pointIntersecting=PTSAB(selectAB);
                npointIntersecting= (numel(pointIntersecting))/3;%count number of intersection points
                if npointIntersecting>1 %if point int. bigger than 1 a line of intersection is totally included in the
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
                    
                    %Calculation of length of the intersection between
                    %discontinuity disc and scan-plane
                    Int_Length= sqrt((X(2,1) - X(1,1))^2 + (Y(2,1) - Y(1,1))^2 + (Z(2,1) - Z(1,1))^2);
                    if Int_Length>disccutoff
                        countCWDisc_inter(i,1)=1;
                    end
                    
                    
                end
                
                
            end
            
        end
        
    end
end
end
