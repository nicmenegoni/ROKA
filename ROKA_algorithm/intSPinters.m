function int_in_out = intSPinters(p1line, p2line,Lline, SP_center, SP_N, SPr)
%INTERSECTLINEPLANE Intersection point between a 3D line and a plane
%%code strongly modified by Niccolò Menegoni 10/06/2020
%   INPUT:
%   p1line = point1 of intersection
%   p2line = point2 of intersection
%   Lline = length of intersection
%   SP_center = Scanplane center
%   SP_N = Scanplane normal
%   SPr = Scanplane radius


%   PT = intersectLinePlane(LINE, PLANE)
%   Returns the intersection point of the given line and the given plane.
%   LINE:  [x0 y0 z0 dx dy dz]
%   PLANE: [x0 y0 z0 dx1 dy1 dz1 dx2 dy2 dz2]
%   PT:    [xi yi zi]
%   If LINE and PLANE are parallel, return [NaN NaN NaN].
%   If LINE (or PLANE) is a matrix with 6 (or 9) columns and N rows, result
%   is an array of points with N rows and 3 columns.
%
%   PT = intersectLinePlane(LINE, PLANE, TOL)
%   Specifies the tolerance factor to test if a line is parallel to a
%   plane. Default is 1e-14.
%
%   Example
%     % define horizontal plane through origin
%     plane = [0 0 0   1 0 0   0 1 0];
%     % intersection with a vertical line
%     line = [2 3 4  0 0 1];
%     intersectLinePlane(line, plane)
%     ans =
%        2   3   0
%     % intersection with a line "parallel" to plane
%     line = [2 3 4  1 2 0];
%     intersectLinePlane(line, plane)
%     ans =
%       NaN  NaN  NaN
%
%   See also:
%   lines3d, planes3d, points3d, clipLine3d
%
%   ---------
%   author : David Legland
%   INRA - TPV URPOI - BIA IMASTE
%   created the 17/02/2005.
%

%   HISTORY
%   24/11/2005 add support for multiple input
%   23/06/2006 correction from Songbai Ji allowing different number of
%       lines or plane if other input has one row
%   14/12/2006 correction for parallel lines and plane normals
%   05/01/2007 fixup for parallel lines and plane normals
%   24/04/2007 rename as 'intersectLinePlane'
%   11/19/2010 Added bsxfun functionality for improved speed (Sven Holcombe)
%   01/02/2011 code cleanup, add option for tolerance, update doc

% extract tolerance if needed
tol = 1e-5;
int_in_out=zeros(length(p1line(:,1)),1);%create an array to store the instersection test results

for i = 1 : length(p1line(:,1))
    line=[p1line(i,:),p2line(i,:)-p1line(i,:)];
    % difference between origins of plane and line
    dp = SP_center - line(1:3);
    
    % dot product of line direction with plane normal
    denom = sum(SP_N .* line(4:6));
    
    % relative position of intersection point on line (can be inf in case of a
    % line parallel to the plane)
    t = sum(SP_N .* dp) ./ denom;
    
    % compute coord of intersection point
    point = (line(1:3) + ([t t t] .* line(4:6)));
    % set indices of line and plane which are parallel to NaN
    par = abs(denom) < tol;
    point(par,:) = NaN;
    
    if ~isnan(point(1))
        if  (Lline(i) > sqrt( (p1line(i,1)-point(1))^2 + (p1line(i,2)-point(2))^2 + (p1line(i,3)-point(3))^2)) &&...
                (Lline(i) > sqrt( (p2line(i,1)-point(1))^2 + (p2line(i,2)-point(2))^2 + (p2line(i,3)-point(3))^2)) &&...
                    (SPr * 1.01 > sqrt( (SP_center(1)-point(1))^2 + (SP_center(2)-point(2))^2 + (SP_center(3)-point(3))^2))% SPr*1.01 is used to add a 1% of uncertanity onto the intersection calculation;
            
            int_in_out(i) = 1;
            
        end
        
    end
    
end

end
