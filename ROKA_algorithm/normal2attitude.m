function [dipdir,dip]= normal2attitude(N)


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
elseif cosAlpha > 0 && cosBeta == 0
    dipdir = 90;
elseif cosAlpha == 0 && cosBeta > 0
    dipdir = 0;
elseif cosAlpha == 0 && cosBeta== 0
    dipdir= 0;
else
    dipdir=0;
end

if dip > 90
    dip = dip - 180;
%     if dipdir<180
%         dipdir=dipdir+180;
%     else
%         dipdir=dipdir-180;
%     end
end
end



% REFERENCES:
%
% - Jones, R. R., Pearce, M. A., Jacquemyn, C., & Watson, F. E. (2016).
%   Robust best-fit planes from geospatial data. Geosphere, 12(1), 196-202.
%
% - Seers, T. D., & Hodgetts, D. (2016). Probabilistic constraints on
%   structural lineament best fit plane precision obtained through
%   numerical analysis. Journal of Structural Geology, 82, 37-47.