function [xyz, xyzdot] = oe2xyz(a,e,incl,w,OM,M)
% convert orbital elements (a,e,incl,w,OM,M) (rad)
% to r,v cartesian elements in inertial reference system

global GM maxiter acc 

% i.
% calculate (E) from Kepler's equation
% M = E- e*sinE
[E, ~, ~] = Kepler(M, e, maxiter, acc);
% calculate true anomaly with e and E (f)
f=acos((cos(E)-e)./(1-e.*cos(E)));

% ii.
% satellite-earth distance (r)
r=(a.*(1-e.^2))./(1+e.*cos(f));

% iii.
% Mean motion (eta):
eta=sqrt(GM./(a.^3));

% iv. 
% calculate (q) and (qdot) PQW frame
q=[(r.*cos(f));(r.*sin(f));0];
qdot=eta.*a./sqrt(1-e.^2).*[(-sin(f));(e+cos(f));0];

% v.
% transformation by the rotation matrices of the axis
rotate=Rot(3,-OM)'*Rot(1,-incl)'*Rot(3,-w)';
xyz=rotate*q;
xyzdot=rotate*qdot;

end

