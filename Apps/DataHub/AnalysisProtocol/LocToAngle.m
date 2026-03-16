%This function is to change location into angle
%i.e. cartesian into polar
%Define the upward direction as zero
function [Angle, Amp]=LocToAngle(X,Y);
Angle=NaN*ones(size(X));
Amp=NaN*ones(size(X));

for i=1:size(X,1)

 AngleCurr=atan2d(Y(i,:),X(i,:));%Transform into angle in degree
 AngleCurr=AngleCurr+90+180;
if sum(AngleCurr>=360)
                 AngleCurr(AngleCurr>=360)=AngleCurr(AngleCurr>=360)-360;
end
            
AmpCurr=sqrt(X(i,:).^2+Y(i,:).^2);
Angle(i,:)=AngleCurr;
Amp(i,:)=AmpCurr;

end


end