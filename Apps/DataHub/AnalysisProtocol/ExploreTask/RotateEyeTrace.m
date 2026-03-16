function [RotatedEye_X,RotatedEye_Y]=RotateEyeTrace(X,Y,Angle)
%Rotate eye trace relative to the target location;
%Align all target location to the up front(0 deg)
%by Xuefei Yu 08252020
RotatedEye_X=[];
RotatedEye_Y=[];
 for i=1:length(Angle)
      TransformMat=[cos(Angle(i)*pi/180),-sin(Angle(i)*pi/180);sin(Angle(i)*pi/180),cos(Angle(i)*pi/180)];
      temp=[X(i,:)',Y(i,:)']*TransformMat;
           RotatedEye_X=[RotatedEye_X,temp(:,1)',NaN];
           RotatedEye_Y=[RotatedEye_Y,temp(:,2)',NaN];
       
 end

end