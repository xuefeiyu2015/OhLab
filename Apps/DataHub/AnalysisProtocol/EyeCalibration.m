function CalibratedEyes=EyeCalibration(EyeDataRaw,EyeData,EyeCalInfo,Ref_Loc,CalPlotCheck);
%Program to calibrateEyes automatically;(from central fixations )
%by Xuefei Yu 08202020
EyeDataRaw_X=EyeDataRaw{1};
EyeDataRaw_Y=EyeDataRaw{2};


EyeData_X=EyeData{1};
EyeData_Y=EyeData{2};
X_Gain=EyeCalInfo(1,1);
Y_Gain=EyeCalInfo(1,3);

X_Bias=EyeCalInfo(1,2);
Y_Bias=EyeCalInfo(1,4);
UniqueLoc=unique(Ref_Loc,'row');





%First calibrate center
CenterLoc=ismembertol(Ref_Loc,[0,0],10^-4,'ByRows',true);
EyeData_X_C=EyeData_X(CenterLoc);
EyeData_Y_C=EyeData_Y(CenterLoc);


Offset_X=mean(cellfun(@(x) ReturnMostFreq(x),EyeData_X_C));
Offset_Y=mean(cellfun(@(x) ReturnMostFreq(x),EyeData_Y_C));


Cal_X=cellfun(@(x) x-Offset_X,EyeData_X,'Uniform',0);
Cal_Y=cellfun(@(x) x-Offset_Y,EyeData_Y,'Uniform',0);


%Calibrate the raw eye data
EyeDataRaw_X_Cal=cellfun(@(x) x-Offset_X,EyeDataRaw_X,'Uniform',0);
EyeDataRaw_Y_Cal=cellfun(@(x) x-Offset_Y,EyeDataRaw_Y,'Uniform',0);


%Next is for the gain
%Then calibrate Horizontal axis
HorizontalLoc_Left=ismembertol(Ref_Loc,UniqueLoc(abs(UniqueLoc(:,2))<10^-4 & UniqueLoc(:,1)<-10^-4,:),10^-4,'ByRows',true);
HorizontalLoc_Right=ismembertol(Ref_Loc,UniqueLoc(abs(UniqueLoc(:,2))<10^-4 & UniqueLoc(:,1)>10^-4,:),10^-4,'ByRows',true);
%HorizontalLoc=ismembertol(Ref_Loc,[10,0;-10,0],10^-4,'ByRows',true);

Num_Hori_Left=sum(HorizontalLoc_Left);
Num_Hori_Right=sum(HorizontalLoc_Right);
%First Left
if Num_Hori_Left<=5
    X_Gain_Cal_Left=X_Gain;
else
%Only start horizontal calibrate whenever there is enough data
EyeData_X_HL=Cal_X(HorizontalLoc_Left);
EyeData_Y_HL=Cal_Y(HorizontalLoc_Left);

X_tmp=mean(cellfun(@(x) ReturnMostFreq(x),EyeData_X_HL)/abs(X_Gain));
Y_tmp=mean(cellfun(@(x) ReturnMostFreq(x),EyeData_Y_HL)/abs(Y_Gain));

X_Gain_Cal_Left=sign(X_Gain)*unique(abs(UniqueLoc((abs(UniqueLoc(:,2))<10^-4 & abs(UniqueLoc(:,1))>10^-4),1)))/abs(X_tmp);

%X_Gain_Cal=sign(X_Gain)*unique(abs(UniqueLoc((abs(UniqueLoc(:,2))<10^-4 & abs(UniqueLoc(:,1))>10^-4),1)))/sqrt(X_tmp.^2+Y_tmp^2);
end

%Then Right
if Num_Hori_Right<=5
    X_Gain_Cal_Right=X_Gain_Cal_Left;
else
%Only start horizontal calibrate whenever there is enough data
EyeData_X_HR=Cal_X(HorizontalLoc_Right);
EyeData_Y_HR=Cal_Y(HorizontalLoc_Right);

X_tmp=mean(cellfun(@(x) ReturnMostFreq(x),EyeData_X_HR)/abs(X_Gain));
Y_tmp=mean(cellfun(@(x) ReturnMostFreq(x),EyeData_Y_HR)/abs(Y_Gain));

X_Gain_Cal_Right=sign(X_Gain)*unique(abs(UniqueLoc((abs(UniqueLoc(:,2))<10^-4 & abs(UniqueLoc(:,1))>10^-4),1)))/abs(X_tmp);

%X_Gain_Cal=sign(X_Gain)*unique(abs(UniqueLoc((abs(UniqueLoc(:,2))<10^-4 & abs(UniqueLoc(:,1))>10^-4),1)))/sqrt(X_tmp.^2+Y_tmp^2);
end



%Next calibrate Vertical axis

VerticalLocUp=ismembertol(Ref_Loc,UniqueLoc((abs(UniqueLoc(:,1))<10^-4 & UniqueLoc(:,2)>10^-4),:),10^-4,'ByRows',true);
VerticalLocDown=ismembertol(Ref_Loc,UniqueLoc((abs(UniqueLoc(:,1))<10^-4 & UniqueLoc(:,2)<-10^-4),:),10^-4,'ByRows',true);



%HorizontalLoc=ismembertol(Ref_Loc,[10,0;-10,0],10^-4,'ByRows',true);
Num_Ver_Up=sum(VerticalLocUp);
Num_Ver_Down=sum(VerticalLocDown);
Num_Ver=Num_Ver_Up+Num_Ver_Down;
if Num_Ver_Up<=5
    Y_Gain_Cal_Up=Y_Gain;
else
%Only start horizontal calibrate whenever there is enough data
EyeData_X_VU=Cal_X(VerticalLocUp);
EyeData_Y_VU=Cal_Y(VerticalLocUp);

X_tmp=mean(cellfun(@(x) ReturnMostFreq(x),EyeData_X_VU)/abs(X_Gain));
Y_tmp=mean(cellfun(@(x) ReturnMostFreq(x),EyeData_Y_VU)/abs(Y_Gain));

Y_Gain_Cal_Up=sign(Y_Gain)*unique(abs(UniqueLoc((abs(UniqueLoc(:,1))<10^-4 & abs(UniqueLoc(:,2))>10^-4),2)))/abs(Y_tmp);
%Y_Gain_Cal=sign(Y_Gain)*unique(abs(UniqueLoc((abs(UniqueLoc(:,1))<10^-4 & abs(UniqueLoc(:,2))>10^-4),2)))/sqrt(X_tmp.^2+Y_tmp^2);
    
end

if Num_Ver_Down<=5
    Y_Gain_Cal_Down=Y_Gain_Cal_Up;
else
%Only start horizontal calibrate whenever there is enough data
EyeData_X_VD=Cal_X(VerticalLocDown);
EyeData_Y_VD=Cal_Y(VerticalLocDown);

X_tmp=mean(cellfun(@(x) ReturnMostFreq(x),EyeData_X_VD)/abs(X_Gain));
Y_tmp=mean(cellfun(@(x) ReturnMostFreq(x),EyeData_Y_VD)/abs(Y_Gain));

Y_Gain_Cal_Down=sign(Y_Gain)*unique(abs(UniqueLoc((abs(UniqueLoc(:,1))<10^-4 & abs(UniqueLoc(:,2))>10^-4),2)))/abs(Y_tmp);
%Y_Gain_Cal=sign(Y_Gain)*unique(abs(UniqueLoc((abs(UniqueLoc(:,1))<10^-4 & abs(UniqueLoc(:,2))>10^-4),2)))/sqrt(X_tmp.^2+Y_tmp^2);
    
end



Cal_X_AfterGain=CalibrateWithGain(Cal_X,X_Gain,[X_Gain_Cal_Left,X_Gain_Cal_Right]);
Cal_Y_AfterGain=CalibrateWithGain(Cal_Y,Y_Gain,[Y_Gain_Cal_Down,Y_Gain_Cal_Up]);


%Calibrate the raw eye data
EyeDataRaw_X_Cal_AfterGain=CalibrateWithGain(EyeDataRaw_X_Cal,X_Gain,[X_Gain_Cal_Left,X_Gain_Cal_Right]);
EyeDataRaw_Y_Cal_AfterGain=CalibrateWithGain(EyeDataRaw_Y_Cal,Y_Gain,[Y_Gain_Cal_Down,Y_Gain_Cal_Up]);


%Finally the coupling
Couple_XY=0;
Couple_YX=0;

%
%First horizontal
Num_Hori=Num_Hori_Right+Num_Hori_Left;
if Num_Hori<=5
    Couple_XY=0;
else
    %Only start horizontal calibrate whenever there is enough data
    Horizontal_Loc_tmp=UniqueLoc((abs(UniqueLoc(:,2))<10^-4 & abs(UniqueLoc(:,1))>10^-4),:);
    %Choose the larger number of data site
    for i=1:size(Horizontal_Loc_tmp,1)
        num_tmp(i)=sum(ismembertol(Ref_Loc,Horizontal_Loc_tmp(i,:),10^-4,'ByRows',true));
        
    end
    
    
    Horizontal_Loc_tmp=Horizontal_Loc_tmp(num_tmp==max(num_tmp),:);%One side is enough
    Horizontal_Loc_tmp=Horizontal_Loc_tmp(1,:);
   %Horizontal_Loc_tmp=Horizontal_Loc_tmp(2,:);
    
    Selection_H=ismembertol(Ref_Loc,Horizontal_Loc_tmp,10^-4,'ByRows',true);

EyeData_X_H=Cal_X_AfterGain(Selection_H);
EyeData_Y_H=Cal_Y_AfterGain(Selection_H);

X_tmp=mean(cellfun(@(x) ReturnMostFreq(x),EyeData_X_H));
Y_tmp=mean(cellfun(@(x) ReturnMostFreq(x),EyeData_Y_H));


%Now the coupling

CosTheta = max(min(dot([X_tmp,Y_tmp],Horizontal_Loc_tmp)/(norm([X_tmp,Y_tmp])*norm(Horizontal_Loc_tmp)),1),-1);
ThetaInDegrees = real(acosd(CosTheta));
%ThetaInDegrees =-10;
Couple_XY=-sign(Horizontal_Loc_tmp(1))*sign(-Y_tmp)*ThetaInDegrees;

   
end
%Cal_X_AfterXCoupling=cellfun(@(x,y) x*cos(Couple_XY*pi/180)+y*sin(Couple_XY*pi/180),Cal_X_AfterGain,Cal_Y_AfterGain,'Uniform',0);
Cal_Y_AfterXCoupling=cellfun(@(x,y) x*sin(-Couple_XY*pi/180)+y*cos(Couple_XY*pi/180),Cal_X_AfterGain,Cal_Y_AfterGain,'Uniform',0);
Cal_X_AfterXCoupling=Cal_X_AfterGain;

%Calibrate the raw eye data
EyeDataRaw_Y_Cal_AfterXCoupling=cellfun(@(x,y) x*sin(-Couple_XY*pi/180)+y*cos(Couple_XY*pi/180),EyeDataRaw_X_Cal_AfterGain,EyeDataRaw_Y_Cal_AfterGain,'Uniform',0);
EyeDataRaw_X_Cal_AfterXCoupling=EyeDataRaw_X_Cal_AfterGain;

%}


%
%Y coupling
%Next calibrate Vertical axis

if Num_Ver<=5
    Couple_YX=0;
else
%Only start vertical calibrate whenever there is enough data
Vertical_Loc_tmp=UniqueLoc((abs(UniqueLoc(:,1))<10^-4 & abs(UniqueLoc(:,2))>10^-4),:);
%Choose the larger number of data site
for i=1:size(Vertical_Loc_tmp,1)
    
        num_tmp(i)=sum(ismembertol(Ref_Loc,Vertical_Loc_tmp(i,:),10^-4,'ByRows',true));
    
        
end
    
    
Vertical_Loc_tmp=Vertical_Loc_tmp(num_tmp==max(num_tmp),:);%One side is enough
Vertical_Loc_tmp=Vertical_Loc_tmp(1,:);

Selection_V=ismembertol(Ref_Loc,Vertical_Loc_tmp,10^-4,'ByRows',true);

EyeData_X_V=Cal_X_AfterXCoupling(Selection_V);
EyeData_Y_V=Cal_Y_AfterXCoupling(Selection_V);

X_tmp=mean(cellfun(@(x) ReturnMostFreq(x),EyeData_X_V));
Y_tmp=mean(cellfun(@(x) ReturnMostFreq(x),EyeData_Y_V));

%Now the coupling

CosTheta = max(min(dot([X_tmp,Y_tmp],Vertical_Loc_tmp)/(norm([X_tmp,Y_tmp])*norm(Vertical_Loc_tmp)),1),-1);
ThetaInDegrees = real(acosd(CosTheta));

Couple_YX=-sign(Vertical_Loc_tmp(2))*sign(X_tmp)*ThetaInDegrees;

    
end
Cal_X_AfterYCoupling=cellfun(@(x,y) x*cos(Couple_YX*pi/180)+y*sin(Couple_YX*pi/180),Cal_X_AfterXCoupling,Cal_Y_AfterXCoupling,'Uniform',0);
%Cal_Y_AfterYCoupling=cellfun(@(x,y) x*sin(-Couple_YX*pi/180)+y*cos(Couple_YX*pi/180),Cal_X_AfterXCoupling,Cal_Y_AfterXCoupling,'Uniform',0);
Cal_Y_AfterYCoupling=Cal_Y_AfterXCoupling;

%Calibrate the raw eye data
EyeDataRaw_X_Cal_AfterYCoupling=cellfun(@(x,y) x*cos(Couple_YX*pi/180)+y*sin(Couple_YX*pi/180),EyeDataRaw_X_Cal_AfterXCoupling,EyeDataRaw_Y_Cal_AfterXCoupling,'Uniform',0);
EyeDataRaw_Y_Cal_AfterYCoupling=EyeDataRaw_Y_Cal_AfterXCoupling;
%}


Cal_X_Final=Cal_X_AfterYCoupling;
Cal_Y_Final=Cal_Y_AfterYCoupling;
%
%Visualize the calibration 
if CalPlotCheck
figure

for i=1:size(UniqueLoc,1)
    subplot(1,2,1)
    hold on
    cellfun(@(x,y) plot(x,y,'.r'),EyeData_X(ismembertol(Ref_Loc,UniqueLoc(i,:),10^-4,'ByRows',true)),EyeData_Y(ismembertol(Ref_Loc,UniqueLoc(i,:),10^-4,'ByRows',true)),'Uniform',0)
    plot(UniqueLoc(i,1),UniqueLoc(i,2),'+b','LineWidth',5);
    title('Before Cali');
    xlim([-20,20]);
    ylim([-20,20]);
    subplot(1,2,2)
    hold on
   
     cellfun(@(x,y) plot(x,y,'.r'),Cal_X_Final(ismembertol(Ref_Loc,UniqueLoc(i,:),10^-4,'ByRows',true)),Cal_Y_Final(ismembertol(Ref_Loc,UniqueLoc(i,:),10^-4,'ByRows',true)),'Uniform',0)
    plot(UniqueLoc(i,1),UniqueLoc(i,2),'+b','LineWidth',5);
    title('After Cali');
    xlim([-20,20]);
    ylim([-20,20]);
end
end

%}
%Return the raw trace and the calibrated eye information
CalibratedEyes.EyeTrace={EyeDataRaw_X_Cal_AfterYCoupling,EyeDataRaw_Y_Cal_AfterYCoupling};
CalibratedEyes.EyeInfo.XGain=[X_Gain_Cal_Left,X_Gain_Cal_Right];
CalibratedEyes.EyeInfo.YGain=[Y_Gain_Cal_Down,Y_Gain_Cal_Up];
CalibratedEyes.EyeInfo.Offset=[Offset_X,Offset_Y];
CalibratedEyes.EyeInfo.Couple=[Couple_XY,Couple_YX];



end
function val=ReturnMostFreq(array);
[val_tmp,idx_tmp]=histcounts(array);
idx=idx_tmp(1:end-1)+diff(idx_tmp)/2;
val=idx(find(val_tmp==max(val_tmp),1,'first'));

end
function cal=CalibrateWithGain(eye,Gain_ori,Gain_Cal);
neg=Gain_Cal(1);%Left on horizontal/Down on verticle
pos=Gain_Cal(2);%Right on horizontal/Up on verticle
ori=Gain_ori;
cal=cellfun(@(x) SeparateAssign(x,neg,pos,ori),eye,'Uniform',0);



end
function x=SeparateAssign(x_ori,neg,pos,ori);
x=x_ori;
x(x_ori<0)=x_ori(x_ori<0)/ori*neg;
x(x_ori>0)=x_ori(x_ori>0)/ori*pos;

end
%{
function F=EyeBasic(x,y,x_loc,x_offset,x_gain,couple_xy)
F=(x+x_offset)*x_gain+y*couple_xy-x_loc;

end
%}