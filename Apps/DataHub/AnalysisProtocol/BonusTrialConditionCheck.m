%Bonus trial condition extra check;
%In some condition with uncouled fractals with between set scenes,
%there might be some coincidence that the uncouled fractals were actually
%coupled fractals associated with the between-set scenes;
%This function is to find them out and make them marked.
%{ 
Here are the need-to-correct conditions
//7.Trained Scene Between---Fractal1: Trained Objects(Between) v.s. Fractal2: Novel Objects 
8.Trained Scene Between---Trained Objects(Within) v.s. Trained(Between)
 %}
function ContainCoupleIndexUpdate=BonusTrialConditionCheck(ContainCoupleIndex,SceneIndexCheck,FractalIndexCheck,TrialConditionCheck);
ECode;
TrialIndex=unique(TrialConditionCheck(:,1));
ContainCoupleIndexUpdate=ContainCoupleIndex;
for i=1:length(TrialIndex)
    CurrTrialIndex=TrialIndex(i);
    ContainCoupleIndexCurr=ContainCoupleIndexUpdate(ContainCoupleIndexUpdate(:,1)==CurrTrialIndex,:);
    
    
    NumberOfRound=sum(FractalIndexCheck(:,1)==CurrTrialIndex);
    %SceneNameCurr=repmat(SceneIndexCheck(SceneIndexCheck(:,1)==CurrTrialIndex,3),NumberOfRound,2);
    SceneNameCurr=SceneIndexCheck(SceneIndexCheck(:,1)==CurrTrialIndex,3);
    FractalNameCurr=FractalIndexCheck(FractalIndexCheck(:,1)==CurrTrialIndex,[3,5]);
    
    %Only one exception for the scene sequence
    if FractalNameCurr<=1037 & FractalNameCurr>=1030
        if FractalNameCurr<=1033
            TheoreticalFractalSet=157;
        else
             TheoreticalFractalSet=156;     
        end
        
        
        
    else
    distance=SceneNameCurr-FORAGESCENE_START;
    Times=floor(distance/10);
    Remain=mod(distance,10);
    TheoreticalFractalSet=FORAGEFRACTAL_START+Times*2;
    if(Remain>3)
        TheoreticalFractalSet=TheoreticalFractalSet+1;
    end
    
    end
    
  
    TheoreticalFractal=[TheoreticalFractalSet*10+Remain,TheoreticalFractalSet*10+Remain+4];%0-3Good Object;4-7:Bad Object;
    ScreenOut=ismember(FractalNameCurr,TheoreticalFractal);
    if sum(ScreenOut,'all')~=0
        ChangeTrialIndex=sum(ScreenOut,2)>0;
        ContainCoupleIndexCurr(ChangeTrialIndex,:)=1;
        ContainCoupleIndexUpdate(ContainCoupleIndexUpdate(:,1)==CurrTrialIndex,:)=ContainCoupleIndexCurr;
        disp('Extra coupled trials found!!')
        
    end
    
end




end