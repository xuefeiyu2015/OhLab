function BonusTrialCondition=BonusTrialConditionRescue(SceneIndexBonusUse,TargetIndexBonusUse)
%{
//1.Trained Scene Within---Coupled Objects v.s. UnCoupled(Within)
	//2.Trained Scene Within---Coupled Objects v.s. UnCoupled(Between)
	//3.Trained Scene Within--- Coupled Objects v.s. Novel Objects
	//4.Trained Scene Within---UnCoupled Objects(Within) v.s. UnCoupled(Between)
	//5.Trained Scene Within---Novel v.s. Novel

	//6.Trained Scene Between---Trained Objects(Within) v.s. Novel Objects
	//7.Trained Scene Between---Trained Objects(Between) v.s. Novel Objects
	//8.Trained Scene Between---Trained Objects(Within) v.s. Trained(Between)
	//9.Trained Scene Between---Novel v.s. Novel

	//10.Novel Scene--Trained Objects(Within) v.s. Novel Objects
	//11.Novel Scene---Trained Objects(Between) v.s. Novel Objects
	//12.Novel Scene---Trained Objects(Between) v.s. Trained Objects(Within)
	//13.Novel Scene---Novel Objects v.s. Novel Objects

%}

%This is to rescue bonus trial condition when there is an error
%Xuefei Yu 02212020

NovelFractalStart=80;

for i=1:size(TargetIndexBonusUse,1)
   
    
    currSceneIndex=SceneIndexBonusUse(ismember(SceneIndexBonusUse(:,1),TargetIndexBonusUse(i,1)),2);
    currFractalIndex=TargetIndexBonusUse(i,[2,4]);
    %Make this hardwired!!
    if currSceneIndex<=12 & currSceneIndex>=7
        %Novel
        if sum(currFractalIndex>=NovelFractalStart & currFractalIndex<93)==2
            %Both are novel
            currCondition(i)=13;
        elseif sum(currFractalIndex>=NovelFractalStart & currFractalIndex<93)==1
            %One novel
            Other=currFractalIndex(~(currFractalIndex>=NovelFractalStart & currFractalIndex<93));
            if Other>=13 & Other<=36
                %Between set
                currCondition(i)=11;
            else
                %Within set
                currCondition(i)=10;
            end
        elseif sum(currFractalIndex>=NovelFractalStart & currFractalIndex<93)==0
            %Trained Objects(Between) v.s. Trained Objects(Within)
            currCondition(i)=12;
        end
      
        
        
        
    elseif currSceneIndex<7 
        %Between set, trained
        
        if sum(currFractalIndex>=NovelFractalStart & currFractalIndex<93)==2
           %Both are novel
            currCondition(i)=9; 
        elseif sum(currFractalIndex>=NovelFractalStart & currFractalIndex<93)==1
            %One novel
            Other=currFractalIndex(~(currFractalIndex>=NovelFractalStart & currFractalIndex<93));
            if Other>=13 & Other<=36
                %Between set
                currCondition(i)=7;
            else
                %Within set
                currCondition(i)=6;
            end
        elseif sum(currFractalIndex>=NovelFractalStart & currFractalIndex<93)==0
            %Trained Objects(Between) v.s. Trained Objects(Within)
            currCondition(i)=8;
        
        end
        
        
    else
        %Within set, trained
        if sum(currFractalIndex>=NovelFractalStart & currFractalIndex<93)==2
           %Both are novel
            currCondition(i)=5; 
        elseif sum(currFractalIndex>=NovelFractalStart & currFractalIndex<93)==1
            %Only one novel
            currCondition(i)=3;
        elseif sum(currFractalIndex>=13 & currFractalIndex<=36)==1
            %One between set
            Other=currFractalIndex(~(currFractalIndex>=13 & currFractalIndex<=36));
            %Within set
            %Hard wired it 
            SceneFractalIndexLib=[[37:48]',[51:2:73]',[50:2:72]'];
            TheoreticalIndex=SceneFractalIndexLib(SceneFractalIndexLib(:,1)==currSceneIndex,2:end);
            if ismember(Other,TheoreticalIndex)
                %Coupled
                currCondition(i)=2;
            else
                %Uncoupled within
                 currCondition(i)=4;
                
            end
        elseif sum(currFractalIndex>=50 & currFractalIndex<=73)==2
            %Both are within
            currCondition(i)=1;
            
            
            
           
            
        end
        
    end
    
    
    
end

BonusTrialCondition=[TargetIndexBonusUse(:,1),currCondition'];

end