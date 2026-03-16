function loc=FindFixLoc(EyeX,EyeY,fixtime);
%Function to find out the fixation location according to fixation times
%Xuefei Yu; 04302020
loc=[];
for i=1:size(fixtime,2)
    
    curr_fixtime=fixtime(:,i);
    loc_x=[];
    loc_y=[];
    
    loc_x=EyeX(curr_fixtime(1):curr_fixtime(2));
    loc_y=EyeY(curr_fixtime(1):curr_fixtime(2));
    
    
    loc{i}=[loc_x;loc_y];
    
    
    
    
end



end