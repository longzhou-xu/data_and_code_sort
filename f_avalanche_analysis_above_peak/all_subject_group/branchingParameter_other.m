function [BranchingParameter]=branchingParameter_other(BinEvent)
%%
% The branchingParameter is design for calculate the branching parameter

      EN=0;
        ratio_1=0;
        for I = 1:length(BinEvent)-1
            if BinEvent(I,1) >0
               EN=EN+1;
               ratio_1(EN,1) = BinEvent(I+1,1) ./ BinEvent(I,1);
            end
        end
        BranchingParameter(1)=mean(ratio_1);
        
        EN=0;
        ratio_2=0;
        for I = 1:length(BinEvent)-1
            if BinEvent(I,1) == 1
               EN=EN+1;
               ratio_2(EN,1) = BinEvent(I+1,1) ./ BinEvent(I,1);
            end
        end
        BranchingParameter(2)=mean(ratio_2);
        
        EN=0;
        ratio_3=0;
        for I = 2:length(BinEvent)-1
            if BinEvent(I,1) >0 && BinEvent(I-1,1) == 0
               EN=EN+1;
               ratio_3(EN,1) = BinEvent(I+1,1) ./ BinEvent(I,1);
            end
        end
        BranchingParameter(3) = mean(ratio_3);
        
        EN=0;
        ratio_4=0;
        for I = 2:length(BinEvent)-1
            if BinEvent(I,1) ==1 && BinEvent(I-1,1) == 0
               EN=EN+1;
               ratio_4(EN,1) = BinEvent(I+1,1) ./ BinEvent(I,1);
            end
        end
        BranchingParameter(4) = mean(ratio_4);
      
end