function [bp]=branchingParameter(binFirings_sumsub)
% The branchingParameter is design for calculate the branching parameter
      J=0;
      P=0;
      for I = 1:length(binFirings_sumsub)-1
          if binFirings_sumsub(I,1) >0
             J=J+1;
             P(J,1) = binFirings_sumsub(I+1,1) ./ binFirings_sumsub(I,1);
          end
      end
      bp=mean(P);
end