classdef partialDCT
   properties
      m;
      n;
      J;
      adjoint;
   end
   methods      
        function obj = partialDCT(n,m,J)
            obj.m = m;
            obj.n = n;
            obj.J = J;
            obj.adjoint = 0;
        end        
   end
end
