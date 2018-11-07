classdef mtxIDCT
   properties
      m;
      n;
      eig;
      adjoint;
   end
   methods      
        function obj = mtxIDCT(n)
            obj.m = n;
            obj.n = n;
            obj.eig = 1;
            obj.adjoint = 0;
        end        
   end
end
