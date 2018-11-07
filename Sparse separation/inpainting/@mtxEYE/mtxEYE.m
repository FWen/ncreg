classdef mtxEYE
   properties
      m;
      n;
      eig;
      adjoint;
   end
   methods      
        function obj = mtxEYE(n)
            obj.m = n;
            obj.n = n;
            obj.eig = 1;
            obj.adjoint = 0;
        end        
   end
end
