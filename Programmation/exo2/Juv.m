function [f,g] = Juv(var,A,b)

    f = var'*A*var - b'*var;
    
    if nargout > 1
       g=2*A*var -b; 
    end


end

