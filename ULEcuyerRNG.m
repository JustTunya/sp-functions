function result = ULEcuyerRNG(varargin)
    optional_input_argument_count = size(varargin,2);
    
    if (optional_input_argument_count == 0)
        row_count = 1;
        column_count = 1;
    else
        if (optional_input_argument_count == 1)
            row_count = 1;
            
            if (isscalar(varargin{1}))
                column_count = round(varargin{1});
            else
                error('Wrong column number!');
            end
        else
            if (optional_input_argument_count == 2)
                if (isscalar(varargin{1}))
                    row_count = round(varargin{1});
                else
                    error('Wrong row number!');
                end
    
                if (isscalar(varargin{2}))
                    column_count = round(varargin{2});
                else
                    error('Wrong column number!');
                end
            else
                error('Too many input arguments!');
            end
        end
    end
    
    persistent seed1 seed2
    
    if (isempty(seed1))
        seed1 = 55555;
    end
    
    if (isempty(seed2))
        seed2 = 99999;
    end
    
    persistent factor
    
    if (isempty(factor))
        factor = 1.0/2147483563.0;
    end
    
    result = zeros(row_count, column_count);
    
    for i = 1:row_count
        for j = 1:column_count
        
            k = seed1 / 53668;
            seed1 = 40014 * mod(seed1, 53668) - k * 12211;
    
            if (seed1 < 0) 
                seed1 = seed1 + 2147483563;
            end
    
            k = seed2 / 52774;
    
            seed2 = 40692 * mod(seed2, 52774) - k * 3791;
    
            if (seed2 < 0) 
                seed2 = seed2 + 2147483399;
            end
    
            z = (seed1 - 2147483563) + seed2;
    
            if (z < 1) 
                z = z + 2147483562;
            end
    
            result(i, j) = z * factor;
        end
    end
end