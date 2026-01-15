function result = UMersenneTwisterRNG(varargin)
    persistent MT
    persistent index

    function initialize_generator(seed)
        index = 0;

        MT = zeros(1, 624);
        MT(1) = uint64(seed);

        for k = 1:623
            MT(k+1) = bitand(...
                uint64(1812433253 * bitxor(MT(k), bitshift(MT(k), -30)) + k), ...
                4294967295);
        end
    end

    function generate_numbers()
        for k = 0:623
            y = bitand(MT(k + 1), 2147483648) + ...
                bitand(2147483647, MT( mod(k + 1, 624) + 1));
                
            MT(k + 1) = bitxor(MT(mod(k + 397, 624) + 1), bitshift(uint32(y), -1));
            
            if (mod(y, 2) ~= 0)
                MT(k + 1) = bitxor(MT(k + 1), 2567483615);
            end
        end
    end

    function y = extract_number() 
        if (index == 0)
            generate_numbers();
        end

        y = MT(index + 1);
        y = bitxor(y, bitshift(y, -11));
        y = bitxor(y, bitand(bitshift(y, 7), 2636928640));
        y = bitxor(y, bitand(bitshift(y, 15), 4022730752));
        y = bitxor(y, bitshift(y, -18));

        index = mod(index + 1, 624);
    end

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
                error("Wrong column number!");
            end
        else
            if (optional_input_argument_count == 2)
                if (isscalar(varargin{1}))
                    row_count = round(varargin{1});
                else
                    error("Wrong row number!");
                end

                if (isscalar(varargin{2}))
                    column_count = round(varargin{2});
                else
                    error("Wrong column number!");
                end
            else
                error("Too many input arguments!");
            end
        end
    end

    if (isempty(MT) && isempty(index))
        initialize_generator(6199);
    end

    result = zeros(row_count, column_count);

    for i = 1:row_count
        for j = 1:column_count
            result(i, j) = extract_number() / 4294967295;
        end
    end
end