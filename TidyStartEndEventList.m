function [I,I2] = TidyStartEndEventList(I,I2,N)

    if ~isempty(I)
        if I(1)>I2(1)
            I=[1;I];
        end
        if I(end)>I2(end)
            I2=[I2;N];
        end
    end
    