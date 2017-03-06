function [TrapDepth,TrapLocation] = GetTrapDepth(Field)
    if max(Field) < 0 
        TrapDepth = 0;
        return;
    end
    dField = diff(Field(50:round(length(Field)/2)));
    [~,TrapLocation] = min(abs(dField));   
    if TrapLocation == length(dField) || dField(TrapLocation-1)*dField(TrapLocation+1) > 0    
        TrapDepth = 0;    
    else
        TrapDepth = abs(Field(TrapLocation+49));
    end    
end