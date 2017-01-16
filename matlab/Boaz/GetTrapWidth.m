function TrapWidth = GetTrapWidth(potential,threshold)

    [C,~] = max(potential);
    first = find(potential>threshold*C,1);
    last = find(potential>threshold*C,1,'last');
    TrapWidth = last-first;
end