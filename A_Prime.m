function specificity = A_Prime(hit, falseAlarm, chance)
    % There are varying ways to implement the A_Prime calculation, but this is the basic way to do it.
    
    numerator = (hit-falseAlarm) .* (1+hit-falseAlarm);
    denominator = (4*hit) .* (1-falseAlarm);
    
    specificity = chance + ( numerator ./ denominator );
    
end

