function specificity = A_Prime(hit, falseAlarm, chance)
    % Note: There are varying ways to implement the A_Prime calculation.
    
    
    % Basic way
%    numerator = (hit-falseAlarm) .* (1+hit-falseAlarm);
%    denominator = (4*hit) .* (1-falseAlarm);
%    
%    specificity = chance + ( numerator ./ denominator );
    
    
    % "Flexible" way
%     if (hit > falseAlarm)
%         biggerValue = hit;
%         smallerValue = falseAlarm;
%     elseif (falseAlarm > hit)
%         biggerValue = falseAlarm;
%         smallerValue = hit;
%         
%     elseif (hit == falseAlarm)
%         if (hit==0 && falseAlarm==0)
%             disp('WARNING: BOTH A_PRIME ARGUMENTS EQUAL ZERO - PERFORMING SHODDY WORKAROUND');
%             hit = 0.01;
%             falseAlarm = 0.01;
%         end
%         biggerValue = hit;
%         smallerValue = falseAlarm;
%         
%     end
%     
%     numerator = (biggerValue-smallerValue) .* (1+biggerValue-smallerValue);
%     denominator = (4*biggerValue) .* (1-smallerValue);
%     
%     specificity = chance + ( numerator ./ denominator );
    
    
    % "Flexible" way; supports vectors as input arguments, but is less readable
    % (Also replicates deBettencourt et al, 2018)
    length_hit = length(hit);
    length_falseAlarm = length(falseAlarm);
    assert(length_hit==length_falseAlarm, 'WARNING: hit and falseAlarm vectors are not the same size.');
    specificity = nan(1,length_hit);
    
    for iSubject=1:length_hit % or length_falseAlarm, doesn't matter
        
        current_hit = hit(iSubject);
        current_falseAlarm = falseAlarm(iSubject);
        
        if (current_hit > current_falseAlarm)
            biggerValue = current_hit;
            smallerValue = current_falseAlarm;
        elseif (current_falseAlarm > current_hit)
            biggerValue = current_falseAlarm;
            smallerValue = current_hit;
        
        elseif (current_hit == current_falseAlarm)
            if (current_hit==0 && current_falseAlarm==0)
                disp(['WARNING: BOTH A_PRIME ARGUMENTS EQUAL ZERO - PERFORMING WORKAROUND (iSubject=', num2str(iSubject), ')']);
                current_hit = 0.01;
                current_falseAlarm = 0.01;
            end
            biggerValue = current_hit;
            smallerValue = current_falseAlarm;

        end
        
        numerator = (biggerValue-smallerValue) .* (1+biggerValue-smallerValue);
        denominator = (4*biggerValue) .* (1-smallerValue);

        specificity(iSubject) = chance + ( numerator ./ denominator );
        
    end

    
end

