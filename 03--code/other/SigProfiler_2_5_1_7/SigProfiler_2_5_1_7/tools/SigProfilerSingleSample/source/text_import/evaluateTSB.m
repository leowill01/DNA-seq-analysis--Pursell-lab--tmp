function [Pboth, Pdir] = evaluateTSB(type2, type1, genome)

    expectedValues = round(sum((genome(type1)+genome(type2))/2));
    [Ppos, Pneg,Pboth] = Fisherextest(sum(genome(type1)), sum(genome(type2)), expectedValues, expectedValues);
    if ( sum(genome(type1)) >= sum(genome(type2)) )
        Pdir = -1;
    else
        Pdir = +1;
    end
    
end