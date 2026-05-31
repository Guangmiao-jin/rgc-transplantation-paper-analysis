function  BI = calculateBiasIndex(peakHeight_on, peakHeight_off)
% Calculates bias index for neurons. Essentailly a modualtion index to give
% a relative strength of ON vs OFF responses, +1 for ON cells through 0 for
% ON-OFF cells to −1 for OFF cells, see:
% 
% Farrow K, Masland RH. Physiological clustering of visual channels in the
% mouse retina. J
% Neurophysiol. 2011 Apr;105(4):1516-30. doi: 10.1152/jn.00331.2010. Epub
% 2011 Jan 27. PMID: 21273316; PMCID: PMC3075295
%%
if isnumeric(peakHeight_on) && isnumeric(peakHeight_off)
    BI = (peakHeight_on - peakHeight_off)./ (peakHeight_on + peakHeight_off);

elseif iscell(peakHeight_on) &&  iscell(peakHeight_off)

    a = cellfun(@minus, peakHeight_on, peakHeight_off, 'UniformOutput',false);
    b = cellfun(@plus, peakHeight_on, peakHeight_off, 'UniformOutput',false);
    BI =cellfun(@(x,y) x./y, a, b, 'UniformOutput',false);
    
else
    error('Please check iput to calculateBiasIndex')
end

end
