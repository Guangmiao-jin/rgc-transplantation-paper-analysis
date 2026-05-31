function ISI_CV = calculateISI_CV(trSpikesStruct)
% Calcuate the interspike interval coeffient of varation based on the event
% aligned spike structure (usually ON time and OFF time)

for bl = 1:length(trSpikesStruct)
    ISI_vec =[];
    for tr = 1:length(trSpikesStruct{bl})
        s = trSpikesStruct{bl}{tr};
        ISI_vec = [ISI_vec; diff(s)];
    end

    ISI_CV(bl) = std(ISI_vec,'omitnan')/mean(ISI_vec, 'omitnan'); % Coefficient of variance
end
end