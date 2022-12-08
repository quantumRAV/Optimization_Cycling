%% Get the Data from Prilutsky to compare agains
PriData = readtable("PrilutskyData.xlsx");

muscNames = unique(PriData.Muscle);
PriData.AdjustedForce = zeros(length(PriData.Muscle),1);
datav = zeros(length(PriData.Muscle),1);

for j = 1:length(muscNames)
    vPD = PriData(strcmp(PriData.Muscle,muscNames{j}),:);
    minF = min(vPD.RawValue);
    maxF = max(vPD.RawValue);
    gradv = (vPD.MaxV-vPD.MinV)./(maxF-minF);
    c = vPD.MaxV-gradv.*maxF;
    adjData = gradv.*(vPD.RawValue) + c;
    datav(strcmp(PriData.Muscle,muscNames{j})) = adjData;

end
PriData.AdjustedActivation = datav;
writetable(PriData,"PrilutskyData_AdjustedActivation.csv");