function [] = Experiment1SylvesterFormat_2Polys_Batch()



arrExamples = {'1','2','3','4','5','6','7','8','9','10','11','12','13','14'};

nExamples = length(arrExamples);

for i = 1:1:nExamples
    
    ex_num = arrExamples{i};
    Experiment1SylvesterFormat_2Polys(ex_num)
end

end
