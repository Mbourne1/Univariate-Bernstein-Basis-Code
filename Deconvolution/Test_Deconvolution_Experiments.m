function [] = Test_Deconvolution_Experiments

arr_ex_num = {'1','2','3','4','5','6'};
arr_noise = {1e-8,1e-10,1e-12};
arr_bool_preproc = {'y','n'};


for i1 = 1:1:length(arr_ex_num)
    ex_num = arr_ex_num{i1};
    for i2 = 1:1:length(arr_noise)
        emin = arr_noise{i2};
        for i3 = 1:1:length(arr_bool_preproc)
            bool_preproc = arr_bool_preproc{i3};
            
            Test_Deconvolution(ex_num,emin,bool_preproc)
            
        end
    end
end


end