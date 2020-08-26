function [dataname] = enterData(data_num)
%% choose a network to load from networks folder
clc
disp('***DISPLAY RESULTS***')
dirName = [pwd,'\simulations\ALL*.mat'];
Allinpnames = dir(dirName);

if isempty(data_num)
    disp(sprintf('\nChoose data file:'))
    for i=1:length(Allinpnames)
        disp([num2str(i),'. ', Allinpnames(i).name])
    end
    x = input(sprintf('\nEnter Data Number: '));
else
    x = data_num;
end
dataname=['\simulations\',Allinpnames(x).name];
end

