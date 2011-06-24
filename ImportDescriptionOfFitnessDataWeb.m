% import description of fitness data from website
str=urlread('http://genomics.lbl.gov/supplemental/MR1fitness2011/MR1_experiments.tab');
a=textscan(str,'%s','delimiter','\n')
a=a{1};
temp=textscan(a{1},'%s','delimiter','\t');
data=cell(length(a),length(temp));
for i = 1:length(a)
    temp=textscan(a{i},'%s','delimiter','\t');
    for j = 1:length(temp{1})
        data{i,j}=temp{1}{j};
    end
end
F.chip=str2double(data(:,1));
F.group=data(:,2);
F.growthmethod=data(:,3);
F.info=data(:,4);

    