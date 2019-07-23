%{
MATLAB's existing SVM (support vector machine) model.
}%
t = 200;
data = squeeze(cat(3,l_eeg(:,t,:),r_eeg(:,t,:)))'; trueLabels = [ones(size(l_eeg,3),1); ...
              2*ones(size(r_eeg,3),1)];
svmModel = fitcsvm(data,trueLabels);
catlabel = predict(svmModel,data);
accu = mean(catlabel==trueLabels);
traindata = data;
traindata(triali,:) = [];
templabels = trueLabels;
templabels(triali) = [];
svmModel = fitcsvm(traindata,templabels); 
catLabel = predict(svmModel,data(triali,:)); 
accu(ti,triali) = catLabel==trueLabels(triali);
