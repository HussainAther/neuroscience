%{
MATLAB's existing SVM (support vector machine) model.
}%
t = 200;
data = squeeze(cat(3,l_eeg(:,t,:),r_eeg(:,t,:)))'; trueLabels = [ones(size(l_eeg,3),1); ...
              2*ones(size(r_eeg,3),1)];
svmModel = fitcsvm(data,trueLabels);
