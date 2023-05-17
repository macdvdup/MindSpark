clc; clear all; close all;

% Step 1: Load and split the dataset
features=readtable('features.csv'); % Replace 'your_dataset.mat' with your actual dataset
X = table2array(features(:,1:end-1)); % Replace 'features' with the appropriate variable name
y = features.Class; % Replace 'labels' with the appropriate variable name

% Step 2: Preprocess the features if necessary
% (e.g., feature scaling using z-score normalization)

kFold = 5; % Set the ratio for the training set
% Step 3: Split the dataset into training and testing sets
% With Cross validation
pt = cvpartition(y, 'KFold', kFold);
accs = zeros(kFold,1);
bestAcc=0;
for fold=1:kFold
    trainIndex = training(pt,fold);
    testIndex = test(pt,fold);
    X_train = X(trainIndex,:);
    y_train  = y(trainIndex);

    [X_train_norm,mu,sigma] = zscore(X_train, 0, 1);
    %X_test = X(testIndex,:);

    X_test = (X(testIndex,:)-mu)./sigma;
    y_test  = y(testIndex);

    % Step 4: Train Model
    classifier = fitcknn(X_train_norm, y_train,'Distance','cityblock','NumNeighbors',5);
    %classifier.ClassNames = {"Natural","Stress"}
    % Step 5: Evaluate the model using the testing data
    predictedY = predict(classifier, X_test);
    
    accs(fold,1) = sum(predictedY == y_test) / numel(y_test);
    if accs(fold,1)>bestAcc
        bestClassifier = classifier;
        bestAcc = accs(fold,1);
    end
end

disp(['Accuracy ' num2str(mean(accs))]);
save('bestModelBoth.mat','bestClassifier')
save('miu.mat','mu')
save('sigma.mat','sigma')