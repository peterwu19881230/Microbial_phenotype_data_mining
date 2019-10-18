def SupportVectorMachinePipeline (NewPheno,Label):    
    
    from sklearn.model_selection import train_test_split
    from sklearn import svm
    from sklearn.metrics import confusion_matrix
    import pandas as pd

    X=NewPheno
    y=Label
    
    metrics=dict.fromkeys(range(Label.shape[1]))
    
    count=0
    for col in Label: #this way all columns in the dataframe can be iterated over (https://stackoverflow.com/questions/28218698/how-to-iterate-over-columns-of-pandas-dataframe-to-run-regression/32558621)
        y=Label[col]
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.25, random_state=42,stratify=y) #not sure if stratify=y works when there are more than 2 classes for the labels
    
        #instantiate the model + fit the data
    
    
        
        model = svm.SVC(kernel='linear') # Linear Kernel (Ref: https://www.datacamp.com/community/tutorials/svm-classification-scikit-learn-python)
    
        model.fit(X_train, y_train)
    
    
    
        #predict the test data and calculate accuracy, sensitivity...etc
    
        y_pred=model.predict(X_test)
    
    
        confusionMat= confusion_matrix(y_test, y_pred) 
        
        TN=confusionMat[0,0]
        FN=confusionMat[1,0]
        FP=confusionMat[0,1]
        TP=confusionMat[1,1]
        
        metrics[count]={
                "Annotation":col,
                "accuracy":(TP+TN)/(TN+FN+FP+TP),
                "sensitivity":TP/(TP+FN),
                "precision":TP/(TP+FP), 
                #when TP+FP=0 -> precision=NaN, there will be RunrimeWarning, but I think it's not a big problem
                
                
                "specificity":TN/(TN+FP),
                }
        
        
        count=count+1
    
    
    
    metrics_df=pd.DataFrame(metrics).transpose()
    metrics_df=metrics_df.sort_values(by=["precision","sensitivity"],ascending=False)
    
    #RuntimeWarning: invalid value encountered in long_scalars => (verified) this happens when TP+FP=0 => -> precision=NaN
    
    return(metrics_df)