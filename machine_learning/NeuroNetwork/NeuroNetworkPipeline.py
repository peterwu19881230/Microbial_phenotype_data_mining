def NeuroNetworkPipeline (NewPheno,Label):    
    
    from sklearn.model_selection import train_test_split
    from sklearn.metrics import confusion_matrix
    import pandas as pd
    from keras.layers import Dense #I did this in terminal first: 1. conda install keras 
    from keras.models import Sequential
    from keras.utils import to_categorical
    import numpy as np

    X=NewPheno
    y=Label
    
    metrics=dict.fromkeys(range(Label.shape[1]))
    
    count=0
    for col in Label: #this way all columns in the dataframe can be iterated over (https://stackoverflow.com/questions/28218698/how-to-iterate-over-columns-of-pandas-dataframe-to-run-regression/32558621)
        y=Label[col]
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.25, random_state=42,stratify=y) #not sure if stratify=y works when there are more than 2 classes for the labels
    
        #instantiate the model + fit the data
    
    
        
        ##The following code is a combination of Ref1 and Ref2
        ##Ref 1 : https://campus.datacamp.com/courses/deep-learning-in-python/building-deep-learning-models-with-keras?ex=9
        ##Ref 2 : https://www.kaggle.com/akashsri99/deep-learning-iris-dataset-keras
        model = Sequential()
        model.add(Dense(10,activation='relu',input_shape=(X.shape[1],))) #X.shape[1] should equal the number of features
        model.add(Dense(10,activation='relu'))
        model.add(Dense(2,activation='softmax')) #the number of the output layer should equal the number of unique outcomes (response variables)
        
        model.compile(optimizer='sgd',loss='categorical_crossentropy',metrics=['accuracy'])
        
        
        
        model.fit(X_train, to_categorical(y_train),epochs=5)
        
        
        y_pred = model.predict(X_test)
        y_pred = np.argmax(y_pred,axis=1)
    
    
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