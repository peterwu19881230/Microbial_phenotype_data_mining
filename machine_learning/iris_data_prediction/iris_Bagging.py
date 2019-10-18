

from sklearn import datasets
from sklearn.model_selection import train_test_split
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import BaggingClassifier
from sklearn.metrics import accuracy_score

iris = datasets.load_iris()
X = iris.data
y = iris.target



#create traning set and test set
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.25, random_state=42,stratify=y) #not sure if stratify=y works when there are more than 2 classes for the labels


#instantiate the model + fit the data



model = DecisionTreeClassifier(random_state=43)


ensemble=BaggingClassifier(base_estimator=model, n_estimators=300, n_jobs=-1) #I tried to increase n_estimators to 3000, but accuracy is still 0.947

ensemble.fit(X_train, y_train)



#predict the test data and calculate accuracy, sensitivity...etc

y_pred=ensemble.predict(X_test)


acc = accuracy_score(y_test, y_pred) #0.947 (If only decision tree is used, acc=0.92)
