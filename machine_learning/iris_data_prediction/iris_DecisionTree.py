

from sklearn import datasets
from sklearn.model_selection import train_test_split
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import accuracy_score

iris = datasets.load_iris()
X = iris.data
y = iris.target



#create traning set and test set
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.25, random_state=42,stratify=y) #not sure if stratify=y works when there are more than 2 classes for the labels


#instantiate the model + fit the data



model = DecisionTreeClassifier(random_state=43)

model.fit(X_train, y_train)



#predict the test data and calculate accuracy, sensitivity...etc

y_pred=model.predict(X_test)


acc = accuracy_score(y_test, y_pred) #0.92
