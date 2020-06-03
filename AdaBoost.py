import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sklearn as sk
from sklearn import datasets 
from sklearn.model_selection import train_test_split
from sklearn.ensemble import AdaBoostRegressor
from sklearn import metrics
import csv

# Load data, select features and ignore target in the X list:
iris = datasets.load_iris()
iris = pd.read_csv("PrimaryUnfiltered.csv")
iris = iris[(iris.ClSize>2)]
X = iris[["t","PathLen","Angle","SecondMom","ClSize","D0","Phi0",
          "Chi2","Emax","EnergySpread","ShowerArea","EnergyDensity","EnergyPerCrystal"]].values
Y = iris[["dE"]].values
X =iris.drop(["dE", "evt","run"], axis=1)

#split data, fit and predict:
X_train , X_test, y_train, y_test = train_test_split(X,Y,test_size = 0.5, random_state = 100)
regr = AdaBoostRegressor(random_state=0, n_estimators=100,learning_rate=3.7)
regr.fit(X_train, y_train)
AdaBoostRegressor(n_estimators=10, random_state=0)
y_pred = regr.predict(X_test)

#Write out the prediction
pd.DataFrame(y_pred).to_csv("predicteddETest.csv", index=False)

# print results and make simple plot:
print('reg score: ',regr.score(X_test,y_test))
print('Mean Absolute Error:', metrics.mean_absolute_error(y_test, y_pred))
print('Mean Squared Error:', metrics.mean_squared_error(y_test, y_pred))
print('Root Mean Squared Error:', np.sqrt(metrics.mean_squared_error(y_test, y_pred)))
print('r2', np.sqrt(metrics.r2_score(y_test, y_pred)))

fig, ax = plt.subplots()
ax.scatter(y_test, y_pred) 

lims = [
    np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
    np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
]

ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
ax.set_aspect('equal')
ax.set_xlim(lims)
ax.set_ylim(lims)
plt.xlabel("Actual dE")
plt.ylabel("Predicted dE")
fig.show()
