import sklearn
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn import tree
import random
import seaborn as sns

# Custom functions
import lib
import importlib
importlib.reload(lib)

np.random.seed(42) # Seeding

## DATA LOADING
# Update these paths with outputs from Chromdet.
Y = pd.read_csv('./data/sample_groups.tab', sep = "\t", header = None)[1]
data = pd.read_csv('./data/samples_collapsed_filtered.tab', sep = '\t')
X = pd.DataFrame([x.split("|") for x in data.iloc[:,4]]).transpose()

############ The data, for now, is treated as if it were linear.
############ This is FALSE.
# Feature exploding like a MFC would do ?
# Or a one-hot encoding regardless.

# Decision tree
clf = tree.DecisionTreeClassifier()
clf.fit(X,Y)
lib.exportTree(clf, './results/tree', range(len(X.columns)))
# Overfits massively, as expected.

# LDA
lda = LinearDiscriminantAnalysis(n_components=3)
lda.fit(X, Y)
lda_result = lda.transform(X)

r = lambda: random.randint(0,255); palette = [('#%02X%02X%02X' % (r(),r(),r())) for i in range(13)]
plt.rcParams["figure.figsize"] = [10,10]
lib.plot_lda(X_r = lda_result, y=np.array(Y),
        target_names=np.unique(Y), title="cell_type_lymph",
        colors=palette, output_path="./results/lda.png")
lib.plot_lda_3d(X_r = lda_result, y=np.array(Y),
        target_names=np.unique(Y), title="cell_type_lymph",
        colors=palette, output_path="./results/lda_3d.png")


## ElasticNet logistic regressor
myEnClas = sklearn.linear_model.SGDClassifier(loss='log', penalty='elasticnet')
myEnClas.fit(X,Y)

## Does it predict well ?
# This is invalid because we did not split the data in training and testing, but gives a good rough first idea
prediction = myEnClas.predict(X)

plt.figure(figsize=(10,6))
sns.heatmap(pd.crosstab(prediction,Y),annot=True,cmap='Blues')
plt.savefig('./results/elastic_net_crosstab.png')

# Weights of the features : one line per class, one column per feature
features_weights = myEnClas.coef_

# TODO Try again with different random seeds, to see if the results are the same.

plt.figure(figsize=(120,15))
sns.heatmap(features_weights,cmap='PiYG')
plt.savefig('./results/elastic_net_features.png')
