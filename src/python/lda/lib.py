from sklearn import tree
import graphviz

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

def exportTree(treeToPrint, output_path, features_names):
    """
    Quick wrapper to export a tree to a given output path with default parameters
    """
    dot_data = tree.export_graphviz(treeToPrint, out_file=None,filled=True,
                                    feature_names=features_names,
                                    rounded=True, special_characters=True)
    graph_decision = graphviz.Source(dot_data)
    graph_decision.render(output_path)


def plot_lda(X_r, y, target_names, title, colors, output_path):
    """Will plot an LDA.
    This function will plot each LDA class on a different subplot."""

    plt.close('all')
    plt.figure()

    target_values = target_names
    plt.title('LDA of '+title+' dataset')

    for color, i, target_name in zip(colors, target_values, target_names):
        plt.scatter(X_r[y == i, 0], X_r[y == i, 1], alpha=.4, color=color, label=target_name)
    plt.legend(loc='best', shadow=False, scatterpoints=1)
    plt.savefig(output_path)


def plot_lda_3d(X_r, y, target_names, title, colors, output_path):
    """ Plots an LDA with at least three components"""
    plt.close('all')

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    target_values = target_names
    plt.title('3D LDA of '+title+' dataset')

    for color, i, target_name in zip(colors, target_values, target_names):
        ax.scatter(xs = X_r[y == i, 0], ys = X_r[y == i, 1],zs = X_r[y == i, 2],
         alpha=.4, c=color, label=target_name)
    plt.legend(loc='best', shadow=False, scatterpoints=1)
    plt.savefig(output_path)