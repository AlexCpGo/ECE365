# Name: Yuhang Chen
# Netid: yuhangc3
import time
import numpy as np
import scipy.spatial.distance as dist
from scipy import stats
from sklearn import neighbors
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

class Question1(object):
    # Without loop
    def bayesClassifier(self,data,pi,means,cov):
        # using dot and sum to make delta_y
        C_inv = np.linalg.inv(cov)
        delta_y = np.log(pi) + np.dot(means, C_inv.dot(data.T)).T - 1 / 2 * np.sum(np.dot(means, C_inv) * means, axis=1)
        labels = np.argmax(delta_y, axis=1)
        return labels

    def classifierError(self,truelabels,estimatedlabels):
        error = np.sum(truelabels != estimatedlabels) / len(truelabels)
        return error


class Question2(object):
    def trainLDA(self,trainfeat,trainlabel):
        nlabels = int(trainlabel.max())+1 # Assuming all labels up to nlabels exist.
        pi = np.zeros(nlabels)            # Store your prior in here
        means = np.zeros((nlabels,trainfeat.shape[1]))            # Store the class means in here
        cov = np.zeros((trainfeat.shape[1],trainfeat.shape[1]))   # Store the covariance matrix in here
        # Put your code below
        # pi: (1, 3), means: (3, 2), cov: (2, 2)
        for i in range(nlabels):
            pi[i] = trainfeat[trainlabel == i].shape[0] / trainfeat.shape[0]
            means[i] = np.sum(trainfeat[trainlabel == i], axis=0) / trainfeat[trainlabel == i].shape[0]
        for i in range(nlabels):
            cov += np.dot((trainfeat[trainlabel == i] - means[i]).T, (trainfeat[trainlabel == i] - means[i]))
        cov /= trainfeat.shape[0] - means.shape[0]  

        # Don't change the output!
        return (pi,means,cov)

    def estTrainingLabelsAndError(self,trainingdata,traininglabels):
        q1 = Question1()
        # You can use results from Question 1 by calling q1.bayesClassifier(...), etc.
        # If you want to refer to functions under the same class, you need to use self.fun(...)
        pi, means, cov = self.trainLDA(trainingdata, traininglabels)
        esttrlabels = q1.bayesClassifier(trainingdata, pi, means, cov)
        trerror = q1.classifierError(traininglabels, esttrlabels)
        # Don't change the output!
        return (esttrlabels, trerror)

    def estValidationLabelsAndError(self,trainingdata,traininglabels,valdata,vallabels):
        q1 = Question1()
        # You can use results from Question 1 by calling q1.bayesClassifier(...), etc.
        # If you want to refer to functions under the same class, you need to use self.fun(...)
        pi, means, cov = self.trainLDA(trainingdata, traininglabels)
        estvallabels = q1.bayesClassifier(valdata, pi, means, cov)
        valerror = q1.classifierError(vallabels, estvallabels)
        # Don't change the output!
        return (estvallabels, valerror)



class Question3(object):
    def kNN(self,trainfeat,trainlabel,testfeat, k):
         # Using trainfeat: (N, d), trainlabel: (N, 1), testfeat: (V, d) to train and test
        distance = dist.cdist(testfeat, trainfeat, metric='euclidean')  
        index = np.argpartition(distance, k-1)[:, :k]  

        # (V, k)--k_labels
        k_labels = np.zeros((distance.shape[0], k))
        for i in range(distance.shape[0]):
            for j in range(k):
                k_labels[i][j] = trainlabel[index[i][j]]

        labels = stats.mode(k_labels, axis=1)[0][:, 0]
        return labels

    def kNN_errors(self,trainingdata, traininglabels, valdata, vallabels):
        q1 = Question1()
        trainingError = np.zeros(4)
        validationError = np.zeros(4)
        k_array = [1,3,4,5]

        for i in range(len(k_array)):
            # Please store the two error arrays in increasing order with k
            # This function should call your previous self.kNN() function.
            # Put your code below
            labels = self.kNN(trainingdata, traininglabels, trainingdata, k_array[i])
            trainingError[i] = q1.classifierError(traininglabels, labels)

            labels = self.kNN(trainingdata, traininglabels, valdata, k_array[i])
            validationError[i] = q1.classifierError(vallabels, labels)


        # Don't change the output!
        return (trainingError, validationError)

class Question4(object):
    def sklearn_kNN(self,traindata,trainlabels,valdata,vallabels):
        # Make record to all of the method used
        # method = {'ball_tree', 'kd_tree', 'brute', 'auto'}
        # ball_tree: 4.6836957931518555 s
        # kd_tree: 5.789940357208252 s
        # brute: 0.18962836265563965 s
        # auto: 5.8512420654296875 s
        # >>> brute is the fastest
        classifier = neighbors.KNeighborsClassifier(n_neighbors=1, algorithm='brute')
        start = time.time()
        classifier.fit(traindata, trainlabels)
        mid = time.time()
        prediction = classifier.predict(valdata)
        end = time.time()
        valerror = np.sum(prediction != vallabels) / vallabels.shape[0]
        fitTime = mid - start
        predTime = end - mid
        # total_time = end - start
        # Don't change the output!
        return (classifier, valerror, fitTime, predTime)

    def sklearn_LDA(self,traindata,trainlabels,valdata,vallabels):
        classifier = LinearDiscriminantAnalysis()
        start = time.time()
        classifier.fit(traindata, trainlabels)
        mid = time.time()
        prediction = classifier.predict(valdata)
        end = time.time()
        valerror = np.sum(prediction != vallabels) / vallabels.shape[0]
        fitTime = mid - start
        predTime = end - mid
        # Don't change the output!
        return (classifier, valerror, fitTime, predTime)

###
