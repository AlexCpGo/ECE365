import numpy as np
from sklearn import neighbors
import scipy.spatial.distance as dist
from sklearn import linear_model
from sklearn.model_selection import train_test_split

class Question1(object):
    def kMeans(self,data,K,niter):
        """ Implement the K-Means algorithm.

        **For grading purposes only:**

        Do NOT change the random seed, otherwise we are not able to grade your code! This is true throughout this script. However, in practice, you never want to set a random seed like this.
        For your own interest, after you have finished implementing this function, you can change the seed to different values and check your results.
        Please use numpy library for initial random choice. This will use the seed above. Scipy library is using a different seeding system, so that would probably result in an error during our grading.

        Parameters:
        1. data     (N, d) numpy ndarray. The unlabelled data with each row as a feature vector.
        2. K        Integer. It indicates the number of clusters.
        3. niter    Integer. It gives the number of iterations. An iteration includes an assignment process and an update process.

        Outputs:
        1. labels   (N,) numpy array. It contains which cluster (0,...,K-1) a feature vector is in. It should be the (niter+1)-th assignment.
        2. centers  (K, d) numpy ndarray. The i-th row should contain the i-th center.
        """
        np.random.seed(12312)
        # Put your code below
        cluster_idx=np.random.choice(data.shape[0] , K, replace=False)
        centers=data[cluster_idx,:]
        labels=np.zeros(data.shape[0])
        # Use loop to produce the new center for each iterating time of the k-Means clustering
        for i in range(niter):
            dis=dist.cdist(data,centers,'euclidean') #Nxk
            labels=np.argmin(dis,axis=1)
            for j in range(K):
                new_center=data[labels==j].mean(axis=0)
                centers[j]=new_center

        # Remember to check your data types: labels should be integers!
        return (labels, centers)

    def calculateJ(self,data):
        """ Calculate the J_k value for K=2,...,10.

        This function should call your self.kMeans() function and set niter=100.

        Parameters:
        1. data     (N, d) numpy ndarray. The unlabelled data with each row as a feature vector.

        Outputs:
        1. err      (9,) numpy array. The i-th element contains the J_k value when k = i+2.
        """
        err = np.zeros(9)
        # Put your code below
        niter=100
        for iK in range(2,11):
            (labels,centers)=self.kMeans(data,iK,niter)
            for j in range(iK):
                err[iK-2]+=((data[labels==j]-centers[j])**2).sum()
        return err

from sklearn.cluster import KMeans
from sklearn.neighbors import NearestNeighbors

class Question2(object):
    def trainVQ(self,image,B,K):
        """ Generate a codebook for vector quantization.

        You can use the KMeans function from the sklearn package.

        **For grading purposes only:**

        Do NOT change the random seed, otherwise we are not able to grade your code!
        Please flatten any matrix in *row-major* order. If you prefer, you can use np.flatten(xxx) to flatten your matrix.

        Parameters:
        1. image        (N, M) numpy ndarray. It represents a grayscale image.
        2. B            Integer. You will use B×B blocks for vector quantization. You may assume that both N and M are divisible by B.
        3. K            Integer. It gives the size of your codebook.

        Outputs:
        1. codebook     (K, B^2) numpy ndarray. It is the codebook you should return.
        """
        np.random.seed(12345)
        # Put your code below
        classifier=KMeans(n_clusters=K,init='k-means++')
        N=image.shape[0]
        M=image.shape[1]
        img_data=np.zeros((N*M//B**2,B**2))
        for i in range(N):
            for j in range (M):
                img_data[i//B*M//B+j//B][i%B*B+j%B]=image[i][j]
        
        classifier.fit(img_data)
        codebook=classifier.cluster_centers_
        return codebook

    def compressImg(self,image,codebook,B):
        """ Compress an image using a given codebook.

        You can use the nearest neighbor classifier from scikit-learn if you want (though it is not necessary) to map blocks to their nearest codeword.

        **For grading purposes only:**

        Please flatten any matrix in *row-major* order. If you prefer, you can use np.flatten(xxx) to flatten your matrix.

        Parameters:
        1. image        (N, M) numpy ndarray. It represents a grayscale image. You may assume that both N and M are divisible by B.
        2. codebook     (K, B^2) numpy ndarray. The codebook used in compression.
        3. B            Integer. Block size.

        Outputs:
        1. cmpimg       (N//B, M//B) numpy ndarray. It consists of the indices in the codebook used to approximate the image.
        """
        # Put your code below
        classifier=NearestNeighbors(n_neighbors=1, algorithm='auto')
        classifier.fit(codebook)
        
        N=image.shape[0]
        M=image.shape[1]
        img_data=np.zeros((N*M//B**2,B**2))
        img_data=np.zeros((N*M//B**2,B**2))
        for i in range(N):
            for j in range (M):
                img_data[i//B*M//B+j//B][i%B*B+j%B]=image[i][j]
        idx=classifier.kneighbors(img_data,return_distance=False)
        cmpimg=np.resize(idx,(N//B,M//B))
        # Check that your indices are integers!
        return cmpimg

    def decompressImg(self,indices,codebook,B):
        """ Reconstruct an image from its codebook.

        You can use np.reshape() to reshape the flattened array.

        Parameters:
        1. indices      (N//B, M//B) numpy ndarray. It contains the indices of the codebook for each block.
        2. codebook     (K, B^2) numpy ndarray. The codebook used in compression.
        3. B            Integer. Block size.

        Outputs:
        1. rctimage     (N, M) numpy ndarray. It consists of the indices in the codebook used to approximate the image.
        """
        # Put your code below
        N=indices.shape[0]
        M=indices.shape[1]
        rctimage=np.zeros((N*B,M*B))
        for i in range(N):
            for j in range(M):
                rctimage[i*B:(i+1)*B,j*B:(j+1)*B]=codebook[indices[i][j]].reshape(B,-1)
        return rctimage

from sklearn.metrics import zero_one_loss
class Question3(object):
    def generatePrototypes(self,traindata,trainlabels,K_list):
        """ Generate prototypes from labeled data.

        You can use the KMeans function from the sklearn package.

        **For grading purposes only:**

        Do NOT change the random seed, otherwise we are not able to grade your code!

        Parameters:
        1. traindata        (Nt, d) numpy ndarray. The features in the training set.
        2. trainlabels      (Nt,) numpy array. The labels in the training set.
        3. K_list           List. A list of integers corresponding to the number of prototypes under each class.

        Outputs:
        1. proto_dat_list   A length len(K_list) list. The K-th element in the list is a (K * num_classes, d) numpy ndarray, representing the prototypes selected if using K prototypes under each class. You should keep the order as in the given K_list.
        2. proto_lab_list   A length len(K_list) list. The K-th element in the list is a (K * num_classes,) numpy array, representing the corresponding labels if using K prototypes under each class. You should keep the order as in the given K_list.
        """
        np.random.seed(56789)   # As stated before, do NOT change this line!
        proto_dat_list = []
        proto_lab_list = []

        # Put your code below
        label_num=len(np.unique(trainlabels))
        label=np.unique(trainlabels)
        for e in label:
            e=int(e)
        for K in K_list:
            proto_dat=np.zeros((K*label_num,traindata.shape[1]))
            proto_lab=np.zeros((K*label_num))
            for i in range(label_num):
                idx=(trainlabels==label[i])
                proto_lab[i*K:(i+1)*K]=np.repeat(label[i],K)
                classifier=KMeans(n_clusters=K,init='k-means++').fit(traindata[idx])
                proto_dat[i*K:(i+1)*K,:]=classifier.cluster_centers_
            proto_lab=proto_lab.astype(int)
            proto_dat_list.append(proto_dat)
            proto_lab_list.append(proto_lab)
        # Check that your proto_lab_list only contains integer arrays!
        return (proto_dat_list, proto_lab_list)

    def protoValError(self,proto_dat_list,proto_lab_list,valdata,vallabels):
        """ Generate prototypes from labeled data.

        You may assume there are at least min(K_list) examples under each class. set(trainlabels) will give you the set of labels.

        Parameters:
        1. proto_dat_list   A list of (K * num_classes, d) numpy ndarray. A list of prototypes selected. This should be one of the outputs from your previous function.
        2. proto_lab_list   A list of (K * num_classes,) numpy array. A list of corresponding labels for the selected prototypes. This should be one of the outputs from your previous function.
        3. valdata          (Nv, d) numpy ndarray. The features in the validation set.
        4. vallabels        (Nv,) numpy array. The labels in the validation set.

        Outputs:
        1. proto_err        (len(proto_dat_list),) numpy ndarray. The validation error for each K value (in the same order as the given K_list).
        """
        proto_err = np.zeros(len(proto_dat_list))
        # Put your code below
        for i in range(len(proto_lab_list)):
            classifier=neighbors.KNeighborsClassifier(n_neighbors=1)
            data=proto_dat_list[i]
            label=proto_lab_list[i]
            classifier.fit(data,label)
            est=classifier.predict(valdata)
            err=zero_one_loss(est,vallabels)
            proto_err[i]=err


        return proto_err

class Question4(object):
    def benchmarkRSS(self,trainfeat,trainresp,valfeat,valresp):
        """ Return the benchmark RSS.

        In particular, always predict the response as zero (mean response on the training data).

        Calculate the validation RSS for this model. Please use the formula as defined in the jupyter notebook.

        Parameters:
        1. trainfeat    (Nt, d) numpy ndarray. The features in the training set.
        2. trainresp    (Nt,) numpy array. The responses in the training set.
        3. valfeat      (Nv, d) numpy ndarray. The features in the validation set.
        4. valresp      (Nv,) numpy array. The responses in the validation set.

        Outputs:
        1. rss          Scalar. The validation RSS.
        """
        # Put your code below
        rss=np.linalg.norm(0-valresp)**2.0
        rss=rss/len(valresp)
        return rss

    def OLSRSS(self,trainfeat,trainresp,valfeat,valresp):
        """ Return the RSS from the ordinary least squares model.

        Use sklearn.linear_model.LinearRegression() with the default parameters.

        Calculate the validation RSS for this model. Please use the formula as defined in the jupyter notebook.

        Note: The .score() method returns an  R^2 value, not the RSS, so you shouldn't use it anywhere in this problem.

        Parameters:
        1. trainfeat    (Nt, d) numpy ndarray. The features in the training set.
        2. trainresp    (Nt,) numpy array. The responses in the training set.
        3. valfeat      (Nv, d) numpy ndarray. The features in the validation set.
        4. valresp      (Nv,) numpy array. The responses in the validation set.

        Outputs:
        1. rss          Scalar. The validation RSS.
        """
        # Put your code below
        model=linear_model.LinearRegression()
        model.fit(trainfeat, trainresp) 
        est=model.predict(valfeat)
        rss=np.linalg.norm(est-valresp)**2.0
        rss=rss/len(valresp)
        return rss

    def RidgeRSS(self,trainfeat,trainresp,valfeat,valresp):
        """ Return the RSS from the ridge regression.

        Apply ridge regression with sklearn.linear_model.Ridge. Sweep the regularization/tuning parameter α = 0,...,100 with 1000 equally spaced values.

        Note: Larger values of α shrink the weights in the model more, and α=0 corresponds to the LS solution.

        Parameters:
        1. trainfeat    (Nt, d) numpy ndarray. The features in the training set.
        2. trainresp    (Nt,) numpy array. The responses in the training set.
        3. valfeat      (Nv, d) numpy ndarray. The features in the validation set.
        4. valresp      (Nv,) numpy array. The responses in the validation set.

        Outputs:
        1. rss_array    (1000,). The validation RSS array. This is used for plotting. This will not be tested by the autograder.
        2. best_a       Scalar. The alpha that minimizes the RSS.
        3. best_rss     Scalar. The corresponding RSS.
        4. coef         (d,) numpy array. The minimizing coefficient. This is for visualization only. This will not be tested by the autograder.
        """
        a = np.linspace(0,100,1000)
        rss_array = np.zeros(a.shape)
        # Put your code below
        for i in range(len(a)):
            model=linear_model.Ridge(alpha=a[i])
            model.fit(trainfeat, trainresp) 
            est=model.predict(valfeat)
            rss=np.linalg.norm(est-valresp)**2
            rss=rss/len(valresp)
            rss_array[i]=rss

        idx=np.argmin(rss_array)
        best_a=a[idx]
        best_rss=rss_array[idx]
        best_model=linear_model.Ridge(alpha=best_a)
        best_model.fit(trainfeat, trainresp)
        coef=best_model.coef_

        return (rss_array, best_a, best_rss, coef)

    def LassoRSS(self,trainfeat,trainresp,valfeat,valresp):
        """ Return the RSS from the Lasso regression.

        Apply lasso regression with sklearn.linear_model.Lasso. Sweep the regularization/tuning parameter α = 0,...,1 with 1000 equally spaced values.

        Note: Larger values of α will lead to sparser solutions (i.e. less features used in the model), with a sufficiently large value of α leading to a constant prediction. Small values of α are closer to the LS solution, with α=0 being the LS solution.

        Parameters:
        1. trainfeat    (Nt, d) numpy ndarray. The features in the training set.
        2. trainresp    (Nt,) numpy array. The responses in the training set.
        3. valfeat      (Nv, d) numpy ndarray. The features in the validation set.
        4. valresp      (Nv,) numpy array. The responses in the validation set.

        Outputs:
        1. rss_array    (1000,). The validation RSS array. This is used for plotting. This will not be tested by the autograder.
        2. best_a       Scalar. The alpha that minimizes the RSS.
        3. best_rss     Scalar. The corresponding RSS.
        4. coef         (d,) numpy array. The minimizing coefficient. This is for visualization only. This will not be tested by the autograder.
        """
        a = np.linspace(0.00001,1,1000)     # Since 0 will give an error, we use 0.00001 instead.
        rss_array = np.zeros(a.shape)
        # Put your code below
        for i in range(len(a)):
            model=linear_model.Lasso(alpha=a[i])
            model.fit(trainfeat, trainresp) 
            est=model.predict(valfeat)
            rss=np.linalg.norm(est-valresp)**2
            rss=rss/len(valresp)
            rss_array[i]=rss

        idx=np.argmin(rss_array)
        best_a=a[idx]
        best_rss=rss_array[idx]
        best_model=linear_model.Lasso(alpha=best_a)
        best_model.fit(trainfeat, trainresp)
        coef=best_model.coef_  
      
        return (rss_array, best_a, best_rss, coef)
