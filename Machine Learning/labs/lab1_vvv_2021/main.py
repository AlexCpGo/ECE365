import numpy as np


class Lab1(object):
    def solver(self, A, b):
        ans = np.linalg.inv(A).dot(b)
        return ans

    def fitting(self, x, y):
        onestock = np.ones((x.size,1))
        x = x.reshape((-1,1))
        stock = np.hstack((x,onestock))
        coeff = np.linalg.pinv(stock).dot(y)
        return coeff

    def naive5(self, X, A, Y):
        # Calculate the matrix with $(i,j$)-th entry as  $\mathbf{x}_i^\top A \mathbf{y}_j$ by looping over the rows of $X,Y$.
        output = np.zeros((X.shape[0],Y.shape[0]),dtype = "int")
        for i in range(X.shape[0]):
            for j in range(Y.shape[0]):
                output[i,j] = np.dot(np.dot(X[i],A),Y[j])
        return output

    def matrix5(self, X, A, Y):
        # Repeat part (a), but using only matrix operations (no loops!).
        return np.dot(np.dot(X,A),Y.T)

    def naive6(self, X, A):
        # Calculate a vector with $i$-th component $\mathbf{x}_i^\top A \mathbf{x}_i$ by looping over the rows of $X$.
        output = np.zeros((X.shape[0]),dtype = "int")
        for i in range(X.shape[0]):
            output[i] = np.dot(np.dot(X[i],A),X[i])
        return output

    def matrix6(self, X, A):
        # Repeat part (a) using matrix operations (no loops!).
        return np.sum(np.multiply(X.dot(A),X), axis=1)
