import numpy as np

class FFT():
    def fft(X, Y):
        dx = X[1]-X[0]
        N  = len(Y)
        K  = np.fft.fftfreq(N, dx)
        F  = np.fft.fft(Y)
        return K, F
    def ifft(K, F): 
        dk = K[1] - K[0]
        N  = len(F)
        X  = np.fft.fftfreq(N, dk)
        Y  = np.fft.ifft(F)
        return X, Y
