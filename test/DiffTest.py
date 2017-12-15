import numdifftools as nd
import numpy as np

def func(x):
    return np.sum(x**2)

if __name__ == '__main__':
    x = np.array([1.0, 2.0, 3.0])
    print func(x)
    df = nd.Gradient(func)
    print df(x)
    df = nd.Hessian(func)
    print df(x)
