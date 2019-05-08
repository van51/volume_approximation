import numpy as np
import volesti


def random_polytope(n, d):
    a = np.random.randn(n, d)
    a = a / np.linalg.norm(a, ord=2, axis=1, keepdims=True)
    b = np.ones(n)
    return volesti.Polytope(a, b)
