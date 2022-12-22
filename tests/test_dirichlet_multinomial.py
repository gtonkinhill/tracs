# test llk gives the same answer as found using Sage math cloud.
from mtran.dirichlet_multinomial import find_dirichlet_priors
import numpy as np

def test_dirichlet_multinomial():
    # compare to the following code in R
    # d=matrix(c(1, 1, 0, 5, 3, 3, 5, 1, 2, 7, 0, 0, 4, 6, 0, 0, 23, 4, 2, 9, 19, 19, 33, 19, 17, 13, 6, 23, 29, 6, 17, 10, 5, 25, 5, 16, 14, 9, 19, 24, 73, 90, 53, 91, 57, 77, 89, 85, 67, 99, 96, 86, 85, 65, 86, 91, 73, 96, 71, 78), ncol=3, byrow = FALSE, dimnames = list(NULL, c('a','b','c')))
    # d=t(apply(d, 1, sort))
    # MGLM::MGLMfit(d, dist='DM')
    r_result = np.array([20.8156311152126,4.38181182238621,0.889048781117318])

    count = np.array([[1, 19, 73] ,[1, 19, 90] ,[0, 33, 53] ,[5, 19, 91] ,[3, 17, 57] ,[3, 13, 77] ,[5, 6, 89] ,[1, 23, 85] ,[2, 29, 67] ,[7, 6, 99] ,[0, 17, 96] ,[0, 10, 86] ,[4, 5, 85] ,[6, 25, 65] ,[0, 5, 86] ,[0, 16, 91] ,[23, 14, 73] ,[4, 9, 96] ,[2, 19, 71] ,[9, 24, 78]])
    
    alphas = find_dirichlet_priors(count, tol=1e-10, method='FP')
    assert np.max(alphas - r_result) < 1e-3

    alphas = find_dirichlet_priors(count, tol=1e-10, method='LOO')
    assert np.max(alphas - r_result) < 1e-3

    return