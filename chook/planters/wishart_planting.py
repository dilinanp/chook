import scipy as sp
from scipy import linalg

def generate_problem(num_nodes, M=1, discretize_bonds=False, num_decimals_R=None):

    # Generate a K_num_nodes (complete graph) Ising weight matrix J with the
    # (+)^num_nodes state as a planted GS.

    # Diagonal of J is zero.

    # Hamiltonian is zero field, i.e:
    # E(s) = -0.5*s'*J*s

    # M specifies the number of columns in W (for M>=num_nodes, FM and
    # easy.)

    # num_decimals_R: number of decimal points to round the uncorrelated
    # Gaussian used to generate the w elements. This is to avoid
    # numerical issues where a spurious state takes over as the GS.

    # Alternatively, can even replace the Gaussian with a bounded
    # range uniform discrete distribution in [-range,+range]...

    # Plants the FM GS
    t = sp.ones( (num_nodes, 1) )
    
    # Sample correlated Gaussian with covariance matrix sigma
    # Note: rank(sigma) = num_nodes-1
    sigma = num_nodes/(num_nodes-1.)*sp.eye( num_nodes ) - 1./(num_nodes-1)*t.dot( t.T )
    sigma_sqrt = sp.sqrt((num_nodes-1.)/num_nodes)*sigma

    if discretize_bonds:
        R = sp.random.choice([-1,1], (num_nodes, M))
    else:
        R = sp.randn(num_nodes, M)

        if num_decimals_R is not None:
                R = R.round(decimals=num_decimals_R)


    W = sigma_sqrt.dot(R)
    J_tilde = -1./num_nodes*W.dot(W.T)
    J = J_tilde - sp.diag(sp.diag(J_tilde))

    if discretize_bonds:
        J *= num_nodes*num_nodes*(num_nodes-1)
        J = J.round().astype(int)

    return J
