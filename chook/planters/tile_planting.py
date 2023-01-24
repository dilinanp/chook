import scipy as sp

def generate_3D_problem(L, p2FP=0.05, p4FP=0.05):

    # Make a voxel problem on an LxLxL lattice with periodic BC.

    # 1. Define three types of voxel having "+" as a G.S. The three
    # types have 2, 4, or 6 frustrated facet plaquettes. There are 2
    # subtypes within the set of 2 and 4 FP voxels, and 1 within the
    # 6FP set.
    
    # 2. Define a probability distribution over the voxel classes,
    # p2FP, p4FP, and of course, p6FP = 1-p2FP-p4FP

    # L must be even to have periodic boundary under this partition
    # scheme
    
    # Current hardest regime seems to be to set p6FP = 1, i.e. 2,4=0.

    
    if L%2 !=0 or L < 4:
        raise Exception("L must be even and >= 4!")

    # For coding clarity
    M = L
    N = L
    num_nodes = L**3

    M_skip = 2
    N_skip = 2
    
    M_max = M
    N_max = N
    L_max = L
    
    J = sp.zeros((num_nodes, num_nodes))

    for l in range(L_max):
        M_offset = l%2
        N_offset = M_offset
        for n in range(N_offset, N_max, N_skip):
            for m in range(M_offset, M_max, M_skip):

                # Voxel vertices
                curr_node    = m + M*n + M*N*l
                nbr_node_M   = (m+1)%M + M*n + M*N*l
                nbr_node_N   = m + M*( (n+1)%N ) + M*N*l
                nbr_node_L   = m + M*n + M*N*((l+1)%L)
                nbr_node_MN  = (m+1)%M + M*( (n+1)%N ) + M*N*l
                nbr_node_ML  = (m+1)%M + M*n + M*N*((l+1)%L)
                nbr_node_NL  = m + M*((n+1)%N) + M*N*((l+1)%L)
                nbr_node_MNL = (m+1)%M + M*((n+1)%N) + M*N*((l+1)%L)

                voxel_vars = [curr_node, nbr_node_M, nbr_node_N, nbr_node_MN, nbr_node_L, nbr_node_ML, nbr_node_NL, nbr_node_MNL]

                # print voxel_vars
                J_vox = sample_voxel(p2FP, p4FP)
    
                # Put into problem
                J[ sp.ix_(voxel_vars, voxel_vars) ] = J_vox
    
    return J


def generate_2D_problem(L, p1, p2, p3):

    # p1,p2,p3: probabilities of generating a plaquette (on the
    # checkerboard sublattice) with 1,2,3 GSs, modulo Z2, and
    # including the FM. p4 = 1-(p1+p2+p3).
    
    if L%2 !=0 or L < 4:
        raise Exception("L must be even and >= 4!")

    # For coding clarity
    M = L
    N = L
    num_nodes = L**2

    M_skip = 2
    
    M_max = M
    N_max = N
    
    J = sp.zeros((num_nodes, num_nodes))
    h = sp.zeros((num_nodes,))

    for n in range(N_max):
        M_offset = n%2
        for m in range(M_offset, M_max, M_skip):

            # Get plaqutte vertices
            curr_node   = m + M*n
            nbr_node_M  = (m+1)%M + M*n
            nbr_node_N  = m + M*( (n+1)%N )
            nbr_node_MN = (m+1)%M + M*( (n+1)%N )

            plaquette_vars = [curr_node, nbr_node_M, nbr_node_N, nbr_node_MN]

            # print plaquetteVars
            J_plaq = sample_plaquette(p1, p2, p3)
            # FIX!! This does *not* break the degeneracy, but it makes
            # things more difficult for SA.
            #J_plaq = sp.rand()*J_plaq
             
            # Put into problem
            J[ sp.ix_(plaquette_vars, plaquette_vars)] = J_plaq
    
    return J


def sample_voxel(p2FP, p4FP):

    # This always allows for reflections via useInversion.

    # These settings replicate the data distribution sent to TAMU
    # group in Feb 2017. Note that with only C22 instances (pC21=0)
    # the GS degeneracy for class 2FP is *higher* (4GS) than class 4FP
    # for any pC41 setting (2GS for all.) That explains Wenlong's
    # results showing "26FP" being harder than "46FP."
    pC21 = 0.
    pC41 = 0.
    #pC41 = 1.
          
    R = sp.rand()

    if R < p2FP:
        # Make "root" voxel with 2 frustrated plaquettes
        J_vox_2FP = make_voxel_adj()
        if sp.rand() < pC21:
        #if sp.rand() < 0.5:
            # C2,1: A single FM bond. 12 distinct elements following
            # 48 transformations in Oh. One bond broken in GS; 1 GS
            # mod Z2.
            J_vox_2FP[0,1] = -1; J_vox_2FP[1,0] = -1;
        else:
            # C2,2: Opposite on same face FM bonds. 12 distinct
            # elements following 48 transformations in Oh. Two bonds
            # broken in GS; 4 GS mod Z2.
            J_vox_2FP[0,1] = -1; J_vox_2FP[1,0] = -1;
            J_vox_2FP[2,3] = -1; J_vox_2FP[3,2] = -1;
        J_vox = transform_voxel(J_vox_2FP)
    elif R < p2FP + p4FP:

        # 4 frustrated plaquettes
        J_vox_4FP = make_voxel_adj()
        if sp.rand() < pC41:
        #if sp.rand() < 0.5:
            # C4,1: Perpendicular/opposite FM bonds. 24 distinct
            # elements following 48 transformations in Oh. Two bonds
            # broken in GS; 2 GS mod Z2.
            J_vox_4FP[0,1] = -1; J_vox_4FP[1,0] = -1;
            J_vox_4FP[2,6] = -1; J_vox_4FP[6,2] = -1;
        else:
            # C4,2: Diagonally opposite FM bonds. 6 distinct elements
            # following 48 transformations in Oh. Two bonds broken in
            # GS; 2 GS mod Z2.
            J_vox_4FP[0,1] = -1; J_vox_4FP[1,0] = -1
            J_vox_4FP[6,7] = -1; J_vox_4FP[7,6] = -1
        
        J_vox = transform_voxel(J_vox_4FP)
    else:
        # 6 frustrated plaquettes
        J_vox_6FP = make_voxel_adj()
        # C6,1: The only element. 8 distinct elements following 48
        # transformations in Oh. Three bonds broken in GS; 8 GS mod
        # Z2.
        J_vox_6FP[0,1] = -1; J_vox_6FP[1,0] = -1
        J_vox_6FP[2,6] = -1; J_vox_6FP[6,2] = -1
        J_vox_6FP[5,7] = -1; J_vox_6FP[7,5] = -1
        
        J_vox = transform_voxel(J_vox_6FP)

    return J_vox


def transform_voxel(J_vox, use_rotation=True, use_inversion=True):

    # Apply a random element of the cube isometry group (Oh) to the J
    # specifying a unit cube (voxel.) This consists of selecting a
    # random rotation (from 24) followed by inversion (reflection
    # through the origin.) If both are active, this implements
    # reflections as well.
    
    # Assumes vertex indices translate to (m,n,l) with m fastest, i.e.
    # 0 -> [0,0,0], 1 -> [1,0,0], 2 -> [0,1,0], 3 -> [1,1,0], etc.

    # All 48 possible transformations appear uniformly, though
    # depending on J_vox there may be fewer distinct outcomes of
    # course.
    
    # This is to pad with zeros and ensure length 3
    format_str = "{0:0"+str(3)+"b}"

    U = sp.zeros((3,8),dtype=sp.int64)
    for i in range(8):
        u_str = format_str.format(i)
        U[:,i] = [int(u) for u in u_str]
    U = sp.flipud(U)

    # Translate to [-1,1] vertices
    V = 2*U-1

    if use_rotation:
        # Random rotation matrix
        R = get_random_rot()
    else:
        # No rotation
        R = sp.eye( len(V) )

    V_prime = R.dot(V)
    U_prime = (V_prime+1)/2

    # This is to apply the full octahedral transformation group
    # (i.e. including reflections!)  Inversion just "reflects through
    # the origin", but then translate back so corner (-1,-1,-1) is at
    # the origin.
    if use_inversion is True:
        if sp.rand() < 0.5:
            U_prime *= -1
            U_prime += sp.ones((3,1)).dot(sp.ones((1,8)))

    # Get the mapped variable indices, i.e. variable i maps to
    # corresponding
    I_prime = sp.array( [[1,2,4]]).dot(U_prime)

    I_prime = sp.int64( I_prime[0,:] )
    
    # Need to invert the map
    I_P_inv = sp.argsort(I_prime)
    J_vox_prime = J_vox[sp.ix_(I_P_inv,I_P_inv)]
        
    return J_vox_prime

        
def get_random_rot():

    # Create a random rotation matrix consisting of uniformly picking
    # one of 24 possible rotations of the cube in multiples of pi/2. 

    # Works as follows:

    # Imagine the facet orthogonal to the x axis, with x = 1.
    # 1. Select a random rotation of this facet (and the whole cube)
    # about the x axis

    # 2. Assign this facet randomly to one of the 6 of the cube. This
    # is done by properly deciding whether to do a y or a z rotation.
    # Rotation about the y axis maps the base facet into one of 4
    # facets (including the base/no rotation) while the z rotation by
    # pi/2, 3pi/2 maps it to the other two. Note that in this order,
    # we do *not* consider z rotations of 0 or pi as this will double
    # count. Also, to be uniform over rotations, choose a y rotation
    # w/Pr 4/6 and a z with 2/6.

    from scipy import cos, sin, pi
    
    kX = int( sp.rand()*4 )
    # This avoids double-counting of the 0 and pi cases.
    if sp.rand() < 4./6:
        kY = int( sp.rand()*4 )
        kZ = 0.
    else:
        kY = 0.
        kZ = 1+2*int( sp.rand()*2 )

    # # print kX,kY,kZ
    
    thetaX = kX*pi/2
    thetaY = kY*pi/2
    thetaZ = kZ*pi/2

    Rx = sp.array( [[1,0,0],[0,cos(thetaX),sin(thetaX)],[0,-sin(thetaX),cos(thetaX)]] )
    Ry = sp.array( [[cos(thetaY),0,-sin(thetaY)],[0,1,0],[sin(thetaY),0,cos(thetaY)]] )
    Rz = sp.array( [[cos(thetaZ),sin(thetaZ),0],[-sin(thetaZ),cos(thetaZ),0],[0,0,1]] )

    # Round is because they should all be integers anyway!
    return ( Rz.dot(Ry.dot(Rx)) ).round()

def make_voxel_adj():

    vox_adj_mat = sp.zeros((8,8),dtype=sp.int32)

    # Set manually
    
    vox_adj_mat[0,1] = 1
    vox_adj_mat[0,2] = 1
    vox_adj_mat[0,4] = 1
    vox_adj_mat[1,3] = 1
    vox_adj_mat[1,5] = 1
    vox_adj_mat[2,3] = 1
    vox_adj_mat[2,6] = 1
    vox_adj_mat[3,7] = 1
    vox_adj_mat[4,5] = 1
    vox_adj_mat[4,6] = 1
    vox_adj_mat[5,7] = 1
    vox_adj_mat[6,7] = 1

    vox_adj_mat = vox_adj_mat + vox_adj_mat.T
    return vox_adj_mat


def sample_plaquette(p1, p2, p3, normalize_GS=False):

    # Generates one of 4 types of plaquette depending on parameters
    # p1..3 (p4=1-p1-p2-p3)
    # Type i is a plaqutte containing i GSs (up to flip) including the
    # FM state. Achieved by generating i weak bonds in the cycle, one
    # of which is AFM.

    # Returned J magnitudes on cycle edges are in {1,2}. One of
    # the mag-1 edges is AFM; all other edges are FM.

    # Note that this (may) return an *integer-valued* matrix, but when put into full problem in generate_2D_problem, gets converted to real.

    R = sp.rand()
    if R < p1:
        num_GS = 1
    elif R < p1+p2:
        num_GS = 2
    elif R < p1+p2+p3:
        num_GS = 3
    else:
        num_GS = 4

    J_plaq = 2*make_plaquette_adj()
    edges = sp.where(sp.triu(J_plaq))
    num_edges = len(edges[0]) # i.e 4

    P = sp.random.permutation(num_edges)
    # Set the weak edges
    for e in P[0:num_GS]:
        (v1,v2) = edges[0][e],edges[1][e]
        J_plaq[v1,v2] *= 0.5
        J_plaq[v2,v1] *= 0.5
    # Arbitrarily select the first of the permutation to be AFM
    e_flip = P[0]
    (v1,v2) = edges[0][e_flip],edges[1][e_flip]
    J_plaq[v1,v2] *= -1
    J_plaq[v2,v1] *= -1

    # Scale such that each cycle has the same GS energy
    if normalize_GS:
        J_plaq = sp.float64(J_plaq)
        minus_EGS = sp.triu(J_plaq).sum()
        J_plaq /= minus_EGS

    return J_plaq


def make_plaquette_adj():

    plaquette_adj_mat = sp.zeros((4,4),dtype=sp.int32)

    plaquette_adj_mat[0,1] = 1
    plaquette_adj_mat[0,2] = 1
    plaquette_adj_mat[1,3] = 1
    plaquette_adj_mat[2,3] = 1

    plaquette_adj_mat = plaquette_adj_mat + plaquette_adj_mat.T
    return plaquette_adj_mat
    


