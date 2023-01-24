class lattice(dict):
    
    def __init__(self, Lx, Ly):
        
        # Get neighbors given a point in the lattice
        def get_nn(x, y):
                
            # Define module
            mod = lambda x,y: (x % Lx, y % Ly)
            
            # Get neighbors
            nn = []
            
            if 0<=(x-1)<Lx: nn.append(mod(x-1, y  )) 
            if 0<=(y-1)<Lx: nn.append(mod(x  , y-1)) 
            if 0<=(x+1)<Lx: nn.append(mod(x+1, y  )) 
            if 0<=(y+1)<Lx: nn.append(mod(x  , y+1)) 
            
            return nn
        
        self.Lx = Lx
        self.Ly = Ly    
        super().__init__({(x,y) : get_nn(x,y) for x in range(Lx) for y in range(Lx)})
    
    # Get a random loop in the lattice
    def get_loop(self, min_size=None, max_size=None, max_attempts=100000):
        
        def gen_loop():
    
            from random import randint, choice
            
            # Loop
            loop = [(randint(0, self.Lx-1), randint(0, self.Ly-1))]
            loop.append(choice(self[loop[-1]]))
            
            # Set of visited positions
            s = set(loop)
            
            while True:
                loop.append(choice([x for x in self[loop[-1]] if x != loop[-2]]))
                if loop[-1] in s: break
                s.add(loop[-1])
            
            return loop
        
        def prune_loop(loop):
            return loop[list(filter(lambda x: x[1] == loop[-1], enumerate(loop)))[0][0]:]
        
        ok = False
        for _ in range(max_attempts):
            loop = prune_loop(gen_loop())
            if (min_size == None or len(loop)-1 >= min_size) and (max_size == None or len(loop)-1 <= max_size): 
                ok = True
                break
                    
        if not ok: raise Exception('Too many attempts.')
                
        return loop
    
# Generate an FCL Ising problem
def gen_fcl_ising(lattice, num_loops, min_size_loop=None, max_size_loop=None, min_rho=None, max_rho=None, max_attempts=100000):
    
    def add_loop(ising):
        
        from random import randint, choice
        get_coupling = lambda i: tuple(sorted((loop[i], loop[i+1])))
        
        ok = False
        for _ in range(max_attempts):
            temp = dict(ising)
            loop = lattice.get_loop(min_size_loop, max_size_loop, max_attempts)
            
            for i in range(len(loop)-1):
                idx = get_coupling(i)
                if idx in temp: temp[idx] -= 1
                else:           temp[idx] = -1
                    
            temp[get_coupling(randint(0, len(loop)-2))] += 2
            
            limits = True
            for i in range(len(loop)-1):
                idx = get_coupling(i)
                if (min_rho != None and temp[idx] < min_rho) or (max_rho != None and temp[idx] > max_rho):
                    limits = False
                    break
                        
            if limits:
                ok = True
                break
            
        if not ok: raise Exception('Too many attempts.')
        
        return temp
    
    ising = {}
    for _ in range(num_loops): ising = add_loop(ising)
        
    # Remove zero terms
    ising = {x:y for x,y in ising.items() if y != 0}
        
    return ising

# Get Chimera Hamiltonian given Ising
def get_dcl_chimera(ising, lambda_par, Mx=16, My=16, KK=4):
    
    from more_itertools import flatten
    
    Lx = max(x for x,y in flatten(ising.keys())) + 1
    Ly = max(y for x,y in flatten(ising.keys())) + 1
    
    if min(x for x,y in flatten(ising.keys())) < 0 or min(y for x,y in flatten(ising.keys())) < 0:
        raise Exception("Error.")
            
    if Lx > Mx or Ly > My:
        raise Exception("Too large problem.")
            
    chimera = {}
    
    index = lambda x,y,k,l: 2*KK*My*x + 2*KK*y + KK*l + k
    coupling = lambda x,y: tuple(sorted((x,y)))
    
    # Generate intra-cell couplings
    for x in range(Lx):
        for y in range(Ly):
            for k1 in range(KK):
                for k2 in range(KK):
                    chimera[coupling(index(x,y,k1,0), index(x,y,k2,1))] = -1
    
    # Generate inter-cell couplings
    for w,J in ising.items():
        ax, ay = w[0]
        bx, by = w[1]
        
        for k in range(KK):
            if ax == bx:
                chimera[coupling(index(ax,ay,k,1), index(bx,by,k,1))] = J * lambda_par
            elif ay == by:
                chimera[coupling(index(ax,ay,k,0), index(bx,by,k,0))] = J * lambda_par
            else:
                raise Exception("Error.")
    
    return sorted([(*w, J) for w,J in chimera.items()])

