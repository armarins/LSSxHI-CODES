import time
import numpy as np
from copy import deepcopy as dcopy

c_light = 299792458/1e3 #km/s
nu_HI   = 1420.405751768 #MHz
STERADIAN = np.radians(1/60)**2 #arcmin to steradian
####################################
def find_active_z_support(W_func, z_min, z_max, n_scan=5000, rtol=1e-5):
    '''
    Find the subinterval [z_active_min, z_active_max] within [z_min, z_max]
    where |W_func(z)|>rtol*max(|W_func|)
    '''
    zscan = np.linspace(z_min, z_max, n_scan)
    Wvals = np.abs(np.nan_to_num(W_func(zscan), nan=0.0))
    Wmax = np.max(Wvals)
    if Wmax==0.0: return z_min,z_max
        
    idx = np.where(Wvals>rtol*Wmax)[0]
    if idx.size==0: return z_min,z_max

    dz = zscan[1]-zscan[0]
    za = max( z_min, zscan[idx[ 0]]-dz )
    zb = min( z_max, zscan[idx[-1]]+dz )    
    return float(za), float(zb)

def kernel_galaxy_nodes(z_nodes, W_overd_func, config_dict, camb_params=None, camb_results=None, include_rsd=False, f_interp=None):
    """
    Evaluate galaxy clustering kernels on quadrature nodes z_nodes:
        overdensity: W_g^(overd)(z) = b_g(z)*n_i(z)
        rsd:         W_g^(rsd)(z)   = f(z)*n_i(z)  
    """
    z_nodes = np.asarray(z_nodes, dtype=float)
    W_overd = np.nan_to_num(W_overd_func(z_nodes), nan=0.0)
    W_rsd   = None

    if include_rsd:
        import cross_functions_theory as cxft
        if f_interp is None:
            f_interp = cxft.get_growth_rate_interp(camb_params, camb_results)
        bg = cxft.get_bias_model(bias_model=config_dict['bias_model'],bias_value=config_dict['bias_value'],
                                 pars=config_dict,
                                 path_to_file=config_dict['path_to_bias_file'],
                                 z_range=z_nodes)
        ni   = np.where(bg!=0.0, W_overd/bg, 0.0)
        W_rsd = f_interp(z_nodes)*ni

    return W_overd, W_rsd
    
def kernel_hi_nodes(z_nodes, z_min, z_max, camb_params, camb_results, config_dict, include_rsd_HI=False, f_interp=None):
    import cross_functions_theory_4forecast as c4ft
    '''
    Evaluate (ith)kernel at z:  
        overdensity: WHI^(overd)(z) = bHI(z)*THI(z)/(zmax-zmin)
        rsd:         WHI^(rsd)(z)   =   f(z)*THI(z)/(zmax-zmin)        
    directly on quadrature nodes znodes inside [zmin,zmax]
    '''
    W_overd,W_rsd = None,None
    z_nodes  = np.asarray(z_nodes, dtype=float)
    w_tophat = 1/(z_max-z_min)  # (normalized) tophat windoe for i-th channel
    THI = c4ft.hi_brightness_temperature(camb_params=camb_params, camb_results=camb_results,
                                         z=z_nodes, 
                                         model=config_dict['omegaHI_model'], fit_OHI=config_dict['fit_omegaHI'],
                                         OHI=config_dict['omegaHI'], use_input=False)
    THI = np.asarray(THI, dtype=float).ravel() #flattens, make the array into a 1Darray. faster and cheaper than flatten

    if (not config_dict['fit_biasHI']) or (config_dict['biasHI'] is None):
        bHI = c4ft.omegaHI_biasHI(model=config_dict['biasHI_model'], z=z_nodes)['biasHI']
        bHI = np.asarray(bHI, dtype=float).ravel()
    else:
        bHI = float(config_dict['biasHI'])*np.ones_like(z_nodes)
    W_overd=bHI*THI*w_tophat
    if include_rsd_HI:
        if f_interp is None:
            import cross_functions_theory as cxft
            f_interp = cxft.get_growth_rate_interp(camb_params, camb_results)        
        W_rsd    = f_interp(z_nodes)*THI*w_tophat
    return W_overd,W_rsd

    
    
def gauss_legendre_z_grid(za,zb,nz):
    '''
    Return sorted nodes z_nodes and weights w(z) for integral over [za,zb]

    idea: int_za^zb dz f(z) ~ 0.5(zb-za) sum w(x)f(x)
    x  = points for Gauss-Legendre quadrature approximation
    nz = provides the degree, that is used to fix the number of sample points and weights

    tip: for delta chi ~35-45 Mpc and k<1/Mpc, as the ones for LOWZ-HI, 
    nz within [150,250] is more than enough to high precision
    For CMBlensing or Galaxy bins ng=10, nz within [400,800] instead.
    '''
    x,w = np.polynomial.legendre.leggauss(nz) #on [-1,1] already ascending
    znodes = 0.5*(zb-za)*x + 0.5*(za+zb)
    wz = 0.5*(zb-za)*w
    return znodes, wz

def build_hybrid_k_grid(chi_max, k_min=1e-5,k_split=1e-2,k_max=0.6,nklog =150,npp=10):
    '''
    it is better to use a linear spacing where there is oscillatory Bessel regime (k>1e-2/Mpc)
    for the other parts, tipically within k in [1e-5,1e-2]/Mpc, the logspacing is more than enough
    k_plit: define the regime where we split from log to linear regime.
    [kmin , kplit]: log regime
    [kplit,  kmax]: linear regime    

    spacing in linear one can be generate use a k-step of delta_k = pi/npp/chi_max; where npp is the number of point per half-period (osc)
    '''
    klog = np.geomspace(k_min,k_split,nklog,endpoint=False)
    dk   = np.pi/(npp*chi_max)
    nlin = int(np.ceil( (k_max-k_split)/dk ))+1
    klin = np.linspace(k_split,k_max,nlin)
    return np.concatenate([klog,klin])

def build_jl_table(ell, x_max, dx=0.02, include_rsd=False):
    """
    Evaluate jl(x) once on a fine 1D x-grid for all tracers at this ell
    """
    from scipy.special import spherical_jn
    ell    = int(ell)
    x_tab  = np.arange(0.0, x_max+2*dx, dx)
    jl_tab = spherical_jn(int(ell), x_tab)
    if not include_rsd:
        return x_tab, jl_tab, None
    #sphBESSEL for recursive relation
    c0 = (2*ell**2+2*ell-1)/((2*ell-1)*(2*ell+3))         #c_(l,0) = (2l^2+2l-1)/((2l-1)*(2l+3))
    cp = -((ell+1)*(ell+2)) / ((2*ell + 1) * (2*ell + 3)) #c_(l,2) = - ((l+1)(l+2))/((2l+1)*(2l+3))
    neg_jldd_tab = c0 * jl_tab + cp * spherical_jn(ell + 2, x_tab)
    if ell >= 2:
        cm = -(ell*(ell-1))/((2*ell-1)*(2*ell+1))  #c_(l,-2) = - (l(l-1))/((2l-1)*(2l+1))
        neg_jldd_tab += cm * spherical_jn(ell - 2, x_tab)
    return x_tab, jl_tab, neg_jldd_tab

    
def compute_delta_ell(ell, k_grid, chi_nodes, weighted_overd, weighted_rsd=None, jl_table=None):
    delta   = np.zeros_like(k_grid)
    chi_min = chi_nodes[0]
    chi_max = chi_nodes[-1]
    dchi    = max(chi_max - chi_min, 1.0)
    klow    = max(k_grid[0], (ell-6*(ell**(1./3.))-15.)/chi_max) if ell>20 else k_grid[0]#lower k
    khigh   = (ell+0.5)/chi_min + 20.0*np.pi/ dchi
    mask    = (k_grid>=klow) & (k_grid<=khigh)  
    
    neg_jldd_tab = None #neg_jldd_tab: neg=negative; jldd=d^2/(dxdx)(jl). -jl''
    x  = chi_nodes[:,None]*k_grid[mask][None,:]  #x = k*chi. shape: [nz,nk]
    if jl_table is not None:
        x_tab, jl_tab, neg_jldd_tab = jl_table
        jl = np.interp(x, x_tab, jl_tab)
    else:
        from scipy.special import spherical_jn
        jl = spherical_jn(int(ell), x)

    delta[mask] = np.sum(weighted_overd[:,mask]*jl, axis=0)
    
    if (weighted_rsd is not None) and (neg_jldd_tab is not None):
        neg_jldd = np.interp(x, x_tab, neg_jldd_tab)
        delta[mask] += np.sum(weighted_rsd[:,mask]*neg_jldd, axis=0)    
    return delta

def cl_from_deltas(k_grid,deltaA, deltaB, int_type='s'):
    integrand = (k_grid**2)*deltaA*deltaB    
    if int_type.lower() in ['simpson', 'simp','s']:
        from scipy.integrate import simpson
        return (2/np.pi)*simpson(integrand,x=k_grid)
    elif int_type.lower() in ['trapezoid', 'trapz','trap', 't']:
        return np.trapezoid(integrand,x=k_grid)*(2/np.pi)
        