import numpy as np
import matplotlib.pyplot as plt


# https://physics.byu.edu/faculty/colton/docs/phy442-winter20/lecture-11-lorentz-oscillator-model.pdf
# https://empossible.net/wp-content/uploads/2018/03/Lecture-2-Lorentz-and-Drude-models.pdf
# https://www.fzu.cz/~dominecf/eps/index.html

def eV_2_w(eV):
    return 2*np.pi*eV*2.41e14 # w, unit: rad

def w_plasma(N, q):
    e = 1.60218e-19 # unit: C
    m_e = 9.109e-31 # mass of an electron, unit: kg
    epsilon_o = 8.854e-12 # vacuum permittivity, unit: F/m

    return np.sqrt( (N*((q*e)**2)) / (m_e*epsilon_o) )

def invL(w, w_o, gamma):
    A = w_o**2 - np.square(w)
    B = w*gamma
    real_part = A / (np.square(A) + np.square(B))
    imag_part = B / (np.square(A) + np.square(B))
    return real_part + 1j*imag_part

def susceptibility_X1(w, f, w_p, w_o, gamma):
    chi_1 = np.zeros_like(w)
    for i in range(len(f)):
        chi_1 = chi_1 + f[i]*np.square(w_p)*invL(w, w_o[i], gamma[i])

    return chi_1

def susceptibility_X2(w, a, N, q, w_o, gamma):
    e = 1.60218e-19 # unit: C
    m_e = 9.109e-31 # mass of an electron, unit: kg
    epsilon_o = 8.854e-12 # vacuum permittivity, unit: F/m

    _, ax = plt.subplots(2,3)
    w_p = w_plasma(N, q)
    chi_1 = susceptibility_X1(w, [1], w_p, [w_o], [gamma])
    f = w/(2*np.pi)
    dual_color_plot(ax[0,0], f, np.real(chi_1), np.imag(chi_1), '$\chi_{real}$', '$\chi_{imag}$')
    ax[0,0].set_title('$\chi^{(2)}$ material')
    ax[0,0].set_xticks([])

    # SHG, X2(2w: w, w)
    chi_2_SHG = (N*a*((q*e)**3)/(epsilon_o*(m_e**2))) * invL(2*w, w_o, gamma)*np.power(invL(w, w_o, gamma),2)
    dual_color_plot(ax[1,0], f, np.real(chi_2_SHG), np.imag(chi_2_SHG), 'SHG $\chi^{(2)}_{real}$', 'SHG $\chi^{(2)}_{imag}$')
    ax[1,0].set_xlabel(' f (Hz) ')

    # wavelength plot
    wavelength = np.divide(c, f) * 1e9
    dual_color_plot(ax[0,1], wavelength, np.real(chi_1), np.imag(chi_1), '$\chi_{real}$', '$\chi_{imag}$')
    ax[0,1].set_title('$\chi^{(2)}$ material')
    ax[0,1].set_xticks([])
    dual_color_plot(ax[1,1], wavelength, np.real(chi_2_SHG), np.imag(chi_2_SHG), 'SHG $\chi^{(2)}_{real}$', 'SHG $\chi^{(2)}_{imag}$')
    ax[1,1].set_xlabel(' $\lambda$ (nm) ')

    # SFG, X2(w3: w1, w2)
    w1, w2 = np.meshgrid(w, w)
    w3 = w1 + w2
    chi_2_2D = (N*a*((q*e)**3)/(epsilon_o*(m_e**2))) * invL(w3, w_o, gamma)*invL(w1, w_o, gamma)*invL(w2, w_o, gamma)
    ax[0,2].set_title('$SFG_{real}$')
    ax[0,2].imshow(np.real(chi_2_2D), extent=[f[0],f[-1],f[-1],f[0]])
    ax[0,2].set_ylabel(' f (Hz) ')
    ax[1,2].set_title('$SFG_{imag}$')
    ax[1,2].imshow(np.imag(chi_2_2D), extent=[f[0],f[-1],f[-1],f[0]])
    ax[1,2].set_xlabel(' f (Hz) ')
    ax[1,2].set_ylabel(' f (Hz) ')

    return chi_2_SHG

def susceptibility_X3(w, b, N, q, w_o, gamma):
    e = 1.60218e-19 # unit: C
    m_e = 9.109e-31 # mass of an electron, unit: kg
    epsilon_o = 8.854e-12 # vacuum permittivity, unit: F/m

    _, ax = plt.subplots(3,2)
    
    w_p = w_plasma(N, q)
    chi_1 = susceptibility_X1(w, [1], w_p, [w_o], [gamma])
    f = w/(2*np.pi)
    dual_color_plot(ax[0,0], f, np.real(chi_1), np.imag(chi_1), '$\chi^{(1)}_{real}$', '$\chi^{(1)}_{imag}$')
    ax[0,0].set_title('$\chi^{(3)}$ material')
    ax[0,0].set_xticks([])

    # SPM and nonlinear absorption, X3(w: -w, w, w)
    chi_3_SPM = (N*b*((q*e)**4)/(epsilon_o*(m_e**3))) * invL(-w, w_o, gamma)*np.power(invL(w, w_o, gamma), 3)
    dual_color_plot(ax[1,0], f, np.real(chi_3_SPM), np.imag(chi_3_SPM), 'SPM $\chi^{(3)}_{real}$', 'NA $\chi^{(3)}_{imag}$')
    ax[1,0].set_xticks([])

    # THG, X3(3w: w, w, w)
    chi_3_THG = (N*b*((q*e)**4)/(epsilon_o*(m_e**3))) * invL(3*w, w_o, gamma)*np.power(invL(w, w_o, gamma), 3)
    dual_color_plot(ax[2,0], f, np.real(chi_3_THG), np.imag(chi_3_THG), 'THG $\chi^{(3)}_{real}$', 'THG $\chi^{(3)}_{imag}$')
    ax[2,0].set_xlabel(' f (Hz) ')


    # wavelength plot
    wavelength = np.divide(c, f) * 1e9
    dual_color_plot(ax[0,1], wavelength, np.real(chi_1), np.imag(chi_1), '$\chi^{(1)}_{real}$', '$\chi^{(1)}_{imag}$')
    ax[0,1].set_title('$\chi^{(3)}$ material')
    ax[0,1].set_xticks([])
    dual_color_plot(ax[1,1], wavelength, np.real(chi_3_SPM), np.imag(chi_3_SPM), 'SPM $\chi^{(3)}_{real}$', 'NA $\chi^{(3)}_{imag}$')
    ax[1,1].set_xticks([])
    dual_color_plot(ax[2,1], wavelength, np.real(chi_3_THG), np.imag(chi_3_THG), 'THG $\chi^{(3)}_{real}$', 'THG $\chi^{(3)}_{imag}$')
    ax[2,1].set_xlabel(' $\lambda$ (nm) ')

    return chi_3_SPM, chi_3_THG

def permittivity(chi):
    epsilon_r = 1 + chi # relative permittivity, unit: none
    epsilon_o = 8.854e-12 # vacuum permittivity, unit: F/m
    epsilon = epsilon_o*epsilon_r

    return epsilon_r, epsilon

def refractive_index(epsilon_r):
    n_complex = np.sqrt(epsilon_r)
    n = np.real(n_complex)
    kappa = np.imag(n_complex)
    return n, kappa, n_complex

def absorption(w, kappa):
    c = 3e8 # unit: m/s
    return 2*w*kappa/c

def reflectivity(n_1, n_2):
    r = (n_2 - n_1) / (n_2 + n_1)
    R = np.square(np.abs(r))
    return R

def conductivity(w, gamma, w_p):
    epsilon_o = 8.854e-12 # vacuum permittivity, unit: F/m
    tau = 1/gamma
    sigma_o = epsilon_o*(w_p**2)*tau
    return sigma_o/(1+np.square(w)*(tau**2))

def light_speed(epsilon):
    mu_o = 1.257e-6 # vacuum permeability, unit: N/A^2
    return 1/np.sqrt(mu_o*epsilon)

def characteristics_plot(w, f, w_p, w_o, gamma, material):
    chi = susceptibility_X1(w, f, w_p, w_o, gamma)
    epsilon_r, _ = permittivity(chi)
    n, kappa, n_complex = refractive_index(epsilon_r)
    alpha = absorption(w, kappa)
    R = reflectivity(1, n_complex)

    # frequency plot
    _, ax = plt.subplots(3,2)
    f = w/(2*np.pi)
    dual_color_plot(ax[0,0], f, np.real(chi), np.imag(chi), '$\chi_{real}$', '$\chi_{imag}$')
    ax[0,0].set_xticks([])
    ax[0,0].set_title(material)
    dual_color_plot(ax[1,0], f, n, kappa, 'n', '$\kappa$')
    ax[1,0].set_xticks([])
    dual_color_plot(ax[2,0], f, alpha, R, r'$\alpha$ (m$^{-1}$)', 'Reflectance')
    ax[2,0].set_xlabel(' f (Hz) ')

    # wavelength plot
    wavelength = np.divide(c, f) * 1e9
    dual_color_plot(ax[0,1], wavelength, np.real(chi), np.imag(chi), '$\chi_{real}$', '$\chi_{imag}$')
    ax[0,1].set_xticks([])
    ax[0,1].set_title(material)
    dual_color_plot(ax[1,1], wavelength, n, kappa, 'n', '$\kappa$')
    ax[1,1].set_xticks([])
    dual_color_plot(ax[2,1], wavelength, alpha, R, r'$\alpha$ (m$^{-1}$)', 'Reflectance')
    ax[2,1].set_xlabel(' $\lambda$ (nm) ')

def dual_color_plot(ax, x, y1, y2, y1_label, y2_label):
    color = 'tab:red'
    lns1 = ax.plot(x, y1, color=color, label=y1_label)
    ax.tick_params(axis ='y', labelcolor = color)

    color = 'tab:green'
    ax_twin = ax.twinx()
    lns2 = ax_twin.plot(x, y2, color=color, label=y2_label)
    ax_twin.tick_params(axis ='y', labelcolor = color)
    
    lns = lns1+lns2
    labs = [l.get_label() for l in lns]
    ax_twin.legend(lns, labs, loc=1)



#%% spectral config.
c = 3e8 # unit: m/s
wavelength = np.linspace(100, 800, 701)*1e-9 # unit: m
f = np.divide(c, wavelength)
w = 2*np.pi*f


#%% gas: 
# http://jungwirth.uochb.cas.cz/assets/papers/paper101.pdf
f = [0.2, 0.8] # molecules fraction
w_p = 2*np.pi*9e11 # f~30 cm^-1
w_o = [3.7e11, 2e11]
gamma = [2.86e10, 1.1e10]

characteristics_plot(w, f, w_p, w_o, gamma, 'N$_2$ + O$_2$ (air)')


#%% solid: UV resonances (oscillation of orbital electrons)
# assuming a single resonance and no damping at all
d = 2.7e-10 # unit: m
N = 1/(d**3) # number of atoms, unit: m^-3
q = 1
f = [1]
w_p = w_plasma(N, q)
w_o = [1.26e16]
gamma = [1e15]

characteristics_plot(w, f, w_p, w_o, gamma, 'SiO$_2$ (dielectric)')


#%% solution: infrared resonances (oscillation of ions)
# The masses are much larger, so the resonant frequencies are much smaller.
# 1M NaCl solution conductivity 8 S/m
f = [1]
w_p = 2*np.pi*2.62e14
w_o = [6.2e13]
gamma = [1e16]

characteristics_plot(w, f, w_p, w_o, gamma, 'NaCl (ionic)')


#%% solid: metals (no oscillation of valence free electrons)
# Copper conductivity 5.96e7 S/m
f = [1]
w_p = 2*np.pi*4.2e15
w_o = [0]
gamma = [1.05e14]
sigma_o = conductivity(0, gamma[0], w_p)

characteristics_plot(w, f, w_p, w_o, gamma, 'Copper (metal)')


#%% nonlinear second-order susceptibility, second harmonic generation (SHG), sum frequency generation (SFG)
# BaBiO3 (BBO)
wavelength = np.linspace(100, 800, 701)*1e-9 # unit: m
f = np.divide(c, wavelength)
w = 2*np.pi*f

d = 2.7e-10 # unit: m
N = 1/(d**3) # number of atoms, unit: m^-3
q = 1
w_p = w_plasma(N, q)
w_o = 1.26e16
gamma = 1e15

r0 = 3e-10
a = (w_o**2)/r0
susceptibility_X2(w, a, N, q, w_o, gamma)


#%% third harmonic generation (THG), self-phase modulation (SPM), and nonlinear absorption (NA)
# https://sci-hub.se/https://www.sciencedirect.com/science/article/abs/pii/B9781782424642000105?via%3Dihub
# hypothetical material
wavelength = np.linspace(100, 1000, 901)*1e-9 # unit: m
f = np.divide(c, wavelength)
w = 2*np.pi*f

N = 1e28 # number of atoms, unit: m^-3
q = 1
w_p = w_plasma(N, q)
w_o = 1e16
gamma = 4.5e13

r0 = 3e-10
b = (w_o**2)/(r0**2)
susceptibility_X3(w, b, N, q, w_o, gamma)

i = 1