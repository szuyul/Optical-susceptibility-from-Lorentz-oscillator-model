# Optical-susceptibility-from-Lorentz-oscillator-model
Computation of linear and nonlinear susceptibility based on Drude-Lorentz oscillator model

The plasma frequency (w_p), oscillation frequency (w_o), and damping coefficient (γ) determine the first order susceptibility (X), which can further imply permittivity (ε), refractive index (n), absorption (imaginary part of n), and reflectivity. For multiple (w_o, gamma) pairs, we use superposition of the models with partition fraction, f. 
We model the first order susceptibility of few common materials:
1. air (80% nitrogen + 20% oxygen)
2. glass (silicon dioxide)
3. ion (NaCl)
4. metal (copper)

Based on the 1st order susceptibility, we calculate the 2nd and 3rd order susceptibility of frictional materials which contribute to:
1. second harmonic generation (SHG), X(2w: w, w)
2. sum frequency generation (SFG), X(w3: w1, w2), a 2D map
3. third harmonic generation (THG), X(3w: w, w, w)
4. self phase modulation (SPM), nonlinear absorption (NA), the real and imaginary part of X(w: -w, w, w)

The Lorentz oscillator model may not be accurate, but provides great insight in the light-matter interaction and the corresponding optical properties of materials.
