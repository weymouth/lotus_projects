# Scaling

The top layer of the ocean and lakes contains an energetic wave-affected thin zone. The propogation of waves along the surface induces oscillating currents at depths of less than fifty metres.

This work focuses on a zero mean oscillating freestream, $ U=U_m\sin(\phi) $, where $\phi=\omega t/T $ is the phase angle and $U_m$ is the velocity amplitude, and T is the time period. We use $U_m=0.4ms^{-1}, T=4s$ which are typical value for the induced current from a surface wave in shallow water (cite). We then define the oscillation frequency of $\omega=2\pi / T = 1.57$.

We define the viscous scales of the flow with a Reynold's number based on the oscillation amplitude, $Re_a=\frac{Re_{\delta_s}^2}{2}=a U_m/\nu = 102,000$, where the amplitude of the freestream motion $a=U_m/\omega=0.255m$. From here we scale all lengths by $a$, using the velocity amplitude $U_m$ to scale velocity and $T$ to scale time.

For the simulation we define a reduced frequency of $St = 2a f /U_m = 0.318$

The Reynold's number based on the Stoke's length ($Re_{\delta_s} = U_m \delta_s / \nu = 2\sqrt{Re_a}$), where $\nu=1\times 10^{-6}$ is the kinematic viscosity of the fluid, and the Stoke's length is $\delta_s=\sqrt{2\nu/\omega}=0.00112m$ (cite), which gives us $Re_{\delta_s} = 451$. The flow will oscillate with an amplitude of $226 \delta_s$ which will probably be significant for scaling.