---
layout: mathjax
---

The mass equation for barotropic model is

$$
\begin{equation}
  \pdt{\phi_d} = - \frac{1}{a \cos{\varphi}} \left(\pd{H U}{\lambda} + \pd{H V \cos{\varphi}}{\varphi}\right)
  \label{eq:mass}
\end{equation}
$$

The polar cap of the South Pole is

$$
\begin{align}
  A_\text{SP} & = \int_{0}^{2 \pi} d\lambda \int_{- \frac{\pi}{2}}^{-\frac{\pi}{2} + \frac{\Delta{\varphi}}{2}} a^2 \cos{\varphi} d\varphi \nonumber \\
  \approx & N_\lambda a^2 \Delta{\lambda} {\color{red}\frac{1}{2} \left(\cos\left(- \frac{\pi}{2} + \frac{\Delta{\varphi}}{2}\right) - \cos\left(- \frac{\pi}{2}\right)\right) \frac{\Delta{\varphi}}{2}} \nonumber \\
  & = N_{\lambda} a^2 \Delta{\lambda} \cos{\varphi_{1\frac{1}{2}}} \frac{1}{4} \Delta{\varphi}
  \label{eq:polar-cap-area}
\end{align}
$$

For the North Pole it is

$$
\begin{align}
  A_\text{NP} & = N_{\lambda} a^2 \Delta{\lambda} \cos{\varphi_{N_\varphi-\frac{1}{2}}} \frac{1}{4} \Delta{\varphi}
\end{align}
$$

# A-grid

Integrate (\ref{eq:mass}) at polar cap with the South Pole as example. For the right hand side

$$
\begin{align}
  & \frac{1}{A_\text{SP}} \int_{- \frac{\pi}{2}}^{- \frac{\pi}{2} + \frac{\Delta{\varphi}}{2}} \left[\int_{0}^{2 \pi} \left(\cancel{\frac{1}{a \cos{\varphi}} \pd{H U}{\lambda}} + \frac{1}{a \cos{\varphi}} \pd{H V \cos{\varphi}}{\varphi}\right) d\lambda\right] a^2 \cos{\varphi} d\varphi \nonumber \\
  = & \frac{1}{A_\text{SP}} \int_{0}^{2 \pi} \left[\int_{- \frac{\pi}{2}}^{- \frac{\pi}{2} + \frac{\Delta{\varphi}}{2}} \left(\pd{H V \cos{\varphi}}{\varphi}\right) d\varphi\right] a d\lambda \nonumber \\
  = & \frac{1}{A_\text{SP}} \int_{0}^{2 \pi} \left[\left.\left(H V \cos{\varphi}\right)\right|_{- \frac{\pi}{2} + \frac{\Delta{\varphi}}{2}} - \left.\left(H V \cos{\varphi}\right)\right|_{- \frac{\pi}{2}}\right] a d\lambda \nonumber \\
  = & \frac{1}{A_\text{SP}} \int_{0}^{2 \pi} \left.\left(H V \cos{\varphi}\right)\right|_{- \frac{\pi}{2} + \frac{\Delta{\varphi}}{2}} a d\lambda \nonumber \\
  = & \frac{1}{A_\text{SP}} a \Delta{\lambda} \sum_{i = 1}^{N_\lambda} \frac{1}{2} \left(\cancel{H_{i,1} V_{i,1} \cos{\varphi_1}} + H_{i,2} V_{i,2} \cos{\varphi_2}\right) \nonumber \\
  = & \frac{1}{A_\text{SP}} a \Delta{\lambda} \frac{\cos{\varphi_2}}{2} \sum_{i = 1}^{N_\lambda} H_{i,2} V_{i,2} \nonumber \\
  = & \frac{2}{N_\lambda a \Delta{\varphi}} \frac{\cos{\varphi_2}}{\cos{\varphi_1\frac{1}{2}}} \sum_{i = 1}^{N_\lambda} H_{i,2} V_{i,2}
  \label{eq:a-grid-mass-polar-integral}
\end{align}
$$

Then we get the tendencies for the Poles

$$
\begin{align}
  \left.\pdt{\phi_d}\right|_\text{SP} & = \frac{2}{N_\lambda a \Delta{\varphi}} \frac{\cos{\varphi_2}}{\cos{\varphi_1\frac{1}{2}}} \sum_{i = 1}^{N_\lambda} H_{i,2} V_{i,2} \\
  \left.\pdt{\phi_d}\right|_\text{NP} & = - \frac{2}{N_\lambda a \Delta{\varphi}} \frac{\cos{\varphi_{N_\varphi-1}}}{\cos{\varphi_{N_\varphi-\frac{1}{2}}}} \sum_{i = 1}^{N_\lambda} H_{i,N_\varphi-1} V_{i,N_\varphi-1}
\end{align}
$$

# C-grid

Integrate (\ref{eq:mass}) at polar cap with the South Pole as example. For the right hand side

$$
\begin{align}
  & \frac{1}{A_\text{SP}} \int_{- \frac{\pi}{2}}^{- \frac{\pi}{2} + \frac{\Delta{\varphi}}{2}} \left[\int_{0}^{2 \pi} \left(\cancel{\frac{1}{a \cos{\varphi}} \pd{H U}{\lambda}} + \frac{1}{a \cos{\varphi}} \pd{H V \cos{\varphi}}{\varphi}\right) d\lambda\right] a^2 \cos{\varphi} d\varphi \nonumber \\
  = & \frac{1}{A_\text{SP}} \int_{0}^{2 \pi} \left[\int_{- \frac{\pi}{2}}^{- \frac{\pi}{2} + \frac{\Delta{\varphi}}{2}} \left(\pd{H V \cos{\varphi}}{\varphi}\right) d\varphi\right] a d\lambda \nonumber \\
  = & \frac{1}{A_\text{SP}} \int_{0}^{2 \pi} \left[\left.\left(H V \cos{\varphi}\right)\right|_{- \frac{\pi}{2} + \frac{\Delta{\varphi}}{2}} - \left.\left(H V \cos{\varphi}\right)\right|_{- \frac{\pi}{2}}\right] a d\lambda \nonumber \\
  = & \frac{1}{A_\text{SP}} \int_{0}^{2 \pi} \left.\left(H V \cos{\varphi}\right)\right|_{- \frac{\pi}{2} + \frac{\Delta{\varphi}}{2}} a d\lambda \nonumber \\
  = & \frac{1}{A_\text{SP}} a \Delta{\lambda} \frac{\cos{\varphi_{1\frac{1}{2}}}}{2} \sum_{i = 1}^{N_\lambda} \left(H_{i,1} + H_{i,2}\right) V_{i,1\frac{1}{2}} \nonumber \\
  = & \frac{2}{N_\lambda a \Delta{\varphi}} \sum_{i = 1}^{N_\lambda} \left(H_{i,1} + H_{i,2}\right) V_{i,1\frac{1}{2}}
  \label{eq:c-grid-mass-polar-integral}
\end{align}
$$

The tendencies for the Poles are

$$
\begin{align}
  \left.\pdt{\phi_d}\right|_\text{SP} & = \frac{2}{N_\lambda a \Delta{\varphi}} \sum_{i = 1}^{N_\lambda} \left(H_{i,1} + H_{i,2}\right) V_{i,1\frac{1}{2}} \\
  \left.\pdt{\phi_d}\right|_\text{NP} & = - \frac{2}{N_\lambda a \Delta{\varphi}} \sum_{i = 1}^{N_\lambda} \left(H_{i,N_\varphi-1} + H_{i,N_\varphi}\right) V_{i,N_\varphi-\frac{1}{2}}
\end{align}
$$