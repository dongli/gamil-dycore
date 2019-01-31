---
layout: mathjax
---

$$
\begin{align}
  & \frac{1}{C_p \rho} \dot{\pi} - \dot{\pi} \d{\overline{T}}{\pi} \\
  = & \left(\frac{1}{C_p \rho} - \d{\overline{T}}{\pi}\right) \dot{\pi} \\
  = & \frac{1}{\pi} \left(\frac{\pi}{C_p \rho} - \d{\overline{T}}{\ln \pi}\right) \dot{\pi} \\
  = & \frac{1}{\pi} \left({\color{red}\frac{\pi}{p}} \frac{R T}{C_p} - \d{\overline{T}}{\ln \pi}\right) \dot{\pi} \\
  = & \frac{1}{\pi} \left({\color{red}\frac{\pi}{p}} \frac{R T}{C_p} - \frac{R \overline{T}}{C_p} + \frac{R \overline{T}}{C_p} - \d{\overline{T}}{\ln \pi}\right) \dot{\pi} \\
\end{align}
$$

Introduce static stability parameter $\overline{c}^2 = R \left(\frac{R \overline{T}}{C_p} - \d{\overline{T}}{\ln \pi}\right)$

$$
\begin{align}
  & \frac{1}{\pi} \left({\color{red}\frac{\pi}{p}} \frac{R T}{C_p} - \frac{R \overline{T}}{C_p} + \frac{R \overline{T}}{C_p} - \d{\overline{T}}{\ln \pi}\right) \dot{\pi} \\
  = & \frac{1}{\pi} \left({\color{red}\frac{\pi}{p}} \frac{R T}{C_p} - \frac{R \overline{T}}{C_p} + \frac{\overline{c}^2}{R}\right) \dot{\pi}
\end{align}
$$