---
layout: mathjax
---

To support hydrostatic and non-hydrostatic modes, the equation set of Laprise (1992) is adopted:

Prognostic equations

$$
\begin{align}
  \dt{\vec{v}} & = - \frac{1}{\rho} \gradpi{p} - \pd{p}{\pi} \gradpi{\phi} - f^* \vec{k} \times \vec{v} \\
  \gamma \dt{w} & = g \left(\pd{p}{\pi} - 1\right) \\
  \dt{T} & = \frac{1}{C_p \rho} \dt{p} + \frac{Q}{C_p} \\
  \dt{\ln{p}} & = - \frac{C_p}{C_V} + \frac{Q}{C_V T}
\end{align}
$$

Diagnostic relations

$$
\begin{align}
	\divpi{\vec{v}} + \pd{\dot{\pi}}{\pi} & = 0 \\
	\pd{\phi}{\pi} & = - \frac{1}{\rho} \\
	p & = \rho R T
\end{align}
$$

Definitions:

$$
\begin{align}
	& \dt{} = \left(\pdt{}\right)_\pi + \vec{v} \cdot \gradpi + \dot{\pi} \pd{}{\pi} \\
	& D_3 = \divpi{\vec{v}} + \rho \left(\gradpi{\phi}\right) \cdot \pd{\vec{v}}{\pi} - \rho g \pd{w}{\pi}
\end{align}
$$