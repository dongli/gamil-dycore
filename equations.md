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
  \dt{\ln{p}} & = - \frac{C_p}{C_V} D_3 + \frac{Q}{C_V T}
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

First we will introduce the standard hydrostatic profile:

$$
\begin{align}
  T(\mathbf{x}, \pi, t) & = \overline{T}(\pi) + T^\prime(\mathbf{x}, \pi, t) \\
  \phi(\mathbf{x}, \pi, t) & = \overline{\phi}(\pi) + \phi^\prime(\mathbf{x}, \pi, t) \\
  p(\mathbf{x}, \pi, t) & = \pi + p^\prime(\mathbf{x}, \pi, t)
\end{align}
$$

where $\mathbf{x}$ is the horizonal coordinates (e.g., $\lambda$, $\varphi$), and $\overline{T}$ and $\overline{\phi}$ are related by hydrostatic relation:

$$
\begin{equation}
	\pd{\overline{\phi}}{\pi} = - \frac{R \overline{T}}{\pi}
\end{equation}
$$

Prognostic equations

$$
\begin{align}
	\dt{\vec{v}} & = - \frac{1}{\rho} \gradpi{p^\prime} - \left(1 + \pd{p^\prime}{\pi}\right) \gradpi{\phi^\prime} - f^* \vec{k} \times \vec{v} \\
	\gamma \dt{w} & = g \pd{p^\prime}{\pi} \\
	\dt{T^\prime} & = \frac{1}{C_p} \left(\frac{R T}{p} - C_p \pd{\overline{T}}{\pi}\right) \dot{\pi} - \frac{1}{C_p \rho} \dt{p^\prime} \\
	\dt{p^\prime} & = \left(- \frac{C_p}{C_V} D_3 + \frac{Q}{C_V T}\right) p - \dot{\pi}
\end{align}
$$

Diagnostic relations

$$
\begin{align}
	\divpi{\vec{v}} + \pd{\dot{\pi}}{\pi} & = 0 \\
	\pd{\phi^\prime}{\pi} & = - \frac{1}{\rho} - \pd{\overline{\phi}}{\pi} \\
	p & = \rho R T
\end{align}
$$