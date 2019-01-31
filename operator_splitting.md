---
layout: mathjax
---

# Problems

$$
\begin{align}
  \pdt{\vec{F}} = - \mathcal{L}_1\left(\vec{F}\right) - \mathcal{L}_2\left(\vec{F}\right)
\end{align}
$$

Where $\mathcal{L}_1$ is fast adaptive processes (e.g. gravity waves) and $\mathcal{L}_2$ is slow evolution processes (e.g. advection).

$$
\begin{align}
  \pdtt{\vec{F}} & = - \mathcal{L}_1^\prime\left(\vec{F}\right) \pdt{F} - \mathcal{L}_2^\prime\left(\vec{F}\right) \pdt{\vec{F}} \\
  & = {\mathcal{L}_1}^\prime \left(\mathcal{L}_1 + \mathcal{L}_2\right) + {\mathcal{L}_2}^\prime \left(\mathcal{L}_1 + \mathcal{L}_2\right)
\end{align}
$$

$$
\begin{align}
  \vec{F}^{n+1} & = \vec{F}^n + \Delta{t} \left(\pdt{\vec{F}}\right)^n + \frac{\Delta{t}^2}{2} \left(\pdtt{\vec{F}}\right)^n + O\left(\Delta{t}^3\right) \\
  & = \vec{F}^n - \Delta{t} \mathcal{L}_1^n - \Delta{t} \mathcal{L}_2^n + \frac{\Delta{t}^2}{2} \left[\left({\mathcal{L}_1}^\prime \mathcal{L}_1\right)^n + \left({\mathcal{L}_1}^\prime \mathcal{L}_2\right)^n + \left({\mathcal{L}_2}^\prime \mathcal{L}_1\right)^n + \left({\mathcal{L}_2}^\prime \mathcal{L}_2\right)^n\right] + O\left(\Delta{t}^3\right)
\end{align}
$$

where $\mathcal{L}_1^n = \mathcal{L}_1 \left(\vec{F}^n\right)$, etc.

# Predicator-corrector splitting scheme

$$
\begin{align}
  \pdt{\vec{P}} & = - L_1 \vec{P} - L_2 \vec{F}^n \\
  \pdt{\vec{Q}} & = - \frac{1}{2} \left(L_2 \vec{Q} - L_2 \vec{F}^n\right)
\end{align}
$$

where $\vec{P}^n = \vec{F}^n$, $\vec{Q}^n = \vec{P}^{n+1}$.

The accuracy of this scheme can be proved as following:

$$
\begin{equation}
  \pdtt{\vec{P}} = - L_1^\prime(\vec{P}) \pdt{\vec{P}} = L_1^\prime L_1 + L_1^\prime L_2^n
\end{equation}
$$

by noting that $\left(L_2^n\right)^\prime = 0$.

$$
\begin{align}
  \vec{P}^{n+1} & = \vec{P}^n - \Delta{t} L_1^n - \Delta{t} L_2^n + \frac{\Delta{t}}{2} \left[\left(L_1^\prime L_1\right)^n + \left(L_1^\prime L_2^n\right)^n\right] + O\left(\Delta{t}^3\right)
\end{align}
$$