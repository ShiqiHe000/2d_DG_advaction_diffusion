# 2d_DG_advaction_diffusion

## Approximation of wave equation
The basic model of wave propgation is the wave equation. The wave equation is:

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial&space;^{2}p}{\partial&space;t^{2}}-c^{2}(p_{xx}&plus;p_{yy})=0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial&space;^{2}p}{\partial&space;t^{2}}-c^{2}(p_{xx}&plus;p_{yy})=0" title="\frac{\partial ^{2}p}{\partial t^{2}}-c^{2}(p_{xx}+p_{yy})=0" /></a>

The variable `p` might represent the acoustic pressure and `c` would be the sound speed. 

Inspite of solving this PDE, we rewrite it as a system of three first order equations. 

In order to convert the wave equation to a systerm of first order equations, let:

<a href="https://www.codecogs.com/eqnedit.php?latex=u_{t}&space;=&space;-p_{x},&space;u_{t}&space;=&space;-p_{y}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?u_{t}&space;=&space;-p_{x},&space;u_{t}&space;=&space;-p_{y}" title="u_{t} = -p_{x}, u_{t} = -p_{y}" /></a>

`u` and `v` correspond to the velocities in a fluid flow. Assuming the order of mixed partial derivatives does not matter, then

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial&space;^2&space;p}{\partial&space;t^2}&space;&plus;&space;c^2&space;((u_{x})_{t}&space;&plus;&space;(v_{y})_{t})&space;=&space;0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial&space;^2&space;p}{\partial&space;t^2}&space;&plus;&space;c^2&space;((u_{x})_{t}&space;&plus;&space;(v_{y})_{t})&space;=&space;0" title="\frac{\partial ^2 p}{\partial t^2} + c^2 ((u_{x})_{t} + (v_{y})_{t}) = 0" /></a>

With proper initial conditions,

<a href="https://www.codecogs.com/eqnedit.php?latex=p_{t}&space;&plus;&space;c^2&space;(u_{x}&plus;v_{y})&space;=&space;0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?p_{t}&space;&plus;&space;c^2&space;(u_{x}&plus;v_{y})&space;=&space;0" title="p_{t} + c^2 (u_{x}+v_{y}) = 0" /></a>

The the system of equations becomes a group of equations for the pressure and two velocities componenets:

<a href="https://www.codecogs.com/eqnedit.php?latex=\begin{bmatrix}&space;p\\&space;u\\&space;v&space;\end{bmatrix}_t&space;&plus;&space;\begin{bmatrix}&space;0&space;&&space;c^2&space;&&space;0&space;\\&space;1&space;&&space;0&space;&&space;0\\&space;0&space;&&space;0&space;&&space;0&space;\end{bmatrix}&space;\begin{bmatrix}&space;p\\&space;u\\&space;v&space;\end{bmatrix}_x&space;&plus;&space;\begin{bmatrix}&space;0&space;&&space;0&space;&&space;c^2&space;\\&space;0&space;&&space;0&space;&&space;0\\&space;1&space;&&space;0&space;&&space;0&space;\end{bmatrix}&space;\begin{bmatrix}&space;p\\&space;u\\&space;v&space;\end{bmatrix}_y&space;=0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\begin{bmatrix}&space;p\\&space;u\\&space;v&space;\end{bmatrix}_t&space;&plus;&space;\begin{bmatrix}&space;0&space;&&space;c^2&space;&&space;0&space;\\&space;1&space;&&space;0&space;&&space;0\\&space;0&space;&&space;0&space;&&space;0&space;\end{bmatrix}&space;\begin{bmatrix}&space;p\\&space;u\\&space;v&space;\end{bmatrix}_x&space;&plus;&space;\begin{bmatrix}&space;0&space;&&space;0&space;&&space;c^2&space;\\&space;0&space;&&space;0&space;&&space;0\\&space;1&space;&&space;0&space;&&space;0&space;\end{bmatrix}&space;\begin{bmatrix}&space;p\\&space;u\\&space;v&space;\end{bmatrix}_y&space;=0" title="\begin{bmatrix} p\\ u\\ v \end{bmatrix}_t + \begin{bmatrix} 0 & c^2 & 0 \\ 1 & 0 & 0\\ 0 & 0 & 0 \end{bmatrix} \begin{bmatrix} p\\ u\\ v \end{bmatrix}_x + \begin{bmatrix} 0 & 0 & c^2 \\ 0 & 0 & 0\\ 1 & 0 & 0 \end{bmatrix} \begin{bmatrix} p\\ u\\ v \end{bmatrix}_y =0" /></a>

Or:

<a href="https://www.codecogs.com/eqnedit.php?latex=\overrightarrow{q}_t&space;&plus;&space;B&space;\overrightarrow{q}_x&space;&plus;&space;C&space;\overrightarrow{q}_y=0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\overrightarrow{q}_t&space;&plus;&space;B&space;\overrightarrow{q}_x&space;&plus;&space;C&space;\overrightarrow{q}_y=0" title="\overrightarrow{q}_t + B \overrightarrow{q}_x + C \overrightarrow{q}_y=0" /></a>

## Numerical flux
We intent to use upwind flux. Which is the upwind side is determined by the sign of the wave speed. Positive wave speed(with respect to x direction) means that the boundary condition is forced on the left. Luckly, the systerm describes the wave equation couples three wave speeds, positive, negative, and zero, with respect to the direction vector <a href="https://www.codecogs.com/eqnedit.php?latex=n_x\widehat{x}&plus;n_y\widehat{y}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?n_x\widehat{x}&plus;n_y\widehat{y}" title="n_x\widehat{x}+n_y\widehat{y}" /></a>. By decoupling the wave components into right going, left going and stationary waves, the outgoing waves are approximated by the interior solution(upwind), and the incoming waves are specified from the external state(also upwind). The derivation of the boundary fluxes is known as **_Riemann problem_**

**Numerical flux** is computed from the internal and external states(with the designation determined relative to the normal at the boundary) as:
<a href="https://www.codecogs.com/eqnedit.php?latex=\overrightarrow{F^*}(Q^L,&space;Q^R;&space;\widehat{n})=A^&plus;Q^L&plus;A^-Q^R" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\overrightarrow{F^*}(Q^L,&space;Q^R;&space;\widehat{n})=A^&plus;Q^L&plus;A^-Q^R" title="\overrightarrow{F^*}(Q^L, Q^R; \widehat{n})=A^+Q^L+A^-Q^R" /></a>

One would be able to derive the numerical flux either by diagonlizing the coefficient matrix A, or by using **_Rankine-Hugonist Condition_**, which yields to:

<a href="https://www.codecogs.com/eqnedit.php?latex=\widehat{A}\overrightarrow{q^*}=\widehat{A}\left&space;\{&space;\left&space;\{&space;q&space;\right&space;\}&space;\right&space;\}&plus;\frac{1}{2}\left&space;|&space;\widehat{A}&space;\right&space;|\left&space;\|&space;q&space;\right&space;\|" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\widehat{A}\overrightarrow{q^*}=\widehat{A}\left&space;\{&space;\left&space;\{&space;q&space;\right&space;\}&space;\right&space;\}&plus;\frac{1}{2}\left&space;|&space;\widehat{A}&space;\right&space;|\left&space;\|&space;q&space;\right&space;\|" title="\widehat{A}\overrightarrow{q^*}=\widehat{A}\left \{ \left \{ q \right \} \right \}+\frac{1}{2}\left | \widehat{A} \right |\left \| q \right \|" /></a>

Where,

<a href="https://www.codecogs.com/eqnedit.php?latex=\widehat{A}=\widehat{n_x}A_x&plus;\widehat{n_y}A_y" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\widehat{A}=\widehat{n_x}A_x&plus;\widehat{n_y}A_y" title="\widehat{A}=\widehat{n_x}A_x+\widehat{n_y}A_y" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=\left&space;|&space;\widehat{A}&space;\right&space;|&space;=&space;S\left&space;|&space;\Lambda&space;\right&space;|S^{-1}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\left&space;|&space;\widehat{A}&space;\right&space;|&space;=&space;S\left&space;|&space;\Lambda&space;\right&space;|S^{-1}" title="\left | \widehat{A} \right | = S\left | \Lambda \right |S^{-1}" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=\left&space;\{&space;\left&space;\{&space;q&space;\right&space;\}&space;\right&space;\}&space;=&space;\frac{q^L&space;&plus;&space;q^R}{2}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\left&space;\{&space;\left&space;\{&space;q&space;\right&space;\}&space;\right&space;\}&space;=&space;\frac{q^L&space;&plus;&space;q^R}{2}" title="\left \{ \left \{ q \right \} \right \} = \frac{q^L + q^R}{2}" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=\left&space;\|&space;q&space;\right&space;\|&space;=&space;q^L\cdot&space;n^-&plus;q^R\cdot&space;n^&plus;&space;=&space;q^L-q^R" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\left&space;\|&space;q&space;\right&space;\|&space;=&space;q^L\cdot&space;n^-&plus;q^R\cdot&space;n^&plus;&space;=&space;q^L-q^R" title="\left \| q \right \| = q^L\cdot n^-+q^R\cdot n^+ = q^L-q^R" /></a>

## Performance
![2d_advection_error](pics/2d_Gaussian_wave.png)

## Test case: 1D wave equation
### Governing equation
<a href="https://www.codecogs.com/eqnedit.php?latex=p_{tt}&plus;c^2p_{xx}=0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?p_{tt}&plus;c^2p_{xx}=0" title="p_{tt}+c^2p_{xx}=0" /></a>

With <a href="https://www.codecogs.com/eqnedit.php?latex=u_t&space;=&space;-&space;p_x" target="_blank"><img src="https://latex.codecogs.com/gif.latex?u_t&space;=&space;-&space;p_x" title="u_t = - p_x" /></a>

The solution is:

<a href="https://www.codecogs.com/eqnedit.php?latex=p&space;=&space;ccos(x-ct)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?p&space;=&space;ccos(x-ct)" title="p = ccos(x-ct)" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=u&space;=&space;cos(x-ct)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?u&space;=&space;cos(x-ct)" title="u = cos(x-ct)" /></a>

We can write the system in conservation law form

<a href="https://www.codecogs.com/eqnedit.php?latex=\overrightarrow{q_t}&plus;\bigtriangledown&space;\cdot&space;\overrightarrow{F}=0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\overrightarrow{q_t}&plus;\bigtriangledown&space;\cdot&space;\overrightarrow{F}=0" title="\overrightarrow{q_t}+\bigtriangledown \cdot \overrightarrow{F}=0" /></a>

With proper initial conditions,

<a href="https://www.codecogs.com/eqnedit.php?latex=p_t&plus;c^2u_x=0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?p_t&plus;c^2u_x=0" title="p_t+c^2u_x=0" /></a>

Then we obtain the system of equations by grouping the equations for the pressure and velocity:

<a href="https://www.codecogs.com/eqnedit.php?latex=\begin{bmatrix}&space;p\\&space;u&space;\end{bmatrix}_t&space;&plus;&space;\begin{bmatrix}&space;0&space;&&space;c^2\\&space;1&space;&&space;0&space;\end{bmatrix}&space;\begin{bmatrix}&space;p\\&space;u&space;\end{bmatrix}_x&space;=0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\begin{bmatrix}&space;p\\&space;u&space;\end{bmatrix}_t&space;&plus;&space;\begin{bmatrix}&space;0&space;&&space;c^2\\&space;1&space;&&space;0&space;\end{bmatrix}&space;\begin{bmatrix}&space;p\\&space;u&space;\end{bmatrix}_x&space;=0" title="\begin{bmatrix} p\\ u \end{bmatrix}_t + \begin{bmatrix} 0 & c^2\\ 1 & 0 \end{bmatrix} \begin{bmatrix} p\\ u \end{bmatrix}_x =0" /></a>

### Boundary condition
Exact solution is imposed on both side of the boundaries.

### Numerical fluxes
For wave equation, the natural choice of flux is Lax-Friedrichs flux
<a href="https://www.codecogs.com/eqnedit.php?latex=(F)^{*}=c\frac{u^{-}&plus;u^{&plus;}}{2}&plus;\left&space;|&space;c&space;\right&space;|\frac{1-\alpha&space;}{2}(u^{-}\cdot&space;\widehat{n}^{-}&plus;u^{&plus;}&space;\cdot&space;\widehat{n}^{&plus;})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?(F)^{*}=c\frac{u^{-}&plus;u^{&plus;}}{2}&plus;\left&space;|&space;c&space;\right&space;|\frac{1-\alpha&space;}{2}(u^{-}\cdot&space;\widehat{n}^{-}&plus;u^{&plus;}&space;\cdot&space;\widehat{n}^{&plus;})" title="(F)^{*}=c\frac{u^{-}+u^{+}}{2}+\left | c \right |\frac{1-\alpha }{2}(u^{-}\cdot \widehat{n}^{-}+u^{+} \cdot \widehat{n}^{+})" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=0\leq&space;\alpha&space;\leq&space;1" target="_blank"><img src="https://latex.codecogs.com/gif.latex?0\leq&space;\alpha&space;\leq&space;1" title="0\leq \alpha \leq 1" /></a>

If <a href="https://www.codecogs.com/eqnedit.php?latex=\alpha&space;=1" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\alpha&space;=1" title="\alpha =1" /></a>, the numerical flux is the average of the two solutions, known as central flux. For <a href="https://www.codecogs.com/eqnedit.php?latex=\alpha&space;=0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\alpha&space;=0" title="\alpha =0" /></a>, we recover a flux which always takes information from where it is coming; that is, it is an upwind flux.

In this case the Lax-Friedrichs flux yelds to
<a href="https://www.codecogs.com/eqnedit.php?latex=(F)^{*}=\widehat{n}&space;\cdot&space;c(\frac{u^{-}&plus;u^{&plus;}}{2}&plus;\frac{1-\alpha&space;}{2}(u^{-}-u^{&plus;}))" target="_blank"><img src="https://latex.codecogs.com/gif.latex?(F)^{*}=\widehat{n}&space;\cdot&space;c(\frac{u^{-}&plus;u^{&plus;}}{2}&plus;\frac{1-\alpha&space;}{2}(u^{-}-u^{&plus;}))" title="(F)^{*}=\widehat{n} \cdot c(\frac{u^{-}+u^{+}}{2}+\frac{1-\alpha }{2}(u^{-}-u^{+}))" /></a>

### Test case Performance
One element, domain:[-1.0, 1.0]
![1d_advection_error](pics/1d_advection_error.png)

Two elements, domain:[0.0, 1.0]
![1d_two_elements](pics/1d_2elem.png)

## Documentation from the source code
[Source code documentation]( https://shiqihe000.github.io/2d_DG_advaction_diffusion/output/html/index.html)

