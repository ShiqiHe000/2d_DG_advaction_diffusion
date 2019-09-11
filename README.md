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



## Documentation from the source code
[Source code documentation]( https://shiqihe000.github.io/2d_DG_advaction_diffusion/output/html/index.html)

