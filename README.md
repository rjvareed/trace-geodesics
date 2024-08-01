# Compiling and running
This program is meant to be run on Linux and was compiled with the GNU C++ compiler. The libraries it uses do have equivalent Windows versions and it will probably work with other compilers but I have not tested it in any other environment.

Make sure you have gtkmm4, ginac, and gsl installed
```
sudo apt install libgtkmm-4.0-dev libginac-dev ginac-tools libgsl-dev
```
Compile and run with
```
make
echo -ne "1+0.01*(x^2+y^2)\n0\n0\n1+0.01*(x^2+y^2)" | ./app
```

# Overview
This program will trace and display geodesics. It does this by solving differential equations and displaying the result in an XY plane. You give it a two dimensional metric tensor for Cartesian coordinates and it will generate the corresponding differential equations for geodesics. A simple example of this would be:

$$
g_{ab}=\\
\begin{pmatrix}
1+x^2+y^2 & 0 \\
0 & 1+x^2+y^2
\end{pmatrix}
\\
$$

Which gives the contravariant metric tensor:

$$
g^{ab}=\\
\begin{pmatrix}
\frac{1}{1+x^2+y^2} & 0 \\
0 & \frac{1}{1+x^2+y^2}
\end{pmatrix}
\\
$$

You can now measure distance with these constructs, in this case

$$
ds^2=g_{ab}dx^adx^b=dx^2(1+x^2+y^2)+dy^2(1+x^2+y^2)
$$

The Christoffel symbols are given by

$$
\Gamma^a_{\\ \\ bc}=\frac{1}{2}g^{ad}(\partial_b g_{cd}+\partial_c g_{bd}-\partial_d g_{bc})
$$

Which can be used to solve the geodesic equation. The geodesic equation is a 2nd order non-linear differential equation whose solution will answer the question: What is the shortest path between two given points? This is known as a geodesic. This equation can be derived using variational calculus to minimize the distance $\int ds$ and is given as the following system:

$$
\frac{\partial^2x^a}{\partial s^2}-\Gamma^a_{\\ \\ bc} \frac{\partial x^b}{\partial s}\frac{\partial x^c}{\partial s}=0
$$

Which in this case is:

$$
\frac{\partial^2 x}{\partial s^2}=\frac{x(\frac{\partial y}{\partial s})^2-x(\frac{\partial x}{\partial s})^2-2y \frac{\partial y}{\partial s} \frac{\partial x}{\partial s}}{1+x^2+y^2}
$$

$$
\frac{\partial^2 y}{\partial s^2}=\frac{y(\frac{\partial x}{\partial s})^2-y(\frac{\partial y}{\partial s})^2-2x \frac{\partial y}{\partial s} \frac{\partial x}{\partial s}}{1+x^2+y^2}
$$

The program will numerically solve these equations for $(x(s),y(s))$ for the following initial conditions:

$$
x(0)=-10.0
$$

$$
\frac{\partial x}{\partial s}|_{s=0}=2.0
$$

$$
y(0)=-10.0+20.0/15*n
$$

$$
\frac{\partial y}{\partial s}|_{s=0}=0.0
$$

For

$$
n \in \\{ 0,1,...,14 \\}
$$

These solutions represent equally spaced lines moving from the negative x direction with different starting y values. The program will also find 15 more lines starting from the negative y direction with equally spaced starting x values. It will then visually display the solution:
