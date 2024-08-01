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

Where

$$
n \in \\{ 0,1,...,14 \\}
$$

These solutions represent equally spaced lines moving from the negative x direction with different starting y values. The program will also find 15 more lines starting from the negative y direction with equally spaced starting x values. It will then visually display the solution:

![Metric for `echo -ne "1+0.01*(x^2+y^2)\n0\n0\n1+0.01*(x^2+y^2)" | ./app`](images/img1.png)

Each line represents a geodesic starting from the edge of the screen (left and lower) and moving to the opposite edge. The program is set to display the coordinates from x=-10 to x=10 and from y=-10 to y=10.

The program is able to calculate these geodesics for arbitrary metrics, which are taken as input from stdin and parsed algebraically using a symbolic library called GiNaC. That said, it is recommended to keep the metric tensor symmetric (give the same entries for g01 and g10) as I am not sure how to properly handle the case where they differ and the resulting calculations could be erroneous. Different coordinate system metrics are not allowed, to display them in the program you must input their equivalent cartesian form which can be calculated by hand using [tensor transformation rules](https://en.wikipedia.org/wiki/Covariant_transformation#With_coordinates).

# More examples

Here are some cool / useful examples of geodesics you can display:

Flat space metric

$$
g_{ab}=\\
\begin{pmatrix}
1 & 0 \\
0 & 1
\end{pmatrix}
\\
$$

`echo -ne "1\n0\n0\n1" | ./app`

![Metric for `echo -ne "1\n0\n0\n1" | ./app`](images/img2.png)

Sinusoidal example 1

$$
g_{ab}=\\
\begin{pmatrix}
1+\frac{1}{x^2+y^2} & sin^2(x) \\
sin^2(x) & 1+\frac{1}{x^2+y^2}
\end{pmatrix}
\\
$$

`echo -ne "1+1/(x^2+y^2)\nsin(x)^2\nsin(x)^2\n1+1/(x^2+y^2)" | ./app`

![Metric for `echo -ne "1+1/(x^2+y^2)\nsin(x)^2\nsin(x)^2\n1+1/(x^2+y^2)" | ./app`](images/img3.png)

Sinusoidal example 2

$$
g_{ab}=\\
\begin{pmatrix}
1+sin^2(y) & sin(xy) \\
sin(xy) & 1+sin^2(x)
\end{pmatrix}
\\
$$

`echo -ne "1+sin(y)^2\nsin(x*y)\nsin(x*y)\n1+sin(x)^2" | ./app`

![Metric for `echo -ne "1+sin(y)^2\nsin(x*y)\nsin(x*y)\n1+sin(x)^2" | ./app`](images/img4.png)

Flipped / inverted shape

$$
g_{ab}=\\
\begin{pmatrix}
1+0.007(x^2-y^2) & 0.007xy \\
0.007xy & 1+0.007(x^2+y^2)
\end{pmatrix}
\\
$$

`echo -ne "1+0.007*(x^2-y^2)\n0.007*x*y\n0.007*x*y\n1+0.007*(x^2+y^2)" | ./app`

![Metric for `echo -ne "1+0.007*(x^2-y^2)\n0.007*x*y\n0.007*x*y\n1+0.007*(x^2+y^2)" | ./app`](images/img5.png)

Spatial Schwarzschild metric (r_s = 2.0)

$$
g_{ab}=\\
\begin{pmatrix}
\frac{y^2}{(x^2+y^2)(1-\frac{\sqrt{x^2+y^2}}{2.0})}+\frac{1}{1-\frac{2.0}{\sqrt{x^2+y^2}}} & \frac{xy}{(x^2+y^2)(\frac{\sqrt{x^2+y^2}}{2.0}-1)} \\
\frac{xy}{(x^2+y^2)(\frac{\sqrt{x^2+y^2}}{2.0}-1)} & \frac{x^2}{(x^2+y^2)(1-\frac{\sqrt{x^2+y^2}}{2.0})}+\frac{1}{1-\frac{2.0}{\sqrt{x^2+y^2}}}
\end{pmatrix}
\\
$$

`echo -ne "y^2/((x^2+y^2)*(1-(x^2+y^2)^(1/2)/2.0))+1/(1-2.0/(x^2+y^2)^(1/2))\nx*y/((x^2+y^2)*((x^2+y^2)^(1/2)/2.0-1))\nx*y/((x^2+y^2)*((x^2+y^2)^(1/2)/2.0-1))\nx^2/((x^2+y^2)*(1-(x^2+y^2)^(1/2)/2.0))+1/(1-2.0/(x^2+y^2)^(1/2))" | ./app`

![Metric for `echo -ne "y^2/((x^2+y^2)*(1-(x^2+y^2)^(1/2)/2.0))+1/(1-2.0/(x^2+y^2)^(1/2))\nx*y/((x^2+y^2)*((x^2+y^2)^(1/2)/2.0-1))\nx*y/((x^2+y^2)*((x^2+y^2)^(1/2)/2.0-1))\nx^2/((x^2+y^2)*(1-(x^2+y^2)^(1/2)/2.0))+1/(1-2.0/(x^2+y^2)^(1/2))" | ./app`](images/img6.png)

(I was really hoping to be able to calculate Schwarzschild orbits or at least lightlike paths with this but it did not match data from actual Schwarzschild metric calculations, likely because time needs to be fully treated as a coordinate for a full relativistic / gravitational interaction. I might code something for this later)

Spatial [FLRW metric](https://en.wikipedia.org/wiki/Friedmann%E2%80%93Lema%C3%AEtre%E2%80%93Robertson%E2%80%93Walker_metric) with k set to 0.1 (takes forever to calculate, errors out on some of the equations, and gives nondeterministic results - probably due to coordinate singularities)

$$
g_{ab}=\\
\begin{pmatrix}
\frac{ky^2-1}{k(x^2+y^2)-1} & \frac{kxy}{1-k(x^2+y^2)} \\
\frac{kxy}{1-k(x^2+y^2)} & \frac{kx^2-1}{k(x^2+y^2)-1}
\end{pmatrix}
\\
$$

`echo -e "(k*y^2-1)/(k*(x^2+y^2)-1)\nk*x*y/(1-k*(x^2+y^2))\nk*x*y/(1-k*(x^2+y^2))\n(k*x^2-1)/(k*(x^2+y^2)-1)" | sed "s/k/0.1/g" | ./app`

![Metric for `echo -e "(k*y^2-1)/(k*(x^2+y^2)-1)\nk*x*y/(1-k*(x^2+y^2))\nk*x*y/(1-k*(x^2+y^2))\n(k*x^2-1)/(k*(x^2+y^2)-1)" | sed "s/k/0.1/g" | ./app`](images/img7.png)

Not sure what to make of this one

Previous metric but with k set to -0.1

`echo -e "(k*y^2-1)/(k*(x^2+y^2)-1)\nk*x*y/(1-k*(x^2+y^2))\nk*x*y/(1-k*(x^2+y^2))\n(k*x^2-1)/(k*(x^2+y^2)-1)" | sed "s/k/(-0.1)/g" | ./app`

![Metric for `echo -e "(k*y^2-1)/(k*(x^2+y^2)-1)\nk*x*y/(1-k*(x^2+y^2))\nk*x*y/(1-k*(x^2+y^2))\n(k*x^2-1)/(k*(x^2+y^2)-1)" | sed "s/k/(-0.1)/g" | ./app`](images/img8.png)

# Internals
The program uses GTK's C++ port gtkmm version 4 to display the graphics in a window. It uses a C++ symbolic library called GiNaC to parse the entries for the metric tensor and algebraically manipulate them to calculate the Christoffel symbols and Christoffel symbol derivatives. It also uses GiNaC to bind these Christoffel symbols and Christoffel symbol derivatives to C++ functions (using GiNaC's [compile](https://www.ginac.de/tutorial/#Compiling-expressions-to-C-function-pointers) feature). Then it uses a C library called GSL to numerically solve the differential equations. In this example it uses GSL's implementation of the Runge-Kutta-Fehlberg method to solve the equations.

Let me know if you have any questions or run into any problems.

