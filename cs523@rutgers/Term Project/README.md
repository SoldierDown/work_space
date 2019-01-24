> *Don't forget to add [Mathjax](https://www.mathjax.org/) plugin for GitHub*
# Abstract
> Smoke simulation is one of the most important parts in physics-based simulation. In this project, our goal is to visualize 2D smoke  
> based on 2D Navier-Stokes equation. It is composed of the following steps:  
> - First of all, we apply numerical methods and solve Navier-Stokes equation by splitting the original equation into three major parts:  
>   - Advection: $\frac{Dq}{dt} = 0$  
>  
> Through advection, the properties of smoke move along with the velocity field  
>   - Diffusion:  $\frac{\partial\vec{u}}{\partial t} = \nu\nabla^2\vec{u}$  
>  
> Through diffusion, density tends to be uniform in the domain.  
>   - Projection: $\frac{\partial\vec{u}}{\partial t} + \frac{1}{\rho}\nabla p = 0$  
>  
> Projection is to make the velocity field divergence-free, which means the smoke is incompressible.  
> - Secondly, we design a GUI using QT. This GUI helps us interact with the smoke. We can add source to the domain by simply using mouse, 
> pause/restart the simulation and run the simulation frame by frame.
>  
> View my [technical report](./Term_Project.pdf)  
> View my [source code](./SmokeSimulation)
