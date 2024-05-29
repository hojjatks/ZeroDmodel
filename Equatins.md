#  Governing equations:


The equilibrium equation is given by:
$m\ddot \delta+ c\dot \delta + k(\delta - v_{pl}t)+\sigma \mu=0$
where:

$\mu=\mu_*+a\ln(\dot \delta/v_*)+b\ln(v_*\theta/d_{rs})$
and 
$\dot \theta=1-\dot \delta\theta/d{rs}$.

By changing the variables ($u=\delta-v_{pl}t$) we have:

$\dot \delta=\dot u+v_{pl}$

Then we have:

$m\ddot u+ c (\dot u+v_{pl}) + ku+\sigma \mu=0$

$\mu=\mu_*+a\ln((\dot u+v_{pl})/v_*)+b\ln(v_*\theta/d_{rs})$

$\dot \theta=1- (\dot u+v_{pl})\theta/d{rs}$. 


Let's use $x_1,x_2,x_3$ for the state space of the dynamical system.

$x_1=u=\delta -v_{pl}t$

$x_2=\dot u=\dot \delta -v_{pl}$

$x_3=\theta$

Now, we can write the equations in terms of $\dot x=f(x)$:


$\dot x_1=x_2$

$\dot x_2=-(1/m)[ c (x_2 +v_{pl}) + kx_1+\sigma \mu]$

$\dot x_3= 1- (x_2+v_{pl})x_3/d{rs}$

where $\mu$ is given by:

$\mu=\mu_*+a\ln((x_2+v_{pl})/v_*)+b\ln(v_*x_3/d_{rs})$







