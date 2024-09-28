## Estimating parameters for epidemic models from viral load data collected when infections are rising at an exponential rate.

Here I present a summary of one of my PhD projects, which was to was to estimate the parameters of epidemic models given one had access to the the rate of exponential growth of infections during the beginning of an epidemic $r$ and data of people's viral loads. We demonstrate with a classic SIR model:
```math
\begin{equation}
\begin{aligned}
\dot{S} &= -\beta S I \\
\dot{I} &= \beta S I - \gamma I \\
\dot{R} &= \gamma I
\end{aligned}
\end{equation}
```
Let $\mathcal{A}(\tau, t)$ be the proportion of the infected population at time $t$ that has age-of-infection $\tau$. Let $i(t)$ be the disease incidence at time $t$, that is the number of people who became infected at time $t$, while $K(\tau)$ is the probability that someone infected $\tau$ units ago is still infected. Then,
```math
\begin{align*}
	\mathcal{A}(\tau, t) &= \frac{i(t - \tau) K(\tau)}{\int_0^t i(t-x)K(x) \mathrm{d}x}
\end{align*}
```
A cool thing about this expression is that during the beginning of an epidemic when disease incidence is growing exponentially with rate $r$, the distribution of the ages-of-infection among the infected quickly approaches a limiting distribution. In the case of the SIR model, this is:
```math
\begin{align*}
	\mathcal{A}(\tau) &= (\gamma + r) \mathrm{e}^{-(\gamma + r ) \tau} \\
	&= \rho \mathrm{e}^{-\rho \tau}
\end{align*}
```
where $\rho$ is the amount of people that an infected person would infect per unit time in a totally susceptible population.
To give a clarifying example, the proportion of the infected population during the exponential growth phase of an epidemic that is well-described by an SIR model that were infected 1-2 days ago would be $\int_1^2 \mathcal{A}(\tau)\mathrm{d}\tau$. 

How can we put this to use?

### Epidemiological model parameters and viral loads

A key question for policymakers when confronted with an epidemic is to determine the basic reproductive number $\mathbb{R}_0$. This is a measure of how infectious a disease is and can be used to estimate the total amount of people that would be infected if the disease was let to propogate assuming nobody changed their behaviour in response to the epidemic as well as the vaccination coverage needed to prevent it spreading in the first place. 

In the SIR model, $\mathbb{R}_0 = \frac{\beta}{\gamma} = \beta t_I$ where $t_I = \frac{1}{\gamma}$ is the expected time an infected person will spend infected, i.e. in the I compartment of the model. The difficulty is that $\mathbb{R}_0$ is hard to measure directly. Specifically $\mathbb{R}_0$ is a theoretical quantity that equals the amount of infections that would be caused by one infectious person in a totally susceptible population. It is difficult to directly measure and thus indirect measures are used. Among one is the $\mathbb{R}_0-r$ relation that connects the basic reproductive number to the exponential growth rate in disease incidence during the beginning of the epidemic. This is useful as $r$ is much easier to directly measure by, for example, fitting an exponential curve to the amount of reported cases. In the SIR model, the $\mathbb{R0}-r$ relation is
```math
\mathbb{R}_0 = \frac{1}{1 + rt_I}
```
However we are still stuck as even if we can measure $r$, we still need $t_I$. 
