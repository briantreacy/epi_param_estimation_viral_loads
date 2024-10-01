## Estimating epidemic parameters from viral load data

I present a summary of one of my PhD projects, which outlined a method to estimate the parameters of a class of epidemic models given one had access to the rate of exponential growth $r$ seen in the infected population during the beginning of an epidemic along with viral load data collected during this same period. We demonstrate with an SIR model:
```math
\begin{equation}
\begin{aligned}
\dot{S} &= -\beta S I \\
\dot{I} &= \beta S I - \gamma I \\
\dot{R} &= \gamma I
\end{aligned}
\end{equation}
```
Let $\mathcal{A}(\tau, t)$ be the proportion of the infected population at time $t$ that have age-of-infection $\tau$. The disease incidence $i(t)$ is the number of people who became infected at time $t$, while $K(\tau)$ is the probability that someone infected $\tau$ units ago is still infected now. Then,
```math
\begin{align*}
	\mathcal{A}(\tau, t) &= \frac{i(t - \tau) K(\tau)}{\int_0^t i(t-x)K(x) \mathrm{d}x}
\end{align*}
```
A cool thing about this expression is that it approaches a limiting distribution during the beginning of an epidemic when disease incidence is growing exponentially with rate $r$ and thus loses its dependence on the epidemic time $t$. For the SIR model, this limiting distribution is:
```math
\begin{align*}
	\mathcal{A}(\tau) &= (\gamma + r) \mathrm{e}^{-(\gamma + r ) \tau} \\
	&= \rho \mathrm{e}^{-\rho \tau}
\end{align*}
```
where $\rho$ is the amount of people that an infected person would infect per unit time in a totally susceptible population.
To give a clarifying example, the proportion of the infected population that were infected 1-2 days ago during the exponential growth phase of an epidemic would be $\int_1^2 \mathcal{A}(\tau)\mathrm{d}\tau$. 

How can we put this to use?

### Shortcomings of the $\mathbb{R}_0 - r$ method

A key challenge for modellers when confronted with an epidemic is to determine the basic reproductive number $\mathbb{R}_0$. This is a measure of infectiousness and can be used to estimate the total amount of people that would be infected if the disease was let to propogate assuming nobody changed their behaviour in response to the outbreak as well as the vaccination coverage needed to prevent the epidemic occurring. 

In the SIR model, $\mathbb{R}_0 = \frac{\beta}{\gamma} = \beta t_I$ where $t_I = \frac{1}{\gamma}$ is the expected time an infected person will spend infected, i.e. in the I compartment of the model. $\mathbb{R}_0$ is a theoretical quantity that equals the amount of infections that would be caused by one infectious person in a totally susceptible population. It is difficult to directly measure and thus indirect approaches are used. Among one is the $\mathbb{R}_0-r$ relation that connects the basic reproductive number to the exponential growth rate in disease incidence during the beginning of the epidemic. This is useful as $r$ is much easier to directly measure by, for example, fitting an exponential curve to the amount of reported cases. In the SIR model, the $\mathbb{R}_0-r$ relation is
```math
\mathbb{R}_0 = \frac{1}{1 + rt_I}
```
However we are still stuck as even if we can measure $r$, we still need $t_I$. Directly measuring the length of people's infectious times is difficult. It can be difficult to estimate exactly when someone was infected though still possible with contact tracing data. Exactly when one stopped being infectious is harder still. 

### Viral load data

Viral load data can be collected via RT-PCR tests where the Ct value is a measure of the viral load in a sample. The idea is to take a set of viral loads $\\{v_i\\}$ that were collected from individuals who were infected during the exponential growth phase of the epidemic and use them to infer the average age of infection $t_I$. To do this, we assume that a person's viral load is function of their age of infection $\nu(\tau)$. For the sake of simplicity, we made $\nu$ deterministic though in reality there is considerable variability between individuals. It is a piecewise function which is composed of two functions, $\nu_G$ and $\nu_D$, which describe the viral load during the growth and decline stage. We assume that upon infection individuals having an initial viral load $\nu_G(0)$ which grows monotonically to a maximum value $\nu_G(\tau^+) = \nu_D(\tau^+)$ whereupon it then decays monotonically towards 0:

![Viral load function $\nu(\tau)$](https://github.com/briantreacy/epi_param_estimation_viral_loads/blob/main/images/each_viral_load_has_two_ages_of_infection.png?raw=true)

The particular form I choose for $\nu_G$ and $\nu_D$ was exponential functions
```math
\nu(\tau) =
\begin{cases}
	\begin{aligned}
		\nu_G(\tau) &= \mathrm{e}^{\alpha \tau} && \text{if } \tau <= \tau^+ \\
		\nu_D(\tau) &= \nu^+ \mathrm{e}^{- \Phi(\tau - \tau^+)} && \text{if } \tau > \tau^+
	\end{aligned}
\end{cases}
```

Importantly, we can derive the inverse function $\nu^{-1}$ that returns the set of ages-of-infection that could have produced a given viral load:
```math
\nu^{-1}(v) =
\begin{cases}
	\{\tau^+ + \frac{1}{\Phi} \ln(\frac{\nu^+}{v})\} & \text{if } 0 < v < 1 \\
	\{\frac{1}{\alpha} \ln(v), \tau^+ + \frac{1}{\Phi} \ln(\frac{\nu^+}{v}) \} & \text{if } 1 \leq v \leq \nu^+ \\
\end{cases}
```
### Viral load distribution $\mathcal{V}$
Using $\mathcal{A}(\tau)$ and $\nu(\tau)$, we can derive the probability distribution $\mathcal{V}$ of the viral loads among the infected during the exponential growth phase of the epidemic:
```math
\begin{align*}
	1 &= \int_{0}^{\infty} \mathcal{A}(\tau) \mathrm{d} \tau \qquad \text{as $\mathcal{A}$ is a probability distribution} \\
	&= \int_{0}^{\tau^+} \mathcal{A}(\tau) \mathrm{d}\tau + \int_{\tau^+}^{\infty} \mathcal{A}(\tau) \mathrm{d} \tau
\end{align*}
```
where $\int_{0}^{\tau^+} \mathcal{A}(\tau) \mathrm{d}\tau$ and $\int_{\tau^+}^{\infty} \mathcal{A}(\tau) \mathrm{d} \tau$ are the probabilities that a person's viral load is in its growth and decay stage respectively. 

A person's age-of-infection during the viral growth stage is given by the inverse of $\nu_G$:
```math
\begin{equation*}
	\tau = \nu_G^{-1}(v) =  \frac{1}{\alpha} \ln(v)
\end{equation*}
```
Thus
```math
\begin{equation*}
	\mathrm{d}\tau = \frac{1}{\alpha} \frac{1}{v} \mathrm{d}v
\end{equation*} 
```
Changing variables: 
```math
\begin{align*}
	\int_{0}^{\tau^+} \mathcal{A}(\tau) \mathrm{d}\tau &=\int_{0}^{\tau^+} \rho \mathrm{e}^{-\rho \tau} \mathrm{d} \tau \\ 
	&= \int_{1}^{v^+} \frac{\beta}{\alpha} v^{- (\frac{\beta}{\alpha} + 1)} \mathrm{d} v \\
	&= 1 - \mathrm{e}^{-\rho \tau^+}
\end{align*}
```
where $\frac{\beta}{\alpha} v^{- (\frac{\beta}{\alpha} + 1)}$ is the probability that a person in their viral growth stage has viral load $v$. 

Likewise the inverse of $\nu_D$ gives a person's age-of-infection during the viral decay stage:
```math
\begin{align*}
	\tau = \nu_D^{-1}(v) = \tau^+ + \frac{1}{\Phi} \ln(\frac{\nu^+}{v})
\end{align*}
```
thus
```math
\begin{align*}
	\mathrm{d}\tau = - \frac{1}{\Phi} \frac{1}{v} \mathrm{d} v
\end{align*}
```
Changing of variables:
```math
\begin{align*}
	\int_{\tau^+}^{\infty} \mathcal{A}(\tau) \mathrm{d}\tau &= \int_{\tau^+}^{\infty} \rho \mathrm{e}^{-\rho \tau} \mathrm{d} \tau \\
	&= \int_{0}^{\nu^+} \frac{\rho}{\Phi} \mathrm{e}^{\rho \tau^+} (\nu^+) ^{\frac{\rho}{\Phi}} v^{\frac{\rho}{\Phi} - 1} \mathrm{d}v \\
	&= \mathrm{e}^{-\rho \tau+}
\end{align*}
```
where $\frac{\rho}{\Phi} \mathrm{e}^{\rho \tau^+} (\nu^+) ^{\frac{\rho}{\Phi}} v^{\frac{\rho}{\Phi} - 1}$ is the probability that a person in their viral decay stage has viral load $v$.  

Therefore our distribution $\mathcal{V}$ of the viral loads among those infected during the exponential rise of the epidemic is
```math
\mathcal{V}(v) =
\begin{cases}
	\frac{\rho}{\Phi} \mathrm{e}^{-\rho \tau^+} (\nu^+)^{- \frac{\rho}{\Phi}} v^{\frac{\rho}{\Phi} - 1} & \text{if } 0 < v < 1 \\
	\frac{\rho}{\alpha} v^{-(\frac{\rho}{\alpha} + 1)} + \frac{\rho}{\Phi} \mathrm{e}^{-\rho \tau^+} (\nu^+)^{-\frac{\rho}{\Phi}} v ^{\frac{\rho}{\Phi} - 1} & \text{if } 1 \leq v \leq v^+ 
\end{cases}
```

### Inferring $t_I$ from viral loads

Given a set of viral loads $[v_1, v_2, ..., v_{N_v} ]$, we can construct a maximum likelihood estimator (MLE) of $\rho$ with the following likelihood function
```math
\begin{align*}
	\mathcal{L}(\rho) = \prod_{i = 1}^{N_v} \mathcal{V}(v_i; \rho)
\end{align*}
```
We tested the MLE on virtual data as seen in figure \ref{fig:Boxplots of R0 estimates}. For a given value of $r$ and $t_I$, we generated viral loads by drawing ages-of-infection from $\mathcal{A}(\tau)$ and applying $\nu$ to them. We then worked backwards by applying our MLE to those viral loads along with $r$, which is assumed to be observed from case incidence data. The goal was to estimate the $\rho$ which is then used to estimate $\mathbb{R}_0$. 

![Viral load function $\nu(\tau)$](https://github.com/briantreacy/epi_param_estimation_viral_loads/blob/main/images/boxplot_R0_estimates_SIR.png?raw=true)
