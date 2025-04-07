# Oxy_Cyclical_Regression

This project describes a Bayesian implementation of a cyclical regression model to describe seasonal variation in stable oxygen isotopic values ($\delta^{18}\text{O}$) derived from serially-sampled animal tooth enamel. By modeling isotopic profiles using linear regression with sinusoidal components, we can summarize a tooth's isotopic profile for comparison and analysis. Further, the Bayesian model's structure allows you to summarize the assemblage using the model parameter values to describe an averaged profile for an assemblage, which can be used as the basis for inter-assemblage comparisons.

## What are isotopic profiles?

Serial samples along a tooth’s growth axis creates a time series of stable isotopic samples that reflect the conditions an animal experienced while their tooth mineralized (Passey and Cerling 2002; Sealy and van der Merwe 1986; Balasse 2003). Unlike bone collagen, tooth enamel does not remodel throughout the animal’s life, providing a reliable snapshot of environmental conditions and diet during an animal’s early life. Variation in these isotopic systems thus allows researchers to explore seasonal variation in animal diet and environment. Variation in stable oxygen isotopic values ($\delta^{18}\text{O}$) reflect variation in evapo-transpiration in an environment; generally, in temperate environments, higher values indicate drier conditions or decreased rainfall, though $\delta^{18}\text{O}$ values can also vary with altitude, distance from the ocean, and changes in utilized water sources (Pederzani and Britton 2019). Interpreting seasonal variation in isotopic time series thus requires an understanding of what local factors are relevant for inducing variation in $\delta^{18}\text{O}$ values at the seasonal scale (for interpretation of a single sequence) and at the annual scale (for interpretation of multiple specimens).

## Cyclical Regression

Time series of isotopic signals are difficult to summarize and interpret directly. Early analyses qualitatively described the shape of these “isotopic profiles” to identify teeth with similar patterns and to investigate individual differences in the range of isotopic values (e.g., Henton, et al. 2010). More recently, researchers developed regression models to quantitatively describe individual tooth profiles, describing the cyclical variation in isotopic values—specifically $\delta^{18}\text{O}$ values—with trigonometric regression terms (Balasse, et al. 2012; Chazin, et al. 2019). Such regression models allow researchers to summarize different aspects of tooth profiles, such as amplitude, periodicity, and the relative location of extreme values in the model (Stolwijk, et al. 1999).

The cyclical regression model estimates the expected $\delta^{18}\text{O}$ value for a given distance from the root-enamel-junction (REJ), $x$. Actual observaitons can then vary from this expectation with normally-distributed error terms.

Equation 1 (error from expected value):
$$\delta^{18}\text{O}\_{\text{observed}} = \text{Normal}(\mu\_{\delta^{18}\text{O}}[x], \sigma\_{\delta^{18}\text{O}})$$

Equation 2 (cyclical regression formula):
$$\mu\_{\delta^{18}\text{O}}[x] = (\alpha + \beta\_{0} \* x) + \beta_{1} \* \text{cos}(2 \pi \frac{x}{\lambda}) + \beta_{2} \* \text{sin}(2 \pi \frac{x}{\lambda})$$

This equation includes several key parameters that are estimated by the model for each tooth (see also Balasse, et al. 2012). Two of the terms describe the "linear trend" of the sinusoidal regression curve, which is the line around which the sinusoidal curve varies:
- Intercept ($\alpha$): the $\delta^{18}\text{O}$ at "distance from REJ (mm)" = 0.
- Slope ($\beta_{0}$): the change in average $\delta^{18}\text{O}$ value for every change of 1 in "distance from REJ (mm)".

The other parameters describe the nature of the sinusoidal curve:
- Periodicity ($\lambda$): The length (in mm) that covers the full sequence of the sinusoidal curve (e.g., from peak-to-peak).
- Regression coefficients ($\beta_{1}$ and $\beta_{2}$): coefficients that help determine the curve's amplitude and relative position of extreme values within a periodic cycle (see Stolwijk, et al. 1999; Chazin, et al. 2019).
- Amplitude ($A$): deviation from the linear trend (in $\delta^{18}\text{O}$) at the extreme values (peaks/troughs of the curve). Amplitude is calculated from the regression coefficients by the formula $A = \sqrt{\beta_{1}^2 + \beta_{2}^2}$

## This Project

This project page currently holds the poster presented at the 2023 ICAZ Stable Isotopes in Zooarchaeology Working Group meeting (2023 SIZWG Berlin) and the code used to produce the poster's analyses and figures.

A brief vignette for applying the Bayesian cyclical regression model to enamel $\delta^{18}\text{O}$ data is also included as a separate .R file.

A description of the model is being planned for publication aimed at _Archaeometry_; the process of building that paper will also be included in this project page.
