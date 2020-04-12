
# Forecasting the impact of covid-19 interventions with SEIR model using Bayesian estimates of age-specific contacts in the UK

## 1. `contact.matrix` 
The code allows to estimate the contact matrix for every country in the **POLYMOD contact survey** (Mossong et al., 2008) using the model proposed by van de Kassteele et al. (2017). 
These countries include BE (Belgium), DE (Germany), FI (Finland), GB (Great Britain), IT (Italy), LU (Luxembourg), NL (the Netherlands) or PL (Poland).
The user can specify the age range and age bands over which they want to obtain the contact matrix. 

The **eurostat 2011 Census database** (https://ec.europa.eu/eurostat/web/population-and-housing-census/census-data/2011-census) is used as demographic data.

The smooth contact surface by one year age bands estimated with a CAR model for every country is:
![Belgium](https://github.com/melodiemonod/social.mixing.ncov/blob/master/contact.matrix/figures/c.smooth_BE.pdf) ![Germany](/contact.matrix/figures/c.smooth_DE.pdf) ![Finland](figures/contact.matrix/c.smooth_FI.pdf) ![Great Britain](figures/contact.matrix/c.smooth_GB.pdf) 
![Italy](figures/contact.matrix/c.smooth_IT.pdf) ![Luxembourg](figures/contact.matrix/c.smooth_LU.pdf) ![Netherlands](figures/contact.matrix/c.smooth_NL.pdf) ![Poland](figures/contact.matrix/c.smooth_PL.pdf) 
<p>
    <img src="/contact.matrix/figures/c.smooth_DE.pdf" width="220" height="240" />
</p>
## 2. `Epidemic model`
We use specific age contact patterns and the latest estimates of the epidemiological parameters to simulate the evolution of the COVID-19 outbreak in the UK under a SEIR model. We compare the impact of the epidemic under different non-pharmaceutical physical distancing interventions (schools closure, social distancing, and lockdown).

## References
- Mossong, J., Hens, N., Jit, M., Beutels, P., Auranen, K., Mikolajczyk, R., Massari, M., Salmaso, S., Tomba, G. S., Wallinga, J., Heijne, J., Sadkowska-Todys, M., Rosinska, M.and Edmunds, W. J. (2008), ‘Social contacts and mixing patterns relevant to the spread of infectious diseases’, PLOS Medicine 5(3), 1–1.
- van de Kassteele, J., van Eijkeren, J. and Wallinga, J. (2017), ‘Efficient estimation of age-specific social contact rates between men and women’, Ann. Appl. Stat. 11(1), 320–339.
- Contact-patterns repository from van de Kassteele, https://github.com/kassteele/Contact-patterns
