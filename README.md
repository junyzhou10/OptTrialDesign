# OptTrialDesign
This package is used to find the clinical trial design which maximizing the the financial revenue and taking regulatory concerns into considerations at the same time. Regulatory agencies and health technology assessment (HTA) bodies usually concern not only the statistical significance of the results, but also the practical importance, and reliability. In particular, we consider data maturity and clinical significance, which can be specified in the App to make sure the output optimal design will meet these requirements.

## Distribution assumptions
The algorithm adopts the time-to-event endpoint under the proportional hazard (PH) assumption, which is common in clinical trial designs.

The event time and loss to follow-up time follow the exponential distributions.

Both uniform and piecewise accrual rate are supported.

## Installation
To install the package: 
```r
devtools::install_github("junyzhou10/OptTrialDesign")
```
## Usage
Users specify the required study design parameters, such as alpha, power, allocation of arms, median survival time and loss to follow-up (dropout) rate of both arms; financial parameters, such as fixed cost per subject, per unit time, and revenue information; regulatory-related parameters, such as data maturity constraints, and clincial meaningful thresholds. 

There is a [shiny app](https://junyzhou.shinyapps.io/OptTrialDesign/) available for easy exploration of the package.

## Other
Please also refer to [here](https://jzhou.org/posts/odt-capm/) to see some interesting thoughts about the idea behind OptTrialDesign.
