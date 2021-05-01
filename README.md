# OptTrialDesign
This package is used to find the clinical trial design which maximizing the the financial revenue and taking regulatory concerns into considerations at the same time. Regulatory agencies and health technology assessment (HTA) bodies usually concern not only the statistical significance of the results, but also the practical importance, and reliability. In particular, we consider data maturity and clinical significance, which can be specified in the App to make sure the output optimal design will meet these requirements.

## Distribution assumptions
The algorithm adopts the time-to-event endpoint under the proportional hazard (PH) assumption, which is common in clinical trial designs.

The event time and loss to follow-up time follow the exponential distributions.

Both uniform and piecewise accrual rate are supported.

## Installation
To install the package: devtools::install_github("junyzhou10/OptTrialDesign")
