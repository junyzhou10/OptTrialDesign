# OptTrialDesign
This package is used to find the optimal clinical trial design which maximizing the the financial revenue, and at meantime, take several regulatory concerns into considerations, including, data maturity and clinical meaningful results. 

## Distribution assumptions
The algorithm adopts the time-to-event endpoint under the proportional hazard (PH) assumption, which is common in clinical trial designs.

The event time and loss to follow-up time are assumed from exponential distributions.

## Installation
To install the package: devtools::install_github("junyzhou10/OptTrialDesign")
