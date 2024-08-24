
-----------------------------------------------------------
Replication package for

"An Analytical Framework for Public Debt Management"

By Jean-Paul Renne, University of Lausanne.

Disclaimer: The views expressed in this paper are those of the author only.

Corresponding author: Jean-Paul.renne@unil.ch
-----------------------------------------------------------

This package contains the codes and data allowing to replicate the above-mentioned paper.


------ How to use the codes?:

To generate all tables and figures of the paper, simply run "main.m" using R.
To re-estimate the model set "indic_estim" to 1.
By default, tables and figures associated with the calibration section are produced.
To produce tables and chart associated with the Demand/Supply section of the model use "indic_DemSup <- 1".
In order to generate the tables and charts associated with the issuance performance analysis, use "indic_run_performances <- 1".


------ Estimation data:

The dataset used in the paper is "Data/data.Rda". It is automatically loaded. Estimation data can be updated by setting "indic_load_data" to 1; data will then be collected from the Internet:
- CPI and real GDP are extracted from the FRED website.
- Nominal and real yields are extracted from the federal Reserve Board website.


------ Outputs:

Figures and tables can be found in folders of the same name.
