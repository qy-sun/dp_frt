# **Code for "Differentially Private Fisher Randomization Test with Binary Outcomes"**

---

This repository contains code for reproducing the results in the paper *"Differentially Private Fisher Randomization Test with Binary Outcomes"*. 

---

## **Description**

The repository is organized into the following main directories:

| Folder | Description |
|:--|:--|
| **DP-FRT.R** | Core implementation of the differentially private Fisher Randomization Test using Bayesian denoising approach. |
| **DP-Decision.R** | Implements Bayesian and frequentist decision-making rules built on the DP-FRT framework. |
| **Illustrative Figs/** | Scripts for generating illustrative figures **Figure 1** and **Figure 2**. |
| **DP-Simu/** | Code for DP simulation studies, including **Table 1**, **Figure 3**, and **Table 2**. |
| **Causal-Simu/** | Code for causal simulation studies, including **Table 3**, **Figure 4**, and **Figure 5**. |
| **Real/** | Real-data application to ADAPTABLE cardiovascular trial for producing **Figure 6** and **Table 4**. |

---

## **Additional Notes**

Some experiments were executed on the Duke Compute Cluster (DCC) using parallel jobs to accelerate replication. Minor numerical discrepancies may occur if one chooose to run the codes locally due to randomization and parallel execution order.
