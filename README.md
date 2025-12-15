# **Code for "Differentially Private Fisher Randomization Tests for Binary Outcomes"**

Available on [arXiv](https://arxiv.org/abs/2511.20884) and under review.

---

This repository contains code for reproducing the results in the paper *"Differentially Private Fisher Randomization Tests for Binary Outcomes"*. 

---

## **Description**

The repository is organized into the following main directories:

| Folder | Description |
|:--|:--|
| **DP-FRT.R** | Core implementation of the differentially private Fisher randomization tests, including DP-FRT-Bayes (recommended), DP-FRT-p, and DP-FRT-t. |
| **DP-Decision.R** | Implements Bayesian and frequentist decision-making rules built on DP-FRT-Bayes. |
| **DP-TopUp.R** | Implements Bayesian sequential decision under additional privacy budget. |
| **Illustrative Figs/** | Code for generating illustrative figures **Figures 1--2**. |
| **DP-Simu/** | Code for assessing DP-FRT-Bayes estimates of p-values, including **Table 1** and **Table 2**. |
| **Causal-Simu/** | Code for evaluating decision rules under DP-FRT-Bayes, including **Table 3** and **Figures 3--5**. |
| **Real/** | Analysis of ADAPTABLE trial for producing **Figure 6**. |
| **Supp/** | Code for additional simulation studies in the supplementary material, including **Tables S1--S2** and **Figures S1--S4**. |

---

## **Additional Notes**

Some experiments were executed on the Duke Compute Cluster (DCC) using parallel jobs to accelerate replication. Minor numerical discrepancies may occur if one chooose to run the codes locally due to randomization and parallel execution order.
