![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)

# CacheTest

This project is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The software and data in this repository are a snapshot of the software and data that were used in the research reported on in the paper  [Fast multinomial logistic regression with group sparsity](https://doi.org/10.1287/ijoc.2024.0796) by Sheng Fu, Shixiang Li, Kai Yu, Piao Chen and  Zhisheng Ye. 

## Cite

To cite the contents of this repository, please cite both the paper and this repo, using the following DOIs.

[https://doi.org/10.1287/ijoc.2024.0796](https://doi.org/10.1287/ijoc.2024.0796)

[https://doi.org/10.1287/ijoc.2024.0796.cd](https://doi.org/10.1287/ijoc.2024.0796.cd)

Below is the BibTex for citing this version of the code.
```latex
@misc{fu2025fast,
  author =    {Sheng Fu, Shixiang Li, Kai Yu, Piao Chen, and Zhisheng Ye},
  publisher =     {INFORMS Journal on Computing},
  title =         {Fast multinomial logistic regression with group sparsity},
  year =          {2025},
  doi =           {10.1287/ijoc.2024.0796},
  url =           {https://github.com/INFORMSJoC/2024.0796},
  note =           {Available for download at https://github.com/INFORMSJoC/2024.0796},
}  
```
## Description

This directory contains the code for the *group simplex-based multinomial logistic regression (GSMLR)* algorithm.

This project contains three folders: `data`, `src`, and `scripts`.
- `data`: include five real datasets (.Rdata file) used in the paper.
- `src`: include the source codes (.R file).
- `scripts`: include codes (.R file) to replicate the experiments for Figure 3 in the paper. All outputs necessary for the study can be reproduced with similar codes.

## Replicating
To get the Figures and Tables in the paper, please run the R codes in the `scripts` folder. 

## Ongoing Development

This code is being developed on an on-going basis at the author's [GitHub site](https://github.com/statcp).

## Support

For support in using this software, submit an
[issue](https://github.com/statcp).
