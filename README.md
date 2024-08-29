# Low-Complexity Distributed XL-MIMO for Multiuser Detection

This is a research-oriented code package that is primarily intended to allow readers to replicate the results of the article mentioned below and also encourage and accelerate further research on this topic:

V. C. Rodrigues, A. Amiri, T. Abrão, E. de Carvalho and P. Popovski, "Low-Complexity Distributed XL-MIMO for Multiuser Detection," 2020 IEEE International Conference on Communications Workshops (ICC Workshops), Dublin, Ireland, 2020, pp. 1-6, doi: 10.1109/ICCWorkshops49005.2020.9145378. Available on: https://ieeexplore.ieee.org/abstract/document/9145378.

The package is based on the Matlab language and can, in fact, reproduce all the numerical results and figures discussed in the article. To contextualize, in the sequel, we present the abstract of the article and other important information.

I hope this content helps in your research and contributes to building the precepts behind open science. Remarkably, in order to boost the idea of open science and further drive the evolution of science, we also motivate you to share your published results to the public.

If you have any questions and if you have encountered any inconsistency, please do not hesitate to contact me via victorcroisfelt@gmail.com.

## Abstract
In this paper, the zero-forcing and regularized zero-forcing schemes operating in crowded extra-large MIMO (XL-MIMO) scenarios with a fixed number of subarrays have been emulated using the randomized Kaczmarz algorithm (rKA). For that, non-stationary properties have been deployed through the concept of visibility regions when considering two different power normalization methods of non-stationary channels. We address the randomness design of rKA based on the exploitation of spatial non-stationary properties. Numerical results show that, in general, the proposed rKA-based combiner applicable to XL-MIMO systems can considerably decrease computational complexity of the signal detector by paying with small performance losses.

## Content
The codes provided herein can be used to simulate Figs. 2a, 2b, and 3. This is done by running the scripts that have "simulation" in their names, while those with "function" in their names are called by the main scripts. Further details about each file can be found inside them.

## Acknowledgments
This research was supported in part by the Coordenação de Aperfeiçoamento de Pessoal de Nível Superior (CAPES) under grant 88887.461434/2019-00.

## Citing this Repository and License
This code is subject to the MIT license. If you use any part of this repository for research, please consider citing our aforementioned work.

```bibtex
@INPROCEEDINGS{9145378,
  author={Rodrigues, Victor Croisfelt and Amiri, Abolfazl, and Abrão, Taufik and de Carvalho, Elisabeth and Popovski, Petar},
  booktitle={2020 IEEE International Conference on Communications Workshops (ICC Workshops)}, 
  title={Low-Complexity Distributed XL-MIMO for Multiuser Detection}, 
  year={2020},
  volume={},
  number={},
  pages={1-6},
  keywords={Antenna arrays;MIMO communication;Computational complexity;Signal to noise ratio;Optimization;Convergence},
  doi={10.1109/ICCWorkshops49005.2020.9145378}
}

