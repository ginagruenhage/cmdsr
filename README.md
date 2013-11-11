# cmdsr

`cmdsr` provides tools to compute and analyze continuous embeddings that are computed using cMDS. Continuous embeddings are a useful way to visualize data that has inherent continuous parameters, such as time or experimental control variables. It is also particularly useful to visualize the effects of changing distance functions. An example could be a family of distance functions that is parameterized by a hyperparameter, such as scale. Applying different distance functions from this family to a dataset might lead to significantly different output. Various effects of such changing distances on data can be visualized using cMDS. 
You can track the development of `cmdsr` at https://github.com/ginagruenhage/cmdsr. `cmdsr` is easy to install using the `devtools` package:
`devtools::install_github("cmdsr","ginagruenhage")`

## Components of cmdsr
The essential function of the package is `cmds()`. It takes a list of distance matrices and gives a list output. The essential component of this output is the list of configuration matrices that contains the coordinates of embedding points. 

The output list can then be summarized using the `summary.cmds()` function and plotted using either `plot.cmds()` for a simple plot of the embedding or `googleVis.cmds()` for an interactive visualization using googleVis charts. 

We recommend taking a look at the *personality* vignette to get started with the package.


## Citation info

The cMDS algorithm is described in: 
Gina Gruenhage & Simon Barthelme,  Visualizing the effects of a changing distance using continuous embeddings, [http://arxiv.org/abs/1311.1911](arxiv.org/abs/1311.1911).
