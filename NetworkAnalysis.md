---
name: NetworkAnalysis
topic: Network Analysis
maintainer: Fabio Ashtar Telarico, Pavel N. Krivitsky, James Hollway
email: Fabio-Ashtar.Telarico@fdv.uni-lj.si
version: 2024-12-30
source: https://github.com/cran-task-views/NetworkAnalysis/
---

This CRAN task view provides a curated list of R packages for analyzing and
modelling networks (also known as _relational data_ or _graphs_). These tools
facilitate the exploration of natural, social, and other phenomena by focusing
on the relationships between entities.

This page lists a number of packages, and sometimes core functions, in several
sections based on their scope and focus:

1. The first section outlines the main ecosystems of R packages that include
basic network-analytic operations such as creating, manipulating, and describing
relational data. Here we also list choices of graphical packages for
visualizing or drawing networks. For those new to network analysis in R, we
recommend starting with the [`igraph`
introduction](https://igraph.org/c/doc/igraph-Introduction.html) or the
[`statnet` tutorial](https://statnet.org/workshop-intro-sna-tools/).

2. Subsequently, packages and functions for advanced network-analytical tasks
are presented. We currently structure these into three subsections: (1)
centrality, (2) community detection, and (3) model-based clustering.

3. Then, packages offering modelling and inferential tools applicable across
disciplines and fields of interest are discussed. A distinction is drawn
between models that are primarily for cross-sectional anddynamic data, with an
extra section on special models for multimodal, multilevel, and multiplex data. 

4. Finally, the focus shifts to packages containing data structures, methods,
and models with a narrower field of application. The list includes some of the
areas where network methods are more widely applied: ecology, bibliometrics,
life and natural sciences, neurosciences, psychology, public health, social
sciences and economics.

The list excludes packages that primarily deal with graph representations of
conditional in/dependence between variables. This includes Bayesian networks
and Markovian graphs, which, despite their relevance to statistical modeling,
are covered under the CRAN task view `r view("GraphicalModels")`. This
distinction keeps the list focused on network analysis to explore broader
relational dynamics.

Some packages could appear under multiple headings because they can perform
multiple tasks (e.g., clustering and visualization). But, for the sake of
brevity, non-core packages are listed only once: in the section that described
each package's main use case.

If you think that a package is missing from the list, please file an issue in
the GitHub repository or contact the maintainer.

## Table of contents

- [Ecosystems](#ecosystems)
   - [Network construction](#network-construction)
   - [Visualization](#visualization)
- [Network analysis](#network-analysis)
   - [Centrality](#centrality)
   - [Community detection](#community-detection)
   - [Model-based clustering ](#model-based-clustering)
- [Statistical modeling](#statistical-modeling)
   - [Cross-sectional networks](#cross-sectional-networks)
   - [Multimodal and multilevel networks](#multimodal-and-multilevel-networks)
   - [Dynamic networks](#dynamic-networks)
- [Field packages](#field-packages)
   - [Ecological networks](#ecological-networks)
   - [Bibliometric networks](#bibliometric-networks)
   - [Networks in the natural and life sciences](#networks-in-the-natural-and-life-sciences)
   - [Neurosciences and psychology](#neurosciences-and-psychology)
   - [Spatial networks](#spatial-networks)
   - [Public-health networks](#public-health-networks)
   - [Social and economic networks](#social-and-economic-networks)
- [References](#references)


## Ecosystems and Data

### Ecosystems

The starting point for analyzing networks in R is to familiarize with the main
package 'families' or ecosystems. Using them, users can access functions to
create, import/export, edit, and otherwise operate on relational data.

- `r pkg("igraph", priority = "core")` provides tools for creating,
  manipulating, and analyzing network structures with a focus on graphical
  representations and fast algorithms to operate on large datasets (particularly
  manipulating data and dealing with vertex attributes). Rather than coming
  bundled with other packages, `r pkg("igraph")` offers an ecosystem of extensions
  and add-ons.

  - *Approach*: This R package is built upon a C library that is shared also by
    implementations in Python and Mathematica. The C code is efficient, and the
    R interface is increasingly consistent and easy to use with a lot of basic
    functionality, including calculating network properties, generating random
    graphs for simulations, etc.

  - *Flexibility*: Many of the functions are provided in two versions: for
    direct assignment and a pipe-able version.
  
  - *Support*: There are tons of online resources answering virtually any
    question concerning _how_ to do almost anything thanks to a large community
    and active maintainers.
  
  - *Extensions*: There are a number of add-on packages that integrate igraph
    with `r pkg("ggplot2")` and other drawing tools or provide sample datasets.
    According to the latest review this is the largest network-analysis
    ecosystem in R by number of extensions (see Kanevsky 2016).

- `r pkg("statnet", priority = "core")`, sometimes used only as `r pkg("sna",
  priority = "core")` and/or `r pkg("network", priority = "core")` is a more
  organic family of packages for statistical network analysis built upon common
  data representations and design choices. But there are also 'add-ons',
  although not as many as in the previous ecosystem. 

  - *Approach*: The approach in this group of packages is more disaggregated and
    involves a number of different packages for different purposes. Not least,
    they provide support to the packages in the `r pkg("ergm")` family.
 
  - *Flexibility*: The packages in `r pkg("statnet")` can be used to analyze
  (social) network data using both direct assignment and a 'piped' approach. 
 
  - *Comprehensiveness*: This ecosystem's biggest advantage is that it allows
    users to carry out most network-analytical operations, especially in dealing
    with social network analysis (SNA). But the trade-off, of course, is that
    the number of possibilities makes life harder for new users.

- `r pkg("manynet", priority = "core")` is built upon many of the other network
  packages in this list and offers interoperability with many different classes
  of objects as well as network visualization and analytic tools.

  - *Approach*: Inspired by the [_tidyverse_](https://www.tidyverse.org/),
    `manynet` offers a piped syntax, with even more structure, consistency, and
    sensible defaults. For example, `manynet::net_*` functions always return a
    single value (scalar), and `manynet::node_*` and `manynet::tie_*` always
    return a vector the length of the nodes or ties in the network,
    respectively. Similarly, `manynet::is_*` functions always return logical
    values (`TRUE`/`FALSE`), for example. This helps new users and those working
    with multiple packages alike.
 
  - *Flexibility*: the package leverages S3 dispatching so that all the
    functions work with many _classes_ of network objects, such as `r pkg
    ("network")` or `r bioc("graph")` objects, but also edge lists, matrices, and
    various other more specific classes.  Its coercion routines retain more
    information than the long-standing alternative `r pkg ("intergraph")`, being
    able to transform between more classes and receiving more frequent updates. 
    It also includes more complete import and export routines than generalist
    packages.
 
  - *Comprehensiveness*: the package wraps most of `r pkg("igraph")`'s offerings
    but extends or corrects them to treat many _types_ of networks, including
    two-mode, multiplex, and dynamic networks. The package also implements a
    number of functions for network analysis that are not available elsewhere.

  - *Documentation*: The package benefits from clear and concise documentation,
    making it accessible to users of all levels. Additionally, it offers
    straightforward tutorials and examples to help users get started with
    network analysis in R.  Further explanation, examples, and references are
    continually being to the documentation to provide user reference and support
    user experience.

### Relational data management and conversion tools

Although the 'core' packages for network analysis in R can create a wide range
of networks from different types of inputs, there are also specialized packages
for constructing more specialized formats or for converting or coercing between
different formats.

- `r pkg("intergraph")` is not a network analysis package _per se_. Rather it
allows to easily convert objects produced by `r pkg("statnet")` packages into
`r pkg("igraph")` objects (or a data frame) and vice versa. Thus, it helps
leveraging multiple packages' functionalities and ensuring compatibility between
several users' workflows too many additional functionalities.

- `r pkg("BoolNet")` provides tools for assembling, analyzing and visualizing
synchronous and asynchronous (probabilistic) Boolean networks as well as simpler
Boolean networks. All the main functions are described in a handy
[vignette](https://cran.r-project.org/web/packages/BoolNet/vignettes/BoolNet_package_vignette.pdf).

- `r pkg("egor")` allows to create ego-centric networks starting from exports
from [_EgoNet_](https://github.cm/egonet/egonet), [_EgoWeb
2.0_](https://www.qualintitative.com/egoweb/) and
[_openeddi_](https://github.cm/jfaganUK/openeddi). It includes a Shiny app and
procedures for creating and visualizing clustered graphs.

- `r pkg("networkDynamic")` from `r pkg("statnet")` allows the management of
dynamic networks.

- `r pkg("ionet")` creates network starting by turning input-output tables into
weighted adjacency matrices.

- `r pkg("rgraph6")` allows to encode relational data (adjacency matrices, edge
lists, `r pkg("network")` and `r pkg("igraph")` objects) as ASCII strings and
vice versa using `graph6`, `sparse6`, and `digraph6`
[formats](http://users.cecs.anu.edu.au/~bdm/data/formats.txt).

- `r pkg("tidygraph")` is designed for handling and manipulating graph data
within the [_tidyverse_](https://www.tidyverse.org/) framework. It does not make
it into "core" packages because it lacks a comprehensive set of tools for
network analysis. Yet, it provides a flexible piped approach to working with
relational data, allowing users to apply familiar data manipulation techniques
from the tidyverse to graphs. Users can easily perform tasks such as filtering,
summarizing, and joining graph data using familiar tidyverse syntax. Given that
both `r pkg ("igraph")` and `r pkg ("sna")` provide piped functions for most
operations, `r pkg ("tidygraph")`'s added value lies mainly in the possibility
of accessing directly either the node data, the edge data or the graph itself
while computing inside verbs.

- `r pkg("backbone")` enables the extraction of the sparse and unweighted
subgraph of a network called 'backbone'.

## Exploratory Data Analysis

`r pkg ("igraph")`, `r pkg ("sna")`, and `r pkg("manynet")` offer functions for
a similar set of network-analytic and visualisation operations, whereas `r pkg ("tidygraph")` is
more limited. However, some algorithms differ from each other and from those
are some specialized packages for their implementation, speed, or defaults.

### Visualization

#### Interactive visualization

More details in the CRAN task view `r view("DynamicVisualizations")`.

- `r pkg("visNetwork")` focuses on interactive network visualization using the
_vis.js_ library. The package allows users to create visually appealing and
interactive network visualizations with features such as zooming, panning, and
node highlighting. It offers a user-friendly interface for creating interactive
network visualizations, making it suitable for un-experienced users.

- `r pkg("networkD3")` provides functions that turns edge lists into a _D3_
JavaScript network, tree, dendrogram, or Sankey plots.

- `r pkg("bipartiteD3")` uses the _D3_ and _viz.js_ libraries for plotting
networks produced with the `r pkg("bipartite")` package.

- `r pkg("ndtv")` renders network objects from the package
`r pkg("networkDynamic")` as videos or interactive animations.

#### Static visualization

- `r pkg("diagram")` was born as a companion to the book _A Practical Guide to
Ecological Modelling_ by K. Soetaert and P.M.J. Herman. But it can visualize as
a flow diagram, a web or grid any network given in the form of a transition
matrix. 

- `r pkg("neatmaps")` tries to simplify the exploratory step of data analysis by
providing function to easily produce hierarchical clustering
(`neatmaps::hierarchy`), consensus clustering (`neatmaps::consClustResTable`)
and heatmaps of multiple networks (`neatmaps::neatmap`).

- `r pkg("manynet")` builds on `r pkg("ggraph")` for graphing networks in
`ggplot2`-style but eases the complicated syntax and offers sensible defaults to
make it easy to visualize and explore network data while retaining the
flexibility to theme. 

  - `manynet::graphr` is for quick, easy network visualization.

  - `manynet::graphs` is for comparing ego networks or subgraphs side by side.

  - `manynet::grapht` is for developing dynamic or longitudinal networks into GIFs.

  - `manynet` includes plotting methods for much of its output, including blockmodels and dendrograms for clustering.

#### Extensions for `ggplot2`

- `r pkg("ggnetwork")` offers geometries to plot `r pkg("network")` objects.

- `r pkg("ggraph")` allows to plot `r pkg("igraph")` objects by building up
plots layer by layer.

- `r pkg("ggsom")` offers functions to plot self-organizing maps (SOMs). 

- `r pkg("snahelper")` is an add-on allowing access to a GUI for visualizing and
analyzing networks. Once the visualization is set, the relevant code is
automatically added to the script. 

- `r pkg("roughnet")` leverages the _rough.js_ library to draw sketchy,
hand-drawn-like networks

- `r pkg("gganimate")` allows to produce GIFs and MP4s version of evolving
`ggplots`, including those representing networks. 

- `r pkg("ggdendro")` makes it easy to make ggplots of dendrograms create using
the functions `tree`, `hclust`, `dendrogram`, and `rpart`.

- `r pkg("multigraph")` is a powerful tool providing easier visualizations of
multigraphs, valued/signed networks, bipartite networks, multilevel networks,
and Cayley graphs with various layout options.

#### Layouts

- `r pkg("ggforce")` offers functions for specialized plots, some of which also
find application in network analysis. Most importantly, alluvial plots can be
used to visualize composition of groups in a dynamic network. 

- `r pkg("graphlayouts")` adds several layout algorithms to `r pkg("igraph")`
and `r pkg("ggraph")` based on the concept of stress majorization (See also
`r pkg("edgebundle")`).

- `r pkg("manynet")` includes a few more layout algorithms for multimodal
networks.

- `r bioc("Rgraphviz")`, available on Bioconductor, creates a direct link
between the `r bioc("graph")` package and the _graphviz_ library. 

- `r pkg("patchwork")` allows for arbitrarily complex composition of plots that
can be used, for example, in visualizing multipartite and other complex
networks.

### Centrality

Both main ecosystems can compute betweenness, eigenvalue, power, and closeness
centrality, but `r pkg ("igraph")` offers more options than `r pkg ("sna")` and
`r pkg ("tidygraph")` overall. In addition:

- `r pkg("centiserve")` adds dozens of centrality measures for `r pkg
("igraph")` objects such as bottleneck, decay, and entropy centrality.

- `r pkg("birankr")` provides optimized functions for estimating various
centrality measures in bipartite/two-mode networks. It can also estimate
efficiently page-rank in one-mode networks, project two-mode networks to
one-mode ones, and convert edge lists and matrices to the `sparseMatrix` format
offered in the package `r pkg("Matrix")`. It supports edge lists (in the
`data.frame`, `data.table::data.table`, or `tidydata::tbl_df` class) and
adjacency matrices (either in the built-in `matrix` class or in `Matrix`'s
`dgCMatrix` class).

- `r pkg("netrankr")` offers index-free centrality rankings via
_neighborhood-inclusion_ or _positional dominance_ and based on probabilistic
methods like computing expected node ranks and relative rank probabilities (see
also Schoch 2018).

## Group detection

### Community detection

- `r pkg("igraph")` is the package of choice for the implementation of most
modularity-based community-detection algorithms. Available approaches include
betweenness, greedy algorithm, infomap, label propagation Leiden, Generalized
Louvain, and walktrap amongst others.

- `r pkg("cencrne")` proposes a regularized network-embedding model to
simultaneously estimate the community structure and the number of communities in
an asymptotically consistent way. The method is mainly used in life sciences but
is applicable across the board.

- `r pkg("linkcomm")` provides functions for generating, visualizing, and
analyzing overlapping communities in networks of arbitrary size and type. Unlike
`r pkg("igraph")` and `r pkg ("network")`, it can compute relatedness
`linkcomm::getClusterRelatedness` and `linkcomm::getCommunityConnectedness` as
well as generate a mesoscopic matrix (`linkcomm::getCommunityMatrix`). Moreover,
it can produce membership for hierarchical communities
(`linkcomm::getNestedHierarchies`).

### Blockmodeling

#### Generalized (structural and/or regular equivalence)

- `r pkg("sna")` implements a simple version of structural-equivalence
blockmodel (`sna::blockmodel`). It can also generate networks with a given
blockmodel as well as print and plot the results.

- `r pkg("concorR")` implements the classical CONCOR (Convergence of iterated
Correlation) algorithm (Breiger, Boorman, and Arabie 1975) for one- and
multi-mode un/directed networks.

- `r pkg("BMconcor")` allows the simultaneous blockmodeling of networks based on
structural and regular equivalence through singular value decomposition (SVD) by
blocks (see Lafosse and Hanafi 1997)

- `r pkg("blockmodeling")`: this package offers and implementation of
generalized blockmodeling (`blockmodeling::optRandomParC`) as well as functions
for computation of (dis)similarities in terms of structural or regular
equivalence and plotting. Furthermore, it includes implementations of the REGE
algorithm (`blockmodeling::REGE`).

- `r pkg("BlockmodelingGUI")` is a Shiny app providing a graphical interface for
generalized blockmodeling of single-relation, one-mode networks from the package
`r pkg("blockmodeling")`. It includes several ways to visualize networks and
partitions using `r pkg("igraph")`, `r pkg("network")`, and more.

- `r pkg("kmBlock")` implements a k-means like approach to the blockmodeling of
one-mode and linked networks as presented in Žiberna (2020).

- `r pkg("dBlockmodeling")` contains functions to apply blockmodeling of signed
(positive and negative weights are assigned to the links), one-mode and valued
one-mode and two-mode.

- `r pkg("signnet")` offers to functions implementing the generalized
blockmodeling with structural equivalence (`signnet::signed_blockmodel`) and
generalized equivalence (`signnet::signed_blockmodel_general`) of signed
networks based on objects from `r pkg("igraph")` 

- `r pkg("oaqc")` enables efficient computation of the orbit-aware quad census
from Ortmann and Brandes (2017, `r doi("10.1007/s41109-017-0027-2")`).

#### Stochastic (SBM)

- `r pkg("igraph")` cannot run SBMs, but it can generate a random graph
according to a specified SBM (`igraph::sample_sbm`) or an arbitrary hierarchical
SBM (`igraph::sample_hierarchical_sbm`)

- `r pkg("blockmodels")` allows to run the SBM or the Latent Block Model (LBM,
an SBM for bipartite networks) of static networks using a Variational
Expectation Maximization algorithm. Various `S4` functions implement three
probability distributions: `blockmodels::BM_bernoulli` for binary data,
`blockmodels::BM_poisson` for discrete/count weights, `blockmodels::BM_gaussian`
for continuous weights. It allows for SBMs and LBM with or without node
covariates and supports multiplex binary networks via
`blockmodels::BM_bernoulli_multiplex`.

- `r pkg("sbm")` is an extension of `blockmodels` for bi- and multi-partite as
well as multiplex networks as proposed by Barbillon and colleagues (2017)
through dedicated `R6` classes. It includes functions to plot the resulting
partition.

- `r pkg ("greed")` leverages a combination of greedy local search and a genetic
algorithm (see Côme et al 2021) to execute (degree-corrected) SBM and LBM.

- `r pkg("dynsbm")`, archived from the CRAN repository on 2023-10-27 due to a
faulty dependence, implements the model for temporal networks presented in
Matias and Miele (2020) which combines a static
SBM with independent Markov chains for the dynamic part. It supports binary and
weighted networks with both discrete and continuous edges. Includes also
functions for plotting (`dynsbm::adjacency.plot`, `dynsbm::alluvial.plot`,
`dynsbm::connectivity.plot`) the partition and automatically constructs matrices
as an array of the right format.

- `r pkg("MLVSBM")` Implements the SBM of multilevel networks where the
different matrices each represent an interaction layer either weighter or
binary. It generalizes the approach proposed by Chabert-Liddell et al. (2021) to
more than two layers.

- `r pkg("StochBlock")` implements the stochastic blockmodeling of one-mode and
linked networks as implemented in Škulj and Žiberna (2022). It includes utilities to plot the results
but cannot choose automatically the 'right' number of clusters and tends to be
very slow according to subsequent reviews (see Cugmas and Žiberna 2023).

- `r pkg("GREMLINS")` implements the SBM of generalized multipartite networks
where the different matrices each involve nodes that can be partitioned into 
a-priori defined _functional groups_ as proposed by Bar-Hen, Barbillon, and 
Donnet (2018).

### Others

- `r pkg("clustNet")` allows to cluster units in a network using a Bayesian
mixture model that can account for node and edge covariates.

- `r pkg("collpcm")` provides Monte-Carlo Markov Chain (MCMC) inference for
collapsed latent space models that allow to search over the model space,
including deciding on the number of clusters.

- `r pkg("graphclust")` implements an agglomerative algorithm to maximize the
integrated classification likelihood criterion and a mixture of stochastic block
models based on `r pkg("igraph")` objects.

- `r pkg("latentnet")` provides functions to fit and simulate latent position
and cluster model using `r pkg("network")` objects and compatibly with
`r pkg("ergm")` approaches.

- `r pkg("latenetwork")` implements Hoshino and Yanagi's (2023) method for the
causal inference with noncompliance and network interference of unknown form on
average causal using instrumental variables.

- `r pkg("netClust")` provides a function to cluster one-layer
(`netClust::netEM_unilayer`) and multilayer (`netClust::netEM_multilayer`)
networks by means of finite mixtures and expectation-maximization.


## Statistical modeling

Statistical modelling in network analysis enables researchers to uncover
patterns, test hypotheses, and make predictions about network structures and
dynamics. This section introduces R packages that support a range of statistical
approaches, from modelling static (cross-sectional) networks to analyzing
dynamic, multimodal, and multilevel networks. These methods provide tools to
infer underlying processes that generate observed network data, assess the
significance of observed patterns, and simulate network structures under various
conditions.

### Cross-sectional networks

- `r pkg("ergm")` from the `r pkg("statnet")` ecosystem provides functions to fit, simulate and analyze
exponential-family random graph models (ERGM). Depending on specific needs,
several specialized extensions are available.

  | Use case                                    | Package               |
  |:--------------------------------------------|:----------------------|
  | Count weights                               | `r pkg("ergm.count")` |
  | Egocentrically sampled networks             | `r pkg("ergm.ego")`   |
  | Multilayer networks and samples of networks | `r pkg("ergm.multi")` |
  | Rank-order networks                         | `r pkg("ergm.rank")`  |
  | Modeling ERGM-generating processes          | `r pkg("ergmgp")`     |
  | Fit ERGM to small networks                  | `r pkg("ergmgp")`     |
  | Small hierarchical ERGMs                    | `r pkg("lightergm")`  |
  | Large hierarchical ERGMs                    | `r pkg("bigergm")`    |

- `r pkg("amen")` offers additive and multiplicative effect (AME) models with
regression terms, covariance structure of the social relations model (Warner,
Kenny and Stoto (1979), and multiplicative factor models (Hoff 2009).
It supports binary networks as well as valued ones (assuming a Gaussian,
zero-inflated/tobit, ordinal, or fixed-rank nomination model)

- `r pkg("bootnet")` implements bootstrap procedures to assess accuracy and
stability of estimated structures and centrality indices on undirected networks
(For an alternative see `r pkg("localboot")`).

- `r pkg("fastnet")` allows to simulate large-scale social networks and retrieve
their most relevant metrics following the approach proposed by Dong, Castro,
and Shaikh (2020).

- `r pkg("nda")` gathers non-parametric dimensionality-reduction functions 
with/out (automated) feature selection and limited plotting capabilities.

- `r pkg("lolog")` implements Latent Order Logistic (LOLOG) models, a network
formation process in which edges are added one at a time drawn from a
distribution conditional on edges already added, with order unknown.

- `r pkg("MoNAn")` implements the method to analyze the structure of weighted
mobility networks or distribution networks outlined in Block et al (2022).

- `r pkg("ERPM")` extends exponential random graph models (ERGMs) to explain
cross-sectional or longitudinal observed partitions, i.e. sets of
non-overlapping groups, such as face-to-face interactions, animal herds, or
political coalitions, through group formation processes based on individual
attributes, relations between individuals, and size-related factors.

### Multimodal and multilevel networks

- `r pkg("migraph")` is an `r pkg("igraph")` extension to analyze multimodal
networks as described in Knoke, Diani, Hollway, and Christopoulos (2021).

- `r pkg("multinets")` is an `r pkg ("igraph")` extension to analyze multilevel
networks as described in Lazega and Snijders (2016).

- `r pkg("multiplex")` makes possible, among other things, to create and
manipulate multiplex, multimode, and multilevel network data with different
formats (Rivero Ostoic 2020, `r doi("10.18637/jss.v092.i11")`).

- `r pkg("dyads")` offers functions for the MCMC simulation of dyadic network
models j2, p2 (also multilevel) and b2 model.

- `r pkg("tnet")` includes functions for analyzing two-mode, weighted, and
longitudinal networks.

- `r pkg("incidentally")` implements methods to generate incidence matrices as
described in Neal (2022: `r doi("10.31219/osf.io/ectms")`).

### Dynamic networks

The following packages focus on modeling and simulation of networks that evolve
over time and network processes that occur over time.

#### Relational events

Relational event data contains information about exact times during which the
nodes interact. This is commonly observed for e-mail, radio, and other
communications.

- `r pkg("rem")` and `r pkg("relevent")` both contain functions to fit and
simulate dyad-oriented relational event models. But only `r pkg("relevent")` can
estimate event sequence data without time stamps.

- `r pkg("goldfish")` offers functions to fit and simulate actor-oriented
dynamic network actor models and dyad-oriented relational event models as
described in Stadtfeld et al. (2017).

#### Discrete observations

The following package are focused on modeling series of networks, also known as
panel data.

- `r pkg("tergm")` a set of extensions for `r pkg("ergm")` for fitting and
simulating discrete-time models for series of networks (or a long-term
equilibrium of a discrete-time network process) where each time step is modeled
as a draw from an ERGM conditional on the prior time steps.

- `r pkg("dnr")` estimation of discrete-time models for series of networks where
each time-step is modeled as a draw from an ERGM conditional on prior time
steps, subject to the constraint that within each time step, edge variables are
independent. Varying node sets are also supported.

- `r pkg("btergm")` bootstrap inference for discrete-time models for series of
networks where each time step is modeled as a draw from an ERGM conditional on
the prior time steps.

- `r pkg("RSiena")` estimation of continuous-time Stochastic Actor-Oriented
Models (SAOMs) for panel network data.

- `r pkg("idopNetwork")` implments the model proposed by Cao et al (2022) to
convert static data into their 'dynamic' form contextually inferring
informative, dynamic, multi-directional networks with clusterable structures.


#### Diffusion on networks

- `r pkg("tsna")` is used to analyze them.

- `r pkg("EpiModel")` allows to simulate mathematical models of infectious
disease dynamics.

- `r pkg("manynet")` can manipulate, visualize, and analyze longitudinal and
network event data, including running contagion/diffusion processes and
compartmental models.

- `r pkg("netdiffuseR")` was developed as part of Valente et al. (2015:
`r doi("10.1016/j.socscimed.2015.10.001")`). The package is for empirical
statistical analysis, visualization and simulation of network diffusion and
contagion processes. It implements algorithms for calculating network diffusion
statistics such as transmission rate, hazard rates, exposure models, network
threshold levels, infectiousness (contagion), and susceptibility.

## Field packages

As an interdisciplinary approach, network analysis is used in a number of
fields, where the specific needs and interests of those fields are addressed by
particular packages.

### Ecological networks

- `r pkg("econetwork")` is a collection of advanced functions to analyze and
models of ecological networks (mainly food webs and host-parasite relations, but
also plant-pollinator and other mutualistic ones) statically and dynamically
(based on Ohlmann _et al._ 2023). 

- `r pkg("AnimalHabitatNetwork")` provides functions for generating and
visualizing networks representing the physical configurations of animal
habitats. It implements an original network-generating algorithm based on
pair-wise Euclidean distances and can output undirected network either weighted
or binary, fully connected or sparse). The package is associated with a PDF on
modelling the physical configurations of animal habitats using networks.

- `r pkg("aniSNA")` allows to obtain network structures from animal GPS
telemetry observations and statistically analyze them to assess their adequacy
for social network analysis. Methods include pre-network data permutations,
bootstrapping techniques to obtain confidence intervals for global and
node-level network metrics, and correlation and regression analysis of the local
network metrics.

- `r pkg("asnipe")` implements several tools that are used in animal social
network analysis to cluster, and generate networks, perform permutation tests,
calculate association rates, and perform multiple regression analysis.

- `r pkg("ATNr")` estimates _allometric trophic models_ (ATN) for the species
biomasses in dynamic food-webs and allows to generate synthetic networks. It
also provides access to the ODE solver deSolve.

- `r pkg("BIEN")` allows to access the _Botanical Information and Ecology
Network Database_ in R

- `r pkg("bipartite")` offers functions to visualize food webs and calculate
some ecological indices on them.

- `r pkg("cassandRa")` deals with under-sampling in ecological networks by
fitting a variety of statistical models and sample coverage estimators to
correct for (likely) missing ties. It works only on bipartite networks.

- `pkg("EcoNetGen")` to simulate and sample from ecological networks.

- `pkg("econullnetr")` to carry out null-model analysis for ecological
networks. 

### Bibliometric networks

- `r pkg("bibliometrix")` includes functions to import bibliographic data from
the main publication databases online ('SCOPUS', 'Clarivate Analytics Web of
Science', 'Digital Science Dimensions', 'Cochrane Library', 'Lens', and
'PubMed'). It can also build networks (`bibliometrix::biblioNetwork`) for
co-citation, coupling, scientific collaboration and co-word analysis including
their dynamic versions (`bibliometrix::histNetwork`). It allows to plot the
data using `VOSviewer.jar`.

- `r pkg("bibliometrixData")` contains example datasets for testing
`bibliometrix`.

- `r pkg("biblionetwork")` proposes functions to identify and weight the edges
in a bibliometric network. All functions are optimized for large datasets. It
implements different methods for different types of relations (based on
Perianes-Rodriguez, Waltman, and Van Eck 2016; Leydesdorff and Park 2017):
Co-authorship supports simple counting, (refined) fractional weight with or with
cosine normalization. Bibliographic coupling supports: coupling strength and
angle. Co-citation supports the cosine normalization of count weights.

- `r pkg("Diderot")` is geared towards the analysis of citation networks using
modularity and heterocitation metrics based on Scopus data.

### Networks in the natural and life sciences

- `r bioc("Rcy3")` provides access to [Cytoscape](https://cytoscape.org/), one
of the most used network tool in the field of molecular biology, allowing to
vizualize, analyze and explore networks using a single function for each
operation executable through Cytoscape's graphical interface.

- `r pkg("WGCNA")` focuses on the analysis of weighted correlation networks. It
has functions for network construction, modularity computation, gene selection,
topological analysis, generating data, plotting, and exports to third-party
software. Notably, the underlying data mining approach has been used beyond the
natural sciences. There are several packages on Bioconductor that
reverse-depend/extend these functionalities.

- `r pkg("c3net")` allows to infer gene-regulation networks with direct physical
interactions using `C3NET` (Altay and Emmert-Streib 2010). Other packages
implement improvements/variants of this algorithm based on the literature, such
as:

- `r pkg("Ac3net")` Infers directional conservative causal core in gene network
based on the algorithm for directional network proposed by Altay (2018).

- `r pkg("bc3net")` implements the BC3NET algorithm for inference on
gene-regulation networks (Simoes and Emmert-Streib 2012). In essence it offers
a Bayesian approach with noninformative prior to the C3NET algorithm.

- `r bioc("BioNAR")` implements a detailed topologically based network analysis
with functions that create networks based on laboratory-produced meta-data. It
includes functions for vertex centrality measure and modularity computation.
Additionally, it provides a robust synaptic proteome network for data
validation.

- `r pkg("BASiNET")` and `r pkg("BASiNETEntropy")` provide functions for
classifying RNA sequences using network algorithms and notions from information
theory.

- `r pkg("bionetdata")` is a collection of relation datasets of biological and
chemical nature.

- `r pkg("Cascade")` includes functions for gene selection, reverse engineering,
and prediction in cascade networks.

### Neurosciences and psychology

- `r pkg("NetworkToolbox")` implements network analysis and graph theory
measures used in neuroscience, cognitive science, and psychology. Methods
include various filtering methods and approaches such as threshold and
dependency. It can also execute some basic operations such as computing
centrality of nodes and community or the network's clustering coefficient.

- `r pkg("qgraph")` provides tools for visualizing and analyzing weighted
networks and a Gaussian graphical model for plotting. It is compatible with
`r pkg ("igraph")` through the `qgraph::as.igraph.qgraph` function. It is mostly
used in psychology and neurosciences.

- `r pkg("HospitalNetwork")` provides functions to construct a one-mode network
of hospitals based on the linked two-mode networks of hospitalized patients'
transfers.


### Spatial networks

- `r pkg("geonetwork")` handles networks or graphs whose nodes are locations.
The functions includes the creation of objects of class `geonetwork` as a graph
with node coordinates, the computation of network measures, the support of
spatial operations (projection to different coordinate reference systems,
handling of bounding boxes, etc.) and the plotting of the `geonetwork` object
combined with supplementary cartography for spatial representation. It is
compatible with `r pkg ("igraph")`.

- `r pkg("sfnetworks")` combines the work-horse spatial-data package in R (`sf`)
and `r pkg("tidygraph")` for tidyverse-friendly classes and routines (shortest
paths, cleaning, and editing) for geospatial networks.

- `r pkg("chessboard")` provides functions to work with un/directed undirected
spatial networks. It allows to create connectivity matrices
(`chessboard::connectivity_matrix`) and exports results to several formats: node
list, neighbor list, edge list, connectivity matrix, Eigenvector maps. It also
implements connectivity for chess pieces via specific functions:
`chessboard::bishop`, `chessboard::knight`, `chessboard::pawn`,
`chessboard::queen`, `chessboard::rook`, besides introducing two sets of
movement rules `chessboard::fool` and `chessboard::wizard`.

- `r pkg("epanet2toolkit")` interfaces R with the
[EPANET](https://github.cm/OpenWaterAnalytics/EPANET) programmer's toolkit to
carry out basic (`epanet2toolkit::ENepanet`) or customized
(`epanet2toolkit::ENopen`) simulations.

- `r pkg("intensitynet")` includes functions to analyze point patterns in space
occurring over planar network structures derived from graph-related intensity
measures for un/directed and mixed networks

### Public-health networks

- `r pkg("epinet")` simulates contact networks to predict the transmission of
contagious diseases through Bayesian inference.

- `r pkg("hybridModels")` offers a meta-population model that assigns nodes to
sub-populations to better model disease spreading through cluster contagion
using stochastic simulation algorithm and an individual-based approach.

- `r pkg("netdiffuseR")` provides functions for calculating network effects such
as transmission rate, hazard rates, exposure models, network threshold levels,
infectiousness (contagion), and susceptibility. 

- `r pkg ("EpiModel")` builds on Statnet's for epidemic modelling. But more on
this field of application can be found in the CRAN task view
`r view("Epidemiology")`.

### Social and economic networks

- `r pkg ("sna")` implements many operations commonly carried out on networks in
the social and economic sciences with the ability of regress a network variable
on others using ordinary least square, linear network autocorrelation models or
a logistic regression (More on this type of applications can be found in the
CRAN task view `r view("GraphicalModels")`.

- `r pkg("FinNet")` provides classes, methods, and functions to deal with
financial networks involving both physical and legal persons. The package
assists in creating various types of financial networks: ownership, board
interlocks, or both. It support different tie-weighting procedures (valued or
binary), and renders them in the most common formats (adjacency matrix,
incidence matrix, edge list, `r pkg ("igraph")`, `r pkg ("statnet")`).

- `r pkg("ITNr")` gathers functions to clean and process international trade
data into an adjacency matrix. It can also extract the network's backbone,
compute centrality, run blockmodels and other clustering procedures, or
highlight regional trade patterns.

- `r pkg("modnets")` models moderator variables in cross-sectional, temporal,
and multi-level networks.

- `r pkg("multinet")` provides functions for the creation/generation and
analysis of multi-layer social networks


## References

- Altay, Gokmen. 2018. "Directed Conservative Causal Core Gene Networks".
_bioRxiv_. `r doi("10.1101/271031")`

- Altay, Gökmen, and Frank Emmert-Streib. 2010. "Inferring the Conservative
Causal Core of Gene Regulatory Networks". _BMC Systems Biology_ 4(1): 132.
`r doi("10.1186/1752-0509-4-132")`

- Barbillon, Pierre, Sophie Donnet, Emmanuel Lazega, and Avner Bar-Hen. 2017.
"Stochastic Block Models for Multiplex Networks: An Application to a Multilevel
Network of Researchers". _Journal of the Royal Statistical Society Series A:
Statistics in Society_ 180 (1): 295-314. `r doi("10.1111/rssa.12193")`

- Bar-Hen, Avner, Pierre Barbillon, and Sophie Donnet. 2022. "Block Models for
Generalized Multipartite Networks: Applications in Ecology and Ethnobiology".
_Statistical Modelling_ 22 (4): 273-96. `r doi("10.1177/1471082X20963254")`

- Batagelj, Vladimir, and Matjaž Zaveršnik. "Fast Algorithms for Determining
(Generalized) Core Groups in Social Networks". _Advances in Data Analysis and
Classification_ 5, no. 2 (1 July 2011): 129-45.
`r doi("10.1007/s11634-010-0079-y")`

- Block, Per, Christoph Stadtfeld, and Garry Robins. 2022. "A Statistical Model
for the Analysis of Mobility Tables as Weighted Networks with an Application to
Faculty Hiring Networks". _Social Networks_ 68: 264-278.
`r doi("10.1016/j.socnet.2021.08.003")`

- Breiger, Ronald L, Scott A Boorman, and Phipps Arabie. 1975. "An Algorithm for
Clustering Relational Data with Applications to Social Network Analysis and
Comparison with Multidimensional Scaling". _Journal of Mathematical Psychology_
12 (3): 328-83. `r doi("10.1016/0022-2496(75)90028-0")`

- Cao, Xiaocang, Ang Dong, Guangbo Kang, Xiaoli Wang, Liyun Duan, Huixing Hou,
Tianming Zhao, et al. 2022. "Modeling Spatial Interaction Networks of the Gut
Microbiota". _Gut Microbes_ 14(1): 2106103.
`r doi("10.1080/19490976.2022.2106103")`

- Chabert-Liddell, Saint-Clair, Pierre Barbillon, Sophie Donnet, and Emmanuel
Lazega. 2021. "A Stochastic Block Model Approach for the Analysis of Multilevel
Networks: An Application to the Sociology of Organizations". _Computational
Statistics & Data Analysis_ 158 (June): 107179.
`r doi("10.1016/j.csda.2021.107179")`

- Côme, E., Jouvin, N., Latouche, P. et al. 2021. "Hierarchical Clustering with
Discrete Latent Variable Models and the Integrated Classification Likelihood".
_Advances in Data Analysis and Classification_ 15: 957-986.
`r doi("10.1007/s11634-021-00440-z")`

- Cugmas, Marjan, and Aleš Žiberna. 2023. "Approaches to Blockmodeling Dynamic
Networks: A Monte Carlo Simulation Study". _Social Networks_ 73 (May): 7-19.
`r doi("10.1016/j.socnet.2022.12.003")`

- Dong, Xu, Luis Castro, and Nazrul Shaikh. 2020. "fastnet: An R Package for
Fast Simulation and Analysis of Large-Scale Social Networks". _Journal of
Statistical Software_ 96: 1-23. `r doi("10.18637/jss.v096.i07")`

- Hoff, Peter D. 2009. "Multiplicative Latent Factor Models for Description and
Prediction of Social Networks". _Computational and Mathematical Organization
Theory_ 15(4): 261-72. `r doi("10.1007/s10588-008-9040-4")`

- Hoffman, Marion, Per Block, Timon Elmer, and Christoph Stadtfeld. 2020. "A
Model for the Dynamics of Face-to-Face Interactions in Social Groups". _Network
Science_ 8(S1): S4-S25. `r doi("10.1017/nws.2020.3")`

- Hoshino, Tadao, and Takahide Yanagi. "Causal Inference with Noncompliance and
Unknown Interference". _arXiv_, 22 Octover 2023.
`r doi("10.48550/arXiv.2108.07455")`.

- Kanevsky, Gregory. 2016. "R Graph Objects: igraph vs. network. _R Bloggers_.
January 30, 2016.
<https://www.R-bloggers.com/2016/01/r-graph-objects-igraph-vs-network/>

- Knoke, David, Mario Diani, James Hollway, and Dimitris Christopoulos. 2021.
_Multimodal Political Networks_. Cambridge: Cambridge University Press.
`r doi("10.1017/9781108985000")`

- Lafosse, R., and Hanafi, M. 1997. "Concordance d'un tableau avec $K$ tableaux:
Définition de $K + 1$ uples synthétiques." _Revue de Statistique Appliquée_
45(4): 111-126. <https://eudml.org/doc/106424>

- Lazega, Emmanuel, and Tom A.B. Snijders, eds. 2016. _Multilevel Network
Analysis for the Social Sciences_. Cham: Springer International Publishing.
`r doi("10.1007/978-3-319-24520-1")`

- Leydesdorff, Loet, and Han Woo Park. 2017. "Full and Fractional Counting in
Bibliometric Networks". _Journal of Informetrics_ 11(1): 117-20.
`r doi("10.1016/j.joi.2016.11.007")`

- Martin, Evan A., and A. Fu. 2019. "Approximate Bayesian Inference of Directed
Acyclic Graphs in Biology with Flexible Priors on Edge States".
`r doi("10.48550/arXiv.1909.10678")`

- Matias, Catherine, and Vincent Miele. 2017. "Statistical Clustering of
Temporal Networks through a Dynamic Stochastic Block Model". _Journal of the
Royal Statistical Society. Series B (Statistical Methodology)_ 79 (4): 1119-41.
<https://www.jstor.org/stable/26773154")`

- Ohlmann, Marc, Catherine Matias, Giovanni Poggiato, Stéphane Dray, Wilfried
Thuiller, and Vincent Miele. 2023. "Quantifying the Overall Effect of Biotic
Interactions on Species Distributions along Environmental Gradients".
_Ecological Modelling_ 483: 110424. `r doi("10.1016/j.ecolmodel.2023.110424")`

- Perianes-Rodriguez, Antonio, Ludo Waltman, and Nees Jan Van Eck. 2016.
"Constructing Bibliometric Networks: A Comparison between Full and Fractional
Counting". _Journal of Informetrics_ 10(4): 1178-95.
`r doi("10.1016/j.joi.2016.10.006")`

- Schoch, David. 2018. "Centrality without Indices: Partial Rankings and Rank
Probabilities in Networks". _Social Networks_ 54: 50-60.
`r doi("10.1016/j.socnet.2017.12.003")`

- Simoes, Ricardo de Matos, and Frank Emmert-Streib. 2012. "Bagging Statistical
Network Inference from Large-Scale Gene Expression Data". _PLOS ONE_ 7(3):
e33624. `r doi("10.1371/journal.pone.0033624")`

- Škulj, Damjan, and Aleš Žiberna. 2022. "Stochastic Blockmodeling of Linked
Networks". _Social Networks_ 70: 240-52. `r doi("10.1016/j.socnet.2022.02.001")`

- Stadtfeld, Christoph, James Hollway, and Per Block. 2017. "Dynamic Network
Actor Models: Investigating Coordination Ties through Time". _Sociological
Methodology_ 47(1): 1-40. `r doi("10.1177/0081175017709295")`

- Warner, Rebecca M., David A. Kenny, and Michael Stoto. 1979. "A New Round
Robin Analysis of Variance for Social Interaction Data." _Journal of Personality
and Social Psychology_ 37(10): 1742-57. `r doi("10.1037/0022-3514.37.10.1742")`

- Žiberna, Aleš. 2020. "K-Means-Based Algorithm for Blockmodeling Linked
Networks". _Social Networks_ 61: 153-69. `r doi("10.1016/j.socnet.2019.10.006")`


### Links

- [`igraph` introduction](https://igraph.org/c/doc/igraph-Introduction.html)

- [`statnet` tutorial](https://statnet.org/workshop-intro-sna-tools/)
