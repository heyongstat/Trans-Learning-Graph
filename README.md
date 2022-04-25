# Trans-Learning-Graph
Trans-Copula-CLIME
The code implements the Trans-Copula-CLIME method in  "Transfer Learning in High-dimensional Semi-parametric Graphical Models with Application to Brain Connectivity Analysis" by He Yong, Li Qiushi, Hu Qinqin and Liu Lei. 



Transfer learning has drawn growing attention with the target of improving statistical efciency of one study (dataset) by digging up information from similar
and related auxiliary studies (datasets). In the article, we consider transfer learning problem in estimating undirected semi-parametric graphical model. We propose
an algorithm called Trans-Copula-CLIME for estimating an undirected graphical
model while uncovering information from similar auxiliary studies, characterizing
the similarity between the target graph and each auxiliary graph by the sparsity of
a divergence matrix. The proposed method relaxes the restrictive Gaussian distribution assumption, which deviates from reality for the fMRI dataset related to Attention
Defcit Hyperactivity Disorder (ADHD) considered here. Nonparametric rank-based
correlation coefcient estimators are utilized in the Trans-Copula-CLIME procedure to achieve robustness against normality. We establish the convergence rate
of the Trans-Copula-CLIME estimator under some mild conditions, which demonstrates that if the similarity between the auxiliary studies and the target study is
sufciently high and the number of informative auxiliary samples is sufciently
large, the Trans-Copula-CLIME estimator shows great advantage over the existing
non-transfer-learning ones. Simulation studies also show that Trans-Copula-CLIME
estimator has better performance especially when data are not from Gaussian distribution. Finally, the proposed method is applied to infer functional brain connectivity
pattern for ADHD patients in the target Beijing site by leveraging the fMRI datasets
from some other sites.
KEYWORDS:



See the arxiv version [arXiv:2112.13356].
