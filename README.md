# A multi-omics temporal atlas during cell transdifferentiation uncovers causal relationships between chromatin and gene expression
Beatrice Borsari<sup>1,@,</sup>\*, Silvia González-López<sup>1,</sup>\*, Amaya Abad<sup>1</sup>, Vasilis F. Ntasis<sup>1</sup>, Cecilia C. Klein<sup>1,2,&</sup>, Diego Garrido-Martín<sup>1,2</sup>, Carme Arnan<sup>1</sup>, Alexandre Esteban<sup>1,$</sup>, Emilio Palumbo<sup>1</sup>, Marina Ruiz-Romero<sup>1</sup>, Raül G. Veiga<sup>1,%</sup>, Maria Sanz<sup>1,£</sup>, Bruna R. Correa<sup>1</sup>, Rory Johnson<sup>1,#</sup>, Sílvia Pérez-Lluch<sup>1,†</sup> and Roderic Guigó<sup>1,3,†</sup>

<sup>1</sup> Centre for Genomic Regulation (CRG), The Barcelona Institute for Science and Technology, Barcelona (BIST), Catalonia, Spain.

<sup>2</sup> Departament de Genètica, Microbiologia i Estadística, Facultat de Biologia and Institut de Biomedicina (IBUB), Universitat de Barcelona, Barcelona 08028, Catalonia, Spain.

<sup>3</sup> Universitat Pompeu Fabra (UPF), Barcelona, Catalonia, Spain.

<sup>@</sup> Current address: Departament de Genètica, Microbiologia i Estadística, Universitat de Barcelona (UB), Barcelona 08028, Catalonia, Spain.

<sup>&</sup> Current address: Discovery and Translational Science Consulting, Life Sciences and Healthcare, Clarivate, Barcelona, Catalonia, Spain.

<sup>$</sup> Current address: “la Caixa” Foundation, Department of Research and Innovation, Barcelona 08028, Catalonia, Spain.

<sup>%</sup> Current address: Vall d’Hebron Research Institute (VHIR), Vall d’Hebron Barcelona Hospital Campus, Barcelona, Catalonia, Spain.

<sup>#</sup> Current address: School of Medicine, University College Dublin, Ireland; Conway Institute, University College Dublin; Systems Biology Ireland, University College Dublin, Ireland.

\* Equal contribution

<sup>†</sup> Corresponding authors


## Abstract
Chromatin marking by post-translational modifications of histone tails is known to be associated with gene expression. However, whether marking is the cause or the consequence of expression remains controversial[^1][^2][^3][^4][^5][^6][^7][^8]. Temporal series are a powerful approach to assess causality[^9]. They are less disruptive than perturbation experiments, preserving the natural state of the system. Here, we generated densely-spaced multi-omics maps in a time-course, cell-homogeneous, human  transdifferentiation system that occurs with dramatic transcriptomic and epigenomic changes. We used these maps first to characterize the temporal dynamics of chromatin marking at candidate _cis_-Regulatory Elements (cCREs) produced by the ENCODE project[^10]. We found that most of these regions are in a limited number of chromatin states (combinations of chromatin marks) and that temporal transitions between states follow a limited number of temporal chromatin trajectories. Then, we used transitions between chromatin states as interventions in causal inference methods to assess the causal impact of chromatin marking on gene expression. We found that the sequential deposition of H3K4me1 and H3K4me2, and eventually of H3K27ac, plays a causal role in triggering gene activation. The subsequent deposition of H3K9ac and H3K4me3 further increases expression, but occurs after gene activation. We also found that gene activation, in turn, has a causal impact promoting the deposition of most canonically activating chromatin marks. Our model serves as an initial framework to integrate apparently contradictory observations in the field.  We believe that our results demonstrate that time-series data is a powerful approach, complementary to perturbation experiments, to untangle the complex causal relationships characterizing biological systems.

[^1]: Dorighi, K. M. et al. Mll3 and Mll4 Facilitate Enhancer RNA Synthesis and Transcription from Promoters Independently of H3K4 Monomethylation. Molecular Cell 66, 568–576 (2017).
[^2]: Rickels, R. et al. Histone H3K4 monomethylation catalyzed by Trr and mammalian COMPASS-like proteins at enhancers is dispensable for development and viability. Nature Genetics 49, 1647–1653 (2017).
[^3]: Douillet, D. et al. Uncoupling histone H3K4 trimethylation from developmental gene expression via an equilibrium of COMPASS, Polycomb and DNA methylation. Nature Genetics 52, 615–625 (2020).
[^4]: Policarpi, C., Munafò, M., Tsagkris, S., Carlini, V. & Hackett, J. A. Systematic epigenome editing captures the context-dependent instructive function of chromatin modifications. Nature genetics 56, 1168–1180 (2024).
[^5]: Hogg, S. J. et al. Targeting histone acetylation dynamics and oncogenic transcription by catalytic P300/CBP inhibition. Molecular cell 81, 2183–2200 (2021).
[^6]: Zhang, T., Zhang, Z., Dong, Q., Xiong, J. & Zhu, B. Histone H3K27 acetylation is dispensable for enhancer activity in mouse embryonic stem cells. Genome Biology 21 (2020).
[^7]: Boileau, R. M., Chen, K. X. & Blelloch, R. Loss of MLL3/4 decouples enhancer H3K4 monomethylation, H3K27 acetylation, and gene activation during embryonic stem cell differentiation. Genome biology 24 (2023).
[^8]: Kubo, N. et al. H3K4me1 facilitates promoter-enhancer interactions and gene activation during embryonic stem cell differentiation. Molecular cell 84, 1742–1752 (2024).
[^9]: Bar-Joseph, Z., Gitter, A. & Simon, I. Studying and modelling dynamic biological processes using time-series gene expression data. Nature reviews. Genetics 13, 552–564 (2012).
[^10]: Moore, J. E. et al. An expanded registry of candidate cis-regulatory elements. Nature (2026).
