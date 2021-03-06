MLTOOLS software
Version 0.138		Friday 22 Jul 2011 at 16:07

This toolbox provides various machine learning tools, either through wrapping other toolboxes (such as NETLAB) or providing the tool directly. It was designed originally as part of splitting the GPLVM and FGPLVM toolboxes.

Version 0.138
-------------

Minor tweak of model write result and model load result to allow specification of the data loading function.

Version 0.137
-------------

Release for latest version of VARGPLVM toolbox.

Version 0.136
-------------

Minor Mods.

Version 0.135
-------------

Minor mods.


Version 0.134
-------------

Added pmvu model.

Version 0.133
-------------

Added functionality for writing model files using modelDeconstruct commands to keep written files smaller.

Version 0.132
-------------

Add click visualise functionality for LVM visualization, Laplacian eigenmaps and wrapper for MVU. 

Version 0.1311
--------------

Minor change to lvmScatterPlot to fix bug caused when minimum values were positive.

Version 0.131
-------------

Minor changes to toolbox to fix reading in of files written by C++ code.

Version 0.13
------------

Added paramNameRegularExpressionLookup.m to regular expression match a parameter name in a model and return the associated indices. paramNameReverseLookup.m does the same thing but for the specific parameter name. Also added multimodel type, which allows for multi-task style learning of existing models. Added linear mapping type of model. 

Version 0.1291
--------------

Changes to modelOutputGrad.m, modelOut.m, kbrOutputGrad.m, kbrExpandParam.m, modelOptimise.m to allow compatibility with SGPLVM and NCCA toolboxes. Added a preliminary coding of LLE.


Version 0.129
-------------

Added dnet type model for GTM and density networks. Added various lvm helper files for doing nearest neighbour and plotting results for latent variable models. Added lmvu and mvu embedding wrapper. Added ppca model type. Added output gradients for model out functions (for magnification factor computation in dnet models). Added helpers for reading various models from FID mapmodel, matrix etc.).
Added rbfOutputGradX and visualisation for spring dampers type.

Version 0.128
-------------

Fixed Viterbi alignment algorithm, thanks to Raquel Urtasun for pointing out the problems with it.

Carl Henrik Ek added embeddings with maximum variance unfolding (landmark and normal) to the toolbox. Also several files modified by Carl to allow a single output dimension of a model to be manipulated.

Version 0.127
-------------

Minor modifications including adding file modelAddDynamics to replace fgplvmAddDynamics.

Version 0.126
-------------

Modified kbrParamInit to scale alpha weights and biases by number of data. Added 'dynamicsSliderChange' to lvmClassVisualise to allow visulisation of models with 'gpTime' style-dynamics.

Version 0.125
-------------

Added multimodel for learning multiple indepedent models with shared parameters.

Version 0.124
-------------

Added model gradient checker and added rbfperiodic function to provide a length scale for the gibbsperiodic kernel.

Version 0.123
-------------

Minor release in line with IVM toolbox 0.4.

Version 0.122
-------------

Added Hessian code for base model type and for MLP. Added Viterbi alignment code, viterbiAlign.

Version 0.121
-------------

Various minor bug fixes and changes which seem to have gone undocumented.

Version 0.12
------------

Extended model type to be a generic container module for optimising any model. Added model test for testing a created model. The code is still in a bit of flux though with some design decisions not made and some code untested.

Version 0.111
-------------

Fixed bug in kbr where bias parameter fields where still being referred to as b.Also acknowledged the fact that the isomap graph may not be fully connected in isomapEmbed, but don't yet deal with it properly. Finally added lleEmbed.m for wrapping the lle code.


Version 0.11
------------

Updated release for operation with FGPLVM toolbox 0.13. Structure of model creation changed and functions of the form modelOptions.m included for setting default options of the various machine learning models.

Version 0.1
-----------

The first release of the toolbox with various wrappers for NETLAB functions. Also latent variable model visualisation code was moved into this toolbox.


MATLAB Files
------------

Matlab files associated with the toolbox are:

dnetExpandParam.m: Update dnet model with new vector of parameters.
modelPointLogLikelihood.m: Compute the log likelihood of a given point.
isomapOptions.m: Options for a isomap.
modelWriteResult.m: Write a model to file.
demOilLle3.m: Demonstrate LLE on the oil data.
dnetUpdateOutputWeights.m: Do an M-step (update parameters) on an Density Network model.
mlpLogLikelihood.m: Multi-layer perceptron log likelihood.
rbfCreate.m: Wrapper for NETLAB's rbf `net'.
spectralUpdateLaplacian.m: Update the Laplacian using graph connections.
mogEstep.m: Do an E-step on an MOG model.
dnetObjective.m: Wrapper function for Density Network objective.
kbrParamInit.m: KBR model parameter initialisation.
smallrandEmbed.m: Embed data set with small random values.
vectorModify.m: Helper code for visualisation of vectorial data.
modelOptions.m: Returns a default options structure for the given model.
mogProject.m: Project a mixture of Gaussians to a low dimensional space.
ppcaOut.m: Output of an PPCA model.
modelLogLikelihood.m: Compute a model log likelihood.
modelWriteToFID.m: Write to a stream a given model.
dnetReconstruct.m: Reconstruct an DNET form component parts.
modelPosteriorVar.m: variances of the posterior at points given by X.
imageVisualise.m: Helper code for showing an image during 2-D visualisation.
mogUpdateCovariance.m: Update the covariances of an MOG model.
modelOutputGrad.m: Compute derivatives with respect to params of model outputs.
rbfperiodicOutputGrad.m: Evaluate derivatives of RBFPERIODIC model outputs with respect to parameters.
demSwissRollLle3.m: Demonstrate LLE on the oil data.
linearCreate.m: Create a linear model.
modelOutputGradX.m: Compute derivatives with respect to model inputs of model outputs.
lvmScatterPlotNeighbours.m: 2-D scatter plot of the latent points with neighbourhood.
leOptimise.m: Optimise an LE model.
kpcaEmbed.m: Embed data set with kernel PCA.
lleReconstruct.m: Reconstruct an LLE form component parts.
distanceWarp.m: Dynamic Time Warping Algorithm
plot3Modify.m: Helper code for visualisation of 3-d data.
kbrDisplay.m: Display parameters of the KBR model.
kbrExpandParam.m: Create model structure from KBR model's parameters.
springDampersModify.m: Helper code for visualisation of springDamper data.
mvuCreate.m: Maximum variance unfolding embedding model.
multimodelParamInit.m: MULTIMODEL model parameter initialisation.
linearDisplay.m: Display a linear model.
mogLowerBound.m: Computes lower bound on log likelihood for an MOG model.
lvmClassVisualise.m: Callback function for visualising data.
rbfOptions.m: Default options for RBF network.
isomapDeconstruct.m: break isomap in pieces for saving.
dnetWriteResult.m: Write a DNET result.
modelExpandParam.m: Update a model structure with parameters.
dnetDeconstruct.m: break DNET in pieces for saving.
dnetOptimise.m: Optimise an DNET model.
lvmVisualise.m: Visualise the manifold.
dnetLogLikelihood.m: Density network log likelihood.
multimodelOptions.m: Create a default options structure for the MULTIMODEL model.
mogUpdatePrior.m: Update the priors of an MOG model.
leOptions.m: Options for a Laplacian eigenmaps.
mogSample.m: Sample from a mixture of Gaussians model.
modelGradientCheck.m: Check gradients of given model.
mlpCreate.m: Multi-layer peceptron model.
kbrExtractParam.m: Extract parameters from the KBR model structure.
ppcaDeconstruct.m: break PPCA in pieces for saving.
modelPosteriorMeanVar.m: Mean and variances of the posterior at points given by X.
mlpLogLikeGradients.m: Multi-layer perceptron gradients.
lleOptimise.m: Optimise an LLE model.
lvmResultsDynamic.m: Load a results file and visualise them.
ppcaCreate.m: Density network model.
modelGradient.m: Gradient of error function to minimise for given model.
mogUpdateMean.m: Update the means of an MOG model.
spectrumVisualise.m: Helper code for showing an spectrum during 2-D visualisation.
mvuOptimise.m: Optimise an MVU model.
linearLogLikeGradients.m: Linear model gradients.
doubleMatrixReadFromFID.m: Read a full matrix from an FID.
modelLoadResult.m: Load a previously saved result.
kbrCreate.m: Create a KBR model.
mlpOutputGrad.m: Evaluate derivatives of mlp model outputs with respect to parameters.
dnetTest.m: Test some settings for the density network.
modelSetOutputWeights.m: Wrapper function to return set output weight and bias matrices.
modelGetOutputWeights.m: Wrapper function to return output weight and bias matrices.
ppcaPosteriorVar.m: Mean and variances of the posterior at points given by X.
lvmScatterPlotColor.m: 2-D scatter plot of the latent points with color.
lvmScatterPlot.m: 2-D scatter plot of the latent points.
lvmClickVisualise.m: Visualise the manifold using clicks.
modelHessian.m: Hessian of error function to minimise for given model.
demSwissRollFullLle1.m: Demonstrate LLE on the oil data.
dnetUpdateBeta.m: Do an M-step (update parameters) on an Density Network model.
modelSamp.m: Give a sample from a model for given X.
matrixReadFromFID.m: Read a matrix from an FID.
multimodelExpandParam.m: Create model structure from MULTIMODEL model's parameters.
rbfperiodicLogLikelihood.m: Log likelihood of RBFPERIODIC model.
demSwissRollLle4.m: Demonstrate LLE on the oil data.
lmvuEmbed.m: Embed data set with landmark MVU
rbfExpandParam.m: Update rbf model with new vector of parameters.
mlpOut.m: Output of an MLP model.
doubleMatrixWriteToFID.m: Writes a double matrix to an FID.
lfmVisualise.m: Visualise the outputs in a latent force model
rbfOutputGrad.m: Evaluate derivatives of rbf model outputs with respect to parameters.
rbfperiodicDisplay.m: Display parameters of the RBFPERIODIC model.
swissRollScatter.m: 3-D scatter plot with colors.
dnetOptions.m: Options for a density network.
ppcaReconstruct.m: Reconstruct an PPCA form component parts.
lvmLoadResult.m: Load a previously saved result.
multimodelExtractParam.m: Extract parameters from the MULTIMODEL model structure.
multimodelLogLikelihood.m: Log likelihood of MULTIMODEL model.
linearLogLikelihood.m: Linear model log likelihood.
modelLogLikeGradients.m: Compute a model's gradients wrt log likelihood.
linearOutputGrad.m: Evaluate derivatives of linear model outputs with respect to parameters.
mltoolsToolboxes.m: Load in the relevant toolboxes for the MLTOOLS.
mlpDisplay.m: Display the multi-layer perceptron model.
isomapOptimise.m: Optimise an ISOMAP model.
isomapReconstruct.m: Reconstruct an isomap form component parts.
linearOutputGradX.m: Evaluate derivatives of linear model outputs with respect to inputs.
lleDeconstruct.m: break LLE in pieces for saving.
mvuDeconstruct.m: break MVU in pieces for saving.
springDampersVisualise.m: Helper code for showing an spring dampers during 2-D visualisation.
lfmResultsDynamic.m: Load a results file and visualise them.
isomapCreate.m: isomap embedding model.
demOilLle2.m: Demonstrate LLE on the oil data.
lfmClassVisualise.m: Callback function to visualize LFM in 2D
kbrOptimise.m: Optimise a KBR model.
linearExtractParam.m: Extract weights from a linear model.
lvmThreeDPlot.m: Helper function for plotting the labels in 3-D.
demSwissRollFullLle4.m: Demonstrate LLE on the oil data.
demSwissRollLle2.m: Demonstrate LLE on the oil data.
mogCreate.m: Create a mixtures of Gaussians model.
mlpExpandParam.m: Update mlp model with new vector of parameters.
lleOptions.m: Options for a locally linear embedding.
lvmNearestNeighbour.m: Give the number of errors in latent space for 1 nearest neighbour.
rbfOptimise.m: Optimise RBF for given inputs and outputs.
rbfperiodicOutputGradX.m: Evaluate derivatives of a RBFPERIODIC model's output with respect to inputs.
demSwissRollFullLle5.m: Demonstrate LLE on the oil data.
dnetPosteriorMeanVar.m: Mean and variances of the posterior at points given by X.
dnetOutputGradX.m: Evaluate derivatives of DNET model outputs with respect to inputs.
findAcyclicNeighbours.m: find the k nearest neighbours for each point in Y preventing cycles in the graph.
modelObjective.m: Objective function to minimise for given model.
paramNameReverseLookup.m: Returns the index of the parameter with the given name.
kbrOptions.m: Create a default options structure for the KBR model.
linearOut.m: Obtain the output of the linear model.
dnetLoadResult.m: Load a previously saved result.
dnetOut.m: Output of an DNET model.
lvmTwoDPlot.m: Helper function for plotting the labels in 2-D.
mlpOptions.m: Options for the multi-layered perceptron.
mapmodelReadFromFID.m: Load from a FID produced by C++ code.
mvuOptions.m: Options for a MVU.
dnetLowerBound.m: Computes lower bound on log likelihood for an DNET model.
rbfperiodicCreate.m: Create a RBFPERIODIC model.
linearParamInit.m: Initialise the parameters of an LINEAR model.
mlpOptimise.m: Optimise MLP for given inputs and outputs.
leDeconstruct.m: break LE in pieces for saving.
dnetLogLikeGradients.m: Density network gradients.
modelDisplay.m: Display a text output of a model.
rbfExtractParam.m: Wrapper for NETLAB's rbfpak.
rbfperiodicExpandParam.m: Create model structure from RBFPERIODIC model's parameters.
mvuReconstruct.m: Reconstruct an MVU form component parts.
rbfOutputGradX.m: Evaluate derivatives of a RBF model's output with respect to inputs.
findDirectedNeighbours.m: find the k nearest neighbours for each point in Y preventing cycles in the graph.
rbfperiodicLogLikeGradients.m: Gradient of RBFPERIODIC model log likelihood with respect to parameters.
lvmClassVisualisePath.m: Latent variable model path drawing in latent space.
lvmSetPlot.m: Sets up the plot for visualization of the latent space.
modelParamInit.m: Initialise the parameters of the model.
rbfDisplay.m: Display an RBF network.
mogOptimise.m: Optimise an MOG model.
modelCreate.m: Create a model of the specified type.
lvmResultsClick.m: Load a results file and visualise them with clicks
lleEmbed.m: Embed data set with LLE.
dnetGradient.m: Density Network gradient wrapper.
modelExtractParam.m: Extract the parameters of a model.
kbrOutputGrad.m: Evaluate derivatives of KBR model outputs with respect to parameters.
mogLogLikelihood.m: Mixture of Gaussian's log likelihood.
modelAddDynamics.m: Add a dynamics kernel to the model.
mappingOptimise.m: Optimise the given model.
rbfperiodicOut.m: Compute the output of a RBFPERIODIC model given the structure and input X.
rbfOut.m: Output of an RBF model.
kbrOut.m: Compute the output of a KBR model given the structure and input X.
linearOptimise.m: Optimise a linear model.
mogPrintPlot.m: Print projection of MOG into two dimensions.
linearOptions.m: Options for learning a linear model.
demOilLle4.m: Demonstrate LLE on the oil data.
multimodelDisplay.m: Display parameters of the MULTIMODEL model.
isomapEmbed.m: Embed data set with Isomap.
modelReadFromFID.m: Load from a FID produced by C++ code.
linearExpandParam.m: Update linear model with vector of parameters.
rbfperiodicOptions.m: Create a default options structure for the RBFPERIODIC model.
spectrumModify.m: Helper code for visualisation of spectrum data.
mlpOutputGradX.m: Evaluate derivatives of mlp model outputs with respect to inputs.
vectorVisualise.m:  Helper code for plotting a vector during 2-D visualisation.
mlpParamInit.m: Initialise the parameters of an MLP model.
findAcyclicNeighbours2.m: find the k nearest neighbours for each point in Y preventing cycles in the graph.
mlpExtractParam.m: Extract weights and biases from an MLP.
demOilLle1.m: Demonstrate LLE on the oil data.
mogMeanCov.m: Project a mixture of Gaussians to a low dimensional space.
paramNameRegularExpressionLookup.m: Returns the indices of the parameter containing the given regular expression.
demSwissRollFullLle2.m: Demonstrate LLE on the oil data.
demSwissRollFullLle3.m: Demonstrate LLE on the oil data.
viterbiAlign.m: Compute the Viterbi alignment.
dnetCreate.m: Density network model.
leCreate.m: Laplacian eigenmap model.
demSwissRollLle1.m: Demonstrate LLE on the oil data.
multimodelLogLikeGradients.m: Gradient of MULTIMODEL model log likelihood with respect to parameters.
mvuEmbed.m: Embed data set with MVU.
modelOut.m: Give the output of a model for given X.
lvmClassClickVisualise.m: Callback function for visualising data in 2-D with clicks.
lvmScoreModel.m: Score model with a GP log likelihood.
demMppca1.m: Demonstrate MPPCA on a artificial dataset.
mogTwoDPlot.m: Helper function for plotting the labels in 2-D.
modelOptimise.m: Optimise the given model.
ppcaOptions.m: Options for probabilistic PCA.
spectralUpdateX.m: Update the latent representation for spectral model.
dnetExtractParam.m: Extract weights and biases from an DNET.
modelTest.m: Run some tests on the specified model.
mlpLogLikeHessian.m: Multi-layer perceptron Hessian.
dnetEstep.m: Do an E-step (update importance weights) on an Density Network model.
modelReadFromFile.m: Read model from a file FID produced by the C++ implementation.
imageModify.m: Helper code for visualisation of image data.
modelTieParam.m: Tie parameters of a model together.
dnetOutputGrad.m: Evaluate derivatives of dnet model outputs with respect to parameters.
ppcaEmbed.m: Embed data set with probabilistic PCA.
rbfperiodicExtractParam.m: Extract parameters from the RBFPERIODIC model structure.
leReconstruct.m: Reconstruct an LE form component parts.
multimodelCreate.m: Create a MULTIMODEL model.
plot3Visualise.m:  Helper code for plotting a plot3 visualisation.
findNeighbours.m: find the k nearest neighbours for each point in Y.
rbfperiodicParamInit.m: RBFPERIODIC model parameter initialisation.
lvmPrintPlot.m: Print latent space for learnt model.
ppcaPosteriorMeanVar.m: Mean and variances of the posterior at points given by X.
lleCreate.m: Locally linear embedding model.
mogOptions.m: Sets the default options structure for MOG models.
demSwissRollLle5.m: Demonstrate LLE on the oil data.
