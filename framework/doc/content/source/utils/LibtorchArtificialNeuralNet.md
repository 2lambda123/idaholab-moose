# LibtorchArtificialNeuralNet

## Overview

This class can be used to generate a simple, feedworward artificial neural network (ANN)
using the underlying objects imported from `libtorch` (C++ API of `pytorch`). For
the instructions on how to enable `libtorch` support in MOOSE, visit [installation/install_libtorch.md].
For a more detailed introduction to neural networks, we refer the reader to [!cite](muller1995neural).
The architecture of a simple feedforward neural network is presented below. The first layer
from the left the right are referred to as input and output layers,
while the layers between them are the hidden layers.

!media stochastic_tools/surrogates/libtorch_nn/neural_net.png style=width:65%; id=structure
      caption=The architecture of the simple feed-forward neural network in MOOSE-STM.


We see that the outputs ($\textbf{y}$) of the neural net can be expressed as function of the
inputs ($\textbf{x}$) and the corresponding model parameters (weights $w_{i,j}$, organized in
weight matrics $\textbf{W}$ and biases $b_i$ organized in the bias vector $\textbf{b}$)
in the following nested form:

!equation id=nn-explicit
\textbf{y} = \sigma(\textbf{W}^{(3)}\sigma(\textbf{W}^{(2)}\sigma(\textbf{W}^{(1)}\textbf{x}+\textbf{b}^{(1)})+\textbf{b}^{(2)})+\textbf{b}^{(3)}),

where $\sigma$ denotes the activation function. At the moment, the Moose implementation
supports `relu`, `elu`, `gelu`, `sigmoid` and `linear` activation functions.
 In this class, no activation function is applied on the
output layer. It is apparent that the real functional dependence between the inputs and outputs
is approximated by the function in [nn-explicit]. As in most cases, the error in this approximation depends on the
smoothness of the function and the values of the model parameters. The weights and
biases in the function are determined by minimizing the error between the
approximate outputs of the neural net corresponding reference (training) values
over a training set.

## Example usage

To be able to use this neural network, we have to construct one using a name,
the number of expected inputs and outputs, an expected hideen-layer-structure and
the activation functions for the layers. If no activation function is given,
`relu` is used for every hidden neuron:

!listing test/src/utils/LibtorchArtificialNeuralNetTest.C start=Define neurons end=activation_functions include-end=true

For training a neural network, we need to initialize an optimizer (ADAM in this case),
then supply a known input-output combinations for the function-to-be-approximated
and let the optimizer set the parameters of the neural network to ensure that the
answer supplied by the neural network is as close to the supplied values as possible.
Once step in this optimization process is shown below:

!listing test/src/utils/LibtorchArtificialNeuralNetTest.C start=Create an Adam optimizer end=optimizer.step(); include-end=true


For more detailed instructions on training a neural network, visit the Stochastic Tools module!
