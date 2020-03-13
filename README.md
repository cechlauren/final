# Write Up Final: Constructing a Neural Network

[![Build
Status](https://travis-ci.org/cechlauren/HW3_skeleton.svg?branch=master)](https://travis-ci.org/cechlauren/HW3_skeleton)

Contains: autoencoder and predictions for RAP1 binding sites.

## STRUCTURE
```
.
├── README.md
├── requirements.txt
├── various .png
│   ...
├── NNfxns
│   ├── neural_network.py
│   ├── auto.py
│   └── train.py
├── test
│    └── test_auto.py
└── Jupyter notebook
    ├── 8x8 reconstruction
    └── other workflow
```

Any dependencies noted in `requirements.txt`. 

## Building an autoencoder:
## Background
Autoencoders are a type of feedforward neural networks where the input is the same thing as the output. Basically, they will take the input and squish it into a low dimension summary code that represents what should be reconstructed in the output. To make an autoencoder (and later on a neural network), one needs at least three things: (1) an encoding method, (2) a decoding method, and (3) a loss function. The encoder and decoder are fully connected feedforward NN, while the code is a single layer of artificial NN that has a hyperparameter-defined number of nodes. Additional hyperparameters that should be set before training are: (1) code size, (2) number of layers, (3) number of nodes per layer, and (4) a loss function like MSE or binary crossentropy. Some of those hyperparameters were already set for us in this final…

## Implementation: reconstructing an 8x8 identity matrix with an autoencoder containing 3 hidden neurons
The 8x3x8 autoencoder used specifically to process the 8bit data can be found here [auto.py](https://github.com/cechlauren/final/blob/master/NNfxns/auto.py).

The input is represented here for ease of testing by others [eightBit.txt](https://github.com/cechlauren/final/blob/master/data/eightBit.txt).

The function that produces the reconstruction is here [test_auto.py](https://github.com/cechlauren/final/blob/master/test/test_auto.py). 


The input for this autoencoder was as an 8x8 identity matrix made using 
<img src="8by8code.png" /><br />

which makes this

<img src="8by8.png" /><br />

If you see nothing, see:
[8by8code.png](https://github.com/cechlauren/final/blob/master/8by8code.png) and
[8by8.png](https://github.com/cechlauren/final/blob/master/8by8.png).

Then the autoencoder vectorizes (ie. reduces the dimension of) the input in a bitwise fashion for the hidden layer(pretty much all of the work is done for the autoencoder already!). 

Then that hidden, reduced encoding is processed by the decoding layer to give an output like the following:

<img src="test_auto.png" /><br />

Also here: [test_auto.png](https://github.com/cechlauren/final/blob/master/test_auto.png) or viewable in Jupyter notebook.

Its not pretty, but this seems to do pretty well given that I've added quite a bit of noise to the autoencoder.
There seems to be some difference in learning activation for the central section, but I'd say its decent enough since this part doesn't need to be perfect.

## Develop a fully connected neural network that predicts transcription factor binding with the training data provided

### Data Preprocessing: Describe and implement a data preprocessing approach
How one prepares the NN data is probably one of the more imporant parts to getting reasonable predictions; "crap in, crap out" as they say. 

There are 4 different nucleotides potentially describing each base in the given sequences. Each nucleotide can be represented in a bitwise fashion using one-hot encoding. For example, from this code: [neural_network.py](https://github.com/cechlauren/final/blob/master/NNfxns/neural_network.py)
```
self.base_binary_conversion = {'A': '0001',
                                       'C': '0010',
                                       'T': '0100',
                                       'G': '1000'
                                       }
```

Not only is this a binary representation, but can be adjusted for wildcard nucleotides, or we can change how we assign each nucleotide if we think there may be a bias for certain nucleotides in the training set (we don't want certain features to overpower the others just because there happens to be randomly more 'G' rich sequences in the training set, for instance.




### DNA Sequence Representation: Describe and implement a way to represent DNA sequence
### Network Architecture: Develop and describe your network architecture
To make something more powerful than the 8x3x8 autoencoder, one can increase the number of layers, nodes per layer, and more importantly the code size. Increasing those hyperparameters allows neural networks to learn more complex codings. Overdoing this part, however, could cause overfitting since the NN will simply learn to copy the inputs as the output without learning anything meaningful. By making more of a sandwich where the code size is small, the NN won’t be able to directly copy the input to the output and so is forced to learn representative features. 
Sometimes we can force an NN to learn useful features by adding random noise so that the NN must determine the meaningful data. 


## Training Regime: Develop a training regime (K-fold cross validation, bagging, etc) to test model 
### Describe and implement the regime and answer question 3 subquestions 

## Cross-Validation: Perform cross-validation experiments to test model hyperparameters
### Develop and describe your choice of model hyperparameters
### Question 4

## Testing: Test model performance on test data



Ref and citations:

