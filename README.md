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

There are ~137 positive sequences each with 17 nucleotides that describe RAP1 binding sites. Its possible that our negatives, which are ~1000 nucleotides long, contain at least partial positives. Thus, we may see some terrible predictions if our negatives aren't really negatives. 

One way to do that might be breaking up the negatives into 17 nucleotide chunks and comparing it with the positive set over and over (ie. recycle them) until one runs through the all the distinct negatives. See line 57 in [train.py](https://github.com/cechlauren/final/blob/master/NNfxns/train.py).

```
partitions = int(args['<partitions>']) #partitions: the number of different parts we should make for cross-validation
    neg_division = len(negative_sites) / float(partitions) # how many partitions we can make from the sites given
    neg_randomly_partitioned_list = [negative_sites[int(round(neg_division * i)): int(round(neg_division * (i + 1)))]
                                     for i in range(partitions)] # makes a list of those partitions thru each site

    pos_division = len(positive_sites) / float(partitions) # ditto ^^
    pos_randomly_partitioned_list = [positive_sites[int(round(pos_division * i)): int(round(pos_division * (i + 1)))]
                                     for i in range(partitions)]

```

One will import the RAP1 binding sites and the negative sites (parsing the fasta for the 1000 bases).
Then for each sequence, iterate over each 17 base chunk and train with an expected of 0 for neg and 1 for pos. 

Clearly, the number of negative instances is going to be much more than that of positives. 

Some ways to deal with this might be:
- increasing our positives (small effect overall)
- making the positive:negative ratio 1:1 (large effect overall)

The latter point supports a decision to train against all positive training sites for every 137 negative training instances.




### DNA Sequence Representation: Describe and implement a way to represent DNA sequence

There are 4 different nucleotides potentially describing each base in the given sequences. Each nucleotide can be represented in a bitwise fashion using one-hot encoding. For example, from this code: [neural_network.py](https://github.com/cechlauren/final/blob/master/NNfxns/neural_network.py)
```
self.base_binary_conversion = {'A': '0001',
                                       'C': '0010',
                                       'T': '0100',
                                       'G': '1000'
                                       }
```

Not only is this a binary (read less memory) representation, but can be adjusted for wildcard nucleotides, or we can change how we assign each nucleotide if we think there may be a bias for certain nucleotides in the training set (we don't want certain features to overpower the others just because there happens to be randomly more 'G' rich sequences in the training set, for instance.

For this NN, the input_DNA string 'GA' would be vectorized like this:

```
self.base_binary_conversion = {'A': '0001',
                                       'C': '0010',
                                       'T': '0100',
                                       'G': '1000'
                                       }
...                                       
                                       
def _construct_input_vector(self, input_DNA):
       
        temp_vector_list = []

        for base in input_DNA:
            for number in self.base_binary_conversion[base]:
                temp_vector_list.append(float(number))

        self.input_vector = np.asarray(temp_vector_list)

```
This may end up looking like 
```
[10000 
0001]

```

or more likely
```
[.99910 0.00001 0.00001 0.00001
 0.00001 0.00001 0.00001 .99910]

```



### Network Architecture: Develop and describe your network architecture
To make something more powerful than the 8x3x8 autoencoder, one can increase the number of layers, nodes per layer, and more importantly the code size. Increasing those hyperparameters allows neural networks to learn more complex codings. Overdoing this part, however, could cause overfitting since the NN will simply learn to copy the inputs as the output without learning anything meaningful. By making more of a sandwich where the code size is small, the NN won’t be able to directly copy the input to the output and so is forced to learn representative features. 
Sometimes we can force an NN to learn useful features by adding random noise so that the NN must determine the meaningful data. 

In this instance we have :
* Input layer with 17*4 nodes (because 17 different nucleotides defining each sequence, and 4 different
          possible nucleotides describing those positions) + bias 
* Hidden layer with 23-35 nodes (merge at least 2 input neurons)+ bias   
* One output layer node (number of neurons in the output layer will equal the number of outputs associated with each input; we want one answer)
* All with sigmoid activation

The network will accept 17 base units and will output a value from around 0 to around 1, where 0 is not a RAP1 binding site.
 

## Training Regime: Develop a training regime (K-fold cross validation, bagging, etc) to test model 
### Describe and implement the regime
### Answer question 3 subquestions 

- How was your training regime designed so as to prevent the negative training data from overwhelming the positive training data?

See answer above in Data Preprocessing.
But in short, I resampled from the positives, while using all 17 nt chunks possible from each negative instance.
May have been interesting to train on negative examples that had higher homology, rather than just pooling all of those negative instances together...that way the NN would have a better chance at learning the features that separate RAP1 sites.


- What was your stop criterion for convergence in your learned parameters? How did you decide this?
In most instances, there are three distinct splits of data. In one we train, in the others we test and validate. When we train, we train across epochs and evaluate loss. If we visualize this, it produces what's called a "loss curve." 
Sometimes we do early stopping if things dont change towards the end of the curve. Like if it looks flat and doesn't change. 
On the test curve we see something similar, but the curve will start to go back up after a certain amount of epochs. 
That's what overfitting looks like. Ideally, we'd stop training at the intersect of the training/test curves. Like this: 

how do you know when you're done training? Three splits of data (at least two) some segment we train on and some that we test/validate on. When we train, we train across epochs and evaluate our loss. Our traibned set we see a loss curve. Sometimes we do early so=toping if things dont change towards the end of the curve, On the test curve we see someting similar but the curve starts go back up, so this is a result of overfitting. Stop training at the intersect. 
We can just train until we think its appropriate, "loss curves"

resample from positives (sample like 3 times), downsample from negs (helps prevent overfitting),  want 1:1. If there was more time, then you'd look at positives and see how many the nucleotides show up. 


## Cross-Validation: Perform cross-validation experiments to test model hyperparameters
### Develop and describe your choice of model hyperparameters

### Question 4

## Testing: Test model performance on test data
* For alpha = 1: [NN_predictions.txt](https://github.com/cechlauren/final/blob/master/NN_predictions.txt)
* For alpha = 



Ref and citations:

