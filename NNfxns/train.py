"""
Training the neural net

What can be used:
     autoencoder
     testme <splits> <sampling>
     test
Arguments:
    autoencoder
        Run ze autoencoder
    testme
        Do the rap1 learning task including cross-validation
    test
        Classify test data and output to tsv file type
    <splits>
        Number of splits to make for cross-valitation (k-fold!)
    <sampling>
        Sampling method for neuralnet training input
        (slide) Go over each sequence in 17 nucleotide sliding frame (the length of the positives/test binding sites)
        (space) Each sequence is cut into 17 nucleotide bits for inputs (the binary bits that are useable by our model)
"""

def testme():
    """
    Train neural network on RAP1 binding sites
        * Input layer with 17*4 nodes (because 17 different nucleotides defining each sequence, and 4 different
          possible nucleotides describing those positions) + bias
        * Hidden layer with 23-35 nodes (merge at least 2 input neurons)+ bias
        * One output layer node (number of neurons in the output layer will equal the number of outputs
          associated with each input; we want one answer)
    Train against negative and positive binding sites
        * Import all negative sequences from .fa file
        * For each sequence, iterate every 17 bases and train with
          expected of 0 (since these are negatives, while positives should be 1)
        * Because it is so important to have a 1:1 ratio of negatives:positives when training our model, 
          for every 137 negative training instances, need to train it against all positive binding sites
          with expected of 1 
        * Will continue until negative binding sites have been fully ran through
    """

    # This part takes care of bringing in all those positive sequences
    
    pos = [pos_sequence.strip() for pos_sequence in open('data/rap1-lieb-positives.txt')]

    
    
    # This part takes care of bringing in all those neg sequences
    
    neg = list(SeqIO.parse('data/yeast-upstream-1k-negative.fa', 'fasta'))
    
    

    # Separate into random sections for the k-fold x-validation for both our positives and negatives 
    # Taken from : http://stackoverflow.com/questions/3352737/python-randomly-partition-a-list-into-n-nearly-equal-parts
    
    
    splits = int(args['<splits>'])
    negative_division = len(neg) / float(splits) # how many splits we can make from the sites given
    neg_split_up = [neg[int(round(negative_division * i)): int(round(negative_division * (i + 1)))]
                                     for i in range(splits)] # makes a list of those splits thru each site

    pos_division = len(pos) / float(splits) # ditto ^^
    pos_split_up = [pos[int(round(pos_division * i)): int(round(pos_division * (i + 1)))]
                                     for i in range(splits)]

    
    
    # Go thru neg sites subsets for x-validation, keep track of how many separations we do based on splits
    
    separation = 0
    for index in range(int(args['<splits>'])):
        # Set up cross-validation sets for the positives and negatives
        neg_site_list_copy = copy.deepcopy(neg_split_up)
        del neg_site_list_copy[index]
        neg_site_training = [seq for partition in neg_site_list_copy for seq in partition]
        neg_cross_validation_set = neg_split_up[index]

        pos_site_list_copy = copy.deepcopy(pos_split_up)
        del pos_site_list_copy[index]
        pos_site_training = [seq for partition in pos_site_list_copy for seq in partition]
        pos_cross_validation_set = pos_split_up[index]

        print("Training on the training set...")

        # Input our hyperparameters = # nodes
        neuralnet = neural_network(68, 23, 1)

        # See neural_net.py to get info on initialization
        neuralnet.initialize_values()

        pos_counter = 0
        counter = 0

        # If we're sampling from tneg then we'll slide over 17 nucleotides 
        
        if args['<sampling>'] == 'slide':
            for site in neg_site_training:

                # Iterate over site in 17 nucleotide sliding frames in negative sites, decide which model to use
                for chunky in range(len(site) - 16):
                    slice = site[chunky:(chunky + 17)].seq
                    if slice not in pos:
                        if all([slice[4] == 'C', slice[5] == 'C', slice[9] == 'C']) == False:
                            neuralnet.setin_n_exp_values(slice, autoencoder=False, negative=True)
                            neuralnet.forward_propogation()
                            neuralnet.backward_propogation()
                            neuralnet.update_weights_and_bias()
                            pos_counter += 1
                        else:
                            print(slice)

                    if pos_counter == len(pos_site_training):
                        for pos_site in pos_site_training:
                            neuralnet.setin_n_exp_values(pos_site, autoencoder=False, negative=False)
                            neuralnet.forward_propogation()
                            neuralnet.backward_propogation()
                            neuralnet.update_weights_and_bias()

                        pos_counter = 0

                # have reset the positives counter and will now say that we've done some training on those
                
                counter += 1

                print("Training set: {}/{} completed...".format(counter, len(neg_cross_validation_set)))

                greatestdelta_1 = neuralnet.matrix_1_errors.max()
                smallestdelta_1 = neuralnet.matrix_1_errors.min()
                greatestdelta_2 = neuralnet.matrix_2_errors.max()
                smallestdelta_2 = neuralnet.matrix_2_errors.min()

                if any([greatestdelta_1 < 0.00000000001 and greatestdelta_1 > 0,
                        smallestdelta_1 > -.00000000001 and smallestdelta_1 < 0]) and any(
                    [greatestdelta_2 < 0.00000000001 and greatestdelta_2 > 0,
                     smallestdelta_2 > -0.00000000001 and smallestdelta_2 < 0]):
                    print("Stop criterion met after {} iterations".format(counter))
                    break

        #when we sample from the negatives we only take 17 nucleotide chunks from each site
        
        if args['<sampling>'] == 'space':
            for site in neg_site_training:
                
                number_of_chunkys = int(len(site) / 17) #length of neg site tells us the amount of 17 length chunks possible

                for chunky in range(number_of_chunkys):
                    slice = site[(chunky * 17):((chunky + 1) * 17)].seq
                    if slice not in pos:
                        if all([slice[4] == 'C', slice[5] == 'C', slice[9] == 'C']) == False:
                            neuralnet.setin_n_exp_values(slice, autoencoder=False, negative=True)
                            neuralnet.forward_propogation()
                            neuralnet.backward_propogation()
                            neuralnet.update_weights_and_bias()
                            pos_counter += 1

                        else:
                            print(slice)

                    #quick check to make sure that we've finished going thru the positives yet
                    if pos_counter == len(pos_site_training):
                        for pos_site in pos_site_training:
                            neuralnet.setin_n_exp_values(pos_site, autoencoder=False, negative=False)
                            neuralnet.forward_propogation()
                            neuralnet.backward_propogation()
                            neuralnet.update_weights_and_bias()

                        pos_counter = 0

                    counter += 1

                greatestdelta_1 = neuralnet.matrix_1_errors.max()
                smallestdelta_1 = neuralnet.matrix_1_errors.min()
                greatestdelta_2 = neuralnet.matrix_2_errors.max()
                smallestdelta_2 = neuralnet.matrix_2_errors.min()

                if any([greatestdelta_1 < 0.00000000001 and greatestdelta_1 > 0,
                        smallestdelta_1 > -.00000000001 and smallestdelta_1 < 0]) and any(
                    [greatestdelta_2 < 0.00000000001 and greatestdelta_2 > 0,
                     smallestdelta_2 > -0.00000000001 and smallestdelta_2 < 0]):
                    print("Stop criterion met after {} iterations".format(counter))
                    break

        # taken each partition and trained model
        print("Performing Cross-validation")

        pos_list = []
        neg_list = []

        
        # Return the sets of positives and negatives from the x-validation
        
        print("Negative cross-validation set...")
        counter = 0
        for site in neg_cross_validation_set:
            for slice in range(len(site) - 16):
                neuralnet.setin_n_exp_values(site[slice:slice + 17].seq, autoencoder=False, negative=True)
                neuralnet.forward_propogation()
                neg_list.append(neuralnet.output_layer_output)
            counter += 1
            print("Negative cross-validation: {}/{} completed...".format(counter, len(neg_cross_validation_set)))
            break

        print("Positive cross-validation set...")
        for site in pos_cross_validation_set:
            neuralnet.setin_n_exp_values(site, autoencoder=False)
            neuralnet.forward_propogation()
            pos_list.append(neuralnet.output_layer_output)

        print('Positive avg: {}'.format(sum(pos_list) / len(pos_list)))
        print('Negative avg: {}'.format(sum(neg_list) / len(neg_list)))
        print(neuralnet.matrix_1_bias)
        print(neuralnet.matrix_2_bias)

        # Output the coneuralnetection matrices with greatest separation between the average positive and negative scores
        if ((sum(pos_list) / len(pos_list)) - (sum(neg_list) / len(neg_list))) > separation:
            np.savetxt('cnx_matrix_1.csv', neuralnet.matrix_1_bias, delimiter=',')
            np.savetxt('cnx_matrix_2.csv', neuralnet.matrix_2_bias, delimiter=',')
            separation = (sum(pos_list) / len(pos_list)) - (sum(neg_list) / len(neg_list))


# A simple definition of the autoencoder that uses those same neuralnet parameters
def autoencoder():
  
    neuralnet = neural_network()
    neuralnet.setin_n_exp_values('GA', autoencoder=True)
    neuralnet.initialize_values()

    # Stop criterion
    finished_working = False

    while finished_working == False:
        neuralnet.forward_propogation()
        neuralnet.backward_propogation()
        neuralnet.update_weights_and_bias()

        greatestdelta_1 = neuralnet.matrix_1_errors.max()
        smallestdelta_1 = neuralnet.matrix_1_errors.min()
        greatestdelta_2 = neuralnet.matrix_2_errors.max()
        smallestdelta_2 = neuralnet.matrix_2_errors.min()

        if any([greatestdelta_1 < 0.00001 and greatestdelta_1 > 0,
                smallestdelta_1 > -.00001 and smallestdelta_1 < 0]) or any(
            [greatestdelta_2 < 0.00001 and greatestdelta_2 > 0,
             smallestdelta_2 > -0.00001 and smallestdelta_2 < 0]):
            finished_working = True

    print(neuralnet.output_layer_output)

def test():
    test_sequences = open('data/rap1-lieb-test.txt')
    neuralnet = neural_network(68, 23, 1)
    neuralnet.matrix_1_bias = np.loadtxt('cnx_matrix_1.csv', delimiter=',')
    neuralnet.matrix_2_bias = np.loadtxt('cnx_matrix_2.csv', delimiter=',')

    neuralnet_outputs = open('neuralnet_predictions.txt', 'w')

    for test_seq in test_sequences:
        neuralnet.setin_n_exp_values(test_seq.strip())
        neuralnet.forward_propogation()
        neuralnet_outputs.write('{}\t{}\n'.format(test_seq.strip(), neuralnet.output_layer_output[0]))

    neuralnet_outputs.close()

if __name__ == 'train':
    import docopt
    import numpy as np
    from Bio import SeqIO
    import copy
    from .neuralnetfxns import neural_network

    args = docopt.docopt(__doc__)

    if args['autoencoder']:
        autoencoder()

    if args['testme']:
        testme()

    if args['test']:
        test()
