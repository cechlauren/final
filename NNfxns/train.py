"""
Training the neural net

What can be used:
     autoencoder
     %%%%%%%%werk <partitions> <sampling>
     test
Arguments:
    autoencoder
        Run ze autoencoder
    %%%%%%%%werk
        Do the rap1 learning task including cross-validation
    test
        Classify test data and output to tsv file type
    <partitions>
        Number of partitions to make for cross-valitation (k-fold!)
    <sampling>
        Sampling method for NN training input
        (slide) Go over each sequence in 17 nucleotide sliding frame (the length of the positives/test binding sites)
        (space) Each sequence is cut into 17 nucleotide bits for inputs (the binary bits that are useable by our model)
"""

def werk():
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
    
    positive_sites = [pos_seq.strip() for pos_seq in open('data/rap1-lieb-positives.txt')]

    
    
    # This part takes care of bringing in all those neg sequences
    
    negative_sites = list(SeqIO.parse('data/yeast-upstream-1k-negative.fa', 'fasta'))
    
    

    # Separate into random sections for the k-fold x-validation for both our positives and negatives 
    # Taken from : http://stackoverflow.com/questions/3352737/python-randomly-partition-a-list-into-n-nearly-equal-parts
    
    
    partitions = int(args['<partitions>'])
    neg_division = len(negative_sites) / float(partitions) # how many partitions we can make from the sites given
    neg_randomly_partitioned_list = [negative_sites[int(round(neg_division * i)): int(round(neg_division * (i + 1)))]
                                     for i in range(partitions)] # makes a list of those partitions thru each site

    pos_division = len(positive_sites) / float(partitions) # ditto ^^
    pos_randomly_partitioned_list = [positive_sites[int(round(pos_division * i)): int(round(pos_division * (i + 1)))]
                                     for i in range(partitions)]

    
    
    # Go thru neg sites subsets for x-validation, keep track of how many separations we do based on partitions
    
    separation = 0
    for index in range(int(args['<partitions>'])):
        # Set up cross-validation sets for the positives and negatives
        neg_site_list_copy = copy.deepcopy(neg_randomly_partitioned_list)
        del neg_site_list_copy[index]
        neg_site_training = [seq for partition in neg_site_list_copy for seq in partition]
        neg_cross_validation_set = neg_randomly_partitioned_list[index]

        pos_site_list_copy = copy.deepcopy(pos_randomly_partitioned_list)
        del pos_site_list_copy[index]
        pos_site_training = [seq for partition in pos_site_list_copy for seq in partition]
        pos_cross_validation_set = pos_randomly_partitioned_list[index]

        print("Training on the training set...")

        # Input our hyperparameters = # nodes
        NN = neural_network(68, 23, 1)

        # See neural_net.py to get info on initialization
        NN.initialize_values()

        pos_counter = 0
        counter = 0

        # If we're sampling from tneg then we'll slide over 17 nucleotides 
        
        if args['<sampling>'] == 'slide':
            for site in neg_site_training:

                # Iterate over site in 17 nucleotide sliding frames in negative sites, decide which model to use
                for block in range(len(site) - 16):
                    slice = site[block:(block + 17)].seq
                    if slice not in positive_sites:
                        if all([slice[4] == 'C', slice[5] == 'C', slice[9] == 'C']) == False:
                            NN.set_input_and_expected_values(slice, autoencoder=False, negative=True)
                            NN.forward_propogation()
                            NN.backward_propogation()
                            NN.update_weights_and_bias()
                            pos_counter += 1
                        else:
                            print(slice)

                    if pos_counter == len(pos_site_training):
                        for pos_site in pos_site_training:
                            NN.set_input_and_expected_values(pos_site, autoencoder=False, negative=False)
                            NN.forward_propogation()
                            NN.backward_propogation()
                            NN.update_weights_and_bias()

                        pos_counter = 0

                # have reset the positives counter and will now say that we've done some training on those
                
                counter += 1

                print("Training set: {}/{} completed...".format(counter, len(neg_cross_validation_set)))

                max_change_1 = NN.matrix_1_errors.max()
                min_change_1 = NN.matrix_1_errors.min()
                max_change_2 = NN.matrix_2_errors.max()
                min_change_2 = NN.matrix_2_errors.min()

                if any([max_change_1 < 0.00000000001 and max_change_1 > 0,
                        min_change_1 > -.00000000001 and min_change_1 < 0]) and any(
                    [max_change_2 < 0.00000000001 and max_change_2 > 0,
                     min_change_2 > -0.00000000001 and min_change_2 < 0]):
                    print("Stop criterion met after {} iterations".format(counter))
                    break

        #when we sample from the negatives we only take 17 nucleotide chunks from each site
        
        if args['<sampling>'] == 'space':
            for site in neg_site_training:
                
                number_of_blocks = int(len(site) / 17) #length of neg site tells us the amount of 17 length chunks possible

                for block in range(number_of_blocks):
                    slice = site[(block * 17):((block + 1) * 17)].seq
                    if slice not in positive_sites:
                        if all([slice[4] == 'C', slice[5] == 'C', slice[9] == 'C']) == False:
                            NN.set_input_and_expected_values(slice, autoencoder=False, negative=True)
                            NN.forward_propogation()
                            NN.backward_propogation()
                            NN.update_weights_and_bias()
                            pos_counter += 1

                        else:
                            print(slice)

                    #quick check to make sure that we've finished going thru the positives yet
                    if pos_counter == len(pos_site_training):
                        for pos_site in pos_site_training:
                            NN.set_input_and_expected_values(pos_site, autoencoder=False, negative=False)
                            NN.forward_propogation()
                            NN.backward_propogation()
                            NN.update_weights_and_bias()

                        pos_counter = 0

                    counter += 1

                max_change_1 = NN.matrix_1_errors.max()
                min_change_1 = NN.matrix_1_errors.min()
                max_change_2 = NN.matrix_2_errors.max()
                min_change_2 = NN.matrix_2_errors.min()

                if any([max_change_1 < 0.00000000001 and max_change_1 > 0,
                        min_change_1 > -.00000000001 and min_change_1 < 0]) and any(
                    [max_change_2 < 0.00000000001 and max_change_2 > 0,
                     min_change_2 > -0.00000000001 and min_change_2 < 0]):
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
                NN.set_input_and_expected_values(site[slice:slice + 17].seq, autoencoder=False, negative=True)
                NN.forward_propogation()
                neg_list.append(NN.output_layer_output)
            counter += 1
            print("Negative cross-validation: {}/{} completed...".format(counter, len(neg_cross_validation_set)))
            break

        print("Positive cross-validation set...")
        for site in pos_cross_validation_set:
            NN.set_input_and_expected_values(site, autoencoder=False)
            NN.forward_propogation()
            pos_list.append(NN.output_layer_output)

        print('Positive avg: {}'.format(sum(pos_list) / len(pos_list)))
        print('Negative avg: {}'.format(sum(neg_list) / len(neg_list)))
        print(NN.matrix_1_bias)
        print(NN.matrix_2_bias)

        # Output the connection matrices with greatest separation between the average positive and negative scores
        if ((sum(pos_list) / len(pos_list)) - (sum(neg_list) / len(neg_list))) > separation:
            np.savetxt('connection_matrix_1.csv', NN.matrix_1_bias, delimiter=',')
            np.savetxt('connection_matrix_2.csv', NN.matrix_2_bias, delimiter=',')
            separation = (sum(pos_list) / len(pos_list)) - (sum(neg_list) / len(neg_list))


# A simple definition of the autoencoder that uses those same NN parameters
def autoencoder():
  
    NN = neural_network()
    NN.set_input_and_expected_values('GA', autoencoder=True)
    NN.initialize_values()

    # Stop criterion
    finished_working = False

    while finished_working == False:
        NN.forward_propogation()
        NN.backward_propogation()
        NN.update_weights_and_bias()

        max_change_1 = NN.matrix_1_errors.max()
        min_change_1 = NN.matrix_1_errors.min()
        max_change_2 = NN.matrix_2_errors.max()
        min_change_2 = NN.matrix_2_errors.min()

        if any([max_change_1 < 0.00001 and max_change_1 > 0,
                min_change_1 > -.00001 and min_change_1 < 0]) or any(
            [max_change_2 < 0.00001 and max_change_2 > 0,
             min_change_2 > -0.00001 and min_change_2 < 0]):
            finished_working = True

    print(NN.output_layer_output)

def test():
    test_sequences = open('data/rap1-lieb-test.txt')
    NN = neural_network(68, 23, 1)
    NN.matrix_1_bias = np.loadtxt('connection_matrix_1.csv', delimiter=',')
    NN.matrix_2_bias = np.loadtxt('connection_matrix_2.csv', delimiter=',')

    NN_outputs = open('NN_predictions.txt', 'w')

    for test_seq in test_sequences:
        NN.set_input_and_expected_values(test_seq.strip())
        NN.forward_propogation()
        NN_outputs.write('{}\t{}\n'.format(test_seq.strip(), NN.output_layer_output[0]))

    NN_outputs.close()

if __name__ == 'train':
    import docopt
    import numpy as np
    from Bio import SeqIO
    import copy
    from .NNfxns import neural_network

    args = docopt.docopt(__doc__)

    if args['autoencoder']:
        autoencoder()

    if args['werk']:
        werk()

    if args['test']:
        test()
