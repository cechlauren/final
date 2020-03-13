"""
Testing the neural net/AE...more comments in the train.py file


"""
from final.NNfxns import neural_network
import numpy as np

#choosing random input to use for eval

def test_neural_network():
  
    NN = neural_network()
    NN.set_input_and_expected_values('AG', autoencoder=True)
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

