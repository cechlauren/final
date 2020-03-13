"""
Testing the AE on the 8bit info


"""
from final.NNfxns import auto
import numpy as np
import sys



def test_auto():
  
    eightBit = open('final/data/eightBit.txt')
    NN = auto()
    for bits in eightBit:
        NN.set_input_and_expected_values(bits.strip(), autoencoder=True)
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
    
