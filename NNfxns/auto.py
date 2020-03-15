import numpy as np

"""
-Make a neuralnet class
-define my functions
  -autoencoder and its parameters
  -initial values/input
  -vecotorization
  -expected values
  -output
  -backwards propagation
  -feedforward propagation
  -weights/bias
  -sigmoid function
  -sigmoid derivative

"""
# Make a class to allow self reference
class neural_network():

# The neuralnet will have a number of input neurons, a hidden layer, and number of an output neurons of the specified numbers, 
# which we apply to the class and can adjust as necessary later

    def __init__(self, input_neur=8, hidden_neur=3, output_neur=8):
        self.input_neur = input_neur
        self.hidden_neur = hidden_neur
        self.output_neur = output_neur

        # We want the autoencoder to learn the best weights for the data so we 
        # initialize matrices with random weights 0.1 > x > -0.1
        # These will link the autoencoder's different layers.
        
        self.cnx_matrix_1 = np.random.randn(self.input_neur, self.hidden_neur) / 10
        self.cnx_matrix_2 = np.random.randn(self.hidden_neur, self.output_neur) / 10

        # We have several vectors that our encoder will make to take in (input) and put out (output) data.
        # Recall the hidden layer will reduce the dimension, in a vector-wise fashion
        
        self.input_vector = None
        self.hidden_neur_output = None
        self.out_neur_output = None

        # We can also start off our autoencoder with a biased vector or matrices, depending on what we need from our model
        
        self.input_with_bias = None
        self.bias_mat1 = None
        self.bias_mat2 = None

        # Recall that we need to take the derivative of the cost function with respect to weight, and develop matrices from this
        
        self.hidden_dx_matrix = np.zeros([self.hidden_neur, self.hidden_neur])
        self.output_dx_matrix = np.zeros([self.output_neur, self.output_neur])

        # Hyperparameter: Learning Rate (controls how quickly a neuralnet learns a problem)
        # Typically in a range between [0.1, 1], configureable
        # A perfect learning rate will make the model learn to best approximate the function given available resources (layers/nodes)
        # in a given number of epochs
        # We'll start with a large one because we've got other stuff to do, dudes. The pitfall of that is that it will likely arrive at a suboptimal
        # set of weights.
        
        self.learning_rate = 2

        
        
        # Bit conversion for DNA into neuralnet inputs because I'm too lazy to think of a different way to encode these. Recall, we want either 0 or 1, not anything continuous.
        
        self.make_mulah = {'0': '0',
                                       '1': '1'
                                       }

        
        
        
        # These are the expected values for our neuralnet to return
        
        self.expected_values = None



    # Now we make moves to fill in and interpret those matrices
    
    # Initialize the values of the bias matrix given the information passing through the hidden layer and output layer
    
    def initialize_values(self):
        bias_ones_1 = np.ones([1, self.hidden_neur])
        self.bias_mat1 = np.append(self.cnx_matrix_1, bias_ones_1, axis=0)

        bias_ones_2 = np.ones([1, self.output_neur])
        self.bias_mat2 = np.append(self.cnx_matrix_2, bias_ones_2, axis=0)

    # We'll make the binding sites into binary bits that the encoder can interpret, and then tell it what to expect 
    
    def setin_n_exp_values(self, dnaIN, autoencoder=True, negative=True):
        # Convert input DNA sequence into binary bits
        self._construct_input_vector(dnaIN)
        # Set expected value depending on autoencoder or testme
        self._set_expected_values(autoencoder, negative)
        # Weight matrices and input/output vectors with bias applied
        self.input_with_bias = np.append(self.input_vector, [1])

    
    
    # This directly handles conversion of the DNA binding site string to the 1/0 vector describing it, and assigns it to the original input class
    
    def _construct_input_vector(self, dnaIN):
       
        temp_vector_list = []

        for base in dnaIN:
            for number in self.make_mulah[base]:
                temp_vector_list.append(float(number))

        self.input_vector = np.asarray(temp_vector_list)

    
    
    # This will set the values we expect from the neuralnet depending on whether or not we use the autoencoder (T/F).
    
    def _set_expected_values(self, autoencoder=True, negative=True):
        if autoencoder == True:
            self.expected_values = self.input_vector

        if autoencoder == False:
            if negative == True:
                self.expected_values = 0
            if negative == False:
                self.expected_values = 1

    # Recall an autoencoder is a feedforward method. We need to convert input to hidden layer reduced dim info to output layer
    # that is the same as the input.
    
    def forward_propogation(self):
        # Generates hidden layer outputs
        
        output_one_list = []

        for element in np.nditer(np.dot(self.input_with_bias, self.bias_mat1)):
            output_one_list.append(self._sigmoid_function(element))
        self.hidden_neur_output = np.asarray(output_one_list)

        # Calculate the square derivate matrix for the hidden layer outputs
        
        for position, value in enumerate(self.hidden_neur_output):
            self.hidden_dx_matrix[position][position] = self._sigmoid_function_derivative(value)

        # The results from the output layer 
        # Add bias to hidden_neur_output
        
        self.hidden_output_bias = np.append(self.hidden_neur_output, [1])

        output_two_list = []
        for element in np.nditer(np.dot(self.hidden_output_bias, self.bias_mat2)):
            output_two_list.append(self._sigmoid_function(element))
        self.out_neur_output = np.asarray(output_two_list)

        # Calculate square derivate matrix for output layer outputs
        
        for position, value in enumerate(self.out_neur_output):
            self.output_dx_matrix[position][position] = self._sigmoid_function_derivative(value)

    
    # Recall that in stochastic gradient descent, we estimate the error gradient for the current state of the model using
    # examples from the training. The backwards propagation is what we use to update the weights of the model. 
    
    def backward_propogation(self):
        # Output Layer error
        
        deviations = self.out_neur_output - self.expected_values

        out_neur_errors = np.dot(self.output_dx_matrix, deviations)

        # Hidden Layer error
        
        hidden_neur_errors = np.dot(np.dot(self.hidden_dx_matrix, self.cnx_matrix_2), out_neur_errors)

        # Matrix 2 Errors (those of our output layer)
        
        output_rated_row_vector = np.asmatrix(-(self.learning_rate * out_neur_errors)).transpose()
        errors_mat2_transposed = np.dot(output_rated_row_vector, np.asmatrix(self.hidden_output_bias))
        self.errors_mat2 = errors_mat2_transposed.transpose()

        # Matrix 1 Errors (those of our hidden layer)
        
        hidden_rated_row_vector = np.asmatrix(-(self.learning_rate * hidden_neur_errors)).transpose()
        errors_mat1_transposed = np.dot(hidden_rated_row_vector, np.asmatrix(self.input_with_bias))
        self.errors_mat1 = errors_mat1_transposed.transpose()

    
    # Basically the sum between the output and hidden layers' bias and errors is what we use...
    
    def weight_bias_renew(self):
        self.bias_mat1 = self.bias_mat1 + self.errors_mat1
        self.bias_mat2 = self.bias_mat2 + self.errors_mat2

        self.cnx_matrix_1 = self.bias_mat1[:-1]
        self.cnx_matrix_2 = self.bias_mat2[:-1]

    
    # The activation function for the layers. Some people use relu. Its more of less simple to take the deriv of sigmoid. 
    
    def _sigmoid_function(self, input):
        return float(1 / (1 + np.exp(-input)))

    
    # Recall we also take the derivative of the activatio function. 
    def _sigmoid_function_derivative(self, input):
        return float(self._sigmoid_function(input) * (1 - self._sigmoid_function(input)))
