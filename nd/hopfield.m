from neurodynex.hopfield_network import network, pattern_tools, plot_tools

%{
Hopefield network model of associative memory
}%

pattern_size = 5

# Create an instance of the class HopfieldNetwork.
hopfield_net = network.HopfieldNetwork(nr_neurons= pattern_size**2)
# Instantiate a pattern factory.
factory = pattern_tools.PatternFactory(pattern_size, pattern_size)
# Create a checkerboard pattern and add it to the pattern list.
checkerboard = factory.create_checkerboard()
pattern_list = [checkerboard]
