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

# Add random patterns to the list.
pattern_list.extend(factory.create_random_pattern_list(nr_patterns=3, on_probability=0.5))
plot_tools.plot_pattern_list(pattern_list)
# How similar are the random patterns and the checkerboard? Check the overlaps.
overlap_matrix = pattern_tools.compute_overlap_matrix(pattern_list)
plot_tools.plot_overlap_matrix(overlap_matrix)

# Let the hopfield network "learn" the patterns. Note: they are not stored
# explicitly but only network weights are updated !
hopfield_net.store_patterns(pattern_list)

