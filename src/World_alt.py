import numpy as np
from Models_alt import *

# Global Variables
global elements, products, emergent_compound_products
elements = ["A","B","C","D"]
products = ["0", "♁", "♆", "🜺", "🜍", "♄", "☽", "☉"]
emergent_compound_products = {
	(): products[0],
	("A",): products[1],
	("B",): products[2],
	("C",): products[3],
	("D",): products[4],
	("A","B"): products[5],
	("A", "B", "D"): products[6],
	("B", "C", "D"): products[7]
}
# products = ["0", "♁", "♆", "🜺", "🜍", "♄", "☽", "☉"]
# emergent_compound_products = {
# 	(): products[0],
# 	("A",): products[1],
# 	("B",): products[2],
# 	("C",): products[3],
# 	("D",): products[4],
# 	("A","B"): products[5],
# 	("A", "B", "D"): None,
# 	("B", "C", "D"): None
# }

# Random Seed
np.random.seed(2023)

def main():
	ground_truth = Compound_Model(elements, products, emergent_compound_products)
	print(ground_truth.find_subsets([("A", "C", "D"), ("A","B","C")], common_subset = True))
	print(ground_truth.find_subsets([("A", "C", "D"), ("A","B","C")], common_subset = True, target = ground_truth._E_compounds))
	return
	partitions = ground_truth.find_subsets(ground_truth._compounds, target = ground_truth._E_compounds)
	all_products = ground_truth.get_products(ground_truth._compounds, target = "i")
	for index in range(len(ground_truth._compounds)):
		print(ground_truth._compounds[index])
		print(partitions[index])
		print(all_products[index])
		print()


	# ground_truth = Compound_Model(elements, emergent_compound_rewards)
	# bayes_model = Bayesian_Model(elements, emergent_compound_rewards)
	# print(bayes_model._uninitialized_E_compounds)
	# print(len(bayes_model._uninitialized_E_compounds))
	# bayes_model.add_emergent_compounds({("A","B"):7, ("C", "E"): -5})
	# print(bayes_model._uninitialized_E_compounds)
	# print(len(bayes_model._uninitialized_E_compounds))
	# for c in bayes_model.compounds: print(c, ":", bayes_model.get_rewards(c)[0])
	# bayes_model.del_emergent_compounds(None)
	# return

	# ground_truth = Compound_Model(elements, emergent_compound_rewards)
	# # The reward values can be accessed like this
	# print("Initial Rewards: ")
	# for c in ground_truth.compounds: print(c, ":", ground_truth.get_rewards(c)[0])

	# # The Compound_Model support dynamic assignment of emergent compounds. For
	# # example:
	# print("New Rewards: ")
	# ground_truth.add_emergent_compounds({("A", "E"): 9, ("A", "B", "E"): 13})
	# ground_truth.del_emergent_compounds([("A", "B"), ("A", "C", "D")])
	# # note that all the related predicted compound rewards are changed as well
	# for c in ground_truth.compounds: print(c, ":", ground_truth.get_rewards(c)[0])
	# # You can reinitialize the rewards
	# ground_truth.initialize_rewards(emergent_compound_rewards)

	# # Another core functionality is finding subset/superset information.
	# # For example:
	# subsets = ground_truth.find_subsets([("A", "C", "D"), ("B", "E")])
	# print("Subset of (A, C, D): ")
	# for compound in subsets[0]: print("\t", compound)
	# print("Subset of (B, E): ")
	# for compound in subsets[1]: print("\t", compound)
	# # You can also constrain the space to search for subsets. For example:
	# subset = ground_truth.find_subsets(("A", "C", "D"), target = ground_truth.E_compounds)
	# print("Subset of (A, C, D) among the Emergent Compounds: ")
	# for compound in subset[0]: print("\t", compound)
	# # Finally, you can also query subset relationships between several
	# # compounds:
	# subset_mat = ground_truth.is_subset([("A", "C", "D"), ("B", "E"), ("C",)], [("A", "B", "C", "D"), ("B", "D", "E")])
	# print("Subset relationship bewteen [(A, C, D), (B, E), (C,)] and [(A, B, C, D), (B, D, E)]: ")
	# print(subset_mat)
	# # These functions are mirrored for superset relations

	return

if __name__ == "__main__":
	main()