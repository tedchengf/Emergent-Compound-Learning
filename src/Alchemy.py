import numpy as np
from Models import *

global elements, products, emergent_compound_products, product_to_symbols, input_to_compound, compound_to_input
# DO NOT MODIFY THE ELEMENTS
elements = ["A","B","C","D"]
# You can modify the products to be whatever they want
products = ["0", "♁", "♆", "🜺", "🜍", "♄", "☽", "☉"]
# A dictionary of emergent compounds. You can also modify the composition here, 
# but be sure to keep the empty product
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
# A dictinoary mapping the input the subject get to choose and the model input.
# DO NOT MODIFY THE VALUES but feel free to modify the keys.
input_to_compound = {
	'': '',
	"A": "A",
	"E": "B",
	"F": "C",
	"W": "D"
}

def main():
	ground_truth = Compound_Model(elements, products, emergent_compound_products)
	past_observations = []
	model_history = []

	# a little while loop for demonstration
	print("-------------------------------------------------------------------")
	while True:
		curr_input = input("Compound to test: ")
		if curr_input == "Q" or curr_input == "q": exit()
		try:
			usr_rep, model_rep = Alchemy_Interface(ground_truth, input_to_compound, curr_input)
			past_observations.append(usr_rep)
			model_history.append(model_rep)
			print()
			print("Tested Compound:", usr_rep[0])
			print("Products:")
			for prod in usr_rep[1]: 
				print("  - " + prod)
			print("-------------------------------------------------------------------")
		except KeyError:
			print("Invalid compound. The compound must be made by using elements in set (A, E, F, W),")
	return

# Throws KeyError if the compound does not exist
# Note that by default the empty product is also present.
def Alchemy_Interface(compound_model, input_to_compound, input_str):
	compound = tuple(sorted([input_to_compound[c] for c in input_str]))
	curr_products = np.array(compound_model.get_products(compound)[0], dtype = int)
	curr_products_symbols = np.array([compound_model._id_to_prod[p] for p in curr_products])
	prod_idx = np.argsort(curr_products)
	return (tuple([c for c in input_str]), curr_products_symbols[prod_idx]), (compound, curr_products[prod_idx])
	
if __name__ == "__main__":
	main()